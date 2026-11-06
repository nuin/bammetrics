//! A gzip file read as `java.util.zip.GZIPInputStream` reads it: one member
//! after another (so BGZF works too, each block being a member), every
//! member's CRC and length checked, and bytes after the last member that do
//! not start a new member ignored. `noodles-bgzf` needs the BGZF extra
//! field, so a plain `gzip` file goes through this reader; it inflates with
//! the `zlib-rs` that `noodles-bgzf` already builds on.

use std::io::{self, BufRead, Read};

use zlib_rs::{Inflate, InflateFlush, Status};

/// RFC 1952 member magic; the first byte starts every member.
const MAGIC: [u8; 2] = [0x1f, 0x8b];

/// `zlib` window bits for a gzip header and a 32 KiB window.
const GZIP_WINDOW_BITS: u8 = 16 + 15;

pub struct GzipReader<R> {
    inner: R,
    /// The member being inflated; `None` between members.
    member: Option<Inflate>,
    members: usize,
    done: bool,
}

impl<R: BufRead> GzipReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            member: None,
            members: 0,
            done: false,
        }
    }

    /// Start the next member, or end the stream: at the end of the input,
    /// or, after the first member, at a byte that cannot start a gzip
    /// header. The inflater checks the rest of the header.
    fn start_member(&mut self) -> io::Result<bool> {
        let first = self.inner.fill_buf()?.first().copied();
        if first != Some(MAGIC[0]) {
            if self.members == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Not in GZIP format",
                ));
            }
            return Ok(false);
        }
        self.member = Some(Inflate::new(true, GZIP_WINDOW_BITS));
        self.members += 1;
        Ok(true)
    }
}

fn inflate_error(error: zlib_rs::InflateError) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, error.as_str())
}

impl<R: BufRead> Read for GzipReader<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        while !self.done && !out.is_empty() {
            if self.member.is_none() && !self.start_member()? {
                self.done = true;
                break;
            }
            let member = self.member.as_mut().expect("started above");
            let input = self.inner.fill_buf()?;
            let at_end = input.is_empty();
            let (in_before, out_before) = (member.total_in(), member.total_out());
            let status = member
                .decompress(input, out, InflateFlush::NoFlush)
                .map_err(inflate_error)?;
            let consumed = usize::try_from(member.total_in() - in_before).expect("fits");
            let produced = usize::try_from(member.total_out() - out_before).expect("fits");
            self.inner.consume(consumed);
            if status == Status::StreamEnd {
                self.member = None;
            } else if at_end && produced == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "Unexpected end of ZLIB input stream",
                ));
            }
            if produced > 0 {
                return Ok(produced);
            }
        }
        Ok(0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `printf 'hello\n' | gzip -n` and `printf 'world\n' | gzip -n`.
    const HELLO: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xcb, 0x48, 0xcd, 0xc9, 0xc9,
        0xe7, 0x02, 0x00, 0x20, 0x30, 0x3a, 0x36, 0x06, 0x00, 0x00, 0x00,
    ];
    const WORLD: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0x2b, 0xcf, 0x2f, 0xca, 0x49,
        0xe1, 0x02, 0x00, 0xa8, 0x61, 0x38, 0xdd, 0x06, 0x00, 0x00, 0x00,
    ];

    fn inflate(bytes: &[u8], capacity: usize) -> io::Result<String> {
        let reader = io::BufReader::with_capacity(capacity, bytes);
        let mut text = String::new();
        GzipReader::new(reader).read_to_string(&mut text)?;
        Ok(text)
    }

    #[test]
    fn members_are_read_in_turn() {
        let both = [HELLO, WORLD].concat();
        for capacity in [1, 2, 3, 7, 64, 8192] {
            assert_eq!(
                inflate(&both, capacity).expect("inflates"),
                "hello\nworld\n"
            );
        }
    }

    #[test]
    fn bytes_after_the_last_member_are_ignored() {
        let trailing = [HELLO, b"junk".as_slice()].concat();
        assert_eq!(inflate(&trailing, 8192).expect("inflates"), "hello\n");
    }

    #[test]
    fn a_file_that_is_not_gzip_or_is_cut_short_fails() {
        assert!(inflate(b">chr1\nACGT\n", 8192).is_err());
        assert!(inflate(b"", 8192).is_err());
        assert!(inflate(&HELLO[..HELLO.len() - 3], 8192).is_err());
        let mut corrupt = HELLO.to_vec();
        corrupt[20] ^= 0xff;
        assert!(inflate(&corrupt, 8192).is_err());
    }
}

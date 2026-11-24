//! The BGZF blocks of a BAM handed on one whole block at a time
//! ([`BlockEnds`]), so a file cut inside a block header fails and bytes
//! after the end-of-file block are not read.

use std::io::{self, Read};

/// A BAM whose last bytes are part of a BGZF block header: the error
/// [`BlockEnds`] gives, which `pass::read_error` calls a truncated BAM.
#[derive(Debug)]
pub(crate) struct CutBlockHeader;

impl std::fmt::Display for CutBlockHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("truncated BAM: the file ends inside a BGZF block header")
    }
}

impl std::error::Error for CutBlockHeader {}

/// Bytes of a BGZF block header; `BSIZE` (the block size less 1) is its
/// last two.
const BGZF_HEADER_BYTES: usize = 18;

/// Bytes of the gzip trailer: CRC-32, then `ISIZE`, the inflated size.
const TRAILER_BYTES: usize = 8;

/// The bytes of a BAM as they come, one BGZF block at a time (at most
/// 65,536 bytes held, two after an empty block).
///
/// noodles-bgzf 0.52 reads a header with `read_exact` and takes an
/// `UnexpectedEof` there as the end of the file, so a file cut 1 to 17
/// bytes into a header would read as whole. Here the end of the file
/// inside a header is [`CutBlockHeader`], of kind `InvalidData`, which
/// noodles passes on. A file that ends on a block boundary (no end-of-file
/// block included) reads as before; one that ends inside a block's data
/// fails in noodles.
///
/// htsjdk's `BlockCompressedInputStream` reads an empty block (`ISIZE` 0)
/// as the end of the stream, so Picard never reads the bytes after the
/// end-of-file block. Here the block after an empty block is read ahead:
/// when it is a whole BGZF block the read goes on (an empty block in the
/// middle of a file, where htsjdk stops; KNOWN_DIFFERENCES.md), and when
/// the file ends or the bytes are anything else (zeros, a cut or foreign
/// gzip member, other bytes) the empty block ends the stream, as in
/// htsjdk.
pub(crate) struct BlockEnds<R> {
    inner: R,
    /// The block being handed on.
    block: Vec<u8>,
    /// Bytes of `block` handed on.
    at: usize,
    /// The whole block read after an empty one, handed on next.
    ahead: Option<Vec<u8>>,
    /// What follows `block`.
    then: Then,
}

/// What [`BlockEnds`] does once its block is handed on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Then {
    /// Read the next block.
    Block,
    /// Hand on the rest of the file as it is: bytes that start no BGZF
    /// block, which noodles refuses.
    Raw,
    /// Nothing more: the file ended, or an empty block ended the stream.
    End,
}

/// One read of a block from the file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Got {
    /// No byte: the file ended on a block boundary.
    End,
    /// A whole block; `empty` when its `ISIZE` is 0.
    Whole { empty: bool },
    /// The file ended inside the block's header.
    CutHeader,
    /// The file ended inside the block's data.
    CutData,
    /// The bytes start no BGZF block as noodles reads one.
    Foreign,
}

/// The header noodles-bgzf 0.52 accepts (`is_valid_header`).
fn is_bgzf_header(header: &[u8]) -> bool {
    header.len() == BGZF_HEADER_BYTES
        && header[..4] == [0x1f, 0x8b, 0x08, 0x04]
        && header[10..16] == [0x06, 0x00, b'B', b'C', 0x02, 0x00]
}

/// Read into `buf` until it holds `want` bytes or the file ends.
fn fill<R: Read>(inner: &mut R, buf: &mut Vec<u8>, want: usize) -> io::Result<()> {
    let have = buf.len();
    (&mut *inner)
        .take((want - have) as u64)
        .read_to_end(buf)
        .map(|_| ())
}

/// Read one block of `inner` into `buf`.
fn read_block<R: Read>(inner: &mut R, buf: &mut Vec<u8>) -> io::Result<Got> {
    buf.clear();
    fill(inner, buf, BGZF_HEADER_BYTES)?;
    if buf.is_empty() {
        return Ok(Got::End);
    }
    if buf.len() < BGZF_HEADER_BYTES {
        return Ok(Got::CutHeader);
    }
    let size = usize::from(u16::from_le_bytes([buf[16], buf[17]])) + 1;
    if !is_bgzf_header(buf) || size < BGZF_HEADER_BYTES + TRAILER_BYTES {
        return Ok(Got::Foreign);
    }
    fill(inner, buf, size)?;
    if buf.len() < size {
        return Ok(Got::CutData);
    }
    Ok(Got::Whole {
        empty: buf[size - 4..] == [0, 0, 0, 0],
    })
}

impl<R: Read> BlockEnds<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            block: Vec::new(),
            at: 0,
            ahead: None,
            then: Then::Block,
        }
    }

    /// Read the next block into `block` and decide what follows it.
    fn next_block(&mut self) -> io::Result<()> {
        self.at = 0;
        let got = match self.ahead.take() {
            Some(block) => {
                self.block = block;
                Got::Whole {
                    empty: self.block.ends_with(&[0, 0, 0, 0]),
                }
            }
            None => read_block(&mut self.inner, &mut self.block)?,
        };
        match got {
            Got::End => self.then = Then::End,
            Got::CutHeader => {
                self.block.clear();
                return Err(io::Error::new(io::ErrorKind::InvalidData, CutBlockHeader));
            }
            // noodles' `read_exact` fails inside the data.
            Got::CutData => self.then = Then::End,
            Got::Foreign => self.then = Then::Raw,
            Got::Whole { empty: false } => {}
            Got::Whole { empty: true } => {
                let mut next = Vec::new();
                match read_block(&mut self.inner, &mut next)? {
                    Got::Whole { .. } => self.ahead = Some(next),
                    _ => self.then = Then::End,
                }
            }
        }
        Ok(())
    }
}

impl<R: Read> Read for BlockEnds<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if buf.is_empty() {
            return Ok(0);
        }
        while self.at == self.block.len() {
            match self.then {
                Then::End => return Ok(0),
                Then::Raw => return self.inner.read(buf),
                Then::Block => self.next_block()?,
            }
        }
        let n = buf.len().min(self.block.len() - self.at);
        buf[..n].copy_from_slice(&self.block[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    /// The 28-byte BGZF end-of-file block.
    const EOF_BLOCK: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0x1b, 0, 0x03,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    /// One BGZF block holding `text`.
    fn data_block(text: &[u8]) -> Vec<u8> {
        let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
        writer.write_all(text).expect("bgzf");
        let mut bytes = writer.finish().expect("bgzf");
        bytes.truncate(bytes.len() - EOF_BLOCK.len());
        bytes
    }

    /// Every byte `bytes` hands on, read `step` bytes at a time.
    fn read_all(bytes: &[u8], step: usize) -> io::Result<Vec<u8>> {
        let mut reader = BlockEnds::new(bytes);
        let mut buf = vec![0; step];
        let mut out = Vec::new();
        loop {
            match reader.read(&mut buf)? {
                0 => return Ok(out),
                n => out.extend_from_slice(&buf[..n]),
            }
        }
    }

    fn join(parts: &[&[u8]]) -> Vec<u8> {
        parts.concat()
    }

    #[test]
    fn a_cut_header_fails_and_a_block_boundary_does_not() {
        let data = data_block(b"some text");
        let whole = join(&[&data, &data, &EOF_BLOCK]);
        for step in [1, 5, 18, 28, 4096] {
            assert_eq!(read_all(&whole, step).expect("whole"), whole);
            assert_eq!(read_all(&data, step).expect("boundary"), data);
            assert_eq!(read_all(&[], step).expect("empty"), b"");
            // Inside a block's data: noodles finds that one itself.
            let cut = &whole[..data.len() + 20];
            assert_eq!(read_all(cut, step).expect("data"), cut);
            for at in [1, 10, 17] {
                let cut = &whole[..data.len() + at];
                let err = read_all(cut, step).expect_err("cut header");
                assert_eq!(err.kind(), io::ErrorKind::InvalidData);
                assert!(err.to_string().contains("truncated BAM"), "{err}");
            }
        }
    }

    #[test]
    fn bytes_after_the_end_of_file_block_are_not_read() {
        let data = data_block(b"records");
        let file = join(&[&data, &EOF_BLOCK]);
        let mut gzip = plain_gzip_member();
        gzip.truncate(gzip.len().min(40));
        let tails: Vec<Vec<u8>> = vec![
            vec![0],
            vec![0; 17],
            vec![0; 18],
            vec![0; 1000],
            vec![0xff; 10],
            EOF_BLOCK[..4].to_vec(),
            EOF_BLOCK[..18].to_vec(),
            data[..data.len() - 1].to_vec(),
            gzip,
        ];
        for tail in tails {
            let bytes = join(&[&file, &tail]);
            for step in [1, 7, 4096] {
                assert_eq!(read_all(&bytes, step).expect("tail"), file, "{tail:?}");
            }
        }
    }

    #[test]
    fn a_whole_block_after_an_empty_one_is_read() {
        let data = data_block(b"more");
        let twice = join(&[&EOF_BLOCK, &EOF_BLOCK]);
        let middle = join(&[&data, &EOF_BLOCK, &data, &EOF_BLOCK]);
        for step in [1, 28, 4096] {
            assert_eq!(read_all(&twice, step).expect("two"), twice);
            assert_eq!(read_all(&middle, step).expect("middle"), middle);
        }
    }

    #[test]
    fn foreign_bytes_after_a_data_block_are_handed_on() {
        let data = data_block(b"x");
        let bytes = join(&[&data, &[0; 30]]);
        assert_eq!(read_all(&bytes, 4096).expect("raw"), bytes);
    }

    /// A gzip member that is no BGZF block: no `FEXTRA` field.
    fn plain_gzip_member() -> Vec<u8> {
        let mut member = vec![0x1f, 0x8b, 0x08, 0, 0, 0, 0, 0, 0, 0xff];
        member.extend_from_slice(&[0x03, 0x00]);
        member.extend_from_slice(&[0; 8]);
        member
    }
}

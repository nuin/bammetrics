//! The BGZF block boundaries of a BAM followed as its bytes are read, so a
//! file cut inside a block header fails ([`BlockEnds`]).

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

/// The bytes of a BAM as they come, with the BGZF block boundaries
/// followed through each header's `BSIZE`. noodles-bgzf 0.52 reads a
/// header with `read_exact` and takes an `UnexpectedEof` there as the end
/// of the file, so a file cut 1 to 17 bytes into a header reads as whole.
/// Here the end of the file inside a header is [`CutBlockHeader`], of kind
/// `InvalidData`, which noodles passes on. A file that ends on a block
/// boundary (no end-of-file block included) reads as before; one that ends
/// inside a block's data fails in noodles, as before.
pub(crate) struct BlockEnds<R> {
    inner: R,
    header: [u8; BGZF_HEADER_BYTES],
    /// Bytes of the current header seen.
    have: usize,
    /// Bytes of the current block after its header still to come.
    rest: u64,
}

impl<R> BlockEnds<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            header: [0; BGZF_HEADER_BYTES],
            have: 0,
            rest: 0,
        }
    }

    fn follow(&mut self, mut bytes: &[u8]) {
        while !bytes.is_empty() {
            if self.rest > 0 {
                let take = bytes
                    .len()
                    .min(usize::try_from(self.rest).unwrap_or(usize::MAX));
                self.rest -= take as u64;
                bytes = &bytes[take..];
                continue;
            }
            let take = bytes.len().min(BGZF_HEADER_BYTES - self.have);
            self.header[self.have..self.have + take].copy_from_slice(&bytes[..take]);
            self.have += take;
            bytes = &bytes[take..];
            if self.have == BGZF_HEADER_BYTES {
                let bsize = u16::from_le_bytes([self.header[16], self.header[17]]);
                // noodles refuses a block shorter than its header and
                // trailer, so a short `BSIZE` never reaches the next block.
                self.rest = (u64::from(bsize) + 1).saturating_sub(BGZF_HEADER_BYTES as u64);
                self.have = 0;
            }
        }
    }
}

impl<R: Read> Read for BlockEnds<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n == 0 && !buf.is_empty() && self.have > 0 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, CutBlockHeader));
        }
        self.follow(&buf[..n]);
        Ok(n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The 28-byte BGZF end-of-file block.
    const EOF_BLOCK: [u8; 28] = [
        0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, b'B', b'C', 0x02, 0, 0x1b, 0, 0x03,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
    ];

    /// Every byte of `bytes` read `step` bytes at a time.
    fn read_all(bytes: &[u8], step: usize) -> io::Result<usize> {
        let mut reader = BlockEnds::new(bytes);
        let mut buf = vec![0; step];
        let mut total = 0;
        loop {
            match reader.read(&mut buf)? {
                0 => return Ok(total),
                n => total += n,
            }
        }
    }

    #[test]
    fn a_cut_header_fails_and_a_block_boundary_does_not() {
        let two: Vec<u8> = EOF_BLOCK.iter().chain(&EOF_BLOCK).copied().collect();
        for step in [1, 5, 18, 28, 4096] {
            assert_eq!(read_all(&two, step).expect("whole"), 56);
            assert_eq!(read_all(&two[..28], step).expect("boundary"), 28);
            assert_eq!(read_all(&[], step).expect("empty"), 0);
            // Inside a block's data: noodles finds that one itself.
            assert_eq!(read_all(&two[..50], step).expect("data"), 50);
            for cut in [29, 36, 45] {
                let err = read_all(&two[..cut], step).expect_err("cut header");
                assert_eq!(err.kind(), io::ErrorKind::InvalidData);
                assert!(err.to_string().contains("truncated BAM"), "{err}");
            }
        }
    }
}

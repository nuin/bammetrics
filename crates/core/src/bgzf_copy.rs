//! BGZF blocks copied as bytes: htsjdk 5.0.0 `BamFileIoUtils.blockCopyBamFile`,
//! the copy under `GatherBamFiles` (`gatherWithBlockCopying`) and under
//! `AddCommentsToBam` and `ReplaceSamHeader` on a BAM (`reheaderBamFile`).
//!
//! htsjdk copies a BAM after its header in three steps, and [`Source::copy`]
//! takes the same three:
//!
//! 1. `findVirtualOffsetOfFirstRecordInBam`: the header is read (magic,
//!    text, reference list) and the virtual offset after it is taken as
//!    `BlockCompressedInputStream.getFilePointer` gives it: the next block
//!    at offset 0 when the header ends at the end of a block.
//! 2. A new stream seeks there, and the rest of the block it stands in is
//!    inflated and deflated again through a `BlockCompressedOutputStream`
//!    at `COMPRESSION_LEVEL`. At offset 0 that is the whole first block of
//!    records; an empty block there is passed over to the next one.
//! 3. The blocks after that one are copied as bytes up to the end of the
//!    file, less the 28-byte end-of-file block when the caller skips it.
//!
//! Without the header skip (the first file of a gather), step 3 starts at
//! byte 0.
//!
//! Step 2 deflates with zlib-rs, which writes other bytes than the JDK's
//! zlib at levels 1 to 9 (docs/concordance.md, "BAM bytes"). Such a block
//! inflates to Picard's bytes; at level 0 it is Picard's block, a stored
//! deflate block. Every block of step 3 is the input's.
//!
//! Every length read from a file is checked. A block header must be a BGZF
//! header (gzip magic with `FEXTRA`, `XLEN` 6 and the `BC` subfield) whose
//! `BSIZE` gives 26 to 65,536 bytes; an inflated block must hold `ISIZE`
//! bytes, at most 65,536. The header text is read through and not held;
//! each reference name is read through [`read_append`], so no buffer is
//! sized from a length before its bytes are read, and the reference list
//! and the text share one budget of [`MAX_HEADER_BYTES`]. The blocks of
//! step 3 are walked header by header: a file cut inside a block, or a byte
//! that starts no block, fails the copy where htsjdk copies it.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use zlib_rs::{DeflateConfig, InflateConfig, Method, ReturnCode, Strategy};

use crate::bam_raw::read_append;
use crate::error::{Error, Result};

/// The BGZF end-of-file block (SAM specification 4.1.2).
pub const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// htsjdk's `GZIP_BLOCK_PREAMBLE`: the first 16 bytes of every block it
/// writes (no modification time, `XFL` 0, `OS` 255).
const PREAMBLE: [u8; 16] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
];

const HEADER_BYTES: usize = 18;
const FOOTER_BYTES: usize = 8;
/// The largest block, compressed or inflated.
const MAX_BLOCK_BYTES: usize = 1 << 16;
/// htsjdk's `DEFAULT_UNCOMPRESSED_BLOCK_SIZE`: the bytes a
/// `BlockCompressedOutputStream` puts in one block.
pub const BLOCK_DATA_BYTES: usize = 65_498;
/// The header text and the reference list of one BAM together: a Java
/// `String` holds 2^31 characters, so htsjdk reads no larger text.
pub const MAX_HEADER_BYTES: u64 = 1 << 31;
/// The bytes of one copy buffer.
const COPY_BUFFER: usize = 1 << 20;

/// `BlockCompressedInputStream.checkTermination`: how a file ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Termination {
    /// The 28-byte end-of-file block.
    Terminator,
    /// A block that ends at the last byte, other than the end-of-file block.
    HealthyLastBlock,
    /// Anything else: fewer than 28 bytes, or a last block cut short.
    Defective,
}

/// How the file `file` of `size` bytes ends, read from its last 64 KiB as
/// htsjdk reads it: the end-of-file block, else the last block header found
/// scanning back from 28 bytes before the end, which must end at the end.
pub fn termination(file: &mut File, size: u64) -> io::Result<Termination> {
    if size < EOF_BLOCK.len() as u64 {
        return Ok(Termination::Defective);
    }
    let tail = size.min(MAX_BLOCK_BYTES as u64);
    file.seek(SeekFrom::Start(size - tail))?;
    let mut buf = Vec::new();
    read_append(file, &mut buf, tail as usize)?;
    if buf.ends_with(&EOF_BLOCK) {
        return Ok(Termination::Terminator);
    }
    for at in (0..=buf.len() - EOF_BLOCK.len()).rev() {
        if buf[at..].starts_with(&PREAMBLE) {
            let bsize = usize::from(u16::from_le_bytes([buf[at + 16], buf[at + 17]]));
            return Ok(if buf.len() - at == bsize + 1 {
                Termination::HealthyLastBlock
            } else {
                Termination::Defective
            });
        }
    }
    Ok(Termination::Defective)
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// The size of the block whose first 18 bytes are `header`, or why they
/// start no BGZF block.
fn block_size(header: &[u8]) -> io::Result<usize> {
    if header[..4] != [0x1f, 0x8b, 0x08, 0x04] {
        return Err(invalid(
            "no BGZF block starts here (gzip magic with FEXTRA)",
        ));
    }
    if header[10..16] != [6, 0, b'B', b'C', 2, 0] {
        return Err(invalid("a gzip block without the BGZF BC subfield"));
    }
    let size = usize::from(u16::from_le_bytes([header[16], header[17]])) + 1;
    if size < HEADER_BYTES + FOOTER_BYTES {
        return Err(invalid(format!("a BGZF block of {size} bytes")));
    }
    Ok(size)
}

/// The next whole block of `reader` into `raw`; false at the end of the
/// file, before any byte of a header.
fn read_raw_block<R: Read>(reader: &mut R, raw: &mut Vec<u8>) -> io::Result<bool> {
    raw.clear();
    raw.resize(HEADER_BYTES, 0);
    let mut filled = 0;
    while filled < HEADER_BYTES {
        match reader.read(&mut raw[filled..]) {
            Ok(0) => break,
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    if filled == 0 {
        raw.clear();
        return Ok(false);
    }
    if filled < HEADER_BYTES {
        return Err(invalid("the file ends inside a BGZF block header"));
    }
    let size = block_size(raw)?;
    raw.resize(size, 0);
    reader
        .read_exact(&mut raw[HEADER_BYTES..])
        .map_err(|e| match e.kind() {
            io::ErrorKind::UnexpectedEof => invalid("the file ends inside a BGZF block"),
            _ => e,
        })?;
    Ok(true)
}

/// The inflated bytes of the block `raw` into `out`.
fn inflate(raw: &[u8], out: &mut Vec<u8>) -> io::Result<()> {
    let footer = &raw[raw.len() - 4..];
    let expected = u32::from_le_bytes([footer[0], footer[1], footer[2], footer[3]]) as usize;
    if expected > MAX_BLOCK_BYTES {
        return Err(invalid(format!(
            "a BGZF block claims {expected} inflated bytes"
        )));
    }
    // One byte more than expected, so a block that inflates to more fails.
    out.clear();
    out.resize(expected + 1, 0);
    let payload = &raw[HEADER_BYTES..raw.len() - FOOTER_BYTES];
    let (done, code) = zlib_rs::decompress_slice(out, payload, InflateConfig { window_bits: -15 });
    let got = done.len();
    if code != ReturnCode::Ok || got != expected {
        return Err(invalid(format!(
            "a BGZF block does not inflate to its {expected} bytes ({code:?})"
        )));
    }
    out.truncate(expected);
    Ok(())
}

/// One inflated block and where it stands.
struct Block {
    address: u64,
    compressed: usize,
    data: Vec<u8>,
}

/// htsjdk's `BlockCompressedInputStream` over a file, as far as the header
/// read and the seek after it use it: a block at a time, an empty block of
/// compressed size 0 at the end of the file.
struct Inflating<R> {
    reader: R,
    /// The address of the next block.
    next: u64,
    raw: Vec<u8>,
    current: Option<Block>,
    offset: usize,
    size: u64,
}

impl<R: Read + Seek> Inflating<R> {
    fn new(reader: R, size: u64) -> Self {
        Self {
            reader,
            next: 0,
            raw: Vec::new(),
            current: None,
            offset: 0,
            size,
        }
    }

    /// `readBlock`: the next block becomes the current one.
    fn read_block(&mut self) -> io::Result<()> {
        let address = self.next;
        let mut data = self.current.take().map(|b| b.data).unwrap_or_default();
        let compressed = if read_raw_block(&mut self.reader, &mut self.raw)? {
            inflate(&self.raw, &mut data)?;
            self.raw.len()
        } else {
            data.clear();
            0
        };
        self.next = address + compressed as u64;
        self.current = Some(Block {
            address,
            compressed,
            data,
        });
        self.offset = 0;
        Ok(())
    }

    /// `available`: the bytes left in the current block, after reading the
    /// next block when the current one is used up.
    fn available(&mut self) -> io::Result<usize> {
        if self
            .current
            .as_ref()
            .is_none_or(|block| self.offset == block.data.len())
        {
            self.read_block()?;
        }
        Ok(self.current.as_ref().map_or(0, |b| b.data.len()) - self.offset)
    }

    /// `read` into `out` until it is full: an empty block reads as the end
    /// of the stream, as in htsjdk.
    fn read_exact(&mut self, mut out: &mut [u8]) -> io::Result<()> {
        while !out.is_empty() {
            let available = self.available()?;
            if available == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "the BAM header ends before its last byte (Premature EOF)",
                ));
            }
            let n = available.min(out.len());
            let block = self.current.as_ref().expect("read above");
            out[..n].copy_from_slice(&block.data[self.offset..self.offset + n]);
            self.offset += n;
            out = &mut out[n..];
        }
        Ok(())
    }

    /// Read `count` bytes and drop them.
    fn skip(&mut self, mut count: u64) -> io::Result<()> {
        let mut scratch = [0u8; 4096];
        while count > 0 {
            let n = count.min(scratch.len() as u64) as usize;
            self.read_exact(&mut scratch[..n])?;
            count -= n as u64;
        }
        Ok(())
    }

    fn read_i32(&mut self) -> io::Result<i32> {
        let mut bytes = [0u8; 4];
        self.read_exact(&mut bytes)?;
        Ok(i32::from_le_bytes(bytes))
    }

    /// `getFilePointer`: the block and offset, the next block at 0 once the
    /// current one is read to its end.
    fn file_pointer(&self) -> (u64, usize) {
        match &self.current {
            None => (0, 0),
            Some(block) if self.offset > 0 && self.offset == block.data.len() => {
                (block.address + block.compressed as u64, 0)
            }
            Some(block) => (block.address, self.offset),
        }
    }

    /// `seek(address, offset)`: the block at `address` read, then, as
    /// `available` does, the next one when it is empty.
    fn seek(&mut self, address: u64, offset: usize) -> io::Result<()> {
        self.reader.seek(SeekFrom::Start(address))?;
        self.next = address;
        self.current = None;
        self.read_block()?;
        let available = self.available()?;
        let at_end = self.current.as_ref().is_some_and(|b| {
            let after = b.address + b.compressed as u64;
            after >= self.size || self.size - after == EOF_BLOCK.len() as u64
        });
        if offset > available || (offset == available && !at_end) {
            return Err(invalid(format!(
                "Invalid file pointer: {}",
                (address << 16) | offset as u64
            )));
        }
        self.offset = offset;
        Ok(())
    }
}

/// Where a BAM's records start, and its binary reference list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BamStart {
    /// The virtual offset of the first record: block address and offset.
    pub first_record: (u64, usize),
    /// Each reference's name (without its NUL) and length.
    pub references: Vec<(Vec<u8>, i32)>,
    /// A byte follows the header before an empty block or the end: the
    /// file holds records.
    pub has_records: bool,
}

/// `BAMFileReader.readHeader` on the inflated stream: the magic, the text
/// (read through), the reference list.
fn read_start<R: Read + Seek>(stream: &mut Inflating<R>) -> io::Result<BamStart> {
    let mut magic = [0u8; 4];
    stream.read_exact(&mut magic)?;
    if magic != *b"BAM\x01" {
        return Err(invalid("Invalid BAM file header"));
    }
    let l_text = stream.read_i32()?;
    let l_text = u64::try_from(l_text)
        .map_err(|_| invalid(format!("the BAM header text claims {l_text} bytes")))?;
    let mut budget = MAX_HEADER_BYTES.saturating_sub(l_text);
    stream.skip(l_text)?;
    let count = stream.read_i32()?;
    let count = u64::try_from(count)
        .map_err(|_| invalid(format!("the BAM reference list claims {count} entries")))?;
    // Each entry is at least its two lengths and a NUL.
    if count > budget / 9 {
        return Err(invalid(format!(
            "the BAM reference list claims {count} entries, over the {MAX_HEADER_BYTES}-byte header budget"
        )));
    }
    let mut references = Vec::new();
    for _ in 0..count {
        let l_name = stream.read_i32()?;
        let l_name = u64::try_from(l_name)
            .ok()
            .filter(|&n| n >= 1 && n + 4 <= budget)
            .ok_or_else(|| invalid(format!("a BAM reference name claims {l_name} bytes")))?;
        budget -= l_name + 4;
        let mut name = Vec::new();
        read_append(&mut ReadAdapter(stream), &mut name, l_name as usize)?;
        name.pop();
        let length = stream.read_i32()?;
        references.push((name, length));
    }
    let first_record = stream.file_pointer();
    Ok(BamStart {
        first_record,
        references,
        has_records: stream.available()? > 0,
    })
}

/// [`Inflating::read_exact`] as a [`Read`], for [`read_append`].
struct ReadAdapter<'a, R>(&'a mut Inflating<R>);

impl<R: Read + Seek> Read for ReadAdapter<'_, R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        self.0.read_exact(out).map(|()| out.len())
    }
}

/// One BAM to copy: its file, size and termination, and where its records
/// start.
pub struct Source {
    path: PathBuf,
    size: u64,
    termination: Termination,
    start: BamStart,
}

impl Source {
    /// Read how `path` ends and its header. A file htsjdk finds
    /// [`Termination::Defective`] fails as `blockCopyBamFile` fails it, so a
    /// gather checks every input before it writes a byte.
    pub fn open(path: &Path) -> Result<Self> {
        let input = |source| Error::Input {
            path: path.to_path_buf(),
            source,
        };
        let mut file = File::open(path).map_err(input)?;
        let meta = file.metadata().map_err(input)?;
        if !meta.is_file() {
            return Err(input(io::Error::new(
                io::ErrorKind::InvalidInput,
                "a pipe or a device: a BAM is block copied from a regular file",
            )));
        }
        let size = meta.len();
        let termination = termination(&mut file, size).map_err(input)?;
        if termination == Termination::Defective {
            return Err(input(invalid(
                "does not have a valid GZIP block at the end of the file",
            )));
        }
        file.seek(SeekFrom::Start(0)).map_err(input)?;
        let mut stream = Inflating::new(BufReader::new(file), size);
        let start = read_start(&mut stream).map_err(input)?;
        Ok(Self {
            path: path.to_path_buf(),
            size,
            termination,
            start,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn start(&self) -> &BamStart {
        &self.start
    }

    pub fn termination(&self) -> Termination {
        self.termination
    }

    /// `blockCopyBamFile`: with `skip_header`, the rest of the block the
    /// first record starts in, deflated again at `level`, then every later
    /// block as bytes; without it, every block from the first. With
    /// `skip_terminator`, an end-of-file block at the end is left out.
    pub fn copy(
        &self,
        out: &mut dyn Write,
        skip_header: bool,
        skip_terminator: bool,
        level: u8,
    ) -> Result<CopyStats> {
        let input = |source| Error::Input {
            path: self.path.clone(),
            source,
        };
        let mut stats = CopyStats::default();
        let mut file = File::open(&self.path).map_err(input)?;
        let from = if skip_header {
            let (address, offset) = self.start.first_record;
            let reader = BufReader::new(file.try_clone().map_err(input)?);
            let mut stream = Inflating::new(reader, self.size);
            stream.seek(address, offset).map_err(input)?;
            let remaining = stream.available().map_err(input)?;
            let block = stream.current.as_ref().expect("seek reads a block");
            let rest = &block.data[stream.offset..stream.offset + remaining];
            stats.deflated_blocks += write_blocks(out, rest, level).map_err(output_error)?;
            stream.offset += remaining;
            stream.file_pointer().0
        } else {
            0
        };
        let end = match self.termination {
            Termination::Terminator if skip_terminator => self.size - EOF_BLOCK.len() as u64,
            _ => self.size,
        };
        if from < end {
            file.seek(SeekFrom::Start(from)).map_err(input)?;
            let mut reader = BufReader::with_capacity(COPY_BUFFER, file).take(end - from);
            let mut raw = Vec::with_capacity(MAX_BLOCK_BYTES);
            let mut at = from;
            while read_raw_block(&mut reader, &mut raw)
                .map_err(|e| input(invalid(format!("at byte {at}: {e}"))))?
            {
                out.write_all(&raw).map_err(output_error)?;
                at += raw.len() as u64;
                stats.copied_blocks += 1;
            }
            stats.ends_with_eof = end == self.size && self.termination == Termination::Terminator;
            if at != end {
                return Err(input(invalid(format!(
                    "the BGZF blocks end at byte {at}, not at byte {end}"
                ))));
            }
        }
        Ok(stats)
    }
}

/// What one [`Source::copy`] wrote.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct CopyStats {
    /// Blocks deflated again (step 2).
    pub deflated_blocks: u64,
    /// Blocks copied as bytes (step 3).
    pub copied_blocks: u64,
    /// The copy ended with the input's end-of-file block.
    pub ends_with_eof: bool,
}

/// A write error; the caller names the output.
fn output_error(source: io::Error) -> Error {
    Error::Record(source)
}

/// `data` as `BlockCompressedOutputStream` writes it and `flush` ends it:
/// blocks of at most [`BLOCK_DATA_BYTES`] bytes deflated at `level`, none
/// for no byte. Returns the number of blocks.
pub fn write_blocks(out: &mut dyn Write, data: &[u8], level: u8) -> io::Result<u64> {
    let mut blocks = 0;
    let mut payload = Vec::with_capacity(MAX_BLOCK_BYTES);
    for chunk in data.chunks(BLOCK_DATA_BYTES) {
        deflate(chunk, level, &mut payload)?;
        let total = HEADER_BYTES + payload.len() + FOOTER_BYTES;
        let bsize =
            u16::try_from(total - 1).map_err(|_| invalid("a deflated block over 64 KiB"))?;
        out.write_all(&PREAMBLE)?;
        out.write_all(&bsize.to_le_bytes())?;
        out.write_all(&payload)?;
        out.write_all(&zlib_rs::crc32::crc32(0, chunk).to_le_bytes())?;
        out.write_all(&(chunk.len() as u32).to_le_bytes())?;
        blocks += 1;
    }
    Ok(blocks)
}

/// `chunk` as a raw deflate stream into `payload`, as
/// `BlockCompressedOutputStream`'s `Deflater(level, true)` writes it: one
/// stored block at level 0 (the bytes zlib writes there), else zlib-rs at
/// `level`, and a stored block when that does not fit in a block, as htsjdk
/// falls back to its no-compression deflater.
fn deflate(chunk: &[u8], level: u8, payload: &mut Vec<u8>) -> io::Result<()> {
    payload.clear();
    if level > 0 {
        let config = DeflateConfig {
            level: i32::from(level.min(9)),
            method: Method::Deflated,
            window_bits: -15,
            mem_level: 8,
            strategy: Strategy::Default,
        };
        payload.resize(MAX_BLOCK_BYTES - HEADER_BYTES - FOOTER_BYTES, 0);
        let (done, code) = zlib_rs::compress_slice(payload, chunk, config);
        let length = done.len();
        if code == ReturnCode::Ok {
            payload.truncate(length);
            return Ok(());
        }
        payload.clear();
    }
    let length = u16::try_from(chunk.len()).map_err(|_| invalid("a stored block over 65535"))?;
    payload.push(1);
    payload.extend_from_slice(&length.to_le_bytes());
    payload.extend_from_slice(&(!length).to_le_bytes());
    payload.extend_from_slice(chunk);
    Ok(())
}

/// `file` behind the 1 MiB buffer a copy writes through.
pub fn writer(file: File) -> BufWriter<File> {
    BufWriter::with_capacity(COPY_BUFFER, file)
}

#[cfg(test)]
mod tests;

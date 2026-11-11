//! The BAM index (`.bai`): a reader, and a writer whose bytes are the ones
//! htsjdk 5.0.0 `BAMIndexer` writes for the same BAM bytes.
//!
//! # The file
//!
//! `BAI\1`, the reference count, then per reference: the bin count, each
//! bin (number, chunk count, chunks as two virtual offsets), the interval
//! count and one virtual offset per 16384-base window; last, the count of
//! records without a coordinate. A virtual offset is the file offset of a
//! BGZF block shifted left 16 bits, plus the offset inside the inflated
//! block.
//!
//! # What htsjdk writes ([`Indexer`])
//!
//! `BAMIndexer.processAlignment`, `BinningIndexBuilder.processFeature` and
//! `BinaryBAMIndexWriter.writeReference`, read from Picard 3.5.0's output on
//! the hand-written BAMs of `tests/concordance/index_bams.py`:
//!
//! - The bin is computed from the record's start and end. The bin field
//!   stored in the record is not read. A record without an end (unmapped,
//!   or mapped with no reference base before position 2) counts as one
//!   base. An end past 2^29 is an error.
//! - A record goes into its bin as the chunk from the offset before it to
//!   the offset after it. It extends the bin's last chunk when that chunk
//!   ends in the BGZF block the record starts in; otherwise it starts a new
//!   chunk. No other merging is done: no bin is folded into its parent.
//! - The linear index holds, per window the record overlaps, the smallest
//!   start offset. A record without an end uses the window of the base
//!   before its start. Windows no record touched take the value of the
//!   window before them; the array ends at the last window touched.
//! - A reference with a bin also gets the pseudo-bin 37450: the first and
//!   last offset of its records, then the counts of records without and
//!   with the unmapped flag. A reference without records is 8 zero bytes.
//! - A record whose start is 0 has no coordinate: it is counted in the
//!   last 8 bytes and nowhere else.
//!
//! [`Indexer`] and [`index_bam`] live in `bai/indexer.rs`, the reader and
//! the file layout here.
//!
//! # The offsets ([`index_bam`])
//!
//! The offset before a record is where its first byte is. The offset after
//! it is where the next byte would be read: when the record ends a block
//! that is the start of the next block, as htsjdk's
//! `BlockCompressedInputStream.getFilePointer` reports it. The index of a
//! BAM is therefore the one Picard `BuildBamIndex` writes for those bytes.
//!
//! The index Picard writes beside a BAM it is writing (`CREATE_INDEX=true`)
//! differs from its own `BuildBamIndex` of that BAM in one offset. The last
//! record of the file ends the last block, which the writer has not closed
//! yet, so the writer reports the offset inside that block (`8374551:1514`)
//! and the reader the EOF block after it (`8375239:0`). That offset is the
//! end of one chunk and the last offset of one metadata bin, and only when
//! the last record has a position. bammetrics writes the reader's value for
//! a BAM it wrote too: the `.bai` beside a BAM is always `BuildBamIndex` of
//! that BAM (`tests/concordance/test_mark_duplicates.py`).
//!
//! htsjdk stops reading at an empty BGZF block in the middle of a file and
//! indexes only the records before it. [`index_bam`] reads on and indexes
//! every record (KNOWN_DIFFERENCES.md).

use std::io::{self, Read, Write};

mod indexer;
#[cfg(test)]
mod tests;

pub use indexer::{Alignment, BamSource, Counts, Indexer, Refusal, index_bam, reference_lengths};

pub const MAGIC: [u8; 4] = *b"BAI\x01";

/// The pseudo-bin that holds a reference's metadata, htsjdk
/// `GenomicIndexUtil.MAX_BINS`.
pub const METADATA_BIN: u32 = 37450;

/// The chunks of the metadata pseudo-bin: the offsets, then the counts.
const METADATA_CHUNKS: i32 = 2;

/// The last position a BAI bin covers: 2^29.
pub const MAX_POSITION: i64 = 1 << 29;

/// A linear index window is 2^14 bases.
pub(crate) const WINDOW_SHIFT: u32 = 14;

/// The first bin of the deepest level, one bin per window.
pub(crate) const FIRST_LEAF_BIN: u32 = 4681;

/// The largest `.bai` [`Index::read`] holds in memory. A 30x human genome
/// indexed by Picard is about 9 MiB.
pub const MAX_INDEX_BYTES: u64 = 1 << 30;

/// The file from virtual offset `start` up to `end`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Chunk {
    pub start: u64,
    pub end: u64,
}

/// The metadata pseudo-bin of one reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Metadata {
    /// The smallest offset a record of the reference starts at.
    pub first_offset: u64,
    /// The largest offset a record of the reference ends at.
    pub last_offset: u64,
    /// Records without the unmapped flag.
    pub aligned: u64,
    /// Records with the unmapped flag, placed on the reference.
    pub unaligned: u64,
}

/// One reference of an index.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Reference {
    /// The bins in file order, without the metadata pseudo-bin.
    pub bins: Vec<(u32, Vec<Chunk>)>,
    pub metadata: Option<Metadata>,
    /// One offset per 16384-base window.
    pub intervals: Vec<u64>,
}

/// A `.bai` in memory.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Index {
    pub references: Vec<Reference>,
    /// The records without a coordinate; `None` for an index that ends
    /// before the count (samtools before 0.1.8).
    pub no_coordinate: Option<u64>,
}

/// What `BamIndexStats` reads of an index: the counts, no bin kept.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// The metadata of the first references, at most the number asked for.
    pub metadata: Vec<Option<Metadata>>,
    /// The references the index holds.
    pub references: usize,
    pub no_coordinate: Option<u64>,
}

fn invalid(what: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, what.into())
}

fn read_array<const N: usize>(reader: &mut impl Read) -> io::Result<[u8; N]> {
    let mut bytes = [0u8; N];
    reader.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn read_u32(reader: &mut impl Read) -> io::Result<u32> {
    read_array(reader).map(u32::from_le_bytes)
}

fn read_u64(reader: &mut impl Read) -> io::Result<u64> {
    read_array(reader).map(u64::from_le_bytes)
}

/// A count field: an `int32_t` that is not negative.
fn read_count(reader: &mut impl Read, what: &str) -> io::Result<u64> {
    let count = read_array(reader).map(i32::from_le_bytes)?;
    u64::try_from(count).map_err(|_| invalid(format!("BAI: {what} count {count}")))
}

fn read_chunk(reader: &mut impl Read) -> io::Result<Chunk> {
    Ok(Chunk {
        start: read_u64(reader)?,
        end: read_u64(reader)?,
    })
}

/// Skip `bytes` of `reader`; a shorter file is an error.
fn skip(reader: &mut impl Read, bytes: u64) -> io::Result<()> {
    let skipped = io::copy(&mut reader.by_ref().take(bytes), &mut io::sink())?;
    if skipped == bytes {
        Ok(())
    } else {
        Err(io::ErrorKind::UnexpectedEof.into())
    }
}

/// The magic and the reference count.
fn read_head(reader: &mut impl Read) -> io::Result<u64> {
    if read_array::<4>(reader)? != MAGIC {
        return Err(invalid(
            "not a BAI index: the file does not start with BAI\\1",
        ));
    }
    read_count(reader, "reference")
}

/// The metadata pseudo-bin after its bin number: two chunks.
fn read_metadata(reader: &mut impl Read) -> io::Result<Metadata> {
    let chunks = read_array(reader).map(i32::from_le_bytes)?;
    if chunks != METADATA_CHUNKS {
        return Err(invalid(format!(
            "BAI: the metadata bin holds {chunks} chunks, not {METADATA_CHUNKS}"
        )));
    }
    let offsets = read_chunk(reader)?;
    let counts = read_chunk(reader)?;
    Ok(Metadata {
        first_offset: offsets.start,
        last_offset: offsets.end,
        aligned: counts.start,
        unaligned: counts.end,
    })
}

/// One reference. Without `keep`, bins and intervals are read past and
/// only the metadata comes back. Nothing is allocated from a count before
/// its items are read, so a count the file does not back is an early end.
fn read_reference(reader: &mut impl Read, keep: bool) -> io::Result<Reference> {
    let mut reference = Reference::default();
    for _ in 0..read_count(reader, "bin")? {
        let bin = read_u32(reader)?;
        if bin == METADATA_BIN {
            if reference.metadata.replace(read_metadata(reader)?).is_some() {
                return Err(invalid("BAI: a reference holds two metadata bins"));
            }
            continue;
        }
        let chunks = read_count(reader, "chunk")?;
        if keep {
            let chunks = (0..chunks)
                .map(|_| read_chunk(reader))
                .collect::<io::Result<_>>()?;
            reference.bins.push((bin, chunks));
        } else {
            skip(reader, chunks.saturating_mul(16))?;
        }
    }
    let intervals = read_count(reader, "interval")?;
    if keep {
        reference.intervals = (0..intervals)
            .map(|_| read_u64(reader))
            .collect::<io::Result<_>>()?;
    } else {
        skip(reader, intervals.saturating_mul(8))?;
    }
    Ok(reference)
}

/// The count after the last reference; an index may end before it.
fn read_no_coordinate(reader: &mut impl Read) -> io::Result<Option<u64>> {
    match read_u64(reader) {
        Ok(count) => Ok(Some(count)),
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => Ok(None),
        Err(e) => Err(e),
    }
}

impl Summary {
    /// The metadata of the first `keep` references of the index in
    /// `reader`, the reference count and the no-coordinate count. The
    /// index is read through once and no bin is held, so memory grows with
    /// `keep` only.
    pub fn read(reader: &mut impl Read, keep: usize) -> io::Result<Self> {
        let references = read_head(reader)?;
        let mut metadata = Vec::new();
        for index in 0..references {
            let reference = read_reference(reader, false)?;
            if usize::try_from(index).is_ok_and(|i| i < keep) {
                metadata.push(reference.metadata);
            }
        }
        Ok(Self {
            metadata,
            references: usize::try_from(references).map_err(|e| invalid(e.to_string()))?,
            no_coordinate: read_no_coordinate(reader)?,
        })
    }
}

impl Index {
    /// The index in `reader`, whole. At most [`MAX_INDEX_BYTES`] are read:
    /// a longer index is an error.
    pub fn read(reader: &mut impl Read) -> io::Result<Self> {
        let mut reader = reader.take(MAX_INDEX_BYTES + 1);
        let index = Self::read_capped(&mut reader);
        if reader.limit() == 0 {
            return Err(invalid(format!(
                "BAI: the index is longer than {MAX_INDEX_BYTES} bytes"
            )));
        }
        index
    }

    fn read_capped(reader: &mut impl Read) -> io::Result<Self> {
        let references = (0..read_head(reader)?)
            .map(|_| read_reference(reader, true))
            .collect::<io::Result<_>>()?;
        Ok(Self {
            references,
            no_coordinate: read_no_coordinate(reader)?,
        })
    }

    /// The index as its file bytes.
    pub fn write(&self, out: &mut impl Write) -> io::Result<()> {
        write_head(out, self.references.len())?;
        for reference in &self.references {
            let bins = reference.bins.len() + usize::from(reference.metadata.is_some());
            write_count(out, bins)?;
            for (bin, chunks) in &reference.bins {
                write_bin(out, *bin, chunks)?;
            }
            if let Some(metadata) = &reference.metadata {
                write_metadata(out, metadata)?;
            }
            write_count(out, reference.intervals.len())?;
            for offset in &reference.intervals {
                out.write_all(&offset.to_le_bytes())?;
            }
        }
        match self.no_coordinate {
            Some(count) => out.write_all(&count.to_le_bytes()),
            None => Ok(()),
        }
    }
}

pub(crate) fn write_count(out: &mut impl Write, count: usize) -> io::Result<()> {
    let count = i32::try_from(count).map_err(|_| invalid("BAI: a count past 2^31"))?;
    out.write_all(&count.to_le_bytes())
}

pub(crate) fn write_head(out: &mut impl Write, references: usize) -> io::Result<()> {
    out.write_all(&MAGIC)?;
    write_count(out, references)
}

pub(crate) fn write_bin(out: &mut impl Write, bin: u32, chunks: &[Chunk]) -> io::Result<()> {
    out.write_all(&bin.to_le_bytes())?;
    write_count(out, chunks.len())?;
    for chunk in chunks {
        out.write_all(&chunk.start.to_le_bytes())?;
        out.write_all(&chunk.end.to_le_bytes())?;
    }
    Ok(())
}

pub(crate) fn write_metadata(out: &mut impl Write, metadata: &Metadata) -> io::Result<()> {
    let offsets = Chunk {
        start: metadata.first_offset,
        end: metadata.last_offset,
    };
    let counts = Chunk {
        start: metadata.aligned,
        end: metadata.unaligned,
    };
    write_bin(out, METADATA_BIN, &[offsets, counts])
}

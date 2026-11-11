//! The index builder: htsjdk's `BAMIndexer` over alignments in file
//! order ([`Indexer`]), and the read of a BAM that feeds it its records
//! with their virtual offsets ([`BamSource`], [`index_bam`]).

use std::fs::File;
use std::io::{self, BufRead, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf, sam};

use super::{
    Chunk, FIRST_LEAF_BIN, MAX_POSITION, Metadata, WINDOW_SHIFT, write_bin, write_count,
    write_head, write_metadata,
};
use crate::bam_raw::RawRecord;
use crate::cram;
use crate::error::{Error, Result};
use crate::pass;

/// `GenomicIndexUtil.regionToBin` for the zero-based half-open `[beg, end)`.
pub(super) fn region_to_bin(beg: i64, end: i64) -> u32 {
    let end = end - 1;
    for (shift, first) in [(14, FIRST_LEAF_BIN), (17, 585), (20, 73), (23, 9), (26, 1)] {
        if beg >> shift == end >> shift {
            return first + (beg >> shift) as u32;
        }
    }
    0
}

/// `LinearIndex.convertToLinearIndexOffset`: the window of a 1-based
/// position.
pub(super) fn window(position: i64) -> i64 {
    (position.max(1) - 1) >> WINDOW_SHIFT
}

/// The compressed file offset of the block a virtual offset points into.
fn block_address(offset: u64) -> u64 {
    offset >> 16
}

/// What the indexer reads of one alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Alignment {
    /// `refID`, -1 for none.
    pub reference: i32,
    /// htsjdk `getAlignmentStart`: 1-based, 0 for none.
    pub start: i32,
    /// htsjdk `getAlignmentEnd`: 1-based and inclusive, 0 for a record
    /// with the unmapped flag.
    pub end: i32,
    pub unmapped: bool,
    /// From the offset before the record to the offset after it.
    pub chunk: Chunk,
}

impl Alignment {
    /// The alignment of `record`, which the file holds as `chunk`. The end
    /// is the start plus the reference bases of the CIGAR, less one. A
    /// CIGAR kept in the `CG` tag has a stand-in with the same reference
    /// length, so the tag is not read.
    pub fn of(record: &RawRecord, chunk: Chunk) -> Self {
        const UNMAPPED: u16 = 0x4;
        // M, D, N, = and X consume the reference.
        const CONSUMES_REFERENCE: u16 = 0b1_1000_1101;
        let unmapped = record.flags() & UNMAPPED != 0;
        let start = record.alignment_start();
        let end = if unmapped {
            0
        } else {
            let length: i64 = record
                .cigar()
                .filter(|(_, op)| CONSUMES_REFERENCE >> op & 1 == 1)
                .map(|(length, _)| i64::from(length))
                .sum();
            (i64::from(start) + length - 1) as i32
        };
        Self {
            reference: record.reference_id(),
            start,
            end,
            unmapped,
            chunk,
        }
    }
}

/// Why a record cannot be indexed. Picard fails on each of these too.
#[derive(Debug, PartialEq, Eq)]
pub enum Refusal {
    /// The record's reference comes before the one being indexed, or it
    /// has a position and no reference.
    UnexpectedReference { reference: i32, indexing: i64 },
    /// The header has no such reference.
    UnknownReference(i32),
    /// The alignment ends past 2^29.
    PositionTooHigh,
    /// The alignment start is negative.
    NegativeStart(i32),
    /// The start lies in a window past the end of the reference.
    PastReference { bin: u32, last: usize },
}

impl std::fmt::Display for Refusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Refusal::UnexpectedReference {
                reference,
                indexing,
            } => write!(
                f,
                "Unexpected reference {reference} when constructing index for {indexing}"
            ),
            Refusal::UnknownReference(reference) => write!(
                f,
                "Reference name for '{reference}' not found in sequence dictionary."
            ),
            Refusal::PositionTooHigh => f.write_str("Read position too high for BAI bin indexing."),
            Refusal::NegativeStart(start) => write!(f, "alignment start {start} is negative"),
            Refusal::PastReference { bin, last } => write!(
                f,
                "bin {bin} is past {last}, the last bin of the reference (the record starts \
                 past its end)"
            ),
        }
    }
}

/// What an [`Indexer`] counted, per reference of the header: the rows
/// `BamIndexStats` prints for the index written.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    /// Records without and with the unmapped flag, per reference.
    pub references: Vec<(u64, u64)>,
    pub no_coordinate: u64,
    /// Every record given to the indexer.
    pub records: u64,
}

/// The index of one reference while its records arrive.
struct Building {
    /// `BAMIndexBuilder.currentReference`: -1 for a header without
    /// references.
    reference: i64,
    /// The chunks of each bin, by bin number: grown to the largest bin met
    /// and emptied through `seen`, so a reference costs what its records
    /// touch and nothing for its length.
    bins: Vec<Vec<Chunk>>,
    /// The bins that hold a chunk, in the order met.
    seen: Vec<u32>,
    /// The bins the reference's length allows: htsjdk sizes its bin array
    /// by it, and a record past it fails.
    slots: usize,
    /// The smallest start offset per window, 0 for a window not touched.
    windows: Vec<u64>,
    /// `largestIndexSeen`.
    last_window: i64,
    metadata: Metadata,
}

const NO_RECORDS: Metadata = Metadata {
    first_offset: u64::MAX,
    last_offset: 0,
    aligned: 0,
    unaligned: 0,
};

impl Building {
    /// `startNewReference`: the next reference, `length` bases long.
    fn start(&mut self, length: usize) {
        self.reference += 1;
        // `getMaxBinNumberForSequenceLength` plus one.
        self.slots = FIRST_LEAF_BIN as usize + (length >> WINDOW_SHIFT) + 1;
        for bin in self.seen.drain(..) {
            self.bins[bin as usize].clear();
        }
        self.windows.clear();
        self.last_window = -1;
        self.metadata = NO_RECORDS;
    }

    /// `BAMIndexMetaData.recordMetaData` for a record with a start.
    fn count(&mut self, alignment: &Alignment) {
        if alignment.unmapped {
            self.metadata.unaligned += 1;
        } else {
            self.metadata.aligned += 1;
        }
        self.metadata.first_offset = self.metadata.first_offset.min(alignment.chunk.start);
        self.metadata.last_offset = self.metadata.last_offset.max(alignment.chunk.end);
    }

    /// `BinningIndexBuilder.processFeature`.
    fn add(&mut self, alignment: &Alignment) -> std::result::Result<(), Refusal> {
        let start = i64::from(alignment.start);
        if start < 0 {
            return Err(Refusal::NegativeStart(alignment.start));
        }
        // `SAMRecord.computeIndexingBin`.
        let end = match i64::from(alignment.end) {
            end if end <= 0 => start,
            end => end,
        };
        if end > MAX_POSITION {
            return Err(Refusal::PositionTooHigh);
        }
        let bin = region_to_bin(start - 1, end);
        let slot = bin as usize;
        if slot >= self.slots {
            return Err(Refusal::PastReference {
                bin,
                last: self.slots.saturating_sub(1),
            });
        }
        if self.bins.len() <= slot {
            self.bins.resize_with(slot + 1, Vec::new);
        }
        let chunks = &mut self.bins[slot];
        let chunk = alignment.chunk;
        match chunks.last_mut() {
            // `areInSameOrAdjacentBlocks`: the same block, or the block one
            // byte further on, which no BGZF file has.
            Some(previous)
                if block_address(chunk.start).wrapping_sub(block_address(previous.end)) <= 1 =>
            {
                previous.end = chunk.end;
            }
            Some(_) => chunks.push(chunk),
            None => {
                chunks.push(chunk);
                self.seen.push(bin);
            }
        }
        let (first, last) = if alignment.end == 0 {
            let only = window(start - 1);
            (only, only)
        } else {
            (window(start), window(i64::from(alignment.end)))
        };
        self.last_window = self.last_window.max(last);
        for index in first..=last {
            let index = index as usize;
            if self.windows.len() <= index {
                self.windows.resize(index + 1, 0);
            }
            let slot = &mut self.windows[index];
            if *slot == 0 || chunk.start < *slot {
                *slot = chunk.start;
            }
        }
        Ok(())
    }

    /// `BinaryBAMIndexWriter.writeReference`; returns the counts written.
    fn write(&mut self, out: &mut impl Write) -> io::Result<(u64, u64)> {
        if self.seen.is_empty() {
            // No bin and no interval.
            out.write_all(&[0u8; 8])?;
            return Ok((0, 0));
        }
        write_count(out, self.seen.len() + 1)?;
        self.seen.sort_unstable();
        for &bin in &self.seen {
            write_bin(out, bin, &self.bins[bin as usize])?;
        }
        write_metadata(out, &self.metadata)?;
        let intervals = (self.last_window + 1) as usize;
        write_count(out, intervals)?;
        let mut previous = 0u64;
        for index in 0..intervals {
            match self.windows.get(index) {
                Some(&offset) if offset != 0 => previous = offset,
                _ => {}
            }
            out.write_all(&previous.to_le_bytes())?;
        }
        Ok((self.metadata.aligned, self.metadata.unaligned))
    }
}

/// The header's sequence dictionary as the lengths [`Indexer::new`] takes.
pub fn reference_lengths(header: &sam::Header) -> Vec<usize> {
    header
        .reference_sequences()
        .values()
        .map(|reference| usize::from(reference.length()))
        .collect()
}

/// htsjdk `BAMIndexer`: alignments in file order in, the `.bai` bytes out.
/// Each reference is written when the first record of a later one arrives,
/// so memory holds one reference's bins.
pub struct Indexer<W: Write> {
    out: W,
    /// The length of each reference of the header.
    lengths: Vec<usize>,
    /// `BAMIndexer.currentReference`: the references before it are written.
    current: usize,
    building: Building,
    counts: Counts,
}

impl<W: Write> Indexer<W> {
    /// Start the index of a BAM whose header lists references of
    /// `lengths` bases: the magic and the reference count are written.
    pub fn new(mut out: W, lengths: Vec<usize>) -> io::Result<Self> {
        write_head(&mut out, lengths.len())?;
        let mut building = Building {
            reference: -1,
            bins: Vec::new(),
            seen: Vec::new(),
            slots: 0,
            windows: Vec::new(),
            last_window: -1,
            metadata: NO_RECORDS,
        };
        if let Some(&first) = lengths.first() {
            building.start(first);
        }
        Ok(Self {
            out,
            lengths,
            current: 0,
            building,
            counts: Counts::default(),
        })
    }

    /// `advanceToReference`: write every reference before `next`.
    fn advance_to(&mut self, next: usize) -> io::Result<()> {
        while self.current < next {
            let written = self.building.write(&mut self.out)?;
            self.counts.references.push(written);
            self.current += 1;
            if let Some(&length) = self.lengths.get(self.current) {
                self.building.start(length);
            }
        }
        Ok(())
    }

    /// `processAlignment`: the next record of the file.
    pub fn add(&mut self, alignment: &Alignment) -> io::Result<std::result::Result<(), Refusal>> {
        self.counts.records += 1;
        let reference = match usize::try_from(alignment.reference) {
            Ok(reference) if reference < self.lengths.len() => Some(reference),
            Err(_) if alignment.reference == -1 => None,
            _ => return Ok(Err(Refusal::UnknownReference(alignment.reference))),
        };
        if let Some(reference) = reference {
            self.advance_to(reference)?;
        }
        if alignment.start == 0 {
            self.counts.no_coordinate += 1;
            return Ok(Ok(()));
        }
        // A position without a reference fails here too: htsjdk compares
        // -1 with the reference it is indexing.
        if reference.is_none() || i64::from(alignment.reference) != self.building.reference {
            return Ok(Err(Refusal::UnexpectedReference {
                reference: alignment.reference,
                indexing: self.building.reference,
            }));
        }
        self.building.count(alignment);
        Ok(self.building.add(alignment))
    }

    /// `finish`: the remaining references and the no-coordinate count.
    pub fn finish(mut self) -> io::Result<(W, Counts)> {
        self.advance_to(self.lengths.len())?;
        self.out
            .write_all(&self.counts.no_coordinate.to_le_bytes())?;
        self.out.flush()?;
        Ok((self.out, self.counts))
    }
}

/// A BAM opened for indexing: its header read, its records not yet.
pub struct BamSource {
    path: PathBuf,
    header: sam::Header,
    reader: bgzf::io::MultithreadedReader<Box<dyn Read + Send>>,
}

impl BamSource {
    /// Open the BAM at `path` and read its header; `threads` workers
    /// inflate its blocks. A CRAM, a SAM and any other file are refused: a
    /// `.bai` indexes a BAM.
    pub fn open(path: &Path, threads: usize) -> Result<Self> {
        let input = |source| Error::Input {
            path: path.to_path_buf(),
            source,
        };
        let sniffed = cram::sniff(File::open(path).map_err(input)?).map_err(input)?;
        if sniffed.cram {
            return Err(Error::Usage(format!(
                "{}: a .bai indexes a BAM and this is a CRAM, whose index is a .crai (Picard: \
                 Input file must be bam file, not sam file.)",
                path.display()
            )));
        }
        let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        let inner: Box<dyn Read + Send> = Box::new(sniffed.reader);
        let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
            workers, inner,
        ));
        let header = reader
            .read_header()
            .map_err(|e| input(cram::bam_header_error(e)))?;
        Ok(Self {
            path: path.to_path_buf(),
            header,
            reader: reader.into_inner(),
        })
    }

    pub fn header(&self) -> &sam::Header {
        &self.header
    }

    /// Read every record and write the index of the file to `out`, the
    /// file `target` names in a write error.
    pub fn index<W: Write>(mut self, out: W, target: &Path) -> Result<(W, Counts)> {
        let path = self.path.clone();
        let read = |source| pass::read_error(&path, false, source);
        let written = |source| Error::Output {
            path: target.to_path_buf(),
            source,
        };
        let mut indexer = Indexer::new(out, reference_lengths(&self.header)).map_err(written)?;
        let mut record = RawRecord::default();
        loop {
            // The offset before a record is where its first byte is: an
            // exhausted block is left for the next one that holds a byte.
            if self.reader.fill_buf().map_err(read)?.is_empty() {
                break;
            }
            let start = u64::from(self.reader.virtual_position());
            if !record.read_from(&mut self.reader).map_err(read)? {
                break;
            }
            let end = u64::from(self.reader.virtual_position());
            let alignment = Alignment::of(&record, Chunk { start, end });
            if let Err(refusal) = indexer.add(&alignment).map_err(written)? {
                return Err(Error::Metric(format!(
                    "{}: Exception creating BAM index for record {}: {refusal}",
                    path.display(),
                    String::from_utf8_lossy(record.name())
                )));
            }
        }
        self.reader.finish().map_err(read)?;
        indexer.finish().map_err(written)
    }
}

/// The `.bai` of the BAM at `path`, written to `out`, the file `target`
/// names in a write error.
pub fn index_bam<W: Write>(
    path: &Path,
    threads: usize,
    out: W,
    target: &Path,
) -> Result<(W, Counts)> {
    BamSource::open(path, threads)?.index(out, target)
}

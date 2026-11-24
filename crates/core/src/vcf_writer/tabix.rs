//! The `.tbi` htsjdk 5.0.0 writes beside a block-compressed VCF
//! (`TabixIndexCreator` with `TabixFormat.VCF`, `BinningIndexBuilder`,
//! `TabixIndex.write`), from the virtual offset of each record line.
//!
//! - References are numbered as their first record arrives; one met again
//!   after another fails (`Sequence ... added out of order`), and so does a
//!   record that starts before the one in front of it on its reference.
//! - Each record is a feature from `POS` to its end (`INFO/END`, else the
//!   reference length), its chunk from its offset to the next record's.
//!   Bins and the 16 kb linear windows are those of the BAM index
//!   ([`crate::bai`]); a chunk joins the one before it in its bin when it
//!   starts in the block that one ends in, or the block after.
//! - With `##contig` lines a reference takes its length from its line, and
//!   a reference the lines lack fails (htsjdk's `NullPointerException`).
//! - The file: `TBI\1`, the reference count, format 2 (VCF), sequence
//!   column 1, start column 2, end column 0, meta character `#`, 0 lines
//!   to skip, the names each ended by a NUL, then per reference its bins
//!   with their chunks and its linear index, with no pseudo-bin and no
//!   count of unplaced records. The whole is BGZF-compressed.

use std::collections::HashMap;
use std::io::{self, Write};

use crate::bai::{Chunk, region_to_bin, window};

/// The first bin of the 16 kb leaves.
const FIRST_LEAF_BIN: usize = 4681;
/// `GenomicIndexUtil.MAX_BINS`.
const MAX_BINS: usize = 37450;

/// One reference while its features arrive.
struct Building {
    bins: Vec<Vec<Chunk>>,
    seen: Vec<u32>,
    /// The bins the reference's length allows.
    slots: usize,
    windows: Vec<u64>,
    last_window: i64,
}

impl Building {
    fn new(length: i64) -> Self {
        let slots = if length <= 0 {
            MAX_BINS + 1
        } else {
            FIRST_LEAF_BIN + (length as usize >> 14) + 1
        };
        Self {
            bins: Vec::new(),
            seen: Vec::new(),
            slots,
            windows: Vec::new(),
            last_window: -1,
        }
    }

    /// `BinningIndexBuilder.processFeature`.
    fn add(&mut self, start: i32, end: i32, chunk: Chunk) -> Result<(), String> {
        let start = i64::from(start);
        let end_or_start = if end <= 0 { start } else { i64::from(end) };
        let bin = region_to_bin(start - 1, end_or_start);
        let slot = bin as usize;
        if slot >= self.slots {
            return Err(format!(
                "the record at {start} is past the length of its ##contig line (bin {bin})"
            ));
        }
        if self.bins.len() <= slot {
            self.bins.resize_with(slot + 1, Vec::new);
        }
        let chunks = &mut self.bins[slot];
        match chunks.last_mut() {
            Some(previous) if (chunk.start >> 16).wrapping_sub(previous.end >> 16) <= 1 => {
                previous.end = chunk.end;
            }
            Some(_) => chunks.push(chunk),
            None => {
                chunks.push(chunk);
                self.seen.push(bin);
            }
        }
        let (first, last) = if end == 0 {
            let only = window(start - 1);
            (only, only)
        } else {
            (window(start), window(i64::from(end)))
        };
        self.last_window = self.last_window.max(last);
        for index in first.max(0)..=last {
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

    /// `TabixIndex.writeSequence`.
    fn write(mut self, out: &mut Vec<u8>) {
        if self.seen.is_empty() {
            out.extend_from_slice(&0i32.to_le_bytes());
            out.extend_from_slice(&0i32.to_le_bytes());
            return;
        }
        out.extend_from_slice(&(self.seen.len() as i32).to_le_bytes());
        self.seen.sort_unstable();
        for &bin in &self.seen {
            let chunks = &self.bins[bin as usize];
            out.extend_from_slice(&bin.to_le_bytes());
            out.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
            for chunk in chunks {
                out.extend_from_slice(&chunk.start.to_le_bytes());
                out.extend_from_slice(&chunk.end.to_le_bytes());
            }
        }
        let intervals = (self.last_window + 1).max(0) as usize;
        out.extend_from_slice(&(intervals as i32).to_le_bytes());
        let mut previous = 0u64;
        for index in 0..intervals {
            if let Some(&offset) = self.windows.get(index)
                && offset != 0
            {
                previous = offset;
            }
            out.extend_from_slice(&previous.to_le_bytes());
        }
    }
}

/// A record waiting for the offset where the next one starts.
struct Pending {
    reference: usize,
    start: i32,
    end: i32,
    offset: u64,
}

/// `TabixIndexCreator`: features in file order in, the index bytes out.
pub(crate) struct TabixBuilder {
    /// The `##contig` lengths by name, `None` without `##contig` lines.
    lengths: Option<HashMap<String, i64>>,
    names: Vec<String>,
    index_of: HashMap<String, usize>,
    /// The references written so far.
    written: Vec<u8>,
    current: Option<Building>,
    pending: Option<Pending>,
}

impl TabixBuilder {
    pub(crate) fn new(dictionary: Option<&crate::vcf::Dictionary>) -> Self {
        Self {
            lengths: dictionary.map(|d| {
                d.iter()
                    .map(|c| (c.name.clone(), c.length.unwrap_or(0)))
                    .collect()
            }),
            names: Vec::new(),
            index_of: HashMap::new(),
            written: Vec::new(),
            current: None,
            pending: None,
        }
    }

    /// `finalizeFeature`: the waiting record's chunk ends at `end`.
    fn finish_pending(&mut self, end: u64) -> Result<(), String> {
        let Some(pending) = self.pending.take() else {
            return Ok(());
        };
        if pending.offset >= end {
            return Err(format!(
                "Feature start position {} >= feature end position {end}",
                pending.offset
            ));
        }
        let chunk = Chunk {
            start: pending.offset,
            end,
        };
        let building = self.current.as_mut().expect("a reference is open");
        building.add(pending.start, pending.end, chunk)
    }

    /// `addFeature`: the record on `contig` from `start` to `end` whose
    /// line starts at the virtual offset `offset`.
    pub(crate) fn add(
        &mut self,
        contig: &str,
        start: i32,
        end: i32,
        offset: u64,
    ) -> Result<(), String> {
        let reference = match self.index_of.get(contig) {
            Some(&index) if index + 1 == self.names.len() => index,
            Some(_) => {
                return Err(format!(
                    "Sequence {contig} added out of order: the records of a block-compressed VCF \
                     with an index must be sorted by contig"
                ));
            }
            None => self.names.len(),
        };
        if let Some(previous) = &self.pending
            && previous.reference == reference
            && previous.start > start
        {
            return Err(format!(
                "Features added out of order: previous ({contig}:{}) > next ({contig}:{start})",
                previous.start
            ));
        }
        self.finish_pending(offset)?;
        if reference == self.names.len() {
            let length = match &self.lengths {
                None => 0,
                Some(lengths) => *lengths.get(contig).ok_or_else(|| {
                    format!(
                        "contig {contig} has no ##contig line, and the index takes each \
                         contig's length from its line (htsjdk fails with a \
                         NullPointerException)"
                    )
                })?,
            };
            if let Some(done) = self.current.take() {
                done.write(&mut self.written);
            }
            self.current = Some(Building::new(length));
            self.index_of.insert(contig.to_owned(), reference);
            self.names.push(contig.to_owned());
        }
        self.pending = Some(Pending {
            reference,
            start,
            end,
            offset,
        });
        Ok(())
    }

    /// `finalizeIndex(end)` and `TabixIndex.write`: the index's bytes,
    /// before compression.
    pub(crate) fn finish(mut self, end: u64) -> Result<Vec<u8>, String> {
        self.finish_pending(end)?;
        if let Some(done) = self.current.take() {
            done.write(&mut self.written);
        }
        let mut out = Vec::with_capacity(self.written.len() + 64);
        out.extend_from_slice(b"TBI\x01");
        let header = [self.names.len() as i32, 2, 1, 2, 0, i32::from(b'#'), 0];
        for value in header {
            out.extend_from_slice(&value.to_le_bytes());
        }
        let names_length: usize = self.names.iter().map(|n| n.len() + 1).sum();
        out.extend_from_slice(&(names_length as i32).to_le_bytes());
        for name in &self.names {
            out.extend_from_slice(name.as_bytes());
            out.push(0);
        }
        out.extend_from_slice(&self.written);
        Ok(out)
    }
}

/// Write `bytes` BGZF-compressed, with the end-of-file block.
pub(crate) fn write_compressed(out: impl Write, bytes: &[u8]) -> io::Result<()> {
    let mut writer = noodles::bgzf::io::Writer::new(out);
    writer.write_all(bytes)?;
    writer.finish()?.flush()
}

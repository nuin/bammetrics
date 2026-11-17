//! The tallies `doctor` keeps over the first records of a file.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet};
use std::io::{self, Read};

use crate::bam_raw::RawRecord;
use crate::sam_order::CoordinateKey;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const SECONDARY: u16 = 0x100;
const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// The most distinct unknown read group IDs kept for the message.
const SHOWN_IDS: usize = 3;

/// A record that breaks an order: its number (1-based) and what it shows.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Break {
    pub record: u64,
    pub what: String,
}

/// What the first records say.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    /// Records read.
    pub read: u64,
    /// The read error that stopped the read before the limit.
    pub error: Option<String>,
    /// Records neither secondary nor supplementary.
    pub primary: u64,
    pub paired: u64,
    pub mapped: u64,
    pub duplicate: u64,
    /// Primary records with `SEQ`, and their bases.
    pub with_bases: u64,
    pub bases: u64,
    /// Primary records with `QUAL`, and the range of their qualities.
    pub with_qualities: u64,
    pub quality_min: u8,
    pub quality_max: u8,
    /// Mapped primary records, and those with an `MD` and an `NM` tag.
    pub mapped_primary: u64,
    pub with_md: u64,
    pub with_nm: u64,
    /// Records without an `RG` tag.
    pub without_rg: u64,
    /// Records whose `RG` the header lacks, by ID (at most [`SHOWN_IDS`]
    /// IDs kept), and their number.
    pub unknown_rg: BTreeMap<String, u64>,
    pub unknown_rg_records: u64,
    /// The first record out of coordinate order (reference, then start;
    /// unmapped records last).
    pub coordinate_break: Option<Break>,
    /// The first record whose name sorts before the previous one, byte by
    /// byte (htsjdk's queryname order) and in samtools' natural order.
    pub name_break: Option<Break>,
    pub natural_break: Option<Break>,
    /// Whether any record has a reference.
    pub placed: u64,
}

/// samtools' `strnum_cmp`: digit runs compared as numbers (leading zeros
/// skipped, then the longer run is larger), other bytes as bytes.
pub fn natural_order(a: &[u8], b: &[u8]) -> Ordering {
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        if a[i].is_ascii_digit() && b[j].is_ascii_digit() {
            while i < a.len() && a[i] == b'0' {
                i += 1;
            }
            while j < b.len() && b[j] == b'0' {
                j += 1;
            }
            let start = (i, j);
            while i < a.len() && a[i].is_ascii_digit() {
                i += 1;
            }
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
            }
            let (x, y) = (&a[start.0..i], &b[start.1..j]);
            let order = x.len().cmp(&y.len()).then_with(|| x.cmp(y));
            if order != Ordering::Equal {
                return order;
            }
        } else {
            match a[i].cmp(&b[j]) {
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                other => return other,
            }
        }
    }
    (a.len() - i).cmp(&(b.len() - j))
}

/// The tally as records come, with the header's read group IDs and
/// reference names.
pub struct Tally<'a> {
    pub stats: Stats,
    read_groups: &'a HashSet<String>,
    references: &'a [String],
    previous: Option<RawRecord>,
}

impl<'a> Tally<'a> {
    pub fn new(read_groups: &'a HashSet<String>, references: &'a [String]) -> Self {
        Self {
            stats: Stats {
                quality_min: u8::MAX,
                ..Stats::default()
            },
            read_groups,
            references,
            previous: None,
        }
    }

    fn place(&self, record: &RawRecord) -> String {
        let id = record.reference_id();
        let name = usize::try_from(id)
            .ok()
            .and_then(|i| self.references.get(i))
            .map_or("*", String::as_str);
        format!("{name}:{}", record.alignment_start())
    }

    /// Read records from `stream` (each with its `block_size`) until
    /// `limit` are read, the stream ends or a read fails.
    pub fn read_from(&mut self, stream: &mut (dyn Read + '_), limit: u64) {
        let mut record = RawRecord::default();
        while self.stats.read < limit {
            match record.read_from(stream) {
                Ok(true) => {
                    if let Err(e) = self.add(&record) {
                        self.stats.error = Some(format!("record {}: {e}", self.stats.read + 1));
                        return;
                    }
                }
                Ok(false) => return,
                Err(e) => {
                    let what = if e.kind() == io::ErrorKind::UnexpectedEof {
                        "the file ends inside the record".to_owned()
                    } else {
                        e.to_string()
                    };
                    self.stats.error = Some(format!("record {}: {what}", self.stats.read + 1));
                    return;
                }
            }
        }
    }

    /// Count one record.
    pub fn add(&mut self, record: &RawRecord) -> io::Result<()> {
        let flags = record.flags();
        let s = &mut self.stats;
        s.read += 1;
        let n = s.read;
        if record.reference_id() >= 0 {
            s.placed += 1;
        }
        if flags & (SECONDARY | SUPPLEMENTARY) == 0 {
            s.primary += 1;
            s.paired += u64::from(flags & PAIRED != 0);
            s.duplicate += u64::from(flags & DUPLICATE != 0);
            let length = record.sequence_length();
            if length > 0 {
                s.with_bases += 1;
                s.bases += length as u64;
            }
            let qualities = record.qualities();
            if qualities.first().is_some_and(|&q| q != 0xff) {
                s.with_qualities += 1;
                for &q in qualities {
                    s.quality_min = s.quality_min.min(q);
                    s.quality_max = s.quality_max.max(q);
                }
            }
            if flags & UNMAPPED == 0 {
                s.mapped += 1;
                s.mapped_primary += 1;
                s.with_md += u64::from(record.tag(*b"MD")?.is_some());
                s.with_nm += u64::from(record.tag(*b"NM")?.is_some());
            }
        }
        match record.tag(*b"RG")? {
            None => s.without_rg += 1,
            Some(value) => {
                let id = String::from_utf8_lossy(value.as_str().unwrap_or_default()).into_owned();
                if !self.read_groups.contains(&id) {
                    s.unknown_rg_records += 1;
                    if s.unknown_rg.len() < SHOWN_IDS || s.unknown_rg.contains_key(&id) {
                        *s.unknown_rg.entry(id).or_default() += 1;
                    }
                }
            }
        }
        if let Some(previous) = &self.previous {
            let (before, here) = (CoordinateKey::of(previous), CoordinateKey::of(record));
            if s.coordinate_break.is_none() && before.file_order(&here) == Ordering::Greater {
                let what = format!(
                    "record {n} ({}) at {} comes after {}",
                    String::from_utf8_lossy(record.name()),
                    self.place(record),
                    self.place(previous)
                );
                self.stats.coordinate_break = Some(Break { record: n, what });
            }
            let names = (previous.name(), record.name());
            let shown = || {
                format!(
                    "record {n} named {} comes after {}",
                    String::from_utf8_lossy(names.1),
                    String::from_utf8_lossy(names.0)
                )
            };
            let s = &mut self.stats;
            if s.name_break.is_none() && names.0 > names.1 {
                s.name_break = Some(Break {
                    record: n,
                    what: shown(),
                });
            }
            if s.natural_break.is_none() && natural_order(names.0, names.1) == Ordering::Greater {
                s.natural_break = Some(Break {
                    record: n,
                    what: shown(),
                });
            }
        }
        match &mut self.previous {
            Some(previous) => previous.clone_from(record),
            None => self.previous = Some(record.clone()),
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn natural_order_compares_digit_runs_as_numbers() {
        assert_eq!(natural_order(b"r2", b"r10"), Ordering::Less);
        assert_eq!(natural_order(b"r10", b"r2"), Ordering::Greater);
        assert_eq!(natural_order(b"r002", b"r2"), Ordering::Equal);
        assert_eq!(natural_order(b"a", b"ab"), Ordering::Less);
        assert_eq!(natural_order(b"b1", b"a9"), Ordering::Greater);
        assert!(b"r10".as_slice() < b"r2".as_slice());
    }
}

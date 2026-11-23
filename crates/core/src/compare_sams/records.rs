//! The records of `SamComparison`: the primary records of one file in file
//! order (`SecondaryOrSupplementarySkippingIterator`), the key two records
//! are paired by (`PrimaryAlignmentKey`) and the map that holds the records
//! still waiting for their partner (a `LinkedHashMap` in Picard).

use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use std::io::{self, Read};
use std::path::Path;

use crate::bam_raw::{Fields, RawRecord, TagIter, Tags};
use crate::error::{Error, Result};
use crate::pass;

const PAIRED: u16 = 0x1;
pub const UNMAPPED: u16 = 0x4;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
pub const DUPLICATE: u16 = 0x400;
const SUPPLEMENTARY: u16 = 0x800;

/// `PrimaryAlignmentKey.PairStatus`, in its declaration order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Pair {
    Unpaired,
    First,
    Second,
}

/// `PrimaryAlignmentKey`: the read name and the end of the pair. A paired
/// record that is not flagged second of its pair is `First`. One
/// allocation holds both, the pair status in the first byte, so a waiting
/// record costs its name and 16 bytes of key.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Key(Box<[u8]>);

impl Key {
    pub fn new(name: &[u8], pair: Pair) -> Self {
        let mut bytes = Vec::with_capacity(name.len() + 1);
        bytes.push(pair as u8);
        bytes.extend_from_slice(name);
        Self(bytes.into_boxed_slice())
    }

    pub fn name(&self) -> &[u8] {
        &self.0[1..]
    }

    fn pair(&self) -> u8 {
        self.0[0]
    }

    /// `compareTo`: the name as `String.compareTo` orders it (a BAM name is
    /// one character per byte), then the pair status.
    pub fn compare(&self, other: &Self) -> Ordering {
        self.name()
            .cmp(other.name())
            .then(self.pair().cmp(&other.pair()))
    }
}

/// What the comparison reads from one primary record, beside its [`Key`].
#[derive(Debug)]
pub struct Facts {
    /// The reference index in the record's own header, -1 without one.
    pub reference: i32,
    /// `getAlignmentStart`: 1-based, 0 without a position.
    pub start: i32,
    pub flags: u16,
    pub mapping_quality: u8,
    /// The record, kept only under `LENIENT_DUP`: a pair whose duplicate
    /// flags differ is sorted by `SAMRecordDuplicateComparator`, which
    /// reads the CIGAR and the `RG` and `MC` tags. A record that waits in
    /// a map holds no more than those ([`Primary::waiting`]).
    pub raw: Option<Box<RawRecord>>,
}

impl Facts {
    pub fn unmapped(&self) -> bool {
        self.flags & UNMAPPED != 0
    }

    pub fn duplicate(&self) -> bool {
        self.flags & DUPLICATE != 0
    }
}

/// One primary record: its key and its facts.
#[derive(Debug)]
pub struct Primary {
    pub key: Key,
    pub facts: Facts,
}

impl Primary {
    /// The record as a map holds it until its partner comes: its raw
    /// record, when kept, cut down by [`slim`].
    pub fn waiting(mut self) -> Self {
        if let Some(raw) = self.facts.raw.as_mut()
            && let Ok(cut) = slim(raw)
        {
            **raw = cut;
        }
        self
    }
}

/// `record` without its bases, its qualities and every tag but `RG` and
/// `MC`: what [`crate::sam_order::DuplicateKey::of`] reads, in about a
/// third of the bytes. An error (tags that do not parse) leaves the caller
/// with the whole record, so the duplicate order reports it.
pub fn slim(record: &RawRecord) -> io::Result<RawRecord> {
    let mut kept = Vec::with_capacity(2);
    for field in TagIter::new(record.tag_bytes()) {
        let (name, value) = field?;
        if name == *b"RG" || name == *b"MC" {
            kept.push((name, value));
        }
    }
    let cigar: Vec<(u32, u8)> = record.cigar().collect();
    RawRecord::from_fields(&Fields {
        name: record.name(),
        flags: record.flags(),
        reference_id: record.reference_id(),
        alignment_start: record.alignment_start(),
        mapping_quality: record.mapping_quality(),
        cigar: &cigar,
        mate_reference_id: record.mate_reference_id(),
        mate_alignment_start: record.mate_alignment_start(),
        template_length: record.template_length(),
        bases: &[],
        qualities: None,
        tags: &Tags::from_unsorted(kept),
    })
}

/// One input: its record stream and the primary record it stands on.
pub struct Side<'r> {
    path: &'r Path,
    reader: &'r mut dyn Read,
    cram: bool,
    /// The dictionary size: a reference index at or past it is an error.
    references: usize,
    keep_raw: bool,
    record: RawRecord,
    current: Option<Primary>,
}

impl<'r> Side<'r> {
    /// The stream of `path`, standing on its first primary record.
    pub fn new(
        path: &'r Path,
        reader: &'r mut dyn Read,
        cram: bool,
        references: usize,
        keep_raw: bool,
    ) -> Result<Self> {
        let mut side = Self {
            path,
            reader,
            cram,
            references,
            keep_raw,
            record: RawRecord::default(),
            current: None,
        };
        side.advance()?;
        Ok(side)
    }

    /// `hasCurrent` and `getCurrent`.
    pub fn current(&self) -> Option<&Primary> {
        self.current.as_ref()
    }

    /// The next primary record becomes the current one; `None` at the end.
    fn advance(&mut self) -> Result<()> {
        self.current = None;
        while self
            .record
            .read_from(self.reader)
            .map_err(|e| pass::read_error(self.path, self.cram, e))?
        {
            self.in_dictionary(self.record.reference_id())?;
            self.in_dictionary(self.record.mate_reference_id())?;
            let flags = self.record.flags();
            if flags & (SECONDARY | SUPPLEMENTARY) != 0 {
                continue;
            }
            self.current = Some(self.primary(flags));
            break;
        }
        Ok(())
    }

    /// htsjdk resolves the reference and mate reference names of every
    /// record it decodes, secondary and supplementary ones too: an index
    /// outside the dictionary fails the read.
    fn in_dictionary(&self, index: i32) -> Result<()> {
        let known = index == -1 || usize::try_from(index).is_ok_and(|at| at < self.references);
        if known {
            return Ok(());
        }
        Err(Error::Metric(format!(
            "{}: record {}: Reference name for '{index}' not found in sequence dictionary.",
            self.path.display(),
            String::from_utf8_lossy(self.record.name())
        )))
    }

    fn primary(&mut self, flags: u16) -> Primary {
        let pair = match (flags & PAIRED != 0, flags & SECOND != 0) {
            (false, _) => Pair::Unpaired,
            (true, true) => Pair::Second,
            (true, false) => Pair::First,
        };
        let key = Key::new(self.record.name(), pair);
        let reference = self.record.reference_id();
        let start = self.record.alignment_start();
        let mapping_quality = self.record.mapping_quality();
        let raw = self
            .keep_raw
            .then(|| Box::new(std::mem::take(&mut self.record)));
        Primary {
            key,
            facts: Facts {
                reference,
                start,
                flags,
                mapping_quality,
                raw,
            },
        }
    }

    /// The current record, and the stream moved to the next: `getCurrent`
    /// then `advance`.
    pub fn take(&mut self) -> Result<Option<Primary>> {
        let Some(primary) = self.current.take() else {
            return Ok(None);
        };
        self.advance()?;
        Ok(Some(primary))
    }

    /// [`Side::take`] when the current record is `wanted`; `None` when it
    /// is not, or at the end.
    pub fn take_if(&mut self, wanted: impl FnOnce(&Primary) -> bool) -> Result<Option<Primary>> {
        if self.current.as_ref().is_some_and(wanted) {
            self.take()
        } else {
            Ok(None)
        }
    }

    /// `countRemaining`: the primary records from the current one on.
    pub fn count_remaining(&mut self) -> Result<i64> {
        let mut count = 0;
        while self.take()?.is_some() {
            count += 1;
        }
        Ok(count)
    }
}

/// A `LinkedHashMap<PrimaryAlignmentKey, SAMRecord>`: `put` of a key the
/// map holds replaces the record and keeps the key's place, and the records
/// come back out in the order their keys went in. One entry is 48 bytes
/// and the key's allocation (the name and one byte).
#[derive(Debug, Default)]
pub struct Pending {
    map: HashMap<Key, Slot>,
    next: u64,
}

/// One waiting record and the order its key came in.
#[derive(Debug)]
struct Slot {
    at: u64,
    facts: Facts,
}

impl Pending {
    pub fn put(&mut self, primary: Primary) {
        let Primary { key, facts } = primary;
        match self.map.entry(key) {
            Entry::Occupied(mut slot) => slot.get_mut().facts = facts,
            Entry::Vacant(slot) => {
                slot.insert(Slot {
                    at: self.next,
                    facts,
                });
                self.next += 1;
            }
        }
    }

    pub fn remove(&mut self, key: &Key) -> Option<Facts> {
        self.map.remove(key).map(|slot| slot.facts)
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// Every record in the order its key went in; the map is left empty.
    pub fn drain_ordered(&mut self) -> Vec<Primary> {
        let mut entries: Vec<(u64, Primary)> = self
            .map
            .drain()
            .map(|(key, slot)| {
                let facts = slot.facts;
                (slot.at, Primary { key, facts })
            })
            .collect();
        entries.sort_unstable_by_key(|(at, _)| *at);
        entries.into_iter().map(|(_, primary)| primary).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam_order::{DuplicateKey, Libraries};

    /// The duplicate key of a slimmed record is the key of the record, on
    /// every record of three fixtures (pairs with and without `MC`, clips,
    /// read groups, unmapped and unpaired reads).
    #[test]
    fn a_slimmed_record_keeps_its_duplicate_key() {
        let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/data");
        let mut seen = 0;
        for name in ["markdup.bam", "revert.bam", "insert_size.bam"] {
            let kept = pass::with_record_stream(&data.join(name), None, 1, |header, reader, _| {
                let libraries = Libraries::new(header);
                let key = |record: &RawRecord| {
                    DuplicateKey::of(record, &libraries)
                        .map(DuplicateKey::into_owned)
                        .map_err(|e| e.to_string())
                };
                let mut record = RawRecord::default();
                let mut count = 0;
                while record.read_from(reader).map_err(Error::Record)? {
                    let cut = slim(&record).map_err(Error::Record)?;
                    assert_eq!(key(&cut), key(&record), "{name}");
                    assert!(cut.as_bytes().len() <= record.as_bytes().len());
                    count += 1;
                }
                Ok(count)
            });
            seen += kept.expect("read the fixture");
        }
        assert!(seen >= 50, "{seen} records");
    }

    #[test]
    fn a_key_orders_by_name_then_pair_and_a_waiting_entry_is_small() {
        let key = |name: &str, pair| Key::new(name.as_bytes(), pair);
        assert_eq!(
            key("a", Pair::Second).compare(&key("b", Pair::Unpaired)),
            Ordering::Less
        );
        assert_eq!(
            key("a", Pair::Unpaired).compare(&key("a", Pair::First)),
            Ordering::Less
        );
        assert_eq!(
            key("a", Pair::First).compare(&key("a", Pair::Second)),
            Ordering::Less
        );
        assert_ne!(key("a", Pair::First), key("a", Pair::Second));
        assert_eq!(key("ab", Pair::First).name(), b"ab");
        assert!(std::mem::size_of::<(Key, Slot)>() <= 48);
    }
}

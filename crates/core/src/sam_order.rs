//! The record orders of htsjdk 5.0.0, the ones Picard 3.5.0 `SortSam`
//! offers as `SORT_ORDER` plus `unsorted`: [`Order`].
//!
//! Each order is a key read from a [`RawRecord`] and a `compare` on two
//! keys, a port of the htsjdk comparator of the same name:
//!
//! - [`CoordinateKey`], `SAMRecordCoordinateComparator`: reference index
//!   with no reference (-1) last, alignment start, forward before reverse
//!   strand, then read name, flags, mapping quality, mate reference index,
//!   mate start and template length. Two records without a reference tie
//!   on the start, whatever it is.
//! - [`QuerynameKey`], `SAMRecordQueryNameComparator`: read name, then
//!   paired before unpaired, first of pair before second, forward before
//!   reverse, primary before secondary, non-supplementary before
//!   supplementary, a record without `HI` before one with it, `HI`
//!   ascending.
//! - [`DuplicateKey`], `SAMRecordDuplicateComparator`: library, reference
//!   index, unclipped 5' coordinate, pair orientation, the mate's reference
//!   and unclipped 5' coordinate (from `MC`, which every paired record with
//!   a mapped mate must carry), then the score
//!   (`TOTAL_MAPPED_REFERENCE_LENGTH`, the higher first), the read group
//!   and name, the name, first of pair before second.
//! - `unsorted` compares nothing: [`Order::compare`] gives `Equal`, so a
//!   stable sort keeps file order.
//!
//! Names are compared as Java compares the `String`s htsjdk holds:
//! `String.compareTo`, by UTF-16 code unit. htsjdk decodes a BAM read name
//! one byte per character, so the order of two names is the order of their
//! bytes, and [`RawRecord`] names compare as byte strings.
//!
//! Where htsjdk subtracts two `int`s the port subtracts with wrap-around
//! and takes the sign, so an overflow orders as it does in Java.
//!
//! Total order. The coordinate order is a total preorder over every
//! record. The queryname order is one over records whose flags htsjdk
//! calls valid: a paired record is first or second of its pair and not
//! both. Outside that domain it is not transitive (two paired records that
//! are both first and second each sort before the other), in htsjdk and
//! here. The property tests of this module draw random keys in those
//! domains. A sort must be stable and must not assume more: ties stay in
//! file order, as htsjdk's `SortingCollection` keeps them.

use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::HashMap;

use noodles::sam;
use noodles::sam::header::record::value::map::read_group::tag::LIBRARY;

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::mark_duplicates::unclipped;
use crate::sam_text::parse_cigar;

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const QC_FAIL: u16 = 0x200;
const SUPPLEMENTARY: u16 = 0x800;

/// `SAMRecordDuplicateComparator.UNKNOWN_LIBRARY_STRING`.
const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// A sort order `SortSam` takes, or `unsorted`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
}

impl Order {
    /// The `SORT_ORDER` value, as the header's `SO` spells it.
    pub fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "unsorted" => Self::Unsorted,
            "queryname" => Self::Queryname,
            "coordinate" => Self::Coordinate,
            "duplicate" => Self::Duplicate,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Unsorted => "unsorted",
            Self::Queryname => "queryname",
            Self::Coordinate => "coordinate",
            Self::Duplicate => "duplicate",
        }
    }

    /// Compare two records in this order. `libraries` is read by
    /// `duplicate` only. An error is the exception htsjdk throws while it
    /// reads the key of a record (an `HI` that is no number, a missing
    /// `MC`).
    pub fn compare(self, a: &RawRecord, b: &RawRecord, libraries: &Libraries) -> Result<Ordering> {
        Ok(match self {
            Self::Unsorted => Ordering::Equal,
            Self::Coordinate => coordinate_order(a, b),
            Self::Queryname => QuerynameKey::of(a)?.compare(&QuerynameKey::of(b)?),
            Self::Duplicate => {
                DuplicateKey::of(a, libraries)?.compare(&DuplicateKey::of(b, libraries)?)
            }
        })
    }
}

/// The sign of the Java expression `a - b` on two `int`s.
fn java_difference(a: i32, b: i32) -> Ordering {
    a.wrapping_sub(b).cmp(&0)
}

/// What `SAMRecordCoordinateComparator` reads from a record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoordinateKey<'a> {
    /// `refID`, -1 for none.
    pub reference_id: i32,
    /// 1-based, 0 for none.
    pub start: i32,
    pub flags: u16,
    pub name: Cow<'a, [u8]>,
    pub mapping_quality: u8,
    pub mate_reference_id: i32,
    pub mate_start: i32,
    pub template_length: i32,
}

impl<'a> CoordinateKey<'a> {
    pub fn of(record: &'a RawRecord) -> Self {
        Self {
            reference_id: record.reference_id(),
            start: record.alignment_start(),
            flags: record.flags(),
            name: Cow::Borrowed(record.name()),
            mapping_quality: record.mapping_quality(),
            mate_reference_id: record.mate_reference_id(),
            mate_start: record.mate_alignment_start(),
            template_length: record.template_length(),
        }
    }

    /// The key without its borrow of the record.
    pub fn into_owned(self) -> CoordinateKey<'static> {
        CoordinateKey {
            name: Cow::Owned(self.name.into_owned()),
            ..self
        }
    }

    /// `fileOrderCompare`: reference index with -1 last, then start.
    pub fn file_order(&self, other: &Self) -> Ordering {
        match (self.reference_id, other.reference_id) {
            (-1, -1) => Ordering::Equal,
            (-1, _) => Ordering::Greater,
            (_, -1) => Ordering::Less,
            (a, b) => java_difference(a, b).then_with(|| java_difference(self.start, other.start)),
        }
    }

    /// `SAMRecordCoordinateComparator.compare`.
    pub fn compare(&self, other: &Self) -> Ordering {
        let reverse = |key: &Self| key.flags & REVERSE != 0;
        self.file_order(other)
            .then_with(|| reverse(self).cmp(&reverse(other)))
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| self.flags.cmp(&other.flags))
            .then_with(|| self.mapping_quality.cmp(&other.mapping_quality))
            .then_with(|| self.mate_reference_id.cmp(&other.mate_reference_id))
            .then_with(|| self.mate_start.cmp(&other.mate_start))
            .then_with(|| self.template_length.cmp(&other.template_length))
    }
}

/// `SAMRecordCoordinateComparator.compare` on two records.
pub fn coordinate_order(a: &RawRecord, b: &RawRecord) -> Ordering {
    CoordinateKey::of(a).compare(&CoordinateKey::of(b))
}

/// `SAMRecord.getIntegerAttribute("HI")`: an integer as it is, a float
/// through `Number.longValue` (truncated), anything else the
/// `RuntimeException` htsjdk throws, naming the Java class; a value
/// outside `Integer` is an error too.
pub fn hit_index(record: &RawRecord) -> Result<Option<i64>> {
    let Some(value) = record.tag(*b"HI").map_err(Error::Record)? else {
        return Ok(None);
    };
    let (number, class) = match value {
        TagValue::Int(v) => (Some(v), ""),
        TagValue::Float(f) => (Some(f as i64), ""),
        TagValue::Char(_) => (None, "java.lang.Character"),
        TagValue::String(_) => (None, "java.lang.String"),
        TagValue::Hex(_) | TagValue::Array(b'c' | b'C', ..) => (None, "[B"),
        TagValue::Array(b's' | b'S', ..) => (None, "[S"),
        TagValue::Array(b'f', ..) => (None, "[F"),
        TagValue::Array(..) => (None, "[I"),
    };
    match number {
        Some(n) if i64::from(i32::MIN) <= n && n <= i64::from(i32::MAX) => Ok(Some(n)),
        Some(_) => Err(Error::Metric(
            "Value for tag HI is not in Integer range".to_owned(),
        )),
        None => Err(Error::Metric(format!(
            "Value for tag HI is not Number: class {class}"
        ))),
    }
}

/// What `SAMRecordQueryNameComparator` reads from a record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuerynameKey<'a> {
    pub name: Cow<'a, [u8]>,
    pub flags: u16,
    /// The `HI` tag, `None` without one.
    pub hit_index: Option<i64>,
}

impl<'a> QuerynameKey<'a> {
    /// The key from its parts, for a caller that holds them already.
    pub fn new(name: &'a [u8], flags: u16, hit_index: Option<i64>) -> Self {
        Self {
            name: Cow::Borrowed(name),
            flags,
            hit_index,
        }
    }

    /// The key of `record`, `HI` read as [`hit_index`] reads it.
    pub fn of(record: &'a RawRecord) -> Result<Self> {
        Ok(Self::new(record.name(), record.flags(), hit_index(record)?))
    }

    /// The key without its borrow of the record.
    pub fn into_owned(self) -> QuerynameKey<'static> {
        QuerynameKey {
            name: Cow::Owned(self.name.into_owned()),
            ..self
        }
    }

    /// `SAMRecordQueryNameComparator.compare`.
    pub fn compare(&self, other: &Self) -> Ordering {
        let by_name = self.name.cmp(&other.name);
        if by_name != Ordering::Equal {
            return by_name;
        }
        let (a, b) = (self.flags, other.flags);
        let has = |flags: u16, bit: u16| flags & bit != 0;
        if has(a, PAIRED) || has(b, PAIRED) {
            if !has(a, PAIRED) {
                return Ordering::Greater;
            }
            if !has(b, PAIRED) || (has(a, FIRST) && has(b, SECOND)) {
                return Ordering::Less;
            }
            if has(a, SECOND) && has(b, FIRST) {
                return Ordering::Greater;
            }
        }
        for bit in [REVERSE, SECONDARY, SUPPLEMENTARY] {
            if has(a, bit) != has(b, bit) {
                return if has(a, bit) {
                    Ordering::Greater
                } else {
                    Ordering::Less
                };
            }
        }
        match (self.hit_index, other.hit_index) {
            (Some(x), Some(y)) => x.cmp(&y),
            (Some(_), None) => Ordering::Greater,
            (None, Some(_)) => Ordering::Less,
            (None, None) => Ordering::Equal,
        }
    }
}

/// `SAMRecordQueryNameComparator.compare` on two records, for a caller
/// that cannot fail: an `HI` that is not an integer counts as absent,
/// where htsjdk throws ([`QuerynameKey::of`] gives that error).
pub fn queryname_order(a: &RawRecord, b: &RawRecord) -> Ordering {
    let hit = |record: &RawRecord| match record.tag(*b"HI") {
        Ok(Some(TagValue::Int(v))) => Some(v),
        _ => None,
    };
    QuerynameKey::new(a.name(), a.flags(), hit(a)).compare(&QuerynameKey::new(
        b.name(),
        b.flags(),
        hit(b),
    ))
}

/// The library IDs of `SAMRecordDuplicateComparator`: the header's
/// library names and `Unknown Library` in `TreeSet` order (by UTF-16 code
/// unit), numbered from 1.
#[derive(Clone, Debug, Default)]
pub struct Libraries {
    /// The library ID of each read group that names a library.
    by_read_group: HashMap<Vec<u8>, i32>,
    unknown: i32,
}

impl Libraries {
    pub fn new(header: &sam::Header) -> Self {
        let library = |map: &sam::header::record::value::Map<_>| {
            map.other_fields()
                .get(&LIBRARY)
                .map(|name| String::from_utf8_lossy(name).into_owned())
        };
        let mut names: Vec<String> = header.read_groups().values().filter_map(library).collect();
        names.push(UNKNOWN_LIBRARY.to_owned());
        names.sort_by(|a, b| a.encode_utf16().cmp(b.encode_utf16()));
        names.dedup();
        let id_of = |name: &str| {
            names
                .binary_search_by(|probe| probe.encode_utf16().cmp(name.encode_utf16()))
                .map_or(0, |index| index as i32 + 1)
        };
        Self {
            by_read_group: header
                .read_groups()
                .iter()
                .filter_map(|(id, map)| Some((id.to_vec(), id_of(&library(map)?))))
                .collect(),
            unknown: id_of(UNKNOWN_LIBRARY),
        }
    }

    /// `getLibraryId`: the library of the record's read group, or the ID
    /// of `Unknown Library`.
    fn id(&self, read_group: Option<&[u8]>) -> i32 {
        read_group
            .and_then(|id| self.by_read_group.get(id))
            .copied()
            .unwrap_or(self.unknown)
    }
}

/// `FF`, `FR`, `F`, `RF`, `RR`, `R` of `SAMRecordDuplicateComparator`.
fn orientation(flags: u16, both_mapped: bool) -> i32 {
    let reverse = flags & REVERSE != 0;
    if !both_mapped {
        return if reverse { 5 } else { 2 };
    }
    match (reverse, flags & MATE_REVERSE != 0) {
        (false, false) => 0,
        (false, true) => 1,
        (true, false) => 3,
        (true, true) => 4,
    }
}

/// What `SAMRecordDuplicateComparator` reads from a record.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DuplicateKey<'a> {
    pub library: i32,
    pub reference_id: i32,
    /// The unclipped end of a reverse-strand record, else the unclipped
    /// start.
    pub read_coordinate: i32,
    pub flags: u16,
    /// The mate's reference index and unclipped 5' coordinate; -1 unless
    /// the record is paired with both ends mapped.
    pub mate_reference_id: i32,
    pub mate_coordinate: i32,
    /// `DuplicateScoringStrategy.computeDuplicateScore` with
    /// `TOTAL_MAPPED_REFERENCE_LENGTH` and the mate's CIGAR.
    pub score: i16,
    /// The `RG` tag: the read name's prefix in `getCanonicalRecordName`.
    pub read_group: Option<Vec<u8>>,
    pub name: Cow<'a, [u8]>,
}

impl<'a> DuplicateKey<'a> {
    /// The key of `record`. A paired record whose mate is mapped needs an
    /// `MC` tag, as in htsjdk.
    pub fn of(record: &'a RawRecord, libraries: &Libraries) -> Result<Self> {
        let flags = record.flags();
        let name = || String::from_utf8_lossy(record.name()).into_owned();
        let string_tag = |tag: [u8; 2]| match record.tag(tag).map_err(Error::Record)? {
            None => Ok(None),
            Some(TagValue::String(text)) => Ok(Some(text)),
            Some(_) => Err(Error::Metric(format!(
                "{}: Value for tag {} is not a String",
                name(),
                String::from_utf8_lossy(&tag)
            ))),
        };
        let read_group = string_tag(*b"RG")?;
        let cigar: Vec<(u32, u8)> = record.cigar().collect();
        let unmapped = flags & UNMAPPED != 0;
        let reverse = flags & REVERSE != 0;
        let (start, end) = unclipped(record.alignment_start(), &cigar);
        let read_coordinate = match (reverse, unmapped) {
            (false, _) => start,
            // `getAlignmentEnd` of an unmapped record is 0, so only the
            // trailing clips count.
            (true, true) => end
                .wrapping_sub(record.alignment_start())
                .wrapping_sub(record.reference_length())
                .wrapping_add(1),
            (true, false) => end,
        };
        let cap = i32::from(i16::MAX / 2);
        let mut score: i16 = if unmapped {
            0
        } else {
            record.reference_length().min(cap) as i16
        };
        let mate_mapped = flags & PAIRED != 0 && flags & MATE_UNMAPPED == 0;
        let mut mate_reference_id = -1;
        let mut mate_coordinate = -1;
        if mate_mapped {
            let mate_cigar = string_tag(*b"MC")?.ok_or_else(|| {
                Error::Metric(format!("Mate CIGAR (Tag MC) not found: {}", name()))
            })?;
            let ops = parse_cigar(&mate_cigar).map_err(Error::Metric)?;
            let reference_length = crate::bam_raw::reference_length(ops.iter().copied());
            score = score.wrapping_add(reference_length.min(cap) as i16);
            if !unmapped {
                let (mate_start, mate_end) = unclipped(record.mate_alignment_start(), &ops);
                mate_reference_id = record.mate_reference_id();
                mate_coordinate = if flags & MATE_REVERSE != 0 {
                    mate_end
                } else {
                    mate_start
                };
            }
        }
        if flags & QC_FAIL != 0 {
            score = score.wrapping_add(i16::MIN / 2);
        }
        Ok(Self {
            library: libraries.id(read_group.as_deref()),
            reference_id: record.reference_id(),
            read_coordinate,
            flags,
            mate_reference_id,
            mate_coordinate,
            score,
            read_group,
            name: Cow::Borrowed(record.name()),
        })
    }

    /// The key without its borrow of the record.
    pub fn into_owned(self) -> DuplicateKey<'static> {
        DuplicateKey {
            name: Cow::Owned(self.name.into_owned()),
            ..self
        }
    }

    fn has(&self, flag: u16) -> bool {
        self.flags & flag != 0
    }

    /// `pairedEndAndBothMapped`.
    fn both_mapped(&self) -> bool {
        self.has(PAIRED) && !self.has(UNMAPPED) && !self.has(MATE_UNMAPPED)
    }

    /// `hasMappedEnd`.
    fn mapped_end(&self) -> bool {
        !self.has(UNMAPPED) || (self.has(PAIRED) && !self.has(MATE_UNMAPPED))
    }

    /// `hasUnmappedEnd`.
    fn unmapped_end(&self) -> bool {
        self.has(UNMAPPED) || (self.has(PAIRED) && self.has(MATE_UNMAPPED))
    }

    /// `SAMUtils.getCanonicalRecordName`: `RG:name`, or the name alone.
    fn canonical(&self) -> impl Iterator<Item = &u8> {
        let prefix = self
            .read_group
            .iter()
            .flat_map(|group| group.iter().chain(b":"));
        prefix.chain(self.name.iter())
    }

    /// `fileOrderCompare(a, b, false, true)`.
    pub fn file_order(&self, other: &Self) -> Ordering {
        let by_reference = match (self.reference_id, other.reference_id) {
            (-1, -1) => Ordering::Equal,
            (-1, _) => Ordering::Greater,
            (_, -1) => Ordering::Less,
            (a, b) => java_difference(a, b),
        };
        let mut order = java_difference(self.library, other.library)
            .then(by_reference)
            .then_with(|| java_difference(self.read_coordinate, other.read_coordinate))
            .then_with(|| {
                java_difference(
                    orientation(self.flags, self.both_mapped()),
                    orientation(other.flags, other.both_mapped()),
                )
            });
        if self.both_mapped() && other.both_mapped() {
            order = order
                .then_with(|| java_difference(self.mate_reference_id, other.mate_reference_id))
                .then_with(|| java_difference(self.mate_coordinate, other.mate_coordinate));
        }
        order
            .then_with(|| other.mapped_end().cmp(&self.mapped_end()))
            .then_with(|| {
                if self.has(PAIRED) == other.has(PAIRED) {
                    self.unmapped_end().cmp(&other.unmapped_end())
                } else if self.has(PAIRED) {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            })
    }

    /// `SAMRecordDuplicateComparator.compare`.
    pub fn compare(&self, other: &Self) -> Ordering {
        self.file_order(other)
            .then_with(|| {
                // `DuplicateScoringStrategy.compare`: paired first, then
                // the higher score, then the canonical name.
                other
                    .has(PAIRED)
                    .cmp(&self.has(PAIRED))
                    .then_with(|| other.score.cmp(&self.score))
                    .then_with(|| self.canonical().cmp(other.canonical()))
            })
            .then_with(|| self.name.cmp(&other.name))
            .then_with(|| {
                if self.has(PAIRED) && other.has(PAIRED) {
                    other.has(FIRST).cmp(&self.has(FIRST))
                } else {
                    Ordering::Equal
                }
            })
    }
}

#[cfg(test)]
mod tests;

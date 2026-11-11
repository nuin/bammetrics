//! `MarkDuplicates.buildSortedReadEndLists`: the first read of the inputs.
//!
//! Rules copied from Picard 3.5.0:
//! * Every record before the stop counts toward the `PG` values seen (for
//!   the `@PG` records written later). In coordinate order the read stops
//!   at the first unmapped record without a reference; that record's `PG`
//!   still counts, the ones after do not.
//! * Unmapped, secondary and supplementary records make no `ReadEnds`.
//!   Every other record makes a fragment: its unclipped 5' end (the
//!   unclipped end for a reverse-strand read), strand, score, library and
//!   read-name location, and, when its mate is mapped, the mate's
//!   reference as `read2` (which marks it paired).
//! * A record whose mate is mapped also waits under its `RG` ID plus read
//!   name (no `RG`, or one the header lacks, is Picard's
//!   `NullPointerException`) for its mate ([`super::mates`]); the second
//!   to arrive completes the pair. The lower end, by reference then coordinate, becomes
//!   `read1`; equal coordinates with opposite strands are `FR` whatever the
//!   file order. The pair keeps the first mate's library and location, its
//!   score is the sum of both scores in 16-bit arithmetic, and its optical
//!   orientation takes the first-of-pair strand first.
//! * In query-name order each record's index is the index of the first
//!   record with its name, so every record of a template shares it.
//! * Scores (`DuplicateScoringStrategy.computeDuplicateScore`): the sum of
//!   base qualities of at least 15, or the reference length of the CIGAR,
//!   each capped at 16383; `RANDOM` is 8192 plus the low 14 bits of the
//!   htsjdk `Murmur3` (seed 1) hash of the read name. A QC-failed record
//!   loses 16384, in 16-bit arithmetic.

use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use noodles::sam;

use super::inputs::{self, SortOrder};
use super::mates::Mates;
use super::read_ends::{F, FR, R, RF, ReadEnds, orientation_byte};
use super::{Options, Scoring};
use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::optical::Location;
use crate::spill::{Sorter, TempDir};

/// `LibraryIdGenerator.UNKNOWN_LIBRARY`.
pub const UNKNOWN_LIBRARY: &str = "Unknown Library";

/// What the first read leaves for [`super::mark`] and [`super::write`].
pub struct Built<'d> {
    pub pairs: Sorter<'d, ReadEnds>,
    pub fragments: Sorter<'d, ReadEnds>,
    /// Library ID of each library name, assigned in order of first use.
    pub library_ids: HashMap<String, i16>,
    /// Each distinct `PG` value seen, in first-seen order (`None` for a
    /// record without one).
    pub pg_ids: Vec<Option<Vec<u8>>>,
}

/// Read groups of the header: index in `@RG` order and library.
pub struct ReadGroups {
    groups: HashMap<Vec<u8>, (i16, String)>,
    count: i16,
}

impl ReadGroups {
    pub fn new(header: &sam::Header) -> Self {
        let mut groups = HashMap::new();
        for (i, (id, map)) in header.read_groups().iter().enumerate() {
            let library = map
                .other_fields()
                .get(&sam::header::record::value::map::read_group::tag::LIBRARY)
                .map_or_else(|| UNKNOWN_LIBRARY.to_owned(), |lb| latin1(lb));
            groups.entry(id.to_vec()).or_insert((i as i16, library));
        }
        Self {
            groups,
            count: header.read_groups().len() as i16,
        }
    }

    /// `LibraryIdGenerator.getLibraryName`.
    pub fn library(&self, rg: Option<&[u8]>) -> &str {
        rg.and_then(|id| self.groups.get(id))
            .map_or(UNKNOWN_LIBRARY, |(_, library)| library.as_str())
    }

    /// The read group number `buildReadEnds` stores: the `@RG` index, the
    /// group count for an ID the header lacks, 0 without `RG`.
    fn number(&self, rg: Option<&[u8]>) -> i16 {
        match rg {
            None => 0,
            Some(id) => self.groups.get(id).map_or(self.count, |(i, _)| *i),
        }
    }

    fn contains(&self, rg: &[u8]) -> bool {
        self.groups.contains_key(rg)
    }
}

/// A Java `String` made from BAM bytes, one char per byte.
pub fn latin1(bytes: &[u8]) -> String {
    bytes.iter().map(|&b| char::from(b)).collect()
}

/// `String.hashCode` of a string decoded one char per byte.
pub fn java_hash(bytes: &[u8]) -> i32 {
    bytes
        .iter()
        .fold(0i32, |h, &b| h.wrapping_mul(31).wrapping_add(i32::from(b)))
}

/// The string value of `tag`: `None` when absent, an error for another
/// type (htsjdk `getStringAttribute`, or the cast in `getLibraryName`).
pub fn string_tag(tags: &TagCache, tag: [u8; 2]) -> Result<Option<&[u8]>> {
    match tags.get(tag) {
        None => Ok(None),
        Some(TagValue::String(text)) => Ok(Some(text)),
        Some(_) => Err(Error::Metric(format!(
            "Value for tag {} is not a String",
            String::from_utf8_lossy(&tag)
        ))),
    }
}

/// A record's tags, decoded once.
#[derive(Default)]
pub struct TagCache {
    fields: Vec<([u8; 2], TagValue)>,
}

impl TagCache {
    pub fn load(&mut self, record: &RawRecord) -> Result<()> {
        self.fields.clear();
        for field in crate::bam_raw::TagIter::new(record.tag_bytes()) {
            let (name, value) = field.map_err(Error::Record)?;
            match self.fields.iter_mut().find(|(n, _)| *n == name) {
                Some(slot) => slot.1 = value,
                None => self.fields.push((name, value)),
            }
        }
        Ok(())
    }

    pub fn get(&self, tag: [u8; 2]) -> Option<&TagValue> {
        self.fields.iter().find(|(n, _)| *n == tag).map(|(_, v)| v)
    }
}

/// A two-letter tag name (`SAMTag.makeBinaryTag` fails on any other length).
pub fn tag_name(name: &str) -> Result<[u8; 2]> {
    match name.as_bytes() {
        [a, b] => Ok([*a, *b]),
        _ => Err(Error::Metric(format!(
            "String tag does not have length() == 2: {name}"
        ))),
    }
}

/// `EstimateLibraryComplexity.getReadBarcodeValue`.
pub fn barcode_value(tags: &TagCache, tag: Option<&str>) -> Result<i32> {
    let Some(tag) = tag else { return Ok(0) };
    Ok(string_tag(tags, tag_name(tag)?)?.map_or(0, java_hash))
}

/// Unclipped start and end of a mapped record (`SAMUtils.getUnclippedStart`,
/// `getUnclippedEnd`), 1-based.
pub fn unclipped_ends(record: &RawRecord) -> (i32, i32) {
    let ops: Vec<(u32, u8)> = record.cigar().collect();
    unclipped(record.alignment_start(), &ops)
}

/// [`unclipped_ends`] for a start and a CIGAR as `(length, BAM op)`.
pub fn unclipped(start: i32, ops: &[(u32, u8)]) -> (i32, i32) {
    let clip = |op: u8| op == 4 || op == 5;
    let reference_length: i64 = ops
        .iter()
        .filter(|(_, op)| matches!(op, 0 | 2 | 3 | 7 | 8))
        .map(|(len, _)| i64::from(*len))
        .sum();
    let end = (i64::from(start) + reference_length - 1) as i32;
    let leading: i64 = ops
        .iter()
        .take_while(|(_, op)| clip(*op))
        .map(|(len, _)| i64::from(*len))
        .sum();
    let trailing: i64 = ops
        .iter()
        .rev()
        .take_while(|(_, op)| clip(*op))
        .map(|(len, _)| i64::from(*len))
        .sum();
    (
        (i64::from(start) - leading) as i32,
        (i64::from(end) + trailing) as i32,
    )
}

/// htsjdk `Murmur3.hashUnencodedChars` with seed 1, over one char per byte.
fn murmur3(name: &[u8]) -> i32 {
    fn mix_k1(k1: i32) -> i32 {
        k1.wrapping_mul(0xcc9e_2d51_u32 as i32)
            .rotate_left(15)
            .wrapping_mul(0x1b87_3593)
    }
    fn mix_h1(h1: i32, k1: i32) -> i32 {
        (h1 ^ k1)
            .rotate_left(13)
            .wrapping_mul(5)
            .wrapping_add(0xe654_6b64_u32 as i32)
    }
    let mut h1 = 1i32;
    let mut i = 1;
    while i < name.len() {
        let k1 = i32::from(name[i - 1]) | (i32::from(name[i]) << 16);
        h1 = mix_h1(h1, mix_k1(k1));
        i += 2;
    }
    if name.len() % 2 == 1 {
        h1 ^= mix_k1(i32::from(name[name.len() - 1]));
    }
    h1 ^= (2 * name.len()) as i32;
    h1 ^= ((h1 as u32) >> 16) as i32;
    h1 = h1.wrapping_mul(0x85eb_ca6b_u32 as i32);
    h1 ^= ((h1 as u32) >> 13) as i32;
    h1 = h1.wrapping_mul(0xc2b2_ae35_u32 as i32);
    h1 ^ ((h1 as u32) >> 16) as i32
}

/// `DuplicateScoringStrategy.computeDuplicateScore`.
pub fn score(record: &RawRecord, scoring: Scoring) -> i16 {
    let cap = i64::from(i16::MAX / 2);
    let mut score: i16 = match scoring {
        Scoring::SumOfBaseQualities => {
            let quals = record.qualities();
            // htsjdk reads a quality string starting with 0xFF as missing.
            let sum: i64 = if quals.first() == Some(&0xff) {
                0
            } else {
                quals
                    .iter()
                    .map(|&q| i64::from(q as i8))
                    .filter(|&q| q >= 15)
                    .sum()
            };
            sum.min(cap) as i16
        }
        Scoring::TotalMappedReferenceLength => {
            if record.flags() & 0x4 != 0 {
                0
            } else {
                let length: i64 = record
                    .cigar()
                    .filter(|(_, op)| matches!(op, 0 | 2 | 3 | 7 | 8))
                    .map(|(len, _)| i64::from(len))
                    .sum();
                length.min(cap) as i16
            }
        }
        Scoring::Random => {
            let low = (murmur3(record.name()) & 0b11_1111_1111_1111) as i16;
            low.wrapping_add(8192)
        }
    };
    if record.flags() & 0x200 != 0 {
        score = score.wrapping_add(i16::MIN / 2);
    }
    score
}

/// `UmiUtil.getStrand`: top, bottom, or unknown (`None`).
pub fn strand(record: &RawRecord, tags: &TagCache) -> Result<Option<bool>> {
    let flags = record.flags();
    if flags & 0x1 == 0 {
        return Err(Error::Metric(format!(
            "Inappropriate call if not paired read: {}",
            String::from_utf8_lossy(record.name())
        )));
    }
    if flags & 0x4 != 0 || flags & 0x8 != 0 {
        return Ok(None);
    }
    let first = flags & 0x40 != 0;
    let (reference, mate_reference) = (record.reference_id(), record.mate_reference_id());
    if reference != mate_reference {
        return Ok(Some(first == (reference < mate_reference)));
    }
    let (start, end) = unclipped_ends(record);
    let read5 = if flags & 0x10 != 0 { end } else { start };
    let mate5 = mate_unclipped(record, tags, flags & 0x20 != 0)?;
    Ok(Some(first == (read5 <= mate5)))
}

/// `SAMUtils.getMateUnclippedEnd` or `getMateUnclippedStart` from `MC`.
fn mate_unclipped(record: &RawRecord, tags: &TagCache, end: bool) -> Result<i32> {
    let missing = || {
        Error::Metric(format!(
            "Mate CIGAR (Tag MC) not found: {}",
            String::from_utf8_lossy(record.name())
        ))
    };
    let Some(TagValue::String(text)) = tags.get(*b"MC") else {
        return Err(missing());
    };
    let ops = parse_cigar(text).ok_or_else(missing)?;
    let (start, stop) = unclipped(record.mate_alignment_start(), &ops);
    Ok(if end { stop } else { start })
}

/// A CIGAR string as `(length, BAM op)`.
pub fn parse_cigar(text: &[u8]) -> Option<Vec<(u32, u8)>> {
    if text == b"*" {
        return Some(Vec::new());
    }
    let mut ops = Vec::new();
    let mut length: u32 = 0;
    for &c in text {
        if c.is_ascii_digit() {
            length = length.checked_mul(10)?.checked_add(u32::from(c - b'0'))?;
        } else {
            let op = b"MIDNSHP=X".iter().position(|&o| o == c)? as u8;
            ops.push((length, op));
            length = 0;
        }
    }
    Some(ops)
}

/// Java `String.split("-")`: trailing empty strings dropped.
fn java_split_dash(text: &[u8]) -> Vec<&[u8]> {
    if text.is_empty() {
        return vec![text];
    }
    let mut parts: Vec<&[u8]> = text.split(|&b| b == b'-').collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// `UmiUtil.getTopStrandNormalizedUmi`, as the `Objects.hash` Picard keeps.
fn top_strand_umi_hash(record: &RawRecord, tags: &TagCache, opts: &Options) -> Result<i32> {
    let Some(tag) = opts.barcode_tag.as_deref() else {
        return Ok(31);
    };
    let Some(umi) = string_tag(tags, tag_name(tag)?)? else {
        return Ok(31);
    };
    if !umi.iter().all(|b| b"ATCGNatcgn-".contains(b)) {
        return Err(Error::Metric(
            "UMI found with illegal characters.  UMIs must match the regular expression \
             ^[ATCGNatcgn-]*$."
                .to_owned(),
        ));
    }
    if !opts.duplex_umi {
        return Ok(31i32.wrapping_add(java_hash(umi)));
    }
    let parts = java_split_dash(umi);
    if parts.len() != 2 {
        return Err(Error::Metric(format!(
            "Duplex UMIs must be of the form X-Y where X and Y are equal length UMIs, for \
             example AT-GA.  Found UMI, {}",
            String::from_utf8_lossy(umi)
        )));
    }
    let normalized = match strand(record, tags)? {
        Some(false) => [parts[1], b"-", parts[0]].concat(),
        _ => umi.to_vec(),
    };
    Ok(31i32.wrapping_add(java_hash(&normalized)))
}

/// The first read: every record of the input in order through
/// [`Builder::push`], then [`Builder::finish`]. The tool reads its inputs
/// into it; in a chain ([`crate::stage`]) the step before pushes its
/// batches.
pub struct Builder<'d> {
    opts: Options,
    groups: ReadGroups,
    order: SortOrder,
    queryname: bool,
    library_ids: HashMap<String, i16>,
    next_library: i16,
    pending: Mates<'d>,
    pairs: Sorter<'d, ReadEnds>,
    fragments: Sorter<'d, ReadEnds>,
    pg_seen: HashSet<Option<Vec<u8>>>,
    pg_ids: Vec<Option<Vec<u8>>>,
    tags: TagCache,
    /// Index of the next record in the input.
    index: i64,
    /// In query-name order, the index of the first record with the
    /// current name.
    name_index: i64,
    last_name: Option<Vec<u8>>,
    /// True after the first unplaced record of a coordinate-sorted input:
    /// the records after it are not read.
    stopped: bool,
}

pub fn run<'d>(opts: &Options, inputs: &[PathBuf], temp: &'d TempDir) -> Result<Built<'d>> {
    inputs::open(
        inputs,
        opts.reference.as_deref(),
        opts.threads,
        opts.stringency,
        opts.assume_sort_order,
        |mut stream| {
            let mut builder = Builder::new(opts, &stream.header, temp)?;
            let mut record = RawRecord::default();
            while !builder.stopped && stream.next(&mut record)? {
                builder.push(&record)?;
            }
            Ok(builder.finish())
        },
    )
}

/// Picard's check that the input is coordinate or query-name sorted
/// (it runs after the first read; the outcome is the same).
pub fn check_order(order: SortOrder, assume: Option<SortOrder>) -> Result<()> {
    if matches!(order, SortOrder::Coordinate | SortOrder::Queryname) {
        return Ok(());
    }
    Err(Error::Metric(format!(
        "This program requires input that are either coordinate or query sorted (according to \
         the header, or at least ASSUME_SORT_ORDER and the content.) Found ASSUME_SORT_ORDER={} \
         and header sortorder={}",
        assume.map_or("null", SortOrder::name),
        order.name()
    )))
}

impl<'d> Builder<'d> {
    /// The first read of an input whose header (with `ASSUME_SORT_ORDER`
    /// applied) is `header`; spilled runs go to `temp`.
    pub fn new(opts: &Options, header: &sam::Header, temp: &'d TempDir) -> Result<Self> {
        let order = SortOrder::of(header);
        check_order(order, opts.assume_sort_order)?;
        let queryname = order == SortOrder::Queryname;
        Ok(Self {
            opts: opts.clone(),
            groups: ReadGroups::new(header),
            order,
            queryname,
            library_ids: HashMap::new(),
            next_library: 1,
            pending: Mates::new(queryname, temp, opts.records_in_ram),
            pairs: Sorter::new(temp, "pairs", opts.records_in_ram),
            fragments: Sorter::new(temp, "fragments", opts.records_in_ram),
            pg_seen: HashSet::new(),
            pg_ids: Vec::new(),
            tags: TagCache::default(),
            index: 0,
            name_index: i64::MAX,
            last_name: None,
            stopped: false,
        })
    }

    /// The next record of the input. After the first unplaced record of a
    /// coordinate-sorted input every further record is ignored.
    pub fn push(&mut self, record: &RawRecord) -> Result<()> {
        if self.stopped {
            return Ok(());
        }
        self.tags.load(record)?;
        if self.opts.program_record_id.is_some() {
            let pg = string_tag(&self.tags, *b"PG")?.map(<[u8]>::to_vec);
            if self.pg_seen.insert(pg.clone()) {
                self.pg_ids.push(pg);
            }
        }
        if self.queryname && self.last_name.as_deref() != Some(record.name()) {
            self.last_name = Some(record.name().to_vec());
            self.name_index = self.index;
        }
        let flags = record.flags();
        if flags & 0x4 != 0 {
            if record.reference_id() == -1 && self.order == SortOrder::Coordinate {
                self.stopped = true;
                return Ok(());
            }
        } else if flags & 0x900 == 0 {
            let index_for_read = if self.queryname {
                self.name_index
            } else {
                self.index
            };
            self.accept(record, index_for_read)?;
        }
        self.index += 1;
        Ok(())
    }

    /// The sorted lists and the IDs the later steps need.
    pub fn finish(self) -> Built<'d> {
        Built {
            pairs: self.pairs,
            fragments: self.fragments,
            library_ids: self.library_ids,
            pg_ids: self.pg_ids,
        }
    }

    fn library_id(&mut self, rg: Option<&[u8]>) -> i16 {
        let library = self.groups.library(rg).to_owned();
        if let Some(id) = self.library_ids.get(&library) {
            return *id;
        }
        let id = self.next_library;
        self.next_library = self.next_library.wrapping_add(1);
        self.library_ids.insert(library, id);
        id
    }

    /// `buildReadEnds`.
    fn read_ends(&mut self, record: &RawRecord, index: i64) -> Result<ReadEnds> {
        let flags = record.flags();
        let negative = flags & 0x10 != 0;
        let (start, end) = unclipped_ends(record);
        let mut ends = ReadEnds {
            read1_reference: record.reference_id(),
            read1_coordinate: if negative { end } else { start },
            orientation: if negative { R } else { F },
            read1_index: index,
            score: score(record, self.opts.scoring),
            ..ReadEnds::default()
        };
        if flags & 0x1 != 0 && flags & 0x8 == 0 {
            ends.read2_reference = record.mate_reference_id();
        }
        let rg = string_tag(&self.tags, *b"RG")?.map(<[u8]>::to_vec);
        ends.library_id = self.library_id(rg.as_deref());
        let mut location = Location::default();
        let name = latin1(record.name());
        if self.opts.finder.parser.read_location(&name, &mut location) {
            ends.read_group = self.groups.number(rg.as_deref());
        }
        ends.tile = location.tile;
        ends.x = location.x;
        ends.y = location.y;
        if self.opts.use_barcodes() {
            ends.barcode = top_strand_umi_hash(record, &self.tags, &self.opts)?;
            if flags & 0x1 == 0 || flags & 0x40 != 0 {
                ends.read_one_barcode =
                    barcode_value(&self.tags, self.opts.read_one_barcode_tag.as_deref())?;
            } else {
                ends.read_two_barcode =
                    barcode_value(&self.tags, self.opts.read_two_barcode_tag.as_deref())?;
            }
        }
        Ok(ends)
    }

    /// One mapped primary record.
    fn accept(&mut self, record: &RawRecord, index: i64) -> Result<()> {
        let fragment = self.read_ends(record, index)?;
        let flags = record.flags();
        let paired = flags & 0x1 != 0 && flags & 0x8 == 0;
        if !paired {
            return self.fragments.push(fragment);
        }
        let rg = string_tag(&self.tags, *b"RG")?
            .filter(|rg| self.groups.contains(rg))
            .ok_or_else(|| {
                Error::Metric(format!(
                    "Cannot invoke \"htsjdk.samtools.SAMReadGroupRecord.getReadGroupId()\" \
                     because the return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" \
                     is null (read {})",
                    String::from_utf8_lossy(record.name())
                ))
            })?;
        let mut key = rg.to_vec();
        key.extend_from_slice(record.name());
        let Some(mut pair) = self.pending.remove(record.reference_id(), &key)? else {
            self.pending
                .put(fragment.read2_reference, key, fragment.clone())?;
            return self.fragments.push(fragment);
        };
        let negative = flags & 0x10 != 0;
        let mate_reverse = pair.orientation == R;
        if flags & 0x40 != 0 {
            pair.orientation_for_optical = orientation_byte(negative, mate_reverse);
            if self.opts.use_barcodes() {
                pair.read_one_barcode =
                    barcode_value(&self.tags, self.opts.read_one_barcode_tag.as_deref())?;
            }
        } else {
            pair.orientation_for_optical = orientation_byte(mate_reverse, negative);
            if self.opts.use_barcodes() {
                pair.read_two_barcode =
                    barcode_value(&self.tags, self.opts.read_two_barcode_tag.as_deref())?;
            }
        }
        let (mate_reference, mate_coordinate) =
            (fragment.read1_reference, fragment.read1_coordinate);
        if mate_reference > pair.read1_reference
            || (mate_reference == pair.read1_reference && mate_coordinate >= pair.read1_coordinate)
        {
            pair.read2_reference = mate_reference;
            pair.read2_coordinate = mate_coordinate;
            pair.read2_index = index;
            pair.orientation = orientation_byte(mate_reverse, negative);
            if pair.read2_reference == pair.read1_reference
                && pair.read2_coordinate == pair.read1_coordinate
                && pair.orientation == RF
            {
                pair.orientation = FR;
            }
        } else {
            pair.read2_reference = pair.read1_reference;
            pair.read2_coordinate = pair.read1_coordinate;
            pair.read2_index = pair.read1_index;
            pair.read1_reference = mate_reference;
            pair.read1_coordinate = mate_coordinate;
            pair.read1_index = index;
            pair.orientation = orientation_byte(negative, mate_reverse);
        }
        pair.score = pair.score.wrapping_add(fragment.score);
        self.fragments.push(fragment)?;
        self.pairs.push(pair)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_raw::build;

    #[test]
    fn scores_follow_the_strategy_and_the_qc_discount() {
        let quals = [10u8, 15, 30, 40];
        let record = build("r", 0, 0, 99, &[(4, 0)], &quals, b"");
        assert_eq!(score(&record, Scoring::SumOfBaseQualities), 85);
        assert_eq!(score(&record, Scoring::TotalMappedReferenceLength), 4);
        let failed = build("r", 0x200, 0, 99, &[(4, 0)], &quals, b"");
        assert_eq!(score(&failed, Scoring::SumOfBaseQualities), 85 - 16384);
        let missing = build("r", 0, 0, 99, &[(2, 0)], &[0xff, 0xff], b"");
        assert_eq!(score(&missing, Scoring::SumOfBaseQualities), 0);
        let random = score(&record, Scoring::Random);
        assert!((8192..8192 + 16384).contains(&i32::from(random)));
    }

    #[test]
    fn murmur3_matches_guava() {
        // Guava Hashing.murmur3_32(1).hashUnencodedChars, as htsjdk computes it.
        assert_eq!(murmur3(b""), 0x514e_28b7);
    }

    #[test]
    fn unclipped_ends_add_the_clips() {
        // 5S10M2D3M4H at 100: aligned 100..114, unclipped 95..118.
        let ops = [(5, 4), (10, 0), (2, 2), (3, 0), (4, 5)];
        assert_eq!(unclipped(100, &ops), (95, 118));
        assert_eq!(parse_cigar(b"5S10M2D3M4H"), Some(ops.to_vec()));
        assert_eq!(parse_cigar(b"5Q"), None);
    }

    #[test]
    fn dashes_split_like_java() {
        assert_eq!(java_split_dash(b"AT-GA"), [b"AT".as_slice(), b"GA"]);
        assert_eq!(java_split_dash(b"A-T-").len(), 2);
        assert_eq!(java_split_dash(b"-").len(), 0);
        assert_eq!(java_split_dash(b"").len(), 1);
        assert_eq!(java_split_dash(b"A--T").len(), 3);
    }
}

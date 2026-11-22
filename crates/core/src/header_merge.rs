//! `SamFileHeaderMerger` and `MergingSamRecordIterator` of htsjdk 5.0.0,
//! shared by `MarkDuplicates` (several `INPUT`s) and the tools that merge
//! inputs: the merged header, the read group and program IDs renamed on a
//! collision, and the records of several sorted inputs taken in coordinate
//! or query-name order.
//!
//! * The sequence dictionaries must be equal (`SequenceUtil.
//!   assertSequenceDictionariesEqual`: names or `AN` names, lengths, `M5`
//!   when both have one), or with `MERGE_SEQUENCE_DICTIONARIES` they are
//!   merged (`mergeSequenceDictionaries`, [`merge_dictionaries`]) and each
//!   record's reference and mate reference are found again by name. The
//!   merged header has `@HD VN:1.6 GO:none SO:<order>`, the dictionary,
//!   the read groups and programs of every input sorted by ID, and every
//!   comment in input order.
//! * A read group or program ID used by two inputs for different records
//!   gets a suffix `.<base 36 counter>`, and the records of that input are
//!   rewritten to the new ID (an `RG` or `PG` the input's header lacks is
//!   removed, as htsjdk sets it to null).
//! * Records are taken in coordinate or query-name order through a port of
//!   `java.util.PriorityQueue`, so equal records from two inputs come in the
//!   order Picard gives them. With a merged dictionary the coordinate order
//!   compares merged reference indexes, and an unplaced record with a placed
//!   mate takes its mate's (`MergedSequenceDictionaryCoordinateOrderComparator`).
//!   [`MergeOrder::InputOrder`] takes every record of one input before the
//!   next: htsjdk's `unsorted` merge compares `System.identityHashCode`, so
//!   its order between inputs changes from run to run.
//!
//! Also here: [`SortOrder`] (`SAMFileHeader.SortOrder`) and
//! [`set_header_tag`] (`setSortOrder`, `setGroupOrder`).

use std::cmp::Ordering;
use std::collections::hash_map::DefaultHasher;
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::io::Read;

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::Version;
use noodles::sam::header::record::value::map::header::tag::{GROUP_ORDER, SORT_ORDER};
use noodles::sam::header::record::value::map::reference_sequence::tag::{
    ALTERNATIVE_NAMES, MD5_CHECKSUM,
};
use noodles::sam::header::record::value::map::{Inner, Program, ReadGroup, ReferenceSequence};

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::sam_order::{CoordinateKey, coordinate_order, queryname_order};

mod dictionary;
mod programs;

pub use dictionary::merge_dictionaries;
use programs::Left;

/// `SAMFileHeader.SortOrder`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    pub fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "unsorted" => Self::Unsorted,
            "queryname" => Self::Queryname,
            "coordinate" => Self::Coordinate,
            "duplicate" => Self::Duplicate,
            "unknown" => Self::Unknown,
            _ => return None,
        })
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Unsorted => "unsorted",
            Self::Queryname => "queryname",
            Self::Coordinate => "coordinate",
            Self::Duplicate => "duplicate",
            Self::Unknown => "unknown",
        }
    }

    /// `SAMFileHeader.getSortOrder`: `unsorted` without `SO`, `unknown` for
    /// a value htsjdk does not know.
    pub fn of(header: &sam::Header) -> Self {
        match header
            .header()
            .and_then(|h| h.other_fields().get(&SORT_ORDER))
        {
            None => Self::Unsorted,
            Some(text) => Self::parse(&text.to_string()).unwrap_or(Self::Unknown),
        }
    }
}

/// `SAMFileHeader.setSortOrder`: an `@HD` line is made (`VN:1.6`) when
/// there is none.
pub fn set_header_tag(header: &mut sam::Header, tag: HeaderTag, value: &str) {
    let map = header.header_mut().get_or_insert_with(|| {
        Map::<sam::header::record::value::map::Header>::new(Version::new(1, 6))
    });
    let key = match tag {
        HeaderTag::SortOrder => SORT_ORDER,
        HeaderTag::GroupOrder => GROUP_ORDER,
    };
    map.other_fields_mut().insert(key, value.into());
}

#[derive(Clone, Copy, Debug)]
pub enum HeaderTag {
    SortOrder,
    GroupOrder,
}

/// `SamFileHeaderMerger.positiveFourDigitBase36Str`.
pub fn base36(mut value: i32) -> String {
    if value == 0 {
        return "0".to_owned();
    }
    let digits = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    let mut out = Vec::new();
    while value > 0 {
        out.push(digits[(value % 36) as usize]);
        value /= 36;
    }
    out.reverse();
    String::from_utf8(out).expect("ASCII")
}

/// A read group or program ID.
type Id = Vec<u8>;

/// The distinct records under one ID, each with the inputs that hold it.
type Distinct<M> = Vec<(M, Vec<usize>)>;

/// One ID's distinct records, with their positions by [`fingerprint`].
type ById<'e, M> = (&'e [u8], Distinct<M>, HashMap<u64, Vec<usize>>);

/// One header record of one input, by its ID.
#[derive(Clone)]
struct Entry<M> {
    id: Id,
    map: M,
    input: usize,
}

/// The merged header and the ID translations per input.
pub struct HeaderMerger {
    pub header: sam::Header,
    read_group_ids: Vec<HashMap<Id, Id>>,
    program_ids: Vec<HashMap<Id, Id>>,
    read_group_collisions: bool,
    program_collisions: bool,
    counter: i32,
    /// Per input, the merged index of each of its references; empty unless
    /// the dictionaries were merged.
    sequence_maps: Vec<Vec<i32>>,
}

impl HeaderMerger {
    /// The merged header of `headers`, in input order, with the first
    /// input's sort order and equal dictionaries.
    pub fn merge(headers: &[sam::Header]) -> Result<Self> {
        Self::merge_with(headers, SortOrder::of(&headers[0]), false)
    }

    /// `new SamFileHeaderMerger(order, headers, merge_dictionaries)`: the
    /// dictionaries are merged only when they are not all equal.
    pub fn merge_with(
        headers: &[sam::Header],
        order: SortOrder,
        merge_dictionaries: bool,
    ) -> Result<Self> {
        let first = &headers[0];
        let equal = headers[1..]
            .iter()
            .try_for_each(|header| same_dictionary(first, header));
        let (dictionary, sequence_maps) = match equal {
            Ok(()) => (first.reference_sequences().clone(), Vec::new()),
            Err(_) if merge_dictionaries => {
                let dictionary = self::merge_dictionaries(headers)?;
                let maps = headers
                    .iter()
                    .map(|header| {
                        header
                            .reference_sequences()
                            .keys()
                            .map(|name| {
                                dictionary
                                    .get_index_of(name)
                                    .and_then(|at| i32::try_from(at).ok())
                                    .unwrap_or(-1)
                            })
                            .collect()
                    })
                    .collect();
                (dictionary, maps)
            }
            Err(error) => return Err(error),
        };
        let mut merger = Self {
            header: sam::Header::default(),
            read_group_ids: vec![HashMap::new(); headers.len()],
            program_ids: vec![HashMap::new(); headers.len()],
            read_group_collisions: false,
            program_collisions: false,
            counter: 0,
            sequence_maps,
        };
        let mut header = sam::Header::default();
        *header.header_mut() = Some(Map::<sam::header::record::value::map::Header>::new(
            Version::new(1, 6),
        ));
        *header.reference_sequences_mut() = dictionary;
        for (id, map) in merger.merge_programs(headers)? {
            header.programs_mut().as_mut().insert(id.into(), map);
        }
        for (id, map) in merger.merge_read_groups(headers)? {
            header.read_groups_mut().insert(id.into(), map);
        }
        set_header_tag(&mut header, HeaderTag::GroupOrder, "none");
        set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
        for input in headers {
            for comment in input.comments() {
                header.comments_mut().push(comment.clone());
            }
        }
        merger.header = header;
        Ok(merger)
    }

    /// True when the dictionaries were merged.
    pub fn has_merged_dictionary(&self) -> bool {
        !self.sequence_maps.is_empty()
    }

    fn merge_read_groups(&mut self, headers: &[sam::Header]) -> Result<Vec<(Id, Map<ReadGroup>)>> {
        let entries: Vec<Entry<Map<ReadGroup>>> = headers
            .iter()
            .enumerate()
            .flat_map(|(input, h)| {
                h.read_groups().iter().map(move |(id, map)| Entry {
                    id: id.to_vec(),
                    map: map.clone(),
                    input,
                })
            })
            .collect();
        self.counter = 0;
        let mut taken = HashSet::new();
        let mut ids = std::mem::take(&mut self.read_group_ids);
        let (mut result, collisions) =
            self.merge_records(&entries, &mut taken, &mut ids, &mut Vec::new());
        self.read_group_ids = ids;
        self.read_group_collisions = collisions;
        result.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(result)
    }

    fn merge_programs(&mut self, headers: &[sam::Header]) -> Result<Vec<(Id, Map<Program>)>> {
        let entries = headers.iter().enumerate().flat_map(|(input, h)| {
            h.programs().as_ref().iter().map(move |(id, map)| Entry {
                id: id.to_vec(),
                map: map.clone(),
                input,
            })
        });
        self.counter = 0;
        let (mut current, left): (Vec<_>, Vec<_>) =
            entries.partition(|e| programs::previous_of(&e.map).is_none());
        let mut left = Left::new(left);
        let mut taken = HashSet::new();
        let mut overall = Vec::new();
        let mut ids = std::mem::take(&mut self.program_ids);
        while !current.is_empty() {
            let mut changed = Vec::new();
            let (result, collisions) =
                self.merge_records(&current, &mut taken, &mut ids, &mut changed);
            self.program_collisions |= collisions;
            overall.extend(result);
            current = translate(current, &ids);
            left.translate(&ids, &changed);
            let mut next: Vec<usize> = current
                .iter()
                .flat_map(|done| left.take_children(done.input, &done.id))
                .collect();
            next.sort_unstable();
            current = next.into_iter().map(|at| left.entry(at).clone()).collect();
        }
        self.program_ids = ids;
        if left.remaining() > 0 {
            return Err(Error::Metric(format!(
                "{} program groups weren't processed. Do their PP ids point to existing PGs?",
                left.remaining()
            )));
        }
        overall.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(overall)
    }

    /// `mergeHeaderRecords`: one record per distinct (ID, record), a new ID
    /// for every one after the first that wants a taken ID. IDs and records
    /// are found by hash, in the order htsjdk's `LinkedHashMap`s keep;
    /// `changed` gets each (input, ID) given another ID.
    fn merge_records<I: Inner>(
        &mut self,
        entries: &[Entry<Map<I>>],
        taken: &mut HashSet<Id>,
        ids: &mut [HashMap<Id, Id>],
        changed: &mut Vec<(usize, Id)>,
    ) -> (Vec<(Id, Map<I>)>, bool)
    where
        Map<I>: Clone + PartialEq,
    {
        let mut slots: HashMap<&[u8], usize> = HashMap::new();
        let mut by_id: Vec<ById<'_, Map<I>>> = Vec::new();
        for entry in entries {
            let slot = *slots.entry(&entry.id).or_insert_with(|| {
                by_id.push((&entry.id, Vec::new(), HashMap::new()));
                by_id.len() - 1
            });
            let (_, records, index) = &mut by_id[slot];
            let bucket = index.entry(fingerprint(&entry.map)).or_default();
            match bucket.iter().find(|&&at| records[at].0 == entry.map) {
                Some(&at) => records[at].1.push(entry.input),
                None => {
                    bucket.push(records.len());
                    records.push((entry.map.clone(), vec![entry.input]));
                }
            }
        }
        let mut collisions = false;
        let mut result = Vec::new();
        for (id, records, _) in by_id {
            for (map, inputs) in records {
                let new_id = if taken.insert(id.to_vec()) {
                    self.counter += 1;
                    id.to_vec()
                } else {
                    collisions = true;
                    loop {
                        let candidate =
                            format!("{}.{}", String::from_utf8_lossy(id), base36(self.counter))
                                .into_bytes();
                        self.counter += 1;
                        if taken.insert(candidate.clone()) {
                            break candidate;
                        }
                    }
                };
                for input in inputs {
                    if new_id != id {
                        changed.push((input, id.to_vec()));
                    }
                    ids[input].insert(id.to_vec(), new_id.clone());
                }
                result.push((new_id, map));
            }
        }
        (result, collisions)
    }
}

/// A hash of a header record's fields that ignores their order, as
/// `AbstractSAMHeaderRecord.equals` compares them as a map.
fn fingerprint<I: Inner>(map: &Map<I>) -> u64 {
    map.other_fields().iter().fold(0u64, |sum, (tag, value)| {
        let mut hasher = DefaultHasher::new();
        tag.as_ref().hash(&mut hasher);
        <[u8]>::hash(value, &mut hasher);
        sum.wrapping_add(hasher.finish())
    })
}

/// `translateIds(programs, table, false)`: each program's ID through its
/// input's table.
fn translate(
    entries: Vec<Entry<Map<Program>>>,
    ids: &[HashMap<Id, Id>],
) -> Vec<Entry<Map<Program>>> {
    entries
        .into_iter()
        .map(|mut entry| {
            if let Some(new_id) = ids[entry.input].get(&entry.id) {
                entry.id = new_id.clone();
            }
            entry
        })
        .collect()
}

/// `SequenceUtil.assertSequenceDictionariesEqual`: as many sequences, and
/// at each index the same sequence (`SAMSequenceRecord.isSameSequence`).
fn same_dictionary(a: &sam::Header, b: &sam::Header) -> Result<()> {
    let (x, y) = (a.reference_sequences(), b.reference_sequences());
    if x.len() != y.len() {
        return Err(Error::Metric(format!(
            "Sequence dictionaries are not the same size ({}, {})",
            x.len(),
            y.len()
        )));
    }
    for (i, ((name_a, map_a), (name_b, map_b))) in x.iter().zip(y).enumerate() {
        if !same_sequence(name_a, map_a, name_b, map_b)? {
            return Err(Error::Metric(format!(
                "Sequences at index {i} don't match: {name_a}/{} {name_b}/{}",
                map_a.length(),
                map_b.length()
            )));
        }
    }
    Ok(())
}

/// `SAMSequenceRecord.isSameSequence` for two records at one index: the
/// lengths equal, the `M5` values equal as hexadecimal numbers when both
/// have one, and the names equal or one named among the other's `AN`.
fn same_sequence(
    name_a: &[u8],
    a: &Map<ReferenceSequence>,
    name_b: &[u8],
    b: &Map<ReferenceSequence>,
) -> Result<bool> {
    if a.length() != b.length() {
        return Ok(false);
    }
    let md5 = |map: &Map<ReferenceSequence>| map.other_fields().get(&MD5_CHECKSUM).cloned();
    if let (Some(x), Some(y)) = (md5(a), md5(b))
        && hex_number(&x)? != hex_number(&y)?
    {
        return Ok(false);
    }
    let named = |map: &Map<ReferenceSequence>, name: &[u8]| {
        map.other_fields()
            .get(&ALTERNATIVE_NAMES)
            .is_some_and(|names| names.split(|b| *b == b',').any(|alt| alt == name))
    };
    Ok(name_a == name_b || named(a, name_b) || named(b, name_a))
}

/// `new BigInteger(text, 16)` as digits without leading zeros: an `M5`
/// that is no hexadecimal number fails, as Java's constructor throws.
fn hex_number(text: &[u8]) -> Result<Vec<u8>> {
    let digits = text.strip_prefix(b"+").unwrap_or(text);
    if digits.is_empty() || !digits.iter().all(u8::is_ascii_hexdigit) {
        return Err(Error::Metric(format!(
            "M5 value {0}: For input string: \"{0}\" under radix 16",
            String::from_utf8_lossy(text)
        )));
    }
    let first = digits
        .iter()
        .position(|b| *b != b'0')
        .unwrap_or(digits.len());
    Ok(digits[first..].to_ascii_lowercase())
}

/// The order a [`Merger`] takes records in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MergeOrder {
    /// `SAMRecordCoordinateComparator`, or with a merged dictionary
    /// `MergedSequenceDictionaryCoordinateOrderComparator`.
    Coordinate,
    /// `SAMRecordQueryNameComparator`.
    Queryname,
    /// Every record of the first input, then of the second, and so on.
    InputOrder,
    /// `SortOrder.unknown`, which has no comparator: htsjdk's priority
    /// queue fails comparing the first records of two inputs, so more than
    /// one input with a record is an error; otherwise input order.
    NoComparator,
}

/// One input of a [`Merger`]: its records in file order.
pub trait RecordSource {
    /// The next record into `record`; false after the last.
    fn next_record(&mut self, record: &mut RawRecord) -> Result<bool>;
}

/// The records of a BAM record stream.
struct Stream<'a>(&'a mut (dyn Read + 'static));

impl RecordSource for Stream<'_> {
    fn next_record(&mut self, record: &mut RawRecord) -> Result<bool> {
        record.read_from(self.0).map_err(Error::Record)
    }
}

/// `MergingSamRecordIterator` over raw records.
pub struct Merger<'a> {
    sources: Vec<Box<dyn RecordSource + 'a>>,
    heads: Vec<Option<RawRecord>>,
    queue: Vec<usize>,
    order: MergeOrder,
    merged: HeaderMerger,
    /// The input [`MergeOrder::InputOrder`] takes records from.
    current: usize,
    /// The last record had its `RG` or `PG` set.
    rewritten: bool,
}

impl<'a> Merger<'a> {
    /// Take the records of `readers` (one per merged header, in order)
    /// in `order`: coordinate, queryname, or unsorted (input after input).
    pub fn new(
        readers: Vec<&'a mut (dyn Read + 'static)>,
        merged: HeaderMerger,
        order: SortOrder,
    ) -> Result<Self> {
        let order = match order {
            SortOrder::Queryname => MergeOrder::Queryname,
            SortOrder::Coordinate => MergeOrder::Coordinate,
            SortOrder::Unsorted => MergeOrder::InputOrder,
            other => {
                return Err(Error::Metric(format!(
                    "merging inputs in {} order is not supported",
                    other.name()
                )));
            }
        };
        let sources = readers
            .into_iter()
            .map(|reader| Box::new(Stream(reader)) as Box<dyn RecordSource + 'a>)
            .collect();
        Self::from_sources(sources, merged, order)
    }

    /// Take the records of `sources` (one per merged header, in order) in
    /// `order`. The first record of every input is read here.
    pub fn from_sources(
        sources: Vec<Box<dyn RecordSource + 'a>>,
        merged: HeaderMerger,
        order: MergeOrder,
    ) -> Result<Self> {
        let mut merger = Self {
            heads: vec![None; sources.len()],
            sources,
            queue: Vec::new(),
            order,
            merged,
            current: 0,
            rewritten: false,
        };
        for i in 0..merger.sources.len() {
            merger.advance(i)?;
        }
        if order == MergeOrder::NoComparator && merger.heads.iter().flatten().count() > 1 {
            return Err(Error::Metric(
                "SORT_ORDER=unknown has no comparator: htsjdk's MergingSamRecordIterator fails \
                 with a NullPointerException comparing the first records of two inputs"
                    .to_owned(),
            ));
        }
        Ok(merger)
    }

    fn queued(&self) -> bool {
        matches!(self.order, MergeOrder::Coordinate | MergeOrder::Queryname)
    }

    /// Read the next record of input `i` and queue the input if it has one.
    fn advance(&mut self, i: usize) -> Result<()> {
        let mut record = RawRecord::default();
        if self.sources[i].next_record(&mut record)? {
            self.heads[i] = Some(record);
            if self.queued() {
                self.offer(i);
            }
        } else {
            self.heads[i] = None;
        }
        Ok(())
    }

    fn compare(&self, a: usize, b: usize) -> Ordering {
        let (x, y) = (
            self.heads[a].as_ref().expect("queued inputs have a record"),
            self.heads[b].as_ref().expect("queued inputs have a record"),
        );
        match self.order {
            MergeOrder::Coordinate if self.merged.has_merged_dictionary() => {
                self.merged_key(a, x).compare(&self.merged_key(b, y))
            }
            MergeOrder::Queryname => queryname_order(x, y),
            _ => coordinate_order(x, y),
        }
    }

    /// The coordinate key of `record` of input `i` with its reference, or
    /// for an unplaced record its mate's, as an index of the merged
    /// dictionary. Its other fields stay the input's.
    fn merged_key<'r>(&self, i: usize, record: &'r RawRecord) -> CoordinateKey<'r> {
        let map = &self.merged.sequence_maps[i];
        let merged = |id: i32| {
            usize::try_from(id)
                .ok()
                .and_then(|at| map.get(at).copied())
                .unwrap_or(-1)
        };
        let mut key = CoordinateKey::of(record);
        key.reference_id = match (key.reference_id, key.mate_reference_id) {
            (-1, -1) => -1,
            (-1, mate) => merged(mate),
            (reference, _) => merged(reference),
        };
        key
    }

    /// `PriorityQueue.offer` (`siftUp`).
    fn offer(&mut self, item: usize) {
        let mut k = self.queue.len();
        self.queue.push(item);
        while k > 0 {
            let parent = (k - 1) >> 1;
            let e = self.queue[parent];
            if self.compare(item, e) != Ordering::Less {
                break;
            }
            self.queue[k] = e;
            k = parent;
        }
        self.queue[k] = item;
    }

    /// `PriorityQueue.poll` (`siftDown`).
    fn poll(&mut self) -> Option<usize> {
        if self.queue.is_empty() {
            return None;
        }
        let result = self.queue[0];
        let x = self.queue.pop().expect("not empty");
        let n = self.queue.len();
        if n > 0 {
            let mut k = 0;
            let half = n >> 1;
            while k < half {
                let mut child = 2 * k + 1;
                let mut c = self.queue[child];
                let right = child + 1;
                if right < n && self.compare(c, self.queue[right]) == Ordering::Greater {
                    child = right;
                    c = self.queue[child];
                }
                if self.compare(x, c) != Ordering::Greater {
                    break;
                }
                self.queue[k] = c;
                k = child;
            }
            self.queue[k] = x;
        }
        Some(result)
    }

    /// The input of the next record, in input order.
    fn next_input(&mut self) -> Option<usize> {
        while self.current < self.heads.len() && self.heads[self.current].is_none() {
            self.current += 1;
        }
        (self.current < self.heads.len()).then_some(self.current)
    }

    /// The next record into `record`, its IDs renamed and, with a merged
    /// dictionary, its references found again by name; false after the
    /// last.
    pub fn next(&mut self, record: &mut RawRecord) -> Result<bool> {
        Ok(self.next_from(record)?.is_some())
    }

    /// [`Merger::next`], with the index of the input the record came from;
    /// `None` after the last.
    pub fn next_from(&mut self, record: &mut RawRecord) -> Result<Option<usize>> {
        let next = if self.queued() {
            self.poll()
        } else {
            self.next_input()
        };
        let Some(i) = next else {
            return Ok(None);
        };
        *record = self.heads[i].take().expect("queued inputs have a record");
        self.advance(i)?;
        self.rewritten = false;
        if self.merged.read_group_collisions {
            self.rewritten |= translate_tag(record, *b"RG", &self.merged.read_group_ids[i])?;
        }
        if self.merged.program_collisions {
            self.rewritten |= translate_tag(record, *b"PG", &self.merged.program_ids[i])?;
        }
        if let Some(map) = self.merged.sequence_maps.get(i) {
            let name = || String::from_utf8_lossy(record.name()).into_owned();
            let merged = |id: i32| -> Result<i32> {
                let Ok(at) = usize::try_from(id) else {
                    return Ok(-1);
                };
                map.get(at).copied().ok_or_else(|| {
                    Error::Metric(format!(
                        "record {}: No mapping for reference index {id} from the header of \
                         input {}",
                        name(),
                        i + 1
                    ))
                })
            };
            let ids = (
                merged(record.reference_id())?,
                merged(record.mate_reference_id())?,
            );
            record.set_reference_ids(ids.0, ids.1);
        }
        Ok(Some(i))
    }

    /// True when the last record had its `RG` or `PG` tag set, which makes
    /// htsjdk encode the record again.
    pub fn rewritten(&self) -> bool {
        self.rewritten
    }
}

/// `record.setAttribute(tag, table.get(old))` for a record that has `tag`;
/// true when it has it.
fn translate_tag(record: &mut RawRecord, tag: [u8; 2], table: &HashMap<Id, Id>) -> Result<bool> {
    let mut tags = record.tags().map_err(Error::Record)?;
    let Some(old) = tags.get(tag) else {
        return Ok(false);
    };
    let TagValue::String(old) = old else {
        return Err(Error::Metric(format!(
            "{}: the {} tag is not a string",
            String::from_utf8_lossy(record.name()),
            String::from_utf8_lossy(&tag)
        )));
    };
    match table.get(old) {
        Some(new) => tags.set(tag, TagValue::String(new.clone())),
        None => tags.remove(tag),
    }
    record.set_tags(&tags).map_err(Error::Record)?;
    Ok(true)
}

#[cfg(test)]
mod tests;

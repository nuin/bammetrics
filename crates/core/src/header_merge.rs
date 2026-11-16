//! `SamFileHeaderMerger` and `MergingSamRecordIterator` of htsjdk 5.0.0,
//! shared by `MarkDuplicates` (several `INPUT`s) and the tools that merge
//! inputs: the merged header, the read group and program IDs renamed on a
//! collision, and the records of several sorted inputs taken in coordinate
//! or query-name order.
//!
//! * The sequence dictionaries must be equal; the merged header has
//!   `@HD VN:1.6 GO:none SO:<first input's order>`, the first input's
//!   `@SQ` lines, the read groups and programs of every input sorted by ID,
//!   and every comment in input order.
//! * A read group or program ID used by two inputs for different records
//!   gets a suffix `.<base 36 counter>`, and the records of that input are
//!   rewritten to the new ID (an `RG` or `PG` the input's header lacks is
//!   removed, as htsjdk sets it to null).
//! * Records are taken in coordinate or query-name order through a port of
//!   `java.util.PriorityQueue`, so equal records from two inputs come in the
//!   order Picard gives them.
//!
//! Also here: [`SortOrder`] (`SAMFileHeader.SortOrder`) and
//! [`set_header_tag`] (`setSortOrder`, `setGroupOrder`).

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::Read;

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::Version;
use noodles::sam::header::record::value::map::header::tag::{GROUP_ORDER, SORT_ORDER};
use noodles::sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;
use noodles::sam::header::record::value::map::{Program, ReadGroup};

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::sam_order::{coordinate_order, queryname_order};

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
}

impl HeaderMerger {
    /// The merged header of `headers`, in input order.
    pub fn merge(headers: &[sam::Header]) -> Result<Self> {
        let first = &headers[0];
        for header in &headers[1..] {
            same_dictionary(first, header)?;
        }
        let mut merger = Self {
            header: sam::Header::default(),
            read_group_ids: vec![HashMap::new(); headers.len()],
            program_ids: vec![HashMap::new(); headers.len()],
            read_group_collisions: false,
            program_collisions: false,
            counter: 0,
        };
        let mut header = sam::Header::default();
        *header.header_mut() = Some(Map::<sam::header::record::value::map::Header>::new(
            Version::new(1, 6),
        ));
        *header.reference_sequences_mut() = first.reference_sequences().clone();
        for (id, map) in merger.merge_programs(headers)? {
            header.programs_mut().as_mut().insert(id.into(), map);
        }
        for (id, map) in merger.merge_read_groups(headers)? {
            header.read_groups_mut().insert(id.into(), map);
        }
        set_header_tag(&mut header, HeaderTag::GroupOrder, "none");
        set_header_tag(
            &mut header,
            HeaderTag::SortOrder,
            SortOrder::of(first).name(),
        );
        for input in headers {
            for comment in input.comments() {
                header.comments_mut().push(comment.clone());
            }
        }
        merger.header = header;
        Ok(merger)
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
        let (mut result, collisions) = self.merge_records(&entries, &mut taken, &mut ids);
        self.read_group_ids = ids;
        self.read_group_collisions = collisions;
        result.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(result)
    }

    fn merge_programs(&mut self, headers: &[sam::Header]) -> Result<Vec<(Id, Map<Program>)>> {
        let mut left: Vec<Entry<Map<Program>>> = headers
            .iter()
            .enumerate()
            .flat_map(|(input, h)| {
                h.programs().as_ref().iter().map(move |(id, map)| Entry {
                    id: id.to_vec(),
                    map: map.clone(),
                    input,
                })
            })
            .collect();
        self.counter = 0;
        let previous = |e: &Entry<Map<Program>>| {
            e.map
                .other_fields()
                .get(&PREVIOUS_PROGRAM_ID)
                .map(|v| v.to_vec())
        };
        let mut current: Vec<Entry<Map<Program>>> = Vec::new();
        left.retain(|e| {
            if previous(e).is_none() {
                current.push(e.clone());
                false
            } else {
                true
            }
        });
        let mut taken = HashSet::new();
        let mut overall = Vec::new();
        let mut ids = std::mem::take(&mut self.program_ids);
        while !current.is_empty() {
            let (result, collisions) = self.merge_records(&current, &mut taken, &mut ids);
            self.program_collisions |= collisions;
            overall.extend(result);
            current = translate(current, &ids, false);
            left = translate(left, &ids, true);
            let mut next = Vec::new();
            let mut i = 0;
            while i < left.len() {
                let pp = previous(&left[i]).expect("only PP records are left");
                if current
                    .iter()
                    .any(|done| done.input == left[i].input && done.id == pp)
                {
                    next.push(left.remove(i));
                } else {
                    i += 1;
                }
            }
            current = next;
        }
        self.program_ids = ids;
        if !left.is_empty() {
            return Err(Error::Metric(format!(
                "{} program groups weren't processed. Do their PP ids point to existing PGs?",
                left.len()
            )));
        }
        overall.sort_by(|a, b| a.0.cmp(&b.0));
        Ok(overall)
    }

    /// `mergeHeaderRecords`: one record per distinct (ID, record), a new ID
    /// for every one after the first that wants a taken ID.
    fn merge_records<M: Clone + PartialEq>(
        &mut self,
        entries: &[Entry<M>],
        taken: &mut HashSet<Id>,
        ids: &mut [HashMap<Id, Id>],
    ) -> (Vec<(Id, M)>, bool) {
        let mut by_id: Vec<(Id, Distinct<M>)> = Vec::new();
        for entry in entries {
            let records = match by_id.iter_mut().find(|(id, _)| *id == entry.id) {
                Some((_, records)) => records,
                None => {
                    by_id.push((entry.id.clone(), Vec::new()));
                    &mut by_id.last_mut().expect("pushed").1
                }
            };
            match records.iter_mut().find(|(m, _)| *m == entry.map) {
                Some((_, inputs)) => inputs.push(entry.input),
                None => records.push((entry.map.clone(), vec![entry.input])),
            }
        }
        let mut collisions = false;
        let mut result = Vec::new();
        for (id, records) in by_id {
            for (map, inputs) in records {
                let new_id = if taken.insert(id.clone()) {
                    self.counter += 1;
                    id.clone()
                } else {
                    collisions = true;
                    loop {
                        let candidate =
                            format!("{}.{}", String::from_utf8_lossy(&id), base36(self.counter))
                                .into_bytes();
                        self.counter += 1;
                        if taken.insert(candidate.clone()) {
                            break candidate;
                        }
                    }
                };
                for input in inputs {
                    ids[input].insert(id.clone(), new_id.clone());
                }
                result.push((new_id, map));
            }
        }
        (result, collisions)
    }
}

/// `translateIds`: each program's ID, and with `previous` its `PP`, through
/// its input's table.
fn translate(
    entries: Vec<Entry<Map<Program>>>,
    ids: &[HashMap<Id, Id>],
    previous: bool,
) -> Vec<Entry<Map<Program>>> {
    entries
        .into_iter()
        .map(|mut entry| {
            let table = &ids[entry.input];
            if let Some(new_id) = table.get(&entry.id) {
                entry.id = new_id.clone();
            }
            if previous {
                let pp = entry
                    .map
                    .other_fields()
                    .get(&PREVIOUS_PROGRAM_ID)
                    .map(|v| v.to_vec());
                if let Some(new_pp) = pp.and_then(|pp| table.get(&pp).cloned()) {
                    entry
                        .map
                        .other_fields_mut()
                        .insert(PREVIOUS_PROGRAM_ID, new_pp.into());
                }
            }
            entry
        })
        .collect()
}

/// `SequenceUtil.assertSequenceDictionariesEqual`: names and lengths.
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
        if name_a != name_b || map_a.length() != map_b.length() {
            return Err(Error::Metric(format!(
                "Sequences at index {i} don't match: {name_a}/{} {name_b}/{}",
                map_a.length(),
                map_b.length()
            )));
        }
    }
    Ok(())
}

/// The order records are merged in.
#[derive(Clone, Copy)]
enum Order {
    Coordinate,
    Queryname,
}

/// `MergingSamRecordIterator` over raw records.
pub struct Merger<'a> {
    readers: Vec<&'a mut (dyn Read + 'static)>,
    heads: Vec<Option<RawRecord>>,
    queue: Vec<usize>,
    order: Order,
    merged: HeaderMerger,
}

impl<'a> Merger<'a> {
    /// Take the records of `readers` (one per merged header, in order)
    /// in `order`, coordinate or queryname.
    pub fn new(
        readers: Vec<&'a mut (dyn Read + 'static)>,
        merged: HeaderMerger,
        order: SortOrder,
    ) -> Result<Self> {
        let order = match order {
            SortOrder::Queryname => Order::Queryname,
            SortOrder::Coordinate => Order::Coordinate,
            other => {
                return Err(Error::Metric(format!(
                    "merging inputs in {} order is not supported",
                    other.name()
                )));
            }
        };
        let mut merger = Self {
            heads: vec![None; readers.len()],
            readers,
            queue: Vec::new(),
            order,
            merged,
        };
        for i in 0..merger.readers.len() {
            merger.advance(i)?;
        }
        Ok(merger)
    }

    /// Read the next record of input `i` and queue the input if it has one.
    fn advance(&mut self, i: usize) -> Result<()> {
        let mut record = RawRecord::default();
        if record.read_from(self.readers[i]).map_err(Error::Record)? {
            self.heads[i] = Some(record);
            self.offer(i);
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
            Order::Coordinate => coordinate_order(x, y),
            Order::Queryname => queryname_order(x, y),
        }
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

    /// The next record into `record`, its IDs renamed; false after the
    /// last.
    pub fn next(&mut self, record: &mut RawRecord) -> Result<bool> {
        let Some(i) = self.poll() else {
            return Ok(false);
        };
        *record = self.heads[i].take().expect("queued inputs have a record");
        self.advance(i)?;
        if self.merged.read_group_collisions {
            translate_tag(record, *b"RG", &self.merged.read_group_ids[i])?;
        }
        if self.merged.program_collisions {
            translate_tag(record, *b"PG", &self.merged.program_ids[i])?;
        }
        Ok(true)
    }
}

/// `record.setAttribute(tag, table.get(old))` for a record that has `tag`.
fn translate_tag(record: &mut RawRecord, tag: [u8; 2], table: &HashMap<Id, Id>) -> Result<()> {
    let mut tags = record.tags().map_err(Error::Record)?;
    let Some(old) = tags.get(tag) else {
        return Ok(());
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
    record.set_tags(&tags).map_err(Error::Record)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base36_counts_like_htsjdk() {
        assert_eq!(base36(0), "0");
        assert_eq!(base36(35), "Z");
        assert_eq!(base36(36), "10");
        assert_eq!(base36(1297), "101");
    }
}

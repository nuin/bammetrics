//! The record sort of `SortSam` and of the sort step of a chain: htsjdk's
//! `SortingCollection` over [`RawRecord`]s in one of the orders of
//! [`crate::sam_order`].
//!
//! Each record goes into the [`Sorter`] with what its order reads that the
//! record bytes do not hold at a fixed offset: the `HI` tag for queryname,
//! the whole [`DuplicateKey`] for duplicate. A run held in memory is sorted
//! by the port of Java's `TimSort` ([`crate::timsort`]), as htsjdk sorts it
//! with `Arrays.sort`, so a queryname sort of records flagged both first
//! and second of a pair (an order that is not total) gives Picard's order
//! and never panics. Runs spill as BGZF level 3 ([`SPILL_LEVEL`]) and merge
//! with ties to the earlier run, the order `SortingCollection` gives.
//!
//! The buffer holds at most a number of records (`MAX_RECORDS_IN_RAM`, the
//! bound Picard has) and at most a number of bytes ([`Limits::bytes`],
//! `--sort-memory` of a chain). [`Keyed::memory`] counts the record bytes,
//! the item and the two index slots the sort uses.

use std::cmp::Ordering;
use std::io::{self, BufRead};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::header_merge::{HeaderTag, set_header_tag};
use crate::sam_order::{DuplicateKey, Libraries, Order, QuerynameKey, coordinate_order, hit_index};
use crate::spill::{
    SortStats, Sorted, Sorter, Spill, TempDir, get_u64, need_i64, need_u64, put_i64, put_u64,
};
use crate::stage::{BATCH_RECORDS, Records};

/// BGZF level of spilled runs and spools. docs/07 section 2.1 says 1.
/// Measured on 8,445,233 exome records (docs/chaining.md section 11): level
/// 1 of the zlib-rs deflater spilled a queryname sort as 1,142,443,761
/// bytes in 16.7 s, level 3 as 871,834,880 bytes in 17.5 s; a coordinate
/// sort 820,851,148 bytes in 8.1 s against 616,683,279 in 8.9 s.
pub const SPILL_LEVEL: u8 = 3;

/// How much a sort holds in memory, and how many threads compress a run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Limits {
    /// Most records held (`MAX_RECORDS_IN_RAM`).
    pub records: usize,
    /// Most bytes held, [`Keyed::memory`] summed.
    pub bytes: usize,
    /// BGZF workers that compress a spilled run.
    pub threads: usize,
}

/// Bytes [`Keyed::memory`] counts for one record beyond its own bytes (and
/// the key of a duplicate sort): the item and the two index slots of the
/// sort.
pub const ITEM_BYTES: usize = std::mem::size_of::<Keyed>() + 2 * std::mem::size_of::<usize>();

/// One record in the sort, with the parts of its key the bytes lack.
#[derive(Debug)]
pub struct Keyed {
    pub record: RawRecord,
    order: Order,
    /// `HI` of a queryname sort.
    hit: Option<i64>,
    /// The key of a duplicate sort.
    duplicate: Option<Box<DuplicateKey<'static>>>,
}

impl Keyed {
    /// `record` keyed for `order`. An error is the exception htsjdk throws
    /// while it reads the key: an `HI` that is no number, a paired record
    /// with a mapped mate and no `MC` under `duplicate`.
    pub fn new(record: RawRecord, order: Order, libraries: &Libraries) -> Result<Self> {
        let (hit, duplicate) = match order {
            Order::Queryname => (hit_index(&record)?, None),
            Order::Duplicate => {
                let key = DuplicateKey::of(&record, libraries)?.into_owned();
                (None, Some(Box::new(key)))
            }
            Order::Coordinate | Order::Unsorted => (None, None),
        };
        Ok(Self {
            record,
            order,
            hit,
            duplicate,
        })
    }

    fn queryname_key(&self) -> QuerynameKey<'_> {
        QuerynameKey::new(self.record.name(), self.record.flags(), self.hit)
    }
}

impl PartialEq for Keyed {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Keyed {}

impl PartialOrd for Keyed {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// The htsjdk comparator of the order. Not a total order for queryname
/// outside valid flags; [`Spill::sort_run`] and the run merge accept that.
impl Ord for Keyed {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.order, &self.duplicate, &other.duplicate) {
            (Order::Coordinate, ..) => coordinate_order(&self.record, &other.record),
            (Order::Queryname, ..) => self.queryname_key().compare(&other.queryname_key()),
            (Order::Duplicate, Some(a), Some(b)) => a.compare(b),
            _ => Ordering::Equal,
        }
    }
}

fn order_code(order: Order) -> u8 {
    match order {
        Order::Unsorted => 0,
        Order::Queryname => 1,
        Order::Coordinate => 2,
        Order::Duplicate => 3,
    }
}

fn invalid(what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("sort run: {what}"))
}

fn need_i32(input: &mut impl BufRead) -> io::Result<i32> {
    i32::try_from(need_i64(input)?).map_err(|_| invalid("a value past 32 bits"))
}

fn put_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    put_u64(out, bytes.len() as u64);
    out.extend_from_slice(bytes);
}

fn need_bytes(input: &mut impl BufRead) -> io::Result<Vec<u8>> {
    let len = usize::try_from(need_u64(input)?).map_err(|_| invalid("a length past memory"))?;
    let mut bytes = vec![0; len];
    input.read_exact(&mut bytes)?;
    Ok(bytes)
}

impl Spill for Keyed {
    fn encode(&self, out: &mut Vec<u8>) {
        out.push(order_code(self.order));
        match self.hit {
            Some(hit) => {
                out.push(1);
                put_i64(out, hit);
            }
            None => out.push(0),
        }
        if let Some(key) = &self.duplicate {
            for value in [
                key.library,
                key.reference_id,
                key.read_coordinate,
                key.mate_reference_id,
                key.mate_coordinate,
                i32::from(key.score),
            ] {
                put_i64(out, i64::from(value));
            }
            put_u64(out, u64::from(key.flags));
            match &key.read_group {
                Some(group) => {
                    out.push(1);
                    put_bytes(out, group);
                }
                None => out.push(0),
            }
        }
        put_bytes(out, self.record.as_bytes());
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        let Some(code) = get_u64(input)? else {
            return Ok(None);
        };
        let order = match code {
            0 => Order::Unsorted,
            1 => Order::Queryname,
            2 => Order::Coordinate,
            3 => Order::Duplicate,
            _ => return Err(invalid("an unknown order")),
        };
        let mut flag = [0u8];
        input.read_exact(&mut flag)?;
        let hit = match flag[0] {
            0 => None,
            _ => Some(need_i64(input)?),
        };
        let mut duplicate = None;
        if order == Order::Duplicate {
            let library = need_i32(input)?;
            let reference_id = need_i32(input)?;
            let read_coordinate = need_i32(input)?;
            let mate_reference_id = need_i32(input)?;
            let mate_coordinate = need_i32(input)?;
            let score = i16::try_from(need_i64(input)?).map_err(|_| invalid("a score"))?;
            let flags = u16::try_from(need_u64(input)?).map_err(|_| invalid("flags"))?;
            input.read_exact(&mut flag)?;
            let read_group = match flag[0] {
                0 => None,
                _ => Some(need_bytes(input)?),
            };
            duplicate = Some((
                library,
                reference_id,
                read_coordinate,
                flags,
                mate_reference_id,
                mate_coordinate,
                score,
                read_group,
            ));
        }
        let record = RawRecord::from_bytes(need_bytes(input)?)?;
        let duplicate = duplicate.map(
            |(
                library,
                reference_id,
                read_coordinate,
                flags,
                mate_reference_id,
                mate_coordinate,
                score,
                read_group,
            )| {
                Box::new(DuplicateKey {
                    library,
                    reference_id,
                    read_coordinate,
                    flags,
                    mate_reference_id,
                    mate_coordinate,
                    score,
                    read_group,
                    name: std::borrow::Cow::Owned(record.name().to_vec()),
                })
            },
        );
        Ok(Some(Self {
            record,
            order,
            hit,
            duplicate,
        }))
    }

    fn memory(&self) -> usize {
        // Two index slots: the index TimSort moves and its merge buffer.
        let fixed = ITEM_BYTES;
        let duplicate = self.duplicate.as_ref().map_or(0, |key| {
            std::mem::size_of::<DuplicateKey<'static>>()
                + key.name.len()
                + key.read_group.as_ref().map_or(0, Vec::len)
        });
        fixed + self.record.as_bytes().len() + duplicate
    }

    fn sort_run(run: &mut Vec<Self>) -> Result<()> {
        crate::timsort::sort_by(run, Ord::cmp).map_err(|v| Error::Metric(v.to_string()))
    }
}

/// The header of a sorted stream: `SAMFileHeader.setSortOrder`, with the
/// version a htsjdk writer gives it. `GO` and `SS` stay as they are, as
/// Picard 3.5.0 `SortSam` keeps them.
pub fn sorted_header(header: &sam::Header, order: Order) -> sam::Header {
    let mut out = header.clone();
    set_header_tag(&mut out, HeaderTag::SortOrder, order.name());
    if let Some(map) = out.header_mut().as_mut() {
        *map.version_mut() = sam::header::record::value::map::header::Version::new(1, 6);
    }
    out
}

/// Records pushed in any order, given back in `order`.
pub struct RecordSorter<'d> {
    order: Order,
    libraries: Libraries,
    sorter: Sorter<'d, Keyed>,
}

impl<'d> RecordSorter<'d> {
    /// A sort of the records of `header` in `order`, spilling to `dir`.
    pub fn new(order: Order, header: &sam::Header, dir: &'d TempDir, limits: Limits) -> Self {
        let libraries = if order == Order::Duplicate {
            Libraries::new(header)
        } else {
            Libraries::default()
        };
        Self {
            order,
            libraries,
            sorter: Sorter::new(dir, "sort", limits.records)
                .byte_limit(limits.bytes)
                .compressed(SPILL_LEVEL, limits.threads),
        }
    }

    pub fn push(&mut self, record: RawRecord) -> Result<()> {
        self.sorter
            .push(Keyed::new(record, self.order, &self.libraries)?)
    }

    pub fn stats(&self) -> SortStats {
        self.sorter.stats()
    }

    /// Every record in order, and what the sort did.
    pub fn finish(self) -> Result<(SortedRecords<'d>, SortStats)> {
        let (sorted, stats) = self.sorter.finish_with_stats()?;
        Ok((SortedRecords { sorted }, stats))
    }
}

/// The sorted records; run files are removed once read.
pub struct SortedRecords<'d> {
    sorted: Sorted<'d, Keyed>,
}

impl SortedRecords<'_> {
    pub fn next_record(&mut self) -> Result<Option<RawRecord>> {
        Ok(self.sorted.next_item()?.map(|keyed| keyed.record))
    }
}

impl Records for SortedRecords<'_> {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        out.clear();
        while out.len() < BATCH_RECORDS {
            match self.next_record()? {
                Some(record) => out.push(record),
                None => break,
            }
        }
        Ok(!out.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_raw::build;

    fn names(records: &[RawRecord]) -> Vec<(String, u16)> {
        records
            .iter()
            .map(|r| (String::from_utf8_lossy(r.name()).into_owned(), r.flags()))
            .collect()
    }

    fn sort(order: Order, records: Vec<RawRecord>, limit: usize) -> (Vec<RawRecord>, SortStats) {
        let dir = TempDir::new(None);
        let header = sam::Header::default();
        let limits = Limits {
            records: limit,
            bytes: usize::MAX,
            threads: 2,
        };
        let mut sorter = RecordSorter::new(order, &header, &dir, limits);
        for record in records {
            sorter.push(record).expect("push");
        }
        let (mut sorted, stats) = sorter.finish().expect("finish");
        let mut out = Vec::new();
        while let Some(record) = sorted.next_record().expect("next") {
            out.push(record);
        }
        (out, stats)
    }

    /// Records whose queryname order is not total (both first and second
    /// of a pair) never panic. Runs of 3 records are sorted by binary
    /// insertion and merged, and the sort completes, as htsjdk's does with
    /// `MAX_RECORDS_IN_RAM=3`. One run of 200 ends in the violation Java's
    /// `TimSort` throws for these records.
    #[test]
    fn a_queryname_sort_of_first_and_second_records_never_panics() {
        let both = 0x1 | 0x40 | 0x80;
        let records: Vec<RawRecord> = (0..200)
            .map(|i| {
                build(
                    &format!("r{}", i % 7),
                    both | (i % 2) * 0x10,
                    0,
                    1,
                    &[],
                    &[],
                    b"",
                )
            })
            .collect();
        let (out, _) = sort(Order::Queryname, records.clone(), 3);
        assert_eq!(out.len(), 200);
        let mut seen: Vec<_> = names(&out).into_iter().map(|(n, _)| n).collect();
        assert!(seen.windows(2).all(|w| w[0] <= w[1]), "names in order");
        seen.dedup();
        assert_eq!(seen.len(), 7);
        let dir = TempDir::new(None);
        let limits = Limits {
            records: 1000,
            bytes: usize::MAX,
            threads: 1,
        };
        let mut sorter = RecordSorter::new(Order::Queryname, &sam::Header::default(), &dir, limits);
        for record in records {
            sorter.push(record).expect("push");
        }
        let error = sorter.finish().err().expect("the violation");
        assert_eq!(
            error.to_string(),
            "Comparison method violates its general contract!"
        );
    }

    /// Equal keys keep the order they came in, through every merge level.
    #[test]
    fn ties_keep_file_order_through_three_merge_levels() {
        let records: Vec<RawRecord> = (0..9000)
            .map(|i| {
                build(
                    &format!("n{:05}", i),
                    0,
                    0,
                    (i % 3) as i32 + 1,
                    &[],
                    &[],
                    b"",
                )
            })
            .collect();
        let (out, stats) = sort(Order::Coordinate, records.clone(), 2);
        assert!(stats.merge_levels >= 3, "{stats:?}");
        let mut expected = records;
        expected.sort_by_key(|r| r.alignment_start());
        assert_eq!(names(&out), names(&expected));
    }

    #[test]
    fn spilled_runs_round_trip_every_key() {
        let mut tagged = build("b", 0x1 | 0x40, 0, 5, &[], &[], b"");
        let mut tags = tagged.tags().expect("tags");
        tags.set(*b"HI", crate::bam_raw::TagValue::Int(2));
        tagged.set_tags(&tags).expect("set");
        let records = vec![
            tagged,
            build("b", 0x1 | 0x40, 0, 5, &[], &[], b""),
            build("a", 0x4, -1, 0, &[], &[], b""),
        ];
        let (out, _) = sort(Order::Queryname, records, 1);
        assert_eq!(
            names(&out),
            [
                ("a".to_owned(), 0x4),
                ("b".to_owned(), 0x41),
                ("b".to_owned(), 0x41)
            ]
        );
        assert_eq!(hit_index(&out[2]).expect("HI"), Some(2));
    }

    #[test]
    fn the_header_takes_the_order_and_version_1_6() {
        let text = "@HD\tVN:1.4\tSO:coordinate\tGO:query\n@SQ\tSN:c\tLN:9\n";
        let header: sam::Header = text.parse().expect("header");
        let sorted = sorted_header(&header, Order::Queryname);
        let out = crate::bam_raw::header_text(&sorted).expect("text");
        assert!(
            out.starts_with("@HD\tVN:1.6\tSO:queryname\tGO:query\n"),
            "{out}"
        );
        let bare = sorted_header(&sam::Header::default(), Order::Duplicate);
        let out = crate::bam_raw::header_text(&bare).expect("text");
        assert_eq!(out, "@HD\tVN:1.6\tSO:duplicate\n");
    }
}

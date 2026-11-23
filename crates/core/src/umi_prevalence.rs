//! `CollectUmiPrevalenceMetrics`, mirroring Picard 3.5.0: how many
//! distinct UMIs each duplicate set holds, as a histogram of duplicate sets
//! by UMI count. No metrics class: the file holds the `## HISTOGRAM`
//! section alone (`numUmis`, `duplicateSets`, `java.lang.Integer` keys).
//!
//! Picard reads every record once, in file order, so neither an index nor
//! a sort order is needed (its help asks for an indexed BAM; the code never
//! opens the index). A record is kept when it passes, in this order:
//! mapped; mapping quality at least `MINIMUM_MQ`; neither secondary nor
//! supplementary; carries the `BARCODE_TAG` tag; the barcode quality
//! filter; and, with `FILTER_UNPAIRED_READS`, paired with its mate mapped.
//!
//! The barcode quality filter is Picard's `BarcodeQualityFilter` as
//! written: a record without the `BARCODE_BQ` tag is kept, and a record
//! with it is kept only when at least one quality (Phred+33, spaces
//! removed) is below `MINIMUM_BARCODE_BQ`. A barcode whose every base is
//! good is dropped. The test is inverted in Picard's source; bammetrics
//! keeps it so the counts equal Picard's (KNOWN_DIFFERENCES.md lists it as
//! Picard's behavior).
//!
//! The kept records are sorted as `DuplicateSetIterator` sorts them,
//! `SAMRecordDuplicateComparator` through a `SortingCollection`
//! ([`crate::spill::Sorter`] with `MAX_RECORDS_IN_RAM` and `TMP_DIR`, a run
//! also spilled at 256 MiB of keys and UMIs), and each run of records
//! equal to the run's first under `duplicateSetCompare` ([`DuplicateKey::set_order`]) is one duplicate
//! set. Each read of a pair is its own record, so a pair's two reads fall
//! in two sets. A set counts the distinct `BARCODE_TAG` strings of its
//! records, held in one hash set as Picard holds them, so peak memory grows
//! with the distinct UMIs of the largest set.
//!
//! The key needs `MC` on a paired record whose mate is mapped, as htsjdk
//! does. htsjdk computes it only when it compares two records, so a run
//! that keeps one record never fails on it; bammetrics holds that error
//! until a second record is kept.

use std::cmp::Ordering;
use std::collections::HashSet;
use std::io::{self, BufRead};
use std::path::PathBuf;

use crate::bam_raw::{HtsjdkTags, RawRecord, TagIter, TagValue};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput, ToolOutput};
use crate::metrics_file::KeyClass;
use crate::pass;
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_order::{DuplicateKey, Libraries, Scoring};
use crate::spill::{Sorter, Spill, TempDir, get_u64, put_bytes};

/// The histogram labels Picard writes.
pub const KEY_LABEL: &str = "numUmis";
pub const SERIES_LABEL: &str = "duplicateSets";

/// Picard's `MAX_RECORDS_IN_RAM` default, the size of a sorted run.
const DEFAULT_RECORDS_IN_RAM: usize = 500_000;

/// A sorted run also spills at this many bytes of keys and UMIs, so long
/// `BARCODE_TAG` strings cannot hold `MAX_RECORDS_IN_RAM` records in memory.
/// Picard bounds the run by records alone.
const SORT_BYTES: usize = 256 << 20;

const UNMAPPED: u16 = 0x4;
const PAIRED: u16 = 0x1;
const MATE_UNMAPPED: u16 = 0x8;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

pub const TOOL: Tool = Tool {
    name: "CollectUmiPrevalenceMetrics",
    summary: "a histogram of duplicate sets by the number of distinct UMIs they hold \
              (numUmis, duplicateSets), over mapped primary records with a UMI; no metrics \
              row.",
    options: &[
        Opt::int(
            &["MINIMUM_MQ", "MQ"],
            "Lowest mapping quality a record needs (default 30).",
        ),
        Opt::text(&["BARCODE_TAG"], "The UMI tag (default RX)."),
        Opt::text(
            &["BARCODE_BQ"],
            "The UMI base quality tag, Phred+33 with spaces ignored (default BQ).",
        ),
        Opt::int(
            &["MINIMUM_BARCODE_BQ", "BQ"],
            "Picard's inverted test: a record with the BARCODE_BQ tag is kept only when a \
             quality is below this (default 30); one without the tag is kept.",
        ),
        Opt::bool(
            &["FILTER_UNPAIRED_READS", "FUR"],
            "Drop records that are not paired or whose mate is unmapped (default true).",
        ),
        Opt::int(
            &["PROGRESS_STEP_INTERVAL"],
            "Progress is not logged; 0 fails once a set is found, as Picard's progress \
             logger divides by it (default 1000000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the sorted runs (default: the system temporary directory).",
        )
        .repeated(),
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more \
             (default 500000).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

/// A tag option: its two bytes, or the error htsjdk's
/// `SAMTag.makeBinaryTag` throws when a record is first tested with it.
#[derive(Clone, Debug, PartialEq, Eq)]
enum TagName {
    Tag([u8; 2]),
    Refused(String),
}

impl TagName {
    fn parse(name: Option<&str>, default: &str) -> Self {
        match name.unwrap_or(default) {
            "null" => Self::Refused("a null tag name".to_owned()),
            text => match text.as_bytes() {
                [a, b] => Self::Tag([*a, *b]),
                _ => Self::Refused(format!("String tag does not have length() == 2: {text}")),
            },
        }
    }

    fn get(&self) -> Result<[u8; 2]> {
        match self {
            Self::Tag(tag) => Ok(*tag),
            Self::Refused(why) => Err(Error::Metric(why.clone())),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Options {
    minimum_mq: i32,
    barcode_tag: TagName,
    barcode_bq: TagName,
    /// `None` for `null`, which fails when the first quality is tested.
    minimum_barcode_bq: Option<i32>,
    filter_unpaired_reads: bool,
    progress_step_interval: i32,
    records_in_ram: usize,
    tmp_dir: Option<PathBuf>,
}

/// An `Integer` option that may be `null`.
fn nullable_int(inv: &Invocation, name: &str, default: i32) -> Result<Option<i32>> {
    match inv.get(name) {
        Some("null") => Ok(None),
        _ => inv.parse(name, default).map(Some),
    }
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let minimum_mq = nullable_int(inv, "MINIMUM_MQ", 30)?.ok_or_else(|| {
            Error::Usage("MINIMUM_MQ=null: a mapping quality is needed".to_owned())
        })?;
        let records: i64 = inv.parse("MAX_RECORDS_IN_RAM", DEFAULT_RECORDS_IN_RAM as i64)?;
        let records_in_ram = usize::try_from(records)
            .ok()
            .filter(|n| *n >= 1)
            .ok_or_else(|| {
                Error::Usage(format!("MAX_RECORDS_IN_RAM={records}: expected 1 or more"))
            })?;
        Ok(Self {
            minimum_mq,
            barcode_tag: TagName::parse(inv.get("BARCODE_TAG"), "RX"),
            barcode_bq: TagName::parse(inv.get("BARCODE_BQ"), "BQ"),
            minimum_barcode_bq: nullable_int(inv, "MINIMUM_BARCODE_BQ", 30)?,
            filter_unpaired_reads: inv.flag("FILTER_UNPAIRED_READS", true)?,
            progress_step_interval: inv.parse("PROGRESS_STEP_INTERVAL", 1_000_000)?,
            records_in_ram,
            tmp_dir: inv.all("TMP_DIR").first().map(PathBuf::from),
        })
    }
}

fn record_name(record: &RawRecord) -> String {
    String::from_utf8_lossy(record.name()).into_owned()
}

/// `getAttribute`: the field htsjdk's tag list keeps ([`HtsjdkTags`],
/// loaded for `record`), `None` when the list holds no such tag.
fn tag_value(record: &RawRecord, tags: &HtsjdkTags, tag: [u8; 2]) -> Result<Option<TagValue>> {
    let Some(index) = tags.position(tag) else {
        return Ok(None);
    };
    let field = TagIter::new(record.tag_bytes()).nth(index).transpose();
    Ok(field.map_err(Error::Record)?.map(|(_, value)| value))
}

/// `getStringAttribute`: the text, `None` without the tag, an error for
/// another type.
fn string_tag(record: &RawRecord, tags: &HtsjdkTags, tag: [u8; 2]) -> Result<Option<Vec<u8>>> {
    match tag_value(record, tags, tag)? {
        None => Ok(None),
        Some(TagValue::String(text)) => Ok(Some(text)),
        Some(_) => Err(Error::Metric(format!(
            "record {}: Value for tag {} is not a String",
            record_name(record),
            String::from_utf8_lossy(&tag)
        ))),
    }
}

impl Options {
    /// `BarcodeQualityFilter.filterOut`: true when the record has the tag
    /// and none of its qualities is below the minimum.
    fn good_barcode(&self, record: &RawRecord, tags: &HtsjdkTags) -> Result<bool> {
        let tag = self.barcode_bq.get()?;
        if tags.position(tag).is_none() {
            return Ok(false);
        }
        let text = string_tag(record, tags, tag)?.unwrap_or_default();
        // `SAMUtils.fastqToPhred` decodes the whole text before any test.
        let mut qualities = Vec::with_capacity(text.len());
        for &byte in text.iter().filter(|&&b| b != b' ') {
            if !(33..=126).contains(&byte) {
                return Err(Error::Metric(format!(
                    "record {}: Invalid fastq character in the {} tag",
                    record_name(record),
                    String::from_utf8_lossy(&tag)
                )));
            }
            qualities.push(i32::from(byte - 33));
        }
        let mut any_low = false;
        for quality in qualities {
            let minimum = self.minimum_barcode_bq.ok_or_else(|| {
                Error::Metric(
                    "MINIMUM_BARCODE_BQ=null: a barcode quality cannot be compared with it"
                        .to_owned(),
                )
            })?;
            if quality < minimum {
                any_low = true;
                break;
            }
        }
        Ok(!any_low)
    }

    /// The UMI of a record the filters keep, `None` for one they drop, in
    /// the order of Picard's `AggregateFilter`.
    fn umi(&self, record: &RawRecord, tags: &mut HtsjdkTags) -> Result<Option<Vec<u8>>> {
        let flags = record.flags();
        if flags & UNMAPPED != 0
            || i32::from(record.mapping_quality()) < self.minimum_mq
            || flags & (SECONDARY | SUPPLEMENTARY) != 0
        {
            return Ok(None);
        }
        let tag = self.barcode_tag.get()?;
        tags.try_load(TagIter::new(record.tag_bytes()).map(|field| field.map(|(name, _)| name)))
            .map_err(Error::Record)?;
        if tags.position(tag).is_none() || self.good_barcode(record, tags)? {
            return Ok(None);
        }
        if self.filter_unpaired_reads && (flags & PAIRED == 0 || flags & MATE_UNMAPPED != 0) {
            return Ok(None);
        }
        // The set's lambda reads the tag as a String.
        string_tag(record, tags, tag)
    }
}

/// One kept record in the sort: its duplicate key and UMI.
#[derive(Debug)]
struct Member {
    key: DuplicateKey<'static>,
    umi: Vec<u8>,
}

impl PartialEq for Member {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Member {}

impl PartialOrd for Member {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// `SAMRecordDuplicateComparator.compare`, the order of the sort.
impl Ord for Member {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.compare(&other.key)
    }
}

impl Spill for Member {
    fn encode(&self, out: &mut Vec<u8>) {
        put_bytes(out, &self.umi);
        self.key.encode(out);
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        let Some(len) = get_u64(input)? else {
            return Ok(None);
        };
        let umi = crate::spill::read_bytes(input, len)?;
        let key = DuplicateKey::decode(input)?;
        Ok(Some(Self { key, umi }))
    }

    fn memory(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.umi.len()
            + self.key.heap_bytes()
            + 2 * std::mem::size_of::<usize>()
    }

    fn sort_run(run: &mut Vec<Self>) -> Result<()> {
        crate::timsort::sort_by(run, Ord::cmp).map_err(|v| Error::Metric(v.to_string()))
    }
}

/// The kept records: the sort, and the key error htsjdk raises only once a
/// second record is kept.
struct Collected<'t> {
    sorter: Sorter<'t, Member>,
    kept: u64,
    held_error: Option<Error>,
}

/// `DuplicateSetIterator.next` over the sorted records: each set is the
/// run of records whose `duplicateSetCompare` against the set's
/// representative is 0. A record after the representative that compares
/// greater is htsjdk's "not sorted in duplicate order" error. The
/// representative becomes a later record that `compare` puts first.
fn count_sets(sorted: &mut crate::spill::Sorted<'_, Member>, interval: i32) -> Result<Histogram> {
    let mut histogram = Histogram::new();
    let Some(first) = sorted.next_item()? else {
        return Ok(histogram);
    };
    let mut scoring = Scoring::default();
    scoring.next(&first.key)?;
    let mut representative = first.key;
    let mut umis: HashSet<Vec<u8>> = HashSet::from([first.umi]);
    let mut close = |umis: &mut HashSet<Vec<u8>>| -> Result<()> {
        if interval == 0 {
            return Err(Error::Metric(
                "PROGRESS_STEP_INTERVAL=0: Picard's progress logger divides by it (/ by zero)"
                    .to_owned(),
            ));
        }
        histogram.increment(umis.len() as i64);
        umis.clear();
        Ok(())
    };
    while let Some(member) = sorted.next_item()? {
        scoring.next(&member.key)?;
        match representative.set_order(&member.key) {
            Ordering::Greater => {
                return Err(Error::Metric(format!(
                    "The input records were not sorted in duplicate order: {} then {}",
                    String::from_utf8_lossy(&representative.name),
                    String::from_utf8_lossy(&member.key.name)
                )));
            }
            Ordering::Less => {
                close(&mut umis)?;
                representative = member.key;
            }
            Ordering::Equal => {
                if representative.compare(&member.key) == Ordering::Greater {
                    representative = member.key;
                }
            }
        }
        umis.insert(member.umi);
    }
    close(&mut umis)?;
    Ok(histogram)
}

/// Read `INPUT`, keep, sort and count: the histogram Picard writes.
pub fn histogram(inv: &Invocation, opts: &Options) -> Result<Histogram> {
    let temp = TempDir::new(opts.tmp_dir.as_deref(), "umi-prevalence");
    pass::with_record_stream_at(
        &inv.input,
        inv.reference(),
        inv.threads,
        inv.stringency()?,
        |header, reader, cram| {
            let libraries = Libraries::new(header);
            let mut collected = Collected {
                sorter: Sorter::new(&temp, "umi", opts.records_in_ram)
                    .byte_limit(inv.sort_memory.unwrap_or(SORT_BYTES)),
                kept: 0,
                held_error: None,
            };
            let mut record = RawRecord::default();
            let mut tags = HtsjdkTags::default();
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&inv.input, cram, e))?
            {
                let Some(umi) = opts.umi(&record, &mut tags)? else {
                    continue;
                };
                collected.kept += 1;
                match DuplicateKey::of(&record, &libraries) {
                    Ok(key) => collected.sorter.push(Member {
                        key: key.into_owned(),
                        umi,
                    })?,
                    Err(error) => {
                        collected.held_error.get_or_insert(error);
                    }
                }
                if collected.kept >= 2
                    && let Some(error) = collected.held_error.take()
                {
                    return Err(error);
                }
            }
            if collected.held_error.is_some() {
                // One record kept: one set of one UMI, never compared.
                return count_one(opts.progress_step_interval);
            }
            let mut sorted = collected.sorter.finish()?;
            count_sets(&mut sorted, opts.progress_step_interval)
        },
    )
}

/// The histogram of a run that kept one record whose key could not be
/// read: one set with one UMI.
fn count_one(interval: i32) -> Result<Histogram> {
    if interval == 0 {
        return Err(Error::Metric(
            "PROGRESS_STEP_INTERVAL=0: Picard's progress logger divides by it (/ by zero)"
                .to_owned(),
        ));
    }
    let mut histogram = Histogram::new();
    histogram.increment(1);
    Ok(histogram)
}

/// The file: the histogram section alone.
pub fn output(histogram: Histogram) -> MetricsOutput {
    MetricsOutput::histogram_only(HistogramSection {
        key_label: KEY_LABEL,
        key_class: KeyClass::Integer,
        series: vec![HistogramSeries {
            label: SERIES_LABEL.to_owned(),
            histogram,
        }],
    })
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    Ok(output(histogram(inv, &opts)?).into())
}

#[cfg(test)]
mod tests;

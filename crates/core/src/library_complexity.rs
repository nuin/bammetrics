//! `EstimateLibraryComplexity`, mirroring Picard 3.5.0
//! (`EstimateLibraryComplexity`, `ElcIdenticalBasesDuplicatesFinder`,
//! `ElcHashBasedDuplicatesFinder`, `DuplicationMetrics`).
//!
//! Duplicates are found from the read sequences alone, so the input may be
//! unaligned. Rules copied from the source:
//! * Only paired primary records with the first or second flag count. The
//!   two mates of a name are joined within one input file; a name seen a
//!   third time starts a new, never finished pair.
//! * A pair is kept when both mates pass the quality check: at least
//!   `MIN_IDENTICAL_BASES` bases, no no-call among the first
//!   `MIN_IDENTICAL_BASES`, and a mean quality (integer division over at
//!   most `MAX_READ_LENGTH` bases) of at least `MIN_MEAN_QUALITY`. A
//!   reverse-strand mate is reverse complemented.
//! * The read group index, and so the library, is set only when the read
//!   name yields a tile, x and y; a pair whose name does not parse, or
//!   whose `RG` is missing from the header, goes to library `Unknown`, as
//!   does a group without `LB`. `READ_NAME_REGEX=null` therefore sends
//!   every pair to `Unknown`.
//! * Pairs are sorted, stably, by the first `MIN_IDENTICAL_BASES` bases of
//!   each mate. A run of pairs equal there is a group; a group larger than
//!   `MAX_GROUP_RATIO` times the expected group size is skipped. Within a
//!   group and library, pairs are taken in order: each untaken pair
//!   collects every later untaken pair that differs from it in at most
//!   `floor(MAX_DIFF_RATE * compared bases)` positions (bases from
//!   `MIN_IDENTICAL_BASES` to the shorter read, capped at
//!   `MAX_READ_LENGTH`; equal barcode tags when any is set). The set size
//!   goes to the library's `duplication_group_count` histogram, and the
//!   set's optical duplicates (the pair itself the keeper) to the optical
//!   count of that size. Picard compares groups of 100 pairs and more
//!   through 32-bit hashes of interleaved bases and only compares bases
//!   where the hashes differ; a hash collision would there hide a
//!   difference. bammetrics compares bases for every group.
//! * Metrics per library: only histogram bins with at least
//!   `MIN_GROUP_COUNT` sets count; `READ_PAIRS_EXAMINED` adds size times
//!   sets, `READ_PAIR_DUPLICATES` size minus one times sets,
//!   `READ_PAIR_OPTICAL_DUPLICATES` the optical count. The other counts
//!   stay 0. Rows and histogram columns come in Java `HashMap` order of
//!   the library names ([`crate::java`]).

use std::collections::HashMap;
use std::path::Path;

use noodles::sam::alignment::record::data::field::{Tag, Value as FieldValue};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag;
use noodles::{bam, sam};

use crate::duplicate_metrics::Row;
use crate::duplication::estimate_library_size;
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::java::{hash_map_order, string_hash_code};
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput, ToolOutput};
use crate::metrics_file::KeyClass;
use crate::optical::{
    DEFAULT_MAX_DUPLICATE_SET_SIZE, DEFAULT_OPTICAL_DUPLICATE_DISTANCE, Location,
    OpticalDuplicateFinder, ReadNameParser,
};
use crate::pass::{self, Accumulator};
use crate::quality::base_qualities;
use crate::record::record_name;
use crate::registry::{Invocation, Opt, Tool};
use crate::sequence::{is_no_call, reverse_complement};

/// The library of a pair without a read group or without `LB`.
pub const UNKNOWN_LIBRARY: &str = "Unknown";

#[derive(Clone, Debug)]
pub struct Options {
    pub min_identical_bases: i32,
    pub max_diff_rate: f64,
    pub min_mean_quality: i32,
    pub max_group_ratio: i32,
    pub barcode_tag: Option<Tag>,
    pub read_one_barcode_tag: Option<Tag>,
    pub read_two_barcode_tag: Option<Tag>,
    /// 0 means no cap.
    pub max_read_length: i32,
    pub min_group_count: i32,
    pub finder: OpticalDuplicateFinder,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let tag = |name: &str| -> Result<Option<Tag>> {
            match inv.get(name) {
                None | Some("null") => Ok(None),
                Some(text) => match text.as_bytes() {
                    [a, b] => Ok(Some(Tag::new(*a, *b))),
                    _ => Err(Error::Usage(format!("{name}={text}: not a two-letter tag"))),
                },
            }
        };
        let opts = Self {
            min_identical_bases: inv.parse("MIN_IDENTICAL_BASES", 5)?,
            max_diff_rate: inv.parse("MAX_DIFF_RATE", 0.03)?,
            min_mean_quality: inv.parse("MIN_MEAN_QUALITY", 20)?,
            max_group_ratio: inv.parse("MAX_GROUP_RATIO", 500)?,
            barcode_tag: tag("BARCODE_TAG")?,
            read_one_barcode_tag: tag("READ_ONE_BARCODE_TAG")?,
            read_two_barcode_tag: tag("READ_TWO_BARCODE_TAG")?,
            max_read_length: inv.parse("MAX_READ_LENGTH", 0)?,
            min_group_count: inv.parse("MIN_GROUP_COUNT", 2)?,
            finder: OpticalDuplicateFinder {
                parser: ReadNameParser::new(inv.get("READ_NAME_REGEX"))?,
                pixel_distance: inv.parse(
                    "OPTICAL_DUPLICATE_PIXEL_DISTANCE",
                    DEFAULT_OPTICAL_DUPLICATE_DISTANCE,
                )?,
                max_set_size: inv.parse(
                    "MAX_OPTICAL_DUPLICATE_SET_SIZE",
                    DEFAULT_MAX_DUPLICATE_SET_SIZE,
                )?,
            },
        };
        // `customCommandLineValidation`.
        if 0 < opts.max_read_length && opts.max_read_length < opts.min_identical_bases {
            return Err(Error::Usage(
                "MAX_READ_LENGTH must be greater than MIN_IDENTICAL_BASES".to_owned(),
            ));
        }
        if opts.min_identical_bases <= 0 {
            return Err(Error::Usage(
                "MIN_IDENTICAL_BASES must be greater than 0".to_owned(),
            ));
        }
        Ok(opts)
    }

    fn use_barcodes(&self) -> bool {
        self.barcode_tag.is_some()
            || self.read_one_barcode_tag.is_some()
            || self.read_two_barcode_tag.is_some()
    }

    /// `Integer.MAX_VALUE` for no cap.
    fn read_length_cap(&self) -> usize {
        if self.max_read_length <= 0 {
            i32::MAX as usize
        } else {
            self.max_read_length as usize
        }
    }
}

pub const TOOL: Tool = Tool {
    name: "EstimateLibraryComplexity",
    options: &[
        Opt::value(&["MIN_IDENTICAL_BASES"]),
        Opt::value(&["MAX_DIFF_RATE"]),
        Opt::value(&["MIN_MEAN_QUALITY"]),
        Opt::value(&["MAX_GROUP_RATIO"]),
        Opt::value(&["BARCODE_TAG"]),
        Opt::value(&["READ_ONE_BARCODE_TAG"]),
        Opt::value(&["READ_TWO_BARCODE_TAG"]),
        Opt::value(&["MAX_READ_LENGTH"]),
        Opt::value(&["MIN_GROUP_COUNT"]),
        Opt::value(&["READ_NAME_REGEX"]),
        Opt::value(&["OPTICAL_DUPLICATE_PIXEL_DISTANCE"]),
        Opt::value(&["MAX_OPTICAL_DUPLICATE_SET_SIZE"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let mut reader = Reader::new(opts);
    for input in inv.all("INPUT") {
        reader = pass::run(Path::new(input), inv.threads, |header| {
            reader.start_file(header);
            Ok(reader)
        })?;
    }
    Ok(reader.finish_all().into())
}

/// `PairedReadSequence`: both mates of a template, original orientation.
#[derive(Clone, Debug)]
struct Pair {
    location: Location,
    quality_ok: bool,
    read1: Option<Vec<u8>>,
    read2: Option<Vec<u8>>,
    barcode: i32,
    read_one_barcode: i32,
    read_two_barcode: i32,
}

impl Pair {
    fn new() -> Self {
        Self {
            location: Location::default(),
            quality_ok: true,
            read1: None,
            read2: None,
            barcode: 0,
            read_one_barcode: 0,
            read_two_barcode: 0,
        }
    }
}

/// A complete pair, as sorted and compared.
#[derive(Clone, Debug)]
struct Complete {
    location: Location,
    read1: Vec<u8>,
    read2: Vec<u8>,
    barcodes: [i32; 3],
}

/// Reads every input into complete pairs.
struct Reader {
    opts: Options,
    /// Every `@RG` line of every input so far, in order.
    read_groups: Vec<(Vec<u8>, Map<ReadGroup>)>,
    /// Read group ID to its index in `read_groups`, for the current file.
    index_by_id: HashMap<Vec<u8>, i16>,
    pending: HashMap<Vec<u8>, Pair>,
    pairs: Vec<Complete>,
    /// `ProgressLogger.getCount`: records that passed the flag filters.
    records: i64,
    quals: Vec<i8>,
}

impl Reader {
    fn new(opts: Options) -> Self {
        Self {
            opts,
            read_groups: Vec::new(),
            index_by_id: HashMap::new(),
            pending: HashMap::new(),
            pairs: Vec::new(),
            records: 0,
            quals: Vec::new(),
        }
    }

    /// `readGroups.addAll(header.getReadGroups())`: a line equal to one
    /// already listed (same ID and fields) keeps the earlier index, as
    /// `List.indexOf` finds it.
    fn start_file(&mut self, header: &sam::Header) {
        self.pending.clear();
        self.index_by_id.clear();
        for (id, map) in header.read_groups() {
            let id = id.to_vec();
            let index = self
                .read_groups
                .iter()
                .position(|(known_id, known)| *known_id == id && known == map)
                .unwrap_or_else(|| {
                    self.read_groups.push((id.clone(), map.clone()));
                    self.read_groups.len() - 1
                });
            self.index_by_id.entry(id).or_insert(index as i16);
        }
    }

    /// `passesQualityCheck`.
    fn passes_quality_check(&mut self, record: &bam::Record, bases: &[u8]) -> Result<bool> {
        let seed = self.opts.min_identical_bases as usize;
        if bases.len() < seed {
            return Ok(false);
        }
        if bases[..seed].iter().any(|&b| is_no_call(b)) {
            return Ok(false);
        }
        base_qualities(record, &mut self.quals);
        let length = bases.len().min(self.opts.read_length_cap());
        if self.quals.len() < length {
            return Err(Error::MalformedRecord {
                name: record_name(record),
                what: "fewer qualities than bases",
            });
        }
        let total: i64 = self.quals[..length].iter().map(|&q| i64::from(q)).sum();
        Ok(total / length as i64 >= i64::from(self.opts.min_mean_quality))
    }

    /// `getReadBarcodeValue`: the Java hash code of the tag's string, 0
    /// without the tag.
    fn barcode(record: &bam::Record, tag: Option<Tag>) -> Result<i32> {
        let Some(tag) = tag else {
            return Ok(0);
        };
        match record.data().get(&tag) {
            None => Ok(0),
            Some(Err(e)) => Err(e.into()),
            Some(Ok(FieldValue::String(text))) => {
                Ok(string_hash_code(&String::from_utf8_lossy(text)))
            }
            Some(Ok(_)) => Err(Error::MalformedRecord {
                name: record_name(record),
                what: "barcode tag is not a string",
            }),
        }
    }

    /// `SAMRecord.getReadGroup`: the header index of the record's `RG`.
    fn read_group_index(&self, record: &bam::Record) -> Result<Option<i16>> {
        match record.data().get(&Tag::READ_GROUP) {
            None => Ok(None),
            Some(Err(e)) => Err(e.into()),
            Some(Ok(FieldValue::String(id))) => {
                Ok(self.index_by_id.get(AsRef::<[u8]>::as_ref(id)).copied())
            }
            Some(Ok(_)) => Err(Error::MalformedRecord {
                name: record_name(record),
                what: "RG tag is not a string",
            }),
        }
    }

    fn library_name(&self, pair: &Complete) -> String {
        if pair.location.read_group == -1 {
            return UNKNOWN_LIBRARY.to_owned();
        }
        self.read_groups[pair.location.read_group as usize]
            .1
            .other_fields()
            .get(&tag::LIBRARY)
            .map_or_else(|| UNKNOWN_LIBRARY.to_owned(), ToString::to_string)
    }

    fn finish_all(self) -> MetricsOutput {
        estimate(self)
    }
}

impl Accumulator for Reader {
    type Output = Reader;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if !flags.is_segmented()
            || (!flags.is_first_segment() && !flags.is_last_segment())
            || flags.is_secondary()
            || flags.is_supplementary()
        {
            return Ok(());
        }
        let name: Vec<u8> = record.name().map(|n| n.to_vec()).unwrap_or_default();
        // A pending pair is taken out of the map for good: a name seen a
        // third time starts a new pair, and the one taken out is lost
        // unless this record completes it (Picard puts only a new pair in
        // `pendingByName`).
        let (mut pair, created) = match self.pending.remove(&name) {
            Some(pair) => (pair, false),
            None => {
                let mut pair = Pair::new();
                let text = String::from_utf8_lossy(&name);
                if let Some((tile, x, y)) = self.opts.finder.parser.location(&text) {
                    pair.location.tile = tile;
                    pair.location.x = x;
                    pair.location.y = y;
                    if let Some(index) = self.read_group_index(record)? {
                        pair.location.read_group = index;
                    }
                }
                (pair, true)
            }
        };
        let mut bases: Vec<u8> = record.sequence().iter().collect();
        let passes = self.passes_quality_check(record, &bases)?;
        pair.quality_ok = pair.quality_ok && passes;
        if flags.is_reverse_complemented() {
            bases = reverse_complement(&bases);
        }
        let barcodes = self.opts.use_barcodes();
        if flags.is_first_segment() {
            if barcodes {
                pair.barcode = Self::barcode(record, self.opts.barcode_tag)?;
                pair.read_one_barcode = Self::barcode(record, self.opts.read_one_barcode_tag)?;
            }
            pair.read1 = Some(bases);
        } else {
            if barcodes {
                pair.read_two_barcode = Self::barcode(record, self.opts.read_two_barcode_tag)?;
            }
            pair.read2 = Some(bases);
        }
        self.records += 1;
        if let (Some(read1), Some(read2), true) = (&pair.read1, &pair.read2, pair.quality_ok) {
            self.pairs.push(Complete {
                location: pair.location,
                read1: read1.clone(),
                read2: read2.clone(),
                barcodes: [pair.barcode, pair.read_one_barcode, pair.read_two_barcode],
            });
        } else if created {
            self.pending.insert(name, pair);
        }
        Ok(())
    }

    fn finish(self) -> Result<Reader> {
        Ok(self)
    }
}

/// One library's two histograms.
#[derive(Default)]
struct LibraryHistograms {
    duplication: Histogram,
    optical: Histogram,
}

/// `getNextGroup`: pairs sharing their first `n` bases on both mates.
fn same_group(a: &Complete, b: &Complete, n: usize) -> bool {
    a.read1[..n] == b.read1[..n] && a.read2[..n] == b.read2[..n]
}

/// `ElcIdenticalBasesDuplicatesFinder.matches`.
fn matches(lhs: &Complete, rhs: &Complete, opts: &Options) -> bool {
    let cap = opts.read_length_cap();
    let read1_length = lhs.read1.len().min(rhs.read1.len()).min(cap);
    let read2_length = lhs.read2.len().min(rhs.read2.len()).min(cap);
    let max_errors = ((read1_length + read2_length) as f64 * opts.max_diff_rate).floor() as i64;
    if opts.use_barcodes() && lhs.barcodes != rhs.barcodes {
        return false;
    }
    let mut errors = 0i64;
    let skip = opts.min_identical_bases as usize;
    for (a, b) in lhs.read1[..read1_length]
        .iter()
        .zip(&rhs.read1[..read1_length])
        .skip(skip)
    {
        if a != b {
            errors += 1;
            if errors > max_errors {
                return false;
            }
        }
    }
    for (a, b) in lhs.read2[..read2_length]
        .iter()
        .zip(&rhs.read2[..read2_length])
        .skip(skip)
    {
        if a != b {
            errors += 1;
            if errors > max_errors {
                return false;
            }
        }
    }
    true
}

/// `searchDuplicates` plus `fillHistogram` for one library's pairs of one
/// group.
fn search_duplicates(pairs: &[&Complete], opts: &Options, histograms: &mut LibraryHistograms) {
    let mut taken = vec![false; pairs.len()];
    for i in 0..pairs.len() {
        if taken[i] {
            continue;
        }
        let mut dupes: Vec<usize> = Vec::new();
        for j in i + 1..pairs.len() {
            if !taken[j] && matches(pairs[i], pairs[j], opts) {
                dupes.push(j);
                taken[j] = true;
            }
        }
        if dupes.is_empty() {
            histograms.duplication.increment(1);
            continue;
        }
        dupes.push(i);
        let count = dupes.len() as i64;
        histograms.duplication.increment(count);
        let locations: Vec<Location> = dupes.iter().map(|&k| pairs[k].location).collect();
        let flags = opts
            .finder
            .find_optical_duplicates(&locations, Some(locations.len() - 1));
        for flagged in flags {
            if flagged {
                histograms.optical.increment(count);
            }
        }
    }
}

/// `doWork` after reading: sort, group, find duplicates, write metrics.
fn estimate(mut reader: Reader) -> MetricsOutput {
    let opts = reader.opts.clone();
    let n = opts.min_identical_bases as usize;
    let mut pairs = std::mem::take(&mut reader.pairs);
    pairs.sort_by(|a, b| {
        a.read1[..n]
            .cmp(&b.read1[..n])
            .then(a.read2[..n].cmp(&b.read2[..n]))
    });
    let expected = 4f64
        .powi(opts.min_identical_bases * 2)
        .min(f64::from(i32::MAX)) as i64;
    let mean_group_size = ((reader.records / 2) / expected).max(1) as i32;
    let too_large = mean_group_size.wrapping_mul(opts.max_group_ratio);
    let mut libraries: Vec<String> = Vec::new();
    let mut by_library: HashMap<String, LibraryHistograms> = HashMap::new();
    let mut start = 0;
    while start < pairs.len() {
        let mut end = start + 1;
        while end < pairs.len() && same_group(&pairs[start], &pairs[end], n) {
            end += 1;
        }
        let group = &pairs[start..end];
        start = end;
        if group.len() as i64 > i64::from(too_large) {
            eprintln!(
                "bammetrics: omitting group with over {} times the expected mean number of read pairs. Mean={mean_group_size}, Actual={}",
                opts.max_group_ratio,
                group.len()
            );
            continue;
        }
        // `splitByLibrary`, a HashMap whose iteration order sets when each
        // library's histograms are first made.
        let mut group_libraries: Vec<String> = Vec::new();
        let mut members: Vec<Vec<&Complete>> = Vec::new();
        for pair in group {
            let library = reader.library_name(pair);
            match group_libraries.iter().position(|l| *l == library) {
                Some(i) => members[i].push(pair),
                None => {
                    group_libraries.push(library);
                    members.push(vec![pair]);
                }
            }
        }
        let names: Vec<&str> = group_libraries.iter().map(String::as_str).collect();
        for i in hash_map_order(&names) {
            let library = &group_libraries[i];
            if !by_library.contains_key(library) {
                libraries.push(library.clone());
                by_library.insert(library.clone(), LibraryHistograms::default());
            }
            let histograms = by_library.get_mut(library).expect("inserted");
            search_duplicates(&members[i], &opts, histograms);
        }
    }
    let names: Vec<&str> = libraries.iter().map(String::as_str).collect();
    let mut rows = Vec::new();
    let mut series = Vec::new();
    for i in hash_map_order(&names) {
        let library = &libraries[i];
        let histograms = &by_library[library];
        let mut row = Row {
            library: library.clone(),
            ..Row::default()
        };
        for bin in histograms.duplication.keys() {
            let groups = histograms.duplication.get(bin).unwrap_or(0.0);
            let optical = histograms.optical.get(bin).unwrap_or(0.0);
            if groups >= f64::from(opts.min_group_count) {
                row.read_pairs_examined += (bin as f64 * groups) as i64;
                row.read_pair_duplicates += ((bin - 1) as f64 * groups) as i64;
                row.read_pair_optical_duplicates += optical as i64;
            }
        }
        row.estimated_library_size = estimate_library_size(
            row.read_pairs_examined - row.read_pair_optical_duplicates,
            row.read_pairs_examined - row.read_pair_duplicates,
        );
        row.percent_duplication = if row.unpaired_reads_examined + row.read_pairs_examined != 0 {
            (row.unpaired_read_duplicates + row.read_pair_duplicates * 2) as f64
                / (row.unpaired_reads_examined + row.read_pairs_examined * 2) as f64
        } else {
            0.0
        };
        rows.push(row);
        series.push(HistogramSeries {
            label: library.clone(),
            histogram: histograms.duplication.clone(),
        });
    }
    let histogram = (!series.is_empty()).then_some(HistogramSection {
        key_label: "duplication_group_count",
        key_class: KeyClass::Integer,
        series,
    });
    MetricsOutput::from_rows(&rows, histogram)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn options() -> Options {
        Options {
            min_identical_bases: 5,
            max_diff_rate: 0.03,
            min_mean_quality: 20,
            max_group_ratio: 500,
            barcode_tag: None,
            read_one_barcode_tag: None,
            read_two_barcode_tag: None,
            max_read_length: 0,
            min_group_count: 2,
            finder: OpticalDuplicateFinder {
                parser: ReadNameParser::Default,
                pixel_distance: 100,
                max_set_size: DEFAULT_MAX_DUPLICATE_SET_SIZE,
            },
        }
    }

    fn pair(read1: &str, read2: &str, x: i32) -> Complete {
        Complete {
            location: Location {
                read_group: 0,
                tile: 1,
                x,
                y: 0,
            },
            read1: read1.as_bytes().to_vec(),
            read2: read2.as_bytes().to_vec(),
            barcodes: [0; 3],
        }
    }

    #[test]
    fn matches_allow_floor_of_rate_times_compared_bases() {
        let opts = options();
        let a = pair("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT", 0);
        // 40 compared bases, 30 after the seeds: floor(40 * 0.03) = 1 error.
        let mut b = a.clone();
        b.read1[10] = b'A';
        assert!(matches(&a, &b, &opts));
        b.read2[19] = b'A';
        assert!(!matches(&a, &b, &opts));
        // Differences in the first five bases are never counted.
        let mut c = a.clone();
        c.read1[0] = b'T';
        c.read2[4] = b'A';
        assert!(matches(&a, &c, &opts));
        // A shorter mate shortens the comparison.
        let short = pair("ACGTACGTAC", "TTTTTTTTTT", 0);
        assert!(matches(&a, &short, &opts));
        // MAX_READ_LENGTH caps the compared bases: 24 bases allow no error,
        // and the differences past 12 are not seen.
        let capped = Options {
            max_read_length: 12,
            ..options()
        };
        let mut d = a.clone();
        d.read1[15] = b'A';
        d.read2[15] = b'A';
        assert!(!matches(&a, &d, &opts));
        assert!(matches(&a, &d, &capped));
        d.read2[8] = b'A';
        assert!(!matches(&a, &d, &capped));
    }

    #[test]
    fn duplicate_sets_fill_the_histograms() {
        let opts = options();
        let a = pair("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT", 0);
        let b = pair("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT", 50);
        let c = pair("ACGTACGTACGTACGTACGT", "TTTTTTTTTTTTTTTTTTTT", 5000);
        let d = pair("ACGTAGGGGGGGGGGGGGGG", "TTTTTTTTTTTTTTTTTTTT", 0);
        let mut histograms = LibraryHistograms::default();
        search_duplicates(&[&a, &b, &c, &d], &opts, &mut histograms);
        assert_eq!(histograms.duplication.get(3), Some(1.0));
        assert_eq!(histograms.duplication.get(1), Some(1.0));
        // a and b are within 100 pixels; a (the keeper, last in the set)
        // stays, b is the optical duplicate.
        assert_eq!(histograms.optical.get(3), Some(1.0));
    }
}

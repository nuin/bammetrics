//! `CollectRrbsMetrics`, as Picard 3.5.0 computes it (`CollectRrbsMetrics`,
//! `RrbsMetricsCollector`, `CpgLocation`).
//!
//! `METRICS_FILE_PREFIX` is a prefix: a `.` is added unless it ends in one,
//! then `rrbs_summary_metrics` and `rrbs_detail_metrics`. No plot is drawn,
//! so `DO_NOT_CREATE_PLOTS` changes nothing.
//!
//! Every record that is not flagged unmapped counts, on a contig named in
//! `SEQUENCE_NAMES` when the list is not empty: secondary, supplementary,
//! duplicate and QC-failed records count, and each end of a pair counts.
//! The reference contig of a record is taken by its index in the FASTA
//! (`ReferenceSequenceFileWalker.get`), not by its name, and only forward.
//! A counted record adds to `READS_ALIGNED`. It is then ignored as short
//! when it has fewer than `MINIMUM_READ_LENGTH` bases, and for mismatches
//! when its aligned bases hold more mismatches than
//! `Math.round(length x MAX_MISMATCH_RATE)` (`SequenceUtil.countMismatches`
//! with bisulfite equality: a read T on a reference C on the forward
//! strand, a read A on a reference G on the reverse strand, is no
//! mismatch).
//!
//! Each alignment block of a kept record is read on the strand of the
//! read: a reverse-strand block is reverse complemented, reference and
//! read, and its qualities reversed. A reference `CG` is a CpG site. It
//! counts when the read has a C or a T at the C, the reference base at the
//! G, the C at `C_QUALITY_THRESHOLD` or above and the base after it at
//! `NEXT_BASE_QUALITY_THRESHOLD` or above; it is converted when the read
//! has a T. Its `POSITION` is the 0-based position of the forward-strand
//! C. A reference C that starts no CpG counts as a non-CpG base when the
//! read has a C or a T there, the base after it is equal or bisulfite
//! equal, and the qualities pass. Picard tests those qualities on the
//! whole read in file order at the index of the block, not on the block,
//! and this is kept. The last base of a block starts nothing. A kept read
//! with no counted CpG adds to `READS_WITH_NO_CPG`; its non-CpG bases
//! still count.
//!
//! One summary row per accumulation unit (`levels.rs`), empty units
//! included: `MEAN_CPG_COVERAGE` of a unit with no CpG prints `?`.
//! Detail rows follow the units, each unit's sites sorted by contig name
//! (Java `String.compareTo`) then position.
//!
//! Memory: a unit holds the sites of the contig it reads in a hash map and
//! packs them at 3 to 4 bytes a site (`sites.rs`) when it reads another
//! contig. The detail rows are made from the packed sites as the file is
//! written ([`RowSource`]), never held as `Value`s by the Picard writer.
//! A unit never comes back to a contig it left: the reference walker
//! refuses a record on an earlier contig, as Picard's does, even under
//! `ASSUME_SORTED=true`. One contig of the reference is held at a time.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::{bam, sam};

use crate::blocks::alignment_blocks;
use crate::error::{Error, Result};
use crate::java;
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{MetricRow, MetricsOutput, RowSource, ToolOutput, Value};
use crate::pass::{self, Accumulator, Collector, Limits};
use crate::quality::base_qualities;
use crate::record::record_name;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool, picard_list};
use crate::sequence::{bases_equal, bisulfite_bases_equal, complement};

mod sites;

use sites::{Site, Sites};

pub const SUMMARY_EXTENSION: &str = "rrbs_summary_metrics";
pub const DETAIL_EXTENSION: &str = "rrbs_detail_metrics";

/// The options of a parsed `CollectRrbsMetrics` command line.
#[derive(Clone, Debug)]
pub struct Options {
    pub prefix: PathBuf,
    pub reference: PathBuf,
    pub c_quality_threshold: i64,
    pub next_base_quality_threshold: i64,
    pub minimum_read_length: i64,
    pub max_mismatch_rate: f64,
    pub sequence_names: Vec<String>,
    pub levels: Vec<Level>,
    pub assume_sorted: bool,
}

impl Options {
    /// The options, checked as `customCommandLineValidation` checks them:
    /// every failing check is listed.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let int = |name: &str, default: i64| -> Result<i64> {
            inv.get(name).map_or(Ok(default), |text| {
                java::parse_integer(text, 32)
                    .map_err(|e| Error::Usage(format!("{name}={text}: {e}")))
            })
        };
        let max_mismatch_rate = inv.get("MAX_MISMATCH_RATE").map_or(Ok(0.1), |text| {
            java::parse_double(text)
                .map_err(|e| Error::Usage(format!("MAX_MISMATCH_RATE={text}: {e}")))
        })?;
        let required = |name: &str| Error::Usage(format!("{name} is required"));
        let opts = Self {
            prefix: inv
                .output
                .clone()
                .ok_or_else(|| required("METRICS_FILE_PREFIX"))?,
            reference: inv
                .get("REFERENCE")
                .map(PathBuf::from)
                .ok_or_else(|| required("REFERENCE"))?,
            c_quality_threshold: int("C_QUALITY_THRESHOLD", 20)?,
            next_base_quality_threshold: int("NEXT_BASE_QUALITY_THRESHOLD", 10)?,
            minimum_read_length: int("MINIMUM_READ_LENGTH", 5)?,
            max_mismatch_rate,
            sequence_names: picard_list(Vec::new(), inv.all("SEQUENCE_NAMES"), |s| {
                Ok(s.to_owned())
            })?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
            assume_sorted: inv.flag("ASSUME_SORTED", false)?,
        };
        let mut problems = Vec::new();
        if opts.max_mismatch_rate < 0.0 || opts.max_mismatch_rate > 1.0 {
            problems.push("MAX_MISMATCH_RATE must be in the range of 0-1");
        }
        if opts.c_quality_threshold < 0 {
            problems.push("C_QUALITY_THRESHOLD must be >= 0");
        }
        if opts.next_base_quality_threshold < 0 {
            problems.push("NEXT_BASE_QUALITY_THRESHOLD must be >= 0");
        }
        if opts.minimum_read_length <= 0 {
            problems.push("MINIMUM_READ_LENGTH must be > 0");
        }
        if !problems.is_empty() {
            return Err(Error::Usage(problems.join("; ")));
        }
        Ok(opts)
    }

    /// `METRICS_FILE_PREFIX`, a `.` unless it ends in one, and `extension`.
    pub fn path(&self, extension: &str) -> PathBuf {
        file_path(&self.prefix, extension)
    }
}

/// `prefix`, a `.` unless it ends in one, and `extension`.
pub fn file_path(prefix: &Path, extension: &str) -> PathBuf {
    let mut name = OsString::from(prefix.as_os_str());
    if !prefix.as_os_str().as_encoded_bytes().ends_with(b".") {
        name.push(".");
    }
    name.push(extension);
    PathBuf::from(name)
}

/// `CollectRrbsMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectRrbsMetrics",
    summary: "METRICS_FILE_PREFIX is a prefix: PREFIX.rrbs_summary_metrics, one row per \
              accumulation unit, and PREFIX.rrbs_detail_metrics, one row per CpG site and unit; \
              Python returns one row per file: file and rows.",
    options: &[
        Opt::output(&["METRICS_FILE_PREFIX", "M"])
            .with_help("The output prefix of the two files.")
            .required(),
        Opt::path(&["REFERENCE", "R"], "Reference FASTA with a .fai index.").required(),
        Opt::int(
            &["C_QUALITY_THRESHOLD"],
            "Lowest quality of a C for it to count (default 20).",
        ),
        Opt::int(
            &["NEXT_BASE_QUALITY_THRESHOLD"],
            "Lowest quality of the base after a C for the C to count (default 10).",
        ),
        Opt::int(
            &["MINIMUM_READ_LENGTH"],
            "Reads with fewer bases are ignored as short (default 5).",
        ),
        Opt::float(
            &["MAX_MISMATCH_RATE"],
            "Reads with more mismatches than this fraction of their length are ignored \
             (default 0.1).",
        ),
        Opt::text(
            &["SEQUENCE_NAMES"],
            "Count only reads on these contigs; each adds one, null clears (default every \
             contig).",
        )
        .repeated(),
        crate::registry::LEVEL_OPTION,
        Opt::bool(
            &["ASSUME_SORTED", "AS"],
            "true reads a file whose header is not coordinate sorted (default false).",
        ),
        Opt::ignored_bool(&["DO_NOT_CREATE_PLOTS"], "no plot is drawn"),
        Opt::rejected(&["REFERENCE_SEQUENCE"]),
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    // Picard asserts both files writable before it reads.
    crate::registry::check_writable(
        [opts.path(SUMMARY_EXTENSION), opts.path(DETAIL_EXTENSION)]
            .iter()
            .map(PathBuf::as_path),
        &inv.inputs,
    )?;
    let limits = Limits {
        stringency: inv.stringency()?,
        ..Limits::default()
    };
    let files = pass::run_limited(
        &inv.input,
        Some(&opts.reference),
        inv.threads,
        limits,
        |header| build(inv, header),
    )?;
    Ok(files)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    build(inv, header)
}

/// The sort check of `doWork`, the reference walker and the collector.
fn build(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    let opts = Options::from_invocation(inv)?;
    if !opts.assume_sorted {
        pass::check_coordinate_sorted(&inv.input, header)?;
    }
    let rrbs = Rrbs::new(&opts, header)?;
    Ok(pass::boxed(rrbs, move |report: Report| {
        Ok(report.files(&opts))
    }))
}

/// `RrbsSummaryMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct SummaryRow {
    pub reads_aligned: i64,
    pub non_cpg_bases: i64,
    pub non_cpg_converted_bases: i64,
    pub pct_non_cpg_bases_converted: f64,
    pub cpg_bases_seen: i64,
    pub cpg_bases_converted: i64,
    pub pct_cpg_bases_converted: f64,
    pub mean_cpg_coverage: f64,
    pub median_cpg_coverage: i64,
    pub reads_with_no_cpg: i64,
    pub reads_ignored_short: i64,
    pub reads_ignored_mismatches: i64,
    pub labels: Labels,
}

impl MetricRow for SummaryRow {
    const CLASS: &'static str = "picard.analysis.RrbsSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "READS_ALIGNED",
        "NON_CPG_BASES",
        "NON_CPG_CONVERTED_BASES",
        "PCT_NON_CPG_BASES_CONVERTED",
        "CPG_BASES_SEEN",
        "CPG_BASES_CONVERTED",
        "PCT_CPG_BASES_CONVERTED",
        "MEAN_CPG_COVERAGE",
        "MEDIAN_CPG_COVERAGE",
        "READS_WITH_NO_CPG",
        "READS_IGNORED_SHORT",
        "READS_IGNORED_MISMATCHES",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.reads_aligned.into(),
            self.non_cpg_bases.into(),
            self.non_cpg_converted_bases.into(),
            self.pct_non_cpg_bases_converted.into(),
            self.cpg_bases_seen.into(),
            self.cpg_bases_converted.into(),
            self.pct_cpg_bases_converted.into(),
            self.mean_cpg_coverage.into(),
            self.median_cpg_coverage.into(),
            self.reads_with_no_cpg.into(),
            self.reads_ignored_short.into(),
            self.reads_ignored_mismatches.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// `RrbsCpgDetailMetrics`. `PCT_CONVERTED` is 0 when no read was
/// converted, else `CONVERTED_SITES / TOTAL_SITES`.
#[derive(Clone, Debug, PartialEq)]
pub struct DetailRow {
    pub sequence_name: String,
    pub position: i64,
    pub total_sites: i64,
    pub converted_sites: i64,
    pub labels: Labels,
}

/// `PCT_CONVERTED`: 0 when no read was converted.
fn pct_converted(total: i64, converted: i64) -> f64 {
    if converted == 0 {
        0.0
    } else {
        converted as f64 / total as f64
    }
}

impl MetricRow for DetailRow {
    const CLASS: &'static str = "picard.analysis.RrbsCpgDetailMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SEQUENCE_NAME",
        "POSITION",
        "TOTAL_SITES",
        "CONVERTED_SITES",
        "PCT_CONVERTED",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sequence_name.clone().into(),
            self.position.into(),
            self.total_sites.into(),
            self.converted_sites.into(),
            pct_converted(self.total_sites, self.converted_sites).into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// Every row of both files: the summary rows, and the detail rows as the
/// site tables they are made from.
#[derive(Debug, Default)]
pub struct Report {
    pub summary: Vec<SummaryRow>,
    pub details: Details,
}

impl Report {
    /// The two files in the order Picard writes them.
    fn files(self, opts: &Options) -> ToolOutput {
        let details = MetricsOutput {
            source: Some(Arc::new(self.details)),
            ..MetricsOutput::from_rows::<DetailRow>(&[], None)
        };
        ToolOutput {
            main: None,
            extra: vec![
                (
                    opts.path(SUMMARY_EXTENSION),
                    MetricsOutput::from_rows(&self.summary, None),
                ),
                (opts.path(DETAIL_EXTENSION), details),
            ],
        }
    }
}

/// The detail rows: each unit's labels and its sites by contig, contigs
/// in the order of their names. A [`RowSource`]: each row is made as the
/// file is written.
#[derive(Debug, Default)]
pub struct Details {
    names: Vec<String>,
    units: Vec<(Labels, Vec<Contig>)>,
    rows: usize,
}

/// One contig's index and its sites.
type Contig = (u32, Sites);

impl RowSource for Details {
    fn len(&self) -> usize {
        self.rows
    }

    fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
        let mut cells = vec![Value::Null; DetailRow::COLUMNS.len()];
        for (labels, contigs) in &self.units {
            cells[5] = labels.sample.clone().into();
            cells[6] = labels.library.clone().into();
            cells[7] = labels.read_group.clone().into();
            for (id, sites) in contigs {
                cells[0] = self.names[*id as usize].as_str().into();
                for site in sites.iter() {
                    let (total, converted) = (java_int(site.total), java_int(site.converted));
                    cells[1] = i64::from(site.position).into();
                    cells[2] = total.into();
                    cells[3] = converted.into();
                    cells[4] = pct_converted(total, converted).into();
                    row(&cells)?;
                }
            }
        }
        Ok(())
    }
}

/// The reads of a site still being counted.
#[derive(Clone, Copy, Debug, Default)]
struct Count {
    total: u32,
    converted: u32,
}

/// The read and base counters of a unit.
#[derive(Clone, Copy, Debug, Default)]
struct Counters {
    reads_aligned: i64,
    short: i64,
    mismatches: i64,
    no_cpg: i64,
    non_cpg: i64,
    non_cpg_converted: i64,
}

/// One accumulation unit (`PerUnitRrbsMetricsCollector`).
#[derive(Clone, Debug, Default)]
struct Unit {
    counters: Counters,
    /// The sites of `open_contig` by position.
    open: HashMap<u32, Count>,
    /// The contig whose sites `open` holds.
    open_contig: Option<u32>,
    /// The sites of each contig the unit has left, in the order read.
    closed: Vec<Contig>,
}

/// What a record adds to every unit it reaches.
#[derive(Debug, Default, PartialEq)]
struct Observation {
    outcome: Outcome,
    non_cpg: i64,
    non_cpg_converted: i64,
    /// Counted CpG sites: 0-based position of the C and whether it was
    /// converted.
    cpgs: Vec<(u32, bool)>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum Outcome {
    #[default]
    Counted,
    Short,
    Mismatched,
}

impl Unit {
    fn add(&mut self, contig: u32, o: &Observation) {
        let c = &mut self.counters;
        c.reads_aligned += 1;
        match o.outcome {
            Outcome::Short => c.short += 1,
            Outcome::Mismatched => c.mismatches += 1,
            Outcome::Counted => {
                c.non_cpg += o.non_cpg;
                c.non_cpg_converted += o.non_cpg_converted;
                if o.cpgs.is_empty() {
                    c.no_cpg += 1;
                    return;
                }
                self.enter(contig);
                for &(position, converted) in &o.cpgs {
                    let site = self.open.entry(position).or_default();
                    site.total = site.total.saturating_add(1);
                    site.converted = site.converted.saturating_add(u32::from(converted));
                }
            }
        }
    }

    /// Make `contig` the open one: the sites of the contig left go to
    /// `closed`. The reference walker refuses a record on an earlier
    /// contig, so a contig left is never read again.
    fn enter(&mut self, contig: u32) {
        if self.open_contig != Some(contig) {
            self.close();
            self.open_contig = Some(contig);
        }
    }

    /// Move the sites of the open contig to `closed`, packed.
    fn close(&mut self) {
        let Some(contig) = self.open_contig.take() else {
            return;
        };
        let mut sites: Vec<Site> = self
            .open
            .drain()
            .map(|(position, count)| Site {
                position,
                total: count.total,
                converted: count.converted,
            })
            .collect();
        self.closed.push((contig, Sites::pack(&mut sites)));
    }

    /// The counters and every contig's sites.
    fn into_contigs(mut self) -> (Counters, Vec<Contig>) {
        self.close();
        (self.counters, self.closed)
    }
}

/// The thresholds a record is judged by.
#[derive(Clone, Copy, Debug)]
struct Thresholds {
    c_quality: i64,
    next_base_quality: i64,
    minimum_read_length: i64,
    max_mismatch_rate: f64,
}

/// Buffers reused from one record to the next.
#[derive(Default)]
struct Scratch {
    bases: Vec<u8>,
    qualities: Vec<i8>,
    reference: Vec<u8>,
    read: Vec<u8>,
    block_qualities: Vec<i8>,
}

/// The pass accumulator.
pub struct Rrbs {
    thresholds: Thresholds,
    /// `SEQUENCE_NAMES`; empty takes every contig.
    sequence_names: HashSet<String>,
    reference: Reference,
    names: Vec<String>,
    levels: Levels<Unit>,
    scratch: Scratch,
    observation: Observation,
}

impl Rrbs {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        // The walker reads the FASTA by index: no dictionary check.
        let reference = Reference::open(&opts.reference, &sam::Header::default())?;
        Ok(Self {
            thresholds: Thresholds {
                c_quality: opts.c_quality_threshold,
                next_base_quality: opts.next_base_quality_threshold,
                minimum_read_length: opts.minimum_read_length,
                max_mismatch_rate: opts.max_mismatch_rate,
            },
            sequence_names: opts.sequence_names.iter().cloned().collect(),
            reference,
            names: header
                .reference_sequences()
                .keys()
                .map(|name| String::from_utf8_lossy(name).into_owned())
                .collect(),
            levels: Levels::new(header, &opts.levels, Unit::default()),
            scratch: Scratch::default(),
            observation: Observation::default(),
        })
    }
}

impl Accumulator for Rrbs {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if record.flags().is_unmapped() {
            return Ok(());
        }
        let id = record.reference_sequence_id().transpose()?;
        let name = id
            .and_then(|i| self.names.get(i))
            .map_or("*", String::as_str);
        if !self.sequence_names.is_empty() && !self.sequence_names.contains(name) {
            return Ok(());
        }
        let id = id.ok_or_else(|| Error::MalformedRecord {
            name: record_name(record),
            what: "mapped without a reference sequence",
        })?;
        let index = u32::try_from(id).map_err(|_| Error::MalformedRecord {
            name: record_name(record),
            what: "a reference sequence index past 2^32",
        })?;
        let contig = self.reference.contig(id)?;
        observe(
            record,
            contig,
            self.thresholds,
            &mut self.scratch,
            &mut self.observation,
        )?;
        let observation = &self.observation;
        self.levels
            .route(record, |unit| unit.add(index, observation));
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        let mut rank: Vec<usize> = (0..self.names.len()).collect();
        rank.sort_by(|&a, &b| java_compare(&self.names[a], &self.names[b]));
        let mut order = vec![0; self.names.len()];
        for (place, &id) in rank.iter().enumerate() {
            order[id] = place;
        }
        let mut summary = Vec::new();
        let mut units = Vec::new();
        let mut rows = 0;
        for bucket in self.levels.into_buckets() {
            let (counters, mut contigs) = bucket.state.into_contigs();
            contigs.sort_unstable_by_key(|(id, _)| order[*id as usize]);
            summary.push(summary_row(counters, &contigs, bucket.labels.clone()));
            rows += contigs.iter().map(|(_, sites)| sites.len()).sum::<usize>();
            units.push((bucket.labels, contigs));
        }
        Ok(Report {
            summary,
            details: Details {
                names: self.names,
                units,
                rows,
            },
        })
    }
}

/// `finish` and `buildSummaryMetrics` of a unit. The site histogram gives
/// `getSumOfValues` (CpG bases), `getMeanBinSize` (the sum over the number
/// of sites, NaN without one) and `getMedianBinSize` (0 without a site).
fn summary_row(c: Counters, contigs: &[Contig], labels: Labels) -> SummaryRow {
    // Sites by their read count: the sums and the median without a copy
    // of every site.
    let mut by_total: BTreeMap<u32, u64> = BTreeMap::new();
    let (mut seen_sum, mut converted_sum, mut n) = (0u64, 0u64, 0u64);
    for site in contigs.iter().flat_map(|(_, sites)| sites.iter()) {
        *by_total.entry(site.total).or_default() += 1;
        seen_sum += u64::from(site.total);
        converted_sum += u64::from(site.converted);
        n += 1;
    }
    let (seen_sum, converted_sum) = (seen_sum as f64, converted_sum as f64);
    // The `k`th smallest read count, from 0.
    let kth = |k: u64| {
        let mut below = 0;
        by_total
            .iter()
            .find(|&(_, &count)| {
                below += count;
                below > k
            })
            .map_or(0.0, |(&total, _)| f64::from(total))
    };
    let median = match n {
        0 => 0.0,
        n if n % 2 == 0 => (kth(n / 2) + kth(n / 2 - 1)) / 2.0,
        n => kth(n / 2),
    };
    let seen = java_d2i(seen_sum);
    let converted = java_d2i(converted_sum);
    let ratio = |part: i64, whole: i64| {
        if whole == 0 {
            0.0
        } else {
            part as f64 / whole as f64
        }
    };
    SummaryRow {
        reads_aligned: c.reads_aligned,
        non_cpg_bases: c.non_cpg,
        non_cpg_converted_bases: c.non_cpg_converted,
        pct_non_cpg_bases_converted: ratio(c.non_cpg_converted, c.non_cpg),
        cpg_bases_seen: seen,
        cpg_bases_converted: converted,
        pct_cpg_bases_converted: ratio(converted, seen),
        mean_cpg_coverage: seen_sum / n as f64,
        median_cpg_coverage: java_d2i(median),
        reads_with_no_cpg: c.no_cpg,
        reads_ignored_short: c.short,
        reads_ignored_mismatches: c.mismatches,
        labels,
    }
}

/// Java's `(int)` of a double: truncated, saturated at the `int` range.
fn java_d2i(value: f64) -> i64 {
    value as i32 as i64
}

/// A count as the `int` Picard keeps it in.
fn java_int(value: u32) -> i64 {
    i64::from(i32::try_from(value).unwrap_or(i32::MAX))
}

/// `String.compareTo`: UTF-16 code units in order.
fn java_compare(a: &str, b: &str) -> std::cmp::Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

/// `isAboveCytoQcThreshold`: the quality at `index` and the one after it.
fn above_threshold(qualities: &[i8], index: usize, t: Thresholds) -> bool {
    index + 1 < qualities.len()
        && i64::from(qualities[index]) >= t.c_quality
        && i64::from(qualities[index + 1]) >= t.next_base_quality
}

/// `isC`: a reference C and a read C or T.
fn is_c(reference: u8, read: u8) -> bool {
    bases_equal(reference, b'C') && bisulfite_bases_equal(false, read, reference)
}

/// `isBisulfiteConverted` on the forward strand: a read T on a reference C.
fn is_converted(read: u8, reference: u8) -> bool {
    bases_equal(reference, b'C') && bases_equal(read, b'T')
}

/// `PerUnitRrbsMetricsCollector.acceptRecord` up to the unit counters.
fn observe(
    record: &bam::Record,
    contig: &[u8],
    t: Thresholds,
    s: &mut Scratch,
    o: &mut Observation,
) -> Result<()> {
    o.outcome = Outcome::Counted;
    o.non_cpg = 0;
    o.non_cpg_converted = 0;
    o.cpgs.clear();
    s.bases.clear();
    s.bases.extend(record.sequence().iter());
    base_qualities(record, &mut s.qualities);
    let length = s.bases.len();
    if (length as i64) < t.minimum_read_length {
        o.outcome = Outcome::Short;
        return Ok(());
    }
    let malformed = |what| Error::MalformedRecord {
        name: record_name(record),
        what,
    };
    let reverse = record.flags().is_reverse_complemented();
    let blocks = alignment_blocks(record)?;
    let mut mismatches = 0i64;
    for block in &blocks {
        let read_start = block.read_start - 1;
        let reference_start = usize::try_from(block.reference_start - 1)
            .map_err(|_| malformed("alignment before the contig start"))?;
        let read = s
            .bases
            .get(read_start..read_start + block.length)
            .ok_or_else(|| malformed("CIGAR longer than the sequence"))?;
        let reference = contig
            .get(reference_start..reference_start + block.length)
            .ok_or_else(|| malformed("alignment past the end of the contig"))?;
        mismatches += read
            .iter()
            .zip(reference)
            .filter(|&(&r, &f)| !bisulfite_bases_equal(reverse, r, f))
            .count() as i64;
    }
    let allowed = (length as f64 * t.max_mismatch_rate).round() as i64;
    if mismatches > allowed {
        o.outcome = Outcome::Mismatched;
        return Ok(());
    }
    for block in &blocks {
        let len = block.length;
        let read_start = block.read_start - 1;
        let reference_start = (block.reference_start - 1) as usize;
        s.reference.clear();
        s.reference
            .extend_from_slice(&contig[reference_start..reference_start + len]);
        s.read.clear();
        s.read
            .extend_from_slice(&s.bases[read_start..read_start + len]);
        s.block_qualities.clear();
        if s.qualities.is_empty() {
            // `Arrays.copyOfRange` of an empty array: zeros from its start,
            // an exception from any later offset.
            if read_start > 0 {
                return Err(malformed(
                    "no base qualities and an alignment block after the first base",
                ));
            }
            s.block_qualities.resize(len, 0);
        } else {
            let qualities = s
                .qualities
                .get(read_start..read_start + len)
                .ok_or_else(|| malformed("fewer qualities than bases"))?;
            s.block_qualities.extend_from_slice(qualities);
        }
        if reverse {
            reverse_complement_in_place(&mut s.reference);
            reverse_complement_in_place(&mut s.read);
            s.block_qualities.reverse();
        }
        let (fragment, read, qualities) = (&s.reference, &s.read, &s.block_qualities);
        let mut i = 0;
        while i + 1 < len {
            if bases_equal(fragment[i], b'C') && bases_equal(fragment[i + 1], b'G') {
                if is_c(fragment[i], read[i])
                    && bases_equal(fragment[i + 1], read[i + 1])
                    && above_threshold(qualities, i, t)
                {
                    let position = if reverse {
                        block.reference_start - 1 + (len as i64 - 1) - i as i64 - 1
                    } else {
                        block.reference_start - 1 + i as i64
                    };
                    let position = u32::try_from(position)
                        .map_err(|_| malformed("a CpG position past 2^32"))?;
                    o.cpgs.push((position, is_converted(read[i], fragment[i])));
                }
                i += 1;
            } else if is_c(fragment[i], read[i])
                && above_threshold(&s.qualities, i, t)
                && bisulfite_bases_equal(false, read[i + 1], fragment[i + 1])
            {
                o.non_cpg += 1;
                o.non_cpg_converted += i64::from(is_converted(read[i], fragment[i]));
            }
            i += 1;
        }
    }
    Ok(())
}

/// `SequenceUtil.reverseComplement` of a byte array, in place.
fn reverse_complement_in_place(bases: &mut [u8]) {
    bases.reverse();
    for base in bases.iter_mut() {
        *base = complement(*base);
    }
}

#[cfg(test)]
mod tests;

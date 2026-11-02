//! `CollectAlignmentSummaryMetrics` (docs/01 section 2.4), mirroring Picard
//! 3.5.0 `AlignmentSummaryMetricsCollector` at every
//! `METRIC_ACCUMULATION_LEVEL`.
//!
//! Each primary record is described once (`Describer::describe`) and the
//! description is added to its categories in every bucket it belongs to:
//! `FIRST_OF_PAIR` or `SECOND_OF_PAIR` plus `PAIR` for a paired read,
//! `UNPAIRED` otherwise. The read-length histograms come from the
//! `ALL_READS` bucket only; without that level Picard writes none.

mod describe;

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use noodles::{bam, sam};

use crate::adapter::{Adapters, DEFAULT_ADAPTERS};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{HistogramSection, HistogramSeries, MetricRow, MetricsOutput, Value};
use crate::pass::{self, Accumulator};
use crate::record::Orientation;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use describe::{Describer, Read};

/// `ChimeraUtil.DEFAULT_INSERT_SIZE_LIMIT`.
pub const DEFAULT_MAX_INSERT_SIZE: i32 = 100_000;

#[derive(Clone, Debug)]
pub struct Options {
    /// `REFERENCE_SEQUENCE`. Without it no base is a mismatch, so the
    /// mismatch and error rates are 0, as Picard prints them.
    pub reference: Option<PathBuf>,
    /// `MAX_INSERT_SIZE`: a mapped pair with a larger |TLEN| is chimeric.
    pub max_insert_size: i32,
    /// `EXPECTED_PAIR_ORIENTATIONS`: a pair in any other orientation is chimeric.
    pub expected_orientations: Vec<Orientation>,
    /// `ADAPTER_SEQUENCE`.
    pub adapter_sequences: Vec<String>,
    /// `COLLECT_ALIGNMENT_INFORMATION`: false counts reads only (and no-calls
    /// for `BAD_CYCLES`); every alignment field stays 0.
    pub collect_alignment_information: bool,
    /// `IS_BISULFITE_SEQUENCED`: a converted base (read T on reference C,
    /// or A on G for a reverse read) is no mismatch and leaves the mismatch
    /// rate denominators.
    pub bisulfite: bool,
    /// `METRIC_ACCUMULATION_LEVEL`: one set of rows per level.
    pub levels: Vec<Level>,
    /// `STOP_AFTER`: read this many records, every flag counted; 0 or less
    /// reads them all.
    pub stop_after: i64,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            reference: None,
            max_insert_size: DEFAULT_MAX_INSERT_SIZE,
            expected_orientations: vec![Orientation::Fr],
            adapter_sequences: DEFAULT_ADAPTERS.iter().map(|s| (*s).to_owned()).collect(),
            collect_alignment_information: true,
            bisulfite: false,
            levels: vec![Level::AllReads],
            stop_after: 0,
        }
    }
}

impl Options {
    /// The options of a parsed `CollectAlignmentSummaryMetrics` command line.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let defaults = Self::default();
        Ok(Self {
            reference: inv.get("REFERENCE_SEQUENCE").map(PathBuf::from),
            max_insert_size: inv.parse("MAX_INSERT_SIZE", defaults.max_insert_size)?,
            expected_orientations: inv.list(
                "EXPECTED_PAIR_ORIENTATIONS",
                defaults.expected_orientations,
                str::parse,
            )?,
            adapter_sequences: inv.list("ADAPTER_SEQUENCE", defaults.adapter_sequences, |s| {
                Ok(s.to_owned())
            })?,
            collect_alignment_information: inv.flag(
                "COLLECT_ALIGNMENT_INFORMATION",
                defaults.collect_alignment_information,
            )?,
            bisulfite: inv.flag("IS_BISULFITE_SEQUENCED", defaults.bisulfite)?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
            stop_after: inv.parse("STOP_AFTER", defaults.stop_after)?,
        })
    }
}

/// `CollectAlignmentSummaryMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectAlignmentSummaryMetrics",
    options: &[
        Opt::ignored(&["HISTOGRAM_FILE", "H"], "no read-length plot is drawn"),
        Opt::value(&["REFERENCE_SEQUENCE", "R"]),
        Opt::value(&["MAX_INSERT_SIZE"]),
        Opt::value(&["EXPECTED_PAIR_ORIENTATIONS"]),
        Opt::value(&["ADAPTER_SEQUENCE"]),
        Opt::value(&["METRIC_ACCUMULATION_LEVEL", "LEVEL"]),
        Opt::value(&["IS_BISULFITE_SEQUENCED", "BS"]),
        Opt::value(&["COLLECT_ALIGNMENT_INFORMATION"]),
        Opt::value(&["STOP_AFTER"]),
    ],
    run: run_tool,
};

fn run_tool(inv: &Invocation) -> Result<Option<MetricsOutput>> {
    let opts = Options::from_invocation(inv)?;
    let report = pass::run(&inv.input, inv.threads, |header| {
        AlignmentSummary::new(&opts, header)
    })?;
    Ok(Some(report.output()))
}

/// `AlignmentSummaryMetrics.Category`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Category {
    FirstOfPair,
    SecondOfPair,
    Pair,
    Unpaired,
}

impl Category {
    pub fn as_str(self) -> &'static str {
        match self {
            Category::FirstOfPair => "FIRST_OF_PAIR",
            Category::SecondOfPair => "SECOND_OF_PAIR",
            Category::Pair => "PAIR",
            Category::Unpaired => "UNPAIRED",
        }
    }
}

/// `AlignmentSummaryMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 36] = [
    "CATEGORY",
    "TOTAL_READS",
    "PF_READS",
    "PCT_PF_READS",
    "PF_NOISE_READS",
    "PF_READS_ALIGNED",
    "PCT_PF_READS_ALIGNED",
    "PF_ALIGNED_BASES",
    "PF_HQ_ALIGNED_READS",
    "PF_HQ_ALIGNED_BASES",
    "PF_HQ_ALIGNED_Q20_BASES",
    "PF_HQ_MEDIAN_MISMATCHES",
    "PF_MISMATCH_RATE",
    "PF_HQ_ERROR_RATE",
    "PF_INDEL_RATE",
    "MEAN_READ_LENGTH",
    "SD_READ_LENGTH",
    "MEDIAN_READ_LENGTH",
    "MAD_READ_LENGTH",
    "MIN_READ_LENGTH",
    "MAX_READ_LENGTH",
    "MEAN_ALIGNED_READ_LENGTH",
    "READS_ALIGNED_IN_PAIRS",
    "PCT_READS_ALIGNED_IN_PAIRS",
    "PF_READS_IMPROPER_PAIRS",
    "PCT_PF_READS_IMPROPER_PAIRS",
    "BAD_CYCLES",
    "STRAND_BALANCE",
    "PCT_CHIMERAS",
    "PCT_ADAPTER",
    "PCT_SOFTCLIP",
    "PCT_HARDCLIP",
    "AVG_POS_3PRIME_SOFTCLIP_LENGTH",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// One `AlignmentSummaryMetrics` row. Java `long` fields are `i64`, Java
/// `double` fields `f64`; a field Picard leaves unset is 0.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub category: Category,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pct_pf_reads: f64,
    pub pf_noise_reads: i64,
    pub pf_reads_aligned: i64,
    pub pct_pf_reads_aligned: f64,
    pub pf_aligned_bases: i64,
    pub pf_hq_aligned_reads: i64,
    pub pf_hq_aligned_bases: i64,
    pub pf_hq_aligned_q20_bases: i64,
    pub pf_hq_median_mismatches: f64,
    pub pf_mismatch_rate: f64,
    pub pf_hq_error_rate: f64,
    pub pf_indel_rate: f64,
    pub mean_read_length: f64,
    pub sd_read_length: f64,
    pub median_read_length: f64,
    pub mad_read_length: f64,
    pub min_read_length: f64,
    pub max_read_length: f64,
    pub mean_aligned_read_length: f64,
    pub reads_aligned_in_pairs: i64,
    pub pct_reads_aligned_in_pairs: f64,
    pub pf_reads_improper_pairs: i64,
    pub pct_pf_reads_improper_pairs: f64,
    pub bad_cycles: i64,
    pub strand_balance: f64,
    pub pct_chimeras: f64,
    pub pct_adapter: f64,
    pub pct_softclip: f64,
    pub pct_hardclip: f64,
    pub avg_pos_3prime_softclip_length: f64,
    /// The `SAMPLE`, `LIBRARY` and `READ_GROUP` columns.
    pub labels: Labels,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.AlignmentSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &COLUMNS;

    fn values(&self) -> Vec<Value> {
        vec![
            self.category.as_str().into(),
            self.total_reads.into(),
            self.pf_reads.into(),
            self.pct_pf_reads.into(),
            self.pf_noise_reads.into(),
            self.pf_reads_aligned.into(),
            self.pct_pf_reads_aligned.into(),
            self.pf_aligned_bases.into(),
            self.pf_hq_aligned_reads.into(),
            self.pf_hq_aligned_bases.into(),
            self.pf_hq_aligned_q20_bases.into(),
            self.pf_hq_median_mismatches.into(),
            self.pf_mismatch_rate.into(),
            self.pf_hq_error_rate.into(),
            self.pf_indel_rate.into(),
            self.mean_read_length.into(),
            self.sd_read_length.into(),
            self.median_read_length.into(),
            self.mad_read_length.into(),
            self.min_read_length.into(),
            self.max_read_length.into(),
            self.mean_aligned_read_length.into(),
            self.reads_aligned_in_pairs.into(),
            self.pct_reads_aligned_in_pairs.into(),
            self.pf_reads_improper_pairs.into(),
            self.pct_pf_reads_improper_pairs.into(),
            self.bad_cycles.into(),
            self.strand_balance.into(),
            self.pct_chimeras.into(),
            self.pct_adapter.into(),
            self.pct_softclip.into(),
            self.pct_hardclip.into(),
            self.avg_pos_3prime_softclip_length.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

/// The read-length histograms of the `ALL_READS` bucket's `PAIR` and
/// `UNPAIRED` categories.
#[derive(Clone, Debug, PartialEq, Default)]
pub struct ReadLengths {
    /// `PAIRED_TOTAL_LENGTH_COUNT`: read lengths of PF paired reads.
    pub paired_total: Histogram,
    /// `PAIRED_ALIGNED_LENGTH_COUNT`: bases in M, I, = and X operations.
    pub paired_aligned: Histogram,
    /// `UNPAIRED_TOTAL_LENGTH_COUNT`.
    pub unpaired_total: Histogram,
    /// `UNPAIRED_ALIGNED_LENGTH_COUNT`.
    pub unpaired_aligned: Histogram,
}

/// The rows Picard writes and, when `ALL_READS` is a level, the read-length
/// histograms (`CollectAlignmentSummaryMetrics.finish`).
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub rows: Vec<Row>,
    pub histograms: Option<ReadLengths>,
}

impl Report {
    /// The metrics file content. Histograms with no reads are left out when
    /// written (`MetricsFile.printHistogram`).
    pub fn output(&self) -> MetricsOutput {
        let histogram = self.histograms.as_ref().map(|h| {
            let series = |label: &str, histogram: &Histogram| HistogramSeries {
                label: label.to_owned(),
                histogram: histogram.clone(),
            };
            HistogramSection {
                key_label: "READ_LENGTH",
                series: vec![
                    series("PAIRED_TOTAL_LENGTH_COUNT", &h.paired_total),
                    series("PAIRED_ALIGNED_LENGTH_COUNT", &h.paired_aligned),
                    series("UNPAIRED_TOTAL_LENGTH_COUNT", &h.unpaired_total),
                    series("UNPAIRED_ALIGNED_LENGTH_COUNT", &h.unpaired_aligned),
                ],
            }
        });
        MetricsOutput::from_rows(&self.rows, histogram)
    }
}

/// Write a report as a Picard `AlignmentSummaryMetrics` file.
pub fn write_metrics(
    report: &Report,
    headers: Vec<String>,
    out: &mut impl Write,
) -> io::Result<()> {
    report.output().write(headers, out)
}

/// `IndividualAlignmentSummaryMetricsCollector` for one category.
#[derive(Clone, Default)]
struct Collector {
    total_reads: i64,
    pf_reads: i64,
    pf_noise_reads: i64,
    pf_reads_aligned: i64,
    pf_reads_improper_pairs: i64,
    reads_aligned_in_pairs: i64,
    positive_strand: i64,
    chimeras: i64,
    chimeras_denominator: i64,
    adapter_reads: i64,
    soft_clipped: i64,
    three_prime_soft_clipped: i64,
    reads_with_three_prime_soft_clips: i64,
    hard_clipped: i64,
    pf_hq_aligned_reads: i64,
    pf_aligned_bases: i64,
    non_bisulfite_aligned_bases: i64,
    pf_hq_aligned_bases: i64,
    hq_non_bisulfite_aligned_bases: i64,
    pf_hq_aligned_q20_bases: i64,
    indels: i64,
    read_lengths: Histogram,
    aligned_lengths: Histogram,
    mismatches: Histogram,
    hq_mismatches: Histogram,
    bad_cycles: Histogram,
}

impl Collector {
    /// `collectReadData` then `collectQualityData`.
    fn add(&mut self, read: &Read) {
        self.total_reads += 1;
        if !read.qc_fail {
            self.pf_reads += 1;
            self.pf_noise_reads += i64::from(read.noise);
            self.read_lengths.increment(read.length);
            self.aligned_lengths.increment(read.aligned_length);
            self.adapter_reads += i64::from(read.adapter);
            self.hard_clipped += read.hard_clipped;
            if let Some((soft, three_prime)) = read.clips {
                self.soft_clipped += soft;
                if three_prime > 0 {
                    self.three_prime_soft_clipped += three_prime;
                    self.reads_with_three_prime_soft_clips += 1;
                }
            }
            if let Some(a) = &read.alignment {
                self.pf_reads_aligned += 1;
                self.pf_reads_improper_pairs += i64::from(a.improper_pair);
                self.positive_strand += i64::from(a.positive_strand);
                self.reads_aligned_in_pairs += i64::from(a.aligned_in_pair);
                if let Some(chimeric) = a.chimera {
                    self.chimeras_denominator += 1;
                    self.chimeras += i64::from(chimeric);
                }
            }
        }
        for &c in &read.bad_cycles {
            self.bad_cycles.increment(c);
        }
        if let Some(b) = &read.bases {
            self.pf_hq_aligned_reads += i64::from(b.hq_read);
            self.pf_aligned_bases += b.aligned;
            self.non_bisulfite_aligned_bases += b.non_bisulfite;
            self.pf_hq_aligned_bases += b.hq_aligned;
            self.hq_non_bisulfite_aligned_bases += b.hq_non_bisulfite;
            self.pf_hq_aligned_q20_bases += b.hq_q20;
            self.mismatches.increment(b.mismatches);
            self.hq_mismatches.increment(b.hq_mismatches);
            self.indels += b.indels;
        }
    }

    /// `IndividualAlignmentSummaryMetricsCollector.finish`. Ratios go
    /// through `MathUtil.divide` (0 for a zero denominator) except
    /// `PCT_PF_READS` and `PCT_ADAPTER`, which divide directly.
    fn finish(&self, category: Category, labels: &Labels, do_ref_metrics: bool) -> Result<Row> {
        let mut row = Row {
            category,
            total_reads: self.total_reads,
            pf_reads: self.pf_reads,
            pct_pf_reads: 0.0,
            pf_noise_reads: self.pf_noise_reads,
            pf_reads_aligned: self.pf_reads_aligned,
            pct_pf_reads_aligned: 0.0,
            pf_aligned_bases: self.pf_aligned_bases,
            pf_hq_aligned_reads: self.pf_hq_aligned_reads,
            pf_hq_aligned_bases: self.pf_hq_aligned_bases,
            pf_hq_aligned_q20_bases: self.pf_hq_aligned_q20_bases,
            pf_hq_median_mismatches: 0.0,
            pf_mismatch_rate: 0.0,
            pf_hq_error_rate: 0.0,
            pf_indel_rate: 0.0,
            mean_read_length: 0.0,
            sd_read_length: 0.0,
            median_read_length: 0.0,
            mad_read_length: 0.0,
            min_read_length: 0.0,
            max_read_length: 0.0,
            mean_aligned_read_length: 0.0,
            reads_aligned_in_pairs: self.reads_aligned_in_pairs,
            pct_reads_aligned_in_pairs: 0.0,
            pf_reads_improper_pairs: self.pf_reads_improper_pairs,
            pct_pf_reads_improper_pairs: 0.0,
            bad_cycles: 0,
            strand_balance: 0.0,
            pct_chimeras: 0.0,
            pct_adapter: 0.0,
            pct_softclip: 0.0,
            pct_hardclip: 0.0,
            avg_pos_3prime_softclip_length: 0.0,
            labels: labels.clone(),
        };
        if self.pf_reads == 0 {
            if self.total_reads > 0 {
                return Err(Error::Usage(format!(
                    "Input file contains no PF_READS ({} reads in {}, all QC-failed)",
                    self.total_reads,
                    category.as_str()
                )));
            }
            return Ok(row);
        }
        let total = self.total_reads as f64;
        let pf = self.pf_reads as f64;
        let h = &self.read_lengths;
        row.pct_pf_reads = pf / total;
        row.pct_adapter = self.adapter_reads as f64 / pf;
        row.mean_read_length = h.mean();
        row.sd_read_length = h.standard_deviation();
        row.median_read_length = h.median();
        row.mad_read_length = h.median_absolute_deviation();
        row.min_read_length = h.min().unwrap_or_default() as f64;
        row.max_read_length = h.max().unwrap_or_default() as f64;
        // A cycle is bad when no-calls and mismatches there reach 80% of
        // TOTAL_READS (QC-failed reads included in both).
        row.bad_cycles = self
            .bad_cycles
            .keys()
            .filter(|&k| self.bad_cycles.get(k).unwrap_or_default() / total >= 0.8)
            .count() as i64;
        if do_ref_metrics {
            let aligned = self.pf_reads_aligned as f64;
            let total_bases = h.sum();
            row.pct_pf_reads_aligned = divide(aligned, pf);
            row.pct_reads_aligned_in_pairs = divide(self.reads_aligned_in_pairs as f64, aligned);
            row.pct_pf_reads_improper_pairs = divide(self.pf_reads_improper_pairs as f64, aligned);
            row.mean_aligned_read_length = self.aligned_lengths.mean();
            row.strand_balance = divide(self.positive_strand as f64, aligned);
            row.pct_chimeras = divide(self.chimeras as f64, self.chimeras_denominator as f64);
            row.pf_indel_rate = divide(self.indels as f64, self.pf_aligned_bases as f64);
            row.pf_mismatch_rate = divide(
                self.mismatches.sum(),
                self.non_bisulfite_aligned_bases as f64,
            );
            row.pf_hq_error_rate = divide(
                self.hq_mismatches.sum(),
                self.hq_non_bisulfite_aligned_bases as f64,
            );
            row.pct_hardclip = divide(self.hard_clipped as f64, total_bases);
            row.pct_softclip = divide(self.soft_clipped as f64, total_bases);
            row.avg_pos_3prime_softclip_length = divide(
                self.three_prime_soft_clipped as f64,
                self.reads_with_three_prime_soft_clips as f64,
            );
            row.pf_hq_median_mismatches = self.hq_mismatches.median();
        }
        Ok(row)
    }
}

/// Picard `MathUtil.divide`: 0 when the denominator is within 1e-6 of 0.
fn divide(numerator: f64, denominator: f64) -> f64 {
    if (0.0 - denominator).abs() > 0.000_001 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `GroupAlignmentSummaryMetricsPerUnitMetricCollector`: the four
/// categories of one accumulation unit.
#[derive(Clone, Default)]
struct Group {
    first: Collector,
    second: Collector,
    pair: Collector,
    unpaired: Collector,
}

impl Group {
    fn add(&mut self, read: &Read, paired: bool, first: bool) {
        if paired {
            if first {
                self.first.add(read);
            } else {
                self.second.add(read);
            }
            self.pair.add(read);
        } else {
            self.unpaired.add(read);
        }
    }

    /// `addMetricsToFile`: paired rows only when FIRST_OF_PAIR has reads,
    /// with PAIR's BAD_CYCLES the sum of the two ends; UNPAIRED when it has
    /// reads or when FIRST_OF_PAIR has none (an all-zero row for a unit
    /// with no reads).
    fn rows(&self, labels: &Labels, do_ref: bool) -> Result<Vec<Row>> {
        let unpaired_row = self.unpaired.finish(Category::Unpaired, labels, do_ref)?;
        let first_row = self.first.finish(Category::FirstOfPair, labels, do_ref)?;
        let second_row = self.second.finish(Category::SecondOfPair, labels, do_ref)?;
        let mut pair_row = self.pair.finish(Category::Pair, labels, do_ref)?;
        let mut rows = Vec::with_capacity(4);
        let has_pairs = first_row.total_reads > 0;
        if has_pairs {
            pair_row.bad_cycles = first_row.bad_cycles + second_row.bad_cycles;
            rows.extend([first_row, second_row, pair_row]);
        }
        if unpaired_row.total_reads > 0 || !has_pairs {
            rows.push(unpaired_row);
        }
        Ok(rows)
    }
}

/// `CollectAlignmentSummaryMetrics` as a pass [`Accumulator`]: one group of
/// category collectors per accumulation unit, plus the reference walker
/// when `R` is given.
pub struct AlignmentSummary {
    reference: Option<Reference>,
    describer: Describer,
    levels: Levels<Group>,
    stop_after: Option<u64>,
    records: u64,
}

impl AlignmentSummary {
    /// Opens the reference, when there is one, and checks it against the
    /// BAM header.
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        let reference = opts
            .reference
            .as_deref()
            .map(|r| Reference::open(r, header))
            .transpose()?;
        Ok(Self {
            reference,
            describer: Describer {
                adapters: Adapters::new(&opts.adapter_sequences),
                max_insert_size: i64::from(opts.max_insert_size),
                expected: opts.expected_orientations.clone(),
                do_ref_metrics: opts.collect_alignment_information,
                bisulfite: opts.bisulfite,
                bases: Vec::new(),
                cigar: Vec::new(),
            },
            levels: Levels::new(header, &opts.levels, Group::default()),
            stop_after: u64::try_from(opts.stop_after).ok().filter(|&n| n > 0),
            records: 0,
        })
    }
}

impl Accumulator for AlignmentSummary {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if self.done() {
            return Ok(());
        }
        // `SinglePassSamProgram` counts every record toward STOP_AFTER and
        // fetches the contig of every placed one, secondary and unmapped
        // ones included, so the sort check covers them.
        self.records += 1;
        let bases = match (
            &mut self.reference,
            record.reference_sequence_id().transpose()?,
        ) {
            (Some(r), Some(id)) => Some(r.contig(id)?),
            _ => None,
        };
        let flags = record.flags();
        // `AlignmentSummaryMetricsCollector.acceptRecord` drops secondary and
        // supplementary records before any counting.
        if flags.is_secondary() || flags.is_supplementary() {
            return Ok(());
        }
        let read = self.describer.describe(record, bases)?;
        let (paired, first) = (flags.is_segmented(), flags.is_first_segment());
        self.levels
            .route(record, |group| group.add(&read, paired, first));
        Ok(())
    }

    fn done(&self) -> bool {
        self.stop_after.is_some_and(|n| self.records >= n)
    }

    /// Every collector is finished before any row is kept, so a unit with
    /// QC-failed reads only fails the run as Picard's `finish` does.
    fn finish(self) -> Result<Report> {
        let do_ref = self.describer.do_ref_metrics;
        let mut rows = Vec::new();
        let mut histograms = None;
        for bucket in self.levels.into_buckets() {
            rows.extend(bucket.state.rows(&bucket.labels, do_ref)?);
            if bucket.level == Level::AllReads {
                let group = bucket.state;
                histograms = Some(ReadLengths {
                    paired_total: group.pair.read_lengths,
                    paired_aligned: group.pair.aligned_lengths,
                    unpaired_total: group.unpaired.read_lengths,
                    unpaired_aligned: group.unpaired.aligned_lengths,
                });
            }
        }
        Ok(Report { rows, histograms })
    }
}

/// `StringHeader` command line in Picard's shape.
pub fn command_header(opts: &Options, bam: &Path, output: Option<&Path>) -> String {
    let mut command = format!("CollectAlignmentSummaryMetrics INPUT={}", bam.display());
    if let Some(output) = output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(reference) = &opts.reference {
        command.push_str(&format!(" REFERENCE_SEQUENCE={}", reference.display()));
    }
    let orientations: Vec<&str> = opts
        .expected_orientations
        .iter()
        .map(|o| o.as_str())
        .collect();
    let levels: Vec<&str> = opts.levels.iter().map(|l| l.name()).collect();
    command.push_str(&format!(
        " MAX_INSERT_SIZE={} EXPECTED_PAIR_ORIENTATIONS=[{}] ADAPTER_SEQUENCE=[{}] METRIC_ACCUMULATION_LEVEL=[{}] IS_BISULFITE_SEQUENCED={} COLLECT_ALIGNMENT_INFORMATION={}",
        opts.max_insert_size,
        orientations.join(", "),
        opts.adapter_sequences.join(", "),
        levels.join(", "),
        opts.bisulfite,
        opts.collect_alignment_information
    ));
    if opts.stop_after > 0 {
        command.push_str(&format!(" STOP_AFTER={}", opts.stop_after));
    }
    command
}

/// Compute the metrics of one BAM with the default BGZF worker count.
pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Report> {
    pass::run(path.as_ref(), pass::default_threads(), |header| {
        AlignmentSummary::new(opts, header)
    })
}

#[cfg(test)]
mod tests {
    use super::describe::BaseCounts;
    use super::*;

    #[test]
    fn math_util_divide_guards_zero() {
        assert_eq!(divide(3.0, 0.0), 0.0);
        assert_eq!(divide(3.0, 0.000_000_5), 0.0);
        assert_eq!(divide(3.0, 4.0), 0.75);
    }

    #[test]
    fn empty_unit_is_an_all_zero_unpaired_row() {
        let labels = Labels {
            sample: Some("S1".to_owned()),
            library: Some("L1".to_owned()),
            read_group: None,
        };
        let rows = Group::default()
            .rows(&labels, true)
            .expect("no reads is not an error");
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].category, Category::Unpaired);
        assert_eq!(rows[0].mean_read_length, 0.0);
        let mut out = Vec::new();
        let report = Report {
            rows,
            histograms: Some(ReadLengths::default()),
        };
        write_metrics(&report, Vec::new(), &mut out).expect("write to Vec");
        let text = String::from_utf8(out).expect("utf8");
        assert!(text.ends_with("\nUNPAIRED\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tS1\tL1\t\n\n\n"));
    }

    #[test]
    fn qc_failed_only_category_is_an_error() {
        let mut collector = Collector::default();
        collector.add(&Read {
            qc_fail: true,
            ..Read::default()
        });
        assert!(matches!(
            collector.finish(Category::Unpaired, &Labels::default(), true),
            Err(Error::Usage(m)) if m.contains("no PF_READS")
        ));
    }

    #[test]
    fn bad_cycle_needs_eighty_percent_of_all_reads() {
        let mut collector = Collector::default();
        for i in 0..5 {
            collector.add(&Read {
                length: 10,
                bad_cycles: if i < 4 { vec![3, 7] } else { vec![3] },
                ..Read::default()
            });
        }
        let row = collector
            .finish(Category::Unpaired, &Labels::default(), true)
            .expect("finish");
        // cycle 3: 5/5, cycle 7: 4/5 = 0.8
        assert_eq!(row.bad_cycles, 2);
    }

    #[test]
    fn bisulfite_bases_leave_the_rate_denominators() {
        let mut collector = Collector::default();
        collector.add(&Read {
            length: 10,
            bases: Some(BaseCounts {
                hq_read: true,
                aligned: 10,
                non_bisulfite: 8,
                hq_aligned: 10,
                hq_non_bisulfite: 8,
                mismatches: 2,
                hq_mismatches: 2,
                ..BaseCounts::default()
            }),
            ..Read::default()
        });
        let row = collector
            .finish(Category::Unpaired, &Labels::default(), true)
            .expect("finish");
        assert_eq!(row.pf_aligned_bases, 10);
        assert_eq!(row.pf_mismatch_rate, 0.25);
        assert_eq!(row.pf_hq_error_rate, 0.25);
    }
}

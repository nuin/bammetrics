//! `CollectRnaSeqMetrics`: Picard's `RnaSeqMetricsCollector` over the
//! refFlat model of [`crate::ref_flat`] (M17), with ribosomal intervals,
//! strand specificity, every accumulation level and the 101-row normalized
//! coverage histogram.
//!
//! Picard 3.5.0 behaviour (`picard.analysis.directed.RnaSeqMetricsCollector`):
//! * A QC-failed record counts nowhere. `PF_BASES` adds the base count of
//!   every other record that is not secondary, unmapped and supplementary
//!   ones included. A mapped non-secondary record on an `IGNORE_SEQUENCE`
//!   contig adds one to `IGNORED_READS` and nothing else. Secondary and
//!   unmapped records stop there.
//! * The fragment of a paired record whose mate is mapped on the same contig
//!   runs from the lower of the two starts over `|TLEN|` bases; an unpaired
//!   record's fragment is its alignment; a record with an unmapped or
//!   cross-contig mate has none. When the longest single ribosomal interval
//!   overlap over the fragment length reaches `RRNA_FRAGMENT_PERCENTAGE`,
//!   every aligned base of the record is ribosomal and the record counts in
//!   no other column (a 0-length fragment divides to NaN and is never
//!   ribosomal). Without `RIBOSOMAL_INTERVALS`, `RIBOSOMAL_BASES` and
//!   `PCT_RIBOSOMAL_BASES` print empty.
//! * Each aligned base takes the strongest `LocusFunction` over every
//!   transcript of every gene overlapping the alignment span;
//!   `PF_ALIGNED_BASES` counts them all. Duplicates and MAPQ are not
//!   filtered. The strand and template columns use a non-supplementary
//!   record with a coding or UTR base that overlaps exactly one gene;
//!   `NUM_R1_`, `NUM_R2_` and `NUM_UNEXPLAINED_READS` count first-of-pair
//!   and unpaired records whose template (mate end from `MC`, else the
//!   read length) is FR and lies within the gene.
//! * Coverage per transcript counts every exonic base of every alignment
//!   block except the block's last (`ref_flat.rs`). The coverage columns
//!   use the best transcript of each gene with coverage: the highest mean
//!   among transcripts at least `max(MINIMUM_LENGTH, END_BIAS_BASES)` long
//!   with a mean of at least 1, ties to the first in `HashMap` order, and
//!   then the transcripts whose mean reaches the 1001st-highest value.
//!   Coverage is reversed for a negative-strand gene. CV is
//!   `sqrt(mean of squares - mean^2) / mean`; the end biases are the means
//!   of the first and last `END_BIAS_BASES` over the mean, their ratio 0
//!   when the 3' mean is within 1e-6 of 0. The medians are over
//!   `Histogram<Double>` bins (NaN sorts last, as `Double.compareTo`).
//!   Each percent bin of the histogram averages the bases within 0.5% of
//!   the position, over the mean, over the transcript count, summed over
//!   the picked transcripts in `HashMap` iteration order ([`crate::java`]).
//! * Rows and histogram columns follow [`crate::levels`]; a unit with no
//!   picked transcript has an empty histogram, which the file leaves out.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam;

use crate::blocks::{AlignmentBlock, alignment_blocks, alignment_end};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::intervals::{self, OverlapIndex};
use crate::java::hash_map_order_by_hash;
use crate::levels::{Level, Levels, parse_levels};
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator, Collector};
use crate::record::{Orientation, mate_reference_length, orientation, position};
use crate::ref_flat::{GeneModel, LocusFunction};
use crate::registry::{Invocation, Opt, Tool};

/// `RnaSeqMetricsCollector.StrandSpecificity`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    None,
    FirstReadTranscriptionStrand,
    SecondReadTranscriptionStrand,
}

impl Strand {
    fn parse(text: &str) -> Result<Self> {
        match text {
            "NONE" => Ok(Strand::None),
            "FIRST_READ_TRANSCRIPTION_STRAND" => Ok(Strand::FirstReadTranscriptionStrand),
            "SECOND_READ_TRANSCRIPTION_STRAND" => Ok(Strand::SecondReadTranscriptionStrand),
            _ => Err(Error::Usage(format!(
                "'{text}' is not a valid value for StrandSpecificity."
            ))),
        }
    }
}

pub const DEFAULT_MINIMUM_LENGTH: i64 = 500;
pub const DEFAULT_RRNA_FRAGMENT_PERCENTAGE: f64 = 0.8;
pub const DEFAULT_END_BIAS_BASES: i64 = 100;

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub ref_flat: PathBuf,
    pub ribosomal_intervals: Option<PathBuf>,
    pub strand: Strand,
    pub minimum_length: i64,
    pub ignore_sequence: Vec<String>,
    pub rrna_fragment_percentage: f64,
    pub levels: Vec<Level>,
    pub end_bias_bases: i64,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let required = |name: &str| {
            inv.get(name)
                .map(str::to_owned)
                .ok_or_else(|| Error::Usage(format!("Option '{name}' is required.")))
        };
        let ref_flat = PathBuf::from(required("REF_FLAT")?);
        let strand = Strand::parse(&required("STRAND_SPECIFICITY")?)?;
        let opts = Self {
            ref_flat,
            ribosomal_intervals: inv.get("RIBOSOMAL_INTERVALS").map(PathBuf::from),
            strand,
            minimum_length: inv.parse("MINIMUM_LENGTH", DEFAULT_MINIMUM_LENGTH)?,
            ignore_sequence: inv.list("IGNORE_SEQUENCE", Vec::new(), |v| Ok(v.to_owned()))?,
            rrna_fragment_percentage: inv
                .parse("RRNA_FRAGMENT_PERCENTAGE", DEFAULT_RRNA_FRAGMENT_PERCENTAGE)?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
            end_bias_bases: inv.parse("END_BIAS_BASES", DEFAULT_END_BIAS_BASES)?,
        };
        if opts.ribosomal_intervals.is_none() && opts.rrna_fragment_percentage == 0.0 {
            return Err(Error::Usage(
                "Must use a RIBOSOMAL_INTERVALS file if RRNA_FRAGMENT_PERCENTAGE = 0.0".to_owned(),
            ));
        }
        Ok(opts)
    }
}

/// `CollectRnaSeqMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectRnaSeqMetrics",
    options: &[
        Opt::value(&["REF_FLAT"]),
        Opt::value(&["RIBOSOMAL_INTERVALS"]),
        Opt::value(&["STRAND_SPECIFICITY", "STRAND"]),
        Opt::value(&["MINIMUM_LENGTH"]),
        Opt::ignored(&["CHART_OUTPUT", "CHART"], "no coverage plot is drawn"),
        Opt::value(&["IGNORE_SEQUENCE"]),
        Opt::value(&["RRNA_FRAGMENT_PERCENTAGE"]),
        Opt::value(&["METRIC_ACCUMULATION_LEVEL", "LEVEL"]),
        Opt::value(&["END_BIAS_BASES"]),
        Opt::value(&["STOP_AFTER"]),
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    let opts = Options::from_invocation(inv)?;
    Ok(pass::boxed(
        RnaSeqMetrics::new(&opts, header, &inv.input)?,
        |report| Ok(report.into()),
    ))
}

/// `RnaSeqMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 30] = [
    "PF_BASES",
    "PF_ALIGNED_BASES",
    "RIBOSOMAL_BASES",
    "CODING_BASES",
    "UTR_BASES",
    "INTRONIC_BASES",
    "INTERGENIC_BASES",
    "IGNORED_READS",
    "CORRECT_STRAND_READS",
    "INCORRECT_STRAND_READS",
    "NUM_R1_TRANSCRIPT_STRAND_READS",
    "NUM_R2_TRANSCRIPT_STRAND_READS",
    "NUM_UNEXPLAINED_READS",
    "PCT_R1_TRANSCRIPT_STRAND_READS",
    "PCT_R2_TRANSCRIPT_STRAND_READS",
    "PCT_RIBOSOMAL_BASES",
    "PCT_CODING_BASES",
    "PCT_UTR_BASES",
    "PCT_INTRONIC_BASES",
    "PCT_INTERGENIC_BASES",
    "PCT_MRNA_BASES",
    "PCT_USABLE_BASES",
    "PCT_CORRECT_STRAND_READS",
    "MEDIAN_CV_COVERAGE",
    "MEDIAN_5PRIME_BIAS",
    "MEDIAN_3PRIME_BIAS",
    "MEDIAN_5PRIME_TO_3PRIME_BIAS",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// One `RnaSeqMetrics` row.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub pf_bases: i64,
    pub pf_aligned_bases: i64,
    pub ribosomal_bases: Option<i64>,
    pub coding_bases: i64,
    pub utr_bases: i64,
    pub intronic_bases: i64,
    pub intergenic_bases: i64,
    pub ignored_reads: i64,
    pub correct_strand_reads: i64,
    pub incorrect_strand_reads: i64,
    pub num_r1_transcript_strand_reads: i64,
    pub num_r2_transcript_strand_reads: i64,
    pub num_unexplained_reads: i64,
    pub pct_r1_transcript_strand_reads: f64,
    pub pct_r2_transcript_strand_reads: f64,
    pub pct_ribosomal_bases: Option<f64>,
    pub pct_coding_bases: f64,
    pub pct_utr_bases: f64,
    pub pct_intronic_bases: f64,
    pub pct_intergenic_bases: f64,
    pub pct_mrna_bases: f64,
    pub pct_usable_bases: f64,
    pub pct_correct_strand_reads: f64,
    pub median_cv_coverage: f64,
    pub median_5prime_bias: f64,
    pub median_3prime_bias: f64,
    pub median_5prime_to_3prime_bias: f64,
    pub sample: Option<String>,
    pub library: Option<String>,
    pub read_group: Option<String>,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.RnaSeqMetrics";
    const COLUMNS: &'static [&'static str] = &COLUMNS;

    fn values(&self) -> Vec<Value> {
        vec![
            self.pf_bases.into(),
            self.pf_aligned_bases.into(),
            self.ribosomal_bases.into(),
            self.coding_bases.into(),
            self.utr_bases.into(),
            self.intronic_bases.into(),
            self.intergenic_bases.into(),
            self.ignored_reads.into(),
            self.correct_strand_reads.into(),
            self.incorrect_strand_reads.into(),
            self.num_r1_transcript_strand_reads.into(),
            self.num_r2_transcript_strand_reads.into(),
            self.num_unexplained_reads.into(),
            self.pct_r1_transcript_strand_reads.into(),
            self.pct_r2_transcript_strand_reads.into(),
            self.pct_ribosomal_bases.into(),
            self.pct_coding_bases.into(),
            self.pct_utr_bases.into(),
            self.pct_intronic_bases.into(),
            self.pct_intergenic_bases.into(),
            self.pct_mrna_bases.into(),
            self.pct_usable_bases.into(),
            self.pct_correct_strand_reads.into(),
            self.median_cv_coverage.into(),
            self.median_5prime_bias.into(),
            self.median_3prime_bias.into(),
            self.median_5prime_to_3prime_bias.into(),
            self.sample.clone().into(),
            self.library.clone().into(),
            self.read_group.clone().into(),
        ]
    }
}

/// A transcript: gene index and transcript index in the model.
type TranscriptKey = (usize, usize);

/// `PerUnitRnaSeqMetricsCollector`: the counters and the coverage array of
/// every transcript of every gene a record overlapped, in first-use order
/// (the insertion order of Picard's `coverageByTranscript` map).
#[derive(Clone, Debug)]
struct Unit {
    pf_bases: i64,
    pf_aligned_bases: i64,
    ribosomal_bases: Option<i64>,
    coding_bases: i64,
    utr_bases: i64,
    intronic_bases: i64,
    intergenic_bases: i64,
    ignored_reads: i64,
    correct_strand_reads: i64,
    incorrect_strand_reads: i64,
    r1_reads: i64,
    r2_reads: i64,
    unexplained_reads: i64,
    coverage_index: HashMap<TranscriptKey, usize>,
    coverage: Vec<(TranscriptKey, Vec<i32>)>,
}

impl Unit {
    fn new(ribosomal_bases: Option<i64>) -> Self {
        Self {
            pf_bases: 0,
            pf_aligned_bases: 0,
            ribosomal_bases,
            coding_bases: 0,
            utr_bases: 0,
            intronic_bases: 0,
            intergenic_bases: 0,
            ignored_reads: 0,
            correct_strand_reads: 0,
            incorrect_strand_reads: 0,
            r1_reads: 0,
            r2_reads: 0,
            unexplained_reads: 0,
            coverage_index: HashMap::new(),
            coverage: Vec::new(),
        }
    }

    fn coverage_of(&mut self, key: TranscriptKey, length: usize) -> &mut Vec<i32> {
        let index = *self.coverage_index.entry(key).or_insert_with(|| {
            self.coverage.push((key, vec![0; length]));
            self.coverage.len() - 1
        });
        &mut self.coverage[index].1
    }

    /// `acceptRecord` after the filters, for one unit.
    fn tally(&mut self, seen: &Seen, model: &GeneModel) -> Result<()> {
        self.pf_bases += seen.read_length;
        match &seen.class {
            Class::Dropped => {}
            Class::Ignored => self.ignored_reads += 1,
            Class::Ribosomal { aligned_bases } => {
                let Some(ribosomal) = &mut self.ribosomal_bases else {
                    return Err(Error::Metric(
                        "a fragment reached RRNA_FRAGMENT_PERCENTAGE without RIBOSOMAL_INTERVALS (Picard throws a NullPointerException)".to_owned(),
                    ));
                };
                *ribosomal += aligned_bases;
                self.pf_aligned_bases += aligned_bases;
            }
            Class::Aligned(aligned) => {
                for (block, functions) in aligned.blocks.iter().zip(&aligned.functions) {
                    for &gene in &aligned.genes {
                        let gene_model = &model.genes()[gene];
                        for (t, transcript) in gene_model.transcripts.iter().enumerate() {
                            let coverage = self.coverage_of((gene, t), transcript.length as usize);
                            transcript.add_coverage(
                                block.reference_start,
                                block.reference_end(),
                                coverage,
                            );
                        }
                    }
                    for function in functions {
                        self.pf_aligned_bases += 1;
                        match function {
                            LocusFunction::Intergenic => self.intergenic_bases += 1,
                            LocusFunction::Intronic => self.intronic_bases += 1,
                            LocusFunction::Utr => self.utr_bases += 1,
                            LocusFunction::Coding => self.coding_bases += 1,
                            LocusFunction::Ribosomal => {
                                if let Some(ribosomal) = &mut self.ribosomal_bases {
                                    *ribosomal += 1;
                                }
                            }
                        }
                    }
                }
                match aligned.strand {
                    Some(true) => self.correct_strand_reads += 1,
                    Some(false) => self.incorrect_strand_reads += 1,
                    None => {}
                }
                match aligned.template {
                    Some(Template::R1) => self.r1_reads += 1,
                    Some(Template::R2) => self.r2_reads += 1,
                    Some(Template::Unexplained) => self.unexplained_reads += 1,
                    None => {}
                }
            }
        }
        Ok(())
    }

    /// `finish` plus `computeCoverageMetrics`: the row and the unit's
    /// normalized coverage histogram.
    fn finish(
        self,
        model: &GeneModel,
        opts: &Options,
        labels: crate::levels::Labels,
    ) -> (Row, Histogram) {
        let mut row = Row {
            pf_bases: self.pf_bases,
            pf_aligned_bases: self.pf_aligned_bases,
            ribosomal_bases: self.ribosomal_bases,
            coding_bases: self.coding_bases,
            utr_bases: self.utr_bases,
            intronic_bases: self.intronic_bases,
            intergenic_bases: self.intergenic_bases,
            ignored_reads: self.ignored_reads,
            correct_strand_reads: self.correct_strand_reads,
            incorrect_strand_reads: self.incorrect_strand_reads,
            num_r1_transcript_strand_reads: self.r1_reads,
            num_r2_transcript_strand_reads: self.r2_reads,
            num_unexplained_reads: self.unexplained_reads,
            pct_r1_transcript_strand_reads: 0.0,
            pct_r2_transcript_strand_reads: 0.0,
            pct_ribosomal_bases: None,
            pct_coding_bases: 0.0,
            pct_utr_bases: 0.0,
            pct_intronic_bases: 0.0,
            pct_intergenic_bases: 0.0,
            pct_mrna_bases: 0.0,
            pct_usable_bases: 0.0,
            pct_correct_strand_reads: 0.0,
            median_cv_coverage: 0.0,
            median_5prime_bias: 0.0,
            median_3prime_bias: 0.0,
            median_5prime_to_3prime_bias: 0.0,
            sample: labels.sample,
            library: labels.library,
            read_group: labels.read_group,
        };
        if self.pf_aligned_bases > 0 {
            let aligned = self.pf_aligned_bases as f64;
            row.pct_ribosomal_bases = self.ribosomal_bases.map(|r| r as f64 / aligned);
            row.pct_coding_bases = self.coding_bases as f64 / aligned;
            row.pct_utr_bases = self.utr_bases as f64 / aligned;
            row.pct_intronic_bases = self.intronic_bases as f64 / aligned;
            row.pct_intergenic_bases = self.intergenic_bases as f64 / aligned;
            row.pct_mrna_bases = row.pct_coding_bases + row.pct_utr_bases;
            row.pct_usable_bases =
                (self.coding_bases + self.utr_bases) as f64 / self.pf_bases as f64;
        }
        if self.correct_strand_reads > 0 || self.incorrect_strand_reads > 0 {
            row.pct_correct_strand_reads = self.correct_strand_reads as f64
                / (self.correct_strand_reads + self.incorrect_strand_reads) as f64;
        }
        let examined = self.r1_reads + self.r2_reads;
        if examined > 0 {
            row.pct_r1_transcript_strand_reads = self.r1_reads as f64 / examined as f64;
            row.pct_r2_transcript_strand_reads = self.r2_reads as f64 / examined as f64;
        }
        let picked = pick_transcripts(&self.coverage, model, opts);
        let coverage_metrics = coverage_metrics(&picked, &self.coverage, model, opts);
        row.median_cv_coverage = coverage_metrics.median_cv;
        row.median_5prime_bias = coverage_metrics.median_5prime;
        row.median_3prime_bias = coverage_metrics.median_3prime;
        row.median_5prime_to_3prime_bias = coverage_metrics.median_5prime_to_3prime;
        (row, coverage_metrics.histogram)
    }
}

/// `pickTranscripts`: indices into `coverage` of the transcripts the
/// coverage metrics use, in the iteration order of the map Picard builds.
/// Each `HashMap` and `HashSet` on the way is iterated in Java order, so
/// the histogram sums run in Picard's order.
fn pick_transcripts(
    coverage: &[(TranscriptKey, Vec<i32>)],
    model: &GeneModel,
    opts: &Options,
) -> Vec<usize> {
    let transcript_hash = |key: TranscriptKey| model.genes()[key.0].transcripts[key.1].hash;
    // coverageByTranscript.keySet(), then the HashSet of their genes.
    let key_hashes: Vec<i32> = coverage.iter().map(|(k, _)| transcript_hash(*k)).collect();
    let mut genes: Vec<usize> = Vec::new();
    for i in hash_map_order_by_hash(&key_hashes) {
        let gene = coverage[i].0.0;
        if !genes.contains(&gene) {
            genes.push(gene);
        }
    }
    let gene_hashes: Vec<i32> = genes.iter().map(|&g| model.genes()[g].hash).collect();
    let by_key: HashMap<TranscriptKey, usize> = coverage
        .iter()
        .enumerate()
        .map(|(i, (k, _))| (*k, i))
        .collect();
    let minimum = opts.minimum_length.max(opts.end_bias_bases);
    // bestPerGene, in the gene set's order.
    let mut best_per_gene: Vec<(usize, f64)> = Vec::new();
    for g in hash_map_order_by_hash(&gene_hashes) {
        let gene = genes[g];
        let mut best: Option<(usize, f64)> = None;
        for (t, transcript) in model.genes()[gene].transcripts.iter().enumerate() {
            let Some(&index) = by_key.get(&(gene, t)) else {
                continue;
            };
            if transcript.length < minimum {
                continue;
            }
            let mean = mean(&coverage[index].1);
            if mean < 1.0 {
                continue;
            }
            if best.is_none_or(|(_, best_mean)| mean > best_mean) {
                best = Some((index, mean));
            }
        }
        if let Some(best) = best {
            best_per_gene.push(best);
        }
    }
    let best_hashes: Vec<i32> = best_per_gene
        .iter()
        .map(|&(i, _)| transcript_hash(coverage[i].0))
        .collect();
    let mut means: Vec<f64> = best_per_gene.iter().map(|&(_, m)| m).collect();
    means.sort_by(f64::total_cmp);
    let min = if means.is_empty() {
        0.0
    } else {
        means[means.len().saturating_sub(1001)]
    };
    let picked: Vec<usize> = hash_map_order_by_hash(&best_hashes)
        .into_iter()
        .filter(|&b| best_per_gene[b].1 >= min)
        .map(|b| best_per_gene[b].0)
        .collect();
    let picked_hashes: Vec<i32> = picked
        .iter()
        .map(|&i| transcript_hash(coverage[i].0))
        .collect();
    hash_map_order_by_hash(&picked_hashes)
        .into_iter()
        .map(|p| picked[p])
        .collect()
}

/// `MathUtil.mean(promote(counts), 0, length)`.
fn mean(counts: &[i32]) -> f64 {
    let total: f64 = counts.iter().fold(0.0, |t, &c| t + f64::from(c));
    total / counts.len() as f64
}

/// `MathUtil.mean(in, start, stop)` over an array of doubles; the bounds
/// may run outside the array for a negative `END_BIAS_BASES`, where Java
/// sums nothing and divides by the negative length.
fn mean_range(values: &[f64], start: i64, stop: i64) -> f64 {
    let mut total = 0.0;
    let mut i = start;
    while i < stop {
        total += values[i as usize];
        i += 1;
    }
    total / (stop - start) as f64
}

/// `MathUtil.stddev(in, mean)`: the population form from the mean of
/// squares.
fn stddev(values: &[f64], mean: f64) -> f64 {
    let total: f64 = values.iter().fold(0.0, |t, &v| t + v * v);
    ((total / values.len() as f64) - (mean * mean)).sqrt()
}

/// `MathUtil.divide`: 0 when the denominator is within 1e-6 of 0.
fn divide(numerator: f64, denominator: f64) -> f64 {
    if (0.0 - denominator).abs() > 0.000001 {
        numerator / denominator
    } else {
        0.0
    }
}

/// `Histogram<Double>.getMedian` over one increment per value: bins in
/// `Double.compareTo` order (NaN last), counts by value.
fn double_median(values: &[f64]) -> f64 {
    let mut bins: Vec<(f64, f64)> = Vec::new();
    for &value in values {
        let key = if value.is_nan() { f64::NAN } else { value };
        match bins
            .iter_mut()
            .find(|(k, _)| k.total_cmp(&key) == std::cmp::Ordering::Equal)
        {
            Some((_, count)) => *count += 1.0,
            None => bins.push((key, 1.0)),
        }
    }
    bins.sort_by(|a, b| a.0.total_cmp(&b.0));
    let count: f64 = bins.iter().map(|(_, c)| c).sum();
    if count == 0.0 {
        return 0.0;
    }
    if count == 1.0 {
        return bins[0].0;
    }
    let (mid_low, mid_high) = if count % 2.0 == 0.0 {
        (count / 2.0, count / 2.0 + 1.0)
    } else {
        ((count / 2.0).ceil(), (count / 2.0).ceil())
    };
    let mut total = 0.0;
    let mut low = None;
    let mut high = None;
    for (key, value) in &bins {
        total += value;
        if low.is_none() && total >= mid_low {
            low = Some(*key);
        }
        if high.is_none() && total >= mid_high {
            high = Some(*key);
        }
        if low.is_some() && high.is_some() {
            break;
        }
    }
    (low.unwrap_or(0.0) + high.unwrap_or(0.0)) / 2.0
}

struct CoverageMetrics {
    median_cv: f64,
    median_5prime: f64,
    median_3prime: f64,
    median_5prime_to_3prime: f64,
    histogram: Histogram,
}

/// The loop of `computeCoverageMetrics` over the picked transcripts.
fn coverage_metrics(
    picked: &[usize],
    coverage: &[(TranscriptKey, Vec<i32>)],
    model: &GeneModel,
    opts: &Options,
) -> CoverageMetrics {
    let transcript_count = picked.len() as f64;
    let mut cvs = Vec::with_capacity(picked.len());
    let mut five_prime = Vec::with_capacity(picked.len());
    let mut three_prime = Vec::with_capacity(picked.len());
    let mut five_to_three = Vec::with_capacity(picked.len());
    let mut histogram = Histogram::new();
    for &index in picked {
        let (key, counts) = &coverage[index];
        let mut values: Vec<f64> = counts.iter().map(|&c| f64::from(c)).collect();
        if model.genes()[key.0].negative_strand {
            values.reverse();
        }
        let mean = mean_range(&values, 0, values.len() as i64);
        let stdev = stddev(&values, mean);
        cvs.push(stdev / mean);
        let end_bias = opts.end_bias_bases;
        let five = mean_range(&values, 0, end_bias);
        let three = mean_range(&values, values.len() as i64 - end_bias, values.len() as i64);
        five_prime.push(five / mean);
        three_prime.push(three / mean);
        five_to_three.push(divide(five, three));
        let last_index = (values.len() - 1) as f64;
        for percent in 0..=100i64 {
            let p = percent as f64 / 100.0;
            let start = 0f64.max(last_index * (p - 0.005)) as usize;
            let end = last_index.min(last_index * (p + 0.005)) as usize;
            let length = (end - start + 1) as f64;
            let sum: f64 = values[start..=end].iter().fold(0.0, |s, &v| s + v);
            let normalized = (sum / length) / mean;
            histogram.increment_by(percent, normalized / transcript_count);
        }
    }
    CoverageMetrics {
        median_cv: double_median(&cvs),
        median_5prime: double_median(&five_prime),
        median_3prime: double_median(&three_prime),
        median_5prime_to_3prime: double_median(&five_to_three),
        histogram,
    }
}

/// Where a record's template falls (`NUM_R1_`, `NUM_R2_`,
/// `NUM_UNEXPLAINED_READS`).
#[derive(Clone, Copy, Debug)]
enum Template {
    R1,
    R2,
    Unexplained,
}

/// The per-record facts of a mapped primary record that is not ribosomal.
struct Aligned {
    blocks: Vec<AlignmentBlock>,
    /// One function per base of each block.
    functions: Vec<Vec<LocusFunction>>,
    /// Indices of the overlapping genes.
    genes: Vec<usize>,
    /// Correct (true) or incorrect strand, when counted.
    strand: Option<bool>,
    template: Option<Template>,
}

enum Class {
    /// Secondary or unmapped: `PF_BASES` only.
    Dropped,
    Ignored,
    Ribosomal {
        aligned_bases: i64,
    },
    Aligned(Aligned),
}

/// What one record adds to every unit it belongs to.
struct Seen {
    read_length: i64,
    class: Class,
}

/// The collector: one [`Unit`] per accumulation bucket over the shared
/// gene model and ribosomal index.
pub struct RnaSeqMetrics {
    opts: Options,
    model: GeneModel,
    ribosomal: OverlapIndex,
    ignored_contigs: HashSet<usize>,
    levels: Levels<Unit>,
}

impl RnaSeqMetrics {
    /// `CollectRnaSeqMetrics.setup`: the refFlat, the ribosomal intervals
    /// (uniqued, dictionary checked) and the ignored contigs.
    pub fn new(opts: &Options, header: &sam::Header, input: &Path) -> Result<Self> {
        let model = crate::ref_flat::load(&opts.ref_flat, header)?;
        let contigs = header.reference_sequences();
        let ribosomal = match &opts.ribosomal_intervals {
            Some(path) => intervals::read_files(&[path], header).map_err(|e| match e {
                Error::Intervals { path, message } => Error::Intervals {
                    path,
                    message: format!(
                        "Sequence dictionaries differ in {} and the interval list: {message}",
                        input.display()
                    ),
                },
                other => other,
            })?,
            None => Vec::new(),
        };
        let ribosomal = OverlapIndex::new(intervals::uniqued(&ribosomal), contigs.len(), 0);
        let mut ignored_contigs = HashSet::new();
        for name in &opts.ignore_sequence {
            let index = contigs.get_index_of(name.as_bytes()).ok_or_else(|| {
                Error::Usage(format!(
                    "Unrecognized sequence {name} passed as argument to IGNORE_SEQUENCE"
                ))
            })?;
            ignored_contigs.insert(index);
        }
        let ribosomal_initial = opts.ribosomal_intervals.as_ref().map(|_| 0);
        Ok(Self {
            opts: opts.clone(),
            model,
            ribosomal,
            ignored_contigs,
            levels: Levels::new(header, &opts.levels, Unit::new(ribosomal_initial)),
        })
    }

    /// The filters and the per-record facts of `acceptRecord`, computed
    /// once for every unit.
    fn classify(&self, record: &bam::Record) -> Result<Seen> {
        let flags = record.flags();
        let read_length = if flags.is_secondary() {
            0
        } else {
            record.sequence().len() as i64
        };
        let contig = record.reference_sequence_id().transpose()?;
        if !flags.is_unmapped()
            && !flags.is_secondary()
            && contig.is_some_and(|c| self.ignored_contigs.contains(&c))
        {
            return Ok(Seen {
                read_length,
                class: Class::Ignored,
            });
        }
        if flags.is_secondary() || flags.is_unmapped() {
            return Ok(Seen {
                read_length,
                class: Class::Dropped,
            });
        }
        let placement = Placement::of(record)?;
        let blocks = alignment_blocks(record)?;
        if let Some((fragment_start, fragment_end)) = placement.fragment(record)
            && self.ribosomal_fraction(placement.contig, fragment_start, fragment_end)
                >= self.opts.rrna_fragment_percentage
        {
            return Ok(Seen {
                read_length,
                class: Class::Ribosomal {
                    aligned_bases: blocks.iter().map(|b| b.length as i64).sum(),
                },
            });
        }
        let genes = self
            .model
            .overlapping(placement.contig, placement.start, placement.end);
        let (functions, overlaps_exon) = self.locus_functions(&genes, &blocks);
        let (strand, template) = if !flags.is_supplementary() && overlaps_exon && genes.len() == 1 {
            self.strand_and_template(record, &placement, &self.model.genes()[genes[0]])?
        } else {
            (None, None)
        };
        Ok(Seen {
            read_length,
            class: Class::Aligned(Aligned {
                blocks,
                functions,
                genes,
                strand,
                template,
            }),
        })
    }

    /// The longest single ribosomal interval overlap over the fragment
    /// length (NaN for a 0-length fragment, which no percentage reaches).
    fn ribosomal_fraction(&self, contig: usize, fragment_start: i64, fragment_end: i64) -> f64 {
        let mut intersection = 0i64;
        for i in self
            .ribosomal
            .overlaps(contig, fragment_start, fragment_end)
        {
            let interval = self.ribosomal.intervals()[i];
            let overlap = interval.end.min(fragment_end) - interval.start.max(fragment_start) + 1;
            intersection = intersection.max(overlap);
        }
        intersection as f64 / (fragment_end - fragment_start + 1) as f64
    }

    /// The function of each base of each block over every transcript of
    /// the overlapping genes, and whether any base is coding or UTR.
    fn locus_functions(
        &self,
        genes: &[usize],
        blocks: &[AlignmentBlock],
    ) -> (Vec<Vec<LocusFunction>>, bool) {
        let mut overlaps_exon = false;
        let functions = blocks
            .iter()
            .map(|block| {
                let mut functions = vec![LocusFunction::Intergenic; block.length];
                for &gene in genes {
                    for transcript in &self.model.genes()[gene].transcripts {
                        transcript.assign_locus_function(block.reference_start, &mut functions);
                    }
                }
                overlaps_exon |= functions
                    .iter()
                    .any(|f| matches!(f, LocusFunction::Utr | LocusFunction::Coding));
                functions
            })
            .collect();
        (functions, overlaps_exon)
    }

    /// The strand-specificity verdict and the template class of a record
    /// over exactly one gene.
    fn strand_and_template(
        &self,
        record: &bam::Record,
        placement: &Placement,
        gene: &crate::ref_flat::Gene,
    ) -> Result<(Option<bool>, Option<Template>)> {
        let flags = record.flags();
        let negative_transcription = gene.negative_strand;
        let read_one_or_unpaired = !placement.paired || flags.is_first_segment();
        let negative_read = flags.is_reverse_complemented();
        let strand = (self.opts.strand != Strand::None).then(|| {
            let agree = negative_read == negative_transcription;
            let first_expected = self.opts.strand == Strand::FirstReadTranscriptionStrand;
            let this_expected = read_one_or_unpaired == first_expected;
            agree == this_expected
        });
        if !read_one_or_unpaired {
            return Ok((strand, None));
        }
        let (proper, left, right) = if !placement.paired {
            (true, placement.start, placement.end)
        } else if flags.is_mate_unmapped() {
            (false, 0, 0)
        } else {
            let mate_start = placement.mate_start.ok_or_else(|| Error::MalformedRecord {
                name: crate::record::record_name(record),
                what: "missing mate alignment start",
            })?;
            let mate_length =
                mate_reference_length(record)?.unwrap_or(record.sequence().len() as i64);
            let proper = placement.mate_contig == Some(placement.contig)
                && orientation(record)? == Orientation::Fr;
            (
                proper,
                placement.start.min(mate_start),
                placement.end.max(mate_start + mate_length - 1),
            )
        };
        let template = if proper && left >= gene.start && right <= gene.end {
            if negative_read == negative_transcription {
                Template::R1
            } else {
                Template::R2
            }
        } else {
            Template::Unexplained
        };
        Ok((strand, Some(template)))
    }
}

/// Where a mapped record and its mate lie.
struct Placement {
    contig: usize,
    start: i64,
    end: i64,
    paired: bool,
    mate_contig: Option<usize>,
    mate_start: Option<i64>,
}

impl Placement {
    fn of(record: &bam::Record) -> Result<Self> {
        let contig =
            record
                .reference_sequence_id()
                .transpose()?
                .ok_or_else(|| Error::MalformedRecord {
                    name: crate::record::record_name(record),
                    what: "mapped record without a reference sequence",
                })?;
        let start = position(record, record.alignment_start(), "missing alignment start")?;
        Ok(Self {
            contig,
            start,
            end: alignment_end(record)?.unwrap_or(start - 1),
            paired: record.flags().is_segmented(),
            mate_contig: record.mate_reference_sequence_id().transpose()?,
            mate_start: match record.mate_alignment_start() {
                Some(p) => Some(p?.get() as i64),
                None => None,
            },
        })
    }

    /// The fragment the rRNA test measures: the alignment of an unpaired
    /// record, `|TLEN|` bases from the lower start of a pair on one contig,
    /// none for a pair with an unmapped or cross-contig mate.
    fn fragment(&self, record: &bam::Record) -> Option<(i64, i64)> {
        if !self.paired {
            return Some((self.start, self.end));
        }
        if record.flags().is_mate_unmapped() || self.mate_contig != Some(self.contig) {
            return None;
        }
        let fragment_start = self.start.min(self.mate_start.unwrap_or(self.start));
        let length = i64::from(record.template_length()).abs();
        Some((fragment_start, fragment_start + length - 1))
    }
}
/// The rows and histogram columns of a run, in file order.
pub struct Report {
    pub rows: Vec<Row>,
    pub series: Vec<HistogramSeries>,
}

impl Report {
    /// The metrics file: every unit's row, then the `normalized_position`
    /// histogram with one `<key>.normalized_coverage` column per unit
    /// (an empty one is left out of the file).
    pub fn metrics_output(&self) -> MetricsOutput {
        MetricsOutput::from_rows(
            &self.rows,
            Some(HistogramSection {
                key_label: "normalized_position",
                key_class: KeyClass::Integer,
                series: self.series.clone(),
            }),
        )
    }
}

impl From<Report> for ToolOutput {
    fn from(report: Report) -> Self {
        report.metrics_output().into()
    }
}

impl Accumulator for RnaSeqMetrics {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if record.flags().is_qc_fail() {
            return Ok(());
        }
        let seen = self.classify(record)?;
        let model = &self.model;
        let mut failure = None;
        self.levels.route(record, |unit| {
            if failure.is_none()
                && let Err(e) = unit.tally(&seen, model)
            {
                failure = Some(e);
            }
        });
        failure.map_or(Ok(()), Err)
    }

    fn finish(self) -> Result<Report> {
        let mut rows = Vec::new();
        let mut series = Vec::new();
        for bucket in self.levels.into_buckets() {
            let (row, histogram) = bucket.state.finish(&self.model, &self.opts, bucket.labels);
            rows.push(row);
            series.push(HistogramSeries {
                label: format!("{}.normalized_coverage", bucket.key),
                histogram,
            });
        }
        Ok(Report { rows, series })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn double_median_follows_histogram_get_median() {
        assert_eq!(double_median(&[]), 0.0);
        assert_eq!(double_median(&[0.5]), 0.5);
        assert_eq!(double_median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(double_median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
        // Equal values share a bin: the bin of 2 holds the 2nd and 3rd.
        assert_eq!(double_median(&[2.0, 2.0, 1.0, 5.0]), 2.0);
        assert!(double_median(&[1.0, f64::NAN]).is_nan());
        assert_eq!(double_median(&[1.0, 1.0, f64::NAN]), 1.0);
    }

    #[test]
    fn math_util_pieces() {
        assert_eq!(mean(&[1, 2, 3, 6]), 3.0);
        let values = [1.0, 2.0, 3.0, 6.0];
        assert_eq!(mean_range(&values, 0, 2), 1.5);
        assert_eq!(mean_range(&values, 2, 4), 4.5);
        assert!((stddev(&values, 3.0) - (12.5f64 - 9.0).sqrt()).abs() < 1e-12);
        assert_eq!(divide(1.0, 0.0), 0.0);
        assert_eq!(divide(1.0, 0.0000005), 0.0);
        assert_eq!(divide(3.0, 2.0), 1.5);
    }

    #[test]
    fn strand_values_are_picards() {
        assert_eq!(Strand::parse("NONE").unwrap(), Strand::None);
        assert!(
            Strand::parse("BOGUS")
                .unwrap_err()
                .to_string()
                .contains("'BOGUS' is not a valid value for StrandSpecificity.")
        );
    }
}

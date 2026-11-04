//! `CollectHsMetrics` and `CollectTargetedPcrMetrics`: Picard 3.5.0
//! `CollectTargetedMetrics.doWork` over `TargetMetricsCollector` (its
//! `PerUnitTargetMetricCollector`, `acceptRecord` and `finish`) with
//! `HsMetricCollector.convertMetric` or
//! `TargetedPcrMetricsCollector.convertMetric`, at every accumulation level.
//!
//! One pass over the BAM in any order. Each record is placed against the
//! uniqued baits and targets with binary searches (`intervals::OverlapIndex`),
//! and coverage is kept only for target bases: one depth array per uniqued
//! target, laid end to end, at high quality and unfiltered, per unit.
//!
//! The four histograms (`highQualityDepthHistogram`,
//! `unfilteredDepthHistogram`, `unfilteredBaseQHistogram`,
//! `unCappedUnfilteredBaseQHistogram`) are fields of the collector, not of
//! the unit, so every unit's `finish` adds to the same four. A unit's
//! `MEDIAN_TARGET_COVERAGE`, `FOLD_80_BASE_PENALTY` and `HET_SNP_*` read
//! what the units before it left there, the file prints the final
//! histograms once per unit, and the check that the capped base-quality
//! counts equal the unfiltered depth total fails for every unit that
//! finishes after a unit with a counted base
//! (`numbers of bases in the base quality histogram and the coverage
//! histogram are not equal`). `Shared` keeps that behaviour.

use std::io::{self, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::{bam, sam};

use crate::adapter::Adapters;
use crate::duplication::{estimate_library_size, estimate_roi};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::intervals::{self, Interval, OverlapIndex};
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::{KeyClass, format_double};
use crate::pass::{self, Accumulator};
use crate::record::{position, record_name};
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool, write_files};
use crate::sensitivity::{
    calculate_sensitivities, het_snp_sensitivity, math_round, normalize_histogram,
    phred_from_error_probability,
};

/// `CollectHsMetrics` overrides of `CollectTargetedMetrics` defaults, and
/// the defaults it keeps.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 20;
pub const DEFAULT_MINIMUM_BASE_QUALITY: i32 = 20;
pub const DEFAULT_COVERAGE_CAP: i32 = 200;
/// `TargetedPcrMetricsCollector.NEAR_PROBE_DISTANCE_DEFAULT`.
pub const DEFAULT_NEAR_DISTANCE: i32 = 250;
pub const DEFAULT_SAMPLE_SIZE: usize = 10_000;
/// `CollectTargetedMetrics.ALLELE_FRACTION` default.
pub const DEFAULT_ALLELE_FRACTIONS: [f64; 9] = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5];

/// `TargetMetricsCollector.LOG_ODDS_THRESHOLD`.
const LOG_ODDS_THRESHOLD: f64 = 3.0;
/// Bases at or below this quality add no coverage at all, even unfiltered.
const MIN_UNFILTERED_QUALITY: i32 = 2;
/// `baseQHistogramArray` has `Byte.MAX_VALUE` slots, qualities 0 to 126.
const QUALITY_SLOTS: usize = 127;
/// `calculateTargetCoverageMetrics.targetBasesDepth`.
const TARGET_BASE_DEPTHS: [i64; 18] = [
    0, 1, 2, 10, 20, 30, 40, 50, 100, 250, 500, 1000, 2500, 5000, 10_000, 25_000, 50_000, 100_000,
];
/// `HsMetricCollector.convertMetric`: coverage goals of the `HS_PENALTY_*X` fields.
const PENALTY_GOALS: [i32; 6] = [10, 20, 30, 40, 50, 100];
/// `calculateGcMetrics`: one bin per whole GC percent.
const GC_BINS: usize = 101;

/// Which `CollectTargetedMetrics` program: they differ in defaults, option
/// names, the metrics class and its probe columns.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Variant {
    Hs,
    TargetedPcr,
}

impl Variant {
    pub fn name(self) -> &'static str {
        match self {
            Variant::Hs => "CollectHsMetrics",
            Variant::TargetedPcr => "CollectTargetedPcrMetrics",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub variant: Variant,
    /// `BAIT_INTERVALS` (`AMPLICON_INTERVALS`, one file): interval lists (or
    /// BED files) of the probes.
    pub bait_intervals: Vec<PathBuf>,
    /// `TARGET_INTERVALS`.
    pub target_intervals: Vec<PathBuf>,
    /// `BAIT_SET_NAME` (`CUSTOM_AMPLICON_SET_NAME`); by default the probe
    /// file names up to their first dot.
    pub bait_set_name: Option<String>,
    /// `REFERENCE_SEQUENCE`: needed for `AT_DROPOUT`, `GC_DROPOUT` and
    /// `PER_TARGET_COVERAGE`.
    pub reference: Option<PathBuf>,
    /// `METRIC_ACCUMULATION_LEVEL`.
    pub levels: Vec<Level>,
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub clip_overlapping_reads: bool,
    /// `INCLUDE_INDELS`: inserted bases count as on target, deleted bases as covered.
    pub include_indels: bool,
    pub coverage_cap: i32,
    pub near_distance: i32,
    /// `SAMPLE_SIZE` of the sensitivity simulations.
    pub sample_size: usize,
    /// `PER_BASE_COVERAGE`: the high-quality depth of every target base of
    /// the `ALL_READS` unit. A name ending in `.gz`, `.gzip` or `.bfq` is
    /// compressed.
    pub per_base_coverage: Option<PathBuf>,
    /// `PER_TARGET_COVERAGE`: one line per uniqued target of the `ALL_READS`
    /// unit, written with the GC metrics.
    pub per_target_coverage: Option<PathBuf>,
    /// `THEORETICAL_SENSITIVITY_OUTPUT`: a second metrics file, one row per
    /// entry of `allele_fractions`.
    pub theoretical_sensitivity_output: Option<PathBuf>,
    /// `ALLELE_FRACTION`.
    pub allele_fractions: Vec<f64>,
}

impl Default for Options {
    /// `CollectHsMetrics` defaults.
    fn default() -> Self {
        Self {
            variant: Variant::Hs,
            bait_intervals: Vec::new(),
            target_intervals: Vec::new(),
            bait_set_name: None,
            reference: None,
            levels: vec![Level::AllReads],
            minimum_mapping_quality: DEFAULT_MINIMUM_MAPPING_QUALITY,
            minimum_base_quality: DEFAULT_MINIMUM_BASE_QUALITY,
            clip_overlapping_reads: true,
            include_indels: false,
            coverage_cap: DEFAULT_COVERAGE_CAP,
            near_distance: DEFAULT_NEAR_DISTANCE,
            sample_size: DEFAULT_SAMPLE_SIZE,
            per_base_coverage: None,
            per_target_coverage: None,
            theoretical_sensitivity_output: None,
            allele_fractions: DEFAULT_ALLELE_FRACTIONS.to_vec(),
        }
    }
}

impl Options {
    /// `CollectTargetedPcrMetrics` defaults: the `CollectTargetedMetrics`
    /// values `CollectHsMetrics` overrides.
    pub fn targeted_pcr() -> Self {
        Self {
            variant: Variant::TargetedPcr,
            minimum_mapping_quality: 1,
            minimum_base_quality: 0,
            clip_overlapping_reads: false,
            ..Self::default()
        }
    }

    /// Picard's argument parser messages for a missing probe or target file:
    /// `BAIT_INTERVALS` and `TARGET_INTERVALS` are lists, `AMPLICON_INTERVALS`
    /// is one required value.
    pub fn check_interval_files(&self) -> Result<()> {
        let missing = |name: &str| {
            Err(Error::Usage(format!(
                "Option '{name}' must be specified at least 1 times."
            )))
        };
        if self.bait_intervals.is_empty() {
            return match self.variant {
                Variant::Hs => missing("BAIT_INTERVALS"),
                Variant::TargetedPcr => Err(Error::Usage(
                    "Option 'AMPLICON_INTERVALS' is required.".to_owned(),
                )),
            };
        }
        if self.target_intervals.is_empty() {
            return missing("TARGET_INTERVALS");
        }
        Ok(())
    }

    fn from_invocation(inv: &Invocation, variant: Variant) -> Result<Self> {
        let defaults = match variant {
            Variant::Hs => Self::default(),
            Variant::TargetedPcr => Self::targeted_pcr(),
        };
        let (bait_intervals, bait_set_name) = match variant {
            Variant::Hs => (
                inv.all("BAIT_INTERVALS")
                    .iter()
                    .map(PathBuf::from)
                    .collect(),
                inv.get("BAIT_SET_NAME").map(str::to_owned),
            ),
            Variant::TargetedPcr => {
                let amplicons = inv.all("AMPLICON_INTERVALS");
                if amplicons.len() > 1 {
                    return Err(Error::Usage(
                        "Option 'AMPLICON_INTERVALS' cannot be specified more than once."
                            .to_owned(),
                    ));
                }
                (
                    amplicons.iter().map(PathBuf::from).collect(),
                    inv.get("CUSTOM_AMPLICON_SET_NAME").map(str::to_owned),
                )
            }
        };
        let opts = Self {
            variant,
            bait_intervals,
            target_intervals: inv
                .all("TARGET_INTERVALS")
                .iter()
                .map(PathBuf::from)
                .collect(),
            bait_set_name,
            reference: inv.get("REFERENCE_SEQUENCE").map(PathBuf::from),
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
            minimum_mapping_quality: inv
                .parse("MINIMUM_MAPPING_QUALITY", defaults.minimum_mapping_quality)?,
            minimum_base_quality: inv
                .parse("MINIMUM_BASE_QUALITY", defaults.minimum_base_quality)?,
            clip_overlapping_reads: inv
                .flag("CLIP_OVERLAPPING_READS", defaults.clip_overlapping_reads)?,
            include_indels: inv.flag("INCLUDE_INDELS", false)?,
            coverage_cap: inv.parse("COVERAGE_CAP", DEFAULT_COVERAGE_CAP)?,
            near_distance: inv.parse("NEAR_DISTANCE", DEFAULT_NEAR_DISTANCE)?,
            sample_size: inv.parse("SAMPLE_SIZE", DEFAULT_SAMPLE_SIZE)?,
            per_base_coverage: inv.get("PER_BASE_COVERAGE").map(PathBuf::from),
            per_target_coverage: inv.get("PER_TARGET_COVERAGE").map(PathBuf::from),
            theoretical_sensitivity_output: inv
                .get("THEORETICAL_SENSITIVITY_OUTPUT")
                .map(PathBuf::from),
            allele_fractions: inv.list(
                "ALLELE_FRACTION",
                DEFAULT_ALLELE_FRACTIONS.to_vec(),
                |v| {
                    v.parse()
                        .map_err(|_| Error::Usage(format!("ALLELE_FRACTION={v}: not a number")))
                },
            )?,
        };
        opts.check_interval_files()?;
        // CollectTargetedMetrics.customCommandLineValidation.
        if opts.per_target_coverage.is_some() && opts.levels != [Level::AllReads] {
            return Err(Error::Usage(
                "PER_TARGET_COVERAGE can be specified only when METRIC_ACCUMULATION_LEVEL is set to ALL_READS."
                    .to_owned(),
            ));
        }
        if opts.per_target_coverage.is_some() && opts.reference.is_none() {
            return Err(Error::Usage(
                "Must supply REFERENCE_SEQUENCE when supplying PER_TARGET_COVERAGE".to_owned(),
            ));
        }
        Ok(opts)
    }
}

const COMMON_OPTIONS: [Opt; 15] = [
    Opt::value(&["TARGET_INTERVALS", "TI"]),
    Opt::value(&["REFERENCE_SEQUENCE", "R"]),
    Opt::value(&["METRIC_ACCUMULATION_LEVEL", "LEVEL"]),
    Opt::value(&["PER_TARGET_COVERAGE"]),
    Opt::value(&["PER_BASE_COVERAGE"]),
    Opt::value(&["NEAR_DISTANCE"]),
    Opt::value(&["MINIMUM_MAPPING_QUALITY", "MQ"]),
    Opt::value(&["MINIMUM_BASE_QUALITY", "Q"]),
    Opt::value(&["CLIP_OVERLAPPING_READS"]),
    Opt::value(&["INCLUDE_INDELS"]),
    Opt::value(&["COVERAGE_CAP", "covMax"]),
    Opt::value(&["SAMPLE_SIZE"]),
    Opt::value(&["THEORETICAL_SENSITIVITY_OUTPUT"]),
    Opt::value(&["ALLELE_FRACTION"]),
    // Not a SinglePassSamProgram: Picard's parser rejects both.
    Opt::rejected(&["ASSUME_SORTED", "AS"]),
];

macro_rules! options {
    ($($own:expr),*) => {{
        const OPTIONS: [Opt; COMMON_OPTIONS.len() + 3] = {
            let mut all = [Opt::rejected(&["STOP_AFTER"]); COMMON_OPTIONS.len() + 3];
            let own = [$($own),*];
            let mut i = 0;
            while i < COMMON_OPTIONS.len() {
                all[i] = COMMON_OPTIONS[i];
                i += 1;
            }
            let mut j = 0;
            while j < own.len() {
                all[COMMON_OPTIONS.len() + j] = own[j];
                j += 1;
            }
            all
        };
        &OPTIONS
    }};
}

pub const TOOL: Tool = Tool {
    name: "CollectHsMetrics",
    options: options!(
        Opt::value(&["BAIT_INTERVALS", "BI"]),
        Opt::value(&["BAIT_SET_NAME", "N"])
    ),
    run: run_hs,
    collector: Some(hs_collector),
};

pub const PCR_TOOL: Tool = Tool {
    name: "CollectTargetedPcrMetrics",
    options: options!(
        Opt::value(&["AMPLICON_INTERVALS", "AI"]),
        Opt::value(&["CUSTOM_AMPLICON_SET_NAME", "N"])
    ),
    run: run_pcr,
    collector: None,
};

fn run_hs(inv: &Invocation) -> Result<ToolOutput> {
    run_variant(inv, Variant::Hs)
}

fn run_pcr(inv: &Invocation) -> Result<ToolOutput> {
    run_variant(inv, Variant::TargetedPcr)
}

/// `CollectHsMetrics` as a one-pass collector, for `qc`: the metrics file
/// and, when asked for, the sensitivity file. An `ALLELE_FRACTION` outside
/// 0 to 1 fails the collector and writes no file. Picard has no
/// `CollectMultipleMetrics` program for it, and `run_hs` keeps its own
/// read.
fn hs_collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv, Variant::Hs)?;
    let capture = CaptureMetrics::new(&opts, header)?;
    Ok(pass::boxed(capture, move |report| {
        let mut output = ToolOutput::file(report.metrics_output());
        if let Some(file) = sensitivity_file(&opts, &report)? {
            output.extra.push(file);
        }
        Ok(output)
    }))
}

/// The `THEORETICAL_SENSITIVITY_OUTPUT` file, when the option is given.
fn sensitivity_file(opts: &Options, report: &Report) -> Result<Option<(PathBuf, MetricsOutput)>> {
    let Some(path) = &opts.theoretical_sensitivity_output else {
        return Ok(None);
    };
    let rows = calculate_sensitivities(
        opts.sample_size,
        &report.unfiltered_depth,
        &report.unfiltered_base_quality,
        &opts.allele_fractions,
    )?;
    Ok(Some((path.clone(), MetricsOutput::from_rows(&rows, None))))
}

/// `CollectTargetedMetrics.doWork`: the pass, the metrics file, then the
/// sensitivity file from the collector's shared histograms. Picard writes
/// `OUTPUT` before it computes the sensitivities, so an `ALLELE_FRACTION`
/// outside 0 to 1 leaves the metrics file behind and exits 1: the file is
/// written here before the error is returned.
fn run_variant(inv: &Invocation, variant: Variant) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv, variant)?;
    let report = pass::run(&inv.input, inv.threads, |header| {
        CaptureMetrics::new(&opts, header)
    })?;
    let mut output = ToolOutput::file(report.metrics_output());
    match sensitivity_file(&opts, &report) {
        Ok(file) => output.extra.extend(file),
        Err(error) => {
            if inv.output.is_some() {
                write_files(inv, &output)?;
            }
            return Err(error);
        }
    }
    Ok(output)
}

/// `HsMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 69] = [
    "BAIT_SET",
    "BAIT_TERRITORY",
    "BAIT_DESIGN_EFFICIENCY",
    "ON_BAIT_BASES",
    "NEAR_BAIT_BASES",
    "OFF_BAIT_BASES",
    "PCT_SELECTED_BASES",
    "PCT_OFF_BAIT",
    "ON_BAIT_VS_SELECTED",
    "MEAN_BAIT_COVERAGE",
    "PCT_USABLE_BASES_ON_BAIT",
    "PCT_USABLE_BASES_ON_TARGET",
    "FOLD_ENRICHMENT",
    "HS_LIBRARY_SIZE",
    "HS_PENALTY_10X",
    "HS_PENALTY_20X",
    "HS_PENALTY_30X",
    "HS_PENALTY_40X",
    "HS_PENALTY_50X",
    "HS_PENALTY_100X",
    "TARGET_TERRITORY",
    "GENOME_SIZE",
    "TOTAL_READS",
    "PF_READS",
    "PF_BASES",
    "PF_UNIQUE_READS",
    "PF_UQ_READS_ALIGNED",
    "PF_BASES_ALIGNED",
    "PF_UQ_BASES_ALIGNED",
    "ON_TARGET_BASES",
    "PCT_PF_READS",
    "PCT_PF_UQ_READS",
    "PCT_PF_UQ_READS_ALIGNED",
    "MEAN_TARGET_COVERAGE",
    "MEDIAN_TARGET_COVERAGE",
    "MAX_TARGET_COVERAGE",
    "MIN_TARGET_COVERAGE",
    "ZERO_CVG_TARGETS_PCT",
    "PCT_EXC_DUPE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_OFF_TARGET",
    "FOLD_80_BASE_PENALTY",
    "PCT_TARGET_BASES_1X",
    "PCT_TARGET_BASES_2X",
    "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "PCT_TARGET_BASES_40X",
    "PCT_TARGET_BASES_50X",
    "PCT_TARGET_BASES_100X",
    "PCT_TARGET_BASES_250X",
    "PCT_TARGET_BASES_500X",
    "PCT_TARGET_BASES_1000X",
    "PCT_TARGET_BASES_2500X",
    "PCT_TARGET_BASES_5000X",
    "PCT_TARGET_BASES_10000X",
    "PCT_TARGET_BASES_25000X",
    "PCT_TARGET_BASES_50000X",
    "PCT_TARGET_BASES_100000X",
    "AT_DROPOUT",
    "GC_DROPOUT",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// `TargetedPcrMetrics` columns in Picard 3.5.0 file order: the amplicon
/// columns, the three `TargetMetricsBase` pair counts, then the
/// `PanelMetricsBase` columns `HsMetrics` shares.
pub const PCR_COLUMNS: [&str; 62] = [
    "CUSTOM_AMPLICON_SET",
    "AMPLICON_TERRITORY",
    "ON_AMPLICON_BASES",
    "NEAR_AMPLICON_BASES",
    "OFF_AMPLICON_BASES",
    "PCT_AMPLIFIED_BASES",
    "PCT_OFF_AMPLICON",
    "ON_AMPLICON_VS_SELECTED",
    "MEAN_AMPLICON_COVERAGE",
    "FOLD_ENRICHMENT",
    "PF_SELECTED_PAIRS",
    "PF_SELECTED_UNIQUE_PAIRS",
    "ON_TARGET_FROM_PAIR_BASES",
    "TARGET_TERRITORY",
    "GENOME_SIZE",
    "TOTAL_READS",
    "PF_READS",
    "PF_BASES",
    "PF_UNIQUE_READS",
    "PF_UQ_READS_ALIGNED",
    "PF_BASES_ALIGNED",
    "PF_UQ_BASES_ALIGNED",
    "ON_TARGET_BASES",
    "PCT_PF_READS",
    "PCT_PF_UQ_READS",
    "PCT_PF_UQ_READS_ALIGNED",
    "MEAN_TARGET_COVERAGE",
    "MEDIAN_TARGET_COVERAGE",
    "MAX_TARGET_COVERAGE",
    "MIN_TARGET_COVERAGE",
    "ZERO_CVG_TARGETS_PCT",
    "PCT_EXC_DUPE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_OFF_TARGET",
    "FOLD_80_BASE_PENALTY",
    "PCT_TARGET_BASES_1X",
    "PCT_TARGET_BASES_2X",
    "PCT_TARGET_BASES_10X",
    "PCT_TARGET_BASES_20X",
    "PCT_TARGET_BASES_30X",
    "PCT_TARGET_BASES_40X",
    "PCT_TARGET_BASES_50X",
    "PCT_TARGET_BASES_100X",
    "PCT_TARGET_BASES_250X",
    "PCT_TARGET_BASES_500X",
    "PCT_TARGET_BASES_1000X",
    "PCT_TARGET_BASES_2500X",
    "PCT_TARGET_BASES_5000X",
    "PCT_TARGET_BASES_10000X",
    "PCT_TARGET_BASES_25000X",
    "PCT_TARGET_BASES_50000X",
    "PCT_TARGET_BASES_100000X",
    "AT_DROPOUT",
    "GC_DROPOUT",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// `PanelMetricsBase`: the fields both metric classes share, with the
/// `MultilevelMetrics` labels.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Panel {
    pub target_territory: i64,
    pub genome_size: i64,
    pub total_reads: i64,
    pub pf_reads: i64,
    pub pf_bases: i64,
    pub pf_unique_reads: i64,
    pub pf_uq_reads_aligned: i64,
    pub pf_bases_aligned: i64,
    pub pf_uq_bases_aligned: i64,
    pub on_target_bases: i64,
    pub pct_pf_reads: f64,
    pub pct_pf_uq_reads: f64,
    pub pct_pf_uq_reads_aligned: f64,
    pub mean_target_coverage: f64,
    pub median_target_coverage: f64,
    pub max_target_coverage: i64,
    pub min_target_coverage: i64,
    pub zero_cvg_targets_pct: f64,
    pub pct_exc_dupe: f64,
    pub pct_exc_adapter: f64,
    pub pct_exc_mapq: f64,
    pub pct_exc_baseq: f64,
    pub pct_exc_overlap: f64,
    pub pct_exc_off_target: f64,
    pub fold_80_base_penalty: f64,
    /// `PCT_TARGET_BASES_1X` to `PCT_TARGET_BASES_100000X`, the depths of
    /// `TARGET_BASE_DEPTHS` after 0.
    pub pct_target_bases: [f64; 17],
    pub at_dropout: f64,
    pub gc_dropout: f64,
    pub het_snp_sensitivity: f64,
    /// A Java `double` holding an `int` Phred score.
    pub het_snp_q: f64,
    pub labels: Labels,
}

impl Panel {
    /// The `TARGET_TERRITORY` to `READ_GROUP` cells.
    fn values(&self) -> Vec<Value> {
        let mut values = vec![
            self.target_territory.into(),
            self.genome_size.into(),
            self.total_reads.into(),
            self.pf_reads.into(),
            self.pf_bases.into(),
            self.pf_unique_reads.into(),
            self.pf_uq_reads_aligned.into(),
            self.pf_bases_aligned.into(),
            self.pf_uq_bases_aligned.into(),
            self.on_target_bases.into(),
            self.pct_pf_reads.into(),
            self.pct_pf_uq_reads.into(),
            self.pct_pf_uq_reads_aligned.into(),
            self.mean_target_coverage.into(),
            self.median_target_coverage.into(),
            self.max_target_coverage.into(),
            self.min_target_coverage.into(),
            self.zero_cvg_targets_pct.into(),
            self.pct_exc_dupe.into(),
            self.pct_exc_adapter.into(),
            self.pct_exc_mapq.into(),
            self.pct_exc_baseq.into(),
            self.pct_exc_overlap.into(),
            self.pct_exc_off_target.into(),
            self.fold_80_base_penalty.into(),
        ];
        values.extend(self.pct_target_bases.iter().map(|&p| Value::Float(p)));
        values.extend([
            self.at_dropout.into(),
            self.gc_dropout.into(),
            self.het_snp_sensitivity.into(),
            self.het_snp_q.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]);
        values
    }
}

/// Picard's `TargetMetrics` after `PerUnitTargetMetricCollector.finish`:
/// the probe columns under their generic names plus the panel.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TargetMetrics {
    pub probe_set: String,
    pub probe_territory: i64,
    pub on_probe_bases: i64,
    pub near_probe_bases: i64,
    pub off_probe_bases: i64,
    pub pct_selected_bases: f64,
    pub pct_off_probe: f64,
    pub on_probe_vs_selected: f64,
    pub mean_probe_coverage: f64,
    pub fold_enrichment: f64,
    pub pf_selected_pairs: i64,
    pub pf_selected_unique_pairs: i64,
    pub on_target_from_pair_bases: i64,
    pub panel: Panel,
}

/// The `HsMetrics` row (`HsMetricCollector.convertMetric`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Row {
    pub bait_set: String,
    pub bait_territory: i64,
    pub bait_design_efficiency: f64,
    pub on_bait_bases: i64,
    pub near_bait_bases: i64,
    pub off_bait_bases: i64,
    pub pct_selected_bases: f64,
    pub pct_off_bait: f64,
    pub on_bait_vs_selected: f64,
    pub mean_bait_coverage: f64,
    pub pct_usable_bases_on_bait: f64,
    pub pct_usable_bases_on_target: f64,
    pub fold_enrichment: f64,
    /// `estimateLibrarySize` returns null when it cannot estimate.
    pub hs_library_size: Option<i64>,
    /// `HS_PENALTY_10X` to `HS_PENALTY_100X`, in `PENALTY_GOALS` order.
    pub hs_penalty: [f64; 6],
    pub panel: Panel,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.directed.HsMetrics";
    const COLUMNS: &'static [&'static str] = &COLUMNS;

    fn values(&self) -> Vec<Value> {
        let mut values = vec![
            self.bait_set.clone().into(),
            self.bait_territory.into(),
            self.bait_design_efficiency.into(),
            self.on_bait_bases.into(),
            self.near_bait_bases.into(),
            self.off_bait_bases.into(),
            self.pct_selected_bases.into(),
            self.pct_off_bait.into(),
            self.on_bait_vs_selected.into(),
            self.mean_bait_coverage.into(),
            self.pct_usable_bases_on_bait.into(),
            self.pct_usable_bases_on_target.into(),
            self.fold_enrichment.into(),
            self.hs_library_size.into(),
        ];
        values.extend(self.hs_penalty.iter().map(|&p| Value::Float(p)));
        values.extend(self.panel.values());
        values
    }
}

/// The `TargetedPcrMetrics` row (`TargetedPcrMetricsCollector.convertMetric`).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PcrRow {
    pub custom_amplicon_set: String,
    pub amplicon_territory: i64,
    pub on_amplicon_bases: i64,
    pub near_amplicon_bases: i64,
    pub off_amplicon_bases: i64,
    pub pct_amplified_bases: f64,
    pub pct_off_amplicon: f64,
    pub on_amplicon_vs_selected: f64,
    pub mean_amplicon_coverage: f64,
    pub fold_enrichment: f64,
    pub pf_selected_pairs: i64,
    pub pf_selected_unique_pairs: i64,
    pub on_target_from_pair_bases: i64,
    pub panel: Panel,
}

impl MetricRow for PcrRow {
    const CLASS: &'static str = "picard.analysis.directed.TargetedPcrMetrics";
    const COLUMNS: &'static [&'static str] = &PCR_COLUMNS;

    fn values(&self) -> Vec<Value> {
        let mut values = vec![
            self.custom_amplicon_set.clone().into(),
            self.amplicon_territory.into(),
            self.on_amplicon_bases.into(),
            self.near_amplicon_bases.into(),
            self.off_amplicon_bases.into(),
            self.pct_amplified_bases.into(),
            self.pct_off_amplicon.into(),
            self.on_amplicon_vs_selected.into(),
            self.mean_amplicon_coverage.into(),
            self.fold_enrichment.into(),
            self.pf_selected_pairs.into(),
            self.pf_selected_unique_pairs.into(),
            self.on_target_from_pair_bases.into(),
        ];
        values.extend(self.panel.values());
        values
    }
}

impl TargetMetrics {
    /// `HsMetricCollector.convertMetric`.
    fn into_hs(self) -> Row {
        let mut row = Row {
            bait_set: self.probe_set.clone(),
            bait_territory: self.probe_territory,
            on_bait_bases: self.on_probe_bases,
            near_bait_bases: self.near_probe_bases,
            off_bait_bases: self.off_probe_bases,
            pct_selected_bases: self.pct_selected_bases,
            pct_off_bait: self.pct_off_probe,
            on_bait_vs_selected: self.on_probe_vs_selected,
            mean_bait_coverage: self.mean_probe_coverage,
            fold_enrichment: self.fold_enrichment,
            ..Row::default()
        };
        row.bait_design_efficiency = self.panel.target_territory as f64 / row.bait_territory as f64;
        row.pct_usable_bases_on_bait = row.on_bait_bases as f64 / self.panel.pf_bases as f64;
        row.pct_usable_bases_on_target =
            self.panel.on_target_bases as f64 / self.panel.pf_bases as f64;
        row.hs_library_size =
            estimate_library_size(self.pf_selected_pairs, self.pf_selected_unique_pairs);
        row.hs_penalty = PENALTY_GOALS.map(|goal| hs_penalty(row.hs_library_size, &self, goal));
        row.panel = self.panel;
        row
    }

    /// `TargetedPcrMetricsCollector.convertMetric`.
    fn into_pcr(self) -> PcrRow {
        PcrRow {
            custom_amplicon_set: self.probe_set,
            amplicon_territory: self.probe_territory,
            on_amplicon_bases: self.on_probe_bases,
            near_amplicon_bases: self.near_probe_bases,
            off_amplicon_bases: self.off_probe_bases,
            pct_amplified_bases: self.pct_selected_bases,
            pct_off_amplicon: self.pct_off_probe,
            on_amplicon_vs_selected: self.on_probe_vs_selected,
            mean_amplicon_coverage: self.mean_probe_coverage,
            fold_enrichment: self.fold_enrichment,
            pf_selected_pairs: self.pf_selected_pairs,
            pf_selected_unique_pairs: self.pf_selected_unique_pairs,
            on_target_from_pair_bases: self.on_target_from_pair_bases,
            panel: self.panel,
        }
    }
}

/// `HsMetricCollector.calculateHsPenalty`: the fold of sequencing needed to
/// bring 80% of target bases to `goal`x, searched by halving steps on the
/// pair multiplier until `estimateRoi` is within 0.1% of the goal. -1 when
/// it does not converge in 10000 steps, 0 without a library size.
fn hs_penalty(library_size: Option<i64>, m: &TargetMetrics, goal: i32) -> f64 {
    let Some(size) = library_size else {
        return 0.0;
    };
    let mean_coverage = m.on_target_from_pair_bases as f64 / m.panel.target_territory as f64;
    let fold80 = m.panel.fold_80_base_penalty;
    let (pairs, unique_pairs) = (m.pf_selected_pairs, m.pf_selected_unique_pairs);
    let on_target_pct = m.panel.on_target_bases as f64 / m.panel.pf_uq_bases_aligned as f64;
    let goal_multiplier = (f64::from(goal) / mean_coverage) * fold80;
    let mut pair_multiplier = goal_multiplier;
    let mut increment = 1.0;
    let mut going_up = goal_multiplier >= 1.0;
    let mut final_multiplier = -1.0;
    for _ in 0..10_000 {
        let unique_multiplier = estimate_roi(size, pair_multiplier, pairs, unique_pairs);
        if (unique_multiplier - goal_multiplier).abs() / goal_multiplier <= 0.001 {
            final_multiplier = pair_multiplier;
            break;
        } else if (unique_multiplier > goal_multiplier && going_up)
            || (unique_multiplier < goal_multiplier && !going_up)
        {
            increment /= 2.0;
            going_up = !going_up;
        }
        pair_multiplier += if going_up { increment } else { -increment };
    }
    if final_multiplier == -1.0 {
        return -1.0;
    }
    let unique_fraction =
        (unique_pairs as f64 * goal_multiplier) / (pairs as f64 * final_multiplier);
    (1.0 / unique_fraction) * fold80 * (1.0 / on_target_pct)
}

/// What one run produces: the finished units in file order and the four
/// shared histograms in their final state.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub variant: Variant,
    pub rows: Vec<TargetMetrics>,
    /// `highQualityDepthHistogram`: target bases per high-quality depth.
    pub depth: Histogram,
    /// `unCappedUnfilteredBaseQHistogram`: on-target bases per quality,
    /// every base above Q2.
    pub base_quality: Histogram,
    /// `unfilteredDepthHistogram`, capped at `COVERAGE_CAP`: the depth
    /// distribution of the sensitivity models.
    pub unfiltered_depth: Histogram,
    /// `unfilteredBaseQHistogram`: the base qualities counted under the cap.
    pub unfiltered_base_quality: Histogram,
}

impl Report {
    /// The `HsMetrics` rows, one per unit.
    pub fn hs_rows(&self) -> Vec<Row> {
        self.rows
            .iter()
            .cloned()
            .map(TargetMetrics::into_hs)
            .collect()
    }

    /// The `TargetedPcrMetrics` rows, one per unit.
    pub fn pcr_rows(&self) -> Vec<PcrRow> {
        self.rows
            .iter()
            .cloned()
            .map(TargetMetrics::into_pcr)
            .collect()
    }

    /// The `## HISTOGRAM` section: every unit adds the same two histograms
    /// (`addMetricsToFile`), so their columns repeat once per unit.
    pub fn histogram_section(&self) -> HistogramSection {
        let series = self
            .rows
            .iter()
            .flat_map(|_| {
                [
                    HistogramSeries {
                        label: "high_quality_coverage_count".to_owned(),
                        histogram: self.depth.clone(),
                    },
                    HistogramSeries {
                        label: "unfiltered_baseq_count".to_owned(),
                        histogram: self.base_quality.clone(),
                    },
                ]
            })
            .collect();
        HistogramSection {
            key_label: "coverage_or_base_quality",
            key_class: KeyClass::Integer,
            series,
        }
    }

    /// The `OUTPUT` file of the run's program.
    pub fn metrics_output(&self) -> MetricsOutput {
        let histogram = Some(self.histogram_section());
        match self.variant {
            Variant::Hs => MetricsOutput::from_rows(&self.hs_rows(), histogram),
            Variant::TargetedPcr => MetricsOutput::from_rows(&self.pcr_rows(), histogram),
        }
    }
}

/// Write the report as the Picard metrics file of its program.
pub fn write_metrics(
    report: &Report,
    headers: Vec<String>,
    out: &mut impl Write,
) -> io::Result<()> {
    report.metrics_output().write(headers, out)
}

/// `CollectHsMetrics.getProbeSetName`: `BAIT_SET_NAME`, or the distinct
/// bait file names cut at their first dot, sorted and joined with dots
/// (`CollectTargetedPcrMetrics.getProbeSetName` is the one-file case).
pub fn bait_set_name(opts: &Options) -> String {
    if let Some(name) = &opts.bait_set_name {
        return name.clone();
    }
    let names: std::collections::BTreeSet<String> = opts
        .bait_intervals
        .iter()
        .map(|p| {
            let file = p
                .file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_default();
            file.split('.').next().unwrap_or_default().to_owned()
        })
        .collect();
    names.into_iter().collect::<Vec<_>>().join(".")
}

/// The `TargetMetrics` accumulators of `PerUnitTargetMetricCollector`.
#[derive(Clone, Default)]
struct Counts {
    total_reads: i64,
    pf_reads: i64,
    pf_unique_reads: i64,
    pf_uq_reads_aligned: i64,
    pf_bases: i64,
    pf_bases_aligned: i64,
    pf_uq_bases_aligned: i64,
    pf_selected_pairs: i64,
    pf_selected_unique_pairs: i64,
    on_probe_bases: i64,
    near_probe_bases: i64,
    off_probe_bases: i64,
    on_target_bases: i64,
    on_target_from_pair_bases: i64,
    /// Aligned bases of duplicate reads (`PCT_EXC_DUPE` before division).
    excluded_dupe: i64,
    /// `CountingAdapterFilter.getFilteredBases`.
    excluded_adapter: i64,
    /// `CountingMapQFilter.getFilteredBases`.
    excluded_mapq: i64,
    excluded_baseq: i64,
    excluded_overlap: i64,
    excluded_off_target: i64,
}

impl Counts {
    /// The read-based and base-based totals of `acceptRecord`, counted
    /// before any filter so they agree with `CollectAlignmentSummaryMetrics`.
    fn add_read(&mut self, flags: Flags, aligned: i64, read_length: i64) {
        let (supplementary, qc_fail, duplicate, unmapped) = (
            flags.is_supplementary(),
            flags.is_qc_fail(),
            flags.is_duplicate(),
            flags.is_unmapped(),
        );
        if !supplementary {
            self.total_reads += 1;
            if !qc_fail {
                self.pf_reads += 1;
                if !duplicate {
                    self.pf_unique_reads += 1;
                    if !unmapped {
                        self.pf_uq_reads_aligned += 1;
                    }
                }
            }
        }
        if qc_fail {
            return;
        }
        // Supplementary records add aligned bases but not read bases.
        if !supplementary {
            self.pf_bases += read_length;
        }
        if !unmapped {
            self.pf_bases_aligned += aligned;
            if !duplicate {
                self.pf_uq_bases_aligned += aligned;
            }
        }
    }
}

/// Depths of every uniqued target, one array after another.
#[derive(Clone)]
struct Coverage {
    /// Start of each target's depths in `high_quality` and `unfiltered`.
    offsets: Vec<usize>,
    high_quality: Vec<i32>,
    unfiltered: Vec<i32>,
    /// `Coverage.readCount`: reads with a high-quality base on the target.
    read_count: Vec<i64>,
}

impl Coverage {
    fn new(targets: &[Interval]) -> Self {
        let mut offsets = Vec::with_capacity(targets.len() + 1);
        let mut total = 0usize;
        for target in targets {
            offsets.push(total);
            total += usize::try_from(target.length()).unwrap_or(0);
        }
        offsets.push(total);
        Self {
            offsets,
            high_quality: vec![0; total],
            unfiltered: vec![0; total],
            read_count: vec![0; targets.len()],
        }
    }

    fn range(&self, target: usize) -> Range<usize> {
        self.offsets[target]..self.offsets[target + 1]
    }
}

/// `Coverage.addBase`: one more read at a target base, unless the depth is
/// about to overflow an `int`.
fn add_depth(depth: &mut i32) -> i32 {
    if *depth < i32::MAX - 1 {
        *depth += 1;
    }
    *depth
}

/// The read's CIGAR as (operation, length) pairs.
type Cigar = Vec<(Kind, usize)>;

fn consumes_read(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match
            | Kind::Insertion
            | Kind::SoftClip
            | Kind::SequenceMatch
            | Kind::SequenceMismatch
    )
}

fn consumes_reference(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

/// `CigarOperator.isAlignment`: M, = and X, the operations of alignment blocks.
fn is_alignment(kind: Kind) -> bool {
    matches!(
        kind,
        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
    )
}

/// `SAMUtils.getNumOverlappingAlignedBasesToClip` for a mapped read with a
/// mapped mate: the read bases at or after the mate's start, counted only
/// on the left-most end (the first of pair when both start together). The
/// mates' contigs are not compared, as in htsjdk.
///
/// Per CIGAR element reaching the mate start: an M element adds the bases
/// from the mate start on (all of it when it starts after); S, H, P and N
/// add nothing; any other element (I, D, and also = and X) adds its whole
/// read length, so a straddling = or X element is clipped entirely.
fn overlapping_bases_to_clip(
    start: i64,
    mate_start: i64,
    first_of_pair: bool,
    cigar: &Cigar,
) -> i64 {
    if mate_start < start || (mate_start == start && first_of_pair) {
        return 0;
    }
    let mut clip = 0i64;
    let mut ref_pos = start;
    for &(kind, len) in cigar {
        let len = len as i64;
        let ref_len = if consumes_reference(kind) { len } else { 0 };
        if mate_start < ref_pos + ref_len {
            match kind {
                Kind::Match => {
                    clip += if mate_start < ref_pos {
                        ref_len
                    } else {
                        ref_pos + ref_len - mate_start
                    };
                }
                Kind::SoftClip | Kind::HardClip | Kind::Pad | Kind::Skip => {}
                _ => {
                    if consumes_read(kind) {
                        clip += len;
                    }
                }
            }
        }
        ref_pos += ref_len;
    }
    clip.max(0)
}

/// `CigarUtil.clipEndOfRead(clipFrom, cigar, SOFT_CLIP)`: soft-clip read
/// bases `clip_from..` (1-based). The element holding `clip_from - 1` ends
/// the copy; M, = and X keep their bases before the clip, I and S fold
/// them into the new S, and trailing hard clips are kept.
fn soft_clip_end(clip_from: i64, cigar: &Cigar) -> Result<Cigar> {
    let read_length: i64 = cigar
        .iter()
        .filter(|(k, _)| consumes_read(*k))
        .map(|&(_, n)| n as i64)
        .sum();
    let clipped = read_length - clip_from + 1;
    let trailing_hard = match cigar.last() {
        Some(&(Kind::HardClip, n)) => n,
        _ => 0,
    };
    let mut out = Vec::with_capacity(cigar.len() + 1);
    let mut pos = 1i64;
    for &(kind, len) in cigar {
        let length = if consumes_read(kind) { len as i64 } else { 0 };
        let end_pos = pos + length - 1;
        if end_pos < clip_from - 1 {
            out.push((kind, len));
            pos = end_pos + 1;
            continue;
        }
        // `mergeClippingCigarElement` with SOFT_CLIP.
        let relative = clip_from - pos;
        let mut amount = clipped;
        if consumes_read(kind) {
            if consumes_reference(kind) && relative > 0 {
                out.push((kind, relative as usize));
            }
            if !consumes_reference(kind) || kind == Kind::SoftClip {
                amount = clipped + relative;
            }
        } else if relative != 0 {
            return Err(Error::Metric(format!(
                "cannot clip a CIGAR at read base {clip_from} inside a {kind:?} element"
            )));
        }
        out.push((Kind::SoftClip, usize::try_from(amount).unwrap_or(0)));
        if trailing_hard > 0 {
            out.push((Kind::HardClip, trailing_hard));
        }
        break;
    }
    Ok(out)
}

/// Bases of alignment blocks that fall inside `baits` (`onBaitBases`).
fn bases_on_baits(start: i64, cigar: &Cigar, baits: &[Interval]) -> i64 {
    let mut on = 0;
    let mut ref_pos = start;
    for &(kind, len) in cigar {
        let len = len as i64;
        if is_alignment(kind) {
            let block_end = ref_pos + len - 1;
            for bait in baits {
                on += (block_end.min(bait.end) - ref_pos.max(bait.start) + 1).max(0);
            }
        }
        if consumes_reference(kind) {
            ref_pos += len;
        }
    }
    on
}

/// The per-unit state of `PerUnitTargetMetricCollector`.
#[derive(Clone)]
struct Unit {
    coverage: Coverage,
    hq_max_depth: i64,
    /// `baseQHistogramArray`: qualities of bases added while the unfiltered depth is within the cap.
    base_q_capped: [i64; QUALITY_SLOTS],
    /// `unCappedBaseQHistogramArray`.
    base_q: [i64; QUALITY_SLOTS],
    counts: Counts,
    /// `coveredTargets` of the record being added: targets that already
    /// counted it as a read.
    covered: Vec<usize>,
}

impl Unit {
    fn new(targets: &[Interval]) -> Self {
        Self {
            coverage: Coverage::new(targets),
            hq_max_depth: 0,
            base_q_capped: [0; QUALITY_SLOTS],
            base_q: [0; QUALITY_SLOTS],
            counts: Counts::default(),
            covered: Vec::new(),
        }
    }
}

/// What the per-base loop of `acceptRecord` needs from the record and the
/// collector, read-only, so the loop can run once per unit.
struct BaseLoop<'a> {
    record: &'a bam::Record,
    cigar: &'a Cigar,
    start: i64,
    /// The record's overlapping targets, as an index range and the slice.
    first_target: usize,
    targets: &'a [Interval],
    mapped_in_pair: bool,
    minimum_base_quality: i32,
    coverage_cap: i32,
    include_indels: bool,
}

impl Unit {
    /// The per-base loop of `acceptRecord` over the (clipped) CIGAR:
    /// summary counts for every aligned base, then target depths and
    /// quality histograms for on-target bases above Q2.
    fn add_coverage(&mut self, b: &BaseLoop<'_>) -> Result<()> {
        let malformed = |what| Error::MalformedRecord {
            name: record_name(b.record),
            what,
        };
        let qualities = b.record.quality_scores();
        let qualities = qualities.as_ref();
        // htsjdk reads a 0xFF quality string as no qualities.
        let qualities = if qualities.first() == Some(&0xFF) {
            &[][..]
        } else {
            qualities
        };
        self.covered.clear();
        let mut cursor = 0usize;
        let (mut read_offset, mut ref_pos) = (0usize, b.start);
        for &(kind, len) in b.cigar {
            let counted = is_alignment(kind)
                || (b.include_indels && matches!(kind, Kind::Insertion | Kind::Deletion));
            let (read_step, ref_step) = (
                usize::from(consumes_read(kind)),
                i64::from(consumes_reference(kind)),
            );
            if !counted {
                read_offset += read_step * len;
                ref_pos += ref_step * len as i64;
                continue;
            }
            for _ in 0..len {
                // A Java byte: qualities above 127 are negative.
                let quality = i32::from(
                    *qualities
                        .get(read_offset)
                        .ok_or_else(|| malformed("aligned base without a quality"))?
                        as i8,
                );
                // Targets are sorted and disjoint, and ref_pos never moves
                // back, so one cursor finds the target holding ref_pos.
                while cursor < b.targets.len() && b.targets[cursor].end < ref_pos {
                    cursor += 1;
                }
                let target = (cursor < b.targets.len() && b.targets[cursor].start <= ref_pos)
                    .then_some(cursor);
                let high_quality = quality >= b.minimum_base_quality;
                let c = &mut self.counts;
                if !high_quality {
                    c.excluded_baseq += 1;
                } else if target.is_none() {
                    c.excluded_off_target += 1;
                } else {
                    c.on_target_bases += 1;
                    if b.mapped_in_pair {
                        c.on_target_from_pair_bases += 1;
                    }
                }
                // An inserted base has no target position.
                if let Some(t) =
                    target.filter(|_| quality > MIN_UNFILTERED_QUALITY && kind != Kind::Insertion)
                {
                    let slot = usize::try_from(quality)
                        .ok()
                        .filter(|&q| q < QUALITY_SLOTS)
                        .ok_or_else(|| malformed("base quality above 126"))?;
                    let index = b.first_target + t;
                    let offset = self.coverage.offsets[index]
                        + usize::try_from(ref_pos - b.targets[t].start).unwrap_or(0);
                    if add_depth(&mut self.coverage.unfiltered[offset]) <= b.coverage_cap {
                        self.base_q_capped[slot] += 1;
                    }
                    self.base_q[slot] += 1;
                    if high_quality {
                        let depth = add_depth(&mut self.coverage.high_quality[offset]);
                        self.hq_max_depth = self.hq_max_depth.max(i64::from(depth));
                        if !self.covered.contains(&index) {
                            self.covered.push(index);
                            self.coverage.read_count[index] += 1;
                        }
                    }
                }
                read_offset += read_step;
                ref_pos += ref_step;
            }
        }
        Ok(())
    }
}

/// What `emitPerBaseCoverageIfRequested` and the per-target writer need
/// besides the depths.
struct FileNames {
    /// `@SQ` names, by contig index.
    contigs: Vec<String>,
    /// Name of each uniqued target, as `IntervalList.uniqued()` merges them.
    target_names: Vec<String>,
}

/// `PerUnitTargetMetricCollector.emitPerBaseCoverageIfRequested`: a
/// `chrom pos target coverage` header, then one line per target base in
/// uniqued target order, with the 1-based position and the high-quality
/// depth. htsjdk compresses when the file name ends in `.gz`, `.gzip` or
/// `.bfq`; this writes BGZF there, which gzip readers take as plain gzip.
fn write_per_base_coverage(
    path: &Path,
    names: &FileNames,
    targets: &[Interval],
    coverage: &Coverage,
) -> Result<()> {
    let output_error = |source| Error::Output {
        path: path.to_path_buf(),
        source,
    };
    let file = std::fs::File::create(path).map_err(output_error)?;
    let name = path.file_name().map(|n| n.to_string_lossy());
    let compressed = name.is_some_and(|n| [".gz", ".gzip", ".bfq"].iter().any(|e| n.ends_with(e)));
    let write = |out: &mut dyn Write| -> io::Result<()> {
        out.write_all(b"chrom\tpos\ttarget\tcoverage\n")?;
        for (i, target) in targets.iter().enumerate() {
            let contig = names.contigs[target.contig].as_bytes();
            let name = names.target_names[i].as_bytes();
            for (offset, depth) in coverage.high_quality[coverage.range(i)].iter().enumerate() {
                out.write_all(contig)?;
                write!(out, "\t{}\t", target.start + offset as i64)?;
                out.write_all(name)?;
                writeln!(out, "\t{depth}")?;
            }
        }
        out.flush()
    };
    if compressed {
        let mut out = noodles::bgzf::io::Writer::new(file);
        write(&mut out).map_err(output_error)?;
        out.finish().map_err(output_error)?;
    } else {
        write(&mut io::BufWriter::with_capacity(1 << 16, file)).map_err(output_error)?;
    }
    Ok(())
}

/// The `PER_TARGET_COVERAGE` writer of `calculateGcMetrics`: one line per
/// uniqued target with its GC fraction, mean and normalised coverage, the
/// extremes, the fraction of bases at 0x and the read count, every number
/// through `FormatUtil`.
struct PerTarget {
    out: io::BufWriter<std::fs::File>,
    path: PathBuf,
}

impl PerTarget {
    fn create(path: &Path) -> Result<Self> {
        let file = std::fs::File::create(path).map_err(|source| Error::Output {
            path: path.to_path_buf(),
            source,
        })?;
        let mut writer = Self {
            out: io::BufWriter::new(file),
            path: path.to_path_buf(),
        };
        writer.write(|out| {
            out.write_all(b"chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage\tnormalized_coverage\tmin_normalized_coverage\tmax_normalized_coverage\tmin_coverage\tmax_coverage\tpct_0x\tread_count\n")
        })?;
        Ok(writer)
    }

    fn write(&mut self, f: impl FnOnce(&mut dyn Write) -> io::Result<()>) -> Result<()> {
        f(&mut self.out).map_err(|source| Error::Output {
            path: self.path.clone(),
            source,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn line(
        &mut self,
        contig: &str,
        target: &Interval,
        name: &str,
        gc: f64,
        depths: &[i32],
        mean_target_coverage: f64,
        read_count: i64,
    ) -> Result<()> {
        let coverage =
            depths.iter().map(|&d| i64::from(d)).sum::<i64>() as f64 / target.length() as f64;
        let mut min = f64::from(i32::MAX);
        let mut max = f64::from(i32::MIN);
        let mut at_0x = 0.0;
        for &d in depths {
            if d == 0 {
                at_0x += 1.0;
            }
            let d = f64::from(d);
            if d < min {
                min = d;
            }
            if max < d {
                max = d;
            }
        }
        let fields = [
            format_double(gc),
            format_double(coverage),
            format_double(coverage / mean_target_coverage),
            format_double(min / mean_target_coverage),
            format_double(max / mean_target_coverage),
            format_double(min),
            format_double(max),
            format_double(at_0x / target.length() as f64),
            read_count.to_string(),
        ];
        self.write(|out| {
            writeln!(
                out,
                "{contig}\t{}\t{}\t{}\t{name}\t{}",
                target.start,
                target.end,
                target.length(),
                fields.join("\t")
            )
        })
    }

    fn finish(mut self) -> Result<()> {
        self.write(|out| out.flush())
    }
}

/// Inputs to `finish` that do not come from the reads.
struct Territory {
    probe_set: String,
    probe_territory: i64,
    target_territory: i64,
    genome_size: i64,
    /// `allTargets.getIntervals().size()`: target intervals before merging,
    /// the denominator of `ZERO_CVG_TARGETS_PCT`.
    raw_target_count: usize,
    /// GC fraction of each uniqued target, when a reference was given.
    target_gc: Option<Vec<f64>>,
    names: FileNames,
}

/// The collector-level histograms every unit's `finish` adds to.
#[derive(Default)]
struct Shared {
    hq_depth: Histogram,
    unfiltered_depth: Histogram,
    unfiltered_base_quality: Histogram,
    uncapped_base_quality: Histogram,
}

/// `SequenceUtil.calculateGc`.
fn gc_fraction(bases: &[u8]) -> f64 {
    let gc = bases
        .iter()
        .filter(|&&b| matches!(b, b'C' | b'G' | b'c' | b'g'))
        .count();
    gc as f64 / bases.len() as f64
}

/// Which files this unit writes: only the `ALL_READS` collector gets the
/// per-target and per-base paths (`makeAllReadCollector`).
struct UnitFiles<'a> {
    per_base: Option<&'a Path>,
    per_target: Option<&'a Path>,
}

impl Unit {
    /// `PerUnitTargetMetricCollector.finish`: the ratios, then the three
    /// calculations in Picard's order, then the per-base file.
    fn finish(
        self,
        opts: &Options,
        territory: &Territory,
        targets: &[Interval],
        labels: Labels,
        files: UnitFiles<'_>,
        shared: &mut Shared,
    ) -> Result<TargetMetrics> {
        let mut m = self.summary_metrics(territory, labels);
        self.coverage_metrics(&mut m.panel, targets, territory, shared)?;
        self.het_sensitivity(&mut m.panel, opts, targets, shared)?;
        self.gc_metrics(&mut m.panel, territory, targets, files.per_target)?;
        // Picard writes the per-base file at the end of `finish`, so a unit
        // that fails above writes none.
        if let Some(path) = files.per_base {
            write_per_base_coverage(path, &territory.names, targets, &self.coverage)?;
        }
        Ok(m)
    }

    /// The counts and the ratios `finish` derives from them before any
    /// histogram is touched.
    fn summary_metrics(&self, territory: &Territory, labels: Labels) -> TargetMetrics {
        let c = &self.counts;
        let mut m = TargetMetrics {
            probe_set: territory.probe_set.clone(),
            probe_territory: territory.probe_territory,
            on_probe_bases: c.on_probe_bases,
            near_probe_bases: c.near_probe_bases,
            off_probe_bases: c.off_probe_bases,
            pf_selected_pairs: c.pf_selected_pairs,
            pf_selected_unique_pairs: c.pf_selected_unique_pairs,
            on_target_from_pair_bases: c.on_target_from_pair_bases,
            panel: Panel {
                target_territory: territory.target_territory,
                genome_size: territory.genome_size,
                total_reads: c.total_reads,
                pf_reads: c.pf_reads,
                pf_bases: c.pf_bases,
                pf_unique_reads: c.pf_unique_reads,
                pf_uq_reads_aligned: c.pf_uq_reads_aligned,
                pf_bases_aligned: c.pf_bases_aligned,
                pf_uq_bases_aligned: c.pf_uq_bases_aligned,
                on_target_bases: c.on_target_bases,
                labels,
                ..Panel::default()
            },
            ..TargetMetrics::default()
        };
        let p = &mut m.panel;
        p.pct_pf_reads = c.pf_reads as f64 / c.total_reads as f64;
        p.pct_pf_uq_reads = c.pf_unique_reads as f64 / c.total_reads as f64;
        p.pct_pf_uq_reads_aligned = c.pf_uq_reads_aligned as f64 / c.pf_unique_reads as f64;
        let denominator = (c.on_probe_bases + c.near_probe_bases + c.off_probe_bases) as f64;
        m.pct_selected_bases = (c.on_probe_bases + c.near_probe_bases) as f64 / denominator;
        m.pct_off_probe = c.off_probe_bases as f64 / denominator;
        m.on_probe_vs_selected =
            c.on_probe_bases as f64 / (c.on_probe_bases + c.near_probe_bases) as f64;
        m.mean_probe_coverage = c.on_probe_bases as f64 / m.probe_territory as f64;
        m.fold_enrichment = (c.on_probe_bases as f64 / denominator)
            / (m.probe_territory as f64 / p.genome_size as f64);
        let aligned = c.pf_bases_aligned as f64;
        p.pct_exc_dupe = c.excluded_dupe as f64 / aligned;
        p.pct_exc_adapter = c.excluded_adapter as f64 / aligned;
        p.pct_exc_mapq = c.excluded_mapq as f64 / aligned;
        p.pct_exc_baseq = c.excluded_baseq as f64 / aligned;
        p.pct_exc_overlap = c.excluded_overlap as f64 / aligned;
        p.pct_exc_off_target = c.excluded_off_target as f64 / aligned;
        m
    }

    /// `calculateTargetCoverageMetrics`. Depths are tallied per value, then
    /// added to the shared histogram once per value: the same integer counts
    /// Picard adds base by base. The median and the 20th percentile come
    /// from the shared histogram as it stands after this unit's additions.
    fn coverage_metrics(
        &self,
        p: &mut Panel,
        targets: &[Interval],
        territory: &Territory,
        shared: &mut Shared,
    ) -> Result<()> {
        for i in 0..self.hq_max_depth {
            shared.hq_depth.increment_by(i, 0.0);
        }
        let hq_max = usize::try_from(self.hq_max_depth).unwrap_or(0);
        let mut depth_counts = vec![0i64; hq_max + 1];
        let mut zero_targets = 0usize;
        let mut total_coverage = 0i64;
        let mut min_depth = i64::MAX;
        for (i, target) in targets.iter().enumerate() {
            let depths = &self.coverage.high_quality[self.coverage.range(i)];
            if !depths.iter().any(|&d| d > 0) {
                zero_targets += 1;
                depth_counts[0] += target.length();
                min_depth = 0;
                continue;
            }
            for &depth in depths {
                depth_counts[depth as usize] += 1;
                total_coverage += i64::from(depth);
                min_depth = min_depth.min(i64::from(depth));
            }
        }
        for (depth, &count) in depth_counts.iter().enumerate() {
            if count > 0 {
                shared.hq_depth.increment_by(depth as i64, count as f64);
            }
        }
        // Bases at or above each TARGET_BASE_DEPTHS entry.
        let mut target_bases = [0i64; TARGET_BASE_DEPTHS.len()];
        for (d, &count) in depth_counts.iter().enumerate() {
            for (bases, &threshold) in target_bases.iter_mut().zip(&TARGET_BASE_DEPTHS) {
                if d as i64 >= threshold {
                    *bases += count;
                }
            }
        }
        p.mean_target_coverage = total_coverage as f64 / p.target_territory as f64;
        p.median_target_coverage = shared.hq_depth.median();
        p.max_target_coverage = self.hq_max_depth;
        p.min_target_coverage = min_depth.min(self.hq_max_depth);
        let percentile_20 = shared.hq_depth.percentile(0.2).ok_or_else(|| {
            Error::Metric("Cannot calculate percentiles when total is zero.".to_owned())
        })?;
        p.fold_80_base_penalty = p.mean_target_coverage / percentile_20;
        p.zero_cvg_targets_pct = zero_targets as f64 / territory.raw_target_count as f64;
        for (pct, &bases) in p.pct_target_bases.iter_mut().zip(&target_bases[1..]) {
            *pct = bases as f64 / target_bases[0] as f64;
        }
        Ok(())
    }

    /// `calculateTheoreticalHetSensitivity`: the capped unfiltered depths
    /// and the base qualities counted under the cap go into the shared
    /// histograms, the two totals must agree, and the model runs on the
    /// shared histograms as they stand.
    fn het_sensitivity(
        &self,
        p: &mut Panel,
        opts: &Options,
        targets: &[Interval],
        shared: &mut Shared,
    ) -> Result<()> {
        let cap = i64::from(opts.coverage_cap);
        for i in 0..cap {
            shared.unfiltered_depth.increment_by(i, 0.0);
        }
        for (i, target) in targets.iter().enumerate() {
            let depths = &self.coverage.unfiltered[self.coverage.range(i)];
            if !depths.iter().any(|&d| d > 0) {
                shared
                    .unfiltered_depth
                    .increment_by(0, target.length() as f64);
                continue;
            }
            for &d in depths {
                shared.unfiltered_depth.increment(i64::from(d).min(cap));
            }
        }
        if self.base_q_capped.iter().sum::<i64>() as f64 != shared.unfiltered_depth.sum() {
            return Err(Error::Metric(
                "numbers of bases in the base quality histogram and the coverage histogram are not equal"
                    .to_owned(),
            ));
        }
        for q in 0..QUALITY_SLOTS {
            shared
                .unfiltered_base_quality
                .increment_by(q as i64, self.base_q_capped[q] as f64);
            shared
                .uncapped_base_quality
                .increment_by(q as i64, self.base_q[q] as f64);
        }
        p.het_snp_sensitivity = het_snp_sensitivity(
            &normalize_histogram(&shared.unfiltered_depth),
            &normalize_histogram(&shared.unfiltered_base_quality),
            opts.sample_size,
            LOG_ODDS_THRESHOLD,
        )?;
        p.het_snp_q = f64::from(phred_from_error_probability(1.0 - p.het_snp_sensitivity));
        Ok(())
    }

    /// `calculateGcMetrics`, only with a reference: `AT_DROPOUT` and
    /// `GC_DROPOUT` over 101 GC bins, with the per-target file written on the
    /// way through.
    fn gc_metrics(
        &self,
        p: &mut Panel,
        territory: &Territory,
        targets: &[Interval],
        per_target: Option<&Path>,
    ) -> Result<()> {
        if let Some(gc) = &territory.target_gc {
            let mut per_target = per_target.map(PerTarget::create).transpose()?;
            let mut target_by_gc = [0i64; GC_BINS];
            let mut aligned_by_gc = [0i64; GC_BINS];
            for (i, target) in targets.iter().enumerate() {
                if target.length() <= 0 {
                    continue;
                }
                let depths = &self.coverage.high_quality[self.coverage.range(i)];
                let bin = usize::try_from(math_round(gc[i] * 100.0)).unwrap_or(0);
                target_by_gc[bin] += target.length();
                aligned_by_gc[bin] += depths.iter().map(|&d| i64::from(d)).sum::<i64>();
                if let Some(out) = &mut per_target {
                    out.line(
                        &territory.names.contigs[target.contig],
                        target,
                        &territory.names.target_names[i],
                        gc[i],
                        depths,
                        p.mean_target_coverage,
                        self.coverage.read_count[i],
                    )?;
                }
            }
            if let Some(out) = per_target {
                out.finish()?;
            }
            let total_target: i64 = target_by_gc.iter().sum();
            let total_aligned: i64 = aligned_by_gc.iter().sum();
            for i in 0..GC_BINS {
                let target_pct = target_by_gc[i] as f64 / total_target as f64;
                let aligned_pct = aligned_by_gc[i] as f64 / total_aligned as f64;
                let dropout = (aligned_pct - target_pct) * 100.0;
                if dropout < 0.0 {
                    if i <= 50 {
                        p.at_dropout += dropout.abs();
                    }
                    if i >= 50 {
                        p.gc_dropout += dropout.abs();
                    }
                }
            }
        }
        Ok(())
    }
}

/// `CollectHsMetrics` or `CollectTargetedPcrMetrics` as a pass [`Accumulator`].
pub struct CaptureMetrics {
    opts: Options,
    levels: Levels<Unit>,
    targets: OverlapIndex,
    probes: OverlapIndex,
    adapters: Adapters,
    territory: Territory,
    cigar: Cigar,
    bases: Vec<u8>,
}

impl CaptureMetrics {
    /// Reads and uniques the interval files against the BAM header, and
    /// with a reference computes the GC fraction of every target.
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        opts.check_interval_files()?;
        let contigs = header.reference_sequences().len();
        let raw_targets = intervals::read_named_files(&opts.target_intervals, header)?;
        let unique_named = intervals::uniqued_named(&raw_targets);
        let unique_targets: Vec<Interval> = unique_named.iter().map(|n| n.interval).collect();
        let unique_baits =
            intervals::uniqued(&intervals::read_files(&opts.bait_intervals, header)?);
        let target_gc = match &opts.reference {
            None => None,
            Some(reference) => {
                let mut reference = Reference::open(reference, header)?;
                let gc = unique_targets
                    .iter()
                    .map(|t| {
                        Ok(gc_fraction(
                            &reference.subsequence(t.contig, t.start, t.end)?,
                        ))
                    })
                    .collect::<Result<Vec<f64>>>()?;
                Some(gc)
            }
        };
        let territory = Territory {
            probe_set: bait_set_name(opts),
            probe_territory: intervals::count_bases(&unique_baits),
            target_territory: intervals::count_bases(&unique_targets),
            genome_size: header
                .reference_sequences()
                .values()
                .map(|c| c.length().get() as i64)
                .sum(),
            raw_target_count: raw_targets.len(),
            target_gc,
            names: FileNames {
                contigs: header
                    .reference_sequences()
                    .keys()
                    .map(|name| name.to_string())
                    .collect(),
                target_names: unique_named.into_iter().map(|n| n.name).collect(),
            },
        };
        Ok(Self {
            opts: opts.clone(),
            levels: Levels::new(header, &opts.levels, Unit::new(&unique_targets)),
            targets: OverlapIndex::new(unique_targets, contigs, 0),
            probes: OverlapIndex::new(unique_baits, contigs, i64::from(opts.near_distance)),
            adapters: Adapters::default(),
            territory,
            cigar: Vec::new(),
            bases: Vec::new(),
        })
    }
}

impl CaptureMetrics {
    /// The bait counts of `acceptRecord`, taken before the duplicate,
    /// adapter, MAPQ and overlap filters so they describe the assay rather
    /// than the usable data: `PF_SELECTED_PAIRS` (first ends of mapped pairs
    /// near a bait, for `HS_LIBRARY_SIZE`), the on, near and off bait bases,
    /// and the duplicate's bases under `PCT_EXC_DUPE`.
    fn bait_counts(
        &mut self,
        record: &bam::Record,
        start: i64,
        aligned: i64,
        probes: Range<usize>,
    ) {
        let flags = record.flags();
        let duplicate = flags.is_duplicate();
        let selected_pair = !flags.is_supplementary()
            && flags.is_segmented()
            && flags.is_first_segment()
            && !flags.is_mate_unmapped()
            && !probes.is_empty();
        let on_bait = if probes.is_empty() {
            None
        } else {
            Some(bases_on_baits(
                start,
                &self.cigar,
                &self.probes.intervals()[probes],
            ))
        };
        self.levels.route(record, |u| {
            let c = &mut u.counts;
            if selected_pair {
                c.pf_selected_pairs += 1;
                if !duplicate {
                    c.pf_selected_unique_pairs += 1;
                }
            }
            match on_bait {
                None => c.off_probe_bases += aligned,
                Some(on) => {
                    c.on_probe_bases += on;
                    c.near_probe_bases += aligned - on;
                }
            }
            if duplicate {
                c.excluded_dupe += aligned;
            }
        });
    }

    /// The overlap clipping of `acceptRecord`: count the clipped bases in
    /// every unit, then soft-clip the CIGAR. False when the read starts at
    /// or after its mate and loses every base (`clipOverlappingAlignedBases`
    /// unmaps it), so it adds no coverage.
    fn clip_overlap(&mut self, record: &bam::Record, start: i64, read_length: i64) -> Result<bool> {
        let mate_start = match record.mate_alignment_start().transpose()? {
            Some(p) => p.get() as i64,
            None => 0,
        };
        let first_of_pair = record.flags().is_first_segment();
        let clip = overlapping_bases_to_clip(start, mate_start, first_of_pair, &self.cigar);
        self.levels
            .route(record, |u| u.counts.excluded_overlap += clip);
        if clip == 0 {
            return Ok(true);
        }
        if mate_start <= start {
            return Ok(false);
        }
        // Clip from the base after the last kept one, before any trailing
        // soft clip (a trailing H is not looked past).
        let mut clip_from = read_length - clip + 1;
        if let Some(&(Kind::SoftClip, n)) = self.cigar.last() {
            clip_from -= n as i64;
        }
        self.cigar = soft_clip_end(clip_from, &self.cigar)?;
        Ok(true)
    }
}

impl Accumulator for CaptureMetrics {
    type Output = Report;

    /// `PerUnitTargetMetricCollector.acceptRecord`, with the record-level
    /// decisions made once and applied to every unit the record belongs to.
    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if flags.is_secondary() {
            return Ok(());
        }
        let unmapped = flags.is_unmapped();
        let supplementary = flags.is_supplementary();
        let duplicate = flags.is_duplicate();
        let paired = flags.is_segmented();
        let mate_unmapped = flags.is_mate_unmapped();
        let qc_fail = flags.is_qc_fail();
        let mapped_in_pair = paired && !unmapped && !mate_unmapped && !supplementary;
        self.cigar.clear();
        for op in record.cigar().iter() {
            let op = op?;
            self.cigar.push((op.kind(), op.len()));
        }
        let aligned: i64 = if unmapped {
            0
        } else {
            self.cigar
                .iter()
                .filter(|(k, _)| is_alignment(*k))
                .map(|&(_, n)| n as i64)
                .sum()
        };
        let read_length = record.sequence().len() as i64;
        self.levels.route(record, |u| {
            u.counts.add_read(flags, aligned, read_length);
        });
        if qc_fail || unmapped {
            return Ok(());
        }

        let contig =
            record
                .reference_sequence_id()
                .transpose()?
                .ok_or_else(|| Error::MalformedRecord {
                    name: record_name(record),
                    what: "mapped read without a reference sequence",
                })?;
        let start = position(record, record.alignment_start(), "missing alignment start")?;
        let reference_length: i64 = self
            .cigar
            .iter()
            .filter(|(k, _)| consumes_reference(*k))
            .map(|&(_, n)| n as i64)
            .sum();
        let end = start + reference_length - 1;
        let targets = self.targets.overlaps(contig, start, end);
        let probes = self.probes.overlaps(contig, start, end);

        self.bait_counts(record, start, aligned, probes);
        if duplicate {
            return Ok(());
        }
        // htsjdk reports a missing MAPQ as 255.
        let mapq = record.mapping_quality().map_or(255, |q| q.get());
        if mapq == 0 {
            self.bases.clear();
            self.bases.extend(record.sequence().iter());
            if self
                .adapters
                .is_adapter(&self.bases, false, flags.is_reverse_complemented(), mapq)
            {
                self.levels
                    .route(record, |u| u.counts.excluded_adapter += aligned);
                return Ok(());
            }
        }
        if i32::from(mapq) < self.opts.minimum_mapping_quality {
            self.levels
                .route(record, |u| u.counts.excluded_mapq += aligned);
            return Ok(());
        }
        if self.opts.clip_overlapping_reads
            && paired
            && !mate_unmapped
            && !self.clip_overlap(record, start, read_length)?
        {
            return Ok(());
        }
        let Self {
            levels,
            cigar,
            targets: index,
            opts,
            ..
        } = self;
        let base_loop = BaseLoop {
            record,
            cigar,
            start,
            first_target: targets.start,
            targets: &index.intervals()[targets],
            mapped_in_pair,
            minimum_base_quality: opts.minimum_base_quality,
            coverage_cap: opts.coverage_cap,
            include_indels: opts.include_indels,
        };
        let mut result = Ok(());
        levels.route(record, |u| {
            if result.is_ok() {
                result = u.add_coverage(&base_loop);
            }
        });
        result
    }

    /// `MultiLevelCollector.finish` then `addAllLevelsToFile`: units in
    /// file order, each adding to the shared histograms.
    fn finish(self) -> Result<Report> {
        let targets = self.targets.intervals();
        let mut shared = Shared::default();
        let mut rows = Vec::new();
        for bucket in self.levels.into_buckets() {
            let all_reads = bucket.level == Level::AllReads;
            let files = UnitFiles {
                per_base: self.opts.per_base_coverage.as_deref().filter(|_| all_reads),
                per_target: self
                    .opts
                    .per_target_coverage
                    .as_deref()
                    .filter(|_| all_reads),
            };
            rows.push(bucket.state.finish(
                &self.opts,
                &self.territory,
                targets,
                bucket.labels,
                files,
                &mut shared,
            )?);
        }
        Ok(Report {
            variant: self.opts.variant,
            rows,
            depth: shared.hq_depth,
            base_quality: shared.uncapped_base_quality,
            unfiltered_depth: shared.unfiltered_depth,
            unfiltered_base_quality: shared.unfiltered_base_quality,
        })
    }
}

/// `StringHeader` command line in Picard's shape, for the native
/// `hs-metrics` subcommand.
pub fn command_header(opts: &Options, bam: &Path, output: Option<&Path>) -> String {
    let list = |paths: &[PathBuf]| {
        paths
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(", ")
    };
    let mut command = format!(
        "CollectHsMetrics BAIT_INTERVALS=[{}] TARGET_INTERVALS=[{}] INPUT={}",
        list(&opts.bait_intervals),
        list(&opts.target_intervals),
        bam.display()
    );
    if let Some(output) = output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(per_base) = &opts.per_base_coverage {
        command.push_str(&format!(" PER_BASE_COVERAGE={}", per_base.display()));
    }
    if let Some(name) = &opts.bait_set_name {
        command.push_str(&format!(" BAIT_SET_NAME={name}"));
    }
    if let Some(reference) = &opts.reference {
        command.push_str(&format!(" REFERENCE_SEQUENCE={}", reference.display()));
    }
    let levels: Vec<&str> = opts.levels.iter().map(|l| l.name()).collect();
    command.push_str(&format!(
        " METRIC_ACCUMULATION_LEVEL=[{}] NEAR_DISTANCE={} MINIMUM_MAPPING_QUALITY={} MINIMUM_BASE_QUALITY={} CLIP_OVERLAPPING_READS={} INCLUDE_INDELS={} COVERAGE_CAP={} SAMPLE_SIZE={}",
        levels.join(", "),
        opts.near_distance,
        opts.minimum_mapping_quality,
        opts.minimum_base_quality,
        opts.clip_overlapping_reads,
        opts.include_indels,
        opts.coverage_cap,
        opts.sample_size
    ));
    command
}

/// Compute the metrics of one BAM with the default BGZF worker count.
pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Report> {
    pass::run(path.as_ref(), pass::default_threads(), |header| {
        CaptureMetrics::new(opts, header)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cigar(text: &str) -> Cigar {
        let mut out = Vec::new();
        let mut n = 0usize;
        for c in text.bytes() {
            if c.is_ascii_digit() {
                n = n * 10 + usize::from(c - b'0');
                continue;
            }
            let kind = match c {
                b'M' => Kind::Match,
                b'I' => Kind::Insertion,
                b'D' => Kind::Deletion,
                b'N' => Kind::Skip,
                b'S' => Kind::SoftClip,
                b'H' => Kind::HardClip,
                b'P' => Kind::Pad,
                b'=' => Kind::SequenceMatch,
                _ => Kind::SequenceMismatch,
            };
            out.push((kind, n));
            n = 0;
        }
        out
    }

    #[test]
    fn overlap_clip_counts_from_the_mate_start() {
        // Read 100..199, mate at 150: 50 bases overlap.
        assert_eq!(
            overlapping_bases_to_clip(100, 150, true, &cigar("100M")),
            50
        );
        // The right-most end is never clipped; equal starts clip the second end.
        assert_eq!(
            overlapping_bases_to_clip(150, 100, false, &cigar("100M")),
            0
        );
        assert_eq!(overlapping_bases_to_clip(100, 100, true, &cigar("100M")), 0);
        assert_eq!(
            overlapping_bases_to_clip(100, 100, false, &cigar("100M")),
            100
        );
        // No overlap.
        assert_eq!(overlapping_bases_to_clip(100, 300, true, &cigar("100M")), 0);
        // An insertion after the mate start counts, a deletion does not.
        assert_eq!(
            overlapping_bases_to_clip(100, 150, true, &cigar("60M5I2D35M")),
            10 + 5 + 35
        );
        // An = element straddling the mate start is counted whole.
        assert_eq!(
            overlapping_bases_to_clip(100, 150, true, &cigar("100=")),
            100
        );
        // Soft clips add nothing.
        assert_eq!(
            overlapping_bases_to_clip(100, 150, true, &cigar("10S90M")),
            40
        );
    }

    #[test]
    fn soft_clip_end_follows_cigar_util() {
        assert_eq!(
            soft_clip_end(51, &cigar("100M")).expect("clips"),
            cigar("50M50S")
        );
        assert_eq!(
            soft_clip_end(51, &cigar("90M10S")).expect("clips"),
            cigar("50M50S")
        );
        // An insertion straddling the clip folds into the soft clip.
        assert_eq!(
            soft_clip_end(48, &cigar("45M5I50M")).expect("clips"),
            cigar("45M55S")
        );
        // Trailing hard clips stay.
        assert_eq!(
            soft_clip_end(51, &cigar("100M5H")).expect("clips"),
            cigar("50M50S5H")
        );
        // A clip right after a deletion lands on the M before it, so the
        // deletion is dropped.
        assert_eq!(
            soft_clip_end(51, &cigar("50M3D50M")).expect("clips"),
            cigar("50M50S")
        );
    }

    #[test]
    fn on_bait_counts_block_bases_inside_baits() {
        let baits = [
            Interval {
                contig: 0,
                start: 90,
                end: 109,
            },
            Interval {
                contig: 0,
                start: 150,
                end: 160,
            },
        ];
        // Blocks 100..149 and 160..209 (after a 10 bp deletion).
        assert_eq!(bases_on_baits(100, &cigar("50M10D50M"), &baits), 10 + 1);
    }

    #[test]
    fn bait_set_name_is_the_file_stem() {
        let opts = Options {
            bait_intervals: vec![
                PathBuf::from("/x/panel.v2.interval_list"),
                PathBuf::from("a.bed"),
            ],
            ..Options::default()
        };
        assert_eq!(bait_set_name(&opts), "a.panel");
        let named = Options {
            bait_set_name: Some("kit".to_owned()),
            ..opts
        };
        assert_eq!(bait_set_name(&named), "kit");
    }

    #[test]
    fn rows_match_their_columns() {
        let hs = Row {
            bait_set: "baits".to_owned(),
            hs_library_size: None,
            panel: Panel {
                het_snp_q: 22.0,
                labels: Labels {
                    sample: Some("S1".to_owned()),
                    ..Labels::default()
                },
                ..Panel::default()
            },
            ..Row::default()
        };
        let values = hs.values();
        assert_eq!(values.len(), COLUMNS.len());
        assert_eq!(values[0].picard_text(), "baits");
        assert_eq!(values[13].picard_text(), "");
        assert_eq!(values[65].picard_text(), "22");
        assert_eq!(values[66].picard_text(), "S1");
        assert_eq!(values[67].picard_text(), "");
        let pcr = PcrRow {
            custom_amplicon_set: "amp".to_owned(),
            pf_selected_pairs: 7,
            ..PcrRow::default()
        };
        let values = pcr.values();
        assert_eq!(values.len(), PCR_COLUMNS.len());
        assert_eq!(values[10].picard_text(), "7");
        // The two classes share the PanelMetricsBase tail.
        assert_eq!(&COLUMNS[20..], &PCR_COLUMNS[13..]);
    }

    #[test]
    fn targeted_pcr_defaults_differ_from_hs() {
        let pcr = Options::targeted_pcr();
        assert_eq!(pcr.variant, Variant::TargetedPcr);
        assert_eq!(
            (
                pcr.minimum_mapping_quality,
                pcr.minimum_base_quality,
                pcr.clip_overlapping_reads
            ),
            (1, 0, false)
        );
        assert_eq!(pcr.coverage_cap, DEFAULT_COVERAGE_CAP);
        let registered: Vec<&str> = [&TOOL, &PCR_TOOL]
            .iter()
            .flat_map(|t| t.options.iter().map(|o| o.names[0]))
            .collect();
        assert!(registered.contains(&"BAIT_INTERVALS"));
        assert!(registered.contains(&"AMPLICON_INTERVALS"));
        assert_eq!(registered.iter().filter(|n| **n == "STOP_AFTER").count(), 2);
    }

    #[test]
    fn missing_interval_files_name_the_option_as_picard_does() {
        let parse = |tool: &Tool, variant, args: &[&str]| {
            let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
            let inv = tool.parse(&args).expect("parses");
            Options::from_invocation(&inv, variant)
                .expect_err("an interval file is missing")
                .to_string()
        };
        assert_eq!(
            parse(&TOOL, Variant::Hs, &["I=x.bam", "BI=b.il"]),
            "Option 'TARGET_INTERVALS' must be specified at least 1 times."
        );
        assert_eq!(
            parse(&TOOL, Variant::Hs, &["I=x.bam", "TI=t.il"]),
            "Option 'BAIT_INTERVALS' must be specified at least 1 times."
        );
        assert_eq!(
            parse(&PCR_TOOL, Variant::TargetedPcr, &["I=x.bam", "AI=a.il"]),
            "Option 'TARGET_INTERVALS' must be specified at least 1 times."
        );
        assert_eq!(
            parse(&PCR_TOOL, Variant::TargetedPcr, &["I=x.bam", "TI=t.il"]),
            "Option 'AMPLICON_INTERVALS' is required."
        );
    }

    #[test]
    fn per_target_file_needs_all_reads_and_a_reference() {
        let parse = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(|s| (*s).to_owned()).collect();
            let inv = TOOL.parse(&args).expect("parses");
            Options::from_invocation(&inv, Variant::Hs)
        };
        let err = parse(&["I=x.bam", "BI=b.il", "TI=t.il", "PER_TARGET_COVERAGE=t.txt"])
            .expect_err("no reference");
        assert!(err.to_string().contains("REFERENCE_SEQUENCE"));
        let err = parse(&[
            "I=x.bam",
            "BI=b.il",
            "TI=t.il",
            "PER_TARGET_COVERAGE=t.txt",
            "R=r.fa",
            "LEVEL=SAMPLE",
        ])
        .expect_err("levels");
        assert!(err.to_string().contains("ALL_READS"));
        let opts = parse(&[
            "I=x.bam",
            "BI=b.il",
            "TI=t.il",
            "PER_TARGET_COVERAGE=t.txt",
            "R=r.fa",
            "ALLELE_FRACTION=null",
            "ALLELE_FRACTION=0.25",
            "NEAR_DISTANCE=-5",
        ])
        .expect("valid");
        assert_eq!(opts.allele_fractions, [0.25]);
        assert_eq!(opts.near_distance, -5);
        let err = TOOL
            .parse(&["I=x.bam".to_owned(), "STOP_AFTER=5".to_owned()])
            .expect_err("rejected");
        assert!(err.to_string().contains("Unrecognized option: STOP_AFTER"));
    }
}

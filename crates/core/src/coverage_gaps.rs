//! `CollectCoverageGaps` (docs/06 section 1): the target bases under a
//! depth, as maximal runs. Not a Picard tool.
//!
//! The depth of a target base is the high-quality depth `CollectHsMetrics`
//! writes to `PER_BASE_COVERAGE`: the tool runs the Hs accumulator
//! ([`CaptureMetrics`]) with the targets as both baits and targets and reads
//! its `ALL_READS` depth arrays ([`CaptureMetrics::into_target_depths`]), so
//! `PCT_REPORTABLE` at `MIN_DEPTH=20` equals `PCT_TARGET_BASES_20X` of
//! `CollectHsMetrics` on the same input and filter options. `COUNT_UNPAIRED=false`
//! is a pre-filter in front of that accumulator, the `CollectWgsMetrics`
//! meaning: unpaired reads and pairs with the mate unmapped add nothing.

use std::fmt::Write as _;
use std::path::PathBuf;

use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::hs_metrics::{self, CaptureMetrics, TargetDepths};
use crate::levels::Level;
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::{KeyClass, format_double};
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

pub const DEFAULT_MIN_DEPTH: i32 = 20;

/// The metrics class string. Not a Picard class: the last name is what
/// readers key on.
pub const CLASS: &str = "bammetrics.qc.CoverageGapMetrics";

/// The columns of the gaps BED, in its header line.
pub const BED_COLUMNS: [&str; 6] = ["chrom", "start", "end", "target", "mean_depth", "min_depth"];

/// The one row of the metrics file.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub target_territory: i64,
    pub target_count: i64,
    pub min_depth: i64,
    pub bases_at_or_above_min_depth: i64,
    pub pct_reportable: f64,
    pub gap_bases: i64,
    pub gap_count: i64,
    pub largest_gap: i64,
    pub targets_with_gap: i64,
    pub pct_targets_complete: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &[
        "TARGET_TERRITORY",
        "TARGET_COUNT",
        "MIN_DEPTH",
        "BASES_AT_OR_ABOVE_MIN_DEPTH",
        "PCT_REPORTABLE",
        "GAP_BASES",
        "GAP_COUNT",
        "LARGEST_GAP",
        "TARGETS_WITH_GAP",
        "PCT_TARGETS_COMPLETE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.target_territory.into(),
            self.target_count.into(),
            self.min_depth.into(),
            self.bases_at_or_above_min_depth.into(),
            self.pct_reportable.into(),
            self.gap_bases.into(),
            self.gap_count.into(),
            self.largest_gap.into(),
            self.targets_with_gap.into(),
            self.pct_targets_complete.into(),
        ]
    }
}

/// One maximal run of target bases below `MIN_DEPTH`.
#[derive(Clone, Debug, PartialEq)]
pub struct Gap {
    /// The merged target the run lies in.
    pub target: usize,
    /// 1-based, closed, as the targets are.
    pub start: i64,
    pub end: i64,
    pub mean_depth: f64,
    pub min_depth: i32,
}

impl Gap {
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// The runs of `depth` (one target's bases from `start`) below `min_depth`,
/// in position order.
pub fn gap_runs(target: usize, start: i64, depth: &[i32], min_depth: i32) -> Vec<Gap> {
    let mut gaps = Vec::new();
    let mut run: Option<(usize, i64, i32)> = None;
    let close = |gaps: &mut Vec<Gap>, from: usize, to: usize, sum: i64, min: i32| {
        gaps.push(Gap {
            target,
            start: start + from as i64,
            end: start + to as i64 - 1,
            mean_depth: sum as f64 / (to - from) as f64,
            min_depth: min,
        });
    };
    for (i, &d) in depth.iter().enumerate() {
        match (&mut run, d < min_depth) {
            (Some((_, sum, min)), true) => {
                *sum += i64::from(d);
                *min = (*min).min(d);
            }
            (None, true) => run = Some((i, i64::from(d), d)),
            (Some((from, sum, min)), false) => {
                close(&mut gaps, *from, i, *sum, *min);
                run = None;
            }
            (None, false) => {}
        }
    }
    if let Some((from, sum, min)) = run {
        close(&mut gaps, from, depth.len(), sum, min);
    }
    gaps
}

/// The metrics and the gaps of one run.
pub struct Report {
    pub row: Row,
    pub gaps: Vec<Gap>,
    pub depths: TargetDepths,
}

impl Report {
    fn new(depths: TargetDepths, min_depth: i32) -> Self {
        let mut gaps = Vec::new();
        let mut above = 0i64;
        let mut territory = 0i64;
        for (index, target) in depths.targets().iter().enumerate() {
            let depth = depths.depth(index);
            territory += depth.len() as i64;
            above += depth.iter().filter(|&&d| d >= min_depth).count() as i64;
            gaps.extend(gap_runs(index, target.start, depth, min_depth));
        }
        let target_count = depths.targets().len() as i64;
        let mut with_gap: Vec<usize> = gaps.iter().map(|g| g.target).collect();
        with_gap.dedup();
        let targets_with_gap = with_gap.len() as i64;
        let row = Row {
            target_territory: territory,
            target_count,
            min_depth: i64::from(min_depth),
            bases_at_or_above_min_depth: above,
            pct_reportable: above as f64 / territory as f64,
            gap_bases: territory - above,
            gap_count: gaps.len() as i64,
            largest_gap: gaps.iter().map(Gap::length).max().unwrap_or(0),
            targets_with_gap,
            pct_targets_complete: (target_count - targets_with_gap) as f64 / target_count as f64,
        };
        Self { row, gaps, depths }
    }

    /// The `## HISTOGRAM` section: gap runs by length.
    pub fn histogram(&self) -> HistogramSection {
        let mut lengths = Histogram::new();
        for gap in &self.gaps {
            lengths.increment(gap.length());
        }
        HistogramSection {
            key_label: "gap_length",
            key_class: KeyClass::Integer,
            series: vec![HistogramSeries {
                label: "gap_count".to_owned(),
                histogram: lengths,
            }],
        }
    }

    /// The metrics file.
    pub fn metrics_output(&self) -> MetricsOutput {
        MetricsOutput::from_rows(std::slice::from_ref(&self.row), Some(self.histogram()))
    }

    /// The `GAPS_FILE` BED (docs/06 section 1.4): a header line, then one
    /// line per gap in target order, 0-based half-open.
    pub fn gaps_output(&self) -> MetricsOutput {
        let mut text = format!("#{}\n", BED_COLUMNS.join("\t"));
        let mut rows = Vec::with_capacity(self.gaps.len());
        for gap in &self.gaps {
            let target = &self.depths.targets()[gap.target];
            let chrom = self.depths.contig_name(target.contig);
            let name = self.depths.target_name(gap.target);
            let _ = writeln!(
                text,
                "{chrom}\t{}\t{}\t{name}\t{}\t{}",
                gap.start - 1,
                gap.end,
                format_double(gap.mean_depth),
                gap.min_depth
            );
            rows.push(vec![
                Value::from(chrom),
                Value::Int(gap.start - 1),
                Value::Int(gap.end),
                Value::from(name),
                Value::Float(gap.mean_depth),
                Value::Int(i64::from(gap.min_depth)),
            ]);
        }
        MetricsOutput::text(&BED_COLUMNS, rows, text)
    }
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub target_intervals: Vec<PathBuf>,
    pub gaps_file: Option<PathBuf>,
    pub min_depth: i32,
    pub count_unpaired: bool,
    /// The Hs filter options: `MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`,
    /// `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`.
    pub hs: hs_metrics::Options,
}

impl Options {
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        let targets: Vec<PathBuf> = inv
            .all("TARGET_INTERVALS")
            .iter()
            .map(PathBuf::from)
            .collect();
        if targets.is_empty() {
            return Err(Error::Usage(
                "Option 'TARGET_INTERVALS' must be specified at least 1 times.".to_owned(),
            ));
        }
        let defaults = hs_metrics::Options::default();
        let hs = hs_metrics::Options {
            bait_intervals: targets.clone(),
            target_intervals: targets.clone(),
            levels: vec![Level::AllReads],
            minimum_mapping_quality: inv
                .parse("MINIMUM_MAPPING_QUALITY", defaults.minimum_mapping_quality)?,
            minimum_base_quality: inv
                .parse("MINIMUM_BASE_QUALITY", defaults.minimum_base_quality)?,
            clip_overlapping_reads: inv
                .flag("CLIP_OVERLAPPING_READS", defaults.clip_overlapping_reads)?,
            include_indels: inv.flag("INCLUDE_INDELS", false)?,
            ..defaults
        };
        let min_depth: i32 = inv.parse("MIN_DEPTH", DEFAULT_MIN_DEPTH)?;
        if min_depth < 1 {
            return Err(Error::Usage(format!(
                "MIN_DEPTH={min_depth}: a depth of 1 or more"
            )));
        }
        Ok(Self {
            target_intervals: targets,
            gaps_file: inv.get("GAPS_FILE").map(PathBuf::from),
            min_depth,
            count_unpaired: inv.flag("COUNT_UNPAIRED", true)?,
            hs,
        })
    }
}

/// The Hs accumulator behind an optional pair filter.
pub struct CoverageGaps {
    capture: CaptureMetrics,
    min_depth: i32,
    count_unpaired: bool,
}

impl CoverageGaps {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        Ok(Self {
            capture: CaptureMetrics::new(&opts.hs, header)?,
            min_depth: opts.min_depth,
            count_unpaired: opts.count_unpaired,
        })
    }
}

impl Accumulator for CoverageGaps {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if !self.count_unpaired && (!flags.is_segmented() || flags.is_mate_unmapped()) {
            return Ok(());
        }
        self.capture.accept(record)
    }

    fn finish(self) -> Result<Report> {
        Ok(Report::new(
            self.capture.into_target_depths()?,
            self.min_depth,
        ))
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectCoverageGaps",
    summary: "one row of the 10 CoverageGapMetrics columns: the target bases under MIN_DEPTH as \
              maximal runs, with the depth CollectHsMetrics writes per base; the metrics file \
              adds a histogram of gap lengths, and GAPS_FILE a BED of the gaps. Not a Picard \
              tool (docs/06).",
    options: &[
        Opt::path(
            &["TARGET_INTERVALS", "TI"],
            "Target intervals: an interval list, or a BED file; repeat for several.",
        )
        .repeated()
        .required(),
        Opt::path(
            &["GAPS_FILE"],
            "Also write a BED of the gaps: chrom, start, end, target, mean_depth, min_depth.",
        )
        .written(),
        Opt::int(
            &["MIN_DEPTH"],
            "A target base at or above this depth is reportable (default 20).",
        ),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Reads below this MAPQ add no coverage (default 20).",
        ),
        Opt::int(
            &["MINIMUM_BASE_QUALITY", "Q"],
            "Bases below this quality add no coverage (default 20).",
        ),
        Opt::bool(
            &["CLIP_OVERLAPPING_READS"],
            "Count the overlap of two mates once (default true).",
        ),
        Opt::bool(
            &["INCLUDE_INDELS"],
            "Count deleted bases as covered (default false).",
        ),
        Opt::bool(
            &["COUNT_UNPAIRED"],
            "false drops unpaired reads and pairs with the mate unmapped (default true, as CollectHsMetrics counts them).",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    let gaps = CoverageGaps::new(&opts, header)?;
    Ok(pass::boxed(gaps, move |report| {
        let mut output = ToolOutput::file(report.metrics_output());
        if let Some(path) = &opts.gaps_file {
            output.extra.push((path.clone(), report.gaps_output()));
        }
        Ok(output)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spans(gaps: &[Gap]) -> Vec<(i64, i64, i32, f64)> {
        gaps.iter()
            .map(|g| (g.start, g.end, g.min_depth, g.mean_depth))
            .collect()
    }

    #[test]
    fn gap_runs_are_maximal_and_stay_inside_the_target() {
        assert_eq!(gap_runs(0, 101, &[20, 25, 20], 20), []);
        assert_eq!(gap_runs(0, 101, &[], 20), []);
        // A whole-target gap.
        assert_eq!(
            spans(&gap_runs(0, 101, &[0, 1, 2], 20)),
            [(101, 103, 0, 1.0)]
        );
        // Gaps at both ends, one base at the threshold between them.
        assert_eq!(
            spans(&gap_runs(3, 1001, &[5, 19, 20, 7, 3], 20)),
            [(1001, 1002, 5, 12.0), (1004, 1005, 3, 5.0)]
        );
        // The threshold itself is reportable; one below is not.
        assert_eq!(
            spans(&gap_runs(0, 1, &[20, 19, 20], 20)),
            [(2, 2, 19, 19.0)]
        );
        assert_eq!(gap_runs(0, 1, &[20, 19, 20], 20)[0].target, 0);
        assert_eq!(gap_runs(0, 1, &[20, 19, 20], 20)[0].length(), 1);
    }

    #[test]
    fn row_columns_match_values() {
        let row = Row {
            target_territory: 10,
            target_count: 2,
            min_depth: 20,
            bases_at_or_above_min_depth: 8,
            pct_reportable: 0.8,
            gap_bases: 2,
            gap_count: 1,
            largest_gap: 2,
            targets_with_gap: 1,
            pct_targets_complete: 0.5,
        };
        assert_eq!(row.values().len(), Row::COLUMNS.len());
        assert_eq!(row.values()[4], Value::Float(0.8));
    }

    #[test]
    fn options_need_targets_and_a_positive_depth() {
        let strings =
            |args: &[&str]| -> Vec<String> { args.iter().map(|s| (*s).to_owned()).collect() };
        let err = TOOL
            .parse(&strings(&["I=x.bam"]))
            .and_then(|inv| Options::from_invocation(&inv))
            .expect_err("no targets");
        assert!(err.to_string().contains("TARGET_INTERVALS"), "{err}");
        let err = TOOL
            .parse(&strings(&["I=x.bam", "TI=t.bed", "MIN_DEPTH=0"]))
            .and_then(|inv| Options::from_invocation(&inv))
            .expect_err("depth 0");
        assert!(err.to_string().contains("MIN_DEPTH=0"), "{err}");
        let opts = Options::from_invocation(
            &TOOL
                .parse(&strings(&[
                    "I=x.bam",
                    "TI=t.bed",
                    "MQ=1",
                    "COUNT_UNPAIRED=false",
                    "GAPS_FILE=g.bed",
                ]))
                .expect("parses"),
        )
        .expect("options");
        assert_eq!(opts.hs.minimum_mapping_quality, 1);
        assert_eq!(opts.hs.bait_intervals, opts.target_intervals);
        assert!(!opts.count_unpaired);
        assert_eq!(
            opts.gaps_file.as_deref(),
            Some(std::path::Path::new("g.bed"))
        );
        assert_eq!(opts.min_depth, DEFAULT_MIN_DEPTH);
    }
}

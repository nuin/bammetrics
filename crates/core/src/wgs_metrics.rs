//! `CollectWgsMetrics`, `CollectRawWgsMetrics` and
//! `CollectWgsMetricsWithNonZeroCoverage`, mirroring Picard 3.5.0
//! (`CollectWgsMetrics.doWork`, `WgsMetricsCollector.addInfo`,
//! `FastWgsMetricsCollector`, `WgsMetrics.calculateDerivedFields`).
//!
//! The pileup comes from [`crate::locus`]. At every emitted locus whose
//! reference base is not `N`, the default collector walks the bases in
//! accumulation order: a quality of 2 or less is excluded as base quality;
//! the next `COVERAGE_CAP` bases above that feed the unfiltered depth and
//! the base-quality histogram; a base under `MINIMUM_BASE_QUALITY` or a
//! no-call is excluded as base quality; a second base from a read name
//! already seen at the locus is excluded as overlap; the rest is the
//! high-quality depth, capped at `COVERAGE_CAP` with the excess counted as
//! capped. `USE_FAST_ALGORITHM=true` swaps in `FastWgsMetricsCollector`,
//! which does the same arithmetic per alignment block when the block's
//! BEGIN event is reached, keeps per-position counters in a sliding window
//! sized by `READ_LENGTH`, and tests overlap against the other blocks of
//! the same name still open.
//!
//! The three tools share the columns: `CollectRawWgsMetrics` changes the
//! defaults and the class name, `CollectWgsMetricsWithNonZeroCoverage`
//! prints a second row computed after the zero-depth bin is cleared, with
//! the `CATEGORY` column first and the two depth histograms side by side.

use std::collections::HashMap;
use std::path::PathBuf;

use noodles::sam;

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::intervals::{self, uniqued};
use crate::locus::{
    BaseEntry, Classic, Edge, EdgeEntry, Engine, Filters, Finished, NameSet, Settings, Sink,
};
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput, ToolOutput, Value};
use crate::metrics_file::KeyClass;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sensitivity::{het_snp_sensitivity, normalize_histogram, phred_from_error_probability};
use crate::sequence::is_no_call;

/// `WgsMetrics.LOG_ODDS_THRESHOLD`.
const LOG_ODDS_THRESHOLD: f64 = 3.0;
/// `unfilteredBaseQHistogramArray` has `Byte.MAX_VALUE` slots.
const QUALITY_SLOTS: usize = 127;
/// `FastWgsMetricsCollector.ARRAY_SIZE_PER_READ_LENGTH`.
const ARRAY_SIZE_PER_READ_LENGTH: i64 = 2000;
/// The depths of the `PCT_nX` columns.
const PCT_DEPTHS: [usize; 14] = [1, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100];

/// `WgsMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 32] = [
    "GENOME_TERRITORY",
    "MEAN_COVERAGE",
    "SD_COVERAGE",
    "MEDIAN_COVERAGE",
    "MAD_COVERAGE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_DUPE",
    "PCT_EXC_UNPAIRED",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_CAPPED",
    "PCT_EXC_TOTAL",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_15X",
    "PCT_20X",
    "PCT_25X",
    "PCT_30X",
    "PCT_40X",
    "PCT_50X",
    "PCT_60X",
    "PCT_70X",
    "PCT_80X",
    "PCT_90X",
    "PCT_100X",
    "FOLD_80_BASE_PENALTY",
    "FOLD_90_BASE_PENALTY",
    "FOLD_95_BASE_PENALTY",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
];

/// `WgsMetricsWithNonZeroCoverage`: `CATEGORY` then the `WgsMetrics` columns.
pub const NON_ZERO_COLUMNS: [&str; 33] = [
    "CATEGORY",
    "GENOME_TERRITORY",
    "MEAN_COVERAGE",
    "SD_COVERAGE",
    "MEDIAN_COVERAGE",
    "MAD_COVERAGE",
    "PCT_EXC_ADAPTER",
    "PCT_EXC_MAPQ",
    "PCT_EXC_DUPE",
    "PCT_EXC_UNPAIRED",
    "PCT_EXC_BASEQ",
    "PCT_EXC_OVERLAP",
    "PCT_EXC_CAPPED",
    "PCT_EXC_TOTAL",
    "PCT_1X",
    "PCT_5X",
    "PCT_10X",
    "PCT_15X",
    "PCT_20X",
    "PCT_25X",
    "PCT_30X",
    "PCT_40X",
    "PCT_50X",
    "PCT_60X",
    "PCT_70X",
    "PCT_80X",
    "PCT_90X",
    "PCT_100X",
    "FOLD_80_BASE_PENALTY",
    "FOLD_90_BASE_PENALTY",
    "FOLD_95_BASE_PENALTY",
    "HET_SNP_SENSITIVITY",
    "HET_SNP_Q",
];

/// Which of the three programs runs: they differ in defaults, class name
/// and row layout.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Variant {
    Wgs,
    Raw,
    NonZero,
}

impl Variant {
    fn class(self) -> &'static str {
        match self {
            Variant::Wgs => "picard.analysis.WgsMetrics",
            Variant::Raw => "picard.analysis.CollectRawWgsMetrics$RawWgsMetrics",
            Variant::NonZero => {
                "picard.analysis.CollectWgsMetricsWithNonZeroCoverage$WgsMetricsWithNonZeroCoverage"
            }
        }
    }

    fn name(self) -> &'static str {
        match self {
            Variant::Wgs => "CollectWgsMetrics",
            Variant::Raw => "CollectRawWgsMetrics",
            Variant::NonZero => "CollectWgsMetricsWithNonZeroCoverage",
        }
    }

    /// `MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`, `COVERAGE_CAP`
    /// and `LOCUS_ACCUMULATION_CAP` defaults (`CollectRawWgsMetrics`'s
    /// constructor overrides the first four).
    fn defaults(self) -> (i32, i32, i32, i32) {
        match self {
            Variant::Raw => (0, 3, 100_000, 200_000),
            Variant::Wgs | Variant::NonZero => (20, 20, 250, 100_000),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub variant: Variant,
    pub reference: PathBuf,
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub coverage_cap: i32,
    pub locus_accumulation_cap: i32,
    /// `STOP_AFTER`: non-N loci to process; -1 (the default) or 0 for all.
    pub stop_after: i64,
    pub include_bq_histogram: bool,
    pub count_unpaired: bool,
    pub sample_size: i32,
    pub intervals: Option<PathBuf>,
    pub use_fast_algorithm: bool,
    pub read_length: i32,
}

impl Options {
    fn from_invocation(inv: &Invocation, variant: Variant) -> Result<Self> {
        let reference = inv.get("REFERENCE_SEQUENCE").ok_or_else(|| {
            Error::Usage(format!(
                "{}: REFERENCE_SEQUENCE is required",
                variant.name()
            ))
        })?;
        let (mq, bq, cap, locus_cap) = variant.defaults();
        let coverage_cap: i32 = inv.parse("COVERAGE_CAP", cap)?;
        let mut locus_accumulation_cap: i32 = inv.parse("LOCUS_ACCUMULATION_CAP", locus_cap)?;
        // "it doesn't make sense for the locus accumulation cap to be lower
        // than the coverage cap"
        if locus_accumulation_cap < coverage_cap {
            locus_accumulation_cap = coverage_cap;
        }
        if coverage_cap <= 0 {
            return Err(Error::Usage("Coverage cap must be positive.".to_owned()));
        }
        Ok(Self {
            variant,
            reference: reference.into(),
            minimum_mapping_quality: inv.parse("MINIMUM_MAPPING_QUALITY", mq)?,
            minimum_base_quality: inv.parse("MINIMUM_BASE_QUALITY", bq)?,
            coverage_cap,
            locus_accumulation_cap,
            stop_after: inv.parse("STOP_AFTER", -1)?,
            include_bq_histogram: inv.flag("INCLUDE_BQ_HISTOGRAM", false)?,
            count_unpaired: inv.flag("COUNT_UNPAIRED", false)?,
            sample_size: inv.parse("SAMPLE_SIZE", 10_000)?,
            intervals: inv.get("INTERVALS").map(PathBuf::from),
            use_fast_algorithm: inv.flag("USE_FAST_ALGORITHM", false)?,
            read_length: inv.parse("READ_LENGTH", 150)?,
        })
    }
}

/// The options of the three tools, the same rows in each; `$extra` adds
/// rows. The defaults of the raw variant differ (`Variant::defaults`).
macro_rules! wgs_options {
    ($($extra:expr),*) => {
        &[
            Opt::path(&["REFERENCE_SEQUENCE", "R"], "Reference FASTA with a .fai index.")
                .required(),
            Opt::int(
                &["MINIMUM_MAPPING_QUALITY", "MQ"],
                "Reads below this MAPQ add no coverage (default 20; raw 0).",
            ),
            Opt::int(
                &["MINIMUM_BASE_QUALITY", "Q"],
                "Bases below this quality add no coverage (default 20; raw 3).",
            ),
            Opt::int(
                &["COVERAGE_CAP", "CAP"],
                "Depths above this count as the cap (default 250; raw 100000).",
            ),
            Opt::int(
                &["LOCUS_ACCUMULATION_CAP"],
                "Reads processed per locus at most (default 100000; raw 200000).",
            ),
            Opt::int(&["STOP_AFTER"], "Stop after this many loci; -1 reads all (default -1)."),
            Opt::bool(
                &["INCLUDE_BQ_HISTOGRAM"],
                "Add the base quality histogram (default false).",
            ),
            Opt::bool(&["COUNT_UNPAIRED"], "Count unpaired reads (default false)."),
            Opt::int(
                &["SAMPLE_SIZE"],
                "Samples of the het SNP sensitivity simulation (default 10000).",
            ),
            Opt::path(&["INTERVALS"], "An interval list: only these loci count."),
            Opt::bool(
                &["USE_FAST_ALGORITHM"],
                "Picard's fast algorithm; the metrics are the same (default false).",
            ),
            Opt::int(
                &["READ_LENGTH"],
                "Window of the fast algorithm (default 150).",
            ),
            Opt::rejected(&["ASSUME_SORTED", "AS"]),
            $($extra),*
        ]
    };
}

const OPTIONS: &[Opt] = wgs_options!();

const NON_ZERO_OPTIONS: &[Opt] = wgs_options!(Opt::ignored(
    &["CHART_OUTPUT", "CHART"],
    "no chart is drawn"
));

pub const TOOL: Tool = Tool {
    name: "CollectWgsMetrics",
    summary: "one row of the 32 WgsMetrics columns; the metrics file adds the coverage \
              histogram.",
    options: OPTIONS,
    run: run_wgs,
    collector: None,
};

pub const RAW_TOOL: Tool = Tool {
    name: "CollectRawWgsMetrics",
    summary: "CollectWgsMetrics with MINIMUM_MAPPING_QUALITY 0, MINIMUM_BASE_QUALITY 3, \
              COVERAGE_CAP 100000 and LOCUS_ACCUMULATION_CAP 200000 as defaults, and the \
              RawWgsMetrics class name.",
    options: OPTIONS,
    run: run_raw,
    collector: None,
};

pub const NON_ZERO_TOOL: Tool = Tool {
    name: "CollectWgsMetricsWithNonZeroCoverage",
    summary: "two rows, CATEGORY WHOLE_GENOME and NON_ZERO_REGIONS (the zero-depth loci left \
              out), with the WgsMetrics columns.",
    options: NON_ZERO_OPTIONS,
    run: run_non_zero,
    collector: None,
};

fn run_wgs(inv: &Invocation) -> Result<ToolOutput> {
    run_variant(inv, Variant::Wgs)
}

fn run_raw(inv: &Invocation) -> Result<ToolOutput> {
    run_variant(inv, Variant::Raw)
}

fn run_non_zero(inv: &Invocation) -> Result<ToolOutput> {
    run_variant(inv, Variant::NonZero)
}

/// The per-base totals of `AbstractWgsMetricsCollector`.
#[derive(Clone, Debug)]
pub struct Counts {
    /// Loci per unfiltered depth, 0 to the cap.
    pub unfiltered_depth: Vec<i64>,
    /// Bases per quality, 0 to 126, counted while under the cap.
    pub unfiltered_baseq: Vec<i64>,
    /// Loci per high-quality depth, 0 to the cap.
    pub high_quality_depth: Vec<i64>,
    pub excluded_baseq: i64,
    pub excluded_overlap: i64,
    pub excluded_capped: i64,
    coverage_cap: usize,
    minimum_base_quality: i32,
}

impl Counts {
    fn new(opts: &Options) -> Self {
        let cap = usize::try_from(opts.coverage_cap).unwrap_or(0);
        Self {
            unfiltered_depth: vec![0; cap + 1],
            unfiltered_baseq: vec![0; QUALITY_SLOTS],
            high_quality_depth: vec![0; cap + 1],
            excluded_baseq: 0,
            excluded_overlap: 0,
            excluded_capped: 0,
            coverage_cap: cap,
            minimum_base_quality: opts.minimum_base_quality,
        }
    }

    /// `unfilteredBaseQHistogramArray[quality]++`: an index error above
    /// 126, as the Java array has `Byte.MAX_VALUE` slots.
    fn count_quality(&mut self, quality: i8) -> Result<()> {
        let slot = self
            .unfiltered_baseq
            .get_mut(usize::try_from(quality).unwrap_or(usize::MAX))
            .ok_or_else(|| Error::Metric(format!("base quality {quality} is above 126")))?;
        *slot += 1;
        Ok(())
    }

    /// The end of `addInfo`: the locus's depths into the histograms.
    fn count_locus(&mut self, pileup: usize, unfiltered: usize) {
        let high_quality = pileup.min(self.coverage_cap);
        if high_quality < pileup {
            self.excluded_capped += (pileup - self.coverage_cap) as i64;
        }
        self.high_quality_depth[high_quality] += 1;
        self.unfiltered_depth[unfiltered] += 1;
    }
}

/// `CollectWgsMetrics.WgsMetricsCollector`.
struct ClassicSink {
    counts: Counts,
    names: NameSet,
}

impl Sink<BaseEntry> for ClassicSink {
    fn locus(
        &mut self,
        _contig: usize,
        _position: i64,
        entries: &[BaseEntry],
        _reference: &[u8],
        reference_n: bool,
        _counter: i64,
    ) -> Result<()> {
        if reference_n {
            return Ok(());
        }
        self.names.reset(entries.len());
        let (mut pileup, mut unfiltered) = (0usize, 0usize);
        for (index, entry) in entries.iter().enumerate() {
            let quality = entry.read.quality(entry.offset)?;
            if quality <= 2 {
                self.counts.excluded_baseq += 1;
                continue;
            }
            if unfiltered < self.counts.coverage_cap {
                self.counts.count_quality(quality)?;
                unfiltered += 1;
            }
            if i32::from(quality) < self.counts.minimum_base_quality
                || is_no_call(entry.read.base(entry.offset)?)
            {
                self.counts.excluded_baseq += 1;
                continue;
            }
            let name_of = |i: u32| entries[i as usize].read.name.as_ref();
            if !self
                .names
                .insert(index as u32, entry.read.name_hash, name_of)
            {
                self.counts.excluded_overlap += 1;
                continue;
            }
            pileup += 1;
        }
        self.counts.count_locus(pileup, unfiltered);
        Ok(())
    }
}

/// `picard.analysis.CounterManager`: two per-position counters over a
/// window of `READ_LENGTH * 2000` positions that slides forward as loci
/// come in. A position outside the window is an index error, which Picard
/// blames on a `READ_LENGTH` far below the actual read length.
#[derive(Debug)]
struct CounterManager {
    array_length: i64,
    read_length: i64,
    offset: i64,
    pileup: Vec<i32>,
    unfiltered: Vec<i32>,
}

impl CounterManager {
    fn new(read_length: i32) -> Result<Self> {
        let array_length = i64::from(read_length) * ARRAY_SIZE_PER_READ_LENGTH;
        let size = usize::try_from(array_length)
            .ok()
            .filter(|&n| n <= i32::MAX as usize)
            .ok_or_else(|| {
                Error::Usage(format!("READ_LENGTH={read_length}: negative array size"))
            })?;
        Ok(Self {
            array_length,
            read_length: i64::from(read_length),
            offset: 0,
            pileup: vec![0; size],
            unfiltered: vec![0; size],
        })
    }

    fn check_out_of_bounds(&mut self, locus: i64) {
        if locus - self.offset + self.read_length >= self.array_length {
            if locus - self.offset < self.array_length {
                let skip = usize::try_from(locus - self.offset).unwrap_or(0);
                for array in [&mut self.pileup, &mut self.unfiltered] {
                    array.copy_within(skip.., 0);
                    let keep = array.len() - skip;
                    array[keep..].fill(0);
                }
            } else {
                self.clear();
            }
            self.offset = locus;
        }
    }

    fn clear(&mut self) {
        self.pileup.fill(0);
        self.unfiltered.fill(0);
        self.offset = 0;
    }

    fn index(&self, position: i64) -> Result<usize> {
        usize::try_from(position - self.offset)
            .ok()
            .filter(|&i| i < self.pileup.len())
            .ok_or_else(|| {
                Error::Metric(format!(
                    "The requested index {position} is out of counter bounds. Possible cause of \
                     exception can be wrong READ_LENGTH parameter (much smaller than actual read length)"
                ))
            })
    }
}

/// `FastWgsMetricsCollector`.
struct FastSink {
    counts: Counts,
    stop_after: i64,
    previous_contig: Option<usize>,
    counters: CounterManager,
    /// `readsNames`: the open BEGIN events of each read name.
    names: HashMap<Box<[u8]>, Vec<EdgeEntry>>,
}

impl FastSink {
    /// `excludeByQuality`: the open blocks of the name that do not reach
    /// `position` or hold a base below the quality threshold there.
    fn exclude_by_quality(&self, open: &[EdgeEntry], position: i64) -> Result<usize> {
        let mut excluded = 0;
        for block in open {
            if position - block.ref_pos >= block.length as i64
                || i32::from(block.quality_at(position)?) < self.counts.minimum_base_quality
            {
                excluded += 1;
            }
        }
        Ok(excluded)
    }

    /// `processRecord`: every base of the block, at its own position. The
    /// block joins the open blocks of its name once its bases are done, so
    /// it never tests overlap against itself.
    fn process_record(&mut self, reference: &[u8], block: &EdgeEntry, counter: i64) -> Result<()> {
        let mut open = self.names.remove(&block.read.name).unwrap_or_default();
        let outcome = self.process_bases(reference, block, &open, counter);
        open.push(block.clone());
        self.names.insert(block.read.name.clone(), open);
        outcome
    }

    fn process_bases(
        &mut self,
        reference: &[u8],
        block: &EdgeEntry,
        open: &[EdgeEntry],
        counter: i64,
    ) -> Result<()> {
        let mut processed = counter;
        for i in 0..block.length {
            let position = block.ref_pos + i as i64;
            let base = usize::try_from(position - 1)
                .ok()
                .and_then(|p| reference.get(p))
                .copied()
                .ok_or_else(|| {
                    Error::Metric(format!("position {position} is past the end of the contig"))
                })?;
            if is_no_call(base) {
                continue;
            }
            let quality = block.read.quality(i + block.offset)?;
            if quality <= 2 {
                self.counts.excluded_baseq += 1;
            } else {
                let index = self.counters.index(position)?;
                if (self.counters.unfiltered[index] as usize) < self.counts.coverage_cap {
                    self.counts.count_quality(quality)?;
                    self.counters.unfiltered[index] += 1;
                }
                if i32::from(quality) < self.counts.minimum_base_quality
                    || is_no_call(block.read.base(i + block.offset)?)
                {
                    self.counts.excluded_baseq += 1;
                } else if open.len() - self.exclude_by_quality(open, position)? > 0 {
                    self.counts.excluded_overlap += 1;
                } else {
                    self.counters.pileup[index] += 1;
                }
            }
            processed += 1;
            if self.stop_after > 0 && processed > self.stop_after - 1 {
                break;
            }
        }
        Ok(())
    }

    /// `removeRecordFromMap`.
    fn close_block(&mut self, block: &EdgeEntry) {
        let Some(open) = self.names.get_mut(&block.read.name) else {
            return;
        };
        if open.len() == 1 {
            self.names.remove(&block.read.name);
        } else {
            open.retain(|b| b.id != block.id);
        }
    }
}

impl Sink<EdgeEntry> for FastSink {
    fn locus(
        &mut self,
        contig: usize,
        position: i64,
        entries: &[EdgeEntry],
        reference: &[u8],
        reference_n: bool,
        counter: i64,
    ) -> Result<()> {
        // prepareCollector
        if self.previous_contig != Some(contig) {
            self.names.clear();
            self.counters.clear();
            self.previous_contig = Some(contig);
        }
        self.counters.check_out_of_bounds(position);
        for entry in entries {
            if entry.begin {
                self.process_record(reference, entry, counter)?;
            } else {
                self.close_block(entry);
            }
        }
        if !reference_n {
            let index = self.counters.index(position)?;
            let pileup = self.counters.pileup[index] as usize;
            let unfiltered = self.counters.unfiltered[index] as usize;
            self.counts.count_locus(pileup, unfiltered);
        }
        Ok(())
    }
}

/// `CollectWgsMetricsWithNonZeroCoverage` with `USE_FAST_ALGORITHM=true`.
/// Picard 3.5.0 keeps `WgsMetricsWithNonZeroCoverageCollector`, a
/// `WgsMetricsCollector`, whatever the flag says, but `doWork` still picks
/// the `EdgeReadIterator`, so the collector's `addInfo` casts each
/// `EdgingRecordAndOffset` to `SamLocusIterator.RecordAndOffset` and the
/// run dies with a `ClassCastException` at the first non-N locus holding
/// an event. A locus with no event counts as depth 0, and an N locus is
/// skipped before the cast, so a BAM with no record reaching a locus
/// completes with every depth at 0.
struct NonZeroFastSink {
    counts: Counts,
}

impl Sink<EdgeEntry> for NonZeroFastSink {
    fn locus(
        &mut self,
        _contig: usize,
        _position: i64,
        entries: &[EdgeEntry],
        _reference: &[u8],
        reference_n: bool,
        _counter: i64,
    ) -> Result<()> {
        if reference_n {
            return Ok(());
        }
        if !entries.is_empty() {
            return Err(Error::Metric(
                "CollectWgsMetricsWithNonZeroCoverage with USE_FAST_ALGORITHM=true: Picard 3.5.0 \
                 fails at the first covered locus (its collector casts the fast iterator's \
                 EdgingRecordAndOffset to SamLocusIterator.RecordAndOffset)"
                    .to_owned(),
            ));
        }
        self.counts.count_locus(0, 0);
        Ok(())
    }
}

/// `WgsMetrics` with every field derived.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub category: Option<&'static str>,
    pub genome_territory: i64,
    pub mean_coverage: f64,
    pub sd_coverage: f64,
    pub median_coverage: f64,
    pub mad_coverage: f64,
    pub pct_exc_adapter: f64,
    pub pct_exc_mapq: f64,
    pub pct_exc_dupe: f64,
    pub pct_exc_unpaired: f64,
    pub pct_exc_baseq: f64,
    pub pct_exc_overlap: f64,
    pub pct_exc_capped: f64,
    pub pct_exc_total: f64,
    pub pct_depth: [f64; 14],
    pub fold_80_base_penalty: f64,
    pub fold_90_base_penalty: f64,
    pub fold_95_base_penalty: f64,
    pub het_snp_sensitivity: f64,
    pub het_snp_q: f64,
}

impl Row {
    fn values(&self) -> Vec<Value> {
        let mut values: Vec<Value> = Vec::with_capacity(33);
        if let Some(category) = self.category {
            values.push(category.into());
        }
        values.extend([
            Value::Int(self.genome_territory),
            self.mean_coverage.into(),
            self.sd_coverage.into(),
            self.median_coverage.into(),
            self.mad_coverage.into(),
            self.pct_exc_adapter.into(),
            self.pct_exc_mapq.into(),
            self.pct_exc_dupe.into(),
            self.pct_exc_unpaired.into(),
            self.pct_exc_baseq.into(),
            self.pct_exc_overlap.into(),
            self.pct_exc_capped.into(),
            self.pct_exc_total.into(),
        ]);
        values.extend(self.pct_depth.iter().map(|&p| Value::Float(p)));
        values.extend([
            self.fold_80_base_penalty.into(),
            self.fold_90_base_penalty.into(),
            self.fold_95_base_penalty.into(),
            self.het_snp_sensitivity.into(),
            self.het_snp_q.into(),
        ]);
        values
    }
}

/// `AbstractWgsMetricsCollector.getHistogram`: a bin for every slot, the
/// empty ones included.
fn histogram(array: &[i64]) -> Histogram {
    let mut histogram = Histogram::new();
    for (i, &count) in array.iter().enumerate() {
        histogram.increment_by(i as i64, count as f64);
    }
    histogram
}

/// `CollectWgsMetrics.generateWgsMetrics` (the fractions) followed by
/// `WgsMetrics.calculateDerivedFields`.
fn derive(counts: &Counts, filters: &Filters, sample_size: i32) -> Result<Row> {
    let high_quality = histogram(&counts.high_quality_depth);
    let unfiltered = histogram(&counts.unfiltered_depth);
    let base_quality = histogram(&counts.unfiltered_baseq);

    let total = high_quality.sum();
    let total_with_excludes = total
        + filters.excluded_dupe as f64
        + filters.excluded_adapter as f64
        + filters.excluded_mapq as f64
        + filters.excluded_unpaired as f64
        + counts.excluded_baseq as f64
        + counts.excluded_overlap as f64
        + counts.excluded_capped as f64;
    let pct = |excluded: i64| {
        if total_with_excludes == 0.0 {
            0.0
        } else {
            excluded as f64 / total_with_excludes
        }
    };

    let genome_territory = high_quality.count() as i64;
    let mean_coverage = high_quality.mean();
    let mut pct_depth = [0.0; 14];
    for (pct, &depth) in pct_depth.iter_mut().zip(&PCT_DEPTHS) {
        let at_least: i64 = counts.high_quality_depth.iter().skip(depth).sum();
        *pct = at_least as f64 / genome_territory as f64;
    }
    let fold = |fraction: f64| {
        if high_quality.count() > 0.0 {
            mean_coverage / high_quality.percentile(fraction).unwrap_or(f64::NAN)
        } else {
            0.0
        }
    };
    if sample_size <= 0 {
        return Err(Error::Usage(
            "Sample size is required when a baseQ histogram is given when deriving metrics."
                .to_owned(),
        ));
    }
    let het_snp_sensitivity = het_snp_sensitivity(
        &normalize_histogram(&unfiltered),
        &normalize_histogram(&base_quality),
        usize::try_from(sample_size).unwrap_or(0),
        LOG_ODDS_THRESHOLD,
    )?;
    Ok(Row {
        category: None,
        genome_territory,
        mean_coverage,
        sd_coverage: high_quality.standard_deviation(),
        median_coverage: high_quality.median(),
        mad_coverage: high_quality.median_absolute_deviation(),
        pct_exc_adapter: pct(filters.excluded_adapter),
        pct_exc_mapq: pct(filters.excluded_mapq),
        pct_exc_dupe: pct(filters.excluded_dupe),
        pct_exc_unpaired: pct(filters.excluded_unpaired),
        pct_exc_baseq: pct(counts.excluded_baseq),
        pct_exc_overlap: pct(counts.excluded_overlap),
        pct_exc_capped: pct(counts.excluded_capped),
        pct_exc_total: if total_with_excludes == 0.0 {
            0.0
        } else {
            (total_with_excludes - total) / total_with_excludes
        },
        pct_depth,
        fold_80_base_penalty: fold(0.2),
        fold_90_base_penalty: fold(0.1),
        fold_95_base_penalty: fold(0.05),
        het_snp_sensitivity,
        het_snp_q: f64::from(phred_from_error_probability(1.0 - het_snp_sensitivity)),
    })
}

/// The file: one row (two for the non-zero variant), the depth histogram
/// (two of them for the non-zero variant) and, when asked for, the
/// base-quality histogram.
fn output(opts: &Options, mut counts: Counts, filters: &Filters) -> Result<MetricsOutput> {
    let series = |label: &str, histogram: Histogram| HistogramSeries {
        label: label.to_owned(),
        histogram,
    };
    let mut rows = Vec::new();
    let mut histograms = Vec::new();
    match opts.variant {
        Variant::Wgs | Variant::Raw => {
            rows.push(derive(&counts, filters, opts.sample_size)?.values());
            histograms.push(series(
                "high_quality_coverage_count",
                histogram(&counts.high_quality_depth),
            ));
        }
        Variant::NonZero => {
            histograms.push(series(
                "count_WHOLE_GENOME",
                histogram(&counts.high_quality_depth),
            ));
            histograms.push(series(
                "count_NON_ZERO_REGIONS",
                non_zero_histogram(&counts.high_quality_depth),
            ));
            let mut whole = derive(&counts, filters, opts.sample_size)?;
            whole.category = Some("WHOLE_GENOME");
            counts.high_quality_depth[0] = 0;
            counts.unfiltered_depth[0] = 0;
            let mut non_zero = derive(&counts, filters, opts.sample_size)?;
            non_zero.category = Some("NON_ZERO_REGIONS");
            rows.push(whole.values());
            rows.push(non_zero.values());
        }
    }
    if opts.include_bq_histogram {
        histograms.push(series(
            "unfiltered_baseq_count",
            histogram(&counts.unfiltered_baseq),
        ));
    }
    Ok(MetricsOutput {
        class: opts.variant.class(),
        columns: match opts.variant {
            Variant::NonZero => &NON_ZERO_COLUMNS,
            Variant::Wgs | Variant::Raw => &COLUMNS,
        },
        rows,
        histogram: Some(HistogramSection {
            key_label: "coverage",
            key_class: KeyClass::Integer,
            series: histograms,
        }),
        with_headers: true,
        text: None,
    })
}

/// `getDepthHistogramNonZero`: the bins from 1 up.
fn non_zero_histogram(array: &[i64]) -> Histogram {
    let mut histogram = Histogram::new();
    for (i, &count) in array.iter().enumerate().skip(1) {
        histogram.increment_by(i as i64, count as f64);
    }
    histogram
}

/// `CollectWgsMetrics.doWork`: the reference, the intervals, the locus
/// iterator with its filters, the collector, then the file.
fn run_variant(inv: &Invocation, variant: Variant) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv, variant)?;
    let settings = |locus_cap: Option<usize>| Settings {
        minimum_mapping_quality: opts.minimum_mapping_quality,
        count_unpaired: opts.count_unpaired,
        locus_accumulation_cap: locus_cap,
        stop_after: opts.stop_after,
    };
    let setup = |header: &sam::Header| -> Result<(Reference, Option<Vec<intervals::Interval>>)> {
        let reference = Reference::open(&opts.reference, header)?;
        let unique = match &opts.intervals {
            Some(path) => Some(uniqued(&intervals::read_files(&[path], header)?)),
            None => None,
        };
        Ok((reference, unique))
    };
    let (counts, filters) = if opts.use_fast_algorithm && variant == Variant::NonZero {
        let Finished { sink, filters } = inv.read_whole(&inv.input, |header| {
            let (reference, unique) = setup(header)?;
            let contigs = header.reference_sequences().len();
            Ok(Engine::new(
                Edge::new(unique.as_deref(), contigs),
                NonZeroFastSink {
                    counts: Counts::new(&opts),
                },
                reference,
                header,
                unique.as_deref(),
                settings(None),
            ))
        })?;
        (sink.counts, filters)
    } else if opts.use_fast_algorithm {
        let Finished { sink, filters } = inv.read_whole(&inv.input, |header| {
            let (reference, unique) = setup(header)?;
            let contigs = header.reference_sequences().len();
            Ok(Engine::new(
                Edge::new(unique.as_deref(), contigs),
                FastSink {
                    counts: Counts::new(&opts),
                    stop_after: opts.stop_after,
                    previous_contig: None,
                    counters: CounterManager::new(opts.read_length)?,
                    names: HashMap::new(),
                },
                reference,
                header,
                unique.as_deref(),
                settings(None),
            ))
        })?;
        (sink.counts, filters)
    } else {
        let Finished { sink, filters } = inv.read_whole(&inv.input, |header| {
            let (reference, unique) = setup(header)?;
            Ok(Engine::new(
                Classic,
                ClassicSink {
                    counts: Counts::new(&opts),
                    names: NameSet::default(),
                },
                reference,
                header,
                unique.as_deref(),
                settings(Some(
                    usize::try_from(opts.locus_accumulation_cap).unwrap_or(0),
                )),
            ))
        })?;
        (sink.counts, filters)
    };
    Ok(output(&opts, counts, &filters)?.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn filters() -> Filters {
        Filters::new(20, false)
    }

    fn options(variant: Variant) -> Options {
        Options {
            variant,
            reference: PathBuf::from("ref.fa"),
            minimum_mapping_quality: 20,
            minimum_base_quality: 20,
            coverage_cap: 4,
            locus_accumulation_cap: 4,
            stop_after: -1,
            include_bq_histogram: false,
            count_unpaired: false,
            sample_size: 10,
            intervals: None,
            use_fast_algorithm: false,
            read_length: 150,
        }
    }

    #[test]
    fn fractions_divide_by_bases_plus_every_exclusion() {
        let mut counts = Counts::new(&options(Variant::Wgs));
        // 2 loci at depth 3, 1 at depth 0: 6 high-quality bases.
        counts.high_quality_depth = vec![1, 0, 0, 2, 0];
        counts.unfiltered_depth = vec![1, 0, 0, 2, 0];
        counts.unfiltered_baseq[30] = 6;
        counts.excluded_baseq = 2;
        counts.excluded_overlap = 1;
        counts.excluded_capped = 1;
        let mut filters = filters();
        filters.excluded_dupe = 4;
        filters.excluded_mapq = 6;
        let row = derive(&counts, &filters, 10).expect("row");
        assert_eq!(row.genome_territory, 3);
        assert_eq!(row.mean_coverage, 2.0);
        assert_eq!(row.pct_exc_dupe, 4.0 / 20.0);
        assert_eq!(row.pct_exc_mapq, 6.0 / 20.0);
        assert_eq!(row.pct_exc_baseq, 2.0 / 20.0);
        assert_eq!(row.pct_exc_overlap, 1.0 / 20.0);
        assert_eq!(row.pct_exc_capped, 1.0 / 20.0);
        assert_eq!(row.pct_exc_total, 14.0 / 20.0);
        assert_eq!(row.pct_depth[0], 2.0 / 3.0);
        assert_eq!(row.pct_depth[1], 0.0);
        assert_eq!(row.median_coverage, 3.0);
        assert_eq!(row.fold_80_base_penalty, f64::INFINITY);
    }

    #[test]
    fn empty_histogram_prints_nan_mean_and_zero_folds() {
        let counts = Counts::new(&options(Variant::Wgs));
        let row = derive(&counts, &filters(), 10).expect("row");
        assert_eq!(row.genome_territory, 0);
        assert!(row.mean_coverage.is_nan());
        assert_eq!(row.fold_80_base_penalty, 0.0);
        assert_eq!(row.pct_exc_total, 0.0);
        // No locus at all: every normalised depth is NaN and so is the sum.
        assert!(row.het_snp_sensitivity.is_nan());
        assert_eq!(row.het_snp_q, 0.0);
        assert!(derive(&counts, &filters(), 0).is_err());
        // Loci without a base: depth 0 everywhere, NaN quality weights,
        // every draw rejected, sensitivity 0.
        let mut counts = Counts::new(&options(Variant::Wgs));
        counts.high_quality_depth[0] = 3;
        counts.unfiltered_depth[0] = 3;
        let row = derive(&counts, &filters(), 10).expect("row");
        assert_eq!(row.het_snp_sensitivity, 0.0);
        assert_eq!(row.mean_coverage, 0.0);
        assert!(row.fold_80_base_penalty.is_nan());
    }

    #[test]
    fn non_zero_variant_clears_the_zero_bin_for_its_second_row() {
        let opts = options(Variant::NonZero);
        let mut counts = Counts::new(&opts);
        counts.high_quality_depth = vec![5, 1, 1, 0, 0];
        counts.unfiltered_depth = vec![5, 1, 1, 0, 0];
        counts.unfiltered_baseq[30] = 3;
        let out = output(&opts, counts, &filters()).expect("output");
        assert_eq!(out.columns[0], "CATEGORY");
        assert_eq!(out.rows.len(), 2);
        assert_eq!(out.rows[0][0], Value::Str("WHOLE_GENOME".to_owned()));
        assert_eq!(out.rows[0][1], Value::Int(7));
        assert_eq!(out.rows[1][1], Value::Int(2));
        let section = out.histogram.expect("histograms");
        assert_eq!(section.series.len(), 2);
        assert_eq!(section.series[1].histogram.min(), Some(1));
        assert_eq!(section.series[1].histogram.get(1), Some(1.0));
    }

    #[test]
    fn non_zero_fast_sink_fails_at_the_first_covered_locus() {
        let mut sink = NonZeroFastSink {
            counts: Counts::new(&options(Variant::NonZero)),
        };
        let read = std::sync::Arc::new(crate::locus::Read {
            name: b"r1".as_slice().into(),
            name_hash: 0,
            bases: b"ACGT".as_slice().into(),
            qualities: [30; 4].as_slice().into(),
            start: 2,
            end: 5,
            blocks: Vec::new(),
        });
        let event = EdgeEntry {
            read,
            offset: 0,
            length: 4,
            ref_pos: 2,
            id: 1,
            begin: true,
        };
        let reference = b"NACGTACGT";
        sink.locus(0, 1, &[], reference, true, 0).expect("N locus");
        sink.locus(0, 2, &[], reference, false, 1)
            .expect("empty locus");
        let events = [event];
        sink.locus(0, 3, &events, reference, true, 1)
            .expect("an event on an N locus is skipped before the cast");
        let error = sink
            .locus(0, 4, &events, reference, false, 2)
            .expect_err("a covered locus");
        assert!(error.to_string().contains("USE_FAST_ALGORITHM=true"));
        assert_eq!(sink.counts.high_quality_depth[0], 1);
        assert_eq!(sink.counts.unfiltered_depth[0], 1);
    }

    #[test]
    fn counter_manager_slides_and_rejects_far_positions() {
        let mut counters = CounterManager::new(2).expect("size");
        assert_eq!(counters.pileup.len(), 4000);
        counters.check_out_of_bounds(1);
        let i = counters.index(1).expect("in window");
        counters.pileup[i] = 7;
        counters.check_out_of_bounds(3999);
        assert_eq!(counters.offset, 3999);
        assert!(counters.index(1).is_err());
        counters.check_out_of_bounds(100_000);
        assert_eq!(counters.offset, 100_000);
        assert!(CounterManager::new(-1).is_err());
    }
}

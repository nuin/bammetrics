//! `CollectHsMetrics`, mirroring Picard 3.5.0 `CollectTargetedMetrics.doWork`,
//! `TargetMetricsCollector.PerUnitTargetMetricCollector` (`acceptRecord`,
//! `finish`) and `HsMetricCollector.convertMetric`, at the `ALL_READS` level.
//!
//! One pass over the BAM in any order. Each record is placed against the
//! uniqued baits and targets with binary searches (`intervals::OverlapIndex`),
//! and coverage is kept only for target bases: one depth array per uniqued
//! target, laid end to end, at high quality and unfiltered.

use std::io::{self, Write};
use std::ops::Range;
use std::path::{Path, PathBuf};

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::{bam, sam};

use crate::adapter::Adapters;
use crate::duplication::{estimate_library_size, estimate_roi};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::intervals::{self, Interval, OverlapIndex};
use crate::metrics_file::{HistogramColumn, HistogramTable, MetricsFile, format_double};
use crate::pass::{self, Accumulator};
use crate::record::{position, record_name};
use crate::reference::Reference;
use crate::sensitivity::{het_snp_sensitivity, math_round, phred_from_error_probability};

/// `CollectHsMetrics` overrides of `CollectTargetedMetrics` defaults, and
/// the defaults it keeps.
pub const DEFAULT_MINIMUM_MAPPING_QUALITY: i32 = 20;
pub const DEFAULT_MINIMUM_BASE_QUALITY: i32 = 20;
pub const DEFAULT_COVERAGE_CAP: i32 = 200;
/// `TargetedPcrMetricsCollector.NEAR_PROBE_DISTANCE_DEFAULT`.
pub const DEFAULT_NEAR_DISTANCE: i32 = 250;
pub const DEFAULT_SAMPLE_SIZE: usize = 10_000;

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

#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    /// `BAIT_INTERVALS`: interval lists (or BED files) of the baits.
    pub bait_intervals: Vec<PathBuf>,
    /// `TARGET_INTERVALS`.
    pub target_intervals: Vec<PathBuf>,
    /// `BAIT_SET_NAME`; by default the bait file names up to their first dot.
    pub bait_set_name: Option<String>,
    /// `REFERENCE_SEQUENCE`: needed for `AT_DROPOUT` and `GC_DROPOUT`.
    pub reference: Option<PathBuf>,
    pub minimum_mapping_quality: i32,
    pub minimum_base_quality: i32,
    pub clip_overlapping_reads: bool,
    /// `INCLUDE_INDELS`: inserted bases count as on target, deleted bases as covered.
    pub include_indels: bool,
    pub coverage_cap: i32,
    pub near_distance: i32,
    /// `SAMPLE_SIZE` of the `HET_SNP_SENSITIVITY` simulation.
    pub sample_size: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            bait_intervals: Vec::new(),
            target_intervals: Vec::new(),
            bait_set_name: None,
            reference: None,
            minimum_mapping_quality: DEFAULT_MINIMUM_MAPPING_QUALITY,
            minimum_base_quality: DEFAULT_MINIMUM_BASE_QUALITY,
            clip_overlapping_reads: true,
            include_indels: false,
            coverage_cap: DEFAULT_COVERAGE_CAP,
            near_distance: DEFAULT_NEAR_DISTANCE,
            sample_size: DEFAULT_SAMPLE_SIZE,
        }
    }
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

/// One field of the row, typed as the Java field is.
#[derive(Clone, Debug, PartialEq)]
pub enum Field {
    Text(String),
    Long(i64),
    Double(f64),
    /// A null `Long` or `String`: an empty cell.
    Null,
}

impl Field {
    /// The cell as htsjdk `FormatUtil` prints it.
    pub fn picard_text(&self) -> String {
        match self {
            Field::Text(s) => s.clone(),
            Field::Long(v) => v.to_string(),
            Field::Double(v) => format_double(*v),
            Field::Null => String::new(),
        }
    }
}

/// The `HsMetrics` row. Java `long` fields are `i64`, `double` fields
/// `f64`, and the nullable `Long HS_LIBRARY_SIZE` an `Option`.
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
    pub hs_library_size: Option<i64>,
    /// `HS_PENALTY_10X` to `HS_PENALTY_100X`, in `PENALTY_GOALS` order.
    pub hs_penalty: [f64; 6],
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
}

impl Row {
    /// Every field in `COLUMNS` order. `SAMPLE`, `LIBRARY` and `READ_GROUP`
    /// are null at the `ALL_READS` level.
    pub fn fields(&self) -> Vec<Field> {
        use Field::{Double, Long, Null, Text};
        let mut fields = vec![
            Text(self.bait_set.clone()),
            Long(self.bait_territory),
            Double(self.bait_design_efficiency),
            Long(self.on_bait_bases),
            Long(self.near_bait_bases),
            Long(self.off_bait_bases),
            Double(self.pct_selected_bases),
            Double(self.pct_off_bait),
            Double(self.on_bait_vs_selected),
            Double(self.mean_bait_coverage),
            Double(self.pct_usable_bases_on_bait),
            Double(self.pct_usable_bases_on_target),
            Double(self.fold_enrichment),
            self.hs_library_size.map_or(Null, Long),
        ];
        fields.extend(self.hs_penalty.iter().map(|&p| Double(p)));
        fields.extend([
            Long(self.target_territory),
            Long(self.genome_size),
            Long(self.total_reads),
            Long(self.pf_reads),
            Long(self.pf_bases),
            Long(self.pf_unique_reads),
            Long(self.pf_uq_reads_aligned),
            Long(self.pf_bases_aligned),
            Long(self.pf_uq_bases_aligned),
            Long(self.on_target_bases),
            Double(self.pct_pf_reads),
            Double(self.pct_pf_uq_reads),
            Double(self.pct_pf_uq_reads_aligned),
            Double(self.mean_target_coverage),
            Double(self.median_target_coverage),
            Long(self.max_target_coverage),
            Long(self.min_target_coverage),
            Double(self.zero_cvg_targets_pct),
            Double(self.pct_exc_dupe),
            Double(self.pct_exc_adapter),
            Double(self.pct_exc_mapq),
            Double(self.pct_exc_baseq),
            Double(self.pct_exc_overlap),
            Double(self.pct_exc_off_target),
            Double(self.fold_80_base_penalty),
        ]);
        fields.extend(self.pct_target_bases.iter().map(|&p| Double(p)));
        fields.extend([
            Double(self.at_dropout),
            Double(self.gc_dropout),
            Double(self.het_snp_sensitivity),
            Double(self.het_snp_q),
            Null,
            Null,
            Null,
        ]);
        fields
    }
}

/// The row and the two histograms Picard writes after it.
#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub row: Row,
    /// `highQualityDepthHistogram`: target bases per high-quality depth.
    pub depth: Histogram,
    /// `unCappedUnfilteredBaseQHistogram`: on-target bases per quality,
    /// every base above Q2.
    pub base_quality: Histogram,
}

/// Write the report as a Picard `HsMetrics` file. The two histograms share
/// the `coverage_or_base_quality` key column, as in Picard.
pub fn write_metrics(
    report: &Report,
    headers: Vec<String>,
    out: &mut impl Write,
) -> io::Result<()> {
    MetricsFile {
        headers,
        class: "picard.analysis.directed.HsMetrics",
        columns: &COLUMNS,
        rows: vec![report.row.fields().iter().map(Field::picard_text).collect()],
        histogram: Some(HistogramTable {
            key_label: "coverage_or_base_quality",
            columns: vec![
                HistogramColumn {
                    label: "high_quality_coverage_count".to_owned(),
                    histogram: &report.depth,
                },
                HistogramColumn {
                    label: "unfiltered_baseq_count".to_owned(),
                    histogram: &report.base_quality,
                },
            ],
        }),
    }
    .write(out)
}

/// `CollectHsMetrics.getProbeSetName`: `BAIT_SET_NAME`, or the distinct
/// bait file names cut at their first dot, sorted and joined with dots.
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
#[derive(Default)]
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

/// Depths of every uniqued target, one array after another.
struct Coverage {
    /// Start of each target's depths in `high_quality` and `unfiltered`.
    offsets: Vec<usize>,
    high_quality: Vec<i32>,
    unfiltered: Vec<i32>,
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

/// Per-record state shared by every call to `accept`.
struct Collector {
    opts: Options,
    targets: OverlapIndex,
    probes: OverlapIndex,
    coverage: Coverage,
    hq_max_depth: i64,
    /// `baseQHistogramArray`: qualities of bases added while the unfiltered depth is within the cap.
    base_q_capped: [i64; QUALITY_SLOTS],
    /// `unCappedBaseQHistogramArray`.
    base_q: [i64; QUALITY_SLOTS],
    counts: Counts,
    adapters: Adapters,
    cigar: Cigar,
    bases: Vec<u8>,
}

impl Collector {
    /// `PerUnitTargetMetricCollector.acceptRecord`.
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

        let c = &mut self.counts;
        if !supplementary {
            c.total_reads += 1;
            if !flags.is_qc_fail() {
                c.pf_reads += 1;
                if !duplicate {
                    c.pf_unique_reads += 1;
                    if !unmapped {
                        c.pf_uq_reads_aligned += 1;
                    }
                }
            }
        }
        if flags.is_qc_fail() {
            return Ok(());
        }
        // Supplementary records add aligned bases but not read bases.
        if !supplementary {
            c.pf_bases += record.sequence().len() as i64;
        }
        if unmapped {
            return Ok(());
        }
        c.pf_bases_aligned += aligned;
        if !duplicate {
            c.pf_uq_bases_aligned += aligned;
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

        // HS_LIBRARY_SIZE counts first ends of mapped pairs near a bait.
        if !supplementary
            && paired
            && flags.is_first_segment()
            && !mate_unmapped
            && !probes.is_empty()
        {
            c.pf_selected_pairs += 1;
            if !duplicate {
                c.pf_selected_unique_pairs += 1;
            }
        }
        // Bait metrics come before the duplicate, adapter, MAPQ and overlap
        // filters, so they describe the assay rather than the usable data.
        if probes.is_empty() {
            c.off_probe_bases += aligned;
        } else {
            let on = bases_on_baits(start, &self.cigar, &self.probes.intervals()[probes]);
            c.on_probe_bases += on;
            c.near_probe_bases += aligned - on;
        }
        if duplicate {
            c.excluded_dupe += aligned;
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
                c.excluded_adapter += aligned;
                return Ok(());
            }
        }
        if i32::from(mapq) < self.opts.minimum_mapping_quality {
            c.excluded_mapq += aligned;
            return Ok(());
        }
        if self.opts.clip_overlapping_reads && paired && !mate_unmapped {
            let mate_start = match record.mate_alignment_start().transpose()? {
                Some(p) => p.get() as i64,
                None => 0,
            };
            let clip =
                overlapping_bases_to_clip(start, mate_start, flags.is_first_segment(), &self.cigar);
            c.excluded_overlap += clip;
            if clip > 0 {
                // `clipOverlappingAlignedBases`: a read starting with its mate
                // is unmapped, and so adds no coverage.
                if mate_start <= start {
                    return Ok(());
                }
                // Clip from the base after the last kept one, before any
                // trailing soft clip (a trailing H is not looked past).
                let mut clip_from = record.sequence().len() as i64 - clip + 1;
                if let Some(&(Kind::SoftClip, n)) = self.cigar.last() {
                    clip_from -= n as i64;
                }
                self.cigar = soft_clip_end(clip_from, &self.cigar)?;
            }
        }
        self.add_coverage(record, start, targets, mapped_in_pair)
    }

    /// The per-base loop of `acceptRecord` over the (clipped) CIGAR:
    /// summary counts for every aligned base, then target depths and
    /// quality histograms for on-target bases above Q2.
    fn add_coverage(
        &mut self,
        record: &bam::Record,
        start: i64,
        targets: Range<usize>,
        mapped_in_pair: bool,
    ) -> Result<()> {
        let malformed = |what| Error::MalformedRecord {
            name: record_name(record),
            what,
        };
        let qualities = record.quality_scores();
        let qualities = qualities.as_ref();
        // htsjdk reads a 0xFF quality string as no qualities.
        let qualities = if qualities.first() == Some(&0xFF) {
            &[][..]
        } else {
            qualities
        };
        let base_index = targets.start;
        let targets = &self.targets.intervals()[targets];
        let (min_bq, cap) = (self.opts.minimum_base_quality, self.opts.coverage_cap);
        let mut cursor = 0usize;
        let (mut read_offset, mut ref_pos) = (0usize, start);
        for &(kind, len) in &self.cigar {
            let counted = is_alignment(kind)
                || (self.opts.include_indels && matches!(kind, Kind::Insertion | Kind::Deletion));
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
                while cursor < targets.len() && targets[cursor].end < ref_pos {
                    cursor += 1;
                }
                let target =
                    (cursor < targets.len() && targets[cursor].start <= ref_pos).then_some(cursor);
                let high_quality = quality >= min_bq;
                let c = &mut self.counts;
                if !high_quality {
                    c.excluded_baseq += 1;
                } else if target.is_none() {
                    c.excluded_off_target += 1;
                } else {
                    c.on_target_bases += 1;
                    if mapped_in_pair {
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
                    let offset = self.coverage.offsets[base_index + t]
                        + usize::try_from(ref_pos - targets[t].start).unwrap_or(0);
                    if add_depth(&mut self.coverage.unfiltered[offset]) <= cap {
                        self.base_q_capped[slot] += 1;
                    }
                    self.base_q[slot] += 1;
                    if high_quality {
                        let depth = add_depth(&mut self.coverage.high_quality[offset]);
                        self.hq_max_depth = self.hq_max_depth.max(i64::from(depth));
                    }
                }
                read_offset += read_step;
                ref_pos += ref_step;
            }
        }
        Ok(())
    }
}

/// Inputs to `finish` that do not come from the reads.
struct Territory {
    bait_set: String,
    bait_territory: i64,
    target_territory: i64,
    genome_size: i64,
    /// `allTargets.getIntervals().size()`: target intervals before merging,
    /// the denominator of `ZERO_CVG_TARGETS_PCT`.
    raw_target_count: usize,
    /// GC fraction of each uniqued target, when a reference was given.
    target_gc: Option<Vec<f64>>,
}

/// Depth counts as Picard's histogram holds them: bins `0..prefill` exist
/// even when empty (`LongStream.range(0, n).forEach(i -> increment(i, 0))`),
/// other bins only when counted.
fn depth_histogram(counts: &[i64], prefill: usize) -> Histogram {
    let mut histogram = Histogram::new();
    for (depth, &count) in counts.iter().enumerate() {
        if depth < prefill || count > 0 {
            histogram.increment_by(depth as i64, count as f64);
        }
    }
    histogram
}

/// `TheoreticalSensitivity.normalizeHistogram`: an array as long as the
/// number of bins, slot i holding bin i over the sum of all bins.
fn normalize(histogram: &Histogram) -> Vec<f64> {
    let total = histogram.count();
    (0..histogram.len())
        .map(|i| histogram.get(i as i64).map_or(0.0, |v| v / total))
        .collect()
}

/// `SequenceUtil.calculateGc`.
fn gc_fraction(bases: &[u8]) -> f64 {
    let gc = bases
        .iter()
        .filter(|&&b| matches!(b, b'C' | b'G' | b'c' | b'g'))
        .count();
    gc as f64 / bases.len() as f64
}

/// `HsMetricCollector.calculateHsPenalty`: the fold of sequencing needed to
/// bring 80% of target bases to `goal`x, searched by halving steps on the
/// pair multiplier until `estimateRoi` is within 0.1% of the goal. -1 when
/// it does not converge in 10000 steps, 0 without a library size.
fn hs_penalty(library_size: Option<i64>, c: &Counts, row: &Row, goal: i32) -> f64 {
    let Some(size) = library_size else {
        return 0.0;
    };
    let mean_coverage = c.on_target_from_pair_bases as f64 / row.target_territory as f64;
    let fold80 = row.fold_80_base_penalty;
    let (pairs, unique_pairs) = (c.pf_selected_pairs, c.pf_selected_unique_pairs);
    let on_target_pct = c.on_target_bases as f64 / c.pf_uq_bases_aligned as f64;
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

impl Collector {
    /// `PerUnitTargetMetricCollector.finish` and `HsMetricCollector.convertMetric`.
    fn finish(self, territory: Territory) -> Result<Report> {
        let c = &self.counts;
        let targets = self.targets.intervals();
        let mut row = Row {
            bait_set: territory.bait_set,
            bait_territory: territory.bait_territory,
            target_territory: territory.target_territory,
            genome_size: territory.genome_size,
            on_bait_bases: c.on_probe_bases,
            near_bait_bases: c.near_probe_bases,
            off_bait_bases: c.off_probe_bases,
            total_reads: c.total_reads,
            pf_reads: c.pf_reads,
            pf_bases: c.pf_bases,
            pf_unique_reads: c.pf_unique_reads,
            pf_uq_reads_aligned: c.pf_uq_reads_aligned,
            pf_bases_aligned: c.pf_bases_aligned,
            pf_uq_bases_aligned: c.pf_uq_bases_aligned,
            on_target_bases: c.on_target_bases,
            ..Row::default()
        };
        row.pct_pf_reads = c.pf_reads as f64 / c.total_reads as f64;
        row.pct_pf_uq_reads = c.pf_unique_reads as f64 / c.total_reads as f64;
        row.pct_pf_uq_reads_aligned = c.pf_uq_reads_aligned as f64 / c.pf_unique_reads as f64;
        let denominator = (c.on_probe_bases + c.near_probe_bases + c.off_probe_bases) as f64;
        row.pct_selected_bases = (c.on_probe_bases + c.near_probe_bases) as f64 / denominator;
        row.pct_off_bait = c.off_probe_bases as f64 / denominator;
        row.on_bait_vs_selected =
            c.on_probe_bases as f64 / (c.on_probe_bases + c.near_probe_bases) as f64;
        row.mean_bait_coverage = c.on_probe_bases as f64 / row.bait_territory as f64;
        row.fold_enrichment = (c.on_probe_bases as f64 / denominator)
            / (row.bait_territory as f64 / row.genome_size as f64);
        let aligned = c.pf_bases_aligned as f64;
        row.pct_exc_dupe = c.excluded_dupe as f64 / aligned;
        row.pct_exc_adapter = c.excluded_adapter as f64 / aligned;
        row.pct_exc_mapq = c.excluded_mapq as f64 / aligned;
        row.pct_exc_baseq = c.excluded_baseq as f64 / aligned;
        row.pct_exc_overlap = c.excluded_overlap as f64 / aligned;
        row.pct_exc_off_target = c.excluded_off_target as f64 / aligned;

        // calculateTargetCoverageMetrics. Depths are tallied per value, then
        // put in the histogram once: the same integer counts Picard adds bin
        // by bin.
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
        let depth = depth_histogram(&depth_counts, hq_max);
        // Bases at or above each TARGET_BASE_DEPTHS entry.
        let mut target_bases = [0i64; TARGET_BASE_DEPTHS.len()];
        for (d, &count) in depth_counts.iter().enumerate() {
            for (bases, &threshold) in target_bases.iter_mut().zip(&TARGET_BASE_DEPTHS) {
                if d as i64 >= threshold {
                    *bases += count;
                }
            }
        }
        row.mean_target_coverage = total_coverage as f64 / row.target_territory as f64;
        row.median_target_coverage = depth.median();
        row.max_target_coverage = self.hq_max_depth;
        row.min_target_coverage = min_depth.min(self.hq_max_depth);
        let percentile_20 = depth.percentile(0.2).ok_or_else(|| {
            Error::Metric(
                "Cannot calculate percentiles when total is zero (no target bases)".to_owned(),
            )
        })?;
        row.fold_80_base_penalty = row.mean_target_coverage / percentile_20;
        row.zero_cvg_targets_pct = zero_targets as f64 / territory.raw_target_count as f64;
        for (pct, &bases) in row.pct_target_bases.iter_mut().zip(&target_bases[1..]) {
            *pct = bases as f64 / target_bases[0] as f64;
        }

        // calculateTheoreticalHetSensitivity: unfiltered depths, capped.
        let cap = usize::try_from(self.opts.coverage_cap).unwrap_or(0);
        let mut capped_counts = vec![0i64; cap + 1];
        for (i, target) in targets.iter().enumerate() {
            let depths = &self.coverage.unfiltered[self.coverage.range(i)];
            if !depths.iter().any(|&d| d > 0) {
                capped_counts[0] += target.length();
                continue;
            }
            for &d in depths {
                capped_counts[(d as usize).min(cap)] += 1;
            }
        }
        let unfiltered_depth = depth_histogram(&capped_counts, cap);
        let mut capped_quality = Histogram::new();
        let mut base_quality = Histogram::new();
        for q in 0..QUALITY_SLOTS {
            capped_quality.increment_by(q as i64, self.base_q_capped[q] as f64);
            base_quality.increment_by(q as i64, self.base_q[q] as f64);
        }
        row.het_snp_sensitivity = het_snp_sensitivity(
            &normalize(&unfiltered_depth),
            &normalize(&capped_quality),
            self.opts.sample_size,
            LOG_ODDS_THRESHOLD,
        )?;
        row.het_snp_q = f64::from(phred_from_error_probability(1.0 - row.het_snp_sensitivity));

        // calculateGcMetrics, only with a reference.
        if let Some(gc) = &territory.target_gc {
            let mut target_by_gc = [0i64; GC_BINS];
            let mut aligned_by_gc = [0i64; GC_BINS];
            for (i, target) in targets.iter().enumerate() {
                if target.length() <= 0 {
                    continue;
                }
                let bin = usize::try_from(math_round(gc[i] * 100.0)).unwrap_or(0);
                target_by_gc[bin] += target.length();
                aligned_by_gc[bin] += self.coverage.high_quality[self.coverage.range(i)]
                    .iter()
                    .map(|&d| i64::from(d))
                    .sum::<i64>();
            }
            let total_target: i64 = target_by_gc.iter().sum();
            let total_aligned: i64 = aligned_by_gc.iter().sum();
            for i in 0..GC_BINS {
                let target_pct = target_by_gc[i] as f64 / total_target as f64;
                let aligned_pct = aligned_by_gc[i] as f64 / total_aligned as f64;
                let dropout = (aligned_pct - target_pct) * 100.0;
                if dropout < 0.0 {
                    if i <= 50 {
                        row.at_dropout += dropout.abs();
                    }
                    if i >= 50 {
                        row.gc_dropout += dropout.abs();
                    }
                }
            }
        }

        // HsMetricCollector.convertMetric.
        row.bait_design_efficiency = row.target_territory as f64 / row.bait_territory as f64;
        row.pct_usable_bases_on_bait = row.on_bait_bases as f64 / c.pf_bases as f64;
        row.pct_usable_bases_on_target = row.on_target_bases as f64 / c.pf_bases as f64;
        row.hs_library_size =
            estimate_library_size(c.pf_selected_pairs, c.pf_selected_unique_pairs);
        row.hs_penalty = PENALTY_GOALS.map(|goal| hs_penalty(row.hs_library_size, c, &row, goal));
        Ok(Report {
            row,
            depth,
            base_quality,
        })
    }
}

/// `CollectHsMetrics` as a pass [`Accumulator`].
pub struct HsMetrics {
    collector: Collector,
    territory: Territory,
}

impl HsMetrics {
    /// Reads and uniques the interval files against the BAM header, and
    /// with a reference computes the GC fraction of every target.
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        if opts.bait_intervals.is_empty() || opts.target_intervals.is_empty() {
            return Err(Error::Usage(
                "CollectHsMetrics needs at least one bait and one target interval file".to_owned(),
            ));
        }
        if opts.coverage_cap < 0 || opts.near_distance < 0 {
            return Err(Error::Usage(
                "COVERAGE_CAP and NEAR_DISTANCE must not be negative".to_owned(),
            ));
        }
        let contigs = header.reference_sequences().len();
        let raw_targets = intervals::read_files(&opts.target_intervals, header)?;
        let unique_targets = intervals::uniqued(&raw_targets);
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
            bait_set: bait_set_name(opts),
            bait_territory: intervals::count_bases(&unique_baits),
            target_territory: intervals::count_bases(&unique_targets),
            genome_size: header
                .reference_sequences()
                .values()
                .map(|c| c.length().get() as i64)
                .sum(),
            raw_target_count: raw_targets.len(),
            target_gc,
        };
        let collector = Collector {
            opts: opts.clone(),
            coverage: Coverage::new(&unique_targets),
            targets: OverlapIndex::new(unique_targets, contigs, 0),
            probes: OverlapIndex::new(unique_baits, contigs, i64::from(opts.near_distance)),
            hq_max_depth: 0,
            base_q_capped: [0; QUALITY_SLOTS],
            base_q: [0; QUALITY_SLOTS],
            counts: Counts::default(),
            adapters: Adapters::default(),
            cigar: Vec::new(),
            bases: Vec::new(),
        };
        Ok(Self {
            collector,
            territory,
        })
    }
}

impl Accumulator for HsMetrics {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        self.collector.accept(record)
    }

    fn finish(self) -> Result<Report> {
        self.collector.finish(self.territory)
    }
}

/// `StringHeader` command line in Picard's shape.
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
    if let Some(name) = &opts.bait_set_name {
        command.push_str(&format!(" BAIT_SET_NAME={name}"));
    }
    if let Some(reference) = &opts.reference {
        command.push_str(&format!(" REFERENCE_SEQUENCE={}", reference.display()));
    }
    command.push_str(&format!(
        " METRIC_ACCUMULATION_LEVEL=[ALL_READS] NEAR_DISTANCE={} MINIMUM_MAPPING_QUALITY={} MINIMUM_BASE_QUALITY={} CLIP_OVERLAPPING_READS={} INCLUDE_INDELS={} COVERAGE_CAP={} SAMPLE_SIZE={}",
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
        HsMetrics::new(opts, header)
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
    fn depth_histogram_prefills_empty_bins() {
        let h = depth_histogram(&[0, 0, 3, 0, 1], 4);
        assert_eq!(h.keys().collect::<Vec<_>>(), [0, 1, 2, 3, 4]);
        let h = depth_histogram(&[5, 0, 0], 0);
        assert_eq!(h.keys().collect::<Vec<_>>(), [0]);
    }

    #[test]
    fn normalize_indexes_by_key() {
        let h = depth_histogram(&[1, 3], 2);
        assert_eq!(normalize(&h), [0.25, 0.75]);
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
    fn row_fields_match_columns() {
        let row = Row {
            bait_set: "baits".to_owned(),
            hs_library_size: None,
            het_snp_q: 22.0,
            ..Row::default()
        };
        let fields = row.fields();
        assert_eq!(fields.len(), COLUMNS.len());
        assert_eq!(fields[0].picard_text(), "baits");
        assert_eq!(fields[13].picard_text(), "");
        assert_eq!(fields[65].picard_text(), "22");
        assert_eq!(Field::Double(f64::NAN).picard_text(), "?");
    }
}

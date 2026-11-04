//! `CollectGcBiasMetrics`, mirroring Picard 3.5.0 (`CollectGcBiasMetrics`,
//! `GcBiasMetricsCollector`, `GcBiasUtils`).
//!
//! Two files: `OUTPUT` holds the detail rows (one per GC percentage, 0 to
//! 100, per accumulation unit and per `READS_USED` value), `SUMMARY_OUTPUT`
//! one summary row per unit and `READS_USED` value.
//!
//! The reference is scanned once at setup: every window of
//! `SCAN_WINDOW_SIZE` bases starting at offsets 1 up to, not including,
//! `length - SCAN_WINDOW_SIZE` of every contig is binned by its GC
//! percentage (integer division, `N` in the denominator), and a window with
//! 5 or more `N` is dropped. The first and the last windows of a contig are
//! never scanned. A read's window starts at its alignment start, or at its
//! alignment end minus the window size for a reverse read; a start at 0 or
//! below is skipped, and a start past the scanned range (the contig end, or
//! a contig shorter than the window) reads bin 0 because the per-contig
//! table is zero there, as Picard's byte array is.
//!
//! Per unit: `TOTAL_CLUSTERS` counts unpaired and first-of-pair records,
//! mapped or not; `ALIGNED_READS` every mapped record, secondary and
//! supplementary ones included; `READ_STARTS` the mapped records whose
//! window has a GC bin; the mean base quality of a bin comes from the
//! mismatches plus inserted and deleted bases over the read length of those
//! records. A record with no bases (`*`) is dropped while a unit has seen
//! fewer than 100 of them, then counted. Rows are written only for a unit
//! that saw an aligned read. With `ALSO_IGNORE_DUPLICATES` each unit also
//! keeps the same counts without duplicate-flagged records, written after
//! its `ALL` rows as `UNIQUE`. The level label is `All Reads`, `Sample`,
//! `Library` or `Read Group`, and only the column of that level is set.
//! `MINIMUM_GENOME_FRACTION` is accepted; Picard 3.5.0 never reads it.

use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::{bam, sam};

use crate::blocks::{alignment_blocks, alignment_end};
use crate::error::{Error, Result};
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator, Collector};
use crate::record::{position, record_name};
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sequence::{bases_equal, bisulfite_bases_equal};

/// GC percentages 0 to 100.
const BINS: usize = 101;
/// `*`-sequence records a unit drops before it starts counting them.
const STAR_RECORDS_DROPPED: u32 = 100;

#[derive(Clone, Debug)]
pub struct Options {
    pub reference: std::path::PathBuf,
    pub scan_window_size: i64,
    pub bisulfite: bool,
    pub also_ignore_duplicates: bool,
    pub levels: Vec<Level>,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let reference = inv.get("REFERENCE_SEQUENCE").ok_or_else(|| {
            Error::Usage("CollectGcBiasMetrics: REFERENCE_SEQUENCE is required".to_owned())
        })?;
        let _: f64 = inv.parse("MINIMUM_GENOME_FRACTION", 0.00001)?;
        Ok(Self {
            reference: reference.into(),
            scan_window_size: inv.parse::<i32>("SCAN_WINDOW_SIZE", 100)?.into(),
            bisulfite: inv.flag("IS_BISULFITE_SEQUENCED", false)?,
            also_ignore_duplicates: inv.flag("ALSO_IGNORE_DUPLICATES", false)?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
        })
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectGcBiasMetrics",
    summary: "the detail rows (Picard's OUTPUT file), one row per GC bin, accumulation unit \
              and READS_USED value; SUMMARY_OUTPUT gets the GcBiasSummaryMetrics file.",
    options: &[
        Opt::ignored(&["CHART_OUTPUT", "CHART"], "no chart is drawn"),
        Opt::path(&["SUMMARY_OUTPUT", "S"], "The GcBiasSummaryMetrics file.").required(),
        Opt::path(
            &["REFERENCE_SEQUENCE", "R"],
            "Reference FASTA with a .fai index.",
        )
        .required(),
        Opt::int(
            &["SCAN_WINDOW_SIZE", "WINDOW_SIZE"],
            "Window length for the GC content (default 100).",
        ),
        Opt::float(
            &["MINIMUM_GENOME_FRACTION", "MGF"],
            "Accepted and unused, as in Picard 3.5.0 (default 0.00001).",
        ),
        Opt::bool(
            &["IS_BISULFITE_SEQUENCED", "BS"],
            "Count C after T as bisulfite converted (default false).",
        ),
        crate::registry::LEVEL_OPTION,
        Opt::bool(
            &["ALSO_IGNORE_DUPLICATES"],
            "Add rows with duplicates left out (default false).",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    let summary_path = inv.get("SUMMARY_OUTPUT").ok_or_else(|| {
        Error::Usage("CollectGcBiasMetrics: SUMMARY_OUTPUT (S) is required".to_owned())
    })?;
    let summary_path = std::path::PathBuf::from(summary_path);
    let opts = Options::from_invocation(inv)?;
    Ok(pass::boxed(GcBias::new(&opts, header)?, move |report| {
        Ok(ToolOutput {
            main: Some(MetricsOutput::from_rows(&report.details, None)),
            extra: vec![(
                summary_path,
                MetricsOutput::from_rows(&report.summaries, None),
            )],
        })
    }))
}

#[derive(Clone, Debug, PartialEq)]
pub struct DetailRow {
    pub accumulation_level: &'static str,
    pub reads_used: &'static str,
    pub gc: i64,
    pub windows: i64,
    pub read_starts: i64,
    pub mean_base_quality: i64,
    pub normalized_coverage: f64,
    pub error_bar_width: f64,
    pub labels: Labels,
}

impl MetricRow for DetailRow {
    const CLASS: &'static str = "picard.analysis.GcBiasDetailMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "ACCUMULATION_LEVEL",
        "READS_USED",
        "GC",
        "WINDOWS",
        "READ_STARTS",
        "MEAN_BASE_QUALITY",
        "NORMALIZED_COVERAGE",
        "ERROR_BAR_WIDTH",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.accumulation_level.into(),
            self.reads_used.into(),
            self.gc.into(),
            self.windows.into(),
            self.read_starts.into(),
            self.mean_base_quality.into(),
            self.normalized_coverage.into(),
            self.error_bar_width.into(),
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct SummaryRow {
    pub accumulation_level: &'static str,
    pub reads_used: &'static str,
    pub window_size: i64,
    pub total_clusters: i64,
    pub aligned_reads: i64,
    pub at_dropout: f64,
    pub gc_dropout: f64,
    pub gc_nc: [f64; 5],
    pub labels: Labels,
}

impl MetricRow for SummaryRow {
    const CLASS: &'static str = "picard.analysis.GcBiasSummaryMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "ACCUMULATION_LEVEL",
        "READS_USED",
        "WINDOW_SIZE",
        "TOTAL_CLUSTERS",
        "ALIGNED_READS",
        "AT_DROPOUT",
        "GC_DROPOUT",
        "GC_NC_0_19",
        "GC_NC_20_39",
        "GC_NC_40_59",
        "GC_NC_60_79",
        "GC_NC_80_100",
        "SAMPLE",
        "LIBRARY",
        "READ_GROUP",
    ];

    fn values(&self) -> Vec<Value> {
        let mut values = vec![
            self.accumulation_level.into(),
            self.reads_used.into(),
            self.window_size.into(),
            self.total_clusters.into(),
            self.aligned_reads.into(),
            self.at_dropout.into(),
            self.gc_dropout.into(),
        ];
        values.extend(self.gc_nc.iter().map(|&v| Value::Float(v)));
        values.extend([
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]);
        values
    }
}

/// Both files' rows.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Report {
    pub details: Vec<DetailRow>,
    pub summaries: Vec<SummaryRow>,
}

/// `GcBiasUtils.calculateGc` over a sliding window: the count of G and C
/// and of N, updated base by base.
struct GcState {
    gc: i64,
    n: i64,
    prior: Option<u8>,
}

impl GcState {
    fn new() -> Self {
        Self {
            gc: 0,
            n: 0,
            prior: None,
        }
    }

    /// The window `bases[start..end]`: -1 with more than 4 `N`, else the GC
    /// percentage by integer division. The first call counts the whole
    /// window; later calls add the last base and remove the base before
    /// the window. The base added tests `N` by exact byte, the base removed
    /// by IUPAC mask, as Picard's two branches do.
    fn gc_percent(&mut self, bases: &[u8], start: usize, end: usize) -> i64 {
        let is_gc = |b: u8| bases_equal(b, b'G') || bases_equal(b, b'C');
        match self.prior {
            None => {
                self.gc = 0;
                self.n = 0;
                for &b in &bases[start..end] {
                    if is_gc(b) {
                        self.gc += 1;
                    } else if bases_equal(b, b'N') {
                        self.n += 1;
                    }
                }
            }
            Some(prior) => {
                let new = bases[end - 1];
                if is_gc(new) {
                    self.gc += 1;
                } else if new == b'N' {
                    self.n += 1;
                }
                if is_gc(prior) {
                    self.gc -= 1;
                } else if bases_equal(prior, b'N') {
                    self.n -= 1;
                }
            }
        }
        self.prior = Some(bases[start]);
        if self.n > 4 {
            -1
        } else {
            self.gc * 100 / (end - start) as i64
        }
    }
}

/// `StringUtil.toUpperCase`: ASCII letters only.
fn upper(bases: &mut [u8]) {
    for b in bases {
        b.make_ascii_uppercase();
    }
}

/// The windows Picard scans: each 1-based offset from 1 up to, not
/// including, `length - window`, with its GC bin (-1 for too many `N`).
fn scan(bases: &[u8], window: i64, mut each: impl FnMut(usize, i64)) {
    let last = bases.len() as i64 - window;
    let mut state = GcState::new();
    let mut i = 1i64;
    while i < last {
        let (start, end) = (i as usize, (i + window) as usize);
        each(start, state.gc_percent(bases, start, end));
        i += 1;
    }
}

/// `calculateAllGcs`: the GC bin of the window starting at each scanned
/// 1-based offset, 0 elsewhere (Picard's zeroed byte array, one longer
/// than the contig).
fn window_bins(bases: &[u8], window: i64) -> Vec<i8> {
    let mut bins = vec![0i8; bases.len() + 1];
    scan(bases, window, |start, bin| bins[start] = bin as i8);
    bins
}

/// `calculateRefWindowsByGc`: windows of the whole reference by GC bin.
fn windows_by_gc(reference: &mut Reference, window: i64) -> Result<[i64; BINS]> {
    let mut windows = [0i64; BINS];
    for id in 0..reference.len() {
        let mut bases = reference.load(id)?;
        upper(&mut bases);
        scan(&bases, window, |_, bin| {
            if bin != -1 {
                windows[bin as usize] += 1;
            }
        });
    }
    Ok(windows)
}

/// `GcObject`: one unit's counts for one `READS_USED` value.
#[derive(Clone, Debug)]
struct Counts {
    total_clusters: i64,
    total_aligned_reads: i64,
    reads_by_gc: [i64; BINS],
    bases_by_gc: [i64; BINS],
    errors_by_gc: [i64; BINS],
}

impl Default for Counts {
    fn default() -> Self {
        Self {
            total_clusters: 0,
            total_aligned_reads: 0,
            reads_by_gc: [0; BINS],
            bases_by_gc: [0; BINS],
            errors_by_gc: [0; BINS],
        }
    }
}

/// One accumulation unit: all reads, and without duplicates when asked.
#[derive(Clone, Debug, Default)]
struct Unit {
    all: Counts,
    unique: Option<Counts>,
    star_records: u32,
}

/// What one record adds to a unit.
struct Observation {
    cluster: bool,
    aligned: bool,
    /// The GC bin of the read's window, with the read length and its
    /// error count, for a mapped read whose window start is positive and
    /// binned.
    window: Option<(usize, i64, i64)>,
    duplicate: bool,
    star: bool,
}

impl Counts {
    fn add(&mut self, o: &Observation) {
        if o.cluster {
            self.total_clusters += 1;
        }
        if o.aligned {
            self.total_aligned_reads += 1;
        }
        if let Some((bin, length, errors)) = o.window {
            self.reads_by_gc[bin] += 1;
            self.bases_by_gc[bin] += length;
            self.errors_by_gc[bin] += errors;
        }
    }
}

/// The pass accumulator.
pub struct GcBias {
    opts: Options,
    reference: Reference,
    windows_by_gc: [i64; BINS],
    /// GC bins of the current contig's windows, by 1-based start.
    bins: Vec<i8>,
    contig: Option<usize>,
    levels: Levels<Unit>,
    bases: Vec<u8>,
}

impl GcBias {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        let mut reference = Reference::open(&opts.reference, header)?;
        let windows_by_gc = windows_by_gc(&mut reference, opts.scan_window_size)?;
        Ok(Self {
            opts: opts.clone(),
            reference,
            windows_by_gc,
            bins: Vec::new(),
            contig: None,
            levels: Levels::new(
                header,
                &opts.levels,
                Unit {
                    unique: opts.also_ignore_duplicates.then(Counts::default),
                    ..Unit::default()
                },
            ),
            bases: Vec::new(),
        })
    }

    /// What the record adds. `counted` is false for a record without
    /// bases that every unit it reaches still drops: its window and errors
    /// are then not computed, since Picard returns before touching them.
    fn observe(&mut self, record: &bam::Record, counted: bool) -> Result<Observation> {
        let flags = record.flags();
        self.bases.clear();
        self.bases.extend(record.sequence().iter());
        let cluster = !flags.is_segmented() || flags.is_first_segment();
        let duplicate = flags.is_duplicate();
        let star = self.bases.is_empty();
        if flags.is_unmapped() || !counted {
            return Ok(Observation {
                cluster,
                aligned: false,
                window: None,
                duplicate,
                star,
            });
        }
        let id =
            record
                .reference_sequence_id()
                .transpose()?
                .ok_or_else(|| Error::MalformedRecord {
                    name: record_name(record),
                    what: "mapped without a reference sequence",
                })?;
        if self.contig != Some(id) {
            let mut contig = self.reference.contig(id)?.to_vec();
            upper(&mut contig);
            self.bins = window_bins(&contig, self.opts.scan_window_size);
            self.contig = Some(id);
        }
        let start = position(
            record,
            record.alignment_start(),
            "mapped without a position",
        )?;
        let pos = if flags.is_reverse_complemented() {
            alignment_end(record)?.unwrap_or(start) - self.opts.scan_window_size
        } else {
            start
        };
        let mut window = None;
        if pos > 0 {
            let bin = *self
                .bins
                .get(pos as usize)
                .ok_or_else(|| Error::MalformedRecord {
                    name: record_name(record),
                    what: "alignment past the end of the contig",
                })?;
            if bin >= 0 {
                let contig = self.reference.contig(id)?;
                let errors = errors(record, &self.bases, contig, self.opts.bisulfite)?;
                window = Some((bin as usize, self.bases.len() as i64, errors));
            }
        }
        Ok(Observation {
            cluster,
            aligned: true,
            window,
            duplicate,
            star,
        })
    }
}

/// Mismatches (`SequenceUtil.countMismatches`) plus inserted and deleted
/// bases of a mapped record against its contig.
fn errors(record: &bam::Record, bases: &[u8], contig: &[u8], bisulfite: bool) -> Result<i64> {
    let reverse = record.flags().is_reverse_complemented();
    let malformed = |what| Error::MalformedRecord {
        name: record_name(record),
        what,
    };
    let mut errors = 0i64;
    for block in alignment_blocks(record)? {
        for i in 0..block.length {
            let read = *bases
                .get(block.read_start - 1 + i)
                .ok_or_else(|| malformed("CIGAR longer than the sequence"))?;
            let reference = *contig
                .get((block.reference_start - 1) as usize + i)
                .ok_or_else(|| malformed("alignment past the end of the contig"))?;
            let equal = if bisulfite {
                bisulfite_bases_equal(reverse, read, reference)
            } else {
                bases_equal(read, reference)
            };
            errors += i64::from(!equal);
        }
    }
    for op in record.cigar().iter() {
        let op = op?;
        if matches!(op.kind(), Kind::Insertion | Kind::Deletion) {
            errors += op.len() as i64;
        }
    }
    Ok(errors)
}

impl Accumulator for GcBias {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let star = record.sequence().is_empty();
        let mut counted = !star;
        if star {
            // The unit decides after its first 100 such records; the
            // `unknown` unit comes into being here, as Picard's does.
            self.levels.route(record, |unit| {
                counted |= unit.star_records >= STAR_RECORDS_DROPPED;
            });
        }
        let o = self.observe(record, counted)?;
        self.levels.route(record, |unit| {
            if o.star && unit.star_records < STAR_RECORDS_DROPPED {
                unit.star_records += 1;
                return;
            }
            unit.all.add(&o);
            if let (Some(unique), false) = (&mut unit.unique, o.duplicate) {
                unique.add(&o);
            }
        });
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        let mut report = Report::default();
        for bucket in self.levels.into_buckets() {
            let (level, labels) = level_labels(bucket.level, &bucket.labels);
            let mut add = |counts: &Counts, reads_used: &'static str| {
                if counts.total_aligned_reads > 0 {
                    add_unit(
                        &mut report,
                        counts,
                        &self.windows_by_gc,
                        self.opts.scan_window_size,
                        level,
                        reads_used,
                        &labels,
                    );
                }
            };
            add(&bucket.state.all, "ALL");
            if let Some(unique) = &bucket.state.unique {
                add(unique, "UNIQUE");
            }
        }
        Ok(report)
    }
}

/// `PerUnitGcBiasMetricsCollector.addReadToGcData`: the level label and the
/// one label column set, chosen by the first of read group, library and
/// sample the unit has.
fn level_labels(level: Level, labels: &Labels) -> (&'static str, Labels) {
    let one = |sample, library, read_group| Labels {
        sample,
        library,
        read_group,
    };
    match (level, labels) {
        (_, l) if l.read_group.is_some() => ("Read Group", one(None, None, l.read_group.clone())),
        (_, l) if l.library.is_some() => ("Library", one(None, l.library.clone(), None)),
        (_, l) if l.sample.is_some() => ("Sample", one(l.sample.clone(), None, None)),
        _ => ("All Reads", Labels::default()),
    }
}

/// `addGcDataToFile` for one unit and one `READS_USED` value.
fn add_unit(
    report: &mut Report,
    counts: &Counts,
    windows_by_gc: &[i64; BINS],
    window_size: i64,
    level: &'static str,
    reads_used: &'static str,
    labels: &Labels,
) {
    let total_windows: f64 = windows_by_gc.iter().map(|&w| w as f64).sum();
    let total_reads: f64 = counts.reads_by_gc.iter().map(|&r| r as f64).sum();
    let mean_reads_per_window = total_reads / total_windows;
    let first = report.details.len();
    for (i, &windows) in windows_by_gc.iter().enumerate() {
        let read_starts = counts.reads_by_gc[i];
        let mean_base_quality = if counts.errors_by_gc[i] > 0 {
            phred_from_errors(counts.bases_by_gc[i] as f64, counts.errors_by_gc[i] as f64)
        } else {
            0
        };
        let (normalized_coverage, error_bar_width) = if windows != 0 {
            (
                (read_starts as f64 / windows as f64) / mean_reads_per_window,
                ((read_starts as f64).sqrt() / windows as f64) / mean_reads_per_window,
            )
        } else {
            (0.0, 0.0)
        };
        report.details.push(DetailRow {
            accumulation_level: level,
            reads_used,
            gc: i as i64,
            windows,
            read_starts,
            mean_base_quality,
            normalized_coverage,
            error_bar_width,
            labels: labels.clone(),
        });
    }
    let gc_nc = |start: usize, end: usize| {
        let mut windows_total = 0i64;
        let mut sum = 0.0;
        for (windows, reads) in windows_by_gc[start..=end]
            .iter()
            .zip(&counts.reads_by_gc[start..=end])
        {
            if *windows != 0 {
                sum += *reads as f64;
                windows_total += windows;
            }
        }
        if windows_total == 0 {
            0.0
        } else {
            sum / (windows_total as f64 * mean_reads_per_window)
        }
    };
    let (at_dropout, gc_dropout) = dropout(&report.details[first..]);
    report.summaries.push(SummaryRow {
        accumulation_level: level,
        reads_used,
        window_size,
        total_clusters: counts.total_clusters,
        aligned_reads: counts.total_aligned_reads,
        at_dropout,
        gc_dropout,
        gc_nc: [
            gc_nc(0, 19),
            gc_nc(20, 39),
            gc_nc(40, 59),
            gc_nc(60, 79),
            gc_nc(80, 100),
        ],
        labels: labels.clone(),
    });
}

/// `calculateDropoutMetrics`: the positive differences between the share of
/// windows and the share of read starts in each bin, times 100, summed over
/// GC at most 50 (AT) and above (GC).
fn dropout(details: &[DetailRow]) -> (f64, f64) {
    let total_reads: f64 = details.iter().map(|d| d.read_starts as f64).sum();
    let total_windows: f64 = details.iter().map(|d| d.windows as f64).sum();
    let (mut at, mut gc) = (0.0, 0.0);
    for d in details {
        let relative_reads = d.read_starts as f64 / total_reads;
        let relative_windows = d.windows as f64 / total_windows;
        let dropout = (relative_windows - relative_reads) * 100.0;
        if dropout > 0.0 {
            if d.gc <= 50 {
                at += dropout;
            } else {
                gc += dropout;
            }
        }
    }
    (at, gc)
}

/// `QualityUtil.getPhredScoreFromObsAndErrors`: `Math.round(-10 log10(e/o))`,
/// Java rounding (half up toward positive infinity).
fn phred_from_errors(observations: f64, errors: f64) -> i64 {
    let value = -10.0 * (errors / observations).log10();
    (value + 0.5).floor() as i64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn windows_skip_the_contig_ends_and_n_runs() {
        // 13 bases, window 4: offsets 1 to 8 are scanned, 0 and 9 are not.
        let bases = b"ACGTNNNNNGGCC";
        let bins = window_bins(bases, 4);
        assert_eq!(bins.len(), 14);
        assert_eq!(&bins[..10], &[0, 50, 25, 0, 0, 0, 25, 50, 75, 0]);
        // Window 5: NNNNN at offset 4 has 5 N and is dropped.
        let bins = window_bins(bases, 5);
        assert_eq!(&bins[..9], &[0, 40, 20, 0, -1, 20, 40, 60, 0]);
        let mut scanned = 0;
        scan(bases, 20, |_, _| scanned += 1);
        assert_eq!(scanned, 0);
    }

    #[test]
    fn gc_state_slides() {
        let bases = b"GGCCAAAAGG";
        let mut state = GcState::new();
        assert_eq!(state.gc_percent(bases, 0, 4), 100);
        assert_eq!(state.gc_percent(bases, 1, 5), 75);
        assert_eq!(state.gc_percent(bases, 2, 6), 50);
        assert_eq!(state.gc_percent(bases, 3, 7), 25);
        assert_eq!(state.gc_percent(bases, 4, 8), 0);
        assert_eq!(state.gc_percent(bases, 5, 9), 25);
        assert_eq!(state.gc_percent(bases, 6, 10), 50);
    }

    #[test]
    fn phred_rounds_half_up() {
        assert_eq!(phred_from_errors(1000.0, 1.0), 30);
        assert_eq!(phred_from_errors(100.0, 100.0), 0);
        assert_eq!(phred_from_errors(10.0, 20.0), -3);
        assert_eq!(phred_from_errors(1.0, 1.778_279_41), -2);
    }

    #[test]
    fn dropout_sums_positive_differences() {
        let row = |gc, windows, read_starts| DetailRow {
            accumulation_level: "All Reads",
            reads_used: "ALL",
            gc,
            windows,
            read_starts,
            mean_base_quality: 0,
            normalized_coverage: 0.0,
            error_bar_width: 0.0,
            labels: Labels::default(),
        };
        let rows = [row(10, 50, 10), row(50, 25, 40), row(90, 25, 0)];
        let (at, gc) = dropout(&rows);
        assert!((at - 30.0).abs() < 1e-9);
        assert!((gc - 25.0).abs() < 1e-9);
    }

    #[test]
    fn labels_follow_the_unit_level() {
        let labels = Labels {
            sample: Some("S".to_owned()),
            library: Some("L".to_owned()),
            read_group: None,
        };
        let (level, l) = level_labels(Level::Library, &labels);
        assert_eq!(level, "Library");
        assert_eq!((l.sample, l.library.as_deref()), (None, Some("L")));
        let (level, l) = level_labels(Level::AllReads, &Labels::default());
        assert_eq!(level, "All Reads");
        assert_eq!(l, Labels::default());
    }
}

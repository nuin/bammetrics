//! `CollectInsertSizeMetrics` (see docs/01 section 2.2).

use std::io::{self, Write};
use std::path::Path;

use noodles::{bam, sam};

use crate::error::Result;
use crate::histogram::Histogram;
use crate::levels::{Labels, Level, Levels, parse_levels};
use crate::metric::{HistogramSection, HistogramSeries, MetricRow, MetricsOutput, Value};
use crate::pass::{self, Accumulator};
pub use crate::record::Orientation;
use crate::record::orientation;
use crate::registry::{Invocation, Opt, Tool};

#[derive(Clone, Debug)]
pub struct Options {
    /// Orientations with fewer pairs than this fraction of all pairs are dropped.
    /// `f32` because `CollectInsertSizeMetrics.MINIMUM_PCT` is a Java `float`;
    /// see `passes_minimum_pct`.
    pub minimum_pct: f32,
    /// Mean and SD use bins up to median + deviations * MAD.
    pub deviations: f64,
    /// `HISTOGRAM_WIDTH`: trim to this key instead of median + deviations * MAD.
    pub histogram_width: Option<i64>,
    /// `MIN_HISTOGRAM_WIDTH`: lower bound on the computed trim width.
    pub min_histogram_width: Option<i64>,
    pub include_duplicates: bool,
    /// `METRIC_ACCUMULATION_LEVEL`: one set of rows per level.
    pub levels: Vec<Level>,
}

impl Options {
    /// The options of a parsed `CollectInsertSizeMetrics` command line.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let defaults = Self::default();
        Ok(Self {
            minimum_pct: inv.parse("MINIMUM_PCT", defaults.minimum_pct)?,
            deviations: inv.parse("DEVIATIONS", defaults.deviations)?,
            histogram_width: inv.parse_opt("HISTOGRAM_WIDTH")?,
            min_histogram_width: inv.parse_opt("MINIMUM_HISTOGRAM_WIDTH")?,
            include_duplicates: inv.flag("INCLUDE_DUPLICATES", false)?,
            levels: parse_levels(inv.all("METRIC_ACCUMULATION_LEVEL"))?,
        })
    }
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_pct: 0.05,
            deviations: 10.0,
            histogram_width: None,
            min_histogram_width: None,
            include_duplicates: false,
            levels: vec![Level::AllReads],
        }
    }
}

/// `CollectInsertSizeMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectInsertSizeMetrics",
    options: &[
        Opt::ignored(&["HISTOGRAM_FILE", "H"], "no histogram plot is drawn"),
        Opt::value(&["DEVIATIONS"]),
        Opt::value(&["HISTOGRAM_WIDTH", "W"]),
        Opt::value(&["MINIMUM_HISTOGRAM_WIDTH", "MW"]),
        Opt::value(&["MINIMUM_PCT", "M"]),
        Opt::value(&["INCLUDE_DUPLICATES"]),
        Opt::value(&["METRIC_ACCUMULATION_LEVEL", "LEVEL"]),
    ],
    run: run_tool,
};

/// The rows Picard writes, or `None` (with Picard's warning on standard
/// error) when it writes no file.
fn run_tool(inv: &Invocation) -> Result<Option<MetricsOutput>> {
    let opts = Options::from_invocation(inv)?;
    let rows = pass::run(&inv.input, inv.threads, |header| {
        Ok(InsertSize::new(&opts, header))
    })?;
    if rows.is_empty() {
        eprintln!("{}", no_rows_warning(opts.minimum_pct));
        return Ok(None);
    }
    Ok(Some(output(&rows)))
}

/// The warning `CollectInsertSizeMetrics.finish` logs when it writes no file.
pub fn no_rows_warning(minimum_pct: f32) -> String {
    format!(
        "bammetrics: all data categories were discarded because they contained < {minimum_pct:?} of the total aligned paired data; no metrics written"
    )
}

/// The X of each `WIDTH_OF_X_PERCENT` column, in column order.
pub const WIDTH_PERCENTS: [u32; 11] = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99];

/// `InsertSizeMetrics` columns in Picard 3.5.0 file order.
pub const COLUMNS: [&str; 23] = [
    "MEDIAN_INSERT_SIZE",
    "MODE_INSERT_SIZE",
    "MEDIAN_ABSOLUTE_DEVIATION",
    "MIN_INSERT_SIZE",
    "MAX_INSERT_SIZE",
    "MEAN_INSERT_SIZE",
    "STANDARD_DEVIATION",
    "READ_PAIRS",
    "PAIR_ORIENTATION",
    "WIDTH_OF_10_PERCENT",
    "WIDTH_OF_20_PERCENT",
    "WIDTH_OF_30_PERCENT",
    "WIDTH_OF_40_PERCENT",
    "WIDTH_OF_50_PERCENT",
    "WIDTH_OF_60_PERCENT",
    "WIDTH_OF_70_PERCENT",
    "WIDTH_OF_80_PERCENT",
    "WIDTH_OF_90_PERCENT",
    "WIDTH_OF_95_PERCENT",
    "WIDTH_OF_99_PERCENT",
    "SAMPLE",
    "LIBRARY",
    "READ_GROUP",
];

/// One `MetricsFile` row for one pair orientation.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    /// Histogram column prefix of the row's accumulation unit.
    pub key: String,
    /// The `SAMPLE`, `LIBRARY` and `READ_GROUP` columns.
    pub labels: Labels,
    pub orientation: Orientation,
    pub read_pairs: i64,
    pub min_insert_size: i64,
    pub max_insert_size: i64,
    pub median_insert_size: f64,
    pub mode_insert_size: i64,
    pub median_absolute_deviation: f64,
    /// Mean of the trimmed histogram; NaN when it has no items.
    pub mean_insert_size: f64,
    /// Sample SD of the trimmed histogram; NaN for a single pair.
    pub standard_deviation: f64,
    /// `WIDTH_OF_X_PERCENT` for each X in `WIDTH_PERCENTS`.
    pub widths: [i64; 11],
    /// The trimmed histogram, as written to the `## HISTOGRAM` section.
    pub histogram: Histogram,
}

impl Row {
    /// `WIDTH_OF_X_PERCENT` for one X of `WIDTH_PERCENTS`.
    pub fn width_of_percent(&self, percent: u32) -> Option<i64> {
        WIDTH_PERCENTS
            .iter()
            .position(|&p| p == percent)
            .map(|i| self.widths[i])
    }
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.InsertSizeMetrics";
    const COLUMNS: &'static [&'static str] = &COLUMNS;

    /// `SAMPLE`, `LIBRARY` and `READ_GROUP` are null at the `ALL_READS` level.
    fn values(&self) -> Vec<Value> {
        let mut values = vec![
            self.median_insert_size.into(),
            self.mode_insert_size.into(),
            self.median_absolute_deviation.into(),
            self.min_insert_size.into(),
            self.max_insert_size.into(),
            self.mean_insert_size.into(),
            self.standard_deviation.into(),
            self.read_pairs.into(),
            self.orientation.as_str().into(),
        ];
        values.extend(self.widths.iter().map(|&w| Value::Int(w)));
        values.extend([
            self.labels.sample.clone().into(),
            self.labels.library.clone().into(),
            self.labels.read_group.clone().into(),
        ]);
        values
    }
}

/// The metrics file content for `rows`. Histogram columns are named
/// `<key>.<orientation>_count`, with `All_Reads` as the `ALL_READS` key.
/// Picard 3.5.0 leaves an empty orientation's histogram out of the table
/// (`MetricsFile.printHistogram`) while still writing its all-zero metrics row.
pub fn output(rows: &[Row]) -> MetricsOutput {
    let series = rows
        .iter()
        .map(|r| HistogramSeries {
            label: format!("{}.{}_count", r.key, r.orientation.as_str().to_lowercase()),
            histogram: r.histogram.clone(),
        })
        .collect();
    MetricsOutput::from_rows(
        rows,
        Some(HistogramSection {
            key_label: "insert_size",
            series,
        }),
    )
}

/// Write rows as a Picard `InsertSizeMetrics` file. `headers` become the
/// `StringHeader` lines.
pub fn write_metrics(rows: &[Row], headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
    output(rows).write(headers, out)
}

/// Whether `CollectInsertSizeMetrics` counts this record (docs/01 section 2.2).
fn is_counted(record: &bam::Record, opts: &Options) -> bool {
    let f = record.flags();
    f.is_segmented()
        && !f.is_unmapped()
        && !f.is_mate_unmapped()
        && !f.is_first_segment()
        && !f.is_secondary()
        && !f.is_supplementary()
        && (opts.include_duplicates || !f.is_duplicate())
        && record.template_length() != 0
}

/// `CollectInsertSizeMetrics` as a pass [`Accumulator`]: one histogram
/// per pair orientation.
pub struct InsertSize {
    opts: Options,
    levels: Levels<Orientations>,
}

/// The three orientation histograms of one accumulation unit.
#[derive(Clone, Default)]
struct Orientations {
    fr: Histogram,
    rf: Histogram,
    tandem: Histogram,
}

impl Orientations {
    fn histogram(&self, orientation: Orientation) -> &Histogram {
        match orientation {
            Orientation::Fr => &self.fr,
            Orientation::Rf => &self.rf,
            Orientation::Tandem => &self.tandem,
        }
    }

    fn histogram_mut(&mut self, orientation: Orientation) -> &mut Histogram {
        match orientation {
            Orientation::Fr => &mut self.fr,
            Orientation::Rf => &mut self.rf,
            Orientation::Tandem => &mut self.tandem,
        }
    }
}

impl InsertSize {
    /// An accumulator with one unit per level of `opts.levels`, units
    /// taken from the header's read groups.
    pub fn new(opts: &Options, header: &sam::Header) -> Self {
        Self {
            opts: opts.clone(),
            levels: Levels::new(header, &opts.levels, Orientations::default()),
        }
    }
}

impl Accumulator for InsertSize {
    /// The rows Picard writes; empty when every orientation is dropped, in
    /// which case Picard writes no file.
    type Output = Vec<Row>;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if !is_counted(record, &self.opts) {
            return Ok(());
        }
        let (Some(Ok(reference)), Some(Ok(mate_reference))) = (
            record.reference_sequence_id(),
            record.mate_reference_sequence_id(),
        ) else {
            return Ok(());
        };
        if reference != mate_reference {
            return Ok(());
        }
        let insert_size = i64::from(record.template_length()).abs();
        let orientation = orientation(record)?;
        self.levels.route(record, |unit| {
            unit.histogram_mut(orientation).increment(insert_size);
        });
        Ok(())
    }

    fn finish(self) -> Result<Vec<Row>> {
        let opts = &self.opts;
        let mut rows = Vec::new();
        for bucket in self.levels.into_buckets() {
            let unit = &bucket.state;
            let total = unit.fr.count() + unit.rf.count() + unit.tandem.count();
            for orientation in [Orientation::Fr, Orientation::Rf, Orientation::Tandem] {
                let h = unit.histogram(orientation);
                if passes_minimum_pct(h.count(), total, opts.minimum_pct) {
                    let mut row = finalize_row(orientation, h, opts);
                    row.key = bucket.key.clone();
                    row.labels = bucket.labels.clone();
                    rows.push(row);
                }
            }
        }
        Ok(rows)
    }
}

/// `StringHeader` command line in Picard's shape (parsers such as MultiQC
/// read `INPUT=` from it).
pub fn command_header(opts: &Options, bam: &Path, output: Option<&Path>) -> String {
    let levels: Vec<&str> = opts.levels.iter().map(|l| l.name()).collect();
    let mut command = format!(
        "CollectInsertSizeMetrics INPUT={} MINIMUM_PCT={:?} DEVIATIONS={:?} INCLUDE_DUPLICATES={} METRIC_ACCUMULATION_LEVEL=[{}]",
        bam.display(),
        opts.minimum_pct,
        opts.deviations,
        opts.include_duplicates,
        levels.join(", ")
    );
    if let Some(output) = output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(w) = opts.histogram_width {
        command.push_str(&format!(" HISTOGRAM_WIDTH={w}"));
    }
    if let Some(w) = opts.min_histogram_width {
        command.push_str(&format!(" MIN_HISTOGRAM_WIDTH={w}"));
    }
    command
}

/// Compute the metrics of one BAM with the default BGZF worker count.
pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Vec<Row>> {
    pass::run(path.as_ref(), pass::default_threads(), |header| {
        Ok(InsertSize::new(opts, header))
    })
}

/// The orientation filter in `InsertSizeMetricsCollector.PerUnitInsertSizeMetricsCollector
/// .addMetricsToFile`: keep when `count >= totalInserts * minimumPct`. The
/// collector holds `minimumPct` as a double widened from the tool's `float`
/// argument, so the default 0.05f is 0.0500000007..., and an orientation at
/// exactly 5% (2 of 40) is dropped.
fn passes_minimum_pct(count: f64, total: f64, minimum_pct: f32) -> bool {
    total > 0.0 && count >= total * f64::from(minimum_pct)
}

/// The `WIDTH_OF_X_PERCENT` loop of `PerUnitInsertSizeMetricsCollector
/// .addMetricsToFile`. `low` and `high` are doubles starting at the median,
/// bins are looked up at `(int) low` and `(int) high` (truncation, like Rust
/// `as`), and the window grows by one on each side while
/// `low >= min - 1 || high <= max + 1`. Each width is the window size
/// `(int) (high - low) + 1` the first time the covered fraction reaches X.
/// `f64::from(p) / 100.0` is the correctly rounded double, the same value as
/// Picard's literals `0.1`, `0.2`, ..., `0.99`.
fn widths(h: &Histogram, median: f64) -> [i64; 11] {
    let mut widths = [0i64; 11];
    let (Some(min), Some(max)) = (h.min(), h.max()) else {
        return widths;
    };
    let total = h.count();
    let (min, max) = (min as f64, max as f64);
    let (mut low, mut high) = (median, median);
    let mut covered = 0.0;
    while low >= min - 1.0 || high <= max + 1.0 {
        covered += h.get(low as i64).unwrap_or(0.0);
        if low != high {
            covered += h.get(high as i64).unwrap_or(0.0);
        }
        let fraction = covered / total;
        let distance = (high - low) as i64 + 1;
        for (width, percent) in widths.iter_mut().zip(WIDTH_PERCENTS) {
            if *width == 0 && fraction >= f64::from(percent) / 100.0 {
                *width = distance;
            }
        }
        low -= 1.0;
        high += 1.0;
    }
    widths
}

/// `getWidthToTrimTo`: `HISTOGRAM_WIDTH` if set, else
/// `(int) (median + deviations * MAD)` raised to `MIN_HISTOGRAM_WIDTH`.
fn trim_width(median: f64, mad: f64, opts: &Options) -> i64 {
    opts.histogram_width.unwrap_or_else(|| {
        let width = (median + opts.deviations * mad) as i64;
        opts.min_histogram_width.map_or(width, |m| m.max(width))
    })
}

fn finalize_row(orientation: Orientation, h: &Histogram, opts: &Options) -> Row {
    let (Some(min), Some(max), Some(mode)) = (h.min(), h.max(), h.mode()) else {
        // An empty orientation passes the cutoff only when MINIMUM_PCT is 0.
        // `addMetricsToFile` then skips every computation (`histogram.isEmpty()`
        // guard), so all numeric fields keep Java's default 0.
        return Row {
            key: String::new(),
            labels: Labels::default(),
            orientation,
            read_pairs: 0,
            min_insert_size: 0,
            max_insert_size: 0,
            median_insert_size: 0.0,
            mode_insert_size: 0,
            median_absolute_deviation: 0.0,
            mean_insert_size: 0.0,
            standard_deviation: 0.0,
            widths: [0; 11],
            histogram: Histogram::new(),
        };
    };
    let median = h.median();
    let mad = h.median_absolute_deviation();
    let widths = widths(h, median);

    // Picard trims the histogram in place after the widths, then takes mean
    // and SD from it and writes it as the HISTOGRAM section. Everything above
    // uses the full histogram (concordance trap #1). The trimmed histogram
    // keeps at least the median bin; a single pair gives SD 0/0 = NaN.
    let mut trimmed = h.clone();
    trimmed.trim_by_width(trim_width(median, mad, opts));
    // A width below the smallest insert size empties the histogram. Picard
    // then reports mean 0 and SD 0 (checked with W=0, 1, 100 on one pair).
    let (mean, sd) = if trimmed.min().is_none() {
        (0.0, 0.0)
    } else {
        (trimmed.mean(), trimmed.standard_deviation())
    };

    Row {
        key: String::new(),
        labels: Labels::default(),
        orientation,
        read_pairs: h.count() as i64,
        min_insert_size: min,
        max_insert_size: max,
        median_insert_size: median,
        mode_insert_size: mode,
        median_absolute_deviation: mad,
        mean_insert_size: mean,
        standard_deviation: sd,
        widths,
        histogram: trimmed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hist(items: &[(i64, f64)]) -> Histogram {
        let mut h = Histogram::new();
        for &(k, c) in items {
            h.increment_by(k, c);
        }
        h
    }

    fn width50(h: &Histogram) -> i64 {
        widths(h, h.median())[4]
    }

    #[test]
    fn width_counts_center_bin_once() {
        // 10 (x2) is the median bin: 2/4 = 0.5 reached by the first step, width 1.
        let h = hist(&[(9, 1.0), (10, 2.0), (11, 1.0)]);
        assert_eq!(width50(&h), 1);
    }

    #[test]
    fn width_grows_symmetrically() {
        let h = hist(&[(9, 1.0), (10, 1.0), (11, 1.0), (12, 1.0), (30, 1.0)]);
        // median 11: 1/5, then +10 and +12 -> 3/5 at width 3
        assert_eq!(width50(&h), 3);
    }

    #[test]
    fn width_half_integer_median_truncates_to_lower_bin() {
        // median 11.5 -> start bin 11; 10 and 12 are covered at the second step.
        let h = hist(&[(10, 1.0), (11, 1.0), (12, 1.0), (13, 1.0)]);
        assert_eq!(h.median(), 11.5);
        assert_eq!(width50(&h), 3);
    }

    #[test]
    fn all_widths_fill_in_order() {
        // 100 items at keys 1..=100, median 50.5. Bins are read at (int) low
        // and (int) high, so the window is 50, then 49..=51, ...: width 2k + 1
        // covers 2k + 1 items, the first odd width at or above X.
        let h = hist(&(1..=100).map(|k| (k, 1.0)).collect::<Vec<_>>());
        assert_eq!(
            widths(&h, h.median()),
            [11, 21, 31, 41, 51, 61, 71, 81, 91, 95, 99]
        );
    }

    #[test]
    fn empty_orientation_writes_zero_row_and_no_histogram_column() {
        let opts = Options {
            minimum_pct: 0.0,
            ..Options::default()
        };
        let mut rows = vec![
            finalize_row(Orientation::Fr, &hist(&[(250, 1.0)]), &opts),
            finalize_row(Orientation::Rf, &Histogram::new(), &opts),
        ];
        for row in &mut rows {
            row.key = "All_Reads".to_owned();
        }
        assert_eq!(rows[1].read_pairs, 0);
        assert_eq!(rows[1].mean_insert_size, 0.0);
        assert!(rows[0].standard_deviation.is_nan());
        let mut out = Vec::new();
        write_metrics(&rows, Vec::new(), &mut out).expect("write to Vec");
        let text = String::from_utf8(out).expect("utf8");
        assert!(text.contains("\n250\t250\t0\t250\t250\t250\t?\t1\tFR\t"));
        assert!(text.contains("\n0\t0\t0\t0\t0\t0\t0\t0\tRF\t0\t"));
        assert!(text.contains("\ninsert_size\tAll_Reads.fr_count\n250\t1\n"));
    }

    #[test]
    fn trim_width_options() {
        let opts = Options::default();
        assert_eq!(trim_width(101.5, 0.5, &opts), 106);
        let floor = Options {
            min_histogram_width: Some(500),
            ..Options::default()
        };
        assert_eq!(trim_width(101.5, 0.5, &floor), 500);
        let fixed = Options {
            histogram_width: Some(50),
            min_histogram_width: Some(500),
            ..Options::default()
        };
        assert_eq!(trim_width(101.5, 0.5, &fixed), 50);
    }

    #[test]
    fn minimum_pct_uses_widened_float() {
        assert!(!passes_minimum_pct(2.0, 40.0, 0.05));
        assert!(passes_minimum_pct(2.0, 39.0, 0.05));
        assert!(passes_minimum_pct(3.0, 40.0, 0.05));
        assert!(!passes_minimum_pct(1.0, 40.0, 0.05));
        assert!(!passes_minimum_pct(0.0, 0.0, 0.0));
    }
}

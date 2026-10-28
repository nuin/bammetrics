//! `CollectInsertSizeMetrics` (see docs/01 section 2.2).

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles::bam;

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metrics_file::{HistogramColumn, HistogramTable, MetricsFile, format_double};
pub use crate::record::Orientation;
use crate::record::orientation;

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
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_pct: 0.05,
            deviations: 10.0,
            histogram_width: None,
            min_histogram_width: None,
            include_duplicates: false,
        }
    }
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

    /// The row's fields as Picard prints them, in `COLUMNS` order. `SAMPLE`,
    /// `LIBRARY` and `READ_GROUP` are empty at the `ALL_READS` level.
    fn picard_fields(&self) -> Vec<String> {
        let mut fields = vec![
            format_double(self.median_insert_size),
            self.mode_insert_size.to_string(),
            format_double(self.median_absolute_deviation),
            self.min_insert_size.to_string(),
            self.max_insert_size.to_string(),
            format_double(self.mean_insert_size),
            format_double(self.standard_deviation),
            self.read_pairs.to_string(),
            self.orientation.as_str().to_owned(),
        ];
        fields.extend(self.widths.iter().map(i64::to_string));
        fields.extend([String::new(), String::new(), String::new()]);
        fields
    }
}

/// Write rows as a Picard `InsertSizeMetrics` file. `headers` become the
/// `StringHeader` lines. Histogram columns are named
/// `All_Reads.<orientation>_count`, as Picard names the `ALL_READS` level.
/// Picard 3.5.0 leaves an empty orientation's histogram out of the table
/// (`MetricsFile.printHistogram`) while still writing its all-zero metrics row.
pub fn write_metrics(rows: &[Row], headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
    let columns = rows
        .iter()
        .map(|r| HistogramColumn {
            label: format!("All_Reads.{}_count", r.orientation.as_str().to_lowercase()),
            histogram: &r.histogram,
        })
        .collect();
    MetricsFile {
        headers,
        class: "picard.analysis.InsertSizeMetrics",
        columns: &COLUMNS,
        rows: rows.iter().map(Row::picard_fields).collect(),
        histogram: Some(HistogramTable {
            key_label: "insert_size",
            columns,
        }),
    }
    .write(out)
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

pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Vec<Row>> {
    let path = path.as_ref();
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let mut reader = bam::io::Reader::new(File::open(path).map_err(input)?);
    reader.read_header().map_err(input)?;
    let (mut fr, mut rf, mut tandem) = (Histogram::new(), Histogram::new(), Histogram::new());

    for result in reader.records() {
        let record = result?;
        if !is_counted(&record, opts) {
            continue;
        }
        let (Some(Ok(reference)), Some(Ok(mate_reference))) = (
            record.reference_sequence_id(),
            record.mate_reference_sequence_id(),
        ) else {
            continue;
        };
        if reference != mate_reference {
            continue;
        }
        let insert_size = i64::from(record.template_length()).abs();
        match orientation(&record)? {
            Orientation::Fr => fr.increment(insert_size),
            Orientation::Rf => rf.increment(insert_size),
            Orientation::Tandem => tandem.increment(insert_size),
        }
    }

    let total = fr.count() + rf.count() + tandem.count();
    let rows = [
        (Orientation::Fr, &fr),
        (Orientation::Rf, &rf),
        (Orientation::Tandem, &tandem),
    ]
    .into_iter()
    .filter(|(_, h)| passes_minimum_pct(h.count(), total, opts.minimum_pct))
    .map(|(orientation, h)| finalize_row(orientation, h, opts))
    .collect();
    Ok(rows)
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

    Row {
        orientation,
        read_pairs: h.count() as i64,
        min_insert_size: min,
        max_insert_size: max,
        median_insert_size: median,
        mode_insert_size: mode,
        median_absolute_deviation: mad,
        mean_insert_size: trimmed.mean(),
        standard_deviation: trimmed.standard_deviation(),
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
        let rows = vec![
            finalize_row(Orientation::Fr, &hist(&[(250, 1.0)]), &opts),
            finalize_row(Orientation::Rf, &Histogram::new(), &opts),
        ];
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

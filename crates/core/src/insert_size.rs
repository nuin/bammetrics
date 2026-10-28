//! `CollectInsertSizeMetrics` (see docs/01 section 2.2).

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metrics_file::{HistogramColumn, HistogramTable, MetricsFile, format_double};

/// `SamPairUtil.PairOrientation`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Orientation {
    Fr,
    Rf,
    Tandem,
}

impl Orientation {
    /// The `PAIR_ORIENTATION` column value.
    pub fn as_str(self) -> &'static str {
        match self {
            Orientation::Fr => "FR",
            Orientation::Rf => "RF",
            Orientation::Tandem => "TANDEM",
        }
    }
}

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
/// while still writing its all-zero metrics row.
pub fn write_metrics(rows: &[Row], headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
    let columns = rows
        .iter()
        .filter(|r| !r.histogram.is_empty())
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

fn reference_length(record: &bam::Record) -> io::Result<i64> {
    let mut length = 0usize;
    for op in record.cigar().iter() {
        let op = op?;
        if matches!(
            op.kind(),
            Kind::Match
                | Kind::Deletion
                | Kind::Skip
                | Kind::SequenceMatch
                | Kind::SequenceMismatch
        ) {
            length += op.len();
        }
    }
    Ok(length as i64)
}

fn record_name(record: &bam::Record) -> String {
    record
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_default()
}

fn position(
    record: &bam::Record,
    value: Option<io::Result<noodles::core::Position>>,
    what: &'static str,
) -> Result<i64> {
    match value {
        Some(Ok(p)) => Ok(p.get() as i64),
        Some(Err(e)) => Err(e.into()),
        None => Err(Error::MalformedRecord {
            name: record_name(record),
            what,
        }),
    }
}

/// `Cigar.getReferenceLength` of a SAM text CIGAR such as an `MC` value:
/// the summed lengths of M, D, N, = and X operations. `*` is an empty CIGAR
/// (`TextCigarCodec.decode`), length 0.
fn text_cigar_reference_length(cigar: &[u8]) -> Option<i64> {
    if cigar == b"*" {
        return Some(0);
    }
    let mut total = 0i64;
    let mut len: Option<i64> = None;
    for &c in cigar {
        if c.is_ascii_digit() {
            let digit = i64::from(c - b'0');
            len = Some(len.unwrap_or(0).checked_mul(10)?.checked_add(digit)?);
            continue;
        }
        let n = len.take()?;
        match c {
            b'M' | b'D' | b'N' | b'=' | b'X' => total += n,
            b'I' | b'S' | b'H' | b'P' => {}
            _ => return None,
        }
    }
    (len.is_none() && !cigar.is_empty()).then_some(total)
}

/// Reference length of the mate from the `MC` tag, `None` when the tag is absent.
fn mate_reference_length(record: &bam::Record) -> Result<Option<i64>> {
    let malformed = |what| Error::MalformedRecord {
        name: record_name(record),
        what,
    };
    match record.data().get(&Tag::MATE_CIGAR) {
        None => Ok(None),
        Some(Err(e)) => Err(e.into()),
        Some(Ok(Value::String(cigar))) => text_cigar_reference_length(cigar)
            .map(Some)
            .ok_or_else(|| malformed("unparseable MC tag")),
        Some(Ok(_)) => Err(malformed("MC tag is not a string")),
    }
}

/// 5' position of the reverse-strand mate, seen from the forward-strand read.
/// Mirrors the `else` branch of `SamPairUtil.getPairOrientation`: with `MC`,
/// `CoordMath.getEnd(mateStart, mateCigar.getReferenceLength())`; without it,
/// `CoordMath.getEnd(alignmentStart, inferredInsertSize)`. `getEnd` is
/// `start + length - 1`. The two disagree for dovetailed pairs, where the
/// reverse mate starts left of the forward read.
fn mate_five_prime(
    start: i64,
    template_length: i64,
    mate_start: i64,
    mate_reference_length: Option<i64>,
) -> i64 {
    match mate_reference_length {
        Some(length) => mate_start + length - 1,
        None => start + template_length - 1,
    }
}

/// `SamPairUtil.getPairOrientation` after the strand test: FR when the
/// forward read's 5' end is at or left of the reverse read's 5' end.
fn orient(forward_5p: i64, reverse_5p: i64) -> Orientation {
    if forward_5p <= reverse_5p {
        Orientation::Fr
    } else {
        Orientation::Rf
    }
}

/// `SamPairUtil.getPairOrientation` (htsjdk 5.0.0). Same strand is TANDEM.
/// The forward 5' end is this read's start, or the mate's start when this read
/// is reverse. The reverse 5' end is this read's alignment end when it is
/// reverse, or `mate_five_prime` when it is forward.
fn orientation(record: &bam::Record) -> Result<Orientation> {
    let flags = record.flags();
    if flags.is_reverse_complemented() == flags.is_mate_reverse_complemented() {
        return Ok(Orientation::Tandem);
    }
    let start = position(record, record.alignment_start(), "missing alignment start")?;
    let mate_start = position(
        record,
        record.mate_alignment_start(),
        "missing mate alignment start",
    )?;
    Ok(if flags.is_reverse_complemented() {
        orient(mate_start, start + reference_length(record)? - 1)
    } else {
        let reverse_5p = mate_five_prime(
            start,
            i64::from(record.template_length()),
            mate_start,
            mate_reference_length(record)?,
        );
        orient(start, reverse_5p)
    })
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

    #[test]
    fn mc_cigar_reference_length() {
        assert_eq!(text_cigar_reference_length(b"100M"), Some(100));
        assert_eq!(text_cigar_reference_length(b"5S90M2I3D1N4=2X5H"), Some(100));
        assert_eq!(text_cigar_reference_length(b"*"), Some(0));
        assert_eq!(text_cigar_reference_length(b""), None);
        assert_eq!(text_cigar_reference_length(b"M"), None);
        assert_eq!(text_cigar_reference_length(b"10"), None);
        assert_eq!(text_cigar_reference_length(b"10Q"), None);
    }

    #[test]
    fn mate_end_uses_mc_when_present() {
        // Forward read at 1000, reverse mate at 1200 with MC 100M: end 1299.
        assert_eq!(mate_five_prime(1000, 300, 1200, Some(100)), 1299);
        // MC wins over a TLEN that disagrees with it.
        assert_eq!(mate_five_prime(1000, 5, 1200, Some(100)), 1299);
    }

    #[test]
    fn mate_end_falls_back_to_tlen_minus_one() {
        // CoordMath.getEnd(1000, 300) = 1299, not 1300.
        assert_eq!(mate_five_prime(1000, 300, 1200, None), 1299);
        assert_eq!(mate_five_prime(1300, -400, 1000, None), 899);
    }

    #[test]
    fn dovetail_flips_between_mc_and_tlen() {
        // Forward read 1000..1099; reverse mate starts at 900 with MC 150M (ends 1049).
        // MC: 1000 <= 1049, FR. TLEN -200 fallback: 1000 > 799, RF.
        assert_eq!(
            orient(1000, mate_five_prime(1000, -200, 900, Some(150))),
            Orientation::Fr
        );
        assert_eq!(
            orient(1000, mate_five_prime(1000, -200, 900, None)),
            Orientation::Rf
        );
    }

    #[test]
    fn orient_boundary_is_fr() {
        assert_eq!(orient(1000, 1000), Orientation::Fr);
        assert_eq!(orient(1001, 1000), Orientation::Rf);
    }
}

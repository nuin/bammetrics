//! `CollectInsertSizeMetrics` (see docs/01 section 2.2).

use std::fs::File;
use std::io;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};
use crate::histogram::Histogram;

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
    pub include_duplicates: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_pct: 0.05,
            deviations: 10.0,
            include_duplicates: false,
        }
    }
}

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
    pub mean_insert_size: f64,
    pub standard_deviation: f64,
    pub width_of_50_percent: i64,
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
    .filter(|(_, h)| !h.is_empty() && passes_minimum_pct(h.count(), total, opts.minimum_pct))
    .filter_map(|(orientation, h)| finalize_row(orientation, h, opts))
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

/// `WIDTH_OF_X_PERCENT`: grow a symmetric window from the median bin and
/// return its width the first time it covers `fraction` of the pairs.
/// Picard truncates the median to a bin and counts that bin once.
fn width_of_fraction(h: &Histogram, median: f64, fraction: f64) -> i64 {
    let (Some(min), Some(max)) = (h.min(), h.max()) else {
        return 0;
    };
    let total = h.count();
    let (mut low, mut high) = (median as i64, median as i64);
    let mut covered = 0.0;
    while low >= min - 1 || high <= max + 1 {
        covered += h.get(low).unwrap_or(0.0);
        if low != high {
            covered += h.get(high).unwrap_or(0.0);
        }
        if covered / total >= fraction {
            return high - low + 1;
        }
        low -= 1;
        high += 1;
    }
    0
}

fn finalize_row(orientation: Orientation, h: &Histogram, opts: &Options) -> Option<Row> {
    let median = h.median();
    let mad = h.median_absolute_deviation();

    // Mean and SD come from the histogram trimmed to median + deviations * MAD;
    // everything else uses the full histogram (concordance trap #1).
    let mut trimmed = h.clone();
    trimmed.trim_by_width((median + opts.deviations * mad) as i64);
    let (mean, sd) = if trimmed.is_empty() {
        (0.0, 0.0)
    } else {
        (trimmed.mean(), trimmed.standard_deviation())
    };

    Some(Row {
        orientation,
        read_pairs: h.count() as i64,
        min_insert_size: h.min()?,
        max_insert_size: h.max()?,
        median_insert_size: median,
        mode_insert_size: h.mode()?,
        median_absolute_deviation: mad,
        mean_insert_size: mean,
        standard_deviation: sd,
        width_of_50_percent: width_of_fraction(h, median, 0.5),
    })
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

    #[test]
    fn width_counts_center_bin_once() {
        // 10 (x2) is the median bin: 2/4 = 0.5 reached by the first step, width 1.
        let h = hist(&[(9, 1.0), (10, 2.0), (11, 1.0)]);
        assert_eq!(width_of_fraction(&h, h.median(), 0.5), 1);
    }

    #[test]
    fn width_grows_symmetrically() {
        let h = hist(&[(9, 1.0), (10, 1.0), (11, 1.0), (12, 1.0), (30, 1.0)]);
        // median 11: 1/5, then +10 and +12 -> 3/5 at width 3
        assert_eq!(width_of_fraction(&h, h.median(), 0.5), 3);
    }

    #[test]
    fn width_half_integer_median_truncates_to_lower_bin() {
        // median 10.5 -> start bin 10; 10 and 11 are covered at the second step.
        let h = hist(&[(10, 1.0), (11, 1.0), (12, 1.0), (13, 1.0)]);
        assert_eq!(h.median(), 11.5);
        assert_eq!(width_of_fraction(&h, h.median(), 0.5), 3);
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

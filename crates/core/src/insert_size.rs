//! `CollectInsertSizeMetrics` (see docs/01 section 2.2).

use std::fs::File;
use std::io;
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;

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
    pub minimum_pct: f64,
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

/// `SamPairUtil.getPairOrientation`. The 5' position of the forward read is
/// its start, or the mate's start when this read is the reverse one. The 5'
/// position of the reverse read is its alignment end, or start + TLEN when
/// this read is the forward one (htsjdk uses the inferred insert size, not MC).
fn orientation(record: &bam::Record) -> Result<Orientation> {
    let flags = record.flags();
    if flags.is_reverse_complemented() == flags.is_mate_reverse_complemented() {
        return Ok(Orientation::Tandem);
    }
    let start = position(record, record.alignment_start(), "missing alignment start")?;
    let (forward_5p, reverse_5p) = if flags.is_reverse_complemented() {
        let mate_start = position(
            record,
            record.mate_alignment_start(),
            "missing mate alignment start",
        )?;
        let end = start + reference_length(record)? - 1;
        (mate_start, end)
    } else {
        (start, start + i64::from(record.template_length()))
    };
    Ok(if forward_5p <= reverse_5p {
        Orientation::Fr
    } else {
        Orientation::Rf
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
    .filter(|(_, h)| total > 0.0 && h.count() >= total * opts.minimum_pct && !h.is_empty())
    .filter_map(|(orientation, h)| finalize_row(orientation, h, opts))
    .collect();
    Ok(rows)
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
}

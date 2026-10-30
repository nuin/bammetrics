//! BAM record helpers shared by the metrics: names, positions, and pair
//! orientation as htsjdk `SamPairUtil.getPairOrientation` computes it.

use std::io;

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};

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

impl std::str::FromStr for Orientation {
    type Err = Error;

    /// `PairOrientation.valueOf`: FR, RF or TANDEM, upper case.
    fn from_str(value: &str) -> Result<Self> {
        match value {
            "FR" => Ok(Orientation::Fr),
            "RF" => Ok(Orientation::Rf),
            "TANDEM" => Ok(Orientation::Tandem),
            _ => Err(Error::Usage(format!(
                "pair orientation {value:?}: expected FR, RF or TANDEM"
            ))),
        }
    }
}

pub(crate) fn reference_length(record: &bam::Record) -> io::Result<i64> {
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

pub(crate) fn record_name(record: &bam::Record) -> String {
    record
        .name()
        .map(|n| String::from_utf8_lossy(n.as_ref()).into_owned())
        .unwrap_or_default()
}

pub(crate) fn position(
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

/// 1-based last aligned reference base of the mate, from the mate's start and
/// the `MC` tag (`SAMRecord.getMateAlignmentEnd` semantics). `None` when the
/// record has no mapped mate position or no `MC` tag.
pub fn mate_alignment_end(record: &bam::Record) -> Result<Option<i64>> {
    let Some(length) = mate_reference_length(record)? else {
        return Ok(None);
    };
    match record.mate_alignment_start() {
        Some(start) => Ok(Some(start?.get() as i64 + length - 1)),
        None => Ok(None),
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
pub(crate) fn orientation(record: &bam::Record) -> Result<Orientation> {
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

#[cfg(test)]
mod tests {
    use super::*;

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

//! `CleanSam`, mirroring Picard 3.5.0 (`picard.sam.CleanSam`): every record
//! of the input, with two fixes, through a presorted `SAMFileWriter`
//! ([`Mode::Rewrite`]).
//!
//! 1. A mapped record whose alignment ends past the end of its reference
//!    is soft-clipped at its end, as `AbstractAlignmentMerger.
//!    createNewCigarsIfMapsOffEndOfReference` clips it: `overhang` is the
//!    alignment end less the reference length, `clipFrom` is the read
//!    length (the bases of `SEQ`, 0 for `*`) less `overhang` plus 1, less
//!    the length of a trailing `S`; then htsjdk 5.0.0 `CigarUtil.
//!    softClipEndOfRead` ([`soft_clip_end_of_read`]). The clip is at the
//!    end of the CIGAR whatever the strand, and the start stays. `overhang`
//!    counts reference bases and `clipFrom` read bases, so a deletion or
//!    an insertion before the clip point moves it, as in Picard. A trailing
//!    `H` stays after the new `S`; an `I` or `S` at the clip point joins the
//!    new `S`; a `D`, `N`, `H` or `P` at a clip point inside it fails
//!    (`Unexpected non-0 relativeClippedPosition`). A read that starts past
//!    the end becomes one `S` longer than its bases.
//! 2. The same clip on the mate's CIGAR in `MC` of a paired record whose
//!    mate is mapped, with the mate's reference, `PNEXT` and the read
//!    length of the `MC` CIGAR. `MC` gets the clipped CIGAR as text.
//! 3. An unmapped record gets `MAPQ` 0.
//!
//! A record whose CIGAR or `MC` changes is encoded again as htsjdk encodes
//! a record whose variable fields changed: tags in [`crate::bam_raw::Tags`]
//! order and widths. A record with only a new `MAPQ` keeps its tag bytes.
//! Every record gets the `bin` htsjdk's writer computes. `NM` and `MD` are
//! left as they were. No `@PG` line is added. `VALIDATION_STRINGENCY=
//! STRICT` (the default) reads at `LENIENT`, as Picard does, so a SAM line
//! past the end of its reference is read.
//!
//! Where Picard writes a record no reader opens, bammetrics fails and
//! writes nothing (docs/07 decision rows): a clip that leaves a CIGAR whose
//! read length differs from the bases of `SEQ` (a read that starts past
//! the end of its reference; samtools refuses Picard's file). A `SEQ` of
//! `*` takes any CIGAR, so such a record is written as Picard writes it.
//!
//! In a chain the step cleans the stream at its position; the header is
//! the stream's.

use noodles::sam;

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Setup};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text::{Stringency, parse_cigar};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub const TOOL: Tool = Tool {
    name: "CleanSam",
    summary: "no rows; OUTPUT gets every record of the input, a mapped record past the end \
              of its reference soft-clipped at its end, the MC of a mate past the end clipped \
              the same way, an unmapped record at MAPQ 0: a BAM, or SAM text for a name \
              ending in .sam.",
    options: &[
        record_edit::OUTPUT_OPTION,
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        record_edit::RECORDS_IN_RAM_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// BAM CIGAR operator codes.
const SOFT_CLIP: u8 = 4;
const HARD_CLIP: u8 = 5;

/// htsjdk `CigarOperator.consumesReadBases`: `M`, `I`, `S`, `=`, `X`.
fn consumes_read(op: u8) -> bool {
    matches!(op, 0 | 1 | 4 | 7 | 8)
}

/// htsjdk `CigarOperator.consumesReferenceBases`: `M`, `D`, `N`, `=`, `X`.
fn consumes_reference(op: u8) -> bool {
    matches!(op, 0 | 2 | 3 | 7 | 8)
}

/// A CIGAR length as a Java `int`.
fn int(len: u32) -> i32 {
    len as i32
}

/// `Cigar.getReadLength`, summed as a Java `int`.
fn read_length(cigar: &[(u32, u8)]) -> i32 {
    cigar
        .iter()
        .filter(|(_, op)| consumes_read(*op))
        .fold(0i32, |sum, (len, _)| sum.wrapping_add(int(*len)))
}

/// `Cigar.getReferenceLength`, summed as a Java `int`.
fn reference_length(cigar: &[(u32, u8)]) -> i32 {
    cigar
        .iter()
        .filter(|(_, op)| consumes_reference(*op))
        .fold(0i32, |sum, (len, _)| sum.wrapping_add(int(*len)))
}

/// `new CigarElement(length, op)`: a negative length is an error.
fn element(len: i32, op: u8) -> std::result::Result<(u32, u8), String> {
    u32::try_from(len).map(|len| (len, op)).map_err(|_| {
        format!(
            "Cigar element being constructed with negative length: {len} and operation: {}",
            op_name(op)
        )
    })
}

/// htsjdk `CigarOperator.name()`.
fn op_name(op: u8) -> &'static str {
    match op {
        0 => "M",
        1 => "I",
        2 => "D",
        3 => "N",
        4 => "S",
        5 => "H",
        6 => "P",
        7 => "EQ",
        _ => "X",
    }
}

/// htsjdk 5.0.0 `CigarUtil.softClipEndOfRead(clipFrom, elements)`: the
/// elements before the one holding read position `clip_from` (1-based),
/// then that element merged with a soft clip of every base from there to
/// the end, then a trailing `H` when the CIGAR ended with one.
pub fn soft_clip_end_of_read(
    clip_from: i32,
    cigar: &[(u32, u8)],
) -> std::result::Result<Vec<(u32, u8)>, String> {
    let clipped = read_length(cigar).wrapping_sub(clip_from).wrapping_add(1);
    let trailing_hard = match cigar.last() {
        Some(&(len, HARD_CLIP)) => int(len),
        Some(_) => 0,
        None => return Err("Index -1 out of bounds for length 0".to_owned()),
    };
    let mut out = Vec::with_capacity(cigar.len() + 2);
    let mut pos = 1i32;
    for &(len, op) in cigar {
        let bases = if consumes_read(op) { int(len) } else { 0 };
        let end = pos.wrapping_add(bases).wrapping_sub(1);
        if end < clip_from.wrapping_sub(1) {
            out.push((len, op));
        } else {
            let relative = clip_from.wrapping_sub(1).wrapping_sub(pos.wrapping_sub(1));
            merge_clip(&mut out, op, relative, clipped, trailing_hard)?;
            break;
        }
        pos = end.wrapping_add(1);
    }
    Ok(out)
}

/// htsjdk `CigarUtil.mergeClippingCigarElement` for a soft clip.
fn merge_clip(
    out: &mut Vec<(u32, u8)>,
    op: u8,
    relative: i32,
    clipped: i32,
    trailing_hard: i32,
) -> std::result::Result<(), String> {
    let mut clip = clipped;
    if consumes_read(op) {
        if consumes_reference(op) {
            if relative > 0 {
                out.push(element(relative, op)?);
            }
        } else {
            clip = clipped.wrapping_add(relative);
        }
    } else if relative != 0 {
        return Err(format!(
            "Unexpected non-0 relativeClippedPosition {relative}"
        ));
    }
    out.push(element(clip, SOFT_CLIP)?);
    if trailing_hard > 0 {
        out.push(element(trailing_hard, HARD_CLIP)?);
    }
    Ok(())
}

/// `createNewCigarIfMapsOffEndOfReference` for a mapped alignment ending at
/// `end` on a reference of `reference_length` bases: the clipped CIGAR, or
/// `None` when nothing hangs off.
pub fn clip_off_end(
    reference_length: i32,
    end: i32,
    read_length: i32,
    cigar: &[(u32, u8)],
) -> std::result::Result<Option<Vec<(u32, u8)>>, String> {
    let overhang = end.wrapping_sub(reference_length);
    if overhang <= 0 {
        return Ok(None);
    }
    let mut clip_from = read_length.wrapping_sub(overhang).wrapping_add(1);
    match cigar.last() {
        Some(&(len, SOFT_CLIP)) => clip_from = clip_from.wrapping_sub(int(len)),
        Some(_) => {}
        None => return Err("Index -1 out of bounds for length 0".to_owned()),
    }
    soft_clip_end_of_read(clip_from, cigar).map(Some)
}

/// htsjdk `BAMRecord.isSentinelCigar`: `kSmN` with `k` the read length,
/// the CIGAR a BAM holds when the real one is in the `CG` tag.
fn is_sentinel(cigar: &[(u32, u8)], bases: usize) -> bool {
    matches!(cigar, [(k, SOFT_CLIP), (_, 3)] if bases == 0 || usize::try_from(*k).ok() == Some(bases))
}

/// htsjdk `TextCigarCodec.encode`: `*` for no operator.
pub fn cigar_text(cigar: &[(u32, u8)]) -> Vec<u8> {
    if cigar.is_empty() {
        return b"*".to_vec();
    }
    let mut text = Vec::new();
    for &(len, op) in cigar {
        text.extend_from_slice(len.to_string().as_bytes());
        text.push(b"MIDNSHP=X"[usize::from(op)]);
    }
    text
}

/// The edit: the reference lengths of the header, as Java `int`s, and the
/// tool that names its errors. `RevertOriginalBaseQualitiesAndAddMateCigar`
/// makes the same clip ([`Clean::clip`]).
pub(crate) struct Clean {
    lengths: Vec<i32>,
    tool: &'static str,
}

impl Clean {
    pub(crate) fn new(header: &sam::Header, tool: &'static str) -> Self {
        let lengths = header
            .reference_sequences()
            .values()
            .map(|sequence| i32::try_from(sequence.length().get()).unwrap_or(i32::MAX))
            .collect();
        Self { lengths, tool }
    }

    /// `createNewCigarsIfMapsOffEndOfReference`: the record's CIGAR and its
    /// `MC`, clipped where they end past their reference. True when either
    /// changed; the tags are then encoded again in htsjdk's order.
    pub(crate) fn clip(&self, record: &mut RawRecord) -> Result<bool> {
        let cigar = if record.flags() & 0x4 != 0 {
            None
        } else {
            self.own_clip(record)?
        };
        let mate = self.mate_clip(record)?;
        if cigar.is_none() && mate.is_none() {
            return Ok(false);
        }
        let mut tags = record.tags().map_err(Error::Record)?;
        if let Some(new) = &cigar {
            record.set_cigar(new).map_err(Error::Record)?;
        }
        if let Some(text) = mate {
            tags.set(*b"MC", TagValue::String(text));
        }
        record.set_tags(&tags).map_err(Error::Record)?;
        Ok(true)
    }

    fn failure(&self, record: &RawRecord, what: &str) -> Error {
        failure(self.tool, record, what)
    }

    /// `header.getSequence(index).getSequenceLength()`: a missing sequence
    /// is Picard's `NullPointerException`.
    fn length(&self, record: &RawRecord, index: i32, what: &str) -> Result<i32> {
        usize::try_from(index)
            .ok()
            .and_then(|at| self.lengths.get(at).copied())
            .ok_or_else(|| {
                self.failure(
                    record,
                    &format!(
                        "{what} reference index {index} is not in the sequence dictionary of {} \
                         sequences",
                        self.lengths.len()
                    ),
                )
            })
    }

    /// The record's own CIGAR, clipped when its alignment ends past its
    /// reference.
    fn own_clip(&self, record: &RawRecord) -> Result<Option<Vec<(u32, u8)>>> {
        let length = self.length(record, record.reference_id(), "the record's")?;
        let cigar: Vec<(u32, u8)> = record.cigar().collect();
        if is_sentinel(&cigar, record.sequence_length())
            && record.tag(*b"CG").ok().flatten().is_some()
        {
            return Err(self.failure(
                record,
                "its CIGAR of more than 65535 operators is in a CG tag; htsjdk clips the CIGAR \
                 of the tag, which bammetrics does not read",
            ));
        }
        let end = record
            .alignment_start()
            .wrapping_add(reference_length(&cigar))
            .wrapping_sub(1);
        let read = i32::try_from(record.sequence_length()).unwrap_or(i32::MAX);
        let clipped =
            clip_off_end(length, end, read, &cigar).map_err(|e| self.failure(record, &e))?;
        if let Some(new) = &clipped {
            let bases = record.sequence_length();
            if bases != 0 && usize::try_from(read_length(new)).ok() != Some(bases) {
                return Err(self.failure(
                    record,
                    &format!(
                        "the clip past the end of its reference gives CIGAR {}, {} read bases, \
                         for a read of {bases} bases; Picard writes this record and no reader \
                         opens the file",
                        String::from_utf8_lossy(&cigar_text(new)),
                        read_length(new)
                    ),
                ));
            }
        }
        Ok(clipped)
    }

    /// The mate's `MC`, clipped when the mate's alignment ends past its
    /// reference (`SAMUtils.hasMateCigar`, `getMateAlignmentEnd`).
    fn mate_clip(&self, record: &RawRecord) -> Result<Option<Vec<u8>>> {
        let flags = record.flags();
        if flags & 0x1 == 0 || flags & 0x8 != 0 {
            return Ok(None);
        }
        let Some(value) = record.tag(*b"MC").map_err(Error::Record)? else {
            return Ok(None);
        };
        let TagValue::String(text) = value else {
            return Err(self.failure(record, "Value for tag MC is not a String"));
        };
        let cigar = parse_cigar(&text).map_err(|e| self.failure(record, &e))?;
        if let Some(&(len, op)) = cigar.iter().find(|(len, _)| int(*len) < 0) {
            return Err(self.failure(record, &element(int(len), op).unwrap_err()));
        }
        let length = self.length(record, record.mate_reference_id(), "the mate's")?;
        let end = record
            .mate_alignment_start()
            .wrapping_add(reference_length(&cigar))
            .wrapping_sub(1);
        let clipped = clip_off_end(length, end, read_length(&cigar), &cigar)
            .map_err(|e| self.failure(record, &e))?;
        Ok(clipped.map(|new| cigar_text(&new)))
    }
}

/// The error of `tool` for `record`.
fn failure(tool: &str, record: &RawRecord, what: &str) -> Error {
    Error::Metric(format!(
        "{tool}: record {}: {what}",
        String::from_utf8_lossy(record.name())
    ))
}

impl RecordEdit for Clean {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(input.clone())
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        let unmapped = record.flags() & 0x4 != 0;
        let mut outcome = Outcome::Same;
        if self.clip(record)? {
            outcome = Outcome::Changed;
        }
        if unmapped && record.mapping_quality() != 0 {
            record.set_mapping_quality(0);
            outcome = Outcome::Changed;
        }
        record.set_codec_bin();
        Ok(outcome)
    }
}

fn build(setup: &Setup) -> Result<Box<dyn RecordEdit>> {
    Ok(Box::new(Clean::new(setup.header, TOOL.name)))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut opts = record_edit::Options::from_invocation(inv)?;
    if opts.stringency == Stringency::Strict {
        opts.stringency = Stringency::Lenient;
    }
    record_edit::run(&opts, build)?;
    Ok(ToolOutput::NONE)
}

/// `CleanSam` as a step of a chain.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    record_edit::stage(inv, frame, build)
}

#[cfg(test)]
mod tests;

//! `SetNmMdAndUqTags` and its deprecated alias `SetNmAndUqTags`, mirroring
//! Picard 3.5.0 (`picard.sam.SetNmMdAndUqTags`; `SetNmAndUqTags` extends it
//! with nothing but its usage text): each mapped record of a
//! coordinate-sorted SAM, BAM or CRAM gets `MD`, `NM` and `UQ` computed
//! against `REFERENCE_SEQUENCE`, every other record is written as read.
//!
//! The header must say `SO:coordinate` (`Input must be coordinate-sorted for
//! this program to run. Found: <order>`). The reference is walked forward as
//! `ReferenceSequenceFileWalker` walks it, by each record's reference index,
//! so a mapped record on an earlier contig than the one before it fails. A
//! record with the unmapped flag is not touched, placed or not.
//!
//! The tags, as htsjdk 5.0.0 `SequenceUtil` computes them:
//! * `MD` and `NM` ([`md_and_nm`], `calculateMdAndNm`): aligned bases (`M`,
//!   `=`, `X`) compare by IUPAC code, ignoring case
//!   ([`crate::sequence::bases_equal`]), and a mismatch writes the reference
//!   base as the FASTA holds it, lower case kept. An insertion adds its
//!   length to `NM`, a deletion writes `^` and its reference bases and adds
//!   its length. At the end of the contig the walk stops: an aligned run
//!   stops where the contig does, a deletion keeps the bases the contig has
//!   and adds nothing to `NM`, and `MD` ends with the match count. The
//!   CIGAR operators `=` and `X` say nothing: the bases decide.
//! * `UQ` ([`sum_qualities_of_mismatches`]): the sum of the base qualities
//!   of the aligned bases that differ from the reference. A record without
//!   qualities (`*`, or a BAM's `0xFF` first quality) gets no `UQ`; when it
//!   gets `MD` and `NM`, every quality byte is written `0xFF`, as htsjdk
//!   encodes a record it decoded without qualities. An
//!   aligned base past the end of the contig fails the run, as Picard's
//!   `ArrayIndexOutOfBoundsException` does; at `STRICT` Picard stops such a
//!   record before, at validation.
//! * `IS_BISULFITE_SEQUENCE`: a read `T` on a reference `C` (forward
//!   strand) or a read `A` on a reference `G` (reverse strand) is no
//!   mismatch for `NM` ([`bisulfite_nm`], `calculateSamNmTag`: the
//!   mismatches over the aligned blocks, then the length of every insertion
//!   and deletion) and `UQ`. The strand flag decides; the read end does
//!   not. `MD` is the same as without the option.
//! * `SET_ONLY_UQ`: `UQ` alone; `MD` and `NM` stay as read.
//!
//! Setting a tag decodes the record's tags, so a record that gets one has
//! its tags in htsjdk's order with each integer in its smallest type
//! ([`Tags`]); a record that gets none keeps its bytes. No `@PG` line is
//! added. The run is a [`record_edit`] rewrite: SAM text for an `OUTPUT`
//! ending in `.sam`, else a BAM, `CREATE_INDEX` as `out.bai`. In a chain the
//! step edits the stream at its position (stream name `tagged`), with the
//! chain's `-R` as its reference.

use std::path::Path;

use noodles::sam;

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::mark_duplicates::SortOrder;
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Setup};
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sequence::{bases_equal, bisulfite_bases_equal, iupac_code};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

const OPTIONS: &[Opt] = &[
    Opt::bool(
        &["IS_BISULFITE_SEQUENCE"],
        "A read T on a reference C (forward strand) or A on G (reverse strand) is no mismatch \
         for NM and UQ (default false).",
    ),
    Opt::bool(
        &["SET_ONLY_UQ"],
        "Set UQ alone; MD and NM stay as read (default false).",
    ),
    record_edit::OUTPUT_OPTION,
    Opt::path(
        &["REFERENCE_SEQUENCE", "R"],
        "Reference FASTA with a .fai index; its contigs are the input header's.",
    )
    .required(),
    record_edit::COMPRESSION_OPTION,
    record_edit::INDEX_OPTION,
    record_edit::MD5_OPTION,
    record_edit::RECORDS_IN_RAM_OPTION,
    Opt::rejected(&["ASSUME_SORTED", "AS"]),
];

pub const TOOL: Tool = Tool {
    name: "SetNmMdAndUqTags",
    summary: "no rows; OUTPUT gets every record of a coordinate-sorted input, each mapped one \
              with MD, NM and UQ computed against REFERENCE_SEQUENCE.",
    options: OPTIONS,
    run: |inv| run_tool(inv, "SetNmMdAndUqTags"),
    collector: None,
};

/// `SetNmAndUqTags`: Picard 3.5.0 keeps the old name as a subclass with the
/// same options and output; only its usage text says it is deprecated.
pub const ALIAS: Tool = Tool {
    name: "SetNmAndUqTags",
    summary: "deprecated name of SetNmMdAndUqTags, the same options and output: no rows; \
              OUTPUT gets every record, each mapped one with MD, NM and UQ.",
    options: OPTIONS,
    run: |inv| run_tool(inv, "SetNmAndUqTags"),
    collector: None,
};

const MD: [u8; 2] = *b"MD";
const NM: [u8; 2] = *b"NM";
const UQ: [u8; 2] = *b"UQ";

/// Why a tag cannot be computed: Picard throws there.
fn fails(record: &RawRecord, what: &str) -> Error {
    Error::Metric(format!(
        "record {}: {what}",
        String::from_utf8_lossy(record.name())
    ))
}

/// The base at `at`, `None` where Java's array access throws.
fn base_at(bases: &[u8], at: i64) -> Option<u8> {
    usize::try_from(at)
        .ok()
        .and_then(|at| bases.get(at))
        .copied()
}

/// htsjdk `SequenceUtil.calculateMdAndNm` with no reference offset: `MD`
/// and `NM` of the alignment `cigar` of `read` starting at the 1-based
/// `alignment_start` of `reference`. An error where htsjdk throws: a read
/// base past the read, a reference base before the contig, or a reference
/// byte its table does not hold (127 and up).
pub fn md_and_nm(
    cigar: impl Iterator<Item = (u32, u8)>,
    read: &[u8],
    reference: &[u8],
    alignment_start: i32,
) -> std::result::Result<(Vec<u8>, i64), &'static str> {
    let end = reference.len() as i64;
    let mut ref_pos = i64::from(alignment_start) - 1;
    let mut read_pos = 0i64;
    let mut matches = 0i64;
    let mut nm = 0i64;
    let mut md = Vec::new();
    let flush = |md: &mut Vec<u8>, matches: &mut i64| {
        md.extend_from_slice(matches.to_string().as_bytes());
        *matches = 0;
    };
    for (length, op) in cigar {
        let length = i64::from(length);
        match op {
            0 | 7 | 8 => {
                let mut i = 0;
                while i < length && ref_pos + i < end {
                    let read_base =
                        base_at(read, read_pos + i).ok_or("a read base past the read")?;
                    let ref_base = base_at(reference, ref_pos + i)
                        .ok_or("a reference base before the contig")?;
                    let read_code = iupac_code(read_base).ok_or("a read base out of range")?;
                    let ref_code =
                        iupac_code(ref_base).ok_or("a reference byte out of the base table")?;
                    if read_code == ref_code || read_base == 0 {
                        matches += 1;
                    } else {
                        flush(&mut md, &mut matches);
                        md.push(ref_base);
                        nm += 1;
                    }
                    i += 1;
                }
                if i < length {
                    break;
                }
                ref_pos += length;
                read_pos += length;
            }
            2 => {
                flush(&mut md, &mut matches);
                md.push(b'^');
                let mut i = 0;
                while i < length && ref_pos + i < end {
                    let base = base_at(reference, ref_pos + i)
                        .ok_or("a reference base before the contig")?;
                    if base == 0 {
                        break;
                    }
                    if base >= 0x80 {
                        return Err("a reference byte that is no code point");
                    }
                    md.push(base);
                    i += 1;
                }
                if i < length {
                    break;
                }
                ref_pos += length;
                nm += length;
            }
            1 => {
                read_pos += length;
                nm += length;
            }
            4 => read_pos += length,
            3 => ref_pos += length,
            _ => {}
        }
    }
    flush(&mut md, &mut matches);
    Ok((md, nm))
}

/// The aligned blocks of `cigar` (htsjdk `getAlignmentBlocks`): read
/// offset, reference offset from the alignment start, length; 0-based.
fn blocks(cigar: impl Iterator<Item = (u32, u8)>) -> Vec<(i64, i64, i64)> {
    let (mut read, mut reference) = (0i64, 0i64);
    let mut out = Vec::new();
    for (length, op) in cigar {
        let length = i64::from(length);
        match op {
            0 | 7 | 8 => {
                out.push((read, reference, length));
                read += length;
                reference += length;
            }
            1 | 4 => read += length,
            2 | 3 => reference += length,
            _ => {}
        }
    }
    out
}

/// The read and reference bases of each aligned position, in order; an
/// error where Java's array access throws.
fn aligned_pairs(
    cigar: impl Iterator<Item = (u32, u8)>,
    read: &[u8],
    reference: &[u8],
    alignment_start: i32,
    mut each: impl FnMut(i64, u8, u8),
) -> std::result::Result<(), &'static str> {
    let start = i64::from(alignment_start) - 1;
    for (read_at, ref_at, length) in blocks(cigar) {
        for i in 0..length {
            let read_base = base_at(read, read_at + i).ok_or("a read base past the read")?;
            let ref_base = base_at(reference, start + ref_at + i)
                .ok_or("an aligned base past the end of the reference contig")?;
            each(read_at + i, read_base, ref_base);
        }
    }
    Ok(())
}

/// The base comparison of `bisulfite` data or of any other.
fn same_base(bisulfite: bool, negative: bool, read: u8, reference: u8) -> bool {
    if bisulfite {
        bisulfite_bases_equal(negative, read, reference)
    } else {
        bases_equal(read, reference)
    }
}

/// htsjdk `SequenceUtil.sumQualitiesOfMismatches` with no reference
/// offset: the qualities (Java bytes) of the aligned bases that differ from
/// the reference.
pub fn sum_qualities_of_mismatches(
    cigar: impl Iterator<Item = (u32, u8)>,
    read: &[u8],
    qualities: &[u8],
    reference: &[u8],
    alignment_start: i32,
    negative: bool,
    bisulfite: bool,
) -> std::result::Result<i64, &'static str> {
    if alignment_start <= 0 {
        return Err("an alignment start of 0");
    }
    let mut sum = 0i64;
    let mut quality_past_the_read = false;
    aligned_pairs(
        cigar,
        read,
        reference,
        alignment_start,
        |at, read_base, ref_base| {
            if !same_base(bisulfite, negative, read_base, ref_base) {
                match base_at(qualities, at) {
                    Some(q) => sum += i64::from(q as i8),
                    None => quality_past_the_read = true,
                }
            }
        },
    )?;
    if quality_past_the_read {
        return Err("a quality past the read");
    }
    Ok(sum)
}

/// htsjdk `SequenceUtil.calculateSamNmTag` of bisulfite data: the
/// mismatches over the aligned blocks, then the length of every insertion
/// and deletion of the CIGAR.
pub fn bisulfite_nm(
    cigar: &[(u32, u8)],
    read: &[u8],
    reference: &[u8],
    alignment_start: i32,
    negative: bool,
) -> std::result::Result<i64, &'static str> {
    let mut nm = 0i64;
    aligned_pairs(
        cigar.iter().copied(),
        read,
        reference,
        alignment_start,
        |_, read_base, ref_base| {
            if !bisulfite_bases_equal(negative, read_base, ref_base) {
                nm += 1;
            }
        },
    )?;
    Ok(nm
        + cigar
            .iter()
            .filter(|(_, op)| matches!(op, 1 | 2))
            .map(|(length, _)| i64::from(*length))
            .sum::<i64>())
}

/// `IS_BISULFITE_SEQUENCE` and `SET_ONLY_UQ`.
#[derive(Clone, Copy, Debug, Default)]
struct Flags {
    bisulfite: bool,
    only_uq: bool,
}

impl Flags {
    fn of(inv: &Invocation) -> Result<Self> {
        Ok(Self {
            bisulfite: inv.flag("IS_BISULFITE_SEQUENCE", false)?,
            only_uq: inv.flag("SET_ONLY_UQ", false)?,
        })
    }
}

/// The edit: the tags of each mapped record from the contig the walker
/// holds.
struct SetTags {
    reference: Reference,
    flags: Flags,
}

/// The header check of `doWork`, before the reference is opened.
fn check_order(header: &sam::Header) -> Result<()> {
    match SortOrder::of(header) {
        SortOrder::Coordinate => Ok(()),
        order => Err(Error::Metric(format!(
            "Input must be coordinate-sorted for this program to run. Found: {}",
            order.name()
        ))),
    }
}

impl SetTags {
    fn new(setup: &Setup, flags: Flags) -> Result<Self> {
        check_order(setup.header)?;
        let fasta = setup.reference.ok_or_else(|| {
            Error::Usage("REFERENCE_SEQUENCE is required: the tags are computed from it".to_owned())
        })?;
        Ok(Self {
            reference: Reference::open(fasta, setup.header)?,
            flags,
        })
    }

    /// The contig of a mapped record, as `ReferenceSequenceFileWalker.get`
    /// gives it: forward only.
    fn contig(&mut self, record: &RawRecord) -> Result<&[u8]> {
        let id = record.reference_id();
        match usize::try_from(id) {
            Ok(id) => self.reference.contig(id),
            Err(_) => Err(fails(
                record,
                "mapped with no reference: Reference sequence (-1) not found",
            )),
        }
    }
}

/// A record's qualities as htsjdk hands them out: none for a record without
/// bases or whose first quality is `0xFF`.
fn qualities(record: &RawRecord) -> Option<&[u8]> {
    let qualities = record.qualities();
    match qualities.first() {
        None | Some(0xFF) => None,
        Some(_) => Some(qualities),
    }
}

/// htsjdk decodes a record whose first quality is `0xFF` as one without
/// qualities and encodes it back with `0xFF` for every base, so a record
/// that gets a tag has all its quality bytes set to `0xFF`. True when a
/// byte changed.
fn clear_missing_qualities(record: &mut RawRecord) -> bool {
    let qualities = record.qualities_mut();
    if qualities.first() != Some(&0xFF) || qualities.iter().all(|&q| q == 0xFF) {
        return false;
    }
    qualities.fill(0xFF);
    true
}

impl RecordEdit for SetTags {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        check_order(input)?;
        Ok(input.clone())
    }

    fn rewrites(&self) -> bool {
        true
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        let flags = self.flags;
        let has_qualities = qualities(record).is_some();
        if record.flags() & 0x4 != 0 || (flags.only_uq && !has_qualities) {
            return Ok(Outcome::Same);
        }
        let reference = self.contig(record)?;
        let (bases, cigar) = (record.bases(), record.cigar().collect::<Vec<_>>());
        let (start, negative) = (record.alignment_start(), record.flags() & 0x10 != 0);
        let mut set = Vec::with_capacity(3);
        if !flags.only_uq {
            let (md, nm) = md_and_nm(cigar.iter().copied(), &bases, reference, start)
                .map_err(|what| fails(record, what))?;
            let nm = if flags.bisulfite {
                bisulfite_nm(&cigar, &bases, reference, start, negative)
                    .map_err(|what| fails(record, what))?
            } else {
                nm
            };
            set.push((MD, TagValue::String(md)));
            set.push((NM, TagValue::Int(nm)));
        }
        if let Some(quals) = qualities(record) {
            let uq = sum_qualities_of_mismatches(
                cigar.iter().copied(),
                &bases,
                quals,
                reference,
                start,
                negative,
                flags.bisulfite,
            )
            .map_err(|what| fails(record, what))?;
            set.push((UQ, TagValue::Int(uq)));
        }
        let before = record.tag_bytes().to_vec();
        let mut tags = record.tags().map_err(Error::Record)?;
        for (name, value) in set {
            tags.set(name, value);
        }
        record.set_tags(&tags).map_err(Error::Record)?;
        let qualities_cleared = clear_missing_qualities(record);
        let same = record.tag_bytes() == before.as_slice() && !qualities_cleared;
        Ok(if same {
            Outcome::Same
        } else {
            Outcome::Changed
        })
    }
}

fn build(flags: Flags) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |setup| Ok(Box::new(SetTags::new(setup, flags)?))
}

/// `REFERENCE_SEQUENCE`, which `tool` cannot run without.
fn check_reference(inv: &Invocation, tool: &str) -> Result<()> {
    match inv.reference() {
        Some(fasta) if fasta != Path::new("") => Ok(()),
        _ => Err(Error::Usage(format!(
            "{tool}: REFERENCE_SEQUENCE is required: the tags are computed from it"
        ))),
    }
}

fn run_tool(inv: &Invocation, tool: &str) -> Result<ToolOutput> {
    let flags = Flags::of(inv)?;
    check_reference(inv, tool)?;
    let opts = record_edit::Options::from_invocation(inv)?;
    record_edit::run(&opts, build(flags))?;
    Ok(ToolOutput::NONE)
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, a flag that is no boolean, no reference
/// (the chain's `-R` is the step's).
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    Flags::of(inv)?;
    check_reference(inv, "the step (give the chain -R FASTA)")
}

/// The tool as a step of a chain: the stream at its position, the chain's
/// `-R` as the reference.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let flags = Flags::of(inv)?;
    record_edit::stage(inv, frame, build(flags))
}

#[cfg(test)]
mod tests;

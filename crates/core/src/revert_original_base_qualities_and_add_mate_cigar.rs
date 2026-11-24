//! `RevertOriginalBaseQualitiesAndAddMateCigar`, mirroring Picard 3.5.0
//! (`picard.sam.RevertOriginalBaseQualitiesAndAddMateCigar`).
//!
//! First the early exit (`canSkipSAMFile`): the records are read in order
//! up to `MAX_RECORDS_TO_EXAMINE`. A record with an `OQ` (with
//! `RESTORE_ORIGINAL_QUALITIES`) decides that the file is reverted; else
//! the first paired record with a mapped mate decides: without `MC` the
//! file is reverted, with `MC` it is skipped. A file that ends before any
//! record decides is skipped; one that holds more records than were
//! examined is reverted. A skipped file writes nothing, not even `OUTPUT`,
//! and exits 0 with Picard's line `Can skip the SAM/BAM/CRAM file` on
//! standard error. `MAX_RECORDS_TO_EXAMINE=0` reverts any file with a
//! record.
//!
//! Then each record, as read: the soft clip of `CleanSam` where the
//! alignment or the mate's `MC` ends past its reference
//! (`AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference`), and
//! with `RESTORE_ORIGINAL_QUALITIES` (default true) `QUAL` replaced by `OQ`
//! (each character less 33; a character outside `!` to `~` or a length
//! other than the read's fails) and `OQ` removed. The records are sorted in
//! queryname order and the records of each read name go through htsjdk's
//! `SamPairUtil.SetMateInfoIterator` with `setMateCigar` true
//! ([`set_mate_info`]): the mate fields, `MQ`, `MC` and `TLEN` of the
//! primary pair, and of each paired supplementary record from its mate. A
//! name with two primary first (or second) records, or a primary first
//! record without a primary second (or the reverse), fails.
//!
//! The output header is the input's with `SO` set to `SORT_ORDER` (default
//! the input's). Picard writes through a writer that sorts: `queryname`,
//! `coordinate` and `duplicate` are sorted again, `unsorted` is written in
//! the queryname order of the mate pass, `unknown` fails on a second
//! record. Both sorts hold `MAX_RECORDS_IN_RAM` records per run, spilled to
//! `TMP_DIR`. `CREATE_INDEX` and `CREATE_MD5_FILE` default to true; for an
//! `OUTPUT` ending in `.cram` the index defaults to false, since no index is
//! written beside a CRAM (Picard writes `<CRAM>.bai` there). Picard
//! reads with `EAGERLY_DECODE`, so every record is encoded again: tags in
//! [`crate::bam_raw::Tags`] order and widths, the `bin` htsjdk computes,
//! qualities whose first byte is `0xFF` written as all `0xFF`. No `@PG`
//! line is added.
//!
//! In a chain the step is a barrier that leaves `reverted`: it holds every
//! record, sorts them by name, sets the mate information and hands them on
//! in `SORT_ORDER` (default the stream's). A chain step never skips its
//! stream, so `MAX_RECORDS_TO_EXAMINE` is a usage error there, as are
//! `SORT_ORDER=unknown` and a stream whose order is unknown.

use noodles::sam;

use crate::add_or_replace_read_groups::htsjdk_qualities;
use crate::bam_raw::{RawRecord, TagValue};
use crate::clean_sam::{Clean, cigar_text};
use crate::error::{Error, Result};
use crate::header_merge::{HeaderTag, SortOrder, set_header_tag};
use crate::metric::ToolOutput;
use crate::record_edit::{self, Gate, Outcome, RecordEdit, Reorder, Setup, WriterOrder};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text::parse_cigar;
use crate::spill::TempDir;
use crate::stage::{Frame, Stage, WriteFormat};

pub const TOOL: Tool = Tool {
    name: "RevertOriginalBaseQualitiesAndAddMateCigar",
    summary: "no rows; OUTPUT gets every record with QUAL restored from OQ, a CIGAR past the \
              end of its reference soft-clipped, and MC, MQ and the mate fields set from each \
              mate, in SORT_ORDER: a BAM, SAM text for a name ending in .sam, or a CRAM (with \
              REFERENCE_SEQUENCE) for one ending in .cram. A file whose first records need \
              nothing is skipped and OUTPUT is not written.",
    options: &[
        Opt::text(
            &["SORT_ORDER", "SO"],
            "The order of OUTPUT: unsorted, queryname, coordinate, duplicate or unknown \
             (default: the input's).",
        ),
        Opt::bool(
            &["RESTORE_ORIGINAL_QUALITIES", "OQ"],
            "Replace QUAL by the OQ tag where a record has one, and remove OQ (default true).",
        ),
        Opt::int(
            &["MAX_RECORDS_TO_EXAMINE"],
            "Records read to decide that the file needs nothing (no OQ to restore, MC on the \
             first paired record with a mapped mate) and is skipped; 0 never skips a file with \
             a record (default 10000).",
        ),
        record_edit::OUTPUT_OPTION,
        record_edit::COMPRESSION_OPTION,
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside a coordinate-sorted BAM OUTPUT (default true).",
        ),
        Opt::bool(&["CREATE_MD5_FILE"], "Write OUTPUT.md5 (default true)."),
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more (default \
             500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs (default: the system temporary directory).",
        )
        .repeated(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

const PAIRED: u16 = 0x1;
const UNMAPPED: u16 = 0x4;
const MATE_UNMAPPED: u16 = 0x8;
const REVERSE: u16 = 0x10;
const MATE_REVERSE: u16 = 0x20;
const FIRST: u16 = 0x40;
const SECOND: u16 = 0x80;
const SECONDARY: u16 = 0x100;
const SUPPLEMENTARY: u16 = 0x800;

fn failure(record: &RawRecord, what: &str) -> Error {
    failure_of(TOOL.name, record, what)
}

/// The error of `tool` on `record`.
fn failure_of(tool: &str, record: &RawRecord, what: &str) -> Error {
    Error::Metric(format!(
        "{tool}: record {}: {what}",
        String::from_utf8_lossy(record.name())
    ))
}

/// `getOriginalBaseQualities`: `OQ` as Phred scores, `None` without it or
/// for an empty one. An `OQ` that is no string, or holds a character
/// outside `!` to `~`, fails as htsjdk fails.
fn original_qualities(record: &RawRecord) -> Result<Option<Vec<u8>>> {
    match record.tag(*b"OQ").map_err(Error::Record)? {
        None => Ok(None),
        Some(TagValue::String(text)) if text.is_empty() => Ok(None),
        Some(TagValue::String(text)) => text
            .iter()
            .map(|&c| match c {
                b'!'..=b'~' => Ok(c - 33),
                _ => Err(failure(
                    record,
                    &format!("Invalid fastq character in OQ: {}", char::from(c)),
                )),
            })
            .collect::<Result<Vec<u8>>>()
            .map(Some),
        Some(_) => Err(failure(record, "OQ is not a String")),
    }
}

/// `SAMUtils.getMateCigar`: true when `MC` holds a CIGAR. An `MC` that is
/// no string or no CIGAR fails.
fn has_mate_cigar(record: &RawRecord) -> Result<bool> {
    match record.tag(*b"MC").map_err(Error::Record)? {
        None => Ok(false),
        Some(TagValue::String(text)) => parse_cigar(&text)
            .map(|_| true)
            .map_err(|e| failure(record, &e)),
        Some(_) => Err(failure(record, "Value for tag MC is not a String")),
    }
}

/// Picard's `CanSkipSamFile`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Decision {
    CanSkip,
    FoundOq,
    FoundNoMc,
    NoEvidence,
}

impl Decision {
    fn message(self, examined: i64) -> String {
        match self {
            Self::CanSkip => "Can skip the SAM/BAM/CRAM file".to_owned(),
            Self::FoundOq => {
                "Cannot skip the SAM/BAM/CRAM as we found a record with an OQ".to_owned()
            }
            Self::FoundNoMc => {
                "Cannot skip the SAM/BAM/CRAM as we found a mate with no mate cigar tag".to_owned()
            }
            Self::NoEvidence => format!(
                "Found no evidence of OQ or mate with no mate cigar in the first {examined} \
                 records.  Will continue..."
            ),
        }
    }
}

/// `canSkipSAMFile` over the records as they are read.
struct Skip {
    examine: i64,
    restore: bool,
    seen: i64,
    decision: Option<Decision>,
}

impl Skip {
    fn decide(&mut self, decision: Decision) -> Option<bool> {
        self.decision = Some(decision);
        Some(decision == Decision::CanSkip)
    }
}

impl Gate for Skip {
    fn look(&mut self, record: &RawRecord) -> Result<Option<bool>> {
        if self.seen >= self.examine {
            return Ok(self.decide(Decision::NoEvidence));
        }
        if self.restore && original_qualities(record)?.is_some() {
            return Ok(self.decide(Decision::FoundOq));
        }
        let flags = record.flags();
        if flags & PAIRED != 0 && flags & MATE_UNMAPPED == 0 {
            let decision = if has_mate_cigar(record)? {
                Decision::CanSkip
            } else {
                Decision::FoundNoMc
            };
            return Ok(self.decide(decision));
        }
        self.seen += 1;
        Ok(None)
    }

    fn at_end(&mut self) -> bool {
        self.decision = Some(Decision::CanSkip);
        true
    }
}

/// The edit of each record as read: the clip past a reference end, then
/// `OQ` into `QUAL`, then the encoding of a record htsjdk decoded.
struct Revert {
    clean: Clean,
    restore: bool,
}

impl RecordEdit for Revert {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(input.clone())
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        let mut outcome = Outcome::Same;
        if self.clean.clip(record)? {
            outcome = Outcome::Changed;
        }
        htsjdk_qualities(record);
        let mut tags = record.tags().map_err(Error::Record)?;
        if self.restore
            && let Some(qualities) = original_qualities(record)?
        {
            if qualities.len() != record.sequence_length() {
                return Err(failure(
                    record,
                    &format!(
                        "Mismatch between read length and quals length writing read; read \
                         length: {}; quals length: {}",
                        record.sequence_length(),
                        qualities.len()
                    ),
                ));
            }
            record.qualities_mut().copy_from_slice(&qualities);
            tags.remove(*b"OQ");
            outcome = Outcome::Changed;
        }
        record.set_tags(&tags).map_err(Error::Record)?;
        Ok(outcome)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

/// The output header: `SO` set to `SORT_ORDER`, or the input's order.
struct Ordered {
    edit: Revert,
    order: Option<SortOrder>,
}

impl RecordEdit for Ordered {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        let order = self.order.unwrap_or_else(|| SortOrder::of(input));
        set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
        Ok(header)
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        self.edit.record(record)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

fn set_tag_int(record: &mut RawRecord, name: [u8; 2], value: Option<i64>) -> Result<()> {
    let mut tags = record.tags().map_err(Error::Record)?;
    match value {
        Some(value) => tags.set(name, TagValue::Int(value)),
        None => tags.remove(name),
    }
    record.set_tags(&tags).map_err(Error::Record)
}

fn set_mate_cigar(record: &mut RawRecord, cigar: Option<Vec<u8>>) -> Result<()> {
    let mut tags = record.tags().map_err(Error::Record)?;
    match cigar {
        Some(text) => tags.set(*b"MC", TagValue::String(text)),
        None => tags.remove(*b"MC"),
    }
    record.set_tags(&tags).map_err(Error::Record)
}

fn set_flag(record: &mut RawRecord, flag: u16, on: bool) {
    let flags = record.flags();
    record.set_flags(if on { flags | flag } else { flags & !flag });
}

/// The mate fields of `record` from `mate`: reference, start, strand, the
/// mate-unmapped flag as given.
fn mate_fields(record: &mut RawRecord, mate: &RawRecord, mate_unmapped: bool) {
    record.set_reference_ids(record.reference_id(), mate.reference_id());
    record.set_mate_alignment_start(mate.alignment_start());
    set_flag(record, MATE_REVERSE, mate.flags() & REVERSE != 0);
    set_flag(record, MATE_UNMAPPED, mate_unmapped);
}

/// `getCigarString`: `*` for no operator.
fn cigar_string(record: &RawRecord) -> Vec<u8> {
    cigar_text(&record.cigar().collect::<Vec<_>>())
}

/// htsjdk `getAlignmentEnd`: the last reference base, 0 when unmapped.
fn alignment_end(record: &RawRecord) -> i32 {
    if record.flags() & UNMAPPED != 0 {
        return 0;
    }
    record
        .alignment_start()
        .wrapping_add(record.reference_length())
        .wrapping_sub(1)
}

/// `SamPairUtil.computeInsertSize`, in Java `int` arithmetic.
fn compute_insert_size(a: &RawRecord, b: &RawRecord) -> i32 {
    if a.flags() & UNMAPPED != 0 || b.flags() & UNMAPPED != 0 {
        return 0;
    }
    if a.reference_id() != b.reference_id() {
        return 0;
    }
    let five_prime = |r: &RawRecord| {
        if r.flags() & REVERSE != 0 {
            alignment_end(r)
        } else {
            r.alignment_start()
        }
    };
    let (first, second) = (five_prime(a), five_prime(b));
    let adjustment = if second >= first { 1 } else { -1 };
    second.wrapping_sub(first).wrapping_add(adjustment)
}

/// `SamPairUtil.setMateInfo(a, b, mate_cigar)`: without `mate_cigar` no
/// record keeps an `MC`.
fn set_mate_info(a: &mut RawRecord, b: &mut RawRecord, mate_cigar: bool) -> Result<()> {
    let a_unmapped = a.flags() & UNMAPPED != 0;
    let b_unmapped = b.flags() & UNMAPPED != 0;
    if !a_unmapped && !b_unmapped {
        mate_fields(a, b, false);
        set_tag_int(a, *b"MQ", Some(i64::from(b.mapping_quality())))?;
        mate_fields(b, a, false);
        set_tag_int(b, *b"MQ", Some(i64::from(a.mapping_quality())))?;
        let (a_cigar, b_cigar) = (cigar_string(a), cigar_string(b));
        set_mate_cigar(a, mate_cigar.then_some(b_cigar))?;
        set_mate_cigar(b, mate_cigar.then_some(a_cigar))?;
    } else if a_unmapped && b_unmapped {
        let (a_reverse, b_reverse) = (a.flags() & REVERSE != 0, b.flags() & REVERSE != 0);
        for (record, mate_reverse) in [(&mut *a, b_reverse), (&mut *b, a_reverse)] {
            record.set_reference_ids(-1, -1);
            record.set_alignment_start(0);
            record.set_mate_alignment_start(0);
            set_flag(record, MATE_REVERSE, mate_reverse);
            set_flag(record, MATE_UNMAPPED, true);
            set_tag_int(record, *b"MQ", None)?;
            set_mate_cigar(record, None)?;
            record.set_template_length(0);
        }
    } else {
        let (mapped, unmapped) = if a_unmapped {
            (&mut *b, &mut *a)
        } else {
            (&mut *a, &mut *b)
        };
        unmapped.set_reference_ids(mapped.reference_id(), unmapped.mate_reference_id());
        unmapped.set_alignment_start(mapped.alignment_start());
        mate_fields(mapped, unmapped, true);
        set_tag_int(mapped, *b"MQ", None)?;
        set_mate_cigar(mapped, None)?;
        mapped.set_template_length(0);
        mate_fields(unmapped, mapped, false);
        set_tag_int(unmapped, *b"MQ", Some(i64::from(mapped.mapping_quality())))?;
        set_mate_cigar(unmapped, mate_cigar.then(|| cigar_string(mapped)))?;
        unmapped.set_template_length(0);
    }
    let insert = compute_insert_size(a, b);
    a.set_template_length(insert);
    b.set_template_length(insert.wrapping_neg());
    Ok(())
}

/// `setMateInformationOnSupplementalAlignment(record, mate, mate_cigar)`.
fn set_supplementary_mate(
    record: &mut RawRecord,
    mate: &RawRecord,
    mate_cigar: bool,
) -> Result<()> {
    let mate_unmapped = mate.flags() & UNMAPPED != 0;
    mate_fields(record, mate, mate_unmapped);
    record.set_template_length(mate.template_length().wrapping_neg());
    set_mate_cigar(
        record,
        (mate_cigar && !mate_unmapped).then(|| cigar_string(mate)),
    )?;
    set_tag_int(record, *b"MQ", Some(i64::from(mate.mapping_quality())))
}

/// Two distinct records of `group` at once.
fn two(group: &mut [RawRecord], i: usize, j: usize) -> (&mut RawRecord, &mut RawRecord) {
    if i < j {
        let (left, right) = group.split_at_mut(j);
        (&mut left[i], &mut right[0])
    } else {
        let (left, right) = group.split_at_mut(i);
        (&mut right[0], &mut left[j])
    }
}

/// `SamPairUtil.SetMateInfoIterator.advance` on the records of one read
/// name, in queryname order, with `setMateCigar` true and missing mates an
/// error. Each record then gets the `bin` htsjdk's writer computes: an
/// unmapped read placed beside its mate has a new start.
pub fn mate_group(group: &mut [RawRecord]) -> Result<()> {
    mate_group_with(group, true, false, TOOL.name)
}

/// [`mate_group`] as `SetMateInfoIterator(iterator, mate_cigar,
/// ignore_missing)` makes it: without `mate_cigar` no record of a pair
/// keeps an `MC`; with `ignore_missing` a primary record whose mate is not
/// in the group is handed on as it came (`FixMateInformation`). An error
/// names `tool`.
pub fn mate_group_with(
    group: &mut [RawRecord],
    mate_cigar: bool,
    ignore_missing: bool,
    tool: &str,
) -> Result<()> {
    set_mates(group, mate_cigar, ignore_missing, tool)?;
    group.iter_mut().for_each(RawRecord::set_codec_bin);
    Ok(())
}

fn set_mates(
    group: &mut [RawRecord],
    mate_cigar: bool,
    ignore_missing: bool,
    tool: &str,
) -> Result<()> {
    let mut first = None;
    let mut second = None;
    let mut supplementary = false;
    for (at, record) in group.iter().enumerate() {
        let flags = record.flags();
        if flags & PAIRED == 0 {
            continue;
        }
        if flags & (SECONDARY | SUPPLEMENTARY) == 0 {
            if flags & FIRST != 0 {
                if first.is_some() {
                    return Err(failure_of(
                        tool,
                        record,
                        "Found two records that are paired, not supplementary, and first of \
                         the pair",
                    ));
                }
                first = Some(at);
            } else if flags & SECOND != 0 {
                if second.is_some() {
                    return Err(failure_of(
                        tool,
                        record,
                        "Found two records that are paired, not supplementary, and second of \
                         the pair",
                    ));
                }
                second = Some(at);
            }
        }
        if flags & SUPPLEMENTARY != 0 {
            supplementary = true;
        }
    }
    let (Some(first), Some(second)) = (first, second) else {
        if ignore_missing {
            return Ok(());
        }
        if let Some(at) = first {
            return Err(failure_of(tool, &group[at], "Missing second read of pair"));
        }
        if let Some(at) = second {
            return Err(failure_of(tool, &group[at], "Missing first read of pair"));
        }
        return Ok(());
    };
    let (a, b) = two(group, first, second);
    set_mate_info(a, b, mate_cigar)?;
    if supplementary {
        let (first_record, second_record) = (group[first].clone(), group[second].clone());
        for record in group.iter_mut() {
            let flags = record.flags();
            if flags & PAIRED == 0 || flags & SUPPLEMENTARY == 0 {
                continue;
            }
            let mate = if flags & FIRST != 0 {
                &second_record
            } else {
                &first_record
            };
            set_supplementary_mate(record, mate, mate_cigar)?;
        }
    }
    Ok(())
}

/// The options of both faces.
struct RevertOptions {
    order: Option<SortOrder>,
    restore: bool,
}

impl RevertOptions {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let order = inv
            .get("SORT_ORDER")
            .filter(|value| *value != "null")
            .map(|value| {
                SortOrder::parse(value).ok_or_else(|| {
                    Error::Usage(format!("'{value}' is not a valid value for SortOrder."))
                })
            })
            .transpose()?;
        Ok(Self {
            order,
            restore: inv.flag("RESTORE_ORIGINAL_QUALITIES", true)?,
        })
    }

    fn build(self) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
        move |setup| {
            Ok(Box::new(Ordered {
                edit: Revert {
                    clean: Clean::new(setup.header, TOOL.name),
                    restore: self.restore,
                },
                order: self.order,
            }))
        }
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let options = RevertOptions::from_invocation(inv)?;
    let examine = inv.parse::<i64>("MAX_RECORDS_TO_EXAMINE", 10_000)?;
    let mut opts = record_edit::Options::from_invocation(inv)?;
    // `CREATE_INDEX` defaults to true, as in Picard, which then writes
    // `<CRAM>.bai` beside a CRAM. bammetrics writes no index of a CRAM, so
    // the default asks for none there; `CREATE_INDEX=true` with a CRAM is
    // refused as for every tool (`cram_writer::index_refusal`).
    let cram = WriteFormat::of(&opts.output) == WriteFormat::Cram;
    opts.create_index = inv.flag("CREATE_INDEX", !cram)?;
    opts.create_md5 = inv.flag("CREATE_MD5_FILE", true)?;
    let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
    opts.reorder = Reorder {
        regroup: Some(mate_group),
        writer: WriterOrder::Header {
            order: options.order,
            presorted_if_same: false,
        },
        limits,
        tmp_dir,
    };
    let mut skip = Skip {
        examine,
        restore: options.restore,
        seen: 0,
        decision: None,
    };
    let ran = record_edit::run_gated(&opts, options.build(), Some(&mut skip))?;
    if ran.is_none() {
        let decision = skip.decision.unwrap_or(Decision::CanSkip);
        if !inv.quiet {
            eprintln!(
                "{}: {}; OUTPUT is not written",
                TOOL.name,
                decision.message(examine)
            );
        }
    }
    Ok(ToolOutput::NONE)
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, `MAX_RECORDS_TO_EXAMINE` (a step never
/// skips its stream), `SORT_ORDER=unknown`.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    if inv.get("MAX_RECORDS_TO_EXAMINE").is_some() {
        return Err(Error::Usage(
            "MAX_RECORDS_TO_EXAMINE: a chain step reverts every record of its stream; the \
             early exit of a file in shape belongs to the tool on its own"
                .to_owned(),
        ));
    }
    let options = RevertOptions::from_invocation(inv)?;
    if options.order == Some(SortOrder::Unknown) {
        return Err(Error::Usage(
            "SORT_ORDER=unknown in a chain: Picard sorts unknown with no comparator and fails"
                .to_owned(),
        ));
    }
    record_edit::stage_limits(inv, 2).map(|_| ())
}

/// `RevertOriginalBaseQualitiesAndAddMateCigar` as a step of a chain: a
/// barrier that sorts twice, by name for the mates, then in `SORT_ORDER`.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let options = RevertOptions::from_invocation(inv)?;
    let reorder = Reorder {
        regroup: Some(mate_group),
        writer: WriterOrder::Header {
            order: options.order,
            presorted_if_same: false,
        },
        limits: record_edit::stage_limits(inv, 2)?,
        tmp_dir: None,
    };
    let Stage::Transform(transform) = record_edit::stage(inv, frame, options.build())? else {
        return Err(Error::Metric(format!(
            "{}: the edit of a step is a transform",
            TOOL.name
        )));
    };
    record_edit::sorting_stage(inv, frame, temp, &reorder, transform)
}

#[cfg(test)]
mod tests;

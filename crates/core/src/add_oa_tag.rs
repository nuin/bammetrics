//! `AddOATag`, mirroring Picard 3.5.0 (`picard.sam.AddOATag`): each record
//! gets its alignment appended to its `OA` tag, then every record is
//! written through a presorted writer ([`Mode::Rewrite`]).
//!
//! The entry is `RNAME,POS,STRAND,CIGAR,MAPQ,NM;` for a mapped record and
//! `*,0,STRAND,*,255,;` for one with the unmapped flag (placed or not).
//! `STRAND` is `-` with the reverse flag, else `+`. `NM` is the record's
//! `NM` value as Java's `String.valueOf` prints it (an integer in decimal,
//! a float through `Float.toString`, a `Z` or `A` value as its text), empty
//! without one. The entry goes after the record's own `OA` value, printed
//! the same way, and `OA` becomes a `Z` tag at its place in htsjdk's tag
//! order ([`crate::bam_raw::Tags`]). Secondary and supplementary records
//! are tagged as any other.
//!
//! `INTERVAL_LIST` (htsjdk `IntervalList.fromFile`: an interval list with
//! its own `@SQ` lines, a BED refused with `Interval list file must contain
//! header.`) limits the tag to the records `OverlapDetector.overlapsAny`
//! finds: a record without the unmapped flag whose reference name is a
//! contig of the list and whose `POS ..= POS + reference length - 1` shares
//! a base with an interval. An unmapped record and a record whose CIGAR
//! covers no reference base never overlap. A record left untagged keeps
//! its tag bytes.
//!
//! Where Picard prints a Java object reference (an `NM` or `OA` value that
//! is an array or `H`, `[B@1a2b3c`), bammetrics fails with the record's
//! name. A reference name with a comma fails as in Picard, for the records
//! the tool tags. No `@PG` line is added.
//!
//! In a chain the step tags the records of the stream at its position.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::{RawRecord, TagValue};
use crate::error::{Error, Result};
use crate::intervals::{self, Interval, OverlapIndex};
use crate::metric::ToolOutput;
use crate::record_edit::{self, Outcome, RecordEdit, Setup};
use crate::registry::{Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub const TOOL: Tool = Tool {
    name: "AddOATag",
    summary: "no rows; OUTPUT gets every record of the input with its alignment appended to \
              its OA tag (RNAME,POS,STRAND,CIGAR,MAPQ,NM;), or only the records over \
              INTERVAL_LIST: a BAM, SAM text for a name ending in .sam, or a CRAM (with \
              REFERENCE_SEQUENCE) for one ending in .cram.",
    options: &[
        Opt::path(
            &["INTERVAL_LIST", "L"],
            "An interval list (with @SQ lines; a BED is refused, as in Picard): only records \
             that overlap one of its intervals get the OA tag.",
        ),
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

const OA: [u8; 2] = *b"OA";
const NM: [u8; 2] = *b"NM";
const REVERSE: u16 = 0x10;
const UNMAPPED: u16 = 0x4;

/// `INTERVAL_LIST`: its contig names and its intervals, uniqued, without
/// the empty ones (`OverlapDetector.addLhs` adds no interval without a base).
struct Intervals {
    names: Vec<String>,
    index: OverlapIndex,
}

impl Intervals {
    fn read(path: &Path) -> Result<Self> {
        let (names, named) = intervals::read_interval_list(path)?;
        let list: Vec<Interval> = named.iter().map(|named| named.interval).collect();
        let unique: Vec<Interval> = intervals::uniqued(&list)
            .into_iter()
            .filter(|interval| interval.start <= interval.end)
            .collect();
        let index = OverlapIndex::new(unique, names.len(), 0);
        Ok(Self { names, index })
    }
}

/// The `INTERVAL_LIST` option, read.
fn interval_list(inv: &Invocation) -> Result<Option<(PathBuf, Intervals)>> {
    inv.get("INTERVAL_LIST")
        .map(|path| {
            let path = PathBuf::from(path);
            Intervals::read(&path).map(|list| (path, list))
        })
        .transpose()
}

/// The edit: `OA` of each record the intervals select.
struct AddOa {
    /// The header's reference names, as `getReferenceName` gives them.
    names: Vec<Vec<u8>>,
    /// For each header reference, its contig in `INTERVAL_LIST`; `None`
    /// without the option.
    contigs: Option<Vec<Option<usize>>>,
    index: Option<OverlapIndex>,
}

impl AddOa {
    fn new(setup: &Setup, list: Option<Intervals>) -> Self {
        let names: Vec<Vec<u8>> = setup
            .header
            .reference_sequences()
            .keys()
            .map(|name| name.to_vec())
            .collect();
        let (contigs, index) = match list {
            Some(list) => {
                let by_name: HashMap<&[u8], usize> = list
                    .names
                    .iter()
                    .enumerate()
                    .map(|(at, name)| (name.as_bytes(), at))
                    .collect();
                let contigs = names
                    .iter()
                    .map(|name| by_name.get(name.as_slice()).copied())
                    .collect();
                (Some(contigs), Some(list.index))
            }
            None => (None, None),
        };
        Self {
            names,
            contigs,
            index,
        }
    }

    /// `overlapsAny(record)`: true without `INTERVAL_LIST`.
    fn selects(&self, record: &RawRecord) -> bool {
        let (Some(contigs), Some(index)) = (&self.contigs, &self.index) else {
            return true;
        };
        if record.flags() & UNMAPPED != 0 {
            return false;
        }
        let contig = usize::try_from(record.reference_id())
            .ok()
            .and_then(|id| contigs.get(id).copied().flatten());
        let Some(contig) = contig else {
            return false;
        };
        let start = i64::from(record.alignment_start());
        let end = start + i64::from(record.reference_length()) - 1;
        !index.overlaps(contig, start, end).is_empty()
    }

    /// `getReferenceName`: `*` for no reference.
    fn reference_name(&self, record: &RawRecord) -> Result<&[u8]> {
        let id = record.reference_id();
        if id == -1 {
            return Ok(b"*");
        }
        usize::try_from(id)
            .ok()
            .and_then(|id| self.names.get(id))
            .map(Vec::as_slice)
            .ok_or_else(|| {
                Error::Metric(format!(
                    "record {}: reference index {id} is not in the header, which has {} \
                     sequences",
                    String::from_utf8_lossy(record.name()),
                    self.names.len()
                ))
            })
    }
}

/// `String.valueOf` of a tag value as htsjdk decodes it; `None` for an
/// array or `H` value, which Java prints as an object reference.
fn java_value(value: &TagValue) -> Option<Vec<u8>> {
    match value {
        TagValue::Int(v) => Some(v.to_string().into_bytes()),
        TagValue::Float(f) => Some(crate::java::float_to_string(*f).into_bytes()),
        TagValue::String(text) => Some(text.clone()),
        TagValue::Char(c) => Some(vec![*c]),
        TagValue::Hex(_) | TagValue::Array(..) => None,
    }
}

/// The text of tag `name` for the `OA` entry: empty without the tag.
fn tag_text(record: &RawRecord, name: [u8; 2], value: Option<&TagValue>) -> Result<Vec<u8>> {
    let Some(value) = value else {
        return Ok(Vec::new());
    };
    java_value(value).ok_or_else(|| {
        Error::Metric(format!(
            "record {}: its {} tag is an array, which Picard prints as a Java object \
             reference in the OA tag",
            String::from_utf8_lossy(record.name()),
            String::from_utf8_lossy(&name)
        ))
    })
}

/// `getCigarString`: `*` for no operator.
fn cigar_text(record: &RawRecord, out: &mut Vec<u8>) -> Result<()> {
    const OPS: &[u8; 9] = b"MIDNSHP=X";
    let mut any = false;
    for (length, op) in record.cigar() {
        let op = OPS.get(usize::from(op)).ok_or_else(|| {
            Error::Metric(format!(
                "record {}: CIGAR operator code {op} is not one of MIDNSHP=X",
                String::from_utf8_lossy(record.name())
            ))
        })?;
        out.extend_from_slice(length.to_string().as_bytes());
        out.push(*op);
        any = true;
    }
    if !any {
        out.push(b'*');
    }
    Ok(())
}

impl RecordEdit for AddOa {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        Ok(input.clone())
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        if !self.selects(record) {
            return Ok(Outcome::Same);
        }
        let name = self.reference_name(record)?;
        if name.contains(&b',') {
            return Err(Error::Metric(format!(
                "Reference name for record {} contains a comma character.",
                String::from_utf8_lossy(record.name())
            )));
        }
        let flags = record.flags();
        let strand = if flags & REVERSE != 0 { b'-' } else { b'+' };
        let mut tags = record.tags().map_err(Error::Record)?;
        let mut entry = tag_text(record, OA, tags.get(OA))?;
        if flags & UNMAPPED != 0 {
            entry.extend_from_slice(b"*,0,");
            entry.push(strand);
            entry.extend_from_slice(b",*,255,;");
        } else {
            entry.extend_from_slice(name);
            entry.push(b',');
            entry.extend_from_slice(record.alignment_start().to_string().as_bytes());
            entry.push(b',');
            entry.push(strand);
            entry.push(b',');
            cigar_text(record, &mut entry)?;
            entry.push(b',');
            entry.extend_from_slice(record.mapping_quality().to_string().as_bytes());
            entry.push(b',');
            entry.extend(tag_text(record, NM, tags.get(NM))?);
            entry.push(b';');
        }
        tags.set(OA, TagValue::String(entry));
        record.set_tags(&tags).map_err(Error::Record)?;
        Ok(Outcome::Changed)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

fn build(list: Option<Intervals>) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |setup| Ok(Box::new(AddOa::new(setup, list)))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut opts = record_edit::Options::from_invocation(inv)?;
    let list = interval_list(inv)?.map(|(path, list)| {
        opts.reads.push(path);
        list
    });
    record_edit::run(&opts, build(list))?;
    Ok(ToolOutput::NONE)
}

/// The errors of the step, asked before the first read of a chain: an
/// option of the written file, an `INTERVAL_LIST` that cannot be read. The
/// list is read here and again when the step is built.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    interval_list(inv).map(|_| ())
}

/// `AddOATag` as a step of a chain.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    record_edit::check_stage(inv)?;
    let list = interval_list(inv)?.map(|(_, list)| list);
    record_edit::stage(inv, frame, build(list))
}

#[cfg(test)]
mod tests;

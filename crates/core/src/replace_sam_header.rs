//! `ReplaceSamHeader`, mirroring Picard 3.5.0 (`picard.sam.ReplaceSamHeader`):
//! the header of `HEADER` (SAM text, a BAM or a CRAM) in place of the
//! input's, then the input's records.
//!
//! Picard takes one of two htsjdk paths, by the `INPUT` name:
//! * a name ending in `.bam` is block copied ([`Mode::BlockCopy`],
//!   `BamFileIoUtils.reheaderBamFile`). The new header's sort order must be
//!   `unsorted` or the input's, else `Sort order of new header does not
//!   match the original file, needs to be <order>`. The record bytes go out
//!   as they are, so each record keeps its reference index and takes the
//!   name the new dictionary gives that index. `CREATE_INDEX` writes
//!   `<OUTPUT>.bai`.
//! * any other name is read record by record at `SILENT`
//!   ([`Mode::Rewrite`]): the sort orders must be equal (`Sort orders of
//!   INPUT (<a>) and HEADER (<b>) do not agree.`), each record keeps its
//!   reference name and takes that name's index in the new dictionary
//!   (`SAMRecord.setHeader`), and a name the dictionary lacks fails with
//!   `Reference index for '<name>' not found in sequence dictionary.` The
//!   presorted writer checks the order; `OUTPUT` may be SAM text.
//!
//! `HEADER` is read at `VALIDATION_STRINGENCY`; its records are not read.
//! No `@PG` line is added, and `SO` is the one `HEADER` says.
//!
//! Where Picard writes a BAM that no reader opens, bammetrics refuses
//! (docs/07 decision rows): a block copy whose record names a reference
//! index the new dictionary does not have, and a block copy into a `.sam`
//! name, which Picard fills with BAM bytes.
//!
//! In a chain the step replaces the header of the stream at its position
//! as a block copy does: the records are BAM records, and each keeps its
//! reference index.

use std::path::PathBuf;

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::mark_duplicates::SortOrder;
use crate::metric::ToolOutput;
use crate::pass;
use crate::record_edit::{self, Mode, Outcome, RecordEdit, Setup};
use crate::registry::{Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub const TOOL: Tool = Tool {
    name: "ReplaceSamHeader",
    summary: "no rows; OUTPUT gets the header of HEADER, then every record of the input: a \
              BAM; when the input name does not end in .bam, SAM text for a name ending in \
              .sam or a CRAM (with REFERENCE_SEQUENCE) for one ending in .cram.",
    options: &[
        Opt::path(
            &["HEADER"],
            "The SAM, BAM or CRAM whose header replaces the input's.",
        )
        .required(),
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

/// `HEADER`, read at `VALIDATION_STRINGENCY`.
fn new_header(inv: &Invocation) -> Result<(PathBuf, sam::Header)> {
    let path = inv
        .get("HEADER")
        .map(PathBuf::from)
        .ok_or_else(|| Error::Usage("ReplaceSamHeader: HEADER is required".to_owned()))?;
    let header = pass::read_header_at(&path, inv.stringency()?)?;
    Ok((path, header))
}

/// The edit: the new header, and each record's reference indexes kept
/// ([`Mode::BlockCopy`], a chain) or found again by name
/// ([`Mode::Rewrite`]).
struct Replace {
    header: sam::Header,
    mode: Mode,
    /// For a rewrite: the new index of each input reference, `None` for a
    /// name the new dictionary lacks.
    indexes: Vec<Option<i32>>,
    /// The input's reference names, for messages.
    names: Vec<String>,
}

impl Replace {
    fn new(header: sam::Header, setup: &Setup) -> Self {
        let names: Vec<String> = setup
            .header
            .reference_sequences()
            .keys()
            .map(ToString::to_string)
            .collect();
        let indexes = match setup.mode {
            Mode::Rewrite => names
                .iter()
                .map(|name| {
                    header
                        .reference_sequences()
                        .get_index_of(name.as_bytes())
                        .and_then(|i| i32::try_from(i).ok())
                })
                .collect(),
            Mode::BlockCopy | Mode::Chain => Vec::new(),
        };
        Self {
            header,
            mode: setup.mode,
            indexes,
            names,
        }
    }

    /// The new index of reference `id` of the input.
    fn rename(&self, id: i32) -> Result<i32> {
        let Ok(at) = usize::try_from(id) else {
            return Ok(-1);
        };
        if let Some(Some(index)) = self.indexes.get(at) {
            return Ok(*index);
        }
        let name = self.names.get(at).map_or("?", String::as_str);
        Err(Error::Metric(format!(
            "Reference index for '{name}' not found in sequence dictionary."
        )))
    }

    /// A reference index the new dictionary does not have.
    fn check_kept(&self, record: &RawRecord) -> Result<()> {
        let count = self.header.reference_sequences().len();
        let fits = |id: i32| usize::try_from(id).map_or(true, |id| id < count);
        for id in [record.reference_id(), record.mate_reference_id()] {
            if !fits(id) {
                return Err(Error::Metric(format!(
                    "record {}: reference index {id} is not in the HEADER sequence dictionary, \
                     which has {count} sequences",
                    String::from_utf8_lossy(record.name())
                )));
            }
        }
        Ok(())
    }
}

impl RecordEdit for Replace {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let (old, new) = (SortOrder::of(input), SortOrder::of(&self.header));
        match self.mode {
            Mode::Rewrite if old != new => Err(Error::Metric(format!(
                "Sort orders of INPUT ({}) and HEADER ({}) do not agree.",
                old.name(),
                new.name()
            ))),
            Mode::BlockCopy | Mode::Chain if new != SortOrder::Unsorted && new != old => {
                Err(Error::Metric(format!(
                    "Sort order of new header does not match the original file, needs to be {}",
                    old.name()
                )))
            }
            _ => Ok(self.header.clone()),
        }
    }

    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome> {
        if self.mode != Mode::Rewrite {
            self.check_kept(record)?;
            return Ok(Outcome::Same);
        }
        let (reference, mate) = (record.reference_id(), record.mate_reference_id());
        let renamed = (self.rename(reference)?, self.rename(mate)?);
        if renamed == (reference, mate) {
            return Ok(Outcome::Same);
        }
        record.set_reference_ids(renamed.0, renamed.1);
        Ok(Outcome::Changed)
    }

    /// A block copy copies the input's blocks as bytes when no record can
    /// name an index past the new dictionary: it is at least as long as the
    /// input's. A shorter one needs each record checked.
    fn copies_bytes(&self, input: &sam::Header) -> bool {
        self.mode == Mode::BlockCopy
            && self.header.reference_sequences().len() >= input.reference_sequences().len()
    }
}

fn build(header: sam::Header) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |setup| Ok(Box::new(Replace::new(header, setup)))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut opts = record_edit::Options::from_invocation(inv)?;
    let (path, header) = new_header(inv)?;
    opts.reads.push(path);
    if opts.input.to_string_lossy().ends_with(".bam") {
        opts.mode = Mode::BlockCopy;
    } else {
        opts.stringency = crate::sam_text::Stringency::Silent;
    }
    record_edit::run(&opts, build(header))?;
    Ok(ToolOutput::NONE)
}

/// `ReplaceSamHeader` as a step of a chain.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    record_edit::check_stage(inv)?;
    let (_, header) = new_header(inv)?;
    record_edit::stage(inv, frame, build(header))
}

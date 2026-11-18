//! `AddCommentsToBam`, mirroring Picard 3.5.0 (`picard.sam.AddCommentsToBam`):
//! the input BAM's header with one `@CO` line per `COMMENT` after the
//! header's own, then every record copied unchanged.
//!
//! Picard writes through htsjdk's `BamFileIoUtils.reheaderBamFile`, a block
//! copy ([`Mode::BlockCopy`]): the input must be a BAM, a name ending in
//! `.sam` is refused with `SAM files are not supported`, and a CRAM or SAM
//! text in another name fails because it is no BAM. The record bytes go out
//! as they were read; `CREATE_INDEX` writes `<OUTPUT>.bai`. A comment
//! starting with `@CO` and a tab keeps it once (`SAMFileHeader.addComment`);
//! a newline is refused, as Picard's command line refuses it. No `@PG` line
//! is added.
//!
//! In a chain the step adds the comments to the header of the stream at
//! its position, whatever the source was read from.

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::record_edit::{self, Mode, Outcome, RecordEdit, Setup};
use crate::registry::{Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::stage::{Frame, Stage};

pub const TOOL: Tool = Tool {
    name: "AddCommentsToBam",
    summary: "no rows; OUTPUT gets the input BAM with one @CO line per COMMENT after the \
              header's own, every record copied unchanged.",
    options: &[
        Opt::text(
            &["COMMENT", "C"],
            "A comment for the header, written as a @CO line; repeat for several. Without \
             one the BAM is copied as it is, as Picard copies it.",
        )
        .repeated(),
        Opt::path(
            &["OUTPUT", "O"],
            "The BAM to write: the input with the comments added.",
        )
        .written()
        .required(),
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

/// `COMMENT`, each as the text after `@CO` and a tab.
fn comments(inv: &Invocation) -> Result<Vec<Vec<u8>>> {
    inv.all("COMMENT")
        .iter()
        .map(|comment| {
            if comment.contains('\n') {
                return Err(Error::Usage(format!(
                    "COMMENT={comment:?}: Comments can not contain a new line"
                )));
            }
            let text = record_edit::java_text(comment);
            Ok(text
                .strip_prefix(b"@CO\t")
                .map_or_else(|| text.clone(), <[u8]>::to_vec))
        })
        .collect()
}

/// The edit: the comments after the header's own, records as they are.
struct Comments(Vec<Vec<u8>>);

impl RecordEdit for Comments {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        for comment in &self.0 {
            header.add_comment(comment.as_slice());
        }
        Ok(header)
    }

    fn record(&mut self, _record: &mut RawRecord) -> Result<Outcome> {
        Ok(Outcome::Same)
    }
}

fn build(comments: Vec<Vec<u8>>) -> impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>> {
    move |_| Ok(Box::new(Comments(comments)))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let comments = comments(inv)?;
    let mut opts = record_edit::Options::from_invocation(inv)?;
    opts.mode = Mode::BlockCopy;
    if opts.input.to_string_lossy().ends_with(".sam") {
        return Err(Error::Input {
            path: opts.input,
            source: std::io::Error::other("SAM files are not supported"),
        });
    }
    record_edit::run(&opts, build(comments))?;
    Ok(ToolOutput::NONE)
}

/// `AddCommentsToBam` as a step of a chain.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    let comments = comments(inv)?;
    record_edit::stage(inv, frame, build(comments))
}

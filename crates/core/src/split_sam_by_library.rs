//! `SplitSamByLibrary`, mirroring Picard 3.5.0
//! (`picard.sam.SplitSamByLibrary`): one file per library of the header's
//! `@RG` lines, in the `OUTPUT` directory (the current directory without
//! one).
//!
//! A library's file is `<LB>.<extension>`: the `LB` value made safe as
//! htsjdk's `IOUtil.makeFileNameSafe` makes it ([`file_name_safe`]), and
//! the extension of the input's kind (`bam`, `sam`). Its header is the
//! input's with only that library's read groups. Every library gets its
//! file before the first record, so a library without records gets a
//! header-only file. A record whose `RG` tag is missing, names no header
//! read group, or names one without `LB` goes to `unknown.<extension>`,
//! made with the first such record, under the input header with only the
//! read groups that have no `LB`. A header with no `LB` on any read group
//! is Picard's exit 2, `No individual libraries are specified in the
//! header of <INPUT>`, with nothing written.
//!
//! Picard writes CRAM files for a CRAM input; bammetrics refuses one
//! before the first read (CRAM output is docs/07 wave W15). Two libraries
//! whose safe names are one file, or a library named `unknown` beside
//! records without a library, make Picard write one file over the other;
//! bammetrics fails and writes none.
//!
//! In a chain the step writes a BAM per library of the stream at its
//! position ([`crate::split_sam::chain_step`]).

use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::registry::{Invocation, Opt, Tool};
use crate::split_sam::{self, Format, Run, SplitKind};

pub use crate::split_sam::file_name_safe;

pub const TOOL: Tool = Tool {
    name: "SplitSamByLibrary",
    summary: "no rows; OUTPUT (a directory) gets one BAM, or SAM text for SAM input, per \
              library of the header, and unknown.bam for the records without a library.",
    options: &[
        Opt::path(
            &["OUTPUT", "O"],
            "The directory the per-library files go to (default the current directory).",
        )
        .written(),
        split_sam::COMPRESSION_OPTION,
        split_sam::INDEX_OPTION,
        split_sam::MD5_OPTION,
        split_sam::RECORDS_IN_RAM_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The format of a library's file: the input's kind. A CRAM input gets
/// CRAM files in Picard, a later wave here.
fn format(extension: &str) -> Result<Format> {
    match extension {
        "sam" => Ok(Format::Sam),
        "cram" => Err(Error::Usage(
            "INPUT is a CRAM, and Picard writes a CRAM per library: bammetrics writes BAM and \
             SAM only: CRAM output is wave W15 of docs/07_full_parity_program.md"
                .to_owned(),
        )),
        _ => Ok(Format::Bam),
    }
}

/// The Picard face: Picard's exit code 2 (`NO_LIBRARIES_SPECIFIED_IN_HEADER`)
/// for a header without a library, with its message; any other error as
/// every tool gives it.
pub fn run_cli(args: &[String]) -> Result<u8> {
    match TOOL.run_cli(args) {
        Ok(()) => Ok(0),
        Err(Error::Usage(message)) if message.starts_with(split_sam::NO_LIBRARIES) => {
            eprintln!("bammetrics: {message}");
            Ok(2)
        }
        Err(e) => Err(e),
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let dir = split_sam::output_dir(inv, Some("."))?;
    if crate::cram::is_cram(&inv.input) {
        format("cram")?;
    }
    split_sam::run(&Run {
        inv,
        kind: SplitKind::Library,
        dir,
        format,
    })?;
    Ok(ToolOutput::NONE)
}

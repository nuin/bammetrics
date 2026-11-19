//! `SplitSamByNumberOfReads`, mirroring Picard 3.5.0
//! (`picard.sam.SplitSamByNumberOfReads`): the records in file order,
//! split over files of about the same number of records in the `OUTPUT`
//! directory.
//!
//! The total is `TOTAL_READS_IN_INPUT`, or the records of a first read of
//! `INPUT` (a pipe without the total is refused, as Picard refuses it).
//! The number of files is `SPLIT_TO_N_FILES`, else `ceil(total /
//! SPLIT_TO_N_READS)`; each file takes `ceil(total / files)` records,
//! then every further record with the name of the last one written, so
//! a query-grouped input keeps each template in one file. Every record
//! counts, secondary and supplementary ones too. Files are
//! `<OUT_PREFIX>_<N>.<extension>`, `N` from 1 printed with at least four
//! digits (`shard_0001.bam`), the extension of the input's kind; the first
//! file is made before the first record, so an empty input gives one
//! header-only file. Picard writes a BAM for a CRAM input, under the
//! name `shard_0001.cram`, and so does bammetrics. A `TOTAL_READS_IN_INPUT`
//! other than the count read changes the sizes and prints Picard's
//! warning.
//!
//! In a chain the step needs `TOTAL_READS_IN_INPUT`: the stream is read
//! once ([`crate::split_sam::chain_options`]).

use crate::error::Result;
use crate::metric::ToolOutput;
use crate::registry::{Invocation, Opt, Tool};
use crate::split_sam::{self, CountOptions, Format, Run, SplitKind};

pub const TOOL: Tool = Tool {
    name: "SplitSamByNumberOfReads",
    summary: "no rows; OUTPUT (a directory) gets the records in file order over files of \
              about the same size, OUT_PREFIX_0001.bam and on, a template never cut.",
    options: &[
        Opt::path(
            &["OUTPUT", "O"],
            "The directory the files go to; it must exist.",
        )
        .written()
        .required(),
        Opt::int(
            &["SPLIT_TO_N_READS", "N_READS"],
            "Records per file, about: a file goes on past it until the record name changes. \
             Not with SPLIT_TO_N_FILES.",
        ),
        Opt::int(
            &["SPLIT_TO_N_FILES", "N_FILES"],
            "The number of files, at least 2. Not with SPLIT_TO_N_READS.",
        ),
        Opt::int(
            &["TOTAL_READS_IN_INPUT", "TOTAL_READS"],
            "The records of INPUT; without it INPUT is read twice, and a pipe is refused.",
        ),
        Opt::text(
            &["OUT_PREFIX"],
            "File names are OUT_PREFIX_N.EXT, EXT the input's kind (default shard).",
        ),
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

/// htsjdk's `makeSAMOrBAMWriter`: SAM text for a `.sam` name, else a BAM.
fn format(extension: &str) -> Result<Format> {
    Ok(if extension == "sam" {
        Format::Sam
    } else {
        Format::Bam
    })
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let count = CountOptions::from_invocation(inv)?;
    let dir = split_sam::output_dir(inv, None)?;
    split_sam::run(&Run {
        inv,
        kind: SplitKind::Count(count),
        dir,
        format,
    })?;
    Ok(ToolOutput::NONE)
}

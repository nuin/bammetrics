//! `BuildBamIndex`, mirroring Picard 3.5.0 `picard.sam.BuildBamIndex` with
//! htsjdk 5.0.0 `BAMIndexer.createIndex`.
//!
//! Picard opens `INPUT` as a BAM, requires `SO:coordinate` in its header
//! (`Input bam file must be sorted by coordinate` for `queryname`,
//! `unsorted`, `unknown` and a header without the tag), reads every record
//! and writes the `.bai` to `OUTPUT`. Without `OUTPUT` the index goes
//! beside the BAM: `x.bai` for `x.bam`, `INPUT.bai` for any other name.
//! [`crate::bai`] holds the format and what htsjdk writes into it; the
//! bytes are Picard's for the same BAM bytes.
//!
//! What bammetrics does differently (KNOWN_DIFFERENCES.md):
//!
//! - Records are not validated. Picard's default `VALIDATION_STRINGENCY`
//!   refuses a record that starts past its reference, has a reference and
//!   no position, or is mapped without a CIGAR; bammetrics indexes it as
//!   Picard does under `VALIDATION_STRINGENCY=SILENT`, and fails only where
//!   the indexer itself fails (a reference out of order, an end past 2^29,
//!   a start in a window past the reference).
//! - A failed run writes nothing: the index is written under a temporary
//!   name and renamed into place. Picard leaves an empty `OUTPUT`.
//! - An empty BGZF block in the middle of the file does not end the read.
//! - `OUTPUT` may not be `INPUT`, and a CRAM is refused by name.
//! - `INPUT` from a pipe needs `OUTPUT`: a pipe has no name to derive one
//!   from.
//!
//! The Python function returns the rows `BamIndexStats` prints for the
//! index written: one per reference, then the `*` row.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::bai::{BamSource, Counts};
use crate::bam_index_stats::{COLUMNS, Row, Stats};
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput};
use crate::pass;
use crate::registry::{self, Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "BuildBamIndex",
    summary: "writes the .bai of a coordinate-sorted BAM, byte for byte the index Picard \
              writes for the same file: OUTPUT, or beside INPUT as x.bai for x.bam. Returns \
              one row per reference of the header (reference, length, aligned, unaligned), \
              then a * row whose unaligned count is the records without a coordinate.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The BAM to index; its header must say SO:coordinate.",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The index file (default: INPUT with .bai in place of .bam, or INPUT.bai for \
             another name).",
        )
        .written(),
        Opt::ignored(
            &["REFERENCE_SEQUENCE", "R"],
            "a BAM is indexed without its reference",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// Where Picard writes the index of `input` without `OUTPUT`: `.bai` in
/// place of a final `.bam`, else appended.
pub fn default_output(input: &Path) -> PathBuf {
    let mut name = input.as_os_str().to_owned();
    if let Some(stem) = input.to_str().and_then(|text| text.strip_suffix(".bam")) {
        name = stem.into();
    }
    name.push(".bai");
    PathBuf::from(name)
}

/// `OUTPUT`, or the default name beside `INPUT`. A device or a pipe has no
/// name to put an index beside.
fn output_path(inv: &Invocation) -> Result<PathBuf> {
    if let Some(output) = inv.get("OUTPUT") {
        return Ok(PathBuf::from(output));
    }
    let regular = std::fs::metadata(&inv.input).is_ok_and(|m| m.is_file());
    if inv.input.starts_with("/dev") || (inv.input.exists() && !regular) {
        return Err(Error::Usage(format!(
            "INPUT={}: a pipe or a device has no name to put the index beside; give OUTPUT",
            inv.input.display()
        )));
    }
    Ok(default_output(&inv.input))
}

/// Index the BAM at `input` into `output`, written under a temporary name
/// and renamed into place. Returns what the index counts.
pub fn build(input: &Path, output: &Path, threads: usize) -> Result<Stats> {
    registry::check_writable([output], &[input.to_path_buf()])?;
    if crate::mark_duplicates::same_file(input, output) {
        return Err(Error::Usage(format!(
            "OUTPUT={}: it is the input {}; writing it would destroy it",
            output.display(),
            input.display()
        )));
    }
    let source = BamSource::open(input, threads)?;
    let order = pass::sort_order(source.header());
    if order != "coordinate" {
        return Err(Error::Metric(format!(
            "{}: Input bam file must be sorted by coordinate (its header says the sort order \
             is {order})",
            input.display()
        )));
    }
    let rows: Vec<Row> = source
        .header()
        .reference_sequences()
        .iter()
        .map(|(name, reference)| Row {
            name: name.to_string(),
            length: usize::from(reference.length()) as i64,
            aligned: 0,
            unaligned: 0,
        })
        .collect();
    let written = |source| Error::Output {
        path: output.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(output);
    let out = BufWriter::with_capacity(1 << 16, file.create().map_err(written)?);
    let (mut out, counts) = source.index(out, output)?;
    out.flush().map_err(written)?;
    drop(out);
    file.commit().map_err(written)?;
    Ok(stats(rows, &counts))
}

/// The rows of the header with the counts the index holds.
fn stats(mut rows: Vec<Row>, counts: &Counts) -> Stats {
    for (row, (aligned, unaligned)) in rows.iter_mut().zip(&counts.references) {
        row.aligned = i64::try_from(*aligned).unwrap_or(i64::MAX);
        row.unaligned = i64::try_from(*unaligned).unwrap_or(i64::MAX);
    }
    Stats {
        rows,
        no_coordinate: i64::try_from(counts.no_coordinate).unwrap_or(i64::MAX),
    }
}

/// The tool writes its own file, so the output holds rows and no text:
/// the registry writer has nothing to print.
fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let output = output_path(inv)?;
    let stats = build(&inv.input, &output, inv.threads)?;
    Ok(MetricsOutput::text(&COLUMNS, stats.values(), String::new()).into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_name_is_picards() {
        for (input, index) in [
            ("x.bam", "x.bai"),
            ("dir/a.b.bam", "dir/a.b.bai"),
            ("/abs/x.BAM", "/abs/x.BAM.bai"),
            ("x.cram", "x.cram.bai"),
            ("x", "x.bai"),
            ("x.bam.gz", "x.bam.gz.bai"),
            (".bam", ".bai"),
        ] {
            assert_eq!(
                default_output(Path::new(input)),
                Path::new(index),
                "{input}"
            );
        }
    }

    #[test]
    fn the_option_table_is_picards() {
        let parse = |args: &[&str]| {
            let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
            TOOL.parse(&args)
        };
        let inv = parse(&[
            "I=x.bam",
            "O=y.bai",
            "VALIDATION_STRINGENCY=SILENT",
            "R=ref.fa",
        ])
        .expect("parses");
        assert_eq!(inv.input, Path::new("x.bam"));
        assert_eq!(inv.get("OUTPUT"), Some("y.bai"));
        assert_eq!(inv.output, None, "the tool writes OUTPUT itself");
        for refused in [
            "ASSUME_SORTED=true",
            "GA4GH_CLIENT_SECRETS=x",
            "STOP_AFTER=1",
        ] {
            assert!(parse(&["I=x.bam", refused]).is_err(), "{refused}");
        }
        assert!(parse(&["O=y.bai"]).is_err(), "INPUT is required");
        assert!(parse(&["I=x.bam", "O=a", "O=b"]).is_err(), "one OUTPUT");
    }
}

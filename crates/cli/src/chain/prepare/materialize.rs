//! A materialized step (docs/07 section 2.3): a tool with neither a stage
//! nor a collector, run on a file of the stream at its position
//! ([`bammetrics_core::stage::materialize`]). Here: the names of its files
//! under the chain's prefix, its file option when the tool has none
//! (`BamIndexStats` and `ViewSam` print to standard output on their own),
//! and the kinds of the files it plans.

use std::path::{Path, PathBuf};

use bammetrics_core::Error;
use bammetrics_core::multiple_metrics::{self, file_extensions};
use bammetrics_core::registry::{Invocation, Kind, Tool};
use bammetrics_core::stage::materialize::{self, Run};

use super::PlannedFile;
use crate::format::{self, Format};

/// The file options of `tool` under a prefix, as (option, extension): the
/// table `qc` names with, else the extension of the tool's materialize row
/// on its `OUTPUT`.
pub(super) fn extensions(tool: &Tool) -> Vec<(&'static str, &'static str)> {
    if let Some(table) = file_extensions(tool.name) {
        return table.to_vec();
    }
    materialize::row_of(tool)
        .and_then(|row| row.extension)
        .map(|extension| vec![("OUTPUT", extension)])
        .unwrap_or_default()
}

/// True for a tool whose own `OUTPUT` is refused: its text goes to
/// standard output when it runs alone, and the chain names its file.
fn has_no_output(tool: &Tool) -> bool {
    !matches!(
        tool.option("OUTPUT").map(|opt| opt.kind),
        Some(Kind::Output | Kind::Value)
    )
}

/// `argv` without the `OUTPUT=FILE` the chain named for a tool that takes
/// none, and that file.
pub(super) fn take_output(tool: &Tool, argv: &mut Vec<String>) -> Option<PathBuf> {
    if !has_no_output(tool) {
        return None;
    }
    let at = argv.iter().position(|arg| arg.starts_with("OUTPUT="))?;
    let arg = argv.remove(at);
    arg.strip_prefix("OUTPUT=").map(PathBuf::from)
}

/// The files of a materialized step: a file the chain names for a tool
/// without `OUTPUT` joins the planned list, and a file the tool writes
/// itself (a checksum, SAM text) has one form, so it is planned as `file`
/// with no sibling per further format.
pub(super) fn plan_files(
    tool: &Tool,
    output: Option<&PathBuf>,
    files: &mut Vec<PlannedFile>,
    formats: &mut Vec<Format>,
) {
    let Some(row) = materialize::row_of(tool) else {
        return;
    };
    if let Some(path) = output {
        files.insert(
            0,
            PlannedFile {
                path: path.clone(),
                kind: "metrics",
            },
        );
    }
    if matches!(row.run, Run::WritesOutput | Run::View) {
        files.retain(|file| file.kind == "metrics");
        for file in files.iter_mut() {
            file.kind = "file";
        }
        formats.truncate(1);
    }
}

/// What a materialized step is checked against before the first read.
pub(super) struct Before<'a> {
    /// The chain's source as given.
    pub source: &'a str,
    /// The step reads the source's stream with no step in front.
    pub reads_source: bool,
    pub has_reference: bool,
    /// The chain's `--tmp-dir`.
    pub tmp_dir: Option<&'a Path>,
}

/// The checks of a materialized step before the first read, on its Picard
/// arguments `argv`: the tool's own option parser
/// ([`materialize::check`]), each further input
/// ([`materialize::check_further_input`]), the index of a source it reads
/// directly ([`materialize::check_source`]), and the `--tmp-dir` its copy
/// goes to.
pub(super) fn check(tool: &Tool, argv: &[String], before: &Before) -> Result<(), Error> {
    materialize::check(tool, argv)?;
    let row = materialize::row_of(tool);
    if row.is_some_and(|row| row.several_inputs) {
        for input in argv.iter().skip(1).filter_map(|a| a.strip_prefix("INPUT=")) {
            // A gather list is read by the tool; only a named BAM, SAM or
            // CRAM is checked here.
            let list = tool.name == "GatherBamFiles"
                && ![".bam", ".sam", ".cram"].iter().any(|e| input.ends_with(e));
            if !list && input != "-" {
                materialize::check_further_input(Path::new(input), before.has_reference)?;
            }
        }
    }
    let source = Path::new(before.source);
    let regular = std::fs::metadata(source).is_ok_and(|m| m.is_file());
    if before.reads_source && regular && !source.starts_with("/dev") {
        materialize::check_source(tool, source)?;
    }
    match before.tmp_dir.map(unusable_tmp_dir) {
        Some(Some(why)) => Err(Error::Usage(format!(
            "--tmp-dir {}: {why}; the step's copy of its stream and its spill go there",
            before.tmp_dir.unwrap_or(Path::new("")).display()
        ))),
        _ => Ok(()),
    }
}

/// Why no temporary directory can be made under `dir`; `None` when one
/// can. A missing `dir` is made by the run (`create_dir_all`), so the
/// nearest directory that exists must be one the user may write.
fn unusable_tmp_dir(dir: &Path) -> Option<String> {
    use std::os::unix::ffi::OsStrExt;
    let existing = dir
        .ancestors()
        .find(|path| !path.as_os_str().is_empty() && std::fs::symlink_metadata(path).is_ok());
    let existing = existing.unwrap_or(Path::new("."));
    if !std::fs::metadata(existing).is_ok_and(|m| m.is_dir()) {
        return Some(format!("{} is not a directory", existing.display()));
    }
    let text = std::ffi::CString::new(existing.as_os_str().as_bytes()).ok()?;
    // SAFETY: `text` is a NUL-terminated path that outlives the call.
    if unsafe { libc::access(text.as_ptr(), libc::W_OK | libc::X_OK) } == 0 {
        return None;
    }
    let cause = std::io::Error::last_os_error();
    Some(format!("{}: {cause}", existing.display()))
}

/// The metrics files `CollectMultipleMetrics` names under its prefix
/// before it reads ([`multiple_metrics::planned_files`]), each with a
/// sibling per further format, so a tap that writes one of them is a
/// usage error before the first read.
pub(super) fn program_files(
    tool: &Tool,
    inv: &Invocation,
    formats: &[Format],
) -> Result<Vec<PlannedFile>, Error> {
    if tool.name != multiple_metrics::TOOL.name {
        return Ok(Vec::new());
    }
    let mut files = Vec::new();
    for path in multiple_metrics::planned_files(inv)? {
        let siblings: Vec<PlannedFile> = formats[1..]
            .iter()
            .flat_map(|f| {
                format::files_of(&format::sibling(&path, f.extension()), *f)
                    .into_iter()
                    .map(|path| PlannedFile {
                        path,
                        kind: f.name(),
                    })
            })
            .collect();
        files.push(PlannedFile {
            path,
            kind: "metrics",
        });
        files.extend(siblings);
    }
    Ok(files)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tool(name: &str) -> &'static Tool {
        bammetrics_core::registry::find(name).expect(name)
    }

    #[test]
    fn a_tool_without_a_prefix_name_takes_its_rows_extension() {
        assert_eq!(
            extensions(tool("CollectRawWgsMetrics")),
            vec![("OUTPUT", ".raw_wgs_metrics")]
        );
        assert_eq!(
            extensions(tool("CollectWgsMetrics")),
            vec![("OUTPUT", ".wgs_metrics")]
        );
        assert!(extensions(tool("CollectInsertSizeMetrics")).len() == 1);
    }

    #[test]
    fn the_output_of_a_tool_that_prints_is_taken_out_of_its_arguments() {
        let mut argv = vec!["INPUT=x.bam".to_owned(), "OUTPUT=s.stats".to_owned()];
        let taken = take_output(tool("BamIndexStats"), &mut argv);
        assert_eq!(taken, Some(PathBuf::from("s.stats")));
        assert_eq!(argv, ["INPUT=x.bam"]);
        let mut argv = vec!["INPUT=x.bam".to_owned(), "OUTPUT=s.txt".to_owned()];
        assert_eq!(take_output(tool("CollectRawWgsMetrics"), &mut argv), None);
        assert_eq!(argv.len(), 2);
    }
}

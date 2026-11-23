//! A materialized step (docs/07 section 2.3): a tool with neither a stage
//! nor a collector, run on a file of the stream at its position
//! ([`bammetrics_core::stage::materialize`]). Here: the names of its files
//! under the chain's prefix, its file option when the tool has none
//! (`BamIndexStats` and `ViewSam` print to standard output on their own),
//! and the kinds of the files it plans.

use std::path::PathBuf;

use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::registry::{Kind, Tool};
use bammetrics_core::stage::materialize::{self, Run};

use super::PlannedFile;
use crate::format::Format;

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

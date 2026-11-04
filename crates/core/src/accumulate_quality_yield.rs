//! `AccumulateQualityYieldMetrics`, mirroring Picard 3.5.0
//! (`picard.util.AccumulateQualityYieldMetrics`).
//!
//! Combines `QualityYieldMetrics` files written for shards of one read group.
//! Each `INPUT` contributes the first row of its metrics block; every field
//! is added except `READ_LENGTH`, which is derived again from the sums
//! (`TOTAL_BASES / TOTAL_READS`, integer division). The two
//! `Q20_EQUIVALENT_YIELD` fields are added as printed, so the result can be
//! lower than a single run over the whole input would print. Picard writes
//! the file without `StringHeader` lines; bammetrics does the same.

use std::collections::HashMap;
use std::path::Path;

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput};
use crate::metrics_file::MetricsDoc;
use crate::quality_yield::{CLASS, COLUMNS, Row};
use crate::registry::{Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "AccumulateQualityYieldMetrics",
    summary: "QualityYieldMetrics files (a list of paths) merged into one dict.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A QualityYieldMetrics file; repeat for each.",
        )
        .repeated()
        .required(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut total = Row::default();
    for input in inv.all("INPUT") {
        total.merge(&read_row(Path::new(input))?);
    }
    let mut output = MetricsOutput::from_rows(&[total], None);
    output.with_headers = false;
    Ok(output.into())
}

/// The first `QualityYieldMetrics` row of a Picard metrics file
/// (`MetricsFile.readBeans(file).get(0)`). A column the class does not have,
/// a value that is not an integer, or a file without a row is an error, as
/// in htsjdk; a column the file lacks stays 0.
pub fn read_row(path: &Path) -> Result<Row> {
    let doc = MetricsDoc::read(path)?;
    let bad = |message: String| Error::Metric(format!("{}: {message}", path.display()));
    let block = doc
        .classes
        .first()
        .ok_or_else(|| bad("no ## METRICS CLASS line".to_owned()))?;
    if block.class != CLASS {
        return Err(bad(format!("metrics class {} is not {CLASS}", block.class)));
    }
    let values = block
        .rows
        .first()
        .ok_or_else(|| bad("no metrics row".to_owned()))?;
    let mut fields = HashMap::new();
    for (column, value) in block.columns.iter().zip(values) {
        if !COLUMNS.contains(&column.as_str()) {
            return Err(bad(format!(
                "column {column} is not a QualityYieldMetrics field"
            )));
        }
        let value = value.as_int().ok_or_else(|| {
            bad(format!(
                "{column}={:?} is not an integer",
                value.picard_text()
            ))
        })?;
        fields.insert(column.as_str(), value);
    }
    let field = |name: &str| fields.get(name).copied().unwrap_or(0);
    Ok(Row {
        total_reads: field("TOTAL_READS"),
        pf_reads: field("PF_READS"),
        read_length: field("READ_LENGTH") as i32,
        total_bases: field("TOTAL_BASES"),
        pf_bases: field("PF_BASES"),
        q20_bases: field("Q20_BASES"),
        pf_q20_bases: field("PF_Q20_BASES"),
        q30_bases: field("Q30_BASES"),
        pf_q30_bases: field("PF_Q30_BASES"),
        q20_equivalent_yield: field("Q20_EQUIVALENT_YIELD"),
        pf_q20_equivalent_yield: field("PF_Q20_EQUIVALENT_YIELD"),
    })
}

#[cfg(test)]
mod tests {
    use std::fs;

    use super::*;

    fn write(dir: &Path, name: &str, body: &str) -> std::path::PathBuf {
        let path = dir.join(name);
        fs::write(&path, body).expect("write");
        path
    }

    fn file(row: &str) -> String {
        format!(
            "## htsjdk.samtools.metrics.StringHeader\n# x\n\n## METRICS CLASS\t{CLASS}\n{}\n{row}\n\n",
            COLUMNS.join("\t")
        )
    }

    #[test]
    fn rows_add_and_read_length_is_derived() {
        let dir = std::env::temp_dir().join(format!("bm_acc_{}", std::process::id()));
        fs::create_dir_all(&dir).expect("dir");
        let a = read_row(&write(
            &dir,
            "a.txt",
            &file("2\t2\t100\t200\t200\t10\t10\t5\t5\t7\t7"),
        ))
        .expect("a");
        let b = read_row(&write(
            &dir,
            "b.txt",
            &file("1\t0\t75\t75\t0\t1\t0\t1\t0\t3\t0"),
        ))
        .expect("b");
        let mut total = Row::default();
        total.merge(&a);
        total.merge(&b);
        assert_eq!(
            (total.total_reads, total.total_bases, total.read_length),
            (3, 275, 91)
        );
        assert_eq!(total.q20_equivalent_yield, 10);
        let other = write(&dir, "c.txt", "## METRICS CLASS\tpicard.Other\nA\n1\n");
        assert!(read_row(&other).is_err());
        let empty = write(&dir, "d.txt", &file("").replace("\n\n\n", "\n"));
        assert!(read_row(&empty).is_err());
        fs::remove_dir_all(&dir).ok();
    }
}

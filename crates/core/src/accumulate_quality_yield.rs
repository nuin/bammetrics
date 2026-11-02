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
use std::fs;
use std::path::Path;

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput};
use crate::quality_yield::{CLASS, COLUMNS, Row};
use crate::registry::{Invocation, Tool};

pub const TOOL: Tool = Tool {
    name: "AccumulateQualityYieldMetrics",
    options: &[],
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
    let text = fs::read_to_string(path).map_err(|source| Error::Input {
        path: path.to_path_buf(),
        source,
    })?;
    let bad = |message: String| Error::Metric(format!("{}: {message}", path.display()));
    let mut lines = text.lines();
    let class = lines
        .by_ref()
        .find_map(|l| l.strip_prefix("## METRICS CLASS\t"))
        .ok_or_else(|| bad("no ## METRICS CLASS line".to_owned()))?;
    if class.trim() != CLASS {
        return Err(bad(format!("metrics class {class} is not {CLASS}")));
    }
    let columns: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let values: Vec<&str> = match lines.next() {
        Some(line) if !line.is_empty() && !line.starts_with('#') => line.split('\t').collect(),
        _ => return Err(bad("no metrics row".to_owned())),
    };
    let mut fields = HashMap::new();
    for (i, &column) in columns.iter().enumerate() {
        if !COLUMNS.contains(&column) {
            return Err(bad(format!(
                "column {column} is not a QualityYieldMetrics field"
            )));
        }
        let text = values.get(i).copied().unwrap_or_default();
        let value: i64 = text
            .parse()
            .map_err(|_| bad(format!("{column}={text:?} is not an integer")))?;
        fields.insert(column, value);
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

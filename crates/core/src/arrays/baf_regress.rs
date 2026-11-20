//! `CreateBafRegressMetricsFile`: the text bafRegress prints, one
//! `BafRegressMetrics` row per sample line, in `OUTPUT.bafregress_metrics`.
//!
//! As Picard 3.5.0 reads it: the first line must be exactly
//! `sample<TAB>estimate<TAB>stderr<TAB>tval<TAB>pval<TAB>callrate<TAB>Nhom`
//! (an empty file fails there). Every later line, an empty one included,
//! is split at runs of Java's `\s` (`split("\\s+")`: a line that starts
//! with a space has an empty first field, trailing spaces add none) and
//! must give exactly 7 fields: the sample as it stands, five numbers read
//! by `Double.parseDouble` (`NaN`, `Infinity`, `1e-5`, `0x1p3`, a `d`
//! suffix) and `Nhom` by `Integer.parseInt`. `LOG10_PVAL` is `log10(pval)`.
//! A header alone gives a file with no row. Any failure writes nothing.

use std::path::{Path, PathBuf};

use super::{content_error, read_lines, read_text, split_java_spaces};
use crate::error::{Error, Result};
use crate::java;
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::variant_calling_metrics::output_path;

/// The file extension Picard adds to `OUTPUT`.
pub const EXTENSION: &str = "bafregress_metrics";

/// The one header line Picard accepts.
pub const HEADER: &str = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom";

/// `BafRegressMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub sample: String,
    pub estimate: f64,
    pub stderr: f64,
    pub tval: f64,
    pub pval: f64,
    pub log10_pval: f64,
    pub call_rate: f64,
    pub nhom: i64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.arrays.BafRegressMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SAMPLE",
        "ESTIMATE",
        "STDERR",
        "TVAL",
        "PVAL",
        "LOG10_PVAL",
        "CALL_RATE",
        "NHOM",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sample.as_str().into(),
            self.estimate.into(),
            self.stderr.into(),
            self.tval.into(),
            self.pval.into(),
            self.log10_pval.into(),
            self.call_rate.into(),
            self.nhom.into(),
        ]
    }
}

/// `CreateBafRegressMetricsFile` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CreateBafRegressMetricsFile",
    summary: "one BafRegressMetrics row per sample line of the bafRegress output INPUT; the \
              file is OUTPUT.bafregress_metrics. Python returns one row per file: file and rows.",
    options: &[
        Opt::path(&["INPUT", "I"], "The text bafRegress printed.").required(),
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix: the file is OUTPUT.bafregress_metrics.")
            .required(),
        crate::registry::STRINGENCY_NO_RECORDS,
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

fn prefix(inv: &Invocation) -> Result<PathBuf> {
    inv.output
        .clone()
        .ok_or_else(|| Error::Usage("Option 'OUTPUT' is required.".to_owned()))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let prefix = prefix(inv)?;
    let rows = read(&inv.input)?;
    Ok(ToolOutput {
        main: None,
        extra: vec![(
            output_path(&prefix, EXTENSION),
            MetricsOutput::from_rows(&rows, None),
        )],
    })
}

/// Every row of the bafRegress output at `path`.
pub fn read(path: &Path) -> Result<Vec<Row>> {
    let text = read_text(path)?;
    let lines = read_lines(path, &text)?;
    parse(&lines).map_err(|message| content_error(path, message))
}

/// The rows of the lines of a bafRegress output; the error names the line.
pub fn parse(lines: &[&str]) -> std::result::Result<Vec<Row>, String> {
    let Some((header, body)) = lines.split_first() else {
        return Err(
            "no header line (Picard: NullPointerException reading the first line)".to_owned(),
        );
    };
    if *header != HEADER {
        return Err(format!(
            "the first line is not the bafRegress header {HEADER:?} (Picard: Unrecognized \
             header line)"
        ));
    }
    body.iter()
        .enumerate()
        .map(|(i, line)| row(line).map_err(|e| format!("line {}: {e}", i + 2)))
        .collect()
}

fn row(line: &str) -> std::result::Result<Row, String> {
    let fields = split_java_spaces(line);
    if fields.len() != 7 {
        return Err(format!(
            "Invalid number of entries ({}) in line: {line} (Picard: Error parsing bafRegress \
             Output)",
            fields.len()
        ));
    }
    let number = |i: usize| java::parse_double(fields[i]);
    let pval = number(4)?;
    Ok(Row {
        sample: fields[0].to_owned(),
        estimate: number(1)?,
        stderr: number(2)?,
        tval: number(3)?,
        pval,
        log10_pval: pval.log10(),
        call_rate: number(5)?,
        nhom: java::parse_integer(fields[6], 32)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_parse_as_picard_reads_them() {
        let lines = [
            HEADER,
            "S1\t0.0012\t0.0003\t4.1\t2.1e-5\t0.995\t123456",
            " 1 2 3 1e2 5 +6  ",
        ];
        let rows = parse(&lines).unwrap();
        assert_eq!(rows[0].sample, "S1");
        assert_eq!(rows[0].nhom, 123_456);
        assert_eq!(rows[1].sample, "");
        assert_eq!(rows[1].log10_pval, 2.0);
        assert_eq!(rows[1].nhom, 6);
    }

    #[test]
    fn bad_lines_fail() {
        assert!(parse(&[]).is_err());
        assert!(parse(&["sample estimate"]).is_err());
        assert!(parse(&[HEADER, ""]).is_err());
        assert!(parse(&[HEADER, "S 1 2 3 4 5 6 7"]).is_err());
        assert!(parse(&[HEADER, "S 1 2 3 4 5 6.5"]).is_err());
        assert!(parse(&[HEADER, "S 1 2 3 x 5 6"]).is_err());
        assert!(parse(&[HEADER]).unwrap().is_empty());
    }
}

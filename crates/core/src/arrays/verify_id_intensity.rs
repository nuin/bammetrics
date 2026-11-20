//! `CreateVerifyIDIntensityContaminationMetricsFile`: the text
//! VerifyIDIntensity prints, one `VerifyIDIntensityContaminationMetrics`
//! row per result line, in `OUTPUT.verifyidintensity_metrics`.
//!
//! As Picard 3.5.0 reads it, each line matched whole against a Java
//! pattern (`\s` is Java's six ASCII spaces, `\d` an ASCII digit):
//!
//! 1. `ID\s+%Mix\s+LLK\s+LLK0\s*`, the header;
//! 2. `-+`, the rule under it;
//! 3. every later line `(\d+)\s+([0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s+([-0-9]*\.?[0-9]+)\s*`:
//!    `ID` by `Integer.parseInt`, `%Mix`, `LLK` and `LLK0` by
//!    `Double.parseDouble`. The last two take any run of `-` and digits the
//!    pattern lets through, so `1-2` passes the pattern and fails the
//!    number, as in Picard.
//!
//! A file shorter than two lines fails (Picard's `NullPointerException`);
//! a header and a rule alone give a file with no row. An empty line, a
//! leading space or any other text fails the pattern. Any failure writes
//! nothing.

use std::path::{Path, PathBuf};

use super::{content_error, is_java_space, read_lines, read_text};
use crate::error::{Error, Result};
use crate::java;
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::variant_calling_metrics::output_path;

/// The file extension Picard adds to `OUTPUT`.
pub const EXTENSION: &str = "verifyidintensity_metrics";

/// `VerifyIDIntensityContaminationMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub id: i64,
    pub pct_mix: f64,
    pub llk: f64,
    pub llk0: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.arrays.VerifyIDIntensityContaminationMetrics";
    const COLUMNS: &'static [&'static str] = &["ID", "PCT_MIX", "LLK", "LLK0"];

    fn values(&self) -> Vec<Value> {
        vec![
            self.id.into(),
            self.pct_mix.into(),
            self.llk.into(),
            self.llk0.into(),
        ]
    }
}

/// `CreateVerifyIDIntensityContaminationMetricsFile` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CreateVerifyIDIntensityContaminationMetricsFile",
    summary: "one VerifyIDIntensityContaminationMetrics row per result line of the \
              VerifyIDIntensity output INPUT; the file is OUTPUT.verifyidintensity_metrics. \
              Python returns one row per file: file and rows.",
    options: &[
        Opt::path(&["INPUT", "I"], "The text VerifyIDIntensity printed.").required(),
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix: the file is OUTPUT.verifyidintensity_metrics.")
            .required(),
        crate::registry::STRINGENCY_NO_RECORDS,
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let prefix: PathBuf = inv
        .output
        .clone()
        .ok_or_else(|| Error::Usage("Option 'OUTPUT' is required.".to_owned()))?;
    let rows = read(&inv.input)?;
    Ok(ToolOutput {
        main: None,
        extra: vec![(
            output_path(&prefix, EXTENSION),
            MetricsOutput::from_rows(&rows, None),
        )],
    })
}

/// Every row of the VerifyIDIntensity output at `path`.
pub fn read(path: &Path) -> Result<Vec<Row>> {
    let text = read_text(path)?;
    let lines = read_lines(path, &text)?;
    parse(&lines).map_err(|message| content_error(path, message))
}

/// The rows of the lines of a VerifyIDIntensity output; the error names
/// the line.
pub fn parse(lines: &[&str]) -> std::result::Result<Vec<Row>, String> {
    let mismatch = |number: usize, line: &str| {
        format!(
            "line {number}: {line:?} does not match the VerifyIDIntensity format (Picard: \
             Unrecognized line)"
        )
    };
    let missing = |number: usize| {
        format!("line {number} is missing (Picard: NullPointerException reading the file)")
    };
    let header = lines.first().ok_or_else(|| missing(1))?;
    if !is_header(header) {
        return Err(mismatch(1, header));
    }
    let rule = lines.get(1).ok_or_else(|| missing(2))?;
    if rule.is_empty() || !rule.bytes().all(|b| b == b'-') {
        return Err(mismatch(2, rule));
    }
    lines[2..]
        .iter()
        .enumerate()
        .map(|(i, line)| {
            let number = i + 3;
            let fields = result_fields(line).ok_or_else(|| mismatch(number, line))?;
            row(&fields).map_err(|e| format!("line {number}: {e}"))
        })
        .collect()
}

fn row(fields: &[&str; 4]) -> std::result::Result<Row, String> {
    Ok(Row {
        id: java::parse_integer(fields[0], 32)?,
        pct_mix: java::parse_double(fields[1])?,
        llk: java::parse_double(fields[2])?,
        llk0: java::parse_double(fields[3])?,
    })
}

/// The fields between runs of Java spaces, with no leading space and any
/// trailing ones.
fn tokens(line: &str) -> Option<Vec<&str>> {
    if line.starts_with(is_java_space) {
        return None;
    }
    Some(
        line.split(is_java_space)
            .filter(|t| !t.is_empty())
            .collect(),
    )
}

/// `^ID\s+%Mix\s+LLK\s+LLK0\s*$`: the four words with spaces between
/// them and after them only.
fn is_header(line: &str) -> bool {
    tokens(line).is_some_and(|t| t == ["ID", "%Mix", "LLK", "LLK0"])
}

/// The four groups of a result line, or `None` when the pattern fails.
fn result_fields(line: &str) -> Option<[&str; 4]> {
    let tokens = tokens(line)?;
    let [id, mix, llk, llk0] = tokens.as_slice() else {
        return None;
    };
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    // `[P]*\.?[0-9]+`: an optional point, digits after it, and before it
    // only characters of the class P.
    let decimal = |s: &str, prefix: fn(u8) -> bool| match s.split_once('.') {
        Some((before, after)) => before.bytes().all(prefix) && digits(after),
        None => s.bytes().all(prefix) && s.bytes().last().is_some_and(|b| b.is_ascii_digit()),
    };
    let unsigned = |b: u8| b.is_ascii_digit();
    let signed = |b: u8| b.is_ascii_digit() || b == b'-';
    (digits(id) && decimal(mix, unsigned) && decimal(llk, signed) && decimal(llk0, signed))
        .then_some([*id, *mix, *llk, *llk0])
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEAD: &str = "ID  %Mix   LLK       LLK0";

    #[test]
    fn rows_parse_as_picard_reads_them() {
        let rows = parse(&[
            HEAD,
            "-----",
            "0   0.0034  -1234.5  -1240.25",
            "12\t10 -.5 5 ",
        ])
        .unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].llk0, -1240.25);
        assert_eq!(rows[1].id, 12);
        assert_eq!(rows[1].llk, -0.5);
        assert!(parse(&[HEAD, "-"]).unwrap().is_empty());
    }

    #[test]
    fn lines_outside_the_patterns_fail() {
        for lines in [
            &[][..],
            &[HEAD][..],
            &["ID %Mix LLK", "--"][..],
            &[HEAD, "--x"][..],
            &[HEAD, ""][..],
            &[HEAD, "--", ""][..],
            &[HEAD, "--", " 1 2 3 4"][..],
            &[HEAD, "--", "1 2. 3 4"][..],
            &[HEAD, "--", "1 -2 3 4"][..],
            &[HEAD, "--", "1 2 3 4 5"][..],
            &[HEAD, "--", "1 2 1-2 4"][..],
            &[HEAD, "--", "1 2 3 1e5"][..],
            &[HEAD, "--", "99999999999 2 3 4"][..],
        ] {
            assert!(parse(lines).is_err(), "{lines:?}");
        }
    }
}

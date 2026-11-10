//! `--format` of the tool verbs (docs/05 section 2.2, docs/formats.md):
//! the Picard metrics file, or the same table in a format of
//! [`bammetrics_formats`].
//!
//! `--format` takes one name or a list (`picard,json`). The first format
//! goes to `-o` (or standard output); each further one is a sibling file
//! named `-o` plus [`Format::extension`]. Every format but `picard` is
//! written from one typed [`Document`], so a column has the same type in
//! each.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use bammetrics_core::metric::MetricsOutput;
use bammetrics_core::metrics_file::format_double;
use bammetrics_core::registry::did_you_mean;
use bammetrics_formats::{CheckLine, Document, ListError, Outcome, VerdictSummary};
use bammetrics_verdict::{Status, Verdict};

pub use bammetrics_formats::Format;

use crate::error::VerbError;

/// What a verb prints for one tool's output: the tool, the `StringHeader`
/// command line and the sample (`--sample`, or its default).
pub struct Context<'a> {
    pub tool: &'a str,
    pub command: &'a str,
    pub sample: Option<&'a str>,
}

/// The formats a `--format` value names, the first for the main output. A
/// name that is no format, or one given twice, is a usage error of `verb`
/// that names the flag as `spelled`.
pub fn parse_list(verb: &str, value: &str, spelled: &str) -> Result<Vec<Format>, VerbError> {
    let names = bammetrics_formats::names();
    bammetrics_formats::parse_list(value).map_err(|e| {
        let error = match e {
            ListError::Unknown(name) => {
                let error = VerbError::usage(verb, format!("no format {name:?}")).from(spelled);
                match did_you_mean(&name, names.iter().copied()) {
                    Some(best) => error.hint(format!("did you mean {best}?")),
                    None => error,
                }
            }
            ListError::Repeated(format) => {
                VerbError::usage(verb, format!("format {} is given twice", format.name()))
                    .from(spelled)
            }
        };
        error.hint(format!("formats: {}", names.join(", ")))
    })
}

/// `path` with `.` and the format's extension added: the sibling a
/// further format is written to.
pub fn sibling(path: &Path, extension: &str) -> PathBuf {
    let mut name = OsString::from(path.as_os_str());
    name.push(".");
    name.push(extension);
    PathBuf::from(name)
}

/// `output` in each of `formats`, with one set of `StringHeader` lines, so
/// the Picard file and its JSON sibling carry the same start time. The
/// Picard form is the file the registry writer writes.
pub fn render_all(
    formats: &[Format],
    context: &Context<'_>,
    output: &MetricsOutput,
) -> Vec<String> {
    let headers = output.headers(context.command);
    let document = formats.iter().any(|f| *f != Format::Picard).then(|| {
        let mut doc = Document::from_output(context.tool, headers.clone(), output);
        doc.sample = context.sample.map(str::to_owned);
        doc
    });
    formats
        .iter()
        .map(|format| match (format, &document) {
            (Format::Picard, _) | (_, None) => {
                let mut bytes = Vec::new();
                output
                    .write(headers.clone(), &mut bytes)
                    .expect("writing to memory");
                String::from_utf8_lossy(&bytes).into_owned()
            }
            (format, Some(doc)) => bammetrics_formats::render(*format, doc, None),
        })
        .collect()
}

/// What the Markdown, JUnit and MultiQC writers need of `verdict`.
pub fn summary(verdict: &Verdict) -> VerdictSummary {
    let preset = &verdict.preset;
    VerdictSummary {
        sample: verdict.sample.clone(),
        status: verdict.status.name().to_owned(),
        preset: format!(
            "{} v{} ({})",
            preset.name,
            preset.version,
            preset.hash.get(..6).unwrap_or(&preset.hash)
        ),
        checks: verdict
            .checks
            .iter()
            .map(|check| CheckLine {
                rule: check.rule.clone(),
                label: check.label.clone(),
                metric: format!(
                    "{} {}",
                    check
                        .class
                        .rsplit(['.', '$'])
                        .next()
                        .unwrap_or(&check.class),
                    check.metric
                ),
                outcome: match check.status {
                    Status::Pass => Outcome::Pass,
                    Status::Warn => Outcome::Warn,
                    Status::Fail => Outcome::Fail,
                    Status::NotMeasured => Outcome::NotMeasured,
                },
                value: check.value,
                shown: shown_value(check.value, check.text.as_deref()),
                warn: check.warn.clone(),
                fail: check.fail.clone(),
                reason: check.reason.clone(),
            })
            .collect(),
    }
}

/// The value of a check as `check` prints it: the number, the cell of a
/// text rule (an empty cell as `empty`, as the page does), `?` when not
/// measured.
pub fn shown_value(value: Option<f64>, text: Option<&str>) -> String {
    match (value, text) {
        (Some(v), _) => format_double(v),
        (None, Some("")) => "empty".to_owned(),
        (None, Some(text)) => text.to_owned(),
        (None, None) => "?".to_owned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bammetrics_core::metric::Value;
    use bammetrics_formats::json::Json;

    fn output() -> MetricsOutput {
        MetricsOutput {
            class: "picard.analysis.InsertSizeMetrics",
            columns: &["MEDIAN_INSERT_SIZE", "PAIR_ORIENTATION", "SAMPLE"],
            rows: vec![vec![
                Value::Float(310.0),
                Value::Str("FR".to_owned()),
                Value::Null,
            ]],
            histogram: None,
            with_headers: true,
            text: None,
        }
    }

    fn render(format: Format, context: &Context<'_>, output: &MetricsOutput) -> String {
        render_all(&[format], context, output).remove(0)
    }

    fn context(sample: Option<&str>) -> Context<'_> {
        Context {
            tool: "CollectInsertSizeMetrics",
            command: "CollectInsertSizeMetrics INPUT=x.bam",
            sample,
        }
    }

    #[test]
    fn json_forms_type_the_cells() {
        let json = render(Format::Json, &context(Some("S1")), &output());
        let parsed = Json::parse(&json).expect("valid json");
        let metrics = parsed
            .get("metrics")
            .and_then(Json::as_array)
            .expect("rows");
        assert_eq!(
            metrics[0].get("MEDIAN_INSERT_SIZE"),
            Some(&Json::Number("310".to_owned()))
        );
        assert_eq!(metrics[0].get("SAMPLE"), Some(&Json::Null));
        assert_eq!(parsed.get("sample"), Some(&Json::text("S1")));
        assert_eq!(
            parsed.get("schema"),
            Some(&Json::text("bammetrics.metrics/1"))
        );
        assert_eq!(
            parsed.get("types"),
            Some(&Json::strings(["float", "string", "string"]))
        );
        assert_eq!(
            render(Format::Jsonl, &context(None), &output()),
            "{\"MEDIAN_INSERT_SIZE\":310,\"PAIR_ORIENTATION\":\"FR\",\"SAMPLE\":null}\n"
        );
    }

    #[test]
    fn text_forms_use_picard_cells() {
        assert_eq!(
            render(Format::Tsv, &context(None), &output()),
            "MEDIAN_INSERT_SIZE\tPAIR_ORIENTATION\tSAMPLE\n310\tFR\t\n"
        );
        assert_eq!(
            render(Format::Human, &context(None), &output()),
            "MEDIAN_INSERT_SIZE  PAIR_ORIENTATION  SAMPLE\n               310                FR        \n"
        );
        assert_eq!(
            render(Format::Csv, &context(None), &output()),
            "MEDIAN_INSERT_SIZE,PAIR_ORIENTATION,SAMPLE\r\n310,FR,\r\n"
        );
    }

    #[test]
    fn a_list_shares_one_start_time() {
        let all = render_all(&[Format::Picard, Format::Json], &context(None), &output());
        let started = all[0]
            .lines()
            .find_map(|l| l.strip_prefix("# Started on: "))
            .expect("a Started on header");
        assert!(all[1].contains(&format!("\"Started on: {started}\"")));
        let back = bammetrics_formats::read(&all[1]).expect("the JSON reads");
        assert_eq!(
            bammetrics_formats::render(Format::Picard, &back, None),
            all[0]
        );
    }

    #[test]
    fn format_lists_name_the_bad_word() {
        assert_eq!(
            parse_list("insert-size", "picard,json", "--format").expect("parses"),
            [Format::Picard, Format::Json]
        );
        let unknown = parse_list("insert-size", "tsv,jsn", "--format").expect_err("jsn");
        assert_eq!(
            unknown.to_string(),
            "bammetrics insert-size: no format \"jsn\" (from --format)\n  did you mean json?\n  \
             formats: picard, json, tsv, jsonl, human, csv, markdown"
        );
        let twice = parse_list("insert-size", "tsv,tsv", "--format").expect_err("twice");
        assert!(twice.to_string().contains("format tsv is given twice"));
        assert_eq!(
            sibling(Path::new("out/x.txt"), "json"),
            PathBuf::from("out/x.txt.json")
        );
    }
}

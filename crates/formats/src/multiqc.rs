//! MultiQC custom content: `PREFIX.verdict_mqc.json`, a table section
//! with one row per sample. MultiQC reads any `*_mqc.json` under the
//! directories it is given, so a `qc --also multiqc` run shows up in its
//! report with no plugin. The Picard-class metrics files need nothing:
//! MultiQC's Picard module parses them as Picard's own.
//!
//! Columns: the verdict, the count of failed, warned and not-measured
//! checks, then one column per rule with the value judged (empty when not
//! measured).

use crate::json::Json;
use crate::verdict::{Outcome, VerdictSummary};

fn object(members: Vec<(&str, Json)>) -> Json {
    Json::Object(
        members
            .into_iter()
            .map(|(key, value)| (key.to_owned(), value))
            .collect(),
    )
}

fn count_header(title: &str, description: &str) -> Json {
    object(vec![
        ("title", Json::text(title)),
        ("description", Json::text(description)),
        ("format", Json::text("{:,d}")),
        ("scale", Json::Bool(false)),
    ])
}

/// The custom content file for `verdicts`.
pub fn write(verdicts: &[VerdictSummary]) -> String {
    let equals = |word: &str| Json::Array(vec![object(vec![("s_eq", Json::text(word))])]);
    let mut headers = vec![
        (
            "verdict".to_owned(),
            object(vec![
                ("title", Json::text("Verdict")),
                (
                    "description",
                    Json::text("pass, warn, fail, or incomplete when a rule could not be measured"),
                ),
                ("scale", Json::Bool(false)),
                (
                    "cond_formatting_rules",
                    object(vec![
                        ("pass", equals("pass")),
                        (
                            "warn",
                            Json::Array(vec![
                                object(vec![("s_eq", Json::text("warn"))]),
                                object(vec![("s_eq", Json::text("incomplete"))]),
                            ]),
                        ),
                        ("fail", equals("fail")),
                    ]),
                ),
            ]),
        ),
        (
            "checks_failed".to_owned(),
            count_header("Failed", "Checks that failed"),
        ),
        (
            "checks_warned".to_owned(),
            count_header("Warned", "Checks that warned"),
        ),
        (
            "checks_not_measured".to_owned(),
            count_header(
                "Not measured",
                "Checks whose file, row or column is missing",
            ),
        ),
    ];
    let mut data = Vec::new();
    for verdict in verdicts {
        let mut row = vec![
            ("verdict".to_owned(), Json::text(&verdict.status)),
            (
                "checks_failed".to_owned(),
                Json::int(verdict.count(Outcome::Fail) as i64),
            ),
            (
                "checks_warned".to_owned(),
                Json::int(verdict.count(Outcome::Warn) as i64),
            ),
            (
                "checks_not_measured".to_owned(),
                Json::int(verdict.count(Outcome::NotMeasured) as i64),
            ),
        ];
        for check in &verdict.checks {
            if !headers.iter().any(|(key, _)| *key == check.rule) {
                let mut description = format!("{}: {}", check.label, check.metric);
                for (name, bound) in [("warn", &check.warn), ("fail", &check.fail)] {
                    if let Some(bound) = bound {
                        description.push_str(&format!(", {name} {bound}"));
                    }
                }
                headers.push((
                    check.rule.clone(),
                    object(vec![
                        ("title", Json::text(&check.rule)),
                        ("description", Json::text(&description)),
                        ("format", Json::text("{:,.6g}")),
                        ("scale", Json::Bool(false)),
                    ]),
                ));
            }
            let value = match (check.outcome, check.value) {
                (Outcome::NotMeasured, _) => continue,
                (_, Some(value)) if value.is_finite() => Json::float(value),
                _ => Json::text(&check.shown),
            };
            row.push((check.rule.clone(), value));
        }
        data.push((verdict.sample.clone(), Json::Object(row)));
    }
    let presets: Vec<&str> =
        verdicts
            .iter()
            .map(|v| v.preset.as_str())
            .fold(Vec::new(), |mut seen, preset| {
                if !seen.contains(&preset) {
                    seen.push(preset);
                }
                seen
            });
    object(vec![
        ("id", Json::text("bammetrics_verdict")),
        ("section_name", Json::text("bammetrics verdict")),
        (
            "description",
            Json::text(&format!(
                "The bammetrics verdict per sample and the value each rule judged. Preset {}.",
                presets.join(", ")
            )),
        ),
        ("plot_type", Json::text("table")),
        (
            "pconfig",
            object(vec![
                ("id", Json::text("bammetrics_verdict_table")),
                ("title", Json::text("bammetrics: verdict")),
                ("namespace", Json::text("bammetrics")),
            ]),
        ),
        ("headers", Json::Object(headers)),
        ("data", Json::Object(data)),
    ])
    .render()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verdict::tests::verdict;

    #[test]
    fn one_row_per_sample_and_one_column_per_rule() {
        let mut second = verdict();
        second.sample = "S2".to_owned();
        second.status = "incomplete".to_owned();
        let text = write(&[verdict(), second]);
        let json = Json::parse(&text).expect("valid json");
        assert_eq!(json.get("plot_type"), Some(&Json::text("table")));
        assert_eq!(json.get("id"), Some(&Json::text("bammetrics_verdict")));
        let headers = json
            .get("headers")
            .and_then(Json::as_object)
            .expect("headers");
        let keys: Vec<&str> = headers.iter().map(|(k, _)| k.as_str()).collect();
        assert_eq!(
            keys,
            [
                "verdict",
                "checks_failed",
                "checks_warned",
                "checks_not_measured",
                "dup.rate",
                "q30",
                "insert.median",
                "wgs.depth"
            ]
        );
        let data = json.get("data").and_then(Json::as_object).expect("data");
        assert_eq!(data.len(), 2);
        let (sample, row) = &data[0];
        assert_eq!(sample, "NA|12878");
        assert_eq!(row.get("verdict"), Some(&Json::text("fail")));
        assert_eq!(
            row.get("checks_failed"),
            Some(&Json::Number("1".to_owned()))
        );
        assert_eq!(row.get("dup.rate"), Some(&Json::Number("0.31".to_owned())));
        // A rule that was not measured has no cell.
        assert_eq!(row.get("wgs.depth"), None);
        assert!(text.contains("Preset generic v3 (0123ab)."));
    }
}

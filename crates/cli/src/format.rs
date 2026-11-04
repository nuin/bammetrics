//! `--format` of the tool verbs (docs/05 section 2.2): the Picard metrics
//! file, or the same table as JSON, JSON lines, TSV or aligned text.
//!
//! The table is the metrics class (one row per metrics row), or the
//! histogram for a tool whose file has no metrics class
//! (`QualityScoreDistribution`, `MeanQualityByCycle`). Cells are the typed
//! values: a JSON number for an int or a float (`null` for NaN and
//! infinity), the Picard text in TSV and aligned text.

use bammetrics_core::metric::{MetricsOutput, Value};
use bammetrics_verdict::json::Json;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Picard,
    Json,
    Tsv,
    Jsonl,
    Human,
}

/// Names `--format` takes, in `--help` order.
pub const NAMES: [&str; 5] = ["picard", "json", "tsv", "jsonl", "human"];

impl Format {
    pub fn parse(name: &str) -> Option<Self> {
        Some(match name {
            "picard" => Format::Picard,
            "json" => Format::Json,
            "tsv" => Format::Tsv,
            "jsonl" => Format::Jsonl,
            "human" => Format::Human,
            _ => return None,
        })
    }
}

fn json_value(value: &Value) -> Json {
    match value {
        Value::Int(v) => Json::Number(*v as f64),
        Value::Float(v) => Json::number(Some(*v)),
        Value::Str(s) => Json::text(s),
        Value::Bool(b) => Json::Bool(*b),
        Value::Null => Json::Null,
    }
}

/// The table `--format` prints: column names and rows.
struct Table {
    columns: Vec<String>,
    rows: Vec<Vec<Value>>,
}

fn table(output: &MetricsOutput) -> Table {
    if output.is_histogram_only() {
        let records: Vec<Vec<(String, Value)>> =
            output.histogram.iter().flat_map(|h| h.records()).collect();
        let columns = records
            .first()
            .map(|r| r.iter().map(|(c, _)| c.clone()).collect())
            .unwrap_or_default();
        let rows = records
            .into_iter()
            .map(|r| r.into_iter().map(|(_, v)| v).collect())
            .collect();
        return Table { columns, rows };
    }
    Table {
        columns: output.columns.iter().map(|c| (*c).to_owned()).collect(),
        rows: output.rows.clone(),
    }
}

fn record(columns: &[String], row: &[Value]) -> Json {
    Json::Object(
        columns
            .iter()
            .zip(row)
            .map(|(c, v)| (c.clone(), json_value(v)))
            .collect(),
    )
}

/// The histogram section as {"key_class", "columns", "rows"}.
fn histogram(output: &MetricsOutput) -> Json {
    let Some(section) = &output.histogram else {
        return Json::Null;
    };
    let records = section.records();
    let columns: Vec<String> = records
        .first()
        .map(|r| r.iter().map(|(c, _)| c.clone()).collect())
        .unwrap_or_default();
    Json::Object(vec![
        (
            "columns".to_owned(),
            Json::Array(columns.iter().map(|c| Json::text(c)).collect()),
        ),
        (
            "rows".to_owned(),
            Json::Array(
                records
                    .iter()
                    .map(|r| {
                        Json::Object(r.iter().map(|(c, v)| (c.clone(), json_value(v))).collect())
                    })
                    .collect(),
            ),
        ),
    ])
}

/// What a verb prints for one tool's output: the tool, the `StringHeader`
/// command line and the sample (`--sample`, or its default).
pub struct Context<'a> {
    pub tool: &'a str,
    pub command: &'a str,
    pub sample: Option<&'a str>,
}

/// `output` as `format`. The Picard form is the file the registry writer
/// writes, `StringHeader` lines included; `sample` is added to the JSON
/// forms when there is one.
pub fn render(format: Format, context: &Context<'_>, output: &MetricsOutput) -> String {
    let table = table(output);
    let (tool, sample) = (context.tool, context.sample);
    match format {
        Format::Picard => {
            let mut bytes = Vec::new();
            output
                .write(output.headers(context.command), &mut bytes)
                .expect("writing to memory");
            String::from_utf8_lossy(&bytes).into_owned()
        }
        Format::Json => {
            let mut members = vec![("tool".to_owned(), Json::text(tool))];
            if let Some(sample) = sample {
                members.push(("sample".to_owned(), Json::text(sample)));
            }
            let metrics = if output.is_histogram_only() {
                Vec::new()
            } else {
                table
                    .rows
                    .iter()
                    .map(|r| record(&table.columns, r))
                    .collect()
            };
            members.extend([
                (
                    "class".to_owned(),
                    Json::option_text(Some(output.class).filter(|c| !c.is_empty())),
                ),
                ("metrics".to_owned(), Json::Array(metrics)),
                ("histogram".to_owned(), histogram(output)),
            ]);
            Json::Object(members).render()
        }
        Format::Jsonl => {
            let mut text = String::new();
            for row in &table.rows {
                let mut line = record(&table.columns, row);
                if let (Some(sample), Json::Object(members)) = (sample, &mut line) {
                    members.insert(0, ("sample".to_owned(), Json::text(sample)));
                }
                text.push_str(&line.render_line());
                text.push('\n');
            }
            text
        }
        Format::Tsv => {
            let mut text = table.columns.join("\t");
            text.push('\n');
            for row in &table.rows {
                let cells: Vec<String> = row.iter().map(Value::picard_text).collect();
                text.push_str(&cells.join("\t"));
                text.push('\n');
            }
            text
        }
        Format::Human => human(&table),
    }
}

/// Aligned text: a table when it has at most 6 columns, else each row as
/// one `COLUMN  value` line per column, rows apart by a blank line.
fn human(table: &Table) -> String {
    let cells: Vec<Vec<String>> = table
        .rows
        .iter()
        .map(|r| r.iter().map(Value::picard_text).collect())
        .collect();
    let mut text = String::new();
    if table.columns.len() <= 6 {
        let widths: Vec<usize> = (0..table.columns.len())
            .map(|i| {
                cells
                    .iter()
                    .map(|r| r[i].len())
                    .chain([table.columns[i].len()])
                    .max()
                    .unwrap_or(0)
            })
            .collect();
        let line = |values: &[String]| {
            let padded: Vec<String> = values
                .iter()
                .zip(&widths)
                .map(|(v, w)| format!("{v:>w$}"))
                .collect();
            format!("{}\n", padded.join("  "))
        };
        text.push_str(&line(&table.columns));
        for row in &cells {
            text.push_str(&line(row));
        }
        return text;
    }
    let width = table.columns.iter().map(String::len).max().unwrap_or(0);
    for (i, row) in cells.iter().enumerate() {
        if i > 0 {
            text.push('\n');
        }
        for (column, value) in table.columns.iter().zip(row) {
            text.push_str(&format!("{column:<width$}  {value}\n"));
        }
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use bammetrics_core::metric::MetricsOutput;

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
            Some(&Json::Number(310.0))
        );
        assert_eq!(metrics[0].get("SAMPLE"), Some(&Json::Null));
        assert_eq!(parsed.get("sample"), Some(&Json::text("S1")));
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
    }
}

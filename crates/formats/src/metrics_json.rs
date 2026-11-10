//! The JSON form of a metrics file, schema `bammetrics.metrics/1`
//! (`docs/schema/metrics-1.json`), written and read.
//!
//! ```json
//! {
//!   "schema": "bammetrics.metrics/1",
//!   "tool": "CollectInsertSizeMetrics",
//!   "sample": "NA12878",
//!   "headers": ["CollectInsertSizeMetrics INPUT=...", "Started on: ..."],
//!   "class": "picard.analysis.InsertSizeMetrics",
//!   "columns": ["MEDIAN_INSERT_SIZE", "READ_PAIRS"],
//!   "types": ["float", "int"],
//!   "metrics": [{"MEDIAN_INSERT_SIZE": 310, "READ_PAIRS": 5}],
//!   "non_finite": [],
//!   "histogram": {"key_class": "java.lang.Integer", "columns": [...],
//!                 "types": [...], "rows": [{...}], "non_finite": []}
//! }
//! ```
//!
//! A cell is a JSON number (an int exactly, a float with the digits that
//! read back as the same `f64`), a string, or `null`. `null` stands for an
//! empty cell and for a float JSON cannot hold; `non_finite` lists the
//! second kind (`{"row": 0, "column": "GC_DROPOUT", "value": "NaN"}`, the
//! value `NaN`, `Infinity` or `-Infinity`), so the Picard text comes back
//! with `?` where it had `?`. `tool` and `sample` are present when known.
//! A file with more than one metrics class, which Picard never writes,
//! carries the further ones under `more_classes`.
//!
//! `columns` holds the members of a row object, so no two are equal. A
//! table that repeats a column name (the histogram of `CollectHsMetrics`
//! with two units) gets [`Table::keys`] there (`v`, `v.1`) and the names
//! of the file under `picard_columns`, which is absent otherwise. A JSON
//! parser that keeps one value per member then loses no series.

use bammetrics_core::metric::Value;

use crate::error::{Error, Result};
use crate::json::Json;
use crate::table::{Document, Table, Ty};

pub const SCHEMA: &str = "bammetrics.metrics/1";

/// A cell as JSON; `null` for an empty cell and for NaN and infinity.
pub fn cell(value: &Value) -> Json {
    match value {
        Value::Int(v) => Json::int(*v),
        Value::Float(v) => Json::float(*v),
        Value::Str(_) | Value::Bool(_) => Json::Str(value.picard_text()),
        Value::Null => Json::Null,
    }
}

/// One row as an object, a member per column in column order under the
/// table's `keys` ([`Table::keys`]).
pub fn record(keys: &[String], row: &[Value]) -> Vec<(String, Json)> {
    keys.iter()
        .zip(row)
        .map(|(key, value)| (key.clone(), cell(value)))
        .collect()
}

fn non_finite_name(value: f64) -> &'static str {
    if value.is_nan() {
        "NaN"
    } else if value > 0.0 {
        "Infinity"
    } else {
        "-Infinity"
    }
}

/// The members every table has: `columns` (the row keys), `picard_columns`
/// when a name repeats, `types`, the rows under `rows_key` and
/// `non_finite`, which names a column by its key.
fn table_members(table: Option<&Table>, rows_key: &str) -> Vec<(String, Json)> {
    let empty = Table {
        class: String::new(),
        columns: Vec::new(),
        rows: Vec::new(),
    };
    let table = table.unwrap_or(&empty);
    let keys = table.keys();
    let mut non_finite = Vec::new();
    for (i, row) in table.rows.iter().enumerate() {
        for (key, value) in keys.iter().zip(row) {
            if let Value::Float(v) = value
                && !v.is_finite()
            {
                non_finite.push(Json::Object(vec![
                    ("row".to_owned(), Json::int(i as i64)),
                    ("column".to_owned(), Json::text(key)),
                    ("value".to_owned(), Json::text(non_finite_name(*v))),
                ]));
            }
        }
    }
    let mut members = vec![(
        "columns".to_owned(),
        Json::strings(keys.iter().map(String::as_str)),
    )];
    if table.names().ne(keys.iter().map(String::as_str)) {
        members.push(("picard_columns".to_owned(), Json::strings(table.names())));
    }
    members.extend([
        (
            "types".to_owned(),
            Json::strings(table.columns.iter().map(|(_, ty)| ty.name())),
        ),
        (
            rows_key.to_owned(),
            Json::Array(
                table
                    .rows
                    .iter()
                    .map(|row| Json::Object(record(&keys, row)))
                    .collect(),
            ),
        ),
        ("non_finite".to_owned(), Json::Array(non_finite)),
    ]);
    members
}

fn class_members(table: Option<&Table>) -> Vec<(String, Json)> {
    let class = table.map(|t| t.class.as_str()).filter(|c| !c.is_empty());
    let mut members = vec![("class".to_owned(), Json::option_text(class))];
    members.extend(table_members(table, "metrics"));
    members
}

/// `doc` as a JSON value.
pub fn value(doc: &Document) -> Json {
    let mut members = vec![("schema".to_owned(), Json::text(SCHEMA))];
    for (key, text) in [("tool", &doc.tool), ("sample", &doc.sample)] {
        if let Some(text) = text {
            members.push((key.to_owned(), Json::text(text)));
        }
    }
    members.push((
        "headers".to_owned(),
        Json::strings(doc.headers.iter().map(String::as_str)),
    ));
    members.extend(class_members(doc.tables.first()));
    let histogram = doc.histogram.as_ref().map_or(Json::Null, |table| {
        let mut members = vec![("key_class".to_owned(), Json::text(&table.class))];
        members.extend(table_members(Some(table), "rows"));
        Json::Object(members)
    });
    members.push(("histogram".to_owned(), histogram));
    if doc.tables.len() > 1 {
        let more = doc.tables[1..]
            .iter()
            .map(|t| Json::Object(class_members(Some(t))))
            .collect();
        members.push(("more_classes".to_owned(), Json::Array(more)));
    }
    Json::Object(members)
}

/// `doc` as pretty-printed JSON with a final newline.
pub fn write(doc: &Document) -> String {
    value(doc).render()
}

fn bad(message: impl Into<String>) -> Error {
    Error::Document(message.into())
}

fn strings(json: &Json, key: &str) -> Result<Option<Vec<String>>> {
    let Some(value) = json.get(key) else {
        return Ok(None);
    };
    value
        .as_array()
        .and_then(|items| {
            items
                .iter()
                .map(|i| i.as_str().map(str::to_owned))
                .collect::<Option<Vec<String>>>()
        })
        .map(Some)
        .ok_or_else(|| bad(format!("{key} is not an array of strings")))
}

/// A JSON cell as the text a metrics file would hold, for a table that
/// names no types.
fn cell_text(json: &Json) -> Result<String> {
    Ok(match json {
        Json::Null => String::new(),
        Json::Number(text) => text.clone(),
        Json::Str(text) => text.clone(),
        Json::Bool(b) => if *b { "Y" } else { "N" }.to_owned(),
        Json::Array(_) | Json::Object(_) => return Err(bad("a cell is an array or an object")),
    })
}

/// A JSON cell as a value of a column of type `ty`.
fn typed_cell(json: &Json, ty: Ty) -> Option<Value> {
    Some(match (json, ty) {
        (Json::Null, _) => Value::Null,
        (Json::Number(text), Ty::Int) => Value::Int(text.parse().ok()?),
        (Json::Number(text), Ty::Float) => Value::Float(text.parse().ok()?),
        (Json::Number(text) | Json::Str(text), Ty::Str) => Value::Str(text.clone()),
        (Json::Bool(b), Ty::Str) => Value::Str(if *b { "Y" } else { "N" }.to_owned()),
        _ => return None,
    })
}

/// The cells of one row object in column order: by position when its
/// members are the columns in order (so a repeated column name keeps each
/// cell), else by name, a missing member as `null`.
fn row_cells<'a>(row: &'a Json, columns: &[String], what: &str, i: usize) -> Result<Vec<&'a Json>> {
    let members = row
        .as_object()
        .ok_or_else(|| bad(format!("{what} row {i} is not an object")))?;
    let in_order =
        members.len() == columns.len() && members.iter().zip(columns).all(|((k, _), c)| k == c);
    if in_order {
        return Ok(members.iter().map(|(_, v)| v).collect());
    }
    if let Some((key, _)) = members.iter().find(|(k, _)| !columns.contains(k)) {
        return Err(bad(format!("{what} row {i} has {key:?}, not a column")));
    }
    Ok(columns
        .iter()
        .map(|c| row.get(c).unwrap_or(&Json::Null))
        .collect())
}

/// Put back the floats `non_finite` lists: each names a row and, by its
/// key in `keys`, a column whose cell is `null`.
fn restore_non_finite(table: &mut Table, keys: &[String], json: &Json, what: &str) -> Result<()> {
    let Some(entries) = json.get("non_finite") else {
        return Ok(());
    };
    let entries = entries
        .as_array()
        .ok_or_else(|| bad(format!("{what} non_finite is not an array")))?;
    for entry in entries {
        let field = |key: &str| {
            entry
                .get(key)
                .ok_or_else(|| bad(format!("{what} non_finite entry has no {key}")))
        };
        let row: usize = match field("row")? {
            Json::Number(text) => text.parse().ok(),
            _ => None,
        }
        .filter(|row| *row < table.rows.len())
        .ok_or_else(|| bad(format!("{what} non_finite names a row the table lacks")))?;
        let column = field("column")?.as_str().unwrap_or_default();
        let value = match field("value")?.as_str() {
            Some("NaN") => f64::NAN,
            Some("Infinity") => f64::INFINITY,
            Some("-Infinity") => f64::NEG_INFINITY,
            _ => {
                return Err(bad(format!(
                    "{what} non_finite value is not NaN, Infinity or -Infinity"
                )));
            }
        };
        // The first float column of that key whose cell is still null (a
        // document written by hand may repeat a key).
        let index = keys
            .iter()
            .zip(&table.columns)
            .enumerate()
            .position(|(j, (key, (_, ty)))| {
                key == column && *ty == Ty::Float && table.rows[row][j].is_null()
            })
            .ok_or_else(|| {
                bad(format!(
                    "{what} non_finite names {column:?} in row {row}, which is no empty float cell"
                ))
            })?;
        table.rows[row][index] = Value::Float(value);
    }
    Ok(())
}

/// The table of `json`: `columns` (the row keys), optional
/// `picard_columns` (the names, when they are not the keys), optional
/// `types`, the rows under `rows_key` and `non_finite`. `untyped` builds
/// the table from text cells when `types` is absent.
fn read_table(
    json: &Json,
    class: &str,
    rows_key: &str,
    what: &str,
    untyped: impl FnOnce(&[String], &[Vec<String>]) -> Table,
) -> Result<Table> {
    let rows = match json.get(rows_key) {
        Some(Json::Array(rows)) => rows.as_slice(),
        Some(_) => return Err(bad(format!("{what} {rows_key} is not an array"))),
        None => &[],
    };
    let columns = match strings(json, "columns")? {
        Some(columns) => columns,
        // Without a column list the first row names them.
        None => rows
            .first()
            .and_then(Json::as_object)
            .map(|members| members.iter().map(|(k, _)| k.clone()).collect())
            .unwrap_or_default(),
    };
    let cells = rows
        .iter()
        .enumerate()
        .map(|(i, row)| row_cells(row, &columns, what, i))
        .collect::<Result<Vec<_>>>()?;
    let keys = columns;
    let columns = match strings(json, "picard_columns")? {
        Some(names) if names.len() == keys.len() => names,
        Some(names) => {
            return Err(bad(format!(
                "{what} has {} picard_columns for {} columns",
                names.len(),
                keys.len()
            )));
        }
        None => keys.clone(),
    };
    let mut table = match strings(json, "types")? {
        None => {
            let text = cells
                .iter()
                .map(|row| row.iter().map(|c| cell_text(c)).collect())
                .collect::<Result<Vec<Vec<String>>>>()?;
            untyped(&columns, &text)
        }
        Some(names) => {
            if names.len() != columns.len() {
                return Err(bad(format!(
                    "{what} has {} types for {} columns",
                    names.len(),
                    columns.len()
                )));
            }
            let types = names
                .iter()
                .map(|n| Ty::parse(n).ok_or_else(|| bad(format!("{what} has no type {n:?}"))))
                .collect::<Result<Vec<Ty>>>()?;
            let rows = cells
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    row.iter()
                        .zip(&types)
                        .zip(&columns)
                        .map(|((cell, ty), column)| {
                            typed_cell(cell, *ty).ok_or_else(|| {
                                bad(format!(
                                    "{what} row {i}: {column} is not {} {}",
                                    if *ty == Ty::Int { "an" } else { "a" },
                                    ty.name()
                                ))
                            })
                        })
                        .collect()
                })
                .collect::<Result<Vec<Vec<Value>>>>()?;
            Table {
                class: class.to_owned(),
                columns: columns.into_iter().zip(types).collect(),
                rows,
            }
        }
    };
    restore_non_finite(&mut table, &keys, json, what)?;
    Ok(table)
}

fn read_class(json: &Json, what: &str) -> Result<Table> {
    let class = match json.get("class") {
        Some(Json::Str(class)) => class.as_str(),
        Some(Json::Null) | None => "",
        Some(_) => return Err(bad(format!("{what} class is not a string"))),
    };
    read_table(json, class, "metrics", what, |columns, rows| {
        Table::class_from_text(class, columns, rows)
    })
}

/// Parse a `bammetrics.metrics/1` document.
pub fn read(text: &str) -> Result<Document> {
    let json = Json::parse(text)?;
    if json.as_object().is_none() {
        return Err(bad("not a metrics document: the JSON is not an object"));
    }
    match json.get("schema").and_then(Json::as_str) {
        Some(SCHEMA) => {}
        Some(other) => return Err(bad(format!("schema {other:?} is not {SCHEMA}"))),
        None => return Err(bad(format!("not a metrics document: no schema {SCHEMA}"))),
    }
    let text_of = |key: &str| match json.get(key) {
        Some(Json::Str(text)) => Ok(Some(text.clone())),
        Some(Json::Null) | None => Ok(None),
        Some(_) => Err(bad(format!("{key} is not a string"))),
    };
    let tool = text_of("tool")?;
    let mut tables = Vec::new();
    let first = read_class(&json, "metrics")?;
    if !first.class.is_empty() || !first.columns.is_empty() {
        tables.push(first);
    }
    if let Some(more) = json.get("more_classes") {
        let more = more
            .as_array()
            .ok_or_else(|| bad("more_classes is not an array"))?;
        for class in more {
            tables.push(read_class(class, "more_classes")?);
        }
    }
    let histogram = match json.get("histogram") {
        Some(Json::Null) | None => None,
        Some(histogram) => {
            let key_class = histogram
                .get("key_class")
                .and_then(Json::as_str)
                .ok_or_else(|| bad("histogram has no key_class"))?;
            Some(read_table(
                histogram,
                key_class,
                "rows",
                "histogram",
                |columns, rows| {
                    Table::histogram_from_text(tool.as_deref(), key_class, columns, rows)
                },
            )?)
        }
    };
    Ok(Document {
        sample: text_of("sample")?,
        headers: strings(&json, "headers")?.unwrap_or_default(),
        tool,
        tables,
        histogram,
        verbatim: None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::picard;

    fn document() -> Document {
        Document {
            tool: Some("CollectGcBiasMetrics".to_owned()),
            sample: Some("S1".to_owned()),
            headers: vec!["CollectGcBiasMetrics INPUT=x.bam".to_owned()],
            tables: vec![Table {
                class: "picard.analysis.GcBiasSummaryMetrics".to_owned(),
                columns: vec![
                    ("ACCUMULATION_LEVEL".to_owned(), Ty::Str),
                    ("TOTAL_CLUSTERS".to_owned(), Ty::Int),
                    ("AT_DROPOUT".to_owned(), Ty::Float),
                    ("GC_DROPOUT".to_owned(), Ty::Float),
                ],
                rows: vec![
                    vec![
                        Value::Str("All Reads".to_owned()),
                        Value::Int(i64::MAX),
                        Value::Float(0.1),
                        Value::Float(f64::NAN),
                    ],
                    vec![
                        Value::Null,
                        Value::Null,
                        Value::Float(f64::NEG_INFINITY),
                        Value::Null,
                    ],
                ],
            }],
            histogram: Some(Table {
                class: "java.lang.Integer".to_owned(),
                columns: vec![("k".to_owned(), Ty::Int), ("v".to_owned(), Ty::Float)],
                rows: vec![vec![Value::Int(1), Value::Float(2.5)]],
            }),
            verbatim: None,
        }
    }

    #[test]
    fn every_kind_of_cell_reads_back() {
        let doc = document();
        let text = write(&doc);
        assert!(text.starts_with("{\n  \"schema\": \"bammetrics.metrics/1\",\n  \"tool\": "));
        assert!(text.contains("\"TOTAL_CLUSTERS\": 9223372036854775807"));
        assert!(text.contains("\"GC_DROPOUT\": null"));
        assert!(text.contains("\"value\": \"NaN\""));
        assert!(text.contains("\"value\": \"-Infinity\""));
        let back = read(&text).expect("reads");
        // NaN is not equal to itself: compare the Picard text.
        assert_eq!(
            picard::write(&back).expect("writes"),
            picard::write(&doc).expect("writes")
        );
        assert_eq!(back.sample, doc.sample);
        assert_eq!(back.tables[0].columns, doc.tables[0].columns);
        assert_eq!(back.tables[0].rows[1], doc.tables[0].rows[1]);
        assert_eq!(back.histogram, doc.histogram);
    }

    #[test]
    fn a_document_without_types_is_typed_by_the_glossary() {
        let text = r#"{"schema": "bammetrics.metrics/1", "class": "picard.analysis.InsertSizeMetrics",
            "metrics": [{"MEDIAN_INSERT_SIZE": 310, "PAIR_ORIENTATION": "FR", "SAMPLE": null}]}"#;
        let doc = read(text).expect("reads");
        assert_eq!(
            doc.tables[0].columns,
            [
                ("MEDIAN_INSERT_SIZE".to_owned(), Ty::Float),
                ("PAIR_ORIENTATION".to_owned(), Ty::Str),
                ("SAMPLE".to_owned(), Ty::Str),
            ]
        );
        assert_eq!(doc.tables[0].rows[0][0], Value::Float(310.0));
        assert!(doc.histogram.is_none());
    }

    #[test]
    fn rows_read_by_name_and_repeated_columns_by_position() {
        let by_name = r#"{"schema": "bammetrics.metrics/1", "class": "X", "columns": ["A", "B"],
            "types": ["int", "string"], "metrics": [{"B": "x"}]}"#;
        let doc = read(by_name).expect("reads");
        assert_eq!(
            doc.tables[0].rows,
            [[Value::Null, Value::Str("x".to_owned())]]
        );
        let repeated = r#"{"schema": "bammetrics.metrics/1", "class": null, "columns": [],
            "types": [], "metrics": [], "histogram": {"key_class": "java.lang.Integer",
            "columns": ["k", "v", "v"], "types": ["int", "float", "float"],
            "rows": [{"k": 1, "v": 2, "v": null}],
            "non_finite": [{"row": 0, "column": "v", "value": "NaN"}]}}"#;
        let doc = read(repeated).expect("reads");
        assert!(doc.tables.is_empty());
        let row = &doc.histogram.expect("histogram").rows[0];
        assert_eq!(row[1], Value::Float(2.0));
        assert!(matches!(row[2], Value::Float(v) if v.is_nan()));
    }

    /// The CollectHsMetrics shape with two units: every series is a member
    /// of its own, so a parser that keeps one value per member loses none.
    #[test]
    fn a_repeated_column_is_a_member_of_its_own() {
        let doc = Document {
            histogram: Some(Table {
                class: "java.lang.Integer".to_owned(),
                columns: ["k", "v", "k", "v"]
                    .iter()
                    .map(|n| ((*n).to_owned(), Ty::Float))
                    .collect(),
                rows: vec![vec![
                    Value::Float(0.0),
                    Value::Float(1.0),
                    Value::Float(0.0),
                    Value::Float(f64::NAN),
                ]],
            }),
            ..Document::default()
        };
        let text = write(&doc);
        let json = Json::parse(&text).expect("valid JSON");
        let histogram = json.get("histogram").expect("histogram");
        assert_eq!(
            histogram.get("columns"),
            Some(&Json::strings(["k", "v", "k.1", "v.1"]))
        );
        assert_eq!(
            histogram.get("picard_columns"),
            Some(&Json::strings(["k", "v", "k", "v"]))
        );
        let rows = histogram
            .get("rows")
            .and_then(Json::as_array)
            .expect("rows");
        let members: Vec<&str> = rows[0]
            .as_object()
            .expect("a row object")
            .iter()
            .map(|(k, _)| k.as_str())
            .collect();
        assert_eq!(members, ["k", "v", "k.1", "v.1"]);
        assert!(text.contains("\"column\": \"v.1\""));
        let back = read(&text).expect("reads");
        assert_eq!(
            picard::write(&back).expect("writes"),
            picard::write(&doc).expect("writes")
        );
        assert_eq!(
            back.histogram.expect("histogram").columns,
            doc.histogram.expect("histogram").columns
        );
        // No table of the other fixture repeats a name: no further member.
        assert!(!write(&document()).contains("picard_columns"));
    }

    #[test]
    fn more_than_one_class_reads_back() {
        let mut doc = document();
        doc.tables.push(Table {
            class: "picard.analysis.GcBiasDetailMetrics".to_owned(),
            columns: vec![("GC".to_owned(), Ty::Int)],
            rows: vec![vec![Value::Int(0)]],
        });
        let back = read(&write(&doc)).expect("reads");
        assert_eq!(back.tables.len(), 2);
        assert_eq!(back.tables[1], doc.tables[1]);
    }

    #[test]
    fn documents_that_do_not_fit_the_schema_are_refused() {
        let cases = [
            ("[]", "not a metrics document"),
            ("{}", "no schema bammetrics.metrics/1"),
            (
                r#"{"schema": "bammetrics.verdict/1"}"#,
                "is not bammetrics.metrics/1",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": []}"#,
                "0 types for 1 columns",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": ["int"], "metrics": [{"A": 1.5}]}"#,
                "row 0: A is not an int",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": ["float"], "metrics": [{"A": "x"}]}"#,
                "row 0: A is not a float",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": ["bool"]}"#,
                "no type \"bool\"",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": ["int"], "metrics": [{"B": 1}]}"#,
                "has \"B\", not a column",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A"], "types": ["float"], "metrics": [{"A": 1}],
                    "non_finite": [{"row": 0, "column": "A", "value": "NaN"}]}"#,
                "no empty float cell",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "histogram": {"columns": []}}"#,
                "histogram has no key_class",
            ),
            (
                r#"{"schema": "bammetrics.metrics/1", "columns": ["A", "A.1"], "picard_columns": ["A"]}"#,
                "1 picard_columns for 2 columns",
            ),
            ("{\"schema\": ", "json: byte"),
        ];
        for (text, expected) in cases {
            let error = read(text).expect_err(text).to_string();
            assert!(error.contains(expected), "{text}: {error}");
        }
    }
}

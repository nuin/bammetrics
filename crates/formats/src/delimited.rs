//! TSV, CSV, JSON Lines and aligned text: one table each.
//!
//! * TSV: the Picard cell text (`?` for NaN, at most 6 fraction digits),
//!   with tab, newline, carriage return and backslash in a column name or
//!   a string cell written as `\t`, `\n`, `\r` and `\\`, so every row is
//!   one line with one cell per column.
//! * CSV: RFC 4180, written by hand. `\r\n` line ends; a field with a
//!   comma, a quote, `\r` or `\n` is quoted, a quote doubled. A float
//!   keeps every digit; NaN is `NaN`, infinity `inf` and `-inf`; an empty
//!   cell is an empty field.
//! * JSON Lines: one object per row, a metrics row as its columns, a
//!   histogram row with `"histogram": true` first. A column name the table
//!   repeats is a member of its own (`v`, `v.1`: `Table::keys`).
//! * Aligned text: a table when it has at most 6 columns, else each row as
//!   one `COLUMN  value` line per column.

use bammetrics_core::metric::Value;

use crate::json::Json;
use crate::metrics_json::record;
use crate::table::{Document, Table};

fn tsv_escape(text: &str, out: &mut String) {
    for c in text.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '\t' => out.push_str("\\t"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c => out.push(c),
        }
    }
}

/// `table` as tab-separated text with a header line.
pub fn tsv(table: &Table) -> String {
    let mut text = String::new();
    for (j, name) in table.names().enumerate() {
        if j > 0 {
            text.push('\t');
        }
        tsv_escape(name, &mut text);
    }
    text.push('\n');
    for row in &table.rows {
        for (j, value) in row.iter().enumerate() {
            if j > 0 {
                text.push('\t');
            }
            match value {
                Value::Str(s) => tsv_escape(s, &mut text),
                other => text.push_str(&other.picard_text()),
            }
        }
        text.push('\n');
    }
    text
}

fn csv_field(text: &str, out: &mut String) {
    if text.contains([',', '"', '\r', '\n']) {
        out.push('"');
        out.push_str(&text.replace('"', "\"\""));
        out.push('"');
    } else {
        out.push_str(text);
    }
}

fn csv_cell(value: &Value) -> String {
    match value {
        Value::Float(v) if v.is_nan() => "NaN".to_owned(),
        Value::Float(v) if v.is_infinite() => if *v > 0.0 { "inf" } else { "-inf" }.to_owned(),
        Value::Float(v) => v.to_string(),
        other => other.picard_text(),
    }
}

/// `table` as RFC 4180 CSV with a header record.
pub fn csv(table: &Table) -> String {
    let mut text = String::new();
    for (j, name) in table.names().enumerate() {
        if j > 0 {
            text.push(',');
        }
        csv_field(name, &mut text);
    }
    text.push_str("\r\n");
    for row in &table.rows {
        for (j, value) in row.iter().enumerate() {
            if j > 0 {
                text.push(',');
            }
            csv_field(&csv_cell(value), &mut text);
        }
        text.push_str("\r\n");
    }
    text
}

/// One line per row of every class table, then one per histogram row with
/// `"histogram": true` first; `sample`, when the document has one, leads
/// the columns.
pub fn jsonl(doc: &Document) -> String {
    let mut text = String::new();
    let tables = doc.tables.iter().map(|t| (t, false));
    for (table, histogram) in tables.chain(doc.histogram.iter().map(|t| (t, true))) {
        let keys = table.keys();
        for row in &table.rows {
            let mut members = Vec::with_capacity(row.len() + 2);
            if histogram {
                members.push(("histogram".to_owned(), Json::Bool(true)));
            }
            if let Some(sample) = &doc.sample {
                members.push(("sample".to_owned(), Json::text(sample)));
            }
            members.extend(record(&keys, row));
            text.push_str(&Json::Object(members).render_line());
            text.push('\n');
        }
    }
    text
}

/// Aligned text: a table when it has at most 6 columns, else each row as
/// one `COLUMN  value` line per column, rows apart by a blank line.
pub fn human(table: &Table) -> String {
    let names: Vec<&str> = table.names().collect();
    let cells: Vec<Vec<String>> = table
        .rows
        .iter()
        .map(|r| r.iter().map(Value::picard_text).collect())
        .collect();
    let mut text = String::new();
    if names.len() <= 6 {
        let widths: Vec<usize> = (0..names.len())
            .map(|i| {
                cells
                    .iter()
                    .map(|r| r[i].len())
                    .chain([names[i].len()])
                    .max()
                    .unwrap_or(0)
            })
            .collect();
        let mut line = |values: &mut dyn Iterator<Item = &str>| {
            let padded: Vec<String> = values
                .zip(&widths)
                .map(|(v, w)| format!("{v:>w$}"))
                .collect();
            text.push_str(&padded.join("  "));
            text.push('\n');
        };
        line(&mut names.iter().copied());
        for row in &cells {
            line(&mut row.iter().map(String::as_str));
        }
        return text;
    }
    let width = names.iter().map(|n| n.len()).max().unwrap_or(0);
    for (i, row) in cells.iter().enumerate() {
        if i > 0 {
            text.push('\n');
        }
        for (column, value) in names.iter().zip(row) {
            text.push_str(&format!("{column:<width$}  {value}\n"));
        }
    }
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::table::Ty;

    fn table() -> Table {
        Table {
            class: "X".to_owned(),
            columns: vec![
                ("NAME".to_owned(), Ty::Str),
                ("N".to_owned(), Ty::Int),
                ("X".to_owned(), Ty::Float),
            ],
            rows: vec![
                vec![
                    Value::Str("a\tb\\c\nd\re, \"q\"".to_owned()),
                    Value::Int(-3),
                    Value::Float(1.0 / 3.0),
                ],
                vec![Value::Null, Value::Null, Value::Float(f64::NAN)],
                vec![
                    Value::Str("?".to_owned()),
                    Value::Int(0),
                    Value::Float(f64::NEG_INFINITY),
                ],
            ],
        }
    }

    #[test]
    fn tsv_escapes_string_cells_and_keeps_picard_numbers() {
        assert_eq!(
            tsv(&table()),
            "NAME\tN\tX\n\
             a\\tb\\\\c\\nd\\re, \"q\"\t-3\t0.333333\n\
             \t\t?\n\
             ?\t0\t-?\n"
        );
        let plain = Table {
            class: String::new(),
            columns: vec![("A\tB".to_owned(), Ty::Str)],
            rows: vec![vec![Value::Str("FR".to_owned())]],
        };
        assert_eq!(tsv(&plain), "A\\tB\nFR\n");
    }

    #[test]
    fn csv_follows_rfc_4180() {
        assert_eq!(
            csv(&table()),
            "NAME,N,X\r\n\
             \"a\tb\\c\nd\re, \"\"q\"\"\",-3,0.3333333333333333\r\n\
             ,,NaN\r\n\
             ?,0,-inf\r\n"
        );
        let header = Table {
            class: String::new(),
            columns: vec![("a,b".to_owned(), Ty::Float)],
            rows: vec![vec![Value::Float(f64::INFINITY)], vec![Value::Float(310.0)]],
        };
        assert_eq!(csv(&header), "\"a,b\"\r\ninf\r\n310\r\n");
    }

    #[test]
    fn jsonl_marks_histogram_rows() {
        let doc = Document {
            sample: Some("S1".to_owned()),
            tables: vec![table()],
            histogram: Some(Table {
                class: "java.lang.Integer".to_owned(),
                columns: vec![("k".to_owned(), Ty::Int), ("v".to_owned(), Ty::Float)],
                rows: vec![vec![Value::Int(1), Value::Float(2.0)]],
            }),
            ..Document::default()
        };
        let text = jsonl(&doc);
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines.len(), 4);
        assert_eq!(
            lines[1],
            "{\"sample\":\"S1\",\"NAME\":null,\"N\":null,\"X\":null}"
        );
        assert_eq!(
            lines[3],
            "{\"histogram\":true,\"sample\":\"S1\",\"k\":1,\"v\":2}"
        );
        for line in lines {
            assert!(Json::parse(line).is_ok(), "{line}");
        }
        let repeated = Document {
            histogram: Some(Table {
                class: "java.lang.Integer".to_owned(),
                columns: ["k", "v", "k", "v"]
                    .iter()
                    .map(|n| ((*n).to_owned(), Ty::Int))
                    .collect(),
                rows: vec![(0..4).map(Value::Int).collect()],
            }),
            ..Document::default()
        };
        assert_eq!(
            jsonl(&repeated),
            "{\"histogram\":true,\"k\":0,\"v\":1,\"k.1\":2,\"v.1\":3}\n"
        );
    }

    #[test]
    fn human_aligns_narrow_tables_and_lists_wide_ones() {
        let narrow = Table {
            class: String::new(),
            columns: vec![("MEDIAN".to_owned(), Ty::Float), ("S".to_owned(), Ty::Str)],
            rows: vec![vec![Value::Float(310.0), Value::Null]],
        };
        assert_eq!(human(&narrow), "MEDIAN  S\n   310   \n");
        let wide = Table {
            class: String::new(),
            columns: (0..7).map(|i| (format!("C{i}"), Ty::Int)).collect(),
            rows: vec![(0..7).map(Value::Int).collect(); 2],
        };
        let text = human(&wide);
        assert!(text.starts_with("C0  0\nC1  1\n"));
        assert_eq!(text.matches("\n\n").count(), 1);
    }
}

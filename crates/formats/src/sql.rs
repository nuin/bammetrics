//! The SQL dump of a metrics file: one text file that `sqlite3 DB <
//! FILE.sql` loads (docs/formats.md, "SQL dump").
//!
//! ```sql
//! -- bammetrics SQL dump of one metrics file: sqlite3 DB < FILE.sql
//! BEGIN TRANSACTION;
//! CREATE TABLE "headers" ("position" BIGINT, "text" TEXT);
//! INSERT INTO "headers" VALUES (1, 'CollectInsertSizeMetrics INPUT=in.bam');
//! CREATE TABLE "tables" ("name" TEXT, "kind" TEXT, "class" TEXT, "tool" TEXT, "sample" TEXT);
//! CREATE TABLE "columns" ("table" TEXT, "position" BIGINT, "column" TEXT, "name" TEXT, "type" TEXT);
//! CREATE TABLE "InsertSizeMetrics" ("MEDIAN_INSERT_SIZE" DOUBLE, "READ_PAIRS" BIGINT);
//! INSERT INTO "InsertSizeMetrics" VALUES (310.0, 5);
//! COMMIT;
//! ```
//!
//! * One table per metrics class, named by the class after its last `.`
//!   or `$` (`InsertSizeMetrics`), and the histogram as `histogram`. The
//!   names `headers`, `tables`, `columns` and `histogram` come first; a
//!   class whose name is taken, in any case, gets `_2`, `_3`. A class
//!   name that starts with `sqlite_`, in any case, gets `_` before it:
//!   SQLite refuses to create such a table.
//! * A column is `BIGINT`, `DOUBLE` or `TEXT` by the glossary type of
//!   [`Table`]. SQLite gives them the affinities INTEGER, REAL and TEXT;
//!   DuckDB reads them as 64-bit integers and doubles (its `INTEGER` is 32
//!   bits and its `REAL` 32-bit float).
//! * A column name the table repeats, or repeats in another case, gets
//!   `.1`, `.2` as [`Table::keys`] gives them. `columns` keeps the name as
//!   the file has it, with its type (`int`, `float`, `string`).
//! * Identifiers are in double quotes, a `"` doubled; strings in single
//!   quotes, a `'` doubled. A carriage return or a NUL in a string is
//!   `char(13)` or `char(0)`: the sqlite3 shell drops a carriage return
//!   before a line end. Such an identifier gets U+240D or U+2400 in its
//!   place, and `columns` keeps the name.
//! * A float keeps every digit (the shortest text that reads back as the
//!   same `f64`). NaN and infinity are `NULL`, as an empty cell is, and
//!   SQLite and DuckDB read `-0.0` back as `0.0`: the losses of the form.
//! * One `INSERT` per row, all in one transaction. The text depends on the
//!   document only.

use std::collections::{HashMap, HashSet};
use std::fmt::Write as _;

use bammetrics_core::metric::Value;

use crate::table::{Document, Table, Ty};

/// The first line of every dump; [`is_sql`] looks for it.
pub const FIRST_LINE: &str = "-- bammetrics SQL dump of one metrics file: sqlite3 DB < FILE.sql";

/// The tables every dump has, in order.
const CATALOGUE: [&str; 3] = ["headers", "tables", "columns"];

/// The name the histogram table takes.
const HISTOGRAM: &str = "histogram";

/// The cell of a row shorter than its table.
static NULL: Value = Value::Null;

/// True when `text` is a dump [`write`] wrote.
pub fn is_sql(text: &str) -> bool {
    text.starts_with(FIRST_LINE)
}

/// The SQL type of a column of type `ty`.
pub fn sql_type(ty: Ty) -> &'static str {
    match ty {
        Ty::Int => "BIGINT",
        Ty::Float => "DOUBLE",
        Ty::Str => "TEXT",
    }
}

/// Names no two equal when compared without ASCII case, as SQLite and
/// DuckDB compare identifiers. A name whose folded form is taken gets
/// `{sep}{n}`, `n` from `first` up, skipping a taken one. `last` keeps
/// the highest `n` given per name, so the work is linear in the names.
struct Names {
    taken: HashSet<String>,
    last: HashMap<String, u64>,
    sep: char,
    first: u64,
}

impl Names {
    fn new(sep: char, first: u64) -> Self {
        Self {
            taken: HashSet::new(),
            last: HashMap::new(),
            sep,
            first,
        }
    }

    fn take(&mut self, name: &str) -> String {
        let folded = name.to_ascii_lowercase();
        if self.taken.insert(folded.clone()) {
            return name.to_owned();
        }
        let n = self.last.entry(folded).or_insert(self.first - 1);
        loop {
            *n += 1;
            let candidate = format!("{name}{}{n}", self.sep);
            if self.taken.insert(candidate.to_ascii_lowercase()) {
                return candidate;
            }
        }
    }
}

/// `name` as an identifier can hold it: a carriage return as U+240D, a
/// NUL as U+2400, an empty name as `_`.
fn identifier_text(name: &str) -> String {
    if name.is_empty() {
        return "_".to_owned();
    }
    name.replace('\r', "\u{240d}").replace('\0', "\u{2400}")
}

/// A class name as a table identifier: [`identifier_text`], with `_` before
/// a name that starts with `sqlite_` in any case, which SQLite keeps for
/// its own tables.
fn table_identifier(name: &str) -> String {
    let text = identifier_text(name);
    let reserved = text
        .as_bytes()
        .get(..7)
        .is_some_and(|b| b.eq_ignore_ascii_case(b"sqlite_"));
    if reserved { format!("_{text}") } else { text }
}

fn push_identifier(name: &str, out: &mut String) {
    out.push('"');
    out.push_str(&name.replace('"', "\"\""));
    out.push('"');
}

/// `text` as an SQL string: in single quotes, a quote doubled, a carriage
/// return or a NUL as `char(13)` or `char(0)` joined with `||`.
fn push_string(text: &str, out: &mut String) {
    let mut open = false;
    let mut first = true;
    for c in text.chars() {
        match c {
            '\r' | '\0' => {
                if open {
                    out.push('\'');
                    open = false;
                }
                if !first {
                    out.push_str("||");
                }
                let _ = write!(out, "char({})", u32::from(c));
            }
            c => {
                if !open {
                    if !first {
                        out.push_str("||");
                    }
                    out.push('\'');
                    open = true;
                }
                if c == '\'' {
                    out.push('\'');
                }
                out.push(c);
            }
        }
        first = false;
    }
    if first {
        out.push_str("''");
    } else if open {
        out.push('\'');
    }
}

/// A cell as an SQL literal: `NULL` for an empty cell and for NaN and
/// infinity.
fn push_value(value: &Value, out: &mut String) {
    match value {
        Value::Null => out.push_str("NULL"),
        Value::Int(v) => {
            let _ = write!(out, "{v}");
        }
        Value::Float(v) if v.is_finite() => {
            // `Debug` gives the shortest digits that read back as `v`, with
            // a point or an exponent, so a reader types it as a double.
            let _ = write!(out, "{v:?}");
        }
        Value::Float(_) => out.push_str("NULL"),
        Value::Str(_) | Value::Bool(_) => push_string(&value.picard_text(), out),
    }
}

fn insert<'v>(table: &str, cells: impl IntoIterator<Item = &'v Value>, out: &mut String) {
    out.push_str("INSERT INTO ");
    push_identifier(table, out);
    out.push_str(" VALUES (");
    for (j, cell) in cells.into_iter().enumerate() {
        if j > 0 {
            out.push_str(", ");
        }
        push_value(cell, out);
    }
    out.push_str(");\n");
}

fn create(table: &str, columns: &[(&str, &str)], out: &mut String) {
    out.push_str("CREATE TABLE ");
    push_identifier(table, out);
    out.push_str(" (");
    for (j, (name, ty)) in columns.iter().enumerate() {
        if j > 0 {
            out.push_str(", ");
        }
        push_identifier(name, out);
        out.push(' ');
        out.push_str(ty);
    }
    out.push_str(");\n");
}

/// One table of the document with its SQL name and column names.
struct Placed<'a> {
    name: String,
    kind: &'static str,
    table: &'a Table,
    columns: Vec<String>,
}

/// The tables of `doc` that get an SQL table (a table with no column gets
/// none), each with its name and its column names: the classes, then the
/// histogram.
fn place(doc: &Document) -> Vec<Placed<'_>> {
    let mut names = Names::new('_', 2);
    for reserved in CATALOGUE {
        names.take(reserved);
    }
    let histogram = doc.histogram.as_ref().filter(|t| !t.columns.is_empty());
    let histogram_name = histogram.map(|_| names.take(HISTOGRAM));
    doc.tables
        .iter()
        .filter(|t| !t.columns.is_empty())
        .map(|table| {
            let short = table.class.rsplit(['.', '$']).next().unwrap_or_default();
            let short = if short.is_empty() { "metrics" } else { short };
            (names.take(&table_identifier(short)), "metrics", table)
        })
        .collect::<Vec<_>>()
        .into_iter()
        .chain(
            histogram
                .zip(histogram_name)
                .map(|(t, n)| (n, "histogram", t)),
        )
        .map(|(name, kind, table)| {
            let mut keys = Names::new('.', 1);
            let columns = table
                .names()
                .map(|n| keys.take(&identifier_text(n)))
                .collect();
            Placed {
                name,
                kind,
                table,
                columns,
            }
        })
        .collect()
}

/// The SQL names of the tables of `doc`: each class, then the histogram.
/// A table with no column has none.
pub fn table_names(doc: &Document) -> Vec<String> {
    place(doc).into_iter().map(|p| p.name).collect()
}

fn text(text: &str) -> Value {
    Value::Str(text.to_owned())
}

fn text_or_null(value: Option<&str>) -> Value {
    value.map_or(Value::Null, text)
}

/// `doc` as an SQL dump.
pub fn write(doc: &Document) -> String {
    let placed = place(doc);
    let mut out = String::new();
    out.push_str(FIRST_LINE);
    out.push_str("\nBEGIN TRANSACTION;\n");
    create(
        "headers",
        &[("position", "BIGINT"), ("text", "TEXT")],
        &mut out,
    );
    for (i, header) in (1..).zip(&doc.headers) {
        insert("headers", &[Value::Int(i), text(header)], &mut out);
    }
    let catalogue = [
        ("name", "TEXT"),
        ("kind", "TEXT"),
        ("class", "TEXT"),
        ("tool", "TEXT"),
        ("sample", "TEXT"),
    ];
    create("tables", &catalogue, &mut out);
    for p in &placed {
        let row = [
            text(&p.name),
            text(p.kind),
            text(&p.table.class),
            text_or_null(doc.tool.as_deref()),
            text_or_null(doc.sample.as_deref()),
        ];
        insert("tables", &row, &mut out);
    }
    let catalogue = [
        ("table", "TEXT"),
        ("position", "BIGINT"),
        ("column", "TEXT"),
        ("name", "TEXT"),
        ("type", "TEXT"),
    ];
    create("columns", &catalogue, &mut out);
    for p in &placed {
        for (i, (key, (name, ty))) in (1..).zip(p.columns.iter().zip(&p.table.columns)) {
            let row = [
                text(&p.name),
                Value::Int(i),
                text(key),
                text(name),
                text(ty.name()),
            ];
            insert("columns", &row, &mut out);
        }
    }
    for p in &placed {
        let columns: Vec<(&str, &str)> = p
            .columns
            .iter()
            .zip(&p.table.columns)
            .map(|(key, (_, ty))| (key.as_str(), sql_type(*ty)))
            .collect();
        create(&p.name, &columns, &mut out);
        for row in &p.table.rows {
            let cells = (0..p.columns.len()).map(|j| row.get(j).unwrap_or(&NULL));
            insert(&p.name, cells, &mut out);
        }
    }
    out.push_str("COMMIT;\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn string(text: &str) -> String {
        let mut out = String::new();
        push_string(text, &mut out);
        out
    }

    fn table(class: &str, columns: &[(&str, Ty)], rows: Vec<Vec<Value>>) -> Table {
        Table {
            class: class.to_owned(),
            columns: columns.iter().map(|(n, t)| ((*n).to_owned(), *t)).collect(),
            rows,
        }
    }

    #[test]
    fn strings_double_a_quote_and_spell_a_carriage_return() {
        assert_eq!(string(""), "''");
        assert_eq!(string("it's"), "'it''s'");
        assert_eq!(string("a\tb\nc"), "'a\tb\nc'");
        assert_eq!(string("a\r\nb"), "'a'||char(13)||'\nb'");
        assert_eq!(string("\r"), "char(13)");
        assert_eq!(string("\0x\r"), "char(0)||'x'||char(13)");
        let mut out = String::new();
        push_identifier("a\"b", &mut out);
        assert_eq!(out, "\"a\"\"b\"");
    }

    #[test]
    fn cells_keep_every_digit_and_null_the_rest() {
        let literal = |value: Value| {
            let mut out = String::new();
            push_value(&value, &mut out);
            out
        };
        assert_eq!(literal(Value::Float(310.0)), "310.0");
        assert_eq!(literal(Value::Float(1.0 / 3.0)), "0.3333333333333333");
        assert_eq!(literal(Value::Float(1e-7)), "1e-7");
        assert_eq!(literal(Value::Float(f64::NAN)), "NULL");
        assert_eq!(literal(Value::Float(f64::NEG_INFINITY)), "NULL");
        assert_eq!(literal(Value::Int(i64::MIN)), "-9223372036854775808");
        assert_eq!(literal(Value::Bool(true)), "'Y'");
        assert_eq!(literal(Value::Null), "NULL");
    }

    #[test]
    fn names_never_collide_in_any_case() {
        let doc = Document {
            tables: vec![
                table("a.b.Histogram", &[("A", Ty::Int)], vec![]),
                table("x.Headers", &[("A", Ty::Int)], vec![]),
                table("y.Headers", &[("A", Ty::Int)], vec![]),
                table("z.Headers_2", &[("A", Ty::Int)], vec![]),
                table("Outer$Inner", &[("A", Ty::Int)], vec![]),
                table("", &[("A", Ty::Int)], vec![]),
                table("no.Columns", &[], vec![]),
                table("x.sqlite_master", &[("A", Ty::Int)], vec![]),
                table("x.SQLite_X", &[("A", Ty::Int)], vec![]),
                table("x._sqlite_x", &[("A", Ty::Int)], vec![]),
                table("x.sqlite", &[("A", Ty::Int)], vec![]),
            ],
            histogram: Some(table("java.lang.Integer", &[("k", Ty::Int)], vec![])),
            ..Document::default()
        };
        assert_eq!(
            table_names(&doc),
            [
                "Histogram_2",
                "Headers_2",
                "Headers_3",
                "Headers_2_2",
                "Inner",
                "metrics",
                "_sqlite_master",
                "_SQLite_X",
                "_sqlite_x_2",
                "sqlite",
                "histogram"
            ]
        );
        let mut keys = Names::new('.', 1);
        let names = ["v", "V", "v.1", "v", "", "", "a\rb"];
        let taken: Vec<String> = names
            .iter()
            .map(|n| keys.take(&identifier_text(n)))
            .collect();
        assert_eq!(
            taken,
            ["v", "V.1", "v.1.1", "v.2", "_", "_.1", "a\u{240d}b"]
        );
    }

    #[test]
    fn a_document_is_one_transaction() {
        let doc = Document {
            tool: Some("CollectX".to_owned()),
            headers: vec![
                "CollectX INPUT=in.bam".to_owned(),
                "Started on: x".to_owned(),
            ],
            tables: vec![table(
                "picard.analysis.XMetrics",
                &[("N", Ty::Int), ("X", Ty::Float), ("S", Ty::Str)],
                vec![
                    vec![
                        Value::Int(5),
                        Value::Float(f64::NAN),
                        Value::Str("FR".to_owned()),
                    ],
                    vec![Value::Null, Value::Float(0.5)],
                ],
            )],
            histogram: Some(table(
                "java.lang.Integer",
                &[("k", Ty::Int), ("v", Ty::Float), ("v", Ty::Float)],
                vec![vec![Value::Int(0), Value::Float(1.0), Value::Float(2.0)]],
            )),
            ..Document::default()
        };
        let sql = write(&doc);
        assert!(is_sql(&sql));
        let expected = [
            FIRST_LINE,
            "BEGIN TRANSACTION;",
            r#"CREATE TABLE "headers" ("position" BIGINT, "text" TEXT);"#,
            r#"INSERT INTO "headers" VALUES (1, 'CollectX INPUT=in.bam');"#,
            r#"INSERT INTO "headers" VALUES (2, 'Started on: x');"#,
            r#"CREATE TABLE "tables" ("name" TEXT, "kind" TEXT, "class" TEXT, "tool" TEXT, "sample" TEXT);"#,
            r#"INSERT INTO "tables" VALUES ('XMetrics', 'metrics', 'picard.analysis.XMetrics', 'CollectX', NULL);"#,
            r#"INSERT INTO "tables" VALUES ('histogram', 'histogram', 'java.lang.Integer', 'CollectX', NULL);"#,
            r#"CREATE TABLE "columns" ("table" TEXT, "position" BIGINT, "column" TEXT, "name" TEXT, "type" TEXT);"#,
            r#"INSERT INTO "columns" VALUES ('XMetrics', 1, 'N', 'N', 'int');"#,
            r#"INSERT INTO "columns" VALUES ('XMetrics', 2, 'X', 'X', 'float');"#,
            r#"INSERT INTO "columns" VALUES ('XMetrics', 3, 'S', 'S', 'string');"#,
            r#"INSERT INTO "columns" VALUES ('histogram', 1, 'k', 'k', 'int');"#,
            r#"INSERT INTO "columns" VALUES ('histogram', 2, 'v', 'v', 'float');"#,
            r#"INSERT INTO "columns" VALUES ('histogram', 3, 'v.1', 'v', 'float');"#,
            r#"CREATE TABLE "XMetrics" ("N" BIGINT, "X" DOUBLE, "S" TEXT);"#,
            r#"INSERT INTO "XMetrics" VALUES (5, NULL, 'FR');"#,
            r#"INSERT INTO "XMetrics" VALUES (NULL, 0.5, NULL);"#,
            r#"CREATE TABLE "histogram" ("k" BIGINT, "v" DOUBLE, "v.1" DOUBLE);"#,
            r#"INSERT INTO "histogram" VALUES (0, 1.0, 2.0);"#,
            "COMMIT;",
        ];
        assert_eq!(sql.lines().collect::<Vec<_>>(), expected);
        let empty = write(&Document::default());
        assert!(empty.ends_with("\"type\" TEXT);\nCOMMIT;\n"), "{empty}");
    }
}

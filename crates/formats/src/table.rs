//! The typed table every format is written from (docs/07 section 3.1).
//!
//! A [`Table`] is one metrics class: its columns, each with a [`Ty`], and
//! its rows. A column's type comes from the glossary
//! ([`glossary::lookup`]); a column the glossary has no row for is typed
//! from its cells. A cell that does not fit its column's type widens the
//! column (int to float, either to string), so every cell of a column has
//! the column's type or is null.
//!
//! A [`Document`] is one metrics file: the `StringHeader` texts, the class
//! tables (Picard writes one) and the histogram as a second table.

use bammetrics_core::glossary;
use bammetrics_core::metric::{MetricsOutput, Value};
use bammetrics_core::metrics_file::parse_value;

/// The type of a column.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Ty {
    Int,
    Float,
    Str,
}

impl Ty {
    /// `int`, `float` or `string`, the glossary's words.
    pub fn name(self) -> &'static str {
        match self {
            Ty::Int => "int",
            Ty::Float => "float",
            Ty::Str => "string",
        }
    }

    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "int" => Some(Ty::Int),
            "float" => Some(Ty::Float),
            "string" => Some(Ty::Str),
            _ => None,
        }
    }

    /// The type of a non-null cell.
    fn of(value: &Value) -> Option<Self> {
        match value {
            Value::Int(_) => Some(Ty::Int),
            Value::Float(_) => Some(Ty::Float),
            Value::Str(_) | Value::Bool(_) => Some(Ty::Str),
            Value::Null => None,
        }
    }
}

/// The glossary's type of `column` in `class`.
pub fn glossary_type(class: &str, column: &str) -> Option<Ty> {
    glossary::lookup(class, column).and_then(|e| Ty::parse(e.kind))
}

/// The types of a histogram's columns: the key by its Java class
/// (`java.lang.Integer` and `Byte` are ints, `Double` a float), else by the
/// glossary; each series by the glossary row of the program that wrote it
/// (found by the key label when `tool` is not known), else float, the type
/// htsjdk reads a histogram value as.
fn histogram_types(tool: Option<&str>, key_class: &str, labels: &[String]) -> Vec<Option<Ty>> {
    let Some(key) = labels.first() else {
        return Vec::new();
    };
    let program = tool.unwrap_or("");
    let lookup = |label: &str| {
        glossary::histogram_lookup(program, key, label).and_then(|e| Ty::parse(e.kind))
    };
    let key_ty = match key_class.rsplit('.').next() {
        Some("Integer" | "Byte" | "Short" | "Long") => Some(Ty::Int),
        Some("Double" | "Float") => Some(Ty::Float),
        _ => lookup(key),
    };
    std::iter::once(key_ty)
        .chain(
            labels[1..]
                .iter()
                .map(|label| lookup(label).or(Some(Ty::Float))),
        )
        .collect()
}

/// A cell of a text file as a number when it is one: [`parse_value`], plus
/// a run of digits too long for an `i64` as a float (a double above 2^63
/// prints with no point).
fn parse_cell(text: &str) -> Value {
    match parse_value(text) {
        Value::Str(s) => {
            let digits = s.strip_prefix('-').unwrap_or(&s);
            if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
                s.parse().map_or(Value::Str(s), Value::Float)
            } else {
                Value::Str(s)
            }
        }
        other => other,
    }
}

/// The type of a column whose cells have the types `seen`: the declared
/// type widened by every cell that does not fit it, or, with none
/// declared, the narrowest type every cell fits (string for no cell).
fn column_type(declared: Option<Ty>, seen: impl Iterator<Item = Ty>) -> Ty {
    let widest = seen.max();
    match (declared, widest) {
        (Some(declared), Some(widest)) => declared.max(widest),
        (Some(declared), None) => declared,
        (None, Some(widest)) => widest,
        (None, None) => Ty::Str,
    }
}

/// `value` as a cell of a column of type `ty`; `text` gives the cell's
/// own text for a string column.
fn coerce(value: Value, ty: Ty, text: impl FnOnce(&Value) -> String) -> Value {
    match (ty, value) {
        (_, Value::Null) => Value::Null,
        (Ty::Float, Value::Int(v)) => Value::Float(v as f64),
        (Ty::Str, Value::Str(s)) => Value::Str(s),
        (Ty::Str, other) => Value::Str(text(&other)),
        (_, other) => other,
    }
}

/// One metrics class, or the histogram of a file.
#[derive(Clone, Debug, PartialEq)]
pub struct Table {
    /// The Picard metrics class; for a histogram the Java class of its
    /// key (`java.lang.Integer`); empty for a file with neither.
    pub class: String,
    pub columns: Vec<(String, Ty)>,
    /// One cell per column per row: null, or a value of the column's type.
    pub rows: Vec<Vec<Value>>,
}

impl Table {
    /// A table from typed rows (what a tool returns): `declared` gives a
    /// column's type, a boolean becomes `Y` or `N` as Picard prints it.
    pub fn from_values(
        class: &str,
        columns: &[&str],
        rows: &[Vec<Value>],
        declared: impl Fn(usize, &str) -> Option<Ty>,
    ) -> Self {
        let types: Vec<Ty> = columns
            .iter()
            .enumerate()
            .map(|(j, name)| {
                column_type(
                    declared(j, name),
                    rows.iter().filter_map(|r| r.get(j).and_then(Ty::of)),
                )
            })
            .collect();
        Self {
            class: class.to_owned(),
            columns: columns
                .iter()
                .zip(&types)
                .map(|(name, ty)| ((*name).to_owned(), *ty))
                .collect(),
            rows: rows
                .iter()
                .map(|row| {
                    types
                        .iter()
                        .enumerate()
                        .map(|(j, ty)| {
                            let value = row.get(j).cloned().unwrap_or(Value::Null);
                            coerce(value, *ty, Value::picard_text)
                        })
                        .collect()
                })
                .collect(),
        }
    }

    /// A table from text cells (what a metrics file holds). An empty cell
    /// is null; a cell of a string column keeps its text.
    pub fn from_text(
        class: &str,
        columns: &[String],
        rows: &[Vec<String>],
        declared: impl Fn(usize, &str) -> Option<Ty>,
    ) -> Self {
        fn cell(row: &[String], j: usize) -> &str {
            row.get(j).map_or("", String::as_str)
        }
        let types: Vec<Ty> = columns
            .iter()
            .enumerate()
            .map(|(j, name)| match declared(j, name) {
                Some(Ty::Str) => Ty::Str,
                declared => column_type(
                    declared,
                    rows.iter().filter_map(|r| Ty::of(&parse_cell(cell(r, j)))),
                ),
            })
            .collect();
        Self {
            class: class.to_owned(),
            columns: columns.iter().cloned().zip(types.iter().copied()).collect(),
            rows: rows
                .iter()
                .map(|row| {
                    types
                        .iter()
                        .enumerate()
                        .map(|(j, ty)| {
                            let text = cell(row, j);
                            match ty {
                                Ty::Str if text.is_empty() => Value::Null,
                                Ty::Str => Value::Str(text.to_owned()),
                                _ => coerce(parse_cell(text), *ty, |_| text.to_owned()),
                            }
                        })
                        .collect()
                })
                .collect(),
        }
    }

    /// A class table from text cells, typed by the glossary.
    pub fn class_from_text(class: &str, columns: &[String], rows: &[Vec<String>]) -> Self {
        Self::from_text(class, columns, rows, |_, name| glossary_type(class, name))
    }

    /// A histogram table from text cells; see [`histogram_types`].
    pub fn histogram_from_text(
        tool: Option<&str>,
        key_class: &str,
        labels: &[String],
        rows: &[Vec<String>],
    ) -> Self {
        let types = histogram_types(tool, key_class, labels);
        Self::from_text(key_class, labels, rows, |j, _| types[j])
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.columns.iter().map(|(name, _)| name.as_str())
    }

    /// One key per column, no two equal, for a format that holds a row as
    /// an object (JSON, JSON Lines, a Python dict): the column's name, and
    /// for a name the table repeats the name plus `.1`, `.2` and so on,
    /// skipping a key that is the name of a column. Picard repeats the
    /// histogram columns of `CollectHsMetrics` once per unit.
    pub fn keys(&self) -> Vec<String> {
        let names: Vec<&str> = self.names().collect();
        let mut keys: Vec<String> = Vec::with_capacity(names.len());
        for (j, name) in names.iter().enumerate() {
            let key = if names[..j].contains(name) {
                (1..)
                    .map(|n| format!("{name}.{n}"))
                    .find(|k| !names.contains(&k.as_str()) && !keys.contains(k))
                    .expect("an unbounded range")
            } else {
                (*name).to_owned()
            };
            keys.push(key);
        }
        keys
    }

    /// True for a histogram whose key prints as a Java double (`1.0`).
    pub fn has_double_key(&self) -> bool {
        self.class == "java.lang.Double"
    }
}

/// One metrics file.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Document {
    /// The Picard tool that wrote the file, when known.
    pub tool: Option<String>,
    /// The sample the caller named; never part of a Picard file.
    pub sample: Option<String>,
    /// The `StringHeader` texts: the command line, then `Started on: ...`.
    pub headers: Vec<String>,
    /// The `## METRICS CLASS` blocks. Picard writes one.
    pub tables: Vec<Table>,
    /// The `## HISTOGRAM` section; its `class` is the key's Java class.
    pub histogram: Option<Table>,
    /// The text of a file that is not a metrics file (an interval list):
    /// its Picard form. `tables` still describes it.
    pub verbatim: Option<String>,
}

impl Document {
    /// What `tool` returned, with `headers` as its `StringHeader` texts.
    pub fn from_output(tool: &str, headers: Vec<String>, output: &MetricsOutput) -> Self {
        let tables = if output.is_histogram_only() {
            Vec::new()
        } else {
            vec![Table::from_values(
                output.class,
                output.columns,
                &output.rows,
                |_, name| glossary_type(output.class, name),
            )]
        };
        let histogram = output.histogram.as_ref().and_then(|section| {
            let records = section.records();
            let labels: Vec<String> = records.first()?.iter().map(|(l, _)| l.clone()).collect();
            let names: Vec<&str> = labels.iter().map(String::as_str).collect();
            let rows: Vec<Vec<Value>> = records
                .into_iter()
                .map(|r| r.into_iter().map(|(_, v)| v).collect())
                .collect();
            let key_class = section.key_class.java_name();
            let types = histogram_types(Some(tool), key_class, &labels);
            Some(Table::from_values(key_class, &names, &rows, |j, _| {
                types[j]
            }))
        });
        Self {
            tool: Some(tool.to_owned()),
            sample: None,
            headers,
            tables,
            histogram,
            verbatim: output.text.clone(),
        }
    }

    /// The table a one-table format prints: the first class, or the
    /// histogram of a file without one.
    pub fn main_table(&self) -> Option<&Table> {
        self.tables.first().or(self.histogram.as_ref())
    }

    /// The table named `name`: `metrics` (the first class), `histogram`,
    /// or a class by its full name or the name after its last `.` or `$`.
    pub fn table(&self, name: &str) -> Option<&Table> {
        match name {
            "metrics" => self.tables.first(),
            "histogram" => self.histogram.as_ref(),
            class => self
                .tables
                .iter()
                .find(|t| t.class == class || t.class.rsplit(['.', '$']).next() == Some(class)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(cells: &[&str]) -> Vec<String> {
        cells.iter().map(|c| (*c).to_owned()).collect()
    }

    /// MEDIAN_INSERT_SIZE prints 310; the glossary says float.
    #[test]
    fn glossary_types_win_over_the_cell_text() {
        let table = Table::class_from_text(
            "picard.analysis.InsertSizeMetrics",
            &strings(&[
                "MEDIAN_INSERT_SIZE",
                "READ_PAIRS",
                "PAIR_ORIENTATION",
                "SAMPLE",
            ]),
            &[strings(&["310", "5", "FR", ""])],
        );
        assert_eq!(
            table.columns,
            [
                ("MEDIAN_INSERT_SIZE".to_owned(), Ty::Float),
                ("READ_PAIRS".to_owned(), Ty::Int),
                ("PAIR_ORIENTATION".to_owned(), Ty::Str),
                ("SAMPLE".to_owned(), Ty::Str),
            ]
        );
        assert_eq!(
            table.rows,
            [[
                Value::Float(310.0),
                Value::Int(5),
                Value::Str("FR".to_owned()),
                Value::Null
            ]]
        );
    }

    #[test]
    fn a_string_column_keeps_number_like_text() {
        let table = Table::class_from_text(
            "picard.analysis.InsertSizeMetrics",
            &strings(&["SAMPLE", "LIBRARY"]),
            &[strings(&["1.50", "012"]), strings(&["1e3", "?"])],
        );
        assert_eq!(table.rows[0][0], Value::Str("1.50".to_owned()));
        assert_eq!(table.rows[1][0], Value::Str("1e3".to_owned()));
        assert_eq!(table.rows[1][1], Value::Str("?".to_owned()));
    }

    #[test]
    fn a_class_the_glossary_lacks_is_inferred() {
        let table = Table::class_from_text(
            "some.other.Metrics",
            &strings(&["N", "X", "S", "EMPTY", "MIXED"]),
            &[
                strings(&["1", "2", "a", "", "1"]),
                strings(&["2", "2.5", "3", "", "x"]),
                strings(&["", "?", "", "", ""]),
            ],
        );
        let types: Vec<Ty> = table.columns.iter().map(|(_, ty)| *ty).collect();
        assert_eq!(types, [Ty::Int, Ty::Float, Ty::Str, Ty::Str, Ty::Str]);
        assert_eq!(table.rows[0][1], Value::Float(2.0));
        assert!(matches!(table.rows[2][1], Value::Float(v) if v.is_nan()));
        assert_eq!(table.rows[1][2], Value::Str("3".to_owned()));
        assert_eq!(table.rows[0][4], Value::Str("1".to_owned()));
        assert_eq!(table.rows[2][0], Value::Null);
    }

    #[test]
    fn a_cell_that_does_not_fit_widens_its_column() {
        // READ_PAIRS is an int in the glossary.
        let float = Table::class_from_text(
            "picard.analysis.InsertSizeMetrics",
            &strings(&["READ_PAIRS"]),
            &[strings(&["5"]), strings(&["2.5"])],
        );
        assert_eq!(float.columns[0].1, Ty::Float);
        assert_eq!(float.rows, [[Value::Float(5.0)], [Value::Float(2.5)]]);
        let text = Table::class_from_text(
            "picard.analysis.InsertSizeMetrics",
            &strings(&["READ_PAIRS"]),
            &[strings(&["5"]), strings(&["many"])],
        );
        assert_eq!(text.columns[0].1, Ty::Str);
        assert_eq!(text.rows[0][0], Value::Str("5".to_owned()));
        // A double above 2^63 prints with no point.
        let big = Table::class_from_text(
            "picard.analysis.InsertSizeMetrics",
            &strings(&["MEAN_INSERT_SIZE"]),
            &[strings(&["12345678901234566000000000"])],
        );
        assert_eq!(big.rows[0][0], Value::Float(1.234_567_890_123_456_6e25));
    }

    #[test]
    fn typed_rows_take_the_glossary_type_and_booleans_print_as_picard() {
        let table = Table::from_values(
            "picard.analysis.InsertSizeMetrics",
            &["MEDIAN_INSERT_SIZE", "FLAG", "SAMPLE"],
            &[vec![Value::Int(310), Value::Bool(true), Value::Null]],
            |_, name| glossary_type("picard.analysis.InsertSizeMetrics", name),
        );
        assert_eq!(
            table.rows,
            [[Value::Float(310.0), Value::Str("Y".to_owned()), Value::Null]]
        );
        assert_eq!(table.columns[1].1, Ty::Str);
    }

    #[test]
    fn a_repeated_column_gets_a_key_of_its_own() {
        let table = |names: &[&str]| Table {
            class: String::new(),
            columns: names.iter().map(|n| ((*n).to_owned(), Ty::Float)).collect(),
            rows: Vec::new(),
        };
        assert_eq!(table(&["k", "v", "w"]).keys(), ["k", "v", "w"]);
        assert_eq!(
            table(&["k", "v", "k", "v", "k", "v"]).keys(),
            ["k", "v", "k.1", "v.1", "k.2", "v.2"]
        );
        // `v.1` is a column, so the second `v` takes the next number.
        assert_eq!(table(&["v", "v", "v.1"]).keys(), ["v", "v.2", "v.1"]);
        assert_eq!(table(&["v.1", "v", "v"]).keys(), ["v.1", "v", "v.2"]);
    }

    #[test]
    fn histogram_keys_follow_their_java_class() {
        let labels = strings(&["BIN", "CoverageMult"]);
        let rows = [strings(&["1.0", "1.5"])];
        let table = Table::histogram_from_text(None, "java.lang.Double", &labels, &rows);
        assert_eq!(table.columns[0].1, Ty::Float);
        assert!(table.has_double_key());
        let labels = strings(&["insert_size", "All_Reads.fr_count", "unknown_series"]);
        let rows = [strings(&["300", "1", "2"])];
        let table = Table::histogram_from_text(None, "java.lang.Integer", &labels, &rows);
        let types: Vec<Ty> = table.columns.iter().map(|(_, ty)| *ty).collect();
        assert_eq!(types, [Ty::Int, Ty::Float, Ty::Float]);
        assert_eq!(table.rows[0][1], Value::Float(1.0));
    }
}

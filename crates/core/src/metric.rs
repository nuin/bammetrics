//! Typed metric rows and the one writer every tool shares.
//!
//! A tool describes its rows with [`MetricRow`] (column names in Picard
//! order, one [`Value`] per column) and returns a [`MetricsOutput`]. The
//! generic [`MetricsOutput::write`] prints it as a Picard metrics file and the
//! Python binding turns the same values into dicts, so a new tool writes
//! neither a formatter nor a converter.
//!
//! Printing rules, checked against Picard 3.5.0 (`MetricsFile.write`):
//! integers print as decimals, doubles through [`format_double`]
//! (at most 6 fraction digits), booleans as `Y` or `N` (htsjdk `FormatUtil`), strings
//! as they are, and null as an empty field.

use std::borrow::Cow;
use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::Arc;

use crate::histogram::Histogram;
use crate::metrics_file::{
    HistogramColumn, HistogramTable, KeyClass, MetricsFile, format_double, headers,
};

/// One field of a metrics row.
#[derive(Clone, Debug, PartialEq)]
pub enum Value {
    Int(i64),
    Float(f64),
    Str(String),
    Bool(bool),
    /// A null field: empty text in the file, `None` in Python.
    Null,
}

impl Value {
    /// The field text as Picard prints it.
    pub fn picard_text(&self) -> String {
        match self {
            Value::Int(v) => v.to_string(),
            Value::Float(v) => format_double(*v),
            Value::Str(v) => v.clone(),
            Value::Bool(v) => if *v { "Y" } else { "N" }.to_owned(),
            Value::Null => String::new(),
        }
    }

    /// The integer, for an `Int` field only.
    pub fn as_int(&self) -> Option<i64> {
        match self {
            Value::Int(v) => Some(*v),
            _ => None,
        }
    }

    /// The number: a `Float`, or an `Int` widened. NaN for a `?` field.
    pub fn as_float(&self) -> Option<f64> {
        match self {
            Value::Float(v) => Some(*v),
            Value::Int(v) => Some(*v as f64),
            _ => None,
        }
    }

    /// The text, for a `Str` field only.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Value::Str(v) => Some(v),
            _ => None,
        }
    }

    pub fn is_null(&self) -> bool {
        matches!(self, Value::Null)
    }
}

impl From<i64> for Value {
    fn from(v: i64) -> Self {
        Value::Int(v)
    }
}

impl From<f64> for Value {
    fn from(v: f64) -> Self {
        Value::Float(v)
    }
}

impl From<bool> for Value {
    fn from(v: bool) -> Self {
        Value::Bool(v)
    }
}

impl From<&str> for Value {
    fn from(v: &str) -> Self {
        Value::Str(v.to_owned())
    }
}

impl From<String> for Value {
    fn from(v: String) -> Self {
        Value::Str(v)
    }
}

impl<T: Into<Value>> From<Option<T>> for Value {
    fn from(v: Option<T>) -> Self {
        v.map_or(Value::Null, Into::into)
    }
}

/// A row of a Picard metrics class.
pub trait MetricRow {
    /// Fully qualified Picard class, e.g. `picard.analysis.InsertSizeMetrics`.
    const CLASS: &'static str;
    /// Column names in the order Picard writes them.
    const COLUMNS: &'static [&'static str];

    /// One value per column, in `COLUMNS` order.
    fn values(&self) -> Vec<Value>;
}

/// One labelled histogram column of the `## HISTOGRAM` section.
#[derive(Clone, Debug)]
pub struct HistogramSeries {
    pub label: String,
    pub histogram: Histogram,
}

/// The `## HISTOGRAM` section: a key column and one column per series.
#[derive(Clone, Debug)]
pub struct HistogramSection {
    pub key_label: &'static str,
    pub key_class: KeyClass,
    pub series: Vec<HistogramSeries>,
}

impl HistogramSection {
    /// The table as printed: one row per key in the union of the non-empty
    /// series, the key then each non-empty series' value (0 where it has no
    /// bin), labelled like the header line.
    pub fn records(&self) -> Vec<Vec<(String, Value)>> {
        let series: Vec<&HistogramSeries> = self
            .series
            .iter()
            .filter(|s| !s.histogram.is_empty())
            .collect();
        let keys: std::collections::BTreeSet<i64> =
            series.iter().flat_map(|s| s.histogram.keys()).collect();
        keys.into_iter()
            .map(|key| {
                let key_value = match self.key_class {
                    KeyClass::Double => Value::Float(key as f64),
                    KeyClass::Integer | KeyClass::Byte => Value::Int(key),
                    KeyClass::MqPair => Value::Str(self.key_class.format(key)),
                };
                std::iter::once((self.key_label.to_owned(), key_value))
                    .chain(series.iter().map(|s| {
                        (
                            s.label.clone(),
                            Value::Float(s.histogram.get(key).unwrap_or(0.0)),
                        )
                    }))
                    .collect()
            })
            .collect()
    }
}

/// Everything one tool run produces: the file for its `OUTPUT` option
/// (`None` where Picard writes none) and any further files with the path
/// each goes to (`CollectGcBiasMetrics` writes `SUMMARY_OUTPUT`,
/// `CollectMultipleMetrics` one file per program).
#[derive(Clone, Debug, Default)]
pub struct ToolOutput {
    pub main: Option<MetricsOutput>,
    pub extra: Vec<(PathBuf, MetricsOutput)>,
}

impl ToolOutput {
    /// No file at all.
    pub const NONE: Self = Self {
        main: None,
        extra: Vec::new(),
    };

    /// The `OUTPUT` file only.
    pub fn file(output: MetricsOutput) -> Self {
        Self {
            main: Some(output),
            extra: Vec::new(),
        }
    }
}

impl From<MetricsOutput> for ToolOutput {
    fn from(output: MetricsOutput) -> Self {
        Self::file(output)
    }
}

/// Rows made one at a time as a file is written, for a file whose rows
/// as [`Value`]s would not fit in memory: `CollectRrbsMetrics` writes one
/// row per CpG site and unit, about 300 bytes each as `Value`s and 3 to 4
/// packed in its own table.
pub trait RowSource: std::fmt::Debug + Send + Sync {
    /// The number of rows.
    fn len(&self) -> usize;

    /// True when there is no row.
    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Call `row` with each row in file order; the first error stops.
    fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()>;
}

/// Everything a tool returns: enough to write the file or build Python rows.
#[derive(Clone, Debug)]
pub struct MetricsOutput {
    pub class: &'static str,
    pub columns: &'static [&'static str],
    pub rows: Vec<Vec<Value>>,
    pub histogram: Option<HistogramSection>,
    /// False for a file Picard writes without `StringHeader` lines
    /// (`AccumulateQualityYieldMetrics`).
    pub with_headers: bool,
    /// A file that is not a metrics file (`BedToIntervalList` writes an
    /// interval list): written as it is, with no header lines. `rows` and
    /// `columns` still describe it for Python.
    pub text: Option<String>,
    /// Rows made as they are written, in place of `rows` (which is then
    /// empty). The Picard writer streams them; a reader of `rows` calls
    /// [`MetricsOutput::materialized`] first.
    pub source: Option<Arc<dyn RowSource>>,
}

impl MetricsOutput {
    pub fn from_rows<R: MetricRow>(rows: &[R], histogram: Option<HistogramSection>) -> Self {
        Self {
            class: R::CLASS,
            columns: R::COLUMNS,
            rows: rows.iter().map(MetricRow::values).collect(),
            histogram,
            with_headers: true,
            text: None,
            source: None,
        }
    }

    /// A file with no metrics class, only the histogram
    /// (`QualityScoreDistribution`, `MeanQualityByCycle`).
    pub fn histogram_only(histogram: HistogramSection) -> Self {
        Self {
            class: "",
            columns: &[],
            rows: Vec::new(),
            histogram: Some(histogram),
            with_headers: true,
            text: None,
            source: None,
        }
    }

    /// A plain text file described as a table for Python: `columns` name
    /// the fields of each row.
    pub fn text(columns: &'static [&'static str], rows: Vec<Vec<Value>>, text: String) -> Self {
        Self {
            class: "",
            columns,
            rows,
            histogram: None,
            with_headers: false,
            text: Some(text),
            source: None,
        }
    }

    /// The `StringHeader` texts for this file: the command line and the
    /// start time, or none.
    pub fn headers(&self, command: &str) -> Vec<String> {
        if self.with_headers {
            headers(command.to_owned())
        } else {
            Vec::new()
        }
    }

    /// True for a file without a metrics class, whose content is the
    /// histogram table.
    pub fn is_histogram_only(&self) -> bool {
        self.columns.is_empty()
    }

    /// This output with `rows` filled from its [`RowSource`], for a reader
    /// of `rows` (the JSON and delimited formats, Python); borrowed when it
    /// has no source.
    pub fn materialized(&self) -> Cow<'_, Self> {
        let Some(source) = &self.source else {
            return Cow::Borrowed(self);
        };
        let mut rows = Vec::with_capacity(source.len());
        source
            .each(&mut |row| {
                rows.push(row.to_vec());
                Ok(())
            })
            .expect("collecting rows in memory");
        Cow::Owned(Self {
            rows,
            source: None,
            ..self.clone()
        })
    }

    /// Rows as `(column, value)` pairs, the shape the Python converter wants.
    pub fn records(&self) -> impl Iterator<Item = impl Iterator<Item = (&'static str, &Value)>> {
        self.rows
            .iter()
            .map(|row| self.columns.iter().copied().zip(row.iter()))
    }

    /// Write the Picard metrics file; `headers` are the `StringHeader` texts.
    pub fn write(&self, headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
        if let Some(text) = &self.text {
            return out.write_all(text.as_bytes());
        }
        let histogram = self.histogram.as_ref().map(|h| HistogramTable {
            key_label: h.key_label,
            key_class: h.key_class,
            columns: h
                .series
                .iter()
                .map(|s| HistogramColumn {
                    label: s.label.clone(),
                    histogram: &s.histogram,
                })
                .collect(),
        });
        let file = MetricsFile {
            headers,
            class: self.class,
            columns: self.columns,
            rows: Vec::new(),
            histogram,
        };
        let Some(source) = &self.source else {
            return file.write_rows(out, &self.rows, |row| line(row));
        };
        file.write_lines(out, !source.is_empty(), |out| {
            source.each(&mut |row| writeln!(out, "{}", line(row)))
        })
    }
}

/// A row as its Picard line, without the line end.
fn line(row: &[Value]) -> String {
    row.iter()
        .map(Value::picard_text)
        .collect::<Vec<_>>()
        .join("\t")
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Demo {
        count: i64,
        mean: f64,
        flag: bool,
        name: Option<String>,
    }

    impl MetricRow for Demo {
        const CLASS: &'static str = "picard.analysis.Demo";
        const COLUMNS: &'static [&'static str] = &["COUNT", "MEAN", "FLAG", "NAME"];

        fn values(&self) -> Vec<Value> {
            vec![
                self.count.into(),
                self.mean.into(),
                self.flag.into(),
                self.name.clone().into(),
            ]
        }
    }

    #[test]
    fn rows_print_like_picard() {
        let rows = [
            Demo {
                count: 3,
                mean: 2.5,
                flag: true,
                name: None,
            },
            Demo {
                count: -1,
                mean: f64::NAN,
                flag: false,
                name: Some("x".to_owned()),
            },
        ];
        let mut out = Vec::new();
        MetricsOutput::from_rows(&rows, None)
            .write(vec!["Demo".to_owned()], &mut out)
            .expect("write to Vec");
        assert_eq!(
            String::from_utf8(out).expect("utf8"),
            "## htsjdk.samtools.metrics.StringHeader\n# Demo\n\n\
             ## METRICS CLASS\tpicard.analysis.Demo\nCOUNT\tMEAN\tFLAG\tNAME\n\
             3\t2.5\tY\t\n-1\t?\tN\tx\n\n\n"
        );
    }

    /// The rows of a [`RowSource`], held here as `Value`s.
    #[derive(Debug)]
    struct Held(Vec<Vec<Value>>);

    impl RowSource for Held {
        fn len(&self) -> usize {
            self.0.len()
        }

        fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
            self.0.iter().try_for_each(|r| row(r))
        }
    }

    #[test]
    fn a_row_source_writes_the_bytes_of_its_rows() {
        let rows = [Demo {
            count: 3,
            mean: 2.5,
            flag: true,
            name: None,
        }];
        let held = MetricsOutput::from_rows(&rows, None);
        for source in [held.rows.clone(), Vec::new()] {
            let streamed = MetricsOutput {
                rows: Vec::new(),
                source: Some(Arc::new(Held(source.clone()))),
                ..held.clone()
            };
            let expected = MetricsOutput {
                rows: source.clone(),
                ..held.clone()
            };
            let write = |output: &MetricsOutput| {
                let mut out = Vec::new();
                output
                    .write(vec!["Demo".to_owned()], &mut out)
                    .expect("write to Vec");
                out
            };
            assert_eq!(write(&streamed), write(&expected));
            assert_eq!(streamed.materialized().rows, source);
        }
    }
}

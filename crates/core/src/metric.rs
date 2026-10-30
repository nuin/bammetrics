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

use std::io::{self, Write};

use crate::histogram::Histogram;
use crate::metrics_file::{HistogramColumn, HistogramTable, MetricsFile, format_double};

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
    pub series: Vec<HistogramSeries>,
}

/// Everything a tool returns: enough to write the file or build Python rows.
#[derive(Clone, Debug)]
pub struct MetricsOutput {
    pub class: &'static str,
    pub columns: &'static [&'static str],
    pub rows: Vec<Vec<Value>>,
    pub histogram: Option<HistogramSection>,
}

impl MetricsOutput {
    pub fn from_rows<R: MetricRow>(rows: &[R], histogram: Option<HistogramSection>) -> Self {
        Self {
            class: R::CLASS,
            columns: R::COLUMNS,
            rows: rows.iter().map(MetricRow::values).collect(),
            histogram,
        }
    }

    /// Rows as `(column, value)` pairs, the shape the Python converter wants.
    pub fn records(&self) -> impl Iterator<Item = impl Iterator<Item = (&'static str, &Value)>> {
        self.rows
            .iter()
            .map(|row| self.columns.iter().copied().zip(row.iter()))
    }

    /// Write the Picard metrics file; `headers` are the `StringHeader` texts.
    pub fn write(&self, headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
        let histogram = self.histogram.as_ref().map(|h| HistogramTable {
            key_label: h.key_label,
            columns: h
                .series
                .iter()
                .map(|s| HistogramColumn {
                    label: s.label.clone(),
                    histogram: &s.histogram,
                })
                .collect(),
        });
        MetricsFile {
            headers,
            class: self.class,
            columns: self.columns,
            rows: self
                .rows
                .iter()
                .map(|r| r.iter().map(Value::picard_text).collect())
                .collect(),
            histogram,
        }
        .write(out)
    }
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
}

//! Picard metrics text, read into a [`Document`] and written from one.
//!
//! Both directions go through the core's [`TextDoc`], the text-level form
//! of the writer every tool uses, so a file Picard or bammetrics wrote
//! comes back byte for byte: headers, class, columns, every cell and the
//! histogram.
//!
//! What does not come back as it was: a header class other than
//! `StringHeader` (written as one), `\r\n` line ends, a row shorter than
//! its column line (filled), and a number in a form `FormatUtil` does not
//! print (`0,5`, `1.50`, `1e3` in a numeric column print as `0.5`, `1.5`,
//! `1000`).
//!
//! Picard text has no escape, so [`write`] refuses a document it cannot
//! hold: a tab, line feed or carriage return in a string cell, a column
//! name or a class, a line feed or carriage return in a header, a row
//! whose line would be blank (one column, an empty cell) and a row whose
//! first cell starts with `## `. Only a JSON document can carry one; a
//! file read from Picard text never does.

use bammetrics_core::java::double_to_string;
use bammetrics_core::metric::Value;
use bammetrics_core::metrics_file::{TextBlock, TextDoc};
use bammetrics_core::registry;

use crate::error::{Error, Result};
use crate::table::{Document, Table};

/// The registry tool a file's first `StringHeader` names: Picard writes
/// the command line there, the class name first.
fn tool_of(headers: &[String]) -> Option<String> {
    let name = headers.first()?.split_whitespace().next()?;
    registry::find(name).map(|tool| tool.name.to_owned())
}

/// Parse the text of a Picard metrics file.
pub fn read(text: &str) -> Result<Document> {
    let doc = read_or_empty(text)?;
    if doc.tables.is_empty() && doc.histogram.is_none() {
        return Err(Error::Metrics(
            "not a metrics file: no ## METRICS CLASS and no ## HISTOGRAM section".to_owned(),
        ));
    }
    Ok(doc)
}

/// [`read`] for a file a tool is known to have written: text with no
/// metrics class and no histogram is a document with the headers and no
/// table. `CollectDuplicateMetrics` writes such a file for a BAM with no
/// read, and its Picard text comes back byte for byte.
pub fn read_or_empty(text: &str) -> Result<Document> {
    let doc = TextDoc::parse(text).map_err(|e| Error::Metrics(e.to_string()))?;
    let tool = tool_of(&doc.headers);
    Ok(Document {
        tables: doc
            .classes
            .iter()
            .map(|b| Table::class_from_text(&b.class, &b.columns, &b.rows))
            .collect(),
        histogram: doc
            .histogram
            .as_ref()
            .map(|b| Table::histogram_from_text(tool.as_deref(), &b.class, &b.columns, &b.rows)),
        tool,
        sample: None,
        headers: doc.headers,
        verbatim: None,
    })
}

/// The cells of `table` as Picard prints them. A histogram's `Double` key
/// prints as `Double.toString` does (`1.0`), every other cell as
/// `FormatUtil` does.
fn block(table: &Table, histogram: bool) -> TextBlock {
    let double_key = histogram && table.has_double_key();
    TextBlock {
        class: table.class.clone(),
        columns: table.names().map(str::to_owned).collect(),
        rows: table
            .rows
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, value)| match value {
                        Value::Float(key) if double_key && j == 0 => double_to_string(*key),
                        other => other.picard_text(),
                    })
                    .collect()
            })
            .collect(),
    }
}

/// The first character of `text` that would end a field or a line of a
/// Picard file, by its name; a tab only when `tab` (a header keeps one).
fn breaks_a_line(text: &str, tab: bool) -> Option<&'static str> {
    text.chars().find_map(|c| match c {
        '\t' if tab => Some("tab"),
        '\n' => Some("line feed"),
        '\r' => Some("carriage return"),
        _ => None,
    })
}

/// An error for the first header, class, column name, string cell or row
/// of `doc` that Picard text cannot hold, naming where it is.
fn check_writable(doc: &Document) -> Result<()> {
    let refuse = |found: Option<&'static str>, place: &dyn Fn() -> String| match found {
        Some(character) => Err(Error::Picard(format!(
            "Picard text cannot hold the {character} in {}",
            place()
        ))),
        None => Ok(()),
    };
    for (i, header) in doc.headers.iter().enumerate() {
        refuse(breaks_a_line(header, false), &|| format!("header {i}"))?;
    }
    let classes = doc.tables.iter().map(|t| ("metrics", t));
    for (what, table) in classes.chain(doc.histogram.iter().map(|t| ("histogram", t))) {
        refuse(breaks_a_line(&table.class, true), &|| {
            format!("the {what} class {:?}", table.class)
        })?;
        for (j, name) in table.names().enumerate() {
            refuse(breaks_a_line(name, true), &|| {
                format!("the name of {what} column {j}, {name:?}")
            })?;
        }
        for (i, row) in table.rows.iter().enumerate() {
            for (name, value) in table.names().zip(row) {
                if let Value::Str(text) = value {
                    refuse(breaks_a_line(text, true), &|| {
                        format!("{what} row {i}, column {name}: {text:?}")
                    })?;
                }
            }
            // A table's rows end at a blank line or at a `## ` line.
            let first = match row.first() {
                Some(Value::Str(text)) => text.as_str(),
                Some(Value::Null) | None => "",
                // A number or a boolean prints neither empty nor as `## `.
                Some(_) => "0",
            };
            let blank = row.len() < 2 && first.is_empty();
            refuse(blank.then_some("blank line"), &|| {
                format!(
                    "{what} row {i}: the one cell of the row is empty and a blank line ends \
                     the table"
                )
            })?;
            refuse(first.starts_with("## ").then_some("section mark"), &|| {
                format!("{what} row {i}: its first cell, {first:?}, starts a new section")
            })?;
        }
    }
    Ok(())
}

/// `doc` as a Picard metrics file; the verbatim text for a file that is
/// not one. An error when a cell, a column name, the class or a header
/// holds a character that ends a field or a line (see the module text):
/// written as it is, the file would read back with other rows.
pub fn write(doc: &Document) -> Result<String> {
    if let Some(text) = &doc.verbatim {
        return Ok(text.clone());
    }
    check_writable(doc)?;
    let text = TextDoc {
        headers: doc.headers.clone(),
        classes: doc.tables.iter().map(|t| block(t, false)).collect(),
        histogram: doc.histogram.as_ref().map(|t| block(t, true)),
    };
    let mut bytes = Vec::new();
    text.write(&mut bytes).expect("writing to memory");
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::table::Ty;
    use bammetrics_core::histogram::Histogram;
    use bammetrics_core::metric::{HistogramSection, HistogramSeries, MetricsOutput};
    use bammetrics_core::metrics_file::KeyClass;

    const INSERT_SIZE: &str = "## htsjdk.samtools.metrics.StringHeader\n\
# CollectInsertSizeMetrics INPUT=x.bam\n\
## htsjdk.samtools.metrics.StringHeader\n\
# Started on: Mon Nov 02 18:48:01 MST 2026\n\
\n\
## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n\
MEDIAN_INSERT_SIZE\tSTANDARD_DEVIATION\tREAD_PAIRS\tPAIR_ORIENTATION\tSAMPLE\n\
310\t8.164966\t5\tFR\t\n\
\n\
## HISTOGRAM\tjava.lang.Integer\n\
insert_size\tAll_Reads.fr_count\n\
300\t1\n\
310\t2\n\
\n";

    #[test]
    fn a_picard_file_reads_typed_and_writes_back() {
        let doc = read(INSERT_SIZE).expect("reads");
        assert_eq!(doc.tool.as_deref(), Some("CollectInsertSizeMetrics"));
        assert_eq!(doc.headers.len(), 2);
        assert_eq!(doc.tables[0].columns[0].1, Ty::Float);
        assert_eq!(doc.tables[0].rows[0][0], Value::Float(310.0));
        let histogram = doc.histogram.as_ref().expect("histogram");
        assert_eq!(histogram.class, "java.lang.Integer");
        assert_eq!(histogram.rows[1], [Value::Int(310), Value::Float(2.0)]);
        assert_eq!(write(&doc).expect("writes"), INSERT_SIZE);
    }

    #[test]
    fn nan_null_and_the_double_key_write_back() {
        let text = "## htsjdk.samtools.metrics.StringHeader\n# CollectDuplicateMetrics I=x\n\n\
                    ## METRICS CLASS\tpicard.sam.DuplicationMetrics\n\
                    LIBRARY\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n\
                    lib\t?\t\n\t-?\t-0\n\n\
                    ## HISTOGRAM\tjava.lang.Double\nBIN\tCoverageMult\n1.0\t1\n2.0\t1.5\n\n";
        let doc = read(text).expect("reads");
        assert!(matches!(doc.tables[0].rows[0][1], Value::Float(v) if v.is_nan()));
        assert_eq!(doc.tables[0].rows[0][2], Value::Null);
        assert_eq!(doc.tables[0].rows[1][0], Value::Null);
        assert_eq!(write(&doc).expect("writes"), text);
    }

    /// A JSON document can carry what Picard text has no form for; written
    /// raw, a tab makes a further cell and a line feed a further row.
    #[test]
    fn a_cell_that_would_break_a_line_is_refused() {
        let doc = |edit: &dyn Fn(&mut Document)| {
            let mut doc = read(INSERT_SIZE).expect("reads");
            edit(&mut doc);
            write(&doc).map_err(|e| e.to_string())
        };
        assert!(doc(&|_| ()).is_ok());
        let cell = |text: &'static str| {
            doc(&move |d| d.tables[0].rows[0][3] = Value::Str(text.to_owned()))
        };
        assert_eq!(
            cell("a\tb").expect_err("tab"),
            "Picard text cannot hold the tab in metrics row 0, column PAIR_ORIENTATION: \"a\\tb\""
        );
        assert!(cell("a\nb").expect_err("lf").contains("the line feed in"));
        assert!(
            cell("a\r")
                .expect_err("cr")
                .contains("the carriage return in")
        );
        assert!(cell("a b\\c, \"d\"").is_ok());
        assert_eq!(
            doc(&|d| d.tables[0].class = "c\nd".to_owned()).expect_err("class"),
            "Picard text cannot hold the line feed in the metrics class \"c\\nd\""
        );
        assert_eq!(
            doc(&|d| d.tables[0].columns[1].0 = "A\tB".to_owned()).expect_err("column"),
            "Picard text cannot hold the tab in the name of metrics column 1, \"A\\tB\""
        );
        assert_eq!(
            doc(&|d| d.headers[1] = "a\nb".to_owned()).expect_err("header"),
            "Picard text cannot hold the line feed in header 1"
        );
        // A row's first cell that starts a section would end the table.
        let section = doc(&|d| d.tables[0].rows[0][0] = Value::Str("## METRICS CLASS".to_owned()));
        assert_eq!(
            section.expect_err("section"),
            "Picard text cannot hold the section mark in metrics row 0: its first cell, \
             \"## METRICS CLASS\", starts a new section"
        );
        // One column and an empty cell: the row would be a blank line.
        let one_column = |cell: Value| {
            write(&Document {
                tables: vec![Table {
                    class: "X".to_owned(),
                    columns: vec![("A".to_owned(), crate::table::Ty::Str)],
                    rows: vec![vec![Value::Str("a".to_owned())], vec![cell]],
                }],
                ..Document::default()
            })
            .map_err(|e| e.to_string())
        };
        assert!(one_column(Value::Str("b".to_owned())).is_ok());
        assert_eq!(
            one_column(Value::Null).expect_err("blank"),
            "Picard text cannot hold the blank line in metrics row 1: the one cell of the row \
             is empty and a blank line ends the table"
        );
        // A header keeps a tab: the text after `# ` is read whole.
        let tabbed = doc(&|d| d.headers[1] = "a\tb".to_owned()).expect("a tab in a header");
        assert_eq!(read(&tabbed).expect("reads").headers[1], "a\tb");
        let histogram = doc(&|d| {
            let table = d.histogram.as_mut().expect("histogram");
            table.columns[1].0 = "x\ry".to_owned();
        });
        assert!(
            histogram
                .expect_err("histogram column")
                .contains("the name of histogram column 1")
        );
    }

    #[test]
    fn text_that_is_no_metrics_file_is_refused() {
        for text in ["", "hello\n", "{\"schema\": 1}", "@HD\tVN:1.6\n"] {
            let error = read(text).expect_err(text);
            assert!(
                error.to_string().starts_with("not a metrics file"),
                "{error}"
            );
        }
        assert!(read("## METRICS CLASS\tX\n").is_err());
    }

    /// What CollectDuplicateMetrics writes for a BAM with no read.
    #[test]
    fn a_file_with_headers_only_reads_as_an_empty_document() {
        let text = "## htsjdk.samtools.metrics.StringHeader\n\
                    # CollectDuplicateMetrics INPUT=zero.bam\n\
                    ## htsjdk.samtools.metrics.StringHeader\n\
                    # Started on: Tue Nov 10 15:35:04 UTC 2026\n\n\n\n";
        assert!(read(text).is_err());
        let doc = read_or_empty(text).expect("reads");
        assert_eq!(doc.tool.as_deref(), Some("CollectDuplicateMetrics"));
        assert_eq!(doc.headers.len(), 2);
        assert!(doc.tables.is_empty() && doc.histogram.is_none());
        assert_eq!(write(&doc).expect("writes"), text);
        let json = crate::metrics_json::write(&doc);
        assert!(json.contains("\"class\": null"), "{json}");
        let back = crate::metrics_json::read(&json).expect("the JSON reads");
        assert_eq!(write(&back).expect("writes"), text);
        assert_eq!(crate::delimited::jsonl(&doc), "");
    }

    /// The Picard form of a tool's output is the bytes the registry writer
    /// writes for it.
    #[test]
    fn a_tool_output_writes_as_the_registry_writer_does() {
        let mut counts = Histogram::new();
        counts.increment_by(1, 1.0);
        counts.increment_by(3, 2.5);
        for key_class in [KeyClass::Integer, KeyClass::Byte, KeyClass::Double] {
            let output = MetricsOutput {
                class: "picard.analysis.InsertSizeMetrics",
                columns: &[
                    "MEDIAN_INSERT_SIZE",
                    "READ_PAIRS",
                    "PAIR_ORIENTATION",
                    "SAMPLE",
                ],
                rows: vec![
                    vec![
                        Value::Float(310.0),
                        Value::Int(5),
                        Value::Str("FR".to_owned()),
                        Value::Null,
                    ],
                    vec![
                        Value::Float(f64::NAN),
                        Value::Int(0),
                        Value::Bool(true),
                        Value::Str("S".to_owned()),
                    ],
                ],
                histogram: Some(HistogramSection {
                    key_label: "insert_size",
                    key_class,
                    series: vec![
                        HistogramSeries {
                            label: "All_Reads.fr_count".to_owned(),
                            histogram: counts.clone(),
                        },
                        HistogramSeries {
                            label: "All_Reads.rf_count".to_owned(),
                            histogram: Histogram::new(),
                        },
                    ],
                }),
                with_headers: true,
                text: None,
            };
            let headers = vec!["CollectInsertSizeMetrics INPUT=x".to_owned()];
            let mut expected = Vec::new();
            output
                .write(headers.clone(), &mut expected)
                .expect("write to Vec");
            let doc = Document::from_output("CollectInsertSizeMetrics", headers, &output);
            assert_eq!(
                write(&doc).expect("writes"),
                String::from_utf8(expected).expect("utf8")
            );
        }
    }
}

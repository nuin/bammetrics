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
    let doc = TextDoc::parse(text).map_err(|e| Error::Metrics(e.to_string()))?;
    if doc.classes.is_empty() && doc.histogram.is_none() {
        return Err(Error::Metrics(
            "not a metrics file: no ## METRICS CLASS and no ## HISTOGRAM section".to_owned(),
        ));
    }
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

/// `doc` as a Picard metrics file; the verbatim text for a file that is
/// not one.
pub fn write(doc: &Document) -> String {
    if let Some(text) = &doc.verbatim {
        return text.clone();
    }
    let text = TextDoc {
        headers: doc.headers.clone(),
        classes: doc.tables.iter().map(|t| block(t, false)).collect(),
        histogram: doc.histogram.as_ref().map(|t| block(t, true)),
    };
    let mut bytes = Vec::new();
    text.write(&mut bytes).expect("writing to memory");
    String::from_utf8_lossy(&bytes).into_owned()
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
        assert_eq!(write(&doc), INSERT_SIZE);
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
        assert_eq!(write(&doc), text);
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
            assert_eq!(write(&doc), String::from_utf8(expected).expect("utf8"));
        }
    }
}

//! Picard `MetricsFile` text output (`htsjdk.samtools.metrics.MetricsFile.write`).
//!
//! Layout: `## htsjdk.samtools.metrics.StringHeader` lines each followed by a
//! `# ` text line, a blank line, `## METRICS CLASS\t<class>`, a tab-separated
//! column row, one line per metric, a blank line, then optionally
//! `## HISTOGRAM\tjava.lang.Integer`, the histogram table and a blank line.

use std::collections::BTreeSet;
use std::io::{self, Write};

use crate::histogram::Histogram;

/// htsjdk `FormatUtil.format(double)`: a `DecimalFormat` with grouping off,
/// at most 6 fraction digits, no trailing zeros, and NaN and infinity both
/// printed as `?`. `DecimalFormat` rounds the exact binary value half-even;
/// Rust's `{:.6}` also rounds the exact value, and a tie at the 7th decimal
/// cannot occur for a binary double, so the two agree digit for digit.
pub fn format_double(value: f64) -> String {
    if !value.is_finite() {
        return "?".to_owned();
    }
    let fixed = format!("{value:.6}");
    let trimmed = fixed.trim_end_matches('0').trim_end_matches('.');
    trimmed.to_owned()
}

/// One histogram column: its header label and its bins.
pub struct HistogramColumn<'a> {
    pub label: String,
    pub histogram: &'a Histogram,
}

/// `MetricsFile.writeHistogram` for several histograms sharing a key column:
/// the rows are the union of all keys in ascending order, and a histogram
/// without a bin at a key prints 0.
pub struct HistogramTable<'a> {
    pub key_label: &'a str,
    pub columns: Vec<HistogramColumn<'a>>,
}

/// A metrics file ready to write.
pub struct MetricsFile<'a> {
    /// Text of each `StringHeader`, written after `# `.
    pub headers: Vec<String>,
    /// Fully qualified Picard metrics class, e.g. `picard.analysis.InsertSizeMetrics`.
    pub class: &'a str,
    pub columns: &'a [&'a str],
    /// Values already formatted; empty strings are Picard's null fields.
    pub rows: Vec<Vec<String>>,
    pub histogram: Option<HistogramTable<'a>>,
}

impl MetricsFile<'_> {
    pub fn write(&self, out: &mut impl Write) -> io::Result<()> {
        for header in &self.headers {
            writeln!(out, "## htsjdk.samtools.metrics.StringHeader")?;
            writeln!(out, "# {header}")?;
        }
        writeln!(out)?;
        writeln!(out, "## METRICS CLASS\t{}", self.class)?;
        writeln!(out, "{}", self.columns.join("\t"))?;
        for row in &self.rows {
            writeln!(out, "{}", row.join("\t"))?;
        }
        writeln!(out)?;
        if let Some(table) = &self.histogram {
            table.write(out)?;
        }
        Ok(())
    }
}

impl HistogramTable<'_> {
    fn write(&self, out: &mut impl Write) -> io::Result<()> {
        writeln!(out, "## HISTOGRAM\tjava.lang.Integer")?;
        let labels: Vec<&str> = self.columns.iter().map(|c| c.label.as_str()).collect();
        writeln!(out, "{}\t{}", self.key_label, labels.join("\t"))?;
        let keys: BTreeSet<i64> = self
            .columns
            .iter()
            .flat_map(|c| c.histogram.keys())
            .collect();
        for key in keys {
            let values: Vec<String> = self
                .columns
                .iter()
                .map(|c| format_double(c.histogram.get(key).unwrap_or(0.0)))
                .collect();
            writeln!(out, "{key}\t{}", values.join("\t"))?;
        }
        writeln!(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn doubles_print_like_format_util() {
        assert_eq!(format_double(300.0), "300");
        assert_eq!(format_double(310.5), "310.5");
        assert_eq!(format_double(297.297_297_297), "297.297297");
        assert_eq!(format_double(16.439_899_4), "16.439899");
        assert_eq!(format_double(0.123_456_500_1), "0.123457");
        assert_eq!(format_double(0.0), "0");
        assert_eq!(format_double(f64::NAN), "?");
        assert_eq!(format_double(f64::INFINITY), "?");
        assert_eq!(format_double(-2.5), "-2.5");
    }

    #[test]
    fn histogram_table_unions_keys() {
        let mut a = Histogram::new();
        a.increment(200);
        a.increment_by(300, 36.0);
        let mut b = Histogram::new();
        b.increment_by(400, 2.0);
        let file = MetricsFile {
            headers: vec!["X INPUT=in.bam".to_owned()],
            class: "picard.analysis.InsertSizeMetrics",
            columns: &["A", "B"],
            rows: vec![vec!["1".to_owned(), String::new()]],
            histogram: Some(HistogramTable {
                key_label: "insert_size",
                columns: vec![
                    HistogramColumn {
                        label: "All_Reads.fr_count".to_owned(),
                        histogram: &a,
                    },
                    HistogramColumn {
                        label: "All_Reads.rf_count".to_owned(),
                        histogram: &b,
                    },
                ],
            }),
        };
        let mut out = Vec::new();
        file.write(&mut out).expect("write to Vec");
        let text = String::from_utf8(out).expect("utf8");
        assert_eq!(
            text,
            "## htsjdk.samtools.metrics.StringHeader\n# X INPUT=in.bam\n\n\
             ## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\nA\tB\n1\t\n\n\
             ## HISTOGRAM\tjava.lang.Integer\n\
             insert_size\tAll_Reads.fr_count\tAll_Reads.rf_count\n\
             200\t1\t0\n300\t36\t0\n400\t0\t2\n\n"
        );
    }
}

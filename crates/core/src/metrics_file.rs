//! Picard `MetricsFile` text output (`htsjdk.samtools.metrics.MetricsFile.write`)
//! and the reader for it ([`MetricsDoc`]).
//!
//! Layout: `## htsjdk.samtools.metrics.StringHeader` lines each followed by a
//! `# ` text line, a blank line, `## METRICS CLASS\t<class>`, a tab-separated
//! column row, one line per metric, a blank line, then optionally
//! `## HISTOGRAM\tjava.lang.Integer`, the histogram table and a blank line.

use std::collections::BTreeSet;
use std::fs;
use std::io::{self, Write};
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::Value;

/// htsjdk `FormatUtil.format(double)`: a `DecimalFormat` with grouping off,
/// at most 6 fraction digits and no trailing zeros. NaN and positive infinity
/// print as `?`, negative infinity as `-?`, and a negative value that rounds
/// to zero as `-0`.
///
/// `DecimalFormat` does not round the exact binary value. It takes the
/// shortest decimal digits that round-trip the double (what `Double.toString`
/// shows), rounds those at the 6th fraction digit, and pads big values with
/// zeros: 1.2345678901234566E25 prints as those 17 digits followed by zeros,
/// where an exact expansion would print other digits. When the shortest digits
/// end in an exact tie (a final 5 right after the kept digits) the JDK breaks
/// it with the exact binary value, rounding up if the double is above the tie
/// and down if below. A double exactly on the tie (123/128 = 0.9609375)
/// rounds down, because htsjdk sets `RoundingMode.HALF_DOWN`.
pub fn format_double(value: f64) -> String {
    if value.is_nan() {
        return "?".to_owned();
    }
    if value.is_infinite() {
        return if value < 0.0 { "-?" } else { "?" }.to_owned();
    }
    let (integer, fraction) = shortest_decimal(value.abs());
    let (integer, fraction) = round_at(value.abs(), integer, fraction, 6);
    let fraction = fraction.trim_end_matches('0');
    let sign = if value.is_sign_negative() { "-" } else { "" };
    if fraction.is_empty() {
        format!("{sign}{integer}")
    } else {
        format!("{sign}{integer}.{fraction}")
    }
}

/// Integer and fraction digits of the shortest round-trip decimal form.
fn shortest_decimal(value: f64) -> (String, String) {
    // `{:e}` prints the shortest round-trip digits as `d.ddde<exp>`.
    let scientific = format!("{value:e}");
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let mut digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    prefer_even_on_tie(value, &mut digits);
    // Rust's own output always has a valid exponent; 0 is only a fallback.
    let exponent: i64 = exponent.parse().unwrap_or(0);
    let point = exponent + 1;
    if point <= 0 {
        let zeros = "0".repeat(point.unsigned_abs() as usize);
        ("0".to_owned(), format!("{zeros}{digits}"))
    } else if point as usize >= digits.len() {
        let zeros = "0".repeat(point as usize - digits.len());
        (format!("{digits}{zeros}"), String::new())
    } else {
        let (integer, fraction) = digits.split_at(point as usize);
        (integer.to_owned(), fraction.to_owned())
    }
}

/// A double whose exact value ends in a 5 one digit past the shortest
/// round-trip digits sits between two candidates that both round-trip. Rust
/// takes the upper one, the JDK (`FloatingDecimal`) the one with the even last
/// digit, so an odd last digit moves down by one when the exact value is that
/// tie. Only large values, with 16 or 17 digits, reach it.
fn prefer_even_on_tie(value: f64, digits: &mut String) {
    let Some(last) = digits.bytes().last() else {
        return;
    };
    let n = digits.len();
    if (last - b'0').is_multiple_of(2) {
        return;
    }
    let exact = format!("{value:.40e}");
    let exact: Vec<u8> = exact
        .split_once('e')
        .map_or(&exact[..], |(m, _)| m)
        .bytes()
        .filter(u8::is_ascii_digit)
        .collect();
    let lower = last - 1;
    if exact[..n - 1] == digits.as_bytes()[..n - 1]
        && exact[n - 1] == lower
        && exact[n] == b'5'
        && exact[n + 1..].iter().all(|&b| b == b'0')
    {
        digits.pop();
        digits.push(char::from(lower));
    }
}

/// Round the shortest decimal `integer.fraction` of `value` to `places`
/// fraction digits as `DecimalFormat` does (see `format_double`).
fn round_at(value: f64, integer: String, fraction: String, places: usize) -> (String, String) {
    if fraction.len() <= places {
        return (integer, fraction);
    }
    let rest = &fraction.as_bytes()[places..];
    let mut kept: Vec<u8> = integer.bytes().chain(fraction[..places].bytes()).collect();
    let tie = rest[0] == b'5' && rest[1..].iter().all(|&b| b == b'0');
    let round_up = if tie {
        // A double's exact expansion ends within 80 places for any value that
        // can sit on a tie at the 7th decimal.
        let exact = format!("{value:.80}");
        let tie_text = format!("{integer}.{fraction:0<80}");
        match exact.as_str().cmp(tie_text.as_str()) {
            std::cmp::Ordering::Greater => true,
            // `RoundingMode.HALF_DOWN`: a double exactly on the tie rounds down.
            std::cmp::Ordering::Less | std::cmp::Ordering::Equal => false,
        }
    } else {
        rest[0] >= b'5'
    };
    if round_up {
        let mut carry = true;
        for digit in kept.iter_mut().rev() {
            if *digit == b'9' {
                *digit = b'0';
            } else {
                *digit += 1;
                carry = false;
                break;
            }
        }
        if carry {
            kept.insert(0, b'1');
        }
    }
    let split = kept.len() - places;
    let text = String::from_utf8_lossy(&kept).into_owned();
    (text[..split].to_owned(), text[split..].to_owned())
}

/// `StringHeader` text for the second header Picard writes, `Started on: ` and
/// the time in `java.util.Date.toString` form (UTC here, where Picard prints
/// the JVM's zone).
pub fn started_on(seconds_since_epoch: u64) -> String {
    const DAYS: [&str; 7] = ["Thu", "Fri", "Sat", "Sun", "Mon", "Tue", "Wed"];
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let days = seconds_since_epoch / 86_400;
    let in_day = seconds_since_epoch % 86_400;
    // Civil date from a day count (Hinnant's days_from_civil, inverted).
    let shifted = days + 719_468;
    let era = shifted / 146_097;
    let day_of_era = shifted % 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let march_month = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * march_month + 2) / 5 + 1;
    let month = if march_month < 10 {
        march_month + 3
    } else {
        march_month - 9
    };
    let year = year_of_era + era * 400 + u64::from(month <= 2);
    format!(
        "Started on: {} {} {day:02} {:02}:{:02}:{:02} UTC {year}",
        DAYS[(days % 7) as usize],
        MONTHS[(month - 1) as usize],
        in_day / 3600,
        in_day % 3600 / 60,
        in_day % 60
    )
}

/// The two `StringHeader` lines Picard writes: the command line, then
/// `Started on` with the current time.
pub fn headers(command: String) -> Vec<String> {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs());
    vec![command, started_on(now)]
}

/// The Java class of a histogram's keys, printed after `## HISTOGRAM`.
/// Keys are stored as integers; `Double` keys print with a `.0`
/// (`Double.toString` of a whole number).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum KeyClass {
    #[default]
    Integer,
    Byte,
    Double,
}

impl KeyClass {
    pub fn java_name(self) -> &'static str {
        match self {
            KeyClass::Integer => "java.lang.Integer",
            KeyClass::Byte => "java.lang.Byte",
            KeyClass::Double => "java.lang.Double",
        }
    }

    pub fn format(self, key: i64) -> String {
        match self {
            KeyClass::Integer | KeyClass::Byte => key.to_string(),
            KeyClass::Double => format!("{key}.0"),
        }
    }
}

/// One histogram column: its header label and its bins.
pub struct HistogramColumn<'a> {
    pub label: String,
    pub histogram: &'a Histogram,
}

/// `MetricsFile.printHistogram` for several histograms sharing a key column:
/// the rows are the union of all keys in ascending order, and a histogram
/// without a bin at a key prints 0.
pub struct HistogramTable<'a> {
    pub key_label: &'a str,
    pub key_class: KeyClass,
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
        // `printBeanMetrics` prints nothing, not even the class line, without
        // a row; a histogram-only file has two blank lines here.
        if !self.rows.is_empty() {
            writeln!(out, "## METRICS CLASS\t{}", self.class)?;
            writeln!(out, "{}", self.columns.join("\t"))?;
            for row in &self.rows {
                writeln!(out, "{}", row.join("\t"))?;
            }
        }
        writeln!(out)?;
        if let Some(table) = &self.histogram {
            table.write(out)?;
        }
        // `MetricsFile.write` ends with a blank line whether or not a
        // histogram was printed.
        writeln!(out)
    }
}

impl HistogramTable<'_> {
    /// `MetricsFile.printHistogram`: empty histograms are left out, and with
    /// none left nothing is printed, not even the `## HISTOGRAM` line.
    fn write(&self, out: &mut impl Write) -> io::Result<()> {
        let columns: Vec<&HistogramColumn> = self
            .columns
            .iter()
            .filter(|c| !c.histogram.is_empty())
            .collect();
        if columns.is_empty() {
            return Ok(());
        }
        writeln!(out, "## HISTOGRAM\t{}", self.key_class.java_name())?;
        let labels: Vec<&str> = columns.iter().map(|c| c.label.as_str()).collect();
        writeln!(out, "{}\t{}", self.key_label, labels.join("\t"))?;
        let keys: BTreeSet<i64> = columns.iter().flat_map(|c| c.histogram.keys()).collect();
        for key in keys {
            let values: Vec<String> = columns
                .iter()
                .map(|c| format_double(c.histogram.get(key).unwrap_or(0.0)))
                .collect();
            writeln!(out, "{}\t{}", self.key_class.format(key), values.join("\t"))?;
        }
        Ok(())
    }
}

const METRICS_CLASS_LINE: &str = "## METRICS CLASS\t";
const HISTOGRAM_LINE: &str = "## HISTOGRAM\t";
const MAJOR_HEADER: &str = "## ";
const MINOR_HEADER: &str = "# ";

/// A parsed metrics file (`MetricsFile.read`, without the Java classes).
///
/// Columns are kept by name, so a class the reader has never seen, or a
/// column Picard added later, parses like any other. Cells are typed from
/// their text with [`parse_value`]: a report asks for a column by name and
/// gets an int, a float, a string or null.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetricsDoc {
    /// The `# ` text under each `## <header class>` line, in file order.
    pub headers: Vec<String>,
    /// Every `## METRICS CLASS` block, in file order.
    pub classes: Vec<ClassBlock>,
    /// The `## HISTOGRAM` section, when the file has one.
    pub histogram: Option<HistogramBlock>,
}

/// One `## METRICS CLASS` block: the class name, its column line and one
/// `Value` per column per row.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClassBlock {
    pub class: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<Value>>,
}

/// The `## HISTOGRAM` section: the key's Java class, the label line (the key
/// label first) and one row per key. Keys are typed by their text; the
/// series cells are floats, as htsjdk reads them, or null when empty.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HistogramBlock {
    pub key_class: String,
    pub columns: Vec<String>,
    pub rows: Vec<Vec<Value>>,
}

impl MetricsDoc {
    /// Parse the text of a metrics file.
    ///
    /// Tolerated: `\r\n` line ends, a missing blank line between sections,
    /// a row shorter than the column line (the rest is null), a histogram
    /// with labels and no rows, and lines outside the known sections.
    /// Errors: a `## METRICS CLASS` or `## HISTOGRAM` line with no column
    /// line after it, a row with more cells than columns, and a second
    /// `## HISTOGRAM` section.
    pub fn parse(text: &str) -> Result<Self> {
        let mut doc = Self::default();
        let mut lines = text.lines().map(|l| l.trim_end_matches('\r')).peekable();
        while let Some(line) = lines.next() {
            if let Some(class) = line.strip_prefix(METRICS_CLASS_LINE) {
                let class = class.trim().to_owned();
                let (columns, rows) = read_table(&mut lines, &class, |_, cell| parse_value(cell))?;
                doc.classes.push(ClassBlock {
                    class,
                    columns,
                    rows,
                });
            } else if let Some(key_class) = line.strip_prefix(HISTOGRAM_LINE) {
                if doc.histogram.is_some() {
                    return Err(bad("a second ## HISTOGRAM section".to_owned()));
                }
                let key_class = key_class.trim().to_owned();
                let (columns, rows) = read_table(&mut lines, "HISTOGRAM", |i, cell| {
                    if i == 0 {
                        parse_value(cell)
                    } else {
                        parse_series_value(cell)
                    }
                })?;
                doc.histogram = Some(HistogramBlock {
                    key_class,
                    columns,
                    rows,
                });
            } else if line.starts_with(MAJOR_HEADER) {
                // `## <header class>` then `# <text>`; a header class without
                // its text line contributes nothing.
                if let Some(text) = lines.next_if(|l| l.starts_with(MINOR_HEADER)) {
                    doc.headers.push(text[MINOR_HEADER.len()..].to_owned());
                }
            }
        }
        Ok(doc)
    }

    /// Read and parse the file at `path`.
    pub fn read(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        Self::parse(&text).map_err(|e| Error::Metric(format!("{}: {e}", path.display())))
    }

    /// The first block of `class`, by the full Java name
    /// (`picard.analysis.InsertSizeMetrics`) or the name after the last `.`
    /// or `$` (`InsertSizeMetrics`, `QualityYieldMetrics`).
    pub fn class(&self, class: &str) -> Option<&ClassBlock> {
        self.classes
            .iter()
            .find(|b| b.class == class || b.class.rsplit(['.', '$']).next() == Some(class))
    }
}

impl ClassBlock {
    /// Index of `column` in the column line.
    pub fn column(&self, column: &str) -> Option<usize> {
        self.columns.iter().position(|c| c == column)
    }

    /// The cell at `row` under `column`; `None` when either is absent.
    pub fn get(&self, row: usize, column: &str) -> Option<&Value> {
        self.rows.get(row)?.get(self.column(column)?)
    }

    /// The first row whose `column` prints as `text` (`CATEGORY=PAIR`).
    pub fn find(&self, column: &str, text: &str) -> Option<&[Value]> {
        let index = self.column(column)?;
        self.rows
            .iter()
            .find(|row| row[index].picard_text() == text)
            .map(Vec::as_slice)
    }

    /// Rows as `(column, value)` pairs.
    pub fn records(&self) -> impl Iterator<Item = impl Iterator<Item = (&str, &Value)>> {
        zip_records(&self.columns, &self.rows)
    }
}

/// Each row zipped with the column names.
fn zip_records<'a>(
    columns: &'a [String],
    rows: &'a [Vec<Value>],
) -> impl Iterator<Item = impl Iterator<Item = (&'a str, &'a Value)>> {
    rows.iter()
        .map(move |row| columns.iter().map(String::as_str).zip(row.iter()))
}

impl HistogramBlock {
    /// Index of `label` in the label line (0 is the key).
    pub fn column(&self, label: &str) -> Option<usize> {
        self.columns.iter().position(|c| c == label)
    }

    /// `(key, value)` down one series, skipping null cells.
    pub fn series(&self, label: &str) -> Option<Vec<(&Value, f64)>> {
        let index = self.column(label)?;
        Some(
            self.rows
                .iter()
                .filter_map(|row| row[index].as_float().map(|v| (&row[0], v)))
                .collect(),
        )
    }

    /// Rows as `(label, value)` pairs.
    pub fn records(&self) -> impl Iterator<Item = impl Iterator<Item = (&str, &Value)>> {
        zip_records(&self.columns, &self.rows)
    }
}

fn bad(message: String) -> Error {
    Error::Metric(message)
}

fn split_cells(line: &str) -> Vec<String> {
    line.split('\t').map(str::to_owned).collect()
}

/// The column line and the rows under a `## METRICS CLASS` or `## HISTOGRAM`
/// line: rows run to a blank line, a `## ` line or the end. `cell` types
/// the text at column index `i`.
fn read_table<'a>(
    lines: &mut std::iter::Peekable<impl Iterator<Item = &'a str>>,
    what: &str,
    cell: impl Fn(usize, &str) -> Value,
) -> Result<(Vec<String>, Vec<Vec<Value>>)> {
    let columns = split_cells(
        lines
            .next()
            .ok_or_else(|| bad(format!("{what} has no column line")))?,
    );
    let mut rows = Vec::new();
    while let Some(row) = lines.next_if(|l| !l.is_empty() && !l.starts_with(MAJOR_HEADER)) {
        let cells = split_cells(row);
        if cells.len() > columns.len() {
            return Err(bad(format!(
                "{what} row {} has {} cells for {} columns",
                rows.len() + 1,
                cells.len(),
                columns.len()
            )));
        }
        let mut values: Vec<Value> = cells.iter().enumerate().map(|(i, c)| cell(i, c)).collect();
        values.resize(columns.len(), Value::Null);
        rows.push(values);
    }
    Ok((columns, rows))
}

/// A cell typed from its text. Empty is null; `?` is NaN and `-?` negative
/// infinity (htsjdk `FormatUtil` prints both NaN and infinity as `?`); an
/// integer in canonical form (`-12`, not `012` or `+12`) is an int; a decimal
/// with optional sign, one `.` or `,` and an optional exponent is a float;
/// anything else is the text. The `,` form is what `DecimalFormat` writes
/// under a locale with a comma decimal separator.
pub fn parse_value(text: &str) -> Value {
    if text.is_empty() {
        return Value::Null;
    }
    match text {
        "?" => return Value::Float(f64::NAN),
        "-?" => return Value::Float(f64::NEG_INFINITY),
        // A negative double that rounds to zero at 6 places.
        "-0" => return Value::Float(-0.0),
        _ => {}
    }
    if let Ok(int) = text.parse::<i64>()
        && int.to_string() == text
    {
        return Value::Int(int);
    }
    parse_decimal(text).map_or_else(|| Value::Str(text.to_owned()), Value::Float)
}

/// A histogram series cell: null when empty, else [`parse_value`] widened to
/// a float, with text that is no number kept as text.
fn parse_series_value(text: &str) -> Value {
    match parse_value(text) {
        Value::Int(v) => Value::Float(v as f64),
        other => other,
    }
}

/// `[+-]digits[(.|,)digits][(e|E)[+-]digits]` with at least one mantissa
/// digit and a separator or an exponent, as a float. `None` for any other
/// text, including `inf`, `nan` and `012`, which `f64::from_str` accepts.
fn parse_decimal(text: &str) -> Option<f64> {
    let digits = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
    let (mantissa, exponent) = match unsigned.split_once(['e', 'E']) {
        Some((m, e)) => (m, Some(e.strip_prefix(['+', '-']).unwrap_or(e))),
        None => (unsigned, None),
    };
    let (whole, fraction) = match mantissa.split_once(['.', ',']) {
        Some(parts) => parts,
        None if exponent.is_some() => (mantissa, ""),
        None => return None,
    };
    if exponent.is_some_and(|e| !digits(e)) {
        return None;
    }
    let whole_ok = whole.is_empty() || digits(whole);
    let fraction_ok = fraction.is_empty() || digits(fraction);
    if !whole_ok || !fraction_ok || (whole.is_empty() && fraction.is_empty()) {
        return None;
    }
    text.replace(',', ".").parse().ok()
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
        assert_eq!(format_double(-2.5), "-2.5");
    }

    /// Each expectation below was printed by htsjdk 4.x `FormatUtil` in Picard 3.5.0.
    #[test]
    fn edge_cases_match_htsjdk() {
        assert_eq!(format_double(0.0), "0");
        assert_eq!(format_double(1e-7), "0");
        assert_eq!(format_double(0.5), "0.5");
        assert_eq!(format_double(1.0 / 3.0), "0.333333");
        assert_eq!(format_double(1e20), "100000000000000000000");
        assert_eq!(
            format_double(1.234_567_890_123_456_7e25),
            "12345678901234566000000000"
        );
        assert_eq!(format_double(f64::NAN), "?");
        assert_eq!(format_double(f64::INFINITY), "?");
        assert_eq!(format_double(f64::NEG_INFINITY), "-?");
        assert_eq!(format_double(-0.0), "-0");
        assert_eq!(format_double(-1e-7), "-0");
    }

    #[test]
    fn rounding_works_on_shortest_digits() {
        assert_eq!(format_double(5e-7), "0");
        assert_eq!(format_double(1.5e-6), "0.000002");
        assert_eq!(format_double(2.5e-7), "0");
        assert_eq!(format_double(1_234_567.000_000_5), "1234567");
        assert_eq!(format_double(1e15 + 0.3), "1000000000000000.2");
        assert_eq!(format_double(0.999_999_6), "1");
        assert_eq!(format_double(9.999_999_5), "9.999999");
    }

    /// Shortest-digit ties broken by the exact binary value, and 16 or 17
    /// digit ties broken toward the even digit; both checked against 200000
    /// htsjdk outputs.
    #[test]
    fn ties_follow_the_exact_value() {
        assert_eq!(format_double(23.533_862_5), "23.533863");
        assert_eq!(format_double(91.412_967_5), "91.412967");
        assert_eq!(
            format_double(-1.224_008_541_635_281_3e13),
            "-12240085416352.812"
        );
        assert_eq!(format_double(6.404_429_851_067_573e14), "640442985106757.2");
    }

    /// Doubles exactly on a 7th-decimal tie round down (HALF_DOWN), seen in
    /// Picard 3.5.0 output as PCT_PF_READS 0.960937 for 123/128.
    #[test]
    fn exact_binary_ties_round_down() {
        assert_eq!(format_double(123.0 / 128.0), "0.960937");
        assert_eq!(format_double(1.0 / 128.0), "0.007812");
        assert_eq!(format_double(3.0 / 128.0), "0.023437");
        assert_eq!(format_double(-123.0 / 128.0), "-0.960937");
    }

    #[test]
    fn started_on_is_java_date_form() {
        assert_eq!(started_on(0), "Started on: Thu Jan 01 00:00:00 UTC 1970");
        assert_eq!(
            started_on(1_709_208_000),
            "Started on: Thu Feb 29 12:00:00 UTC 2024"
        );
        assert_eq!(
            started_on(1_793_218_309),
            "Started on: Wed Oct 28 20:11:49 UTC 2026"
        );
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
                key_class: KeyClass::Integer,
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
        let doc = MetricsDoc::parse(&text).expect("parses");
        assert_eq!(doc.headers, ["X INPUT=in.bam"]);
        let block = doc.class("InsertSizeMetrics").expect("block");
        assert_eq!(block.columns, ["A", "B"]);
        assert_eq!(block.rows, [[Value::Int(1), Value::Null]]);
        let histogram = doc.histogram.expect("histogram");
        assert_eq!(histogram.key_class, "java.lang.Integer");
        assert_eq!(
            histogram.series("All_Reads.rf_count").expect("series"),
            [
                (&Value::Int(200), 0.0),
                (&Value::Int(300), 0.0),
                (&Value::Int(400), 2.0)
            ]
        );
    }

    /// The first lines of Picard 3.5.0's own output for the committed fixture.
    const PICARD_INSERT_SIZE: &str = "## htsjdk.samtools.metrics.StringHeader\n\
# CollectInsertSizeMetrics HISTOGRAM_FILE=/dev/null INPUT=tests/data/insert_size.bam OUTPUT=t0_fixture.insert_size.txt VERBOSITY=ERROR QUIET=true    DEVIATIONS=10.0 MINIMUM_PCT=0.05 METRIC_ACCUMULATION_LEVEL=[ALL_READS] INCLUDE_DUPLICATES=false ASSUME_SORTED=true STOP_AFTER=0 VALIDATION_STRINGENCY=STRICT COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=500000 CREATE_INDEX=false CREATE_MD5_FILE=false USE_JDK_DEFLATER=false USE_JDK_INFLATER=false\n\
## htsjdk.samtools.metrics.StringHeader\n\
# Started on: Mon Nov 02 18:48:01 MST 2026\n\
\n\
## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n\
MEDIAN_INSERT_SIZE\tMODE_INSERT_SIZE\tMEDIAN_ABSOLUTE_DEVIATION\tMIN_INSERT_SIZE\tMAX_INSERT_SIZE\tMEAN_INSERT_SIZE\tSTANDARD_DEVIATION\tREAD_PAIRS\tPAIR_ORIENTATION\tWIDTH_OF_10_PERCENT\tWIDTH_OF_20_PERCENT\tWIDTH_OF_30_PERCENT\tWIDTH_OF_40_PERCENT\tWIDTH_OF_50_PERCENT\tWIDTH_OF_60_PERCENT\tWIDTH_OF_70_PERCENT\tWIDTH_OF_80_PERCENT\tWIDTH_OF_90_PERCENT\tWIDTH_OF_95_PERCENT\tWIDTH_OF_99_PERCENT\tSAMPLE\tLIBRARY\tREAD_GROUP\n\
310\t310\t10\t300\t500\t310\t8.164966\t5\tFR\t1\t1\t1\t1\t21\t21\t21\t21\t381\t381\t381\t\t\t\n\
\n\
## HISTOGRAM\tjava.lang.Integer\n\
insert_size\tAll_Reads.fr_count\n\
300\t1\n\
310\t2\n\
320\t1\n\
\n";

    #[test]
    fn picard_insert_size_output_parses() {
        let doc = MetricsDoc::parse(PICARD_INSERT_SIZE).expect("parses");
        assert_eq!(doc.headers.len(), 2);
        assert!(doc.headers[0].starts_with("CollectInsertSizeMetrics HISTOGRAM_FILE=/dev/null"));
        assert_eq!(doc.headers[1], "Started on: Mon Nov 02 18:48:01 MST 2026");
        assert_eq!(doc.classes.len(), 1);
        let block = doc
            .class("picard.analysis.InsertSizeMetrics")
            .expect("block");
        assert_eq!(block.columns.len(), 23);
        assert_eq!(block.rows.len(), 1);
        assert_eq!(
            block.get(0, "MEDIAN_INSERT_SIZE").and_then(Value::as_float),
            Some(310.0)
        );
        assert_eq!(block.get(0, "READ_PAIRS"), Some(&Value::Int(5)));
        assert_eq!(
            block.get(0, "STANDARD_DEVIATION"),
            Some(&Value::Float(8.164966))
        );
        assert_eq!(
            block.get(0, "PAIR_ORIENTATION").and_then(Value::as_str),
            Some("FR")
        );
        assert_eq!(block.get(0, "SAMPLE"), Some(&Value::Null));
        assert_eq!(block.get(0, "READ_GROUP"), Some(&Value::Null));
        assert_eq!(block.get(0, "NO_SUCH_COLUMN"), None);
        assert_eq!(block.get(1, "READ_PAIRS"), None);
        assert_eq!(
            block.find("PAIR_ORIENTATION", "FR").map(|r| r[7].clone()),
            Some(Value::Int(5))
        );
        let histogram = doc.histogram.as_ref().expect("histogram");
        assert_eq!(histogram.columns, ["insert_size", "All_Reads.fr_count"]);
        assert_eq!(histogram.rows.len(), 3);
        assert_eq!(histogram.rows[1], [Value::Int(310), Value::Float(2.0)]);
        let records: Vec<Vec<(&str, &Value)>> =
            histogram.records().map(Iterator::collect).collect();
        assert_eq!(records[0][0], ("insert_size", &Value::Int(300)));
    }

    /// `QualityScoreDistribution`: no metrics class, a `Byte`-keyed histogram.
    #[test]
    fn histogram_only_file_has_no_class() {
        let text = "## htsjdk.samtools.metrics.StringHeader\n# QualityScoreDistribution INPUT=x.bam\n\
                    ## htsjdk.samtools.metrics.StringHeader\n# Started on: Mon Nov 02 18:48:01 MST 2026\n\n\n\
                    ## HISTOGRAM\tjava.lang.Byte\nQUALITY\tCOUNT_OF_Q\n2\t15\n30\t985\n\n";
        let doc = MetricsDoc::parse(text).expect("parses");
        assert!(doc.classes.is_empty());
        let histogram = doc.histogram.expect("histogram");
        assert_eq!(histogram.key_class, "java.lang.Byte");
        assert_eq!(
            histogram.series("COUNT_OF_Q").expect("series"),
            [(&Value::Int(2), 15.0), (&Value::Int(30), 985.0)]
        );
        assert_eq!(histogram.series("COUNT_OF_OQ"), None);
    }

    /// `AccumulateQualityYieldMetrics` writes no `StringHeader` lines.
    #[test]
    fn file_without_headers_parses() {
        let text = "\n## METRICS CLASS\tpicard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics\n\
                    TOTAL_READS\tPF_READS\n10\t10\n\n\n";
        let doc = MetricsDoc::parse(text).expect("parses");
        assert!(doc.headers.is_empty());
        assert_eq!(
            doc.class("QualityYieldMetrics")
                .and_then(|b| b.get(0, "PF_READS"))
                .and_then(Value::as_int),
            Some(10)
        );
    }

    #[test]
    fn several_blocks_extra_columns_empty_cells_and_locale_decimals() {
        let text = "## METRICS CLASS\tpicard.analysis.GcBiasSummaryMetrics\r\n\
                    ACCUMULATION_LEVEL\tAT_DROPOUT\tGC_DROPOUT\tNEW_COLUMN\tSAMPLE\r\n\
                    All Reads\t0,5\t?\tfoo\r\n\
                    Sample\t1.25\t-?\t\tS1\r\n\
                    ## METRICS CLASS\tpicard.analysis.GcBiasDetailMetrics\r\n\
                    GC\tWINDOWS\r\n0\t12\r\n1\t\r\n\r\n";
        let doc = MetricsDoc::parse(text).expect("parses");
        assert_eq!(doc.classes.len(), 2);
        let summary = doc.class("GcBiasSummaryMetrics").expect("summary");
        assert_eq!(summary.column("NEW_COLUMN"), Some(3));
        assert_eq!(
            summary.get(0, "ACCUMULATION_LEVEL"),
            Some(&Value::Str("All Reads".to_owned()))
        );
        assert_eq!(summary.get(0, "AT_DROPOUT"), Some(&Value::Float(0.5)));
        assert!(
            summary
                .get(0, "GC_DROPOUT")
                .and_then(Value::as_float)
                .expect("float")
                .is_nan()
        );
        assert_eq!(summary.get(0, "SAMPLE"), Some(&Value::Null));
        assert_eq!(summary.get(1, "AT_DROPOUT"), Some(&Value::Float(1.25)));
        assert_eq!(
            summary.get(1, "GC_DROPOUT"),
            Some(&Value::Float(f64::NEG_INFINITY))
        );
        assert_eq!(summary.get(1, "NEW_COLUMN"), Some(&Value::Null));
        assert_eq!(summary.find("SAMPLE", "S1").map(<[Value]>::len), Some(5));
        let detail = doc.class("GcBiasDetailMetrics").expect("detail");
        assert_eq!(
            detail.rows,
            [
                [Value::Int(0), Value::Int(12)],
                [Value::Int(1), Value::Null]
            ]
        );
        assert!(doc.histogram.is_none());
    }

    #[test]
    fn empty_tables_and_malformed_files() {
        let doc = MetricsDoc::parse(
            "## METRICS CLASS\tX\nA\tB\n\n## HISTOGRAM\tjava.lang.Integer\nk\tv\n",
        )
        .expect("parses");
        assert_eq!(doc.classes[0].rows, Vec::<Vec<Value>>::new());
        assert_eq!(doc.histogram.as_ref().map(|h| h.rows.len()), Some(0));
        assert_eq!(
            doc.histogram.as_ref().and_then(|h| h.series("v")),
            Some(Vec::new())
        );
        assert!(MetricsDoc::parse("## METRICS CLASS\tX\n").is_err());
        assert!(MetricsDoc::parse("## HISTOGRAM\tjava.lang.Integer\n").is_err());
        assert!(MetricsDoc::parse("## METRICS CLASS\tX\nA\n1\t2\n").is_err());
        assert!(
            MetricsDoc::parse(
                "## HISTOGRAM\tjava.lang.Integer\nk\tv\n\n## HISTOGRAM\tjava.lang.Integer\nk\tv\n"
            )
            .is_err()
        );
        let missing = MetricsDoc::read(Path::new("/nonexistent/metrics.txt")).expect_err("missing");
        assert!(missing.to_string().contains("metrics.txt"));
    }

    #[test]
    fn cells_are_typed_from_their_text() {
        assert_eq!(parse_value(""), Value::Null);
        assert_eq!(parse_value("0"), Value::Int(0));
        assert_eq!(parse_value("-12"), Value::Int(-12));
        assert_eq!(parse_value("310.5"), Value::Float(310.5));
        assert_eq!(parse_value("-0"), Value::Float(-0.0));
        assert_eq!(parse_value(".5"), Value::Float(0.5));
        assert_eq!(parse_value("1,5e3"), Value::Float(1500.0));
        assert_eq!(parse_value("1e3"), Value::Float(1000.0));
        assert!(parse_value("?").as_float().expect("float").is_nan());
        assert_eq!(parse_value("-?"), Value::Float(f64::NEG_INFINITY));
        for text in [
            "012", "+5", "1.2.3", "1,234,5", "inf", "nan", "NaN", "1e", "FR", "Y", "S1.LibA", ".",
        ] {
            assert_eq!(parse_value(text), Value::Str(text.to_owned()), "{text}");
        }
        assert_eq!(Value::Int(3).as_float(), Some(3.0));
        assert_eq!(Value::Float(3.0).as_int(), None);
        assert_eq!(Value::Str("x".to_owned()).as_float(), None);
        assert!(Value::Null.is_null());
    }
}

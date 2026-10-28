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
/// and down if below, and half-even only if it is exactly on it.
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
    if (last - b'0') % 2 == 0 {
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
            std::cmp::Ordering::Less => false,
            std::cmp::Ordering::Equal => kept.last().is_some_and(|b| (b - b'0') % 2 == 1),
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
        writeln!(out, "## HISTOGRAM\tjava.lang.Integer")?;
        let labels: Vec<&str> = columns.iter().map(|c| c.label.as_str()).collect();
        writeln!(out, "{}\t{}", self.key_label, labels.join("\t"))?;
        let keys: BTreeSet<i64> = columns.iter().flat_map(|c| c.histogram.keys()).collect();
        for key in keys {
            let values: Vec<String> = columns
                .iter()
                .map(|c| format_double(c.histogram.get(key).unwrap_or(0.0)))
                .collect();
            writeln!(out, "{key}\t{}", values.join("\t"))?;
        }
        Ok(())
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
    fn rounding_is_half_even_on_shortest_digits() {
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

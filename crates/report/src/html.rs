//! HTML escaping and number formatting for the page.

use bammetrics_core::metrics_file::format_double;

/// `text` with `&`, `<`, `>`, `"` and `'` replaced by entities. Every string
/// that comes from a BAM header, a file name or a metrics cell goes through
/// it before it reaches the page.
pub fn escape(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(c),
        }
    }
    out
}

/// An integer with thousands separators: `1,234,567`.
pub fn thousands(value: i64) -> String {
    let digits = value.unsigned_abs().to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3 + 1);
    if value < 0 {
        out.push('-');
    }
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// A number as Picard prints it (`format_double`), with thousands
/// separators when it is a whole number of 10000 or more.
pub fn number(value: f64) -> String {
    if value.fract() == 0.0 && value.abs() >= 10_000.0 && value.abs() < 9.0e15 {
        thousands(value as i64)
    } else {
        format_double(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn escapes_the_five_html_characters() {
        assert_eq!(
            escape("<script>alert('x & \"y\"')</script>"),
            "&lt;script&gt;alert(&#39;x &amp; &quot;y&quot;&#39;)&lt;/script&gt;"
        );
    }

    #[test]
    fn thousands_groups_digits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1000), "1,000");
        assert_eq!(thousands(1234567), "1,234,567");
        assert_eq!(thousands(-1234), "-1,234");
    }

    #[test]
    fn number_keeps_picard_text_below_ten_thousand() {
        assert_eq!(number(0.377778), "0.377778");
        assert_eq!(number(310.0), "310");
        assert_eq!(number(12345.0), "12,345");
        assert_eq!(number(12345.5), "12345.5");
    }
}

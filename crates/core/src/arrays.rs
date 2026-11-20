//! The genotyping-array tools of Picard 3.5.0 (`picard.arrays`), wave W11:
//!
//! - [`baf_regress`] is `CreateBafRegressMetricsFile`, the text bafRegress
//!   prints as `BafRegressMetrics` rows.
//! - [`verify_id_intensity`] is `CreateVerifyIDIntensityContaminationMetricsFile`,
//!   the text VerifyIDIntensity prints as
//!   `VerifyIDIntensityContaminationMetrics` rows.
//! - [`variant_calling`] is `CollectArraysVariantCallingMetrics`, the detail,
//!   summary and control-code files of an arrays VCF as `GtcToVcf` writes it.
//!
//! The two text tools read their input as `BufferedReader.readLine` does,
//! through [`read_lines`]: the file decoded as UTF-8 (a malformed sequence
//! becomes U+FFFD, as Java's decoder replaces it), lines ended by `\n`, `\r` or
//! `\r\n`. The file is read whole under three caps: [`MAX_TEXT_BYTES`]
//! bytes, [`MAX_LINE_BYTES`] bytes per line and [`MAX_LINES`] lines. Each
//! cap is an error that names it; bafRegress and VerifyIDIntensity print
//! one short line per sample.

use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::error::{Error, Result};

pub mod baf_regress;
pub mod header_date;
pub mod variant_calling;
pub mod verify_id_intensity;

/// The largest text input read, in bytes.
pub const MAX_TEXT_BYTES: u64 = 64 << 20;
/// The longest line of a text input, in bytes.
pub const MAX_LINE_BYTES: usize = 1 << 20;
/// The most lines of a text input.
pub const MAX_LINES: usize = 100_000;

/// An error about the content of `path`.
pub(crate) fn content_error(path: &Path, message: String) -> Error {
    Error::Intervals {
        path: path.to_path_buf(),
        message,
    }
}

/// The whole of `path` as text, at most [`MAX_TEXT_BYTES`] bytes, decoded
/// as `FileReader` decodes it.
pub(crate) fn read_text(path: &Path) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    let mut bytes = Vec::new();
    file.take(MAX_TEXT_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(input)?;
    if bytes.len() as u64 > MAX_TEXT_BYTES {
        return Err(content_error(
            path,
            format!("over {MAX_TEXT_BYTES} bytes, the most a text input may hold"),
        ));
    }
    Ok(java_utf8(&bytes))
}

/// `bytes` decoded as Java's UTF-8 decoder decodes them: each malformed
/// sequence becomes one U+FFFD, as in `String::from_utf8_lossy`, except an
/// encoded surrogate (`ED A0..BF 80..BF`), which Java reads as one
/// malformed sequence of three bytes and Rust as three.
pub(crate) fn java_utf8(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    let mut rest = bytes;
    loop {
        match std::str::from_utf8(rest) {
            Ok(text) => {
                out.push_str(text);
                return out;
            }
            Err(error) => {
                let (valid, bad) = rest.split_at(error.valid_up_to());
                out.push_str(std::str::from_utf8(valid).unwrap_or_default());
                out.push('\u{FFFD}');
                let skip = match bad {
                    [0xED, 0xA0..=0xBF, 0x80..=0xBF, ..] => 3,
                    _ => error.error_len().unwrap_or(bad.len()),
                };
                rest = &bad[skip..];
            }
        }
    }
}

/// The lines of `text` as `BufferedReader.readLine` returns them: ended by
/// `\n`, `\r` or `\r\n`, no empty line after a final line end. A line over
/// [`MAX_LINE_BYTES`] or a line past [`MAX_LINES`] is an error.
pub(crate) fn read_lines<'a>(path: &Path, text: &'a str) -> Result<Vec<&'a str>> {
    let mut lines = Vec::new();
    let mut rest = text;
    while !rest.is_empty() {
        let (line, skip) = match rest.find(['\n', '\r']) {
            Some(end) => (
                &rest[..end],
                end + 1 + usize::from(rest[end..].starts_with("\r\n")),
            ),
            None => (rest, rest.len()),
        };
        if line.len() > MAX_LINE_BYTES {
            return Err(content_error(
                path,
                format!("line {} is over {MAX_LINE_BYTES} bytes", lines.len() + 1),
            ));
        }
        if lines.len() == MAX_LINES {
            return Err(content_error(
                path,
                format!("over {MAX_LINES} lines, the most a text input may hold"),
            ));
        }
        lines.push(line);
        rest = &rest[skip..];
    }
    Ok(lines)
}

/// Java's `\s`: space, tab, line feed, vertical tab, form feed, carriage
/// return.
pub(crate) fn is_java_space(c: char) -> bool {
    matches!(c, ' ' | '\t' | '\n' | '\u{0B}' | '\u{0C}' | '\r')
}

/// `line.split("\\s+")`: the text between runs of [`is_java_space`], an
/// empty first field when the line starts with a space, no empty field at
/// the end, and the whole line when it holds no space.
pub(crate) fn split_java_spaces(line: &str) -> Vec<&str> {
    if !line.contains(is_java_space) {
        return vec![line];
    }
    let mut fields = Vec::new();
    if line.starts_with(is_java_space) {
        fields.push("");
    }
    fields.extend(line.split(is_java_space).filter(|f| !f.is_empty()));
    // Only the leading empty field can be left at the end: a line of
    // spaces splits into no field.
    if fields.len() == 1 && fields[0].is_empty() {
        fields.clear();
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lines_end_as_read_line_ends_them() {
        let path = Path::new("x");
        assert_eq!(read_lines(path, "a\r\nb\rc\n").unwrap(), ["a", "b", "c"]);
        assert_eq!(read_lines(path, "a\n\nb").unwrap(), ["a", "", "b"]);
        assert!(read_lines(path, "").unwrap().is_empty());
        let long = "x".repeat(MAX_LINE_BYTES + 1);
        assert!(read_lines(path, &long).is_err());
    }

    #[test]
    fn malformed_bytes_decode_as_java_decodes_them() {
        assert_eq!(java_utf8(b"S\xED\xA0\x80"), "S\u{FFFD}");
        assert_eq!(java_utf8(b"S\xED\xBF\xBFx"), "S\u{FFFD}x");
        assert_eq!(java_utf8(b"\xC0\x80"), "\u{FFFD}\u{FFFD}");
        assert_eq!(java_utf8(b"a\xE2\x82"), "a\u{FFFD}");
        assert_eq!(java_utf8(b"\xFFb\xED\x9F\xBF"), "\u{FFFD}b\u{D7FF}");
        assert_eq!(java_utf8("\u{1F600}".as_bytes()), "\u{1F600}");
        assert_eq!(java_utf8(b"\xF0\x9F\x98"), "\u{FFFD}");
    }

    #[test]
    fn spaces_split_as_java_splits() {
        assert_eq!(split_java_spaces("a\tb  c "), ["a", "b", "c"]);
        assert_eq!(split_java_spaces(" a b"), ["", "a", "b"]);
        assert_eq!(split_java_spaces(""), [""]);
        assert!(split_java_spaces("   ").is_empty());
        assert_eq!(split_java_spaces("abc"), ["abc"]);
        assert_eq!(split_java_spaces("a\u{0B}b\u{A0}c"), ["a", "b\u{A0}c"]);
    }
}

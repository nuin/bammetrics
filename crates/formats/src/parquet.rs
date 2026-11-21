//! The Parquet form of a metrics file (docs/formats.md, docs/07 section
//! 3.2), written and read through the `parquet` crate (feature `parquet`).
//!
//! One metrics file becomes up to three Parquet files beside each other:
//!
//! | File | Holds |
//! |---|---|
//! | `PATH` | the metrics class table, or the histogram of a file without one |
//! | `STEM.histogram.parquet` | the histogram, when the file has a class and a histogram |
//! | `STEM.long.parquet` | every cell as a row: `sample`, `class`, `row_key`, `column`, `value` |
//!
//! `STEM` is `PATH` without a last `.parquet` ([`companion`]). A column is
//! `INT64` for an int, `DOUBLE` for a float and `BYTE_ARRAY` with the
//! `STRING` annotation for a string, the type of [`Table`]: the glossary's.
//! Every column is optional: an empty cell is null, `?` is NaN, and a
//! float keeps all 64 bits. The column names are [`Table::keys`], so a
//! name Picard repeats gets `.1`, `.2`.
//!
//! The key-value metadata key `bammetrics` holds a JSON object: `schema`
//! (`bammetrics.parquet/1`), `part` (`metrics`, `histogram` or `long`),
//! `tool` and `sample` when known, `headers` (the `StringHeader` texts),
//! `class` (the metrics class, or the histogram key's Java class),
//! `picard_columns` when a name repeats, and in the main file `histogram`,
//! true when the histogram companion belongs to it. Picard text is rebuilt
//! from it byte for byte. The key `bammetrics.md5` holds the MD5 of every
//! byte before the footer, checked before a page is decoded, so a damaged
//! file is an error before the `parquet` crate reads it. A document with
//! headers and no table gets a main file of part `empty`: the long file's
//! columns and no row.
//!
//! In the long file `class` is the metrics class (for a histogram cell
//! too; the tool name in a file without a class), `row_key` is
//! `metrics/<i>` or `histogram/<i>` with `i` the row's index from 0 in its
//! table, `column` is the Picard column name and `value` is the cell as
//! text: an int in decimal, a float in the shortest form that reads back as
//! the same double (`310.0`, `1e-7`, `NaN`, `Infinity`, `-Infinity`), a
//! string as it is, null for an empty cell.
//!
//! Files are written with snappy, one data page version 1 per column chunk
//! up to 1 MiB, no dictionary, no statistics, no page index and
//! `created_by` `bammetrics`, so the same document gives the same bytes.

use std::path::{Path, PathBuf};

#[cfg(feature = "parquet")]
use crate::Error;
use crate::table::Ty;
use bammetrics_core::metric::Value;

#[cfg(feature = "parquet")]
mod read;
#[cfg(feature = "parquet")]
mod scan;
#[cfg(feature = "parquet")]
mod write;

#[cfg(feature = "parquet")]
pub(crate) use read::read_capped;
#[cfg(feature = "parquet")]
pub use read::{describe, read};
#[cfg(feature = "parquet")]
pub use write::write;

/// The `schema` member of the metadata.
pub const SCHEMA: &str = "bammetrics.parquet/1";
/// The key-value metadata key.
pub const KEY: &str = "bammetrics";
/// The key of the MD5 of every byte before the footer.
pub const DIGEST_KEY: &str = "bammetrics.md5";
/// The first four bytes (and the last four) of a Parquet file.
pub const MAGIC: &[u8; 4] = b"PAR1";
/// The companions a main file may have, in write order.
pub const COMPANIONS: [&str; 2] = ["histogram", "long"];
/// The largest Parquet file read: 256 MiB.
pub const MAX_BYTES: u64 = 256 << 20;
/// The most cells read from one file, whatever its size: 8,388,608, about
/// 1 GiB of decoded table.
pub const MAX_CELLS: u64 = 1 << 23;
/// The most cells read per byte of the file. A row bammetrics writes
/// holds at least one value of 8 bytes, which snappy shrinks at most about
/// 21 times, so a written file holds well under this; a file of nulls
/// that says millions of rows in a few kilobytes does not.
pub const CELLS_PER_BYTE: u64 = 64;
/// The most bytes the pages of one file may decompress to: 512 MiB.
pub const MAX_DECODED_BYTES: u64 = 512 << 20;

/// True when `bytes` start as a Parquet file does.
pub fn is_parquet(bytes: &[u8]) -> bool {
    bytes.starts_with(MAGIC)
}

/// The companion `part` of the main file `path`: `out.parquet` gives
/// `out.histogram.parquet`, `out.txt` gives `out.txt.histogram.parquet`.
pub fn companion(path: &Path, part: &str) -> PathBuf {
    if path.extension().is_some_and(|e| e == "parquet") {
        return path.with_extension(format!("{part}.parquet"));
    }
    let mut name = path.as_os_str().to_owned();
    name.push(format!(".{part}.parquet"));
    PathBuf::from(name)
}

/// A cell as the long file's `value`: lossless text, `None` for null.
pub fn cell_text(value: &Value) -> Option<String> {
    match value {
        Value::Null => None,
        Value::Int(v) => Some(v.to_string()),
        Value::Float(v) if v.is_nan() => Some("NaN".to_owned()),
        Value::Float(v) if v.is_infinite() => {
            Some(if *v > 0.0 { "Infinity" } else { "-Infinity" }.to_owned())
        }
        Value::Float(v) => Some(format!("{v:?}")),
        Value::Str(s) => Some(s.clone()),
        Value::Bool(b) => Some(if *b { "Y" } else { "N" }.to_owned()),
    }
}

/// The MD5 of `bytes` in hex.
#[cfg(feature = "parquet")]
fn digest(bytes: &[u8]) -> String {
    let mut md5 = bammetrics_core::md5::Md5::new();
    md5.update(bytes);
    bammetrics_core::md5::hex(&md5.finish())
}

/// The Parquet type of a column: `INT64`, `DOUBLE` or `BYTE_ARRAY`
/// (`STRING`).
pub fn type_name(ty: Ty) -> &'static str {
    match ty {
        Ty::Int => "INT64",
        Ty::Float => "DOUBLE",
        Ty::Str => "BYTE_ARRAY (STRING)",
    }
}

/// True when `path` is a histogram companion bammetrics wrote: a Parquet
/// file whose metadata says part `histogram` and no main file of its own.
/// A writer removes such a file beside a new main file that has no
/// histogram, so no stale histogram stays beside it.
pub fn is_histogram_file(path: &Path) -> bool {
    #[cfg(feature = "parquet")]
    return read::companion_part(path).as_deref() == Some("histogram");
    #[cfg(not(feature = "parquet"))]
    {
        let _ = path;
        false
    }
}

/// Without the feature: the error every Parquet entry point gives.
#[cfg(not(feature = "parquet"))]
pub(crate) fn unavailable() -> crate::Error {
    crate::Error::Parquet(
        "this build has no Parquet: build bammetrics with the feature parquet".to_owned(),
    )
}

#[cfg(feature = "parquet")]
fn bad(message: impl Into<String>) -> Error {
    Error::Parquet(message.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn companions_replace_a_last_parquet_extension() {
        let name = |path: &str, part| companion(Path::new(path), part).display().to_string();
        assert_eq!(name("out/x.parquet", "long"), "out/x.long.parquet");
        assert_eq!(
            name("out/x.txt", "histogram"),
            "out/x.txt.histogram.parquet"
        );
        assert_eq!(name("x.parquet.parquet", "long"), "x.parquet.long.parquet");
        assert_eq!(name("x", "long"), "x.long.parquet");
    }

    #[test]
    fn long_values_read_back_as_the_same_cell() {
        let cases = [
            (Value::Int(-3), "-3"),
            (Value::Float(310.0), "310.0"),
            (Value::Float(0.1), "0.1"),
            (Value::Float(1e-7), "1e-7"),
            (Value::Float(f64::NAN), "NaN"),
            (Value::Float(f64::NEG_INFINITY), "-Infinity"),
            (Value::Str("FR".to_owned()), "FR"),
        ];
        for (value, text) in cases {
            assert_eq!(cell_text(&value).as_deref(), Some(text));
            if let Value::Float(v) = value
                && v.is_finite()
            {
                assert_eq!(text.parse::<f64>().expect("a float"), v);
            }
        }
        assert_eq!(cell_text(&Value::Null), None);
    }
}

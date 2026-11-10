//! Output formats of a metrics file (docs/07 section 3, docs/formats.md).
//!
//! Every format is written from one model: a [`Document`] holds a typed
//! [`Table`] per metrics class and the histogram as a second table, with
//! column types from the glossary. Picard text stays the default and is
//! the file the registry writer writes; the others are views of the same
//! cells.
//!
//! | Format | Module | Read back |
//! |---|---|---|
//! | Picard text | [`picard`] | yes |
//! | JSON, schema `bammetrics.metrics/1` | [`metrics_json`] | yes |
//! | JSON Lines, TSV, CSV, aligned text | [`delimited`] | no |
//! | Markdown | [`markdown`] | no |
//! | JUnit XML, MultiQC custom content | [`junit`], [`multiqc`] | no |
//!
//! [`read`] takes either readable form, so Picard text to JSON to Picard
//! text is byte for byte the first file (`tests/round_trip.rs` runs it on
//! every fixture). The crate uses `std` and the core only: Parquet, XLSX
//! and Arrow, which bring dependencies, come behind features (docs/07
//! W21, W46) and write from the same [`Table`].

pub mod delimited;
pub mod error;
pub mod json;
pub mod junit;
pub mod markdown;
pub mod metrics_json;
pub mod multiqc;
pub mod picard;
pub mod table;
pub mod verdict;

pub use error::{Error, Result};
pub use table::{Document, Table, Ty, glossary_type};
pub use verdict::{CheckLine, Outcome, VerdictSummary};

/// A format of one metrics file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Picard,
    Json,
    Tsv,
    Jsonl,
    Human,
    Csv,
    Markdown,
}

impl Format {
    /// Every format, in `--help` order.
    pub const ALL: [Format; 7] = [
        Format::Picard,
        Format::Json,
        Format::Tsv,
        Format::Jsonl,
        Format::Human,
        Format::Csv,
        Format::Markdown,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Format::Picard => "picard",
            Format::Json => "json",
            Format::Tsv => "tsv",
            Format::Jsonl => "jsonl",
            Format::Human => "human",
            Format::Csv => "csv",
            Format::Markdown => "markdown",
        }
    }

    /// The format named `name`; `md` is `markdown`.
    pub fn parse(name: &str) -> Option<Self> {
        match name {
            "md" => Some(Format::Markdown),
            name => Self::ALL.into_iter().find(|f| f.name() == name),
        }
    }

    /// The extension a sibling file gets: `out.txt` in JSON is
    /// `out.txt.json`.
    pub fn extension(self) -> &'static str {
        match self {
            Format::Picard => "txt",
            Format::Json => "json",
            Format::Tsv => "tsv",
            Format::Jsonl => "jsonl",
            Format::Human => "human.txt",
            Format::Csv => "csv",
            Format::Markdown => "md",
        }
    }
}

/// The names `--format` takes, in `--help` order.
pub fn names() -> Vec<&'static str> {
    Format::ALL.into_iter().map(Format::name).collect()
}

/// A `--format` value: one name or several joined by commas. The error is
/// the first name that is no format, or an empty text for a repeated one.
pub fn parse_list(text: &str) -> std::result::Result<Vec<Format>, ListError> {
    let mut formats = Vec::new();
    for name in text.split(',') {
        let format = Format::parse(name).ok_or_else(|| ListError::Unknown(name.to_owned()))?;
        if formats.contains(&format) {
            return Err(ListError::Repeated(format));
        }
        formats.push(format);
    }
    Ok(formats)
}

/// Why a `--format` list does not parse.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ListError {
    Unknown(String),
    Repeated(Format),
}

/// `doc` in `format`. A one-table format prints `table` (the main table
/// without one).
pub fn render(format: Format, doc: &Document, table: Option<&Table>) -> String {
    let one = || table.or_else(|| doc.main_table());
    match format {
        Format::Picard => picard::write(doc),
        Format::Json => metrics_json::write(doc),
        Format::Jsonl => delimited::jsonl(doc),
        Format::Tsv => one().map(delimited::tsv).unwrap_or_default(),
        Format::Csv => one().map(delimited::csv).unwrap_or_default(),
        Format::Human => one().map(delimited::human).unwrap_or_default(),
        Format::Markdown => markdown::write(&[doc], &[]),
    }
}

/// Read a metrics file in either form that reads back: JSON when the
/// first character that is not whitespace is `{`, else Picard text.
pub fn read(text: &str) -> Result<Document> {
    if text.trim_start().starts_with('{') {
        metrics_json::read(text)
    } else {
        picard::read(text)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_lists_parse() {
        assert_eq!(parse_list("json"), Ok(vec![Format::Json]));
        assert_eq!(
            parse_list("picard,json,md"),
            Ok(vec![Format::Picard, Format::Json, Format::Markdown])
        );
        assert_eq!(parse_list("jsn"), Err(ListError::Unknown("jsn".to_owned())));
        assert_eq!(parse_list("json,"), Err(ListError::Unknown(String::new())));
        assert_eq!(
            parse_list("tsv,json,tsv"),
            Err(ListError::Repeated(Format::Tsv))
        );
        for format in Format::ALL {
            assert_eq!(Format::parse(format.name()), Some(format));
        }
        assert_eq!(names()[..5], ["picard", "json", "tsv", "jsonl", "human"]);
    }

    #[test]
    fn read_tells_json_from_picard_text() {
        let picard = "\n## METRICS CLASS\tX\nA\n1\n\n\n";
        let doc = read(picard).expect("picard");
        assert_eq!(render(Format::Picard, &doc, None), picard);
        let json = render(Format::Json, &doc, None);
        let back = read(&format!("  \n{json}")).expect("json");
        assert_eq!(back, doc);
        assert_eq!(render(Format::Tsv, &doc, None), "A\n1\n");
        assert_eq!(render(Format::Csv, &doc, None), "A\r\n1\r\n");
        assert!(read("nothing").is_err());
    }
}

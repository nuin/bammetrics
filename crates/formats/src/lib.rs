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
//! | Parquet (feature `parquet`) | [`parquet`] | yes |
//! | XLSX, DOCX (stored zip of OOXML parts) | [`xlsx`], [`docx`] | no |
//!
//! [`read`] takes either text form, so Picard text to JSON to Picard
//! text is byte for byte the first file (`tests/round_trip.rs` runs it on
//! every fixture); [`read_path`] takes Parquet too. The crate uses `std`
//! and the core only by default; XLSX and DOCX are hand-written zips of
//! XML ([`ooxml`]). Parquet brings the `parquet` crate behind the feature
//! `parquet`; Arrow comes behind a feature (docs/07 W46). All write from
//! the same [`Table`].

pub mod delimited;
pub mod docx;
pub mod error;
pub mod json;
pub mod junit;
pub mod markdown;
pub mod metrics_json;
pub mod multiqc;
pub mod parquet;
pub mod ooxml;
pub mod picard;
pub mod sheet;
pub mod table;
pub mod verdict;
pub mod xlsx;

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
    Parquet,
    Xlsx,
    Docx,
}

impl Format {
    /// Every format, in `--help` order.
    pub const ALL: [Format; 10] = [
        Format::Picard,
        Format::Json,
        Format::Tsv,
        Format::Jsonl,
        Format::Human,
        Format::Csv,
        Format::Markdown,
        Format::Parquet,
        Format::Xlsx,
        Format::Docx,
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
            Format::Parquet => "parquet",
            Format::Xlsx => "xlsx",
            Format::Docx => "docx",
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
            Format::Parquet => "parquet",
            Format::Xlsx => "xlsx",
            Format::Docx => "docx",
        }
    }

    /// The companion parts a format may write beside its main file
    /// ([`parquet::companion`] names them); empty for a one-file format.
    pub fn companions(self) -> &'static [&'static str] {
        match self {
            Format::Parquet => &parquet::COMPANIONS,
            _ => &[],
        }
    }

    /// True for a format whose files are bytes, not text: [`render_files`]
    /// writes it.
    pub fn is_binary(self) -> bool {
        matches!(self, Format::Parquet | Format::Xlsx | Format::Docx)
    }
}

/// One file of a rendered document: `part` is `None` for the main file,
/// else the companion it is ([`Format::companions`]).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct File {
    pub part: Option<&'static str>,
    pub bytes: Vec<u8>,
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

/// `doc` in a text format. A one-table format prints `table` (the main
/// table without one). `picard` fails on a cell Picard text cannot hold
/// ([`picard::write`]); a binary format is an error here, since
/// [`render_files`] writes it.
pub fn render(format: Format, doc: &Document, table: Option<&Table>) -> Result<String> {
    let one = || table.or_else(|| doc.main_table());
    Ok(match format {
        Format::Parquet => {
            return Err(Error::Parquet(
                "Parquet is binary: render_files writes it".to_owned(),
            ));
        }
        Format::Picard => picard::write(doc)?,
        Format::Json => metrics_json::write(doc),
        Format::Jsonl => delimited::jsonl(doc),
        Format::Tsv => one().map(delimited::tsv).unwrap_or_default(),
        Format::Csv => one().map(delimited::csv).unwrap_or_default(),
        Format::Human => one().map(delimited::human).unwrap_or_default(),
        Format::Markdown => markdown::write(&[doc], &[]),
        Format::Xlsx | Format::Docx => return Err(Error::Binary(format.name())),
    })
}

/// Every file of `doc` in `format`, the main file first. A text format,
/// XLSX and DOCX give one file ([`render_bytes`]); `parquet` gives the
/// main file, the histogram companion when the document has a class and a
/// histogram, and the long file. `table` picks the table of a one-table
/// text format.
pub fn render_files(format: Format, doc: &Document, table: Option<&Table>) -> Result<Vec<File>> {
    if format == Format::Parquet {
        #[cfg(feature = "parquet")]
        return parquet::write(doc);
        #[cfg(not(feature = "parquet"))]
        return Err(parquet::unavailable());
    }
    let bytes = render_bytes(format, doc, table)?;
    Ok(vec![File { part: None, bytes }])
}

/// The metrics file at `path` in any form that reads back: Parquet by
/// its first four bytes (with its histogram companion beside it), else
/// text through [`read`]. A zip (XLSX or DOCX) is refused after its first
/// four bytes.
pub fn read_path(path: &std::path::Path) -> Result<Document> {
    use std::io::Read as _;
    let io = |e: std::io::Error| Error::Io(format!("{}: {e}", path.display()));
    let mut file = std::fs::File::open(path).map_err(io)?;
    let mut head = Vec::with_capacity(4);
    (&mut file).take(4).read_to_end(&mut head).map_err(io)?;
    if parquet::is_parquet(&head) {
        #[cfg(feature = "parquet")]
        return parquet::read(parquet::read_capped(path)?, Some(path));
        #[cfg(not(feature = "parquet"))]
        return Err(parquet::unavailable());
    }
    if head == ZIP_MAGIC {
        return Err(Error::Io(format!("{}: {}", path.display(), TextError::Zip)));
    }
    let mut rest = Vec::new();
    file.read_to_end(&mut rest).map_err(io)?;
    head.extend(rest);
    let text = String::from_utf8(head).map_err(|_| {
        Error::Io(format!(
            "{}: not UTF-8 text and not a Parquet file",
            path.display()
        ))
    })?;
    read(&text)
}

/// `doc` in a one-file format as bytes: the text of [`render`], or the
/// XLSX or DOCX file. XLSX and DOCX fail past their limits
/// ([`Error::Limit`]); `parquet` fails here ([`render_files`] writes it).
pub fn render_bytes(format: Format, doc: &Document, table: Option<&Table>) -> Result<Vec<u8>> {
    match format {
        Format::Xlsx => xlsx::write(doc),
        Format::Docx => docx::write(doc),
        text => render(text, doc, table).map(String::into_bytes),
    }
}

/// `text` without a leading byte order mark (U+FEFF), which some editors
/// and exporters put before UTF-8 text.
fn without_bom(text: &str) -> &str {
    text.strip_prefix('\u{feff}').unwrap_or(text)
}

/// True when `text` is read as JSON: after a byte order mark and
/// whitespace its first character is `{` or `[`. No Picard file starts
/// with either.
pub fn is_json(text: &str) -> bool {
    without_bom(text).trim_start().starts_with(['{', '['])
}

/// Read a metrics file in either form that reads back: JSON by
/// [`is_json`], else Picard text. A leading byte order mark is dropped.
pub fn read(text: &str) -> Result<Document> {
    let text = without_bom(text);
    if is_json(text) {
        metrics_json::read(text)
    } else {
        picard::read(text)
    }
}

/// The zip magic an XLSX or DOCX file starts with.
pub const ZIP_MAGIC: &[u8] = b"PK\x03\x04";

/// Why [`read_text`] gave no text.
#[derive(Debug)]
pub enum TextError {
    /// The input starts as a zip: XLSX and DOCX do not read back.
    Zip,
    /// A read error, or bytes that are not UTF-8.
    Io(std::io::Error),
}

impl std::fmt::Display for TextError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TextError::Zip => f.write_str("a zip (XLSX or DOCX) does not read back"),
            TextError::Io(e) => e.fmt(f),
        }
    }
}

/// The UTF-8 text of `reader`. Its first four bytes are read first, so a
/// zip is refused before the rest of it is read.
pub fn read_text(mut reader: impl std::io::Read) -> std::result::Result<String, TextError> {
    use std::io::Read as _;
    let mut bytes = Vec::new();
    reader
        .by_ref()
        .take(ZIP_MAGIC.len() as u64)
        .read_to_end(&mut bytes)
        .map_err(TextError::Io)?;
    if bytes == ZIP_MAGIC {
        return Err(TextError::Zip);
    }
    reader.read_to_end(&mut bytes).map_err(TextError::Io)?;
    String::from_utf8(bytes).map_err(|_| {
        TextError::Io(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "stream did not contain valid UTF-8",
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_zip_is_refused_after_four_bytes() {
        let mut zip = std::io::Cursor::new(b"PK\x03\x04 and more".to_vec());
        assert!(matches!(read_text(&mut zip), Err(TextError::Zip)));
        assert_eq!(zip.position(), 4);
        assert_eq!(read_text(&b"PK"[..]).expect("text"), "PK");
        let error = read_text(&b"\xff\xfe"[..]).expect_err("not UTF-8");
        assert_eq!(error.to_string(), "stream did not contain valid UTF-8");
    }

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
        let render = |format| render(format, &doc, None).expect("renders");
        assert_eq!(render(Format::Picard), picard);
        let json = render(Format::Json);
        let back = read(&format!("  \n{json}")).expect("json");
        assert_eq!(back, doc);
        assert_eq!(render(Format::Tsv), "A\n1\n");
        assert_eq!(render(Format::Csv), "A\r\n1\r\n");
        assert!(read("nothing").is_err());
    }

    /// A byte order mark before either form is dropped, and text that
    /// starts as a JSON array is reported as JSON.
    #[test]
    fn a_byte_order_mark_and_a_json_array_are_told_apart() {
        let picard = "\n## METRICS CLASS\tX\nA\n1\n\n\n";
        let doc = read(picard).expect("picard");
        let json = render(Format::Json, &doc, None).expect("renders");
        assert_eq!(read(&format!("\u{feff}{json}")).expect("json"), doc);
        assert_eq!(read(&format!("\u{feff}{picard}")).expect("picard"), doc);
        assert!(is_json("\u{feff} \n{"));
        assert!(is_json("[1]"));
        assert!(!is_json("## METRICS CLASS\tX"));
        assert_eq!(
            read("[]").expect_err("an array").to_string(),
            "not a metrics document: the JSON is not an object"
        );
        let deep = read(&"[".repeat(200_000))
            .expect_err("too deep")
            .to_string();
        assert!(deep.starts_with("json: "), "{deep}");
    }

    #[test]
    fn only_picard_text_refuses_a_cell() {
        let cell = r#"{"schema": "bammetrics.metrics/1", "class": "X", "columns": ["A"],
            "types": ["string"], "metrics": [{"A": "a\tb"}]}"#;
        let doc = read(cell).expect("reads");
        assert_eq!(render(Format::Tsv, &doc, None).expect("tsv"), "A\na\\tb\n");
        for format in Format::ALL {
            let rendered = render_bytes(format, &doc, None);
            let fails = matches!(format, Format::Picard | Format::Parquet);
            assert_eq!(rendered.is_err(), fails, "{format:?}");
            let text = render(format, &doc, None);
            assert_eq!(
                text.is_err(),
                format == Format::Picard || format.is_binary()
            );
        }
        let error = render(Format::Picard, &doc, None).expect_err("a tab");
        assert!(matches!(error, Error::Picard(_)), "{error}");
    }
}

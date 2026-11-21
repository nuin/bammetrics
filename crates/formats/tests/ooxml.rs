//! XLSX and DOCX read back, on every file Picard 3.5.0 wrote under
//! `tests/data/picard` and on every metrics file a golden record holds.
//!
//! * The round trip: the zip is read entry by entry (each CRC-32 and the
//!   central directory checked), and every cell of every sheet and table
//!   is compared with the document it was written from. XLSX keeps every
//!   double; DOCX holds the text Picard prints.
//! * Outside readers, each skipped when it is not installed: openpyxl,
//!   pandas and python-docx through `ooxml_readers.py` (the Python of
//!   `OOXML_PYTHON`, else `python3`), pandoc (DOCX to Markdown, and the
//!   tables of its JSON), and DuckDB's `excel` extension.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use bammetrics_core::metric::Value;
use bammetrics_formats::json::Json;
use bammetrics_formats::ooxml::crc32;
use bammetrics_formats::{Document, Format, Table, Ty, read, render_bytes};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn files_under(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).expect("list directory") {
        let path = entry.expect("directory entry").path();
        if path.is_dir() {
            files_under(&path, out);
        } else {
            out.push(path);
        }
    }
}

/// True for the text of a Picard metrics file (as `round_trip.rs` tells
/// it): a `StringHeader` line or a blank first line, and a metrics class
/// or a histogram.
fn is_metrics_file(text: &str) -> bool {
    (text.starts_with("## htsjdk.samtools.metrics.StringHeader\n") || text.starts_with('\n'))
        && text
            .lines()
            .any(|l| l.starts_with("## METRICS CLASS\t") || l.starts_with("## HISTOGRAM\t"))
}

/// Standard output and each written file of a golden record.
fn golden_sections(record: &str) -> Vec<(String, String)> {
    let marker = |line: &str| {
        matches!(line, "## argv" | "## exit" | "## stdout" | "## stderr")
            || line.starts_with("## file ")
    };
    let mut sections: Vec<(String, String)> = Vec::new();
    for line in record.split_inclusive('\n') {
        let bare = line.trim_end_matches('\n');
        if marker(bare) {
            sections.push((bare.trim_start_matches("## ").to_owned(), String::new()));
        } else if let Some((_, text)) = sections.last_mut() {
            text.push_str(line);
        }
    }
    sections
        .into_iter()
        .filter(|(name, _)| name == "stdout" || name.starts_with("file "))
        .collect()
}

/// Every metrics document of the fixtures and the golden records, each
/// with a name.
fn documents() -> Vec<(String, Document)> {
    let mut docs = Vec::new();
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    files.sort();
    for path in files {
        let text = fs::read_to_string(&path).expect("read fixture");
        if is_metrics_file(&text) {
            docs.push((
                path.display().to_string(),
                read(&text).expect("a fixture reads"),
            ));
        }
    }
    let mut records = Vec::new();
    files_under(&root().join("tests/golden"), &mut records);
    records.sort();
    for path in records {
        let record = fs::read_to_string(&path).expect("read golden record");
        for (name, text) in golden_sections(&record) {
            if is_metrics_file(&text) {
                let doc = read(&text).expect("a golden metrics file reads");
                docs.push((format!("{} ({name})", path.display()), doc));
            }
        }
    }
    docs
}

/// The entries of a stored zip, in order. Each local header must match
/// its central directory entry and its CRC-32.
fn unzip(bytes: &[u8]) -> Result<Vec<(String, Vec<u8>)>, String> {
    let u16_at = |i: usize| -> Result<usize, String> {
        bytes
            .get(i..i + 2)
            .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
            .ok_or_else(|| format!("cut at {i}"))
    };
    let u32_at = |i: usize| -> Result<usize, String> {
        bytes
            .get(i..i + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
            .ok_or_else(|| format!("cut at {i}"))
    };
    let end = bytes.len().checked_sub(22).ok_or("no end record")?;
    if u32_at(end)? != 0x0605_4b50 {
        return Err("no end record".to_owned());
    }
    let count = u16_at(end + 10)?;
    let mut at = u32_at(end + 16)?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        if u32_at(at)? != 0x0201_4b50 {
            return Err(format!("no directory entry at {at}"));
        }
        let name_length = u16_at(at + 28)?;
        let offset = u32_at(at + 42)?;
        let central_name = &bytes[at + 46..at + 46 + name_length];
        if u32_at(offset)? != 0x0403_4b50 || u16_at(offset + 8)? != 0 {
            return Err(format!("no stored local entry at {offset}"));
        }
        let size = u32_at(offset + 18)?;
        let local_name_length = u16_at(offset + 26)?;
        let start = offset + 30 + local_name_length + u16_at(offset + 28)?;
        let name = &bytes[offset + 30..offset + 30 + local_name_length];
        let data = bytes.get(start..start + size).ok_or("entry cut")?;
        if name != central_name || crc32(data) as usize != u32_at(offset + 14)? {
            return Err(format!("{} does not check", String::from_utf8_lossy(name)));
        }
        entries.push((String::from_utf8_lossy(name).into_owned(), data.to_vec()));
        at += 46 + name_length + u16_at(at + 30)? + u16_at(at + 32)?;
    }
    Ok(entries)
}

/// XML text with its five escapes undone; `ooxml` also undoes `_xHHHH_`.
fn unescape(text: &str, ooxml: bool) -> String {
    let text = text
        .replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&#13;", "\r")
        .replace("&amp;", "&");
    if !ooxml {
        return text;
    }
    let mut out = String::new();
    let mut rest = text.as_str();
    while let Some(i) = rest.find("_x") {
        out.push_str(&rest[..i]);
        let code = rest
            .get(i + 2..i + 7)
            .filter(|s| s.ends_with('_'))
            .and_then(|s| u32::from_str_radix(&s[..4], 16).ok())
            .and_then(char::from_u32);
        match code {
            Some(c) => {
                out.push(c);
                rest = &rest[i + 7..];
            }
            None => {
                out.push_str("_x");
                rest = &rest[i + 2..];
            }
        }
    }
    out.push_str(rest);
    out
}

/// The text between `<tag ...>` (or `<tag>`) and `</tag>` of each element
/// in `xml`.
fn elements<'a>(xml: &'a str, tag: &str) -> Vec<(&'a str, &'a str)> {
    let open = format!("<{tag}");
    let close = format!("</{tag}>");
    let mut found = Vec::new();
    let mut rest = xml;
    while let Some(i) = rest.find(&open) {
        let after = &rest[i + open.len()..];
        if !after.starts_with(['>', ' ']) {
            rest = after;
            continue;
        }
        let head_end = after.find('>').expect("a closed start tag");
        let head = &after[..head_end];
        if head.ends_with('/') {
            found.push((head, ""));
            rest = &after[head_end + 1..];
            continue;
        }
        let body = &after[head_end + 1..];
        let end = body.find(&close).expect("a closing tag");
        found.push((head, &body[..end]));
        rest = &body[end + close.len()..];
    }
    found
}

/// A cell as XLSX holds it.
#[derive(Debug, PartialEq)]
enum Cell {
    Number(String),
    Text(String),
}

/// The cells of a sheet, by row and column; a missing cell is `None`.
fn sheet_cells(xml: &str) -> Vec<Vec<Option<Cell>>> {
    elements(xml, "row")
        .into_iter()
        .map(|(_, row)| {
            let mut cells = Vec::new();
            for (head, body) in elements(row, "c") {
                let reference = head
                    .split("r=\"")
                    .nth(1)
                    .and_then(|r| r.split('"').next())
                    .expect("a cell reference");
                let letters: String = reference
                    .chars()
                    .take_while(char::is_ascii_uppercase)
                    .collect();
                let column = letters
                    .bytes()
                    .fold(0, |n, b| n * 26 + usize::from(b - b'A') + 1)
                    - 1;
                cells.resize_with(column, || None);
                let cell = if head.contains("t=\"inlineStr\"") {
                    let text = elements(body, "t").first().map_or("", |(_, t)| *t);
                    Cell::Text(unescape(text, true))
                } else {
                    Cell::Number(elements(body, "v")[0].1.to_owned())
                };
                cells.push(Some(cell));
            }
            cells
        })
        .collect()
}

/// What XLSX must hold for `value` in a column of type `ty`: an empty
/// text for a null or empty cell of a text column, no cell for one of a
/// number column.
fn expected_cell(value: &Value, ty: Ty) -> Option<Cell> {
    match value {
        Value::Null if ty == Ty::Str => Some(Cell::Text(String::new())),
        Value::Null => None,
        Value::Str(s) if s.is_empty() && ty != Ty::Str => None,
        Value::Int(v) => Some(Cell::Number(v.to_string())),
        Value::Float(v) if v.is_finite() => Some(Cell::Number(format!("{v:?}"))),
        other => Some(Cell::Text(other.picard_text())),
    }
}

/// True when the cell XLSX holds is `value`: an int with its digits, a
/// float as the same double, every other cell as its text.
fn same_cell(held: Option<&Cell>, value: &Value, ty: Ty) -> bool {
    match (held, value) {
        (Some(Cell::Number(text)), Value::Int(v)) => text.parse::<i64>() == Ok(*v),
        (Some(Cell::Number(text)), Value::Float(v)) => {
            (text.contains(['.', 'e']))
                && text
                    .parse::<f64>()
                    .is_ok_and(|p| p.to_bits() == v.to_bits())
        }
        (held, value) => held == expected_cell(value, ty).as_ref(),
    }
}

fn sheet_names(workbook: &str) -> Vec<String> {
    elements(workbook, "sheet")
        .into_iter()
        .map(|(head, _)| {
            let name = head
                .split("name=\"")
                .nth(1)
                .and_then(|n| n.split('"').next());
            unescape(name.expect("a sheet name"), false)
        })
        .collect()
}

/// The tables of a document in sheet order, with the histogram last.
fn tables(doc: &Document) -> Vec<&Table> {
    doc.tables.iter().chain(doc.histogram.as_ref()).collect()
}

/// Every cell of the XLSX of `doc` equal to the document's.
fn xlsx_round_trip(doc: &Document) -> Result<(), String> {
    let bytes = render_bytes(Format::Xlsx, doc, None).map_err(|e| e.to_string())?;
    let entries = unzip(&bytes)?;
    let part = |name: &str| {
        entries
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, data)| String::from_utf8(data.clone()).expect("UTF-8 XML"))
            .ok_or_else(|| format!("no part {name}"))
    };
    let names = sheet_names(&part("xl/workbook.xml")?);
    let tables = tables(doc);
    if names.len() != tables.len() + 1 || names[0] != "Headers" {
        return Err(format!("sheets {names:?}"));
    }
    let lowered: HashSet<String> = names.iter().map(|n| n.to_lowercase()).collect();
    if lowered.len() != names.len() || names.iter().any(|n| n.encode_utf16().count() > 31) {
        return Err(format!("sheet names {names:?}"));
    }
    let headers = sheet_cells(&part("xl/worksheets/sheet1.xml")?);
    let held: Vec<Option<&Cell>> = headers.iter().map(|row| row[0].as_ref()).collect();
    let want: Vec<Cell> = doc.headers.iter().map(|h| Cell::Text(h.clone())).collect();
    if held != want.iter().map(Some).collect::<Vec<_>>() {
        return Err("the Headers sheet differs".to_owned());
    }
    for (i, table) in tables.iter().enumerate() {
        let cells = sheet_cells(&part(&format!("xl/worksheets/sheet{}.xml", i + 2))?);
        let names: Vec<Option<Cell>> = table
            .names()
            .map(|n| Some(Cell::Text(n.to_owned())))
            .collect();
        if cells.first() != Some(&names) || cells.len() != table.rows.len() + 1 {
            return Err(format!(
                "sheet {} has another header row or row count",
                i + 2
            ));
        }
        for (r, row) in table.rows.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                let held = cells[r + 1].get(j).and_then(Option::as_ref);
                if !same_cell(held, value, table.columns[j].1) {
                    return Err(format!(
                        "sheet {} row {} column {}: {held:?} for {value:?}",
                        i + 2,
                        r + 1,
                        table.columns[j].0
                    ));
                }
            }
        }
    }
    Ok(())
}

/// The text of each cell of each table of a DOCX body.
fn docx_tables(document: &str) -> Vec<Vec<Vec<String>>> {
    elements(document, "w:tbl")
        .into_iter()
        .map(|(_, table)| {
            elements(table, "w:tr")
                .into_iter()
                .map(|(_, row)| {
                    elements(row, "w:tc")
                        .into_iter()
                        .map(|(_, cell)| cell_text(cell))
                        .collect()
                })
                .collect()
        })
        .collect()
}

/// A cell's text: its `w:t` runs, `w:tab` as a tab and `w:br` as a line
/// feed, in document order.
fn cell_text(xml: &str) -> String {
    let mut text = String::new();
    let mut rest = xml;
    while let Some(i) = rest.find('<') {
        rest = &rest[i..];
        if let Some(after) = rest.strip_prefix("<w:tab/>") {
            text.push('\t');
            rest = after;
        } else if let Some(after) = rest.strip_prefix("<w:br/>") {
            text.push('\n');
            rest = after;
        } else if rest.starts_with("<w:t>") || rest.starts_with("<w:t ") {
            let start = rest.find('>').expect("closed") + 1;
            let end = rest.find("</w:t>").expect("closed run text");
            text.push_str(&unescape(&rest[start..end], false));
            rest = &rest[end..];
        } else {
            rest = &rest[1..];
        }
    }
    text
}

/// The tables DOCX must hold for `table`: the header row, then the rows,
/// or one row per column for a table of more than 8 columns and at most
/// 8 rows; a wider table cut into parts of at most 63 columns.
fn expected_docx_table(table: &Table) -> Vec<Vec<Vec<String>>> {
    // A line end of any form is one Word line break.
    let text = |v: &Value| v.picard_text().replace("\r\n", "\n").replace('\r', "\n");
    if table.columns.len() > 8 && table.rows.len() <= 8 {
        let mut rows = vec![
            std::iter::once("Column".to_owned())
                .chain((1..=table.rows.len()).map(|i| i.to_string()))
                .collect(),
        ];
        for (j, name) in table.names().enumerate() {
            rows.push(
                std::iter::once(name.to_owned())
                    .chain(table.rows.iter().map(|r| text(&r[j])))
                    .collect(),
            );
        }
        return vec![rows];
    }
    let names: Vec<String> = table.names().map(str::to_owned).collect();
    (0..names.len())
        .step_by(63)
        .map(|first| {
            let last = (first + 63).min(names.len());
            std::iter::once(names[first..last].to_vec())
                .chain(
                    table
                        .rows
                        .iter()
                        .map(|r| r[first..last].iter().map(text).collect()),
                )
                .collect()
        })
        .collect()
}

fn expected_docx(doc: &Document) -> Vec<Vec<Vec<String>>> {
    tables(doc)
        .into_iter()
        .filter(|t| !t.rows.is_empty())
        .flat_map(expected_docx_table)
        .collect()
}

fn docx_round_trip(doc: &Document) -> Result<(), String> {
    let bytes = render_bytes(Format::Docx, doc, None).map_err(|e| e.to_string())?;
    let entries = unzip(&bytes)?;
    let document = entries
        .iter()
        .find(|(n, _)| n == "word/document.xml")
        .map(|(_, d)| String::from_utf8(d.clone()).expect("UTF-8"))
        .ok_or("no word/document.xml")?;
    let held = docx_tables(&document);
    let want = expected_docx(doc);
    if held != want {
        let first = held
            .iter()
            .zip(&want)
            .position(|(h, w)| h != w)
            .unwrap_or(held.len().min(want.len()));
        return Err(format!(
            "table {first} differs ({} held, {} wanted)",
            held.len(),
            want.len()
        ));
    }
    for header in &doc.headers {
        let mut escaped = String::new();
        for c in header.chars() {
            match c {
                '&' => escaped.push_str("&amp;"),
                '<' => escaped.push_str("&lt;"),
                '>' => escaped.push_str("&gt;"),
                '"' => escaped.push_str("&quot;"),
                c => escaped.push(c),
            }
        }
        if !document.contains(&escaped) {
            return Err(format!("no header {header:?}"));
        }
    }
    Ok(())
}

#[test]
fn every_metrics_file_round_trips_through_xlsx_and_docx() {
    let docs = documents();
    let mut failures = Vec::new();
    for (name, doc) in &docs {
        if let Err(e) = xlsx_round_trip(doc) {
            failures.push(format!("{name} xlsx: {e}"));
        }
        if let Err(e) = docx_round_trip(doc) {
            failures.push(format!("{name} docx: {e}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!(
        "{} metrics files round trip through XLSX and DOCX",
        docs.len()
    );
    assert!(docs.len() >= 150, "only {} metrics files", docs.len());
}

/// A document with every edge a cell can have: NaN and both infinities,
/// the extreme ints, a float with 17 digits, markup, tabs and line ends,
/// non-ASCII text, an empty string and a sheet name that collides.
fn planted() -> Document {
    let s = |text: &str| Value::Str(text.to_owned());
    Document {
        tool: Some("CollectX".to_owned()),
        sample: Some("S <1>".to_owned()),
        headers: vec![
            "CollectX INPUT=a&b.bam".to_owned(),
            "Started on: today".to_owned(),
        ],
        tables: vec![
            Table {
                class: "picard.analysis.Headers".to_owned(),
                columns: vec![
                    ("N".to_owned(), Ty::Int),
                    ("F".to_owned(), Ty::Float),
                    ("S".to_owned(), Ty::Str),
                ],
                rows: vec![
                    vec![
                        Value::Int(i64::MAX),
                        Value::Float(f64::NAN),
                        s("a<b>&\"c\""),
                    ],
                    vec![
                        Value::Int(i64::MIN),
                        Value::Float(f64::INFINITY),
                        s("tab\there"),
                    ],
                    vec![
                        Value::Int(0),
                        Value::Float(f64::NEG_INFINITY),
                        s("line\nend\r\nx"),
                    ],
                    vec![Value::Int(-1), Value::Float(0.1 + 0.2), s("Grüße, 😀")],
                    vec![Value::Null, Value::Float(310.0), s("")],
                    vec![Value::Int(5), Value::Float(-1.5e-300), Value::Null],
                ],
            },
            Table {
                class: "x.Y$History".to_owned(),
                columns: vec![("ONLY".to_owned(), Ty::Str)],
                rows: Vec::new(),
            },
            // Wider than the 63 columns of a Word table and long, so DOCX
            // cuts it into parts.
            Table {
                class: "picard.analysis.WideMetrics".to_owned(),
                columns: (0..70).map(|j| (format!("C{j}"), Ty::Int)).collect(),
                rows: (0..10)
                    .map(|i| (0..70).map(|j| Value::Int(i * 100 + j)).collect())
                    .collect(),
            },
        ],
        histogram: Some(Table {
            class: "java.lang.Integer".to_owned(),
            columns: vec![("key".to_owned(), Ty::Int), ("value".to_owned(), Ty::Float)],
            rows: (0..20)
                .map(|k| vec![Value::Int(k), Value::Float(k as f64 / 7.0)])
                .collect(),
        }),
        verbatim: None,
    }
}

#[test]
fn a_planted_document_round_trips() {
    let doc = planted();
    xlsx_round_trip(&doc).expect("xlsx");
    docx_round_trip(&doc).expect("docx");
    let bytes = render_bytes(Format::Xlsx, &doc, None).expect("xlsx");
    let workbook = unzip(&bytes).expect("unzips")[2].1.clone();
    assert_eq!(
        sheet_names(&String::from_utf8(workbook).expect("UTF-8")),
        [
            "Headers",
            "Headers (2)",
            "History (2)",
            "WideMetrics",
            "Histogram"
        ]
    );
    let docx = render_bytes(Format::Docx, &doc, None).expect("docx");
    let body = String::from_utf8_lossy(&docx);
    assert!(body.contains("<w:t>Columns 64 to 70 of 70.</w:t>"));
}

/// A cell an outside reader must give: null, an int, a float or a text
/// (an empty one for a null cell of a text column).
fn reader_cell(value: &Value, ty: Ty) -> Json {
    let tagged = |tag: &str, text: String| Json::Object(vec![(tag.to_owned(), Json::Str(text))]);
    match value {
        Value::Null if ty == Ty::Str => tagged("text", String::new()),
        Value::Null => Json::Null,
        Value::Str(s) if s.is_empty() && ty != Ty::Str => Json::Null,
        Value::Int(v) => tagged("int", v.to_string()),
        Value::Float(v) if v.is_finite() => tagged("float", format!("{v:?}")),
        other => tagged("text", other.picard_text()),
    }
}

/// What `ooxml_readers.py` checks for one document: its files, the sheet
/// names, every sheet's cells, and the DOCX tables.
fn reader_case(xlsx: &Path, docx: &Path, doc: &Document) -> Json {
    let text_cell = |t: &str| Json::Object(vec![("text".to_owned(), Json::text(t))]);
    let mut sheets = vec![Json::Array(
        doc.headers
            .iter()
            .map(|h| Json::Array(vec![text_cell(h)]))
            .collect(),
    )];
    for table in tables(doc) {
        let mut rows = vec![Json::Array(table.names().map(text_cell).collect())];
        rows.extend(table.rows.iter().map(|row| {
            Json::Array(
                row.iter()
                    .zip(&table.columns)
                    .map(|(value, (_, ty))| reader_cell(value, *ty))
                    .collect(),
            )
        }));
        sheets.push(Json::Array(rows));
    }
    let docx_tables = expected_docx(doc)
        .into_iter()
        .map(|t| {
            Json::Array(
                t.into_iter()
                    .map(|row| Json::Array(row.iter().map(|c| Json::text(c)).collect()))
                    .collect(),
            )
        })
        .collect();
    Json::Object(vec![
        ("xlsx".to_owned(), Json::text(&xlsx.display().to_string())),
        ("docx".to_owned(), Json::text(&docx.display().to_string())),
        ("sheets".to_owned(), Json::Array(sheets)),
        ("docx_tables".to_owned(), Json::Array(docx_tables)),
    ])
}

/// The documents the outside readers get: the Picard fixtures, every
/// metrics file of the golden records and the planted document, written
/// into `dir`.
fn write_reader_files(dir: &Path) -> Vec<(PathBuf, PathBuf, Document)> {
    let mut docs = documents();
    docs.push(("planted".to_owned(), planted()));
    docs.into_iter()
        .enumerate()
        .map(|(i, (_, doc))| {
            let xlsx = dir.join(format!("doc{i}.xlsx"));
            let docx = dir.join(format!("doc{i}.docx"));
            fs::write(&xlsx, render_bytes(Format::Xlsx, &doc, None).expect("xlsx")).expect("write");
            fs::write(&docx, render_bytes(Format::Docx, &doc, None).expect("docx")).expect("write");
            (xlsx, docx, doc)
        })
        .collect()
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-ooxml-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

#[test]
fn openpyxl_pandas_and_python_docx_read_every_cell() {
    let python = std::env::var("OOXML_PYTHON").unwrap_or_else(|_| "python3".to_owned());
    let has_modules = Command::new(&python)
        .args(["-c", "import openpyxl, pandas, docx"])
        .output()
        .is_ok_and(|o| o.status.success());
    if !has_modules {
        eprintln!("skipped: {python} lacks openpyxl, pandas or python-docx (set OOXML_PYTHON)");
        return;
    }
    let dir = scratch("python");
    let cases: Vec<Json> = write_reader_files(&dir)
        .iter()
        .map(|(xlsx, docx, doc)| reader_case(xlsx, docx, doc))
        .collect();
    let expect = dir.join("expect.json");
    fs::write(&expect, Json::Array(cases).render()).expect("write expect.json");
    let script = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/ooxml_readers.py");
    let output = Command::new(&python)
        .arg(&script)
        .arg(&expect)
        .output()
        .expect("run python");
    let _ = fs::remove_dir_all(&dir);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(
        output.status.success(),
        "{stdout}{}",
        String::from_utf8_lossy(&output.stderr)
    );
    eprintln!("{}", stdout.trim_end());
}

#[test]
fn pandoc_reads_every_docx_table() {
    let Ok(version) = Command::new("pandoc").arg("--version").output() else {
        eprintln!("skipped: no pandoc on PATH");
        return;
    };
    assert!(version.status.success());
    let dir = scratch("pandoc");
    let mut checked = 0;
    for (_, docx, doc) in write_reader_files(&dir) {
        let markdown = Command::new("pandoc")
            .args(["-f", "docx", "-t", "markdown"])
            .arg(&docx)
            .output()
            .expect("run pandoc");
        assert!(
            markdown.status.success(),
            "{}",
            String::from_utf8_lossy(&markdown.stderr)
        );
        let text = String::from_utf8_lossy(&markdown.stdout);
        assert!(text.starts_with("# "), "{text}");
        let json = Command::new("pandoc")
            .args(["-f", "docx", "-t", "json"])
            .arg(&docx)
            .output()
            .expect("run pandoc");
        let ast = Json::parse(&String::from_utf8_lossy(&json.stdout)).expect("pandoc JSON");
        let blocks = ast.get("blocks").and_then(Json::as_array).expect("blocks");
        let kind = |b: &Json| b.get("t").and_then(Json::as_str).map(str::to_owned);
        let count = |t: &str| {
            blocks
                .iter()
                .filter(|b| kind(b).as_deref() == Some(t))
                .count()
        };
        assert_eq!(
            count("Table"),
            expected_docx(&doc).len(),
            "{}",
            docx.display()
        );
        assert_eq!(
            count("Header"),
            1 + tables(&doc).len(),
            "{}",
            docx.display()
        );
        checked += 1;
    }
    let _ = fs::remove_dir_all(&dir);
    eprintln!("pandoc read {checked} DOCX files");
}

/// The rows DuckDB's `read_xlsx` gives for `table` with its default
/// types, or `None` when it may fail: a NaN or infinity is the text `?` in a
/// number column. It stops at the first row with no cell, a row whose
/// every cell is a null number.
fn duckdb_default_rows(table: &Table) -> Option<usize> {
    let non_finite = table.rows.iter().any(|row| {
        row.iter()
            .any(|v| matches!(v, Value::Float(f) if !f.is_finite()))
    });
    if non_finite {
        return None;
    }
    let empty = |row: &Vec<Value>| {
        row.iter()
            .zip(&table.columns)
            .all(|(v, (_, ty))| *ty != Ty::Str && matches!(v, Value::Null))
    };
    Some(
        table
            .rows
            .iter()
            .position(empty)
            .unwrap_or(table.rows.len()),
    )
}

/// DuckDB's `excel` extension reads each class sheet with its header row.
/// With its default types: the rows before the first row with no cell,
/// and the sum of an int column equal to the sum of the document's cells.
/// A sheet with a `?` in a number column may fail there; it reads with
/// `all_varchar=true`. Skipped unless `LOAD excel` works offline.
#[test]
fn duckdb_reads_every_class_sheet() {
    let duckdb = |sql: &str| {
        Command::new("duckdb")
            .env_remove("AI_AGENT")
            .args(["-csv", "-noheader", "-c", sql])
            .output()
    };
    if !duckdb("LOAD excel;").is_ok_and(|o| o.status.success()) {
        eprintln!("skipped: duckdb or its excel extension is not installed");
        return;
    }
    let dir = scratch("duckdb");
    let (mut typed, mut varchar) = (0, 0);
    for (xlsx, _, doc) in write_reader_files(&dir) {
        let bytes = fs::read(&xlsx).expect("read back");
        let entries = unzip(&bytes).expect("unzips");
        let workbook = String::from_utf8(entries[2].1.clone()).expect("UTF-8");
        let names = sheet_names(&workbook);
        for (table, sheet) in tables(&doc).into_iter().zip(&names[1..]) {
            if table.rows.is_empty()
                || table.names().collect::<HashSet<_>>().len() != table.columns.len()
            {
                continue;
            }
            let source = format!(
                "read_xlsx('{}', sheet='{}'",
                xlsx.display(),
                sheet.replace('\'', "''")
            );
            let Some(rows) = duckdb_default_rows(table) else {
                let sql = format!("LOAD excel; SELECT count(*) FROM {source}, all_varchar=true);");
                let output = duckdb(&sql).expect("run duckdb");
                let line = String::from_utf8_lossy(&output.stdout).trim().to_owned();
                assert_eq!(line.parse::<usize>(), Ok(table.rows.len()), "{sql}");
                varchar += 1;
                continue;
            };
            // An int column whose every cell a double holds exactly.
            let exact = |j: usize| {
                table.rows[..rows]
                    .iter()
                    .all(|r| r[j].as_int().is_none_or(|v| v.unsigned_abs() < 1 << 53))
            };
            let number =
                (0..table.columns.len()).find(|&j| table.columns[j].1 == Ty::Int && exact(j));
            let sum = number.map_or_else(
                || "0".to_owned(),
                |j| format!("coalesce(sum(\"{}\"), 0)", table.columns[j].0),
            );
            let sql = format!("LOAD excel; SELECT count(*), {sum} FROM {source});");
            let output = duckdb(&sql).expect("run duckdb");
            assert!(
                output.status.success(),
                "{sql}: {}",
                String::from_utf8_lossy(&output.stderr)
            );
            let line = String::from_utf8_lossy(&output.stdout).trim().to_owned();
            let (count, total) = line.split_once(',').expect("two columns");
            assert_eq!(count.parse::<usize>(), Ok(rows), "{sql}");
            if let Some(j) = number {
                let want: i64 = table.rows[..rows]
                    .iter()
                    .filter_map(|r| r[j].as_int())
                    .sum();
                let got: f64 = total.parse().expect("a number");
                assert!(
                    (got - want as f64).abs() <= want.unsigned_abs() as f64 * 1e-12,
                    "{sql}: {got} for {want}"
                );
            }
            typed += 1;
        }
    }
    let _ = fs::remove_dir_all(&dir);
    eprintln!("duckdb read {typed} sheets with its types and {varchar} as text");
    assert!(typed > 0);
}

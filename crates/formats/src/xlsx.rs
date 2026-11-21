//! XLSX: an Office Open XML workbook, hand-written into a stored zip
//! ([`crate::ooxml`]).
//!
//! * The first sheet, `Headers`, holds the `StringHeader` texts, one per
//!   row in column A.
//! * Then one sheet per metrics class, named by the class after its last
//!   `.` or `$`, and a sheet `Histogram` for the histogram. Row 1 holds
//!   the column names in bold and stays in view when the sheet scrolls.
//! * An `int` cell is a number with every digit; a `float` cell is a
//!   number with the shortest digits that read back as the same double,
//!   and always a point or an exponent, so a reader gives a float. NaN and
//!   infinity have no spreadsheet number: they are the text Picard prints,
//!   `?` and `-?`. A null cell and an empty string are no cell.
//! * A string is an inline string. A character XML 1.0 cannot hold is
//!   `_xHHHH_`, the escape of the OOXML string type that Excel decodes,
//!   and a `_x` that would read as one is written `_x005F_x`. A string
//!   is cut at 32,767 UTF-16 units, the most an Excel cell holds.
//! * A sheet name has at most 31 UTF-16 units and none of `: \ / ? * [ ]`;
//!   a name equal to an earlier one in any case gets ` (2)`, ` (3)`.
//!
//! A sheet holds at most 1,048,576 rows and 16,384 columns, as Excel
//! does; a larger table is an error, as is a workbook over the zip cap.
//! The bytes depend on the document only: no date or name of the run.

use std::collections::HashSet;

use bammetrics_core::metric::Value;

use crate::error::{Error, Result};
use crate::ooxml::{DECLARATION, MAX_ZIP_BYTES, Zip, escape_into, package_rels};
use crate::table::{Document, Table};

/// The most rows of a sheet.
pub const MAX_ROWS: usize = 1 << 20;
/// The most columns of a sheet.
pub const MAX_COLUMNS: usize = 1 << 14;
/// The most UTF-16 units of a cell.
pub const MAX_CELL_UNITS: usize = 32_767;
/// The most UTF-16 units of a sheet name.
pub const MAX_NAME_UNITS: usize = 31;

const MAIN: &str = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
const RELATIONSHIPS: &str = "http://schemas.openxmlformats.org/officeDocument/2006/relationships";

const STYLES: &str = "<styleSheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">\
<fonts count=\"2\"><font><sz val=\"11\"/><name val=\"Calibri\"/></font>\
<font><b/><sz val=\"11\"/><name val=\"Calibri\"/></font></fonts>\
<fills count=\"2\"><fill><patternFill patternType=\"none\"/></fill>\
<fill><patternFill patternType=\"gray125\"/></fill></fills>\
<borders count=\"1\"><border><left/><right/><top/><bottom/><diagonal/></border></borders>\
<cellStyleXfs count=\"1\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\"/></cellStyleXfs>\
<cellXfs count=\"2\"><xf numFmtId=\"0\" fontId=\"0\" fillId=\"0\" borderId=\"0\" xfId=\"0\"/>\
<xf numFmtId=\"0\" fontId=\"1\" fillId=\"0\" borderId=\"0\" xfId=\"0\" applyFont=\"1\"/></cellXfs>\
<cellStyles count=\"1\"><cellStyle name=\"Normal\" xfId=\"0\" builtinId=\"0\"/></cellStyles>\
</styleSheet>\n";

/// The letters of column `j` (0 is `A`, 26 is `AA`).
fn column_letters(j: usize) -> String {
    let mut letters = Vec::new();
    let mut n = j + 1;
    while n > 0 {
        n -= 1;
        letters.push(b'A' + (n % 26) as u8);
        n /= 26;
    }
    letters.reverse();
    String::from_utf8(letters).expect("ASCII letters")
}

/// `value` as the text of a float cell: the shortest digits that read
/// back as `value`, with a point or an exponent.
fn float_text(value: f64) -> String {
    let magnitude = value.abs();
    let text = if magnitude != 0.0 && !(1e-5..1e15).contains(&magnitude) {
        format!("{value:e}")
    } else {
        value.to_string()
    };
    if text.contains(['.', 'e']) {
        text
    } else {
        text + ".0"
    }
}

/// `text` cut to at most [`MAX_CELL_UNITS`] UTF-16 units.
fn cut(text: &str) -> &str {
    let mut units = 0;
    for (i, c) in text.char_indices() {
        units += c.len_utf16();
        if units > MAX_CELL_UNITS {
            return &text[..i];
        }
    }
    text
}

/// `text` as the content of an OOXML string: XML escapes, `_xHHHH_` for
/// a character XML cannot hold, and `_x005F_` before a `_` that would
/// start such an escape.
fn string_text(out: &mut String, text: &str) {
    let bytes = text.as_bytes();
    let mut start = 0;
    for (i, _) in text.match_indices("_x") {
        let escape_like = bytes.len() >= i + 7
            && bytes[i + 2..i + 6].iter().all(u8::is_ascii_hexdigit)
            && bytes[i + 6] == b'_';
        if escape_like {
            escape_into(out, &text[start..i], hex_escape);
            out.push_str("_x005F");
            start = i;
        }
    }
    escape_into(out, &text[start..], hex_escape);
}

fn hex_escape(out: &mut String, c: char) {
    out.push_str(&format!("_x{:04X}_", u32::from(c)));
}

/// A sheet's XML being written row by row.
struct Sheet {
    xml: String,
    rows: usize,
}

impl Sheet {
    fn new(rows: usize, columns: usize, frozen_header: bool) -> Self {
        let last = if rows == 0 || columns == 0 {
            "A1".to_owned()
        } else {
            format!("{}{rows}", column_letters(columns - 1))
        };
        let view = if frozen_header {
            "<sheetViews><sheetView workbookViewId=\"0\"><pane ySplit=\"1\" topLeftCell=\"A2\" \
             activePane=\"bottomLeft\" state=\"frozen\"/></sheetView></sheetViews>"
        } else {
            "<sheetViews><sheetView workbookViewId=\"0\"/></sheetViews>"
        };
        let range = if last == "A1" {
            last
        } else {
            format!("A1:{last}")
        };
        Self {
            xml: format!(
                "{DECLARATION}<worksheet xmlns=\"{MAIN}\"><dimension ref=\"{range}\"/>{view}<sheetData>\n"
            ),
            rows: 0,
        }
    }

    /// Add one row of cells; a `None` cell is left out.
    fn row<'a>(&mut self, cells: impl Iterator<Item = Option<Cell<'a>>>) -> Result<()> {
        self.rows += 1;
        let r = self.rows;
        self.xml.push_str(&format!("<row r=\"{r}\">"));
        for (j, cell) in cells.enumerate() {
            let Some(cell) = cell else { continue };
            let at = format!("{}{r}", column_letters(j));
            match cell {
                Cell::Number(text) => {
                    self.xml
                        .push_str(&format!("<c r=\"{at}\"><v>{text}</v></c>"));
                }
                Cell::Text(text, bold) => {
                    let style = if bold { " s=\"1\"" } else { "" };
                    let text = cut(text);
                    let space = if text.starts_with(char::is_whitespace)
                        || text.ends_with(char::is_whitespace)
                    {
                        " xml:space=\"preserve\""
                    } else {
                        ""
                    };
                    self.xml.push_str(&format!(
                        "<c r=\"{at}\"{style} t=\"inlineStr\"><is><t{space}>"
                    ));
                    string_text(&mut self.xml, text);
                    self.xml.push_str("</t></is></c>");
                }
            }
        }
        self.xml.push_str("</row>\n");
        if self.xml.len() > MAX_ZIP_BYTES {
            return Err(Error::Limit(format!(
                "a sheet is over the {MAX_ZIP_BYTES}-byte cap of an XLSX file"
            )));
        }
        Ok(())
    }

    fn finish(mut self) -> String {
        self.xml.push_str("</sheetData></worksheet>\n");
        self.xml
    }
}

/// One cell: a number's text, or a string and whether it is bold.
enum Cell<'a> {
    Number(String),
    Text(&'a str, bool),
}

fn value_cell(value: &Value) -> Option<Cell<'_>> {
    match value {
        Value::Null => None,
        Value::Int(v) => Some(Cell::Number(v.to_string())),
        Value::Float(v) if v.is_finite() => Some(Cell::Number(float_text(*v))),
        Value::Float(_) | Value::Bool(_) => Some(Cell::Text(value_text(value), false)),
        Value::Str(s) if s.is_empty() => None,
        Value::Str(s) => Some(Cell::Text(s, false)),
    }
}

/// The text of a cell with no spreadsheet number: NaN and infinity as
/// Picard prints them, a boolean as `Y` or `N`.
fn value_text(value: &Value) -> &'static str {
    match value {
        Value::Float(v) if v.is_nan() => "?",
        Value::Float(v) if *v < 0.0 => "-?",
        Value::Float(_) => "?",
        Value::Bool(true) => "Y",
        _ => "N",
    }
}

/// The sheet of `table`: the column names, then the rows.
fn table_sheet(table: &Table) -> Result<String> {
    let rows = table.rows.len() + 1;
    let columns = table.columns.len();
    if rows > MAX_ROWS || columns > MAX_COLUMNS {
        return Err(Error::Limit(format!(
            "an XLSX sheet holds at most {MAX_ROWS} rows and {MAX_COLUMNS} columns; \
             {} has {rows} rows and {columns} columns",
            if table.class.is_empty() {
                "the table"
            } else {
                &table.class
            }
        )));
    }
    let mut sheet = Sheet::new(rows, columns, true);
    sheet.row(table.names().map(|name| Some(Cell::Text(name, true))))?;
    for row in &table.rows {
        sheet.row(row.iter().map(value_cell))?;
    }
    Ok(sheet.finish())
}

/// The sheet of the `StringHeader` texts, one per row.
fn headers_sheet(headers: &[String]) -> Result<String> {
    if headers.len() > MAX_ROWS {
        return Err(Error::Limit(format!(
            "an XLSX sheet holds at most {MAX_ROWS} rows; the file has {} header lines",
            headers.len()
        )));
    }
    let mut sheet = Sheet::new(headers.len(), 1, false);
    for header in headers {
        sheet.row(std::iter::once(Some(Cell::Text(header, false))))?;
    }
    Ok(sheet.finish())
}

/// Sheet names as Excel takes them: forbidden characters as `_`, at most
/// [`MAX_NAME_UNITS`] UTF-16 units, no two equal in any case.
struct Names {
    taken: HashSet<String>,
}

impl Names {
    fn new() -> Self {
        // Excel keeps `History` for itself.
        Self {
            taken: HashSet::from(["history".to_owned()]),
        }
    }

    fn fit(base: &str, suffix: &str) -> String {
        let room = MAX_NAME_UNITS - suffix.len();
        let mut name = String::new();
        let mut units = 0;
        for c in base.chars() {
            units += c.len_utf16();
            if units > room {
                break;
            }
            name.push(c);
        }
        name.push_str(suffix);
        name
    }

    fn name(&mut self, wanted: &str) -> String {
        let clean: String = wanted
            .chars()
            .map(|c| {
                if matches!(c, ':' | '\\' | '/' | '?' | '*' | '[' | ']')
                    || !crate::ooxml::xml_char(c)
                    || c.is_control()
                {
                    '_'
                } else {
                    c
                }
            })
            .collect();
        let clean = clean.trim_matches('\'');
        let base = if clean.trim().is_empty() {
            "Sheet"
        } else {
            clean
        };
        let mut n = 1;
        loop {
            let suffix = if n == 1 {
                String::new()
            } else {
                format!(" ({n})")
            };
            let name = Self::fit(base, &suffix);
            if self.taken.insert(name.to_lowercase()) {
                return name;
            }
            n += 1;
        }
    }
}

/// The name of a class's sheet: the class after its last `.` or `$`.
fn class_sheet_name(class: &str) -> &str {
    match class.rsplit(['.', '$']).next() {
        Some(short) if !short.is_empty() => short,
        _ => "Metrics",
    }
}

/// `doc` as an XLSX workbook.
pub fn write(doc: &Document) -> Result<Vec<u8>> {
    let mut names = Names::new();
    let mut sheets = vec![(names.name("Headers"), headers_sheet(&doc.headers)?)];
    for table in &doc.tables {
        let name = names.name(class_sheet_name(&table.class));
        sheets.push((name, table_sheet(table)?));
    }
    if let Some(histogram) = &doc.histogram {
        sheets.push((names.name("Histogram"), table_sheet(histogram)?));
    }

    let mut types = format!(
        "{DECLARATION}<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\
<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\
<Default Extension=\"xml\" ContentType=\"application/xml\"/>\
<Override PartName=\"/xl/workbook.xml\" \
ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>\
<Override PartName=\"/xl/styles.xml\" \
ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml\"/>"
    );
    let mut workbook = format!(
        "{DECLARATION}<workbook xmlns=\"{MAIN}\" xmlns:r=\"{RELATIONSHIPS}\">\
<bookViews><workbookView/></bookViews><sheets>"
    );
    let mut rels = format!(
        "{DECLARATION}<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
    );
    for (i, (name, _)) in sheets.iter().enumerate() {
        let n = i + 1;
        types.push_str(&format!(
            "<Override PartName=\"/xl/worksheets/sheet{n}.xml\" \
ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
        ));
        workbook.push_str("<sheet name=\"");
        escape_into(&mut workbook, name, |_, _| {});
        workbook.push_str(&format!("\" sheetId=\"{n}\" r:id=\"rId{n}\"/>"));
        rels.push_str(&format!(
            "<Relationship Id=\"rId{n}\" Type=\"{RELATIONSHIPS}/worksheet\" \
Target=\"worksheets/sheet{n}.xml\"/>"
        ));
    }
    types.push_str("</Types>\n");
    workbook.push_str("</sheets></workbook>\n");
    rels.push_str(&format!(
        "<Relationship Id=\"rId{}\" Type=\"{RELATIONSHIPS}/styles\" Target=\"styles.xml\"/>\
</Relationships>\n",
        sheets.len() + 1
    ));

    let mut zip = Zip::new();
    zip.add("[Content_Types].xml", types.as_bytes())?;
    zip.add("_rels/.rels", package_rels("xl/workbook.xml").as_bytes())?;
    zip.add("xl/workbook.xml", workbook.as_bytes())?;
    zip.add("xl/_rels/workbook.xml.rels", rels.as_bytes())?;
    zip.add("xl/styles.xml", format!("{DECLARATION}{STYLES}").as_bytes())?;
    for (i, (_, sheet)) in sheets.iter().enumerate() {
        zip.add(
            &format!("xl/worksheets/sheet{}.xml", i + 1),
            sheet.as_bytes(),
        )?;
    }
    Ok(zip.finish())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::table::Ty;

    #[test]
    fn columns_are_lettered_as_excel_letters_them() {
        assert_eq!(column_letters(0), "A");
        assert_eq!(column_letters(25), "Z");
        assert_eq!(column_letters(26), "AA");
        assert_eq!(column_letters(701), "ZZ");
        assert_eq!(column_letters(702), "AAA");
        assert_eq!(column_letters(MAX_COLUMNS - 1), "XFD");
    }

    #[test]
    fn floats_keep_every_digit_and_read_as_floats() {
        assert_eq!(float_text(310.0), "310.0");
        assert_eq!(float_text(8.164_965_809_277_26), "8.16496580927726");
        assert_eq!(float_text(-0.5), "-0.5");
        assert_eq!(float_text(0.0), "0.0");
        assert_eq!(float_text(1e300), "1e300");
        assert_eq!(float_text(1.5e-7), "1.5e-7");
        assert_eq!(float_text(f64::MAX).parse::<f64>(), Ok(f64::MAX));
        assert_eq!(float_text(5e-324).parse::<f64>(), Ok(5e-324));
    }

    #[test]
    fn strings_escape_what_xml_cannot_hold() {
        let mut out = String::new();
        string_text(&mut out, "a<b\u{1}_x0041_ _x12_ _xZZZZ_\r");
        assert_eq!(out, "a&lt;b_x0001__x005F_x0041_ _x12_ _xZZZZ_&#13;");
        let long = "é".repeat(MAX_CELL_UNITS + 5);
        assert_eq!(cut(&long).chars().count(), MAX_CELL_UNITS);
        let wide = "😀".repeat(MAX_CELL_UNITS);
        assert_eq!(cut(&wide).encode_utf16().count(), MAX_CELL_UNITS - 1);
    }

    #[test]
    fn sheet_names_fit_and_differ() {
        let mut names = Names::new();
        assert_eq!(names.name("Headers"), "Headers");
        assert_eq!(names.name("headers"), "headers (2)");
        assert_eq!(names.name("History"), "History (2)");
        assert_eq!(names.name("a/b:c?[d]*"), "a_b_c__d__");
        assert_eq!(names.name("'quoted'"), "quoted");
        assert_eq!(names.name("''"), "Sheet");
        let long = "X".repeat(40);
        assert_eq!(names.name(&long), "X".repeat(31));
        assert_eq!(names.name(&long), format!("{} (2)", "X".repeat(27)));
        assert_eq!(
            class_sheet_name("picard.analysis.InsertSizeMetrics"),
            "InsertSizeMetrics"
        );
        assert_eq!(class_sheet_name("a.B$Inner"), "Inner");
        assert_eq!(class_sheet_name(""), "Metrics");
    }

    fn doc() -> Document {
        Document {
            tool: Some("CollectX".to_owned()),
            sample: None,
            headers: vec![
                "CollectX INPUT=x.bam".to_owned(),
                "Started on: now".to_owned(),
            ],
            tables: vec![Table {
                class: "picard.analysis.XMetrics".to_owned(),
                columns: vec![
                    ("N".to_owned(), Ty::Int),
                    ("F".to_owned(), Ty::Float),
                    ("S".to_owned(), Ty::Str),
                ],
                rows: vec![
                    vec![Value::Int(i64::MAX), Value::Float(f64::NAN), Value::Null],
                    vec![
                        Value::Int(-1),
                        Value::Float(f64::NEG_INFINITY),
                        Value::Str(" a&b ".to_owned()),
                    ],
                ],
            }],
            histogram: Some(Table {
                class: "java.lang.Integer".to_owned(),
                columns: vec![("insert_size".to_owned(), Ty::Int)],
                rows: Vec::new(),
            }),
            verbatim: None,
        }
    }

    #[test]
    fn a_sheet_holds_typed_cells() {
        let sheet = table_sheet(&doc().tables[0]).expect("fits");
        assert!(sheet.contains("<dimension ref=\"A1:C3\"/>"), "{sheet}");
        assert!(
            sheet
                .contains("<row r=\"1\"><c r=\"A1\" s=\"1\" t=\"inlineStr\"><is><t>N</t></is></c>")
        );
        assert!(sheet.contains(
            "<row r=\"2\"><c r=\"A2\"><v>9223372036854775807</v></c>\
             <c r=\"B2\" t=\"inlineStr\"><is><t>?</t></is></c></row>"
        ));
        assert!(sheet.contains("<c r=\"B3\" t=\"inlineStr\"><is><t>-?</t></is></c>"));
        assert!(sheet.contains("<t xml:space=\"preserve\"> a&amp;b </t>"));
    }

    #[test]
    fn a_workbook_is_the_same_bytes_each_time() {
        let bytes = write(&doc()).expect("writes");
        assert_eq!(bytes, write(&doc()).expect("writes"));
        assert_eq!(&bytes[..4], b"PK\x03\x04");
        let text = String::from_utf8_lossy(&bytes);
        for name in ["Headers", "XMetrics", "Histogram"] {
            assert!(text.contains(&format!("<sheet name=\"{name}\"")), "{name}");
        }
        assert!(text.contains("xl/worksheets/sheet3.xml"));
    }

    #[test]
    fn a_table_past_the_sheet_limits_is_an_error() {
        let table = Table {
            class: "Wide".to_owned(),
            columns: (0..=MAX_COLUMNS)
                .map(|j| (j.to_string(), Ty::Int))
                .collect(),
            rows: Vec::new(),
        };
        let error = table_sheet(&table).expect_err("too wide").to_string();
        assert!(
            error.contains("16384 columns; Wide has 1 rows and 16385"),
            "{error}"
        );
    }
}

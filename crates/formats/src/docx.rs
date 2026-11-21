//! DOCX: a WordprocessingML document of text and tables, hand-written
//! into a stored zip ([`crate::ooxml`]). No plot.
//!
//! * A `Heading 1` with the tool (or `Metrics`) and the sample when known,
//!   then one paragraph per `StringHeader` text.
//! * Per metrics class a `Heading 2` with the class, then a table whose
//!   first row holds the column names in bold and repeats on each page. A
//!   class with no row is the line `No rows.`
//! * The histogram under the Markdown heading `Histogram of CLASS, by KEY`.
//! * A table with more than 8 columns and at most 8 rows is turned on its
//!   side, as the Markdown form turns it: a row per column, the column
//!   name first, then one cell per row.
//! * Every cell holds the text Picard prints (at most 6 fraction digits,
//!   `?` for NaN); numbers are right-aligned. A tab is a Word tab, a line
//!   end a line break, and a character XML 1.0 cannot hold is U+FFFD.
//! * The page is A4 landscape.
//!
//! The bytes depend on the document only: no date or name of the run.

use bammetrics_core::metric::Value;

use crate::error::Result;
use crate::markdown::{histogram_heading, sideways};
use crate::ooxml::{DECLARATION, Zip, escape_into, package_rels};
use crate::table::{Document, Table, Ty};

const W: &str = "http://schemas.openxmlformats.org/wordprocessingml/2006/main";

const CONTENT_TYPES: &str = "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\
<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\
<Default Extension=\"xml\" ContentType=\"application/xml\"/>\
<Override PartName=\"/word/document.xml\" \
ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml\"/>\
<Override PartName=\"/word/styles.xml\" \
ContentType=\"application/vnd.openxmlformats-officedocument.wordprocessingml.styles+xml\"/>\
</Types>\n";

const DOCUMENT_RELS: &str = "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\
<Relationship Id=\"rId1\" \
Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles\" \
Target=\"styles.xml\"/></Relationships>\n";

/// Normal text at 9 points, two headings, and a table style with a grid.
const STYLES: &str = "<w:styles xmlns:w=\"http://schemas.openxmlformats.org/wordprocessingml/2006/main\">\
<w:docDefaults><w:rPrDefault><w:rPr><w:rFonts w:ascii=\"Calibri\" w:hAnsi=\"Calibri\" w:cs=\"Calibri\"/>\
<w:sz w:val=\"18\"/></w:rPr></w:rPrDefault>\
<w:pPrDefault><w:pPr><w:spacing w:after=\"60\"/></w:pPr></w:pPrDefault></w:docDefaults>\
<w:style w:type=\"paragraph\" w:default=\"1\" w:styleId=\"Normal\"><w:name w:val=\"Normal\"/></w:style>\
<w:style w:type=\"paragraph\" w:styleId=\"Heading1\"><w:name w:val=\"heading 1\"/>\
<w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:pPr><w:keepNext/><w:spacing w:before=\"240\" \
w:after=\"120\"/><w:outlineLvl w:val=\"0\"/></w:pPr><w:rPr><w:b/><w:sz w:val=\"32\"/></w:rPr></w:style>\
<w:style w:type=\"paragraph\" w:styleId=\"Heading2\"><w:name w:val=\"heading 2\"/>\
<w:basedOn w:val=\"Normal\"/><w:next w:val=\"Normal\"/><w:pPr><w:keepNext/><w:spacing w:before=\"240\" \
w:after=\"120\"/><w:outlineLvl w:val=\"1\"/></w:pPr><w:rPr><w:b/><w:sz w:val=\"26\"/></w:rPr></w:style>\
<w:style w:type=\"table\" w:default=\"1\" w:styleId=\"TableNormal\"><w:name w:val=\"Normal Table\"/>\
<w:tblPr><w:tblCellMar><w:left w:w=\"80\" w:type=\"dxa\"/><w:right w:w=\"80\" w:type=\"dxa\"/>\
</w:tblCellMar></w:tblPr></w:style>\
<w:style w:type=\"table\" w:styleId=\"TableGrid\"><w:name w:val=\"Table Grid\"/>\
<w:basedOn w:val=\"TableNormal\"/><w:pPr><w:spacing w:after=\"0\"/></w:pPr><w:tblPr><w:tblBorders>\
<w:top w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
<w:left w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
<w:bottom w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
<w:right w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
<w:insideH w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
<w:insideV w:val=\"single\" w:sz=\"4\" w:space=\"0\" w:color=\"808080\"/>\
</w:tblBorders></w:tblPr></w:style></w:styles>\n";

/// A4 landscape with 1.5 cm margins, in twentieths of a point.
const SECTION: &str = "<w:sectPr><w:pgSz w:w=\"16838\" w:h=\"11906\" w:orient=\"landscape\"/>\
<w:pgMar w:top=\"851\" w:right=\"851\" w:bottom=\"851\" w:left=\"851\" w:header=\"0\" w:footer=\"0\" \
w:gutter=\"0\"/></w:sectPr>";

/// `text` as runs of one paragraph: tabs as `w:tab`, a line end (`\n`,
/// `\r\n` or `\r`) as `w:br`, a character XML cannot hold as U+FFFD.
fn runs(out: &mut String, text: &str, bold: bool) {
    out.push_str("<w:r>");
    if bold {
        out.push_str("<w:rPr><w:b/></w:rPr>");
    }
    let text = text.replace("\r\n", "\n");
    let mut pieces = text.split(['\t', '\n', '\r']);
    let mut breaks = text.matches(['\t', '\n', '\r']);
    if let Some(first) = pieces.next() {
        piece(out, first);
    }
    for next in pieces {
        out.push_str(match breaks.next() {
            Some("\t") => "<w:tab/>",
            _ => "<w:br/>",
        });
        piece(out, next);
    }
    out.push_str("</w:r>");
}

fn piece(out: &mut String, text: &str) {
    if text.is_empty() {
        return;
    }
    let space = if text.starts_with(' ') || text.ends_with(' ') {
        " xml:space=\"preserve\""
    } else {
        ""
    };
    out.push_str(&format!("<w:t{space}>"));
    escape_into(out, text, |out, _| out.push('\u{fffd}'));
    out.push_str("</w:t>");
}

/// A paragraph of `text` in `style` (`None` for Normal).
fn paragraph(out: &mut String, style: Option<&str>, text: &str) {
    out.push_str("<w:p>");
    if let Some(style) = style {
        out.push_str(&format!("<w:pPr><w:pStyle w:val=\"{style}\"/></w:pPr>"));
    }
    runs(out, text, false);
    out.push_str("</w:p>\n");
}

/// One table cell: its text, bold, right-aligned.
fn cell(out: &mut String, text: &str, bold: bool, right: bool) {
    out.push_str("<w:tc><w:p>");
    if right {
        out.push_str("<w:pPr><w:jc w:val=\"right\"/></w:pPr>");
    }
    if !text.is_empty() {
        runs(out, text, bold);
    }
    out.push_str("</w:p></w:tc>");
}

/// A table of `rows`, the first the header row. `right(i, j)` tells a cell
/// to align right.
fn table_xml(
    out: &mut String,
    columns: usize,
    rows: impl Iterator<Item = Vec<String>>,
    right: impl Fn(usize, usize) -> bool,
) {
    out.push_str(
        "<w:tbl><w:tblPr><w:tblStyle w:val=\"TableGrid\"/><w:tblW w:w=\"0\" w:type=\"auto\"/>\
<w:tblLook w:val=\"04A0\" w:firstRow=\"1\" w:lastRow=\"0\" w:firstColumn=\"0\" w:lastColumn=\"0\" \
w:noHBand=\"0\" w:noVBand=\"1\"/></w:tblPr><w:tblGrid>",
    );
    for _ in 0..columns {
        out.push_str("<w:gridCol/>");
    }
    out.push_str("</w:tblGrid>\n");
    for (i, row) in rows.enumerate() {
        out.push_str("<w:tr>");
        if i == 0 {
            out.push_str("<w:trPr><w:tblHeader/></w:trPr>");
        }
        for (j, text) in row.iter().enumerate() {
            cell(out, text, i == 0, i > 0 && right(i - 1, j));
        }
        out.push_str("</w:tr>\n");
    }
    out.push_str("</w:tbl>\n");
}

fn numeric(ty: Ty) -> bool {
    ty != Ty::Str
}

/// `table` under `heading`, turned on its side when it is wide and short.
fn table_section(out: &mut String, heading: &str, table: &Table) {
    paragraph(out, Some("Heading2"), heading);
    if table.rows.is_empty() {
        paragraph(out, None, "No rows.");
        return;
    }
    let text = |value: &Value| value.picard_text();
    if sideways(table) {
        let header = std::iter::once("Column".to_owned())
            .chain((1..=table.rows.len()).map(|i| i.to_string()))
            .collect();
        let lines = table.names().enumerate().map(|(j, name)| {
            std::iter::once(name.to_owned())
                .chain(table.rows.iter().map(|row| text(&row[j])))
                .collect()
        });
        let right = |j: usize, k: usize| k > 0 && numeric(table.columns[j].1);
        table_xml(
            out,
            table.rows.len() + 1,
            std::iter::once(header).chain(lines),
            right,
        );
    } else {
        let header = table.names().map(str::to_owned).collect();
        let lines = table.rows.iter().map(|row| row.iter().map(text).collect());
        let right = |_: usize, j: usize| numeric(table.columns[j].1);
        table_xml(
            out,
            table.columns.len(),
            std::iter::once(header).chain(lines),
            right,
        );
    }
}

/// `doc` as a DOCX document.
pub fn write(doc: &Document) -> Result<Vec<u8>> {
    let mut body = String::new();
    let title = match (&doc.tool, &doc.sample) {
        (Some(tool), Some(sample)) => format!("{tool}, sample {sample}"),
        (Some(tool), None) => tool.clone(),
        (None, Some(sample)) => format!("Metrics, sample {sample}"),
        (None, None) => "Metrics".to_owned(),
    };
    paragraph(&mut body, Some("Heading1"), &title);
    for header in &doc.headers {
        paragraph(&mut body, None, header);
    }
    for table in &doc.tables {
        let heading = if table.class.is_empty() {
            "Rows"
        } else {
            &table.class
        };
        table_section(&mut body, heading, table);
    }
    if let Some(histogram) = &doc.histogram {
        table_section(&mut body, &histogram_heading(doc, histogram), histogram);
    }
    let document = format!(
        "{DECLARATION}<w:document xmlns:w=\"{W}\"><w:body>\n{body}{SECTION}</w:body></w:document>\n"
    );
    let mut zip = Zip::new();
    zip.add(
        "[Content_Types].xml",
        format!("{DECLARATION}{CONTENT_TYPES}").as_bytes(),
    )?;
    zip.add("_rels/.rels", package_rels("word/document.xml").as_bytes())?;
    zip.add("word/document.xml", document.as_bytes())?;
    zip.add(
        "word/_rels/document.xml.rels",
        format!("{DECLARATION}{DOCUMENT_RELS}").as_bytes(),
    )?;
    zip.add(
        "word/styles.xml",
        format!("{DECLARATION}{STYLES}").as_bytes(),
    )?;
    Ok(zip.finish())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runs_turn_tabs_and_line_ends_into_word_marks() {
        let mut out = String::new();
        runs(&mut out, "a\tb\r\nc\rd\u{1} e ", false);
        assert_eq!(
            out,
            "<w:r><w:t>a</w:t><w:tab/><w:t>b</w:t><w:br/><w:t>c</w:t><w:br/>\
             <w:t xml:space=\"preserve\">d\u{fffd} e </w:t></w:r>"
        );
    }

    fn table(columns: usize, rows: usize) -> Table {
        Table {
            class: "picard.analysis.Demo".to_owned(),
            columns: (0..columns)
                .map(|j| (format!("C{j}"), if j == 0 { Ty::Str } else { Ty::Float }))
                .collect(),
            rows: (0..rows)
                .map(|i| {
                    (0..columns)
                        .map(|j| match j {
                            0 => Value::Str(format!("r{i}<")),
                            _ => Value::Float(i as f64 / 3.0),
                        })
                        .collect()
                })
                .collect(),
        }
    }

    #[test]
    fn narrow_tables_keep_their_rows_and_wide_short_ones_turn() {
        let mut out = String::new();
        table_section(&mut out, "Narrow", &table(2, 2));
        assert_eq!(out.matches("<w:tr>").count(), 3);
        assert!(out.contains("<w:t>r1&lt;</w:t>"));
        assert!(out.contains("<w:jc w:val=\"right\"/></w:pPr><w:r><w:t>0.333333</w:t>"));
        let mut out = String::new();
        table_section(&mut out, "Wide", &table(9, 2));
        assert_eq!(out.matches("<w:tr>").count(), 10);
        assert_eq!(out.matches("<w:gridCol/>").count(), 3);
        let mut out = String::new();
        table_section(&mut out, "Empty", &table(2, 0));
        assert!(out.contains("No rows."));
        assert!(!out.contains("<w:tbl>"));
    }

    #[test]
    fn a_document_is_the_same_bytes_each_time() {
        let doc = Document {
            tool: Some("CollectX".to_owned()),
            sample: Some("S&1".to_owned()),
            headers: vec!["CollectX INPUT=x.bam".to_owned()],
            tables: vec![table(3, 1)],
            histogram: Some(table(2, 3)),
            verbatim: None,
        };
        let bytes = write(&doc).expect("writes");
        assert_eq!(bytes, write(&doc).expect("writes"));
        let text = String::from_utf8_lossy(&bytes);
        assert!(text.contains("<w:t>CollectX, sample S&amp;1</w:t>"));
        assert!(text.contains("<w:t>Histogram of picard.analysis.Demo, by C0</w:t>"));
    }
}

//! A sample sheet read back: CSV (RFC 4180) or TSV, a header row naming
//! the columns, one row per sample. `doctor --sample-sheet` reads it.
//!
//! * The form is told by the first line: a tab in it makes the file TSV,
//!   else CSV.
//! * CSV: a field in double quotes may hold a comma, a line end or a
//!   doubled quote. `\r\n` and `\n` end a row.
//! * TSV: fields split at tabs; `\t`, `\n`, `\r` and `\\` are read back as
//!   [`crate::delimited`] writes them.
//! * A byte order mark is dropped, blank lines are skipped, and a cell is
//!   trimmed of spaces.
//! * Column names are matched without case. A row with fewer cells than
//!   the header has empty cells for the rest; a row with more is an error.

use crate::error::{Error, Result};

/// The largest sheet read: 16 MiB.
pub const MAX_SHEET_BYTES: usize = 16 << 20;

/// The cells of a sheet.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sheet {
    /// The header row, as written.
    pub columns: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Sheet {
    /// The index of the column named `name`, any case.
    pub fn column(&self, name: &str) -> Option<usize> {
        self.columns
            .iter()
            .position(|c| c.eq_ignore_ascii_case(name))
    }

    /// The cell of `row` in column `name`; `None` without such a column.
    pub fn cell<'a>(&self, row: &'a [String], name: &str) -> Option<&'a str> {
        self.column(name)
            .map(|i| row.get(i).map_or("", String::as_str))
    }
}

fn error(line: usize, what: impl std::fmt::Display) -> Error {
    Error::Document(format!("sample sheet line {line}: {what}"))
}

/// One TSV cell with the escapes of [`crate::delimited::tsv`] undone.
fn tsv_cell(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('t') => out.push('\t'),
            Some('n') => out.push('\n'),
            Some('r') => out.push('\r'),
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

/// The rows of CSV text, each with the line it starts on.
fn csv_rows(text: &str) -> Result<Vec<(usize, Vec<String>)>> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut line = 1;
    let mut start = 1;
    let mut chars = text.chars().peekable();
    while let Some(c) = chars.next() {
        if quoted {
            match c {
                '"' if chars.peek() == Some(&'"') => {
                    chars.next();
                    field.push('"');
                }
                '"' => quoted = false,
                '\n' => {
                    line += 1;
                    field.push(c);
                }
                c => field.push(c),
            }
            continue;
        }
        match c {
            '"' if field.trim().is_empty() => {
                field.clear();
                quoted = true;
            }
            '"' => return Err(error(line, "a quote inside an unquoted field")),
            ',' => row.push(std::mem::take(&mut field)),
            '\r' if chars.peek() == Some(&'\n') => {}
            '\n' => {
                row.push(std::mem::take(&mut field));
                rows.push((start, std::mem::take(&mut row)));
                line += 1;
                start = line;
            }
            c => field.push(c),
        }
    }
    if quoted {
        return Err(error(start, "a quoted field is not closed"));
    }
    if !field.is_empty() || !row.is_empty() {
        row.push(field);
        rows.push((start, row));
    }
    Ok(rows)
}

/// The sheet in `text`.
pub fn parse(text: &str) -> Result<Sheet> {
    if text.len() > MAX_SHEET_BYTES {
        return Err(Error::Document(format!(
            "the sample sheet is over the {MAX_SHEET_BYTES}-byte cap"
        )));
    }
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let tab = text.lines().next().is_some_and(|l| l.contains('\t'));
    let raw: Vec<(usize, Vec<String>)> = if tab {
        text.lines()
            .enumerate()
            .map(|(i, l)| {
                let l = l.strip_suffix('\r').unwrap_or(l);
                (i + 1, l.split('\t').map(tsv_cell).collect())
            })
            .collect()
    } else {
        csv_rows(text)?
    };
    let mut rows = raw.into_iter().filter_map(|(line, cells)| {
        let cells: Vec<String> = cells.iter().map(|c| c.trim().to_owned()).collect();
        (!cells.iter().all(String::is_empty)).then_some((line, cells))
    });
    let Some((_, columns)) = rows.next() else {
        return Err(Error::Document(
            "the sample sheet has no header row".to_owned(),
        ));
    };
    let width = columns.len();
    let mut kept = Vec::new();
    for (line, mut cells) in rows {
        if cells.len() > width {
            return Err(error(
                line,
                format!("{} cells, the header names {width} columns", cells.len()),
            ));
        }
        cells.resize(width, String::new());
        kept.push(cells);
    }
    Ok(Sheet {
        columns,
        rows: kept,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_with_quotes_and_crlf() {
        let sheet =
            parse("\u{feff}Sample,BAM,sex\r\n\"NA,1\",\"a \"\"b\"\".bam\",F\r\n\r\nS2,,\r\n")
                .expect("parse");
        assert_eq!(sheet.columns, ["Sample", "BAM", "sex"]);
        assert_eq!(sheet.rows.len(), 2);
        assert_eq!(sheet.rows[0], ["NA,1", "a \"b\".bam", "F"]);
        assert_eq!(sheet.cell(&sheet.rows[1], "bam"), Some(""));
        assert_eq!(sheet.cell(&sheet.rows[1], "build"), None);
    }

    #[test]
    fn tsv_is_told_by_a_tab_in_the_first_line() {
        let sheet = parse("sample\tbuild\nS1\thg38\nS2\n").expect("parse");
        assert_eq!(sheet.rows, [["S1", "hg38"], ["S2", ""]]);
        assert_eq!(
            parse("sample\tbam\nx\\ty\tz\n").expect("tsv").rows[0][0],
            "x\ty"
        );
    }

    #[test]
    fn malformed_sheets_name_their_line() {
        let long = parse("sample\nS1,extra\n").expect_err("too many cells");
        assert_eq!(
            long.to_string(),
            "sample sheet line 2: 2 cells, the header names 1 columns"
        );
        assert!(parse("sample\n\"S1\n").is_err());
        assert!(parse("sample\nS\"1\n").is_err());
        assert!(parse("\n\n").is_err());
    }
}

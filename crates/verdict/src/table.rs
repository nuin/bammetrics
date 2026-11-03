//! The preset file syntax: the subset of TOML the presets use, read by
//! [`parse_tables`].
//!
//! A preset file is valid TOML. This reader accepts the part of TOML the
//! presets need and rejects the rest with a line number: full-line and
//! trailing `#` comments, a `[name]` table, a `[[name]]` array of tables,
//! and `key = value` lines where the value is a basic string in double
//! quotes (escapes `\"`, `\\`, `\n`, `\t`), an integer, a decimal, `true`
//! or `false`. Keys are bare (`[A-Za-z0-9_-]`). Not accepted: dotted or
//! quoted keys, inline tables, arrays, multi-line or literal strings, dates.
//!
//! Decision 4 of docs/05 (the `toml` crate or a hand-written format) stays
//! open. [`parse_tables`] is the one function the choice touches: swapping
//! to the `toml` crate replaces its 90 lines with a `toml::Table` walk and
//! adds `toml`, `serde`, `serde_spanned`, `toml_datetime`, `toml_parser`,
//! `toml_writer` and `winnow` to the build (`toml_edit` and `serde` are
//! already build dependencies of pyo3). The preset files need no change
//! either way.

use crate::error::{Error, Result};

/// One scalar value of a `key = value` line.
#[derive(Clone, Debug, PartialEq)]
pub enum Scalar {
    Str(String),
    Int(i64),
    Float(f64),
    Bool(bool),
}

impl Scalar {
    /// What the value looks like, for an error message.
    pub fn kind(&self) -> &'static str {
        match self {
            Scalar::Str(_) => "a string",
            Scalar::Int(_) => "an integer",
            Scalar::Float(_) => "a decimal",
            Scalar::Bool(_) => "a boolean",
        }
    }
}

/// A `[name]` table or one `[[name]]` element, with its entries in file
/// order and the line of its header.
#[derive(Clone, Debug, PartialEq)]
pub struct Table {
    pub name: String,
    pub array: bool,
    pub line: usize,
    pub entries: Vec<Entry>,
}

/// One `key = value` line.
#[derive(Clone, Debug, PartialEq)]
pub struct Entry {
    pub key: String,
    pub value: Scalar,
    pub line: usize,
}

/// Read `text` into its tables, in file order. `source` names the file in
/// errors. A `key = value` line before any table header is an error, since
/// the presets keep every key under `[preset]` or a `[[rule]]`.
pub fn parse_tables(text: &str, source: &str) -> Result<Vec<Table>> {
    let mut tables: Vec<Table> = Vec::new();
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        let bad = |message: String| Error::Preset {
            file: source.to_owned(),
            line: Some(line),
            message,
        };
        let content = strip_comment(raw).trim();
        if content.is_empty() {
            continue;
        }
        if let Some(inner) = content.strip_prefix("[[") {
            let name = inner
                .strip_suffix("]]")
                .ok_or_else(|| bad("expected [[name]]".to_owned()))?;
            tables.push(Table {
                name: bare_key(name.trim())
                    .ok_or_else(|| bad(format!("bad table name {name:?}")))?,
                array: true,
                line,
                entries: Vec::new(),
            });
            continue;
        }
        if let Some(inner) = content.strip_prefix('[') {
            let name = inner
                .strip_suffix(']')
                .ok_or_else(|| bad("expected [name]".to_owned()))?;
            tables.push(Table {
                name: bare_key(name.trim())
                    .ok_or_else(|| bad(format!("bad table name {name:?}")))?,
                array: false,
                line,
                entries: Vec::new(),
            });
            continue;
        }
        let (key, value) = content
            .split_once('=')
            .ok_or_else(|| bad("expected [table], [[table]] or key = value".to_owned()))?;
        let key = bare_key(key.trim()).ok_or_else(|| bad(format!("bad key {:?}", key.trim())))?;
        let value = parse_scalar(value.trim()).ok_or_else(|| {
            bad(format!(
                "{key}: expected a quoted string, a number, true or false, got {:?}",
                value.trim()
            ))
        })?;
        let table = tables
            .last_mut()
            .ok_or_else(|| bad(format!("{key} is outside any table")))?;
        table.entries.push(Entry { key, value, line });
    }
    Ok(tables)
}

/// The line without a trailing `#` comment; a `#` inside a quoted string
/// is kept.
fn strip_comment(line: &str) -> &str {
    let mut quoted = false;
    let mut escaped = false;
    for (i, c) in line.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        match c {
            '\\' if quoted => escaped = true,
            '"' => quoted = !quoted,
            '#' if !quoted => return &line[..i],
            _ => {}
        }
    }
    line
}

/// A bare key: letters, digits, `_` and `-`, at least one character.
fn bare_key(text: &str) -> Option<String> {
    let valid = !text.is_empty()
        && text
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-');
    valid.then(|| text.to_owned())
}

fn parse_scalar(text: &str) -> Option<Scalar> {
    if let Some(inner) = text.strip_prefix('"') {
        return unescape(inner.strip_suffix('"')?).map(Scalar::Str);
    }
    match text {
        "true" => return Some(Scalar::Bool(true)),
        "false" => return Some(Scalar::Bool(false)),
        _ => {}
    }
    let unsigned = text.strip_prefix(['+', '-']).unwrap_or(text);
    if unsigned.is_empty() || !unsigned.starts_with(|c: char| c.is_ascii_digit()) {
        return None;
    }
    if let Ok(int) = text.parse::<i64>() {
        return Some(Scalar::Int(int));
    }
    let plain = unsigned
        .bytes()
        .all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'e' | b'E' | b'+' | b'-'));
    if !plain {
        return None;
    }
    text.parse::<f64>().ok().map(Scalar::Float)
}

/// The text between the quotes with TOML's basic escapes applied; `None`
/// for an unknown escape or an unescaped `"` inside.
fn unescape(inner: &str) -> Option<String> {
    let mut out = String::with_capacity(inner.len());
    let mut chars = inner.chars();
    while let Some(c) = chars.next() {
        match c {
            '\\' => out.push(match chars.next()? {
                '"' => '"',
                '\\' => '\\',
                'n' => '\n',
                't' => '\t',
                _ => return None,
            }),
            '"' => return None,
            _ => out.push(c),
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_entries_and_comments_parse() {
        let text = "# a comment\n[preset]\nname = \"clinical-panel\" # trailing\nversion = 1\n\n\
                    [[rule]]\nid = \"a\"\nwarn = \"< 0.98\"\nratio = 1.5\nenabled = false\n\
                    note = \"has # inside\"\n[[rule]]\nid = \"b\\\"q\\\"\"\n";
        let tables = parse_tables(text, "t").expect("parses");
        assert_eq!(tables.len(), 3);
        assert_eq!(tables[0].name, "preset");
        assert!(!tables[0].array);
        assert_eq!(tables[0].line, 2);
        assert_eq!(
            tables[0].entries,
            [
                Entry {
                    key: "name".into(),
                    value: Scalar::Str("clinical-panel".into()),
                    line: 3
                },
                Entry {
                    key: "version".into(),
                    value: Scalar::Int(1),
                    line: 4
                }
            ]
        );
        assert!(tables[1].array);
        let values: Vec<&Scalar> = tables[1].entries.iter().map(|e| &e.value).collect();
        assert_eq!(
            values,
            [
                &Scalar::Str("a".into()),
                &Scalar::Str("< 0.98".into()),
                &Scalar::Float(1.5),
                &Scalar::Bool(false),
                &Scalar::Str("has # inside".into())
            ]
        );
        assert_eq!(tables[2].entries[0].value, Scalar::Str("b\"q\"".into()));
    }

    #[test]
    fn errors_carry_the_line() {
        let cases = [
            ("key = 1\n", 1, "outside any table"),
            ("[preset]\nname = unquoted\n", 2, "quoted string"),
            ("[preset]\nname = \"open\n", 2, "quoted string"),
            ("[preset\n", 1, "expected [name]"),
            ("[[rule]\n", 1, "expected [[name]]"),
            ("[preset]\njust words\n", 2, "key = value"),
            ("[preset]\na.b = 1\n", 2, "bad key"),
            ("[preset]\nx = \"\\q\"\n", 2, "quoted string"),
            ("[preset]\nx = 1.2.3\n", 2, "quoted string"),
            ("[preset]\nx = -\n", 2, "quoted string"),
        ];
        for (text, line, needle) in cases {
            let error = parse_tables(text, "my.toml").expect_err(text).to_string();
            assert!(error.contains(needle), "{text:?}: {error}");
            assert!(error.starts_with(&format!("my.toml:{line}:")), "{error}");
        }
    }

    #[test]
    fn scalars() {
        assert_eq!(parse_scalar("-12"), Some(Scalar::Int(-12)));
        assert_eq!(parse_scalar("+3"), Some(Scalar::Int(3)));
        assert_eq!(parse_scalar("0.5"), Some(Scalar::Float(0.5)));
        assert_eq!(parse_scalar("1e3"), Some(Scalar::Float(1000.0)));
        assert_eq!(parse_scalar("\"\""), Some(Scalar::Str(String::new())));
        assert_eq!(parse_scalar("inf"), None);
        assert_eq!(parse_scalar("nan"), None);
        assert_eq!(parse_scalar("'single'"), None);
        assert_eq!(parse_scalar("\"a\" \"b\""), None);
    }
}

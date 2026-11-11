//! The subset of TOML bammetrics reads, by [`parse_tables`]: the verdict
//! presets (`crates/verdict`) and the plan file of `bammetrics chain
//! --plan` (`crates/cli`).
//!
//! A file in this subset is valid TOML. The reader accepts the part the
//! two formats need and rejects the rest with a line number: full-line and
//! trailing `#` comments, a `[name]` table, a `[[name]]` array of tables,
//! and `key = value` lines where the value is a basic string in double
//! quotes (every escape of a TOML basic string: `\"`, `\\`, `\b`, `\t`,
//! `\n`, `\f`, `\r`, `\uXXXX`, `\UXXXXXXXX`), an integer, a decimal,
//! `true` or `false`, or a one-line array of those (`["a", "b"]`). Keys
//! are bare (`[A-Za-z0-9_-]`). A UTF-8 byte order mark in front of the
//! first line is skipped. Not accepted: dotted or quoted keys, inline
//! tables, nested or multi-line arrays, multi-line or literal strings,
//! dates, `inf` and `nan`. The error for a value in one of those forms
//! names the form (`unread`).
//!
//! Decision 4 of docs/05 (the `toml` crate or a hand-written format) is
//! hand-written. [`parse_tables`] is the one function a switch to the
//! `toml` crate replaces: a `toml::Table` walk in its place adds `toml`,
//! `serde`, `serde_spanned`, `toml_datetime`, `toml_parser`, `toml_writer`
//! and `winnow` to the build. No file needs a change either way.

use std::fmt;

/// One value of a `key = value` line.
#[derive(Clone, Debug, PartialEq)]
pub enum Scalar {
    Str(String),
    Int(i64),
    Float(f64),
    Bool(bool),
    /// A one-line array; its items are never lists themselves.
    List(Vec<Scalar>),
}

impl Scalar {
    /// What the value looks like, for an error message.
    pub fn kind(&self) -> &'static str {
        match self {
            Scalar::Str(_) => "a string",
            Scalar::Int(_) => "an integer",
            Scalar::Float(_) => "a decimal",
            Scalar::Bool(_) => "a boolean",
            Scalar::List(_) => "a list",
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

/// A line [`parse_tables`] cannot read: the 1-based line and why.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParseError {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for ParseError {}

/// `[[name]]` gives the name and true, `[name]` the name and false.
fn table_header(content: &str) -> Option<Result<(String, bool), String>> {
    let (inner, close, array) = match content.strip_prefix("[[") {
        Some(inner) => (inner, "]]", true),
        None => (content.strip_prefix('[')?, "]", false),
    };
    let shape = if array { "[[name]]" } else { "[name]" };
    let Some(name) = inner.strip_suffix(close) else {
        return Some(Err(format!("expected {shape}")));
    };
    Some(
        bare_key(name.trim())
            .map(|name| (name, array))
            .ok_or_else(|| format!("bad table name {name:?}")),
    )
}

/// Read `text` into its tables, in file order. A `key = value` line before
/// any table header is an error: both formats keep every key under a table.
pub fn parse_tables(text: &str) -> Result<Vec<Table>, ParseError> {
    let mut tables: Vec<Table> = Vec::new();
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        let bad = |message: String| ParseError { line, message };
        let content = strip_comment(raw).trim();
        if content.is_empty() {
            continue;
        }
        if let Some(header) = table_header(content) {
            let (name, array) = header.map_err(bad)?;
            tables.push(Table {
                name,
                array,
                line,
                entries: Vec::new(),
            });
            continue;
        }
        let (key, value) = content
            .split_once('=')
            .ok_or_else(|| bad("expected [table], [[table]] or key = value".to_owned()))?;
        let key = bare_key(key.trim()).ok_or_else(|| bad(format!("bad key {:?}", key.trim())))?;
        let value = parse_value(value.trim()).ok_or_else(|| {
            let form = unread(value.trim()).map_or_else(String::new, |why| format!(" ({why})"));
            bad(format!(
                "{key}: expected a quoted string, a number, true or false, got {:?}{form}",
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

/// The opening brace of an inline table.
const INLINE_TABLE: char = '\u{7b}';

/// The TOML form of a value [`parse_value`] rejects, when it is one this
/// reader leaves out: what the error says after the value.
fn unread(value: &str) -> Option<&'static str> {
    let array = value.starts_with('[');
    let item = value.trim_start_matches('[').trim_start();
    if item.starts_with("\"\"\"") || item.starts_with("'''") {
        Some("multi-line strings are not read: write one line in double quotes")
    } else if item.starts_with('\'') {
        Some("literal strings in single quotes are not read: use double quotes")
    } else if array && !value.ends_with(']') {
        Some("multi-line arrays are not read: write the array on one line")
    } else if array && value[1..].trim_start().starts_with('[') {
        Some("nested arrays are not read")
    } else if item.starts_with(INLINE_TABLE) {
        Some("inline tables are not read")
    } else if matches!(item.trim_start_matches(['+', '-']), "inf" | "nan") {
        Some("inf and nan are not read")
    } else {
        None
    }
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

/// A scalar, or a one-line array of scalars.
fn parse_value(text: &str) -> Option<Scalar> {
    match text.strip_prefix('[') {
        Some(inner) => parse_list(inner.strip_suffix(']')?).map(Scalar::List),
        None => parse_scalar(text),
    }
}

/// The items between the brackets of an array, split at the commas outside
/// quoted strings. A trailing comma is allowed, as in TOML.
fn parse_list(inner: &str) -> Option<Vec<Scalar>> {
    let mut items = Vec::new();
    let mut start = 0;
    let mut quoted = false;
    let mut escaped = false;
    for (i, c) in inner.char_indices() {
        if escaped {
            escaped = false;
            continue;
        }
        match c {
            '\\' if quoted => escaped = true,
            '"' => quoted = !quoted,
            ',' if !quoted => {
                items.push(parse_scalar(inner[start..i].trim())?);
                start = i + 1;
            }
            _ => {}
        }
    }
    let last = inner[start..].trim();
    if !last.is_empty() {
        items.push(parse_scalar(last)?);
    }
    Some(items)
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
    // `1e400` overflows to infinity: no value of an option.
    let float = text.parse::<f64>().ok().filter(|x| x.is_finite());
    float.map(Scalar::Float)
}

/// The character of `\uXXXX` (`digits` 4) or `\UXXXXXXXX` (`digits` 8):
/// exactly that many hexadecimal digits naming a Unicode scalar value.
fn unicode_escape(chars: &mut std::str::Chars, digits: usize) -> Option<char> {
    let hex: String = chars.take(digits).collect();
    let valid = hex.len() == digits && hex.bytes().all(|b| b.is_ascii_hexdigit());
    let code = valid.then(|| u32::from_str_radix(&hex, 16).ok()).flatten();
    code.and_then(char::from_u32)
}

/// The text between the quotes with the escapes of a TOML basic string
/// applied; `None` for an unknown escape or an unescaped `"` inside.
fn unescape(inner: &str) -> Option<String> {
    let mut out = String::with_capacity(inner.len());
    let mut chars = inner.chars();
    while let Some(c) = chars.next() {
        match c {
            '\\' => out.push(match chars.next()? {
                '"' => '"',
                '\\' => '\\',
                'b' => '\u{8}',
                't' => '\t',
                'n' => '\n',
                'f' => '\u{c}',
                'r' => '\r',
                'u' => unicode_escape(&mut chars, 4)?,
                'U' => unicode_escape(&mut chars, 8)?,
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
        let tables = parse_tables(text).expect("parses");
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
            ("[preset]\nx = [1, [2]]\n", 2, "quoted string"),
            ("[preset]\nx = [\"a\"\n", 2, "quoted string"),
            ("[preset]\nx = [a]\n", 2, "quoted string"),
        ];
        for (text, line, needle) in cases {
            let error = parse_tables(text).expect_err(text);
            assert!(error.message.contains(needle), "{text:?}: {error}");
            assert_eq!(error.line, line, "{text:?}");
            assert!(error.to_string().starts_with(&format!("line {line}: ")));
        }
    }

    #[test]
    fn scalars() {
        assert_eq!(parse_scalar("-12"), Some(Scalar::Int(-12)));
        assert_eq!(parse_scalar("+3"), Some(Scalar::Int(3)));
        assert_eq!(parse_scalar("0.5"), Some(Scalar::Float(0.5)));
        assert_eq!(parse_scalar("1e3"), Some(Scalar::Float(1000.0)));
        assert_eq!(parse_scalar("\"\""), Some(Scalar::Str(String::new())));
        assert_eq!(parse_scalar("1e400"), None);
        assert_eq!(parse_scalar("inf"), None);
        assert_eq!(parse_scalar("nan"), None);
        assert_eq!(parse_scalar("'single'"), None);
        assert_eq!(parse_scalar("\"a\" \"b\""), None);
    }

    #[test]
    fn every_escape_of_a_basic_string_is_read() {
        assert_eq!(
            unescape(r#"a\"b\\c\td\ne\rf\bg\fh"#).as_deref(),
            Some("a\"b\\c\td\ne\rf\u{8}g\u{c}h")
        );
        assert_eq!(unescape(r"café").as_deref(), Some("caf\u{e9}"));
        assert_eq!(unescape(r"\U0001F600!").as_deref(), Some("\u{1F600}!"));
        // Too few digits, a digit that is not hexadecimal, a surrogate, a
        // value past the last scalar, an escape TOML does not have, and a
        // backslash at the end.
        for bad in [r"\u00e", r"\u00zz", r"\uD800", r"\U00110000", r"\q", r"\"] {
            assert_eq!(unescape(bad), None, "{bad}");
        }
    }

    #[test]
    fn a_byte_order_mark_in_front_is_skipped() {
        let tables = parse_tables("\u{feff}[chain]\nsource = \"a\"\n").expect("parses");
        assert_eq!(tables[0].name, "chain");
        assert_eq!(tables[0].entries[0].value, Scalar::Str("a".into()));
        // Only in front of the first line.
        assert!(parse_tables("[chain]\n\u{feff}source = \"a\"\n").is_err());
    }

    #[test]
    fn a_rejected_value_names_the_toml_form_it_has() {
        let cases = [
            (
                "x = [",
                "(multi-line arrays are not read: write the array on one line)",
            ),
            (
                "x = [\"a\",",
                "(multi-line arrays are not read: write the array on one line)",
            ),
            (
                "x = 'a'",
                "(literal strings in single quotes are not read: use double quotes)",
            ),
            (
                "x = ['a']",
                "(literal strings in single quotes are not read: use double quotes)",
            ),
            (
                "x = \"\"\"",
                "(multi-line strings are not read: write one line in double quotes)",
            ),
            (
                "x = '''a",
                "(multi-line strings are not read: write one line in double quotes)",
            ),
            ("x = [[1], [2]]", "(nested arrays are not read)"),
            ("x = { a = 1 }", "(inline tables are not read)"),
            ("x = -inf", "(inf and nan are not read)"),
            ("x = nan", "(inf and nan are not read)"),
        ];
        for (line, why) in cases {
            let error = parse_tables(&format!("[t]\n{line}\n")).expect_err(line);
            assert!(
                error.message.contains("expected a quoted string"),
                "{error}"
            );
            assert!(error.message.ends_with(why), "{line}: {error}");
        }
        // A value of no TOML form gets no such note.
        let plain = parse_tables("[t]\nx = unquoted\n").expect_err("unquoted");
        assert!(plain.message.ends_with("got \"unquoted\""), "{plain}");
    }

    #[test]
    fn one_line_arrays_hold_scalars() {
        let text = "[step]\nlevel = [\"SAMPLE\", \"a, b\", \"q\\\"]\"]\nn = [1, 2.5, true,]\n\
                    empty = []\n";
        let tables = parse_tables(text).expect("parses");
        let values: Vec<&Scalar> = tables[0].entries.iter().map(|e| &e.value).collect();
        assert_eq!(
            values,
            [
                &Scalar::List(vec![
                    Scalar::Str("SAMPLE".into()),
                    Scalar::Str("a, b".into()),
                    Scalar::Str("q\"]".into())
                ]),
                &Scalar::List(vec![Scalar::Int(1), Scalar::Float(2.5), Scalar::Bool(true)]),
                &Scalar::List(Vec::new()),
            ]
        );
        assert_eq!(values[0].kind(), "a list");
    }
}

//! A JSON value with a writer and a reader, hand-written so the crate
//! depends on nothing but the core.
//!
//! A number keeps its text, so an `i64` above 2^53 and `-0` read back as
//! they were written. NaN and infinity have no JSON form; the callers write
//! `null` for them.

use crate::error::{Error, Result};

/// Arrays and objects nested deeper than this are refused, so a hostile
/// file cannot exhaust the stack.
const MAX_DEPTH: usize = 64;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Json {
    Null,
    Bool(bool),
    /// The number as JSON text.
    Number(String),
    Str(String),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    pub fn text(value: &str) -> Self {
        Json::Str(value.to_owned())
    }

    /// `Str` for `Some`, `Null` for `None`.
    pub fn option_text(value: Option<&str>) -> Self {
        value.map_or(Json::Null, Json::text)
    }

    pub fn int(value: i64) -> Self {
        Json::Number(value.to_string())
    }

    /// The shortest text that reads back as `value`; `Null` for NaN and
    /// infinity.
    pub fn float(value: f64) -> Self {
        if value.is_finite() {
            Json::Number(value.to_string())
        } else {
            Json::Null
        }
    }

    pub fn strings<'a>(values: impl IntoIterator<Item = &'a str>) -> Self {
        Json::Array(values.into_iter().map(Json::text).collect())
    }

    /// The member `key` of an object; `None` for a missing key or a non-object.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Object(members) => members.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }

    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Array(items) => Some(items),
            _ => None,
        }
    }

    pub fn as_object(&self) -> Option<&[(String, Json)]> {
        match self {
            Json::Object(members) => Some(members),
            _ => None,
        }
    }

    /// Pretty-printed with two-space indentation and a final newline.
    pub fn render(&self) -> String {
        let mut out = String::new();
        self.write(&mut out, 0);
        out.push('\n');
        out
    }

    /// On one line with no whitespace, for JSON Lines.
    pub fn render_line(&self) -> String {
        let mut out = String::new();
        self.write_line(&mut out);
        out
    }

    fn write_scalar(&self, out: &mut String) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Number(n) => out.push_str(n),
            Json::Str(s) => write_string(out, s),
            Json::Array(_) => out.push_str("[]"),
            Json::Object(_) => out.push_str("{}"),
        }
    }

    fn write_line(&self, out: &mut String) {
        match self {
            Json::Array(items) => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    item.write_line(out);
                }
                out.push(']');
            }
            Json::Object(members) => {
                out.push('{');
                for (i, (key, value)) in members.iter().enumerate() {
                    if i > 0 {
                        out.push(',');
                    }
                    write_string(out, key);
                    out.push(':');
                    value.write_line(out);
                }
                out.push('}');
            }
            scalar => scalar.write_scalar(out),
        }
    }

    fn write(&self, out: &mut String, depth: usize) {
        match self {
            Json::Array(items) if !items.is_empty() => {
                out.push('[');
                for (i, item) in items.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    indent(out, depth + 1);
                    item.write(out, depth + 1);
                }
                out.push('\n');
                indent(out, depth);
                out.push(']');
            }
            Json::Object(members) if !members.is_empty() => {
                out.push('{');
                for (i, (key, value)) in members.iter().enumerate() {
                    out.push_str(if i == 0 { "\n" } else { ",\n" });
                    indent(out, depth + 1);
                    write_string(out, key);
                    out.push_str(": ");
                    value.write(out, depth + 1);
                }
                out.push('\n');
                indent(out, depth);
                out.push('}');
            }
            scalar => scalar.write_scalar(out),
        }
    }

    /// Parse JSON text. Errors name the byte offset.
    pub fn parse(text: &str) -> Result<Self> {
        let mut parser = Parser {
            bytes: text.as_bytes(),
            pos: 0,
            depth: 0,
        };
        parser.skip_space();
        let value = parser.value()?;
        parser.skip_space();
        if parser.pos != parser.bytes.len() {
            return Err(parser.error("trailing text"));
        }
        Ok(value)
    }
}

fn indent(out: &mut String, depth: usize) {
    for _ in 0..depth {
        out.push_str("  ");
    }
}

fn write_string(out: &mut String, text: &str) {
    out.push('"');
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
}

struct Parser<'a> {
    bytes: &'a [u8],
    pos: usize,
    depth: usize,
}

impl Parser<'_> {
    fn error(&self, message: &str) -> Error {
        Error::Json(format!("byte {}: {message}", self.pos))
    }

    fn skip_space(&mut self) {
        while self
            .peek()
            .is_some_and(|b| matches!(b, b' ' | b'\t' | b'\n' | b'\r'))
        {
            self.pos += 1;
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn expect(&mut self, byte: u8) -> Result<()> {
        if self.peek() == Some(byte) {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.error(&format!("expected {:?}", byte as char)))
        }
    }

    fn literal(&mut self, word: &str, value: Json) -> Result<Json> {
        if self.bytes[self.pos..].starts_with(word.as_bytes()) {
            self.pos += word.len();
            Ok(value)
        } else {
            Err(self.error(&format!("expected {word}")))
        }
    }

    fn value(&mut self) -> Result<Json> {
        match self.peek() {
            Some(b'{') => self.nested(Self::object),
            Some(b'[') => self.nested(Self::array),
            Some(b'"') => self.string().map(Json::Str),
            Some(b't') => self.literal("true", Json::Bool(true)),
            Some(b'f') => self.literal("false", Json::Bool(false)),
            Some(b'n') => self.literal("null", Json::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => Err(self.error("expected a value")),
        }
    }

    fn nested(&mut self, read: fn(&mut Self) -> Result<Json>) -> Result<Json> {
        if self.depth == MAX_DEPTH {
            return Err(self.error("nested too deep"));
        }
        self.depth += 1;
        let value = read(self);
        self.depth -= 1;
        value
    }

    fn object(&mut self) -> Result<Json> {
        self.expect(b'{')?;
        let mut members = Vec::new();
        self.skip_space();
        if self.peek() == Some(b'}') {
            self.pos += 1;
            return Ok(Json::Object(members));
        }
        loop {
            self.skip_space();
            let key = self.string()?;
            self.skip_space();
            self.expect(b':')?;
            self.skip_space();
            let value = self.value()?;
            members.push((key, value));
            self.skip_space();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b'}') => {
                    self.pos += 1;
                    return Ok(Json::Object(members));
                }
                _ => return Err(self.error("expected , or }")),
            }
        }
    }

    fn array(&mut self) -> Result<Json> {
        self.expect(b'[')?;
        let mut items = Vec::new();
        self.skip_space();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(Json::Array(items));
        }
        loop {
            self.skip_space();
            items.push(self.value()?);
            self.skip_space();
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(Json::Array(items));
                }
                _ => return Err(self.error("expected , or ]")),
            }
        }
    }

    fn digits(&mut self) -> usize {
        let start = self.pos;
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            self.pos += 1;
        }
        self.pos - start
    }

    /// `-? (0 | [1-9][0-9]*) (. [0-9]+)? ([eE] [+-]? [0-9]+)?`, kept as text.
    fn number(&mut self) -> Result<Json> {
        let start = self.pos;
        if self.peek() == Some(b'-') {
            self.pos += 1;
        }
        let leading_zero = self.peek() == Some(b'0');
        let whole = self.digits();
        if whole == 0 || (leading_zero && whole > 1) {
            return Err(self.error("bad number"));
        }
        if self.peek() == Some(b'.') {
            self.pos += 1;
            if self.digits() == 0 {
                return Err(self.error("bad number"));
            }
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            self.pos += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.pos += 1;
            }
            if self.digits() == 0 {
                return Err(self.error("bad number"));
            }
        }
        let text = String::from_utf8_lossy(&self.bytes[start..self.pos]).into_owned();
        Ok(Json::Number(text))
    }

    fn string(&mut self) -> Result<String> {
        self.expect(b'"')?;
        let mut out = String::new();
        loop {
            let start = self.pos;
            while self.peek().is_some_and(|b| b != b'"' && b != b'\\') {
                self.pos += 1;
            }
            out.push_str(
                std::str::from_utf8(&self.bytes[start..self.pos])
                    .map_err(|_| self.error("bad utf8"))?,
            );
            match self.peek() {
                Some(b'"') => {
                    self.pos += 1;
                    return Ok(out);
                }
                Some(_) => {
                    self.pos += 1;
                    let escaped = self.peek().ok_or_else(|| self.error("open escape"))?;
                    self.pos += 1;
                    match escaped {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'b' => out.push('\u{8}'),
                        b'f' => out.push('\u{c}'),
                        b'u' => out.push(self.unicode_escape()?),
                        _ => return Err(self.error("unknown escape")),
                    }
                }
                None => return Err(self.error("open string")),
            }
        }
    }

    /// The four hex digits after `\u`, joined with a following low
    /// surrogate escape when the first is a high surrogate.
    fn unicode_escape(&mut self) -> Result<char> {
        let first = self.hex4()?;
        let code = if (0xD800..0xDC00).contains(&first) {
            if !self.bytes[self.pos..].starts_with(b"\\u") {
                return Err(self.error("lone surrogate"));
            }
            self.pos += 2;
            let second = self.hex4()?;
            if !(0xDC00..0xE000).contains(&second) {
                return Err(self.error("lone surrogate"));
            }
            0x10000 + ((first - 0xD800) << 10) + (second - 0xDC00)
        } else {
            first
        };
        char::from_u32(code).ok_or_else(|| self.error("bad unicode escape"))
    }

    fn hex4(&mut self) -> Result<u32> {
        let digits = self
            .bytes
            .get(self.pos..self.pos + 4)
            .and_then(|d| std::str::from_utf8(d).ok())
            .filter(|d| d.bytes().all(|b| b.is_ascii_hexdigit()))
            .ok_or_else(|| self.error("bad unicode escape"))?;
        let code = u32::from_str_radix(digits, 16).map_err(|_| self.error("bad unicode escape"))?;
        self.pos += 4;
        Ok(code)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_and_parses_back() {
        let value = Json::Object(vec![
            ("n".into(), Json::int(i64::MAX)),
            ("f".into(), Json::float(0.1)),
            ("zero".into(), Json::float(-0.0)),
            ("nan".into(), Json::float(f64::NAN)),
            ("none".into(), Json::option_text(None)),
            ("ok".into(), Json::Bool(true)),
            ("empty".into(), Json::Array(vec![])),
            ("nothing".into(), Json::Object(vec![])),
            (
                "list".into(),
                Json::Array(vec![Json::text("a\"\\\n\t\u{1}é"), Json::float(-1.5e-7)]),
            ),
        ]);
        let text = value.render();
        assert!(text.starts_with(
            "{\n  \"n\": 9223372036854775807,\n  \"f\": 0.1,\n  \"zero\": -0,\n  \"nan\": null,"
        ));
        assert!(text.contains("\"a\\\"\\\\\\n\\t\\u0001é\""));
        assert!(text.ends_with("}\n"));
        assert_eq!(Json::parse(&text).expect("parses"), value);
        assert_eq!(Json::parse(&value.render_line()).expect("parses"), value);
        assert_eq!(value.get("missing"), None);
    }

    #[test]
    fn parses_escapes_and_rejects_bad_text() {
        assert_eq!(
            Json::parse(r#""\u00e9\ud83d\ude00\/""#).expect("escapes"),
            Json::text("é😀/")
        );
        assert_eq!(
            Json::parse(" [ 1 , 2.50e+3 ] ").expect("spaces"),
            Json::Array(vec![
                Json::Number("1".into()),
                Json::Number("2.50e+3".into())
            ])
        );
        for bad in [
            "",
            "{",
            "[1,]",
            "{\"a\" 1}",
            "tru",
            "\"open",
            "\"\\q\"",
            "1 2",
            "-",
            "01",
            "1.",
            "1e",
            "+1",
            "NaN",
            "\"\\ud83d\"",
            "\"\\ud83d\\u0041\"",
            "\"\\u+123\"",
            "nul",
        ] {
            let error = Json::parse(bad).expect_err(bad).to_string();
            assert!(error.starts_with("json: byte "), "{bad:?}: {error}");
        }
        let deep = format!("{}{}", "[".repeat(100), "]".repeat(100));
        assert!(Json::parse(&deep).is_err());
        let fits = format!("{}{}", "[".repeat(64), "]".repeat(64));
        assert!(Json::parse(&fits).is_ok());
    }
}

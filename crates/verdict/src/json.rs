//! A JSON value with a writer and a reader, enough for `verdict.json`.
//!
//! Hand-written (docs/05 section 4.7) so the crate adds no dependency.
//! Numbers are `f64`; a whole number prints without a fraction. NaN and
//! infinity have no JSON form, so the verdict stores them as `null`.

use crate::error::{Error, Result};

#[derive(Clone, Debug, PartialEq)]
pub enum Json {
    Null,
    Bool(bool),
    Number(f64),
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

    /// `Number` for a finite value, `Null` otherwise.
    pub fn number(value: Option<f64>) -> Self {
        match value {
            Some(v) if v.is_finite() => Json::Number(v),
            _ => Json::Null,
        }
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

    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Number(n) => Some(*n),
            _ => None,
        }
    }

    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Array(items) => Some(items),
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

    fn write(&self, out: &mut String, depth: usize) {
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Number(n) => out.push_str(&n.to_string()),
            Json::Str(s) => write_string(out, s),
            Json::Array(items) => {
                if items.is_empty() {
                    out.push_str("[]");
                    return;
                }
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
            Json::Object(members) => {
                if members.is_empty() {
                    out.push_str("{}");
                    return;
                }
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
        }
    }

    /// Parse JSON text. Errors name the byte offset.
    pub fn parse(text: &str) -> Result<Self> {
        let mut parser = Parser {
            bytes: text.as_bytes(),
            pos: 0,
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
}

impl Parser<'_> {
    fn error(&self, message: &str) -> Error {
        Error::Json(format!("byte {}: {message}", self.pos))
    }

    fn skip_space(&mut self) {
        while self.pos < self.bytes.len() && self.bytes[self.pos].is_ascii_whitespace() {
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
            Some(b'{') => self.object(),
            Some(b'[') => self.array(),
            Some(b'"') => self.string().map(Json::Str),
            Some(b't') => self.literal("true", Json::Bool(true)),
            Some(b'f') => self.literal("false", Json::Bool(false)),
            Some(b'n') => self.literal("null", Json::Null),
            Some(b'-' | b'0'..=b'9') => self.number(),
            _ => Err(self.error("expected a value")),
        }
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

    fn number(&mut self) -> Result<Json> {
        let start = self.pos;
        while self
            .peek()
            .is_some_and(|b| b.is_ascii_digit() || matches!(b, b'-' | b'+' | b'.' | b'e' | b'E'))
        {
            self.pos += 1;
        }
        let text = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| self.error("bad number"))?;
        text.parse::<f64>()
            .map(Json::Number)
            .map_err(|_| self.error(&format!("bad number {text:?}")))
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
                Some(b'\\') => {
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
                Some(_) => unreachable!("loop stops at a quote or a backslash"),
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
            0x10000 + ((first - 0xD800) << 10) + second.wrapping_sub(0xDC00)
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
            .ok_or_else(|| self.error("short unicode escape"))?;
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
            ("schema".into(), Json::text("bammetrics.verdict/1")),
            ("n".into(), Json::Number(12.0)),
            ("f".into(), Json::Number(0.1)),
            ("nan".into(), Json::number(Some(f64::NAN))),
            ("none".into(), Json::option_text(None)),
            ("ok".into(), Json::Bool(true)),
            ("empty".into(), Json::Array(vec![])),
            ("nothing".into(), Json::Object(vec![])),
            (
                "list".into(),
                Json::Array(vec![Json::text("a\"\\\n\t\u{1}é"), Json::Number(-1.5e-7)]),
            ),
        ]);
        let text = value.render();
        assert!(text.starts_with("{\n  \"schema\": \"bammetrics.verdict/1\",\n  \"n\": 12,\n  \"f\": 0.1,\n  \"nan\": null,"));
        assert!(text.contains("\"a\\\"\\\\\\n\\t\\u0001é\""));
        assert!(text.ends_with("}\n"));
        assert_eq!(Json::parse(&text).expect("parses"), value);
        assert_eq!(value.get("n").and_then(Json::as_f64), Some(12.0));
        assert_eq!(value.get("missing"), None);
        assert_eq!(
            value
                .get("list")
                .and_then(Json::as_array)
                .map(<[Json]>::len),
            Some(2)
        );
    }

    #[test]
    fn parses_escapes_and_rejects_bad_text() {
        assert_eq!(
            Json::parse(r#""\u00e9\ud83d\ude00\/""#).expect("escapes"),
            Json::text("é😀/")
        );
        assert_eq!(
            Json::parse(" [ 1 , 2 ] ").expect("spaces"),
            Json::Array(vec![Json::Number(1.0), Json::Number(2.0)])
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
            "\"\\ud83d\"",
            "nul",
        ] {
            let error = Json::parse(bad).expect_err(bad).to_string();
            assert!(error.starts_with("verdict.json: byte "), "{bad:?}: {error}");
        }
    }
}

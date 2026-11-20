//! `RGDT`: htsjdk 5.0.0 `DateParser.parse` into an `Iso8601Date`, written
//! as `SimpleDateFormat("yyyy-MM-dd'T'HH:mm:ssZ")` writes it.
//!
//! The text is cut by `StringTokenizer` at `-`, `T`, `:`, `.`, `+` and `Z`
//! (each of those a token of its own), then read as `year[-month[-day[Thour
//! [:minute[:second[.fraction]]]]]][Z|+hh[:mm]|-hh[:mm]]` into a lenient UTC
//! `GregorianCalendar`: a field left out takes its first value (January,
//! day 1, 00:00:00), a field past its range carries into the next (month 13
//! is January of the next year), a fraction is cut or padded to 3 digits,
//! an offset of 100 or more without `:` is `hhmm`, and text after the
//! offset is not read. Every number goes through `Integer.parseInt`, so a
//! digit of any script counts. `Iso8601Date` drops the milliseconds,
//! toward zero.
//!
//! The fields go through Java's hybrid calendar
//! ([`crate::java_date`], the reader of `CalculateReadGroupChecksum`):
//! Julian before 15 October 1582, so `1582-10-10` is the Gregorian
//! `1582-10-20`, `1500-03-04` stays as written and year 0 prints as
//! `0001`, the year of its era.
//!
//! The JVM writes the instant in its own time zone (`-0700` in Edmonton in
//! winter). bammetrics writes it in UTC, `+0000`, so the file does not
//! depend on the machine; Picard under `TZ=UTC` writes the same bytes
//! (docs/07 decision rows). A year past [`MAX_YEAR`] is refused: the
//! milliseconds would overflow a 64-bit integer once the other fields
//! carry into it. A year cannot be negative: `-` is a delimiter.

use crate::java::parse_integer;

const DELIMITERS: [char; 6] = ['-', 'T', ':', '.', '+', 'Z'];

/// The largest year read: with every other field at `Integer.MAX_VALUE`
/// carried into it, the instant still fits an `i64` of milliseconds.
pub const MAX_YEAR: i64 = 100_000_000;

/// `java.util.StringTokenizer(text, "-T:.+Z", true)`.
struct Tokens<'a> {
    rest: &'a str,
}

impl<'a> Tokens<'a> {
    fn has_more(&self) -> bool {
        !self.rest.is_empty()
    }

    fn next(&mut self) -> Option<&'a str> {
        let first = self.rest.chars().next()?;
        let end = if DELIMITERS.contains(&first) {
            first.len_utf8()
        } else {
            self.rest.find(DELIMITERS).unwrap_or(self.rest.len())
        };
        let (token, rest) = self.rest.split_at(end);
        self.rest = rest;
        Some(token)
    }

    /// `DateParser.check`: false at the end, true when the next token is
    /// `expected`, else an error.
    fn check(&mut self, expected: &str) -> Result<bool, String> {
        match self.next() {
            None => Ok(false),
            Some(token) if token == expected => Ok(true),
            Some(token) => Err(format!("Missing [{expected}], found [{token}]")),
        }
    }

    fn number(&mut self) -> Result<i64, String> {
        let token = self.next().ok_or("no more tokens")?;
        parse_integer(token, 32)
    }
}

/// The calendar fields as `Calendar.set` leaves them, before they carry.
struct Fields {
    year: i64,
    /// From 0, as `Calendar.MONTH`.
    month: i64,
    day: i64,
    hour: i64,
    minute: i64,
    second: i64,
    millisecond: i64,
    /// Milliseconds the offset adds (`Calendar.add`).
    offset: i64,
}

/// `DateParser.getCalendar` up to the fields.
fn fields(text: &str) -> Result<Fields, String> {
    let mut f = Fields {
        year: 1970,
        month: 0,
        day: 1,
        hour: 0,
        minute: 0,
        second: 0,
        millisecond: 0,
        offset: 0,
    };
    let mut t = Tokens { rest: text };
    if !t.has_more() {
        return Ok(f);
    }
    f.year = t.number()?;
    if !(t.check("-")? && t.has_more()) {
        return Ok(f);
    }
    f.month = t.number()? - 1;
    if !(t.check("-")? && t.has_more()) {
        return Ok(f);
    }
    f.day = t.number()?;
    if !(t.check("T")? && t.has_more()) {
        return Ok(f);
    }
    f.hour = t.number()?;
    if !(t.check(":")? && t.has_more()) {
        return Ok(f);
    }
    f.minute = t.number()?;
    let Some(mut token) = t.next() else {
        return Ok(f);
    };
    if token == ":" {
        if !t.has_more() {
            return Err("No secondes specified".to_owned());
        }
        f.second = t.number()?;
        let Some(next) = t.next() else {
            return Ok(f);
        };
        token = next;
        if token == "." {
            f.millisecond = fraction(t.next().ok_or("no more tokens")?)?;
            let Some(next) = t.next() else {
                return Ok(f);
            };
            token = next;
        }
    }
    if token == "Z" {
        return Ok(f);
    }
    if token != "+" && token != "-" {
        return Err("only Z, + or - allowed".to_owned());
    }
    let plus = token == "+";
    if !t.has_more() {
        return Err("Missing hour field".to_owned());
    }
    let mut hours = t.number()?;
    let mut minutes = 0;
    if t.check(":")? && t.has_more() {
        minutes = t.number()?;
    } else if hours >= 100 {
        minutes = hours % 100;
        hours /= 100;
    }
    let shift = (hours * 60 + minutes) * 60_000;
    f.offset = if plus { -shift } else { shift };
    Ok(f)
}

/// The milliseconds of a fraction: its first three UTF-16 units, padded
/// with `0`.
fn fraction(token: &str) -> Result<i64, String> {
    let mut units: Vec<u16> = token.encode_utf16().collect();
    while units.len() < 3 {
        units.push(u16::from(b'0'));
    }
    let text = String::from_utf16(&units[..3]).map_err(|e| e.to_string())?;
    parse_integer(&text, 32)
}

/// The instant of `text` in milliseconds since 1970, milliseconds dropped
/// toward zero as `Iso8601Date.truncateMilliseconds` drops them.
pub fn parse(text: &str) -> Result<i64, String> {
    let f = fields(text)?;
    if f.year > MAX_YEAR {
        return Err(format!(
            "year {}: bammetrics reads years up to {MAX_YEAR}",
            f.year
        ));
    }
    let time_of_day = ((f.hour * 60 + f.minute) * 60 + f.second) * 1000 + f.millisecond;
    let millis = crate::java_date::calendar_time(f.year, f.month, f.day, time_of_day) + f.offset;
    Ok(millis - millis % 1000)
}

/// `yyyy-MM-dd'T'HH:mm:ssZ` of `millis` in UTC.
pub fn format_utc(millis: i64) -> String {
    crate::java_date::iso8601_utc(millis)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn written(text: &str) -> String {
        format_utc(parse(text).expect(text))
    }

    /// The values Picard 3.5.0 writes under `TZ=UTC` (probed by hand,
    /// `tests/concordance/test_add_or_replace_read_groups.py` holds them).
    #[test]
    fn each_form_gives_the_instant_picard_writes() {
        for (text, expected) in [
            ("2020-01-02", "2020-01-02T00:00:00+0000"),
            ("2020-01-02T03:04:05", "2020-01-02T03:04:05+0000"),
            ("2020-01-02T03:04:05Z", "2020-01-02T03:04:05+0000"),
            ("2020-01-02T03:04:05.123+05:30", "2020-01-01T21:34:05+0000"),
            ("2020-01", "2020-01-01T00:00:00+0000"),
            ("2020", "2020-01-01T00:00:00+0000"),
            ("2020-01-02T03:04", "2020-01-02T03:04:00+0000"),
            ("2020-01-02T03:04:05-0700", "2020-01-02T10:04:05+0000"),
            ("2020-13-45", "2021-02-14T00:00:00+0000"),
            ("1969-12-31T23:59:59Z", "1969-12-31T23:59:59+0000"),
            ("2020-01-02T03:04:05.5Z", "2020-01-02T03:04:05+0000"),
            ("2020-01-02T03:04:05Zjunk", "2020-01-02T03:04:05+0000"),
            ("2020-01-02T25:61:61", "2020-01-03T02:02:01+0000"),
            ("2020-01-02T03:04+05", "2020-01-01T22:04:00+0000"),
            ("١٩٩٩-٠١-٠١", "1999-01-01T00:00:00+0000"),
            ("1500-03-04T05:06:07", "1500-03-04T05:06:07+0000"),
            ("0001-01-01", "0001-01-01T00:00:00+0000"),
            ("1582-10-10", "1582-10-20T00:00:00+0000"),
            ("1582-10-15", "1582-10-15T00:00:00+0000"),
            ("1582-10-04", "1582-10-04T00:00:00+0000"),
            ("20200102", "20200102-01-01T00:00:00+0000"),
        ] {
            assert_eq!(written(text), expected, "{text}");
        }
    }

    #[test]
    fn a_fraction_before_1970_drops_toward_zero() {
        assert_eq!(
            written("1969-12-31T23:59:59.500Z"),
            "1970-01-01T00:00:00+0000"
        );
        assert_eq!(parse("1969-12-31T23:58:58.999Z"), Ok(-61_000));
    }

    /// The largest year with every other field at `Integer.MAX_VALUE`.
    #[test]
    fn the_largest_fields_do_not_overflow() {
        let max = i32::MAX;
        let text = format!("{MAX_YEAR}-{max}-{max}T{max}:{max}:{max}.999-{max}");
        assert!(parse(&text).is_ok(), "{text}");
    }

    /// Text Picard refuses, then years past [`MAX_YEAR`].
    #[test]
    fn text_outside_the_grammar_or_the_years_is_an_error() {
        for text in [
            "garbage",
            "2020-01-02 03:04",
            "2020:01",
            "2020-01-02T03Z",
            "2020-01-02T03:04:",
            "2020-01-02T03:04:05.",
            "2020-01-02T03:04:05X",
            "2020-01-02T03:04:05+",
            "2020-01-02T03:04:05+05Z",
            "100000001",
            "2147483647-13-01",
        ] {
            assert!(parse(text).is_err(), "{text}");
        }
    }
}

//! The two dates of an arrays VCF header, as Picard 3.5.0 reads them with
//! a lenient `SimpleDateFormat` and prints them as `Iso8601Date`:
//!
//! - `##autocallDate`, `MM/dd/yyyy HH:mm` in the JVM's zone, which is the
//!   process's zone here (`TZ` or `/etc/localtime`);
//! - `##imagingDate`, `MM/dd/yyyy hh:mm:ss a` in `America/New_York`.
//!
//! The parse follows `SimpleDateFormat.parse(String)`: each number skips
//! the spaces and tabs before it and takes every ASCII digit after them;
//! each literal (`/`, `:`, the space) must be the next character; text
//! after the last field is ignored. The calendar is lenient, so a month 13
//! or a day 32 rolls over. `hh` reads 12 as 0 and `a` adds 12 hours for
//! `PM`. The `AM` and `PM` markers are matched without case; Picard takes
//! the markers of the JVM's locale (`AM` and `PM` in `en_US`, `a.m.` and
//! `p.m.` in `en_CA`), so the markers here are those `GtcToVcf` writes.
//!
//! `America/New_York` is computed from the US daylight-saving rules from
//! 1967 on (2007: second Sunday of March to first Sunday of November;
//! 1987: first Sunday of April to last Sunday of October; 1976 and
//! 1967: last Sunday of April; 1974 from 6 January, 1975 from 23
//! February); an imaging date before 1967 is refused. A number with a sign
//! or over 2,147,483,647 is refused too, where Java's number parser reads
//! it.

use crate::java_date;

/// The pattern of a header date.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pattern {
    /// `MM/dd/yyyy HH:mm`.
    Autocall,
    /// `MM/dd/yyyy hh:mm:ss a`.
    Imaging,
}

const HOUR: i64 = 3_600_000;
const DAY: i64 = 24 * HOUR;

/// The text of `##autocallDate` or `##imagingDate` as `Iso8601Date`
/// prints it, `yyyy-MM-dd'T'HH:mm:ssZ` in the process's zone.
pub fn iso_date(text: &str, pattern: Pattern) -> Result<String, String> {
    let fields = parse(text, pattern)?;
    let wall = java_date::wall_millis(
        fields.year,
        fields.month - 1,
        fields.day,
        ((fields.hour * 60 + fields.minute) * 60 + fields.second) * 1000,
    );
    let instant = match pattern {
        Pattern::Autocall => java_date::local_instant(wall),
        Pattern::Imaging => Some(new_york_instant(wall)?),
    };
    instant
        .and_then(java_date::iso8601)
        .ok_or_else(|| format!("{text:?}: a date the C library cannot convert"))
}

/// The calendar fields a pattern sets, before the lenient roll-over.
#[derive(Debug, PartialEq, Eq)]
struct Fields {
    year: i64,
    month: i64,
    day: i64,
    hour: i64,
    minute: i64,
    second: i64,
}

struct Cursor<'a> {
    text: &'a str,
    at: usize,
}

impl Cursor<'_> {
    /// `subParse` of a numeric field: spaces and tabs skipped, then the
    /// digits.
    fn number(&mut self) -> Option<i64> {
        self.skip_blanks()?;
        let digits = self.text[self.at..]
            .bytes()
            .take_while(u8::is_ascii_digit)
            .count();
        let value = self.text[self.at..self.at + digits].parse::<i64>().ok()?;
        self.at += digits;
        (value <= i64::from(i32::MAX)).then_some(value)
    }

    fn skip_blanks(&mut self) -> Option<()> {
        let blanks = self.text[self.at..]
            .bytes()
            .take_while(|b| matches!(b, b' ' | b'\t'))
            .count();
        self.at += blanks;
        (self.at < self.text.len()).then_some(())
    }

    fn literal(&mut self, c: u8) -> Option<()> {
        (self.text.as_bytes().get(self.at) == Some(&c)).then(|| self.at += 1)
    }

    /// `a`: 0 for `AM`, 1 for `PM`, without case.
    fn am_pm(&mut self) -> Option<i64> {
        self.skip_blanks()?;
        let marker = self.text.get(self.at..self.at + 2)?;
        let value = if marker.eq_ignore_ascii_case("AM") {
            0
        } else if marker.eq_ignore_ascii_case("PM") {
            1
        } else {
            return None;
        };
        self.at += 2;
        Some(value)
    }
}

fn parse(text: &str, pattern: Pattern) -> Result<Fields, String> {
    let mut c = Cursor { text, at: 0 };
    let fields = (|| {
        let month = c.number()?;
        c.literal(b'/')?;
        let day = c.number()?;
        c.literal(b'/')?;
        let year = c.number()?;
        c.literal(b' ')?;
        let hour = c.number()?;
        c.literal(b':')?;
        let minute = c.number()?;
        let (hour, second) = match pattern {
            Pattern::Autocall => (hour, 0),
            Pattern::Imaging => {
                c.literal(b':')?;
                let second = c.number()?;
                c.literal(b' ')?;
                let pm = c.am_pm()?;
                // `hh` is 1 to 12: 12 is hour 0 of the half day.
                let hour = if hour == 12 { 0 } else { hour };
                (pm * 12 + hour, second)
            }
        };
        Some(Fields {
            year,
            month,
            day,
            hour,
            minute,
            second,
        })
    })();
    let format = match pattern {
        Pattern::Autocall => "MM/dd/yyyy HH:mm",
        Pattern::Imaging => "MM/dd/yyyy hh:mm:ss a",
    };
    fields.ok_or_else(|| {
        format!("{text:?} is not a date of the form {format} (Picard: ParseException)")
    })
}

/// Days from 1970-01-01 of a proleptic Gregorian date (month 1 to 12).
fn days_from_civil(year: i64, month: i64, day: i64) -> i64 {
    let y = if month <= 2 { year - 1 } else { year };
    let era = y.div_euclid(400);
    let yoe = y - era * 400;
    let mp = (month + 9) % 12;
    let doy = (153 * mp + 2) / 5 + day - 1;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468
}

/// The year of a day count from 1970-01-01.
fn year_of(days: i64) -> i64 {
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    yoe + era * 400 + i64::from(month <= 2)
}

/// The day of the `n`th Sunday of a month (`n` 1 to 5; 5 is the last).
fn sunday(year: i64, month: i64, n: i64) -> i64 {
    let first = days_from_civil(year, month, 1);
    // 1970-01-01 was a Thursday: day 3 is the first Sunday.
    let first_sunday = first + (3 - first).rem_euclid(7);
    if n < 5 {
        return first_sunday + 7 * (n - 1);
    }
    let next = days_from_civil(year + i64::from(month == 12), month % 12 + 1, 1);
    let mut last = first_sunday;
    while last + 7 < next {
        last += 7;
    }
    last
}

/// The first and last day of daylight saving in New York in `year`, from
/// 1967 on.
fn new_york_dst_days(year: i64) -> (i64, i64) {
    let october_end = sunday(year, 10, 5);
    match year {
        2007.. => (sunday(year, 3, 2), sunday(year, 11, 1)),
        1987.. => (sunday(year, 4, 1), october_end),
        1976.. => (sunday(year, 4, 5), october_end),
        1975 => (days_from_civil(1975, 2, 23), october_end),
        1974 => (days_from_civil(1974, 1, 6), october_end),
        _ => (sunday(year, 4, 5), october_end),
    }
}

/// The offset of `America/New_York` at the instant `millis`: daylight time
/// from 02:00 standard time (07:00 UTC) of its first day to 02:00 daylight
/// time (06:00 UTC) of its last.
fn new_york_offset(millis: i64) -> i64 {
    let year = year_of(millis.div_euclid(DAY));
    let (start, end) = new_york_dst_days(year);
    if (start * DAY + 7 * HOUR..end * DAY + 6 * HOUR).contains(&millis) {
        -4 * HOUR
    } else {
        -5 * HOUR
    }
}

/// The instant of a New York wall-clock time: in the hour the clocks repeat
/// the earlier one, in the hour they skip the offset before it, as
/// `ZoneInfo.getOffsetsByWall` does.
fn new_york_instant(wall: i64) -> Result<i64, String> {
    if year_of(wall.div_euclid(DAY)) < 1967 {
        return Err(
            "an imagingDate before 1967 (only the US daylight-saving rules from 1967 on are \
             known here)"
                .to_owned(),
        );
    }
    let before = new_york_offset(wall - 14 * HOUR);
    let early = wall - before;
    if new_york_offset(early) == before {
        return Ok(early);
    }
    let after = new_york_offset(wall + 14 * HOUR);
    let late = wall - after;
    Ok(if new_york_offset(late) == after {
        late
    } else {
        early
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn utc(text: &str) -> i64 {
        let f = parse(text, Pattern::Imaging).unwrap();
        new_york_instant(java_date::wall_millis(
            f.year,
            f.month - 1,
            f.day,
            ((f.hour * 60 + f.minute) * 60 + f.second) * 1000,
        ))
        .unwrap()
    }

    #[test]
    fn fields_parse_leniently() {
        let f = parse("04/18/2019 20:57", Pattern::Autocall).unwrap();
        assert_eq!(
            (f.year, f.month, f.day, f.hour, f.minute),
            (2019, 4, 18, 20, 57)
        );
        let f = parse(" 4/ 1/2019  9:05:07 pm and more", Pattern::Imaging).unwrap();
        assert_eq!((f.month, f.day, f.hour, f.second), (4, 1, 21, 7));
        let f = parse("13/45/2019 12:00:00 AM", Pattern::Imaging).unwrap();
        assert_eq!((f.month, f.day, f.hour), (13, 45, 0));
        for bad in [
            "04-18-2019 20:57",
            "04/18/2019",
            "",
            "04/18/2019 20:57:00 XM",
            "-4/18/2019 20:57:00 AM",
        ] {
            assert!(parse(bad, Pattern::Imaging).is_err(), "{bad}");
        }
    }

    #[test]
    fn new_york_follows_the_us_rules() {
        // 21:15:01 EDT is 01:15:01 UTC the next day.
        assert_eq!(
            utc("04/17/2019 09:15:01 PM"),
            days_from_civil(2019, 4, 18) * DAY + 4_501_000
        );
        // January is standard time.
        assert_eq!(
            utc("01/02/2020 12:00:00 PM"),
            days_from_civil(2020, 1, 2) * DAY + 17 * HOUR
        );
        // 2006 ended daylight time on 29 October, 2007 on 4 November.
        assert_eq!(
            new_york_offset(days_from_civil(2006, 10, 30) * DAY),
            -5 * HOUR
        );
        assert_eq!(
            new_york_offset(days_from_civil(2007, 11, 1) * DAY),
            -4 * HOUR
        );
        assert_eq!(sunday(2019, 3, 2), days_from_civil(2019, 3, 10));
        assert_eq!(sunday(2019, 10, 5), days_from_civil(2019, 10, 27));
        assert!(new_york_instant(java_date::wall_millis(1966, 0, 1, 0)).is_err());
    }
}

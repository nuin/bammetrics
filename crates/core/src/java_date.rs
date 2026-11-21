//! The dates htsjdk 5.0.0 reads from an `@RG DT` field and writes back:
//! `TextTagCodec.decodeDate`, then `Iso8601Date.toString`.
//!
//! `decodeDate` tries three readers in turn. [`decode_date`] ports two:
//!
//! 1. `DateParser.parse` (W3C ISO 8601: `YYYY[-MM[-DD[Thh[:mm[:ss[.s]]]
//!    [Z|+hh[:mm]|-hhmm]]]]`), on a lenient `GregorianCalendar` in UTC:
//!    fields out of range roll over, and a date before 15 October 1582 is
//!    Julian.
//! 2. `java.util.Date.parse`, the old free-form reader (`2014/03/24`,
//!    `Mar 24 2014 10:00 EST`), in the JVM's time zone when the text names
//!    none.
//!
//! The reader in between, `DateFormat.getDateTimeInstance().parse`, reads
//! the JVM locale's medium date and time (`Dec 31, 1969, 5:00:00 p.m.` in
//! `en_CA`); it is not ported, so such a text is a date that cannot be
//! read here (docs/07 decision rows). The milliseconds are dropped
//! (`Iso8601Date.truncateMilliseconds`).
//!
//! [`iso8601`] is `Iso8601Date.toString`: `yyyy-MM-dd'T'HH:mm:ssZ` in the
//! time zone of the process, which is the JVM's when both read `TZ` or
//! `/etc/localtime`. The zone offset of an instant comes from the C
//! library (`localtime_r`); the date fields from Java's hybrid calendar.

/// Milliseconds in a day.
const DAY: i64 = 86_400_000;
/// The fixed date (1 January 1 Gregorian is day 1) of 1 January 1970.
const EPOCH_OFFSET: i64 = 719_163;
/// The fixed date of 15 October 1582 Gregorian, the first day after the
/// Julian calendar in `GregorianCalendar`.
const CUTOVER: i64 = 577_736;
/// `gregorianCutoverYear`.
const CUTOVER_YEAR: i64 = 1582;

fn floor_div(a: i64, b: i64) -> i64 {
    a.div_euclid(b)
}

/// The fixed date of day `day` of month `month` (1 to 12, or beyond: it
/// rolls into the years) of `year`, proleptic Gregorian.
fn gregorian_fixed(year: i64, month: i64, day: i64) -> i64 {
    let year = year + floor_div(month - 1, 12);
    let month = (month - 1).rem_euclid(12) + 1;
    // Howard Hinnant's days_from_civil, from 1 March of year 0.
    let y = if month <= 2 { year - 1 } else { year };
    let era = floor_div(y, 400);
    let yoe = y - era * 400;
    let mp = (month + 9) % 12;
    let doy = (153 * mp + 2) / 5;
    let doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;
    era * 146_097 + doe - 719_468 + EPOCH_OFFSET + day - 1
}

/// [`gregorian_fixed`] in the Julian calendar.
fn julian_fixed(year: i64, month: i64, day: i64) -> i64 {
    let year = year + floor_div(month - 1, 12);
    let month = (month - 1).rem_euclid(12) + 1;
    let a = (14 - month) / 12;
    let y = year + 4800 - a;
    let m = month + 12 * a - 3;
    // The Julian day number of the first of the month, then the fixed date.
    let jdn = 1 + (153 * m + 2) / 5 + 365 * y + floor_div(y, 4) - 32_083;
    jdn - 1_721_425 + day - 1
}

/// Year, month and day of `fixed`, proleptic Gregorian.
fn gregorian_date(fixed: i64) -> (i64, i64, i64) {
    let z = fixed - EPOCH_OFFSET + 719_468;
    let era = floor_div(z, 146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    (year, month, day)
}

/// Year, month and day of `fixed` in the Julian calendar.
fn julian_date(fixed: i64) -> (i64, i64, i64) {
    let c = fixed + 1_721_425 + 32_082;
    let d = floor_div(4 * c + 3, 1461);
    let e = c - floor_div(1461 * d, 4);
    let m = (5 * e + 2) / 153;
    let day = e - (153 * m + 2) / 5 + 1;
    let month = m + 3 - 12 * (m / 10);
    let year = d - 4800 + m / 10;
    (year, month, day)
}

/// `GregorianCalendar.computeTime` for the fields `DateParser` sets: the
/// calendar of `year` before its month rolls over, Gregorian from 1582 on
/// when that date is on or after the cutover, else Julian.
pub(crate) fn calendar_time(year: i64, month0: i64, day: i64, time_of_day: i64) -> i64 {
    let (days, time_of_day) = calendar_days(year, month0, day, time_of_day);
    days * DAY + time_of_day
}

/// [`calendar_time`] as days from 1970 and the milliseconds into the day.
fn calendar_days(year: i64, month0: i64, day: i64, time_of_day: i64) -> (i64, i64) {
    let days = floor_div(time_of_day, DAY);
    let time_of_day = time_of_day.rem_euclid(DAY);
    let gregorian = days + gregorian_fixed(year, month0 + 1, day);
    let fixed = if year >= CUTOVER_YEAR && gregorian >= CUTOVER {
        gregorian
    } else {
        days + julian_fixed(year, month0 + 1, day)
    };
    (fixed - EPOCH_OFFSET, time_of_day)
}

/// `java.util.StringTokenizer(text, "-T:.+Z", true)`: runs of other
/// characters, each delimiter a token of its own.
fn iso_tokens(text: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let mut start = 0;
    for (at, c) in text.char_indices() {
        if "-T:.+Z".contains(c) {
            if start < at {
                tokens.push(&text[start..at]);
            }
            tokens.push(&text[at..at + 1]);
            start = at + 1;
        }
    }
    if start < text.len() {
        tokens.push(&text[start..]);
    }
    tokens
}

/// Why `DateParser` refuses a text.
#[derive(Debug, PartialEq, Eq)]
enum Invalid {
    /// `InvalidDateException`: `decodeDate` tries the next reader.
    Format,
    /// `NoSuchElementException` (a `.` with no fraction after it): it
    /// leaves `decodeDate` and fails the run at any stringency.
    Crash,
}
use Invalid::{Crash, Format};

/// `Integer.parseInt`.
fn int(token: &str) -> Result<i64, Invalid> {
    crate::java::parse_integer(token, 32).map_err(|_| Format)
}

/// The calendar fields `DateParser.getCalendar` sets, then the
/// milliseconds after the zone, as `getTime` gives them.
fn date_parser(text: &str) -> Result<i64, Invalid> {
    let tokens = iso_tokens(text);
    let mut at = 0;
    let mut next = || {
        let token = tokens.get(at).copied();
        at += 1;
        token
    };
    // `check`: false with no token left, else the token must be `wanted`.
    let mut fields = Fields::default();
    let Some(year) = next() else {
        return Ok(fields.time());
    };
    fields.year = int(year)?;
    macro_rules! check {
        ($wanted:expr) => {
            match next() {
                None => false,
                Some(token) if token == $wanted => true,
                Some(_) => return Err(Format),
            }
        };
    }
    // `check(st, "-") && st.hasMoreTokens()`: the second call sees the
    // token after the delimiter.
    macro_rules! field {
        ($wanted:expr) => {{ if check!($wanted) { next() } else { None } }};
    }
    let Some(month) = field!("-") else {
        return Ok(fields.time());
    };
    fields.month0 = int(month)? - 1;
    let Some(day) = field!("-") else {
        return Ok(fields.time());
    };
    fields.day = int(day)?;
    let Some(hour) = field!("T") else {
        return Ok(fields.time());
    };
    fields.hour = int(hour)?;
    let Some(minute) = field!(":") else {
        return Ok(fields.time());
    };
    fields.minute = int(minute)?;
    let Some(mut token) = next() else {
        return Ok(fields.time());
    };
    if token == ":" {
        let Some(second) = next() else {
            return Err(Format);
        };
        fields.second = int(second)?;
        let Some(after) = next() else {
            return Ok(fields.time());
        };
        token = after;
        if token == "." {
            // `nextToken` without a check: none left fails the read.
            let fraction = next().ok_or(Crash)?;
            let mut digits: String = fraction.chars().take(3).collect();
            while digits.chars().count() < 3 {
                digits.push('0');
            }
            fields.millisecond = int(&digits)?;
            let Some(after) = next() else {
                return Ok(fields.time());
            };
            token = after;
        }
    }
    let mut time = fields.time();
    if token != "Z" {
        let plus = match token {
            "+" => true,
            "-" => false,
            _ => return Err(Format),
        };
        let mut hours = int(next().ok_or(Format)?)?;
        let mut minutes = 0;
        if let Some(text) = field!(":") {
            minutes = int(text)?;
        } else if hours >= 100 {
            minutes = hours % 100;
            hours /= 100;
        }
        let sign = if plus { -1 } else { 1 };
        time += sign * (hours * 3_600_000 + minutes * 60_000);
    }
    Ok(time)
}

/// The fields of a cleared `GregorianCalendar` (all 0, the year 1970)
/// that `DateParser` sets.
struct Fields {
    year: i64,
    month0: i64,
    day: i64,
    hour: i64,
    minute: i64,
    second: i64,
    millisecond: i64,
}

impl Default for Fields {
    fn default() -> Self {
        Self {
            year: 1970,
            month0: 0,
            day: 1,
            hour: 0,
            minute: 0,
            second: 0,
            millisecond: 0,
        }
    }
}

impl Fields {
    fn time(&self) -> i64 {
        let time_of_day =
            ((self.hour * 60 + self.minute) * 60 + self.second) * 1000 + self.millisecond;
        calendar_time(self.year, self.month0, self.day, time_of_day)
    }
}

/// `Date.parse`'s words: am and pm, the days, the months, and the zones.
const WORDS: [(&str, i32); 32] = [
    ("am", 14),
    ("pm", 1),
    ("monday", 0),
    ("tuesday", 0),
    ("wednesday", 0),
    ("thursday", 0),
    ("friday", 0),
    ("saturday", 0),
    ("sunday", 0),
    ("january", 2),
    ("february", 3),
    ("march", 4),
    ("april", 5),
    ("may", 6),
    ("june", 7),
    ("july", 8),
    ("august", 9),
    ("september", 10),
    ("october", 11),
    ("november", 12),
    ("december", 13),
    ("gmt", 10_000),
    ("ut", 10_000),
    ("utc", 10_000),
    ("est", 10_000 + 5 * 60),
    ("edt", 10_000 + 4 * 60),
    ("cst", 10_000 + 6 * 60),
    ("cdt", 10_000 + 5 * 60),
    ("mst", 10_000 + 7 * 60),
    ("mdt", 10_000 + 6 * 60),
    ("pst", 10_000 + 8 * 60),
    ("pdt", 10_000 + 7 * 60),
];

/// `java.util.Date.parse`: `None` where it throws.
fn date_parse(text: &str) -> Option<i64> {
    let s: Vec<u16> = text.encode_utf16().collect();
    let unit = |c: u8| u16::from(c);
    let (mut year, mut mon, mut mday) = (i32::MIN, -1i32, -1i32);
    let (mut hour, mut min, mut sec) = (-1i32, -1i32, -1i32);
    let mut tzoffset = -1i32;
    let mut prevc = 0u16;
    let limit = s.len();
    let mut i = 0;
    let mut c: u16;
    while i < limit {
        c = s[i];
        i += 1;
        if c <= unit(b' ') || c == unit(b',') {
            continue;
        }
        if c == unit(b'(') {
            let mut depth = 1;
            while i < limit {
                c = s[i];
                i += 1;
                if c == unit(b'(') {
                    depth += 1;
                } else if c == unit(b')') {
                    depth -= 1;
                    if depth <= 0 {
                        break;
                    }
                }
            }
            continue;
        }
        let digit = |c: u16| (unit(b'0')..=unit(b'9')).contains(&c);
        if digit(c) {
            let mut n = i32::from(c - unit(b'0'));
            while i < limit && digit(s[i]) {
                c = s[i];
                n = n.wrapping_mul(10).wrapping_add(i32::from(c - unit(b'0')));
                i += 1;
            }
            // `c` is the character after the number, or its last digit.
            if i < limit {
                c = s[i];
            }
            if prevc == unit(b'+') || prevc == unit(b'-') && year != i32::MIN {
                n = if n < 24 {
                    n.wrapping_mul(60)
                } else {
                    n % 100 + n / 100 * 60
                };
                if prevc == unit(b'+') {
                    n = n.wrapping_neg();
                }
                if tzoffset != 0 && tzoffset != -1 {
                    return None;
                }
                tzoffset = n;
            } else if n >= 70 {
                if year != i32::MIN {
                    return None;
                } else if c <= unit(b' ') || c == unit(b',') || c == unit(b'/') || i >= limit {
                    year = n;
                } else {
                    return None;
                }
            } else if c == unit(b':') {
                if hour < 0 {
                    hour = n;
                } else if min < 0 {
                    min = n;
                } else {
                    return None;
                }
            } else if c == unit(b'/') {
                if mon < 0 {
                    mon = n - 1;
                } else if mday < 0 {
                    mday = n;
                } else {
                    return None;
                }
            } else if i < limit && c != unit(b',') && c > unit(b' ') && c != unit(b'-') {
                return None;
            } else if hour >= 0 && min < 0 {
                min = n;
            } else if min >= 0 && sec < 0 {
                sec = n;
            } else if mday < 0 {
                mday = n;
            } else if year == i32::MIN && mon >= 0 && mday >= 0 {
                year = n;
            } else {
                return None;
            }
            prevc = 0;
        } else if b"/:+-".map(unit).contains(&c) {
            prevc = c;
        } else {
            let start = i - 1;
            let letter = |c: u16| {
                (unit(b'A')..=unit(b'Z')).contains(&c) || (unit(b'a')..=unit(b'z')).contains(&c)
            };
            while i < limit && letter(s[i]) {
                i += 1;
            }
            if i <= start + 1 {
                return None;
            }
            // `regionMatches(true, 0, s, st, i - st)`: the word is a prefix
            // of the table entry, case ignored as Java ignores it; the last
            // entry that matches wins.
            let word = &s[start..i];
            let action = WORDS
                .iter()
                .rev()
                .find(|(entry, _)| {
                    entry.len() >= word.len()
                        && entry
                            .bytes()
                            .zip(word)
                            .all(|(e, &w)| same_ignoring_case(char::from(e), w))
                })
                .map(|(_, action)| *action)?;
            if action == 1 {
                if !(1..=12).contains(&hour) {
                    return None;
                } else if hour < 12 {
                    hour += 12;
                }
            } else if action == 14 {
                if !(1..=12).contains(&hour) {
                    return None;
                } else if hour == 12 {
                    hour = 0;
                }
            } else if (2..=13).contains(&action) {
                if mon < 0 {
                    mon = action - 2;
                } else {
                    return None;
                }
            } else if action >= 10_000 {
                tzoffset = action - 10_000;
            }
            prevc = 0;
        }
    }
    if year == i32::MIN || mon < 0 || mday < 0 {
        return None;
    }
    let mut year = i64::from(year);
    if year < 100 {
        let start = current_year()? - 80;
        year += floor_div(start, 100) * 100;
        if year < start {
            year += 100;
        }
    }
    let (sec, min, hour) = (sec.max(0), min.max(0), hour.max(0));
    // `getCalendarSystem(year)`, then `getFixedDate` with the month
    // rolled over: Gregorian from 1582 on, else Julian.
    let month = i64::from(mon) + 1;
    let fixed = if year >= CUTOVER_YEAR {
        gregorian_fixed(year, month, i64::from(mday))
    } else {
        julian_fixed(year, month, i64::from(mday))
    };
    let wall = (fixed - EPOCH_OFFSET) * DAY
        + ((i64::from(hour) * 60 + i64::from(min)) * 60 + i64::from(sec)) * 1000;
    if tzoffset == -1 {
        return from_wall(wall);
    }
    Some(wall + i64::from(tzoffset) * 60_000)
}

/// `String.regionMatches` with `ignoreCase` on one character pair: equal,
/// equal upper cases, or equal lower cases of the upper cases (Java's
/// one-character case mappings).
fn same_ignoring_case(entry: char, unit: u16) -> bool {
    let Some(c) = char::from_u32(u32::from(unit)) else {
        return false;
    };
    let upper = |c: char| {
        let mut it = c.to_uppercase();
        match (it.next(), it.next()) {
            (Some(u), None) => u,
            _ => c,
        }
    };
    let lower = |c: char| {
        let mut it = c.to_lowercase();
        match (it.next(), it.next()) {
            (Some(l), None) => l,
            _ => c,
        }
    };
    c == entry || upper(c) == upper(entry) || lower(upper(c)) == lower(upper(entry))
}

/// 1900-01-01T00:00:00Z: `ZoneInfo` holds no transition before it.
const UTC1900: i64 = -2_208_988_800_000;

/// The offset of the process's time zone at `millis`, in milliseconds east
/// of UTC, as `sun.util.calendar.ZoneInfo` gives it: `tm_gmtoff` of
/// `localtime_r` from 1900 on, and before 1900 the zone's standard offset
/// today (`getLastRawOffset`; no local mean time).
fn zone_offset(millis: i64) -> Option<i64> {
    if millis < UTC1900 {
        return raw_offset();
    }
    libc_offset(millis)
}

/// The zone's last standard offset (`ZoneInfo.rawOffset`, the last offset
/// of the zone's history, a rule already decided for later years
/// included): the offset without daylight saving at one of two instants
/// half a year apart, ten years from now, else the smaller.
fn raw_offset() -> Option<i64> {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?;
    let later = i64::try_from(now.as_millis()).ok()? + 3650 * DAY;
    let mut offsets = Vec::new();
    for at in [later, later + 182 * DAY] {
        let tm = local_tm(at)?;
        offsets.push((gmtoff_millis(&tm), tm.tm_isdst > 0));
    }
    offsets
        .iter()
        .find(|(_, dst)| !dst)
        .map(|(offset, _)| *offset)
        .or_else(|| offsets.iter().map(|(offset, _)| *offset).min())
}

/// `tm_gmtoff` of `localtime_r` at `millis`, in milliseconds.
fn libc_offset(millis: i64) -> Option<i64> {
    local_tm(millis).map(|tm| gmtoff_millis(&tm))
}

/// `tm_gmtoff` in milliseconds: a C `long`, 32 bits on some targets.
#[allow(clippy::useless_conversion)]
fn gmtoff_millis(tm: &libc::tm) -> i64 {
    i64::from(tm.tm_gmtoff) * 1000
}

/// `localtime_r` at `millis`.
fn local_tm(millis: i64) -> Option<libc::tm> {
    let seconds = libc::time_t::try_from(floor_div(millis, 1000)).ok()?;
    // SAFETY: `localtime_r` writes the zeroed `tm` it is given and reads
    // only `seconds`; a null result is checked.
    unsafe {
        let mut tm: libc::tm = std::mem::zeroed();
        if libc::localtime_r(&seconds, &mut tm).is_null() {
            return None;
        }
        Some(tm)
    }
}

/// The instant of the wall-clock time `wall` in the process's zone, as
/// `ZoneInfo.getOffsetsByWall` takes it: in an overlap the earlier
/// instant, in a gap the offset before it.
fn from_wall(wall: i64) -> Option<i64> {
    let before = zone_offset(wall - 14 * 3_600_000)?;
    let after = zone_offset(wall + 14 * 3_600_000)?;
    let early = wall - before;
    if zone_offset(early)? == before {
        return Some(early);
    }
    let late = wall - after;
    if zone_offset(late)? == after {
        return Some(late);
    }
    Some(early)
}

/// The year now in the process's zone, for `Date.parse`'s two-digit years.
fn current_year() -> Option<i64> {
    let now = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .ok()?;
    let millis = i64::try_from(now.as_millis()).ok()?;
    Some(local_fields(millis)?.0)
}

/// Year, month, day, hour, minute, second and zone offset of `millis` in
/// the process's zone, the date from Java's hybrid calendar.
fn local_fields(millis: i64) -> Option<(i64, i64, i64, i64, i64, i64, i64)> {
    Some(fields_at(millis, zone_offset(millis)?))
}

/// [`local_fields`] at a given zone `offset` in milliseconds.
fn fields_at(millis: i64, offset: i64) -> (i64, i64, i64, i64, i64, i64, i64) {
    let local = millis + offset;
    let fixed = floor_div(local, DAY) + EPOCH_OFFSET;
    let (year, month, day) = if fixed >= CUTOVER {
        gregorian_date(fixed)
    } else {
        julian_date(fixed)
    };
    let in_day = local.rem_euclid(DAY) / 1000;
    (
        year,
        month,
        day,
        in_day / 3600,
        in_day / 60 % 60,
        in_day % 60,
        offset,
    )
}

/// `decodeDate` without the `DateFormat` reader: the instant `text` names,
/// to the second; `Ok(None)` when neither ported reader reads it. An error
/// is the exception `DateParser` lets out (a `.` ending the text), which
/// fails Picard at any stringency.
pub fn decode_date(text: &str) -> Result<Option<i64>, String> {
    let millis = match date_parser(text) {
        Ok(millis) => Some(millis),
        Err(Crash) => {
            return Err(format!(
                "DT tag value '{text}': java.util.NoSuchElementException (no fraction after \
                 the point)"
            ));
        }
        Err(Format) => date_parse(text),
    };
    // `truncateMilliseconds`: Java's remainder takes the sign of `millis`.
    Ok(millis.map(|millis| millis - millis % 1000))
}

/// `Iso8601Date.toString`: `yyyy-MM-dd'T'HH:mm:ssZ` in the process's zone;
/// `None` past what the C library converts.
pub fn iso8601(millis: i64) -> Option<String> {
    local_fields(millis).map(iso8601_text)
}

/// [`iso8601`] in UTC, as a JVM under `TZ=UTC` writes it: the file does not
/// depend on the machine's zone.
pub fn iso8601_utc(millis: i64) -> String {
    iso8601_text(fields_at(millis, 0))
}

fn iso8601_text(
    (year, month, day, hour, minute, second, offset): (i64, i64, i64, i64, i64, i64, i64),
) -> String {
    // `yyyy` prints the year of the era: 1 BC for year 0.
    let year = if year <= 0 { 1 - year } else { year };
    let minutes = offset / 60_000;
    let sign = if minutes < 0 { '-' } else { '+' };
    let minutes = minutes.abs();
    format!(
        "{year:04}-{month:02}-{day:02}T{hour:02}:{minute:02}:{second:02}{sign}{:04}",
        minutes / 60 * 100 + minutes % 60
    )
}

/// The wall-clock time of lenient calendar fields (`month0` from 0, every
/// field free to roll over), as milliseconds from 1970 with no zone:
/// `GregorianCalendar.computeTime` before the zone offset. `None` more
/// than [`MAX_WALL_MILLIS`] from 1970 (about 292,000,000 years), where
/// Java's `long` arithmetic wraps.
pub fn wall_millis(year: i64, month0: i64, day: i64, time_of_day: i64) -> Option<i64> {
    let (days, time_of_day) = calendar_days(year, month0, day, time_of_day);
    days.checked_mul(DAY)?
        .checked_add(time_of_day)
        .filter(|wall| wall.unsigned_abs() <= MAX_WALL_MILLIS as u64)
}

/// The furthest wall time from 1970 [`wall_millis`] gives, a day inside the range
/// of a `long`, so a zone offset added to it cannot overflow.
pub const MAX_WALL_MILLIS: i64 = i64::MAX - DAY;

/// The instant of `wall` (from [`wall_millis`]) in the process's zone, as
/// `GregorianCalendar.computeTime` resolves it through
/// `ZoneInfo.getOffsetsByWall`; `None` past what the C library converts.
pub fn local_instant(wall: i64) -> Option<i64> {
    let before = zone_offset(wall - 14 * 3_600_000)?;
    let after = zone_offset(wall + 14 * 3_600_000)?;
    Some(wall_instant(
        wall,
        before,
        after,
        zone_offset(wall - after)?,
    ))
}

/// `getOffsetsByWall` around one transition from the offset `before` to
/// `after`: a transition takes effect at its instant plus the new offset
/// in wall time. So a wall time the clocks skip takes the offset before,
/// and one they repeat takes the offset after (the later instant).
/// `at_late` is the zone's offset at `wall - after`.
pub fn wall_instant(wall: i64, before: i64, after: i64, at_late: i64) -> i64 {
    if at_late == after {
        wall - after
    } else {
        wall - before
    }
}

#[cfg(test)]
mod tests;

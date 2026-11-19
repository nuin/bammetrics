use super::*;

/// Milliseconds since the epoch of a UTC date and time, proleptic Gregorian.
fn utc(year: i64, month: i64, day: i64, hour: i64, minute: i64, second: i64) -> i64 {
    (gregorian_fixed(year, month, day) - EPOCH_OFFSET) * DAY
        + ((hour * 60 + minute) * 60 + second) * 1000
}

#[test]
fn fixed_dates_are_javas() {
    assert_eq!(gregorian_fixed(1, 1, 1), 1);
    assert_eq!(gregorian_fixed(1970, 1, 1), EPOCH_OFFSET);
    assert_eq!(gregorian_fixed(1582, 10, 15), CUTOVER);
    assert_eq!(julian_fixed(1582, 10, 5), CUTOVER);
    assert_eq!(julian_fixed(1, 1, 3), 1);
    // Months and days roll over.
    assert_eq!(gregorian_fixed(2020, 13, 1), gregorian_fixed(2021, 1, 1));
    assert_eq!(gregorian_fixed(2020, 0, 1), gregorian_fixed(2019, 12, 1));
    assert_eq!(gregorian_fixed(2020, 2, 30), gregorian_fixed(2020, 3, 1));
    for fixed in [
        -800_000,
        -1,
        0,
        1,
        577_735,
        CUTOVER,
        EPOCH_OFFSET,
        3_000_000,
    ] {
        let (y, m, d) = gregorian_date(fixed);
        assert_eq!(gregorian_fixed(y, m, d), fixed);
        let (y, m, d) = julian_date(fixed);
        assert_eq!(julian_fixed(y, m, d), fixed);
    }
}

#[test]
fn iso_dates_read_as_date_parser_reads_them() {
    let read = |text: &str| decode_date(text).expect("no crash");
    assert_eq!(read("2020-01-02T03:04:05Z"), Some(utc(2020, 1, 2, 3, 4, 5)));
    assert_eq!(
        read("2020-01-02T03:04:05.678Z"),
        Some(utc(2020, 1, 2, 3, 4, 5))
    );
    assert_eq!(
        read("2020-01-02T03:04:05.6+0530"),
        Some(utc(2020, 1, 1, 21, 34, 5))
    );
    assert_eq!(
        read("2020-01-02T03:04:05-05:30"),
        Some(utc(2020, 1, 2, 8, 34, 5))
    );
    assert_eq!(read("2020-01-02T03:04+05"), Some(utc(2020, 1, 1, 22, 4, 0)));
    assert_eq!(read("2020-01-02T03"), Some(utc(2020, 1, 2, 3, 0, 0)));
    assert_eq!(read("2020"), Some(utc(2020, 1, 1, 0, 0, 0)));
    assert_eq!(read("2020-13-40"), Some(utc(2021, 2, 9, 0, 0, 0)));
    assert_eq!(read("2020-01-01T25:61:61Z"), Some(utc(2020, 1, 2, 2, 2, 1)));
    // Before the cutover the date is Julian: 10 October 1582 is the 20th.
    assert_eq!(read("1582-10-10"), Some(utc(1582, 10, 20, 0, 0, 0)));
    assert_eq!(
        read("1000-01-01"),
        Some((julian_fixed(1000, 1, 1) - EPOCH_OFFSET) * DAY)
    );
    // Arabic-Indic digits are digits to Integer.parseInt.
    assert_eq!(read("٢٠٢٠-01-02"), Some(utc(2020, 1, 2, 0, 0, 0)));
    assert!(decode_date("2020-01-01T10:00:00.").is_err());
}

#[test]
fn other_dates_read_as_date_parse_reads_them() {
    let read = |text: &str| decode_date(text).expect("no crash");
    assert_eq!(
        read("Mar 24 2014 10:00 GMT"),
        Some(utc(2014, 3, 24, 10, 0, 0))
    );
    assert_eq!(
        read("Tue, 24 Mar 2014 10:00:00 EST"),
        Some(utc(2014, 3, 24, 15, 0, 0))
    );
    assert_eq!(
        read("2014/03/24 10:30 UTC"),
        Some(utc(2014, 3, 24, 10, 30, 0))
    );
    assert_eq!(
        read("24 MARCH 2014 1:05 pm GMT+0200"),
        Some(utc(2014, 3, 24, 11, 5, 0))
    );
    assert_eq!(read("3/24/14 0:00 GMT"), Some(utc(2014, 3, 24, 0, 0, 0)));
    for bad in [
        "notadate",
        "x",
        "Dec 31, 1969, 5:00:00 p.m.",
        "2014 2015 3/4",
        "13:00 pm",
    ] {
        assert_eq!(read(bad), None, "{bad}");
    }
}

#[test]
fn the_iso_text_carries_the_zone_offset_of_the_instant() {
    for millis in [0, utc(2014, 3, 24, 10, 0, 0), utc(1900, 7, 1, 0, 0, 0)] {
        let text = iso8601(millis).expect("in range");
        let offset = zone_offset(millis).expect("offset") / 60_000;
        let sign = if offset < 0 { '-' } else { '+' };
        let tail = format!("{sign}{:04}", offset.abs() / 60 * 100 + offset.abs() % 60);
        assert!(text.ends_with(&tail), "{text} {tail}");
        assert_eq!(text.len(), 24, "{text}");
    }
    // A wall-clock time read back is the instant it names.
    let instant = utc(2014, 3, 24, 10, 0, 0);
    let wall = instant + zone_offset(instant).expect("offset");
    assert_eq!(from_wall(wall), Some(instant));
}

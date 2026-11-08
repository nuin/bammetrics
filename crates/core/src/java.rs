//! Java behaviour some Picard output depends on (plan item M22):
//! `String.hashCode` and the iteration order of a `HashMap<String, _>`.
//!
//! `EstimateLibraryComplexity` keeps its per-library histograms in a
//! `HashMap` and writes rows in its iteration order. That order is the
//! bucket index of each key's spread hash in the table, ties in insertion
//! order (Java 8 keeps a bucket's list order through a resize). The table
//! starts at 16 slots and doubles whenever the entry count passes three
//! quarters of it. A bucket turns into a tree only past 8 entries in one
//! bucket with 64 slots or more, which no library list reaches.

/// `String.hashCode`: `31 * h + c` over UTF-16 code units, wrapping.
pub fn string_hash_code(text: &str) -> i32 {
    text.encode_utf16()
        .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(i32::from(c)))
}

/// `HashMap.hash`: the hash code with its upper half folded in.
fn spread(hash: i32) -> i32 {
    hash ^ ((hash as u32) >> 16) as i32
}

/// Slots of a `HashMap` holding `size` entries, grown from 16 by doubling
/// while `size` exceeds the load threshold.
fn capacity(size: usize) -> usize {
    let mut capacity = 16usize;
    while size > capacity * 3 / 4 {
        capacity *= 2;
    }
    capacity
}

/// Indices of `keys` (given in insertion order, no repeats) in the order a
/// `HashMap<String, _>` iterates them.
pub fn hash_map_order(keys: &[&str]) -> Vec<usize> {
    let hashes: Vec<i32> = keys.iter().map(|k| string_hash_code(k)).collect();
    hash_map_order_by_hash(&hashes)
}

/// [`hash_map_order`] for keys of any class: `hashes` are their `hashCode`
/// values in insertion order. A `HashSet` iterates like the `HashMap` it
/// wraps.
pub fn hash_map_order_by_hash(hashes: &[i32]) -> Vec<usize> {
    let mask = (capacity(hashes.len()) - 1) as i32;
    let mut order: Vec<usize> = (0..hashes.len()).collect();
    order.sort_by_key(|&i| spread(hashes[i]) & mask);
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_code_matches_java() {
        assert_eq!(string_hash_code(""), 0);
        assert_eq!(string_hash_code("Unknown"), 1_379_812_394);
        assert_eq!(string_hash_code("lib1"), 3_321_420);
        assert_eq!(string_hash_code("LibA"), 2_368_124);
        assert_eq!(string_hash_code("polygenelubricants"), i32::MIN);
    }

    /// Each expectation is the `keySet()` order a JVM (OpenJDK 27) printed
    /// for a `HashMap` filled in the given order.
    #[test]
    fn iteration_order_is_by_bucket_then_insertion() {
        let order = hash_map_order(&["L2", "L1", "Unknown", "L3"]);
        assert_eq!(order, [2, 1, 0, 3]);
        let order = hash_map_order(&["LibD", "LibA", "Unknown", "LibC", "LibB"]);
        assert_eq!(order, [2, 1, 4, 3, 0]);
        // "Aa" and "BB" share a hash code: insertion order within the bucket.
        assert_eq!(hash_map_order(&["Aa", "BB"]), [0, 1]);
        assert_eq!(hash_map_order(&["BB", "Aa"]), [0, 1]);
        // 14 keys: the table has grown to 32 slots.
        let keys: Vec<String> = (1..=14).map(|i| format!("a{i}")).collect();
        let keys: Vec<&str> = keys.iter().map(String::as_str).collect();
        let order: Vec<&str> = hash_map_order(&keys).into_iter().map(|i| keys[i]).collect();
        assert_eq!(
            order,
            [
                "a11", "a10", "a13", "a12", "a14", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8",
                "a9"
            ]
        );
        assert_eq!(capacity(12), 16);
        assert_eq!(capacity(13), 32);
        assert_eq!(capacity(25), 64);
    }
}

/// A `java.util.HashSet<String>` reduced to its iteration order, for text
/// Picard builds by joining one (`CompareMetrics` names the columns two
/// files do not share in this order). Each key keeps its insertion number;
/// iteration is by bucket index at the current table size, then insertion.
/// Removing a key leaves the others in place.
#[derive(Clone, Debug)]
pub struct StringHashSet {
    /// Table slots, a power of two.
    capacity: usize,
    keys: Vec<(String, u64)>,
    inserted: u64,
}

impl Default for StringHashSet {
    fn default() -> Self {
        Self::new()
    }
}

impl StringHashSet {
    /// `new HashSet<>()`: 16 slots.
    pub fn new() -> Self {
        Self {
            capacity: 16,
            keys: Vec::new(),
            inserted: 0,
        }
    }

    /// `new HashSet<>(collection)`: the JDK sizes the table for
    /// `max(n, 12)` entries (`HashMap.newHashMap`, `ceil(n / 0.75)` slots
    /// rounded up to a power of two), then adds the keys in the order given.
    pub fn copy_of<'a>(keys: impl IntoIterator<Item = &'a str>) -> Self {
        let keys: Vec<&str> = keys.into_iter().collect();
        let entries = keys.len().max(12);
        let mut set = Self {
            capacity: (entries * 4).div_ceil(3).next_power_of_two(),
            ..Self::new()
        };
        for key in keys {
            set.insert(key);
        }
        set
    }

    /// `add`: false when the key is there already. The table doubles once
    /// the entry count passes three quarters of it.
    pub fn insert(&mut self, key: &str) -> bool {
        if self.contains(key) {
            return false;
        }
        self.keys.push((key.to_owned(), self.inserted));
        self.inserted += 1;
        if self.keys.len() > self.capacity * 3 / 4 {
            self.capacity *= 2;
        }
        true
    }

    pub fn remove(&mut self, key: &str) {
        self.keys.retain(|(k, _)| k != key);
    }

    pub fn contains(&self, key: &str) -> bool {
        self.keys.iter().any(|(k, _)| k == key)
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// `Set.equals`: the same keys, in any order.
    pub fn same_keys(&self, other: &Self) -> bool {
        self.len() == other.len() && self.keys.iter().all(|(k, _)| other.contains(k))
    }

    /// The keys in iteration order.
    pub fn iter(&self) -> Vec<&str> {
        let mask = (self.capacity - 1) as i32;
        let mut keys: Vec<&(String, u64)> = self.keys.iter().collect();
        keys.sort_by_key(|(k, n)| (spread(string_hash_code(k)) & mask, *n));
        keys.into_iter().map(|(k, _)| k.as_str()).collect()
    }
}

/// `String.trim`: every leading and trailing char up to U+0020.
pub fn trim(text: &str) -> &str {
    text.trim_matches(|c: char| c <= ' ')
}

/// `For input string: "<text>"`, the message of most `NumberFormatException`s.
fn bad_number(text: &str) -> String {
    format!("For input string: \"{text}\"")
}

/// The message `FloatingDecimal` gives a mantissa with a second `.`.
const MULTIPLE_POINTS: &str = "multiple points";

/// A Java floating-point literal after its sign: decimal text Rust parses
/// (`NaN` and `inf` included), or a hexadecimal value, which is rounded
/// once to the target type.
enum Literal {
    Decimal(String),
    Hex(HexLiteral),
}

/// `Double.parseDouble` up to the rounding: the string trimmed, an optional
/// sign, `NaN`, `Infinity`, or digits with an optional `.` and exponent and
/// an optional `f`, `F`, `d` or `D` suffix. The flag is the minus sign.
fn java_literal(text: &str) -> Result<(bool, Literal), String> {
    let trimmed = trim(text);
    if trimmed.is_empty() {
        return Err("empty String".to_owned());
    }
    let (negative, body) = match trimmed.strip_prefix(['+', '-']) {
        Some(rest) => (trimmed.starts_with('-'), rest),
        None => (false, trimmed),
    };
    match body {
        "NaN" => return Ok((negative, Literal::Decimal("NaN".to_owned()))),
        "Infinity" => return Ok((negative, Literal::Decimal("inf".to_owned()))),
        _ => {}
    }
    if let Some(hex) = body.strip_prefix("0x").or_else(|| body.strip_prefix("0X")) {
        if second_point(hex, |b| b.is_ascii_hexdigit()) {
            return Err(MULTIPLE_POINTS.to_owned());
        }
        let hex = hex.strip_suffix(['f', 'F', 'd', 'D']).unwrap_or(hex);
        let value = hex_literal(hex).ok_or_else(|| bad_number(text))?;
        return Ok((negative, Literal::Hex(value)));
    }
    if second_point(body, |b| b.is_ascii_digit()) {
        return Err(MULTIPLE_POINTS.to_owned());
    }
    let body = body.strip_suffix(['f', 'F', 'd', 'D']).unwrap_or(body);
    let (mantissa, exponent) = match body.split_once(['e', 'E']) {
        Some((m, e)) => (m, Some(e.strip_prefix(['+', '-']).unwrap_or(e))),
        None => (body, None),
    };
    let digits = |s: &str| s.bytes().all(|b| b.is_ascii_digit());
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let mantissa_ok =
        digits(whole) && digits(fraction) && !(whole.is_empty() && fraction.is_empty());
    let exponent_ok = exponent.is_none_or(|e| !e.is_empty() && digits(e));
    if !mantissa_ok || !exponent_ok {
        return Err(bad_number(text));
    }
    Ok((negative, Literal::Decimal(body.to_owned())))
}

/// Whether the leading run of digits and points holds two points, which
/// `FloatingDecimal` reports before it checks the rest (`1..5x` and
/// `0x1..p0` are `multiple points`, `1a.5.` is not).
fn second_point(text: &str, digit: impl Fn(u8) -> bool) -> bool {
    let mut seen = false;
    for byte in text.bytes() {
        if byte == b'.' {
            if seen {
                return true;
            }
            seen = true;
        } else if !digit(byte) {
            break;
        }
    }
    false
}

/// A hexadecimal literal: its first 30 significant digits, whether a digit
/// past them is nonzero, and the power of two that scales the digits.
struct HexLiteral {
    digits: u128,
    sticky: bool,
    shift: i64,
}

impl HexLiteral {
    /// The value rounded once, half to even, to a binary format of
    /// `precision` significand bits whose least subnormal is 2^`min_lsb`
    /// (53 and -1074 for a double, 24 and -149 for a float). The result is
    /// exact in a double; past the format's range it is infinite.
    fn round(&self, precision: i64, min_lsb: i64) -> f64 {
        if self.digits == 0 {
            return 0.0;
        }
        let bits = 128 - i64::from(self.digits.leading_zeros());
        let top = self.shift + bits - 1;
        let lsb = (top - (precision - 1)).max(min_lsb);
        let dropped = lsb - self.shift;
        let units = if dropped <= 0 {
            // Fewer than `precision` bits, so the shift cannot overflow.
            self.digits << dropped.unsigned_abs()
        } else if dropped >= 128 {
            // At most 120 bits, so the value is below half a unit.
            0
        } else {
            let kept = self.digits >> dropped;
            let rest = self.digits & ((1u128 << dropped) - 1);
            let half = 1u128 << (dropped - 1);
            let up = rest > half || (rest == half && (self.sticky || kept & 1 == 1));
            kept + u128::from(up)
        };
        units as f64 * pow2(lsb)
    }
}

/// 2^`exponent` as a double: infinite above 1023, zero below -1074.
fn pow2(exponent: i64) -> f64 {
    if exponent > 1023 {
        f64::INFINITY
    } else if exponent >= -1022 {
        f64::from_bits(((exponent + 1023) as u64) << 52)
    } else if exponent >= -1074 {
        f64::from_bits(1u64 << (exponent + 1074))
    } else {
        0.0
    }
}

/// `<hex digits>[.<hex digits>](p|P)[sign]<digits>`, after the `0x`.
fn hex_literal(hex: &str) -> Option<HexLiteral> {
    let (mantissa, exponent) = hex.split_once(['p', 'P'])?;
    let unsigned = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
    if unsigned.is_empty() || !unsigned.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    // Past 2^40 every literal is zero or infinite; the clamp keeps the shift
    // arithmetic in range.
    let limit = 1i64 << 40;
    let fallback = if exponent.starts_with('-') {
        -limit
    } else {
        limit
    };
    let exponent = exponent
        .parse::<i64>()
        .map_or(fallback, |e| e.clamp(-limit, limit));
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = format!("{whole}{fraction}");
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_hexdigit()) {
        return None;
    }
    let significant = digits.trim_start_matches('0');
    let (kept, rest) = significant.split_at(significant.len().min(30));
    let value = if kept.is_empty() {
        0
    } else {
        u128::from_str_radix(kept, 16).ok()?
    };
    Some(HexLiteral {
        digits: value,
        sticky: rest.bytes().any(|b| b != b'0'),
        shift: exponent + 4 * rest.len() as i64 - 4 * fraction.len() as i64,
    })
}

/// `Double.parseDouble`; the error is the `NumberFormatException` message.
pub fn parse_double(text: &str) -> Result<f64, String> {
    let (negative, literal) = java_literal(text)?;
    let value = match literal {
        Literal::Decimal(body) => body.parse::<f64>().map_err(|_| bad_number(text))?,
        Literal::Hex(hex) => hex.round(53, -1074),
    };
    Ok(if negative { -value } else { value })
}

/// `Float.parseFloat`, rounded once from the literal.
pub fn parse_float(text: &str) -> Result<f32, String> {
    let (negative, literal) = java_literal(text)?;
    let value = match literal {
        Literal::Decimal(body) => body.parse::<f32>().map_err(|_| bad_number(text))?,
        // Exact in a double, so the conversion does not round again.
        Literal::Hex(hex) => hex.round(24, -149) as f32,
    };
    Ok(if negative { -value } else { value })
}

/// The first code point of each run of ten decimal digits (Unicode category
/// Nd) in the Basic Multilingual Plane, as OpenJDK 27 `Character.digit`
/// reads them (370 digits). A digit outside the plane is a surrogate pair
/// in Java, and neither half is a digit.
const DIGIT_ZEROS: [u32; 37] = [
    0x0030, 0x0660, 0x06F0, 0x07C0, 0x0966, 0x09E6, 0x0A66, 0x0AE6, 0x0B66, 0x0BE6, 0x0C66, 0x0CE6,
    0x0D66, 0x0DE6, 0x0E50, 0x0ED0, 0x0F20, 0x1040, 0x1090, 0x17E0, 0x1810, 0x1946, 0x19D0, 0x1A80,
    0x1A90, 0x1B50, 0x1BB0, 0x1C40, 0x1C50, 0xA620, 0xA8D0, 0xA900, 0xA9D0, 0xA9F0, 0xAA50, 0xABF0,
    0xFF10,
];

/// `Character.digit(c, 10)`.
fn decimal_digit(c: char) -> Option<i128> {
    let code = u32::from(c);
    DIGIT_ZEROS
        .iter()
        .find(|&&zero| (zero..zero + 10).contains(&code))
        .map(|zero| i128::from(code - zero))
}

/// `Long.parseLong` (`bits` 64) or `Integer.parseInt` (32): an optional
/// ASCII sign, then decimal digits of any script (`١٢` is 12), inside the
/// type's range. Nothing is trimmed.
pub fn parse_integer(text: &str, bits: u32) -> Result<i64, String> {
    let negative = text.starts_with('-');
    let body = text.strip_prefix(['+', '-']).unwrap_or(text);
    if body.is_empty() {
        return Err(bad_number(text));
    }
    let limit = 1i128 << (bits - 1);
    let mut magnitude: i128 = 0;
    for c in body.chars() {
        let digit = decimal_digit(c).ok_or_else(|| bad_number(text))?;
        magnitude = magnitude * 10 + digit;
        if magnitude > limit {
            return Err(bad_number(text));
        }
    }
    let value = if negative { -magnitude } else { magnitude };
    if value >= limit {
        return Err(bad_number(text));
    }
    Ok(value as i64)
}

/// `Short.parseShort`: `Integer.parseInt`, then its own range message.
pub fn parse_short(text: &str) -> Result<i64, String> {
    let value = parse_integer(text, 32)?;
    if !(-32768..=32767).contains(&value) {
        return Err(format!("Value out of range. Value:\"{text}\" Radix:10"));
    }
    Ok(value)
}

/// `Double.toString` (JDK 19 and later, the form string concatenation
/// prints): the shortest digits that round-trip, as a plain decimal with at
/// least one fraction digit from 10^-3 up to 10^7 (`0.001`, `123.0`), else
/// as `d.dddE<exponent>` (`1.0E7`, `9.99E-4`). Of two shortest candidates
/// equally close, the one with the even last digit is printed. Where the
/// shortest is one digit of a subnormal, the closest two-digit decimal that
/// round-trips is printed (`Double.MIN_VALUE` is `4.9E-324`, not `5.0E-324`).
pub fn double_to_string(value: f64) -> String {
    match special_text(
        value.is_nan(),
        value.is_infinite(),
        value == 0.0,
        value.is_sign_negative(),
    ) {
        Some(text) => text,
        None => {
            let (mut digits, mut exponent) = crate::metrics_file::shortest_digits(value.abs());
            if digits.len() == 1 && value.abs() < f64::MIN_POSITIVE {
                let two = format!("{:.1e}", value.abs());
                if two.parse::<f64>() == Ok(value.abs()) {
                    (digits, exponent) = scientific_digits(&two);
                }
            }
            signed_layout(value.is_sign_negative(), &digits, exponent)
        }
    }
}

/// `Float.toString`: the same rules over the digits of the float
/// (`1234567.25f` is `1234567.2`, `Float.MIN_VALUE` is `1.4E-45`).
pub fn float_to_string(value: f32) -> String {
    match special_text(
        value.is_nan(),
        value.is_infinite(),
        value == 0.0,
        value.is_sign_negative(),
    ) {
        Some(text) => text,
        None => {
            let (mut digits, mut exponent) = scientific_digits(&format!("{:e}", value.abs()));
            // A float converts to a double exactly, so its exact digits are
            // the same. Below a power of two the lower candidate can fall
            // outside the float's rounding interval; it is kept only when it
            // reads back as the same float.
            let mut even = digits.clone();
            crate::metrics_file::prefer_even_on_tie(f64::from(value.abs()), &mut even);
            let reads_back = |d: &str| {
                format!("{}.{}e{exponent}", &d[..1], &d[1..]).parse::<f32>() == Ok(value.abs())
            };
            if even != digits && reads_back(&even) {
                digits = even;
            }
            if digits.len() == 1 && value.abs() < f32::MIN_POSITIVE {
                let two = format!("{:.1e}", value.abs());
                if two.parse::<f32>() == Ok(value.abs()) {
                    (digits, exponent) = scientific_digits(&two);
                }
            }
            signed_layout(value.is_sign_negative(), &digits, exponent)
        }
    }
}

/// The digits and exponent of Rust's `d.ddde<exp>` text.
fn scientific_digits(scientific: &str) -> (String, i64) {
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((scientific, "0"));
    let digits = mantissa.chars().filter(char::is_ascii_digit).collect();
    (digits, exponent.parse().unwrap_or(0))
}

/// `NaN`, `Infinity` and the zeros, which have no digits to lay out.
fn special_text(nan: bool, infinite: bool, zero: bool, negative: bool) -> Option<String> {
    let sign = if negative { "-" } else { "" };
    if nan {
        Some("NaN".to_owned())
    } else if infinite {
        Some(format!("{sign}Infinity"))
    } else if zero {
        Some(format!("{sign}0.0"))
    } else {
        None
    }
}

/// The JDK layout of significant `digits` whose first digit has decimal
/// `exponent`.
fn signed_layout(negative: bool, digits: &str, exponent: i64) -> String {
    let sign = if negative { "-" } else { "" };
    if (-3..7).contains(&exponent) {
        if exponent < 0 {
            let zeros = "0".repeat(exponent.unsigned_abs() as usize - 1);
            return format!("{sign}0.{zeros}{digits}");
        }
        let point = exponent as usize + 1;
        if digits.len() <= point {
            let zeros = "0".repeat(point - digits.len());
            return format!("{sign}{digits}{zeros}.0");
        }
        return format!("{sign}{}.{}", &digits[..point], &digits[point..]);
    }
    let (first, rest) = digits.split_at(1);
    let rest = if rest.is_empty() { "0" } else { rest };
    format!("{sign}{first}.{rest}E{exponent}")
}

#[cfg(test)]
mod lang_tests {
    use super::*;

    /// Each expectation is what OpenJDK 27 printed for the same value.
    #[test]
    fn doubles_print_like_double_to_string() {
        let cases = [
            (0.012_000_000_000_000_01, "0.01200000000000001"),
            (1e-3, "0.001"),
            (9.99e-4, "9.99E-4"),
            (1e7, "1.0E7"),
            (9_999_999.0, "9999999.0"),
            (123_456_789.0, "1.23456789E8"),
            (1e-7, "1.0E-7"),
            (1e300, "1.0E300"),
            (f64::MAX, "1.7976931348623157E308"),
            (100.04 - 100.0, "0.04000000000000625"),
            (0.1 + 0.2, "0.30000000000000004"),
            (-0.0, "-0.0"),
            (0.0, "0.0"),
            (2e23, "2.0E23"),
            (1.0 / 3.0, "0.3333333333333333"),
            (4.35, "4.35"),
            (1e21, "1.0E21"),
            (1e16, "1.0E16"),
            (123.0, "123.0"),
            (-500.0, "-500.0"),
            (f64::from_bits(1), "4.9E-324"),
            (f64::from_bits(2), "9.9E-324"),
            (f64::from_bits(3), "1.5E-323"),
            (f64::from_bits(10), "4.9E-323"),
            (-1e-320, "-1.0E-320"),
            (2.5e-323, "2.5E-323"),
            (f64::NAN, "NaN"),
            (f64::NEG_INFINITY, "-Infinity"),
        ];
        for (value, text) in cases {
            assert_eq!(double_to_string(value), text, "{value:e}");
        }
    }

    #[test]
    fn floats_print_like_float_to_string() {
        let cases = [
            (0.1f32, "0.1"),
            (1.0f32 / 3.0, "0.33333334"),
            (1e-4f32, "1.0E-4"),
            (1e7f32, "1.0E7"),
            (16_777_216f32, "1.6777216E7"),
            (f32::MAX, "3.4028235E38"),
            (0.3f32, "0.3"),
            (1_234_567f32 + 0.25, "1234567.2"),
            (3_897_416f32 + 0.25, "3897416.2"),
            (-1_234_567f32 - 0.75, "-1234567.8"),
            (16_777_215f32, "1.6777215E7"),
            (9.5f32, "9.5"),
            (f32::from_bits(1), "1.4E-45"),
            (f32::from_bits(3), "4.2E-45"),
            (f32::from_bits(4), "5.6E-45"),
            (f32::from_bits(7), "9.8E-45"),
            (1e-40f32, "1.0E-40"),
            (1e-39f32, "1.0E-39"),
        ];
        for (value, text) in cases {
            assert_eq!(float_to_string(value), text, "{value:e}");
        }
    }

    /// The `NumberFormatException` messages OpenJDK 27 gives.
    #[test]
    fn numbers_parse_like_java() {
        assert_eq!(
            parse_integer("", 64),
            Err("For input string: \"\"".to_owned())
        );
        assert_eq!(
            parse_integer("-", 64),
            Err("For input string: \"-\"".to_owned())
        );
        assert!(parse_integer(" 5", 64).is_err());
        assert!(parse_integer("9223372036854775808", 64).is_err());
        assert_eq!(parse_integer("9223372036854775807", 64), Ok(i64::MAX));
        assert!(parse_integer("2147483648", 32).is_err());
        assert_eq!(parse_integer("+5", 32), Ok(5));
        assert_eq!(
            parse_short("70000"),
            Err("Value out of range. Value:\"70000\" Radix:10".to_owned())
        );
        assert_eq!(parse_double(""), Err("empty String".to_owned()));
        assert_eq!(parse_double("   "), Err("empty String".to_owned()));
        assert_eq!(
            parse_double("abc"),
            Err("For input string: \"abc\"".to_owned())
        );
        assert_eq!(
            parse_double("1e"),
            Err("For input string: \"1e\"".to_owned())
        );
        assert_eq!(parse_double("."), Err("For input string: \".\"".to_owned()));
        assert_eq!(
            parse_double("inf"),
            Err("For input string: \"inf\"".to_owned())
        );
        assert!(parse_double("+NaN").is_ok_and(f64::is_nan));
        assert_eq!(parse_double("-Infinity"), Ok(f64::NEG_INFINITY));
        assert_eq!(parse_double("1."), Ok(1.0));
        assert_eq!(parse_double(".5f"), Ok(0.5));
        assert_eq!(parse_double(" 1 "), Ok(1.0));
        assert_eq!(parse_double("1e400"), Ok(f64::INFINITY));
        assert_eq!(parse_float("1e40"), Ok(f32::INFINITY));
        assert_eq!(parse_double("0x1p3"), Ok(8.0));
        assert_eq!(parse_double("-0x1.8P-1d"), Ok(-0.75));
        assert_eq!(parse_double("0x0p9999"), Ok(0.0));
        assert!(parse_double("0x1").is_err());
        for text in [
            "1..5", "..5", "0.5.", "1.2.3e4", ".5.5", "-1..5", " 1..5", "0x1..p0",
        ] {
            assert_eq!(
                parse_double(text),
                Err("multiple points".to_owned()),
                "{text}"
            );
        }
        for text in ["1a.5.", "1.5e3.2", "1e5..", "0x1.p0.", "0xg..p0"] {
            assert_eq!(parse_double(text), Err(bad_number(text)), "{text}");
        }
        assert_eq!(parse_float("1..5"), Err("multiple points".to_owned()));
    }

    /// Hexadecimal literals round once, half to even, as OpenJDK 27 rounds
    /// them (the raw bits Java printed).
    #[test]
    fn hex_literals_round_like_java() {
        let doubles = [
            ("0x0.80000001p-1074", 0x1),
            ("0x1.00000000000008p0", 0x3ff0_0000_0000_0000),
            ("0x1.000000000000081p0", 0x3ff0_0000_0000_0001),
            (
                "0x1.0000000000000800000000000000000001p0",
                0x3ff0_0000_0000_0001,
            ),
            ("0x.0000000000000000000000000000000000000000001p-900", 0x4),
            ("0xfffffffffffffffffffffffffffffffffffffp-1200", 0x40_0000),
            ("0x1.fffffffffffff8p1023", 0x7ff0_0000_0000_0000),
            ("0x1.fffffffffffff7ffffp1023", 0x7fef_ffff_ffff_ffff),
            ("0x0.0000000000001p-1022", 0x1),
            ("0x0.00000000000008p-1022", 0x0),
            ("0x0.000000000000080001p-1022", 0x1),
            ("0x0.00000000000018p-1022", 0x2),
        ];
        for (text, bits) in doubles {
            assert_eq!(parse_double(text).map(f64::to_bits), Ok(bits), "{text}");
        }
        let floats = [
            ("0x1.000001p0", 0x3f80_0000),
            ("0x1.0000010000001p0", 0x3f80_0001),
            ("0x1.00000100000000000000001p0", 0x3f80_0001),
            ("0x1p-150", 0x0),
            ("0x1.0000000001p-150", 0x1),
        ];
        for (text, bits) in floats {
            assert_eq!(parse_float(text).map(f32::to_bits), Ok(bits), "{text}");
        }
    }

    /// `Long.parseLong` reads a decimal digit of any script in the Basic
    /// Multilingual Plane, as OpenJDK 27 does.
    #[test]
    fn integers_read_unicode_digits() {
        assert_eq!(parse_integer("\u{661}\u{662}", 64), Ok(12));
        assert_eq!(parse_integer("\u{96b}", 32), Ok(5));
        assert_eq!(parse_integer("\u{ff11}\u{ff12}", 64), Ok(12));
        assert_eq!(parse_integer("-\u{ff11}", 64), Ok(-1));
        assert_eq!(parse_integer("+\u{663}", 32), Ok(3));
        assert_eq!(parse_integer("1\u{662}", 64), Ok(12));
        assert_eq!(
            parse_integer("922337203685477580\u{ff17}", 64),
            Ok(i64::MAX)
        );
        assert_eq!(
            parse_integer("-922337203685477580\u{ff18}", 64),
            Ok(i64::MIN)
        );
        assert!(parse_integer("922337203685477580\u{ff18}", 64).is_err());
        assert!(parse_integer("214748364\u{ff18}", 32).is_err());
        assert_eq!(parse_integer("-214748364\u{ff18}", 32), Ok(-2_147_483_648));
        // Mathematical bold digits lie outside the plane.
        assert!(parse_integer("\u{1d7cf}\u{1d7d0}", 64).is_err());
        assert_eq!(
            parse_integer(&format!("{}1", "\u{ff10}".repeat(19)), 64),
            Ok(1)
        );
    }

    /// Orders an OpenJDK 27 `HashSet` printed for the same insertions.
    #[test]
    fn hash_set_iterates_like_java() {
        let keys =
            |n: usize| -> Vec<String> { (0..n).map(|i| format!("K{}", i * 7919 % 1000)).collect() };
        let expected = [
            (3, "K838,K919,K0"),
            (
                12,
                "K838,K271,K919,K595,K0,K352,K676,K433,K757,K109,K514,K190",
            ),
            (
                13,
                "K271,K595,K0,K352,K676,K433,K757,K109,K514,K838,K28,K919,K190",
            ),
            (
                24,
                "K380,K271,K595,K0,K299,K461,K352,K676,K137,K542,K433,K785,K757,K218,K109,K623,\
                 K514,K866,K838,K28,K704,K947,K919,K190",
            ),
            (
                25,
                "K271,K595,K0,K299,K352,K676,K137,K433,K757,K218,K514,K838,K28,K919,K380,K461,\
                 K542,K785,K109,K56,K623,K866,K704,K947,K190",
            ),
        ];
        for (n, order) in expected {
            let keys = keys(n);
            let set = StringHashSet::copy_of(keys.iter().map(String::as_str));
            assert_eq!(set.iter().join(","), order, "{n} keys");
        }
        let base: Vec<String> = (0..10).map(|i| format!("A{i}")).collect();
        let mut grown = StringHashSet::copy_of(base.iter().map(String::as_str));
        for i in 0..20 {
            grown.insert(&format!("B{i}"));
        }
        assert_eq!(
            grown.iter().join(","),
            "A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B0,B1,B2,B3,\
             B4,B5,B6,B7,B8,B9"
        );
    }
}

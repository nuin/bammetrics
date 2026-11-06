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

/// The text Rust parses for a Java floating-point literal
/// (`Double.parseDouble`): the string trimmed, an optional sign, `NaN`,
/// `Infinity`, or digits with an optional `.` and exponent and an optional
/// `f`, `F`, `d` or `D` suffix. A hexadecimal literal (`0x1.8p1`) is
/// converted through a double, so a float parsed from one can round twice.
fn java_decimal(text: &str) -> Result<String, String> {
    let trimmed = trim(text);
    if trimmed.is_empty() {
        return Err("empty String".to_owned());
    }
    let (sign, body) = match trimmed.strip_prefix(['+', '-']) {
        Some(rest) => (&trimmed[..1], rest),
        None => ("", trimmed),
    };
    match body {
        "NaN" => return Ok("NaN".to_owned()),
        "Infinity" => return Ok(format!("{sign}inf")),
        _ => {}
    }
    let body = body.strip_suffix(['f', 'F', 'd', 'D']).unwrap_or(body);
    if let Some(hex) = body.strip_prefix("0x").or_else(|| body.strip_prefix("0X")) {
        let value = hex_value(hex).ok_or_else(|| bad_number(text))?;
        return Ok(format!("{sign}{value:e}"));
    }
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
    Ok(format!("{sign}{body}"))
}

/// `<hex digits>[.<hex digits>](p|P)[sign]<digits>`, after the `0x`.
fn hex_value(hex: &str) -> Option<f64> {
    let (mantissa, exponent) = hex.split_once(['p', 'P'])?;
    let unsigned = exponent.strip_prefix(['+', '-']).unwrap_or(exponent);
    if unsigned.is_empty() || !unsigned.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let exponent: i32 = exponent.parse().unwrap_or(if exponent.starts_with('-') {
        i32::MIN / 2
    } else {
        i32::MAX / 2
    });
    let (whole, fraction) = mantissa.split_once('.').unwrap_or((mantissa, ""));
    let digits = format!("{whole}{fraction}");
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_hexdigit()) {
        return None;
    }
    let significant = digits.trim_start_matches('0');
    let kept = &significant[..significant.len().min(30)];
    let dropped = (significant.len() - kept.len()) as i32;
    let value = u128::from_str_radix(if kept.is_empty() { "0" } else { kept }, 16).ok()? as f64;
    let shift = exponent
        .saturating_add(4 * dropped)
        .saturating_sub(4 * fraction.len() as i32);
    if value == 0.0 {
        return Some(0.0);
    }
    // Two steps keep 2^shift finite where the product is.
    let half = shift.clamp(-1100, 1100) / 2;
    Some(value * 2f64.powi(half) * 2f64.powi(shift.clamp(-1100, 1100) - half))
}

/// `Double.parseDouble`; the error is the `NumberFormatException` message.
pub fn parse_double(text: &str) -> Result<f64, String> {
    java_decimal(text)?.parse().map_err(|_| bad_number(text))
}

/// `Float.parseFloat`, rounded once from the decimal text.
pub fn parse_float(text: &str) -> Result<f32, String> {
    java_decimal(text)?.parse().map_err(|_| bad_number(text))
}

/// `Long.parseLong` (`bits` 64) or `Integer.parseInt` (32): an optional
/// sign, then ASCII digits, inside the type's range. Nothing is trimmed.
pub fn parse_integer(text: &str, bits: u32) -> Result<i64, String> {
    let body = text.strip_prefix(['+', '-']).unwrap_or(text);
    if body.is_empty() || !body.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad_number(text));
    }
    let value: i128 = text.parse().map_err(|_| bad_number(text))?;
    let limit = 1i128 << (bits - 1);
    if value < -limit || value >= limit {
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
/// as `d.dddE<exponent>` (`1.0E7`, `9.99E-4`). A subnormal value can differ
/// in its last digit (`Double.MIN_VALUE` is `5.0E-324` here, `4.9E-324` in
/// Java).
pub fn double_to_string(value: f64) -> String {
    match special_text(
        value.is_nan(),
        value.is_infinite(),
        value == 0.0,
        value.is_sign_negative(),
    ) {
        Some(text) => text,
        None => {
            let (digits, exponent) = crate::metrics_file::shortest_digits(value.abs());
            signed_layout(value.is_sign_negative(), &digits, exponent)
        }
    }
}

/// `Float.toString`: the same layout over the shortest digits of the float.
pub fn float_to_string(value: f32) -> String {
    match special_text(
        value.is_nan(),
        value.is_infinite(),
        value == 0.0,
        value.is_sign_negative(),
    ) {
        Some(text) => text,
        None => {
            let scientific = format!("{:e}", value.abs());
            let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
            let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
            signed_layout(
                value.is_sign_negative(),
                &digits,
                exponent.parse().unwrap_or(0),
            )
        }
    }
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

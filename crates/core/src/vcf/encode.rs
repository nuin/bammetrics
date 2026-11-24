//! A record line read as htsjdk 5.0.0's `VCFCodec` reads it and written
//! again as its `VCFEncoder` writes it (`VCFWriter.add`), one line at a
//! time.
//!
//! The line is first checked as [`super::site::parse_line`] checks it, so
//! a line the codec refuses is refused here with the codec's reason. Then:
//!
//! - `CHROM` and `ID` as given; `POS` as a Java `int` prints (`+10` and
//!   `010` print `10`).
//! - `REF` upper-cased; each `ALT` upper-cased unless symbolic, a `.`
//!   dropped, trailing empty ones dropped, `.` when none is left.
//! - `QUAL` from `log10PError` (`QUAL / -10`): missing prints `.`, any
//!   other value `String.format("%.2f")` of `log10PError * -10` with a
//!   trailing `.00` cut (`30.0` prints `30`, `12.345` prints `12.35`).
//!   Java rounds half up on the shortest decimal form of the double. A
//!   `QUAL` of `-0` prints `0`, as Picard 3.5.0 writes it.
//! - `FILTER`: `.` and `PASS` as given; names split on `;`, each once,
//!   sorted (`String.compareTo`) and joined by `;`. A name with no
//!   `##FILTER` line fails.
//! - `INFO`: `.` stays. Each `;` entry is `KEY=VALUE` or a key alone; a
//!   key alone is a flag when its `##INFO` line declares one or none, else
//!   it gets the value `.`; an empty value is `.`; `FLAG=0` for a declared
//!   flag is dropped; a later entry of a key replaces an earlier one. The
//!   keys print in sorted order, a key whose `Number` is 0 without its
//!   value. A key with no
//!   `##INFO` line fails. From VCF 4.3 each value is percent-decoded.
//! - `FORMAT` and the samples are written as they were read
//!   (`LazyGenotypesContext`), unless the header's sample names are not in
//!   sorted order: then the codec decodes them and the encoder writes them
//!   anew ([`super::encode_genotypes`]).

use std::collections::BTreeMap;

use super::full_header::{Count, VcfHeader, java_compare};
use super::site::{self, Rules, java_split, percent_decode, split_condensed};
use super::{Site, encode_genotypes};

/// One record line in, one record line out.
#[derive(Clone, Debug)]
pub struct Encoder {
    rules: Rules,
}

/// A key of the INFO map, ordered as Java's `TreeMap<String, _>` orders it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct JavaKey(String);

impl PartialOrd for JavaKey {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for JavaKey {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        java_compare(&self.0, &other.0)
    }
}

/// An INFO value as `parseInfo` holds it.
enum InfoValue {
    /// A key alone, read as a flag (`Boolean.TRUE`).
    Flag,
    Text(String),
}

/// `VCFEncoder.fieldIsMissingFromHeaderError`.
fn missing_from_header(key: &str, field: &str, site: &Site) -> String {
    format!(
        "Key {key} found in VariantContext field {field} at {}:{} but this key isn't defined in \
         the VCFHeader.  We require all VCFs to have complete VCF headers by default.",
        site.contig, site.start
    )
}

impl Encoder {
    pub fn new(header: &VcfHeader) -> Self {
        Self {
            rules: Rules::of(header.site()),
        }
    }

    /// Append the encoded form of `line` (no line end) to `out`; `None`
    /// for a line that starts with `#`, which the codec skips.
    pub fn encode(
        &self,
        header: &VcfHeader,
        line: &str,
        out: &mut String,
    ) -> Result<Option<Site>, String> {
        let Some((site, allele_count, genotype_text)) = site::parse_fixed(line, &self.rules)?
        else {
            return Ok(None);
        };
        // htsjdk decodes the genotypes as it reads the record, so their
        // errors come before the encoder's.
        let decoded = if self.rules.eager_genotypes {
            Some(encode_genotypes::decode(
                header,
                genotype_text,
                allele_count,
                &site,
            )?)
        } else {
            None
        };
        let max = if self.rules.samples > 0 { 9 } else { 8 };
        let parts = split_condensed(line, max);
        out.push_str(parts[0]);
        out.push('\t');
        out.push_str(&site.start.to_string());
        out.push('\t');
        out.push_str(parts[2]);
        out.push('\t');
        let reference = parts[3].to_ascii_uppercase();
        out.push_str(&reference);
        out.push('\t');
        let alternates = alternates(parts[4]);
        if alternates.is_empty() {
            out.push('.');
        } else {
            out.push_str(&alternates.join(","));
        }
        out.push('\t');
        match site::log10_p_error(parts[5])? {
            None => out.push('.'),
            // `-0.0 + 0.0` is `0.0`: a `QUAL` of `-0` prints `0`.
            Some(log10) => out.push_str(&format_quality(log10 * -10.0 + 0.0)),
        }
        out.push('\t');
        out.push_str(&filters(parts[6], header, &site)?);
        out.push('\t');
        self.info(parts[7], header, &site, out)?;
        if max == 9 {
            out.push('\t');
            match &decoded {
                Some(decoded) => encode_genotypes::write_decoded(header, decoded, &site, out)?,
                None => out.push_str(parts[8]),
            }
        }
        Ok(Some(site))
    }

    /// `parseInfo`, then `VCFEncoder.writeInfoString`.
    fn info(
        &self,
        text: &str,
        header: &VcfHeader,
        site: &Site,
        out: &mut String,
    ) -> Result<(), String> {
        if text == "." {
            out.push('.');
            return Ok(());
        }
        let decode = |value: &str| {
            if self.rules.percent_decoded {
                percent_decode(value)
            } else {
                value.to_owned()
            }
        };
        let mut map: BTreeMap<JavaKey, InfoValue> = BTreeMap::new();
        for entry in text.split(';') {
            let (key, value) = match entry.split_once('=') {
                Some((key, raw)) => {
                    let value = if raw.contains(',') {
                        raw.split(',').map(decode).collect::<Vec<_>>().join(",")
                    } else {
                        let value = decode(raw);
                        if value == "0" && header.info(key).is_some_and(|d| d.flag) {
                            continue;
                        }
                        value
                    };
                    let value = if value.is_empty() {
                        ".".to_owned()
                    } else {
                        value
                    };
                    (key, InfoValue::Text(value))
                }
                None => match header.info(entry) {
                    Some(definition) if !definition.flag => {
                        (entry, InfoValue::Text(".".to_owned()))
                    }
                    _ => (entry, InfoValue::Flag),
                },
            };
            map.insert(JavaKey(key.to_owned()), value);
        }
        if map.is_empty() {
            out.push('.');
            return Ok(());
        }
        for (i, (JavaKey(key), value)) in map.iter().enumerate() {
            let Some(definition) = header.info(key) else {
                return Err(missing_from_header(key, "INFO", site));
            };
            if i > 0 {
                out.push(';');
            }
            out.push_str(key);
            if let InfoValue::Text(value) = value
                && definition.count != Count::Fixed(0)
            {
                out.push('=');
                out.push_str(value);
            }
        }
        Ok(())
    }
}

/// The alternates `parseAlleles` keeps, as `Allele.getDisplayString`
/// prints them.
fn alternates(text: &str) -> Vec<String> {
    let alts = if text.contains(',') {
        java_split(text, ',')
    } else {
        vec![text]
    };
    alts.into_iter()
        .filter(|alt| *alt != ".")
        .map(|alt| {
            if site::is_symbolic(alt.as_bytes()) {
                alt.to_owned()
            } else {
                alt.to_ascii_uppercase()
            }
        })
        .collect()
}

/// `VCFEncoder.getFilterString` over the codec's filters.
fn filters(text: &str, header: &VcfHeader, site: &Site) -> Result<String, String> {
    if matches!(text, "." | "PASS") {
        return Ok(text.to_owned());
    }
    let mut names: Vec<&str> = java_split(text, ';');
    names.sort_by(|a, b| java_compare(a, b));
    names.dedup();
    if names.is_empty() {
        return Ok("PASS".to_owned());
    }
    if let Some(missing) = names.iter().find(|name| !header.has_filter(name)) {
        return Err(missing_from_header(missing, "FILTER", site));
    }
    Ok(names.join(";"))
}

/// `VCFEncoder.formatQualValue`: `String.format(Locale.US, "%.2f", q)`
/// with a trailing `.00` cut.
pub fn format_quality(quality: f64) -> String {
    let text = format_fixed2(quality);
    match text.strip_suffix(".00") {
        Some(cut) => cut.to_owned(),
        None => text,
    }
}

/// Java's `%.2f`: the shortest decimal digits of the double rounded half
/// up at the second decimal (`FormattedFloatingDecimal.applyPrecision`
/// looks at the first digit dropped).
pub fn format_fixed2(value: f64) -> String {
    if value.is_nan() {
        return "NaN".to_owned();
    }
    if value.is_infinite() {
        return if value > 0.0 { "Infinity" } else { "-Infinity" }.to_owned();
    }
    let scientific = format!("{:e}", value.abs());
    let (mantissa, exponent) = scientific.split_once('e').unwrap_or((&scientific, "0"));
    let exponent: i64 = exponent.parse().unwrap_or(0);
    let digits: Vec<u8> = mantissa
        .bytes()
        .filter(u8::is_ascii_digit)
        .map(|b| b - b'0')
        .collect();
    // The digits before the point of `value * 100`.
    let whole = exponent + 1 + 2;
    let mut kept: Vec<u8> = if whole <= 0 {
        Vec::new()
    } else {
        let whole = whole as usize;
        let mut kept: Vec<u8> = digits.iter().take(whole).copied().collect();
        kept.resize(whole, 0);
        kept
    };
    let next = if whole < 0 {
        None
    } else {
        digits.get(whole as usize).copied()
    };
    if next.is_some_and(|d| d >= 5) {
        let mut carry = true;
        for digit in kept.iter_mut().rev() {
            if *digit == 9 {
                *digit = 0;
            } else {
                *digit += 1;
                carry = false;
                break;
            }
        }
        if carry {
            kept.insert(0, 1);
        }
    }
    while kept.len() < 3 {
        kept.insert(0, 0);
    }
    while kept.len() > 3 && kept[0] == 0 {
        kept.remove(0);
    }
    let split = kept.len() - 2;
    let mut out = String::with_capacity(kept.len() + 2);
    if value.is_sign_negative() {
        out.push('-');
    }
    out.extend(kept[..split].iter().map(|d| char::from(b'0' + d)));
    out.push('.');
    out.extend(kept[split..].iter().map(|d| char::from(b'0' + d)));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qualities_print_as_java_formats_them() {
        let cases = [
            (30.0, "30"),
            (12.345, "12.35"),
            (0.001, "0"),
            (100.0, "100"),
            (1_000_000.5, "1000000.50"),
            (45.678, "45.68"),
            (30.5, "30.50"),
            (1.005, "1.01"),
            (2.675, "2.68"),
            (0.005, "0.01"),
            (0.004, "0"),
            (9.995, "10"),
            (99.999, "100"),
            (-0.0, "-0"),
            (f64::INFINITY, "Infinity"),
            (f64::NAN, "NaN"),
            (1e20, "100000000000000000000"),
        ];
        for (value, expected) in cases {
            assert_eq!(format_quality(value), expected, "{value}");
        }
    }

    fn header(meta: &str, samples: &str) -> VcfHeader {
        let text = format!(
            "##fileformat=VCFv4.2\n{meta}#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO{samples}\n"
        );
        VcfHeader::read(&mut text.as_bytes()).expect("a header")
    }

    fn encode(header: &VcfHeader, line: &str) -> Result<String, String> {
        let mut out = String::new();
        Encoder::new(header).encode(header, line, &mut out)?;
        Ok(out)
    }

    #[test]
    fn sites_print_as_the_encoder_prints_them() {
        let h = header(
            "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             ##INFO=<ID=FL,Number=0,Type=Flag,Description=\"f\">\n\
             ##INFO=<ID=AF,Number=A,Type=Float,Description=\"a\">\n\
             ##FILTER=<ID=q10,Description=\"q\">\n##FILTER=<ID=LowQual,Description=\"l\">\n",
            "\tFORMAT\tS1",
        );
        let cases = [
            (
                "c1\t+010\t.\ta\tc\t30\tPASS\tFL;DP=10;AF=0.5\tGT:GQ\t0/1:3.5",
                "c1\t10\t.\tA\tC\t30\tPASS\tAF=0.5;DP=10;FL\tGT:GQ\t0/1:3.5",
            ),
            (
                "c1\t11\trs1;rs2\tAC\tA,act,.\t30.0\tq10;LowQual;q10\tDP=10;DP=11;FL=1;AF=\tGT\t0",
                "c1\t11\trs1;rs2\tAC\tA,ACT\t30\tLowQual;q10\tAF=.;DP=11;FL\tGT\t0",
            ),
            (
                "c1\t12\t.\tG\t.\t-1\t.\tDP;FL=0\tGT\t.",
                "c1\t12\t.\tG\t.\t.\t.\tDP=.\tGT\t.",
            ),
            (
                "c1\t13\t.\tG\t<del>\t12.345\t.\t.\tGT\t.",
                "c1\t13\t.\tG\t<del>\t12.35\t.\t.\tGT\t.",
            ),
            (
                "c1\t14\t.\tG\tC\t-0\t.\t.\tGT\t.",
                "c1\t14\t.\tG\tC\t0\t.\t.\tGT\t.",
            ),
            (
                "c1\t15\t.\tG\tC\t-0.0\t.\t.\tGT\t.",
                "c1\t15\t.\tG\tC\t0\t.\t.\tGT\t.",
            ),
        ];
        for (line, expected) in cases {
            assert_eq!(encode(&h, line).as_deref(), Ok(expected), "{line}");
        }
        for bad in [
            "c1\t1\t.\tA\tC\t.\tother\t.\tGT\t0",
            "c1\t1\t.\tA\tC\t.\t.\tXX=1\tGT\t0",
        ] {
            assert!(encode(&h, bad).is_err(), "{bad}");
        }
    }
}

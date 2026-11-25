//! A record line read for its site as htsjdk 5.0.0's `VCFCodec.decode`
//! reads it (`AbstractVCFCodec.parseVCFLine`), the genotypes left undecoded
//! as `LazyGenotypesContext` leaves them.
//!
//! The rules, in the codec's order:
//!
//! - A line starting with `#` is skipped. The line is cut at tabs into 8
//!   fields (9 when the header names samples), the last holding the rest
//!   (`ParsingUtils.split` with `condenseTrailingTokens`); a leading tab is
//!   dropped. Any other count is an error, so an empty line fails.
//! - `POS` is a Java `int` (`Integer.parseInt`: `+10` and `-5` are read).
//! - `ID` may not be empty; `.` is no ID, any other text is kept whole.
//! - `QUAL` is `.` or a Java double, or `inf`, `infinity` or `nan` in any
//!   case with a sign (`VCFUtils.parseVcfDouble`). A value within 5e-5 of
//!   -1 is missing; any other negative value but -10 fails
//!   `VariantContext`'s check that `log10PError` is not above 0.
//! - `FILTER`: `.` and `PASS` are unfiltered, `0` and an empty field are
//!   errors, else the names split on `;` (trailing empty names dropped)
//!   must each match `^[!-:<-~]+$`. A record with names is filtered.
//! - `INFO`: not empty, no tab or space; `.` is no attribute. Each
//!   `;`-separated entry is `KEY=VALUE` or a flag; a later `END` replaces
//!   an earlier one. `END=0` under an `END` declared as `Flag` is dropped.
//!   An `END` whose value is not one Java `int` (a list, a flag, `.`) is
//!   an error. In VCF 4.3 and later the value is percent-decoded first.
//!   Without `END` the end is `POS` plus the reference length less 1 in
//!   Java `int` arithmetic.
//! - `REF` is upper-cased and checked as `checkAllele` checks it: not
//!   empty, not symbolic, not `.`, not `*`, no leading `D` or `I`, only
//!   `ACGTN`. Each `ALT` (the field split on `,`, trailing empty ones
//!   dropped) is not empty, no leading `D` or `I`, symbolic or `ACGTN`
//!   in any case, `*` or `.` (a `.` is dropped). Two alleles with the same
//!   bases are an error (`Duplicate allele added to VariantContext`).
//! - Without `END` and without a symbolic allele, an end that wrapped past
//!   2^31 fails `validateStop`.
//! - When the header's sample names are not in sorted order, htsjdk
//!   decodes the genotypes at once (`samplesWereAlreadySorted`), so their
//!   errors are the record's.

use std::collections::HashSet;

use super::Site;
use super::site_header::SiteHeader;
use crate::error::excerpt;

/// What the header tells the line parser.
#[derive(Clone, Debug, Default)]
pub struct Rules {
    /// The header names samples: a record has 9 fields, `FORMAT` and the
    /// samples in the last.
    pub samples: usize,
    /// The sample names are not in sorted order: genotypes are decoded.
    pub eager_genotypes: bool,
    /// `##INFO=<ID=END,...,Type=Flag>`.
    pub end_is_flag: bool,
    /// VCF 4.3 or later: `INFO` values are percent-decoded
    /// (`VCFPercentEncodedTextTransformer`).
    pub percent_decoded: bool,
}

impl Rules {
    /// What [`VcfReader::next_site`] needs from the header: the sample count,
    /// whether the samples are in sorted order (`String.compareTo`), and
    /// whether `END` is declared as a flag.
    pub fn of(header: &SiteHeader) -> Self {
        let names = &header.samples;
        let sorted = names
            .windows(2)
            .all(|w| w[0].encode_utf16().le(w[1].encode_utf16()));
        Self {
            samples: names.len(),
            eager_genotypes: !sorted,
            end_is_flag: header.end_is_flag,
            percent_decoded: header.minor >= 3,
        }
    }
}

/// The parsed line: `None` for a line starting with `#`.
pub fn parse_line(line: &str, rules: &Rules) -> Result<Option<Site>, String> {
    let Some((site, alleles, genotypes)) = parse_fixed(line, rules)? else {
        return Ok(None);
    };
    if rules.eager_genotypes {
        super::encode_genotypes::check(genotypes, rules.samples, alleles, rules.percent_decoded)?;
    }
    Ok(Some(site))
}

/// [`parse_line`] without the genotypes decoded: the site, the record's
/// allele count and the `FORMAT` and sample text (empty without samples),
/// for a caller that decodes the genotypes itself.
pub(super) fn parse_fixed<'l>(
    line: &'l str,
    rules: &Rules,
) -> Result<Option<(Site, usize, &'l str)>, String> {
    if line.starts_with('#') {
        return Ok(None);
    }
    let max = if rules.samples > 0 { 9 } else { 8 };
    let parts = split_condensed(line, max);
    if parts.len() != max {
        return Err(format!(
            "there aren't enough columns for line {} (we expected {max} tokens, and saw {} )",
            excerpt(line),
            parts.len()
        ));
    }
    let start = crate::java::parse_integer(parts[1], 32).map_err(|_| {
        format!(
            "{} is not a valid start position in the VCF format",
            excerpt(parts[1])
        )
    })? as i32;
    let id = match parts[2] {
        "" => return Err("The VCF specification requires a valid ID field".to_owned()),
        "." => None,
        id => Some(id.to_owned()),
    };
    let reference = parts[3].to_uppercase();
    check_quality(parts[5])?;
    let filtered = parse_filters(parts[6])?;
    let end_value = parse_info_end(parts[7], rules)?;
    let end = match end_value {
        Some(text) => Some(
            crate::java::parse_integer(&text, 32)
                .map_err(|_| "the END value in the INFO field is not valid".to_owned())?
                as i32,
        ),
        None => None,
    };
    let alleles = parse_alleles(&reference, parts[4])?;
    let end = match end {
        Some(end) => end,
        None => {
            let wide = i64::from(start) + reference.encode_utf16().count() as i64 - 1;
            if wide > i64::from(i32::MAX) && !alleles.symbolic {
                return Err(format!(
                    "BUG: GenomeLoc {}:{start}-{} has a size == {} but the variation reference \
                     allele has length {}",
                    parts[0],
                    wide as i32,
                    i64::from(wide as i32) - i64::from(start) + 1,
                    reference.len()
                ));
            }
            wide as i32
        }
    };
    let site = Site {
        contig: parts[0].to_owned(),
        start,
        end,
        id,
        filtered,
    };
    Ok(Some((
        site,
        alleles.count,
        parts.get(8).copied().unwrap_or(""),
    )))
}

/// `ParsingUtils.split(line, parts, '\t', true)` with `max` parts: a
/// leading tab skipped, the last part holding the rest of the line.
pub(super) fn split_condensed(line: &str, max: usize) -> Vec<&str> {
    let mut parts = Vec::with_capacity(max);
    let mut rest = line;
    if rest.starts_with('\t') {
        if rest.len() == 1 {
            return parts;
        }
        rest = &rest[1..];
    }
    while parts.len() < max {
        match rest.find('\t') {
            Some(tab) if parts.len() + 1 < max => {
                parts.push(&rest[..tab]);
                rest = &rest[tab + 1..];
            }
            _ => {
                parts.push(rest);
                break;
            }
        }
    }
    parts
}

/// `AbstractVCFCodec.parseQual`, then `CommonInfo.setLog10PError`.
fn check_quality(text: &str) -> Result<(), String> {
    log10_p_error(text).map(|_| ())
}

/// `VCFConstants.VCF_ENCODING_EPSILON`.
const ENCODING_EPSILON: f64 = 5e-5;

/// `AbstractVCFCodec.parseQual`, then `CommonInfo.setLog10PError`: the
/// `log10PError` of `QUAL`, `None` for a missing quality. `.` is missing,
/// and so is a negative value within 5e-5 of -1 (the VCF 3 spelling). Any
/// other value is divided by -10; a result above 0 other than exactly 1
/// fails (the message prints the builder's value before the set, 1.0), and
/// exactly 1 (`QUAL` -10) is `NO_LOG10_PERROR`, missing too.
pub(crate) fn log10_p_error(text: &str) -> Result<Option<f64>, String> {
    if text == "." {
        return Ok(None);
    }
    let quality = match crate::java::parse_double(text) {
        Ok(q) => q,
        Err(message) => special_double(text).ok_or(message)?,
    };
    if quality < 0.0 && (quality + 1.0).abs() < ENCODING_EPSILON {
        return Ok(None);
    }
    let log10 = quality / -10.0;
    if log10 > 0.0 && log10 != 1.0 {
        return Err("BUG: log10PError cannot be > 0 : 1.0".to_owned());
    }
    Ok((log10 != 1.0).then_some(log10))
}

/// `VCFUtils.parseVcfDouble`'s `^[-+]?(INF|INFINITY|NAN)$` without regard
/// to case.
pub(super) fn special_double(text: &str) -> Option<f64> {
    let (negative, body) = match text.as_bytes().first() {
        Some(b'-') => (true, &text[1..]),
        Some(b'+') => (false, &text[1..]),
        _ => (false, text),
    };
    let value = if body.eq_ignore_ascii_case("nan") {
        f64::NAN
    } else if body.eq_ignore_ascii_case("inf") || body.eq_ignore_ascii_case("infinity") {
        f64::INFINITY
    } else {
        return None;
    };
    Some(if negative { -value } else { value })
}

/// `VCFCodec.parseFilters` and the filter check of `VariantContext`: true
/// when the record is filtered.
fn parse_filters(text: &str) -> Result<bool, String> {
    match text {
        "." | "PASS" => return Ok(false),
        "0" => return Err("0 is an invalid filter name in vcf4".to_owned()),
        "" => {
            return Err(format!(
                "The VCF specification requires a valid filter status: filter was {text}"
            ));
        }
        _ => {}
    }
    let names = java_split(text, ';');
    for name in &names {
        if name.is_empty() || !name.bytes().all(|b| matches!(b, b'!'..=b':' | b'<'..=b'~')) {
            return Err(format!(
                "Filter '{}' contains an illegal character. It must conform to the regex \
                 ;'^[!-:<-~]+$",
                excerpt(name)
            ));
        }
    }
    Ok(!names.is_empty())
}

/// `VariantContext.determineType` from the `REF` and `ALT` text of a line
/// the codec read: no alternate is `NO_VARIATION`; each alternate against
/// the reference is `SYMBOLIC`, a `SNP` or `MNP` (same length) or an
/// `INDEL`; two alternates of different types make `MIXED`.
pub(super) fn kind_of(reference: &str, alternates: &str) -> super::VariantType {
    use super::VariantType;
    let alts = if alternates.contains(',') {
        java_split(alternates, ',')
    } else {
        vec![alternates]
    };
    let mut kind = None;
    for alt in alts.into_iter().filter(|alt| *alt != ".") {
        let this = if is_symbolic(alt.as_bytes()) {
            VariantType::Symbolic
        } else if alt.len() == reference.len() {
            if alt.len() == 1 {
                VariantType::Snp
            } else {
                VariantType::Mnp
            }
        } else {
            VariantType::Indel
        };
        match kind {
            None => kind = Some(this),
            Some(k) if k != this => return VariantType::Mixed,
            Some(_) => {}
        }
    }
    kind.unwrap_or(VariantType::NoVariation)
}

/// `String.split` on one character: trailing empty parts dropped.
pub(super) fn java_split(text: &str, separator: char) -> Vec<&str> {
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// `AbstractVCFCodec.parseInfo` as far as `END` goes: the text of the
/// last `END` entry as `Object.toString` gives it (`true` for a flag, a
/// list in brackets), `None` without one.
fn parse_info_end(info: &str, rules: &Rules) -> Result<Option<String>, String> {
    let end_is_flag = rules.end_is_flag;
    if info.is_empty() {
        return Err("The VCF specification requires a valid (non-zero length) info field".into());
    }
    if info == "." {
        return Ok(None);
    }
    if info.contains(['\t', ' ']) {
        return Err(format!(
            "The VCF specification does not allow for whitespace in the INFO field. Offending \
             field value was \"{}\"",
            excerpt(info)
        ));
    }
    let mut end = None;
    for entry in info.split(';') {
        match entry.split_once('=') {
            Some(("END", value)) => {
                let values: Vec<&str> = value.split(',').collect();
                if values.len() == 1 {
                    let value = if rules.percent_decoded {
                        percent_decode(value)
                    } else {
                        value.to_owned()
                    };
                    if end_is_flag && value == "0" {
                        continue;
                    }
                    end = Some(if value.is_empty() {
                        ".".to_owned()
                    } else {
                        value
                    });
                } else {
                    end = Some(format!("[{}]", values.join(", ")));
                }
            }
            None if entry == "END" => {
                end = Some(if end_is_flag { "true" } else { "." }.to_owned());
            }
            _ => {}
        }
    }
    Ok(end)
}

/// `VCFPercentEncodedTextTransformer.decodeText`: each `%` with two hex
/// digits after it becomes the character of that code; any other `%`
/// stays.
pub(super) fn percent_decode(text: &str) -> String {
    if !text.contains('%') {
        return text.to_owned();
    }
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    let mut i = 0;
    while i < chars.len() {
        let code = (chars[i] == '%' && i + 2 < chars.len())
            .then(|| chars[i + 1].to_digit(16).zip(chars[i + 2].to_digit(16)))
            .flatten();
        match code {
            Some((high, low)) => {
                out.push(char::from((high * 16 + low) as u8));
                i += 3;
            }
            None => {
                out.push(chars[i]);
                i += 1;
            }
        }
    }
    out
}

/// The alleles a record holds after `parseAlleles`.
struct Alleles {
    /// The reference and the alternates other than `.`.
    count: usize,
    symbolic: bool,
}

/// `Allele.wouldBeSymbolicAllele`.
pub(super) fn is_symbolic(bytes: &[u8]) -> bool {
    bytes.len() > 1
        && (bytes[0] == b'<'
            || bytes[bytes.len() - 1] == b'>'
            || bytes.contains(&b'[')
            || bytes.contains(&b']')
            || bytes[0] == b'.'
            || bytes[bytes.len() - 1] == b'.')
}

/// `AbstractVCFCodec.checkAllele` and `Allele.acceptableAlleleBases`.
fn check_allele(allele: &str, is_ref: bool) -> Result<(), String> {
    let bytes = allele.as_bytes();
    if bytes.is_empty() {
        return Err("empty alleles are not permitted in VCF records".to_owned());
    }
    if is_symbolic(bytes) {
        if is_ref {
            return Err(format!(
                "Symbolic alleles not allowed as reference allele: {}",
                excerpt(allele)
            ));
        }
        return Ok(());
    }
    if matches!(bytes[0], b'D' | b'I') {
        return Err(
            "Insertions/Deletions are not supported when reading 3.x VCF's. Please convert your \
             file to VCF4 using VCFTools, available at http://vcftools.sourceforge.net/index.html"
                .to_owned(),
        );
    }
    let acceptable = match bytes {
        [b'-'] => false,
        [b'.'] => true,
        [b'*'] => !is_ref,
        _ => bytes
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N')),
    };
    if !acceptable {
        return Err(if allele.contains(['[', ']', ':', '.']) {
            "VCF support for complex rearrangements with breakends has not yet been implemented"
                .to_owned()
        } else {
            format!("unparsable vcf record with allele {}", excerpt(allele))
        });
    }
    if is_ref && allele == "." {
        return Err("The reference allele cannot be missing".to_owned());
    }
    Ok(())
}

/// `AbstractVCFCodec.parseAlleles`, then the duplicate check of
/// `VariantContext.makeAlleles` (bases compared, the reference flag
/// ignored), one hash lookup per allele.
fn parse_alleles(reference: &str, alternates: &str) -> Result<Alleles, String> {
    check_allele(reference, true)?;
    let mut seen: HashSet<String> = HashSet::from([reference.to_owned()]);
    let mut symbolic = false;
    let alts = if alternates.contains(',') {
        java_split(alternates, ',')
    } else {
        vec![alternates]
    };
    for alt in alts {
        check_allele(alt, false)?;
        if alt == "." {
            continue;
        }
        let bases = if is_symbolic(alt.as_bytes()) {
            symbolic = true;
            alt.to_owned()
        } else {
            alt.to_ascii_uppercase()
        };
        if seen.contains(&bases) {
            return Err(format!(
                "Duplicate allele added to VariantContext: {}",
                excerpt(&bases)
            ));
        }
        seen.insert(bases);
    }
    Ok(Alleles {
        count: seen.len(),
        symbolic,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(line: &str) -> Result<Option<Site>, String> {
        parse_line(line, &Rules::default())
    }

    fn span(line: &str) -> (i32, i32) {
        let s = site(line).expect("parses").expect("a site");
        (s.start, s.end)
    }

    #[test]
    fn fields_split_as_parsing_utils_splits_them() {
        assert_eq!(split_condensed("a\tb\tc", 2), ["a", "b\tc"]);
        assert_eq!(split_condensed("\ta\tb", 8), ["a", "b"]);
        assert_eq!(split_condensed("", 8), [""]);
        assert!(split_condensed("\t", 8).is_empty());
        assert_eq!(split_condensed("a\t\tb", 3), ["a", "", "b"]);
    }

    #[test]
    fn ends_come_from_end_or_the_reference() {
        assert_eq!(span("c1\t10\t.\tACGT\tA\t.\t.\t."), (10, 13));
        assert_eq!(span("c1\t10\t.\tacgt\tA\t.\t.\t."), (10, 13));
        assert_eq!(span("c1\t10\t.\tA\t<DEL>\t.\t.\tEND=5"), (10, 5));
        assert_eq!(span("c1\t10\t.\tA\t<DEL>\t.\t.\tEND=20;END=30"), (10, 30));
        assert_eq!(span("c1\t+10\t.\tA\t<DEL>\t.\t.\tEND=+20"), (10, 20));
        assert_eq!(span("c1\t-5\t.\tA\tC\t.\t.\t."), (-5, -5));
        for bad in ["END=.", "END", "END=20,30", "END= 20", "END=abc", "END="] {
            assert!(
                site(&format!("c1\t10\t.\tA\t<DEL>\t.\t.\t{bad}")).is_err(),
                "{bad}"
            );
        }
        let flag = Rules {
            end_is_flag: true,
            ..Rules::default()
        };
        let s = parse_line("c1\t10\t.\tA\tC\t.\t.\tEND=0", &flag)
            .expect("ok")
            .expect("site");
        assert_eq!(s.end, 10);
        let v43 = Rules {
            percent_decoded: true,
            ..Rules::default()
        };
        let s = parse_line("c1\t1\t.\tA\t<DEL>\t.\t.\tEND=%32%30", &v43)
            .expect("ok")
            .expect("site");
        assert_eq!(s.end, 20);
        assert!(site("c1\t1\t.\tA\t<DEL>\t.\t.\tEND=%32%30").is_err());
        assert_eq!(percent_decode("a%3Db%2"), "a=b%2");
    }

    #[test]
    fn an_end_past_two_to_the_31_fails_without_a_symbolic_allele() {
        assert!(site("c1\t2147483647\t.\tAC\tC\t.\t.\t.").is_err());
        assert_eq!(
            span("c1\t2147483647\t.\tAC\t<DEL>\t.\t.\t."),
            (i32::MAX, i32::MIN)
        );
    }

    #[test]
    fn filters_follow_the_codec() {
        let filtered =
            |f: &str| site(&format!("c1\t1\t.\tA\tC\t.\t{f}\t.")).map(|s| s.unwrap().filtered);
        assert_eq!(filtered("."), Ok(false));
        assert_eq!(filtered("PASS"), Ok(false));
        assert_eq!(filtered("q10"), Ok(true));
        assert_eq!(filtered("q10;PASS"), Ok(true));
        assert_eq!(filtered("q10;"), Ok(true));
        assert!(filtered(";q10").is_err());
        assert!(filtered("").is_err());
        assert!(filtered("0").is_err());
        assert!(filtered("a b").is_err());
    }

    #[test]
    fn qualities_follow_the_codec() {
        let ok = |q: &str| site(&format!("c1\t1\t.\tA\tC\t{q}\t.\t.")).is_ok();
        for good in [
            ".", "5", " 5 ", "5d", "0x1p3", "-10", "nan", "Infinity", "INF", "+inf", "-0", "-1",
            "-1.00001", "-0.99996",
        ] {
            assert!(ok(good), "{good}");
        }
        for bad in [
            "-5",
            "abc",
            "-inf",
            "-Infinity",
            "",
            "-0.5",
            "-10.000001",
            "-1.0001",
        ] {
            assert!(!ok(bad), "{bad}");
        }
    }

    #[test]
    fn alleles_follow_the_codec() {
        let ok = |r: &str, a: &str| site(&format!("c1\t1\t.\t{r}\t{a}\t.\t.\t.")).is_ok();
        for (r, a) in [
            ("A", "*"),
            ("A", "."),
            ("A", "C,."),
            ("A", "C,"),
            ("A", "A[c1:20["),
            ("n", "<X>"),
        ] {
            assert!(ok(r, a), "{r} {a}");
        }
        for (r, a) in [
            ("A", "DA"),
            ("A", "C,c"),
            ("A", "A"),
            ("*", "C"),
            ("<X>", "C"),
            (".", "C"),
            ("AXZ", "C"),
            ("A", "Q"),
            ("A", ",C"),
            ("", "C"),
        ] {
            assert!(!ok(r, a), "{r} {a}");
        }
    }

    /// 100,000 alternates, then the first again: one hash lookup each, so
    /// the record is read in well under a second (a scan of the alleles
    /// seen took 10 s here).
    #[test]
    fn many_alternates_are_checked_in_linear_time() {
        let alternate = |i: usize| -> String {
            (0..9)
                .map(|k| ['A', 'C', 'G', 'T'][(i >> (2 * k)) & 3])
                .collect()
        };
        let alts: Vec<String> = (0..100_000).map(alternate).collect();
        let started = std::time::Instant::now();
        let line = |alts: &str| format!("c1\t1\t.\tAAAAAAAAAA\t{alts}\t.\t.\t.");
        let joined = alts.join(",");
        assert!(site(&line(&joined)).is_ok());
        assert!(site(&line(&format!("{joined},{}", alts[0]))).is_err());
        assert!(started.elapsed() < std::time::Duration::from_secs(5));
    }

    #[test]
    fn records_need_their_columns() {
        assert!(site("").is_err());
        assert!(site("c1\t10\t.\tA\tC\t.\t.").is_err());
        assert!(site("c1\t10\t.\tA\tC\t.\t.\t.\tGT").is_err());
        assert!(site("c1\t10\t\tA\tC\t.\t.\t.").is_err());
        assert!(site("c1\t10\t.\tA\tC\t.\t.\t").is_err());
        assert!(site("#c1\t10").expect("skipped").is_none());
        let with_samples = Rules {
            samples: 2,
            ..Rules::default()
        };
        assert!(parse_line("c1\t10\t.\tA\tC\t.\t.\t.", &with_samples).is_err());
        assert!(parse_line("c1\t10\t.\tA\tC\t.\t.\t.\tGT\tzz", &with_samples).is_ok());
        let eager = Rules {
            eager_genotypes: true,
            ..with_samples
        };
        assert!(parse_line("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1", &eager).is_err());
        assert!(parse_line("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/5\t0/1", &eager).is_err());
        assert!(parse_line("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1|1", &eager).is_ok());
    }
}

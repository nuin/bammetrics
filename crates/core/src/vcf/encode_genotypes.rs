//! The genotypes of a record decoded as htsjdk 5.0.0's
//! `AbstractVCFCodec.createGenotypeMap` decodes them and written again as
//! `VCFEncoder.appendGenotypeData` writes them. htsjdk takes this path when
//! the header's sample names are not in sorted order; with sorted names
//! the text goes through unread ([`super::encode`]).
//!
//! Decoding, per sample: the values split on `:` (no more than the keys);
//! `GT` first when present (`/`, `|` and `\` separate allele indexes, empty
//! tokens skipped, `.` a no-call; any `|` makes the genotype phased); `FT`
//! `.` and `PASS` unfiltered, else the names sorted and joined by `;`; any
//! other `.` value absent; `GQ` `Math.round` of the number, absent when
//! that is `-1`; `DP` an `int`, `-1` absent; `AD` and `PL` a list of `int`
//! (absent when one is not); `GL` turned into `PL`
//! (`GenotypeLikelihoods.fromGLField(...).getAsPLs()`) until a `PL` value
//! of any sample of the record was read, after which `GL` is skipped;
//! every other key kept as its text. From VCF 4.3 each value is
//! percent-decoded.
//!
//! Encoding: the keys are `GT` when a genotype has alleles, then in sorted
//! order every key a genotype holds (`GQ`, `DP`, `AD`, `PL` when one has a
//! value, `FT` when one is filtered). A key without a `##FORMAT` line
//! fails. A genotype prints `PASS` for an unfiltered `FT`, a `GQ` over 99
//! as 99, `.` for an absent value whatever the key's `Number`, and drops
//! its trailing missing values.

use std::borrow::Cow;
use std::collections::HashSet;

use super::Site;
use super::full_header::{VcfHeader, java_compare};
use super::site::percent_decode;

/// A run of values in one of the arrays of [`Genotypes`].
#[derive(Clone, Copy, Debug, Default)]
struct Span {
    start: u32,
    end: u32,
}

impl Span {
    fn of<T>(values: &[T], start: usize) -> Self {
        Self {
            start: start as u32,
            end: values.len() as u32,
        }
    }

    fn in_<T>(self, values: &[T]) -> &[T] {
        &values[self.start as usize..self.end as usize]
    }
}

/// One sample as the codec decodes it. Its alleles and integer lists are
/// spans of the record's arrays, so a sample costs no allocation of its
/// own unless it holds a key outside the standard ones.
#[derive(Default)]
struct Decoded<'t> {
    /// Allele indexes, `None` for a no-call; empty when there is no `GT`.
    alleles: Span,
    phased: bool,
    gq: Option<i32>,
    dp: Option<i32>,
    ad: Option<Span>,
    pl: Option<Span>,
    /// The filters joined by `;`, `None` when unfiltered.
    filters: Option<String>,
    /// Every other key with its value, in `FORMAT` order.
    attributes: Vec<(&'t str, Cow<'t, str>)>,
}

impl Decoded<'_> {
    fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .rev()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.as_ref())
    }
}

/// The genotypes of one record: the keys of `FORMAT`, each sample, and
/// the arrays their spans point into.
#[derive(Default)]
struct Genotypes<'t> {
    keys: Vec<&'t str>,
    samples: Vec<Decoded<'t>>,
    alleles: Vec<Option<usize>>,
    ints: Vec<i32>,
}

/// `Math.round(double)`: the nearest `long`, a half rounded up, NaN 0,
/// saturated at the `long` range.
fn java_round(value: f64) -> i64 {
    if value.is_nan() {
        return 0;
    }
    let floor = value.floor();
    // `value - floor` is exact for every finite double.
    let rounded = if value - floor >= 0.5 {
        floor + 1.0
    } else {
        floor
    };
    if rounded >= 9.223_372_036_854_776e18 {
        i64::MAX
    } else if rounded <= -9.223_372_036_854_776e18 {
        i64::MIN
    } else {
        rounded as i64
    }
}

/// `Math.round(double)` cast to an `int`, as `GenotypeBuilder.GQ` takes it.
fn java_round_int(value: f64) -> i32 {
    java_round(value) as i32
}

/// `VCFUtils.parseVcfDouble`.
fn vcf_double(text: &str) -> Result<f64, String> {
    crate::java::parse_double(text).or_else(|e| super::site::special_double(text).ok_or(e))
}

/// `String.split(",")`: no comma gives the text whole, trailing empty
/// parts are dropped.
fn java_split_commas(text: &str) -> Vec<&str> {
    if !text.contains(',') {
        return vec![text];
    }
    let mut parts: Vec<&str> = text.split(',').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `Math.max(double, double)` and `Math.min`: NaN when either is.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// `GenotypeLikelihoods.fromGLField(text).getAsPLs()`: `None` when every
/// value is `.`, an error when some are. Each PL is
/// `(int) Math.round(Math.min(-10 * (GL - max GL), 2^31 - 1))`.
fn gl_to_pl(text: &str, ints: &mut Vec<i32>) -> Result<Option<Span>, String> {
    if text == "." {
        return Ok(None);
    }
    let parts = java_split_commas(text);
    let mut missing = 0usize;
    let mut likelihoods = Vec::with_capacity(parts.len());
    for part in &parts {
        if *part == "." {
            missing += 1;
        } else {
            likelihoods.push(vcf_double(part)?);
        }
    }
    if missing > 0 {
        return if missing == parts.len() {
            Ok(None)
        } else {
            Err("partial missing values for GL field".to_owned())
        };
    }
    let max = likelihoods
        .iter()
        .fold(f64::NEG_INFINITY, |max, &gl| java_max(max, gl));
    let start = ints.len();
    ints.extend(
        likelihoods
            .iter()
            .map(|&gl| java_round(java_min(-10.0 * (gl - max), 2_147_483_647.0)) as i32),
    );
    Ok(Some(Span::of(ints, start)))
}

/// `decodeInts` into `ints`: `None` (and `ints` as it was) when a value is
/// no Java `int`.
fn decode_ints(text: &str, ints: &mut Vec<i32>) -> Option<Span> {
    let start = ints.len();
    for value in text.split(',') {
        match crate::java::parse_integer(value, 32) {
            Ok(n) => ints.push(n as i32),
            Err(_) => {
                ints.truncate(start);
                return None;
            }
        }
    }
    Some(Span::of(ints, start))
}

/// `parseFilters` for `FT`: the filters, `None` for `.` and `PASS`.
fn genotype_filters(text: &str) -> Result<Option<String>, String> {
    match text {
        "." | "PASS" => Ok(None),
        "0" => Err("0 is an invalid filter name in vcf4".to_owned()),
        "" => Err(format!(
            "The VCF specification requires a valid filter status: filter was {text}"
        )),
        _ => {
            let mut names: Vec<&str> = if text.contains(';') {
                super::site::java_split(text, ';')
            } else {
                vec![text]
            };
            if names.is_empty() {
                return Ok(None);
            }
            names.sort_by(|a, b| java_compare(a, b));
            Ok(Some(names.join(";")))
        }
    }
}

/// `parseGenotypeAlleles` into `out`: `StringTokenizer` over `/|\`.
fn genotype_alleles(
    text: &str,
    alleles: usize,
    out: &mut Vec<Option<usize>>,
) -> Result<Span, String> {
    let start = out.len();
    for token in text
        .split(['/', '|', '\\'])
        .filter(|token| !token.is_empty())
    {
        if token == "." {
            out.push(None);
            continue;
        }
        let index = crate::java::parse_integer(token, 32).map_err(|_| {
            format!("The following invalid GT allele index was encountered in the file: {token}")
        })?;
        match usize::try_from(index) {
            Ok(i) if i < alleles => out.push(Some(i)),
            _ => {
                return Err(format!(
                    "The allele with index {token} is not defined in the REF/ALT columns in the \
                     record"
                ));
            }
        }
    }
    Ok(Span::of(out, start))
}

/// What [`decode_sample`] reads one sample with.
struct Reading<'k, 't> {
    keys: &'k [&'t str],
    alleles: usize,
    percent: bool,
    /// The codec's flag for the whole record: set once any sample's `PL`
    /// was read, after which `GL` is skipped.
    pl_read: bool,
    /// The values of the sample being read, reused from one to the next.
    values: Vec<Cow<'t, str>>,
}

/// One sample, its alleles and integer lists pushed onto `genotypes`.
fn decode_sample<'t>(
    reading: &mut Reading<'_, 't>,
    text: &'t str,
    genotypes: &mut Genotypes<'t>,
) -> Result<Decoded<'t>, String> {
    let keys = reading.keys;
    let mut values = std::mem::take(&mut reading.values);
    values.clear();
    values.extend(text.split(':').map(|v| {
        if reading.percent {
            Cow::Owned(percent_decode(v))
        } else {
            Cow::Borrowed(v)
        }
    }));
    if keys.len() < values.len() {
        return Err(format!(
            "There are too many keys for the sample, keys = {}, values = {text}",
            keys.join(":")
        ));
    }
    let mut decoded = Decoded::default();
    let mut gt_at = None;
    for (i, key) in keys.iter().enumerate() {
        if *key == "GT" {
            gt_at = Some(i);
            continue;
        }
        let Some(value) = values.get(i) else {
            continue;
        };
        if *key == "FT" {
            decoded.filters = genotype_filters(value)?;
            continue;
        }
        if value == "." {
            continue;
        }
        match *key {
            "GQ" => {
                // `GenotypeBuilder.GQ(-1)` is no GQ.
                let q = java_round_int(vcf_double(value)?);
                decoded.gq = (q != -1).then_some(q);
            }
            "DP" => {
                let dp = crate::java::parse_integer(value, 32)?;
                decoded.dp = (dp != -1).then_some(dp as i32);
            }
            "AD" => decoded.ad = decode_ints(value, &mut genotypes.ints),
            "PL" => {
                decoded.pl = decode_ints(value, &mut genotypes.ints);
                reading.pl_read = true;
            }
            "GL" => {
                if !reading.pl_read {
                    decoded.pl = gl_to_pl(value, &mut genotypes.ints)?;
                }
            }
            _ => decoded.attributes.push((*key, value.clone())),
        }
    }
    if let Some(at) = gt_at {
        if at > 0 {
            return Err(format!(
                "Saw GT field at position {at}, but it must be at the first position for \
                 genotypes when present"
            ));
        }
        let gt = values.get(at).map_or("", |v| v.as_ref());
        decoded.alleles = genotype_alleles(gt, reading.alleles, &mut genotypes.alleles)?;
        decoded.phased = gt.contains('|');
    }
    reading.values = values;
    Ok(decoded)
}

/// `isMissingValue`: only `.` and `,`.
fn is_missing(value: &str) -> bool {
    value.bytes().all(|b| b == b'.' || b == b',')
}

/// `createGenotypeMap` over `text` (`FORMAT` and the samples) for a header
/// of `samples` samples and a record of `alleles` alleles: the keys and
/// each sample decoded. `place` names the record in the count error.
fn decode_all<'t>(
    text: &'t str,
    samples: usize,
    alleles: usize,
    percent: bool,
    place: &str,
) -> Result<Genotypes<'t>, String> {
    let mut fields = text.split('\t');
    let format = fields.next().unwrap_or_default();
    let found = text.bytes().filter(|&b| b == b'\t').count();
    if found != samples {
        return Err(format!(
            "there are {found} genotypes while the header requires that {samples} genotypes be \
             present for all records{place}"
        ));
    }
    let keys: Vec<&str> = format.split(':').collect();
    let mut genotypes = Genotypes {
        samples: Vec::with_capacity(samples),
        ..Genotypes::default()
    };
    let mut reading = Reading {
        keys: &keys,
        alleles,
        percent,
        pl_read: false,
        values: Vec::new(),
    };
    for sample in fields {
        let decoded = decode_sample(&mut reading, sample, &mut genotypes)?;
        genotypes.samples.push(decoded);
    }
    genotypes.keys = keys;
    Ok(genotypes)
}

/// The errors the codec gives decoding the genotypes of a record whose
/// header names its samples out of sorted order.
pub(super) fn check(
    text: &str,
    samples: usize,
    alleles: usize,
    percent: bool,
) -> Result<(), String> {
    decode_all(text, samples, alleles, percent, "").map(|_| ())
}

/// `calcVCFGenotypeKeys`: `GT` first when a genotype has alleles, then
/// every other key a genotype holds in sorted order; `GT` alone when none.
fn genotype_keys<'t>(decoded: &[Decoded<'t>]) -> Vec<&'t str> {
    let mut seen: HashSet<&str> = HashSet::new();
    for genotype in decoded {
        for (key, _) in &genotype.attributes {
            seen.insert(key);
        }
    }
    let any = |test: fn(&Decoded) -> bool| decoded.iter().any(test);
    for (name, present) in [
        ("GQ", any(|g| g.gq.is_some())),
        ("DP", any(|g| g.dp.is_some())),
        ("AD", any(|g| g.ad.is_some())),
        ("PL", any(|g| g.pl.is_some())),
        ("FT", any(|g| g.filters.is_some())),
    ] {
        if present {
            seen.insert(name);
        }
    }
    let mut names: Vec<&str> = seen.into_iter().collect();
    names.sort_by(|a, b| java_compare(a, b));
    if any(|g| g.alleles.end > g.alleles.start) {
        names.insert(0, "GT");
    }
    if names.is_empty() {
        names.push("GT");
    }
    names
}

/// Append one genotype's value of `name` to `out`, as
/// `appendGenotypeData` prints it: an absent value is `.` whatever the
/// key's `Number`.
fn push_value(out: &mut String, genotypes: &Genotypes<'_>, genotype: &Decoded<'_>, name: &str) {
    use std::fmt::Write as _;
    let mut ints = |values: Option<Span>| match values {
        None => out.push('.'),
        Some(span) => {
            for (i, value) in span.in_(&genotypes.ints).iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                let _ = write!(out, "{value}");
            }
        }
    };
    match name {
        "AD" => ints(genotype.ad),
        "PL" => ints(genotype.pl),
        "FT" => out.push_str(genotype.filters.as_deref().unwrap_or("PASS")),
        // `GQAccessor.getValue`: `Math.min(GQ, 99)`.
        "GQ" => match genotype.gq {
            None => out.push('.'),
            Some(q) => {
                let _ = write!(out, "{}", q.min(99));
            }
        },
        "DP" => match genotype.dp {
            None => out.push('.'),
            Some(d) => {
                let _ = write!(out, "{d}");
            }
        },
        _ => out.push_str(genotype.attribute(name).unwrap_or(".")),
    }
}

/// The genotypes of one record, decoded by [`decode`] before the fixed
/// columns are written, so their errors come first as in htsjdk.
pub(super) struct Decoding<'t>(Genotypes<'t>);

/// Decode `text` (`FORMAT` and the samples) of a record with `alleles`
/// alleles.
pub(super) fn decode<'t>(
    header: &VcfHeader,
    text: &'t str,
    alleles: usize,
    site: &Site,
) -> Result<Decoding<'t>, String> {
    let place = format!(" at {}:{}", site.contig, site.start);
    decode_all(
        text,
        header.samples().len(),
        alleles,
        header.minor() >= 3,
        &place,
    )
    .map(Decoding)
}

/// Append the genotypes [`decode`] read, encoded anew.
pub(super) fn write_decoded(
    header: &VcfHeader,
    decoding: &Decoding<'_>,
    site: &Site,
    out: &mut String,
) -> Result<(), String> {
    write(header, &decoding.0, site, out)
}

fn write(
    header: &VcfHeader,
    genotypes: &Genotypes<'_>,
    site: &Site,
    out: &mut String,
) -> Result<(), String> {
    use std::fmt::Write as _;
    if header.minor() < 1 && !genotypes.keys.contains(&"GT") {
        return Err(
            "Unable to find the GT field for the record; the GT field is required before VCF4.1"
                .to_owned(),
        );
    }
    let names = genotype_keys(&genotypes.samples);
    if let Some(name) = names
        .iter()
        .find(|n| **n != "GT" && header.format(n).is_none())
    {
        return Err(format!(
            "Key {name} found in VariantContext field FORMAT at {}:{} but this key isn't defined \
             in the VCFHeader.  We require all VCFs to have complete VCF headers by default.",
            site.contig, site.start
        ));
    }
    out.push_str(&names.join(":"));
    let with_gt = names.first().is_some_and(|n| *n == "GT");
    let others = if with_gt { &names[1..] } else { &names[..] };
    // One sample's values, each ended at an offset of `values`.
    let mut values = String::new();
    let mut ends: Vec<usize> = Vec::with_capacity(others.len());
    for genotype in &genotypes.samples {
        out.push('\t');
        if with_gt {
            let alleles = genotype.alleles.in_(&genotypes.alleles);
            if alleles.is_empty() {
                return Err(
                    "GTs cannot be missing for some samples if they are available for others in \
                     the record"
                        .to_owned(),
                );
            }
            let separator = if genotype.phased { '|' } else { '/' };
            for (i, allele) in alleles.iter().enumerate() {
                if i > 0 {
                    out.push(separator);
                }
                match allele {
                    Some(index) => {
                        let _ = write!(out, "{index}");
                    }
                    None => out.push('.'),
                }
            }
        }
        values.clear();
        ends.clear();
        for name in others {
            push_value(&mut values, genotypes, genotype, name);
            ends.push(values.len());
        }
        // Trailing missing values are dropped.
        let mut kept = ends.len();
        while kept > 0 {
            let start = if kept > 1 { ends[kept - 2] } else { 0 };
            if !is_missing(&values[start..ends[kept - 1]]) {
                break;
            }
            kept -= 1;
        }
        let mut start = 0;
        for (i, &end) in ends[..kept].iter().enumerate() {
            if i > 0 || with_gt {
                out.push(':');
            }
            out.push_str(&values[start..end]);
            start = end;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::encode::Encoder;
    use super::*;

    fn header() -> VcfHeader {
        let text = "##fileformat=VCFv4.2\n\
                    ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
                    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
                    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
                    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
                    ##FORMAT=<ID=FT,Number=1,Type=String,Description=\"f\">\n\
                    ##FORMAT=<ID=XX,Number=2,Type=String,Description=\"x\">\n\
                    ##FORMAT=<ID=XA,Number=A,Type=String,Description=\"x\">\n\
                    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"p\">\n\
                    ##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"l\">\n\
                    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n";
        VcfHeader::read(&mut text.as_bytes()).expect("a header")
    }

    fn encode(line: &str) -> Result<String, String> {
        let h = header();
        let mut out = String::new();
        Encoder::new(&h).encode(&h, line, &mut out)?;
        Ok(out)
    }

    #[test]
    fn genotypes_are_decoded_and_written_anew() {
        let cases = [
            (
                "GT:GQ:DP:AD\t0/1:30.5:10:3,7\t1|1:.:.:.",
                "GT:AD:DP:GQ\t0/1:3,7:10:31\t1|1",
            ),
            ("GT:FT\t0/1:PASS\t0/1:b;a", "GT:FT\t0/1:PASS\t0/1:a;b"),
            ("GT:FT\t0/1:PASS\t0/1:.", "GT\t0/1\t0/1"),
            ("GT:XX\t0/1:.\t0/1:z", "GT:XX\t0/1\t0/1:z"),
            ("GQ:XX\t30:.\t.:z", "GQ:XX\t30\t.:z"),
            ("GT:AD\t0/1:1,x\t./.:2,3", "GT:AD\t0/1\t./.:2,3"),
            ("GT\t0//1\t1/", "GT\t0/1\t1"),
            // Picard 3.5.0 VcfFormatConverter on the same records.
            ("GT:XA:XX\t0/1:b\t0/1:.:c", "GT:XA:XX\t0/1:b\t0/1:.:c"),
            ("GT:XA:DP\t0/1:.:5\t0/1:.:6", "GT:DP\t0/1:5\t0/1:6"),
            ("GT:XX:XA\t0/1:.:b\t0/1", "GT:XA\t0/1:b\t0/1"),
            ("GT:XA:XX\t0/1::b\t0/1:z", "GT:XA:XX\t0/1::b\t0/1:z"),
            ("GT:XA:XX\t0/1:.:b\t0/1:z", "GT:XA:XX\t0/1:.:b\t0/1:z"),
            ("GT:GQ\t0/1:-1.5\t0/1:-0.6", "GT\t0/1\t0/1"),
            ("GT:GQ\t0/1:-1.0\t0/1:-0.4", "GT:GQ\t0/1\t0/1:0"),
            (
                "GT:GL\t0/1:-1.5,-0.1,-100\t1/1:0,-1,-2",
                "GT:PL\t0/1:14,0,999\t1/1:0,10,20",
            ),
            (
                "GT:GL:PL\t0/1:-1.5,-0.1,-100:5,6,7\t1/1:0,-1,-2",
                "GT:PL\t0/1:5,6,7\t1/1",
            ),
            (
                "GT:PL:GL\t0/1:5,6,7:-1,0,-2\t1/1:0,-1,-2",
                "GT:PL\t0/1:5,6,7\t1/1:0,-1,-2",
            ),
            ("GT:GL\t0/1:.,.,.\t1/1:.", "GT\t0/1\t1/1"),
            (
                "GT:GL\t0/1:1.5,0.2,-3\t1/1:-0.05,-0.15,-0.25",
                "GT:PL\t0/1:0,13,45\t1/1:0,1,2",
            ),
            (
                "GT:GL\t0/1:NaN,0,-1\t1/1:-Infinity,0,-1",
                "GT:PL\t0/1:0,0,0\t1/1:2147483647,0,10",
            ),
        ];
        for (genotypes, expected) in cases {
            let line = format!("c1\t1\t.\tA\tC\t.\t.\t.\t{genotypes}");
            assert_eq!(
                encode(&line).as_deref(),
                Ok(format!("c1\t1\t.\tA\tC\t.\t.\t.\t{expected}").as_str()),
                "{genotypes}"
            );
        }
        for bad in [
            "GT\t0/1\t0/1:5",
            "GQ:GT\t1:0\t1:0",
            "GT:ZZ\t0:1\t0:1",
            "GT\t0/3\t0",
            "GT:GL\t0/1:x,1,2\t0",
            "GT:GL\t0/1:.,1,2\t0",
        ] {
            let line = format!("c1\t1\t.\tA\tC\t.\t.\t.\t{bad}");
            assert!(encode(&line).is_err(), "{bad}");
        }
    }
}

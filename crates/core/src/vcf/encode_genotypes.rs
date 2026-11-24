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
//! other `.` value absent; `GQ` `-1` absent, else `Math.round` of the
//! number; `DP` an `int`, `-1` absent; `AD` and `PL` a list of `int`
//! (absent when one is not); every other key kept as its text. From VCF
//! 4.3 each value is percent-decoded.
//!
//! Encoding: the keys are `GT` when a genotype has alleles, then in sorted
//! order every key a genotype holds (`GQ`, `DP`, `AD`, `PL` when one has a
//! value, `FT` when one is filtered). A key without a `##FORMAT` line
//! fails. A genotype prints `PASS` for an unfiltered `FT`, a `GQ` over 99
//! as 99, `.` for an
//! absent value (`.,.` for a key whose `Number` asks for more than one),
//! and drops its trailing missing values.

use std::collections::HashSet;

use super::Site;
use super::full_header::{Count, VcfHeader, java_compare};
use super::site::percent_decode;

/// One sample as the codec decodes it.
#[derive(Default)]
struct Decoded {
    /// Allele indexes, `None` for a no-call; empty when there is no `GT`.
    alleles: Vec<Option<usize>>,
    phased: bool,
    gq: Option<i32>,
    dp: Option<i32>,
    ad: Option<Vec<i32>>,
    pl: Option<Vec<i32>>,
    /// The filters joined by `;`, `None` when unfiltered.
    filters: Option<String>,
    /// Every other key with its value, in `FORMAT` order.
    attributes: Vec<(String, String)>,
}

impl Decoded {
    fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .rev()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }
}

/// `Math.round(double)` cast to an `int`, as `GenotypeBuilder.GQ` takes it.
fn java_round_int(value: f64) -> i32 {
    if value.is_nan() {
        return 0;
    }
    let rounded = (value + 0.5).floor();
    let long = if rounded >= 9.223_372_036_854_776e18 {
        i64::MAX
    } else if rounded <= -9.223_372_036_854_776e18 {
        i64::MIN
    } else {
        rounded as i64
    };
    long as i32
}

/// `decodeInts`: `None` when a value is no Java `int`.
fn decode_ints(text: &str) -> Option<Vec<i32>> {
    text.split(',')
        .map(|v| crate::java::parse_integer(v, 32).ok().map(|n| n as i32))
        .collect()
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

/// `parseGenotypeAlleles`: `StringTokenizer` over `/|\`.
fn genotype_alleles(text: &str, alleles: usize) -> Result<Vec<Option<usize>>, String> {
    text.split(['/', '|', '\\'])
        .filter(|token| !token.is_empty())
        .map(|token| {
            if token == "." {
                return Ok(None);
            }
            let index = crate::java::parse_integer(token, 32).map_err(|_| {
                format!(
                    "The following invalid GT allele index was encountered in the file: {token}"
                )
            })?;
            match usize::try_from(index) {
                Ok(i) if i < alleles => Ok(Some(i)),
                _ => Err(format!(
                    "The allele with index {token} is not defined in the REF/ALT columns in the \
                     record"
                )),
            }
        })
        .collect()
}

fn decode_sample(
    keys: &[&str],
    text: &str,
    alleles: usize,
    percent: bool,
) -> Result<Decoded, String> {
    let values: Vec<String> = text
        .split(':')
        .map(|v| {
            if percent {
                percent_decode(v)
            } else {
                v.to_owned()
            }
        })
        .collect();
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
                decoded.gq = if value == "-1" {
                    None
                } else {
                    let q = crate::java::parse_double(value)
                        .or_else(|e| super::site::special_double(value).ok_or(e))?;
                    Some(java_round_int(q))
                };
            }
            "DP" => {
                let dp = crate::java::parse_integer(value, 32)?;
                decoded.dp = (dp != -1).then_some(dp as i32);
            }
            "AD" => decoded.ad = decode_ints(value),
            "PL" => decoded.pl = decode_ints(value),
            _ => decoded.attributes.push(((*key).to_owned(), value.clone())),
        }
    }
    if let Some(at) = gt_at {
        if at > 0 {
            return Err(format!(
                "Saw GT field at position {at}, but it must be at the first position for \
                 genotypes when present"
            ));
        }
        let gt = values.get(at).map_or("", String::as_str);
        decoded.alleles = genotype_alleles(gt, alleles)?;
        decoded.phased = gt.contains('|');
    }
    Ok(decoded)
}

/// `isMissingValue`: only `.` and `,`.
fn is_missing(value: &str) -> bool {
    value.bytes().all(|b| b == b'.' || b == b',')
}

fn join_ints(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
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
) -> Result<(Vec<&'t str>, Vec<Decoded>), String> {
    let fields: Vec<&str> = text.split('\t').collect();
    if fields.len() != samples + 1 {
        return Err(format!(
            "there are {} genotypes while the header requires that {samples} genotypes be \
             present for all records{place}",
            fields.len().saturating_sub(1),
        ));
    }
    let keys: Vec<&str> = fields[0].split(':').collect();
    let decoded = fields[1..]
        .iter()
        .map(|sample| decode_sample(&keys, sample, alleles, percent))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((keys, decoded))
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
fn genotype_keys(decoded: &[Decoded]) -> Vec<String> {
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
    let mut names: Vec<String> = seen.into_iter().map(str::to_owned).collect();
    names.sort_by(|a, b| java_compare(a, b));
    if any(|g| !g.alleles.is_empty()) {
        names.insert(0, "GT".to_owned());
    }
    if names.is_empty() {
        names.push("GT".to_owned());
    }
    names
}

/// One genotype's value of `name`, as `appendGenotypeData` prints it.
fn value_of(
    genotype: &Decoded,
    name: &str,
    header: &VcfHeader,
    alleles: usize,
    ploidy: usize,
) -> String {
    let ints = |values: Option<&[i32]>| values.map_or_else(|| ".".to_owned(), join_ints);
    match name {
        "FT" => genotype
            .filters
            .clone()
            .unwrap_or_else(|| "PASS".to_owned()),
        // `GQAccessor.getValue`: `Math.min(GQ, 99)`.
        "GQ" => genotype
            .gq
            .map_or_else(|| ".".to_owned(), |q| q.min(99).to_string()),
        "DP" => genotype
            .dp
            .map_or_else(|| ".".to_owned(), |d| d.to_string()),
        "AD" => ints(genotype.ad.as_deref()),
        "PL" => ints(genotype.pl.as_deref()),
        _ => match genotype.attribute(name) {
            Some(value) => value.to_owned(),
            None => {
                let count = header
                    .format(name)
                    .map_or(Count::Unbounded, |d| d.count)
                    .values(alleles, ploidy);
                let count = usize::try_from(count).unwrap_or(1).max(1);
                vec!["."; count].join(",")
            }
        },
    }
}

/// Decode `text` (`FORMAT` and the samples) and append it encoded anew.
pub(super) fn encode(
    header: &VcfHeader,
    text: &str,
    alleles: &[String],
    site: &Site,
    out: &mut String,
) -> Result<(), String> {
    let place = format!(" at {}:{}", site.contig, site.start);
    let (keys, decoded) = decode_all(
        text,
        header.samples().len(),
        alleles.len(),
        header.minor() >= 3,
        &place,
    )?;
    if header.minor() < 1 && !keys.contains(&"GT") {
        return Err(
            "Unable to find the GT field for the record; the GT field is required before VCF4.1"
                .to_owned(),
        );
    }
    let names = genotype_keys(&decoded);
    if let Some(name) = names
        .iter()
        .find(|n| *n != "GT" && header.format(n).is_none())
    {
        return Err(format!(
            "Key {name} found in VariantContext field FORMAT at {}:{} but this key isn't defined \
             in the VCFHeader.  We require all VCFs to have complete VCF headers by default.",
            site.contig, site.start
        ));
    }
    out.push_str(&names.join(":"));
    let ploidy = decoded
        .iter()
        .map(|g| g.alleles.len())
        .max()
        .filter(|&p| p > 0)
        .unwrap_or(2);
    let with_gt = names.first().is_some_and(|n| n == "GT");
    for genotype in &decoded {
        out.push('\t');
        if with_gt {
            if genotype.alleles.is_empty() {
                return Err(
                    "GTs cannot be missing for some samples if they are available for others in \
                     the record"
                        .to_owned(),
                );
            }
            let separator = if genotype.phased { '|' } else { '/' };
            for (i, allele) in genotype.alleles.iter().enumerate() {
                if i > 0 {
                    out.push(separator);
                }
                match allele {
                    Some(index) => out.push_str(&index.to_string()),
                    None => out.push('.'),
                }
            }
        }
        let mut values: Vec<String> = names
            .iter()
            .filter(|n| *n != "GT")
            .map(|name| value_of(genotype, name, header, alleles.len(), ploidy))
            .collect();
        while values.last().is_some_and(|v| is_missing(v)) {
            values.pop();
        }
        for (i, value) in values.iter().enumerate() {
            if i > 0 || with_gt {
                out.push(':');
            }
            out.push_str(value);
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
        ] {
            let line = format!("c1\t1\t.\tA\tC\t.\t.\t.\t{bad}");
            assert!(encode(&line).is_err(), "{bad}");
        }
    }
}

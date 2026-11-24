//! The whole header of a VCF as htsjdk 5.0.0's `VCFCodec` reads it and its
//! `VCFWriter` writes it again, for the tools that write a VCF
//! ([`crate::vcf_writer`]).
//!
//! Reading (`AbstractVCFCodec.parseHeaderFromLines`, checked by
//! [`super::site_header::parse`]):
//!
//! - `##INFO` and `##FORMAT` lines keep `ID`, `Number`, `Type`,
//!   `Description`, and from VCF 4.2 `Source` and `Version`; any other key
//!   of the line is dropped. A missing `Description` reads `Not provided in original VCF
//!   header`. A `Flag` with a fixed `Number` gets `Number=0`.
//! - A line of the standard set (`VCFStandardHeaderLines`: `FORMAT` `GT`,
//!   `GQ`, `DP`, `PL`, `AD`, `FT`, `PS`, `PQ`; `INFO` `END`, `DB`, `DP`,
//!   `SB`, `AF`, `AC`, `AN`, `MQ0`, `MQ`, `SOMATIC`) whose count kind,
//!   fixed count or type differs from the standard one is replaced by it
//!   (`repairStandardHeaderLines`). `FT` with `Number=1` becomes
//!   `Number=.,...,Description="Genotype-level filter"`.
//! - `##FILTER`, `##ALT`, `##contig`, `##META`, `##SAMPLE` and from VCF
//!   4.3 `##PEDIGREE` keep every key in order (`VCF4Parser.parseLine`: a
//!   repeated key replaces the earlier value in its place). A `##contig`
//!   whose `ID` an earlier one gave stays a line; the dictionary takes the
//!   first.
//! - Any other `##` line with a `=` is a key and a value; one without is
//!   dropped. The `##fileformat` line is not kept.
//!
//! Writing (`VCFWriter.writeHeader`): `##fileformat=VCFv4.2` first, then
//! the lines in sorted order, each once: contig lines by their place among
//! the contig lines, every other pair of lines by their text
//! (`VCFHeaderLine.compareTo`, `String.compareTo`). A value is quoted when
//! it holds a comma or a space or its key is `Description`, `Source` or
//! `Version`; a quote inside not after a backslash gets one
//! (`escapeQuotes`). Then the column line, with `FORMAT` and the samples
//! when the header names samples, each once.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;

use super::site_header::{self, RawHeader, SiteHeader, after, column_samples, parse_fields};
use super::{Contig, Dictionary};

/// The `Number` of an `##INFO` or `##FORMAT` line (`VCFHeaderLineCount`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Count {
    Fixed(i32),
    /// One per alternate allele.
    A,
    /// One per allele.
    R,
    /// One per genotype.
    G,
    /// `.`: any number.
    Unbounded,
}

impl Count {
    fn parse(text: &str) -> Option<Self> {
        match text {
            "A" => Some(Self::A),
            "R" => Some(Self::R),
            "G" => Some(Self::G),
            "." => Some(Self::Unbounded),
            _ => crate::java::parse_integer(text, 32)
                .ok()
                .map(|n| Self::Fixed(n as i32)),
        }
    }

    fn text(self) -> String {
        match self {
            Self::Fixed(n) => n.to_string(),
            Self::A => "A".to_owned(),
            Self::R => "R".to_owned(),
            Self::G => "G".to_owned(),
            Self::Unbounded => ".".to_owned(),
        }
    }

    /// `VCFCompoundHeaderLine.getCount(VariantContext)` for a record with
    /// `alleles` alleles (the reference included) and a largest ploidy of
    /// `ploidy`: -1 for `.`.
    pub fn values(self, alleles: usize, ploidy: usize) -> i64 {
        match self {
            Self::Fixed(n) => i64::from(n),
            Self::A => alleles as i64 - 1,
            Self::R => alleles as i64,
            Self::G => likelihoods(alleles, ploidy),
            Self::Unbounded => -1,
        }
    }
}

/// `GenotypeLikelihoods.numLikelihoods`: the genotypes of `ploidy` alleles
/// drawn from `alleles`, `C(alleles + ploidy - 1, ploidy)`.
fn likelihoods(alleles: usize, ploidy: usize) -> i64 {
    let n = (alleles + ploidy).saturating_sub(1) as i64;
    let k = ploidy as i64;
    let mut result: i64 = 1;
    for i in 0..k {
        result = result.saturating_mul(n - i) / (i + 1);
    }
    result
}

/// What an `##INFO` or `##FORMAT` line declares.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Definition {
    pub count: Count,
    pub flag: bool,
}

/// One kept header line.
#[derive(Clone, Debug, PartialEq, Eq)]
enum Line {
    /// `##INFO` or `##FORMAT`.
    Compound {
        key: &'static str,
        id: String,
        count: Count,
        ty: String,
        description: String,
        source: Option<String>,
        version: Option<String>,
    },
    /// A line of `<key=value,...>` pairs; `contig` is its place among the
    /// contig lines.
    Structured {
        key: String,
        fields: Vec<(String, String)>,
        contig: Option<usize>,
    },
    /// `##key=value`.
    Other { key: String, value: String },
}

/// `VCFHeaderLine.escapeQuotes`: `replaceAll("([^\\\\])\"", "$1\\\\\"")`,
/// a quote after any character but a backslash gets a backslash; matches
/// do not overlap, so of two quotes in a row only the first gets one.
fn escape_quotes(value: &str) -> String {
    let chars: Vec<char> = value.chars().collect();
    let mut out = String::with_capacity(value.len());
    let mut i = 0;
    while i < chars.len() {
        if i + 1 < chars.len() && chars[i] != '\\' && chars[i + 1] == '"' {
            out.push(chars[i]);
            out.push('\\');
            out.push('"');
            i += 2;
        } else {
            out.push(chars[i]);
            i += 1;
        }
    }
    out
}

/// `VCFHeaderLine.toStringEncoding(Map)`.
fn encode_fields(fields: &[(&str, &str)]) -> String {
    let mut out = String::from("<");
    for (i, (key, value)) in fields.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(key);
        out.push('=');
        let quoted = value.contains(',')
            || value.contains(' ')
            || matches!(*key, "Description" | "Source" | "Version");
        if quoted {
            out.push('"');
            out.push_str(&escape_quotes(value));
            out.push('"');
        } else {
            out.push_str(value);
        }
    }
    out.push('>');
    out
}

impl Line {
    /// `toString`: the line without its `##`.
    fn text(&self) -> String {
        match self {
            Self::Compound {
                key,
                id,
                count,
                ty,
                description,
                source,
                version,
            } => {
                let count = count.text();
                let mut fields = vec![
                    ("ID", id.as_str()),
                    ("Number", count.as_str()),
                    ("Type", ty.as_str()),
                    ("Description", description.as_str()),
                ];
                if let Some(source) = source {
                    fields.push(("Source", source));
                }
                if let Some(version) = version {
                    fields.push(("Version", version));
                }
                format!("{key}={}", encode_fields(&fields))
            }
            Self::Structured { key, fields, .. } => {
                let pairs: Vec<(&str, &str)> = fields
                    .iter()
                    .map(|(k, v)| (k.as_str(), v.as_str()))
                    .collect();
                format!("{key}={}", encode_fields(&pairs))
            }
            Self::Other { key, value } => format!("{key}={value}"),
        }
    }
}

/// The standard lines of `VCFStandardHeaderLines`: key, ID, count, type,
/// description.
const STANDARD: [(&str, &str, Count, &str, &str); 18] = [
    ("FORMAT", "GT", Count::Fixed(1), "String", "Genotype"),
    (
        "FORMAT",
        "GQ",
        Count::Fixed(1),
        "Integer",
        "Genotype Quality",
    ),
    (
        "FORMAT",
        "DP",
        Count::Fixed(1),
        "Integer",
        "Approximate read depth (reads with MQ=255 or with bad mates are filtered)",
    ),
    (
        "FORMAT",
        "PL",
        Count::G,
        "Integer",
        "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification",
    ),
    (
        "FORMAT",
        "AD",
        Count::R,
        "Integer",
        "Allelic depths for the ref and alt alleles in the order listed",
    ),
    (
        "FORMAT",
        "FT",
        Count::Unbounded,
        "String",
        "Genotype-level filter",
    ),
    (
        "FORMAT",
        "PS",
        Count::Fixed(1),
        "Integer",
        "Phasing set (typically the position of the first variant in the set)",
    ),
    (
        "FORMAT",
        "PQ",
        Count::Fixed(1),
        "Float",
        "Read-backed phasing quality",
    ),
    (
        "INFO",
        "END",
        Count::Fixed(1),
        "Integer",
        "Stop position of the interval",
    ),
    ("INFO", "DB", Count::Fixed(0), "Flag", "dbSNP Membership"),
    (
        "INFO",
        "DP",
        Count::Fixed(1),
        "Integer",
        "Approximate read depth; some reads may have been filtered",
    ),
    ("INFO", "SB", Count::Fixed(1), "Float", "Strand Bias"),
    (
        "INFO",
        "AF",
        Count::A,
        "Float",
        "Allele Frequency, for each ALT allele, in the same order as listed",
    ),
    (
        "INFO",
        "AC",
        Count::A,
        "Integer",
        "Allele count in genotypes, for each ALT allele, in the same order as listed",
    ),
    (
        "INFO",
        "AN",
        Count::Fixed(1),
        "Integer",
        "Total number of alleles in called genotypes",
    ),
    (
        "INFO",
        "MQ0",
        Count::Fixed(1),
        "Integer",
        "Total Mapping Quality Zero Reads",
    ),
    (
        "INFO",
        "MQ",
        Count::Fixed(1),
        "Float",
        "RMS Mapping Quality",
    ),
    ("INFO", "SOMATIC", Count::Fixed(0), "Flag", "Somatic event"),
];

/// `Standards.repair`: the standard line in place of `line` when their
/// count kind, fixed count or type differ.
fn repair(line: Line) -> Line {
    let Line::Compound {
        key, id, count, ty, ..
    } = &line
    else {
        return line;
    };
    let Some((_, _, standard_count, standard_ty, description)) = STANDARD
        .iter()
        .find(|(k, i, ..)| k == key && i == id)
        .copied()
    else {
        return line;
    };
    let same_kind = std::mem::discriminant(count) == std::mem::discriminant(&standard_count);
    let bad_count = same_kind && *count != standard_count;
    if same_kind && !bad_count && ty == standard_ty {
        return line;
    }
    Line::Compound {
        key,
        id: id.clone(),
        count: standard_count,
        ty: standard_ty.to_owned(),
        description: description.to_owned(),
        source: None,
        version: None,
    }
}

/// The text htsjdk puts in place of a missing `Description`.
const NO_DESCRIPTION: &str = "Not provided in original VCF header";

/// `##INFO` or `##FORMAT` (`VCFCompoundHeaderLine`): the checks are
/// [`site_header::parse`]'s, made before.
fn compound(key: &'static str, text: &str, minor: u8) -> Result<Line, String> {
    let recommended: &[&str] = if minor >= 2 {
        &["Source", "Version"]
    } else {
        &[]
    };
    let fields = parse_fields(
        text,
        Some(&["ID", "Number", "Type", "Description"]),
        recommended,
    )?;
    let get = |name: &str| {
        fields
            .iter()
            .find(|(k, _)| k == name)
            .map(|(_, v)| v.clone())
    };
    let ty = get("Type").unwrap_or_default();
    let mut count = get("Number")
        .as_deref()
        .and_then(Count::parse)
        .unwrap_or(Count::Unbounded);
    if ty == "Flag" && matches!(count, Count::Fixed(_)) {
        count = Count::Fixed(0);
    }
    Ok(repair(Line::Compound {
        key,
        id: get("ID").unwrap_or_default(),
        count,
        ty,
        description: get("Description").unwrap_or_else(|| NO_DESCRIPTION.to_owned()),
        // htsjdk keeps `Source` and `Version` from VCF 4.2 on.
        source: get("Source").filter(|_| minor >= 2),
        version: get("Version").filter(|_| minor >= 2),
    }))
}

/// The header a VCF writer writes, with what its encoder looks up.
#[derive(Clone, Debug)]
pub struct VcfHeader {
    lines: Vec<Line>,
    samples: Vec<String>,
    minor: u8,
    site: SiteHeader,
    infos: HashMap<String, Definition>,
    formats: HashMap<String, Definition>,
    filters: HashSet<String>,
    dictionary: Dictionary,
}

/// Java's `String.compareTo`: UTF-16 code units in order.
pub fn java_compare(a: &str, b: &str) -> Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

impl VcfHeader {
    /// Read the header of `reader` as `VCFCodec` reads it, leaving the
    /// reader at the first record.
    pub fn read(reader: &mut impl BufRead) -> Result<Self, String> {
        Self::from_raw(&site_header::read_raw(reader)?)
    }

    /// The header of `raw`, checked line by line.
    pub fn from_raw(raw: &RawHeader) -> Result<Self, String> {
        let site = site_header::parse(raw)?;
        let mut lines = Vec::new();
        let mut contigs: HashSet<String> = HashSet::new();
        let mut dictionary = Vec::new();
        let mut contig_lines = 0usize;
        for line in &raw.lines {
            let kept = if line.starts_with("##INFO") {
                Some(compound("INFO", after(line, 7)?, raw.minor)?)
            } else if line.starts_with("##FORMAT") {
                Some(compound("FORMAT", after(line, 9)?, raw.minor)?)
            } else if line.starts_with("##contig") {
                let fields = parse_fields(after(line, 9)?, None, &[])?;
                let id = value(&fields, "ID").unwrap_or_default();
                if contigs.insert(id.clone()) {
                    dictionary.push(Contig {
                        name: id,
                        length: value(&fields, "length").and_then(|l| l.parse().ok()),
                    });
                }
                contig_lines += 1;
                Some(Line::Structured {
                    key: "contig".to_owned(),
                    fields,
                    contig: Some(contig_lines),
                })
            } else if let Some((key, skip)) = structured_key(line, raw.minor) {
                Some(Line::Structured {
                    key: key.to_owned(),
                    fields: parse_fields(after(line, skip)?, None, &[])?,
                    contig: None,
                })
            } else {
                line.find('=').and_then(|equals| {
                    let key = &line[2..equals];
                    (!matches!(key, "fileformat" | "format")).then(|| Line::Other {
                        key: key.to_owned(),
                        value: line[equals + 1..].to_owned(),
                    })
                })
            };
            lines.extend(kept);
        }
        let mut header = Self {
            lines,
            samples: column_samples(&raw.columns)?,
            minor: raw.minor,
            site,
            infos: HashMap::new(),
            formats: HashMap::new(),
            filters: HashSet::new(),
            dictionary,
        };
        header.index();
        Ok(header)
    }

    /// The lookups of the encoder: the first line of each `ID`.
    fn index(&mut self) {
        for line in &self.lines {
            match line {
                Line::Compound {
                    key, id, count, ty, ..
                } => {
                    let map = if *key == "INFO" {
                        &mut self.infos
                    } else {
                        &mut self.formats
                    };
                    map.entry(id.clone()).or_insert(Definition {
                        count: *count,
                        flag: ty == "Flag",
                    });
                }
                Line::Structured { key, fields, .. } if key == "FILTER" => {
                    if let Some(id) = value(fields, "ID") {
                        self.filters.insert(id);
                    }
                }
                _ => {}
            }
        }
    }

    /// The sample names, each once, in column order.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// The minor version of `VCFv4.N` the input named.
    pub fn minor(&self) -> u8 {
        self.minor
    }

    /// What the site reader of [`super::site`] needs.
    pub fn site(&self) -> &SiteHeader {
        &self.site
    }

    /// The `##contig` lines as a dictionary, `None` without any.
    pub fn dictionary(&self) -> Option<&Dictionary> {
        (!self.dictionary.is_empty()).then_some(&self.dictionary)
    }

    pub fn info(&self, id: &str) -> Option<&Definition> {
        self.infos.get(id)
    }

    pub fn format(&self, id: &str) -> Option<&Definition> {
        self.formats.get(id)
    }

    pub fn has_filter(&self, id: &str) -> bool {
        self.filters.contains(id)
    }

    /// The header text as `VCFWriter.writeHeader` writes it, each line
    /// ended by `\n`.
    pub fn text(&self) -> String {
        let mut lines: Vec<(Option<usize>, String)> = self
            .lines
            .iter()
            .map(|line| {
                let contig = match line {
                    Line::Structured { contig, .. } => *contig,
                    _ => None,
                };
                (contig, line.text())
            })
            .collect();
        lines.sort_by(|(a_contig, a), (b_contig, b)| match (a_contig, b_contig) {
            (Some(x), Some(y)) => x.cmp(y),
            _ => java_compare(a, b),
        });
        lines.dedup_by(|(a_contig, a), (b_contig, b)| {
            a_contig.is_none() && b_contig.is_none() && a == b
        });
        let mut out = String::from("##fileformat=VCFv4.2\n");
        for (_, line) in &lines {
            out.push_str("##");
            out.push_str(line);
            out.push('\n');
        }
        out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
        if !self.samples.is_empty() {
            out.push_str("\tFORMAT");
            for sample in &self.samples {
                out.push('\t');
                out.push_str(sample);
            }
        }
        out.push('\n');
        out
    }
}

/// The key and text start of a line of `<...>` pairs other than `INFO`,
/// `FORMAT` and `contig`, in the codec's order of tests.
fn structured_key(line: &str, minor: u8) -> Option<(&'static str, usize)> {
    if line.starts_with("##FILTER") {
        Some(("FILTER", 9))
    } else if line.starts_with("##ALT") {
        Some(("ALT", 6))
    } else if line.starts_with("##PEDIGREE") && minor >= 3 {
        Some(("PEDIGREE", 11))
    } else if line.starts_with("##META") {
        Some(("META", 7))
    } else if line.starts_with("##SAMPLE") {
        Some(("SAMPLE", 9))
    } else {
        None
    }
}

fn value(fields: &[(String, String)], key: &str) -> Option<String> {
    fields
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.clone())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(text: &str) -> VcfHeader {
        VcfHeader::read(&mut text.as_bytes()).expect("a header")
    }

    const COLUMNS: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

    #[test]
    fn lines_are_sorted_repaired_and_kept_once() {
        let text = format!(
            "##fileformat=VCFv4.1\n##source=x y\n##contig=<ID=c2,length=3000>\n\
             ##contig=<ID=c1,length=5000>\n##contig=<ID=c2,length=1>\n\
             ##INFO=<ID=XT,Number=1,Type=String,Description=\"a\",Foo=bar>\n\
             ##INFO=<ID=NODESC,Number=1,Type=Integer>\n\
             ##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n\
             ##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"mq\">\n\
             ##META=<ID=Assay,Values=[A, B]>\n##noequals\n##k=v\n##k=v\n{COLUMNS}"
        );
        let written = header(&text).text();
        assert_eq!(
            written,
            "##fileformat=VCFv4.2\n\
             ##FORMAT=<ID=FT,Number=.,Type=String,Description=\"Genotype-level filter\">\n\
             ##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">\n\
             ##INFO=<ID=NODESC,Number=1,Type=Integer,Description=\"Not provided in original VCF header\">\n\
             ##INFO=<ID=XT,Number=1,Type=String,Description=\"a\">\n\
             ##META=<ID=Assay,Values=B]>\n\
             ##contig=<ID=c2,length=3000>\n##contig=<ID=c1,length=5000>\n\
             ##contig=<ID=c2,length=1>\n##k=v\n##source=x y\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
        );
    }

    #[test]
    fn quotes_are_escaped_as_htsjdk_escapes_them() {
        assert_eq!(escape_quotes("a\"b"), "a\\\"b");
        assert_eq!(escape_quotes("a\\\"b"), "a\\\"b");
        assert_eq!(escape_quotes("\"a"), "\"a");
        assert_eq!(escape_quotes("a\"\"b"), "a\\\"\"b");
    }

    #[test]
    fn genotype_counts_follow_the_alleles() {
        assert_eq!(Count::G.values(2, 2), 3);
        assert_eq!(Count::G.values(3, 2), 6);
        assert_eq!(Count::G.values(2, 1), 2);
        assert_eq!(Count::A.values(3, 2), 2);
        assert_eq!(Count::R.values(3, 2), 3);
    }
}

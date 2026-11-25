//! The whole header of a VCF as htsjdk 5.0.0's `VCFCodec` reads it and its
//! `VCFWriter` writes it again, for the tools that write a VCF
//! ([`crate::vcf_writer`]).
//!
//! Reading (`AbstractVCFCodec.parseHeaderFromLines`, checked by
//! [`super::site_header::parse_without_contigs`]):
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
/// [`site_header::parse_without_contigs`]'s, made before.
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

/// One kept line as written, without its `##`: its bytes in
/// [`VcfHeader::text`], and `contig`, its place among the `##contig` lines
/// from 1, 0 for any other line.
#[derive(Clone, Copy, Debug)]
struct Written {
    contig: u32,
    start: u32,
    end: u32,
}

/// The header a VCF writer writes, with what its encoder looks up. Each
/// kept line is held once, as the text it is written as, in one buffer
/// with 12 bytes of offsets per line; the site rules hold no contig line
/// (the dictionary does). The raw lines are freed as they are read, so a
/// header of `n` lines of `b` bytes costs about `2b + 100n` bytes at its
/// peak (docs/chaining.md section 33).
#[derive(Clone, Debug)]
pub struct VcfHeader {
    text: String,
    lines: Vec<Written>,
    samples: Vec<String>,
    minor: u8,
    site: SiteHeader,
    infos: HashMap<String, Definition>,
    formats: HashMap<String, Definition>,
    filters: HashSet<String>,
    dictionary: Dictionary,
}

/// Keep the first contig of each name, in order: the dictionary takes the
/// first `##contig` line of each `ID`. A sort of 4-byte indexes finds the
/// repeats, so no second copy of the names is held.
fn first_of_each_name(dictionary: &mut Dictionary) {
    let mut order: Vec<u32> = (0..dictionary.len() as u32).collect();
    order.sort_by(|&a, &b| {
        dictionary[a as usize]
            .name
            .cmp(&dictionary[b as usize].name)
            .then(a.cmp(&b))
    });
    let mut keep = vec![true; dictionary.len()];
    for pair in order.windows(2) {
        if dictionary[pair[0] as usize].name == dictionary[pair[1] as usize].name {
            keep[pair[1] as usize] = false;
        }
    }
    let mut at = 0;
    dictionary.retain(|_| {
        at += 1;
        keep[at - 1]
    });
}

/// Java's `String.compareTo`: UTF-16 code units in order.
pub fn java_compare(a: &str, b: &str) -> Ordering {
    a.encode_utf16().cmp(b.encode_utf16())
}

impl VcfHeader {
    /// Read the header of `reader` as `VCFCodec` reads it, leaving the
    /// reader at the first record.
    pub fn read(reader: &mut impl BufRead) -> Result<Self, String> {
        Self::from_raw(site_header::read_raw(reader)?)
    }

    /// The header of `raw`, checked line by line; each raw line is freed
    /// once it is kept as written.
    pub fn from_raw(raw: RawHeader) -> Result<Self, String> {
        let site = site_header::parse_without_contigs(&raw)?;
        let RawHeader {
            lines: raw_lines,
            columns,
            minor,
        } = raw;
        // Sized from the lines read, so no buffer doubles past its need.
        let mut lines = Vec::with_capacity(raw_lines.len());
        let mut text = String::with_capacity(raw_lines.iter().map(String::len).sum());
        let mut dictionary = Vec::with_capacity(
            raw_lines
                .iter()
                .filter(|line| line.starts_with("##contig"))
                .count(),
        );
        let mut contig_lines = 0usize;
        let mut header = Self {
            text: String::new(),
            lines: Vec::new(),
            samples: column_samples(&columns)?,
            minor,
            site,
            infos: HashMap::new(),
            formats: HashMap::new(),
            filters: HashSet::new(),
            dictionary: Vec::new(),
        };
        for line in raw_lines {
            let line = line.as_str();
            let kept = if line.starts_with("##INFO") {
                Some(compound("INFO", after(line, 7)?, minor)?)
            } else if line.starts_with("##FORMAT") {
                Some(compound("FORMAT", after(line, 9)?, minor)?)
            } else if line.starts_with("##contig") {
                let fields = parse_fields(after(line, 9)?, None, &[])?;
                dictionary.push(Contig {
                    name: value(&fields, "ID").unwrap_or_default(),
                    length: value(&fields, "length").and_then(|l| l.parse().ok()),
                });
                contig_lines += 1;
                Some(Line::Structured {
                    key: "contig".to_owned(),
                    fields,
                    contig: Some(contig_lines),
                })
            } else if let Some((key, skip)) = structured_key(line, minor) {
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
            if let Some(line) = kept {
                header.index(&line);
                let contig = match &line {
                    Line::Structured {
                        contig: Some(rank), ..
                    } => u32::try_from(*rank).unwrap_or(u32::MAX),
                    _ => 0,
                };
                let start = text.len();
                text.push_str(&line.text());
                let offset = |at: usize| {
                    u32::try_from(at).map_err(|_| "the VCF header is over 4 GiB".to_owned())
                };
                lines.push(Written {
                    contig,
                    start: offset(start)?,
                    end: offset(text.len())?,
                });
            }
        }
        first_of_each_name(&mut dictionary);
        lines.shrink_to_fit();
        text.shrink_to_fit();
        dictionary.shrink_to_fit();
        header.text = text;
        header.lines = lines;
        header.dictionary = dictionary;
        Ok(header)
    }

    /// The text of a kept line, without its `##`.
    fn line(&self, line: Written) -> &str {
        &self.text[line.start as usize..line.end as usize]
    }

    /// The lookups of the encoder: the first line of each `ID`.
    fn index(&mut self, line: &Line) {
        match line {
            Line::Compound {
                key, id, count, ty, ..
            } => {
                let map = if *key == "INFO" {
                    &mut self.infos
                } else {
                    &mut self.formats
                };
                if !map.contains_key(id) {
                    map.insert(
                        id.clone(),
                        Definition {
                            count: *count,
                            flag: ty == "Flag",
                        },
                    );
                }
            }
            Line::Structured { key, fields, .. } if key == "FILTER" => {
                if let Some(id) = value(fields, "ID") {
                    self.filters.insert(id);
                }
            }
            _ => {}
        }
    }

    /// The sample names, each once, in column order.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// The same lines with `samples` as the sample columns:
    /// `new VCFHeader(getMetaDataInInputOrder(), samples)`.
    pub fn with_samples(&self, samples: Vec<String>) -> Self {
        let mut header = self.clone();
        header.set_samples(samples);
        header
    }

    /// [`VcfHeader::with_samples`] in place.
    pub fn set_samples(&mut self, samples: Vec<String>) {
        self.site.samples.clone_from(&samples);
        self.samples = samples;
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

    /// The header a reader of the file [`crate::vcf_writer::VcfWriter`]
    /// writes under `self` gets: the text as ISO-8859-1 bytes read back as
    /// UTF-8, `VCFv4.2`, the lines in written order. A chain hands this on
    /// after a transform, so the steps behind it read what they would read
    /// from the transform's file.
    ///
    /// The lines go straight to [`VcfHeader::from_raw`], each split where
    /// the reader splits it (`\n`, `\r`, `\r\n`), so no copy of the whole
    /// text is held beside them.
    pub fn reread(&self) -> Result<Self, String> {
        // The lines this header writes: its own, the version and the columns.
        let mut lines: Vec<String> = Vec::with_capacity(self.lines.len() + 2);
        let mut current = String::new();
        let mut after_cr = false;
        self.each_line(&mut |piece| {
            let piece = crate::java::iso_8859_1(piece);
            let piece = crate::java::decode_utf8(&piece);
            for c in piece.chars() {
                match c {
                    '\n' if after_cr => {}
                    '\n' | '\r' => {
                        // Each line at its own length: the buffer is reused.
                        lines.push(current.as_str().to_owned());
                        current.clear();
                    }
                    _ => current.push(c),
                }
                after_cr = c == '\r';
            }
            Ok(())
        })
        .map_err(|e| e.to_string())?;
        let Some(at) = lines.iter().position(|line| !line.starts_with("##")) else {
            return Err("the written header has no column line".to_owned());
        };
        lines.truncate(at + 1);
        let columns = lines.pop().unwrap_or_default();
        Self::from_raw(RawHeader {
            lines,
            columns,
            minor: 2,
        })
    }

    /// What a reader of the file [`crate::vcf_writer::VcfWriter`] writes
    /// under `self` reads: [`VcfHeader::reread`], with no second copy of
    /// the lines for a header of ASCII text. Writing such a header changes
    /// no character, so its kept lines read back as written; only the
    /// version (4.2) and what the reader takes from the first line of each
    /// `ID` in written order change: the `INFO`, `FORMAT` and `FILTER`
    /// lookups and the `END` flag of the site reader. Those are read again
    /// from the written lines.
    ///
    /// `None` when that is `self` (a `VCFv4.2` header whose lookups read
    /// back the same), so the caller keeps its copy.
    pub fn read_back(&self) -> Result<Option<Self>, String> {
        if !self.is_ascii() {
            return self.reread().map(Some);
        }
        let mut lookups = Self {
            text: String::new(),
            lines: Vec::new(),
            samples: Vec::new(),
            minor: 2,
            site: SiteHeader::default(),
            infos: HashMap::new(),
            formats: HashMap::new(),
            filters: HashSet::new(),
            dictionary: Vec::new(),
        };
        for written in self.written_lines() {
            let text = self.line(written);
            let parsed = if text.starts_with("INFO") || text.starts_with("FORMAT") {
                let line = format!("##{text}");
                if text.starts_with("INFO") {
                    compound("INFO", after(&line, 7)?, 2)?
                } else {
                    compound("FORMAT", after(&line, 9)?, 2)?
                }
            } else if text.starts_with("FILTER") {
                let line = format!("##{text}");
                let Some((key, skip)) = structured_key(&line, 2) else {
                    continue;
                };
                Line::Structured {
                    key: key.to_owned(),
                    fields: parse_fields(after(&line, skip)?, None, &[])?,
                    contig: None,
                }
            } else {
                continue;
            };
            lookups.index(&parsed);
        }
        let end_is_flag = lookups.info("END").is_some_and(|d| d.flag);
        if self.minor == 2
            && end_is_flag == self.site.end_is_flag
            && lookups.infos == self.infos
            && lookups.formats == self.formats
            && lookups.filters == self.filters
        {
            return Ok(None);
        }
        let mut header = self.clone();
        header.minor = 2;
        header.site.minor = 2;
        header.site.end_is_flag = end_is_flag;
        header.infos = lookups.infos;
        header.formats = lookups.formats;
        header.filters = lookups.filters;
        Ok(Some(header))
    }

    /// True when the header's text is ASCII, so writing it as ISO-8859-1
    /// and reading it back as UTF-8 changes no character.
    pub fn is_ascii(&self) -> bool {
        self.text.is_ascii() && self.samples.iter().all(|s| s.is_ascii())
    }

    /// `VCFHeader.getSequenceDictionary` builds a `SAMSequenceRecord` per
    /// `##contig` line, which checks its name: the first `ID` outside the
    /// SAM name pattern fails with htsjdk's message. `MakeSitesOnlyVcf`
    /// always asks for the dictionary, `SplitVcfs` without
    /// `SEQUENCE_DICTIONARY`.
    pub fn check_contig_names(&self) -> Result<(), String> {
        match self
            .dictionary
            .iter()
            .find(|contig| !crate::create_sequence_dictionary::legal_name(&contig.name))
        {
            Some(contig) => Err(format!(
                "Sequence name '{}' doesn't match regex: '{}'",
                contig.name,
                crate::create_sequence_dictionary::NAME_PATTERN
            )),
            None => Ok(()),
        }
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
        let mut out = String::new();
        self.each_line(&mut |line| {
            out.push_str(line);
            Ok(())
        })
        .unwrap_or_default();
        out
    }

    /// The kept lines in the order `VCFWriter.writeHeader` writes them,
    /// each once: contig lines by their place among the contig lines,
    /// every other pair by its text, two equal lines that are no contig
    /// written once (`TreeSet`). An index of 4 bytes per line is sorted.
    fn written_lines(&self) -> Vec<Written> {
        let mut order: Vec<u32> = (0..self.lines.len() as u32).collect();
        order.sort_by(|&a, &b| {
            let (a, b) = (self.lines[a as usize], self.lines[b as usize]);
            if a.contig > 0 && b.contig > 0 {
                a.contig.cmp(&b.contig)
            } else {
                java_compare(self.line(a), self.line(b))
            }
        });
        let mut written: Vec<Written> = Vec::with_capacity(order.len());
        for at in order {
            let line = self.lines[at as usize];
            if let Some(before) = written.last()
                && before.contig == 0
                && line.contig == 0
                && self.line(*before) == self.line(line)
            {
                continue;
            }
            written.push(line);
        }
        written
    }

    /// Hand `put` each piece of [`VcfHeader::text`] in order, so a writer
    /// never holds the whole text: the lines are sorted by an index of
    /// 4 bytes per line.
    pub fn each_line(
        &self,
        put: &mut dyn FnMut(&str) -> std::io::Result<()>,
    ) -> std::io::Result<()> {
        put("##fileformat=VCFv4.2\n")?;
        for line in self.written_lines() {
            put("##")?;
            put(self.line(line))?;
            put("\n")?;
        }
        put("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
        if !self.samples.is_empty() {
            put("\tFORMAT")?;
            for sample in &self.samples {
                put("\t")?;
                put(sample)?;
            }
        }
        put("\n")
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

    #[test]
    fn read_back_equals_a_reread() {
        let body = "##source=x y\n##contig=<ID=c2,length=3000>\n##contig=<ID=c1,length=5000>\n\
             ##contig=<ID=c2,length=1>\n\
             ##INFO=<ID=XT,Number=1,Type=String,Description=\"a\",Foo=bar,Source=s,Version=3>\n\
             ##INFO=<ID=NODESC,Number=1,Type=Integer>\n\
             ##INFO=<ID=END,Number=1,Type=Flag,Description=\"e\">\n\
             ##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filter\">\n\
             ##FORMAT=<ID=ZZ,Number=1,Type=String,Description=\"z\">\n\
             ##FORMAT=<ID=ZZ,Number=A,Type=Integer,Description=\"a\">\n\
             ##FILTER=<ID=q10,Description=\"low, \\\"q\\\"\">\n\
             ##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"mq\">\n\
             ##META=<ID=Assay,Values=[A, B]>\n##PEDIGREE=<ID=c,Father=f,Mother=m>\n\
             ##SAMPLE=<ID=S1,Assay=WholeGenome>\n##ALT=<ID=DEL,Description=\"d\">\n\
             ##INFO=<ID=AA,Number=A,Type=Flag,Description=\"a\">\n\
             ##INFO=<ID=AA,Number=2,Type=Integer,Description=\"z\">\n\
             ##FILTER=<ID=q10,Description=\"again\">\n\
             ##noequals\n##k=v\n##k=v\n";
        for version in 0..=3 {
            for body in [body, ""] {
                let read = header(&format!("##fileformat=VCFv4.{version}\n{body}{COLUMNS}"));
                let fast = read
                    .read_back()
                    .expect("reads back")
                    .unwrap_or_else(|| read.clone());
                let reread = read.reread().expect("reads back");
                assert_eq!(fast.text(), reread.text(), "4.{version}");
                assert_eq!(fast.site(), reread.site(), "4.{version}");
                assert_eq!(fast.minor(), reread.minor());
                assert_eq!(fast.samples(), reread.samples());
                assert_eq!(fast.dictionary(), reread.dictionary());
                assert_eq!(fast.filters, reread.filters, "4.{version}");
                let keys: Vec<&String> = fast
                    .infos
                    .keys()
                    .chain(fast.formats.keys())
                    .chain(reread.infos.keys())
                    .chain(reread.formats.keys())
                    .collect();
                for key in keys {
                    assert_eq!(fast.info(key), reread.info(key), "4.{version} {key}");
                    assert_eq!(fast.format(key), reread.format(key), "4.{version} {key}");
                }
            }
        }
    }

    #[test]
    fn a_header_read_back_is_at_version_4_2_with_its_text_as_written() {
        let text = format!(
            "##fileformat=VCFv4.0\n##INFO=<ID=END,Number=1,Type=Flag,Description=\"e\">\n\
             ##INFO=<ID=AA,Number=A,Type=Flag,Description=\"a\">\n\
             ##INFO=<ID=AA,Number=2,Type=Integer,Description=\"z\">\n{COLUMNS}"
        );
        let read = header(&text);
        let reread = read.reread().expect("reads back");
        assert_eq!(reread.minor(), 2);
        assert_eq!(reread.text(), read.text());
        // The first AA line of the written order is the one a reader keeps.
        assert_eq!(read.info("AA").map(|d| d.flag), Some(true));
        assert_eq!(reread.info("AA").map(|d| d.flag), Some(false));
        assert!(!reread.site().end_is_flag);
        // Text a VCF cannot hold as ISO-8859-1 changes on the way back.
        let text = format!("##fileformat=VCFv4.2\n##source=caller \u{e9} \u{65e5}\n{COLUMNS}");
        let read = header(&text);
        assert!(!read.is_ascii());
        assert!(
            read.reread()
                .expect("reads back")
                .text()
                .contains("##source=caller \u{fffd} ?\n")
        );
    }
}

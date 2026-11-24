//! M20: the VCF layer the variant tools share, as htsjdk 5.0.0 reads a
//! file for Picard 3.5.0.
//!
//! [`VcfReader`] opens a plain or BGZF-compressed VCF (detected by the gzip
//! magic bytes, whatever the extension) through noodles and turns each line
//! into a [`Variant`]: the alleles as `VCFCodec.parseAlleles` builds them
//! (reference first, a `.` alternate dropped, bases upper-cased), the record
//! end from `INFO/END` or the reference length, the site filter status, and
//! one [`Genotype`] per sample as `AbstractVCFCodec.createGenotypeMap`
//! decodes `GT`, `GQ`, `DP`, `AD` and `FT`. The lazy noodles record gives the
//! raw text of each field; the field rules here are the codec's, because the
//! metrics depend on them (a `GQ` of `20.5` rounds to 21, an `AD` with a
//! non-integer is null, `FT=PASS` is unfiltered).
//!
//! [`VariantType`] is `VariantContext.determineType`, [`GenotypeType`] is
//! `Genotype.determineType`, and [`Variant::sub_context`] is
//! `VariantContext.subContextFromSample`, which re-derives the allele list
//! from the one genotype and so can turn a SNP record into `NO_VARIATION`.
//!
//! [`DbSnpSets`] holds `DbSnpBitSetUtil.createSnpAndIndelBitSets`: one bit
//! per reference position, set for every position a SNP (or an indel)
//! record spans, in pages made only where a bit is set, at most
//! [`MAX_DBSNP_PAGES`]. A dbSNP or dictionary VCF is read under
//! [`MAX_LINE_BYTES`] per line and [`MAX_DBSNP_LINES`] lines. [`Regions`] is the interval list support the VCF tools need
//! (`IntervalList.fromPath`, `uniqued`, `intersection`, `union`, base count
//! and overlap queries by contig name).

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};

use noodles::vcf;

mod bit_set;
mod encode;
mod encode_genotypes;
mod full_header;
mod header;
mod line_cap;
mod site;
mod site_header;
mod subset;

pub use bit_set::MAX_DBSNP_PAGES;
use bit_set::PAGE_BITS;
use bit_set::{BitSet, OverBudget};
pub use encode::{Encoder, format_quality};
pub use full_header::{Count, Definition, VcfHeader, java_compare};
use line_cap::LineCap;
pub use site_header::{ContigLine, MAX_HEADER_BYTES, RawHeader, SiteHeader};
pub use subset::{SampleSubset, tree_set};

use crate::bed_to_interval_list;
use crate::error::{Error, Result};

/// One `@SQ` or `##contig` entry: a name and a length (a VCF contig line may
/// omit the length).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Contig {
    pub name: String,
    pub length: Option<i64>,
}

/// A sequence dictionary in file order.
pub type Dictionary = Vec<Contig>;

/// `SequenceDictionaryUtils.assertSequenceDictionariesEqual`: both absent is
/// fine, one absent is Picard's `NullPointerException`, both present must
/// list the same names and lengths in the same order.
pub fn assert_dictionaries_equal(
    first: Option<&Dictionary>,
    first_source: &str,
    second: Option<&Dictionary>,
    second_source: &str,
) -> Result<()> {
    let differ = || {
        Error::Usage(format!(
            "Sequence dictionary for ({first_source}) does not match sequence dictionary for ({second_source})"
        ))
    };
    match (first, second) {
        (None, None) => Ok(()),
        (Some(a), Some(b)) => {
            if a.len() != b.len() {
                return Err(differ());
            }
            for (x, y) in a.iter().zip(b) {
                let same_length = match (x.length, y.length) {
                    (Some(p), Some(q)) => p == q,
                    _ => true,
                };
                if x.name != y.name || !same_length {
                    return Err(differ());
                }
            }
            Ok(())
        }
        _ => Err(Error::Usage(format!(
            "{}: no sequence dictionary (Picard fails with a NullPointerException comparing it with {})",
            if first.is_none() {
                first_source
            } else {
                second_source
            },
            if first.is_none() {
                second_source
            } else {
                first_source
            }
        ))),
    }
}

/// `SAMSequenceDictionaryExtractor.extractDictionary`: a VCF's `##contig`
/// lines, or the `@SQ` lines of a `.dict`, BAM, SAM, interval list or FASTA
/// (its `.dict`, else the `.fai`).
pub fn read_dictionary(path: &Path) -> Result<Dictionary> {
    let name = path.to_string_lossy().to_ascii_lowercase();
    if [".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"]
        .iter()
        .any(|e| name.ends_with(e))
    {
        let reader = VcfReader::open_capped(path, MAX_LINE_BYTES, MAX_DBSNP_LINES)?;
        return reader.dictionary().ok_or_else(|| Error::Intervals {
            path: path.to_path_buf(),
            message: "no ##contig line: cannot extract a sequence dictionary".to_owned(),
        });
    }
    Ok(bed_to_interval_list::read_dictionary(path)?
        .into_iter()
        .map(|s| Contig {
            name: s.name,
            length: Some(s.length),
        })
        .collect())
}

/// One allele of a record: `ACGTN` bases (upper case), or a symbolic allele
/// such as `<NON_REF>`, kept as written.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Allele {
    pub bases: String,
    pub symbolic: bool,
}

impl Allele {
    /// `Allele.create`: symbolic when bracketed or a breakend
    /// (`Allele.wouldBeSymbolicAllele`), else upper-cased bases.
    fn parse(text: &str, is_ref: bool) -> std::result::Result<Self, String> {
        if text.is_empty() {
            return Err("empty alleles are not permitted in VCF records".to_owned());
        }
        let bytes = text.as_bytes();
        let symbolic = bytes.len() > 1
            && (bytes[0] == b'<'
                || bytes[bytes.len() - 1] == b'>'
                || bytes.contains(&b'[')
                || bytes.contains(&b']')
                || (bytes[0] == b'.' || bytes[bytes.len() - 1] == b'.'));
        if symbolic {
            if is_ref {
                return Err(format!(
                    "Symbolic alleles not allowed as reference allele: {text}"
                ));
            }
            return Ok(Self {
                bases: text.to_owned(),
                symbolic: true,
            });
        }
        if text == "*" {
            if is_ref {
                return Err("The reference allele cannot be a spanning deletion".to_owned());
            }
            return Ok(Self {
                bases: text.to_owned(),
                symbolic: false,
            });
        }
        if !bytes
            .iter()
            .all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T' | b'N'))
        {
            return Err(format!("unparsable vcf record with allele {text}"));
        }
        Ok(Self {
            bases: text.to_ascii_uppercase(),
            symbolic: false,
        })
    }

    pub fn len(&self) -> usize {
        self.bases.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bases.is_empty()
    }
}

/// `VariantContext.Type`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VariantType {
    NoVariation,
    Snp,
    Mnp,
    Indel,
    Symbolic,
    Mixed,
}

impl VariantType {
    /// The name `MetricBase` prints for the enum.
    pub fn as_str(self) -> &'static str {
        match self {
            VariantType::NoVariation => "NO_VARIATION",
            VariantType::Snp => "SNP",
            VariantType::Mnp => "MNP",
            VariantType::Indel => "INDEL",
            VariantType::Symbolic => "SYMBOLIC",
            VariantType::Mixed => "MIXED",
        }
    }
}

/// `GenotypeType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenotypeType {
    NoCall,
    HomRef,
    Het,
    HomVar,
    Unavailable,
    Mixed,
}

/// One sample's genotype as the codec builds it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Genotype {
    /// One entry per `GT` allele: the index into [`Variant::alleles`], or
    /// `None` for `.`. Empty when the record has no `GT` key.
    pub alleles: Vec<Option<usize>>,
    /// `GQ`, rounded as `Math.round` does; -1 when absent or `.`.
    pub gq: i64,
    /// `DP`; -1 when absent or `.`.
    pub dp: i64,
    /// `AD`; `None` when absent, `.` or not all integers.
    pub ad: Option<Vec<i64>>,
    /// `FT` present and neither `.` nor `PASS`.
    pub filtered: bool,
}

impl Genotype {
    /// `Genotype.determineType`.
    pub fn kind(&self) -> GenotypeType {
        if self.alleles.is_empty() {
            return GenotypeType::Unavailable;
        }
        let mut saw_no_call = false;
        let mut first: Option<usize> = None;
        let mut multiple = false;
        for allele in &self.alleles {
            match (*allele, first) {
                (None, _) => saw_no_call = true,
                (Some(a), None) => first = Some(a),
                (Some(a), Some(f)) => multiple |= a != f,
            }
        }
        match (saw_no_call, first) {
            (true, None) => GenotypeType::NoCall,
            (true, Some(_)) => GenotypeType::Mixed,
            (false, None) => GenotypeType::Unavailable,
            (false, Some(f)) => {
                if multiple {
                    GenotypeType::Het
                } else if f == 0 {
                    GenotypeType::HomRef
                } else {
                    GenotypeType::HomVar
                }
            }
        }
    }

    pub fn is_hom_ref(&self) -> bool {
        self.kind() == GenotypeType::HomRef
    }

    pub fn is_hom_var(&self) -> bool {
        self.kind() == GenotypeType::HomVar
    }

    pub fn is_het(&self) -> bool {
        self.kind() == GenotypeType::Het
    }

    pub fn is_no_call(&self) -> bool {
        self.kind() == GenotypeType::NoCall
    }

    pub fn is_mixed(&self) -> bool {
        self.kind() == GenotypeType::Mixed
    }
}

/// One VCF record with its genotypes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Variant {
    pub contig: String,
    /// 1-based.
    pub start: i64,
    /// `INFO/END`, else `start + reference length - 1`.
    pub end: i64,
    /// The reference allele first, then the alternates in file order.
    pub alleles: Vec<Allele>,
    /// `FILTER` is neither `.` nor `PASS`.
    pub filtered: bool,
    /// One per sample, in header order.
    pub genotypes: Vec<Genotype>,
}

impl Variant {
    pub fn reference(&self) -> &Allele {
        &self.alleles[0]
    }

    pub fn alternates(&self) -> &[Allele] {
        &self.alleles[1..]
    }

    /// `VariantContext.determineType`: one allele is `NO_VARIATION`; each
    /// alternate against the reference is symbolic, a SNP or MNP (same
    /// length) or an indel; two alternates of different kinds make `MIXED`.
    pub fn kind(&self) -> VariantType {
        let reference = self.reference();
        let mut kind = None;
        for alt in self.alternates() {
            let this = if alt.symbolic {
                VariantType::Symbolic
            } else if reference.len() == alt.len() {
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

    pub fn is_biallelic(&self) -> bool {
        self.alleles.len() == 2
    }

    /// `VariantContext.isSimpleIndel`: a biallelic indel whose alleles
    /// share their first base and one of which is one base long.
    fn is_simple_indel(&self) -> bool {
        let (reference, alt) = (self.reference(), &self.alleles[1]);
        self.kind() == VariantType::Indel
            && self.is_biallelic()
            && !reference.is_empty()
            && !alt.is_empty()
            && reference.bases.as_bytes()[0] == alt.bases.as_bytes()[0]
            && (reference.len() == 1 || alt.len() == 1)
    }

    pub fn is_simple_insertion(&self) -> bool {
        self.is_simple_indel() && self.reference().len() == 1
    }

    pub fn is_simple_deletion(&self) -> bool {
        self.is_simple_indel() && self.alleles[1].len() == 1
    }

    pub fn is_complex_indel(&self) -> bool {
        self.kind() == VariantType::Indel
            && !self.is_simple_deletion()
            && !self.is_simple_insertion()
    }

    /// `VariantContextUtils.isTransition` for a biallelic SNP: A<>G or C<>T.
    pub fn is_transition(&self) -> bool {
        let r = self.reference().bases.as_bytes()[0];
        let a = self.alleles[1].bases.as_bytes()[0];
        matches!(
            (r, a),
            (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C')
        )
    }

    /// `VariantContext.subContextFromSample`: the record with only this
    /// sample's genotype and the alleles re-derived from it (the reference
    /// always, then every called allele, in the record's order).
    pub fn sub_context(&self, sample: usize) -> Variant {
        let genotype = &self.genotypes[sample];
        let mut keep = vec![false; self.alleles.len()];
        keep[0] = true;
        for allele in genotype.alleles.iter().flatten() {
            keep[*allele] = true;
        }
        let mut new_index = vec![None; self.alleles.len()];
        let mut alleles = Vec::new();
        for (i, allele) in self.alleles.iter().enumerate() {
            if keep[i] {
                new_index[i] = Some(alleles.len());
                alleles.push(allele.clone());
            }
        }
        let mut sub = genotype.clone();
        sub.alleles = genotype
            .alleles
            .iter()
            .map(|a| a.and_then(|i| new_index[i]))
            .collect();
        Variant {
            contig: self.contig.clone(),
            start: self.start,
            end: self.end,
            alleles,
            filtered: self.filtered,
            genotypes: vec![sub],
        }
    }
}

/// A VCF file open for reading, plain or BGZF.
pub struct VcfReader {
    path: PathBuf,
    header: vcf::Header,
    reader: vcf::io::Reader<Box<dyn BufRead + Send>>,
    record: vcf::Record,
    line: u64,
}

/// `path` as text lines under the caps, more than `max_header_lines`
/// lines starting with `#` included: plain, or gzip when its first two
/// bytes are the gzip magic, as `GZIPInputStream` reads a `.vcf.gz` for
/// htsjdk. A gzip member whose header carries the BGZF `BC` extra field is
/// read as BGZF; any other gzip member by member. One open: the first
/// bytes are handed back in front of the rest, so a named pipe is read
/// once.
fn open_text(
    path: &Path,
    max_line_bytes: usize,
    max_lines: u64,
    max_header_lines: u64,
) -> Result<Box<dyn BufRead + Send>> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let mut file = File::open(path).map_err(input)?;
    let mut magic = [0u8; 14];
    let mut read = 0;
    while read < magic.len() {
        match file.read(&mut magic[read..]) {
            Ok(0) => break,
            Ok(n) => read += n,
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => return Err(input(error)),
        }
    }
    let gzip = read >= 2 && magic[..2] == [0x1f, 0x8b];
    let bgzf = read == magic.len() && magic[3] & 0x04 != 0 && magic[12..14] == *b"BC";
    let file = io::Cursor::new(magic[..read].to_vec()).chain(file);
    let text: Box<dyn BufRead + Send> = if gzip && bgzf {
        Box::new(noodles::bgzf::io::Reader::new(file))
    } else if gzip {
        Box::new(io::BufReader::new(crate::gzip::GzipReader::new(
            io::BufReader::new(file),
        )))
    } else {
        Box::new(io::BufReader::new(file))
    };
    Ok(Box::new(
        LineCap::new(text, max_line_bytes, max_lines).with_header_lines(max_header_lines),
    ))
}

impl VcfReader {
    pub fn open(path: &Path) -> Result<Self> {
        Self::open_capped(path, usize::MAX, u64::MAX)
    }

    /// [`VcfReader::open`] where a line (header or record) over
    /// `max_line_bytes`, or more than `max_lines` lines, is a read error.
    pub fn open_capped(path: &Path, max_line_bytes: usize, max_lines: u64) -> Result<Self> {
        Self::open_capped_header(path, max_line_bytes, max_lines, u64::MAX)
    }

    /// [`VcfReader::open_capped`] where more than `max_header_lines` lines
    /// starting with `#` is a read error too, so a header of millions of
    /// `##contig` lines is refused before noodles holds it.
    pub fn open_capped_header(
        path: &Path,
        max_line_bytes: usize,
        max_lines: u64,
        max_header_lines: u64,
    ) -> Result<Self> {
        let input = |source| Error::Input {
            path: path.to_path_buf(),
            source,
        };
        let capped = open_text(path, max_line_bytes, max_lines, max_header_lines)?;
        let mut reader = vcf::io::Reader::new(capped);
        let header = header::read_header(reader.get_mut()).map_err(input)?;
        Ok(Self {
            path: path.to_path_buf(),
            header,
            reader,
            record: vcf::Record::default(),
            line: 0,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Sample names in header order (`VCFHeader.getGenotypeSamples`).
    pub fn samples(&self) -> Vec<String> {
        self.header.sample_names().iter().cloned().collect()
    }

    /// The `##contig` lines, `None` without any (`getSequenceDictionary`).
    pub fn dictionary(&self) -> Option<Dictionary> {
        let contigs = self.header.contigs();
        if contigs.is_empty() {
            return None;
        }
        Some(
            contigs
                .iter()
                .map(|(name, contig)| Contig {
                    name: name.clone(),
                    length: contig.length().map(|l| l as i64),
                })
                .collect(),
        )
    }

    /// The next record, `None` at the end of the file.
    pub fn next_variant(&mut self) -> Result<Option<Variant>> {
        let read = self
            .reader
            .read_record(&mut self.record)
            .map_err(|source| Error::Input {
                path: self.path.clone(),
                source,
            })?;
        if read == 0 {
            return Ok(None);
        }
        self.line += 1;
        let samples = self.header.sample_names().len();
        parse_variant(&self.record, samples)
            .map(Some)
            .map_err(|message| Error::Intervals {
                path: self.path.clone(),
                message: format!("record {}: {message}", self.line),
            })
    }

    /// The header as noodles read it.
    pub fn header(&self) -> &vcf::Header {
        &self.header
    }

    /// The raw record [`VcfReader::next_variant`] read last, for the
    /// fields a [`Variant`] does not keep (each `FILTER` name, a `FORMAT`
    /// value such as `GTA`).
    pub fn record(&self) -> &vcf::Record {
        &self.record
    }

    /// Every remaining record.
    pub fn variants(mut self) -> Result<Vec<Variant>> {
        let mut out = Vec::new();
        while let Some(v) = self.next_variant()? {
            out.push(v);
        }
        Ok(out)
    }
}

/// A VCF read for its sites as htsjdk's `VCFCodec` reads them for
/// `VcfToIntervalList`: the header through `vcf/site_header.rs`, then each
/// record line through `vcf/site.rs`, without noodles.
pub struct SiteReader {
    path: PathBuf,
    reader: Box<dyn BufRead + Send>,
    header: SiteHeader,
    rules: site::Rules,
    line: u64,
    /// The line read last, as bytes.
    raw: Vec<u8>,
}

impl SiteReader {
    /// Open `path` (plain or gzip by its first bytes) and read its header,
    /// a line over `max_line_bytes`, more than `max_lines` lines or more
    /// than `max_header_lines` header lines being a read error.
    pub fn open(
        path: &Path,
        max_line_bytes: usize,
        max_lines: u64,
        max_header_lines: u64,
    ) -> Result<Self> {
        let mut reader = open_text(path, max_line_bytes, max_lines, max_header_lines)?;
        let header = site_header::read(&mut reader).map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message,
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
            rules: site::Rules::of(&header),
            header,
            line: 0,
            raw: Vec::new(),
        })
    }

    pub fn header(&self) -> &SiteHeader {
        &self.header
    }

    /// The next record's position, `ID` and filter status, read from its
    /// line as htsjdk's `VCFCodec` reads it for a site (`vcf/site.rs`): the
    /// fixed columns checked, the genotypes left undecoded unless the
    /// header's samples are not in sorted order; a line starting with `#`
    /// skipped; `None` at the end of the file. A line is decoded as UTF-8,
    /// a byte that is not becoming U+FFFD.
    pub fn next_site(&mut self) -> Result<Option<Site>> {
        loop {
            let read = line_cap::read_line(&mut self.reader, &mut self.raw).map_err(|source| {
                Error::Input {
                    path: self.path.clone(),
                    source,
                }
            })?;
            if read == 0 {
                return Ok(None);
            }
            self.line += 1;
            let text = crate::java::decode_utf8(&self.raw);
            match site::parse_line(&text, &self.rules) {
                Ok(Some(site)) => return Ok(Some(site)),
                Ok(None) => {}
                Err(message) => {
                    return Err(Error::Intervals {
                        path: self.path.clone(),
                        message: format!("record {}: {message}", self.line),
                    });
                }
            }
        }
    }
}

/// A VCF read line by line for a writer: the whole header as htsjdk's
/// `VCFCodec` reads it ([`VcfHeader`]), then each record line as text, at
/// most [`MAX_LINE_BYTES`] each, one line in memory at a time.
pub struct VcfLines {
    path: PathBuf,
    reader: Box<dyn BufRead + Send>,
    header: std::sync::Arc<VcfHeader>,
    raw: Vec<u8>,
    line: u64,
}

impl VcfLines {
    /// Open `path` (plain, gzip or BGZF by its first bytes) and read its
    /// header.
    pub fn open(path: &Path) -> Result<Self> {
        let mut reader = open_text(path, MAX_LINE_BYTES, u64::MAX, u64::MAX)?;
        let header = VcfHeader::read(&mut reader).map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message,
        })?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
            header: std::sync::Arc::new(header),
            raw: Vec::new(),
            line: 0,
        })
    }

    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// The header, one copy shared by every step that holds it.
    pub fn shared_header(&self) -> std::sync::Arc<VcfHeader> {
        std::sync::Arc::clone(&self.header)
    }

    /// The header as the reader holds it, for a caller that rewrites it in
    /// place ([`SampleSubset::new_shared`]); the lines read after are not
    /// read by it.
    pub fn header_mut(&mut self) -> &mut std::sync::Arc<VcfHeader> {
        &mut self.header
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The record lines read so far.
    pub fn lines(&self) -> u64 {
        self.line
    }

    /// Replace `out` by the next line without its line end (`\n`, `\r` or
    /// `\r\n`), decoded as UTF-8 (a byte that is not becomes U+FFFD, as
    /// Java decodes it). False at the end of the file.
    pub fn next_line(&mut self, out: &mut String) -> Result<bool> {
        let read = line_cap::read_line(&mut self.reader, &mut self.raw).map_err(|source| {
            Error::Input {
                path: self.path.clone(),
                source,
            }
        })?;
        if read == 0 {
            return Ok(false);
        }
        self.line += 1;
        out.clear();
        out.push_str(&crate::java::decode_utf8(&self.raw));
        Ok(true)
    }
}

/// One record line (no line end) decoded into a [`Variant`] as
/// [`VcfReader::next_variant`] decodes a record of a file whose header
/// names `samples` samples: the same noodles fields and the same rules.
pub fn variant_from_line(line: &str, samples: usize) -> std::result::Result<Variant, String> {
    let mut reader = vcf::io::Reader::new(line.as_bytes());
    let mut record = vcf::Record::default();
    reader.read_record(&mut record).map_err(|e| e.to_string())?;
    parse_variant(&record, samples)
}

/// The first bytes of the text of the regular file at `path` (plain, gzip
/// or BGZF), at most 4; `None` for a pipe, a device or a file that cannot
/// be read. A pipe or a device is never opened here, so it is never read
/// ahead of its reader.
fn text_start(path: &Path) -> Option<Vec<u8>> {
    if !fs::metadata(path).is_ok_and(|m| m.is_file()) || path.starts_with("/dev") {
        return None;
    }
    let mut text = open_text(path, MAX_LINE_BYTES, u64::MAX, u64::MAX).ok()?;
    let mut start = [0u8; 4];
    let mut filled = 0;
    while filled < start.len() {
        match text.read(&mut start[filled..]) {
            Ok(0) | Err(_) => break,
            Ok(n) => filled += n,
        }
    }
    Some(start[..filled].to_vec())
}

/// True when `path` is a regular file whose text (plain, gzip or BGZF)
/// starts with `#`: a VCF to the chain verb. htsjdk reads the `##` lines
/// in any order, so `##fileformat` need not be first; a text that names no
/// version fails when its header is read, with the codec's message.
pub fn is_vcf(path: &Path) -> bool {
    text_start(path).is_some_and(|start| start.first() == Some(&b'#'))
}

/// Why the file at `path`, which [`is_vcf`] refuses, is no VCF.
pub fn why_not_vcf(path: &Path) -> String {
    match text_start(path) {
        None => "it is not a regular file that can be read".to_owned(),
        Some(start) if start.is_empty() => "it is empty".to_owned(),
        Some(start) if start.starts_with(b"BCF") => {
            "it is a BCF, and bammetrics reads VCF text only".to_owned()
        }
        Some(_) => "its text does not start with #, as a VCF header does".to_owned(),
    }
}

/// One record as [`SiteReader::next_site`] reads it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Site {
    pub contig: String,
    /// `POS`, a Java `int`.
    pub start: i32,
    /// `INFO/END`, else `start + reference length - 1`.
    pub end: i32,
    /// `ID`, `None` for `.`.
    pub id: Option<String>,
    /// `FILTER` is neither `.` nor `PASS`.
    pub filtered: bool,
}

/// `VCFCodec.parseFilters`: `.` and `PASS` are unfiltered. noodles hands
/// a missing (`.`) `FILTER` column over as an empty string.
fn is_filtered(text: &str) -> bool {
    !matches!(text, "" | "." | "PASS")
}

fn parse_variant(record: &vcf::Record, samples: usize) -> std::result::Result<Variant, String> {
    let contig = record.reference_sequence_name().to_owned();
    let start = match record.variant_start() {
        Some(Ok(p)) => usize::from(p) as i64,
        Some(Err(e)) => return Err(format!("bad position: {e}")),
        // noodles reads `POS` 0 (a telomere) as no position; htsjdk keeps 0.
        None => 0,
    };
    let reference = Allele::parse(record.reference_bases(), true)?;
    let mut alleles = vec![reference];
    let alternates = record.alternate_bases();
    let alts: &str = alternates.as_ref();
    // noodles hands a missing (`.`) `ALT` column over as an empty string;
    // a `.` among the alternates is `Allele.NO_CALL`, which the codec drops.
    for alt in alts.split(',') {
        if alt.is_empty() || alt == "." {
            continue;
        }
        alleles.push(Allele::parse(alt, false)?);
    }
    let filtered = is_filtered(record.filters().as_ref());
    let info_field = record.info();
    let info: &str = info_field.as_ref();
    let mut end = start + alleles[0].len() as i64 - 1;
    for field in info.split(';') {
        if let Some(value) = field.strip_prefix("END=") {
            end = value
                .parse()
                .map_err(|_| "the END value in the INFO field is not valid".to_owned())?;
        }
    }
    let raw = record.samples();
    let keys: Vec<&str> = raw.keys().iter().collect();
    let mut genotypes = Vec::with_capacity(samples);
    for sample in raw.iter() {
        genotypes.push(parse_genotype(&keys, sample.as_ref(), alleles.len())?);
    }
    if genotypes.len() != samples {
        return Err(format!(
            "there are {} genotypes while the header requires that {samples} genotypes be present for all records",
            genotypes.len()
        ));
    }
    Ok(Variant {
        contig,
        start,
        end,
        alleles,
        filtered,
        genotypes,
    })
}

/// `AbstractVCFCodec.createGenotypeMap` for one sample.
fn parse_genotype(
    keys: &[&str],
    text: &str,
    n_alleles: usize,
) -> std::result::Result<Genotype, String> {
    let values: Vec<&str> = text.split(':').collect();
    if values.len() > keys.len() {
        return Err(format!(
            "There are too many keys for the sample, keys = {}, values = {text}",
            keys.join(":")
        ));
    }
    let mut genotype = Genotype {
        alleles: Vec::new(),
        gq: -1,
        dp: -1,
        ad: None,
        filtered: false,
    };
    for (i, key) in keys.iter().enumerate() {
        let value = values.get(i).copied();
        if *key == "GT" {
            if i > 0 {
                return Err(format!(
                    "Saw GT field at position {i}, but it must be at the first position for genotypes when present"
                ));
            }
            let gt = value.unwrap_or(".");
            // `StringTokenizer` over `/|\`: empty tokens are skipped.
            for token in gt.split(['/', '|', '\\']).filter(|t| !t.is_empty()) {
                if token == "." {
                    genotype.alleles.push(None);
                } else {
                    let index: usize = token
                        .parse()
                        .map_err(|_| format!("Invalid genotype allele {token:?} in {gt}"))?;
                    if index >= n_alleles {
                        return Err(format!(
                            "Invalid genotype index {index} in {gt}: the record has {n_alleles} alleles"
                        ));
                    }
                    genotype.alleles.push(Some(index));
                }
            }
            continue;
        }
        let Some(value) = value else { continue };
        if *key == "FT" {
            if value.is_empty() {
                return Err("The VCF specification requires a valid filter status".to_owned());
            }
            genotype.filtered = is_filtered(value);
            continue;
        }
        if value == "." {
            continue;
        }
        match *key {
            "GQ" => {
                let q: f64 = value
                    .parse()
                    .map_err(|_| format!("GQ {value:?} is not a number"))?;
                genotype.gq = (q + 0.5).floor() as i64;
            }
            "DP" => {
                genotype.dp = value
                    .parse()
                    .map_err(|_| format!("DP {value:?} is not an integer"))?;
            }
            "AD" => {
                genotype.ad = value
                    .split(',')
                    .map(|v| v.parse::<i64>().ok())
                    .collect::<Option<Vec<i64>>>();
            }
            _ => {}
        }
    }
    Ok(genotype)
}

/// The longest line of a dbSNP or dictionary VCF, in bytes.
pub const MAX_LINE_BYTES: usize = 16 << 20;
/// The most lines of a dbSNP or dictionary VCF (dbSNP 156 has about 1.2
/// billion records).
pub const MAX_DBSNP_LINES: u64 = 4_000_000_000;

/// The error for bit sets over [`MAX_DBSNP_PAGES`].
fn over_budget(path: &Path) -> Error {
    Error::Intervals {
        path: path.to_path_buf(),
        message: format!(
            "the sites set over {MAX_DBSNP_PAGES} pages of {PAGE_BITS} positions (1 GiB of bit \
             set), the most a dbSNP file may set"
        ),
    }
}

/// `DbSnpBitSetUtil`: one bit set per contig.
#[derive(Clone, Debug, Default)]
pub struct DbSnpBitSet {
    by_contig: HashMap<String, BitSet>,
}

impl DbSnpBitSet {
    /// `isDbSnpSite`.
    pub fn is_site(&self, contig: &str, position: i64) -> bool {
        self.by_contig
            .get(contig)
            .is_some_and(|bits| bits.get(position))
    }
}

/// `DbSnpBitSetUtil.createSnpAndIndelBitSets`: SNP records set the SNP
/// bits, indel records the indel bits, each from the record start to its
/// end; every other record type sets nothing but still creates the contig's
/// (empty) set.
#[derive(Clone, Debug, Default)]
pub struct DbSnpSets {
    pub snps: DbSnpBitSet,
    pub indels: DbSnpBitSet,
}

impl DbSnpSets {
    /// Load `path`, whose `##contig` lines must equal `dictionary`
    /// (`assertSequenceDictionariesEqual`, the input's dictionary against
    /// the dbSNP's). With `regions`, only records overlapping a region
    /// count (`ByIntervalListVariantContextIterator`).
    pub fn load(
        path: &Path,
        dictionary: Option<&Dictionary>,
        regions: Option<&Regions>,
    ) -> Result<Self> {
        Self::load_with(path, dictionary, regions, true)
    }

    /// [`DbSnpSets::load`] where a dbSNP file without `##contig` lines is
    /// read without the dictionary check (`CollectArraysVariantCallingMetrics`
    /// completes on one).
    pub fn load_lenient(path: &Path, dictionary: &Dictionary) -> Result<Self> {
        Self::load_with(path, Some(dictionary), None, false)
    }

    fn load_with(
        path: &Path,
        dictionary: Option<&Dictionary>,
        regions: Option<&Regions>,
        needs_dictionary: bool,
    ) -> Result<Self> {
        let mut sets = Self::default();
        let mut budget = MAX_DBSNP_PAGES;
        load_db_snp(path, dictionary, regions, needs_dictionary, |variant| {
            let kind = variant.kind();
            sets.snps
                .add(variant, kind == VariantType::Snp, &mut budget)?;
            sets.indels
                .add(variant, kind == VariantType::Indel, &mut budget)
        })
        .map_err(|e| e.into_error(path))?;
        Ok(sets)
    }
}

impl DbSnpBitSet {
    /// `new DbSnpBitSetUtil(file, dictionary, EnumSet.noneOf(VariantType.class),
    /// intervals)`: every record sets its bits, whatever its type. The
    /// dictionary and region rules are those of [`DbSnpSets::load`].
    pub fn load_all(
        path: &Path,
        dictionary: Option<&Dictionary>,
        regions: Option<&Regions>,
    ) -> Result<Self> {
        let mut set = Self::default();
        let mut budget = MAX_DBSNP_PAGES;
        load_db_snp(path, dictionary, regions, false, |variant| {
            set.add(variant, true, &mut budget)
        })
        .map_err(|e| e.into_error(path))?;
        Ok(set)
    }

    /// The record's contig gets a set (empty when not `wanted`), and a
    /// wanted record sets every position from its start to its end.
    fn add(
        &mut self,
        variant: &Variant,
        wanted: bool,
        budget: &mut usize,
    ) -> std::result::Result<(), LoadError> {
        let bits = match self.by_contig.get_mut(&variant.contig) {
            Some(bits) => bits,
            None => self.by_contig.entry(variant.contig.clone()).or_default(),
        };
        if wanted {
            bits.set_range(variant.start, variant.end, budget)
                .map_err(|OverBudget| LoadError::Budget)?;
        }
        Ok(())
    }
}

/// Why a dbSNP load stopped.
enum LoadError {
    Read(Error),
    Budget,
}

impl LoadError {
    fn into_error(self, path: &Path) -> Error {
        match self {
            Self::Read(error) => error,
            Self::Budget => over_budget(path),
        }
    }
}

impl From<Error> for LoadError {
    fn from(error: Error) -> Self {
        Self::Read(error)
    }
}

/// `DbSnpBitSetUtil.loadVcf`: check the dictionaries, then hand `each`
/// every record of `path` (only those overlapping `regions`, when given).
/// Without `needs_dictionary`, a VCF without `##contig` lines skips the
/// check: `CollectSequencingArtifactMetrics` completes on one, while
/// `CollectVariantCallingMetrics` fails.
fn load_db_snp(
    path: &Path,
    dictionary: Option<&Dictionary>,
    regions: Option<&Regions>,
    needs_dictionary: bool,
    mut each: impl FnMut(&Variant) -> std::result::Result<(), LoadError>,
) -> std::result::Result<(), LoadError> {
    let mut reader = VcfReader::open_capped(path, MAX_LINE_BYTES, MAX_DBSNP_LINES)?;
    let own = reader.dictionary();
    if needs_dictionary || own.is_some() {
        assert_dictionaries_equal(
            own.as_ref(),
            &format!("DBSNP: {}", path.display()),
            dictionary,
            "INPUT",
        )?;
    }
    while let Some(variant) = reader.next_variant()? {
        if let Some(regions) = regions
            && !regions.overlaps(&variant.contig, variant.start, variant.end)
        {
            continue;
        }
        if let Some(dictionary) = dictionary
            && !dictionary.iter().any(|c| c.name == variant.contig)
        {
            return Err(LoadError::Read(Error::Intervals {
                path: path.to_path_buf(),
                message: format!(
                    "contig {} is not in the sequence dictionary (Picard fails with a NullPointerException)",
                    variant.contig
                ),
            }));
        }
        each(&variant)?;
    }
    Ok(())
}

/// The interval list support of the VCF tools, over `interval_ops`.
pub use crate::interval_ops::Regions;

/// `Tribble.indexPath` or `Tribble.tabixIndexPath` exists next to `path`.
pub fn index_exists(path: &Path) -> bool {
    let name = path.to_string_lossy();
    Path::new(&format!("{name}.idx")).exists() || Path::new(&format!("{name}.tbi")).exists()
}

/// Reads a whole Picard metrics file's `## METRICS CLASS` section as rows
/// of text (`MetricsFile.read` without the typed fields), for the
/// `Accumulate*` tools.
pub fn read_metrics_rows(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let text = fs::read_to_string(path).map_err(|source| Error::Input {
        path: path.to_path_buf(),
        source,
    })?;
    let mut lines = text.lines();
    let mut columns = Vec::new();
    let mut rows = Vec::new();
    while let Some(line) = lines.next() {
        if !line.starts_with("## METRICS CLASS") {
            continue;
        }
        columns = lines
            .next()
            .unwrap_or_default()
            .split('\t')
            .map(str::to_owned)
            .collect();
        for line in lines.by_ref() {
            if line.trim().is_empty() {
                break;
            }
            let row: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if row.len() != columns.len() {
                return Err(Error::Intervals {
                    path: path.to_path_buf(),
                    message: format!(
                        "metrics row with {} fields for {} columns",
                        row.len(),
                        columns.len()
                    ),
                });
            }
            rows.push(row);
        }
        break;
    }
    Ok((columns, rows))
}

/// `FormatUtil.parseDouble`: `?` is NaN.
pub fn parse_metric_double(text: &str) -> Result<f64> {
    if text == "?" || text == "-?" {
        return Ok(f64::NAN);
    }
    text.parse()
        .map_err(|_| Error::Metric(format!("{text:?} is not a number")))
}

pub fn parse_metric_long(text: &str) -> Result<i64> {
    if text.is_empty() {
        return Ok(0);
    }
    text.parse()
        .map_err(|_| Error::Metric(format!("{text:?} is not an integer")))
}

/// `io::Error` to the crate error for a path.
pub fn input_error(path: &Path, source: io::Error) -> Error {
    Error::Input {
        path: path.to_path_buf(),
        source,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn allele(bases: &str) -> Allele {
        Allele::parse(bases, false).expect("allele")
    }

    fn variant(reference: &str, alts: &[&str], gt: &[Option<usize>]) -> Variant {
        let mut alleles = vec![Allele::parse(reference, true).expect("ref")];
        alleles.extend(alts.iter().map(|a| allele(a)));
        Variant {
            contig: "c1".to_owned(),
            start: 10,
            end: 10 + reference.len() as i64 - 1,
            alleles,
            filtered: false,
            genotypes: vec![Genotype {
                alleles: gt.to_vec(),
                gq: -1,
                dp: -1,
                ad: None,
                filtered: false,
            }],
        }
    }

    #[test]
    fn variant_types_follow_determine_type() {
        assert_eq!(variant("A", &[], &[]).kind(), VariantType::NoVariation);
        assert_eq!(variant("A", &["G"], &[]).kind(), VariantType::Snp);
        assert_eq!(variant("A", &["G", "T"], &[]).kind(), VariantType::Snp);
        assert_eq!(variant("AT", &["GC"], &[]).kind(), VariantType::Mnp);
        assert_eq!(variant("A", &["AT"], &[]).kind(), VariantType::Indel);
        assert_eq!(variant("AT", &["A", "ATT"], &[]).kind(), VariantType::Indel);
        assert_eq!(variant("A", &["G", "AT"], &[]).kind(), VariantType::Mixed);
        assert_eq!(variant("A", &["<DEL>"], &[]).kind(), VariantType::Symbolic);
        let ins = variant("A", &["AT"], &[]);
        assert!(ins.is_simple_insertion() && !ins.is_simple_deletion());
        let del = variant("AT", &["A"], &[]);
        assert!(del.is_simple_deletion() && !del.is_complex_indel());
        assert!(variant("AT", &["GCC"], &[]).is_complex_indel());
        assert!(variant("A", &["AT", "ATT"], &[]).is_complex_indel());
        assert!(variant("A", &["G"], &[]).is_transition());
        assert!(!variant("A", &["C"], &[]).is_transition());
    }

    #[test]
    fn genotype_types_and_sub_context() {
        let v = variant("A", &["G", "T"], &[Some(0), Some(2)]);
        assert_eq!(v.genotypes[0].kind(), GenotypeType::Het);
        let sub = v.sub_context(0);
        assert_eq!(sub.alleles, vec![allele("A"), allele("T")]);
        assert_eq!(sub.genotypes[0].alleles, vec![Some(0), Some(1)]);
        let hom_ref = variant("A", &["G"], &[Some(0), Some(0)]).sub_context(0);
        assert_eq!(hom_ref.kind(), VariantType::NoVariation);
        let mixed = variant("A", &["G"], &[None, Some(1)]);
        assert_eq!(mixed.genotypes[0].kind(), GenotypeType::Mixed);
        assert_eq!(mixed.sub_context(0).kind(), VariantType::Snp);
        let no_call = variant("A", &["G"], &[None, None]);
        assert_eq!(no_call.genotypes[0].kind(), GenotypeType::NoCall);
        assert_eq!(no_call.sub_context(0).kind(), VariantType::NoVariation);
        assert_eq!(
            variant("A", &["G"], &[]).genotypes[0].kind(),
            GenotypeType::Unavailable
        );
        assert_eq!(
            variant("A", &["G"], &[Some(1), Some(1)]).genotypes[0].kind(),
            GenotypeType::HomVar
        );
    }

    #[test]
    fn genotype_fields_follow_the_codec() {
        let keys = ["GT", "GQ", "DP", "AD", "FT"];
        let g = parse_genotype(&keys, "0/1:20.5:7:3,4:PASS", 2).expect("parses");
        assert_eq!(g.alleles, vec![Some(0), Some(1)]);
        assert_eq!((g.gq, g.dp), (21, 7));
        assert_eq!(g.ad, Some(vec![3, 4]));
        assert!(!g.filtered);
        let g = parse_genotype(&keys, "1|.:.:.:3,.:lowq", 2).expect("parses");
        assert_eq!(g.alleles, vec![Some(1), None]);
        assert_eq!((g.gq, g.dp, g.ad.clone(), g.filtered), (-1, -1, None, true));
        let g = parse_genotype(&keys, "0/0", 2).expect("short");
        assert_eq!((g.gq, g.dp), (-1, -1));
        let g = parse_genotype(&keys, ".", 2).expect("dot");
        assert_eq!(g.alleles, vec![None]);
        assert!(parse_genotype(&keys, "0/2", 2).is_err());
        assert!(parse_genotype(&["GQ", "GT"], "20:0/1", 2).is_err());
        assert!(parse_genotype(&["GT"], "0/1:20", 2).is_err());
    }

    #[test]
    fn dictionaries_compare_by_name_and_length() {
        let a = vec![Contig {
            name: "c1".to_owned(),
            length: Some(10),
        }];
        let b = vec![Contig {
            name: "c1".to_owned(),
            length: Some(11),
        }];
        assert!(assert_dictionaries_equal(Some(&a), "x", Some(&a), "y").is_ok());
        assert!(assert_dictionaries_equal(Some(&a), "x", Some(&b), "y").is_err());
        assert!(assert_dictionaries_equal(None, "x", None, "y").is_ok());
        assert!(assert_dictionaries_equal(Some(&a), "x", None, "y").is_err());
    }

    #[test]
    fn bit_sets_grow() {
        let mut bits = BitSet::default();
        bits.set_range(130, 130, &mut 1).unwrap();
        assert!(bits.get(130) && !bits.get(129) && !bits.get(1000));
    }
}

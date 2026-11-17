//! The line parsers of [`super::GeneModel::read`]: GTF, GFF3 (its parent
//! chain resolved, parents read after their children included), BED12 and
//! BED.

use std::collections::HashMap;
use std::io::{BufRead, Read};
use std::path::Path;

use super::{
    Exon, Feature, Format, Gene, GeneModel, MAX_BLOCKS, MAX_FEATURES, MAX_LINE_BYTES, Strand, Tag,
    Transcript,
};
use crate::error::{Error, Result};

/// One feature line, its names resolved.
#[derive(Clone)]
struct Record {
    kind: Feature,
    contig: String,
    start: i64,
    end: i64,
    strand: Strand,
    gene_id: String,
    gene_name: Option<String>,
    transcript_id: String,
    number: Option<i64>,
    tag: Tag,
}

/// A GFF3 feature with an `ID` that is not an exon or a CDS: a transcript
/// or a gene, as far as its children need it.
#[derive(Default)]
struct GffNode {
    parent: Option<String>,
    transcript_id: Option<String>,
    gene_id: Option<String>,
    gene_name: Option<String>,
    name: Option<String>,
    tag: Tag,
}

/// A GFF3 exon or CDS whose parent was not read yet.
struct Pending {
    line: usize,
    parent: String,
    record: Record,
}

pub(super) struct Parser<'a> {
    path: &'a Path,
    feature: Feature,
    format: Option<Format>,
    gff3_directive: bool,
    contigs: Vec<String>,
    contig_index: HashMap<String, usize>,
    genes: Vec<Gene>,
    gene_index: HashMap<String, usize>,
    transcripts: Vec<Transcript>,
    transcript_index: HashMap<String, usize>,
    exons: Vec<Exon>,
    nodes: HashMap<String, GffNode>,
    pending: Vec<Pending>,
    /// BED regions per name so far, for their numbers.
    regions: HashMap<usize, i64>,
    /// [`MAX_FEATURES`], lower in a test.
    limit: usize,
}

impl<'a> Parser<'a> {
    pub(super) fn new(path: &'a Path, feature: Feature) -> Self {
        Self {
            path,
            feature,
            format: None,
            gff3_directive: false,
            contigs: Vec::new(),
            contig_index: HashMap::new(),
            genes: Vec::new(),
            gene_index: HashMap::new(),
            transcripts: Vec::new(),
            transcript_index: HashMap::new(),
            exons: Vec::new(),
            nodes: HashMap::new(),
            pending: Vec::new(),
            regions: HashMap::new(),
            limit: MAX_FEATURES,
        }
    }

    #[cfg(test)]
    pub(super) fn with_limit(mut self, limit: usize) -> Self {
        self.limit = limit;
        self
    }

    fn error(&self, line: usize, message: impl std::fmt::Display) -> Error {
        Error::Intervals {
            path: self.path.to_path_buf(),
            message: format!("line {line}: {message}"),
        }
    }

    pub(super) fn read(&mut self, mut reader: Box<dyn BufRead>) -> Result<()> {
        let mut buffer = Vec::new();
        let mut number = 0usize;
        loop {
            buffer.clear();
            let read = (&mut reader)
                .take(MAX_LINE_BYTES as u64 + 1)
                .read_until(b'\n', &mut buffer)
                .map_err(|source| Error::Input {
                    path: self.path.to_path_buf(),
                    source,
                })?;
            if read == 0 {
                return Ok(());
            }
            number += 1;
            if buffer.len() > MAX_LINE_BYTES {
                return Err(self.error(number, format!("longer than {MAX_LINE_BYTES} bytes")));
            }
            let mut line = std::str::from_utf8(&buffer)
                .map_err(|_| self.error(number, "not UTF-8 text"))?
                .trim_end_matches(['\n', '\r']);
            if number == 1 {
                line = line.strip_prefix('\u{feff}').unwrap_or(line);
            }
            if line.starts_with("##FASTA") {
                return Ok(());
            }
            self.line(number, line)?;
        }
    }

    fn line(&mut self, number: usize, line: &str) -> Result<()> {
        if line.starts_with("##gff-version 3") {
            self.gff3_directive = true;
        }
        if line.trim().is_empty() || line.starts_with('#') || is_bed_header(line) {
            return Ok(());
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let format = match self.format {
            Some(format) => format,
            None => {
                let format = self.sniff(&fields).ok_or_else(|| {
                    self.error(
                        number,
                        "neither GTF, GFF3 nor BED: a GTF or GFF3 line has 9 tab-separated \
                         fields, a BED line 3 to 12 with integer start and end",
                    )
                })?;
                self.format = Some(format);
                format
            }
        };
        let outcome = match format {
            Format::Gtf | Format::Gff3 => self.feature_line(format, number, &fields),
            Format::Bed12 => self.bed12_line(&fields),
            Format::Bed => self.bed_line(&fields),
        };
        outcome.map_err(|message| self.error(number, message))?;
        // Genes and transcripts count too: a GTF line of the other feature
        // kind adds both and no exon.
        let limit = self.limit;
        if self.exons.len() + self.pending.len() > limit
            || self.nodes.len() > limit
            || self.genes.len() > limit
            || self.transcripts.len() > limit
        {
            return Err(self.error(
                number,
                format!("more than {limit} features, genes or transcripts"),
            ));
        }
        Ok(())
    }

    fn sniff(&self, fields: &[&str]) -> Option<Format> {
        let int = |i: usize| fields.get(i).is_some_and(|f| f.parse::<i64>().is_ok());
        // A GTF or GFF3 feature column is text; a BED9 line with integer
        // name and score has integers in columns 3 and 4 too.
        if fields.len() == 9 && !int(2) && int(3) && int(4) {
            let attributes = fields[8];
            let gff3 =
                self.gff3_directive || (!attributes.contains('"') && attributes.contains('='));
            return Some(if gff3 { Format::Gff3 } else { Format::Gtf });
        }
        if (3..=12).contains(&fields.len()) && int(1) && int(2) {
            return Some(if fields.len() == 12 {
                Format::Bed12
            } else {
                Format::Bed
            });
        }
        None
    }

    /// A GTF or GFF3 line: `exon` and `CDS` features are kept, and in GFF3
    /// any other feature with an `ID` is remembered for its children.
    fn feature_line(
        &mut self,
        format: Format,
        number: usize,
        fields: &[&str],
    ) -> std::result::Result<(), String> {
        let [contig, _, kind, start, end, _, strand, _, attributes] = fields[..] else {
            return Err(format!(
                "{} fields; a {} line has 9",
                fields.len(),
                format.name()
            ));
        };
        let kind = match kind {
            "exon" => Some(Feature::Exon),
            "CDS" => Some(Feature::Cds),
            _ => None,
        };
        if format == Format::Gtf {
            let Some(kind) = kind else { return Ok(()) };
            let (start, end, strand) = coordinates(start, end, strand)?;
            let record = gtf_record(kind, contig, start, end, strand, attributes)?;
            return self.add(record);
        }
        let attributes = gff3_attributes(attributes)?;
        let first = |key: &str| {
            attributes
                .iter()
                .find(|(k, _)| *k == key)
                .and_then(|(_, values)| values.first().cloned())
        };
        let tag = attributes
            .iter()
            .filter(|(k, _)| *k == "tag")
            .flat_map(|(_, values)| values.iter())
            .map(|t| Tag::of(t))
            .max()
            .unwrap_or_default();
        let Some(kind) = kind else {
            if let Some(id) = first("ID") {
                self.nodes.insert(
                    id,
                    GffNode {
                        parent: first("Parent"),
                        transcript_id: first("transcript_id"),
                        gene_id: first("gene_id"),
                        gene_name: first("gene_name"),
                        name: first("Name"),
                        tag,
                    },
                );
            }
            return Ok(());
        };
        let (start, end, strand) = coordinates(start, end, strand)?;
        let parents: Vec<String> = attributes
            .iter()
            .find(|(k, _)| *k == "Parent")
            .map(|(_, v)| v.clone())
            .unwrap_or_default();
        if parents.is_empty() {
            return Err(format!("{} feature without a Parent", kind.name()));
        }
        let exon_number = match first("exon_number").or_else(|| first("rank")) {
            None => None,
            Some(n) => Some(n.parse::<i64>().map_err(|_| format!("exon number {n:?}"))?),
        };
        for parent in parents {
            let record = Record {
                kind,
                contig: contig.to_owned(),
                start,
                end,
                strand,
                gene_id: first("gene_id").unwrap_or_default(),
                gene_name: first("gene_name"),
                transcript_id: first("transcript_id").unwrap_or_default(),
                number: exon_number,
                tag,
            };
            if self.resolvable(&parent) {
                let record = self.resolve(&parent, record);
                self.add(record)?;
            } else {
                self.pending.push(Pending {
                    line: number,
                    parent,
                    record,
                });
            }
        }
        Ok(())
    }

    /// The parent and its own parent, when it has one, are both read.
    fn resolvable(&self, parent: &str) -> bool {
        self.nodes.get(parent).is_some_and(|node| {
            node.parent
                .as_ref()
                .is_none_or(|g| self.nodes.contains_key(g))
        })
    }

    /// A GFF3 record's transcript and gene from its parent chain: the
    /// feature's own attributes first, then the transcript's, then the
    /// gene's, then the IDs without their `transcript:` or `gene:` prefix.
    fn resolve(&self, parent: &str, mut record: Record) -> Record {
        let node = self.nodes.get(parent);
        let gene_node = node
            .and_then(|n| n.parent.as_deref())
            .and_then(|g| self.nodes.get(g));
        if record.transcript_id.is_empty() {
            record.transcript_id = node
                .and_then(|n| n.transcript_id.clone())
                .unwrap_or_else(|| strip(parent, "transcript:"));
        }
        if record.gene_id.is_empty() {
            record.gene_id = node
                .and_then(|n| n.gene_id.clone())
                .or_else(|| gene_node.and_then(|g| g.gene_id.clone()))
                .or_else(|| {
                    node.and_then(|n| n.parent.as_deref())
                        .map(|g| strip(g, "gene:"))
                })
                .unwrap_or_else(|| record.transcript_id.clone());
        }
        if record.gene_name.is_none() {
            record.gene_name = node
                .and_then(|n| n.gene_name.clone())
                .or_else(|| gene_node.and_then(|g| g.gene_name.clone().or_else(|| g.name.clone())));
        }
        record.tag = record.tag.max(node.map_or(Tag::None, |n| n.tag));
        record
    }

    fn bed12_line(&mut self, fields: &[&str]) -> std::result::Result<(), String> {
        if fields.len() != 12 {
            return Err(format!("{} fields; a BED12 line has 12", fields.len()));
        }
        let (start, end) = bed_span(fields[1], fields[2])?;
        let name = fields[3];
        if name.is_empty() {
            return Err("BED12 line without a name".to_owned());
        }
        let strand = Strand::parse(fields[5]).ok_or_else(|| format!("strand {:?}", fields[5]))?;
        let thick_start = bounded(fields[6])?;
        let thick_end = bounded(fields[7])?;
        let count = int(fields[9])?;
        let sizes = int_list(fields[10])?;
        let offsets = int_list(fields[11])?;
        if count < 1 || sizes.len() as i64 != count || offsets.len() as i64 != count {
            return Err(format!(
                "blockCount {count} with {} sizes and {} starts",
                sizes.len(),
                offsets.len()
            ));
        }
        // 1-based closed blocks, sorted by start.
        let mut blocks: Vec<(i64, i64)> = Vec::with_capacity(sizes.len());
        for (&size, &offset) in sizes.iter().zip(&offsets) {
            // Both under 2^31 before they are added, so the sums stay far
            // inside i64.
            let in_range =
                (1..=MAX_COORDINATE).contains(&size) && (0..=MAX_COORDINATE).contains(&offset);
            let first = start + if in_range { offset } else { 0 };
            let last = first + if in_range { size } else { 0 } - 1;
            if !in_range || last > end {
                return Err(format!(
                    "block of {size} bases at offset {offset} outside {}..{end}",
                    start
                ));
            }
            blocks.push((first, last));
        }
        blocks.sort_unstable();
        let thick = (thick_start + 1, thick_end);
        let count = blocks.len() as i64;
        for (rank, &(first, last)) in blocks.iter().enumerate() {
            let number = if strand == Strand::Reverse {
                count - rank as i64
            } else {
                rank as i64 + 1
            };
            let mut record = Record {
                kind: Feature::Exon,
                contig: fields[0].to_owned(),
                start: first,
                end: last,
                strand,
                gene_id: name.to_owned(),
                gene_name: None,
                transcript_id: name.to_owned(),
                number: Some(number),
                tag: Tag::None,
            };
            let (cds_first, cds_last) = (first.max(thick.0), last.min(thick.1));
            let cds = (cds_first <= cds_last).then_some((cds_first, cds_last));
            self.add(record.clone())?;
            if let Some((cds_first, cds_last)) = cds {
                record.kind = Feature::Cds;
                record.start = cds_first;
                record.end = cds_last;
                self.add(record)?;
            }
        }
        Ok(())
    }

    fn bed_line(&mut self, fields: &[&str]) -> std::result::Result<(), String> {
        if !(3..=11).contains(&fields.len()) {
            return Err(format!(
                "{} fields; a BED line has 3 to 11, or 12 for BED12 in every line",
                fields.len()
            ));
        }
        if self.feature == Feature::Cds {
            return Err("FEATURE=CDS needs a GTF, GFF3 or BED12 annotation".to_owned());
        }
        let (start, end) = bed_span(fields[1], fields[2])?;
        let name = match fields
            .get(3)
            .copied()
            .filter(|n| !n.is_empty() && *n != ".")
        {
            Some(name) => name.to_owned(),
            None => format!("{}:{}-{}", fields[0], start, end),
        };
        let strand = match fields.get(5) {
            Some(s) => Strand::parse(s).ok_or_else(|| format!("strand {s:?}"))?,
            None => Strand::Unknown,
        };
        let gene = self.gene(&name, None);
        let transcript = self.transcript(&name, gene);
        let count = self.regions.entry(transcript).or_insert(0);
        *count += 1;
        let number = *count;
        self.add(Record {
            kind: Feature::Exon,
            contig: fields[0].to_owned(),
            start,
            end,
            strand,
            gene_id: name.clone(),
            gene_name: None,
            transcript_id: name,
            number: Some(number),
            tag: Tag::None,
        })
    }

    fn gene(&mut self, id: &str, name: Option<&str>) -> usize {
        if let Some(&i) = self.gene_index.get(id) {
            if let Some(name) = name
                && self.genes[i].name == self.genes[i].id
            {
                self.genes[i].name = name.to_owned();
            }
            return i;
        }
        self.genes.push(Gene {
            id: id.to_owned(),
            name: name.unwrap_or(id).to_owned(),
        });
        self.gene_index.insert(id.to_owned(), self.genes.len() - 1);
        self.genes.len() - 1
    }

    /// The transcript `id`; a transcript keeps the gene of its first line.
    fn transcript(&mut self, id: &str, gene: usize) -> usize {
        if let Some(&i) = self.transcript_index.get(id) {
            return i;
        }
        self.transcripts.push(Transcript {
            id: id.to_owned(),
            gene,
            tag: Tag::None,
            cds_bases: 0,
            exon_bases: 0,
        });
        self.transcript_index
            .insert(id.to_owned(), self.transcripts.len() - 1);
        self.transcripts.len() - 1
    }

    fn contig(&mut self, name: &str) -> usize {
        if let Some(&i) = self.contig_index.get(name) {
            return i;
        }
        self.contigs.push(name.to_owned());
        self.contig_index
            .insert(name.to_owned(), self.contigs.len() - 1);
        self.contigs.len() - 1
    }

    fn add(&mut self, record: Record) -> std::result::Result<(), String> {
        if record.gene_id.is_empty() {
            return Err("feature without a gene_id".to_owned());
        }
        if record.transcript_id.is_empty() {
            return Err("feature without a transcript_id".to_owned());
        }
        let gene = self.gene(&record.gene_id, record.gene_name.as_deref());
        let transcript = self.transcript(&record.transcript_id, gene);
        let t = &mut self.transcripts[transcript];
        t.tag = t.tag.max(record.tag);
        let length = record.end - record.start + 1;
        match record.kind {
            Feature::Exon => t.exon_bases += length,
            Feature::Cds => t.cds_bases += length,
        }
        if record.kind == self.feature {
            let contig = self.contig(&record.contig);
            self.exons.push(Exon {
                transcript,
                contig,
                start: record.start,
                end: record.end,
                strand: record.strand,
                number: record.number,
            });
        }
        Ok(())
    }

    pub(super) fn finish(mut self) -> Result<GeneModel> {
        // Parents read after their children: resolve in file order now.
        let pending = std::mem::take(&mut self.pending);
        for Pending {
            line,
            parent,
            record,
        } in pending
        {
            if !self.nodes.contains_key(&parent) {
                return Err(self.error(
                    line,
                    format!("Parent {parent} is not the ID of any feature"),
                ));
            }
            let record = self.resolve(&parent, record);
            self.add(record)
                .map_err(|message| self.error(line, message))?;
        }
        let format = self.format.ok_or_else(|| Error::Intervals {
            path: self.path.to_path_buf(),
            message: "no feature line".to_owned(),
        })?;
        Ok(GeneModel {
            format,
            feature: self.feature,
            contigs: self.contigs,
            genes: self.genes,
            transcripts: self.transcripts,
            exons: self.exons,
        })
    }
}

/// A BED `track` or `browser` line: the word, then a space, a tab or
/// nothing.
fn is_bed_header(line: &str) -> bool {
    ["track", "browser"].iter().any(|word| {
        line.strip_prefix(word)
            .is_some_and(|rest| rest.is_empty() || rest.starts_with([' ', '\t']))
    })
}

fn strip(id: &str, prefix: &str) -> String {
    id.strip_prefix(prefix).unwrap_or(id).to_owned()
}

/// Start and end of a GTF or GFF3 line, 1-based closed, and its strand.
fn coordinates(
    start: &str,
    end: &str,
    strand: &str,
) -> std::result::Result<(i64, i64, Strand), String> {
    let (start, end) = (int(start)?, int(end)?);
    if start < 1 || end < start || end > MAX_COORDINATE {
        return Err(format!("start {start} and end {end}: 1 <= start <= end"));
    }
    let strand = Strand::parse(strand).ok_or_else(|| format!("strand {strand:?}"))?;
    Ok((start, end, strand))
}

/// A BED start and end as 1-based closed coordinates. A region of no base
/// is an error.
fn bed_span(start: &str, end: &str) -> std::result::Result<(i64, i64), String> {
    let (start, end) = (int(start)?, int(end)?);
    if start < 0 || end <= start || end > MAX_COORDINATE {
        return Err(format!(
            "start {start} and end {end}: 0 <= start < end (a BED region holds a base)"
        ));
    }
    Ok((start + 1, end))
}

/// The largest coordinate, size or offset read: htsjdk keeps positions in
/// an `int`.
const MAX_COORDINATE: i64 = i32::MAX as i64;

/// An integer from 0 to [`MAX_COORDINATE`].
fn bounded(field: &str) -> std::result::Result<i64, String> {
    let value = int(field)?;
    if (0..=MAX_COORDINATE).contains(&value) {
        Ok(value)
    } else {
        Err(format!("{value} is outside 0..{MAX_COORDINATE}"))
    }
}

fn int(field: &str) -> std::result::Result<i64, String> {
    field
        .trim()
        .parse::<i64>()
        .map_err(|_| format!("{field:?} is not an integer"))
}

/// A BED12 list such as `10,20,`: the trailing comma is allowed.
fn int_list(field: &str) -> std::result::Result<Vec<i64>, String> {
    let items: Vec<&str> = field
        .strip_suffix(',')
        .unwrap_or(field)
        .split(',')
        .collect();
    if items.len() > MAX_BLOCKS {
        return Err(format!("more than {MAX_BLOCKS} blocks"));
    }
    items.into_iter().map(int).collect()
}

/// A GTF line's names: `gene_id` and `transcript_id` are required (GTF
/// 2.2); `gene_name`, `exon_number` and any `tag` are read when present.
fn gtf_record(
    kind: Feature,
    contig: &str,
    start: i64,
    end: i64,
    strand: Strand,
    attributes: &str,
) -> std::result::Result<Record, String> {
    let mut record = Record {
        kind,
        contig: contig.to_owned(),
        start,
        end,
        strand,
        gene_id: String::new(),
        gene_name: None,
        transcript_id: String::new(),
        number: None,
        tag: Tag::None,
    };
    for (key, value) in gtf_attributes(attributes)? {
        match key {
            "gene_id" => record.gene_id = value.to_owned(),
            "transcript_id" => record.transcript_id = value.to_owned(),
            "gene_name" => record.gene_name = Some(value.to_owned()),
            "exon_number" => {
                record.number = Some(
                    value
                        .parse::<i64>()
                        .map_err(|_| format!("exon_number {value:?}"))?,
                );
            }
            "tag" => record.tag = record.tag.max(Tag::of(value)),
            _ => {}
        }
    }
    Ok(record)
}

/// GTF attributes: `key "value";` or `key value;` pairs. A `;` inside
/// quotes belongs to the value.
fn gtf_attributes(text: &str) -> std::result::Result<Vec<(&str, &str)>, String> {
    let mut pairs = Vec::new();
    let mut rest = text;
    loop {
        rest = rest.trim_start_matches(|c: char| c == ';' || c.is_whitespace());
        if rest.is_empty() {
            return Ok(pairs);
        }
        let key_end = rest
            .find(|c: char| c == ';' || c.is_whitespace())
            .unwrap_or(rest.len());
        let key = &rest[..key_end];
        rest = rest[key_end..].trim_start_matches([' ', '\t']);
        let value;
        if let Some(quoted) = rest.strip_prefix('"') {
            let close = quoted
                .find('"')
                .ok_or_else(|| format!("attribute {key}: no closing quote"))?;
            value = &quoted[..close];
            rest = &quoted[close + 1..];
        } else {
            let end = rest.find(';').unwrap_or(rest.len());
            value = rest[..end].trim_end();
            rest = &rest[end..];
        }
        pairs.push((key, value));
    }
}

/// GFF3 attributes: `key=value` pairs separated by `;`, each value a
/// comma-separated list with `%XX` escapes decoded.
fn gff3_attributes(text: &str) -> std::result::Result<Vec<(&str, Vec<String>)>, String> {
    let mut pairs = Vec::new();
    for item in text
        .split(';')
        .map(str::trim)
        .filter(|i| !i.is_empty() && *i != ".")
    {
        let (key, value) = item
            .split_once('=')
            .ok_or_else(|| format!("attribute {item:?} without ="))?;
        let values = value
            .split(',')
            .map(percent_decode)
            .collect::<std::result::Result<Vec<String>, String>>()?;
        pairs.push((key.trim(), values));
    }
    Ok(pairs)
}

fn percent_decode(text: &str) -> std::result::Result<String, String> {
    if !text.contains('%') {
        return Ok(text.to_owned());
    }
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' {
            let hex = text
                .get(i + 1..i + 3)
                .and_then(|h| u8::from_str_radix(h, 16).ok())
                .ok_or_else(|| format!("bad escape in {text:?}"))?;
            out.push(hex);
            i += 3;
        } else {
            out.push(bytes[i]);
            i += 1;
        }
    }
    String::from_utf8(out).map_err(|_| format!("{text:?} decodes to bytes that are not UTF-8"))
}

//! Gene annotations for per-exon coverage (docs/07 section 4, wave W07):
//! the `exon` or `CDS` features of a GTF 2.2 or GFF3 file, the blocks of a
//! BED12 file, or the lines of a BED file with 3 to 11 columns as plain
//! regions. A file is plain text or gzip (one or more members, BGZF too),
//! told by its first two bytes. The format is told by the first line that
//! is not a comment, a `track` or `browser` line or blank.
//!
//! Coordinates are 1-based and closed, as in GTF and htsjdk `Interval`; a
//! BED `start` is read as `start + 1`.
//!
//! The file is streamed, never read whole: a line holds at most
//! [`MAX_LINE_BYTES`] and at most [`MAX_FEATURES`] features are kept. Gene,
//! transcript and contig names are stored once each.
//!
//! [`GeneModel::select`] applies the filters of a run (genes, transcripts,
//! one canonical transcript per gene) and merges the overlapping features of
//! one transcript, so a transcript's bases count once. Features of
//! different transcripts are never merged: each transcript reports its own
//! exons.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::{Error, Result};
use crate::gzip::GzipReader;

mod parse;

use parse::Parser;

/// The bytes of one line, line end included.
pub const MAX_LINE_BYTES: usize = 1 << 20;

/// The features of the chosen kind one annotation may hold (GENCODE 46
/// holds 1.65 million exon lines), and the GFF3 genes and transcripts.
pub const MAX_FEATURES: usize = 5_000_000;

/// The bytes of a gene list file, which is read whole.
pub const MAX_LIST_BYTES: u64 = 16 << 20;

/// The blocks of one BED12 line.
pub const MAX_BLOCKS: usize = 100_000;

/// The kind of annotation file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Gtf,
    Gff3,
    Bed12,
    /// BED with 3 to 11 columns: each line one region.
    Bed,
}

impl Format {
    pub fn name(self) -> &'static str {
        match self {
            Format::Gtf => "GTF",
            Format::Gff3 => "GFF3",
            Format::Bed12 => "BED12",
            Format::Bed => "BED",
        }
    }

    /// GTF, GFF3 and BED12 describe transcripts; BED lines are regions.
    pub fn has_transcripts(self) -> bool {
        self != Format::Bed
    }
}

/// The feature a run reports.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Feature {
    #[default]
    Exon,
    Cds,
}

impl Feature {
    /// `exon` or `CDS`, in any case.
    pub fn parse(text: &str) -> Result<Self> {
        if text.eq_ignore_ascii_case("exon") {
            Ok(Feature::Exon)
        } else if text.eq_ignore_ascii_case("cds") {
            Ok(Feature::Cds)
        } else {
            Err(Error::Usage(format!("FEATURE={text}: exon or CDS")))
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Feature::Exon => "exon",
            Feature::Cds => "CDS",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    fn parse(text: &str) -> Option<Self> {
        match text {
            "+" => Some(Strand::Forward),
            "-" => Some(Strand::Reverse),
            "." | "?" => Some(Strand::Unknown),
            _ => None,
        }
    }

    pub fn symbol(self) -> &'static str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        }
    }
}

/// The canonical tag of a transcript, ordered by preference.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum Tag {
    #[default]
    None,
    ManeSelect,
    EnsemblCanonical,
}

impl Tag {
    fn of(tag: &str) -> Self {
        match tag {
            "Ensembl_canonical" => Tag::EnsemblCanonical,
            // GENCODE writes MANE_Select, NCBI RefSeq GTF "MANE Select".
            "MANE_Select" | "MANE Select" => Tag::ManeSelect,
            _ => Tag::None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    /// Index into [`GeneModel::genes`].
    pub gene: usize,
    pub tag: Tag,
    /// Bases of the transcript's CDS features (BED12: its blocks inside
    /// the thick part).
    pub cds_bases: i64,
    /// Bases of its exon features.
    pub exon_bases: i64,
}

/// One feature of the chosen kind.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Exon {
    /// Index into [`GeneModel::transcripts`].
    pub transcript: usize,
    /// Index into [`GeneModel::contigs`].
    pub contig: usize,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    /// `exon_number` (GFF3: also `rank`) as written; [`GeneModel::select`]
    /// numbers the rest along the strand.
    pub number: Option<i64>,
}

impl Exon {
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// The filters of one run.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Selection {
    /// Gene names or IDs (an ID matches with or without its `.version`).
    pub genes: Vec<String>,
    /// Transcript IDs, with or without the `.version`.
    pub transcripts: Vec<String>,
    /// One transcript per gene: `Ensembl_canonical`, else `MANE_Select`,
    /// else the longest CDS.
    pub canonical: bool,
}

/// An annotation file as read: its features of one kind and the genes and
/// transcripts they belong to.
#[derive(Clone, Debug)]
pub struct GeneModel {
    pub format: Format,
    pub feature: Feature,
    pub contigs: Vec<String>,
    pub genes: Vec<Gene>,
    pub transcripts: Vec<Transcript>,
    pub exons: Vec<Exon>,
}

/// One selected feature: merged within its transcript and numbered.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Selected {
    pub transcript: usize,
    pub contig: usize,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub number: i64,
}

impl Selected {
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

impl GeneModel {
    /// Read the features of kind `feature` from `path`. CDS lengths are
    /// summed per transcript whatever the kind, for the canonical choice.
    pub fn read(path: &Path, feature: Feature) -> Result<Self> {
        let file = File::open(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        let mut plain = BufReader::with_capacity(1 << 16, file);
        let gzip = plain
            .fill_buf()
            .map_err(|source| Error::Input {
                path: path.to_path_buf(),
                source,
            })?
            .starts_with(&[0x1f, 0x8b]);
        let reader: Box<dyn BufRead> = if gzip {
            Box::new(BufReader::with_capacity(1 << 16, GzipReader::new(plain)))
        } else {
            Box::new(plain)
        };
        let mut parser = Parser::new(path, feature);
        parser.read(reader)?;
        parser.finish()
    }

    /// The features `selection` keeps, each transcript's overlapping
    /// features merged, numbered and in file order. An error names every
    /// requested gene or transcript the file lacks, and an empty result is
    /// an error.
    pub fn select(&self, selection: &Selection) -> Result<Vec<Selected>> {
        if selection.canonical && !self.format.has_transcripts() {
            return Err(Error::Usage(
                "CANONICAL needs a GTF, GFF3 or BED12 annotation; a BED file has no transcripts"
                    .to_owned(),
            ));
        }
        let genes = self.requested_genes(&selection.genes)?;
        let transcripts = self.requested_transcripts(&selection.transcripts)?;
        let canonical = selection.canonical.then(|| self.canonical_transcripts());
        let keep: Vec<bool> = self
            .transcripts
            .iter()
            .enumerate()
            .map(|(i, t)| {
                genes.as_ref().is_none_or(|g| g[t.gene])
                    && transcripts.as_ref().is_none_or(|k| k[i])
                    && canonical.as_ref().is_none_or(|c| c[t.gene] == Some(i))
            })
            .collect();
        let mut by_transcript: Vec<Vec<usize>> = vec![Vec::new(); self.transcripts.len()];
        for (i, exon) in self.exons.iter().enumerate() {
            if keep[exon.transcript] {
                by_transcript[exon.transcript].push(i);
            }
        }
        let mut selected: Vec<(usize, Selected)> = Vec::new();
        for members in by_transcript.iter().filter(|m| !m.is_empty()) {
            selected.extend(self.merge_transcript(members));
        }
        if selected.is_empty() {
            return Err(Error::Usage(format!(
                "no {} of the annotation passes the filters",
                self.feature.name()
            )));
        }
        selected.sort_by_key(|(first, _)| *first);
        Ok(selected.into_iter().map(|(_, s)| s).collect())
    }

    /// One transcript's features (indices into `exons`) merged where they
    /// overlap on one contig, each with the index of its first feature in
    /// the file. A merged feature keeps the lowest number written; features
    /// without one are numbered along the strand.
    fn merge_transcript(&self, members: &[usize]) -> Vec<(usize, Selected)> {
        let mut sorted = members.to_vec();
        sorted.sort_by_key(|&i| (self.exons[i].contig, self.exons[i].start, i));
        let mut merged: Vec<(usize, Selected, Option<i64>)> = Vec::new();
        for i in sorted {
            let exon = &self.exons[i];
            match merged.last_mut() {
                Some((first, current, number))
                    if current.contig == exon.contig && exon.start <= current.end =>
                {
                    current.end = current.end.max(exon.end);
                    *first = (*first).min(i);
                    *number = match (*number, exon.number) {
                        (Some(a), Some(b)) => Some(a.min(b)),
                        (a, b) => a.or(b),
                    };
                }
                _ => merged.push((
                    i,
                    Selected {
                        transcript: exon.transcript,
                        contig: exon.contig,
                        start: exon.start,
                        end: exon.end,
                        strand: exon.strand,
                        number: 0,
                    },
                    exon.number,
                )),
            }
        }
        let count = merged.len() as i64;
        let reverse = merged
            .first()
            .is_some_and(|(_, s, _)| s.strand == Strand::Reverse);
        merged
            .into_iter()
            .enumerate()
            .map(|(rank, (first, mut selected, number))| {
                let along = if reverse {
                    count - rank as i64
                } else {
                    rank as i64 + 1
                };
                selected.number = number.unwrap_or(along);
                (first, selected)
            })
            .collect()
    }

    /// One flag per gene: requested by name or ID. `None` without a request.
    fn requested_genes(&self, requests: &[String]) -> Result<Option<Vec<bool>>> {
        if requests.is_empty() {
            return Ok(None);
        }
        let mut index: HashMap<&str, Vec<usize>> = HashMap::new();
        for (i, gene) in self.genes.iter().enumerate() {
            index.entry(&gene.name).or_default().push(i);
            index.entry(&gene.id).or_default().push(i);
            if let Some(base) = without_version(&gene.id) {
                index.entry(base).or_default().push(i);
            }
        }
        let names = self.genes.iter().map(|g| g.name.as_str());
        mark(requests, &index, self.genes.len(), "gene", names).map(Some)
    }

    /// One flag per transcript: requested by ID. `None` without a request.
    fn requested_transcripts(&self, requests: &[String]) -> Result<Option<Vec<bool>>> {
        if requests.is_empty() {
            return Ok(None);
        }
        let mut index: HashMap<&str, Vec<usize>> = HashMap::new();
        for (i, transcript) in self.transcripts.iter().enumerate() {
            index.entry(&transcript.id).or_default().push(i);
            if let Some(base) = without_version(&transcript.id) {
                index.entry(base).or_default().push(i);
            }
        }
        let ids = self.transcripts.iter().map(|t| t.id.as_str());
        mark(requests, &index, self.transcripts.len(), "transcript", ids).map(Some)
    }

    /// The canonical transcript of each gene: the best tag
    /// (`Ensembl_canonical`, then `MANE_Select`), then the most CDS bases,
    /// then the most exon bases, then the lowest ID.
    pub fn canonical_transcripts(&self) -> Vec<Option<usize>> {
        let mut best: Vec<Option<usize>> = vec![None; self.genes.len()];
        for (i, t) in self.transcripts.iter().enumerate() {
            let better = match best[t.gene] {
                None => true,
                Some(b) => {
                    let b = &self.transcripts[b];
                    (t.tag, t.cds_bases, t.exon_bases, std::cmp::Reverse(&t.id))
                        > (b.tag, b.cds_bases, b.exon_bases, std::cmp::Reverse(&b.id))
                }
            };
            if better {
                best[t.gene] = Some(i);
            }
        }
        best
    }
}

/// `ENSG00000012048.23` without `.23`; `None` for an ID without a numeric
/// version.
fn without_version(id: &str) -> Option<&str> {
    let (base, version) = id.rsplit_once('.')?;
    (!base.is_empty() && !version.is_empty() && version.bytes().all(|b| b.is_ascii_digit()))
        .then_some(base)
}

/// Flags for the items `requests` name through `index`; an error lists the
/// requests no item answers, the first with the closest name.
fn mark<'a>(
    requests: &[String],
    index: &HashMap<&str, Vec<usize>>,
    len: usize,
    what: &str,
    names: impl Iterator<Item = &'a str>,
) -> Result<Vec<bool>> {
    let mut flags = vec![false; len];
    let mut missing = Vec::new();
    for request in requests {
        match index.get(request.as_str()) {
            Some(items) => items.iter().for_each(|&i| flags[i] = true),
            None => missing.push(request.as_str()),
        }
    }
    if missing.is_empty() {
        return Ok(flags);
    }
    let shown: Vec<&str> = missing.iter().copied().take(10).collect();
    let more = missing.len().saturating_sub(shown.len());
    let mut message = format!(
        "{} {what}{} not in the annotation: {}{}",
        missing.len(),
        if missing.len() == 1 { " is" } else { "s are" },
        shown.join(", "),
        if more > 0 {
            format!(" and {more} more")
        } else {
            String::new()
        }
    );
    if let Some(close) = crate::registry::did_you_mean(missing[0], names) {
        message.push_str(&format!(" (did you mean {close}?)"));
    }
    Err(Error::Usage(message))
}

/// Gene names from a list file: one per line, the first word of the line;
/// blank lines and `#` lines skipped. The file is read whole, at most
/// [`MAX_LIST_BYTES`].
pub fn read_gene_list(path: &Path) -> Result<Vec<String>> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let size = std::fs::metadata(path).map_err(input)?.len();
    if size > MAX_LIST_BYTES {
        return Err(Error::Intervals {
            path: path.to_path_buf(),
            message: format!("the gene list holds {size} bytes; at most {MAX_LIST_BYTES} are read"),
        });
    }
    let mut text = String::new();
    File::open(path)
        .map_err(input)?
        .take(MAX_LIST_BYTES)
        .read_to_string(&mut text)
        .map_err(input)?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| l.split_whitespace().next())
        .map(str::to_owned)
        .collect())
}

#[cfg(test)]
mod tests;

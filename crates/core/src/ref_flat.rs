//! M17: the refFlat gene model as Picard's `RefFlatReader` builds it and
//! `Gene.Transcript` reads it, with the gene overlap detector.
//!
//! Picard 3.5.0 behaviour (`picard.annotation.RefFlatReader`, `Gene`):
//! * A line is 11 tab-separated fields (runs of tabs count as one, blank
//!   lines and `#` lines are skipped): gene name, transcript name,
//!   chromosome, strand, transcription start and end, coding start and end,
//!   exon count, exon starts, exon ends. Starts are 0-based and are moved
//!   to 1-based here; ends stay. Any other field count is an error.
//! * A line on a chromosome the BAM header lacks is dropped. Lines are
//!   grouped by gene name; a gene whose lines disagree on strand or
//!   chromosome, repeat a transcript name, list an exon count that is not
//!   the number of starts or ends, hold an exon of zero or negative extent
//!   or two exons that overlap, is dropped whole (`AnnotationException`,
//!   logged at debug level). A field that is not an integer is an error.
//! * The gene spans the lowest transcription start to the highest end and
//!   takes the strand of its first line. Genes go into an `OverlapDetector`
//!   whose tree merges equal spans into one set, and `Gene.equals` is the
//!   span plus the strand, so a second gene with the span and strand of an
//!   earlier one is dropped. Genes are added in `HashMap` order of their
//!   names, and a gene iterates its transcripts in `HashMap` order of their
//!   names ([`crate::java`]).
//! * A locus is `CODING` inside an exon and the coding range, `UTR` inside
//!   an exon outside it, `INTRONIC` inside the transcription range
//!   otherwise, `INTERGENIC` elsewhere; the strongest function over every
//!   transcript of every overlapping gene wins (`LocusFunction` order). A
//!   non-coding transcript (coding start one past the coding end) is all
//!   `UTR`.
//! * `Transcript.addCoverageCounts(start, end, coverage)` runs `start`
//!   inclusive to `end` exclusive, and `RnaSeqMetricsCollector` passes the
//!   block's last base as `end`, so the last base of every alignment block
//!   adds no coverage. [`Transcript::add_coverage`] keeps that.

use std::fs;
use std::path::Path;

use noodles::sam;

use crate::error::{Error, Result};
use crate::java::{hash_map_order, hash_map_order_by_hash, string_hash_code};

/// `picard.annotation.LocusFunction`, in its strength order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum LocusFunction {
    Intergenic,
    Intronic,
    Utr,
    Coding,
    Ribosomal,
}

/// `Gene.Transcript.Exon`: 1-based, closed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Exon {
    pub start: i64,
    pub end: i64,
}

/// `Gene.Transcript`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Transcript {
    pub name: String,
    pub transcription_start: i64,
    pub transcription_end: i64,
    pub coding_start: i64,
    pub coding_end: i64,
    /// In file order, which the reader requires to be ascending.
    pub exons: Vec<Exon>,
    /// Summed exon lengths: the coverage array length.
    pub length: i64,
    /// `Transcript.hashCode`: the name's hash folded with the four
    /// coordinates, which `HashMap` keyed by transcripts orders by.
    pub hash: i32,
}

impl Transcript {
    /// `assignLocusFunctionForRange`: `functions[i]` describes the locus
    /// `start + i`; a transcript only raises a function, never lowers it.
    pub fn assign_locus_function(&self, start: i64, functions: &mut [LocusFunction]) {
        let last = start + functions.len() as i64 - 1;
        let from = start.max(self.transcription_start);
        let to = last.min(self.transcription_end);
        let mut locus = from;
        while locus <= to {
            let slot = &mut functions[(locus - start) as usize];
            if *slot <= LocusFunction::Coding {
                let function = if self.in_exon(locus) {
                    if self.utr(locus) {
                        LocusFunction::Utr
                    } else {
                        LocusFunction::Coding
                    }
                } else {
                    LocusFunction::Intronic
                };
                if function > *slot {
                    *slot = function;
                }
            }
            locus += 1;
        }
    }

    /// `addCoverageCounts(genomeStart, genomeEnd, coverage)`: one count
    /// per exonic locus from `genome_start` up to, not including,
    /// `genome_end`.
    pub fn add_coverage(&self, genome_start: i64, genome_end: i64, coverage: &mut [i32]) {
        let mut locus = genome_start;
        while locus < genome_end {
            let base = self.transcript_coordinate(locus);
            if base > 0 {
                coverage[(base - 1) as usize] += 1;
            }
            locus += 1;
        }
    }

    /// `getTranscriptCoordinate`: the 1-based position of a genomic locus
    /// in the spliced transcript, -1 outside its exons.
    pub fn transcript_coordinate(&self, locus: i64) -> i64 {
        let mut offset = 0;
        for exon in &self.exons {
            if locus >= exon.start && locus <= exon.end {
                return locus - exon.start + 1 + offset;
            }
            offset += exon.end - exon.start + 1;
        }
        -1
    }

    fn utr(&self, locus: i64) -> bool {
        locus < self.coding_start || locus > self.coding_end
    }

    /// `inExon`: stops at the first exon starting past the locus.
    fn in_exon(&self, locus: i64) -> bool {
        for exon in &self.exons {
            if exon.start > locus {
                return false;
            }
            if locus >= exon.start && locus <= exon.end {
                return true;
            }
        }
        false
    }
}

/// `picard.annotation.Gene`: an interval with its transcripts.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gene {
    pub name: String,
    /// Index into the BAM header's `@SQ` lines.
    pub contig: usize,
    pub start: i64,
    pub end: i64,
    pub negative_strand: bool,
    /// In `Gene.transcripts` (`HashMap<String, Transcript>`) iteration order.
    pub transcripts: Vec<Transcript>,
    /// `Interval.hashCode`: contig name, start and end.
    pub hash: i32,
}

/// One contig's genes in `IntervalTree` order (start, then end) with the
/// running maximum end for overlap queries.
#[derive(Clone, Debug, Default)]
struct ContigIndex {
    genes: Vec<usize>,
    max_end: Vec<i64>,
}

/// The loaded gene model: `OverlapDetector<Gene>` plus `getAll()`.
#[derive(Clone, Debug, Default)]
pub struct GeneModel {
    genes: Vec<Gene>,
    by_contig: Vec<ContigIndex>,
}

impl GeneModel {
    /// Every gene kept, in the order they were added to the detector.
    pub fn genes(&self) -> &[Gene] {
        &self.genes
    }

    /// `getOverlaps(new Interval(contig, start, end))`: indices into
    /// [`GeneModel::genes`] of the genes sharing a base with the query, in
    /// the iteration order of the `HashSet` Picard returns (bucket order of
    /// `Gene.hashCode`, ties in tree order). Empty when `end < start`.
    pub fn overlapping(&self, contig: usize, start: i64, end: i64) -> Vec<usize> {
        let Some(index) = self.by_contig.get(contig) else {
            return Vec::new();
        };
        if end < start {
            return Vec::new();
        }
        let last = index.genes.partition_point(|&g| self.genes[g].start <= end);
        let mut found = Vec::new();
        for i in (0..last).rev() {
            if index.max_end[i] < start {
                break;
            }
            let gene = index.genes[i];
            if self.genes[gene].end >= start {
                found.push(gene);
            }
        }
        found.reverse();
        let hashes: Vec<i32> = found.iter().map(|&g| self.genes[g].hash).collect();
        hash_map_order_by_hash(&hashes)
            .into_iter()
            .map(|i| found[i])
            .collect()
    }

    fn from_genes(genes: Vec<Gene>, contigs: usize) -> Self {
        let mut by_contig = vec![ContigIndex::default(); contigs];
        for (i, gene) in genes.iter().enumerate() {
            if let Some(index) = by_contig.get_mut(gene.contig) {
                index.genes.push(i);
            }
        }
        for index in &mut by_contig {
            index.genes.sort_by_key(|&g| (genes[g].start, genes[g].end));
            let mut max_end = i64::MIN;
            index.max_end = index
                .genes
                .iter()
                .map(|&g| {
                    max_end = max_end.max(genes[g].end);
                    max_end
                })
                .collect();
        }
        Self { genes, by_contig }
    }
}

/// One refFlat line after the field split, coordinates still 0-based.
struct Line<'a> {
    fields: Vec<&'a str>,
}

impl Line<'_> {
    fn field(&self, column: usize) -> &str {
        self.fields[column]
    }

    fn int(&self, column: usize) -> std::result::Result<i64, String> {
        parse_int(self.field(column))
    }
}

const GENE_NAME: usize = 0;
const TRANSCRIPT_NAME: usize = 1;
const CHROMOSOME: usize = 2;
const STRAND: usize = 3;
const TX_START: usize = 4;
const TX_END: usize = 5;
const CDS_START: usize = 6;
const CDS_END: usize = 7;
const EXON_COUNT: usize = 8;
const EXON_STARTS: usize = 9;
const EXON_ENDS: usize = 10;
const COLUMNS: usize = 11;

/// `Integer.parseInt`.
fn parse_int(text: &str) -> std::result::Result<i64, String> {
    text.parse::<i32>()
        .map(i64::from)
        .map_err(|_| format!("For input string: \"{text}\""))
}

/// `String.split(",")`: trailing empty strings are dropped, and an empty
/// string splits into one empty field.
fn split_list(text: &str) -> Vec<&str> {
    if text.is_empty() {
        return vec![""];
    }
    let mut parts: Vec<&str> = text.split(',').collect();
    while parts.last() == Some(&"") {
        parts.pop();
    }
    parts
}

/// `RefFlatReader.load`: the model of `path` against the BAM header's
/// dictionary.
pub fn load(path: &Path, header: &sam::Header) -> Result<GeneModel> {
    let text = fs::read_to_string(path).map_err(|source| Error::Input {
        path: path.to_path_buf(),
        source,
    })?;
    let contigs = header.reference_sequences();
    let failure = |message: String| Error::Intervals {
        path: path.to_path_buf(),
        message,
    };
    let mut groups: Vec<(&str, Vec<Line>)> = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').filter(|f| !f.is_empty()).collect();
        if fields.len() != COLUMNS {
            return Err(failure(format!(
                "Wrong number of fields in refFlat file at line {}. Found '{}' but expected '{COLUMNS}'.",
                i + 1,
                fields.len()
            )));
        }
        if contigs
            .get_index_of(fields[CHROMOSOME].as_bytes())
            .is_none()
        {
            continue;
        }
        let name = fields[GENE_NAME];
        let line = Line { fields };
        match groups.iter_mut().find(|(n, _)| *n == name) {
            Some((_, lines)) => lines.push(line),
            None => groups.push((name, vec![line])),
        }
    }
    let names: Vec<&str> = groups.iter().map(|(n, _)| *n).collect();
    let mut genes: Vec<Gene> = Vec::new();
    for i in hash_map_order(&names) {
        let (name, lines) = &groups[i];
        let gene = match make_gene(name, lines, header) {
            Ok(Some(gene)) => gene,
            Ok(None) => continue,
            Err(message) => return Err(failure(message)),
        };
        let same_span = genes.iter().any(|g| {
            g.contig == gene.contig
                && g.start == gene.start
                && g.end == gene.end
                && g.negative_strand == gene.negative_strand
        });
        if !same_span {
            genes.push(gene);
        }
    }
    Ok(GeneModel::from_genes(genes, contigs.len()))
}

/// `makeGeneFromRefFlatLines`: `Ok(None)` for a gene an
/// `AnnotationException` drops, `Err` for what Picard throws out of the
/// load.
fn make_gene(
    name: &str,
    lines: &[Line],
    header: &sam::Header,
) -> std::result::Result<Option<Gene>, String> {
    let first = &lines[0];
    let strand = first.field(STRAND);
    let chromosome = first.field(CHROMOSOME);
    let mut start = i64::MAX;
    let mut end = i64::MIN;
    for line in lines {
        start = start.min(line.int(TX_START)? + 1);
        end = end.max(line.int(TX_END)?);
    }
    let mut transcripts: Vec<Transcript> = Vec::new();
    for line in lines {
        if line.field(STRAND) != strand || line.field(CHROMOSOME) != chromosome {
            return Ok(None);
        }
        let Some(transcript) = make_transcript(line)? else {
            return Ok(None);
        };
        if transcripts.iter().any(|t| t.name == transcript.name) {
            return Ok(None);
        }
        transcripts.push(transcript);
    }
    let transcript_names: Vec<&str> = transcripts.iter().map(|t| t.name.as_str()).collect();
    let order = hash_map_order(&transcript_names);
    let mut ordered: Vec<Option<Transcript>> = transcripts.into_iter().map(Some).collect();
    let transcripts = order
        .into_iter()
        .map(|i| ordered[i].take().expect("each index once"))
        .collect();
    let contig = header
        .reference_sequences()
        .get_index_of(chromosome.as_bytes())
        .expect("checked while reading");
    let hash = string_hash_code(chromosome)
        .wrapping_mul(31)
        .wrapping_add(start as i32)
        .wrapping_mul(31)
        .wrapping_add(end as i32);
    Ok(Some(Gene {
        name: name.to_owned(),
        contig,
        start,
        end,
        negative_strand: strand == "-",
        transcripts,
        hash,
    }))
}

/// `makeTranscriptFromRefFlatLine`: `Ok(None)` where Picard drops the gene.
fn make_transcript(line: &Line) -> std::result::Result<Option<Transcript>, String> {
    let exon_count = line.int(EXON_COUNT)?;
    let starts = split_list(line.field(EXON_STARTS));
    let ends = split_list(line.field(EXON_ENDS));
    if exon_count != starts.len() as i64 || exon_count != ends.len() as i64 {
        return Ok(None);
    }
    let transcription_start = line.int(TX_START)? + 1;
    let transcription_end = line.int(TX_END)?;
    let coding_start = line.int(CDS_START)? + 1;
    let coding_end = line.int(CDS_END)?;
    let mut exons: Vec<Exon> = Vec::with_capacity(starts.len());
    let mut length = 0;
    for (start, end) in starts.iter().zip(&ends) {
        let exon = Exon {
            start: parse_int(start)? + 1,
            end: parse_int(end)?,
        };
        length += exon.end - exon.start + 1;
        if exon.start > exon.end {
            return Ok(None);
        }
        if exons.last().is_some_and(|last| last.end >= exon.start) {
            return Ok(None);
        }
        exons.push(exon);
    }
    let name = line.field(TRANSCRIPT_NAME);
    let hash = [
        transcription_start,
        transcription_end,
        coding_start,
        coding_end,
    ]
    .into_iter()
    .fold(string_hash_code(name), |h, v| {
        h.wrapping_mul(31).wrapping_add(v as i32)
    });
    Ok(Some(Transcript {
        name: name.to_owned(),
        transcription_start,
        transcription_end,
        coding_start,
        coding_end,
        exons,
        length,
        hash,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header() -> sam::Header {
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10000\n@SQ\tSN:c2\tLN:5000\n"
            .parse()
            .expect("header")
    }

    /// `text` loaded from a file of its own (tests run in parallel).
    fn load_text(text: &str) -> Result<GeneModel> {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "bammetrics_ref_flat_{}_{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::write(&path, text).expect("write");
        let model = load(&path, &header());
        fs::remove_file(&path).expect("remove");
        model
    }

    const TWO_GENES: &str = "G1\tT1\tc1\t+\t999\t2000\t1199\t1800\t2\t999,1499,\t1200,2000,\n\
        G1\tT2\tc1\t+\t999\t2200\t1199\t1800\t2\t999,1899,\t1200,2200,\n\
        G2\tT3\tc1\t-\t3000\t3500\t3500\t3500\t1\t3000,\t3500,\n\
        G3\tT4\tchrUn\t+\t1\t100\t1\t100\t1\t1,\t100,\n";

    #[test]
    fn genes_span_their_transcripts_and_unknown_contigs_are_dropped() {
        let model = load_text(TWO_GENES).expect("loads");
        let mut names: Vec<&str> = model.genes().iter().map(|g| g.name.as_str()).collect();
        names.sort_unstable();
        assert_eq!(names, ["G1", "G2"]);
        let g1 = model.genes().iter().find(|g| g.name == "G1").expect("G1");
        assert_eq!((g1.start, g1.end, g1.negative_strand), (1000, 2200, false));
        assert_eq!(g1.transcripts.len(), 2);
        let t1 = g1.transcripts.iter().find(|t| t.name == "T1").expect("T1");
        assert_eq!(t1.length, 201 + 501);
        assert_eq!(
            t1.exons[1],
            Exon {
                start: 1500,
                end: 2000
            }
        );
    }

    #[test]
    fn locus_functions_follow_the_strongest_transcript() {
        let model = load_text(TWO_GENES).expect("loads");
        let g1 = model.genes().iter().find(|g| g.name == "G1").expect("G1");
        let mut functions = vec![LocusFunction::Intergenic; 1300];
        for t in &g1.transcripts {
            t.assign_locus_function(900, &mut functions);
        }
        let at = |locus: i64| functions[(locus - 900) as usize];
        assert_eq!(at(999), LocusFunction::Intergenic);
        assert_eq!(at(1000), LocusFunction::Utr);
        assert_eq!(at(1200), LocusFunction::Coding);
        assert_eq!(at(1300), LocusFunction::Intronic);
        assert_eq!(at(1500), LocusFunction::Coding);
        assert_eq!(at(1801), LocusFunction::Utr);
        // 2100 is intronic in T1 (past its end: nothing) and exonic in T2.
        assert_eq!(at(2100), LocusFunction::Utr);
        // The non-coding G2 is all UTR.
        let g2 = model.genes().iter().find(|g| g.name == "G2").expect("G2");
        let mut functions = vec![LocusFunction::Intergenic; 10];
        g2.transcripts[0].assign_locus_function(3001, &mut functions);
        assert!(functions.iter().all(|f| *f == LocusFunction::Utr));
    }

    #[test]
    fn coverage_skips_the_end_locus_and_maps_through_exons() {
        let model = load_text(TWO_GENES).expect("loads");
        let g1 = model.genes().iter().find(|g| g.name == "G1").expect("G1");
        let t1 = g1.transcripts.iter().find(|t| t.name == "T1").expect("T1");
        let mut coverage = vec![0; t1.length as usize];
        // A block 1195..=1504: 1195 to 1503 count; 1201..1499 are intronic.
        t1.add_coverage(1195, 1504, &mut coverage);
        assert_eq!(coverage[194..201].iter().sum::<i32>(), 6);
        assert_eq!(coverage[201..205], [1, 1, 1, 1]);
        assert_eq!(coverage.iter().sum::<i32>(), 10);
        assert_eq!(t1.transcript_coordinate(1500), 202);
        assert_eq!(t1.transcript_coordinate(1300), -1);
    }

    #[test]
    fn overlap_queries_find_nested_genes() {
        let model = load_text(TWO_GENES).expect("loads");
        assert_eq!(model.overlapping(0, 2150, 2300).len(), 1);
        assert_eq!(model.overlapping(0, 2201, 2999).len(), 0);
        assert_eq!(model.overlapping(0, 500, 5000).len(), 2);
        assert_eq!(model.overlapping(0, 3500, 3400).len(), 0);
        assert_eq!(model.overlapping(1, 1, 5000).len(), 0);
        assert_eq!(model.overlapping(7, 1, 5000).len(), 0);
    }

    #[test]
    fn bad_genes_are_dropped_and_bad_fields_are_errors() {
        let text = "A\tt1\tc1\t+\t10\t20\t10\t20\t1\t10,\t20,\n\
            A\tt2\tc1\t-\t10\t20\t10\t20\t1\t10,\t20,\n\
            B\tt3\tc1\t+\t10\t20\t10\t20\t2\t10,\t20,\n\
            C\tt4\tc1\t+\t10\t20\t10\t20\t2\t10,15,\t16,20,\n\
            D\tt5\tc1\t+\t10\t20\t10\t20\t1\t15,\t15,\n\
            E\tt6\tc1\t+\t10\t20\t10\t20\t1\t10,\t20,\n\
            E\tt6\tc1\t+\t10\t20\t10\t20\t1\t10,\t20,\n\
            F\tt7\tc1\t+\t30\t40\t30\t40\t1\t30,\t40,\n\
            G\tt8\tc1\t+\t30\t40\t30\t40\t1\t30,\t40,\n";
        let model = load_text(text).expect("loads");
        let names: Vec<&str> = model.genes().iter().map(|g| g.name.as_str()).collect();
        assert_eq!(names.len(), 1);
        assert!(names == ["F"] || names == ["G"]);
        assert!(load_text("A\tt1\tc1\t+\t10\t20\t10\t20\t1\t10,\n").is_err());
        assert!(load_text("A\tt1\tc1\t+\tx\t20\t10\t20\t1\t10,\t20,\n").is_err());
        assert!(
            load_text("# comment\n\n")
                .expect("empty")
                .genes()
                .is_empty()
        );
    }

    #[test]
    fn transcript_hash_matches_java() {
        // "T1".hashCode() = 2653, then 31 * h + coordinate four times.
        let model = load_text(TWO_GENES).expect("loads");
        let g1 = model.genes().iter().find(|g| g.name == "G1").expect("G1");
        let t1 = g1.transcripts.iter().find(|t| t.name == "T1").expect("T1");
        let expected = [1000i32, 2000, 1200, 1800]
            .into_iter()
            .fold(2653i32, |h, v| h.wrapping_mul(31).wrapping_add(v));
        assert_eq!(t1.hash, expected);
    }
}

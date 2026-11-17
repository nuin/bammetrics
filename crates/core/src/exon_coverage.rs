//! `CollectExonCoverageMetrics` (docs/07 wave W07): the depth of every
//! exon of an annotation, one row per exon of each selected transcript,
//! and in `SUMMARY_OUTPUT` one row per gene and per transcript. Not a
//! Picard tool.
//!
//! The depth of a base is the high-quality depth `CollectHsMetrics` writes
//! to `PER_BASE_COVERAGE`: the tool runs the Hs accumulator
//! ([`CaptureMetrics`]) with the selected exons as baits and targets and
//! reads its `ALL_READS` depth arrays ([`CaptureMetrics::into_target_depths`]),
//! as `CollectCoverageGaps` does. So an exon's `MEAN_DEPTH` and
//! `MINIMUM_DEPTH` are the `mean_coverage` and `min_coverage` that
//! `CollectHsMetrics PER_TARGET_COVERAGE` prints for the same interval
//! (`tests/concordance/test_exon_coverage.py` holds both tools to that).
//! The Hs filters apply: duplicates, secondary and QC-failed reads, MAPQ
//! under `MINIMUM_MAPPING_QUALITY`, bases under `MINIMUM_BASE_QUALITY` or
//! at quality 2 or less, the overlap of two mates counted once.
//! `COVERAGE_CAP` caps a Hs histogram only, never these depths, so it is
//! not an option.
//!
//! Exons of different transcripts may overlap; the depth array is per
//! merged interval, and each exon reads its own bases from it. An exon on a
//! contig the BAM header lacks, under its own name and under the `chr`
//! alias, gets a row with empty depth cells and a note on standard error.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;

use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::gene_model::{self, Feature, GeneModel, Selected, Selection};
use crate::hs_metrics::{self, CaptureMetrics, TargetDepths};
use crate::intervals::{Interval, NamedInterval};
use crate::levels::Level;
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::pass::{self, Accumulator};
use crate::registry::{Invocation, Opt, Tool};

/// The metrics class of the exon file. Not a Picard class.
pub const CLASS: &str = "bammetrics.qc.ExonCoverageMetrics";

/// The metrics class of the gene and transcript file.
pub const SUMMARY_CLASS: &str = "bammetrics.qc.GeneCoverageMetrics";

pub const DEFAULT_MIN_DEPTHS: [i32; 3] = [10, 20, 30];

/// The `MIN_DEPTH` values one run takes.
pub const MAX_MIN_DEPTHS: usize = 16;

/// The exon columns before the `BASES_<N>X` and `PCT_<N>X` pairs.
pub const EXON_COLUMNS: [&str; 13] = [
    "GENE_NAME",
    "GENE_ID",
    "TRANSCRIPT_ID",
    "EXON_NUMBER",
    "CHROM",
    "START",
    "END",
    "STRAND",
    "LENGTH",
    "MEAN_DEPTH",
    "MEDIAN_DEPTH",
    "MINIMUM_DEPTH",
    "MAXIMUM_DEPTH",
];

/// The gene and transcript columns before the `BASES_<N>X` and `PCT_<N>X`
/// columns.
pub const SUMMARY_COLUMNS: [&str; 14] = [
    "LEVEL",
    "GENE_NAME",
    "GENE_ID",
    "TRANSCRIPT_ID",
    "CHROM",
    "START",
    "END",
    "STRAND",
    "EXON_COUNT",
    "LENGTH",
    "MEAN_DEPTH",
    "MEDIAN_DEPTH",
    "MINIMUM_DEPTH",
    "MAXIMUM_DEPTH",
];

/// Column lists by class and depths. A [`MetricsOutput`] names its columns
/// with `'static` strings; each distinct list is built once per process, so
/// at most one list per `MIN_DEPTH` set a process uses is kept.
type ColumnCache = HashMap<(&'static str, Vec<i32>), &'static [&'static str]>;
static COLUMNS: Mutex<Option<ColumnCache>> = Mutex::new(None);

/// The columns of `class` for `depths`: the fixed ones, then `BASES_<N>X`
/// for each depth, then `PCT_<N>X` for each.
pub fn columns(class: &'static str, depths: &[i32]) -> &'static [&'static str] {
    let fixed: &[&'static str] = if class == CLASS {
        &EXON_COLUMNS
    } else {
        &SUMMARY_COLUMNS
    };
    let mut cache = COLUMNS
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let cache = cache.get_or_insert_with(HashMap::new);
    cache.entry((class, depths.to_vec())).or_insert_with(|| {
        let mut list: Vec<&'static str> = fixed.to_vec();
        for prefix in ["BASES", "PCT"] {
            for depth in depths {
                list.push(Box::leak(format!("{prefix}_{depth}X").into_boxed_str()));
            }
        }
        Box::leak(list.into_boxed_slice())
    })
}

/// The options of one run.
#[derive(Clone, Debug, PartialEq)]
pub struct Options {
    pub annotation: PathBuf,
    pub feature: Feature,
    pub genes: Vec<String>,
    pub gene_list: Option<PathBuf>,
    pub transcripts: Vec<String>,
    pub canonical: bool,
    /// Ascending, distinct, each at least 1.
    pub min_depths: Vec<i32>,
    pub summary_output: Option<PathBuf>,
    /// The Hs filter options: `MINIMUM_MAPPING_QUALITY`,
    /// `MINIMUM_BASE_QUALITY`, `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`.
    pub hs: hs_metrics::Options,
}

impl Options {
    /// The options alone: no file is opened, so a chain checks them before
    /// its first read.
    pub(crate) fn from_invocation(inv: &Invocation) -> Result<Self> {
        let annotation = inv.get("ANNOTATION").map(PathBuf::from).ok_or_else(|| {
            Error::Usage("Option 'ANNOTATION' must be specified at least 1 times.".to_owned())
        })?;
        let mut min_depths =
            inv.all("MIN_DEPTH")
                .iter()
                .map(|text| {
                    text.parse::<i32>().ok().filter(|&d| d >= 1).ok_or_else(|| {
                        Error::Usage(format!("MIN_DEPTH={text}: a depth of 1 or more"))
                    })
                })
                .collect::<Result<Vec<i32>>>()?;
        if min_depths.is_empty() {
            min_depths = DEFAULT_MIN_DEPTHS.to_vec();
        }
        min_depths.sort_unstable();
        min_depths.dedup();
        if min_depths.len() > MAX_MIN_DEPTHS {
            return Err(Error::Usage(format!(
                "{} MIN_DEPTH values; at most {MAX_MIN_DEPTHS}",
                min_depths.len()
            )));
        }
        let defaults = hs_metrics::Options::default();
        let hs = hs_metrics::Options {
            levels: vec![Level::AllReads],
            minimum_mapping_quality: inv
                .parse("MINIMUM_MAPPING_QUALITY", defaults.minimum_mapping_quality)?,
            minimum_base_quality: inv
                .parse("MINIMUM_BASE_QUALITY", defaults.minimum_base_quality)?,
            clip_overlapping_reads: inv
                .flag("CLIP_OVERLAPPING_READS", defaults.clip_overlapping_reads)?,
            include_indels: inv.flag("INCLUDE_INDELS", false)?,
            ..defaults
        };
        Ok(Self {
            annotation,
            feature: Feature::parse(inv.get("FEATURE").unwrap_or("exon"))?,
            genes: inv.all("GENE").to_vec(),
            gene_list: inv.get("GENE_LIST").map(PathBuf::from),
            transcripts: inv.all("TRANSCRIPT").to_vec(),
            canonical: inv.flag("CANONICAL", false)?,
            min_depths,
            summary_output: inv.get("SUMMARY_OUTPUT").map(PathBuf::from),
            hs,
        })
    }

    /// The selection: `GENE` and the names of `GENE_LIST`, `TRANSCRIPT`,
    /// `CANONICAL`. Reads the gene list.
    fn selection(&self) -> Result<Selection> {
        let mut genes = self.genes.clone();
        if let Some(list) = &self.gene_list {
            let names = gene_model::read_gene_list(list)?;
            if names.is_empty() {
                return Err(Error::Usage(format!(
                    "GENE_LIST={}: no gene name in the file",
                    list.display()
                )));
            }
            genes.extend(names);
        }
        Ok(Selection {
            genes,
            transcripts: self.transcripts.clone(),
            canonical: self.canonical,
        })
    }
}

/// The BAM contig of each annotation contig: its own name, else the name
/// with `chr` added or removed (`chrM` and `MT` too).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ContigMap {
    /// By annotation contig index.
    pub index: Vec<Option<usize>>,
    /// (annotation contig index, BAM name) of each contig read under an
    /// alias.
    pub aliased: Vec<(usize, String)>,
}

impl ContigMap {
    pub fn new(names: &[String], header: &sam::Header) -> Self {
        let contigs = header.reference_sequences();
        let mut map = Self::default();
        for (index, name) in names.iter().enumerate() {
            let found = contigs.get_index_of(name.as_bytes()).or_else(|| {
                aliases(name).into_iter().find_map(|alias| {
                    let i = contigs.get_index_of(alias.as_bytes())?;
                    map.aliased.push((index, alias));
                    Some(i)
                })
            });
            map.index.push(found);
        }
        map
    }
}

/// The other spellings of a contig name: `chr` added or removed, and the
/// mitochondrion as `chrM` and `MT`.
fn aliases(name: &str) -> Vec<String> {
    let mut names = match name.strip_prefix("chr") {
        Some(rest) if !rest.is_empty() => vec![rest.to_owned()],
        _ => vec![format!("chr{name}")],
    };
    match name {
        "chrM" => names.push("MT".to_owned()),
        "MT" => names.push("chrM".to_owned()),
        _ => {}
    }
    names
}

/// The selected exons placed on the BAM header.
struct Plan {
    model: GeneModel,
    exons: Vec<Selected>,
    /// BAM contig of each exon, `None` when the header lacks it.
    placed: Vec<Option<usize>>,
    min_depths: Vec<i32>,
}

impl Plan {
    fn new(opts: &Options, header: &sam::Header) -> Result<(Self, Vec<String>)> {
        let model = GeneModel::read(&opts.annotation, opts.feature)?;
        let exons = model.select(&opts.selection()?)?;
        let map = ContigMap::new(&model.contigs, header);
        let contigs = header.reference_sequences();
        let mut placed = Vec::with_capacity(exons.len());
        let mut used = vec![false; model.contigs.len()];
        for exon in &exons {
            used[exon.contig] = true;
            let bam = map.index[exon.contig];
            if let Some(i) = bam {
                let length = contigs[i].length().get() as i64;
                if exon.end > length {
                    return Err(Error::Intervals {
                        path: opts.annotation.clone(),
                        message: format!(
                            "{} {}:{}-{} of {} extends beyond the end of {} ({length} bp in \
                             the BAM header); is the annotation for another build?",
                            opts.feature.name(),
                            model.contigs[exon.contig],
                            exon.start,
                            exon.end,
                            model.transcripts[exon.transcript].id,
                            String::from_utf8_lossy(
                                contigs.get_index(i).map_or(&[][..], |(n, _)| n.as_ref())
                            ),
                        ),
                    });
                }
            }
            placed.push(bam);
        }
        let mut notes = Vec::new();
        let aliased: Vec<String> = map
            .aliased
            .iter()
            .filter(|(index, _)| used[*index])
            .map(|(index, bam)| format!("{} as {bam}", model.contigs[*index]))
            .collect();
        if !aliased.is_empty() {
            notes.push(format!(
                "CollectExonCoverageMetrics: {} annotation contig{} read under the BAM header's \
                 name: {}",
                aliased.len(),
                if aliased.len() == 1 { "" } else { "s" },
                shortened(&aliased)
            ));
        }
        let missing: Vec<String> = (0..model.contigs.len())
            .filter(|&i| used[i] && map.index[i].is_none())
            .map(|i| model.contigs[i].clone())
            .collect();
        let unplaced = placed.iter().filter(|p| p.is_none()).count();
        if unplaced == exons.len() {
            return Err(Error::Intervals {
                path: opts.annotation.clone(),
                message: format!(
                    "none of the {} selected {}s is on a contig of the BAM header (annotation \
                     contigs {}); is the annotation for another build?",
                    exons.len(),
                    opts.feature.name(),
                    shortened(&missing)
                ),
            });
        }
        if unplaced > 0 {
            notes.push(format!(
                "CollectExonCoverageMetrics: {unplaced} {}{} on {} contig{} the BAM header \
                 lacks ({}): written with empty depth cells, not counted",
                opts.feature.name(),
                if unplaced == 1 { "" } else { "s" },
                missing.len(),
                if missing.len() == 1 { "" } else { "s" },
                shortened(&missing)
            ));
        }
        let plan = Self {
            model,
            exons,
            placed,
            min_depths: opts.min_depths.clone(),
        };
        Ok((plan, notes))
    }

    /// The placed exons as Hs targets, each named by its transcript.
    fn targets(&self) -> Vec<NamedInterval> {
        self.exons
            .iter()
            .zip(&self.placed)
            .filter_map(|(exon, placed)| {
                Some(NamedInterval {
                    interval: Interval {
                        contig: (*placed)?,
                        start: exon.start,
                        end: exon.end,
                    },
                    negative_strand: exon.strand == gene_model::Strand::Reverse,
                    name: self.model.transcripts[exon.transcript].id.clone(),
                })
            })
            .collect()
    }
}

/// Up to five names, then how many more.
fn shortened(names: &[String]) -> String {
    let shown = names.iter().take(5).cloned().collect::<Vec<_>>().join(", ");
    match names.len().saturating_sub(5) {
        0 => shown,
        more => format!("{shown} and {more} more"),
    }
}

/// The depth summary of a set of bases.
#[derive(Clone, Debug, PartialEq)]
pub struct Stats {
    pub length: i64,
    pub mean: f64,
    /// The middle depth; with an even count, the mean of the two middle
    /// depths.
    pub median: f64,
    pub minimum: i32,
    pub maximum: i32,
    /// Bases at or above each `MIN_DEPTH`.
    pub at: Vec<i64>,
}

impl Stats {
    /// The summary of `depths` (at least one base). `depths` is reordered.
    pub fn of(depths: &mut [i32], min_depths: &[i32]) -> Self {
        let length = depths.len() as i64;
        let sum: i64 = depths.iter().map(|&d| i64::from(d)).sum();
        let minimum = depths.iter().copied().min().unwrap_or(0);
        let maximum = depths.iter().copied().max().unwrap_or(0);
        let at = min_depths
            .iter()
            .map(|&m| depths.iter().filter(|&&d| d >= m).count() as i64)
            .collect();
        let half = depths.len() / 2;
        let odd = depths.len() % 2 == 1;
        let median = if depths.is_empty() {
            f64::NAN
        } else {
            let (low, &mut upper, _) = depths.select_nth_unstable(half);
            if odd {
                f64::from(upper)
            } else {
                let lower = low.iter().copied().max().unwrap_or(upper);
                (f64::from(lower) + f64::from(upper)) / 2.0
            }
        };
        Self {
            length,
            // `sum / length` in double, as the PER_TARGET_COVERAGE line divides.
            mean: sum as f64 / length as f64,
            median,
            minimum,
            maximum,
            at,
        }
    }

    fn values(&self, out: &mut Vec<Value>) {
        out.push(Value::Float(self.mean));
        out.push(Value::Float(self.median));
        out.push(Value::Int(i64::from(self.minimum)));
        out.push(Value::Int(i64::from(self.maximum)));
        out.extend(self.at.iter().map(|&n| Value::Int(n)));
        out.extend(
            self.at
                .iter()
                .map(|&n| Value::Float(n as f64 / self.length as f64)),
        );
    }

    /// The depth cells of an interval with no depth: empty.
    fn empty(depths: usize, out: &mut Vec<Value>) {
        out.extend(std::iter::repeat_n(Value::Null, 4 + 2 * depths));
    }
}

/// The depth of `[start, end]` on BAM contig `contig`, from the merged
/// target that holds it.
fn depth_of(depths: &TargetDepths, contig: usize, start: i64, end: i64) -> Result<&[i32]> {
    let missing = || Error::Metric(format!("no merged target holds {contig}:{start}-{end}"));
    let targets = depths.targets();
    let index = targets.partition_point(|t| (t.contig, t.end) < (contig, start));
    let target = targets
        .get(index)
        .filter(|t| t.contig == contig && t.start <= start && end <= t.end)
        .ok_or_else(missing)?;
    let from = usize::try_from(start - target.start).map_err(|_| missing())?;
    let to = usize::try_from(end - target.start + 1).map_err(|_| missing())?;
    depths.depth(index).get(from..to).ok_or_else(missing)
}

/// One gene or transcript row: its exons (indices into the plan).
struct Group {
    gene: usize,
    transcript: Option<usize>,
    contig: usize,
    exons: Vec<usize>,
}

/// The two files of one run.
pub struct Report {
    pub exons: MetricsOutput,
    pub summary: MetricsOutput,
}

impl Report {
    fn new(plan: &Plan, depths: &TargetDepths) -> Result<Self> {
        let depth_count = plan.min_depths.len();
        let mut scratch: Vec<i32> = Vec::new();
        let mut rows = Vec::with_capacity(plan.exons.len());
        for (exon, placed) in plan.exons.iter().zip(&plan.placed) {
            let transcript = &plan.model.transcripts[exon.transcript];
            let gene = &plan.model.genes[transcript.gene];
            let mut row = vec![
                Value::from(gene.name.as_str()),
                Value::from(gene.id.as_str()),
                Value::from(transcript.id.as_str()),
                Value::Int(exon.number),
                Value::from(chrom(plan, depths, exon.contig, *placed)),
                Value::Int(exon.start),
                Value::Int(exon.end),
                Value::from(exon.strand.symbol()),
                Value::Int(exon.length()),
            ];
            match placed {
                Some(contig) => {
                    scratch.clear();
                    scratch.extend_from_slice(depth_of(depths, *contig, exon.start, exon.end)?);
                    Stats::of(&mut scratch, &plan.min_depths).values(&mut row);
                }
                None => Stats::empty(depth_count, &mut row),
            }
            rows.push(row);
        }
        let exons = MetricsOutput {
            rows,
            ..output(CLASS, &plan.min_depths)
        };
        let mut summary_rows = Vec::new();
        for group in groups(plan) {
            summary_rows.push(summary_row(plan, depths, &group, &mut scratch)?);
        }
        let summary = MetricsOutput {
            rows: summary_rows,
            ..output(SUMMARY_CLASS, &plan.min_depths)
        };
        Ok(Self { exons, summary })
    }
}

fn output(class: &'static str, depths: &[i32]) -> MetricsOutput {
    MetricsOutput {
        class,
        columns: columns(class, depths),
        rows: Vec::new(),
        histogram: None,
        with_headers: true,
        text: None,
    }
}

/// The BAM's name of a placed contig, else the annotation's.
fn chrom<'a>(
    plan: &'a Plan,
    depths: &'a TargetDepths,
    contig: usize,
    placed: Option<usize>,
) -> &'a str {
    match placed {
        Some(i) => depths.contig_name(i),
        None => &plan.model.contigs[contig],
    }
}

/// The gene rows in order of their first exon, each followed by its
/// transcript rows in the same order. A gene or transcript on two contigs
/// has a row per contig.
fn groups(plan: &Plan) -> Vec<Group> {
    let mut genes: Vec<(Group, Vec<Group>)> = Vec::new();
    let mut gene_index: HashMap<(usize, usize), usize> = HashMap::new();
    let mut transcript_index: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
    for (i, exon) in plan.exons.iter().enumerate() {
        let gene = plan.model.transcripts[exon.transcript].gene;
        let g = *gene_index.entry((gene, exon.contig)).or_insert_with(|| {
            genes.push((
                Group {
                    gene,
                    transcript: None,
                    contig: exon.contig,
                    exons: Vec::new(),
                },
                Vec::new(),
            ));
            genes.len() - 1
        });
        genes[g].0.exons.push(i);
        let (g, t) = *transcript_index
            .entry((exon.transcript, exon.contig))
            .or_insert_with(|| {
                genes[g].1.push(Group {
                    gene,
                    transcript: Some(exon.transcript),
                    contig: exon.contig,
                    exons: Vec::new(),
                });
                (g, genes[g].1.len() - 1)
            });
        genes[g].1[t].exons.push(i);
    }
    genes
        .into_iter()
        .flat_map(|(gene, transcripts)| std::iter::once(gene).chain(transcripts))
        .collect()
}

/// A gene or transcript row: the union of its exons' bases, so a base two
/// transcripts share counts once in the gene row.
fn summary_row(
    plan: &Plan,
    depths: &TargetDepths,
    group: &Group,
    scratch: &mut Vec<i32>,
) -> Result<Vec<Value>> {
    let gene = &plan.model.genes[group.gene];
    let mut spans: Vec<(i64, i64)> = group
        .exons
        .iter()
        .map(|&i| (plan.exons[i].start, plan.exons[i].end))
        .collect();
    spans.sort_unstable();
    let exon_count = {
        let mut distinct = spans.clone();
        distinct.dedup();
        distinct.len() as i64
    };
    let mut union: Vec<(i64, i64)> = Vec::with_capacity(spans.len());
    for (start, end) in spans {
        match union.last_mut() {
            Some(last) if start <= last.1 => last.1 = last.1.max(end),
            _ => union.push((start, end)),
        }
    }
    let first = &plan.exons[group.exons[0]];
    let placed = plan.placed[group.exons[0]];
    let length: i64 = union.iter().map(|(s, e)| e - s + 1).sum();
    let mut row = vec![
        Value::from(if group.transcript.is_some() {
            "TRANSCRIPT"
        } else {
            "GENE"
        }),
        Value::from(gene.name.as_str()),
        Value::from(gene.id.as_str()),
        Value::from(
            group
                .transcript
                .map_or("", |t| plan.model.transcripts[t].id.as_str()),
        ),
        Value::from(chrom(plan, depths, group.contig, placed)),
        Value::Int(union.first().map_or(0, |u| u.0)),
        Value::Int(union.last().map_or(0, |u| u.1)),
        Value::from(first.strand.symbol()),
        Value::Int(exon_count),
        Value::Int(length),
    ];
    match placed {
        Some(contig) => {
            scratch.clear();
            for &(start, end) in &union {
                scratch.extend_from_slice(depth_of(depths, contig, start, end)?);
            }
            Stats::of(scratch, &plan.min_depths).values(&mut row);
        }
        None => Stats::empty(plan.min_depths.len(), &mut row),
    }
    Ok(row)
}

/// The Hs accumulator over the selected exons.
pub struct ExonCoverage {
    capture: CaptureMetrics,
    plan: Plan,
}

impl ExonCoverage {
    /// Reads the annotation and the gene list and places the exons on the
    /// header. The notes are for standard error.
    pub fn new(opts: &Options, header: &sam::Header) -> Result<(Self, Vec<String>)> {
        let (plan, notes) = Plan::new(opts, header)?;
        let targets = plan.targets();
        let baits: Vec<Interval> = targets.iter().map(|t| t.interval).collect();
        let capture = CaptureMetrics::with_intervals(&opts.hs, header, &targets, &baits)?;
        Ok((Self { capture, plan }, notes))
    }
}

impl Accumulator for ExonCoverage {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        self.capture.accept(record)
    }

    fn finish(self) -> Result<Report> {
        let depths = self.capture.into_target_depths()?;
        Report::new(&self.plan, &depths)
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectExonCoverageMetrics",
    summary: "one ExonCoverageMetrics row per exon of each selected transcript: mean, median, \
              minimum and maximum depth and the bases and fraction at or above each MIN_DEPTH, \
              with the depth CollectHsMetrics writes per base; SUMMARY_OUTPUT adds one \
              GeneCoverageMetrics row per gene and per transcript. Not a Picard tool \
              (docs/07 W07).",
    options: &[
        Opt::path(
            &["ANNOTATION"],
            "The exons: a GTF or GFF3 file (its exon or CDS features), a BED12 file (blocks as \
             exons, the name as the transcript) or a BED file (each line a region); plain or gzip.",
        )
        .required(),
        Opt::text(
            &["GENE"],
            "Only this gene, by gene_name or gene_id (the version may be left out); repeat for several.",
        )
        .repeated(),
        Opt::path(
            &["GENE_LIST"],
            "A file of gene names or IDs, one per line, added to GENE.",
        ),
        Opt::text(
            &["TRANSCRIPT"],
            "Only this transcript ID (the version may be left out); repeat for several.",
        )
        .repeated(),
        Opt::bool(
            &["CANONICAL"],
            "One transcript per gene: tagged Ensembl_canonical, else MANE_Select, else the \
             longest CDS (default false).",
        ),
        Opt::text(
            &["FEATURE"],
            "exon (default) or CDS: the GTF or GFF3 feature reported; CDS of a BED12 file is \
             the thick part of its blocks.",
        ),
        Opt::int(
            &["MIN_DEPTH"],
            "Count the bases at or above this depth; repeat for several, which replace the \
             default 10, 20 and 30.",
        )
        .repeated(),
        Opt::path(
            &["SUMMARY_OUTPUT"],
            "Also write the GeneCoverageMetrics file: one row per gene, then one per transcript.",
        )
        .written(),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Reads below this MAPQ add no depth (default 20).",
        ),
        Opt::int(
            &["MINIMUM_BASE_QUALITY", "Q"],
            "Bases below this quality add no depth (default 20).",
        ),
        Opt::bool(
            &["CLIP_OVERLAPPING_READS"],
            "Count the overlap of two mates once (default true).",
        ),
        Opt::bool(
            &["INCLUDE_INDELS"],
            "Count deleted bases as covered (default false).",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    let (coverage, notes) = ExonCoverage::new(&opts, header)?;
    for note in notes.iter().filter(|_| !inv.quiet) {
        eprintln!("bammetrics: {note}");
    }
    Ok(pass::boxed(coverage, move |report| {
        let mut output = ToolOutput::file(report.exons);
        if let Some(path) = &opts.summary_output {
            output.extra.push((path.clone(), report.summary));
        }
        Ok(output)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn options(args: &[&str]) -> Result<Options> {
        TOOL.parse(&strings(args))
            .and_then(|inv| Options::from_invocation(&inv))
    }

    #[test]
    fn stats_count_mean_median_and_thresholds() {
        let mut depths = vec![5, 30, 10, 20];
        let stats = Stats::of(&mut depths, &[10, 20, 31]);
        assert_eq!(stats.length, 4);
        assert_eq!(stats.mean, 16.25);
        assert_eq!(stats.median, 15.0);
        assert_eq!((stats.minimum, stats.maximum), (5, 30));
        assert_eq!(stats.at, [3, 2, 0]);
        let mut one = vec![7];
        let stats = Stats::of(&mut one, &[7]);
        assert_eq!((stats.median, stats.at[0]), (7.0, 1));
        let mut odd = vec![3, 1, 2];
        assert_eq!(Stats::of(&mut odd, &[]).median, 2.0);
    }

    #[test]
    fn columns_follow_the_depths_and_are_built_once() {
        let a = columns(CLASS, &[10, 20]);
        assert_eq!(&a[13..], ["BASES_10X", "BASES_20X", "PCT_10X", "PCT_20X"]);
        assert!(std::ptr::eq(a, columns(CLASS, &[10, 20])));
        let s = columns(SUMMARY_CLASS, &[5]);
        assert_eq!(s[0], "LEVEL");
        assert_eq!(&s[14..], ["BASES_5X", "PCT_5X"]);
    }

    #[test]
    fn aliases_add_or_remove_chr() {
        assert_eq!(aliases("chr1"), ["1"]);
        assert_eq!(aliases("1"), ["chr1"]);
        assert_eq!(aliases("chrM"), ["M", "MT"]);
        assert_eq!(aliases("MT"), ["chrMT", "chrM"]);
        assert_eq!(aliases("chr"), ["chrchr"]);
    }

    #[test]
    fn options_need_an_annotation_and_take_depths_in_place_of_the_default() {
        let err = options(&["I=x.bam"]).expect_err("no annotation");
        assert!(err.to_string().contains("ANNOTATION"), "{err}");
        let opts = options(&["I=x.bam", "ANNOTATION=a.gtf"]).expect("defaults");
        assert_eq!(opts.min_depths, DEFAULT_MIN_DEPTHS);
        assert_eq!(opts.feature, Feature::Exon);
        assert_eq!(opts.hs.minimum_mapping_quality, 20);
        let opts = options(&[
            "I=x.bam",
            "ANNOTATION=a.gtf",
            "MIN_DEPTH=50",
            "MIN_DEPTH=5",
            "MIN_DEPTH=50",
            "FEATURE=cds",
            "GENE=BRCA1",
            "GENE=TP53",
            "CANONICAL=true",
            "MQ=1",
        ])
        .expect("options");
        assert_eq!(opts.min_depths, [5, 50]);
        assert_eq!(opts.feature, Feature::Cds);
        assert_eq!(opts.genes, ["BRCA1", "TP53"]);
        assert!(opts.canonical);
        assert_eq!(opts.hs.minimum_mapping_quality, 1);
        for (bad, text) in [
            ("MIN_DEPTH=0", "MIN_DEPTH=0"),
            ("MIN_DEPTH=x", "MIN_DEPTH=x"),
            ("FEATURE=gene", "FEATURE=gene"),
        ] {
            let err = options(&["I=x.bam", "ANNOTATION=a.gtf", bad]).expect_err(bad);
            assert!(matches!(err, Error::Usage(_)), "{err}");
            assert!(err.to_string().contains(text), "{err}");
        }
        let many: Vec<String> = (1..=17).map(|d| format!("MIN_DEPTH={d}")).collect();
        let mut args = strings(&["I=x.bam", "ANNOTATION=a.gtf"]);
        args.extend(many);
        let err = TOOL
            .parse(&args)
            .and_then(|inv| Options::from_invocation(&inv))
            .expect_err("17 depths");
        assert!(err.to_string().contains("at most 16"), "{err}");
    }
}

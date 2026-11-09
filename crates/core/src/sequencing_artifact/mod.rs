//! `CollectSequencingArtifactMetrics`, mirroring Picard 3.5.0
//! (`CollectSequencingArtifactMetrics`, `ArtifactCounter`,
//! `ContextAccumulator`, `SequencingArtifactMetrics`, `ErrorSummaryMetrics`).
//!
//! `OUTPUT` is a prefix. Five files go to `OUTPUT` plus an extension plus
//! `FILE_EXTENSION`: `.pre_adapter_detail_metrics`,
//! `.pre_adapter_summary_metrics`, `.bait_bias_detail_metrics`,
//! `.bait_bias_summary_metrics` and `.error_summary_metrics`.
//!
//! A record counts when it passes vendor QC (unless
//! `INCLUDE_NON_PF_READS`), is not secondary (supplementary records count),
//! is not a duplicate (unless `INCLUDE_DUPLICATES`), is mapped, has a
//! mapping quality of at least `MINIMUM_MAPPING_QUALITY` (255 when absent)
//! and, unless `INCLUDE_UNPAIRED`, is paired with `|TLEN|` from
//! `MINIMUM_INSERT_SIZE` to `MAXIMUM_INSERT_SIZE` (0 is no maximum). Each
//! aligned base of a counted record is then tested in Picard's order: base
//! quality (`OQ` first with `USE_OQ`) at least `MINIMUM_QUALITY_SCORE`,
//! read base not `N`, inside `INTERVALS`, not a `DB_SNP` position, the
//! reference context of `CONTEXT_SIZE` bases on each side inside the contig
//! and all `ACGT`, read base `ACGT`. A base that passes counts once for
//! its library, its reference context and its called base, by read end
//! (second of pair, or not) and strand.
//!
//! Libraries come from the `@RG` header lines (`UnknownLibrary` for a line
//! without `LB`); a record whose library is not among them fails the run.
//! Every library gets rows, in Java `HashMap` order of the names, and every
//! row's `SAMPLE_ALIAS` is the header's samples joined by `,` in `HashSet`
//! order.
//!
//! Per library and per substitution `REF>ALT` (`ALT` not `REF`), the
//! detail rows list every context with `REF` at its centre in sorted order.
//! Pre-adapter counts are pro and con by read orientation: without
//! `TANDEM_READS`, first-of-pair forward and second-of-pair reverse reads
//! over the context, and the opposite pair over its reverse complement, are
//! pro; the other four are con. `ERROR_RATE` is
//! `(PRO_ALT - CON_ALT) / total`, at least 1e-10. Bait-bias counts take
//! every read over the context (`FWD`) and over its reverse complement
//! (`REV`); `ERROR_RATE` is `FWD_ERROR_RATE - REV_ERROR_RATE`, at least
//! 1e-10. `QSCORE` is the Phred score of the error rate rounded to an
//! integer (100 at 1e-10).
//!
//! A summary row takes `TOTAL_QSCORE` from the context-free counts, and
//! the worst full, leading (`XYN`) and trailing (`NYZ`) contexts: the
//! lowest `QSCORE`, ties to the lowest context. With `CONTEXT_SIZE=0` the
//! leading and trailing contexts are the centre base alone.
//! `CONTEXTS_TO_PRINT` filters the detail rows only. The error summary sums
//! the printed pre-adapter rows over every library and context, with a
//! `G` or `T` reference folded to its complement, one row per substitution
//! in sorted order.

mod rows;

use std::collections::HashMap;
use std::path::PathBuf;

use noodles::sam::alignment::record::data::field::{Tag, Value as DataValue};
use noodles::{bam, sam};

use crate::blocks::alignment_blocks;
use crate::error::{Error, Result};
use crate::java::hash_map_order;
use crate::levels::read_groups;
use crate::metric::{MetricsOutput, ToolOutput};
use crate::pass::{self, Accumulator, Collector};
use crate::quality::qualities;
use crate::record::record_name;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool, picard_list};
use crate::vcf::{self, Contig, DbSnpBitSet, Dictionary, Regions};

use rows::{
    ALT_TRANSITIONS, Calls, LibraryRows, R1_NEG, R1_POS, R2_NEG, R2_POS, base_index,
    error_summaries,
};
pub use rows::{
    BaitBiasDetail, BaitBiasSummary, ErrorSummary, PreAdapterDetail, PreAdapterSummary, Report,
    Summary,
};

pub const PRE_ADAPTER_SUMMARY_EXT: &str = ".pre_adapter_summary_metrics";
pub const PRE_ADAPTER_DETAILS_EXT: &str = ".pre_adapter_detail_metrics";
pub const BAIT_BIAS_SUMMARY_EXT: &str = ".bait_bias_summary_metrics";
pub const BAIT_BIAS_DETAILS_EXT: &str = ".bait_bias_detail_metrics";
pub const ERROR_SUMMARY_EXT: &str = ".error_summary_metrics";

const UNKNOWN_LIBRARY: &str = "UnknownLibrary";
const UNKNOWN_SAMPLE: &str = "UnknownSample";
/// The most full contexts counted, 4^11 (`CONTEXT_SIZE` 5): 512 MB of
/// counts per library.
const MAX_CONTEXTS: usize = 1 << 22;

/// `CollectSequencingArtifactMetrics` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "CollectSequencingArtifactMetrics",
    summary: "OUTPUT is a prefix: PREFIX.pre_adapter_detail_metrics, \
              .pre_adapter_summary_metrics, .bait_bias_detail_metrics, \
              .bait_bias_summary_metrics and .error_summary_metrics (plus FILE_EXTENSION); \
              Python returns one row per file: file and rows.",
    options: &[
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix of the five files.")
            .required(),
        Opt::path(
            &["REFERENCE_SEQUENCE", "R"],
            "Reference FASTA with a .fai index.",
        )
        .required(),
        Opt::path(
            &["INTERVALS"],
            "An interval list: only bases inside it count.",
        ),
        Opt::path(
            &["DB_SNP"],
            "A VCF: bases at its records' positions do not count.",
        ),
        Opt::int(
            &["MINIMUM_QUALITY_SCORE", "Q"],
            "Lowest base quality counted (default 20).",
        ),
        Opt::int(
            &["MINIMUM_MAPPING_QUALITY", "MQ"],
            "Lowest mapping quality counted (default 30).",
        ),
        Opt::int(
            &["MINIMUM_INSERT_SIZE", "MIN_INS"],
            "Smallest |TLEN| counted (default 60).",
        ),
        Opt::int(
            &["MAXIMUM_INSERT_SIZE", "MAX_INS"],
            "Largest |TLEN| counted; 0 is no maximum (default 600).",
        ),
        Opt::bool(
            &["INCLUDE_UNPAIRED", "UNPAIRED"],
            "Count unpaired reads, and every pair whatever its insert size (default false).",
        ),
        Opt::bool(
            &["INCLUDE_DUPLICATES", "DUPES"],
            "Count duplicate reads (default false).",
        ),
        Opt::bool(
            &["INCLUDE_NON_PF_READS", "NON_PF"],
            "Count reads that fail vendor quality checks (default false).",
        ),
        Opt::bool(
            &["TANDEM_READS", "TANDEM"],
            "Both ends are read from the same strand (default false).",
        ),
        Opt::bool(
            &["USE_OQ"],
            "Test the OQ qualities when a read has them (default true).",
        ),
        Opt::int(
            &["CONTEXT_SIZE"],
            "Reference bases on each side of the counted base (default 1).",
        ),
        Opt::text(
            &["CONTEXTS_TO_PRINT"],
            "Write only these contexts in the detail files; each adds one, null clears.",
        )
        .repeated(),
        Opt::text(
            &["FILE_EXTENSION", "EXT"],
            "Appended to every file name after its extension.",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

/// The parsed options.
#[derive(Clone, Debug)]
pub struct Options {
    pub prefix: PathBuf,
    pub reference: PathBuf,
    pub intervals: Option<PathBuf>,
    pub db_snp: Option<PathBuf>,
    pub minimum_quality_score: i32,
    pub minimum_mapping_quality: i32,
    pub minimum_insert_size: i32,
    pub maximum_insert_size: i32,
    pub include_unpaired: bool,
    pub include_duplicates: bool,
    pub include_non_pf_reads: bool,
    pub tandem_reads: bool,
    pub use_oq: bool,
    pub context_size: i32,
    pub contexts_to_print: Vec<String>,
    pub file_extension: String,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let prefix = inv
            .output
            .clone()
            .ok_or_else(|| Error::Usage("Option 'OUTPUT' is required.".to_owned()))?;
        let reference = inv.get("REFERENCE_SEQUENCE").ok_or_else(|| {
            Error::Usage("Option 'REFERENCE_SEQUENCE' is required.".to_owned())
        })?;
        Ok(Self {
            prefix,
            reference: reference.into(),
            intervals: inv.get("INTERVALS").map(PathBuf::from),
            db_snp: inv.get("DB_SNP").map(PathBuf::from),
            minimum_quality_score: inv.parse("MINIMUM_QUALITY_SCORE", 20)?,
            minimum_mapping_quality: inv.parse("MINIMUM_MAPPING_QUALITY", 30)?,
            minimum_insert_size: inv.parse("MINIMUM_INSERT_SIZE", 60)?,
            maximum_insert_size: inv.parse("MAXIMUM_INSERT_SIZE", 600)?,
            include_unpaired: inv.flag("INCLUDE_UNPAIRED", false)?,
            include_duplicates: inv.flag("INCLUDE_DUPLICATES", false)?,
            include_non_pf_reads: inv.flag("INCLUDE_NON_PF_READS", false)?,
            tandem_reads: inv.flag("TANDEM_READS", false)?,
            use_oq: inv.flag("USE_OQ", true)?,
            context_size: inv.parse("CONTEXT_SIZE", 1)?,
            contexts_to_print: picard_list(Vec::new(), inv.all("CONTEXTS_TO_PRINT"), |s| {
                Ok(s.to_owned())
            })?,
            file_extension: inv.get("FILE_EXTENSION").unwrap_or("").to_owned(),
        })
    }

    /// `customCommandLineValidation`, which the stand-alone tool runs and
    /// `CollectMultipleMetrics` does not: every message at once.
    fn validate(&self) -> Result<()> {
        let mut messages = Vec::new();
        let full_length = 2 * i64::from(self.context_size) + 1;
        if self.context_size < 0 {
            messages.push("CONTEXT_SIZE cannot be negative".to_owned());
        }
        // `CONTEXTS_TO_PRINT` is a `HashSet`: its messages come in that order.
        let mut contexts: Vec<&str> = Vec::new();
        for context in &self.contexts_to_print {
            if !contexts.contains(&context.as_str()) {
                contexts.push(context);
            }
        }
        for context in hash_map_order(&contexts).into_iter().map(|i| contexts[i]) {
            if context.encode_utf16().count() as i64 != full_length {
                messages.push(format!(
                    "Context {context} is not the length implied by CONTEXT_SIZE: {full_length}"
                ));
            }
        }
        if self.minimum_insert_size < 0 {
            messages.push("MINIMUM_INSERT_SIZE cannot be negative".to_owned());
        }
        if self.maximum_insert_size < 0 {
            messages.push("MAXIMUM_INSERT_SIZE cannot be negative".to_owned());
        }
        if self.maximum_insert_size > 0 && self.maximum_insert_size < self.minimum_insert_size {
            messages.push(
                "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0"
                    .to_owned(),
            );
        }
        if messages.is_empty() {
            Ok(())
        } else {
            Err(Error::Usage(messages.join("\n")))
        }
    }

    /// `OUTPUT` plus `extension` plus `FILE_EXTENSION`.
    fn path(&self, extension: &str) -> PathBuf {
        let mut name = self.prefix.clone().into_os_string();
        name.push(extension);
        name.push(&self.file_extension);
        PathBuf::from(name)
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    opts.validate()?;
    let report = inv.single_pass(|header| Artifacts::new(&opts, header))?;
    Ok(report.files(&opts))
}

fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    let opts = Options::from_invocation(inv)?;
    let artifacts = Artifacts::new(&opts, header)?;
    Ok(pass::boxed(artifacts, move |report: Report| {
        Ok(report.files(&opts))
    }))
}

impl Report {
    /// The five files in the order Picard writes them.
    fn files(&self, opts: &Options) -> ToolOutput {
        ToolOutput {
            main: None,
            extra: vec![
                (
                    opts.path(PRE_ADAPTER_DETAILS_EXT),
                    MetricsOutput::from_rows(&self.pre_adapter_details, None),
                ),
                (
                    opts.path(PRE_ADAPTER_SUMMARY_EXT),
                    MetricsOutput::from_rows(&self.pre_adapter_summaries, None),
                ),
                (
                    opts.path(BAIT_BIAS_DETAILS_EXT),
                    MetricsOutput::from_rows(&self.bait_bias_details, None),
                ),
                (
                    opts.path(BAIT_BIAS_SUMMARY_EXT),
                    MetricsOutput::from_rows(&self.bait_bias_summaries, None),
                ),
                (
                    opts.path(ERROR_SUMMARY_EXT),
                    MetricsOutput::from_rows(&self.error_summaries, None),
                ),
            ],
        }
    }
}

/// One library's counts by full context, indexed by the context read as a
/// base-4 number (`A` 0 to `T` 3, first base most significant), so the
/// index order is the sorted order of the context strings.
struct Library {
    name: String,
    counts: Vec<Calls>,
}

/// The pass accumulator: Picard's `setup`, `acceptRead` and `finish`.
pub struct Artifacts {
    opts: Options,
    context_size: usize,
    sample_alias: String,
    libraries: Vec<Library>,
    /// Library index by `@RG` ID.
    library_of_group: HashMap<Vec<u8>, usize>,
    /// Index of `UnknownLibrary` among the libraries, when a header line
    /// gives it.
    unknown_library: Option<usize>,
    reference: Reference,
    /// The current contig, upper case, and its index.
    contig: Option<(usize, Vec<u8>)>,
    /// `INTERVALS` by the interval list's own contig index, and the mask
    /// of the current contig (bit per 1-based position).
    intervals: Option<Vec<Vec<(i64, i64)>>>,
    mask: Vec<u64>,
    db_snp: Option<DbSnpBitSet>,
    /// Contig names by index, for the dbSNP lookup.
    contig_names: Vec<String>,
    /// A record without a reference index was seen: Picard stops reading
    /// (`usesNoRefReads` is false).
    stopped: bool,
    bases: Vec<u8>,
    quals: Vec<i8>,
}

impl Artifacts {
    pub fn new(opts: &Options, header: &sam::Header) -> Result<Self> {
        if opts.context_size < 0 {
            return Err(Error::Usage("CONTEXT_SIZE cannot be negative".to_owned()));
        }
        let context_size = opts.context_size as usize;
        let contexts = 4usize
            .checked_pow(2 * context_size as u32 + 1)
            .filter(|&n| n <= MAX_CONTEXTS)
            .ok_or_else(|| {
                Error::Usage(format!(
                    "CONTEXT_SIZE={context_size}: too many contexts to count"
                ))
            })?;
        let groups = read_groups(header);
        let mut samples: Vec<String> = Vec::new();
        let mut names: Vec<String> = Vec::new();
        for group in &groups {
            let sample = group
                .sample
                .clone()
                .unwrap_or_else(|| UNKNOWN_SAMPLE.to_owned());
            if !samples.contains(&sample) {
                samples.push(sample);
            }
            let library = group
                .library
                .clone()
                .unwrap_or_else(|| UNKNOWN_LIBRARY.to_owned());
            if !names.contains(&library) {
                names.push(library);
            }
        }
        let sample_keys: Vec<&str> = samples.iter().map(String::as_str).collect();
        let sample_alias = hash_map_order(&sample_keys)
            .into_iter()
            .map(|i| sample_keys[i])
            .collect::<Vec<_>>()
            .join(",");
        let name_keys: Vec<&str> = names.iter().map(String::as_str).collect();
        let libraries: Vec<Library> = hash_map_order(&name_keys)
            .into_iter()
            .map(|i| Library {
                name: names[i].clone(),
                counts: vec![[[0; 4]; 4]; contexts],
            })
            .collect();
        let index_of = |name: &str| libraries.iter().position(|l| l.name == name);
        let library_of_group = groups
            .iter()
            .filter_map(|g| {
                let name = g.library.as_deref().unwrap_or(UNKNOWN_LIBRARY);
                index_of(name).map(|i| (g.id.as_bytes().to_vec(), i))
            })
            .collect();
        let unknown_library = index_of(UNKNOWN_LIBRARY);
        let dictionary: Dictionary = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| Contig {
                name: name.to_string(),
                length: Some(map.length().get() as i64),
            })
            .collect();
        let contig_names = dictionary.iter().map(|c| c.name.clone()).collect();
        let regions = opts
            .intervals
            .as_deref()
            .map(|p| Regions::read(p).map(|r| r.uniqued()))
            .transpose()?;
        let db_snp = match &opts.db_snp {
            Some(path) => {
                if regions.is_some() && !vcf::index_exists(path) {
                    // htsjdk opens the VCF for indexed queries over INTERVALS.
                    return Err(Error::Usage(format!(
                        "An index is required, but none found with file ending .idx, for input source: {}",
                        path.display()
                    )));
                }
                Some(DbSnpBitSet::load_all(
                    path,
                    Some(&dictionary),
                    regions.as_ref(),
                )?)
            }
            None => None,
        };
        let intervals = regions.map(|r| {
            let mut by_contig = vec![Vec::new(); r.dictionary.len()];
            for (contig, start, end) in r.intervals {
                by_contig[contig].push((start, end));
            }
            by_contig
        });
        Ok(Self {
            opts: opts.clone(),
            context_size,
            sample_alias,
            libraries,
            library_of_group,
            unknown_library,
            reference: Reference::open(&opts.reference, header)?,
            contig: None,
            intervals,
            mask: Vec::new(),
            db_snp,
            contig_names,
            stopped: false,
            bases: Vec::new(),
            quals: Vec::new(),
        })
    }

    /// `ReferenceSequenceFileWalker.get` for every record with a reference
    /// index, and the interval mask of that contig.
    fn load_contig(&mut self, id: usize) -> Result<()> {
        if self.contig.as_ref().is_some_and(|(i, _)| *i == id) {
            return Ok(());
        }
        let mut bases = self.reference.contig(id)?.to_vec();
        bases.make_ascii_uppercase();
        self.contig = Some((id, bases));
        if let Some(intervals) = &self.intervals {
            self.mask.clear();
            for &(start, end) in intervals.get(id).map_or(&[][..], Vec::as_slice) {
                for position in start..=end {
                    let (word, bit) = ((position / 64) as usize, position % 64);
                    if word >= self.mask.len() {
                        self.mask.resize(word + 1, 0);
                    }
                    self.mask[word] |= 1 << bit;
                }
            }
        }
        Ok(())
    }

    fn in_mask(&self, position: i64) -> bool {
        let (word, bit) = ((position / 64) as usize, position % 64);
        self.mask.get(word).is_some_and(|w| w & (1 << bit) != 0)
    }

    /// The record-level filters (`AggregateFilter`).
    fn passes(&self, record: &bam::Record) -> bool {
        let flags = record.flags();
        if (!self.opts.include_non_pf_reads && flags.is_qc_fail())
            || flags.is_secondary()
            || (!self.opts.include_duplicates && flags.is_duplicate())
            || flags.is_unmapped()
        {
            return false;
        }
        let mapq = record
            .mapping_quality()
            .map_or(255, |q| i32::from(u8::from(q)));
        if mapq < self.opts.minimum_mapping_quality {
            return false;
        }
        if !self.opts.include_unpaired {
            if !flags.is_segmented() {
                return false;
            }
            let maximum = if self.opts.maximum_insert_size == 0 {
                i32::MAX
            } else {
                self.opts.maximum_insert_size
            };
            let insert = record.template_length().wrapping_abs();
            if insert < self.opts.minimum_insert_size || insert > maximum {
                return false;
            }
        }
        true
    }

    /// The library the record counts for.
    fn library(&self, record: &bam::Record) -> Result<usize> {
        let group = match record.data().get(&Tag::READ_GROUP) {
            Some(Ok(DataValue::String(id))) => self
                .library_of_group
                .get(AsRef::<[u8]>::as_ref(id))
                .copied(),
            _ => None,
        };
        group.or(self.unknown_library).ok_or_else(|| {
            Error::Metric(format!(
                "Record contains library that is missing from header: {UNKNOWN_LIBRARY}"
            ))
        })
    }

    /// The context code of the `2k + 1` reference bases centred on 1-based
    /// `position`, when they lie inside the contig and are all `ACGT`.
    fn context_code(&self, contig: &[u8], position: i64) -> Option<usize> {
        let k = self.context_size as i64;
        let start = position - k - 1;
        let length = 2 * k + 1;
        if start < 0 || start + length > contig.len() as i64 {
            return None;
        }
        let mut code = 0usize;
        for &b in &contig[start as usize..(start + length) as usize] {
            code = code * 4 + base_index(b)?;
        }
        Some(code)
    }
}

impl Accumulator for Artifacts {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let Some(id) = record.reference_sequence_id().transpose()? else {
            self.stopped = true;
            return Ok(());
        };
        self.load_contig(id)?;
        if !self.passes(record) {
            return Ok(());
        }
        let library = self.library(record)?;
        let flags = record.flags();
        let read_two = flags.is_segmented() && flags.is_last_segment();
        let slot = match (read_two, flags.is_reverse_complemented()) {
            (false, false) => R1_POS,
            (false, true) => R1_NEG,
            (true, false) => R2_POS,
            (true, true) => R2_NEG,
        };
        let mut bases = std::mem::take(&mut self.bases);
        bases.clear();
        bases.extend(record.sequence().iter());
        let mut quals = std::mem::take(&mut self.quals);
        qualities(record, self.opts.use_oq, &mut quals)?;
        let malformed = |what| Error::MalformedRecord {
            name: record_name(record),
            what,
        };
        let contig = std::mem::take(&mut self.contig);
        let bases_of_contig = contig.as_ref().map_or(&[][..], |(_, b)| b.as_slice());
        let mut result = Ok(());
        'blocks: for block in alignment_blocks(record)? {
            for offset in 0..block.length {
                let read_pos = block.read_start + offset;
                let ref_pos = block.reference_start + offset as i64;
                let Some(&qual) = quals.get(read_pos - 1) else {
                    result = Err(malformed("fewer base qualities than aligned bases"));
                    break 'blocks;
                };
                if i32::from(qual) < self.opts.minimum_quality_score {
                    continue;
                }
                let Some(&base) = bases.get(read_pos - 1) else {
                    result = Err(malformed("fewer bases than the CIGAR aligns"));
                    break 'blocks;
                };
                let base = base.to_ascii_uppercase();
                if base == b'N' {
                    continue;
                }
                if self.intervals.is_some() && !self.in_mask(ref_pos) {
                    continue;
                }
                if let Some(db_snp) = &self.db_snp
                    && db_snp.is_site(&self.contig_names[id], ref_pos)
                {
                    continue;
                }
                let Some(code) = self.context_code(bases_of_contig, ref_pos) else {
                    continue;
                };
                let Some(called) = base_index(base) else {
                    continue;
                };
                self.libraries[library].counts[code][called][slot] += 1;
            }
        }
        self.contig = contig;
        self.bases = bases;
        self.quals = quals;
        result
    }

    fn done(&self) -> bool {
        self.stopped
    }

    fn finish(self) -> Result<Report> {
        let mut report = Report::default();
        let k = self.context_size;
        for library in &self.libraries {
            let rows = LibraryRows::new(
                &self.sample_alias,
                &library.name,
                k,
                self.opts.tandem_reads,
                &library.counts,
            );
            let print = |context: &str| {
                self.opts.contexts_to_print.is_empty()
                    || self.opts.contexts_to_print.iter().any(|c| c == context)
            };
            for (transition, details) in ALT_TRANSITIONS.iter().zip(&rows.details) {
                let (pre, bait) = rows.summary(*transition);
                report.pre_adapter_summaries.push(PreAdapterSummary(pre));
                report.bait_bias_summaries.push(BaitBiasSummary(bait));
                for (pre, bait) in details {
                    if print(&pre.context) {
                        report.pre_adapter_details.push(pre.clone());
                    }
                    if print(&bait.context) {
                        report.bait_bias_details.push(bait.clone());
                    }
                }
            }
        }
        report.error_summaries = error_summaries(&report.pre_adapter_details);
        Ok(report)
    }
}

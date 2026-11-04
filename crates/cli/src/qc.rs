//! `bammetrics qc` and `bammetrics plan` (docs/05 phase 7): the collectors a
//! preset needs on one read of a BAM, then the verdict and the page.
//!
//! Stages, in order (docs/05 section 3):
//! * prepare: the BAM header, the sample name, the input fingerprint, the
//!   preset, every file an option names, the manifest of an earlier run;
//!   then a `--panel` BED written as `PREFIX.panel.interval_list`
//!   (`BedToIntervalList` with the BAM header as the dictionary).
//! * pass: every one-pass collector on one read through
//!   [`pass::run_isolated`]; a collector that fails fails alone. A tool
//!   with no one-pass collector (`CollectWgsMetrics` walks loci) reads the
//!   BAM on its own.
//! * post: the verdict over the files the run holds.
//! * report: `PREFIX.verdict.json` and `PREFIX.report.html`, the files
//!   `bammetrics report` writes for the same metrics files.
//!
//! Each metrics file is the one the Picard face writes for the same
//! arguments (`bammetrics <Tool> INPUT=... OUTPUT=...`), named by
//! [`multiple_metrics::file_extensions`] and written whole, then renamed.

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::collect;
use bammetrics_core::hs_metrics;
use bammetrics_core::intervals;
use bammetrics_core::levels::{self, Level};
use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::pass;
use bammetrics_core::registry::{self, Invocation, Tool};
use bammetrics_verdict::{Preset, bytes_hash};
use clap::{ArgAction, Args, ValueEnum};

use crate::check::{EXIT_CODES, FailOn};
use crate::error::{Error, Result};
use crate::manifest::{Fingerprint, Key, Manifest, Output, Status, StepRecord};
use crate::report::PRESET_HELP;

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// Progress on stderr (on a terminal), a summary of at most 12 lines on stdout.
    Text,
    /// One JSON object per line on stdout: start, progress, file, error,
    /// check and done events.
    Jsonl,
}

/// `STRAND_SPECIFICITY` of `CollectRnaSeqMetrics`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Strand {
    #[value(name = "NONE")]
    Unstranded,
    #[value(name = "FIRST_READ_TRANSCRIPTION_STRAND")]
    FirstRead,
    #[value(name = "SECOND_READ_TRANSCRIPTION_STRAND")]
    SecondRead,
}

impl Strand {
    fn picard(self) -> &'static str {
        match self {
            Strand::Unstranded => "NONE",
            Strand::FirstRead => "FIRST_READ_TRANSCRIPTION_STRAND",
            Strand::SecondRead => "SECOND_READ_TRANSCRIPTION_STRAND",
        }
    }
}

fn after_help() -> String {
    format!(
        "Files: PREFIX plus the CollectMultipleMetrics extensions (.alignment_summary_metrics, \
.insert_size_metrics, .quality_yield_metrics, .gc_bias.detail_metrics and \
.gc_bias.summary_metrics), .hs_metrics (.BAIT_SET.hs_metrics with several --hs), \
.duplicate_metrics, .wgs_metrics, .rna_metrics; .panel.interval_list when --panel is a BED (its \
@SQ lines from the BAM header); then PREFIX.verdict.json, PREFIX.report.html \
and PREFIX.manifest.json. Each metrics file is byte for byte the one `bammetrics <Tool> \
INPUT=BAM OUTPUT=FILE` writes, apart from the Started on line.\n\n\
Without --fail-on, qc exits 0 on any verdict; bammetrics check always gates.\n\n{EXIT_CODES}\n  \
qc also exits 1 when a collector failed, after it writes the files the others made."
    )
}

#[derive(Args)]
#[command(after_help = after_help())]
pub struct QcArgs {
    /// Input BAM.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,
    /// Output prefix (default: ./ plus the BAM name without .bam). A
    /// directory (an existing one, or a path ending in /) gets the BAM name
    /// as the file part. Missing directories are created.
    #[arg(short = 'o', short_alias = 'O', long, value_name = "PREFIX")]
    pub output: Option<PathBuf>,
    /// Reference FASTA with a .fai index: alignment summary mismatches, hs
    /// dropouts, and the collectors that need one (gc bias, wgs).
    #[arg(short = 'R', long, value_name = "FASTA")]
    pub reference: Option<PathBuf>,
    #[arg(long, default_value = "generic", value_name = "NAME|FILE", help = PRESET_HELP)]
    pub preset: String,
    /// One CollectHsMetrics run: bait and target intervals (interval list
    /// or BED). Repeat for more runs.
    #[arg(long = "hs", num_args = 2, value_names = ["BAITS", "TARGETS"], action = ArgAction::Append)]
    pub hs: Vec<PathBuf>,
    /// Baits and targets in one file: the same as --hs FILE FILE.
    #[arg(long, value_name = "FILE")]
    pub panel: Option<PathBuf>,
    /// Gene annotation (refFlat) for CollectRnaSeqMetrics.
    #[arg(long, value_name = "FILE", requires = "strand")]
    pub ref_flat: Option<PathBuf>,
    /// CollectRnaSeqMetrics STRAND_SPECIFICITY.
    #[arg(long, value_enum, requires = "ref_flat")]
    pub strand: Option<Strand>,
    /// CollectRnaSeqMetrics RIBOSOMAL_INTERVALS.
    #[arg(long, value_name = "FILE", requires = "ref_flat")]
    pub ribosomal_intervals: Option<PathBuf>,
    /// Sample name for the verdict, the page and the manifest (default: the
    /// first @RG SM, else the BAM name). Never written into a metrics file.
    #[arg(long, value_name = "NAME")]
    pub sample: Option<String>,
    /// Exit 10 when the verdict reaches this band, 11 when it is incomplete.
    #[arg(long, value_enum)]
    pub fail_on: Option<FailOn>,
    /// Skip every collector whose manifest record matches (argv, tool
    /// version, preset, and the size, mtime and end hashes of the BAM and of
    /// every FASTA, .fai, interval or annotation file it reads) and whose
    /// files exist. The verdict and the page are always rewritten.
    #[arg(long, conflicts_with = "force")]
    pub resume: bool,
    /// Rerun every step over an earlier run's files.
    #[arg(long)]
    pub force: bool,
    /// BGZF decompression threads (default: available cores, at most 8).
    #[arg(long, default_value_t = pass::default_threads())]
    pub threads: usize,
    #[arg(long, value_enum, default_value_t = Format::Text)]
    pub format: Format,
    /// No progress line on stderr.
    #[arg(short = 'q', long)]
    pub quiet: bool,
    /// Draw every plot family the files have data for.
    #[arg(long)]
    pub all_plots: bool,
    /// The date printed on the page (default: now, as 2026-01-02 03:04:05 UTC).
    #[arg(long, value_name = "ISO")]
    pub date: Option<String>,
}

/// What a collector needs beyond the BAM.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Needs {
    Bam,
    Reference,
    Intervals,
    Annotation,
}

impl Needs {
    fn missing(self) -> &'static str {
        match self {
            Needs::Bam => "",
            Needs::Reference => "needs -R FASTA",
            Needs::Intervals => "needs --hs BAITS TARGETS or --panel FILE",
            Needs::Annotation => "needs --ref-flat FILE and --strand",
        }
    }
}

/// The tool that writes each metrics class a rule can name, in the order
/// the collectors run.
const SOURCES: &[(&str, &str, Needs)] = &[
    (
        "QualityYieldMetrics",
        "CollectQualityYieldMetrics",
        Needs::Bam,
    ),
    (
        "AlignmentSummaryMetrics",
        "CollectAlignmentSummaryMetrics",
        Needs::Bam,
    ),
    ("InsertSizeMetrics", "CollectInsertSizeMetrics", Needs::Bam),
    ("DuplicationMetrics", "CollectDuplicateMetrics", Needs::Bam),
    ("HsMetrics", "CollectHsMetrics", Needs::Intervals),
    (
        "GcBiasSummaryMetrics",
        "CollectGcBiasMetrics",
        Needs::Reference,
    ),
    (
        "GcBiasDetailMetrics",
        "CollectGcBiasMetrics",
        Needs::Reference,
    ),
    ("RnaSeqMetrics", "CollectRnaSeqMetrics", Needs::Annotation),
    ("WgsMetrics", "CollectWgsMetrics", Needs::Reference),
];

/// `picard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics` and
/// `QualityYieldMetrics` both name `QualityYieldMetrics`.
fn short_class(class: &str) -> &str {
    let class = class.rsplit('.').next().unwrap_or(class);
    class.rsplit('$').next().unwrap_or(class)
}

/// How a step reads the BAM.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Read {
    /// On the shared read.
    Pass,
    /// On a read of its own, for the reason given.
    Own(&'static str),
    /// The header only: a prepare step (`--panel` BED to interval list).
    Header,
}

impl Read {
    /// The manifest's `stage` of a step that reads the BAM this way.
    fn stage(self) -> &'static str {
        match self {
            Read::Pass => "pass",
            Read::Own(_) => "own-read",
            Read::Header => "prepare",
        }
    }
}

/// The interval list `qc --panel BED` keeps, `PREFIX.panel.interval_list`
/// (docs/05 example 7).
pub const PANEL_SUFFIX: &str = ".panel.interval_list";

/// One collector of the run.
pub struct Step {
    pub tool: &'static Tool,
    pub inv: Invocation,
    /// The record this run writes before it runs: key and planned outputs.
    pub record: StepRecord,
    /// The rules it serves; empty for an hs or rna run asked for by option.
    pub rules: Vec<String>,
    pub read: Read,
    /// The earlier run's record matches (see [`StepRecord::resumes`]).
    pub resumable: bool,
}

impl Step {
    /// True when this run skips the step.
    pub fn skipped(&self, resume: bool) -> bool {
        resume && self.resumable
    }
}

/// A metrics class some rule reads that no collector of this run writes.
pub struct Unmet {
    pub class: String,
    pub rules: Vec<String>,
    pub why: &'static str,
}

/// Where a run name came from, for the plan.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Origin {
    Option,
    ReadGroup,
    FileName,
    /// `-o` named a directory; the BAM name is the file part.
    Directory,
}

/// Everything prepare decides; nothing has been written.
pub struct Plan {
    pub bam: PathBuf,
    pub prefix: PathBuf,
    pub prefix_origin: Origin,
    pub sample: String,
    pub sample_origin: Origin,
    pub preset: Preset,
    pub fingerprint: Fingerprint,
    /// The prepare stage: steps that read no BAM record, run first.
    pub prepare: Vec<Step>,
    /// The collectors.
    pub steps: Vec<Step>,
    pub unmet: Vec<Unmet>,
    /// The earlier run's manifest, or why it could not be read.
    pub previous: std::result::Result<Option<Manifest>, String>,
}

/// `PREFIX` plus `suffix`.
pub fn with_suffix(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = OsString::from(prefix.as_os_str());
    name.push(suffix);
    PathBuf::from(name)
}

/// The BAM's file name without `.bam`.
fn bam_stem(bam: &Path) -> String {
    let name = bam
        .file_name()
        .map_or_else(|| "sample".to_owned(), |n| n.to_string_lossy().into_owned());
    name.strip_suffix(".bam").unwrap_or(&name).to_owned()
}

/// An input file an option names must exist before anything runs.
fn need_file(path: &Path, option: &str) -> Result<()> {
    if path.is_file() {
        return Ok(());
    }
    Err(Error::Input(format!(
        "{} not found (from {option})",
        path.display()
    )))
}

/// The bait and target files of every hs run: `--hs` pairs, then `--panel`.
fn hs_runs(args: &QcArgs) -> Vec<(PathBuf, PathBuf)> {
    let mut runs: Vec<(PathBuf, PathBuf)> = args
        .hs
        .chunks(2)
        .map(|pair| (pair[0].clone(), pair[1].clone()))
        .collect();
    if let Some(panel) = &args.panel {
        runs.push((panel.clone(), panel.clone()));
    }
    runs
}

fn check_inputs(args: &QcArgs) -> Result<()> {
    need_file(&args.bam, "BAM")?;
    if let Some(reference) = &args.reference {
        need_file(reference, "-R")?;
        need_file(&with_suffix(reference, ".fai"), "-R, its .fai index")?;
    }
    for (baits, targets) in hs_runs(args) {
        need_file(&baits, "--hs or --panel")?;
        need_file(&targets, "--hs or --panel")?;
    }
    if let Some(ref_flat) = &args.ref_flat {
        need_file(ref_flat, "--ref-flat")?;
    }
    if let Some(ribosomal) = &args.ribosomal_intervals {
        need_file(ribosomal, "--ribosomal-intervals")?;
    }
    Ok(())
}

/// The run's parameters every step is built from.
struct Context<'a> {
    args: &'a QcArgs,
    prefix: &'a Path,
    preset: &'a Preset,
    fingerprint: &'a Fingerprint,
}

impl Context<'_> {
    fn available(&self, needs: Needs) -> bool {
        match needs {
            Needs::Bam => true,
            Needs::Reference => self.args.reference.is_some(),
            Needs::Intervals => !hs_runs(self.args).is_empty(),
            Needs::Annotation => self.args.ref_flat.is_some(),
        }
    }

    /// Rule ids reading one of `classes`, and the accumulation levels they
    /// read beyond `ALL_READS`.
    fn rules_for(&self, tool: &str) -> (Vec<String>, Vec<Level>) {
        let mut ids = Vec::new();
        let mut levels = Vec::new();
        for rule in &self.preset.rules {
            let class = short_class(&rule.class);
            if SOURCES.iter().any(|(c, t, _)| *c == class && *t == tool) {
                ids.push(rule.id.clone());
                if rule.level != Level::AllReads && !levels.contains(&rule.level) {
                    levels.push(rule.level);
                }
            }
        }
        (ids, levels)
    }

    /// `INPUT`, the file options, then what the tool reads beyond the BAM.
    fn argv(&self, tool: &Tool, files: &[(String, PathBuf)], extra: &[String]) -> Vec<String> {
        let mut argv = vec![format!("INPUT={}", self.args.bam.display())];
        argv.extend(
            files
                .iter()
                .map(|(option, path)| format!("{option}={}", path.display())),
        );
        let reference = match tool.name {
            "CollectAlignmentSummaryMetrics"
            | "CollectHsMetrics"
            | "CollectGcBiasMetrics"
            | "CollectWgsMetrics" => self.args.reference.as_ref(),
            _ => None,
        };
        if let Some(r) = reference {
            argv.push(format!("REFERENCE_SEQUENCE={}", r.display()));
        }
        argv.extend(extra.iter().cloned());
        argv
    }

    /// The file options of `tool` under the prefix.
    fn files(&self, tool: &str) -> Vec<(String, PathBuf)> {
        file_extensions(tool)
            .unwrap_or_default()
            .iter()
            .map(|(option, ext)| ((*option).to_owned(), with_suffix(self.prefix, ext)))
            .collect()
    }

    fn rna_options(&self) -> Vec<String> {
        let mut extra = Vec::new();
        if let Some(ref_flat) = &self.args.ref_flat {
            extra.push(format!("REF_FLAT={}", ref_flat.display()));
        }
        if let Some(strand) = self.args.strand {
            extra.push(format!("STRAND_SPECIFICITY={}", strand.picard()));
        }
        if let Some(ribosomal) = &self.args.ribosomal_intervals {
            extra.push(format!("RIBOSOMAL_INTERVALS={}", ribosomal.display()));
        }
        extra
    }
}

/// The argv keys that name a file a step reads (the BAM itself is the
/// step's input and is left out).
const READ_FILES: &[&str] = &[
    "INPUT",
    "SEQUENCE_DICTIONARY",
    "REFERENCE_SEQUENCE",
    "BAIT_INTERVALS",
    "TARGET_INTERVALS",
    "REF_FLAT",
    "RIBOSOMAL_INTERVALS",
];

/// The fingerprint of every file `argv` names besides `bam` (a FASTA
/// with its `.fai`), each once, in argv order.
fn read_files(argv: &[String], bam: &Path) -> Result<Vec<(PathBuf, Fingerprint)>> {
    let mut paths: Vec<PathBuf> = Vec::new();
    for arg in argv {
        let Some((key, value)) = arg.split_once('=') else {
            continue;
        };
        let path = PathBuf::from(value);
        if !READ_FILES.contains(&key) || path == bam {
            continue;
        }
        if key == "REFERENCE_SEQUENCE" {
            paths.push(path.clone());
            paths.push(with_suffix(&path, ".fai"));
        } else {
            paths.push(path);
        }
    }
    let mut files: Vec<(PathBuf, Fingerprint)> = Vec::new();
    for path in paths {
        if files.iter().all(|(p, _)| *p != path) {
            let print = Fingerprint::of(&path)?;
            files.push((path, print));
        }
    }
    Ok(files)
}

/// Why a tool with no one-pass collector reads the BAM on its own.
fn own_read_reason(tool: &str) -> &'static str {
    match tool {
        "CollectWgsMetrics" => "it walks reference loci, not one record at a time",
        _ => "it has no one-pass collector",
    }
}

/// What one step runs: its id, tool, files, further arguments and rules.
struct Spec {
    id: String,
    tool: &'static Tool,
    files: Vec<(String, PathBuf)>,
    extra: Vec<String>,
    rules: Vec<String>,
}

impl Context<'_> {
    /// The collector step for `spec`, its record keyed on this run.
    fn step(&self, spec: Spec) -> Result<Step> {
        let Spec {
            id,
            tool,
            files,
            extra,
            rules,
        } = spec;
        let argv = self.argv(tool, &files, &extra);
        let read = if tool.collector.is_some() {
            Read::Pass
        } else {
            Read::Own(own_read_reason(tool.name))
        };
        let outputs = files
            .into_iter()
            .map(|(_, path)| (path, "metrics"))
            .collect();
        self.planned(id, tool, argv, outputs, rules, read)
    }

    /// The step running `tool` on `argv` and writing `outputs` (path and
    /// manifest kind), its record keyed on this run.
    fn planned(
        &self,
        id: String,
        tool: &'static Tool,
        argv: Vec<String>,
        outputs: Vec<(PathBuf, &str)>,
        rules: Vec<String>,
        read: Read,
    ) -> Result<Step> {
        let mut inv = tool
            .parse(&argv)
            .map_err(|e| Error::Usage(format!("{}: {e}", tool.name)))?;
        inv.threads = self.args.threads;
        let record = StepRecord {
            id,
            stage: read.stage().to_owned(),
            tool: tool.name.to_owned(),
            key: Key {
                argv_hash: bytes_hash(argv.join("\0").into_bytes()),
                tool_version: env!("CARGO_PKG_VERSION").to_owned(),
                preset_hash: self.preset.hash.clone(),
                input: self.args.bam.clone(),
                fingerprint: self.fingerprint.clone(),
                files: read_files(&argv, &self.args.bam)?,
            },
            argv,
            outputs: outputs
                .into_iter()
                .map(|(path, kind)| Output {
                    path,
                    tool: tool.name.to_owned(),
                    kind: kind.to_owned(),
                })
                .collect(),
            status: Status::Ok,
            error: None,
        };
        Ok(Step {
            tool,
            inv,
            record,
            rules,
            read,
            resumable: false,
        })
    }
}

/// The hs steps: one per run, each named by its `BAIT_SET` when there are
/// several, as `collect` names them.
fn hs_steps(cx: &Context, rules: &[String], levels: &[String]) -> Result<Vec<Step>> {
    let runs = hs_runs(cx.args);
    let request = collect::Request {
        hs_metrics: runs
            .iter()
            .map(|(baits, targets)| hs_metrics::Options {
                bait_intervals: vec![baits.clone()],
                target_intervals: vec![targets.clone()],
                ..hs_metrics::Options::default()
            })
            .collect(),
        ..collect::Request::default()
    };
    let paths = collect::Outputs::new(cx.prefix, &request)
        .map_err(|e| Error::Usage(e.to_string()))?
        .hs_metrics;
    let tool = registry::find("CollectHsMetrics").expect("registered");
    runs.iter()
        .zip(paths)
        .zip(&request.hs_metrics)
        .map(|(((baits, targets), path), opts)| {
            let mut extra = vec![
                format!("BAIT_INTERVALS={}", baits.display()),
                format!("TARGET_INTERVALS={}", targets.display()),
            ];
            extra.extend(levels.iter().cloned());
            let id = if runs.len() == 1 {
                tool.name.to_owned()
            } else {
                format!("{}.{}", tool.name, hs_metrics::bait_set_name(opts))
            };
            cx.step(Spec {
                id,
                tool,
                files: vec![("OUTPUT".to_owned(), path)],
                extra,
                rules: rules.to_vec(),
            })
        })
        .collect()
}

/// The prepare stage: a `--panel` BED as `PREFIX.panel.interval_list`,
/// written by `BedToIntervalList` with the BAM header as the dictionary,
/// for the Picard tools that read only interval lists. The hs step reads
/// the BED itself, so its file and `BAIT_SET` do not depend on this one.
/// A panel that is an interval list already needs no step.
fn prepare_steps(cx: &Context) -> Result<Vec<Step>> {
    let Some(panel) = &cx.args.panel else {
        return Ok(Vec::new());
    };
    let text = fs::read_to_string(panel)
        .map_err(|e| Error::Input(format!("{}: {e} (from --panel)", panel.display())))?;
    if intervals::is_interval_list(&text) {
        return Ok(Vec::new());
    }
    let tool = registry::find("BedToIntervalList").expect("registered");
    let output = with_suffix(cx.prefix, PANEL_SUFFIX);
    let argv = vec![
        format!("INPUT={}", panel.display()),
        format!("SEQUENCE_DICTIONARY={}", cx.args.bam.display()),
        format!("OUTPUT={}", output.display()),
    ];
    let step = cx.planned(
        tool.name.to_owned(),
        tool,
        argv,
        vec![(output, "intervals")],
        Vec::new(),
        Read::Header,
    )?;
    Ok(vec![step])
}

/// Every collector the preset's rules need (plus hs for `--hs`/`--panel`
/// and rna for `--ref-flat`), and the classes no collector can write.
fn collectors(cx: &Context) -> Result<(Vec<Step>, Vec<Unmet>)> {
    let mut steps = Vec::new();
    let mut tools: Vec<&str> = Vec::new();
    for (_, tool, needs) in SOURCES {
        if tools.contains(tool) {
            continue;
        }
        let (rules, levels) = cx.rules_for(tool);
        let asked = match *tool {
            "CollectHsMetrics" => cx.available(Needs::Intervals),
            "CollectRnaSeqMetrics" => cx.available(Needs::Annotation),
            _ => false,
        };
        if (rules.is_empty() && !asked) || !cx.available(*needs) {
            continue;
        }
        tools.push(tool);
        let registered = registry::find(tool).expect("registered");
        let takes_levels = registered
            .options
            .iter()
            .any(|o| o.names[0] == "METRIC_ACCUMULATION_LEVEL");
        let levels: Vec<String> = if takes_levels {
            levels
                .iter()
                .map(|l| format!("METRIC_ACCUMULATION_LEVEL={}", l.name()))
                .collect()
        } else {
            Vec::new()
        };
        if *tool == "CollectHsMetrics" {
            steps.extend(hs_steps(cx, &rules, &levels)?);
            continue;
        }
        let mut extra = levels;
        if *tool == "CollectRnaSeqMetrics" {
            extra.extend(cx.rna_options());
        }
        steps.push(cx.step(Spec {
            id: (*tool).to_owned(),
            tool: registered,
            files: cx.files(tool),
            extra,
            rules,
        })?);
    }
    Ok((steps, unmet(cx, &tools)))
}

/// The classes some rule reads that no step of this run writes, with why.
fn unmet(cx: &Context, tools: &[&str]) -> Vec<Unmet> {
    let mut unmet: Vec<Unmet> = Vec::new();
    for rule in &cx.preset.rules {
        let class = short_class(&rule.class);
        let source = SOURCES.iter().find(|(c, _, _)| *c == class);
        if source.is_some_and(|(_, tool, _)| tools.contains(tool)) {
            continue;
        }
        let why = source.map_or("no collector writes this class", |(_, _, needs)| {
            needs.missing()
        });
        match unmet.iter_mut().find(|u| u.class == class) {
            Some(u) => u.rules.push(rule.id.clone()),
            None => unmet.push(Unmet {
                class: class.to_owned(),
                rules: vec![rule.id.clone()],
                why,
            }),
        }
    }
    unmet
}

/// The sample: `--sample`, else the first `@RG` with an `SM`, else the
/// BAM name. The header is read either way, so a file that is not a BAM
/// stops the run here.
fn sample(args: &QcArgs) -> Result<(String, Origin)> {
    let header =
        pass::read_header(&args.bam).map_err(|e| Error::Input(format!("{e} (from BAM)")))?;
    if let Some(name) = &args.sample {
        return Ok((name.clone(), Origin::Option));
    }
    Ok(levels::read_groups(&header)
        .into_iter()
        .find_map(|rg| rg.sample)
        .map_or_else(
            || (bam_stem(&args.bam), Origin::FileName),
            |sm| (sm, Origin::ReadGroup),
        ))
}

/// `1 collector`, `4 collectors`.
pub fn plural(n: usize, word: &str) -> String {
    if n == 1 {
        format!("1 {word}")
    } else {
        format!("{n} {word}s")
    }
}

/// True for an `-o` that ends in a path separator or is a directory: a
/// prefix there would make hidden files (`out/.manifest.json`) that
/// `report` skips.
fn names_directory(output: &Path) -> bool {
    let text = output.as_os_str().to_string_lossy();
    text.ends_with(std::path::is_separator) || output.is_dir()
}

/// Prepare: everything checked and decided, nothing written.
pub fn prepare(args: &QcArgs) -> Result<Plan> {
    check_inputs(args)?;
    let preset = Preset::resolve(&args.preset)?;
    let (sample, sample_origin) = sample(args)?;
    let (prefix, prefix_origin) = match &args.output {
        Some(dir) if names_directory(dir) => (dir.join(bam_stem(&args.bam)), Origin::Directory),
        Some(prefix) => (prefix.clone(), Origin::Option),
        None => (PathBuf::from(bam_stem(&args.bam)), Origin::FileName),
    };
    if prefix.file_name().is_none() {
        return Err(Error::Usage(format!(
            "-o {}: a prefix needs a file name part, as in out/NA12878",
            prefix.display()
        )));
    }
    let fingerprint = Fingerprint::of(&args.bam)?;
    let cx = Context {
        args,
        prefix: &prefix,
        preset: &preset,
        fingerprint: &fingerprint,
    };
    let mut prepare = prepare_steps(&cx)?;
    let (mut steps, unmet) = collectors(&cx)?;
    let previous = Manifest::read(&manifest_path(&prefix));
    if let Ok(Some(manifest)) = &previous {
        for step in prepare.iter_mut().chain(&mut steps) {
            step.resumable = manifest
                .step(&step.record.id)
                .is_some_and(|old| step.record.resumes(old));
        }
    }
    Ok(Plan {
        bam: args.bam.clone(),
        prefix,
        prefix_origin,
        sample,
        sample_origin,
        preset,
        fingerprint,
        prepare,
        steps,
        unmet,
        previous,
    })
}

pub fn manifest_path(prefix: &Path) -> PathBuf {
    with_suffix(prefix, ".manifest.json")
}

pub fn verdict_path(prefix: &Path) -> PathBuf {
    with_suffix(prefix, ".verdict.json")
}

pub fn page_path(prefix: &Path) -> PathBuf {
    with_suffix(prefix, ".report.html")
}

/// Create the prefix's directory.
pub fn make_prefix_dir(prefix: &Path) -> Result<()> {
    match prefix.parent().filter(|d| !d.as_os_str().is_empty()) {
        Some(dir) => fs::create_dir_all(dir).map_err(|e| Error::io(dir, e)),
        None => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_names_shorten_to_the_rule_form() {
        assert_eq!(
            short_class("picard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics"),
            "QualityYieldMetrics"
        );
        assert_eq!(
            short_class("picard.analysis.directed.HsMetrics"),
            "HsMetrics"
        );
        assert_eq!(short_class("InsertSizeMetrics"), "InsertSizeMetrics");
    }

    #[test]
    fn every_source_tool_is_registered_and_named_by_the_table() {
        for (_, tool, _) in SOURCES {
            assert!(registry::find(tool).is_some(), "{tool}");
            assert!(file_extensions(tool).is_some(), "{tool}");
        }
    }

    #[test]
    fn a_directory_output_is_not_a_prefix() {
        assert!(names_directory(Path::new("out/")));
        assert!(names_directory(Path::new(".")));
        assert!(names_directory(&std::env::temp_dir()));
        assert!(!names_directory(Path::new("no-such-dir/s")));
    }

    #[test]
    fn the_default_prefix_is_the_bam_name() {
        assert_eq!(bam_stem(Path::new("runs/NA12878.bam")), "NA12878");
        assert_eq!(bam_stem(Path::new("x.cram")), "x.cram");
        assert_eq!(
            with_suffix(Path::new("out/s"), ".manifest.json"),
            PathBuf::from("out/s.manifest.json")
        );
    }
}

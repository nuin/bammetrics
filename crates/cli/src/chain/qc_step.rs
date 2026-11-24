//! `qc` as the last step of a chain (docs/chaining.md section 29).
//!
//! The step evaluates a preset over the metrics files of the taps that
//! read its stream and writes `PREFIX.verdict.json`, `PREFIX.report.html`
//! (and `PREFIX.report.pdf`, `--also markdown,junit,multiqc`) as the `qc`
//! verb writes them for the same taps on the same records. The collectors
//! the preset needs and no tap of that stream runs are added to the chain
//! as taps, numbered after the qc step, each built as `qc` builds it. A tool
//! with no one-pass collector (`CollectWgsMetrics`) reads a file of the
//! stream on its own after the chain: the source when no transform is in
//! front, else the stream the last transform kept (`--keep`), else the BAM
//! a `write` step writes of that stream (`-o out.bam`). With no such file
//! it is a materialized step of the chain's read (docs/chaining.md
//! section 31): it runs on a copy of the stream in the temporary directory,
//! as `wgs` given as a step does. A tap or a materialized tool the chain
//! already runs on the stream is not added twice, and the verdict reads its
//! file. The step's exit code is the chain's when every step finished: 10
//! or 11 under `--fail-on`, as `check` gives them.

use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::CivilTime;
use bammetrics_core::registry::Ty;
use bammetrics_report::Options;
use bammetrics_verdict::{Input, Verdict, evaluate};
use clap::ValueEnum;

use super::keep::Kept;
use super::prepare::{Planned, PlannedFile, Role, materialize, planned_files};
use super::run::{StepOutcome, StepStatus, Written};
use super::{QC, VERB, usage};
use crate::also::Also;
use crate::check::{self, FailOn};
use crate::error::VerbError;
use crate::format::Format;
use crate::grammar::{Flag, Key, Parsed};
use crate::pdf::{PageSize, Pdf, PdfArgs};
use crate::qc::{self, QcArgs, Read, ReportForm, Sex, Strand, Unmet, with_suffix};
use crate::report::write_atomic;
use crate::verbs::{self, option_lines};

const PRESET: Key = Key::Own("preset");
const REPORT: Key = Key::Own("report");
const ALSO: Key = Key::Own("also");
const FAIL_ON: Key = Key::Own("fail-on");
const HS: Key = Key::Own("hs");
const PANEL: Key = Key::Own("panel");
const REF_FLAT: Key = Key::Own("ref-flat");
const STRAND: Key = Key::Own("strand");
const RIBOSOMAL: Key = Key::Own("ribosomal-intervals");
const SEX: Key = Key::Own("sex");
const ALL_PLOTS: Key = Key::Own("all-plots");
const DATE: Key = Key::Own("date");
const PDF_FONT: Key = Key::Own("pdf-font");
const PAGE_SIZE: Key = Key::Own("page-size");

/// The options of the qc step: those of the `qc` verb that are not the
/// chain's (the BAM, `-o`, `-R`, `--threads`, `--sample`, `--resume`).
pub fn flags() -> Vec<Flag> {
    let mut hs = Flag::new(
        HS,
        "hs",
        &[],
        Ty::Path,
        "One CollectHsMetrics run: bait and target intervals (interval list or BED).",
    )
    .placeholder("BAITS TARGETS")
    .repeated();
    hs.values = 2;
    vec![
        Flag::new(
            PRESET,
            "preset",
            &[],
            Ty::Text,
            "The rules: generic (default), wgs, clinical-panel, rna-seq, or a preset file.",
        )
        .placeholder("NAME|FILE"),
        Flag::new(
            REPORT,
            "report",
            &[],
            Ty::Text,
            "The sample report: html (PREFIX.report.html, the default), pdf \
             (PREFIX.report.pdf, cargo feature pdf), or html,pdf.",
        )
        .placeholder("FORM[,FORM]")
        .repeated(),
        Flag::new(
            ALSO,
            "also",
            &[],
            Ty::Text,
            "The forms of the verdict: markdown (PREFIX.report.md), junit (PREFIX.junit.xml), \
             multiqc (PREFIX.verdict_mqc.json). The per-file formats are the chain's --also.",
        )
        .placeholder("FORMAT[,FORMAT]")
        .repeated(),
        Flag::new(
            FAIL_ON,
            "fail-on",
            &[],
            Ty::Text,
            "warn or fail: the chain exits 10 when the verdict reaches it, 11 when it is \
             incomplete. Without it the verdict does not change the exit code.",
        )
        .placeholder("BAND"),
        hs,
        Flag::new(
            PANEL,
            "panel",
            &[],
            Ty::Path,
            "Baits and targets in one file: the same as --hs FILE FILE.",
        )
        .placeholder("FILE"),
        Flag::new(
            REF_FLAT,
            "ref-flat",
            &[],
            Ty::Path,
            "Gene annotation (refFlat) for CollectRnaSeqMetrics; needs --strand.",
        )
        .placeholder("FILE"),
        Flag::new(
            STRAND,
            "strand",
            &[],
            Ty::Text,
            "CollectRnaSeqMetrics STRAND_SPECIFICITY: NONE, FIRST_READ_TRANSCRIPTION_STRAND or \
             SECOND_READ_TRANSCRIPTION_STRAND.",
        )
        .placeholder("STRAND"),
        Flag::new(
            RIBOSOMAL,
            "ribosomal-intervals",
            &[],
            Ty::Path,
            "CollectRnaSeqMetrics RIBOSOMAL_INTERVALS; needs --ref-flat.",
        )
        .placeholder("FILE"),
        Flag::new(
            SEX,
            "sex",
            &[],
            Ty::Text,
            "male or female: CollectSexCheckMetrics EXPECTED_SEX, for a preset with a \
             SexCheckMetrics rule.",
        )
        .placeholder("SEX"),
        Flag::new(
            ALL_PLOTS,
            "all-plots",
            &[],
            Ty::Bool,
            "Draw every plot family the files have data for.",
        ),
        Flag::new(
            DATE,
            "date",
            &[],
            Ty::Text,
            "The date printed on the page (default: now, as 2026-01-02 03:04:05 UTC).",
        )
        .placeholder("ISO"),
        Flag::new(
            PDF_FONT,
            "pdf-font",
            &[],
            Ty::Path,
            "A TrueType or OpenType font of at most 32 MiB for the text of the PDF.",
        )
        .placeholder("FONT"),
        Flag::new(
            PAGE_SIZE,
            "page-size",
            &[],
            Ty::Text,
            "The page size of the PDF: a4 (default) or letter.",
        )
        .placeholder("SIZE"),
        Flag::new(Key::Help, "help", &["h"], Ty::Bool, "This help."),
    ]
}

/// `bammetrics chain ... + qc --help`.
pub fn help() -> String {
    format!(
        "bammetrics chain ... + qc [OPTIONS]: the last step. It evaluates the preset over the \
         metrics files of the taps on its stream and writes PREFIX.verdict.json and \
         PREFIX.report.html as bammetrics qc writes them for the same taps on the same \
         records. A collector the preset needs that no tap of the stream runs is added as a \
         tap. A tool with no collector (CollectWgsMetrics) reads a file of the stream after \
         the chain: the source with no transform in front, the last kept stream (--keep) or \
         the BAM a write step wrote (-o out.bam); with no such file it runs on a copy of the \
         stream, as the wgs step does. The chain writes PREFIX.manifest.json (schema \
         bammetrics.manifest/2).\n\n\
         Every option also reads as KEY=VALUE. The source, -o, -R, --threads, --sample and \
         --resume belong to the chain: give them before the first step, or after qc in the \
         form without +. qc may follow a recipe: bammetrics prep BAM -o out.bam qc.\n\n\
         {}\n\nOptions:\n{}",
        check::EXIT_CODES,
        option_lines(&flags())
    )
}

/// What the qc step is planned against.
pub(super) struct Context<'a> {
    pub source_path: &'a Path,
    /// The source is a regular file: a tool that reads on its own and the
    /// dictionary of a `--panel` BED can read it.
    pub replayable: bool,
    pub reference: Option<&'a str>,
    pub prefix: &'a Path,
    /// The formats of every metrics file, `picard` first.
    pub formats: &'a [Format],
    pub threads: usize,
    pub sample: Option<&'a str>,
    /// The chain's source as given, and its `--tmp-dir`: what a
    /// materialized step is checked against before the first read.
    pub source: &'a str,
    pub tmp_dir: Option<&'a Path>,
}

/// Where the qc step stands: its number and label, the stream it reads,
/// its segment and the transforms in front (the `CHAIN=` of a tap).
pub(super) struct At<'a> {
    pub number: usize,
    pub label: &'a str,
    pub reads: &'a str,
    pub segment: usize,
    pub behind: &'a str,
}

/// The qc step after prepare.
pub struct QcPlan {
    /// Its index in `Chain::steps`.
    pub index: usize,
    pub args: QcArgs,
    pub preset: bammetrics_verdict::Preset,
    pub pdf: Option<Pdf>,
    /// The taps it added, by index in `Chain::steps`.
    pub added: Vec<usize>,
    /// Those of them that are materialized steps: a tool with no collector
    /// on a stream no file holds.
    pub materialized: Vec<usize>,
    /// The steps that read a file of the stream on their own after the
    /// chain (`CollectWgsMetrics`).
    pub own: Vec<usize>,
    /// The `--panel` BED as `PREFIX.panel.interval_list`, after the chain.
    pub prepare: Vec<usize>,
    /// The classes some rule reads that no step writes, with why.
    pub unmet: Vec<Unmet>,
    /// Its options as `--KEY=VALUE`, the canonical argv of the manifest.
    pub argv: Vec<String>,
    /// The files its options name, for the manifest.
    pub inputs: Vec<PathBuf>,
}

impl QcPlan {
    /// True for a step that runs after the chain's reads (an own read, the
    /// panel's interval list).
    pub fn after_reads(&self, index: usize) -> bool {
        self.own.contains(&index) || self.prepare.contains(&index)
    }
}

/// The values of one option, each comma list split.
fn listed(parsed: &Parsed, key: Key) -> impl Iterator<Item = (&str, &str)> {
    parsed.all(key).flat_map(|given| {
        given.values.iter().flat_map(move |value| {
            value
                .split(',')
                .map(str::trim)
                .map(move |name| (name, given.spelled.as_str()))
        })
    })
}

fn last(parsed: &Parsed, key: Key) -> Option<(&str, &str)> {
    parsed
        .last(key)
        .map(|given| (given.values[0].as_str(), given.spelled.as_str()))
}

/// `value` as a value of the clap enum `T`, or the usage error naming the
/// values it takes.
fn choice<T: ValueEnum>(place: &str, value: &str, spelled: &str) -> Result<T, VerbError> {
    T::from_str(value, true).map_err(|_| {
        let names: Vec<String> = T::value_variants()
            .iter()
            .filter_map(|v| v.to_possible_value().map(|p| p.get_name().to_owned()))
            .collect();
        usage(format!(
            "{place}: {spelled} {value:?}: expected {}",
            names.join(", ")
        ))
        .from(spelled)
    })
}

/// The `qc` arguments of the step: the chain's source, prefix, reference,
/// threads and sample, then the step's own options.
fn qc_args(cx: &Context, parsed: &Parsed, place: &str) -> Result<QcArgs, VerbError> {
    let mut report = Vec::new();
    for (name, spelled) in listed(parsed, REPORT) {
        let form: ReportForm = choice(place, name, spelled)?;
        if !report.contains(&form) {
            report.push(form);
        }
    }
    if report.is_empty() {
        report.push(ReportForm::Html);
    }
    let mut also = Vec::new();
    for (name, spelled) in listed(parsed, ALSO) {
        let form: Also = choice(place, name, spelled)?;
        if form.per_file().is_some() {
            return Err(usage(format!(
                "{place}: --also {name} is a format of every metrics file and belongs to the \
                 chain"
            ))
            .from(spelled)
            .hint(format!(
                "give --also {name} before the first +; qc takes markdown, junit and multiqc"
            )));
        }
        if !also.contains(&form) {
            also.push(form);
        }
    }
    let path = |key: Key| last(parsed, key).map(|(value, _)| PathBuf::from(value));
    let fail_on = last(parsed, FAIL_ON)
        .map(|(value, spelled)| choice::<FailOn>(place, value, spelled))
        .transpose()?;
    let strand = last(parsed, STRAND)
        .map(|(value, spelled)| choice::<Strand>(place, value, spelled))
        .transpose()?;
    let sex = last(parsed, SEX)
        .map(|(value, spelled)| choice::<Sex>(place, value, spelled))
        .transpose()?;
    let page_size = last(parsed, PAGE_SIZE)
        .map(|(value, spelled)| choice::<PageSize>(place, value, spelled))
        .transpose()?;
    let ref_flat = path(REF_FLAT);
    let ribosomal = path(RIBOSOMAL);
    // The verb's clap rules: an annotation needs its strand, and both go
    // with --ref-flat.
    let needs = [
        (
            ref_flat.is_some() && strand.is_none(),
            "--ref-flat",
            "--strand",
        ),
        (
            strand.is_some() && ref_flat.is_none(),
            "--strand",
            "--ref-flat",
        ),
        (
            ribosomal.is_some() && ref_flat.is_none(),
            "--ribosomal-intervals",
            "--ref-flat",
        ),
    ];
    if let Some((_, given, missing)) = needs.iter().find(|(wanted, _, _)| *wanted) {
        return Err(usage(format!("{place}: {given} needs {missing}")).from(*given));
    }
    Ok(QcArgs {
        bam: cx.source_path.to_path_buf(),
        output: Some(cx.prefix.to_path_buf()),
        reference: cx.reference.map(PathBuf::from),
        preset: last(parsed, PRESET).map_or_else(|| "generic".to_owned(), |(v, _)| v.to_owned()),
        hs: parsed
            .all(HS)
            .flat_map(|given| given.values.iter().map(PathBuf::from))
            .collect(),
        panel: path(PANEL),
        ref_flat,
        strand,
        ribosomal_intervals: ribosomal,
        sex,
        sample: cx.sample.map(str::to_owned),
        fail_on,
        resume: false,
        force: false,
        threads: cx.threads,
        validation_stringency: None,
        format: qc::Format::Text,
        also,
        quiet: true,
        all_plots: parsed.flag(ALL_PLOTS),
        date: last(parsed, DATE).map(|(value, _)| value.to_owned()),
        report,
        pdf: PdfArgs {
            pdf_font: path(PDF_FONT),
            page_size,
        },
        write_page: true,
    })
}

/// The step's options as `--KEY=VALUE`, in the order of [`flags`].
fn canonical(parsed: &Parsed) -> Vec<String> {
    let mut argv = Vec::new();
    for flag in flags() {
        for given in parsed.all(flag.key) {
            argv.push(format!("--{}={}", flag.long, given.values.join(" ")));
        }
    }
    argv
}

/// A qc error as an error of the step at `place`.
fn step_error(place: &str, error: crate::error::Error) -> VerbError {
    match error {
        crate::error::Error::Verb(mut e) => {
            e.what = format!("{place}: {}", e.what);
            VERB.clone_into(&mut e.verb);
            e
        }
        crate::error::Error::Usage(what) => usage(format!("{place}: {what}")),
        other => VerbError::input(VERB, format!("{place}: {other}")),
    }
}

/// The steps the qc step adds after itself, from the verb's plan: a tap per
/// collector no tap of `at.reads` runs, an own read per tool without a
/// collector, the panel's interval list.
struct Added {
    steps: Vec<Planned>,
    own: Vec<usize>,
    /// The tools that run on a copy of the stream in the chain's read.
    materialized: Vec<usize>,
    prepare: Vec<usize>,
    unmet: Vec<Unmet>,
}

/// The file a tool that reads on its own reads for the qc step: the source
/// with no transform in front, else the stream the last transform keeps
/// (`--keep`), else the BAM a `write` step writes of that stream (the BAM
/// `-o out.bam` names). `None` when no file holds the stream: the tool is
/// then a materialized step of the chain's read ([`materialized`]).
fn own_read_file(cx: &Context, at: &At, before: &[Planned], kept: &[Kept]) -> Option<PathBuf> {
    let last = before.iter().rposition(Planned::hands_on);
    match last {
        None => cx.replayable.then(|| cx.source_path.to_path_buf()),
        Some(index) => kept
            .iter()
            .find(|k| k.index == index)
            .map(|k| k.path().to_path_buf())
            .or_else(|| written_bam(at, before)),
    }
}

/// The BAM a `write` step writes of the qc step's stream, when one does:
/// the records the taps of that stream read, in a file.
fn written_bam(at: &At, before: &[Planned]) -> Option<PathBuf> {
    before
        .iter()
        .filter(|step| step.reads == at.reads && !step.stdout)
        .filter_map(|step| step.write.as_ref())
        .find(|options| options.format == bammetrics_core::stage::WriteFormat::Bam)
        .map(|options| options.path.clone())
}

/// A step of the qc verb's plan as a step of the chain, numbered `number`.
fn planned(
    cx: &Context,
    at: &At,
    step: qc::Step,
    number: usize,
    input: Option<&Path>,
) -> Result<Planned, VerbError> {
    let tool = step.tool;
    let verb = verbs::for_class(tool.name).map_or(tool.name, |v| v.name);
    let label = format!("step {number} ({verb})");
    let mut argv = step.record.argv;
    let mut inv = step.inv;
    if let Some(input) = input {
        argv[0] = format!("INPUT={}", input.display());
        inv = tool
            .parse(&argv)
            .map_err(|e| usage(format!("{label}: {e}")))?;
        // The kept BAM is a file the chain writes, not one it reads before
        // the run: no claim may see it as an input.
        inv.inputs.retain(|path| path != input);
    }
    let (formats, files) = if step.read == Read::Header {
        let files = step
            .record
            .outputs
            .iter()
            .map(|o| PlannedFile {
                path: o.path.clone(),
                kind: "intervals",
            })
            .collect();
        (vec![Format::Picard], files)
    } else {
        (cx.formats.to_vec(), planned_files(tool, &argv, cx.formats))
    };
    let mut command = inv.command;
    if step.read == Read::Pass && !at.behind.is_empty() {
        command.push_str(&format!(" CHAIN={}", at.behind));
    }
    Ok(Planned {
        number,
        verb,
        tool: Some(tool),
        role: Role::Tap,
        label,
        reads: at.reads.to_owned(),
        leaves: None,
        segment: at.segment,
        files,
        late_prefix: None,
        stdout: false,
        command: Some(command),
        argv,
        output: inv.output,
        notes: inv.notes,
        formats,
        write: None,
        vcf: None,
        inputs: inv.inputs,
    })
}

/// A tool that reads on its own as a materialized step of the chain's
/// read, for a stream no file holds (a transform in front without `--keep`,
/// or a piped source): the step `VERB` given on the command line would be,
/// with the checks it gets before the first read.
fn materialized(
    cx: &Context,
    at: &At,
    step: qc::Step,
    number: usize,
) -> Result<Planned, VerbError> {
    let mut planned = planned(cx, at, step, number, None)?;
    planned.role = Role::Materialize;
    if let Some(tool) = planned.tool {
        let before = materialize::Before {
            source: cx.source,
            reads_source: false,
            has_reference: cx.reference.is_some(),
            tmp_dir: cx.tmp_dir,
        };
        materialize::check(tool, &planned.argv, &before)
            .map_err(|e| usage(format!("{}: {e}", planned.label)))?;
    }
    if let (false, Some(command)) = (at.behind.is_empty(), planned.command.as_mut()) {
        command.push_str(&format!(" CHAIN={}", at.behind));
    }
    Ok(planned)
}

/// The taps, own reads and prepare step of the verb's plan that the chain
/// does not run already.
fn added(
    cx: &Context,
    at: &At,
    plan: qc::ChainQc,
    before: &[Planned],
    kept: &[Kept],
) -> Result<(Added, bammetrics_verdict::Preset, Option<Pdf>), VerbError> {
    let qc::ChainQc {
        preset,
        prepare,
        steps,
        mut unmet,
        pdf,
    } = plan;
    let mut added = Added {
        steps: Vec::new(),
        own: Vec::new(),
        materialized: Vec::new(),
        prepare: Vec::new(),
        unmet: Vec::new(),
    };
    // A tap of the stream, or a materialized tool on it (`wgs` as a step).
    let covered = |name: &str| {
        before.iter().any(|s| {
            s.reads_only() && s.reads == at.reads && s.tool.is_some_and(|t| t.name == name)
        })
    };
    // The qc step is number `before.len() + 1`; its steps follow it.
    let first = before.len() + 2;
    for step in prepare {
        if !cx.replayable {
            return Err(usage(format!(
                "{}: --panel names a BED, and its interval list takes the dictionary of the \
                 source, which a pipe cannot give twice",
                at.label
            ))
            .hint("give the panel as an interval list, or the source as a file"));
        }
        let number = first + added.steps.len();
        added.prepare.push(number - 1);
        added.steps.push(planned(cx, at, step, number, None)?);
    }
    for step in steps {
        if covered(step.tool.name) {
            continue;
        }
        let number = first + added.steps.len();
        match step.read {
            Read::Own(_) => match own_read_file(cx, at, before, kept) {
                Some(file) => {
                    let input = (file != cx.source_path).then_some(file.as_path());
                    added.own.push(number - 1);
                    added.steps.push(planned(cx, at, step, number, input)?);
                }
                None if bammetrics_core::stage::materialize::materializes(step.tool) => {
                    added.materialized.push(number - 1);
                    added.steps.push(materialized(cx, at, step, number)?);
                }
                None => {
                    // The class stays unmeasured, as for the verb without
                    // the file its tool needs.
                    for rule in &step.rules {
                        let class = preset
                            .rules
                            .iter()
                            .find(|r| &r.id == rule)
                            .map(|r| r.class.rsplit(['.', '$']).next().unwrap_or(&r.class));
                        if let Some(class) = class {
                            match unmet.iter_mut().find(|u| u.class == class) {
                                Some(u) => u.rules.push(rule.clone()),
                                None => unmet.push(Unmet {
                                    class: class.to_owned(),
                                    rules: vec![rule.clone()],
                                    why: "reads a file of the stream on its own: --keep on the \
                                          last transform, or a write step",
                                }),
                            }
                        }
                    }
                }
            },
            Read::Pass | Read::Header => {
                added.steps.push(planned(cx, at, step, number, None)?);
            }
        }
    }
    // A class a tap of the stream writes is measured: its tool's file is
    // among the metrics files the verdict reads.
    unmet.retain(|u| !qc::source_tool(&u.class).is_some_and(covered));
    added.unmet = unmet;
    Ok((added, preset, pdf))
}

/// The files the qc step writes under the prefix.
fn qc_files(prefix: &Path, args: &QcArgs) -> Vec<PlannedFile> {
    let mut files = vec![PlannedFile {
        path: qc::verdict_path(prefix),
        kind: "verdict",
    }];
    if args.writes(ReportForm::Html) {
        files.push(PlannedFile {
            path: qc::page_path(prefix),
            kind: "report",
        });
    }
    if args.writes(ReportForm::Pdf) {
        files.push(PlannedFile {
            path: qc::pdf_path(prefix),
            kind: "pdf",
        });
    }
    for also in &args.also {
        files.push(PlannedFile {
            path: with_suffix(prefix, also.run_suffix()),
            kind: also.name(),
        });
    }
    files
}

/// Plan the qc step at `at`, after the steps `before` (the qc step is the
/// last one asked for): the step itself, then the steps it adds.
pub(super) fn plan(
    cx: &Context,
    at: &At,
    parsed: &Parsed,
    before: &[Planned],
    kept: &[Kept],
) -> Result<(Vec<Planned>, QcPlan), VerbError> {
    let args = qc_args(cx, parsed, at.label)?;
    let verb_plan = qc::plan_chain_step(&args, cx.prefix).map_err(|e| step_error(at.label, e))?;
    let (added, preset, pdf) = added(cx, at, verb_plan, before, kept)?;
    let mut inputs: Vec<PathBuf> = args.hs.clone();
    inputs.extend(
        [
            &args.panel,
            &args.ref_flat,
            &args.ribosomal_intervals,
            &args.pdf.pdf_font,
        ]
        .into_iter()
        .flatten()
        .cloned(),
    );
    if Path::new(&args.preset).is_file() {
        inputs.push(PathBuf::from(&args.preset));
    }
    let notes = added
        .unmet
        .iter()
        .map(|u| format!("{} not measured: {}", u.class, u.why))
        .collect();
    let qc_step = Planned {
        number: at.number,
        verb: QC,
        tool: None,
        role: Role::Qc,
        label: at.label.to_owned(),
        reads: at.reads.to_owned(),
        leaves: None,
        segment: at.segment,
        files: qc_files(cx.prefix, &args),
        late_prefix: None,
        stdout: false,
        command: None,
        argv: canonical(parsed),
        output: None,
        notes,
        formats: Vec::new(),
        write: None,
        vcf: None,
        inputs: inputs.clone(),
    };
    let index = before.len();
    let mut steps = vec![qc_step];
    steps.extend(added.steps);
    let plan = QcPlan {
        index,
        argv: canonical(parsed),
        args,
        preset,
        pdf,
        added: (index + 1..index + steps.len())
            .filter(|i| !added.own.contains(i) && !added.prepare.contains(i))
            .collect(),
        materialized: added.materialized,
        own: added.own,
        prepare: added.prepare,
        unmet: added.unmet,
        inputs,
    };
    Ok((steps, plan))
}

/// What the qc step did.
pub(super) struct Ended {
    pub files: Vec<Written>,
    /// 0, or 10 and 11 under `--fail-on`.
    pub code: u8,
    /// The verdict in one line.
    pub note: String,
}

/// The verdict over the metrics files the taps of the qc step's stream
/// wrote (or kept from the run before), then the report files, as the
/// `qc` verb writes them. `outcomes` holds every other step.
pub(super) fn finish(
    plan: &QcPlan,
    steps: &[Planned],
    outcomes: &[Option<StepOutcome>],
    prefix: &Path,
    sample: &str,
) -> crate::error::Result<Ended> {
    let reads = &steps[plan.index].reads;
    let mut paths: Vec<&Path> = steps
        .iter()
        .zip(outcomes)
        .filter(|(step, _)| step.reads_only() && &step.reads == reads)
        .filter_map(|(_, outcome)| outcome.as_ref())
        .filter(|outcome| matches!(outcome.status, StepStatus::Ok | StepStatus::Skipped))
        .flat_map(|outcome| &outcome.files)
        .filter(|file| file.kind == "metrics")
        .map(|file| file.path.as_path())
        .collect();
    paths.sort();
    let inputs = paths
        .iter()
        .map(|p| Input::read(p))
        .collect::<bammetrics_verdict::Result<Vec<_>>>()?;
    let args = &plan.args;
    let verdict = evaluate(&plan.preset, &inputs, Some(sample));
    let date = args.date.clone().unwrap_or_else(|| CivilTime::now().iso());
    // The PDF is rendered first: a sample name its fonts cannot set fails
    // the step with nothing written.
    let pdf = match &plan.pdf {
        Some(settings) => Some(settings.render(&verdict, &inputs, &date, args.all_plots)?),
        None => None,
    };
    let mut files = Vec::new();
    let mut write = |path: PathBuf, text: &[u8], kind: &'static str| {
        write_atomic(&path, text)?;
        files.push(Written::new(path, kind));
        crate::error::Result::Ok(())
    };
    write(
        qc::verdict_path(prefix),
        verdict.to_json().as_bytes(),
        "verdict",
    )?;
    if args.writes(ReportForm::Html) {
        let page = bammetrics_report::render(
            &verdict,
            &inputs,
            &Options {
                all_plots: args.all_plots,
                date: &date,
            },
        );
        write(qc::page_path(prefix), page.as_bytes(), "report")?;
    }
    for (target, text) in
        crate::also::texts(prefix, &paths, &args.also, sample, &verdict, args.fail_on)?
    {
        write(target.path, &text, target.also.name())?;
    }
    if let (Some(settings), Some(bytes)) = (&plan.pdf, pdf) {
        let path = qc::pdf_path(prefix);
        settings.save(&path, &bytes, &inputs)?;
        files.push(Written::new(path, "pdf"));
    }
    let code = args.fail_on.map_or(0, |gate| {
        check::exit_code(std::slice::from_ref(&verdict), gate)
    });
    Ok(Ended {
        files,
        code,
        note: verdict_line(&verdict),
    })
}

/// `verdict pass: 6 checks, 6 pass, 0 warn, 0 fail, 0 not measured`.
fn verdict_line(verdict: &Verdict) -> String {
    let c = &verdict.counts;
    format!(
        "verdict {}: {} checks, {} pass, {} warn, {} fail, {} not measured",
        verdict.status.name(),
        c.pass + c.warn + c.fail + c.not_measured,
        c.pass,
        c.warn,
        c.fail,
        c.not_measured
    )
}

/// The outcome of the qc step: its files and its verdict line, or its
/// error.
pub(super) fn outcome(number: usize, ended: crate::error::Result<Ended>) -> (StepOutcome, u8) {
    match ended {
        Ok(ended) => (
            StepOutcome {
                number,
                status: StepStatus::Ok,
                error: None,
                files: ended.files,
                counters: None,
                notes: vec![ended.note],
            },
            ended.code,
        ),
        Err(e) => (
            StepOutcome {
                number,
                status: StepStatus::Failed,
                error: Some(e.to_string()),
                files: Vec::new(),
                counters: None,
                notes: Vec::new(),
            },
            0,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_flag_has_one_long_spelling() {
        let flags = flags();
        for flag in &flags {
            assert_eq!(
                flags.iter().filter(|f| f.long == flag.long).count(),
                1,
                "--{}",
                flag.long
            );
        }
        assert!(help().contains("--preset NAME|FILE"));
        assert!(help().contains("--hs BAITS TARGETS"));
    }

    #[test]
    fn the_canonical_argv_follows_the_flag_order() {
        let args: Vec<String> = ["--also", "junit", "--preset", "wgs", "--all-plots"]
            .iter()
            .map(|a| (*a).to_owned())
            .collect();
        let parsed = crate::grammar::parse(QC, &flags(), &args).expect("parses");
        assert_eq!(
            canonical(&parsed),
            ["--preset=wgs", "--also=junit", "--all-plots=true"]
        );
    }
}

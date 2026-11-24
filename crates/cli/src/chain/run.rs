//! A [`Chain`] run: the steps through
//! [`bammetrics_core::stage::run_chain`], every file written, one
//! [`StepOutcome`] per step.
//!
//! Failure (docs/chaining.md section 7). A tap that fails fails alone:
//! every other step finishes and writes its files. A transform or `write`
//! that fails stops its segment: the other steps of that segment are
//! `stopped`, every later step is `not run`. The exit code is 1 when any
//! step did not finish.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use bammetrics_core::cleanup::MadeDirs;
use bammetrics_core::metric::MetricsOutput;
use bammetrics_core::metrics_file::CivilTime;
use bammetrics_core::registry;
use bammetrics_core::split_sam::{self, Claims};
use bammetrics_core::stage::{
    self, ChainPlan, Counters, FileSource, Halt, Products, Step, materialize,
};
use bammetrics_formats::{Document, markdown};

use super::VERB;
use super::lock::PrefixLock;
use super::prepare::{Chain, Planned, PrefixOrigin, Role};
use super::resume::{self, Resume};
use crate::error::VerbError;
use crate::format::{self, Context, Format};
use crate::manifest::Fingerprint;
use crate::verbs;

/// How a step ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepStatus {
    /// It finished and wrote its files.
    Ok,
    /// It failed by itself.
    Failed,
    /// Another step of its segment failed while it ran: it saw part of
    /// the stream and wrote nothing.
    Stopped,
    /// A step in front of it failed before it started.
    NotRun,
    /// `--resume`: the manifest of the run before proves it done, and its
    /// files stay as that run wrote them.
    Skipped,
}

impl StepStatus {
    pub fn name(self) -> &'static str {
        match self {
            StepStatus::Ok => "ok",
            StepStatus::Failed => "failed",
            StepStatus::Stopped => "stopped",
            StepStatus::NotRun => "not run",
            StepStatus::Skipped => "skipped",
        }
    }
}

/// A file a run wrote: `kind` is `metrics`, `gaps`, `file`, `bam`, `index`,
/// `md5`, or the name of a further format (`json`, `markdown`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Written {
    pub path: PathBuf,
    pub kind: &'static str,
    /// Its fingerprint, taken as soon as the run wrote it, for the
    /// manifest; `None` until then and in a run that writes no manifest.
    pub print: Option<Fingerprint>,
}

impl Written {
    pub fn new(path: PathBuf, kind: &'static str) -> Self {
        Self {
            path,
            kind,
            print: None,
        }
    }
}

/// One step after the run.
#[derive(Debug)]
pub struct StepOutcome {
    /// Its place in the chain, from 1: `Chain::steps[number - 1]`.
    pub number: usize,
    pub status: StepStatus,
    /// Why a step that is not `Ok` ended as it did.
    pub error: Option<String>,
    /// The files it wrote, in write order.
    pub files: Vec<Written>,
    /// A transform's records in, out and changed.
    pub counters: Option<Counters>,
    /// What the step did not do and why (an index asked for on a stream
    /// that is not coordinate sorted).
    pub notes: Vec<String>,
}

/// A finished run.
pub struct Outcome {
    /// The chain as planned: each step's verb, role, stream and files.
    pub chain: Chain,
    /// One per step of [`Chain::steps`], in step order.
    pub steps: Vec<StepOutcome>,
    /// The files of the chain itself (`PREFIX.report.md`).
    pub files: Vec<Written>,
    /// The exit code `bammetrics chain` gives the same run: 0, or 1 when
    /// a step did not finish, or the qc step's 10 or 11 under `--fail-on`.
    pub code: u8,
    /// Lines about the run as a whole: why `--resume` reruns every step,
    /// the files of the run before it removed.
    pub notes: Vec<String>,
}

/// An error of step `label` found before the run.
fn step_error(label: &str, e: &bammetrics_core::Error) -> VerbError {
    let mut error = VerbError::core(VERB, e);
    error.what = format!("{label}: {}", error.what);
    error
}

/// Create the directories above every planned file, then fail when a file
/// cannot be written, as the Picard face fails for its `OUTPUT`: before
/// the first read, exit 1. A file the write itself would refuse (a link
/// another user owns) fails here with that refusal. `made` gets every
/// directory created, so a run that fails removes the empty ones.
fn check_writable(chain: &Chain, made: &mut MadeDirs) -> Result<(), VerbError> {
    let mut create_parent = |label: &str, path: &Path| {
        made.create_parent(path).map_err(|e| {
            let parent = path.parent().unwrap_or(path).display();
            VerbError::input(VERB, format!("{label}: {parent}: {e}"))
        })
    };
    for step in &chain.steps {
        let paths = step
            .files
            .iter()
            .map(|file| file.path.as_path())
            .chain(step.late_prefix.as_deref());
        for path in paths {
            create_parent(&step.label, path)?;
        }
        let files = step.files.iter().map(|file| file.path.as_path());
        registry::check_writable(files, &[]).map_err(|e| step_error(&step.label, &e))?;
    }
    if let Some(report) = &chain.report {
        create_parent("--also markdown", report)?;
        registry::check_writable([report.as_path()], &[])
            .map_err(|e| step_error("--also markdown", &e))?;
    }
    if let Some(manifest) = &chain.manifest {
        create_parent("the manifest", manifest)?;
        registry::check_writable([manifest.as_path()], &[])
            .map_err(|e| step_error("the manifest", &e))?;
        // The manifest is written last, under a temporary name beside it:
        // a directory that takes no new file fails here, before the read.
        bammetrics_core::cleanup::AtomicFile::probe(manifest)
            .map_err(|e| manifest_error(manifest, &e))?;
    }
    Ok(())
}

/// The error of a manifest that cannot be written, naming its directory
/// (the temporary name it is written under means nothing to a user).
fn manifest_error(manifest: &Path, e: &std::io::Error) -> VerbError {
    let dir = match manifest.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir,
        _ => Path::new("."),
    };
    VerbError::input(
        VERB,
        format!(
            "the manifest {}: cannot write in {}: {}",
            manifest.display(),
            dir.display(),
            e
        ),
    )
}

/// A split step: its directory and options, every file the chain reads,
/// and the claims every split of the chain shares.
fn split_step(
    chain: &Chain,
    step: &Planned,
    tool: &'static registry::Tool,
    inv: &registry::Invocation,
    claims: &Claims,
) -> bammetrics_core::Result<Step> {
    let (kind, dir) = split_sam::chain_options(tool, inv)?;
    let mut reads = vec![chain.source_path.clone()];
    if let Some(reference) = &chain.reference {
        reads.extend(bammetrics_core::reference::files_read(reference));
    }
    reads.extend(chain.steps.iter().flat_map(|s| s.inputs.iter().cloned()));
    let split = split_sam::ChainSplit {
        kind,
        dir,
        writer: split_sam::WriterOptions::from_invocation(inv)?,
        reads,
        claims: Some((step.label.clone(), Arc::clone(claims))),
    };
    Ok(split_sam::chain_step(&step.label, split))
}

/// The step of the stage layer for `step`, under its label.
fn core_step(
    chain: &Chain,
    index: usize,
    step: &Planned,
    claims: &Claims,
) -> Result<Step, VerbError> {
    let failed = |e: bammetrics_core::Error| step_error(&step.label, &e);
    let unfit = || super::usage(format!("{}: it cannot be a step", step.label));
    let mut built = match (&step.write, step.tool) {
        (Some(options), _) => Step::write(options.clone(), chain.threads),
        (None, Some(tool)) => {
            let mut inv = tool.parse(&step.argv).map_err(failed)?;
            inv.threads = chain.threads;
            inv.quiet = chain.quiet;
            let sort_bytes = super::planner::sort_budget(chain, index);
            inv.sort_memory = Some(usize::try_from(sort_bytes).unwrap_or(usize::MAX));
            inv.read_dash_from_stdin();
            match step.role {
                Role::Transform => Step::tool_stage(tool, inv).ok_or_else(unfit)?,
                Role::Split => split_step(chain, step, tool, &inv, claims).map_err(failed)?,
                Role::Materialize => Step::materialize(
                    &step.label,
                    materialize::Spec {
                        tool,
                        argv: step.argv.clone(),
                        threads: chain.threads,
                        quiet: chain.quiet,
                        sort_memory: inv.sort_memory,
                        tmp_dir: chain.tmp_dir.clone(),
                        stringency: chain.stringency,
                        output: step.output.clone().filter(|_| !step.stdout),
                    },
                ),
                _ => Step::tool_tap(tool, inv)
                    .map_err(failed)?
                    .ok_or_else(unfit)?,
            }
        }
        (None, None) => return Err(unfit()),
    };
    built.label.clone_from(&step.label);
    Ok(built)
}

/// `file` in each of `formats`: the first to `path`, each further one
/// beside it as `FILE.<extension>`, all from one set of header lines.
fn write_metrics(
    path: &Path,
    file: &MetricsOutput,
    formats: &[Format],
    context: &Context<'_>,
) -> bammetrics_core::Result<Vec<Written>> {
    let rendered =
        format::render_all(formats, context, file).map_err(|e| format::unrendered(path, &e))?;
    // The file and the further formats are one set: a failed write leaves
    // none of them.
    let mut set = format::FileSet::default();
    set.add(path, formats, rendered);
    Ok(set
        .write_with_formats()?
        .into_iter()
        .map(|(target, format)| {
            let kind = if target == path {
                "metrics"
            } else {
                format.name()
            };
            Written::new(target, kind)
        })
        .collect())
}

/// The files of a tool step that finished: the main file (to its path, or
/// to standard output), then every further file. A file the tool named
/// only now (under an `OUTPUT` prefix) is claimed first; when another step
/// writes it, or the chain reads it, this step fails and replaces nothing.
fn write_tool_files(
    chain: &mut Chain,
    index: usize,
    products: &Products,
    sample: Option<&str>,
) -> bammetrics_core::Result<Vec<Written>> {
    let Chain {
        steps,
        outputs,
        inputs,
        ..
    } = chain;
    let step = &steps[index];
    let context = Context {
        tool: step.tool.map_or(step.verb, |tool| tool.name),
        command: step.command.as_deref().unwrap_or_default(),
        sample,
    };
    let picard = [Format::Picard];
    // Each further file with its formats and its planned kind. A file that
    // is not a metrics file (a BED of gaps, a per-target table) has one
    // form.
    let mut extras = Vec::with_capacity(products.output.extra.len());
    for (path, file) in &products.output.extra {
        let planned = step.files.iter().find(|planned| planned.path == *path);
        let metrics = file.text.is_none() && planned.is_none_or(|p| p.kind == "metrics");
        let formats = if metrics { &step.formats[..] } else { &picard };
        if planned.is_none() {
            for late in format::all_files(path, formats) {
                if inputs.owner(&late).is_some() {
                    return Err(bammetrics_core::Error::Usage(format!(
                        "{} is an input; writing it would destroy it",
                        late.display()
                    )));
                }
                outputs.claim(&step.label, &late)?;
            }
        }
        extras.push((path, file, formats, planned.map(|p| p.kind)));
    }
    let mut written = Vec::new();
    match (&products.output.main, &step.output, step.stdout) {
        (Some(main), _, true) => {
            let text = format::render_all(&picard, &context, main)
                .map_err(|e| format::unrendered(Path::new("-"), &e))?
                .remove(0)
                .remove(0);
            registry::write_to(None, |out| out.write_all(&text.bytes))?;
        }
        (Some(main), Some(path), false) => {
            written.extend(write_metrics(path, main, &step.formats, &context)?);
        }
        _ => {}
    }
    for (path, file, formats, kind) in extras {
        let mut files = write_metrics(path, file, formats, &context)?;
        if let (Some(kind), Some(first)) = (kind, files.first_mut()) {
            first.kind = kind;
        }
        written.extend(files);
    }
    Ok(written)
}

/// The kind of a file `write` renamed into place.
fn bam_kind(path: &Path) -> &'static str {
    match path.extension().and_then(|e| e.to_str()) {
        Some("bai") => "index",
        Some("md5") => "md5",
        Some("sam") => "sam",
        _ => "bam",
    }
}

/// What step `index` left, written: its files, its counters, its notes.
fn finish_step(
    chain: &mut Chain,
    index: usize,
    products: &Products,
    sample: Option<&str>,
) -> bammetrics_core::Result<StepOutcome> {
    // `write -` names /dev/stdout: the BAM went to the pipe, no file was
    // written.
    let piped = chain.steps[index].stdout;
    // A materialized tool's own file (a checksum, SAM text) has one form.
    let own = chain.steps[index].role == Role::Materialize;
    let mut files: Vec<Written> = products
        .files
        .iter()
        .filter(|_| !piped)
        .map(|path| Written::new(path.clone(), if own { "file" } else { bam_kind(path) }))
        .collect();
    if chain.steps[index].tool.is_some() {
        files.extend(write_tool_files(chain, index, products, sample)?);
    }
    // The tool wrote its files and fails, as it exits 1 on its own.
    let status = match products.failure {
        Some(_) => StepStatus::Failed,
        None => StepStatus::Ok,
    };
    Ok(StepOutcome {
        number: index + 1,
        status,
        error: products.failure.clone(),
        files,
        counters: products.counters,
        notes: products.notes.clone(),
    })
}

/// The outcome of a step with an error: its own, or the halt another
/// step caused.
fn failed_step(step: &Planned, error: &bammetrics_core::Error) -> StepOutcome {
    let status = match stage::halt(error) {
        Some(Halt::Stopped) => StepStatus::Stopped,
        Some(Halt::NotRun) => StepStatus::NotRun,
        None => StepStatus::Failed,
    };
    let text = error.to_string();
    let own = step
        .tool
        .and_then(|tool| text.strip_prefix(tool.name))
        .and_then(|rest| rest.strip_prefix(": "))
        .unwrap_or(&text);
    StepOutcome {
        number: step.number,
        status,
        error: Some(own.to_owned()),
        files: Vec::new(),
        counters: None,
        notes: Vec::new(),
    }
}

/// The sample a chain is about: `--sample`, else the first `@RG SM` of a
/// file source, else the source name (`stdin` for a pipe, whose header
/// cannot be read a second time). The Python `Run.sample`.
pub fn sample_name(chain: &Chain) -> String {
    if let Some(sample) = &chain.sample {
        return sample.clone();
    }
    let from_header = chain
        .replayable()
        .then(|| bammetrics_core::pass::read_header(&chain.source_path).ok())
        .flatten()
        .and_then(|header| {
            bammetrics_core::levels::read_groups(&header)
                .into_iter()
                .find_map(|group| group.sample)
        });
    from_header.unwrap_or_else(|| verbs::stem(&chain.source))
}

/// The sample of the json, jsonl and parquet forms: [`sample_name`], and
/// none for a piped source without `--sample`. Asked only when a step writes
/// one of those forms, as the verb asks it.
fn sample(chain: &Chain) -> Option<String> {
    let wanted = chain.steps.iter().any(|step| {
        step.formats
            .iter()
            .any(|f| matches!(f, Format::Json | Format::Jsonl | Format::Parquet))
    });
    (wanted && (chain.sample.is_some() || chain.replayable())).then(|| sample_name(chain))
}

/// `PREFIX.report.md` of `--also markdown`: every metrics table the run
/// wrote, read back from its Picard file, in step order.
fn write_report(
    report: &Path,
    steps: &[StepOutcome],
    sample: Option<&str>,
) -> bammetrics_core::Result<Written> {
    let mut docs = Vec::new();
    let metrics = steps
        .iter()
        .flat_map(|step| &step.files)
        .filter(|file| file.kind == "metrics");
    for file in metrics {
        let unreadable = |source| bammetrics_core::Error::Input {
            path: file.path.clone(),
            source,
        };
        let text = super::read_capped(&file.path, super::REPORT_FILE_CAP).map_err(unreadable)?;
        let mut doc: Document = bammetrics_formats::picard::read_or_empty(&text)
            .map_err(|e| unreadable(std::io::Error::other(e.to_string())))?;
        doc.sample = sample.map(str::to_owned);
        docs.push(doc);
    }
    let text = markdown::write(&docs.iter().collect::<Vec<_>>(), &[]);
    registry::write_to(Some(report), |out| out.write_all(text.as_bytes()))?;
    Ok(Written::new(report.to_path_buf(), "markdown"))
}

/// Run `chain` and write every file. `Err` for what stops the chain
/// before its first record: an output that cannot be written, a source
/// that cannot be read. A step that fails is in the [`Outcome`]: the
/// other steps still write their files. Nothing is printed.
///
/// A run that does not end with every step finished removes the output
/// directories it made that are still empty, and no other directory.
pub fn execute(chain: Chain) -> Result<Outcome, VerbError> {
    let mut made = MadeDirs::default();
    let outcome = run_steps(chain, &mut made)?;
    if outcome.code == 0 {
        made.keep();
    }
    Ok(outcome)
}

/// What a step of the stage layer stands for: a step of the chain, or the
/// write of the stream a transform keeps.
#[derive(Clone, Copy)]
enum Slot {
    Step(usize),
    Keep(usize),
}

/// The steps of the stage layer for the steps of `chain` the read runs
/// (all of them without `--resume`), each kept stream's write right after
/// its transform.
fn core_steps(
    chain: &Chain,
    resume: &Resume,
    claims: &Claims,
) -> Result<(Vec<Step>, Vec<Slot>), VerbError> {
    let mut steps = Vec::new();
    let mut slots = Vec::new();
    for (index, step) in chain.steps.iter().enumerate() {
        if resume.skipped(index) || !resume::in_read(chain, index) {
            continue;
        }
        steps.push(core_step(chain, index, step, claims)?);
        slots.push(Slot::Step(index));
        if let Some(kept) = chain.kept.iter().find(|kept| kept.index == index) {
            steps.push(super::keep::step(kept, chain.threads, &step.label));
            slots.push(Slot::Keep(index));
        }
    }
    Ok((steps, slots))
}

/// The steps of the qc step that run after the reads (an own read, the
/// panel's interval list): the tool on its own, its files written as a
/// tap's.
fn after_reads(
    chain: &mut Chain,
    resume: &Resume,
    outcomes: &mut [Option<StepOutcome>],
    sample: Option<&str>,
) {
    let indices: Vec<usize> = (0..chain.steps.len())
        .filter(|&i| !resume::in_read(chain, i) && chain.steps[i].role != Role::Qc)
        .filter(|&i| !resume.skipped(i))
        .collect();
    for index in indices {
        let step = &chain.steps[index];
        let result = step.tool.map_or_else(
            || Err(bammetrics_core::Error::Usage("no tool".to_owned())),
            |tool| {
                let mut inv = tool.parse(&step.argv)?;
                inv.threads = chain.threads;
                inv.quiet = true;
                (tool.run)(&inv).map(|output| Products {
                    output,
                    ..Products::default()
                })
            },
        );
        let finished = result.and_then(|products| finish_step(chain, index, &products, sample));
        outcomes[index] = Some(finished.unwrap_or_else(|e| failed_step(&chain.steps[index], &e)));
    }
}

/// The one read of the chain: the steps the run reads with, from the source
/// or the kept stream `--resume` starts at. No step, no read.
fn read_stream(
    chain: &Chain,
    resume: &Resume,
    claims: &Claims,
) -> Result<(Vec<bammetrics_core::Result<Products>>, Vec<Slot>), VerbError> {
    let (steps, slots) = core_steps(chain, resume, claims)?;
    if steps.is_empty() {
        return Ok((Vec::new(), slots));
    }
    let source = match &resume.start {
        Some((_, kept)) => FileSource::new(kept, chain.reference.as_deref())
            .decoded_from_cram(bammetrics_core::cram::is_cram(&chain.source_path)),
        None => FileSource::new(&chain.source_path, chain.reference.as_deref()),
    };
    let plan = ChainPlan {
        source: Box::new(source.with_stringency(chain.stringency)),
        steps,
        threads: chain.threads,
        tmp_dir: chain.tmp_dir.clone(),
        lower: !chain.stage_executor,
    };
    let results = stage::run_chain(plan, None).map_err(|e| VerbError::core(VERB, &e))?;
    Ok((results, slots))
}

/// The outcome of each step the read ran, its files written; a kept
/// stream's files join its transform's.
fn gather(
    chain: &mut Chain,
    slots: Vec<Slot>,
    results: Vec<bammetrics_core::Result<Products>>,
    sample: Option<&str>,
) -> Vec<Option<StepOutcome>> {
    let mut outcomes: Vec<Option<StepOutcome>> = chain.steps.iter().map(|_| None).collect();
    for (slot, result) in slots.into_iter().zip(results) {
        match slot {
            Slot::Step(index) => {
                let finished =
                    result.and_then(|products| finish_step(chain, index, &products, sample));
                outcomes[index] =
                    Some(finished.unwrap_or_else(|e| failed_step(&chain.steps[index], &e)));
            }
            Slot::Keep(index) => {
                let Some(outcome) = outcomes[index].as_mut() else {
                    continue;
                };
                match result {
                    Ok(products) => outcome.files.extend(super::keep::written(&products)),
                    Err(e) if outcome.status == StepStatus::Ok => {
                        outcome.status = StepStatus::Failed;
                        outcome.error = Some(format!("keeping its stream: {e}"));
                    }
                    Err(_) => {}
                }
            }
        }
    }
    outcomes
}

/// Write the manifest of `run` and, after a resumed run that finished,
/// remove the files of the run before that this run does not keep.
fn write_manifest(
    chain: &Chain,
    run: &resume::Run,
    notes: &mut Vec<String>,
) -> Result<Option<Written>, VerbError> {
    let Some(path) = &chain.manifest else {
        return Ok(None);
    };
    let manifest = resume::manifest(chain, run);
    crate::report::write_atomic(path, manifest.to_json()).map_err(|e| match e {
        crate::error::Error::Io { source, .. } => manifest_error(path, &source),
        other => VerbError::input(VERB, format!("the manifest: {other}")),
    })?;
    if let (0, Some(previous)) = (run.code, &run.resume.previous) {
        notes.extend(resume::remove_stale(chain, previous, &manifest).notes());
    }
    if chain.resume {
        notes.extend(resume::orphans(chain).into_iter().map(|path| {
            format!(
                "--resume: left in place: {} (a temporary file of a run that was killed; \
                 remove it)",
                path.display()
            )
        }));
    }
    Ok(Some(Written::new(path.clone(), "manifest")))
}

/// Fingerprint each file a step wrote that has none yet: right after the
/// read or step that wrote it, so the manifest records the bytes this run
/// wrote.
fn stamp(outcomes: &mut [Option<StepOutcome>]) {
    for file in outcomes
        .iter_mut()
        .flatten()
        .flat_map(|o| o.files.iter_mut())
    {
        if file.print.is_none() {
            file.print = crate::manifest::v2::fingerprint(&file.path);
        }
    }
}

/// [`execute`], with the directories it makes recorded in `made`.
fn run_steps(mut chain: Chain, made: &mut MadeDirs) -> Result<Outcome, VerbError> {
    super::planner::check_space(&chain)?;
    check_writable(&chain, made)?;
    // One run per manifest, from before the read to the end.
    let _lock = chain
        .manifest
        .as_deref()
        .map(PrefixLock::take)
        .transpose()?;
    let manifested = chain.manifest.is_some();
    // A chain that writes no manifest fingerprints nothing: it reads what
    // it read before this wave.
    let keys = if chain.manifest.is_some() {
        resume::keys(&chain)
    } else {
        Vec::new()
    };
    let resume = resume::decide(&chain, &keys);
    let mut notes: Vec<String> = resume
        .note
        .iter()
        .map(|note| format!("--resume: {note}"))
        .collect();
    let started = CivilTime::now().iso();
    let clock = Instant::now();
    if let Some(qc) = &chain.qc
        && let Some(pdf) = &qc.pdf
    {
        let sample = sample_name(&chain);
        pdf.check_sample(&sample)
            .map_err(|e| super::usage(format!("{}: {e}", chain.steps[qc.index].label)))?;
    }
    // The planned files again, for the split steps to claim the files they
    // name once the stream is read.
    let claims: Claims = Arc::new(Mutex::new(super::claims::claim_outputs(
        &chain.steps,
        chain.report.as_deref(),
    )?));
    let (results, slots) = read_stream(&chain, &resume, &claims)?;
    let sample = sample(&chain);
    let mut outcomes = gather(&mut chain, slots, results, sample.as_deref());
    if manifested {
        stamp(&mut outcomes);
    }
    for (index, key) in keys.iter().enumerate() {
        if resume.skipped(index) {
            outcomes[index] = Some(resume::skipped_outcome(
                chain.steps[index].number,
                key,
                &resume,
            ));
        }
    }
    after_reads(&mut chain, &resume, &mut outcomes, sample.as_deref());
    if manifested {
        stamp(&mut outcomes);
    }
    let mut gate = 0;
    if let Some(qc) = &chain.qc {
        let ended = super::qc_step::finish(
            qc,
            &chain.steps,
            &outcomes,
            &chain.prefix,
            &sample_name(&chain),
        );
        let (outcome, code) = super::qc_step::outcome(chain.steps[qc.index].number, ended);
        gate = code;
        outcomes[qc.index] = Some(outcome);
        stamp(&mut outcomes);
    }
    let outcomes: Vec<StepOutcome> = outcomes
        .into_iter()
        .zip(&chain.steps)
        .map(|(outcome, step)| {
            outcome.unwrap_or_else(|| StepOutcome {
                number: step.number,
                status: StepStatus::NotRun,
                error: Some("not run: the read stopped before it".to_owned()),
                files: Vec::new(),
                counters: None,
                notes: Vec::new(),
            })
        })
        .collect();
    let failed = outcomes
        .iter()
        .any(|step| !matches!(step.status, StepStatus::Ok | StepStatus::Skipped));
    let code = if failed { 1 } else { gate };
    let mut files = Vec::new();
    if let Some(report) = &chain.report {
        let written = write_report(report, &outcomes, sample.as_deref())
            .map_err(|e| step_error("--also markdown", &e))?;
        files.push(written);
    }
    if chain.manifest.is_some() {
        let run = resume::Run {
            keys: &keys,
            resume: &resume,
            outcomes: &outcomes,
            sample: sample_name(&chain),
            started,
            millis: u64::try_from(clock.elapsed().as_millis()).unwrap_or(u64::MAX),
            code,
        };
        files.extend(write_manifest(&chain, &run, &mut notes)?);
    }
    Ok(Outcome {
        chain,
        steps: outcomes,
        files,
        code,
        notes,
    })
}

/// `1 step`, `3 steps`.
fn count(n: usize, word: &str) -> String {
    if n == 1 {
        format!("1 {word}")
    } else {
        format!("{n} {word}s")
    }
}

/// How a step that failed by itself reads in the summary: a tap fails
/// alone, a transform or `write` stops its segment when another step was
/// halted by it.
fn failure_words(step: &Planned, halted: bool) -> &'static str {
    let alone = step.role == Role::Tap || (step.role == Role::Materialize && step.leaves.is_none());
    match (alone, halted) {
        (true, _) => "failed alone",
        (_, true) => "failed and stopped its segment",
        (_, false) => "failed",
    }
}

/// The lines `bammetrics chain` prints on standard error after a run: a
/// transform's record counts and each step's notes (not with `quiet`),
/// each file written (with `verbose`), one line per step that did not
/// finish saying how it ended, then the count line.
pub fn summary(outcome: &Outcome) -> Vec<String> {
    let chain = &outcome.chain;
    let halted = outcome
        .steps
        .iter()
        .any(|step| matches!(step.status, StepStatus::Stopped | StepStatus::NotRun));
    let mut lines = Vec::new();
    if !chain.quiet {
        lines.extend(outcome.notes.iter().cloned());
    }
    let mut ended = [0usize; 4];
    let mut skipped = 0;
    for (step, planned) in outcome.steps.iter().zip(&chain.steps) {
        let label = &planned.label;
        let error = step.error.as_deref().unwrap_or_default();
        match step.status {
            StepStatus::Ok => {
                ended[0] += 1;
                if let (Some(counters), false) = (step.counters, chain.quiet) {
                    lines.push(format!(
                        "{label}: {} records in, {} out, {} changed",
                        counters.records_in, counters.records_out, counters.changed
                    ));
                }
                if !chain.quiet {
                    lines.extend(step.notes.iter().map(|note| format!("{label}: {note}")));
                }
                if chain.verbose {
                    let wrote =
                        |file: &super::Written| format!("{label}: wrote {}", file.path.display());
                    lines.extend(step.files.iter().map(wrote));
                }
            }
            StepStatus::Failed => {
                ended[1] += 1;
                let words = failure_words(planned, halted);
                lines.push(format!("{label}: {words}: {error}"));
            }
            StepStatus::Stopped => {
                ended[2] += 1;
                lines.push(format!("{label}: {error}"));
            }
            StepStatus::NotRun => {
                ended[3] += 1;
                lines.push(format!("{label}: {error}"));
            }
            StepStatus::Skipped => {
                ended[0] += 1;
                skipped += 1;
                if !chain.quiet {
                    lines.push(format!(
                        "{label}: skipped: the manifest proves it done, its files stay"
                    ));
                }
            }
        }
    }
    let files = outcome
        .steps
        .iter()
        .filter(|step| step.status != StepStatus::Skipped)
        .map(|step| step.files.len())
        .sum::<usize>()
        + outcome.files.len();
    let total = outcome.steps.len();
    let [finished, failed, stopped, not_run] = ended;
    if finished == total {
        if chain.quiet {
            return lines;
        }
        if skipped > 0 {
            lines.push(format!(
                "{} finished, {skipped} of them skipped by --resume, {} written",
                count(total, "step"),
                count(files, "file")
            ));
        } else {
            lines.push(format!(
                "{} finished, {} written",
                count(total, "step"),
                count(files, "file")
            ));
        }
        return lines;
    }
    let mut parts = vec![format!("{failed} failed")];
    if stopped > 0 {
        parts.push(format!("{stopped} stopped"));
    }
    if not_run > 0 {
        parts.push(format!("{not_run} not run"));
    }
    lines.push(format!(
        "{finished} of {} finished, {}; {} written",
        count(total, "step"),
        parts.join(", "),
        count(files, "file")
    ));
    lines
}

/// The command: run `chain`, print what it says on standard error, return
/// the exit code.
pub(super) fn command(chain: Chain) -> Result<u8, VerbError> {
    let say = |line: &str| eprintln!("bammetrics {VERB}: {line}");
    if chain.prefix_origin == PrefixOrigin::SourceName && chain.prefix_used && !chain.quiet {
        say(&format!("no -o, files go to {}.*", chain.prefix.display()));
    }
    for step in &chain.steps {
        if let (true, Some(command)) = (chain.verbose, &step.command) {
            say(&format!("{}: {command}", step.label));
        }
        if !chain.quiet {
            for note in &step.notes {
                say(&format!("{}: {note}", step.label));
            }
        }
    }
    let outcome = execute(chain)?;
    for line in summary(&outcome) {
        say(&line);
    }
    Ok(outcome.code)
}

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

use bammetrics_core::metric::MetricsOutput;
use bammetrics_core::registry;
use bammetrics_core::stage::{self, ChainPlan, Counters, FileSource, Halt, Products, Step};
use bammetrics_formats::{Document, markdown};

use super::VERB;
use super::prepare::{Chain, Planned, PrefixOrigin, Role};
use crate::error::VerbError;
use crate::format::{self, Context, Format};
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
}

impl StepStatus {
    pub fn name(self) -> &'static str {
        match self {
            StepStatus::Ok => "ok",
            StepStatus::Failed => "failed",
            StepStatus::Stopped => "stopped",
            StepStatus::NotRun => "not run",
        }
    }
}

/// A file a run wrote: `kind` is `metrics`, `gaps`, `file`, `bam`, `index`,
/// `md5`, or the name of a further format (`json`, `markdown`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Written {
    pub path: PathBuf,
    pub kind: &'static str,
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
    /// a step did not finish.
    pub code: u8,
}

/// An error of step `label` found before the run.
fn step_error(label: &str, e: &bammetrics_core::Error) -> VerbError {
    let mut error = VerbError::core(VERB, e);
    error.what = format!("{label}: {}", error.what);
    error
}

/// Create the directories above every planned file, then fail when a file
/// cannot be written, as the Picard face fails for its `OUTPUT`: before
/// the first read, exit 1.
fn check_writable(chain: &Chain) -> Result<(), VerbError> {
    let labelled = |label: &str, mut error: VerbError| {
        error.what = format!("{label}: {}", error.what);
        VERB.clone_into(&mut error.verb);
        error
    };
    for step in &chain.steps {
        let paths = step
            .files
            .iter()
            .map(|file| file.path.as_path())
            .chain(step.late_prefix.as_deref());
        for path in paths {
            verbs::create_parent(VERB, path).map_err(|e| labelled(&step.label, e))?;
        }
        let files = step.files.iter().map(|file| file.path.as_path());
        registry::check_writable(files, &[]).map_err(|e| step_error(&step.label, &e))?;
    }
    if let Some(report) = &chain.report {
        verbs::create_parent(VERB, report)?;
        registry::check_writable([report.as_path()], &[])
            .map_err(|e| step_error("--also markdown", &e))?;
    }
    Ok(())
}

/// The step of the stage layer for `step`, under its label.
fn core_step(chain: &Chain, step: &Planned) -> Result<Step, VerbError> {
    let failed = |e: bammetrics_core::Error| step_error(&step.label, &e);
    let unfit = || super::usage(format!("{}: it cannot be a step", step.label));
    let mut built = match (&step.write, step.tool) {
        (Some(options), _) => Step::write(options.clone(), chain.threads),
        (None, Some(tool)) => {
            let mut inv = tool.parse(&step.argv).map_err(failed)?;
            inv.threads = chain.threads;
            inv.read_dash_from_stdin();
            match step.role {
                Role::Transform => Step::tool_stage(tool, inv).ok_or_else(unfit)?,
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
    let texts = format::render_all(formats, context, file);
    let mut written = Vec::with_capacity(formats.len());
    for (index, (format, text)) in formats.iter().zip(texts).enumerate() {
        let (target, kind) = if index == 0 {
            (path.to_path_buf(), "metrics")
        } else {
            (format::sibling(path, format.extension()), format.name())
        };
        registry::write_to(Some(&target), |out| out.write_all(text.as_bytes()))?;
        written.push(Written { path: target, kind });
    }
    Ok(written)
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
            let siblings = formats[1..]
                .iter()
                .map(|format| format::sibling(path, format.extension()));
            for late in std::iter::once(path.clone()).chain(siblings) {
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
            let text = format::render_all(&picard, &context, main).remove(0);
            registry::write_to(None, |out| out.write_all(text.as_bytes()))?;
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
    let mut files: Vec<Written> = products
        .files
        .iter()
        .filter(|_| !piped)
        .map(|path| Written {
            path: path.clone(),
            kind: bam_kind(path),
        })
        .collect();
    if chain.steps[index].tool.is_some() {
        files.extend(write_tool_files(chain, index, products, sample)?);
    }
    Ok(StepOutcome {
        number: index + 1,
        status: StepStatus::Ok,
        error: None,
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

/// The sample of the json and jsonl forms: [`sample_name`], and none for a
/// piped source without `--sample`. Asked only when a step writes one of
/// those forms.
fn sample(chain: &Chain) -> Option<String> {
    let wanted = chain.steps.iter().any(|step| {
        step.formats
            .iter()
            .any(|f| matches!(f, Format::Json | Format::Jsonl))
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
        let text = std::fs::read_to_string(&file.path).map_err(unreadable)?;
        let mut doc: Document = bammetrics_formats::picard::read_or_empty(&text)
            .map_err(|e| unreadable(std::io::Error::other(e.to_string())))?;
        doc.sample = sample.map(str::to_owned);
        docs.push(doc);
    }
    let text = markdown::write(&docs.iter().collect::<Vec<_>>(), &[]);
    registry::write_to(Some(report), |out| out.write_all(text.as_bytes()))?;
    Ok(Written {
        path: report.to_path_buf(),
        kind: "markdown",
    })
}

/// Run `chain` and write every file. `Err` for what stops the chain
/// before its first record: an output that cannot be written, a source
/// that cannot be read. A step that fails is in the [`Outcome`]: the
/// other steps still write their files. Nothing is printed.
pub fn execute(mut chain: Chain) -> Result<Outcome, VerbError> {
    check_writable(&chain)?;
    let steps = chain
        .steps
        .iter()
        .map(|step| core_step(&chain, step))
        .collect::<Result<Vec<_>, _>>()?;
    let plan = ChainPlan {
        source: Box::new(FileSource::new(
            &chain.source_path,
            chain.reference.as_deref(),
        )),
        steps,
        threads: chain.threads,
        tmp_dir: chain.tmp_dir.clone(),
        lower: !chain.stage_executor,
    };
    let results = stage::run_chain(plan, None).map_err(|e| VerbError::core(VERB, &e))?;
    let sample = sample(&chain);
    let mut outcomes = Vec::with_capacity(results.len());
    for (index, result) in results.into_iter().enumerate() {
        let finished = result
            .and_then(|products| finish_step(&mut chain, index, &products, sample.as_deref()));
        outcomes.push(finished.unwrap_or_else(|e| failed_step(&chain.steps[index], &e)));
    }
    let code = u8::from(outcomes.iter().any(|step| step.status != StepStatus::Ok));
    let mut files = Vec::new();
    if let Some(report) = &chain.report {
        let written = write_report(report, &outcomes, sample.as_deref())
            .map_err(|e| step_error("--also markdown", &e))?;
        files.push(written);
    }
    Ok(Outcome {
        chain,
        steps: outcomes,
        files,
        code,
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
fn failure_words(role: Role, halted: bool) -> &'static str {
    match (role, halted) {
        (Role::Tap, _) => "failed alone",
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
    let mut ended = [0usize; 4];
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
                let words = failure_words(planned.role, halted);
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
        }
    }
    let files = outcome
        .steps
        .iter()
        .map(|step| step.files.len())
        .sum::<usize>()
        + outcome.files.len();
    let total = outcome.steps.len();
    let [finished, failed, stopped, not_run] = ended;
    if finished == total {
        if !chain.quiet {
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

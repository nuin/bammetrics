//! `chain --resume` and the chain's manifest (docs/chaining.md section 27).
//!
//! Each step has a key: an id (a transform's stream, else `VERB@STREAM`),
//! the ids of the transforms in front of it, its canonical argv, the hash
//! of that argv with the tool version, the formats, the sample and the
//! source's validation stringency, and the fingerprints of the files it
//! reads by name. A step is skipped when the record of the run before
//! under its id has the same verb and hash, ended `ok`, read files with the
//! same fingerprints, and every file it wrote is still there with its
//! fingerprint, and when every transform in front of it is skipped too.
//! A step that writes standard output and the qc step always run. A
//! materialized tool that writes a BAM (`merge-sam-files`) counts as a
//! transform here: its id is its stream, and the steps behind it run when
//! it does.
//!
//! The steps that run share reads. The read starts at the last transform
//! in front of the first step that runs whose stream is in a BAM the run
//! before wrote and that is still there: the stream it kept (`--keep`), or
//! the BAM a `write` step of that stream wrote when that step is skipped
//! (the BAM `-o out.bam` names). Else it starts at the source. Every
//! transform between that start and the last step that runs runs too, to
//! feed the steps after it.
//!
//! A file is one file whatever its spelling: paths are compared through
//! [`resolve`], the paths of the run before against the working directory
//! that run recorded. The source and the reference are fingerprinted before
//! the first read, each output as soon as the run wrote it.

use std::collections::{HashMap, HashSet};
use std::path::{Component, Path, PathBuf};

use bammetrics_core::registry::{Kind, Ty};
use bammetrics_verdict::bytes_hash;

use super::prepare::{Chain, Planned, Role};
use super::run::{StepOutcome, StepStatus};
use crate::manifest::v2::{self, AnyManifest, ChainManifest, FileRecord, Kept, StepRecord};

/// `path` as an absolute path with every link of the part that exists
/// followed and `.` and `..` taken out: one name per file whatever the
/// spelling (`m.bam`, `./m.bam`, `$PWD/m.bam`, a link to it, `../d/m.bam`
/// from a subdirectory). `-` and a path under `/dev` stay as they are.
pub(super) fn resolve(path: &Path) -> PathBuf {
    match std::env::current_dir() {
        Ok(dir) => resolve_in(&dir, path),
        Err(_) => path.to_path_buf(),
    }
}

/// [`resolve`] with `dir` as the working directory.
pub(super) fn resolve_in(dir: &Path, path: &Path) -> PathBuf {
    if path == Path::new("-") || path.starts_with("/dev") {
        return path.to_path_buf();
    }
    let mut resolved = PathBuf::new();
    let mut exists = true;
    for component in dir.join(path).components() {
        match component {
            Component::Prefix(_) | Component::RootDir => resolved.push(component.as_os_str()),
            Component::CurDir => {}
            Component::ParentDir => {
                resolved.pop();
            }
            Component::Normal(name) => {
                resolved.push(name);
                if exists {
                    match std::fs::canonicalize(&resolved) {
                        Ok(real) => resolved = real,
                        Err(_) => exists = false,
                    }
                }
            }
        }
    }
    resolved
}

/// The working directory of this run, resolved.
fn directory() -> PathBuf {
    resolve(Path::new("."))
}

/// What the run reads first, fingerprinted before the first read: a
/// source swapped during the run then shows as changed to the next one.
#[derive(Clone, Debug)]
pub struct Snapshot {
    /// The working directory, resolved: the paths of the manifest are
    /// relative to it.
    pub directory: PathBuf,
    pub source: FileRecord,
    pub reference: Option<FileRecord>,
}

impl Snapshot {
    pub fn take(chain: &Chain) -> Self {
        Self {
            directory: directory(),
            source: FileRecord::input(&chain.source_path),
            reference: chain.reference.as_deref().map(FileRecord::input),
        }
    }
}

/// What identifies a step in the manifest.
#[derive(Clone, Debug)]
pub struct Key {
    pub id: String,
    /// The ids of the transforms in front of it.
    pub upstream: Vec<String>,
    pub argv: Vec<String>,
    pub options_hash: String,
    /// The files it reads by name, fingerprinted before the run.
    pub inputs: Vec<FileRecord>,
}

/// `STRICT`, `LENIENT` or `SILENT`: how the chain reads a SAM text source.
fn stringency_name(chain: &Chain) -> &'static str {
    use bammetrics_core::sam_text::Stringency;
    match chain.stringency {
        Stringency::Strict => "STRICT",
        Stringency::Lenient => "LENIENT",
        Stringency::Silent => "SILENT",
    }
}

/// The id a step has in the manifest, before a `#2` for a repeat.
fn base_id(step: &Planned) -> String {
    match (step.role, &step.leaves) {
        (Role::Transform | Role::Materialize, Some(stream)) => stream.clone(),
        (Role::Qc, _) => super::QC.to_owned(),
        _ => format!("{}@{}", step.verb, step.reads),
    }
}

/// The canonical argv of step `index`.
fn argv(chain: &Chain, index: usize) -> Vec<String> {
    let step = &chain.steps[index];
    if let Some(options) = &step.write {
        return vec![
            format!("OUTPUT={}", options.path.display()),
            format!("FORMAT={}", options.format.name()),
            format!("COMPRESSION_LEVEL={}", options.compression_level),
            format!("CREATE_INDEX={}", options.create_index),
            format!("CREATE_MD5_FILE={}", options.create_md5),
        ];
    }
    if step.role == Role::Qc {
        return chain
            .qc
            .as_ref()
            .map(|qc| qc.argv.clone())
            .unwrap_or_default();
    }
    let mut argv = step.argv.clone();
    argv.extend(
        chain
            .kept
            .iter()
            .filter(|kept| kept.index == index)
            .map(super::keep::argv),
    );
    argv
}

/// The names of the path options of `step`'s tool: their values are
/// files, hashed by [`resolve`] so a respelled path hashes the same.
fn path_names(step: &Planned) -> Vec<&'static str> {
    if step.write.is_some() {
        return vec!["OUTPUT"];
    }
    step.tool.map_or_else(Vec::new, |tool| {
        tool.table()
            .into_iter()
            .filter(|o| o.ty == Ty::Path || o.kind == Kind::Output)
            .flat_map(|o| o.names.iter().copied())
            .collect()
    })
}

/// `argv` as the options hash reads it: each file value resolved.
fn hashed_argv(step: &Planned, argv: &[String]) -> Vec<String> {
    let names = path_names(step);
    argv.iter()
        .map(|arg| match arg.split_once('=') {
            Some((key, value))
                if (key == super::keep::ARGV_KEY || names.contains(&key))
                    && value != "null"
                    && !value.is_empty() =>
            {
                format!("{key}={}", resolve(Path::new(value)).display())
            }
            _ => arg.clone(),
        })
        .collect()
}

/// The files step `index` reads by name: the source is the manifest's own
/// record, and a pipe has no fingerprint.
fn inputs(chain: &Chain, index: usize) -> Vec<FileRecord> {
    let step = &chain.steps[index];
    let named = match (step.role, &chain.qc) {
        (Role::Qc, Some(qc)) => &qc.inputs,
        _ => &step.inputs,
    };
    let mut files: Vec<FileRecord> = Vec::new();
    for path in named {
        if path == &chain.source_path || files.iter().any(|f| &f.path == path) {
            continue;
        }
        files.push(FileRecord::input(path));
    }
    files
}

/// The key of every step, its input files fingerprinted now.
pub fn keys(chain: &Chain) -> Vec<Key> {
    let mut seen: HashMap<String, usize> = HashMap::new();
    let mut transforms: Vec<String> = Vec::new();
    let mut keys = Vec::with_capacity(chain.steps.len());
    let formats = |step: &Planned| {
        step.formats
            .iter()
            .map(|f| f.name())
            .collect::<Vec<_>>()
            .join(",")
    };
    for (index, step) in chain.steps.iter().enumerate() {
        let base = base_id(step);
        let n = seen.entry(base.clone()).or_insert(0);
        *n += 1;
        let id = if *n == 1 { base } else { format!("{base}#{n}") };
        let argv = argv(chain, index);
        let hashed = [
            format!("version={}", env!("CARGO_PKG_VERSION")),
            format!("formats={}", formats(step)),
            format!("sample={}", chain.sample.as_deref().unwrap_or_default()),
            format!("stringency={}", stringency_name(chain)),
        ];
        let text = hashed_argv(step, &argv)
            .iter()
            .chain(&hashed)
            .map(String::as_str)
            .collect::<Vec<_>>()
            .join("\0");
        keys.push(Key {
            id: id.clone(),
            upstream: transforms.clone(),
            options_hash: bytes_hash(text.into_bytes()),
            argv,
            inputs: inputs(chain, index),
        });
        if step.hands_on() {
            transforms.push(id);
        }
    }
    keys
}

/// Why a step runs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Why {
    /// No `--resume`, or no manifest to resume from.
    Run,
    /// The manifest has no record under its id.
    NoRecord,
    /// Its verb, options, formats or the tool version changed.
    Options,
    /// It did not end `ok` in the run before.
    Status(String),
    /// A file it reads changed.
    Input(PathBuf),
    /// A file it wrote is missing or changed.
    Output(PathBuf),
    /// It writes standard output, which cannot be compared.
    Stdout,
    /// The transform with this number, in front of it, runs.
    Upstream(usize),
    /// It runs to make the stream of the step with this number.
    Feeds(usize),
    /// The qc step: the verdict and report are always written.
    Verdict,
}

impl Why {
    pub fn text(&self) -> String {
        match self {
            Why::Run => "runs".to_owned(),
            Why::NoRecord => "runs: the manifest has no record of it".to_owned(),
            Why::Options => "runs: its options changed".to_owned(),
            Why::Status(status) => format!("runs: it ended {status} in the run before"),
            Why::Input(path) => format!("runs: {} changed", path.display()),
            Why::Output(path) => format!("runs: {} is missing or changed", path.display()),
            Why::Stdout => "runs: it writes standard output".to_owned(),
            Why::Upstream(n) => format!("runs: step {n} in front of it runs"),
            Why::Feeds(n) => format!("runs: its stream feeds step {n}"),
            Why::Verdict => "runs: the verdict and report are written on every run".to_owned(),
        }
    }
}

/// What a run resumes from.
pub struct Resume {
    /// The source and reference as this run found them before its read;
    /// `None` for a run that writes no manifest.
    pub snapshot: Option<Snapshot>,
    /// The manifest of the run before, when it is one of a chain.
    pub previous: Option<ChainManifest>,
    /// Why every step runs, when that is so for the whole run.
    pub note: Option<String>,
    /// One per step: `None` for a step that is skipped, else why it runs.
    pub runs: Vec<Option<Why>>,
    /// The BAM the read starts from, a kept stream or the BAM a skipped
    /// `write` step wrote of it: the transform's index and its
    /// BAM. `None` reads the source.
    pub start: Option<(usize, PathBuf)>,
}

impl Resume {
    fn everything(
        chain: &Chain,
        snapshot: Option<Snapshot>,
        previous: Option<ChainManifest>,
        note: Option<String>,
    ) -> Self {
        Self {
            snapshot,
            previous,
            note,
            runs: vec![Some(Why::Run); chain.steps.len()],
            start: None,
        }
    }

    pub fn skipped(&self, index: usize) -> bool {
        self.runs[index].is_none()
    }
}

/// Why the whole run reruns: the manifest is of another version, the
/// source, the reference or the stringency changed.
fn whole_run(snapshot: &Snapshot, previous: &ChainManifest, stringency: &str) -> Option<String> {
    if previous.tool_version != env!("CARGO_PKG_VERSION") {
        return Some(format!(
            "the manifest was written by bammetrics {}",
            previous.tool_version
        ));
    }
    let source = &snapshot.source;
    if source.print.is_none() {
        return Some(format!(
            "the source {} is not a regular file",
            source.path.display()
        ));
    }
    let then = resolve_in(&previous.directory, &previous.source.path);
    let now = resolve(&source.path);
    if then != now {
        let (then, now) = if previous.source.path == source.path {
            (then, now)
        } else {
            (previous.source.path.clone(), source.path.clone())
        };
        return Some(format!(
            "the run before read {}; this run reads {}",
            then.display(),
            now.display()
        ));
    }
    if previous.source.print != source.print {
        return Some(format!("the source {} changed", source.path.display()));
    }
    let same = match (&snapshot.reference, &previous.reference) {
        (None, None) => true,
        (Some(new), Some(old)) => {
            resolve(&new.path) == resolve_in(&previous.directory, &old.path)
                && new.print == old.print
        }
        _ => false,
    };
    if !same {
        return Some("the reference changed".to_owned());
    }
    if previous.stringency != stringency {
        return Some(format!(
            "the validation stringency was {}",
            previous.stringency
        ));
    }
    None
}

/// Why step `index` cannot keep its record of the run before; `None` when
/// it can.
fn stale(chain: &Chain, index: usize, key: &Key, previous: &ChainManifest) -> Option<Why> {
    let step = &chain.steps[index];
    if step.role == Role::Qc {
        return Some(Why::Verdict);
    }
    if step.stdout {
        return Some(Why::Stdout);
    }
    let Some(old) = previous.step(&key.id) else {
        return Some(Why::NoRecord);
    };
    if old.verb != step.verb || old.options_hash != key.options_hash {
        return Some(Why::Options);
    }
    if old.status != StepStatus::Ok.name() {
        return Some(Why::Status(old.status.clone()));
    }
    let then = &previous.directory;
    for input in &key.inputs {
        let now = resolve(&input.path);
        let same = old.inputs.iter().any(|o| {
            input.print.is_some() && o.print == input.print && resolve_in(then, &o.path) == now
        });
        if !same {
            return Some(Why::Input(input.path.clone()));
        }
    }
    if old.inputs.len() != key.inputs.len() {
        return Some(Why::Options);
    }
    old.outputs
        .iter()
        .find(|output| !output.on_disk_in(then))
        .map(|output| Why::Output(shown(then, &output.path)))
}

/// True for a step the chain's read runs: every step but the qc step and
/// the steps it runs after the reads.
pub fn in_read(chain: &Chain, index: usize) -> bool {
    chain.steps[index].role != Role::Qc
        && !chain.qc.as_ref().is_some_and(|qc| qc.after_reads(index))
}

/// The kept BAM of transform `index` the run before wrote, when the chain
/// keeps it under the same name and it is there with its fingerprint.
fn kept_on_disk(
    chain: &Chain,
    keys: &[Key],
    previous: &ChainManifest,
    index: usize,
) -> Option<PathBuf> {
    let kept = chain.kept.iter().find(|k| k.index == index)?;
    let now = resolve(kept.path());
    let then = &previous.directory;
    previous
        .kept
        .iter()
        .find(|k| {
            k.step == keys[index].id
                && resolve_in(then, &k.file.path) == now
                && k.file.on_disk_in(then)
        })
        .map(|_| kept.path().to_path_buf())
}

/// The BAM a skipped `write` step holds of the stream transform `index`
/// leaves: the records every step of that stream reads, written by the
/// writer `--keep` uses. A `write` that runs in this run, one that writes
/// SAM text and one that writes standard output give none.
fn written_on_disk(chain: &Chain, runs: &[Option<Why>], index: usize) -> Option<PathBuf> {
    let stream = chain.steps[index].leaves.as_deref()?;
    chain
        .steps
        .iter()
        .zip(runs)
        .skip(index + 1)
        .filter(|(step, runs)| runs.is_none() && step.reads == stream && !step.stdout)
        .filter_map(|(step, _)| step.write.as_ref())
        .find(|options| options.format == bammetrics_core::stage::WriteFormat::Bam)
        .map(|options| options.path.clone())
}

/// What the run skips, what it runs and why, and where its read starts.
pub fn decide(chain: &Chain, keys: &[Key]) -> Resume {
    // A chain that writes no manifest fingerprints nothing.
    let snapshot = chain.manifest.as_ref().map(|_| Snapshot::take(chain));
    let (true, Some(path), Some(snap)) = (chain.resume, &chain.manifest, &snapshot) else {
        return Resume::everything(chain, snapshot, None, None);
    };
    if let Some(uid) = foreign_owner(path) {
        let note = format!(
            "{} is owned by user {uid}, not by this user: every step runs and no file is removed",
            path.display()
        );
        return Resume::everything(chain, snapshot, None, Some(note));
    }
    let previous = match v2::read_any(path) {
        Ok(Some(AnyManifest::Chain(manifest))) => *manifest,
        Ok(Some(AnyManifest::Qc(_))) => {
            let note = format!(
                "{} holds a qc run ({}): every step runs",
                path.display(),
                crate::manifest::SCHEMA
            );
            return Resume::everything(chain, snapshot, None, Some(note));
        }
        Ok(None) => {
            let note = format!("no manifest at {}: every step runs", path.display());
            return Resume::everything(chain, snapshot, None, Some(note));
        }
        Err(why) => {
            let note = format!("{why}: every step runs");
            return Resume::everything(chain, snapshot, None, Some(note));
        }
    };
    if let Some(why) = whole_run(snap, &previous, stringency_name(chain)) {
        let note = format!("{why}: every step runs");
        return Resume::everything(chain, snapshot, Some(previous), Some(note));
    }
    let mut runs: Vec<Option<Why>> = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        let why = stale(chain, index, key, &previous).or_else(|| {
            chain.steps[..index]
                .iter()
                .zip(&runs)
                .find(|(step, runs)| step.hands_on() && runs.is_some())
                .map(|(step, _)| Why::Upstream(step.number))
        });
        runs.push(why);
    }
    let reading: Vec<usize> = (0..runs.len())
        .filter(|&i| runs[i].is_some() && in_read(chain, i))
        .collect();
    let mut start = None;
    if let (Some(&first), Some(&last)) = (reading.first(), reading.last()) {
        start = (0..first)
            .rev()
            .filter(|&t| chain.steps[t].hands_on())
            .find_map(|t| {
                kept_on_disk(chain, keys, &previous, t)
                    .or_else(|| written_on_disk(chain, &runs, t))
                    .map(|path| (t, path))
            });
        let from = start.as_ref().map_or(0, |(t, _)| t + 1);
        for t in from..last {
            if chain.steps[t].hands_on() && runs[t].is_none() {
                let fed = reading.iter().find(|&&r| r > t).copied().unwrap_or(last);
                runs[t] = Some(Why::Feeds(chain.steps[fed].number));
            }
        }
    }
    Resume {
        snapshot,
        previous: Some(previous),
        note: None,
        runs,
        start,
    }
}

/// The owner of `path`, or of the file a link at `path` points to, when
/// that owner is another user. The paths a resumed run removes come from
/// the manifest, so a manifest another user can write is never trusted.
fn foreign_owner(path: &Path) -> Option<u32> {
    use std::os::unix::fs::MetadataExt as _;
    // SAFETY: geteuid has no preconditions and cannot fail.
    let me = unsafe { libc::geteuid() };
    [std::fs::symlink_metadata(path), std::fs::metadata(path)]
        .into_iter()
        .filter_map(Result::ok)
        .map(|meta| meta.uid())
        .find(|&uid| uid != me)
}

/// The lines of the plan text about `--resume`.
pub fn plan_lines(chain: &Chain, resume: &Resume) -> Vec<String> {
    let Some(path) = chain.manifest.as_ref().filter(|_| chain.resume) else {
        return Vec::new();
    };
    let mut lines = vec![format!("resume     from {}", path.display())];
    if let Some(note) = &resume.note {
        lines.push(format!("resume     {note}"));
        return lines;
    }
    for (index, step) in chain.steps.iter().enumerate() {
        let what = match (&resume.runs[index], &resume.start) {
            (None, Some((start, kept))) if *start == index => {
                let is_kept = chain
                    .kept
                    .iter()
                    .any(|k| k.index == index && k.path() == kept);
                if is_kept {
                    format!(
                        "skipped; its kept stream {} starts the read",
                        kept.display()
                    )
                } else {
                    format!(
                        "skipped; the BAM a write step wrote of its stream, {}, starts the read",
                        kept.display()
                    )
                }
            }
            (None, _) => "skipped: its record and files match".to_owned(),
            (Some(why), _) => why.text(),
        };
        lines.push(format!(
            "resume     step {} {}: {what}",
            step.number, step.verb
        ));
    }
    if resume.runs.iter().all(Option::is_none)
        || (0..chain.steps.len()).all(|i| resume.skipped(i) || !in_read(chain, i))
    {
        lines.push("resume     no step reads the stream: the source is not read".to_owned());
    }
    lines
}

/// The record of a step that ran in this run, or did not finish.
fn record(
    chain: &Chain,
    index: usize,
    key: &Key,
    outcome: &StepOutcome,
    started: &str,
    millis: u64,
) -> StepRecord {
    let step = &chain.steps[index];
    StepRecord {
        id: key.id.clone(),
        number: step.number,
        verb: step.verb.to_owned(),
        tool: step.tool.map(|tool| tool.name.to_owned()),
        role: step.role.name().to_owned(),
        reads: step.reads.clone(),
        leaves: step.leaves.clone(),
        upstream: key.upstream.clone(),
        argv: key.argv.clone(),
        options_hash: key.options_hash.clone(),
        inputs: key.inputs.clone(),
        outputs: outcome
            .files
            .iter()
            .map(|file| match &file.print {
                Some(print) => FileRecord {
                    path: file.path.clone(),
                    kind: Some(file.kind.to_owned()),
                    print: Some(print.clone()),
                },
                None => FileRecord::output(&file.path, file.kind),
            })
            .collect(),
        status: outcome.status.name().to_owned(),
        error: outcome.error.clone(),
        skipped: false,
        started: started.to_owned(),
        millis,
    }
}

/// The run's manifest.
pub struct Run<'a> {
    pub keys: &'a [Key],
    pub resume: &'a Resume,
    pub outcomes: &'a [StepOutcome],
    pub sample: String,
    pub started: String,
    pub millis: u64,
    pub code: u8,
}

/// The manifest of `run`: a skipped step keeps its record of the run before
/// (with `skipped` true and its number now), every other step gets the
/// record of this run.
pub fn manifest(chain: &Chain, run: &Run) -> ChainManifest {
    let snapshot = run
        .resume
        .snapshot
        .clone()
        .unwrap_or_else(|| Snapshot::take(chain));
    let mut steps = Vec::with_capacity(chain.steps.len());
    for (index, (key, outcome)) in run.keys.iter().zip(run.outcomes).enumerate() {
        let kept = run
            .resume
            .previous
            .as_ref()
            .filter(|_| run.resume.skipped(index))
            .and_then(|previous| Some((previous.step(&key.id)?, &previous.directory)));
        steps.push(match kept {
            Some((old, then)) => {
                let mut old = StepRecord {
                    number: chain.steps[index].number,
                    skipped: true,
                    ..old.clone()
                };
                // The paths of a record of another working directory are
                // kept as that directory names them.
                for file in old.inputs.iter_mut().chain(old.outputs.iter_mut()) {
                    file.path = rebase(then, &snapshot.directory, &file.path);
                }
                old
            }
            None => record(chain, index, key, outcome, &run.started, run.millis),
        });
    }
    let kept = chain
        .kept
        .iter()
        .filter_map(|kept| {
            let id = &run.keys[kept.index].id;
            let file = steps[kept.index]
                .outputs
                .iter()
                .find(|o| o.kind.as_deref() == Some(super::keep::KIND))?;
            Some(Kept {
                step: id.clone(),
                stream: kept.stream.clone(),
                file: file.clone(),
            })
        })
        .collect();
    ChainManifest {
        tool_version: env!("CARGO_PKG_VERSION").to_owned(),
        sample: run.sample.clone(),
        directory: snapshot.directory,
        source: snapshot.source,
        reference: snapshot.reference,
        prefix: chain.prefix.clone(),
        stringency: stringency_name(chain).to_owned(),
        started: run.started.clone(),
        millis: run.millis,
        exit: run.code,
        steps,
        kept,
    }
}

/// The outcome of a step the run skips: its files as the run before wrote
/// them.
pub fn skipped_outcome(number: usize, key: &Key, resume: &Resume) -> StepOutcome {
    let here = resume
        .snapshot
        .as_ref()
        .map_or_else(directory, |s| s.directory.clone());
    let files = resume
        .previous
        .as_ref()
        .and_then(|previous| Some((previous.step(&key.id)?, &previous.directory)))
        .map(|(old, then)| {
            old.outputs
                .iter()
                .map(|o| super::run::Written {
                    path: rebase(then, &here, &o.path),
                    kind: known_kind(o.kind.as_deref().unwrap_or_default()),
                    print: o.print.clone(),
                })
                .collect()
        })
        .unwrap_or_default();
    StepOutcome {
        number,
        status: StepStatus::Skipped,
        error: None,
        files,
        counters: None,
        notes: Vec::new(),
    }
}

/// The kinds a run writes, as the `&'static str` of [`super::run::Written`].
const KINDS: [&str; 25] = [
    "metrics",
    "gaps",
    "file",
    "bam",
    "sam",
    "index",
    "md5",
    "kept",
    "intervals",
    "verdict",
    "report",
    "pdf",
    "manifest",
    "json",
    "jsonl",
    "tsv",
    "csv",
    "parquet",
    "xlsx",
    "docx",
    "sql",
    "markdown",
    "junit",
    "multiqc",
    "human",
];

/// [`KINDS`] entry `kind`, or `file` for a kind no run writes.
fn known_kind(kind: &str) -> &'static str {
    KINDS
        .iter()
        .find(|k| **k == kind)
        .copied()
        .unwrap_or("file")
}

/// True when a manifest may name `path` with `kind`: a kind a chain
/// writes, and for a BAM, SAM, index or MD5 the extension that kind has.
/// A record that calls some other file a `bam` is not trusted.
fn kind_fits(kind: &str, path: &Path) -> bool {
    let extension = path.extension().and_then(|e| e.to_str());
    let wanted = match kind {
        "bam" | "kept" => Some("bam"),
        "sam" => Some("sam"),
        "index" => Some("bai"),
        "md5" => Some("md5"),
        _ => None,
    };
    KINDS.contains(&kind) && wanted.is_none_or(|w| extension == Some(w))
}

/// A path of a run in directory `then` as a run in `here` names it: as
/// recorded when the two are one directory, else relative to `here` when
/// the file is under it, else absolute.
fn rebase(then: &Path, here: &Path, path: &Path) -> PathBuf {
    if then == here {
        return path.to_path_buf();
    }
    let absolute = resolve_in(then, path);
    absolute
        .strip_prefix(here)
        .map_or_else(|_| absolute.clone(), Path::to_path_buf)
}

/// [`rebase`] into this run's working directory.
fn shown(then: &Path, path: &Path) -> PathBuf {
    rebase(then, &directory(), path)
}

/// What a resumed run did with the files of the run before that it neither
/// wrote nor kept.
#[derive(Debug, Default)]
pub struct Stale {
    pub removed: Vec<PathBuf>,
    /// Each file left in place, with why.
    pub left: Vec<(PathBuf, &'static str)>,
}

impl Stale {
    /// One line per file, for standard error.
    pub fn notes(&self) -> Vec<String> {
        let removed = self.removed.iter().map(|path| {
            format!(
                "--resume: removed {}: the run before wrote it and this run does not",
                path.display()
            )
        });
        let left = self
            .left
            .iter()
            .map(|(path, why)| format!("--resume: left in place: {} ({why})", path.display()));
        removed.chain(left).collect()
    }
}

/// The `.bai` and `.md5` names htsjdk gives the companions of BAM or SAM
/// file `bam`.
fn companions(bam: &Path) -> [PathBuf; 3] {
    [
        bam.with_extension("bai"),
        crate::qc::with_suffix(bam, ".bai"),
        crate::qc::with_suffix(bam, ".md5"),
    ]
}

/// After a resumed run that finished: remove the files the run before
/// wrote that this run neither wrote nor kept. A file goes only when the
/// manifest of the run before lists it, it is still a regular file (no
/// link is followed) with the fingerprint recorded, the chain does not
/// read it, and either it sits in the directory of `PREFIX` with a name
/// that starts with `PREFIX`'s and its recorded kind is one this run
/// writes, or it is an `index` or `md5` record naming the `.bai` or `.md5`
/// of a BAM the same step wrote in this run. [`decide`] drops a manifest
/// another user owns, so this never runs on one. The files of a `write` or split step this
/// run does not have, a file elsewhere and a file changed since are left in
/// place and named.
pub fn remove_stale(chain: &Chain, previous: &ChainManifest, now: &ChainManifest) -> Stale {
    let here = directory();
    let written: HashSet<PathBuf> = now
        .steps
        .iter()
        .flat_map(|s| &s.outputs)
        .map(|o| resolve(&o.path))
        .collect();
    // The companion names of each BAM a step wrote in this run, by step id.
    let mut beside: HashMap<&str, HashSet<PathBuf>> = HashMap::new();
    for step in now.steps.iter().filter(|s| !s.skipped) {
        let names = step
            .outputs
            .iter()
            .filter(|o| o.path.extension().is_some_and(|e| e == "bam" || e == "sam"))
            .flat_map(|o| companions(&resolve(&o.path)));
        beside.entry(step.id.as_str()).or_default().extend(names);
    }
    let prefix = resolve(&chain.prefix);
    let prefix_dir = prefix.parent().map(Path::to_path_buf);
    let starts = prefix
        .file_name()
        .map(|name| format!("{}.", name.to_string_lossy()));
    let in_prefix = |path: &Path| {
        let name = path.file_name().map(|n| n.to_string_lossy());
        path.parent() == prefix_dir.as_deref()
            && starts
                .as_deref()
                .zip(name)
                .is_some_and(|(starts, name)| name.starts_with(starts))
    };
    let then = &previous.directory;
    let mut stale = Stale::default();
    let mut seen = HashSet::new();
    for step in &previous.steps {
        for output in &step.outputs {
            let path = rebase(then, &here, &output.path);
            let resolved = resolve(&path);
            let plain = std::fs::symlink_metadata(&path).is_ok_and(|m| m.file_type().is_file());
            if !plain
                || written.contains(&resolved)
                || chain.inputs.owner(&path).is_some()
                || !seen.insert(resolved.clone())
            {
                continue;
            }
            let kind = output.kind.as_deref();
            let companion = matches!(kind, Some("index" | "md5"))
                && beside
                    .get(step.id.as_str())
                    .is_some_and(|names| names.contains(&resolved));
            let why = if matches!(step.role.as_str(), "write" | "split")
                && now.step(&step.id).is_none()
            {
                Some("a write or split step's file; this run has no such step")
            } else if !in_prefix(&resolved) && !companion {
                Some("outside the files of the output prefix")
            } else if !kind.is_some_and(|k| kind_fits(k, &resolved)) {
                Some("its kind is none a chain writes under that name")
            } else if !output.on_disk_in(then) {
                Some("changed since the run before wrote it")
            } else {
                None
            };
            match why {
                Some(why) => stale.left.push((path, why)),
                None => {
                    if std::fs::remove_file(&path).is_ok() {
                        stale.removed.push(path);
                    }
                }
            }
        }
    }
    stale
}

/// The temporary files a run of this prefix left when it was killed
/// (`.NAME.PID.N.HASH.tmp` in the prefix's directory, its process gone):
/// named, never removed, since another program may own a name of that
/// form.
pub fn orphans(chain: &Chain) -> Vec<PathBuf> {
    let Some(name) = chain.prefix.file_name().map(|n| n.to_string_lossy()) else {
        return Vec::new();
    };
    let dir = match chain.prefix.parent() {
        Some(dir) if !dir.as_os_str().is_empty() => dir.to_path_buf(),
        _ => PathBuf::from("."),
    };
    let Ok(entries) = std::fs::read_dir(&dir) else {
        return Vec::new();
    };
    let starts = format!(".{name}.");
    let mut found: Vec<PathBuf> = entries
        .filter_map(Result::ok)
        .filter_map(|entry| {
            let file = entry.file_name().to_string_lossy().into_owned();
            let rest = file.strip_prefix(&starts)?.strip_suffix(".tmp")?;
            // NAME.PID.N.HASH: the process id is the third field from the end.
            let pid: i32 = rest.rsplit('.').nth(2)?.parse().ok()?;
            (!alive(pid)).then(|| dir.join(&file))
        })
        .collect();
    found.sort();
    found
}

/// True while process `pid` runs (or another user's process has it).
fn alive(pid: i32) -> bool {
    if pid <= 0 {
        return false;
    }
    // SAFETY: signal 0 checks for the process and sends nothing.
    let sent = unsafe { libc::kill(pid, 0) };
    sent == 0 || std::io::Error::last_os_error().raw_os_error() == Some(libc::EPERM)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_reason_reads_as_a_clause() {
        let reasons = [
            Why::NoRecord,
            Why::Options,
            Why::Status("failed".to_owned()),
            Why::Input(PathBuf::from("ref.fa")),
            Why::Output(PathBuf::from("out/s.insert_size_metrics")),
            Why::Stdout,
            Why::Upstream(3),
            Why::Feeds(4),
            Why::Verdict,
        ];
        for why in reasons {
            assert!(why.text().starts_with("runs: "), "{why:?}");
        }
        assert_eq!(Why::Upstream(3).text(), "runs: step 3 in front of it runs");
    }

    #[test]
    fn a_written_kind_reads_back_as_itself() {
        for kind in ["metrics", "kept", "index", "json", "verdict", "markdown"] {
            assert_eq!(known_kind(kind), kind);
        }
        assert_eq!(known_kind("nope"), "file");
    }

    #[test]
    fn every_spelling_of_one_file_resolves_to_one_path() {
        let dir = std::env::temp_dir().join(format!("bammetrics-resolve-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("a/b")).expect("dirs");
        let dir = dir.canonicalize().expect("dir");
        std::fs::write(dir.join("a/m.bam"), "x").expect("file");
        let want = dir.join("a/m.bam");
        for spelled in [
            "a/m.bam",
            "./a/m.bam",
            "a/b/../m.bam",
            "a/./b/../../a/m.bam",
        ] {
            assert_eq!(resolve_in(&dir, Path::new(spelled)), want, "{spelled}");
        }
        assert_eq!(resolve_in(&dir.join("a/b"), Path::new("../m.bam")), want);
        assert_eq!(resolve_in(Path::new("/x"), &want), want);
        // A path that does not exist yet: the part that does is resolved.
        assert_eq!(
            resolve_in(&dir, Path::new("a/b/../new/s.txt")),
            dir.join("a/new/s.txt")
        );
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(dir.join("a"), dir.join("link")).expect("link");
            assert_eq!(resolve_in(&dir, Path::new("link/m.bam")), want);
        }
        assert_eq!(resolve_in(&dir, Path::new("-")), PathBuf::from("-"));
        // The paths of a run in a/b, read by a run in `dir`.
        assert_eq!(
            rebase(&dir.join("a/b"), &dir, Path::new("../m.bam")),
            PathBuf::from("a/m.bam")
        );
        assert_eq!(rebase(&dir, &dir, Path::new("./x")), PathBuf::from("./x"));
        assert_eq!(
            rebase(&dir.join("a"), Path::new("/elsewhere"), Path::new("m.bam")),
            want
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_kind_must_be_one_a_chain_writes_with_its_extension() {
        assert!(kind_fits("metrics", Path::new("o/s.insert_size_metrics")));
        assert!(kind_fits("kept", Path::new("o/s.marked.bam")));
        assert!(kind_fits("index", Path::new("o/s.marked.bai")));
        assert!(kind_fits("md5", Path::new("o/s.marked.bam.md5")));
        assert!(!kind_fits("bam", Path::new("o/s.insert_size_metrics")));
        assert!(!kind_fits("index", Path::new("o/s.marked.bam")));
        assert!(!kind_fits("kept", Path::new("o/s.vcf")));
        assert!(!kind_fits("anything", Path::new("o/s.txt")));
        assert!(!kind_fits("", Path::new("o/s.txt")));
    }

    #[test]
    fn a_file_of_another_user_is_foreign() {
        let dir = std::env::temp_dir().join(format!("bammetrics-owner-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let mine = dir.join("m.manifest.json");
        std::fs::write(&mine, "{}").expect("write");
        assert_eq!(foreign_owner(&mine), None);
        // SAFETY: geteuid has no preconditions.
        if unsafe { libc::geteuid() } != 0 {
            assert_eq!(foreign_owner(Path::new("/")), Some(0));
            let link = dir.join("link.manifest.json");
            std::os::unix::fs::symlink("/", &link).expect("link");
            assert_eq!(foreign_owner(&link), Some(0));
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_stale_note_names_each_file() {
        let stale = Stale {
            removed: vec![PathBuf::from("out/s.sorted.bai")],
            left: vec![(
                PathBuf::from("deliver/final.bam"),
                "a write or split step's file",
            )],
        };
        assert_eq!(
            stale.notes(),
            [
                "--resume: removed out/s.sorted.bai: the run before wrote it and this run does not",
                "--resume: left in place: deliver/final.bam (a write or split step's file)",
            ]
        );
    }
}

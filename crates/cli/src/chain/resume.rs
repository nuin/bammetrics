//! `chain --resume` and the chain's manifest (docs/chaining.md section 22).
//!
//! Each step has a key: an id (a transform's stream, else `VERB@STREAM`),
//! the ids of the transforms in front of it, its canonical argv, the hash
//! of that argv with the tool version, the formats, the sample and the
//! source's validation stringency, and the fingerprints of the files it
//! reads by name. A step is skipped when the record of the run before
//! under its id has the same verb and hash, ended `ok`, read files with the
//! same fingerprints, and every file it wrote is still there with its
//! fingerprint, and when every transform in front of it is skipped too.
//! A step that writes standard output and the qc step always run.
//!
//! The steps that run share reads. The read starts at the last transform
//! in front of the first step that runs whose stream the run before kept
//! (`--keep`) and whose BAM is still there, else at the source; every
//! transform between that start and the last step that runs runs too, to
//! feed the steps after it.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use bammetrics_verdict::bytes_hash;

use super::prepare::{Chain, Planned, Role};
use super::run::{StepOutcome, StepStatus};
use crate::manifest::v2::{self, AnyManifest, ChainManifest, FileRecord, Kept, StepRecord};

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
        (Role::Transform, Some(stream)) => stream.clone(),
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
        let text = argv
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
        if step.role == Role::Transform {
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
    /// The manifest of the run before, when it is one of a chain.
    pub previous: Option<ChainManifest>,
    /// Why every step runs, when that is so for the whole run.
    pub note: Option<String>,
    /// One per step: `None` for a step that is skipped, else why it runs.
    pub runs: Vec<Option<Why>>,
    /// The kept stream the read starts from: the transform's index and its
    /// BAM. `None` reads the source.
    pub start: Option<(usize, PathBuf)>,
}

impl Resume {
    fn everything(chain: &Chain, previous: Option<ChainManifest>, note: Option<String>) -> Self {
        Self {
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
fn whole_run(chain: &Chain, previous: &ChainManifest) -> Option<String> {
    if previous.tool_version != env!("CARGO_PKG_VERSION") {
        return Some(format!(
            "the manifest was written by bammetrics {}",
            previous.tool_version
        ));
    }
    let source = FileRecord::input(&chain.source_path);
    if previous.source.path != source.path || source.print.is_none() {
        return Some(format!(
            "the source is not {}",
            previous.source.path.display()
        ));
    }
    if previous.source.print != source.print {
        return Some(format!("the source {} changed", source.path.display()));
    }
    let reference = chain.reference.as_deref().map(FileRecord::input);
    let same = match (&reference, &previous.reference) {
        (None, None) => true,
        (Some(new), Some(old)) => new.path == old.path && new.print == old.print,
        _ => false,
    };
    if !same {
        return Some("the reference changed".to_owned());
    }
    if previous.stringency != stringency_name(chain) {
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
    for input in &key.inputs {
        let same = old
            .inputs
            .iter()
            .any(|o| o.path == input.path && input.print.is_some() && o.print == input.print);
        if !same {
            return Some(Why::Input(input.path.clone()));
        }
    }
    if old.inputs.len() != key.inputs.len() {
        return Some(Why::Options);
    }
    old.outputs
        .iter()
        .find(|output| !output.on_disk())
        .map(|output| Why::Output(output.path.clone()))
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
    previous
        .kept
        .iter()
        .find(|k| k.step == keys[index].id && k.file.path == kept.path() && k.file.on_disk())
        .map(|k| k.file.path.clone())
}

/// What the run skips, what it runs and why, and where its read starts.
pub fn decide(chain: &Chain, keys: &[Key]) -> Resume {
    let (true, Some(path)) = (chain.resume, &chain.manifest) else {
        return Resume::everything(chain, None, None);
    };
    let previous = match v2::read_any(path) {
        Ok(Some(AnyManifest::Chain(manifest))) => *manifest,
        Ok(Some(AnyManifest::Qc(_))) => {
            let note = format!(
                "{} holds a qc run ({}): every step runs",
                path.display(),
                crate::manifest::SCHEMA
            );
            return Resume::everything(chain, None, Some(note));
        }
        Ok(None) => {
            let note = format!("no manifest at {}: every step runs", path.display());
            return Resume::everything(chain, None, Some(note));
        }
        Err(why) => {
            return Resume::everything(chain, None, Some(format!("{why}: every step runs")));
        }
    };
    if let Some(why) = whole_run(chain, &previous) {
        let note = format!("{why}: every step runs");
        return Resume::everything(chain, Some(previous), Some(note));
    }
    let mut runs: Vec<Option<Why>> = Vec::with_capacity(keys.len());
    for (index, key) in keys.iter().enumerate() {
        let why = stale(chain, index, key, &previous).or_else(|| {
            chain.steps[..index]
                .iter()
                .zip(&runs)
                .find(|(step, runs)| step.role == Role::Transform && runs.is_some())
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
            .filter(|&t| chain.steps[t].role == Role::Transform)
            .find_map(|t| kept_on_disk(chain, keys, &previous, t).map(|path| (t, path)));
        let from = start.as_ref().map_or(0, |(t, _)| t + 1);
        for t in from..last {
            if chain.steps[t].role == Role::Transform && runs[t].is_none() {
                let fed = reading.iter().find(|&&r| r > t).copied().unwrap_or(last);
                runs[t] = Some(Why::Feeds(chain.steps[fed].number));
            }
        }
    }
    Resume {
        previous: Some(previous),
        note: None,
        runs,
        start,
    }
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
                format!(
                    "skipped; its kept stream {} starts the read",
                    kept.display()
                )
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
            .map(|file| FileRecord::output(&file.path, file.kind))
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
    let mut steps = Vec::with_capacity(chain.steps.len());
    for (index, (key, outcome)) in run.keys.iter().zip(run.outcomes).enumerate() {
        let kept = run
            .resume
            .previous
            .as_ref()
            .filter(|_| run.resume.skipped(index))
            .and_then(|previous| previous.step(&key.id));
        steps.push(match kept {
            Some(old) => StepRecord {
                number: chain.steps[index].number,
                skipped: true,
                ..old.clone()
            },
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
        source: FileRecord::input(&chain.source_path),
        reference: chain.reference.as_deref().map(FileRecord::input),
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
    let files = resume
        .previous
        .as_ref()
        .and_then(|previous| previous.step(&key.id))
        .map(|old| {
            old.outputs
                .iter()
                .map(|o| super::run::Written {
                    path: o.path.clone(),
                    kind: known_kind(o.kind.as_deref().unwrap_or_default()),
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
fn known_kind(kind: &str) -> &'static str {
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
    KINDS
        .iter()
        .find(|k| **k == kind)
        .copied()
        .unwrap_or("file")
}

/// After a resumed run that finished: remove the files the run before
/// wrote that this run neither wrote nor kept. Only a file the manifest of
/// the run before lists, still a regular file (no link is followed) with
/// the fingerprint recorded, and no file the chain reads or plans to write.
pub fn remove_stale(chain: &Chain, previous: &ChainManifest, now: &ChainManifest) -> usize {
    let listed = |path: &Path| {
        now.steps
            .iter()
            .flat_map(|s| &s.outputs)
            .any(|o| o.path == path)
    };
    let mut removed = 0;
    for output in previous.steps.iter().flat_map(|s| &s.outputs) {
        let path = output.path.as_path();
        let plain = std::fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_file());
        if !plain
            || listed(path)
            || !output.on_disk()
            || chain.inputs.owner(path).is_some()
            || chain.outputs.owner(path).is_some()
        {
            continue;
        }
        if std::fs::remove_file(path).is_ok() {
            removed += 1;
        }
    }
    removed
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
}

//! `bammetrics qc`: the stages of a [`Plan`], run. See [`crate::qc`] for
//! the stage graph and [`crate::manifest`] for what `--resume` compares.

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metric::{MetricsOutput, ToolOutput};
use bammetrics_core::metrics_file::CivilTime;
use bammetrics_core::pass::{self, Limits};
use bammetrics_report::Options;
use bammetrics_verdict::{Band, Input, Status as CheckStatus, Verdict, bytes_hash, evaluate};

use crate::also;
use crate::check::{self, status_word};
use crate::error::{Error, Result};
use crate::events::{self, Reporter, print_stdout};
use crate::manifest::{Key, Manifest, Output, Status, StepRecord};
use crate::qc::{
    Origin, Plan, QcArgs, Read, Step, make_prefix_dir, manifest_path, page_path, plural, prepare,
    verdict_path,
};
use crate::report::write_atomic;

/// What the run did to the files under the prefix.
#[derive(Default)]
struct Files {
    written: usize,
    replaced: usize,
    kept: usize,
    removed: usize,
}

/// One collector after the pass: its final record and whether it ran now.
struct Done {
    record: StepRecord,
    ran: bool,
}

/// `file` rendered as the Picard face writes it for `step`'s arguments.
fn render(step: &Step, file: &MetricsOutput) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    file.write(file.headers(&step.inv.command), &mut bytes)
        .map_err(|e| Error::Input(format!("{}: {e}", step.tool.name)))?;
    Ok(bytes)
}

/// Write `output`'s files atomically, and remove a planned file the tool
/// did not write this time (Picard writes no insert size file when every
/// orientation is dropped), so the prefix holds no file of an older run.
fn write_output(
    step: &Step,
    output: &ToolOutput,
    files: &mut Files,
    out: &Reporter,
) -> Result<Vec<Output>> {
    let mut pending: Vec<(PathBuf, &MetricsOutput)> = Vec::new();
    if let (Some(main), Some(path)) = (&output.main, &step.inv.output) {
        pending.push((path.clone(), main));
    }
    pending.extend(output.extra.iter().map(|(p, f)| (p.clone(), f)));
    let mut written = Vec::new();
    for (path, file) in pending {
        let existed = path.exists();
        write_atomic(&path, render(step, file)?)?;
        files.written += 1;
        files.replaced += usize::from(existed);
        let kind = step
            .record
            .outputs
            .iter()
            .find(|o| o.path == path)
            .map_or("metrics", |o| o.kind.as_str());
        let record = Output {
            path,
            tool: step.tool.name.to_owned(),
            kind: kind.to_owned(),
        };
        out.file(&step.record.id, &record, true);
        written.push(record);
    }
    remove_stale(step, &written, files);
    Ok(written)
}

/// Remove the files `step` plans to write and did not.
fn remove_stale(step: &Step, written: &[Output], files: &mut Files) {
    for planned in &step.record.outputs {
        if !written.iter().any(|w| w.path == planned.path) && fs::remove_file(&planned.path).is_ok()
        {
            files.removed += 1;
        }
    }
}

/// The manifest kinds of the files a step writes: [`remove_dropped`]
/// removes only these, never the verdict, page or manifest.
const STEP_KINDS: &[&str] = &["metrics", "intervals", "gaps"];

/// Remove the step files the earlier run's manifest lists that no step of
/// this run wrote or kept (a collector the new options drop, the interval
/// list of a `--panel` BED no longer given), so the
/// prefix holds what this run's manifest lists and `report DIR` agrees
/// with it. The manifest is a file on disk, so its paths are not trusted:
/// only [`owned_by_prefix`] paths are removed, and a symlink is never followed.
fn remove_dropped(plan: &Plan, done: &[Done], files: &mut Files) {
    let Ok(Some(old)) = &plan.previous else {
        return;
    };
    let current = |path: &Path| {
        done.iter()
            .any(|d| d.record.outputs.iter().any(|o| o.path == path))
    };
    for output in old.steps.iter().flat_map(|s| &s.outputs) {
        if STEP_KINDS.contains(&output.kind.as_str())
            && owned_by_prefix(&plan.prefix, &output.path)
            && !current(&output.path)
            && fs::symlink_metadata(&output.path).is_ok_and(|m| m.file_type().is_file())
            && fs::remove_file(&output.path).is_ok()
        {
            files.removed += 1;
        }
    }
}

/// True when `path` is `<prefix>.<suffix>`: the same directory as `prefix`,
/// a file name of the prefix's file name, a dot and a suffix without a path
/// separator, and no `..` component anywhere.
fn owned_by_prefix(prefix: &Path, path: &Path) -> bool {
    use std::path::Component;
    let plain = |p: &Path| !p.components().any(|c| matches!(c, Component::ParentDir));
    if !plain(prefix) || !plain(path) || path.parent() != prefix.parent() {
        return false;
    }
    let (Some(stem), Some(name)) = (
        prefix.file_name().and_then(|n| n.to_str()),
        path.file_name().and_then(|n| n.to_str()),
    ) else {
        return false;
    };
    name.strip_prefix(stem)
        .and_then(|rest| rest.strip_prefix('.'))
        .is_some_and(|suffix| !suffix.is_empty())
}

/// The record of a step that ran: its files, or its error.
fn finished(
    step: &Step,
    result: bammetrics_core::Result<ToolOutput>,
    files: &mut Files,
    out: &Reporter,
) -> Result<Done> {
    let mut record = step.record.clone();
    match result {
        Ok(output) => record.outputs = write_output(step, &output, files, out)?,
        Err(e) => {
            out.error(&record.id, &e.to_string());
            remove_stale(step, &[], files);
            record.outputs.clear();
            record.status = Status::Failed;
            record.error = Some(e.to_string());
        }
    }
    Ok(Done { record, ran: true })
}

/// The earlier run's record of a step `--resume` skips, its files counted
/// as kept; `None` for a step that runs.
fn kept(
    plan: &Plan,
    step: &Step,
    args: &QcArgs,
    files: &mut Files,
    out: &Reporter,
) -> Option<Done> {
    if !step.skipped(args.resume) {
        return None;
    }
    let old = plan
        .previous
        .as_ref()
        .ok()
        .and_then(Option::as_ref)
        .and_then(|m| m.step(&step.record.id))
        .expect("a resumable step has a record");
    for output in &old.outputs {
        out.file(&old.id, output, false);
    }
    files.kept += old.outputs.len();
    Some(Done {
        record: old.clone(),
        ran: false,
    })
}

/// Prepare, the pass and the own reads: every step not skipped by
/// `--resume`, prepare steps first.
fn collect(plan: &Plan, args: &QcArgs, files: &mut Files, out: &mut Reporter) -> Result<Vec<Done>> {
    let mut done = Vec::with_capacity(plan.prepare.len() + plan.steps.len());
    for step in &plan.prepare {
        match kept(plan, step, args, files, out) {
            Some(old) => done.push(old),
            None => done.push(finished(step, (step.tool.run)(&step.inv), files, out)?),
        }
    }
    let running: Vec<&Step> = plan
        .steps
        .iter()
        .filter(|s| s.read == Read::Pass && !s.skipped(args.resume))
        .collect();
    let mut results = if running.is_empty() {
        Vec::new()
    } else {
        let mut on_progress = |p| out.progress(p);
        let results = pass::run_isolated(
            &plan.bam,
            args.reference.as_deref(),
            args.threads,
            Limits::default(),
            Some(&mut on_progress),
            |header| {
                running
                    .iter()
                    .map(|s| (s.tool.collector.expect("a pass step"))(&s.inv, header))
                    .collect()
            },
        );
        out.end_progress();
        results?
    }
    .into_iter();
    for step in &plan.steps {
        if let Some(old) = kept(plan, step, args, files, out) {
            done.push(old);
            continue;
        }
        let result = match step.read {
            Read::Pass => results.next().expect("one result per pass step"),
            Read::Own(_) | Read::Header => (step.tool.run)(&step.inv),
        };
        done.push(finished(step, result, files, out)?);
    }
    Ok(done)
}

/// A step `qc` runs itself (the verdict, the page), keyed on the run.
fn own_step(
    plan: &Plan,
    id: &str,
    stage: &str,
    argv: Vec<String>,
    outputs: Vec<Output>,
) -> StepRecord {
    StepRecord {
        id: id.to_owned(),
        stage: stage.to_owned(),
        tool: "bammetrics".to_owned(),
        key: Key {
            argv_hash: bytes_hash(argv.join("\0").into_bytes()),
            tool_version: env!("CARGO_PKG_VERSION").to_owned(),
            preset_hash: plan.preset.hash.clone(),
            input: plan.bam.clone(),
            fingerprint: plan.fingerprint.clone(),
            files: Vec::new(),
        },
        argv,
        outputs,
        status: Status::Ok,
        error: None,
    }
}

fn manifest(plan: &Plan, steps: Vec<StepRecord>) -> Manifest {
    Manifest {
        tool_version: env!("CARGO_PKG_VERSION").to_owned(),
        sample: plan.sample.clone(),
        bam: plan.bam.clone(),
        prefix: plan.prefix.clone(),
        preset_name: plan.preset.name.clone(),
        preset_hash: plan.preset.hash.clone(),
        steps,
    }
}

/// An earlier run's files are replaced only on request. A manifest that
/// cannot be read (`plan.previous` is `Err`) is reported with its reason.
fn guard_previous(plan: &Plan, args: &QcArgs) -> Result<()> {
    if matches!(plan.previous, Ok(None)) || args.force {
        return Ok(());
    }
    match (&plan.previous, args.resume) {
        (Err(why), false) => Err(Error::Usage(format!(
            "{why}: --resume or --force runs every step"
        ))),
        (Err(why), true) => {
            eprintln!("bammetrics qc: {why}; running every step");
            Ok(())
        }
        (_, false) => Err(Error::Usage(format!(
            "{} holds an earlier run: --resume skips its finished steps, --force reruns them",
            manifest_path(&plan.prefix).display()
        ))),
        (_, true) => Ok(()),
    }
}

/// The metrics files the run holds, in path order, the order `report`
/// discovers the same files in.
fn metrics_paths(done: &[Done]) -> Vec<&Path> {
    let mut paths: Vec<&Path> = done
        .iter()
        .filter(|d| d.record.status == Status::Ok)
        .flat_map(|d| &d.record.outputs)
        .filter(|o| o.kind == "metrics")
        .map(|o| o.path.as_path())
        .collect();
    paths.sort();
    paths
}

/// Post and report: the verdict over the run's files, then verdict.json and
/// the page, as `bammetrics report` writes them for the same files.
fn post_and_report(
    plan: &Plan,
    args: &QcArgs,
    done: &[Done],
    files: &mut Files,
    out: &Reporter,
) -> Result<(Verdict, Vec<Input>, Vec<StepRecord>)> {
    let paths = metrics_paths(done);
    let inputs = paths
        .iter()
        .map(|p| Input::read(p))
        .collect::<bammetrics_verdict::Result<Vec<_>>>()?;
    let verdict = evaluate(&plan.preset, &inputs, Some(&plan.sample));
    let mut argv = vec![format!("--preset={}", args.preset)];
    argv.extend(paths.iter().map(|p| p.display().to_string()));
    let post = own_step(plan, "verdict", "post", argv, Vec::new());

    let mut written = vec![(verdict_path(&plan.prefix), verdict.to_json(), "verdict")];
    if args.write_page {
        let date = args.date.clone().unwrap_or_else(|| CivilTime::now().iso());
        let page = bammetrics_report::render(
            &verdict,
            &inputs,
            &Options {
                all_plots: args.all_plots,
                date: &date,
            },
        );
        written.push((page_path(&plan.prefix), page, "report"));
    } else {
        remove_stale_page(plan, files);
    }
    let mut outputs = Vec::new();
    for (path, text, kind) in written {
        files.replaced += usize::from(path.exists());
        write_atomic(&path, text)?;
        files.written += 1;
        let output = Output {
            path,
            tool: "bammetrics".to_owned(),
            kind: kind.to_owned(),
        };
        out.file("report", &output, true);
        outputs.push(output);
    }
    let mut report_argv = vec![format!("--all-plots={}", args.all_plots)];
    report_argv.extend(args.date.iter().map(|d| format!("--date={d}")));
    let report = own_step(plan, "report", "report", report_argv, outputs);
    Ok((verdict, inputs, vec![post, report]))
}

/// `--also`: write the files of [`also::texts`] and return their manifest
/// step, `also`; `None` without the option. An `--also` file the earlier
/// run's manifest lists and this run does not write is removed, under the
/// rules of [`remove_dropped`]: a path owned by the prefix, a regular
/// file, never through a symlink.
fn write_also(
    plan: &Plan,
    args: &QcArgs,
    verdict: &Verdict,
    done: &[Done],
    files: &mut Files,
    out: &Reporter,
) -> Result<Option<StepRecord>> {
    let texts = also::texts(
        &plan.prefix,
        &metrics_paths(done),
        &args.also,
        &plan.sample,
        verdict,
        args.fail_on,
    )?;
    let mut outputs = Vec::new();
    for (target, text) in texts {
        files.replaced += usize::from(target.path.exists());
        write_atomic(&target.path, text)?;
        files.written += 1;
        let output = Output {
            path: target.path,
            tool: "bammetrics".to_owned(),
            kind: target.also.name().to_owned(),
        };
        out.file("also", &output, true);
        outputs.push(output);
    }
    if let Ok(Some(old)) = &plan.previous {
        let stale = old
            .steps
            .iter()
            .filter(|s| s.id == "also")
            .flat_map(|s| &s.outputs)
            .filter(|o| also::is_kind(&o.kind))
            .filter(|o| !outputs.iter().any(|new| new.path == o.path));
        for output in stale {
            if owned_by_prefix(&plan.prefix, &output.path)
                && fs::symlink_metadata(&output.path).is_ok_and(|m| m.file_type().is_file())
                && fs::remove_file(&output.path).is_ok()
            {
                files.removed += 1;
            }
        }
    }
    if outputs.is_empty() {
        return Ok(None);
    }
    let names: Vec<&str> = args.also.iter().map(|a| a.name()).collect();
    let argv = vec![format!("--also={}", names.join(","))];
    Ok(Some(own_step(plan, "also", "report", argv, outputs)))
}

/// Without the page (Python `qc(report=False)`), an earlier run's page
/// would show an older verdict: remove it. Only the path this run computes
/// from its prefix, only when the earlier run's manifest lists it as its
/// `report` file, and only a regular file (a symlink is never followed).
/// A page no run of this prefix wrote stays.
fn remove_stale_page(plan: &Plan, files: &mut Files) {
    let Ok(Some(old)) = &plan.previous else {
        return;
    };
    let page = page_path(&plan.prefix);
    let written_before = old
        .steps
        .iter()
        .flat_map(|s| &s.outputs)
        .any(|o| o.kind == "report" && o.path == page);
    if written_before
        && fs::symlink_metadata(&page).is_ok_and(|m| m.file_type().is_file())
        && fs::remove_file(&page).is_ok()
    {
        files.removed += 1;
    }
}

/// `qc`: the exit code, after a fatal error too (printed as `bammetrics qc:
/// ...`, and as `error` and `done` events with `--format jsonl`).
pub fn run(args: &QcArgs) -> u8 {
    qc(args).unwrap_or_else(|e| events::fatal("qc", args.format, &e))
}

/// A finished `qc` run, as [`execute`] returns it.
pub struct Outcome {
    pub plan: Plan,
    pub verdict: Verdict,
    /// The metrics files the verdict was evaluated over, parsed, in path order.
    pub inputs: Vec<Input>,
    /// The manifest's step records: prepare steps and collectors, then the
    /// verdict and the report.
    pub records: Vec<StepRecord>,
    /// The exit code `bammetrics qc` gives the same run: 0, 1 when a step
    /// failed, 10 or 11 under `--fail-on`.
    pub code: u8,
}

/// The stages of `qc`, writing the files the command writes, with no
/// progress line, no event, no summary and no step error on stderr (the
/// Python `qc`).
pub fn execute(args: &QcArgs) -> Result<Outcome> {
    stages(args, true).map(|s| s.outcome)
}

/// The Python `collect(metrics=[...])`: the collectors of `tools` on the
/// prepare stage and the one pass `qc` runs, their files written as `qc`
/// writes them; no manifest, verdict or page. Returns the plan and one
/// record per step, a failed step with its error.
pub fn collect_tools(args: &QcArgs, tools: &[&str]) -> Result<(Plan, Vec<StepRecord>)> {
    let plan = crate::qc::prepare_tools(args, tools)?;
    let mut out = Reporter::silent();
    make_prefix_dir(&plan.prefix)?;
    let mut files = Files::default();
    let done = collect(&plan, args, &mut files, &mut out)?;
    let records = done.into_iter().map(|d| d.record).collect();
    Ok((plan, records))
}

/// A run's outcome plus what the summary needs.
struct Stages {
    outcome: Outcome,
    done: Vec<Done>,
    files: Files,
    out: Reporter,
}

fn stages(args: &QcArgs, silent: bool) -> Result<Stages> {
    let plan = prepare(args)?;
    guard_previous(&plan, args)?;
    let mut out = if silent {
        Reporter::silent()
    } else {
        Reporter::new(args.format, args.quiet, plan.fingerprint.size)
    };
    out.start(&plan, args.resume);
    make_prefix_dir(&plan.prefix)?;
    let mut files = Files::default();
    let done = collect(&plan, args, &mut files, &mut out)?;
    remove_dropped(&plan, &done, &mut files);
    let records: Vec<StepRecord> = done.iter().map(|d| d.record.clone()).collect();
    write_atomic(
        &manifest_path(&plan.prefix),
        manifest(&plan, records.clone()).to_json(),
    )?;
    let (verdict, inputs, own) = post_and_report(&plan, args, &done, &mut files, &out)?;
    let mut all = records;
    all.extend(own);
    all.extend(write_also(&plan, args, &verdict, &done, &mut files, &out)?);
    write_atomic(
        &manifest_path(&plan.prefix),
        manifest(&plan, all.clone()).to_json(),
    )?;
    files.written += 1;
    let code = if done.iter().any(|d| d.record.status == Status::Failed) {
        1
    } else {
        args.fail_on.map_or(0, |gate| {
            check::exit_code(std::slice::from_ref(&verdict), gate)
        })
    };
    Ok(Stages {
        outcome: Outcome {
            plan,
            verdict,
            inputs,
            records: all,
            code,
        },
        done,
        files,
        out,
    })
}

fn qc(args: &QcArgs) -> Result<u8> {
    let Stages {
        outcome,
        done,
        files,
        out,
    } = stages(args, false)?;
    let Outcome {
        plan,
        verdict,
        code,
        ..
    } = &outcome;
    let failed: Vec<&StepRecord> = done
        .iter()
        .map(|d| &d.record)
        .filter(|r| r.status == Status::Failed)
        .collect();
    out.checks(verdict);
    out.done(verdict, *code, &page_path(&plan.prefix));
    if args.format == crate::qc::Format::Text {
        let reads = usize::from(done.iter().any(|d| d.ran && d.record.stage == "pass"))
            + done
                .iter()
                .filter(|d| d.ran && d.record.stage == "own-read")
                .count();
        print_stdout(&summary(&Summary {
            plan,
            verdict,
            failed: &failed,
            files: &files,
            reads,
            records: out.records(),
            code: *code,
            gate: args.fail_on.is_some(),
        }));
    }
    Ok(*code)
}

/// The most lines the summary prints.
pub const SUMMARY_LINES: usize = 12;

struct Summary<'a> {
    plan: &'a Plan,
    verdict: &'a Verdict,
    failed: &'a [&'a StepRecord],
    files: &'a Files,
    /// BAM reads made: the pass, plus each own read.
    reads: usize,
    /// Records the shared pass read.
    records: u64,
    code: u8,
    gate: bool,
}

/// `FAIL`, then the rule and the reason, one line per issue: errors, fails,
/// warns, then not-measured rules (grouped by the class this run could not
/// write).
fn issue_lines(s: &Summary) -> Vec<String> {
    let rules_of = |id: &str| {
        s.plan
            .steps
            .iter()
            .find(|step| step.record.id == id)
            .map_or(&[][..], |step| step.rules.as_slice())
    };
    let mut lines: Vec<(String, String, String)> = s
        .failed
        .iter()
        .map(|r| {
            let n = rules_of(&r.id).len();
            let lost = match n {
                0 => String::new(),
                1 => " (1 rule not measured)".to_owned(),
                n => format!(" ({n} rules not measured)"),
            };
            (
                "ERROR".to_owned(),
                r.id.clone(),
                format!("{}{lost}", r.error.clone().unwrap_or_default()),
            )
        })
        .collect();
    let grouped = |rule: &str| {
        s.plan
            .unmet
            .iter()
            .any(|u| u.rules.iter().any(|r| r == rule))
            || s.failed
                .iter()
                .any(|f| rules_of(&f.id).iter().any(|r| r == rule))
    };
    for status in [
        CheckStatus::Fail,
        CheckStatus::Warn,
        CheckStatus::NotMeasured,
    ] {
        for check in s.verdict.checks.iter().filter(|c| c.status == status) {
            if status == CheckStatus::NotMeasured && grouped(&check.rule) {
                continue;
            }
            lines.push((
                status_word(status).to_owned(),
                check.rule.clone(),
                check.reason.clone().unwrap_or_else(|| check.label.clone()),
            ));
        }
    }
    for unmet in &s.plan.unmet {
        let n = unmet.rules.len();
        lines.push((
            "NOT MEASURED".to_owned(),
            unmet.class.clone(),
            format!("{n} rule{}, {}", if n == 1 { "" } else { "s" }, unmet.why),
        ));
    }
    let width = lines.iter().map(|l| l.1.len()).max().unwrap_or(0);
    lines
        .into_iter()
        .map(|(word, what, why)| {
            format!("{word:<12}  {what:<width$}  {why}")
                .trim_end()
                .to_owned()
        })
        .collect()
}

/// `./NA12878`: a bare prefix reads as a path in the current directory.
fn shown_prefix(prefix: &Path) -> String {
    let text = prefix.display().to_string();
    if prefix.parent().is_some_and(|p| p.as_os_str().is_empty()) {
        format!("./{text}")
    } else {
        text
    }
}

/// The summary of docs/05 section 4.10, at most [`SUMMARY_LINES`] lines: a
/// header, every error, fail, warn and not-measured line (collapsed into
/// `+N more` past the budget), the pass count, the files, and the gate
/// hint or exit reason.
fn summary(s: &Summary) -> String {
    let v = s.verdict;
    let p = &v.preset;
    let mut text = String::new();
    let records = match s.records {
        0 => String::new(),
        1 => " (1 record)".to_owned(),
        n => format!(" ({} records)", events::count(n as f64)),
    };
    let _ = writeln!(
        text,
        "{}  {}  preset {} v{} ({})  {}, {}{records}",
        v.sample,
        v.status.name().to_uppercase(),
        p.name,
        p.version,
        p.hash.get(..6).unwrap_or(&p.hash),
        plural(s.plan.steps.len(), "collector"),
        plural(s.reads, "BAM read"),
    );
    let tail = if s.code != 0 || (!s.gate && v.status != Band::Pass) {
        1
    } else {
        0
    };
    let budget = SUMMARY_LINES - 3 - tail;
    let issues = issue_lines(s);
    let page = page_path(&s.plan.prefix);
    if issues.len() > budget {
        for line in &issues[..budget - 1] {
            let _ = writeln!(text, "{line}");
        }
        let _ = writeln!(
            text,
            "+{} more, see {}",
            issues.len() - (budget - 1),
            page.display()
        );
    } else {
        for line in &issues {
            let _ = writeln!(text, "{line}");
        }
    }
    let _ = writeln!(text, "{:<12}  {} checks", "PASS", v.counts.pass);
    let mut notes = Vec::new();
    for (n, word) in [
        (s.files.kept, "kept by --resume"),
        (s.files.replaced, "replaced"),
        (s.files.removed, "stale removed"),
    ] {
        if n > 0 {
            notes.push(format!("{n} {word}"));
        }
    }
    match s.plan.prefix_origin {
        Origin::FileName => notes.push("no -o: named after the BAM".to_owned()),
        Origin::Directory => notes.push("-o is a directory: named after the BAM".to_owned()),
        Origin::Option | Origin::ReadGroup => {}
    }
    let notes = if notes.is_empty() {
        String::new()
    } else {
        format!(" ({})", notes.join(", "))
    };
    let _ = writeln!(
        text,
        "wrote {} files to {}.*{notes}, report {}",
        s.files.written,
        shown_prefix(&s.plan.prefix),
        page.display()
    );
    match s.code {
        0 if tail == 1 => text.push_str("use --fail-on fail or bammetrics check to gate\n"),
        0 => {}
        1 => {
            let prepare = s.failed.iter().filter(|r| r.stage == "prepare").count();
            let collectors = s.failed.len() - prepare;
            let what = match (collectors, prepare) {
                (n, 0) => format!("{n} collector(s)"),
                (0, p) => format!("{p} prepare step(s)"),
                (n, p) => format!("{n} collector(s) and {p} prepare step(s)"),
            };
            let _ = writeln!(text, "exit 1: {what} failed");
        }
        10 => text.push_str("exit 10: the verdict reaches --fail-on\n"),
        _ => text.push_str("exit 11: the verdict is incomplete\n"),
    }
    text
}

#[cfg(test)]
mod tests {
    use super::owned_by_prefix;
    use std::path::Path;

    #[test]
    fn owned_by_prefix_rejects_paths_outside_the_prefix() {
        let prefix = Path::new("out/s");
        assert!(owned_by_prefix(
            prefix,
            Path::new("out/s.insert_size_metrics")
        ));
        assert!(!owned_by_prefix(prefix, Path::new("out/s/../../victim")));
        assert!(!owned_by_prefix(prefix, Path::new("out/s.x/../../victim")));
        assert!(!owned_by_prefix(
            prefix,
            Path::new("out/sample2.insert_size_metrics")
        ));
        assert!(!owned_by_prefix(prefix, Path::new("/etc/s.passwd")));
        assert!(!owned_by_prefix(prefix, Path::new("out/s.")));
        assert!(!owned_by_prefix(
            prefix,
            Path::new("out/other/s.insert_size_metrics")
        ));
        let abs = Path::new("/data/run/s");
        assert!(owned_by_prefix(abs, Path::new("/data/run/s.hs_metrics")));
        assert!(!owned_by_prefix(
            abs,
            Path::new("/data/run/../run/s.hs_metrics")
        ));
    }
}

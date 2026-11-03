//! `bammetrics plan`: what `qc` with the same options would do, read from
//! the [`Plan`] prepare builds. Nothing is written and no record is read:
//! the stages, the collectors, which join the single pass and why any
//! does not, which steps `--resume` would skip, and the classes no
//! collector can write.

use std::fmt::Write as _;

use bammetrics_verdict::json::Json;

use crate::error::Result;
use crate::qc::{
    Format, Origin, Plan, QcArgs, Read, Step, manifest_path, page_path, prepare, verdict_path,
};

fn origin(o: Origin, option: &str, other: &str) -> String {
    match o {
        Origin::Option => format!("from {option}"),
        Origin::ReadGroup => "from @RG SM".to_owned(),
        Origin::FileName => format!("from {other}"),
    }
}

/// `run`, `skip`, or `run` with the note that `--resume` would skip it.
fn action(step: &Step, resume: bool) -> (&'static str, &'static str) {
    match (step.resumable, resume) {
        (true, true) => ("skip", "  (unchanged since the last run)"),
        (true, false) => ("run", "  (--resume would skip it)"),
        _ => ("run", ""),
    }
}

fn rules_text(step: &Step) -> String {
    match step.rules.len() {
        0 => "asked for by option".to_owned(),
        1 => "1 rule".to_owned(),
        n => format!("{n} rules"),
    }
}

fn manifest_line(plan: &Plan, resume: bool) -> String {
    let path = manifest_path(&plan.prefix);
    match &plan.previous {
        Ok(None) => format!("{}: none yet", path.display()),
        Err(why) => format!("{why}: every step runs"),
        Ok(Some(_)) => {
            let n = plan.steps.iter().filter(|s| s.resumable).count();
            let verb = if resume { "skips" } else { "would skip" };
            format!(
                "{}: an earlier run; --resume {verb} {n} of {} collectors{}",
                path.display(),
                plan.steps.len(),
                if resume {
                    ""
                } else {
                    " (qc needs --resume or --force)"
                }
            )
        }
    }
}

fn outputs(step: &Step) -> String {
    step.record
        .outputs
        .iter()
        .map(|o| o.path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

/// The plan as text.
pub fn text(plan: &Plan, resume: bool) -> String {
    let mut out = String::new();
    let p = &plan.preset;
    let _ = writeln!(out, "qc plan for {}", plan.bam.display());
    let _ = writeln!(
        out,
        "  sample    {} ({})",
        plan.sample,
        origin(plan.sample_origin, "--sample", "the BAM name")
    );
    let _ = writeln!(
        out,
        "  prefix    {} ({})",
        plan.prefix.display(),
        origin(plan.prefix_origin, "-o", "the BAM name, no -o")
    );
    let _ = writeln!(
        out,
        "  preset    {} v{} ({}), {} rules",
        p.name,
        p.version,
        p.hash.get(..6).unwrap_or(&p.hash),
        p.rules.len()
    );
    let _ = writeln!(out, "  manifest  {}", manifest_line(plan, resume));
    let _ = writeln!(
        out,
        "prepare   inputs checked: the BAM header, -R and its .fai, every file an option names"
    );
    let shared: Vec<&Step> = plan.steps.iter().filter(|s| s.read == Read::Pass).collect();
    let reading = shared.iter().any(|s| !s.skipped(resume));
    let _ = writeln!(
        out,
        "pass      {} BAM read, {} collector{}",
        usize::from(reading),
        shared.len(),
        if shared.len() == 1 { "" } else { "s" }
    );
    let width = plan
        .steps
        .iter()
        .map(|s| s.record.id.len())
        .max()
        .unwrap_or(0);
    for step in &shared {
        let (verb, note) = action(step, resume);
        let id = &step.record.id;
        let _ = writeln!(
            out,
            "  {verb:<4}  {id:<width$}  {}  {}{note}",
            outputs(step),
            rules_text(step)
        );
    }
    for step in plan.steps.iter().filter(|s| s.read != Read::Pass) {
        let (verb, note) = action(step, resume);
        let Read::Own(why) = step.read else { continue };
        let _ = writeln!(
            out,
            "own read  {verb} {}  {}  {}: not in the single pass, {why}{note}",
            step.record.id,
            outputs(step),
            rules_text(step)
        );
    }
    for unmet in &plan.unmet {
        let _ = writeln!(
            out,
            "not run   {}: {} ({})",
            unmet.class,
            unmet.why,
            unmet.rules.join(", ")
        );
    }
    let _ = writeln!(
        out,
        "post      verdict: preset {} over the files of the pass",
        p.name
    );
    let _ = writeln!(
        out,
        "report    {}, {} (always rewritten)",
        verdict_path(&plan.prefix).display(),
        page_path(&plan.prefix).display()
    );
    out.push_str("nothing ran: this is a plan\n");
    out
}

fn strings<'a>(items: impl Iterator<Item = &'a str>) -> Json {
    Json::Array(items.map(Json::text).collect())
}

/// The plan as JSON lines: one `plan` event, a `step` event per
/// collector, an `unmet` event per class no collector writes.
pub fn jsonl(plan: &Plan, resume: bool) -> String {
    let mut out = String::new();
    let field = |k: &str, v: Json| (k.to_owned(), v);
    let head = Json::Object(vec![
        field("event", Json::text("plan")),
        field("bam", Json::text(&plan.bam.display().to_string())),
        field("sample", Json::text(&plan.sample)),
        field("prefix", Json::text(&plan.prefix.display().to_string())),
        field("preset", Json::text(&plan.preset.name)),
        field("preset_hash", Json::text(&plan.preset.hash)),
    ]);
    let _ = writeln!(out, "{}", head.render_line());
    for step in &plan.steps {
        let (verb, _) = action(step, resume);
        let read = match step.read {
            Read::Pass => Json::text("pass"),
            Read::Own(_) => Json::text("own"),
        };
        let line = Json::Object(vec![
            field("event", Json::text("step")),
            field("id", Json::text(&step.record.id)),
            field("tool", Json::text(step.tool.name)),
            field("read", read),
            field(
                "reason",
                match step.read {
                    Read::Own(why) => Json::text(why),
                    Read::Pass => Json::Null,
                },
            ),
            field("action", Json::text(verb)),
            field("resumable", Json::Bool(step.resumable)),
            field("argv", strings(step.record.argv.iter().map(String::as_str))),
            field(
                "outputs",
                Json::Array(
                    step.record
                        .outputs
                        .iter()
                        .map(|o| Json::text(&o.path.display().to_string()))
                        .collect(),
                ),
            ),
            field("rules", strings(step.rules.iter().map(String::as_str))),
        ]);
        let _ = writeln!(out, "{}", line.render_line());
    }
    for unmet in &plan.unmet {
        let line = Json::Object(vec![
            field("event", Json::text("unmet")),
            field("class", Json::text(&unmet.class)),
            field("why", Json::text(unmet.why)),
            field("rules", strings(unmet.rules.iter().map(String::as_str))),
        ]);
        let _ = writeln!(out, "{}", line.render_line());
    }
    out
}

pub fn run(args: &QcArgs) -> Result<u8> {
    let plan = prepare(args)?;
    match args.format {
        Format::Text => print!("{}", text(&plan, args.resume)),
        Format::Jsonl => print!("{}", jsonl(&plan, args.resume)),
    }
    Ok(0)
}

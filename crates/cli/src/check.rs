//! `bammetrics check`: the gate. Verdicts from directories of metrics
//! files or from `verdict.json` files, printed as text, JSON or JSON
//! lines, with the exit code of [`EXIT_CODES`].

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use bammetrics_verdict::json::Json;
use bammetrics_verdict::{Band, Preset, Status, Verdict, evaluate};
use clap::{Args, ValueEnum};

use crate::discover::SAMPLE_RULE;
use crate::error::{Error, Result};
use crate::report::{FILE_RULE, PRESET_HELP, samples};

/// The exit code table, as `--help` and docs/cli.md print it.
pub const EXIT_CODES: &str = "Exit codes:\n  \
  0   every sample passes, or none reaches --fail-on\n  \
  1   an input could not be read, or the preset does not parse\n  \
  2   usage error; nothing ran\n  \
  10  a sample reached --fail-on: its verdict is fail (or warn with --fail-on warn)\n  \
  11  a sample is incomplete (a rule's file, row or column is missing) and none reached --fail-on";

/// The band that makes `check` exit 10.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum FailOn {
    /// Exit 10 on a warn or a fail verdict.
    Warn,
    /// Exit 10 on a fail verdict only.
    Fail,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// One block per sample: the band line, then every check.
    Text,
    /// A JSON array of verdicts, each as verdict.json.
    Json,
    /// One JSON object per line: a `check` event per check, then a
    /// `verdict` event per sample.
    Jsonl,
}

fn after_help() -> String {
    format!("{EXIT_CODES}\n\n{FILE_RULE}\n\n{SAMPLE_RULE}")
}

#[derive(Args)]
#[command(after_help = after_help())]
pub struct CheckArgs {
    /// Directories of Picard metrics files (searched at any depth), or
    /// verdict.json files written by report, or a directory report wrote
    /// (its *.verdict.json files are read, no preset needed).
    #[arg(required = true, value_name = "DIR|VERDICT.json")]
    pub inputs: Vec<PathBuf>,
    #[arg(long, default_value = "generic", value_name = "NAME|FILE", help = PRESET_HELP)]
    pub preset: String,
    /// The verdict that exits 10.
    #[arg(long, value_enum, default_value_t = FailOn::Fail)]
    pub fail_on: FailOn,
    #[arg(long, value_enum, default_value_t = Format::Text)]
    pub format: Format,
}

/// The `*.verdict.json` files directly in `dir`, in name order. A directory
/// `report` wrote holds those and no metrics file.
fn verdict_files(dir: &PathBuf) -> Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)
        .map_err(|e| Error::io(dir, e))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|p| {
            p.is_file()
                && p.file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.ends_with(".verdict.json"))
        })
        .collect();
    files.sort();
    Ok(files)
}

/// The verdicts of `inputs`: a directory is discovered and evaluated with
/// `preset`, or, when it holds `*.verdict.json` files and no metrics file,
/// those are read; a file is read as verdict.json.
pub fn verdicts(inputs: &[PathBuf], preset: &str) -> Result<Vec<Verdict>> {
    let mut verdicts = Vec::new();
    let mut preset_cache = None;
    for input in inputs {
        if input.is_dir() && crate::discover::find(input)?.is_empty() {
            let files = verdict_files(input)?;
            if files.is_empty() {
                return Err(Error::Input(format!(
                    "no metrics file under {} (a file with a ## METRICS CLASS block) and no *.verdict.json",
                    input.display()
                )));
            }
            verdicts.extend(self::verdicts(&files, preset)?);
        } else if input.is_dir() {
            let preset = match &preset_cache {
                Some(p) => p,
                None => preset_cache.insert(Preset::resolve(preset)?),
            };
            for sample in samples(std::slice::from_ref(input))? {
                verdicts.push(evaluate(preset, &sample.inputs, Some(&sample.name)));
            }
        } else {
            let text = fs::read_to_string(input).map_err(|e| Error::io(input, e))?;
            let verdict = Verdict::from_json(&text)
                .map_err(|e| Error::Input(format!("{}: {e}", input.display())))?;
            verdicts.push(verdict);
        }
    }
    Ok(verdicts)
}

/// The exit code for `verdicts` under `fail_on`.
pub fn exit_code(verdicts: &[Verdict], fail_on: FailOn) -> u8 {
    let gated = |band: Band| {
        matches!(
            (fail_on, band),
            (_, Band::Fail) | (FailOn::Warn, Band::Warn)
        )
    };
    if verdicts.iter().any(|v| gated(v.status)) {
        10
    } else if verdicts.iter().any(|v| v.status == Band::Incomplete) {
        11
    } else {
        0
    }
}

pub fn status_word(status: Status) -> &'static str {
    match status {
        Status::Pass => "PASS",
        Status::Warn => "WARN",
        Status::Fail => "FAIL",
        Status::NotMeasured => "NOT MEASURED",
    }
}

/// The text form: a band line per sample, every check under it, then the
/// exit line.
pub fn text(verdicts: &[Verdict], code: u8, fail_on: FailOn) -> String {
    let mut out = String::new();
    for v in verdicts {
        let c = v.counts;
        let _ = writeln!(
            out,
            "{}  {}  {} fail, {} warn, {} pass, {} not measured  preset {} v{} ({})",
            v.sample,
            v.status.name().to_uppercase(),
            c.fail,
            c.warn,
            c.pass,
            c.not_measured,
            v.preset.name,
            v.preset.version,
            v.preset.hash.get(..6).unwrap_or(&v.preset.hash)
        );
        let rows: Vec<(&'static str, &str, String, String, &str)> = v
            .checks
            .iter()
            .map(|check| {
                // A text check shows its cell (an empty cell as `empty`,
                // as the page does), a not-measured check `?`.
                let value = match (check.value, &check.text) {
                    (Some(v), _) => bammetrics_core::metrics_file::format_double(v),
                    (None, Some(text)) if text.is_empty() => "empty".to_owned(),
                    (None, Some(text)) => text.clone(),
                    (None, None) => "?".to_owned(),
                };
                let bounds: Vec<String> = [("warn", &check.warn), ("fail", &check.fail)]
                    .into_iter()
                    .filter_map(|(name, b)| b.as_ref().map(|b| format!("{name} {b}")))
                    .collect();
                (
                    status_word(check.status),
                    check.rule.as_str(),
                    value,
                    bounds.join(", "),
                    check.reason.as_deref().unwrap_or(""),
                )
            })
            .collect();
        let column = |pick: fn(&(&str, &str, String, String, &str)) -> usize| {
            rows.iter().map(pick).max().unwrap_or(0)
        };
        let (rule_w, value_w, bounds_w) = (
            column(|r| r.1.len()),
            column(|r| r.2.len()),
            column(|r| r.3.len()),
        );
        for (status, rule, value, bounds, reason) in &rows {
            let line = format!(
                "  {status:<12}  {rule:<rule_w$}  {value:>value_w$}  {bounds:<bounds_w$}  {reason}"
            );
            let _ = writeln!(out, "{}", line.trim_end());
        }
    }
    let gate = match fail_on {
        FailOn::Warn => "warn",
        FailOn::Fail => "fail",
    };
    let why = match code {
        0 => format!("no sample reaches --fail-on {gate}"),
        10 => format!("a sample reaches --fail-on {gate}"),
        _ => "a sample is incomplete".to_owned(),
    };
    let _ = writeln!(out, "exit {code}: {why}");
    out
}

/// One `check` event line per check of `v`: the sample, then the check as
/// verdict.json holds it.
pub fn check_events(v: &Verdict) -> impl Iterator<Item = String> + '_ {
    v.checks.iter().map(|check| {
        let mut fields = vec![
            ("event".to_owned(), Json::text("check")),
            ("sample".to_owned(), Json::text(&v.sample)),
        ];
        if let Json::Object(check_fields) = check.to_json() {
            fields.extend(check_fields);
        }
        Json::Object(fields).render_line()
    })
}

/// The JSON lines form: `check` events, then one `verdict` event per sample.
pub fn jsonl(verdicts: &[Verdict]) -> String {
    let mut out = String::new();
    for v in verdicts {
        for line in check_events(v) {
            let _ = writeln!(out, "{line}");
        }
        let line = Json::Object(vec![
            ("event".to_owned(), Json::text("verdict")),
            ("sample".to_owned(), Json::text(&v.sample)),
            ("status".to_owned(), Json::text(v.status.name())),
            ("counts".to_owned(), v.counts.to_json()),
            (
                "preset".to_owned(),
                Json::Object(vec![
                    ("name".to_owned(), Json::text(&v.preset.name)),
                    ("version".to_owned(), Json::Number(v.preset.version as f64)),
                    ("hash".to_owned(), Json::text(&v.preset.hash)),
                ]),
            ),
        ]);
        let _ = writeln!(out, "{}", line.render_line());
    }
    out
}

pub fn run(args: CheckArgs) -> Result<u8> {
    let verdicts = verdicts(&args.inputs, &args.preset)?;
    let code = exit_code(&verdicts, args.fail_on);
    match args.format {
        Format::Text => print!("{}", text(&verdicts, code, args.fail_on)),
        Format::Json => {
            let items: Vec<String> = verdicts
                .iter()
                .map(|v| v.to_json().trim_end().to_owned())
                .collect();
            println!("[{}]", items.join(",\n"));
        }
        Format::Jsonl => print!("{}", jsonl(&verdicts)),
    }
    Ok(code)
}

#[cfg(test)]
mod tests {
    use super::*;
    use bammetrics_verdict::{Counts, PresetInfo, SCHEMA};

    fn verdict(sample: &str, status: Band) -> Verdict {
        Verdict {
            schema: SCHEMA.to_owned(),
            tool_version: "0.7.0".to_owned(),
            sample: sample.to_owned(),
            preset: PresetInfo {
                name: "generic".to_owned(),
                version: 1,
                hash: "0123456789abcdef".to_owned(),
                source: "built-in".to_owned(),
                note: "starting values".to_owned(),
            },
            status,
            counts: Counts::default(),
            files: Vec::new(),
            checks: Vec::new(),
        }
    }

    #[test]
    fn exit_codes_follow_the_table() {
        let pass = verdict("a", Band::Pass);
        let warn = verdict("b", Band::Warn);
        let fail = verdict("c", Band::Fail);
        let incomplete = verdict("d", Band::Incomplete);
        assert_eq!(exit_code(std::slice::from_ref(&pass), FailOn::Fail), 0);
        assert_eq!(exit_code(&[pass.clone(), warn.clone()], FailOn::Fail), 0);
        assert_eq!(exit_code(&[pass.clone(), warn.clone()], FailOn::Warn), 10);
        assert_eq!(exit_code(&[pass.clone(), fail.clone()], FailOn::Fail), 10);
        assert_eq!(
            exit_code(&[incomplete.clone(), warn.clone()], FailOn::Fail),
            11
        );
        assert_eq!(
            exit_code(&[incomplete.clone(), warn.clone()], FailOn::Warn),
            10
        );
        assert_eq!(exit_code(&[incomplete.clone(), fail], FailOn::Fail), 10);
        assert_eq!(exit_code(&[incomplete], FailOn::Warn), 11);
        assert_eq!(exit_code(&[], FailOn::Fail), 0);
    }

    #[test]
    fn jsonl_has_one_verdict_event_per_sample() {
        let lines = jsonl(&[verdict("a", Band::Pass), verdict("b", Band::Fail)]);
        let lines: Vec<&str> = lines.lines().collect();
        assert_eq!(lines.len(), 2);
        assert!(
            lines[0].starts_with("{\"event\":\"verdict\",\"sample\":\"a\",\"status\":\"pass\"")
        );
        assert!(lines[1].contains("\"status\":\"fail\""));
        assert!(lines[1].contains("\"hash\":\"0123456789abcdef\""));
    }
}

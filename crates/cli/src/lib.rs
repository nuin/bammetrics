//! The `bammetrics` command as a library: every verb's code path, so the
//! binary (`src/main.rs`) and the Python module (`crates/py`) run the same
//! engine, discovery, verdict and report code. [`main`] is the whole binary.

use std::process::ExitCode;

use bammetrics_core::{create_sequence_dictionary, registry};
use clap::error::{ContextKind, ContextValue, ErrorKind};
use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};

pub mod check;
pub mod collect;
pub mod discover;
pub mod engine;
pub mod error;
pub mod events;
pub mod format;
pub mod grammar;
pub mod info;
pub mod manifest;
pub mod plan;
pub mod python;
pub mod qc;
pub mod report;
pub mod verbs;

use error::VerbError;

#[derive(Parser)]
#[command(
    name = "bammetrics",
    version,
    about = "Picard-concordant BAM QC metrics",
    after_help = "Each tool verb runs one Picard class (bammetrics tools lists them); bammetrics \
                  help VERB lists its options in every spelling.\n\n\
                  Picard syntax also works: bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt [M=0.05 ...]\n\
                  bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt [R=ref.fa ...]\n\
                  bammetrics CollectHsMetrics I=in.bam O=out.txt BAIT_INTERVALS=b.interval_list TARGET_INTERVALS=t.interval_list"
)]
pub struct Cli {
    #[command(subcommand)]
    command: Command,
}

/// The verbs clap parses; the tool verbs, `collect` and the info verbs
/// go through the verb grammar ([`verbs`]).
#[derive(Subcommand)]
pub enum Command {
    /// Verdicts, sample pages and a run index from directories of Picard
    /// metrics files: SAMPLE.verdict.json, SAMPLE.report.html, index.html.
    Report(report::ReportArgs),
    /// Evaluate a preset over directories of metrics files, or read
    /// verdict.json files, and gate: print the checks and exit 10 or 11 on
    /// a bad verdict.
    Check(check::CheckArgs),
    /// The collectors a preset needs on one read of a BAM, then the
    /// verdict and the page: PREFIX.<metrics files>, PREFIX.verdict.json,
    /// PREFIX.report.html, PREFIX.manifest.json.
    Qc(qc::QcArgs),
    /// What qc with the same options would run: stages, collectors, which
    /// share the BAM read and why any does not, which steps --resume would
    /// skip. Nothing runs.
    Plan(qc::QcArgs),
}

/// The text of one clap context value.
fn context_text(value: &ContextValue) -> String {
    match value {
        ContextValue::String(s) => s.clone(),
        ContextValue::Strings(v) => v.join(", "),
        other => other.to_string(),
    }
}

/// A clap usage error on `verb` in the style of the grammar verbs:
/// `bammetrics <verb>: <what> (from <option>)`, a suggestion or the valid
/// values, then the usage line.
pub fn clap_error(verb: &str, e: &clap::Error) -> VerbError {
    let context = |kind| e.get(kind).map(context_text);
    let rendered = e.render().to_string();
    let first = rendered
        .lines()
        .next()
        .unwrap_or_default()
        .trim_start_matches("error: ")
        .to_owned();
    // `--fail-on <FAIL_ON>` is written as the option alone.
    let arg = context(ContextKind::InvalidArg)
        .map(|a| a.split_whitespace().next().unwrap_or_default().to_owned());
    let value = context(ContextKind::InvalidValue).unwrap_or_default();
    let mut error = match (e.kind(), &arg) {
        (ErrorKind::MissingRequiredArgument, Some(arg)) => {
            VerbError::usage(verb, format!("missing {arg}"))
        }
        (ErrorKind::InvalidValue, Some(arg)) if value.is_empty() => {
            VerbError::usage(verb, "missing value").from(arg.clone())
        }
        // `--threads x`: clap's reason follows the quoted option.
        (ErrorKind::ValueValidation, Some(arg)) => {
            let reason = first
                .rsplit_once("': ")
                .map_or_else(String::new, |(_, r)| format!(": {r}"));
            VerbError::usage(verb, format!("invalid value {value:?}{reason}")).from(arg.clone())
        }
        (ErrorKind::InvalidValue, Some(arg)) => VerbError::usage(
            verb,
            format!(
                "invalid value {:?}",
                context(ContextKind::InvalidValue).unwrap_or_default()
            ),
        )
        .from(arg.clone()),
        (ErrorKind::UnknownArgument, Some(arg)) => {
            VerbError::usage(verb, format!("unknown option {arg}"))
        }
        (ErrorKind::ArgumentConflict, Some(arg)) => VerbError::usage(
            verb,
            format!(
                "{arg} cannot be used with {}",
                context(ContextKind::PriorArg).unwrap_or_default()
            ),
        ),
        _ => VerbError::usage(verb, first),
    };
    if let Some(suggested) = context(ContextKind::SuggestedArg) {
        error = error.hint(format!("did you mean {suggested}?"));
    } else if let Some(valid) = context(ContextKind::ValidValue) {
        error = error.hint(format!("one of: {valid}"));
    }
    match rendered.lines().find(|l| l.starts_with("Usage: ")) {
        Some(usage) => error.hint(format!("usage: {}", &usage["Usage: ".len()..])),
        None => error.hint(format!("every option: bammetrics help {verb}")),
    }
}

/// Dispatch, in the order of docs/05 section 2.1: (1) a Picard class name
/// runs the Picard face, untouched; (2) a grammar verb (a tool verb,
/// `collect`, `tools`, `presets`, `help`, `completions`); (3) the clap
/// verbs `report`, `check`, `qc` and `plan`. The exit code of the command:
/// the Picard face exits 0 or 1, the grammar verbs 0, 1 or 2, `check` and
/// `qc --fail-on` 10 or 11 on a bad verdict (`bammetrics help check`).
pub fn run() -> error::Result<u8> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Some((name, rest)) = args.split_first() {
        if let Some(tool) = registry::find(name) {
            if let [flag] = rest
                && (flag == "--help" || flag == "-h")
                && let Some(verb) = verbs::for_class(name)
            {
                print!("{}", verb.help());
                return Ok(0);
            }
            // Picard's CreateSequenceDictionary names its .dict when O= is absent;
            // the verb and Python write standard output and rows instead.
            if tool.name == create_sequence_dictionary::TOOL.name {
                tool.run_cli(&create_sequence_dictionary::with_default_output(rest)?)?;
                return Ok(0);
            }
            tool.run_cli(rest)?;
            return Ok(0);
        }
        if let Some(verb) = verbs::find(name) {
            return Ok(verb.run(rest)?);
        }
        match name.as_str() {
            "collect" => return Ok(collect::run(rest)?),
            "tools" => {
                print!("{}", info::tools(&Cli::command()));
                return Ok(0);
            }
            "presets" => {
                print!("{}", info::presets()?);
                return Ok(0);
            }
            "help" => {
                print!("{}", info::help(Cli::command(), rest)?);
                return Ok(0);
            }
            "completions" => {
                info::completions(Cli::command(), rest)?;
                return Ok(0);
            }
            other if !other.starts_with('-') && !info::CLAP_VERBS.contains(&other) => {
                return Err(info::unknown_verb(other).into());
            }
            _ => {}
        }
    }
    let matches = info::command(Cli::command())
        .try_get_matches_from(std::env::args_os())
        .and_then(|m| Cli::from_arg_matches(&m));
    let cli = match matches {
        Ok(cli) => cli,
        Err(e) => match args.first().map(String::as_str) {
            Some(verb)
                if info::CLAP_VERBS.contains(&verb)
                    && !matches!(
                        e.kind(),
                        ErrorKind::DisplayHelp
                            | ErrorKind::DisplayVersion
                            | ErrorKind::DisplayHelpOnMissingArgumentOrSubcommand
                    ) =>
            {
                return Err(clap_error(verb, &e).into());
            }
            _ => e.exit(),
        },
    };
    match cli.command {
        Command::Report(args) => report::run(args),
        Command::Check(args) => check::run(args),
        Command::Qc(args) => Ok(engine::run(&args)),
        Command::Plan(args) => Ok(plan::run(&args)),
    }
}

/// The binary: [`run`], its error printed, its exit code returned.
pub fn main() -> ExitCode {
    match run() {
        Ok(code) => ExitCode::from(code),
        Err(error::Error::Verb(e)) => {
            eprintln!("{e}");
            ExitCode::from(e.code)
        }
        Err(e) => {
            eprintln!("bammetrics: {e}");
            ExitCode::from(e.exit_code())
        }
    }
}

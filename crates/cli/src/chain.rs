//! `bammetrics chain` (docs/chaining.md): several steps on one BAM in one
//! process. Records move from step to step in memory, with no shell pipe
//! and no intermediate file.
//!
//! ```text
//! bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]
//! bammetrics chain --plan FILE [OPTIONS] [BAM]
//! ```
//!
//! The chain's options and its source come first. Every `+` starts a
//! step. A step is a tool verb with a stage (a transform: it changes the
//! records later steps see), a one-pass tool verb (a tap: it reads the
//! records of its position and changes nothing), or `write [PATH]` (a tap
//! that writes those records as a BAM). Each step is parsed by its verb's
//! own flags, so kebab flags, `--PICARD_KEY`, `KEY=VALUE` and did-you-mean
//! work per step.
//!
//! The pieces, in the order a run uses them:
//! * [`args`]: the command line into a [`Spec`].
//! * [`planfile`]: a TOML plan file into the same [`Spec`].
//! * [`prepare`]: a [`Spec`] into a [`Chain`]. Every usage error is found
//!   here, before the first read: nothing is opened for records and
//!   nothing is written.
//! * [`describe`]: the plan text of `--dry-run` and of the `plan` verb.
//! * [`run`]: [`execute`], which runs a [`Chain`] through
//!   [`bammetrics_core::stage::run_chain`], writes every file and returns
//!   the [`Outcome`] with one [`StepOutcome`] per step.
//! * [`python`]: the keyword tables of the Python `chain` and `step`, and
//!   their `bammetrics.pyi` entries.
//!
//! A caller that is not the command line builds a [`Spec`], calls
//! [`prepare`] and then [`execute`]. The Python module spells the argv of
//! its keywords, so [`parse`] builds its [`Spec`].

mod args;
mod describe;
mod planfile;
mod prepare;
pub mod python;
mod run;

use crate::error::VerbError;

pub use args::{Asked, global_flags, help, parse, step_flags, step_help, write_flags};
pub use describe::text as plan_text;
pub use planfile::{parse as parse_plan, read as read_plan};
pub use prepare::{Chain, Planned, PlannedFile, PrefixOrigin, Role, prepare};
pub use run::{Outcome, StepOutcome, StepStatus, Written, execute, sample_name, summary};

/// The verb.
pub const VERB: &str = "chain";
/// The argument that starts a step.
pub const STEP: &str = "+";
/// The step that is not a tool verb.
pub const WRITE: &str = "write";

/// A chain as asked for: the source, the chain's options and the steps,
/// each still as text. The command line and a plan file both make one.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Spec {
    /// The BAM or CRAM, or `-` for standard input.
    pub source: Option<String>,
    /// The FASTA a CRAM is decoded with, handed to every step that reads
    /// a reference and names none.
    pub reference: Option<String>,
    /// The output prefix of every step that names no file of its own.
    pub prefix: Option<String>,
    pub threads: Option<usize>,
    /// Where a transform spills and a piped source is copied.
    pub tmp_dir: Option<String>,
    /// The `--format` list (`picard,json`) of every step without its own.
    pub formats: Option<String>,
    /// The `--also` list (`csv,markdown`).
    pub also: Option<String>,
    /// The sample name of the json and jsonl forms.
    pub sample: Option<String>,
    pub quiet: bool,
    pub verbose: bool,
    /// Run a chain of taps through the stage executor in place of the
    /// shared read of `qc`, to measure the executor.
    pub stage_executor: bool,
    pub steps: Vec<StepSpec>,
}

/// One step as asked for: its verb and its own arguments.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StepSpec {
    /// A tool verb (`insert-size`), its Picard class, or `write`.
    pub verb: String,
    /// The step's arguments, in any spelling its verb takes.
    pub args: Vec<String>,
    /// Where a plan file names the step (`plan.toml:12`), for errors.
    pub origin: Option<String>,
}

impl StepSpec {
    pub fn new(verb: &str, args: &[&str]) -> Self {
        Self {
            verb: verb.to_owned(),
            args: args.iter().map(|a| (*a).to_owned()).collect(),
            origin: None,
        }
    }
}

/// A usage error of the chain itself.
pub(crate) fn usage(what: impl Into<String>) -> VerbError {
    VerbError::usage(VERB, what)
}

/// True when the arguments of `plan` describe a chain: they hold a `+`
/// step or name a plan file. `plan` with `qc` arguments has neither.
pub fn is_chain(args: &[String]) -> bool {
    args.iter()
        .any(|arg| arg == STEP || arg == "--plan" || arg.starts_with("--plan="))
}

/// `bammetrics chain ARGS`: the exit code. 0 when every step finished, 1
/// when a step failed (the others still write their files), 2 for a usage
/// error (found before the first read).
pub fn run(args: &[String]) -> Result<u8, VerbError> {
    dispatch(args, false)
}

/// `bammetrics plan ARGS` for a chain: what `chain ARGS` would run.
pub fn run_plan(args: &[String]) -> Result<u8, VerbError> {
    dispatch(args, true)
}

fn dispatch(args: &[String], plan_only: bool) -> Result<u8, VerbError> {
    let (spec, dry_run) = match parse(args)? {
        Asked::Help(text) => {
            crate::events::print_stdout(&text);
            return Ok(0);
        }
        Asked::Chain { spec, dry_run } => (spec, dry_run),
    };
    let chain = prepare(&spec)?;
    if dry_run || plan_only {
        crate::events::print_stdout(&plan_text(&chain)?);
        return Ok(0);
    }
    run::command(chain)
}

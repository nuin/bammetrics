//! `bammetrics __chain`: the test entry of the stage executor
//! ([`bammetrics_core::stage`], docs/chaining.md). It is hidden: no help,
//! no completion, no `tools` line. The `chain` verb with its `+` grammar,
//! plan file and Python function replaces it in the next wave and keeps
//! [`run`]'s plan building.
//!
//! Picard syntax, one process, no pipe and no intermediate file:
//!
//! ```text
//! bammetrics __chain I=in.bam [R=ref.fa] [threads=N] [TMP_DIR=dir] \
//!     STEP=MarkDuplicates [M=dup.txt] [KEY=VALUE ...] \
//!     STEP=CollectInsertSizeMetrics O=insert.txt [KEY=VALUE ...] \
//!     STEP=write O=marked.bam [COMPRESSION_LEVEL=5] [CREATE_INDEX=true] [CREATE_MD5_FILE=true]
//! ```
//!
//! Every `STEP=` starts a step; the options after it are that step's, in
//! the tool's own Picard spelling. A step is a registry tool with a stage
//! (`MarkDuplicates`), a one-pass registry tool (a tap, which needs `O=`),
//! or `write`. Each step sees the records the step before it emits. Exit 0,
//! 1 when a step failed (the others still write their files), 2 for a usage
//! error.

use std::path::{Path, PathBuf};

use bammetrics_core::pass::default_threads;
use bammetrics_core::registry::{self, Invocation, Kind, Tool};
use bammetrics_core::stage::{self, ChainPlan, FileSource, Products, Step, WriteOptions};

use crate::error::VerbError;

/// The first argument that selects this entry.
pub const ENTRY: &str = "__chain";

type Result<T> = std::result::Result<T, VerbError>;

fn usage(what: impl Into<String>) -> VerbError {
    VerbError::usage(ENTRY, what)
}

/// The options in front of the first `STEP=`.
struct Globals {
    input: String,
    reference: Option<String>,
    threads: usize,
    tmp_dir: Option<PathBuf>,
    lower: bool,
}

/// `NAME=VALUE`, as the Picard face takes it.
fn key_value(arg: &str) -> Result<(&str, &str)> {
    arg.split_once('=')
        .ok_or_else(|| usage(format!("expected NAME=VALUE, got {arg:?}")))
}

fn boolean(name: &str, value: &str) -> Result<bool> {
    match value.to_ascii_lowercase().as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(usage(format!("{name}={value}: expected true or false"))),
    }
}

/// One step as given: its name and its own arguments.
type Named = (String, Vec<String>);

/// The arguments before the first `STEP=` and each step after it.
fn split(args: &[String]) -> Result<(Globals, Vec<Named>)> {
    let mut input = None;
    let mut globals = Globals {
        input: String::new(),
        reference: None,
        threads: default_threads(),
        tmp_dir: None,
        lower: true,
    };
    let mut steps: Vec<Named> = Vec::new();
    for arg in args {
        let (name, value) = key_value(arg)?;
        if name == "STEP" {
            steps.push((value.to_owned(), Vec::new()));
            continue;
        }
        if let Some((_, own)) = steps.last_mut() {
            own.push(arg.clone());
            continue;
        }
        match name {
            "INPUT" | "I" => input = Some(value.to_owned()),
            "REFERENCE_SEQUENCE" | "R" => globals.reference = Some(value.to_owned()),
            "TMP_DIR" => globals.tmp_dir = Some(PathBuf::from(value)),
            "threads" => {
                globals.threads = value
                    .parse()
                    .map_err(|_| usage(format!("threads={value}: not a number")))?;
            }
            // False runs a chain of taps through the stage executor in
            // place of the shared read of qc, to measure the executor.
            "LOWER" => globals.lower = boolean(name, value)?,
            _ => {
                return Err(usage(format!(
                    "unknown option {name} before the first STEP="
                )));
            }
        }
    }
    globals.input = input.ok_or_else(|| usage("INPUT is required"))?;
    if steps.is_empty() {
        return Err(usage("no STEP= given"));
    }
    Ok((globals, steps))
}

/// The options of the `write` step.
fn write_options(args: &[String]) -> Result<WriteOptions> {
    let mut path = None;
    let mut options = WriteOptions::new(Path::new(""));
    for arg in args {
        let (name, value) = key_value(arg)?;
        match name {
            "OUTPUT" | "O" => path = Some(PathBuf::from(value)),
            "COMPRESSION_LEVEL" => {
                options.compression_level = value
                    .parse()
                    .ok()
                    .filter(|level| *level <= 9)
                    .ok_or_else(|| usage(format!("COMPRESSION_LEVEL={value}: expected 0 to 9")))?;
            }
            "CREATE_INDEX" => options.create_index = boolean(name, value)?,
            "CREATE_MD5_FILE" => options.create_md5 = boolean(name, value)?,
            _ => return Err(usage(format!("write: unknown option {name}"))),
        }
    }
    options.path = path.ok_or_else(|| usage("write: OUTPUT is required"))?;
    Ok(options)
}

/// What a step leaves for after the run: the tool's parsed command line,
/// which names its files.
enum Role {
    Tap(Invocation),
    Stage(Invocation),
    Write,
}

/// The Picard arguments of a tool step: the chain's `INPUT`, its
/// `REFERENCE_SEQUENCE` when the tool reads one and the step names none,
/// then the step's own.
fn tool_args(tool: &Tool, globals: &Globals, own: &[String]) -> Vec<String> {
    let mut args = vec![
        format!("INPUT={}", globals.input),
        format!("threads={}", globals.threads),
    ];
    let reads_reference = tool
        .option("REFERENCE_SEQUENCE")
        .is_some_and(|opt| opt.kind == Kind::Value);
    let names_reference = own.iter().any(|arg| {
        let name = arg.trim_start_matches('-');
        name.starts_with("R=") || name.starts_with("REFERENCE_SEQUENCE=")
    });
    if let (Some(reference), true, false) = (&globals.reference, reads_reference, names_reference) {
        args.push(format!("REFERENCE_SEQUENCE={reference}"));
    }
    args.extend(own.iter().cloned());
    args
}

/// One tool step: the step for the plan and its role. `stages` are the
/// stages in front of it, which a tap's `StringHeader` line names.
fn tool_step(
    name: &str,
    globals: &Globals,
    own: &[String],
    stages: &[String],
) -> Result<(Step, Role)> {
    let core = |e: bammetrics_core::Error| VerbError::core(ENTRY, &e);
    let tool = registry::find(name)
        .ok_or_else(|| usage(format!("unknown step {name:?}: a Picard class or write")))?;
    let args = tool_args(tool, globals, own);
    let parse = || tool.parse(&args).map_err(core);
    if let Some(step) = stage::Step::tool_stage(tool, parse()?) {
        let inv = parse()?;
        registry::check_writable(inv.output.as_deref(), &inv.inputs).map_err(core)?;
        return Ok((step, Role::Stage(inv)));
    }
    let Some(step) = stage::Step::tool_tap(tool, parse()?).map_err(core)? else {
        return Err(usage(format!(
            "{name} does not read its input in one pass and has no stage: it cannot be a step yet"
        )));
    };
    let mut inv = parse()?;
    let Some(output) = inv.output.as_deref() else {
        return Err(usage(format!("{name}: a tap needs OUTPUT")));
    };
    registry::check_writable([output], &inv.inputs).map_err(core)?;
    if !stages.is_empty() {
        inv.command
            .push_str(&format!(" CHAIN={}", stages.join("+")));
    }
    Ok((step, Role::Tap(inv)))
}

/// Write what an `Ok` step left. A stage without `METRICS_FILE` writes
/// nothing; its counts go to standard error.
fn write_products(label: &str, role: &Role, products: &Products) -> bammetrics_core::Result<()> {
    if let Some(counters) = products.counters {
        eprintln!(
            "bammetrics {ENTRY}: {label}: {} records in, {} out, {} changed",
            counters.records_in, counters.records_out, counters.changed
        );
    }
    let inv = match role {
        Role::Tap(inv) => inv,
        Role::Stage(inv) if inv.output.is_some() => inv,
        Role::Stage(_) | Role::Write => return Ok(()),
    };
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    registry::write_files(inv, &products.output)
}

/// Run the chain `args` describe. The exit code: 0, or 1 when a step
/// failed.
pub fn run(args: &[String]) -> Result<u8> {
    let (globals, named) = split(args)?;
    let mut steps = Vec::with_capacity(named.len());
    let mut roles = Vec::with_capacity(named.len());
    let mut stages: Vec<String> = Vec::new();
    for (name, own) in &named {
        let (step, role) = if name == "write" {
            let options = write_options(own)?;
            // Found before the first read, which a barrier in front of
            // the step would otherwise finish first.
            let reads: Vec<PathBuf> = std::iter::once(&globals.input)
                .chain(&globals.reference)
                .map(PathBuf::from)
                .collect();
            registry::check_writable([options.path.as_path()], &reads)
                .map_err(|e| VerbError::core(ENTRY, &e))?;
            (Step::write(options, globals.threads), Role::Write)
        } else {
            tool_step(name, &globals, own, &stages)?
        };
        if matches!(role, Role::Stage(_)) {
            stages.push(name.clone());
        }
        steps.push(step);
        roles.push(role);
    }
    let labels: Vec<String> = steps.iter().map(|step| step.label.clone()).collect();
    let plan = ChainPlan {
        source: Box::new(FileSource::new(
            Path::new(&globals.input),
            globals.reference.as_deref().map(Path::new),
        )),
        steps,
        threads: globals.threads,
        tmp_dir: globals.tmp_dir,
        lower: globals.lower,
    };
    let results = stage::run_chain(plan, None).map_err(|e| VerbError::core(ENTRY, &e))?;
    let mut failed = false;
    for ((label, role), result) in labels.iter().zip(&roles).zip(results) {
        let written = result.and_then(|products| write_products(label, role, &products));
        if let Err(e) = written {
            eprintln!("bammetrics {ENTRY}: {label}: {e}");
            failed = true;
        }
    }
    Ok(u8::from(failed))
}

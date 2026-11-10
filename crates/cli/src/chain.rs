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
//! or `write`. Each step sees the records the step before it emits. A step
//! is named by its tool, with `#1`, `#2` after a name given more than once.
//!
//! Exit 0; 2 for a usage error, found before the first read; 1 when a file
//! cannot be read or written or when a step failed (the other steps still
//! write their files). Two steps that name one output file are a usage
//! error: every output of a chain needs its own file.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use bammetrics_core::pass::default_threads;
use bammetrics_core::registry::{self, Invocation, Kind, Tool, Ty};
use bammetrics_core::stage::{self, ChainPlan, FileSource, Products, Step, WriteOptions};

use crate::error::VerbError;

/// The first argument that selects this entry.
pub const ENTRY: &str = "__chain";
/// The step that is not a registry tool.
const WRITE: &str = "write";

type Result<T> = std::result::Result<T, VerbError>;

fn usage(what: impl Into<String>) -> VerbError {
    VerbError::usage(ENTRY, what)
}

/// A step's message under its label. A tool's own errors start with the
/// tool's name, which the label already says.
fn labelled(label: &str, tool: &str, message: &str) -> String {
    let message = message
        .strip_prefix(tool)
        .and_then(|rest| rest.strip_prefix(": "))
        .unwrap_or(message);
    format!("{label}: {message}")
}

/// A core error of step `label` found before the run: exit 2 for a usage
/// error, 1 for a file that cannot be read or written.
fn step_error(label: &str, tool: &str, e: &bammetrics_core::Error) -> VerbError {
    let mut error = VerbError::core(ENTRY, e);
    error.what = labelled(label, tool, &error.what);
    error
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

/// One label per step: the step's name, and `name#k` for the k-th step of
/// a name given more than once, so every message names one step.
fn labels(named: &[Named]) -> Vec<String> {
    let mut total: HashMap<&str, usize> = HashMap::new();
    for (name, _) in named {
        *total.entry(name).or_default() += 1;
    }
    let mut seen: HashMap<&str, usize> = HashMap::new();
    named
        .iter()
        .map(|(name, _)| {
            let k = seen.entry(name).or_default();
            *k += 1;
            if total[name.as_str()] == 1 {
                name.clone()
            } else {
                format!("{name}#{k}")
            }
        })
        .collect()
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

/// The `write` step and the files it may write. Its name, its directory
/// and the inputs are checked here, before the first read, which a barrier
/// in front of the step would otherwise finish first.
fn write_step(label: &str, globals: &Globals, own: &[String]) -> Result<(Step, Vec<PathBuf>)> {
    let options = write_options(own)?;
    let failed = |e: bammetrics_core::Error| step_error(label, WRITE, &e);
    options.check().map_err(failed)?;
    let files = options.files();
    let reads: Vec<PathBuf> = std::iter::once(&globals.input)
        .chain(&globals.reference)
        .map(PathBuf::from)
        .collect();
    registry::check_writable(files.iter().map(PathBuf::as_path), &reads).map_err(failed)?;
    Ok((Step::write(options, globals.threads), files))
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

/// The files the command line of a tap names as written: its file option
/// and every path option marked [`registry::Opt::written`]. A tool whose
/// `OUTPUT` is a prefix names its files when it finishes;
/// [`claim_late_files`] checks those.
fn declared_outputs(tool: &Tool, inv: &Invocation) -> Vec<PathBuf> {
    let written = tool
        .table()
        .into_iter()
        .filter(|opt| opt.kind == Kind::Value && opt.ty == Ty::Path && opt.written)
        .flat_map(|opt| inv.all(opt.picard()))
        .map(PathBuf::from);
    inv.output.iter().cloned().chain(written).collect()
}

/// One tool step: the step for the plan, its role and the files it names
/// as written. `stages` are the stages in front of it, which a tap's
/// `StringHeader` line names.
fn tool_step(
    (name, label): (&str, &str),
    globals: &Globals,
    own: &[String],
    stages: &[String],
) -> Result<(Step, Role, Vec<PathBuf>)> {
    let failed = |e: bammetrics_core::Error| step_error(label, name, &e);
    let tool = registry::find(name)
        .ok_or_else(|| usage(format!("unknown step {name:?}: a Picard class or write")))?;
    let args = tool_args(tool, globals, own);
    let parse = || tool.parse(&args).map_err(failed);
    let writable = |inv: &Invocation, files: Vec<PathBuf>| {
        registry::check_writable(files.iter().map(PathBuf::as_path), &inv.inputs)
            .map_err(failed)?;
        Ok::<_, VerbError>(files)
    };
    if let Some(step) = stage::Step::tool_stage(tool, parse()?) {
        let inv = parse()?;
        // A stage writes its metrics file alone: the step refuses the
        // options of the BAM when it is built, and `write` takes them.
        let files = writable(&inv, inv.output.iter().cloned().collect())?;
        return Ok((step, Role::Stage(inv), files));
    }
    let Some(step) = stage::Step::tool_tap(tool, parse()?).map_err(failed)? else {
        return Err(usage(format!(
            "{label}: {name} does not read its input in one pass and has no stage: it cannot be \
             a step yet"
        )));
    };
    let mut inv = parse()?;
    if inv.output.is_none() {
        return Err(usage(format!("{label}: a tap needs OUTPUT")));
    }
    let files = writable(&inv, declared_outputs(tool, &inv))?;
    if !stages.is_empty() {
        inv.command
            .push_str(&format!(" CHAIN={}", stages.join("+")));
    }
    Ok((step, Role::Tap(inv), files))
}

/// The files of every step so far, each with the label of its step.
type Outputs = Vec<(String, PathBuf)>;

/// Add the files `products` names that the command line of step `label`
/// did not (the files under an `OUTPUT` prefix). One that another step
/// writes fails this step, which then replaces nothing.
fn claim_late_files(
    label: &str,
    products: &Products,
    outputs: &mut Outputs,
) -> bammetrics_core::Result<()> {
    for (path, _) in &products.output.extra {
        let own = outputs
            .iter()
            .any(|(owner, file)| owner == label && file == path);
        if own {
            continue;
        }
        if let Some(other) = stage::writer_of(outputs, path) {
            return Err(bammetrics_core::Error::Usage(format!(
                "{} is a file of {other}: every output of a chain needs its own file",
                path.display()
            )));
        }
        outputs.push((label.to_owned(), path.clone()));
    }
    Ok(())
}

/// Write what an `Ok` step left. A stage without `METRICS_FILE` writes
/// nothing; its counts go to standard error.
fn write_products(
    label: &str,
    role: &Role,
    products: &Products,
    outputs: &mut Outputs,
) -> bammetrics_core::Result<()> {
    if let Some(counters) = products.counters {
        eprintln!(
            "bammetrics {ENTRY}: {label}: {} records in, {} out, {} changed",
            counters.records_in, counters.records_out, counters.changed
        );
    }
    for note in &products.notes {
        eprintln!("bammetrics {ENTRY}: {label}: {note}");
    }
    let inv = match role {
        Role::Tap(inv) => inv,
        Role::Stage(inv) if inv.output.is_some() => inv,
        Role::Stage(_) | Role::Write => return Ok(()),
    };
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    claim_late_files(label, products, outputs)?;
    registry::write_files(inv, &products.output)
}

/// Run the chain `args` describe. The exit code: 0, or 1 when a step
/// failed.
pub fn run(args: &[String]) -> Result<u8> {
    let (globals, named) = split(args)?;
    let labels = labels(&named);
    let mut steps = Vec::with_capacity(named.len());
    let mut roles = Vec::with_capacity(named.len());
    let mut stages: Vec<String> = Vec::new();
    let mut outputs = Outputs::new();
    for ((name, own), label) in named.iter().zip(&labels) {
        let (mut step, role, files) = if name == WRITE {
            let (step, files) = write_step(label, &globals, own)?;
            (step, Role::Write, files)
        } else {
            tool_step((name, label), &globals, own, &stages)?
        };
        if matches!(role, Role::Stage(_)) {
            stages.push(name.clone());
        }
        step.label.clone_from(label);
        outputs.extend(files.into_iter().map(|file| (label.clone(), file)));
        steps.push(step);
        roles.push(role);
    }
    let core = |e: bammetrics_core::Error| VerbError::core(ENTRY, &e);
    stage::check_distinct_outputs(&outputs).map_err(core)?;
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
    let results = stage::run_chain(plan, None).map_err(core)?;
    let mut failed = false;
    let steps = named.iter().zip(&labels).zip(&roles).zip(results);
    for ((((name, _), label), role), result) in steps {
        let written =
            result.and_then(|products| write_products(label, role, &products, &mut outputs));
        if let Err(e) = written {
            let message = labelled(label, name, &e.to_string());
            eprintln!("bammetrics {ENTRY}: {message}");
            failed = true;
        }
    }
    Ok(u8::from(failed))
}

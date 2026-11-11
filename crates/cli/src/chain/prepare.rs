//! A [`Spec`] into a [`Chain`]: every step resolved to its verb, its role,
//! the stream it reads, its Picard command line and the files it writes.
//! Every usage error of a chain is found here, before the first read: no
//! record is read and nothing is written or created.
//!
//! Names (docs/07 section 2.2). The source stream is `in`. A transform
//! leaves a stream under its label (`marked` after `mark-duplicates`, or
//! `--as NAME`). A tap on the last stream writes `PREFIX<extension>` with
//! the extensions of [`file_extensions`], the names `qc` writes. A tap on
//! an earlier stream writes `PREFIX.<stream><extension>`. A transform's
//! own metrics file is `PREFIX.<stream>.duplication_metrics`. `write`
//! without a path writes `PREFIX.<stream>.bam`. A step's own `-o` replaces
//! the name.

use std::path::{Path, PathBuf};

use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::registry::{Kind, Tool, Ty};
use bammetrics_core::stage::{self, Outputs, WriteOptions};
use bammetrics_core::{cram, java, pass};
use clap::ValueEnum;

use super::args::{
    AS, CHAIN_KEYS, COMPRESSION_LEVEL, INDEX, MD5, STEPS_HINT, belongs_to_chain,
    chain_option_on_step, find_step, is_transform, not_a_step, step_parser_flags, unknown_step,
    write_flags,
};
use super::claims::{
    Read, check_one_stdout, claim_outputs, inputs_of, refuse_overwriting_inputs, refuse_unresolved,
};
use super::{EMPTY_PREFIX_HINT, Spec, StepSpec, VERB, WRITE, usage};
use crate::also::Also;
use crate::error::VerbError;
use crate::format::{self, Format};
use crate::grammar::{self, Flag, Key, Parsed};
use crate::qc::with_suffix;
use crate::verbs::{self, Files, Verb};

/// The stream of the source.
const SOURCE_STREAM: &str = "in";

/// The metrics file a transform writes under the prefix: the tool, its
/// file option and the extension after `PREFIX.<stream>`.
const TRANSFORM_FILES: [(&str, &str, &str); 1] =
    [("MarkDuplicates", "METRICS_FILE", ".duplication_metrics")];

/// What a step is in the chain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    /// It changes the records later steps see (`mark-duplicates`).
    Transform,
    /// A one-pass tool on the records of its position.
    Tap,
    /// `write`: the records of its position as a BAM.
    Write,
}

impl Role {
    pub fn name(self) -> &'static str {
        match self {
            Role::Transform => "transform",
            Role::Tap => "tap",
            Role::Write => "write",
        }
    }
}

/// A file a step plans to write: `kind` is `metrics`, `gaps`, `bam`,
/// `index`, `md5`, or the name of a further format (`json`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PlannedFile {
    pub path: PathBuf,
    pub kind: &'static str,
}

/// Where the chain's prefix came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrefixOrigin {
    /// `-o PREFIX`.
    Option,
    /// `-o` named a directory; the source name is the file part.
    Directory,
    /// No `-o`: `./` plus the source name.
    SourceName,
}

/// One step, ready to run.
pub struct Planned {
    /// Its place in the chain, from 1.
    pub number: usize,
    /// The verb (`insert-size`, `mark-duplicates`) or `write`.
    pub verb: &'static str,
    /// The registry tool; `None` for `write`.
    pub tool: Option<&'static Tool>,
    pub role: Role,
    /// `step 2 (insert-size)`: how every message names the step.
    pub label: String,
    /// The stream it reads: `in`, or the label of the transform in front.
    pub reads: String,
    /// The stream a transform leaves.
    pub leaves: Option<String>,
    /// The segment it runs in, from 1. A transform takes its records in
    /// this segment and emits them in the next.
    pub segment: usize,
    /// The files it writes under names known before it runs.
    pub files: Vec<PlannedFile>,
    /// A prefix whose files the tool names when it finishes
    /// (`sequencing-artifact`).
    pub late_prefix: Option<PathBuf>,
    /// It writes standard output (`-o -`, `write -`).
    pub stdout: bool,
    /// The `StringHeader` command line of its metrics files.
    pub command: Option<String>,
    pub(crate) argv: Vec<String>,
    /// The file of the tool's main output.
    pub(crate) output: Option<PathBuf>,
    /// One line per option the tool ignores.
    pub(crate) notes: Vec<String>,
    /// The formats of each metrics file, `picard` first.
    pub(crate) formats: Vec<Format>,
    pub(crate) write: Option<WriteOptions>,
    /// The files the tool reads by name.
    pub(crate) inputs: Vec<PathBuf>,
}

/// A chain, ready to run: nothing has been read or written.
pub struct Chain {
    /// The source as given: a path, or `-`.
    pub source: String,
    /// The file the source is read from (`/dev/stdin` for `-`).
    pub source_path: PathBuf,
    pub reference: Option<PathBuf>,
    pub prefix: PathBuf,
    pub prefix_origin: PrefixOrigin,
    /// Some step is named from the prefix.
    pub prefix_used: bool,
    pub threads: usize,
    pub tmp_dir: Option<PathBuf>,
    pub sample: Option<String>,
    pub quiet: bool,
    pub verbose: bool,
    pub stage_executor: bool,
    pub steps: Vec<Planned>,
    /// `PREFIX.report.md` of `--also markdown`.
    pub report: Option<PathBuf>,
    /// Every planned file with the label of its step.
    pub(crate) outputs: Outputs,
    /// Every file the chain reads: no step may write one.
    pub(crate) inputs: Outputs,
}

impl Chain {
    /// True when the source can be read a second time: a regular file. A
    /// pipe is copied to the temporary directory once in front of a
    /// transform that reads its input twice.
    pub fn replayable(&self) -> bool {
        std::fs::metadata(&self.source_path).is_ok_and(|m| m.is_file())
            && !self.source_path.starts_with("/dev")
    }

    /// The stream after the last transform: `in` for a chain without one.
    /// A tap on it writes `PREFIX<extension>`, the names `qc` writes.
    pub fn last_stream(&self) -> &str {
        self.steps
            .iter()
            .rev()
            .find_map(|step| step.leaves.as_deref())
            .unwrap_or(SOURCE_STREAM)
    }

    /// The segments: one, plus one per transform.
    pub fn segments(&self) -> usize {
        1 + self
            .steps
            .iter()
            .filter(|step| step.role == Role::Transform)
            .count()
    }
}

/// What a step runs.
#[derive(Clone, Copy)]
enum What {
    Tool(&'static Verb, Role),
    Write,
}

/// A step with its arguments parsed by its verb's flags.
struct Draft<'a> {
    spec: &'a StepSpec,
    number: usize,
    label: String,
    what: What,
    flags: Vec<Flag>,
    parsed: Parsed,
}

impl Draft<'_> {
    /// Where an error of this step is: its plan-file line, then its label.
    fn place(&self) -> String {
        place(self.spec, &self.label)
    }

    /// An error of the step's own verb as an error of the chain.
    fn error(&self, mut error: VerbError) -> VerbError {
        let what = match self.what {
            What::Tool(verb, _) => own_message(verb.tool().name, &error.what),
            What::Write => &error.what,
        };
        error.what = format!("{}: {what}", self.place());
        VERB.clone_into(&mut error.verb);
        error
    }

    fn usage(&self, what: impl AsRef<str>) -> VerbError {
        usage(format!("{}: {}", self.place(), what.as_ref()))
    }

    fn core(&self, verb: &Verb, e: &bammetrics_core::Error) -> VerbError {
        self.error(verb.core_error(&self.parsed, e))
    }
}

fn place(spec: &StepSpec, label: &str) -> String {
    match &spec.origin {
        Some(origin) => format!("{origin}: {label}"),
        None => label.to_owned(),
    }
}

/// A tool's message without the tool's name in front, which the step's
/// label already says.
fn own_message<'m>(tool: &str, message: &'m str) -> &'m str {
    message
        .strip_prefix(tool)
        .and_then(|rest| rest.strip_prefix(": "))
        .unwrap_or(message)
}

/// The verb and role of a step that is not `write`.
fn tool_role(spec: &StepSpec, number: usize) -> Result<(&'static Verb, Role), VerbError> {
    let label = format!("step {number}");
    let verb = find_step(&spec.verb)
        .ok_or_else(|| unknown_step(&spec.verb, Some(&place(spec, &label))))?;
    if is_transform(verb) {
        return Ok((verb, Role::Transform));
    }
    if verb.tool().collector.is_some() {
        return Ok((verb, Role::Tap));
    }
    let label = format!("step {number} ({})", verb.name);
    Err(not_a_step(verb, &place(spec, &label)).hint(STEPS_HINT))
}

/// The spelling of a file a step names by an empty value (`-o ''`, a
/// `PATH` of `write` that is `''`). As a prefix it would name hidden files
/// (`.gc_bias.detail_metrics`); as a file it fails when the step writes.
fn empty_path<'a>(draft: &'a Draft) -> Option<&'a str> {
    let is_path = |key: Key| {
        let mut flags = draft.flags.iter();
        flags.any(|f| f.key == key && f.ty == Ty::Path)
    };
    let empty = |given: &&grammar::Given| given.values.iter().any(String::is_empty);
    let mut options = draft.parsed.given.iter().filter(|g| is_path(g.key));
    // Only `write` takes a positional argument, its PATH.
    let path =
        matches!(draft.what, What::Write) && draft.parsed.positional.iter().any(String::is_empty);
    options
        .find(empty)
        .map(|given| given.spelled.as_str())
        .or_else(|| path.then_some("PATH"))
}

/// Step `number` with its arguments parsed, and the arguments that belong
/// to the chain refused.
fn draft(spec: &StepSpec, number: usize) -> Result<Draft<'_>, VerbError> {
    let (what, name, flags) = if spec.verb == WRITE {
        (What::Write, WRITE, write_flags())
    } else {
        let (verb, role) = tool_role(spec, number)?;
        (What::Tool(verb, role), verb.name, step_parser_flags(verb))
    };
    let label = format!("step {number} ({name})");
    let parsed = grammar::parse(name, &flags, &spec.args).map_err(|mut error| {
        if let Some(belongs) = chain_option_on_step(&error, &place(spec, &label)) {
            return belongs;
        }
        error.what = format!("{}: {}", place(spec, &label), error.what);
        VERB.clone_into(&mut error.verb);
        // `write` is a step, not a verb: the chain's help lists its options.
        let help = "every option: bammetrics help write";
        for hint in error.hints.iter_mut().filter(|hint| *hint == help) {
            "every option: bammetrics chain BAM + write --help".clone_into(hint);
        }
        error
    })?;
    let draft = Draft {
        spec,
        number,
        label,
        what,
        flags,
        parsed,
    };
    if let Some(spelled) = empty_path(&draft) {
        return Err(draft
            .usage(format!("{spelled} names no file: its value is empty"))
            .from(spelled));
    }
    let What::Tool(verb, role) = draft.what else {
        return Ok(draft);
    };
    if let Some(extra) = draft.parsed.positional.first() {
        let mut error = draft
            .usage("a step takes no input, the chain's source feeds every step")
            .from(extra.clone());
        if let Some(hint) = draft.parsed.bool_value_hint(0) {
            error = error.hint(hint);
        }
        return Err(error);
    }
    if let Some(given) = draft
        .parsed
        .given
        .iter()
        .find(|given| CHAIN_KEYS.contains(&given.key))
    {
        return Err(belongs_to_chain(&draft.place(), &given.spelled));
    }
    if let (Role::Transform, Files::Bam, Some(given)) =
        (role, verb.files, draft.parsed.last(Key::Output))
    {
        return Err(draft
            .usage(format!(
                "{} names a BAM, and a transform writes none in a chain",
                given.spelled
            ))
            .from(given.spelled.clone())
            .hint(format!(
                "add + write {} after it; --metrics-file names its metrics",
                given.values[0]
            )));
    }
    Ok(draft)
}

/// Where a step stands in the chain.
struct Position {
    reads: String,
    leaves: Option<String>,
    segment: usize,
    /// The Picard classes of the transforms in front, joined by `+`: the
    /// `CHAIN=` of a tap's `StringHeader` line.
    behind: String,
}

/// The name of the stream a transform leaves: `--as`, or its default.
fn stream_name(draft: &Draft, verb: &Verb) -> Result<String, VerbError> {
    let Some(given) = draft.parsed.last(AS) else {
        return Ok(stage::label_of(verb.tool()).unwrap_or(verb.name).to_owned());
    };
    let name = &given.values[0];
    let plain = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    if !plain {
        return Err(draft
            .usage(format!(
                "{name:?} cannot name a stream: letters, digits, _ and - only"
            ))
            .from(given.spelled.clone()));
    }
    if name == SOURCE_STREAM {
        return Err(draft
            .usage(format!(
                "the stream name {name} is taken: it names the source"
            ))
            .from(given.spelled.clone()));
    }
    Ok(name.clone())
}

/// The stream each step reads and the segment it runs in. A second
/// transform is refused: it reads its input twice, and only the source
/// can be read again in this version.
fn positions(drafts: &[Draft]) -> Result<Vec<Position>, VerbError> {
    let mut stream = SOURCE_STREAM.to_owned();
    let mut segment = 1;
    let mut behind: Vec<&str> = Vec::new();
    let mut first: Option<&Draft> = None;
    let mut positions = Vec::with_capacity(drafts.len());
    for draft in drafts {
        let mut position = Position {
            reads: stream.clone(),
            leaves: None,
            segment,
            behind: behind.join("+"),
        };
        if let What::Tool(verb, Role::Transform) = draft.what {
            if let Some(first) = first {
                return Err(draft
                    .usage(format!(
                        "it reads its input twice and cannot follow {}: only a file source \
                         can be read again in this version",
                        first.label
                    ))
                    .hint("wave W05 of docs/07 spools the stream for a second transform"));
            }
            first = Some(draft);
            stream = stream_name(draft, verb)?;
            position.leaves = Some(stream.clone());
            behind.push(verb.class);
            segment += 1;
        }
        positions.push(position);
    }
    Ok(positions)
}

/// What every step is planned against.
struct Context<'a> {
    source: &'a str,
    is_cram: bool,
    reference: Option<&'a str>,
    prefix: &'a Path,
    /// The stream after the last transform.
    last_stream: &'a str,
    /// The chain's `--format` list.
    formats: &'a [Format],
    /// The per-file formats of `--also`.
    also: &'a [Format],
    threads: usize,
}

impl Context<'_> {
    /// The prefix of a step that reads `stream`: the chain's on the last
    /// stream, `PREFIX.<stream>` in front of a transform.
    fn base(&self, stream: &str) -> PathBuf {
        if stream == self.last_stream {
            self.prefix.to_path_buf()
        } else {
            with_suffix(self.prefix, &format!(".{stream}"))
        }
    }
}

/// `-o` when it names a directory (it ends in `/` or is one): the source
/// name becomes the file part.
fn prefix_under(text: &str, source: &str) -> (PathBuf, bool) {
    let path = PathBuf::from(text);
    if text.ends_with('/') || path.is_dir() {
        (path.join(verbs::stem(source)), true)
    } else {
        (path, false)
    }
}

/// A format list that does not start with `picard` is refused: a chain
/// names Picard files, and the further formats go beside them.
fn picard_first(formats: &[Format], value: &str, spelled: &str) -> Result<(), VerbError> {
    let first = formats[0];
    if first == Format::Picard {
        return Ok(());
    }
    Err(VerbError::usage(
        VERB,
        format!(
            "--format {value} starts with {}: a chain names Picard files",
            first.name()
        ),
    )
    .from(spelled)
    .hint(format!(
        "picard,{} writes the Picard files and each as FILE.{} beside it",
        first.name(),
        first.extension()
    )))
}

/// The formats of a step's metrics files: its own `--format` or the
/// chain's, then the per-file formats of `--also`.
fn step_formats(cx: &Context, draft: &Draft) -> Result<Vec<Format>, VerbError> {
    let mut formats = match draft.parsed.last(Key::Format) {
        None => cx.formats.to_vec(),
        Some(given) => format::parse_list(VERB, &given.values[0], &given.spelled)
            .and_then(|list| picard_first(&list, &given.values[0], &given.spelled).map(|()| list))
            .map_err(|mut error| {
                error.what = format!("{}: {}", draft.place(), error.what);
                error
            })?,
    };
    for format in cx.also {
        if !formats.contains(format) {
            formats.push(*format);
        }
    }
    Ok(formats)
}

/// A value an option of type int or float cannot hold is a usage error
/// here, before the first read; the tool itself would find it only when
/// its collector is built.
fn check_types(draft: &Draft) -> Result<(), VerbError> {
    for given in &draft.parsed.given {
        let Key::Picard(name) = given.key else {
            continue;
        };
        let Some(flag) = draft.flags.iter().find(|flag| flag.key == given.key) else {
            continue;
        };
        for value in &given.values {
            let valid = value == "null"
                || match flag.ty {
                    Ty::Int => value.parse::<i64>().is_ok(),
                    Ty::Float => value.parse::<f64>().is_ok() || java::parse_double(value).is_ok(),
                    Ty::Bool | Ty::Text | Ty::Path => true,
                };
            if !valid {
                return Err(draft
                    .usage(format!("{name}={value}: not a valid value"))
                    .from(given.spelled.clone()));
            }
        }
    }
    Ok(())
}

/// The file options a step gets from a prefix or from its own `-o`.
#[derive(Default)]
struct Named {
    /// `OPTION=PATH`, in the order of the tool's file table.
    args: Vec<String>,
    stdout: bool,
    uses_prefix: bool,
}

/// True when `own` (Picard arguments) sets `option`.
fn sets(own: &[String], option: &str) -> bool {
    own.iter().any(|arg| {
        arg.strip_prefix(option)
            .is_some_and(|rest| rest.starts_with('='))
    })
}

/// The file options of a tap: its own `-o` (a file, `-`, or a prefix for a
/// tool that writes several), else the chain's prefix with the names `qc`
/// writes. A file option the step gives itself keeps its value.
fn tap_names(
    cx: &Context,
    draft: &Draft,
    verb: &Verb,
    position: &Position,
    own: &[String],
) -> Result<Named, VerbError> {
    let tool = verb.tool();
    let output = tool.output_name();
    let own_output = draft.parsed.last(Key::Output);
    let mut named = Named::default();
    let base = match (own_output.map(|given| given.values[0].as_str()), verb.files) {
        (Some("-"), Files::Single) => {
            named.stdout = true;
            return Ok(named);
        }
        (Some("-"), _) => {
            return Err(draft
                .usage("standard output takes one file; this step writes several")
                .hint("give the step -o PREFIX"));
        }
        (Some(path), Files::Single) => {
            named.args.push(format!("{output}={path}"));
            return Ok(named);
        }
        (Some(path), _) => prefix_under(path, cx.source).0,
        (None, _) => {
            named.uses_prefix = true;
            cx.base(&position.reads)
        }
    };
    let table = file_extensions(tool.name).unwrap_or_default();
    if table.is_empty() {
        return Err(draft
            .usage(format!("{} has no name under a prefix", verb.class))
            .hint("give the step -o FILE"));
    }
    for (option, extension) in table {
        let main = *option == "OUTPUT";
        let option = if main { output } else { option };
        if main || !sets(own, option) {
            let path = with_suffix(&base, extension);
            named.args.push(format!("{option}={}", path.display()));
        }
    }
    Ok(named)
}

/// The metrics file of a transform: its own option, else
/// `PREFIX.<stream>.duplication_metrics`.
fn transform_names(cx: &Context, tool: &Tool, position: &Position, own: &[String]) -> Named {
    let mut named = Named::default();
    let stream = position.leaves.as_deref().unwrap_or_default();
    for (_, option, extension) in TRANSFORM_FILES
        .iter()
        .filter(|(name, option, _)| *name == tool.name && !sets(own, option))
    {
        let path = with_suffix(cx.prefix, &format!(".{stream}{extension}"));
        named.args.push(format!("{option}={}", path.display()));
        named.uses_prefix = true;
    }
    named
}

/// `REFERENCE_SEQUENCE` for a tool step: the chain's reference when the
/// tool reads one, or decodes a CRAM source, and the step names none.
fn reference_arg(cx: &Context, tool: &Tool, own: &[String]) -> Option<String> {
    let reference = cx.reference?;
    let reads = tool
        .option("REFERENCE_SEQUENCE")
        .is_some_and(|opt| opt.kind == Kind::Value);
    ((reads || cx.is_cram) && !sets(own, "REFERENCE_SEQUENCE"))
        .then(|| format!("REFERENCE_SEQUENCE={reference}"))
}

/// The kind of the file option `option` of `tool` names: `metrics` for the
/// tool's main file and every file of its prefix table, `gaps` for the BED
/// of `CollectCoverageGaps`, `file` for any other file it writes (a
/// per-target table), which has one form.
fn file_kind(tool: &Tool, option: &str) -> &'static str {
    let in_table = file_extensions(tool.name)
        .unwrap_or_default()
        .iter()
        .any(|(name, _)| *name == option);
    match option {
        "GAPS_FILE" => "gaps",
        _ if in_table || option == tool.output_name() => "metrics",
        _ => "file",
    }
}

/// True for a tool whose `OUTPUT` is a prefix: it names its files when it
/// finishes.
fn output_is_prefix(tool: &Tool) -> bool {
    matches!(file_extensions(tool.name), Some([("OUTPUT", "")]))
}

/// The files `argv` names as written, each with a sibling per further
/// format when it is a metrics file.
fn planned_files(tool: &Tool, argv: &[String], formats: &[Format]) -> Vec<PlannedFile> {
    let mut files = Vec::new();
    for arg in argv {
        let Some((name, value)) = arg.split_once('=') else {
            continue;
        };
        let Some(opt) = tool.option(name) else {
            continue;
        };
        let written = match opt.kind {
            Kind::Output => !output_is_prefix(tool),
            Kind::Value => opt.ty == Ty::Path && opt.written,
            Kind::Ignored(_) | Kind::Rejected => false,
        };
        if !written || value == "null" {
            continue;
        }
        let kind = file_kind(tool, opt.picard());
        let path = PathBuf::from(value);
        let siblings = formats[1..]
            .iter()
            .filter(|_| kind == "metrics")
            .map(|f| PlannedFile {
                path: format::sibling(&path, f.extension()),
                kind: f.name(),
            })
            .collect::<Vec<_>>();
        files.push(PlannedFile { path, kind });
        files.extend(siblings);
    }
    files
}

/// A tap or a transform: its Picard command line, checked as the tool's
/// own verb checks it, and its files.
fn tool_step(
    cx: &Context,
    mut draft: Draft,
    (verb, role): (&'static Verb, Role),
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let tool = verb.tool();
    let mut formats = step_formats(cx, &draft)?;
    verb.check_one_spelling(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    check_types(&draft)?;
    draft.parsed.positional.push(cx.source.to_owned());
    let picard = verb
        .picard_arguments(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    let (input, own) = picard.split_first().expect("the source is the first input");
    let named = match role {
        Role::Transform => transform_names(cx, tool, &position, own),
        _ => tap_names(cx, &draft, verb, &position, own)?,
    };
    if named.stdout {
        // Standard output holds the Picard file alone: the chain's further
        // formats do not apply, and the step's own list is an error.
        if let Some(given) = draft.parsed.last(Key::Format).filter(|_| formats.len() > 1) {
            return Err(draft
                .usage(format!(
                    "--format {} writes {} files; standard output holds one",
                    given.values[0],
                    formats.len()
                ))
                .from(given.spelled.clone()));
        }
        formats.truncate(1);
    }
    let mut argv = vec![input.clone()];
    argv.extend(named.args.iter().cloned());
    argv.extend(reference_arg(cx, tool, own));
    argv.extend(own.iter().cloned());
    verb.check_required(&draft.flags, &argv)
        .and_then(|()| verb.check_inputs(&draft.parsed, &argv))
        .map_err(|e| draft.error(e))?;
    let mut inv = tool.parse(&argv).map_err(|e| draft.core(verb, &e))?;
    inv.threads = cx.threads;
    match role {
        Role::Transform => stage::check_stage_options(tool, &inv),
        _ => inv
            .limits()
            .and_then(|_| bammetrics_core::tap_options::check(tool, &inv)),
    }
    .map_err(|e| draft.core(verb, &e))?;
    if role == Role::Tap && !position.behind.is_empty() {
        inv.command.push_str(&format!(" CHAIN={}", position.behind));
    }
    let late_prefix = inv.output.clone().filter(|_| output_is_prefix(tool));
    let planned = Planned {
        number: draft.number,
        verb: verb.name,
        tool: Some(tool),
        role,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: position.leaves,
        segment: position.segment,
        files: planned_files(tool, &argv, &formats),
        late_prefix,
        stdout: named.stdout,
        command: Some(inv.command),
        argv,
        output: inv.output,
        notes: inv.notes,
        formats,
        write: None,
        inputs: inv.inputs,
    };
    Ok((planned, named.uses_prefix))
}

/// The kind of a file `write` makes, by its place in
/// [`WriteOptions::files`].
fn write_kinds(options: &WriteOptions) -> Vec<&'static str> {
    let mut kinds = vec!["bam"];
    if options.create_index {
        kinds.push("index");
    }
    if options.create_md5 {
        kinds.push("md5");
    }
    kinds
}

/// The `write` step: its path (its own, `-`, or `PREFIX.<stream>.bam`) and
/// its options, with a name of a later wave refused.
fn write_step(
    cx: &Context,
    draft: &Draft,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let parsed = &draft.parsed;
    let paths: Vec<&String> = parsed
        .positional
        .iter()
        .chain(parsed.all(Key::Output).flat_map(|given| &given.values))
        .collect();
    if let [_, second, ..] = paths.as_slice() {
        return Err(draft
            .usage(format!("one PATH expected, got {}", paths.len()))
            .from((*second).clone()));
    }
    let (path, stdout, uses_prefix) = match paths.first().map(|path| path.as_str()) {
        Some("-") => (PathBuf::from("/dev/stdout"), true, false),
        Some(path) => (PathBuf::from(path), false, false),
        None => {
            let name = format!(".{}.bam", position.reads);
            (with_suffix(cx.prefix, &name), false, true)
        }
    };
    let mut options = WriteOptions::new(&path);
    if let Some(given) = parsed.last(COMPRESSION_LEVEL) {
        options.compression_level = given.values[0]
            .parse()
            .ok()
            .filter(|level| *level <= 9)
            .ok_or_else(|| {
                draft
                    .usage(format!("expected 0 to 9, got {:?}", given.values[0]))
                    .from(given.spelled.clone())
            })?;
    }
    options.create_index = parsed.flag(INDEX);
    options.create_md5 = parsed.flag(MD5);
    options
        .check()
        .map_err(|e| draft.error(VerbError::core(VERB, &e)))?;
    if let (true, Some(given)) = (stdout, parsed.last(INDEX).or_else(|| parsed.last(MD5))) {
        return Err(draft
            .usage(format!(
                "{} needs a file: standard output takes the BAM alone",
                given.spelled
            ))
            .from(given.spelled.clone()));
    }
    let files = options
        .files()
        .into_iter()
        .zip(write_kinds(&options))
        .map(|(path, kind)| PlannedFile { path, kind })
        .collect();
    let planned = Planned {
        number: draft.number,
        verb: WRITE,
        tool: None,
        role: Role::Write,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: None,
        segment: position.segment,
        files,
        late_prefix: None,
        stdout,
        command: None,
        argv: Vec::new(),
        output: None,
        notes: Vec::new(),
        formats: Vec::new(),
        write: Some(options),
        inputs: Vec::new(),
    };
    Ok((planned, uses_prefix))
}

/// The source: it must be given and, as a path, exist.
fn source(spec: &Spec) -> Result<(String, PathBuf), VerbError> {
    let given = spec.source.clone().ok_or_else(|| {
        usage("missing BAM")
            .hint("bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]")
            .hint("- reads standard input")
    })?;
    if given == "-" {
        return Ok((given, PathBuf::from("/dev/stdin")));
    }
    let path = PathBuf::from(&given);
    if !path.exists() {
        return Err(VerbError::input(VERB, format!("{given} not found")).from("BAM"));
    }
    Ok((given, path))
}

/// The reference: it must exist, and a CRAM source needs one.
fn reference(spec: &Spec, source: &Path) -> Result<bool, VerbError> {
    let is_cram = cram::is_cram(source);
    match &spec.reference {
        Some(reference) if !Path::new(reference).is_file() => {
            Err(VerbError::input(VERB, format!("{reference} not found")).from("-R"))
        }
        None if is_cram => Err(usage(format!(
            "{} is a CRAM: give the reference it was written with as -R FASTA",
            source.display()
        ))),
        _ => Ok(is_cram),
    }
}

/// The chain's prefix: `-o`, with the source name as the file part when
/// it names a directory, else `./` plus the source name.
fn prefix(spec: &Spec, source: &str) -> Result<(PathBuf, PrefixOrigin), VerbError> {
    match spec.prefix.as_deref() {
        // It would name hidden files: `.insert_size_metrics`.
        Some("") => Err(usage("the output prefix is empty")
            .from("-o")
            .hint(EMPTY_PREFIX_HINT)),
        Some("-") => Err(
            usage("standard output takes one file; a chain writes several")
                .from("-o")
                .hint("a single-file step takes -o -, and write - writes the BAM there"),
        ),
        Some(text) => {
            let (path, directory) = prefix_under(text, source);
            let origin = if directory {
                PrefixOrigin::Directory
            } else {
                PrefixOrigin::Option
            };
            Ok((path, origin))
        }
        None => Ok((
            Path::new(".").join(verbs::stem(source)),
            PrefixOrigin::SourceName,
        )),
    }
}

/// The chain's `--format` list: `picard` first.
fn chain_formats(spec: &Spec) -> Result<Vec<Format>, VerbError> {
    let Some(value) = &spec.formats else {
        return Ok(vec![Format::Picard]);
    };
    let formats = format::parse_list(VERB, value, "--format")?;
    picard_first(&formats, value, "--format")?;
    Ok(formats)
}

/// The formats of `--also`. `junit` and `multiqc` print a verdict, which a
/// chain does not have.
fn also(spec: &Spec) -> Result<Vec<Also>, VerbError> {
    let Some(value) = &spec.also else {
        return Ok(Vec::new());
    };
    let mut formats = Vec::new();
    for name in value.split(',').map(str::trim) {
        let format = Also::from_str(name, true).map_err(|_| {
            let names: Vec<&str> = Also::value_variants().iter().map(|a| a.name()).collect();
            let mut error = usage(format!("no format {name:?}")).from("--also");
            if let Some(best) = bammetrics_core::registry::did_you_mean(name, names) {
                error = error.hint(format!("did you mean {best}?"));
            }
            error.hint("formats: json, jsonl, tsv, csv, markdown")
        })?;
        if matches!(format, Also::Junit | Also::Multiqc) {
            return Err(usage(format!(
                "--also {} prints a verdict, and a chain has none",
                format.name()
            ))
            .from("--also")
            .hint(format!(
                "bammetrics check DIR --format {} reads the files a chain wrote",
                format.name()
            ))
            .hint("qc as the last step of a chain is wave W14 of docs/07"));
        }
        if !formats.contains(&format) {
            formats.push(format);
        }
    }
    Ok(formats)
}

/// Prepare: everything checked and named, nothing read or written.
pub fn prepare(spec: &Spec) -> Result<Chain, VerbError> {
    let (source, source_path) = source(spec)?;
    let is_cram = reference(spec, &source_path)?;
    let (prefix, prefix_origin) = prefix(spec, &source)?;
    let formats = chain_formats(spec)?;
    let also = also(spec)?;
    if spec.steps.is_empty() {
        return Err(usage("no step: a chain needs at least one + STEP")
            .hint("bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]")
            .hint("the steps: bammetrics chain --help"));
    }
    let drafts = spec
        .steps
        .iter()
        .enumerate()
        .map(|(index, step)| draft(step, index + 1))
        .collect::<Result<Vec<_>, _>>()?;
    let positions = positions(&drafts)?;
    let last_stream = positions
        .iter()
        .rev()
        .find_map(|position| position.leaves.clone())
        .unwrap_or_else(|| SOURCE_STREAM.to_owned());
    let per_file: Vec<Format> = also.iter().filter_map(|a| a.per_file()).collect();
    // 0 counts as 1, as on every verb. A `Spec` a caller built is checked
    // here; the command line and a plan file checked theirs with its place.
    let threads = spec.threads.unwrap_or_else(pass::default_threads);
    let threads = pass::check_threads(threads)
        .map_err(|e| usage(format!("{threads} threads is {e}")))?
        .max(1);
    let cx = Context {
        source: &source,
        is_cram,
        reference: spec.reference.as_deref(),
        prefix: &prefix,
        last_stream: &last_stream,
        formats: &formats,
        also: &per_file,
        threads,
    };
    let mut steps = Vec::with_capacity(drafts.len());
    let mut prefix_used = false;
    for (draft, position) in drafts.into_iter().zip(positions) {
        let (step, uses_prefix) = match draft.what {
            What::Tool(verb, role) => tool_step(&cx, draft, (verb, role), position)?,
            What::Write => write_step(&cx, &draft, position)?,
        };
        prefix_used |= uses_prefix;
        steps.push(step);
    }
    let report = also
        .contains(&Also::Markdown)
        .then(|| with_suffix(&prefix, Also::Markdown.run_suffix()));
    prefix_used |= report.is_some();
    check_one_stdout(&steps)?;
    refuse_unresolved(&steps, report.as_deref())?;
    let outputs = claim_outputs(&steps, report.as_deref())?;
    let reference = spec.reference.as_ref().map(PathBuf::from);
    let read = Read {
        source: &source_path,
        reference: reference.as_deref(),
        plan_file: spec.plan_file.as_deref().map(Path::new),
    };
    let inputs = inputs_of(&read, &steps);
    refuse_overwriting_inputs(&inputs, &steps, report.as_deref())?;
    Ok(Chain {
        source,
        source_path,
        reference,
        prefix,
        prefix_origin,
        prefix_used,
        threads,
        tmp_dir: spec.tmp_dir.as_ref().map(PathBuf::from),
        sample: spec.sample.clone(),
        quiet: spec.quiet,
        verbose: spec.verbose,
        stage_executor: spec.stage_executor,
        steps,
        report,
        outputs,
        inputs,
    })
}

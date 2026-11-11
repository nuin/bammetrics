//! The Python faces of `bammetrics chain` (docs/python.md): `chain`,
//! `step`, and the `Step` and `Stage` classes. A call runs the command's
//! own code: the keywords spell the argv of the command
//! ([`bammetrics_cli::chain::python`]), the command's parser and
//! `prepare` check it, and `execute` writes the files the command writes.

use std::hash::{DefaultHasher, Hash, Hasher};

use bammetrics_cli::chain::python::{self as table, Keyword, Many};
use bammetrics_cli::chain::{self as cli, Asked, Outcome, StepStatus, VERB};
use bammetrics_cli::error::VerbError;
use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_core::registry::Ty;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyFloat, PyInt, PyList, PyString, PyTuple};

use crate::doc_dict;
use crate::verbs::{self, Run, input_error, usage_error};

/// True for a list or a tuple: several values of one keyword.
fn is_listed(value: &Bound<'_, PyAny>) -> bool {
    value.is_instance_of::<PyList>() || value.is_instance_of::<PyTuple>()
}

/// One value of `keyword` as the text the command takes. `label` is empty
/// for a chain option and `step VERB: ` for a step's; `listed` says the
/// value is an item of a list. A value of the wrong type is a `UsageError`
/// naming the keyword: a float that is NaN or infinite, and an item of a
/// list of texts that is not a str (one number alone stands for its
/// text). A path takes a str or an `os.PathLike` only (`TypeError`, as
/// `qc`).
fn value_text(
    label: &str,
    keyword: &Keyword,
    value: &Bound<'_, PyAny>,
    listed: bool,
) -> PyResult<String> {
    let is_bool = value.is_instance_of::<PyBool>();
    let is_int = value.is_instance_of::<PyInt>() && !is_bool;
    let is_float = value.is_instance_of::<PyFloat>();
    let is_number = is_int || is_float;
    let is_text = value.is_instance_of::<PyString>();
    let finite = !is_float || value.extract::<f64>()?.is_finite();
    let wanted = match keyword.ty {
        Ty::Bool if is_bool => None,
        Ty::Bool => Some("True or False"),
        Ty::Int if is_int => None,
        Ty::Int => Some("an int"),
        Ty::Float | Ty::Text if !finite => Some("a finite number"),
        Ty::Float if is_number => None,
        Ty::Float => Some("a number"),
        Ty::Text if listed && !is_text => Some("a list of str"),
        Ty::Text if is_bool || is_number || is_text => None,
        Ty::Text => Some("a str or a number"),
        Ty::Path if is_text || value.hasattr("__fspath__")? => None,
        Ty::Path => {
            return Err(PyTypeError::new_err(format!(
                "bammetrics {VERB}: {label}{} takes a str or os.PathLike, not {}",
                keyword.name,
                value.get_type().name()?
            )));
        }
    };
    if let Some(wanted) = wanted {
        // `nan` and `inf` say more than `float` does.
        let given = if finite {
            value.get_type().name()?.to_string()
        } else {
            value.repr()?.to_string()
        };
        let holding = if listed && finite { "one holding " } else { "" };
        return Err(usage_error(format!(
            "bammetrics {VERB}: {label}{} takes {wanted}, not {holding}{given}",
            keyword.name
        )));
    }
    let text = if is_bool {
        value.extract::<bool>()?.to_string()
    } else {
        verbs::text(value)?
    };
    if text.contains('\0') {
        return Err(usage_error(format!(
            "bammetrics {VERB}: {label}{} cannot hold a NUL byte",
            keyword.name
        )));
    }
    Ok(text)
}

/// The arguments `keyword=value` spells: `--long=VALUE`, once per item of
/// a list for a repeatable option, once with the items joined by commas
/// for `format` and `also`.
fn keyword_args(label: &str, keyword: &Keyword, value: &Bound<'_, PyAny>) -> PyResult<Vec<String>> {
    let texts = if is_listed(value) {
        if keyword.many == Many::One {
            return Err(usage_error(format!(
                "bammetrics {VERB}: {label}{} takes one value, not a {}",
                keyword.name,
                value.get_type().name()?
            )));
        }
        value
            .try_iter()?
            .map(|item| value_text(label, keyword, &item?, true))
            .collect::<PyResult<Vec<String>>>()?
    } else {
        vec![value_text(label, keyword, value, false)?]
    };
    let flag = |text: &str| format!("--{}={text}", keyword.long);
    Ok(match keyword.many {
        Many::Joined if texts.is_empty() => Vec::new(),
        Many::Joined => vec![flag(&texts.join(","))],
        Many::One | Many::Repeated => texts.iter().map(|text| flag(text)).collect(),
    })
}

/// One step of a chain: a verb and its options, checked when it is made.
#[pyclass(module = "bammetrics", frozen)]
pub struct Step {
    /// The verb as the command spells it: `insert-size`, or `write`.
    #[pyo3(get)]
    verb: String,
    /// `transform`, `tap` or `write`.
    #[pyo3(get)]
    role: String,
    /// The options as the command takes them: `["--minimum-pct=0.1"]`.
    #[pyo3(get)]
    args: Vec<String>,
    /// The keywords as given, without those that were `None`; a list is
    /// kept as a tuple.
    options: Py<PyDict>,
}

/// A [`Step`] as pickle stores it: the verb and the keywords.
type StepState<'py> = (String, Bound<'py, PyDict>);

impl Step {
    /// `name` and `options` checked against the verb's keywords.
    fn build(py: Python<'_>, name: &str, options: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        let table = table::step_table(name).map_err(|e| verbs::verb_error(&e))?;
        let label = format!("step {}: ", table.verb);
        let kept = PyDict::new(py);
        let mut args = Vec::new();
        for (key, value) in options.into_iter().flat_map(|o| o.iter()) {
            let name: String = key.extract()?;
            let keyword = table.keyword(&name).map_err(|e| verbs::verb_error(&e))?;
            if value.is_none() {
                continue;
            }
            args.extend(keyword_args(&label, keyword, &value)?);
            if value.is_instance_of::<PyList>() {
                kept.set_item(
                    key,
                    PyTuple::new(py, value.try_iter()?.collect::<PyResult<Vec<_>>>()?)?,
                )?;
            } else {
                kept.set_item(key, value)?;
            }
        }
        Ok(Step {
            verb: table.verb.to_owned(),
            role: table.role.name().to_owned(),
            args,
            options: kept.unbind(),
        })
    }
}

#[pymethods]
impl Step {
    #[new]
    #[pyo3(signature = (verb, /, **options))]
    fn new(py: Python<'_>, verb: &str, options: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {
        Step::build(py, verb, options)
    }

    /// What pickle and copy rebuild from: the verb and the keywords.
    #[staticmethod]
    fn _restore(py: Python<'_>, verb: &str, options: &Bound<'_, PyDict>) -> PyResult<Self> {
        Step::build(py, verb, Some(options))
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, StepState<'py>)> {
        let step = slf.get();
        Ok((
            slf.get_type().getattr("_restore")?,
            (step.verb.clone(), step.options.bind(slf.py()).copy()?),
        ))
    }

    /// The keywords the step was made with, as a new dict.
    #[getter]
    fn options<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.options.bind(py).copy()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<Step>()
            .is_ok_and(|other| other.get().verb == self.verb && other.get().args == self.args)
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        (&self.verb, &self.args).hash(&mut hasher);
        hasher.finish()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let mut text = format!("Step({:?}", self.verb);
        for (key, value) in self.options.bind(py).iter() {
            text.push_str(&format!(", {}={}", key.str()?, value.repr()?));
        }
        text.push(')');
        Ok(text)
    }
}

/// One step of chain(): a tool verb in snake case (insert_size), its Picard
/// class, or write. The keywords are those of the tool's function
/// (minimum_pct=0.1) without its input, plus output (the step's own file,
/// or prefix for a tool that writes several), format, panel on hs_metrics
/// (baits and targets in one file), and as_ on a transform (the name of the
/// stream it leaves). write takes output, index, md5, compression_level and
/// output_format ("bam" or "sam").
/// An unknown verb or keyword raises UsageError here, with the closest
/// name; so does a float that is NaN or infinite, and a list of texts with
/// an item that is not a str.
//
// The docstring above is the stub's (tests/python/test_mirror.py compares
// them).
#[pyfunction]
#[pyo3(name = "step", signature = (verb, /, **options), text_signature = "(verb, /, **options)")]
pub fn py_step(py: Python<'_>, verb: &str, options: Option<&Bound<'_, PyDict>>) -> PyResult<Step> {
    Step::build(py, verb, options)
}

/// One step of a chain after the run.
#[pyclass(module = "bammetrics", get_all, frozen)]
#[derive(Clone)]
pub struct Stage {
    /// Its place in the chain, from 1.
    number: usize,
    /// The verb as the command spells it: `insert-size`, or `write`.
    verb: String,
    /// `transform`, `tap` or `write`.
    role: String,
    /// `ok`, `failed`, `stopped` (another step of its segment failed) or
    /// `not_run` (a step in front of it failed).
    status: String,
    /// Why a step that is not `ok` ended as it did.
    error: Option<String>,
    /// The stream it reads: `in`, or the name a transform in front left.
    reads: String,
    /// The stream a transform leaves (`marked`); `None` otherwise.
    leaves: Option<String>,
    /// The segment it ran in, from 1.
    segment: usize,
    /// The files it wrote, in write order.
    files: Vec<String>,
    /// A transform's record counts; `None` for a tap and for `write`.
    records_in: Option<u64>,
    records_out: Option<u64>,
    changed: Option<u64>,
}

impl Stage {
    /// `step 2 (insert-size)`, as every message of the command names it.
    fn label(&self) -> String {
        format!("step {} ({})", self.number, self.verb)
    }

    fn cause(&self) -> &str {
        self.error.as_deref().unwrap_or_default()
    }
}

/// The fields of a [`Stage`] in declaration order, as pickle stores them.
type StageState = (
    usize,
    String,
    String,
    String,
    Option<String>,
    String,
    Option<String>,
    usize,
    Vec<String>,
    Option<u64>,
    Option<u64>,
    Option<u64>,
);

#[pymethods]
impl Stage {
    #[staticmethod]
    fn _restore(state: StageState) -> Self {
        let (
            number,
            verb,
            role,
            status,
            error,
            reads,
            leaves,
            segment,
            files,
            records_in,
            records_out,
            changed,
        ) = state;
        Stage {
            number,
            verb,
            role,
            status,
            error,
            reads,
            leaves,
            segment,
            files,
            records_in,
            records_out,
            changed,
        }
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, (StageState,))> {
        let s = slf.get().clone();
        let state = (
            s.number,
            s.verb,
            s.role,
            s.status,
            s.error,
            s.reads,
            s.leaves,
            s.segment,
            s.files,
            s.records_in,
            s.records_out,
            s.changed,
        );
        Ok((slf.get_type().getattr("_restore")?, (state,)))
    }

    fn __repr__(&self) -> String {
        let mut text = format!(
            "Stage(number={}, verb={:?}, role={:?}, status={:?}",
            self.number, self.verb, self.role, self.status
        );
        if let (Some(read), Some(out), Some(changed)) =
            (self.records_in, self.records_out, self.changed)
        {
            text.push_str(&format!(
                ", records_in={read}, records_out={out}, changed={changed}"
            ));
        }
        text.push(')');
        text
    }
}

/// How a step ended, as `Stage.status` names it.
fn status_name(status: StepStatus) -> &'static str {
    match status {
        StepStatus::Ok => "ok",
        StepStatus::Failed => "failed",
        StepStatus::Stopped => "stopped",
        StepStatus::NotRun => "not_run",
    }
}

/// A finished chain as plain values: everything the `Run` holds, made
/// without the GIL.
struct Finished {
    sample: String,
    prefix: String,
    /// Every file written: each step's in step order, then the chain's.
    files: Vec<String>,
    code: u8,
    stages: Vec<Stage>,
    /// `(key, parsed file)` of every metrics file written.
    metrics: Vec<(String, MetricsDoc)>,
}

/// What a call asked for, done.
enum Done {
    /// `dry_run=True`: the plan text.
    Plan(String),
    Ran(Finished),
}

/// `key` when no file has it, else `key#N` with the step's number: two
/// taps of one tool on one stream keep both files.
fn free_key(taken: &[(String, MetricsDoc)], key: String, number: usize) -> String {
    if taken.iter().any(|(other, _)| *other == key) {
        format!("{key}#{number}")
    } else {
        key
    }
}

/// The values of a finished run, with every metrics file read back.
fn finished(outcome: &Outcome) -> Result<Finished, VerbError> {
    let chain = &outcome.chain;
    let mut stages = Vec::with_capacity(outcome.steps.len());
    let mut metrics = Vec::new();
    let mut files = Vec::new();
    for (step, planned) in outcome.steps.iter().zip(&chain.steps) {
        let written: Vec<String> = step
            .files
            .iter()
            .map(|file| file.path.display().to_string())
            .collect();
        for file in step.files.iter().filter(|file| file.kind == "metrics") {
            let doc = MetricsDoc::read(&file.path).map_err(|e| VerbError::core(VERB, &e))?;
            let key = table::metric_key(chain, planned, &file.path);
            metrics.push((free_key(&metrics, key, step.number), doc));
        }
        files.extend(written.iter().cloned());
        stages.push(Stage {
            number: step.number,
            verb: planned.verb.to_owned(),
            role: planned.role.name().to_owned(),
            status: status_name(step.status).to_owned(),
            error: step.error.clone(),
            reads: planned.reads.clone(),
            leaves: planned.leaves.clone(),
            segment: planned.segment,
            files: written,
            records_in: step.counters.map(|c| c.records_in),
            records_out: step.counters.map(|c| c.records_out),
            changed: step.counters.map(|c| c.changed),
        });
    }
    files.extend(outcome.files.iter().map(|f| f.path.display().to_string()));
    Ok(Finished {
        sample: cli::sample_name(chain),
        prefix: chain.prefix.display().to_string(),
        files,
        code: outcome.code,
        stages,
        metrics,
    })
}

/// The command on `argv`: parse, prepare, then the plan text or the run.
fn run(argv: &[String]) -> Result<Done, VerbError> {
    let (spec, dry_run) = match cli::parse(argv)? {
        Asked::Chain { spec, dry_run } => (spec, dry_run),
        // No keyword spells `--help`.
        Asked::Help(text) => return Ok(Done::Plan(text)),
    };
    let chain = cli::prepare(&spec)?;
    if dry_run {
        return cli::plan_text(&chain).map(Done::Plan);
    }
    let outcome = cli::execute(chain)?;
    finished(&outcome).map(Done::Ran)
}

/// The steps of a call as `Step`s: each item a `Step` or a verb name.
fn steps_of<'py>(
    py: Python<'py>,
    steps: Option<&Bound<'py, PyAny>>,
) -> PyResult<Vec<Bound<'py, Step>>> {
    let Some(steps) = steps.filter(|steps| !steps.is_none()) else {
        return Ok(Vec::new());
    };
    let items: Vec<Bound<'py, PyAny>> = if is_listed(steps) {
        steps.try_iter()?.collect::<PyResult<_>>()?
    } else {
        vec![steps.clone()]
    };
    items
        .into_iter()
        .map(|item| {
            if let Ok(step) = item.cast::<Step>() {
                return Ok(step.clone());
            }
            if let Ok(name) = item.cast::<PyString>() {
                return Bound::new(py, Step::build(py, name.to_str()?, None)?);
            }
            Err(PyTypeError::new_err(format!(
                "bammetrics {VERB}: steps takes verb names and step(...) objects, not {}",
                item.get_type().name()?
            )))
        })
        .collect()
}

/// The argv the call spells: the chain's options, the source, then
/// `+ VERB OPTIONS` per step.
fn argv<'py>(
    py: Python<'py>,
    source: Option<&Bound<'py, PyAny>>,
    steps: Option<&Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Vec<String>> {
    let mut argv = Vec::new();
    for (key, value) in options.into_iter().flat_map(|o| o.iter()) {
        let name: String = key.extract()?;
        let keyword = table::chain_keyword(&name).map_err(|e| verbs::verb_error(&e))?;
        if !value.is_none() {
            argv.extend(keyword_args("", &keyword, &value)?);
        }
    }
    if let Some(source) = source.filter(|source| !source.is_none()) {
        argv.extend(keyword_args("", &table::source_keyword(), source)?);
    }
    for step in steps_of(py, steps)? {
        let step = step.get();
        argv.push(cli::STEP.to_owned());
        argv.push(step.verb.clone());
        argv.extend(step.args.iter().cloned());
    }
    Ok(argv)
}

/// `InputError` for a run with a step that did not finish, with `.run`
/// set. The message names each step that failed by itself with its cause,
/// then counts the steps it halted; `Run.errors` holds every message.
fn failed_steps(py: Python<'_>, run: Py<Run>, done: &Finished) -> PyErr {
    let with_status = |status: StepStatus| {
        let name = status_name(status);
        done.stages.iter().filter(move |stage| stage.status == name)
    };
    let failed: Vec<String> = with_status(StepStatus::Failed)
        .map(|stage| format!("{}: {}", stage.label(), stage.cause()))
        .collect();
    let mut message = format!(
        "bammetrics {VERB}: {} step(s) failed: {}",
        failed.len(),
        failed.join("; ")
    );
    for (status, words) in [
        (StepStatus::Stopped, "stopped"),
        (StepStatus::NotRun, "not run"),
    ] {
        let count = with_status(status).count();
        if count > 0 {
            message.push_str(&format!("; {count} {words}"));
        }
    }
    let err = input_error(message);
    if let Err(e) = err.value(py).setattr("run", run) {
        return e;
    }
    err
}

/// The `Run` of a finished chain: no manifest and no verdict.
fn run_object(py: Python<'_>, done: &Finished) -> PyResult<Py<Run>> {
    let errors = PyDict::new(py);
    for stage in done.stages.iter().filter(|stage| stage.error.is_some()) {
        errors.set_item(stage.label(), stage.cause())?;
    }
    let metrics = PyDict::new(py);
    for (key, doc) in &done.metrics {
        metrics.set_item(key, doc_dict(py, doc)?)?;
    }
    let stages = done
        .stages
        .iter()
        .map(|stage| Py::new(py, stage.clone()))
        .collect::<PyResult<_>>()?;
    Py::new(
        py,
        Run {
            sample: done.sample.clone(),
            prefix: done.prefix.clone(),
            files: done.files.clone(),
            manifest: None,
            exit_code: done.code,
            errors: errors.unbind(),
            metrics: metrics.unbind(),
            verdict: None,
            stages,
            inputs: Vec::new(),
            rust_verdict: None,
            all_plots: false,
            date: String::new(),
        },
    )
}

/// bammetrics chain: several steps on one BAM in one process, the records
/// handed from step to step in memory. steps holds verb names (insert_size)
/// or step(...) objects, in order: mark_duplicates and sort_sam are
/// transforms, a one-pass tool is a tap on the records of its position,
/// write writes them as a BAM, or as SAM text for a .sam name or
/// output_format "sam". plan names a TOML plan file in place of steps; a keyword given
/// replaces the file's. Returns a Run: files, metrics (keyed by verb,
/// STREAM.verb in front of a transform) and stages (one Stage per step).
/// dry_run=True returns the plan text and runs nothing. A failed step raises
/// InputError after the other steps wrote their files (its .run holds the
/// Run).
//
// Runs `chain::prepare` and `chain::execute` without the GIL. The text
// signature and the docstring above are the stub's
// (tests/python/test_mirror.py compares them).
#[pyfunction]
#[pyo3(
    name = "chain",
    signature = (source = None, steps = None, **options),
    text_signature = "(source=None, steps=None, *, output=None, reference=None, threads=None, \
                      tmp_dir=None, format=None, also=None, sample=None, sort_memory=None, \
                      force_space=None, validation_stringency=None, plan=None, dry_run=False)"
)]
pub fn py_chain<'py>(
    py: Python<'py>,
    source: Option<&Bound<'py, PyAny>>,
    steps: Option<&Bound<'py, PyAny>>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Py<PyAny>> {
    let argv = argv(py, source, steps, options)?;
    let done = py
        .detach(|| run(&argv))
        .map_err(|e| verbs::verb_error(&e))?;
    let done = match done {
        Done::Plan(text) => return Ok(PyString::new(py, &text).into_any().unbind()),
        Done::Ran(done) => done,
    };
    let run = run_object(py, &done)?;
    if done.code != 0 {
        return Err(failed_steps(py, run, &done));
    }
    Ok(run.into_any())
}

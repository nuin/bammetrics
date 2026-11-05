//! The Python faces of `bammetrics qc`, `check`, `report`, `tools` and
//! `presets`, and the `collect(metrics=[...])` form (docs/05 section 5).
//! Every function runs the command's own code in `bammetrics-cli`:
//! keywords become the command's argv ([`bammetrics_cli::python`]), the
//! command's parser checks them, and the engine writes the same files.

use std::ffi::CString;
use std::path::{Path, PathBuf};

use bammetrics_cli::check::{self, CheckArgs, FailOn};
use bammetrics_cli::engine::{self, Outcome};
use bammetrics_cli::error::{Error as CliError, VerbError};
use bammetrics_cli::manifest::{Output, Status as StepStatus, StepRecord};
use bammetrics_cli::python::{self as table, Keyword, Shape};
use bammetrics_cli::qc::{self, QcArgs};
use bammetrics_cli::report::{self, ReportArgs};
use bammetrics_cli::{info, verbs};
use bammetrics_core::metrics_file::{CivilTime, MetricsDoc};
use bammetrics_core::registry::did_you_mean;
use bammetrics_report::{Input, Options};
use bammetrics_verdict::{Status, Verdict as RustVerdict};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyBool, PyDict, PyList, PyString, PyTuple, PyType};

use crate::doc_dict;

/// The classes `types.py` defines, kept for raising and building.
struct Types {
    usage: Py<PyType>,
    input: Py<PyType>,
    verdict: Py<PyType>,
    tool: Py<PyType>,
    preset: Py<PyType>,
}

static TYPES: PyOnceLock<Types> = PyOnceLock::new();

/// The class names `types.py` defines, all exported.
const TYPE_NAMES: [&str; 6] = [
    "BammetricsError",
    "UsageError",
    "InputError",
    "VerdictError",
    "Tool",
    "PresetInfo",
];

/// Run `types.py` and add its classes to `m`, as members of `bammetrics`.
pub fn add_types(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    let code = CString::new(include_str!("types.py")).expect("no NUL in types.py");
    let types = PyModule::from_code(py, &code, c"bammetrics/_types.py", c"bammetrics._types")?;
    for name in TYPE_NAMES {
        let class = types.getattr(name)?;
        class.setattr("__module__", "bammetrics")?;
        m.add(name, class)?;
    }
    let class = |name: &str| -> PyResult<Py<PyType>> {
        Ok(types.getattr(name)?.cast_into::<PyType>()?.unbind())
    };
    let _ = TYPES.set(
        py,
        Types {
            usage: class("UsageError")?,
            input: class("InputError")?,
            verdict: class("VerdictError")?,
            tool: class("Tool")?,
            preset: class("PresetInfo")?,
        },
    );
    Ok(())
}

fn types(py: Python<'_>) -> &Types {
    TYPES.get(py).expect("add_types ran at import")
}

/// `UsageError(message)`.
pub fn usage_error(message: impl Into<String>) -> PyErr {
    let message = message.into();
    Python::attach(|py| PyErr::from_type(types(py).usage.bind(py).clone(), message))
}

/// `InputError(message)`.
pub fn input_error(message: impl Into<String>) -> PyErr {
    let message = message.into();
    Python::attach(|py| PyErr::from_type(types(py).input.bind(py).clone(), message))
}

/// A command's error as the exception of its exit code: 2 is a
/// `UsageError`, 1 an `InputError`. The message is the line the command
/// prints, `bammetrics qc: ...`.
fn verb_error(e: &VerbError) -> PyErr {
    if e.code == 2 {
        usage_error(e.to_string())
    } else {
        input_error(e.to_string())
    }
}

fn cli_error(verb: &str, e: &CliError) -> PyErr {
    let message = match e {
        CliError::Verb(inner) => return verb_error(inner),
        other => format!("bammetrics {verb}: {other}"),
    };
    if e.exit_code() == 2 {
        usage_error(message)
    } else {
        input_error(message)
    }
}

/// A keyword value as argv text: `os.fspath` for a path, `str` otherwise.
fn text(value: &Bound<'_, PyAny>) -> PyResult<String> {
    if value.is_instance_of::<PyString>() {
        return value.extract();
    }
    if value.hasattr("__fspath__")? {
        let path: PathBuf = value.extract()?;
        return Ok(path.to_string_lossy().into_owned());
    }
    Ok(value.str()?.to_string())
}

/// True for a str or a path: one value, not a sequence of them.
fn is_scalar(value: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(value.is_instance_of::<PyString>() || value.hasattr("__fspath__")?)
}

/// The items of a list or tuple, or `[value]`.
fn items<'py>(value: &Bound<'py, PyAny>) -> PyResult<Vec<Bound<'py, PyAny>>> {
    if !is_scalar(value)? && (value.is_instance_of::<PyList>() || value.is_instance_of::<PyTuple>())
    {
        return value.try_iter()?.collect();
    }
    Ok(vec![value.clone()])
}

/// Python-only keywords (`report`) by name, with their values.
type Extra<'py> = Vec<(String, Bound<'py, PyAny>)>;

/// The argv `options` spell for `verb`: `--long=VALUE`, `--flag`, `--hs A
/// B`, then `--` and the positional values, and the Python-only keywords.
/// `allowed`, when given, limits the keywords.
fn argv<'py>(
    verb: &table::Verb,
    positional: &Bound<'py, PyAny>,
    options: Option<&Bound<'py, PyDict>>,
    allowed: Option<&[&str]>,
) -> PyResult<(Vec<String>, Extra<'py>)> {
    let keywords = verb.keywords();
    let mut args = Vec::new();
    let mut extra = Vec::new();
    for (key, value) in options.into_iter().flat_map(|o| o.iter()) {
        let name: String = key.extract()?;
        let known = keywords
            .iter()
            .filter(|k| !matches!(k.shape, Shape::Positional { .. }))
            .filter(|k| allowed.is_none_or(|a| a.contains(&k.name.as_str())));
        let found: Vec<&Keyword> = known.collect();
        let Some(keyword) = found.iter().find(|k| k.name == name) else {
            let mut message = format!("bammetrics {}: unknown keyword {name:?}", verb.name);
            if let Some(best) = did_you_mean(&name, found.iter().map(|k| k.name.as_str())) {
                message.push_str(&format!("; did you mean {best}?"));
            }
            return Err(usage_error(message));
        };
        if value.is_none() {
            continue;
        }
        if keyword.long.is_empty() {
            extra.push((name, value));
            continue;
        }
        let flag = format!("--{}", keyword.long);
        match keyword.shape {
            Shape::Flag => {
                if !value.is_instance_of::<PyBool>() {
                    return Err(usage_error(format!(
                        "bammetrics {}: {name} takes True or False",
                        verb.name
                    )));
                }
                if value.extract::<bool>()? {
                    args.push(flag);
                }
            }
            Shape::Value => args.push(format!("{flag}={}", text(&value)?)),
            Shape::Pairs => {
                let first = items(&value)?;
                let pairs = match first.first() {
                    Some(item) if is_scalar(item)? => vec![value.clone()],
                    _ => first,
                };
                for pair in pairs {
                    let sides = items(&pair)?;
                    let [a, b] = sides.as_slice() else {
                        return Err(usage_error(format!(
                            "bammetrics {}: {name} takes (baits, targets) pairs",
                            verb.name
                        )));
                    };
                    args.extend([flag.clone(), text(a)?, text(b)?]);
                }
            }
            Shape::Positional { .. } => {}
        }
    }
    args.push("--".to_owned());
    for value in items(positional)? {
        args.push(text(&value)?);
    }
    Ok((args, extra))
}

/// One check of a verdict.
#[pyclass(module = "bammetrics", get_all, frozen)]
pub struct Check {
    rule: String,
    /// `pass`, `warn`, `fail` or `not_measured`.
    status: String,
    value: Option<f64>,
    /// The bound the status reports: the fail bound of a fail, else the
    /// warn bound (the fail bound when the rule has no warn bound).
    limit: Option<String>,
    /// The reason in plain words; `None` for a pass.
    why: Option<String>,
    label: String,
    metric: String,
    level: String,
    unit: Option<String>,
    warn: Option<String>,
    fail: Option<String>,
}

#[pymethods]
impl Check {
    fn __repr__(&self) -> String {
        format!(
            "Check(rule={:?}, status={:?}, value={}, limit={})",
            self.rule,
            self.status,
            self.value
                .map_or_else(|| "None".to_owned(), |v| v.to_string()),
            self.limit
                .as_ref()
                .map_or_else(|| "None".to_owned(), |l| format!("{l:?}"))
        )
    }
}

/// A sample's verdict, as `PREFIX.verdict.json` holds it.
#[pyclass(module = "bammetrics", frozen)]
pub struct Verdict {
    inner: RustVerdict,
}

#[pymethods]
impl Verdict {
    #[getter]
    fn sample(&self) -> &str {
        &self.inner.sample
    }

    /// `pass`, `warn`, `fail` or `incomplete`.
    #[getter]
    fn status(&self) -> &'static str {
        self.inner.status.name()
    }

    #[getter]
    fn checks(&self) -> Vec<Check> {
        self.inner.checks.iter().map(check_of).collect()
    }

    /// `{"pass": n, "warn": n, "fail": n, "not_measured": n}`.
    #[getter]
    fn counts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let c = self.inner.counts;
        let dict = PyDict::new(py);
        dict.set_item("pass", c.pass)?;
        dict.set_item("warn", c.warn)?;
        dict.set_item("fail", c.fail)?;
        dict.set_item("not_measured", c.not_measured)?;
        Ok(dict)
    }

    /// `{"name", "version", "hash", "source", "note"}`.
    #[getter]
    fn preset<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let p = &self.inner.preset;
        let dict = PyDict::new(py);
        dict.set_item("name", &p.name)?;
        dict.set_item("version", p.version)?;
        dict.set_item("hash", &p.hash)?;
        dict.set_item("source", &p.source)?;
        dict.set_item("note", &p.note)?;
        Ok(dict)
    }

    /// The text of `verdict.json`.
    fn to_json(&self) -> String {
        self.inner.to_json()
    }

    /// `verdict.json` loaded with the json module, as `evaluate` returns it.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        py.import("json")?
            .call_method1("loads", (self.inner.to_json(),))
    }

    /// The exit code `bammetrics check --fail-on` gives this verdict: 0, 10
    /// or 11.
    #[pyo3(signature = (fail_on = "fail"))]
    fn exit_code(&self, fail_on: &str) -> PyResult<u8> {
        Ok(check::exit_code(
            std::slice::from_ref(&self.inner),
            fail_on_of(fail_on)?,
        ))
    }

    fn __repr__(&self) -> String {
        let c = self.inner.counts;
        format!(
            "Verdict(sample={:?}, status={:?}, fail={}, warn={}, pass={}, not_measured={})",
            self.inner.sample,
            self.inner.status.name(),
            c.fail,
            c.warn,
            c.pass,
            c.not_measured
        )
    }
}

fn fail_on_of(text: &str) -> PyResult<FailOn> {
    match text {
        "fail" => Ok(FailOn::Fail),
        "warn" => Ok(FailOn::Warn),
        other => Err(usage_error(format!(
            "fail_on {other:?}: one of \"warn\", \"fail\""
        ))),
    }
}

fn check_of(c: &bammetrics_verdict::Check) -> Check {
    let limit = match c.status {
        Status::Fail => c.fail.clone(),
        _ => c.warn.clone().or_else(|| c.fail.clone()),
    };
    Check {
        rule: c.rule.clone(),
        status: c.status.name().to_owned(),
        value: c.value,
        limit,
        why: c.reason.clone(),
        label: c.label.clone(),
        metric: c.metric.clone(),
        level: c.level.clone(),
        unit: c.unit.clone(),
        warn: c.warn.clone(),
        fail: c.fail.clone(),
    }
}

/// A `qc` run, or a `collect(metrics=[...])` run (no verdict).
#[pyclass(module = "bammetrics", frozen)]
pub struct Run {
    #[pyo3(get)]
    sample: String,
    #[pyo3(get)]
    prefix: String,
    /// Every file the run wrote or kept, the manifest last.
    #[pyo3(get)]
    files: Vec<String>,
    /// `PREFIX.manifest.json`; `None` for a collect run.
    #[pyo3(get)]
    manifest: Option<String>,
    /// The exit code `bammetrics qc` gives the run: 0, 1, 10 or 11.
    #[pyo3(get)]
    exit_code: u8,
    /// `{step: message}` for every step that failed.
    #[pyo3(get)]
    errors: Py<PyDict>,
    metrics: Py<PyDict>,
    verdict: Option<Py<Verdict>>,
    /// What the notebook view renders.
    inputs: Vec<Input>,
    rust_verdict: Option<RustVerdict>,
    all_plots: bool,
    date: String,
}

#[pymethods]
impl Run {
    /// `{name: parsed metrics file}`, each as `read_metrics` returns it,
    /// keyed by the verb of the tool that wrote it (`insert-size`; see
    /// docs/python.md for several hs runs and further files).
    #[getter]
    fn metrics<'py>(&self, py: Python<'py>) -> Bound<'py, PyDict> {
        self.metrics.bind(py).clone()
    }

    #[getter]
    fn verdict<'py>(&self, py: Python<'py>) -> Option<Bound<'py, Verdict>> {
        self.verdict.as_ref().map(|v| v.bind(py).clone())
    }

    /// The report body for a notebook: one `<div class="bm">` with its
    /// stylesheet scoped to the div and inline SVG, series cut to 200
    /// points. `None` for a collect run.
    fn _repr_html_(&self) -> Option<String> {
        let verdict = self.rust_verdict.as_ref()?;
        Some(bammetrics_report::fragment(
            verdict,
            &self.inputs,
            &Options {
                all_plots: self.all_plots,
                date: &self.date,
            },
        ))
    }

    fn __repr__(&self) -> String {
        let status = self
            .rust_verdict
            .as_ref()
            .map_or_else(String::new, |v| format!(", status={:?}", v.status.name()));
        format!(
            "Run(sample={:?}{status}, prefix={:?}, files={})",
            self.sample,
            self.prefix,
            self.files.len()
        )
    }
}

/// The metrics files of the steps that finished, parsed.
fn read_outputs(records: &[StepRecord]) -> Result<Vec<Input>, CliError> {
    metrics_outputs(records)
        .map(|(_, output)| {
            Ok(Input {
                path: path_text(&output.path),
                doc: MetricsDoc::read(&output.path)?,
            })
        })
        .collect()
}

/// Every metrics file of the steps that finished, with its step.
fn metrics_outputs(records: &[StepRecord]) -> impl Iterator<Item = (&StepRecord, &Output)> {
    records
        .iter()
        .filter(|r| r.status == StepStatus::Ok)
        .flat_map(|r| r.outputs.iter().map(move |o| (r, o)))
        .filter(|(_, o)| o.kind == "metrics")
}

/// The parsed metrics files of `records`, keyed by [`table::metric_key`];
/// `inputs` holds them parsed.
fn metrics_of<'py>(
    py: Python<'py>,
    records: &[StepRecord],
    inputs: &[Input],
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for (record, output) in metrics_outputs(records) {
        let shown = path_text(&output.path);
        if let Some(input) = inputs.iter().find(|i| i.path == shown) {
            dict.set_item(table::metric_key(record, output), doc_dict(py, &input.doc)?)?;
        }
    }
    Ok(dict)
}

fn path_text(path: &Path) -> String {
    path.display().to_string()
}

/// The `{step: message}` of the failed steps.
fn errors_of<'py>(py: Python<'py>, records: &[StepRecord]) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for record in records.iter().filter(|r| r.status == StepStatus::Failed) {
        dict.set_item(&record.id, record.error.clone().unwrap_or_default())?;
    }
    Ok(dict)
}

/// `InputError` for the failed steps of a run, with `.run` set.
fn failed_steps(py: Python<'_>, verb: &str, run: Py<Run>, records: &[StepRecord]) -> PyErr {
    let failed: Vec<String> = records
        .iter()
        .filter(|r| r.status == StepStatus::Failed)
        .map(|r| format!("{}: {}", r.id, r.error.clone().unwrap_or_default()))
        .collect();
    let err = input_error(format!(
        "bammetrics {verb}: {} step(s) failed: {}",
        failed.len(),
        failed.join("; ")
    ));
    if let Err(e) = err.value(py).setattr("run", run) {
        return e;
    }
    err
}

fn verdict_error(
    py: Python<'_>,
    message: String,
    verdicts: Vec<Py<Verdict>>,
    run: Option<Py<Run>>,
) -> PyErr {
    match types(py).verdict.bind(py).call1((message, verdicts, run)) {
        Ok(instance) => PyErr::from_value(instance),
        Err(e) => e,
    }
}

/// The Python verb named `name`.
fn verb(name: &str) -> &'static table::Verb {
    table::verb(name).expect("a Python verb")
}

/// `qc(bam, **options)`: `bammetrics qc` through [`engine::execute`],
/// without the GIL. Returns a [`Run`].
#[pyfunction]
#[pyo3(name = "qc", signature = (bam, **options))]
pub fn py_qc<'py>(
    py: Python<'py>,
    bam: &Bound<'py, PyAny>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Py<Run>> {
    let qc_verb = verb("qc");
    let (argv, extra) = argv(qc_verb, bam, options, None)?;
    let mut args: QcArgs = qc_verb.parse(&argv).map_err(|e| verb_error(&e))?;
    for (name, value) in extra {
        if name == "report" {
            if !value.is_instance_of::<PyBool>() {
                return Err(usage_error("bammetrics qc: report takes True or False"));
            }
            args.write_page = value.extract()?;
        }
    }
    args.quiet = true;
    let date = args.date.clone().unwrap_or_else(|| CivilTime::now().iso());
    let outcome = py
        .detach(|| engine::execute(&args))
        .map_err(|e| cli_error("qc", &e))?;
    let Outcome {
        plan,
        verdict,
        inputs,
        records,
        code,
    } = outcome;
    let mut files: Vec<String> = records
        .iter()
        .flat_map(|r| &r.outputs)
        .map(|o| path_text(&o.path))
        .collect();
    let manifest = path_text(&qc::manifest_path(&plan.prefix));
    files.push(manifest.clone());
    let py_verdict = Py::new(
        py,
        Verdict {
            inner: verdict.clone(),
        },
    )?;
    let run = Py::new(
        py,
        Run {
            sample: plan.sample.clone(),
            prefix: path_text(&plan.prefix),
            files,
            manifest: Some(manifest),
            exit_code: code,
            errors: errors_of(py, &records)?.unbind(),
            metrics: metrics_of(py, &records, &inputs)?.unbind(),
            verdict: Some(py_verdict.clone_ref(py)),
            inputs,
            rust_verdict: Some(verdict.clone()),
            all_plots: args.all_plots,
            date,
        },
    )?;
    match code {
        1 => Err(failed_steps(py, "qc", run, &records)),
        10 => Err(verdict_error(
            py,
            format!(
                "bammetrics qc: {} is {}, which reaches fail_on",
                plan.sample,
                verdict.status.name()
            ),
            vec![py_verdict],
            Some(run),
        )),
        _ => Ok(run),
    }
}

/// `check(inputs, **options)`: the verdicts `bammetrics check` gates on.
#[pyfunction]
#[pyo3(name = "check", signature = (inputs, **options))]
pub fn py_check<'py>(
    py: Python<'py>,
    inputs: &Bound<'py, PyAny>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Py<PyAny>> {
    let check_verb = verb("check");
    let gate = match options
        .map(|o| o.get_item("fail_on"))
        .transpose()?
        .flatten()
    {
        Some(value) => !value.is_none(),
        None => true,
    };
    let (argv, _) = argv(check_verb, inputs, options, None)?;
    let mut args: CheckArgs = check_verb.parse(&argv).map_err(|e| verb_error(&e))?;
    args.inputs = args.inputs.iter().map(|p| table::check_input(p)).collect();
    let verdicts = py
        .detach(|| check::verdicts(&args.inputs, &args.preset))
        .map_err(|e| cli_error("check", &e))?;
    let code = check::exit_code(&verdicts, args.fail_on);
    let objects: Vec<Py<Verdict>> = verdicts
        .iter()
        .map(|v| Py::new(py, Verdict { inner: v.clone() }))
        .collect::<PyResult<_>>()?;
    if gate && code == 10 {
        let reached: Vec<&str> = verdicts
            .iter()
            .filter(|v| check::exit_code(std::slice::from_ref(*v), args.fail_on) == 10)
            .map(|v| v.sample.as_str())
            .collect();
        return Err(verdict_error(
            py,
            format!("bammetrics check: {} reach fail_on", reached.join(", ")),
            objects,
            None,
        ));
    }
    if is_scalar(inputs)? && objects.len() == 1 {
        let one = objects.into_iter().next().expect("one verdict");
        return Ok(one.into_any());
    }
    Ok(PyList::new(py, objects)?.into_any().unbind())
}

/// `report(dirs, **options)`: the pages, verdicts and index `bammetrics
/// report` writes, without the GIL. Returns the paths written.
#[pyfunction]
#[pyo3(name = "report", signature = (dirs, **options))]
pub fn py_report(
    py: Python<'_>,
    dirs: &Bound<'_, PyAny>,
    options: Option<&Bound<'_, PyDict>>,
) -> PyResult<Vec<String>> {
    let report_verb = verb("report");
    let (argv, _) = argv(report_verb, dirs, options, None)?;
    let args: ReportArgs = report_verb.parse(&argv).map_err(|e| verb_error(&e))?;
    let written = py
        .detach(|| report::write(&args))
        .map_err(|e| cli_error("report", &e))?;
    Ok(written.files.iter().map(|p| path_text(p)).collect())
}

/// The `qc` keywords `collect(metrics=[...])` takes.
pub const COLLECT_KEYWORDS: [&str; 8] = [
    "output",
    "reference",
    "hs",
    "panel",
    "ref_flat",
    "strand",
    "ribosomal_intervals",
    "threads",
];

/// The Picard class of a `collect` metric name: a tool verb `qc` can run,
/// or `hs` as `bammetrics collect --metric hs` spells it.
fn collector_of(name: &str) -> PyResult<&'static str> {
    let tools = qc::collector_tools();
    let verb_names: Vec<&str> = tools
        .iter()
        .filter_map(|t| verbs::for_class(t).map(|v| v.name))
        .collect();
    let class = match name {
        "hs" => Some("CollectHsMetrics"),
        other => verbs::find(other).map(|v| v.class),
    };
    if let Some(class) = class.filter(|c| tools.contains(c)) {
        return Ok(class);
    }
    let mut message = format!("bammetrics collect: no collector {name:?}");
    if let Some(best) = did_you_mean(name, verb_names.iter().copied()) {
        message.push_str(&format!("; did you mean {best}?"));
    }
    message.push_str(&format!("; one of: {}", verb_names.join(", ")));
    Err(usage_error(message))
}

/// `collect(path, metrics=[verb names], output=..., ...)`: the collectors
/// `qc` would run for those metrics, on one read, their files under the
/// output prefix (default `./<bam stem>`), without the GIL.
pub fn collect_run<'py>(
    py: Python<'py>,
    path: &Bound<'py, PyAny>,
    metrics: Option<Vec<String>>,
    options: &Bound<'py, PyDict>,
) -> PyResult<Py<Run>> {
    let has_hs = options.get_item("hs")?.is_some_and(|v| !v.is_none());
    let names = metrics.unwrap_or_else(|| {
        let mut names = vec!["quality-yield", "alignment-summary", "insert-size"];
        if has_hs {
            names.push("hs-metrics");
        }
        names.into_iter().map(str::to_owned).collect()
    });
    let mut tools: Vec<&str> = Vec::new();
    for name in &names {
        let tool = collector_of(name)?;
        if !tools.contains(&tool) {
            tools.push(tool);
        }
    }
    if tools.contains(&"CollectHsMetrics") != has_hs {
        return Err(usage_error(
            "bammetrics collect: hs-metrics and hs=(baits, targets) go together",
        ));
    }
    let qc_verb = verb("qc");
    let (argv, _) = argv(qc_verb, path, Some(options), Some(&COLLECT_KEYWORDS))?;
    let mut args: QcArgs = qc_verb.parse(&argv).map_err(|e| verb_error(&e))?;
    args.quiet = true;
    let (plan, records, inputs) = py
        .detach(|| {
            let (plan, records) = engine::collect_tools(&args, &tools)?;
            let inputs = read_outputs(&records)?;
            Ok::<_, CliError>((plan, records, inputs))
        })
        .map_err(|e| cli_error("collect", &e))?;
    let files = records
        .iter()
        .flat_map(|r| &r.outputs)
        .map(|o| path_text(&o.path))
        .collect();
    let run = Py::new(
        py,
        Run {
            sample: plan.sample.clone(),
            prefix: path_text(&plan.prefix),
            files,
            manifest: None,
            exit_code: u8::from(records.iter().any(|r| r.status == StepStatus::Failed)),
            errors: errors_of(py, &records)?.unbind(),
            metrics: metrics_of(py, &records, &inputs)?.unbind(),
            verdict: None,
            inputs: Vec::new(),
            rust_verdict: None,
            all_plots: false,
            date: String::new(),
        },
    )?;
    if records.iter().any(|r| r.status == StepStatus::Failed) {
        return Err(failed_steps(py, "collect", run, &records));
    }
    Ok(run)
}

/// The tool verbs as `bammetrics tools` lists them: one `Tool` (a str
/// equal to the Picard class) per verb, with `verb` and `summary`.
#[pyfunction]
#[pyo3(name = "tools")]
pub fn py_tools(py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
    let class = types(py).tool.bind(py);
    verbs::VERBS
        .iter()
        .map(|v| Ok(class.call1((v.class, v.name, v.tool().summary))?.unbind()))
        .collect()
}

/// The shipped presets as `bammetrics presets` lists them: one
/// `PresetInfo` (a str equal to the name) each, with `version`, `hash`,
/// `rules` (the rule count) and `note`.
#[pyfunction]
#[pyo3(name = "presets")]
pub fn py_presets(py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
    let class = types(py).preset.bind(py);
    info::shipped_presets()
        .map_err(|e| verb_error(&e))?
        .iter()
        .map(|p| {
            Ok(class
                .call1((&p.name, p.version, &p.hash, p.rules.len(), &p.note))?
                .unbind())
        })
        .collect()
}

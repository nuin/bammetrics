//! The Python face of `bammetrics doctor` (docs/doctor.md, docs/python.md):
//! `doctor(path, ...)` returns a [`Doctor`] with one [`DoctorCheck`] per
//! check. The keywords are the command's flags, parsed by its parser; the
//! checks run in `bammetrics_cli::doctor::execute` without the GIL.

use bammetrics_cli::doctor::{self as cli_doctor, DoctorArgs};
use bammetrics_cli::python as table;
use bammetrics_core::doctor::Report;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::verbs::{argv, cli_error, verb_error};

/// One check: `name` (as the command prints it), `status` (`pass`, `warn`,
/// `fail` or `skip`), `value` (what was found) and `reason`.
#[pyclass(module = "bammetrics", get_all, frozen)]
#[derive(Clone)]
pub struct DoctorCheck {
    name: String,
    status: String,
    value: String,
    reason: String,
}

type CheckState = (String, String, String, String);

#[pymethods]
impl DoctorCheck {
    #[staticmethod]
    fn _restore(state: CheckState) -> Self {
        let (name, status, value, reason) = state;
        DoctorCheck {
            name,
            status,
            value,
            reason,
        }
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, (CheckState,))> {
        let c = slf.get();
        let state = (
            c.name.clone(),
            c.status.clone(),
            c.value.clone(),
            c.reason.clone(),
        );
        Ok((slf.get_type().getattr("_restore")?, (state,)))
    }

    fn __repr__(&self) -> String {
        format!(
            "DoctorCheck(name={:?}, status={:?}, value={:?})",
            self.name, self.status, self.value
        )
    }
}

/// What `doctor` found in one file: `path`, `kind` (`BAM`, `CRAM 3.1`,
/// `SAM`; empty when the file is none of them), `status` (the worst
/// status of any check: `pass`, `warn` or `fail`), `exit_code` (the
/// command's: 0, or 10 when a check reaches fail_on), `checks`, and
/// `to_json()`, the text of `--format json`.
#[pyclass(module = "bammetrics", get_all, frozen)]
pub struct Doctor {
    path: String,
    kind: String,
    status: String,
    exit_code: u8,
    checks: Vec<DoctorCheck>,
    json: String,
}

type DoctorState = (String, String, String, u8, Vec<CheckState>, String);

#[pymethods]
impl Doctor {
    #[staticmethod]
    fn _restore(state: DoctorState) -> Self {
        let (path, kind, status, exit_code, checks, json) = state;
        Doctor {
            path,
            kind,
            status,
            exit_code,
            checks: checks.into_iter().map(DoctorCheck::_restore).collect(),
            json,
        }
    }

    fn __reduce__<'py>(slf: &Bound<'py, Self>) -> PyResult<(Bound<'py, PyAny>, (DoctorState,))> {
        let d = slf.get();
        let checks = d
            .checks
            .iter()
            .map(|c| {
                (
                    c.name.clone(),
                    c.status.clone(),
                    c.value.clone(),
                    c.reason.clone(),
                )
            })
            .collect();
        let state = (
            d.path.clone(),
            d.kind.clone(),
            d.status.clone(),
            d.exit_code,
            checks,
            d.json.clone(),
        );
        Ok((slf.get_type().getattr("_restore")?, (state,)))
    }

    /// The check named `name`; `KeyError` when there is none.
    fn __getitem__(&self, name: &str) -> PyResult<DoctorCheck> {
        self.checks
            .iter()
            .find(|c| c.name == name)
            .cloned()
            .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(name.to_owned()))
    }

    /// The text `bammetrics doctor --format json` prints.
    fn to_json(&self) -> String {
        self.json.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "Doctor(path={:?}, kind={:?}, status={:?}, checks={})",
            self.path,
            self.kind,
            self.status,
            self.checks.len()
        )
    }
}

fn doctor_of(report: &Report, code: u8, records: u64) -> Doctor {
    Doctor {
        path: report.path.to_string_lossy().into_owned(),
        kind: report.kind.clone(),
        status: report.status().name().to_owned(),
        exit_code: code,
        checks: report
            .checks
            .iter()
            .map(|c| DoctorCheck {
                name: c.name.to_owned(),
                status: c.status.name().to_owned(),
                value: c.value.clone(),
                reason: c.reason.clone(),
            })
            .collect(),
        json: cli_doctor::json(report, code, records),
    }
}

/// bammetrics doctor: a pre-flight check of one BAM, CRAM or SAM file (its
/// header, its end-of-file marker, its first records). Returns a Doctor:
/// .checks, one DoctorCheck per check with status pass, warn, fail or skip,
/// .status, the worst, and .exit_code under fail_on. A failed check is a
/// value, never an exception; InputError when the file cannot be opened.
#[pyfunction]
#[pyo3(
    name = "doctor",
    signature = (path, **options),
    text_signature = "(path, *, reference=None, records=None, check_md5=False, \
                      expect_build=None, expect_sex=None, sample_sheet=None, sample=None, \
                      fail_on='fail', threads=None)"
)]
pub fn py_doctor<'py>(
    py: Python<'py>,
    path: &Bound<'py, PyAny>,
    options: Option<&Bound<'py, PyDict>>,
) -> PyResult<Py<Doctor>> {
    let verb = table::verb(cli_doctor::VERB).expect("doctor is a Python verb");
    let (argv, _) = argv(verb, path, options, None)?;
    let args: DoctorArgs = verb.parse(&argv).map_err(|e| verb_error(&e))?;
    let (report, code) = py
        .detach(|| cli_doctor::execute(&args))
        .map_err(|e| cli_error(cli_doctor::VERB, &e))?;
    Py::new(py, doctor_of(&report, code, args.records))
}

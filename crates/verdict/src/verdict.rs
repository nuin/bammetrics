//! Evaluation of a [`Preset`] over one sample's metrics files, and the
//! `verdict.json` contract (schema `bammetrics.verdict/1`, docs/verdict.md).

use std::path::Path;

use bammetrics_core::levels::Level;
use bammetrics_core::metric::Value;
use bammetrics_core::metrics_file::{ClassBlock, MetricsDoc, format_double};

use crate::error::{Error, Result};
use crate::json::Json;
use crate::preset::{Bound, Preset, Rule};

/// The schema name written into `verdict.json`.
pub const SCHEMA: &str = "bammetrics.verdict/1";

/// One metrics file, parsed, with the path it came from.
#[derive(Clone, Debug, PartialEq)]
pub struct Input {
    pub path: String,
    pub doc: MetricsDoc,
}

impl Input {
    pub fn read(path: &Path) -> Result<Self> {
        Ok(Input {
            path: path.display().to_string(),
            doc: MetricsDoc::read(path)?,
        })
    }
}

/// The status of one check.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Status {
    Pass,
    /// The metrics file, row, column or value is absent; never a pass.
    NotMeasured,
    Warn,
    Fail,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Pass => "pass",
            Status::Warn => "warn",
            Status::Fail => "fail",
            Status::NotMeasured => "not_measured",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        match text {
            "pass" => Some(Status::Pass),
            "warn" => Some(Status::Warn),
            "fail" => Some(Status::Fail),
            "not_measured" => Some(Status::NotMeasured),
            _ => None,
        }
    }
}

/// The verdict band over every check: fail when any check fails, else
/// incomplete when any is not measured, else warn when any warns, else pass.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Band {
    Pass,
    Warn,
    Fail,
    Incomplete,
}

impl Band {
    pub fn name(self) -> &'static str {
        match self {
            Band::Pass => "pass",
            Band::Warn => "warn",
            Band::Fail => "fail",
            Band::Incomplete => "incomplete",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        match text {
            "pass" => Some(Band::Pass),
            "warn" => Some(Band::Warn),
            "fail" => Some(Band::Fail),
            "incomplete" => Some(Band::Incomplete),
            _ => None,
        }
    }

    fn of(checks: &[Check]) -> Self {
        let has = |status: Status| checks.iter().any(|c| c.status == status);
        if has(Status::Fail) {
            Band::Fail
        } else if has(Status::NotMeasured) {
            Band::Incomplete
        } else if has(Status::Warn) {
            Band::Warn
        } else {
            Band::Pass
        }
    }
}

/// How many checks ended in each status.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub pass: usize,
    pub warn: usize,
    pub fail: usize,
    pub not_measured: usize,
}

impl Counts {
    fn of(checks: &[Check]) -> Self {
        let mut counts = Counts::default();
        for check in checks {
            match check.status {
                Status::Pass => counts.pass += 1,
                Status::Warn => counts.warn += 1,
                Status::Fail => counts.fail += 1,
                Status::NotMeasured => counts.not_measured += 1,
            }
        }
        counts
    }
}

/// One rule's outcome on this sample.
#[derive(Clone, Debug, PartialEq)]
pub struct Check {
    pub rule: String,
    pub label: String,
    pub class: String,
    /// `PCT_TARGET_BASES_20X` or `PF_Q30_BASES / PF_BASES`.
    pub metric: String,
    /// `ALL_READS`, `SAMPLE`, `LIBRARY` or `READ_GROUP`.
    pub level: String,
    /// `CATEGORY=PAIR`, when the rule selects rows.
    pub select: Option<String>,
    /// The sample, library or read group of the row judged; `None` at `ALL_READS`.
    pub unit: Option<String>,
    /// The value judged; `None` when not measured or not finite.
    pub value: Option<f64>,
    /// The bounds as written: `< 0.98`.
    pub warn: Option<String>,
    pub fail: Option<String>,
    pub status: Status,
    /// One line in plain words for a warn, fail or not-measured check.
    pub reason: Option<String>,
    pub plot: Option<String>,
}

/// The preset a verdict was made with.
#[derive(Clone, Debug, PartialEq)]
pub struct PresetInfo {
    pub name: String,
    pub version: i64,
    pub hash: String,
    pub source: String,
    pub note: String,
}

/// One metrics file the verdict read, with the classes it holds.
#[derive(Clone, Debug, PartialEq)]
pub struct FileInfo {
    pub path: String,
    pub classes: Vec<String>,
}

/// The verdict for one sample.
#[derive(Clone, Debug, PartialEq)]
pub struct Verdict {
    pub schema: String,
    pub tool_version: String,
    pub sample: String,
    pub preset: PresetInfo,
    pub status: Band,
    pub counts: Counts,
    pub files: Vec<FileInfo>,
    pub checks: Vec<Check>,
}

/// Evaluate `preset` over `inputs`, one sample's metrics files. `sample`
/// names the sample; without it the first non-empty `SAMPLE` cell in any
/// class block is used, else the first input's file name up to its first `.`.
pub fn evaluate(preset: &Preset, inputs: &[Input], sample: Option<&str>) -> Verdict {
    let checks: Vec<Check> = preset
        .rules
        .iter()
        .map(|rule| check(rule, inputs))
        .collect();
    Verdict {
        schema: SCHEMA.to_owned(),
        tool_version: env!("CARGO_PKG_VERSION").to_owned(),
        sample: sample.map_or_else(|| sample_name(inputs), str::to_owned),
        preset: PresetInfo {
            name: preset.name.clone(),
            version: preset.version,
            hash: preset.hash.clone(),
            source: preset.source.clone(),
            note: preset.note.clone(),
        },
        status: Band::of(&checks),
        counts: Counts::of(&checks),
        files: inputs
            .iter()
            .map(|input| FileInfo {
                path: input.path.clone(),
                classes: input.doc.classes.iter().map(|b| b.class.clone()).collect(),
            })
            .collect(),
        checks,
    }
}

/// [`evaluate`] over the metrics files at `paths`.
pub fn evaluate_files(preset: &Preset, paths: &[&Path], sample: Option<&str>) -> Result<Verdict> {
    let inputs = paths
        .iter()
        .map(|p| Input::read(p))
        .collect::<Result<Vec<_>>>()?;
    Ok(evaluate(preset, &inputs, sample))
}

fn sample_name(inputs: &[Input]) -> String {
    let from_cells = inputs
        .iter()
        .flat_map(|input| &input.doc.classes)
        .flat_map(|block| (0..block.rows.len()).map(move |row| block.get(row, "SAMPLE")))
        .find_map(|cell| cell.and_then(Value::as_str).map(str::to_owned));
    from_cells.unwrap_or_else(|| {
        inputs.first().map_or_else(String::new, |input| {
            let name = Path::new(&input.path)
                .file_name()
                .map_or(String::new(), |n| n.to_string_lossy().into_owned());
            name.split('.').next().unwrap_or_default().to_owned()
        })
    })
}

/// The outcome of one rule on one row.
struct Outcome {
    status: Status,
    value: Option<f64>,
    unit: Option<String>,
    reason: Option<String>,
}

impl Outcome {
    fn not_measured(reason: String) -> Self {
        Outcome {
            status: Status::NotMeasured,
            value: None,
            unit: None,
            reason: Some(reason),
        }
    }
}

fn check(rule: &Rule, inputs: &[Input]) -> Check {
    let outcome = judge(rule, inputs);
    Check {
        rule: rule.id.clone(),
        label: rule.label.clone(),
        class: rule.class.clone(),
        metric: rule.metric_text(),
        level: rule.level.name().to_owned(),
        select: rule.select_text(),
        unit: outcome.unit,
        value: outcome.value,
        warn: rule.warn.map(|b| b.to_string()),
        fail: rule.fail.map(|b| b.to_string()),
        status: outcome.status,
        reason: outcome.reason,
        plot: rule.plot.clone(),
    }
}

/// The worst outcome over the rows the rule's scope selects.
fn judge(rule: &Rule, inputs: &[Input]) -> Outcome {
    let Some(block) = inputs.iter().find_map(|input| input.doc.class(&rule.class)) else {
        return Outcome::not_measured(format!("no {} file", rule.class));
    };
    if block.column(&rule.metric).is_none() {
        return Outcome::not_measured(format!("{} has no column {}", rule.class, rule.metric));
    }
    if let Some(divisor) = &rule.divisor
        && block.column(divisor).is_none()
    {
        return Outcome::not_measured(format!("{} has no column {divisor}", rule.class));
    }
    let rows: Vec<usize> = (0..block.rows.len())
        .filter(|&row| at_level(block, row, rule.level) && selected(block, row, rule))
        .collect();
    if rows.is_empty() {
        let scope = match rule.select_text() {
            Some(select) => format!("{select} at {}", rule.level.name()),
            None => rule.level.name().to_owned(),
        };
        return Outcome::not_measured(format!("no {} row with {scope}", rule.class));
    }
    rows.into_iter()
        .map(|row| judge_row(rule, block, row))
        .max_by_key(|outcome| outcome.status)
        .expect("at least one row")
}

/// Rows at `level`: `ALL_READS` rows leave `SAMPLE`, `LIBRARY` and
/// `READ_GROUP` empty; a `SAMPLE` row fills only `SAMPLE`; a `LIBRARY` row
/// fills `LIBRARY` and not `READ_GROUP`; a `READ_GROUP` row fills
/// `READ_GROUP`. A column the class lacks counts as empty.
fn at_level(block: &ClassBlock, row: usize, level: Level) -> bool {
    let filled = |column: &str| block.get(row, column).is_some_and(|v| !v.is_null());
    match level {
        Level::AllReads => !filled("SAMPLE") && !filled("LIBRARY") && !filled("READ_GROUP"),
        Level::Sample => filled("SAMPLE") && !filled("LIBRARY") && !filled("READ_GROUP"),
        Level::Library => filled("LIBRARY") && !filled("READ_GROUP"),
        Level::ReadGroup => filled("READ_GROUP"),
    }
}

fn selected(block: &ClassBlock, row: usize, rule: &Rule) -> bool {
    match &rule.select {
        None => true,
        Some((column, text)) => block
            .get(row, column)
            .is_some_and(|v| v.picard_text() == *text),
    }
}

/// The unit a row belongs to: its `READ_GROUP`, `LIBRARY` or `SAMPLE` cell,
/// whichever the level fills.
fn unit_of(block: &ClassBlock, row: usize, level: Level) -> Option<String> {
    let column = match level {
        Level::AllReads => return None,
        Level::Sample => "SAMPLE",
        Level::Library => "LIBRARY",
        Level::ReadGroup => "READ_GROUP",
    };
    block.get(row, column).map(Value::picard_text)
}

fn judge_row(rule: &Rule, block: &ClassBlock, row: usize) -> Outcome {
    let unit = unit_of(block, row, rule.level);
    let number = |column: &str| -> std::result::Result<f64, String> {
        let cell = block.get(row, column).expect("column checked above");
        match cell {
            Value::Null => Err(format!("{column} is empty")),
            Value::Str(text) => Err(format!("{column} is not a number ({text})")),
            Value::Bool(_) => Err(format!("{column} is a yes/no field")),
            Value::Int(_) | Value::Float(_) => {
                let v = cell.as_float().expect("a number");
                if v.is_nan() {
                    Err(format!("{column} is ?"))
                } else {
                    Ok(v)
                }
            }
        }
    };
    let value = number(&rule.metric).and_then(|v| match &rule.divisor {
        None => Ok(v),
        Some(divisor) => {
            let d = number(divisor)?;
            if d == 0.0 {
                Err(format!("{divisor} is 0"))
            } else {
                Ok(v / d)
            }
        }
    });
    let value = match value {
        Ok(v) => v,
        Err(reason) => {
            return Outcome {
                unit,
                ..Outcome::not_measured(reason)
            };
        }
    };
    let (status, broken) = match (rule.fail, rule.warn) {
        (Some(fail), _) if fail.breaks(value) => (Status::Fail, Some(fail)),
        (_, Some(warn)) if warn.breaks(value) => (Status::Warn, Some(warn)),
        _ => (Status::Pass, None),
    };
    Outcome {
        status,
        value: Some(value),
        reason: broken.map(|bound| reason_line(rule, value, status, bound, unit.as_deref())),
        unit,
    }
}

/// `Target bases at 20x is 0.91, below the fail limit of 0.95` (with
/// ` in LibA` after the value when the row belongs to a unit).
fn reason_line(
    rule: &Rule,
    value: f64,
    status: Status,
    bound: Bound,
    unit: Option<&str>,
) -> String {
    let unit = unit.map(|u| format!(" in {u}")).unwrap_or_default();
    format!(
        "{} is {}{unit}, {}",
        rule.label,
        format_double(value),
        bound.describe(status.name())
    )
}

impl Check {
    /// This check as the object `verdict.json` holds under `checks`.
    pub fn to_json(&self) -> Json {
        Json::Object(vec![
            ("rule".into(), Json::text(&self.rule)),
            ("label".into(), Json::text(&self.label)),
            ("class".into(), Json::text(&self.class)),
            ("metric".into(), Json::text(&self.metric)),
            ("level".into(), Json::text(&self.level)),
            ("select".into(), Json::option_text(self.select.as_deref())),
            ("unit".into(), Json::option_text(self.unit.as_deref())),
            ("value".into(), Json::number(self.value)),
            ("warn".into(), Json::option_text(self.warn.as_deref())),
            ("fail".into(), Json::option_text(self.fail.as_deref())),
            ("status".into(), Json::text(self.status.name())),
            ("reason".into(), Json::option_text(self.reason.as_deref())),
            ("plot".into(), Json::option_text(self.plot.as_deref())),
        ])
    }
}

impl Counts {
    /// The `counts` object of `verdict.json`.
    pub fn to_json(&self) -> Json {
        let count = |n: usize| Json::Number(n as f64);
        Json::Object(vec![
            ("pass".into(), count(self.pass)),
            ("warn".into(), count(self.warn)),
            ("fail".into(), count(self.fail)),
            ("not_measured".into(), count(self.not_measured)),
        ])
    }
}

impl Verdict {
    /// The `verdict.json` text.
    pub fn to_json(&self) -> String {
        let text = |s: &String| Json::text(s);
        let files = self
            .files
            .iter()
            .map(|f| {
                Json::Object(vec![
                    ("path".into(), text(&f.path)),
                    (
                        "classes".into(),
                        Json::Array(f.classes.iter().map(text).collect()),
                    ),
                ])
            })
            .collect();
        Json::Object(vec![
            ("schema".into(), text(&self.schema)),
            ("tool_version".into(), text(&self.tool_version)),
            ("sample".into(), text(&self.sample)),
            (
                "preset".into(),
                Json::Object(vec![
                    ("name".into(), text(&self.preset.name)),
                    ("version".into(), Json::Number(self.preset.version as f64)),
                    ("hash".into(), text(&self.preset.hash)),
                    ("source".into(), text(&self.preset.source)),
                    ("note".into(), text(&self.preset.note)),
                ]),
            ),
            ("status".into(), Json::text(self.status.name())),
            ("counts".into(), self.counts.to_json()),
            ("files".into(), Json::Array(files)),
            (
                "checks".into(),
                Json::Array(self.checks.iter().map(Check::to_json).collect()),
            ),
        ])
        .render()
    }

    /// A verdict read back from `verdict.json` text. Only the schema this
    /// crate writes is accepted.
    pub fn from_json(text: &str) -> Result<Self> {
        let json = Json::parse(text)?;
        let schema = string(&json, "schema")?;
        if schema != SCHEMA {
            return Err(Error::Json(format!(
                "schema {schema:?}, expected {SCHEMA:?}"
            )));
        }
        let preset = json.get("preset").ok_or_else(|| missing("preset"))?;
        let counts = json.get("counts").ok_or_else(|| missing("counts"))?;
        let count = |key: &str| {
            counts
                .get(key)
                .and_then(Json::as_f64)
                .ok_or_else(|| missing(&format!("counts.{key}")))
                .map(|n| n as usize)
        };
        let files = array(&json, "files")?
            .iter()
            .map(|f| {
                Ok(FileInfo {
                    path: string(f, "path")?,
                    classes: array(f, "classes")?
                        .iter()
                        .map(|c| {
                            c.as_str()
                                .map(str::to_owned)
                                .ok_or_else(|| missing("files.classes"))
                        })
                        .collect::<Result<Vec<_>>>()?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let checks = array(&json, "checks")?
            .iter()
            .map(|c| {
                let status_text = string(c, "status")?;
                Ok(Check {
                    rule: string(c, "rule")?,
                    label: string(c, "label")?,
                    class: string(c, "class")?,
                    metric: string(c, "metric")?,
                    level: string(c, "level")?,
                    select: optional_string(c, "select")?,
                    unit: optional_string(c, "unit")?,
                    value: c.get("value").and_then(Json::as_f64),
                    warn: optional_string(c, "warn")?,
                    fail: optional_string(c, "fail")?,
                    status: Status::parse(&status_text).ok_or_else(|| {
                        Error::Json(format!("unknown check status {status_text:?}"))
                    })?,
                    reason: optional_string(c, "reason")?,
                    plot: optional_string(c, "plot")?,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let status_text = string(&json, "status")?;
        Ok(Verdict {
            schema,
            tool_version: string(&json, "tool_version")?,
            sample: string(&json, "sample")?,
            preset: PresetInfo {
                name: string(preset, "name")?,
                version: preset
                    .get("version")
                    .and_then(Json::as_f64)
                    .ok_or_else(|| missing("preset.version"))? as i64,
                hash: string(preset, "hash")?,
                source: string(preset, "source")?,
                note: string(preset, "note")?,
            },
            status: Band::parse(&status_text)
                .ok_or_else(|| Error::Json(format!("unknown status {status_text:?}")))?,
            counts: Counts {
                pass: count("pass")?,
                warn: count("warn")?,
                fail: count("fail")?,
                not_measured: count("not_measured")?,
            },
            files,
            checks,
        })
    }
}

fn missing(key: &str) -> Error {
    Error::Json(format!("missing or mistyped field {key}"))
}

fn string(json: &Json, key: &str) -> Result<String> {
    json.get(key)
        .and_then(Json::as_str)
        .map(str::to_owned)
        .ok_or_else(|| missing(key))
}

fn optional_string(json: &Json, key: &str) -> Result<Option<String>> {
    match json.get(key) {
        None | Some(Json::Null) => Ok(None),
        Some(Json::Str(s)) => Ok(Some(s.clone())),
        Some(_) => Err(missing(key)),
    }
}

fn array<'a>(json: &'a Json, key: &str) -> Result<&'a [Json]> {
    json.get(key)
        .and_then(Json::as_array)
        .ok_or_else(|| missing(key))
}

//! `CompareMetrics`, mirroring Picard 3.5.0 (`picard.analysis.CompareMetrics`).
//!
//! Each `INPUT` is read as htsjdk `MetricsFile.read` reads it: the class on
//! the `## METRICS CLASS` line is looked up in [`picard_classes`], every
//! column must be one of its fields, and each cell is parsed by the field's
//! Java type. A column the file lacks keeps the value a new instance holds.
//! [`MetricsDoc`](crate::metrics_file::MetricsDoc) types cells by their text
//! and tolerates short rows, so this module has its own reader: htsjdk fails
//! a short row and ignores the cells past the last column.
//!
//! The comparison, in Picard's order:
//!
//! 1. Rows equal (each field printed by `FormatUtil` is the same text, so
//!    `0.5` and `0.50` agree and a change below the sixth decimal does
//!    not count) and histograms equal: equal. Rows equal and histograms
//!    not: `Metrics Histograms differ`, equal only with
//!    `IGNORE_HISTOGRAM_DIFFERENCES`.
//! 2. Row counts differ: one line, NOT equal. Classes differ: an error.
//! 3. Column sets differ once `METRICS_NOT_REQUIRED` is taken out of both:
//!    one line naming the columns in `HashSet` order, NOT equal.
//! 4. `METRICS_TO_IGNORE` and the `METRIC_ALLOWABLE_RELATIVE_CHANGE` names
//!    must be columns of the first file.
//! 5. Rows are paired by position, or by the `KEY` columns' values. Every
//!    field of the class except the ignored ones is compared: numbers by
//!    their double value (`d2 - d1`; both NaN is no change; a relative
//!    change below the allowed fraction, or a NaN relative change, is no
//!    difference), everything else by `equals`. Each difference is one
//!    line.
//! 6. A histogram difference adds `Metric Histograms differ`, unless
//!    ignored.
//!
//! Picard logs `<class> Metric files A and B are equal` (or `NOT equal`) at
//! INFO and each difference at ERROR, then writes `OUTPUT` (the same lines
//! under a heading with both absolute paths) and `OUTPUT_TABLE` (the value
//! differences as a metrics file). The exit status is 1 when the files are
//! NOT equal; [`run_cli`] prints the log lines without the timestamp and
//! returns it.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

mod reader;

use crate::error::{Error, Result};
use crate::java::{self, StringHashSet};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::picard_classes::MetricClass;
use crate::registry::{Invocation, Opt, Tool, write_files};
use reader::{Cell, CellKey, Doc, java_split, read_doc};

pub const TOOL: Tool = Tool {
    name: "CompareMetrics",
    summary: "one row per difference Picard reports between two metrics files: KEY, METRIC, \
              VALUE1 and VALUE2 for a value, DIFFERENCE (the line Picard logs) for every one; \
              no row when the files are equal. The command line exits 1 when they are not.",
    options: &[
        Opt::path(&["INPUT", "I"], "A metrics file; give exactly two.")
            .repeated()
            .required(),
        Opt::output(&["OUTPUT", "O"]).with_help(
            "The comparison: both files, equal or NOT equal, then one line per difference.",
        ),
        Opt::path(
            &["OUTPUT_TABLE"],
            "The value differences as a metrics file (KEY, METRIC, VALUE1, VALUE2).",
        ),
        Opt::text(
            &["METRICS_TO_IGNORE", "MI"],
            "A column left out of the comparison; both files must still have it.",
        )
        .repeated(),
        Opt::text(
            &["METRICS_NOT_REQUIRED", "MNR"],
            "A column left out of the comparison that either file may lack.",
        )
        .repeated(),
        Opt::text(
            &["METRIC_ALLOWABLE_RELATIVE_CHANGE", "MARC"],
            "NAME:FRACTION: a relative change of NAME below FRACTION is no difference.",
        )
        .repeated(),
        Opt::bool(
            &["IGNORE_HISTOGRAM_DIFFERENCES", "IHD"],
            "Histogram differences do not make the files differ (default false).",
        ),
        Opt::text(
            &["KEY"],
            "A column whose values pair the rows, in place of their order.",
        )
        .repeated(),
        Opt::text(
            &["VERBOSITY"],
            "ERROR or WARNING drops the INFO summary line on standard error (default INFO).",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_NOT_DECODED,
    ],
    run: run_tool,
    collector: None,
};

/// The columns of the rows [`TOOL`] returns.
pub const COLUMNS: &[&str] = &["KEY", "METRIC", "VALUE1", "VALUE2", "DIFFERENCE"];

/// Picard's `OUTPUT_TABLE` class.
pub const TABLE_CLASS: &str = "picard.analysis.CompareMetrics$MetricComparisonDifferences";

const OPTION_MARC: &str = "METRIC_ALLOWABLE_RELATIVE_CHANGE";
/// The NPE `StringUtil.join` raises on a null key value.
const NULL_KEY: &str = "Cannot invoke \"Object.toString()\" because \"obj\" is null";

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    Ok(compare(inv)?.output)
}

/// The Picard face: parse, compare, print the log lines on standard error,
/// write the files. The status is Picard's exit code: 0 equal, 1 NOT equal.
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(args)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    let comparison = compare(&inv)?;
    comparison.print_log();
    let status = comparison.status();
    let ToolOutput { main, extra } = comparison.output;
    write_files(
        &inv,
        &ToolOutput {
            main,
            extra: Vec::new(),
        },
    )?;
    for (path, table) in extra {
        let one = ToolOutput {
            main: None,
            extra: vec![(path.clone(), table)],
        };
        write_files(&inv, &one).map_err(|_| table_write_error(&path))?;
    }
    Ok(status)
}

/// `MetricsFile.write`'s message when `OUTPUT_TABLE` cannot be written.
fn table_write_error(path: &Path) -> Error {
    let file = JavaFile::new(&path.to_string_lossy());
    Error::Metric(format!(
        "Could not write metrics to file: {}",
        file.absolute
    ))
}

/// One comparison's result.
#[derive(Clone, Debug)]
pub struct Comparison {
    pub equal: bool,
    /// `<class> Metric files A and B are equal` (or `NOT equal`).
    pub summary: String,
    /// Picard's log lines: the level (`INFO`, `ERROR`) and the message,
    /// already filtered by `VERBOSITY`.
    pub log: Vec<(&'static str, String)>,
    /// The `OUTPUT` text, whether or not `OUTPUT` was given.
    pub report: String,
    /// The rows of [`COLUMNS`] with the `OUTPUT` text (empty without
    /// `OUTPUT`, so standard output stays empty as Picard's does) as the
    /// main file, and `OUTPUT_TABLE` as a further file.
    pub output: ToolOutput,
}

impl Comparison {
    /// Picard's exit code.
    pub fn status(&self) -> u8 {
        u8::from(!self.equal)
    }

    /// Each log line as Picard prints it, without the timestamp field:
    /// `INFO\tCompareMetrics\t<message>`.
    pub fn print_log(&self) {
        for (level, message) in &self.log {
            eprintln!("{level}\tCompareMetrics\t{message}");
        }
    }
}

/// The options, checked as Picard's parser and
/// `customCommandLineValidation` check them.
#[derive(Clone, Debug, Default)]
struct Options {
    inputs: Vec<String>,
    output: Option<String>,
    ignore: Vec<String>,
    not_required: Vec<String>,
    /// `MetricToAllowableRelativeChange`: names in `HashMap` insertion
    /// order (a repeated name keeps its first place, the last value wins).
    allowed: Vec<(String, f64)>,
    ignore_histograms: bool,
    keys: Vec<String>,
    info: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let text = |v: &str| Ok(v.to_owned());
        let inputs = inv.list("INPUT", Vec::new(), text)?;
        if inputs.len() < 2 {
            return Err(Error::Usage(
                "Option 'INPUT' must be specified at least 2 times.".to_owned(),
            ));
        }
        if inputs.len() > 2 {
            return Err(Error::Usage(
                "Option 'INPUT' cannot be used more than 2 times.".to_owned(),
            ));
        }
        let info = match inv.get("VERBOSITY") {
            None | Some("INFO" | "DEBUG") => true,
            Some("ERROR" | "WARNING") => false,
            Some(other) => {
                return Err(Error::Usage(format!(
                    "VERBOSITY={other}: one of ERROR, WARNING, INFO, DEBUG"
                )));
            }
        };
        let options = Self {
            inputs,
            output: inv
                .output
                .as_ref()
                .map(|p| p.to_string_lossy().into_owned()),
            ignore: inv.list("METRICS_TO_IGNORE", Vec::new(), text)?,
            not_required: inv.list("METRICS_NOT_REQUIRED", Vec::new(), text)?,
            allowed: Vec::new(),
            ignore_histograms: inv.flag("IGNORE_HISTOGRAM_DIFFERENCES", false)?,
            keys: inv.list("KEY", Vec::new(), text)?,
            info,
        };
        if let Some(output) = &options.output {
            check_writable(output)?;
        }
        let allowed = allowed_changes(&inv.list(OPTION_MARC, Vec::new(), text)?)?;
        Ok(Self { allowed, ..options })
    }

    fn allowance(&self, name: &str) -> Option<f64> {
        self.allowed
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| *v)
    }
}

/// `METRIC_ALLOWABLE_RELATIVE_CHANGE` values as Picard validates them: one
/// error line per bad value, all reported together.
fn allowed_changes(values: &[String]) -> Result<Vec<(String, f64)>> {
    let mut allowed: Vec<(String, f64)> = Vec::new();
    let mut errors = Vec::new();
    for value in values {
        let parts = java_split(value, ':');
        let [name, number] = parts.as_slice() else {
            errors.push(format!("Invalid value for Argument '{OPTION_MARC}'"));
            continue;
        };
        match java::parse_double(number) {
            Ok(v) if v > 0.0 => match allowed.iter_mut().find(|(n, _)| n == name) {
                Some(entry) => entry.1 = v,
                None => allowed.push(((*name).to_owned(), v)),
            },
            Ok(_) => errors.push(format!(
                "Value for numeric component of Argument '{OPTION_MARC}' must be > 0.0"
            )),
            Err(_) => errors.push(format!(
                "Invalid value for numeric component of Argument '{OPTION_MARC}'"
            )),
        }
    }
    if errors.is_empty() {
        Ok(allowed)
    } else {
        Err(Error::Usage(errors.join("\n")))
    }
}

/// A path as `java.io.File` holds it: `toString` (repeated `/` collapsed,
/// a trailing `/` dropped) and `getAbsolutePath` (the working directory
/// prefixed to a relative path, nothing resolved).
#[derive(Clone, Debug)]
struct JavaFile {
    given: String,
    absolute: String,
}

impl JavaFile {
    fn new(path: &str) -> Self {
        let mut given = String::new();
        for c in path.chars() {
            if !(c == '/' && given.ends_with('/')) {
                given.push(c);
            }
        }
        if given.len() > 1 && given.ends_with('/') {
            given.pop();
        }
        let absolute = if given.starts_with('/') {
            given.clone()
        } else {
            let cwd = std::env::current_dir()
                .map(|d| d.to_string_lossy().into_owned())
                .unwrap_or_default();
            if given.is_empty() {
                cwd
            } else {
                format!("{cwd}/{given}")
            }
        };
        Self { given, absolute }
    }

    /// `Path.toUri`: `file://`, the absolute path percent-encoded, `/` after
    /// a directory.
    fn uri(&self) -> String {
        let mut uri = String::from("file://");
        for byte in self.absolute.bytes() {
            if byte.is_ascii_alphanumeric() || b"/-_.!~*'()@:&=+$,;".contains(&byte) {
                uri.push(char::from(byte));
            } else {
                uri.push_str(&format!("%{byte:02X}"));
            }
        }
        if Path::new(&self.absolute).is_dir() && !uri.ends_with('/') {
            uri.push('/');
        }
        uri
    }
}

/// `IOUtil.assertFileIsWritable` on `OUTPUT`.
pub(crate) fn check_writable(output: &str) -> Result<()> {
    let file = JavaFile::new(output);
    let path = Path::new(&file.absolute);
    if path.is_dir() {
        return Err(Error::Metric(format!(
            "Cannot write file because it is a directory: {}",
            file.absolute
        )));
    }
    let parent_exists = path.parent().is_some_and(Path::is_dir);
    if !path.exists() && !parent_exists {
        return Err(Error::Metric(format!(
            "Cannot write file: {}. Neither file nor parent directory exist.",
            file.absolute
        )));
    }
    Ok(())
}

/// `IOUtil.assertFileIsReadable` on an `INPUT`, then its text (bytes that
/// are not UTF-8 read as U+FFFD, as `FileReader` reads them).
fn read_text(file: &JavaFile) -> Result<String> {
    let path = Path::new(&file.absolute);
    if !path.exists() {
        return Err(Error::Metric(format!(
            "Cannot read non-existent file: {}",
            file.uri()
        )));
    }
    if path.is_dir() {
        return Err(Error::Metric(format!(
            "Cannot read file because it is a directory: {}",
            file.uri()
        )));
    }
    let bytes = fs::read(path)
        .map_err(|_| Error::Metric(format!("File exists but is not readable: {}", file.uri())))?;
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// `MetricBase.equals` over two lists of rows of one class each.
fn rows_equal(a: &Doc, b: &Doc) -> bool {
    if a.rows.len() != b.rows.len() {
        return false;
    }
    if a.rows.is_empty() {
        return true;
    }
    let same_class = matches!((a.class, b.class), (Some(x), Some(y)) if x.name == y.name);
    same_class
        && a.rows.iter().zip(&b.rows).all(|(x, y)| {
            x.iter().zip(y).all(|(v1, v2)| match v1 {
                Cell::Null => matches!(v2, Cell::Null),
                _ => v1.format() == v2.format(),
            })
        })
}

/// One line of Picard's `differences` list; a value difference also goes
/// to `OUTPUT_TABLE`.
#[derive(Clone, Debug)]
struct Difference {
    message: String,
    value: Option<(String, &'static str, Cell, Cell)>,
}

impl Difference {
    fn line(message: String) -> Self {
        Self {
            message,
            value: None,
        }
    }
}

/// The two inputs with what the comparison prints about them.
struct Inputs {
    files: [JavaFile; 2],
    docs: [Doc; 2],
}

/// `compareMetricsFiles`: whether the files are equal, with each line added
/// to `differences`.
fn compare_docs(
    inputs: &Inputs,
    options: &Options,
    differences: &mut Vec<Difference>,
) -> Result<bool> {
    let [d1, d2] = &inputs.docs;
    let [f1, f2] = &inputs.files;
    let histograms_equal = d1.histograms.len() == d2.histograms.len()
        && d1
            .histograms
            .iter()
            .zip(&d2.histograms)
            .all(|(a, b)| a.equals(b));
    if rows_equal(d1, d2) {
        if histograms_equal {
            return Ok(true);
        }
        if options.ignore_histograms {
            differences.push(Difference::line(
                "Metrics Histograms differ, but the 'IGNORE_HISTOGRAM_DIFFERENCES' flag is set."
                    .to_owned(),
            ));
            return Ok(true);
        }
        differences.push(Difference::line("Metrics Histograms differ".to_owned()));
        return Ok(false);
    }
    if d1.rows.len() != d2.rows.len() {
        differences.push(Difference::line(format!(
            "Number of metric rows differ between {} and {}",
            f1.absolute, f2.absolute
        )));
        return Ok(false);
    }
    // Equal row counts that are not equal rows: both files have a class.
    let (Some(class), Some(class2)) = (d1.class, d2.class) else {
        return Ok(false);
    };
    if class.name != class2.name {
        return Err(Error::Metric(format!(
            "Metrics are of differing class between {} and {}",
            f1.absolute, f2.absolute
        )));
    }
    if let Some(columns) = column_difference(d1, d2, &options.not_required) {
        differences.push(Difference::line(format!(
            "Metric columns differ between {} and {} ({columns})",
            f1.absolute, f2.absolute
        )));
        return Ok(false);
    }
    validate_names(d1, f1, options.ignore.iter().map(String::as_str))?;
    let mut allowed = StringHashSet::new();
    for (name, _) in &options.allowed {
        allowed.insert(name);
    }
    validate_names(d1, f1, allowed.iter())?;
    let ignore: Vec<&str> = options
        .ignore
        .iter()
        .chain(&options.not_required)
        .map(String::as_str)
        .collect();
    let entry = Entry {
        class,
        ignore: &ignore,
        options,
    };
    let mut differ = false;
    if options.keys.is_empty() {
        for (i, (r1, r2)) in d1.rows.iter().zip(&d2.rows).enumerate() {
            differ |= entry.compare(r1, r2, &i.to_string(), differences);
        }
    } else {
        let map1 = KeyedRows::new(d1, class, &options.keys)?;
        let mut map2 = KeyedRows::new(d2, class, &options.keys)?;
        for r1 in map1.remaining() {
            let text = map1.text(r1)?;
            match map2.take(r1) {
                Some(r2) => {
                    differ |= entry.compare(r1, r2, &text, differences);
                }
                None => {
                    differences.push(Difference::line(format!(
                        "KEY {text} found in {} but not in {}",
                        f1.given, f2.given
                    )));
                    differ = true;
                }
            }
        }
        for r2 in map2.remaining() {
            differences.push(Difference::line(format!(
                "KEY {} found in {} but not in {}",
                map2.text(r2)?,
                f2.given,
                f1.given
            )));
            differ = true;
        }
    }
    if !options.ignore_histograms && !histograms_equal {
        differences.push(Difference::line("Metric Histograms differ".to_owned()));
        differ = true;
    }
    Ok(!differ)
}

/// The columns in one file only, after both drop `not_required`, joined by
/// `,` in the order of the `HashSet` Picard builds; `None` when the sets
/// are equal.
fn column_difference(d1: &Doc, d2: &Doc, not_required: &[String]) -> Option<String> {
    let mut c1 = StringHashSet::copy_of(d1.columns.iter());
    let mut c2 = StringHashSet::copy_of(d2.columns.iter());
    for name in not_required {
        c1.remove(name);
        c2.remove(name);
    }
    if c1.same_keys(&c2) {
        return None;
    }
    let mut only = StringHashSet::copy_of(c1.iter());
    for name in c2.iter() {
        only.remove(name);
    }
    let mut only2 = StringHashSet::copy_of(c2.iter());
    for name in c1.iter() {
        only2.remove(name);
    }
    for name in only2.iter() {
        only.insert(name);
    }
    Some(only.iter().join(","))
}

/// `validateMetricNames`: every name must be a column of the first file.
fn validate_names<'a>(
    doc: &Doc,
    file: &JavaFile,
    names: impl IntoIterator<Item = &'a str>,
) -> Result<()> {
    let mut missing = StringHashSet::copy_of(names);
    for column in doc.columns.iter() {
        missing.remove(column);
    }
    if missing.is_empty() {
        return Ok(());
    }
    Err(Error::Metric(format!(
        "Metric(s) of the name: {} were not found in {}",
        missing.iter().join(", "),
        file.absolute
    )))
}

/// `buildMetricsMap`: rows by the values of the `KEY` fields, in first
/// appearance order; a repeated key keeps its place and takes the later row.
/// `index` finds a key's slot in one hash lookup; a slot taken by the
/// pairing is `None`.
struct KeyedRows<'d> {
    rows: Vec<Option<&'d [Cell]>>,
    index: HashMap<Vec<CellKey<'d>>, usize>,
    /// The `KEY` fields' positions in a row.
    slots: Vec<usize>,
}

impl<'d> KeyedRows<'d> {
    fn new(doc: &'d Doc, class: &MetricClass, keys: &[String]) -> Result<Self> {
        let slots = keys
            .iter()
            .map(|k| class.field(k).ok_or_else(|| Error::Metric(k.clone())))
            .collect::<Result<Vec<usize>>>()?;
        let mut map = Self {
            rows: Vec::new(),
            index: HashMap::new(),
            slots,
        };
        for row in &doc.rows {
            let key = map.key(row);
            match map.index.get(&key) {
                Some(&slot) => map.rows[slot] = Some(row),
                None => {
                    map.index.insert(key, map.rows.len());
                    map.rows.push(Some(row));
                }
            }
        }
        Ok(map)
    }

    fn key(&self, row: &'d [Cell]) -> Vec<CellKey<'d>> {
        self.slots.iter().map(|&s| row[s].key()).collect()
    }

    /// The row whose key equals the key of `row` (a row of the other file),
    /// taken out of the map.
    fn take(&mut self, row: &'d [Cell]) -> Option<&'d [Cell]> {
        let slot = *self.index.get(&self.key(row))?;
        self.rows[slot].take()
    }

    /// The rows not taken, in first appearance order.
    fn remaining(&self) -> impl Iterator<Item = &'d [Cell]> + '_ {
        self.rows.iter().filter_map(|r| *r)
    }

    /// `StringUtil.join(",", key)` of `row`'s key.
    fn text(&self, row: &[Cell]) -> Result<String> {
        key_text(self.slots.iter().map(|&s| &row[s]))
    }
}

/// `StringUtil.join(",", key)`: each value's `toString`, an error on null.
fn key_text<'c>(key: impl Iterator<Item = &'c Cell>) -> Result<String> {
    let parts = key
        .map(|cell| match cell {
            Cell::Null => Err(Error::Metric(NULL_KEY.to_owned())),
            other => Ok(other.java_string()),
        })
        .collect::<Result<Vec<String>>>()?;
    Ok(parts.join(","))
}

/// What `compareMetricsForEntry` needs beside the two rows.
struct Entry<'a> {
    class: &'static MetricClass,
    ignore: &'a [&'a str],
    options: &'a Options,
}

impl Entry<'_> {
    /// Every field not ignored; true when one differs.
    fn compare(&self, r1: &[Cell], r2: &[Cell], key: &str, out: &mut Vec<Difference>) -> bool {
        let mut differ = false;
        for (i, field) in self.class.fields.iter().enumerate() {
            if self.ignore.contains(&field.name) {
                continue;
            }
            let (v1, v2) = (&r1[i], &r2[i]);
            let Some(description) = self.difference(v1, v2, field.name) else {
                continue;
            };
            differ = true;
            out.push(Difference {
                message: format!(
                    "Key: {key} Metric: {} values differ. Value1: {} Value2: {} {description}",
                    field.name,
                    v1.java_string(),
                    v2.java_string()
                ),
                value: Some((key.to_owned(), field.name, v1.clone(), v2.clone())),
            });
        }
        differ
    }

    /// `compareMetricValues`: `None` when equal, else the description that
    /// ends the difference line.
    fn difference(&self, v1: &Cell, v2: &Cell, name: &str) -> Option<String> {
        if matches!(v1, Cell::Null) || matches!(v2, Cell::Null) {
            let both = matches!((v1, v2), (Cell::Null, Cell::Null));
            return (!both).then(|| "One of the values is null".to_owned());
        }
        let (Some(d1), Some(d2)) = (v1.number(), v2.number()) else {
            return (!v1.same(v2)).then(String::new);
        };
        let change = if d1.is_nan() && d2.is_nan() {
            0.0
        } else {
            d2 - d1
        };
        if change == 0.0 {
            return None;
        }
        let relative = if d1 == 0.0 { f64::MAX } else { change / d1 };
        let changed = format!(
            "Changed by {} (relative change of {})",
            java::double_to_string(change),
            java::double_to_string(relative)
        );
        match self.options.allowance(name) {
            // `dcmpl` puts a NaN relative change below the tolerance.
            Some(tolerance) if relative.is_nan() || relative.abs() < tolerance => None,
            Some(tolerance) => Some(format!(
                "{changed} which is outside of the allowable relative change tolerance of {}",
                java::double_to_string(tolerance)
            )),
            None => Some(changed),
        }
    }
}

/// One `OUTPUT_TABLE` row.
struct TableRow(Vec<Value>);

impl MetricRow for TableRow {
    const CLASS: &'static str = TABLE_CLASS;
    const COLUMNS: &'static [&'static str] = &["KEY", "METRIC", "VALUE1", "VALUE2"];

    fn values(&self) -> Vec<Value> {
        self.0.clone()
    }
}

/// Run the comparison of the two `INPUT` files.
pub fn compare(inv: &Invocation) -> Result<Comparison> {
    let options = Options::from_invocation(inv)?;
    let files = [
        JavaFile::new(&options.inputs[0]),
        JavaFile::new(&options.inputs[1]),
    ];
    let texts = [read_text(&files[0])?, read_text(&files[1])?];
    let docs = [
        read_doc(&texts[0]).map_err(Error::Metric)?,
        read_doc(&texts[1]).map_err(Error::Metric)?,
    ];
    let class = docs
        .iter()
        .find(|d| !d.rows.is_empty())
        .and_then(|d| d.class)
        .map_or("Unknown", |c| c.name);
    let inputs = Inputs { files, docs };
    let mut differences = Vec::new();
    let equal = compare_docs(&inputs, &options, &mut differences)?;
    let [f1, f2] = &inputs.files;
    let status = if equal { "equal" } else { "NOT equal" };
    let summary = format!(
        "{class} Metric files {} and {} are {status}",
        f1.given, f2.given
    );
    let mut log = Vec::new();
    if options.info {
        log.push(("INFO", summary.clone()));
    }
    log.extend(differences.iter().map(|d| ("ERROR", d.message.clone())));
    let lines: Vec<&str> = differences.iter().map(|d| d.message.as_str()).collect();
    let report = format!(
        "Comparison of {class} metrics between files {} and {}\n\nMetrics are {status}\n\n{}",
        f1.absolute,
        f2.absolute,
        lines.join("\n")
    );
    let rows = differences
        .iter()
        .map(|d| {
            let mut row = match &d.value {
                Some((key, metric, v1, v2)) => vec![
                    Value::Str(key.clone()),
                    Value::Str((*metric).to_owned()),
                    v1.value(),
                    v2.value(),
                ],
                None => vec![Value::Null; 4],
            };
            row.push(Value::Str(d.message.clone()));
            row
        })
        .collect();
    let text = if inv.output.is_some() {
        report.clone()
    } else {
        String::new()
    };
    let mut output = ToolOutput::file(MetricsOutput::text(COLUMNS, rows, text));
    if let Some(path) = inv.get("OUTPUT_TABLE") {
        let table: Vec<TableRow> = differences
            .iter()
            .filter_map(|d| d.value.as_ref())
            .map(|(key, metric, v1, v2)| {
                TableRow(vec![
                    Value::Str(key.clone()),
                    Value::Str((*metric).to_owned()),
                    v1.value(),
                    v2.value(),
                ])
            })
            .collect();
        output
            .extra
            .push((path.into(), MetricsOutput::from_rows(&table, None)));
    }
    Ok(Comparison {
        equal,
        summary,
        log,
        report,
        output,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn invocation(args: &[&str]) -> Invocation {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        TOOL.parse(&args).expect("parses")
    }

    #[test]
    fn allowable_changes_validate_like_picard() {
        let values: Vec<String> = ["A:0.5", "B", "C:0", "D:x", "E:NaN"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let Err(Error::Usage(message)) = allowed_changes(&values) else {
            panic!("usage error");
        };
        assert_eq!(
            message.lines().collect::<Vec<_>>(),
            [
                "Invalid value for Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE'",
                "Value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE' must be > 0.0",
                "Invalid value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE'",
                "Value for numeric component of Argument 'METRIC_ALLOWABLE_RELATIVE_CHANGE' must be > 0.0",
            ]
        );
        let values: Vec<String> = ["A:0.5", "A:0.25d", "B:Infinity"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        assert_eq!(
            allowed_changes(&values).expect("valid"),
            [("A".to_owned(), 0.25), ("B".to_owned(), f64::INFINITY)]
        );
    }

    #[test]
    fn value_differences_read_like_picard() {
        let options = Options {
            allowed: vec![("MEAN_READ_LENGTH".to_owned(), 0.001)],
            ..Options::default()
        };
        let class = crate::picard_classes::metric_class("picard.analysis.AlignmentSummaryMetrics")
            .expect("class");
        let entry = Entry {
            class,
            ignore: &[],
            options: &options,
        };
        let d = |a: Cell, b: Cell, name: &str| entry.difference(&a, &b, name);
        assert_eq!(
            d(Cell::Int(5), Cell::Int(6), "TOTAL_READS").as_deref(),
            Some("Changed by 1.0 (relative change of 0.2)")
        );
        assert_eq!(
            d(Cell::Double(0.0), Cell::Double(0.25), "PCT_CHIMERAS").as_deref(),
            Some("Changed by 0.25 (relative change of 1.7976931348623157E308)")
        );
        assert_eq!(
            d(Cell::Double(f64::NAN), Cell::Double(0.0), "PCT_CHIMERAS").as_deref(),
            Some("Changed by NaN (relative change of NaN)")
        );
        assert_eq!(
            d(
                Cell::Double(f64::NAN),
                Cell::Double(f64::NAN),
                "PCT_CHIMERAS"
            ),
            None
        );
        assert_eq!(
            d(Cell::Double(-0.0), Cell::Double(0.0), "PCT_CHIMERAS"),
            None
        );
        assert_eq!(
            d(
                Cell::Double(100.0),
                Cell::Double(100.04),
                "MEAN_READ_LENGTH"
            ),
            None
        );
        assert_eq!(
            d(Cell::Double(100.0), Cell::Double(101.0), "MEAN_READ_LENGTH").as_deref(),
            Some(
                "Changed by 1.0 (relative change of 0.01) which is outside of the allowable \
                 relative change tolerance of 0.001"
            )
        );
        assert_eq!(
            d(Cell::Null, Cell::Text("S1".into()), "SAMPLE").as_deref(),
            Some("One of the values is null")
        );
        assert_eq!(
            d(
                Cell::Text("PAIR".into()),
                Cell::Text("UNPAIRED".into()),
                "CATEGORY"
            )
            .as_deref(),
            Some("")
        );
    }

    #[test]
    fn missing_columns_are_named_in_hash_set_order() {
        let text =
            include_str!("../../../tests/data/picard/NA12878/NA12878.alignment_summary_metrics");
        let full = read_doc(text).expect("reads");
        let mut lines: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
        let start = lines
            .iter()
            .position(|l| l[0] == "CATEGORY")
            .expect("column line");
        let dropped: Vec<usize> = ["TOTAL_READS", "PF_READS", "MEAN_READ_LENGTH"]
            .iter()
            .map(|c| lines[start].iter().position(|x| x == c).expect("a column"))
            .collect();
        for line in &mut lines[start..start + 4] {
            *line = line
                .iter()
                .enumerate()
                .filter(|(i, _)| !dropped.contains(i))
                .map(|(_, c)| *c)
                .collect();
        }
        let shorter: Vec<String> = lines.iter().map(|l| l.join("\t")).collect();
        let short = read_doc(&shorter.join("\n")).expect("reads");
        // Picard 3.5.0 printed this order for the same files.
        assert_eq!(
            column_difference(&full, &short, &[]).as_deref(),
            Some("MEAN_READ_LENGTH,TOTAL_READS,PF_READS")
        );
        let not_required = [
            "TOTAL_READS".to_owned(),
            "PF_READS".to_owned(),
            "MEAN_READ_LENGTH".to_owned(),
        ];
        assert_eq!(column_difference(&full, &short, &not_required), None);
    }

    #[test]
    fn input_count_and_verbosity_are_checked() {
        let error = Options::from_invocation(&invocation(&["I=a.txt"])).expect_err("one input");
        assert_eq!(
            error.to_string(),
            "Option 'INPUT' must be specified at least 2 times."
        );
        let error = Options::from_invocation(&invocation(&["I=a", "I=b", "I=c"]))
            .expect_err("three inputs");
        assert_eq!(
            error.to_string(),
            "Option 'INPUT' cannot be used more than 2 times."
        );
        let options = Options::from_invocation(&invocation(&["I=a", "I=b", "VERBOSITY=ERROR"]))
            .expect("valid");
        assert!(!options.info);
        assert!(Options::from_invocation(&invocation(&["I=a", "I=b", "VERBOSITY=LOUD"])).is_err());
        assert!(
            TOOL.parse(&["I=a".to_owned(), "AS=true".to_owned()])
                .is_err()
        );
    }

    #[test]
    fn java_file_paths() {
        let file = JavaFile::new("/tmp//a/b.txt/");
        assert_eq!(file.given, "/tmp/a/b.txt");
        assert_eq!(file.absolute, "/tmp/a/b.txt");
        let relative = JavaFile::new(".//x.txt");
        assert_eq!(relative.given, "./x.txt");
        assert!(relative.absolute.ends_with("/./x.txt"));
        assert_eq!(JavaFile::new("/a b").uri(), "file:///a%20b");
    }
}

//! `IntervalListTools`, mirroring Picard 3.5.0 (`picard.util.IntervalListTools`)
//! over the interval layer [`crate::interval_ops`].
//!
//! The steps of `doWork`, in order: each `INPUT` read and padded by
//! `PADDING`, the lists reduced by the action (`INTERSECT` intersects them,
//! every other action concatenates them), the same for `SECOND_INPUT`, then
//! the action. `INVERT` turns `SORT` off and `UNIQUE` on, `UNION` turns both
//! on. The result is sorted, inverted, merged (`DONT_MERGE_ABUTTING` keeps
//! abutting intervals apart) and broken at `BREAK_BANDS_AT_MULTIPLES_OF`,
//! in that order. The output header is the action's, with one `@PG` line
//! (Picard's command line, [`picard_command_line`]) and one `@CO` line per
//! `COMMENT`. `SCATTER_CONTENT` sets `SCATTER_COUNT`; with `OUTPUT` and a
//! count above 1 the list is scattered into
//! `OUTPUT/temp_NNNN_of_N/scattered.interval_list`, else written to
//! `OUTPUT`. `OUTPUT_VALUE` prints the base or interval count of what was
//! written (of the whole list without `OUTPUT`) to standard output or
//! `COUNT_OUTPUT`.
//!
//! An input is an interval list unless its file name is exactly `.vcf`,
//! `.vcf.gz` or `.vcf.bgz`: Picard compares the name with
//! `Path.endsWith(String)`, which matches whole names (D95), and reads such
//! a file as a VCF with `INCLUDE_FILTERED`. A file named `.bcf` is refused.
//!
//! Every output is checked before the first is written: no output is an
//! input, the directory of `OUTPUT` and `COUNT_OUTPUT` exists (a scatter
//! makes its directories, as Picard does), and `COUNT_OUTPUT` is not
//! `OUTPUT`. A scatter is written as one set ([`interval_ops::write_scatter`]).

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::interval_ops::{
    Header, Interval, IntervalList, ScatterMode, Scattered, break_at_band_multiples, concatenate,
    difference, intersection, invert, overlaps, scatter, scatter_count_for_content, subtract,
    union, unique_intervals, write_list, write_scatter,
};
use crate::metric::{MetricsOutput, RowSource, ToolOutput, Value};
use crate::registry::{self, Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "IntervalListTools",
    summary: "one row per interval written (list, contig, start, end, strand, name; list counts \
              the scatter lists from 1); OUTPUT gets the interval list, or the scatter directory \
              when SCATTER_COUNT is above 1.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "An interval list; a file named .vcf, .vcf.gz or .vcf.bgz is read as a VCF, as \
             Picard reads it. Each adds a list.",
        )
        .required()
        .repeated(),
        Opt::path(
            &["OUTPUT", "O"],
            "The interval list to write, or the directory of the scatter lists when \
             SCATTER_COUNT is above 1; without it no list is written.",
        )
        .written(),
        Opt::int(
            &["PADDING"],
            "Bases added to both ends of every input interval before the action (default 0; \
             negative is an error, as in htsjdk 5.0.0).",
        ),
        Opt::bool(
            &["UNIQUE"],
            "Merge overlapping and abutting intervals; implies SORT (default false).",
        ),
        Opt::bool(
            &["DONT_MERGE_ABUTTING"],
            "Keep abutting intervals apart when merging (default false).",
        ),
        Opt::bool(&["SORT"], "Sort the result by coordinate (default true)."),
        Opt::text(
            &["ACTION"],
            "CONCAT (default), UNION, INTERSECT, SUBTRACT, SYMDIFF or OVERLAPS.",
        ),
        Opt::path(
            &["SECOND_INPUT", "SI"],
            "The second lists of SUBTRACT, SYMDIFF and OVERLAPS. Each adds a list.",
        )
        .repeated(),
        Opt::text(
            &["COMMENT"],
            "A @CO line for the output header. Each adds a line.",
        )
        .repeated(),
        Opt::int(
            &["SCATTER_COUNT"],
            "The lists to scatter the result into under OUTPUT (default 1; a scatter that \
             makes more than 100,000 lists is refused).",
        ),
        Opt::int(
            &["SCATTER_CONTENT"],
            "Bases (intervals in the two INTERVAL_COUNT modes) per scatter list; sets \
             SCATTER_COUNT.",
        ),
        Opt::bool(
            &["INCLUDE_FILTERED"],
            "Read the filtered records of a VCF input too (default false).",
        ),
        Opt::int(
            &["BREAK_BANDS_AT_MULTIPLES_OF", "BRK"],
            "Break intervals at multiples of N; 0 breaks none (default 0).",
        ),
        Opt::text(
            &["SUBDIVISION_MODE", "M"],
            "INTERVAL_SUBDIVISION (default), BALANCING_WITHOUT_INTERVAL_SUBDIVISION, \
             BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW, INTERVAL_COUNT or \
             INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER.",
        ),
        Opt::bool(
            &["INVERT"],
            "Write the regions of the dictionary no interval covers (default false).",
        ),
        Opt::text(
            &["OUTPUT_VALUE"],
            "NONE (default), BASES or INTERVALS: the count printed to standard output or \
             COUNT_OUTPUT.",
        ),
        Opt::path(
            &["COUNT_OUTPUT"],
            "The file the OUTPUT_VALUE count goes to; standard output without one.",
        )
        .written(),
        Opt::ignored_int(
            &["COMPRESSION_LEVEL"],
            "a .gz list is written as BGZF at the default level",
        ),
        Opt::ignored_bool(&["USE_JDK_DEFLATER"], "no JDK codec"),
        Opt::ignored_bool(&["USE_JDK_INFLATER"], "no JDK codec"),
        registry::REFERENCE_NOT_DECODED,
        registry::STRINGENCY_NO_RECORDS,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The columns of the rows Python gets.
pub const COLUMNS: [&str; 6] = ["list", "contig", "start", "end", "strand", "name"];

/// `IntervalListTools.Action`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Action {
    Concat,
    Union,
    Intersect,
    Subtract,
    Symdiff,
    Overlaps,
}

impl Action {
    const ALL: [(Action, &'static str); 6] = [
        (Action::Concat, "CONCAT"),
        (Action::Union, "UNION"),
        (Action::Intersect, "INTERSECT"),
        (Action::Subtract, "SUBTRACT"),
        (Action::Symdiff, "SYMDIFF"),
        (Action::Overlaps, "OVERLAPS"),
    ];

    fn parse(text: &str) -> Option<Action> {
        Self::ALL.iter().find(|(_, n)| *n == text).map(|(a, _)| *a)
    }

    fn name(self) -> &'static str {
        Self::ALL
            .iter()
            .find(|(a, _)| *a == self)
            .map_or("", |(_, n)| n)
    }

    fn takes_second_input(self) -> bool {
        matches!(self, Action::Subtract | Action::Symdiff | Action::Overlaps)
    }

    /// `Action.reduceEach`: how the lists of one input option combine.
    fn reduce_each(self, a: IntervalList, b: IntervalList) -> Result<IntervalList> {
        match self {
            Action::Intersect => intersection(&a, &b),
            _ => concatenate(vec![a, b]),
        }
    }

    /// `Action.act`.
    fn act(self, a: IntervalList, b: Option<IntervalList>) -> Result<IntervalList> {
        let second = || b.ok_or_else(|| Error::Usage("SECOND_INPUT is required".to_owned()));
        match self {
            Action::Concat | Action::Union | Action::Intersect => Ok(a),
            Action::Subtract => subtract(&a, &second()?),
            Action::Symdiff => difference(&union(vec![a])?, &union(vec![second()?])?),
            Action::Overlaps => overlaps(&a, &second()?),
        }
    }
}

/// `IntervalListTools.Output`: what is printed after the list is written.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Count {
    None,
    Bases,
    Intervals,
}

struct Options {
    inputs: Vec<PathBuf>,
    second_inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    action: Action,
    padding: i32,
    unique: bool,
    dont_merge_abutting: bool,
    sort: bool,
    comments: Vec<String>,
    scatter_count: i32,
    scatter_content: Option<i32>,
    include_filtered: bool,
    break_bands: i32,
    mode: ScatterMode,
    invert: bool,
    count: Count,
    count_output: Option<PathBuf>,
}

fn usage(message: String) -> Error {
    Error::Usage(message)
}

impl Options {
    /// The options and `customCommandLineValidation`, every problem in one
    /// error as Picard lists them.
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let paths = |name: &str| inv.list(name, Vec::new(), |v| Ok(java_path(v)));
        let path = |name: &str| inv.get(name).filter(|v| *v != "null").map(java_path);
        let action = match inv.get("ACTION") {
            None => Action::Concat,
            Some(text) => Action::parse(text).ok_or_else(|| {
                usage(format!(
                    "ACTION={text}: expected CONCAT, UNION, INTERSECT, SUBTRACT, SYMDIFF or \
                     OVERLAPS"
                ))
            })?,
        };
        let mode = match inv.get("SUBDIVISION_MODE") {
            None => ScatterMode::IntervalSubdivision,
            Some(text) => ScatterMode::parse(text).ok_or_else(|| {
                let names: Vec<&str> = ScatterMode::ALL.iter().map(|m| m.name()).collect();
                usage(format!(
                    "SUBDIVISION_MODE={text}: expected {}",
                    names.join(", ")
                ))
            })?,
        };
        let count = match inv.get("OUTPUT_VALUE") {
            None | Some("NONE") => Count::None,
            Some("BASES") => Count::Bases,
            Some("INTERVALS") => Count::Intervals,
            Some(text) => {
                return Err(usage(format!(
                    "OUTPUT_VALUE={text}: expected NONE, BASES or INTERVALS"
                )));
            }
        };
        let opts = Self {
            inputs: paths("INPUT")?,
            second_inputs: paths("SECOND_INPUT")?,
            output: path("OUTPUT"),
            action,
            padding: inv.parse("PADDING", 0)?,
            unique: inv.flag("UNIQUE", false)?,
            dont_merge_abutting: inv.flag("DONT_MERGE_ABUTTING", false)?,
            sort: inv.flag("SORT", true)?,
            comments: inv.list("COMMENT", Vec::new(), |v| Ok(v.to_owned()))?,
            scatter_count: inv.parse("SCATTER_COUNT", 1)?,
            scatter_content: match inv.get("SCATTER_CONTENT") {
                None | Some("null") => None,
                Some(_) => inv.parse_opt("SCATTER_CONTENT")?,
            },
            include_filtered: inv.flag("INCLUDE_FILTERED", false)?,
            break_bands: inv.parse("BREAK_BANDS_AT_MULTIPLES_OF", 0)?,
            mode,
            invert: inv.flag("INVERT", false)?,
            count,
            count_output: path("COUNT_OUTPUT"),
        };
        let mut problems = Vec::new();
        if opts.inputs.is_empty() {
            problems.push("INPUT is required".to_owned());
        }
        if opts.scatter_count < 1 {
            problems.push("SCATTER_COUNT must be greater than 0.".to_owned());
        }
        if opts.break_bands < 0 {
            problems
                .push("BREAK_BANDS_AT_MULTIPLES_OF must be greater than or equal to 0.".to_owned());
        }
        let name = opts.action.name();
        match (
            opts.second_inputs.is_empty(),
            opts.action.takes_second_input(),
        ) {
            (true, true) => problems.push(format!(
                "SECOND_INPUT was not provided but action {name} requires a second input."
            )),
            (false, false) => problems.push(format!(
                "SECOND_INPUT was provided but action {name} doesn't take a second input."
            )),
            _ => {}
        }
        if opts.count_output.is_some() && opts.count == Count::None {
            problems.push("COUNT_OUTPUT was provided but OUTPUT_VALUE is set to NONE.".to_owned());
        }
        if !problems.is_empty() {
            return Err(usage(problems.join("\n")));
        }
        Ok(opts)
    }
}

/// A path as Java's `Path` holds it: trailing `/` dropped, so `O=sc/`
/// names the file `sc`.
fn java_path(text: &str) -> PathBuf {
    match text.trim_end_matches('/') {
        "" if text.starts_with('/') => PathBuf::from("/"),
        trimmed => PathBuf::from(trimmed),
    }
}

/// True for a file name Picard's `IntervalListInputType.forFile` reads as a
/// VCF: `Path.endsWith(String)` on the file name matches the whole name.
fn vcf_name(path: &Path) -> Option<&'static str> {
    let name = path.file_name()?;
    [".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"]
        .into_iter()
        .find(|e| name == *e)
}

/// One input as `IntervalListInputType.getIntervalList` reads it.
fn read_input(path: &Path, include_filtered: bool) -> Result<IntervalList> {
    match vcf_name(path) {
        Some(".bcf") => Err(usage(format!(
            "{}: Picard reads a file named .bcf as BCF; bammetrics reads no BCF",
            path.display()
        ))),
        Some(_) => IntervalList::from_vcf(path, include_filtered),
        None => IntervalList::read(path),
    }
}

/// `openIntervalLists`: each file read and padded, then reduced by the
/// action; `None` for no file.
fn open(paths: &[PathBuf], opts: &Options) -> Result<Option<IntervalList>> {
    let mut combined: Option<IntervalList> = None;
    for path in paths {
        let list = read_input(path, opts.include_filtered)?.padded(opts.padding, opts.padding)?;
        combined = Some(match combined {
            None => list,
            Some(acc) => opts.action.reduce_each(acc, list)?,
        });
    }
    Ok(combined)
}

/// The output list of `doWork` before it is written.
fn compute(inv: &Invocation, opts: &Options) -> Result<IntervalList> {
    let first = open(&opts.inputs, opts)?.ok_or_else(|| usage("INPUT is required".to_owned()))?;
    let second = open(&opts.second_inputs, opts)?;
    let result = opts.action.act(first, second)?;
    let (mut sort, mut unique) = (opts.sort, opts.unique);
    if opts.invert {
        sort = false;
        unique = true;
    }
    if opts.action == Action::Union {
        sort = true;
        unique = true;
    }
    // The output header is the action's, before the sort sets its order.
    let mut header = result.header.clone();
    let possibly_sorted = if sort { result.into_sorted() } else { result };
    let possibly_inverted = if opts.invert {
        invert(&possibly_sorted)
    } else {
        possibly_sorted
    };
    let mut intervals = if unique {
        unique_intervals(&possibly_inverted, !opts.dont_merge_abutting, true)
    } else {
        possibly_inverted.intervals
    };
    if opts.break_bands > 0 {
        intervals = break_at_band_multiples(&intervals, opts.break_bands)?;
    }
    header.add_program(&picard_command_line(inv), TOOL.name);
    for comment in &opts.comments {
        header.add_comment(comment);
    }
    Ok(IntervalList { header, intervals })
}

/// Refuse an output that cannot be written before any is: a missing
/// directory above `OUTPUT` (a scatter makes its own) or `COUNT_OUTPUT`, an
/// input, or `COUNT_OUTPUT` naming `OUTPUT`.
fn check_outputs(inv: &Invocation, opts: &Options, scattered: bool) -> Result<()> {
    let output = opts.output.as_deref();
    let count = opts.count_output.as_deref();
    let single = output.filter(|_| !scattered);
    registry::check_writable(single.into_iter().chain(count), &inv.inputs)?;
    if let (Some(list), Some(count)) = (output, count)
        && crate::mark_duplicates::same_file(list, count)
    {
        return Err(usage(format!(
            "COUNT_OUTPUT {} is OUTPUT too; one file cannot hold both",
            count.display()
        )));
    }
    Ok(())
}

/// Refuse a `COUNT_OUTPUT` that names a file or directory of the scatter
/// under `dir`: the count would replace a list after the set is written.
fn check_count_beside_scatter(dir: &Path, count: &Path, lists: usize, of: i32) -> Result<()> {
    use crate::mark_duplicates::same_file_key;
    let key = same_file_key(count);
    let root = same_file_key(dir);
    let Some(first) = key.strip_prefix(&root).ok().and_then(|r| r.iter().next()) else {
        return Ok(());
    };
    let first = first.to_string_lossy();
    let index = first
        .strip_prefix("temp_")
        .and_then(|rest| rest.strip_suffix(&format!("_of_{of}")))
        .and_then(|digits| digits.parse::<usize>().ok().map(|n| (n, digits)));
    match index {
        Some((n, digits)) if (1..=lists).contains(&n) && format!("{n:04}") == digits => {
            Err(usage(format!(
                "COUNT_OUTPUT {} is inside list {n} of the scatter under OUTPUT {}; the count \
                 would replace it",
                count.display(),
                dir.display()
            )))
        }
        _ => Ok(()),
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    // `CommandLineProgram.instanceMain` puts the command line in a
    // `StringHeader`, which refuses a line feed; Picard fails before it
    // reads anything.
    if picard_command_line(inv).contains('\n') {
        return Err(usage(
            "an option value holds a line feed, which no header line can hold (Picard fails in \
             StringHeader)"
                .to_owned(),
        ));
    }
    let opts = Options::from_invocation(inv)?;
    let output = compute(inv, &opts)?;
    let scatter_count = match opts.scatter_content {
        Some(content) => scatter_count_for_content(opts.mode, &output, content),
        None => opts.scatter_count,
    };
    let written = match &opts.output {
        Some(dir) if scatter_count > 1 => {
            let lists = scatter(&output, opts.mode, scatter_count)?;
            check_outputs(inv, &opts, true)?;
            if let Some(count) = &opts.count_output {
                check_count_beside_scatter(dir, count, lists.len(), scatter_count)?;
            }
            write_scatter(dir, &lists, scatter_count, &inv.inputs)?;
            lists
        }
        Some(path) => {
            check_outputs(inv, &opts, false)?;
            write_list(&output, path, &inv.inputs)?;
            whole(output)
        }
        None => {
            check_outputs(inv, &opts, false)?;
            whole(output)
        }
    };
    let mut text = match opts.count {
        Count::None => String::new(),
        Count::Bases => format!("{}\n", written.base_count()),
        Count::Intervals => format!("{}\n", written.interval_count()),
    };
    // `COUNT_OUTPUT` gets the count; without it the registry prints it.
    if let Some(path) = &opts.count_output {
        registry::write_to(Some(path), |out| out.write_all(text.as_bytes()))?;
        text.clear();
    }
    let mut file = MetricsOutput::text(&COLUMNS, Vec::new(), text);
    file.source = Some(Arc::new(Rows(written)));
    Ok(file.into())
}

/// A list as a scatter of one.
fn whole(list: IntervalList) -> Scattered {
    Scattered {
        header: list.header,
        parts: vec![list.intervals],
    }
}

/// The rows Python gets: every interval written, with its list number.
#[derive(Debug)]
struct Rows(Scattered);

impl RowSource for Rows {
    fn len(&self) -> usize {
        self.0.interval_count() as usize
    }

    fn each(&self, row: &mut dyn FnMut(&[Value]) -> std::io::Result<()>) -> std::io::Result<()> {
        let header: &Header = &self.0.header;
        for (number, part) in self.0.parts.iter().enumerate() {
            for interval in part {
                row(&values(header, number + 1, interval))?;
            }
        }
        Ok(())
    }
}

fn values(header: &Header, list: usize, interval: &Interval) -> [Value; 6] {
    [
        Value::Int(list as i64),
        Value::Str(header.sequences[interval.contig].name.clone()),
        Value::Int(i64::from(interval.start)),
        Value::Int(i64::from(interval.end)),
        Value::Str(if interval.negative_strand { "-" } else { "+" }.to_owned()),
        interval.name.clone().map_or(Value::Null, Value::Str),
    ]
}

/// How Picard prints an option's value in its command line.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Shape {
    /// A collection: `[a, b]` in the legacy form, one `--KEY v` per value
    /// in the other.
    List,
    /// `true` or `false` (`Boolean.parseBoolean`).
    Bool,
    /// The parsed number.
    Int,
    /// The value as given.
    Text,
    /// Barclay's own `help`, `version` and `showHidden`: printed by the
    /// `--KEY VALUE` parser only.
    Barclay,
}

/// The options of `IntervalListTools` in Picard's declaration order (its
/// own fields, then `CommandLineProgram`'s), with the default printed for
/// an option not given (`None`: printed only when given).
const FIELDS: &[(&str, Shape, Option<&str>)] = &[
    ("INPUT", Shape::List, None),
    ("OUTPUT", Shape::Text, None),
    ("PADDING", Shape::Int, Some("0")),
    ("UNIQUE", Shape::Bool, Some("false")),
    ("DONT_MERGE_ABUTTING", Shape::Bool, Some("false")),
    ("SORT", Shape::Bool, Some("true")),
    ("ACTION", Shape::Text, Some("CONCAT")),
    ("SECOND_INPUT", Shape::List, None),
    ("COMMENT", Shape::List, None),
    ("SCATTER_COUNT", Shape::Int, Some("1")),
    ("SCATTER_CONTENT", Shape::Int, None),
    ("INCLUDE_FILTERED", Shape::Bool, Some("false")),
    ("BREAK_BANDS_AT_MULTIPLES_OF", Shape::Int, Some("0")),
    (
        "SUBDIVISION_MODE",
        Shape::Text,
        Some("INTERVAL_SUBDIVISION"),
    ),
    ("INVERT", Shape::Bool, Some("false")),
    ("OUTPUT_VALUE", Shape::Text, Some("NONE")),
    ("COUNT_OUTPUT", Shape::Text, None),
    ("TMP_DIR", Shape::List, None),
    ("VERBOSITY", Shape::Text, Some("INFO")),
    ("QUIET", Shape::Bool, Some("false")),
    ("VALIDATION_STRINGENCY", Shape::Text, Some("STRICT")),
    ("COMPRESSION_LEVEL", Shape::Int, Some("5")),
    ("MAX_RECORDS_IN_RAM", Shape::Int, Some("500000")),
    ("CREATE_INDEX", Shape::Bool, Some("false")),
    ("CREATE_MD5_FILE", Shape::Bool, Some("false")),
    ("REFERENCE_SEQUENCE", Shape::Text, None),
    ("help", Shape::Barclay, Some("false")),
    ("version", Shape::Barclay, Some("false")),
    ("showHidden", Shape::Barclay, Some("false")),
    ("USE_JDK_DEFLATER", Shape::Bool, Some("false")),
    ("USE_JDK_INFLATER", Shape::Bool, Some("false")),
];

/// One value as Picard prints it once parsed.
fn shown(shape: Shape, value: &str) -> String {
    match shape {
        Shape::Bool => value.eq_ignore_ascii_case("true").to_string(),
        Shape::Int => value
            .parse::<i64>()
            .map_or_else(|_| value.to_owned(), |n| n.to_string()),
        Shape::List | Shape::Text | Shape::Barclay => value.to_owned(),
    }
}

/// `CommandLineProgram.getCommandLine`, the `CL` of the `@PG` line: the
/// tool name, every option given in declaration order, then every other
/// option with a default. The legacy parser prints `KEY=VALUE` (a
/// collection as `[a, b]`, a null as `null`) and four spaces before the
/// defaults; the `--KEY VALUE` parser prints `--KEY VALUE` (one per
/// collection value, a null as nothing) and its own three options.
pub fn picard_command_line(inv: &Invocation) -> String {
    let legacy = inv.legacy_syntax;
    let mut given = Vec::new();
    let mut defaults = Vec::new();
    for &(name, shape, default) in FIELDS {
        let values: Vec<&str> = inv
            .given
            .iter()
            .filter(|(n, _)| *n == name)
            .map(|(_, v)| v.as_str())
            .collect();
        if values.is_empty() {
            match default {
                Some(_) if shape == Shape::Barclay && legacy => {}
                Some(value) if legacy => defaults.push(format!("{name}={value}")),
                Some(value) => defaults.push(format!("--{name} {value}")),
                None => {}
            }
            continue;
        }
        if shape == Shape::List {
            let items = registry::picard_list(Vec::new(), &values, |v| Ok(shown(shape, v)))
                .unwrap_or_default();
            if legacy {
                given.push(format!("{name}=[{}]", items.join(", ")));
            } else if items.is_empty() {
                // A collection a `null` cleared prints as nothing.
                given.push(String::new());
            } else {
                given.extend(items.iter().map(|v| format!("--{name} {v}")));
            }
            continue;
        }
        let value = values[values.len() - 1];
        given.push(match (value, legacy) {
            ("null", true) => format!("{name}=null"),
            ("null", false) => String::new(),
            (_, true) => format!("{name}={}", shown(shape, value)),
            (_, false) => format!("--{name} {}", shown(shape, value)),
        });
    }
    let mut line = TOOL.name.to_owned();
    for piece in &given {
        line.push(' ');
        line.push_str(piece);
    }
    if legacy {
        line.push_str("   ");
    }
    for piece in &defaults {
        line.push(' ');
        line.push_str(piece);
    }
    line
}

#[cfg(test)]
mod tests;

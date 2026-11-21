//! `IntervalListTools` of Picard 3.5.0 over `bammetrics_core::interval_ops`,
//! for the concordance module `tests/concordance/test_interval_ops.py`:
//!
//!     cargo run -p bammetrics-core --example interval_ops -- \
//!         ACTION=SUBTRACT I=a.interval_list SI=b.interval_list O=out.interval_list
//!
//! Arguments are Picard's `KEY=VALUE` (long or short names): `INPUT`/`I`,
//! `SECOND_INPUT`/`SI`, `OUTPUT`/`O`, `ACTION`, `PADDING`, `UNIQUE`, `SORT`,
//! `DONT_MERGE_ABUTTING`, `INVERT`, `BREAK_BANDS_AT_MULTIPLES_OF`/`BRK`,
//! `SCATTER_COUNT`, `SCATTER_CONTENT`, `SUBDIVISION_MODE`/`M`, `COMMENT`,
//! `INCLUDE_FILTERED`, `OUTPUT_VALUE`, `COUNT_OUTPUT`. The `@PG` line it
//! adds holds its own command line. Exit 0, or 1 with the error.

use std::path::PathBuf;

use bammetrics_core::error::{Error, Result};
use bammetrics_core::interval_ops::{
    IntervalList, ScatterMode, break_at_band_multiples, concatenate, difference, intersection,
    invert, overlaps, scatter, scatter_count_for_content, subtract, union, unique_intervals,
    write_list, write_scatter,
};

#[derive(Clone, Copy, PartialEq, Eq)]
enum Action {
    Concat,
    Union,
    Intersect,
    Subtract,
    Symdiff,
    Overlaps,
}

impl Action {
    fn parse(text: &str) -> Option<Action> {
        Some(match text {
            "CONCAT" => Action::Concat,
            "UNION" => Action::Union,
            "INTERSECT" => Action::Intersect,
            "SUBTRACT" => Action::Subtract,
            "SYMDIFF" => Action::Symdiff,
            "OVERLAPS" => Action::Overlaps,
            _ => return None,
        })
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
        let second = || {
            b.clone()
                .ok_or_else(|| Error::Usage("SECOND_INPUT is required".to_owned()))
        };
        match self {
            Action::Concat | Action::Union | Action::Intersect => Ok(a),
            Action::Subtract => subtract(&a, &second()?),
            Action::Symdiff => difference(&union(vec![a])?, &union(vec![second()?])?),
            Action::Overlaps => overlaps(&a, &second()?),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Action::Concat => "CONCAT",
            Action::Union => "UNION",
            Action::Intersect => "INTERSECT",
            Action::Subtract => "SUBTRACT",
            Action::Symdiff => "SYMDIFF",
            Action::Overlaps => "OVERLAPS",
        }
    }
}

struct Options {
    input: Vec<PathBuf>,
    second_input: Vec<PathBuf>,
    output: Option<PathBuf>,
    action: Action,
    padding: i32,
    unique: bool,
    sort: bool,
    dont_merge_abutting: bool,
    invert: bool,
    break_bands: i32,
    scatter_count: i32,
    scatter_content: Option<i32>,
    mode: ScatterMode,
    comments: Vec<String>,
    output_value: String,
    count_output: Option<PathBuf>,
}

fn usage(message: String) -> Error {
    Error::Usage(message)
}

fn int(key: &str, value: &str) -> Result<i32> {
    value
        .parse()
        .map_err(|_| usage(format!("{key}: not an integer: {value}")))
}

fn boolean(key: &str, value: &str) -> Result<bool> {
    match value.to_ascii_lowercase().as_str() {
        "true" | "t" => Ok(true),
        "false" | "f" => Ok(false),
        _ => Err(usage(format!("{key}: not a boolean: {value}"))),
    }
}

fn parse(args: &[String]) -> Result<Options> {
    let mut opts = Options {
        input: Vec::new(),
        second_input: Vec::new(),
        output: None,
        action: Action::Concat,
        padding: 0,
        unique: false,
        sort: true,
        dont_merge_abutting: false,
        invert: false,
        break_bands: 0,
        scatter_count: 1,
        scatter_content: None,
        mode: ScatterMode::IntervalSubdivision,
        comments: Vec::new(),
        output_value: "NONE".to_owned(),
        count_output: None,
    };
    for arg in args {
        let (key, value) = arg
            .split_once('=')
            .ok_or_else(|| usage(format!("not KEY=VALUE: {arg}")))?;
        match key {
            "INPUT" | "I" => opts.input.push(value.into()),
            "SECOND_INPUT" | "SI" => opts.second_input.push(value.into()),
            "OUTPUT" | "O" => opts.output = Some(value.into()),
            "ACTION" => {
                opts.action =
                    Action::parse(value).ok_or_else(|| usage(format!("ACTION: {value}")))?;
            }
            "PADDING" => opts.padding = int(key, value)?,
            "UNIQUE" => opts.unique = boolean(key, value)?,
            "SORT" => opts.sort = boolean(key, value)?,
            "DONT_MERGE_ABUTTING" => opts.dont_merge_abutting = boolean(key, value)?,
            "INVERT" => opts.invert = boolean(key, value)?,
            "INCLUDE_FILTERED" => {
                boolean(key, value)?;
            }
            "BREAK_BANDS_AT_MULTIPLES_OF" | "BRK" => opts.break_bands = int(key, value)?,
            "SCATTER_COUNT" => opts.scatter_count = int(key, value)?,
            "SCATTER_CONTENT" => opts.scatter_content = Some(int(key, value)?),
            "SUBDIVISION_MODE" | "M" => {
                opts.mode = ScatterMode::parse(value)
                    .ok_or_else(|| usage(format!("SUBDIVISION_MODE: {value}")))?;
            }
            "COMMENT" => opts.comments.push(value.to_owned()),
            "OUTPUT_VALUE" => {
                if !["NONE", "BASES", "INTERVALS"].contains(&value) {
                    return Err(usage(format!("OUTPUT_VALUE: {value}")));
                }
                opts.output_value = value.to_owned();
            }
            "COUNT_OUTPUT" => opts.count_output = Some(value.into()),
            _ => return Err(usage(format!("unknown option {key}"))),
        }
    }
    // `customCommandLineValidation`.
    let mut problems = Vec::new();
    if opts.input.is_empty() {
        problems.push("INPUT is required".to_owned());
    }
    if opts.scatter_count < 1 {
        problems.push("SCATTER_COUNT must be greater than 0.".to_owned());
    }
    if opts.break_bands < 0 {
        problems.push("BREAK_BANDS_AT_MULTIPLES_OF must be greater than or equal to 0.".to_owned());
    }
    if opts.second_input.is_empty() && opts.action.takes_second_input() {
        problems.push(format!(
            "SECOND_INPUT was not provided but action {} requires a second input.",
            opts.action.name()
        ));
    }
    if !opts.second_input.is_empty() && !opts.action.takes_second_input() {
        problems.push(format!(
            "SECOND_INPUT was provided but action {} doesn't take a second input.",
            opts.action.name()
        ));
    }
    if opts.count_output.is_some() && opts.output_value == "NONE" {
        problems.push("COUNT_OUTPUT was provided but OUTPUT_VALUE is set to NONE.".to_owned());
    }
    if !problems.is_empty() {
        return Err(usage(problems.join("\n")));
    }
    Ok(opts)
}

/// `openIntervalLists`: each file read and padded, then reduced.
fn open(paths: &[PathBuf], opts: &Options) -> Result<Option<IntervalList>> {
    let mut combined: Option<IntervalList> = None;
    for path in paths {
        let list = IntervalList::read(path)?.padded(opts.padding, opts.padding)?;
        combined = Some(match combined {
            None => list,
            Some(acc) => opts.action.reduce_each(acc, list)?,
        });
    }
    Ok(combined)
}

fn run(args: &[String]) -> Result<()> {
    let mut opts = parse(args)?;
    let inputs: Vec<PathBuf> = opts
        .input
        .iter()
        .chain(&opts.second_input)
        .cloned()
        .collect();
    let first = open(&opts.input, &opts)?.expect("INPUT is required");
    let second = open(&opts.second_input, &opts)?;
    if opts.unique && !opts.sort {
        eprintln!(
            "UNIQUE=true requires sorting but SORT=false was specified.  Results will be sorted."
        );
    }
    let result = opts.action.act(first, second)?;
    if opts.invert {
        opts.sort = false;
        opts.unique = true;
    }
    if opts.action == Action::Union {
        opts.sort = true;
        opts.unique = true;
    }
    // The output header is the action's, before the sort sets its order.
    let mut header = result.header.clone();
    let possibly_sorted = if opts.sort {
        result.into_sorted()
    } else {
        result
    };
    let possibly_inverted = if opts.invert {
        invert(&possibly_sorted)
    } else {
        possibly_sorted
    };
    let mut intervals = if opts.unique {
        unique_intervals(&possibly_inverted, !opts.dont_merge_abutting, true)
    } else {
        possibly_inverted.intervals
    };
    if opts.break_bands > 0 {
        intervals = break_at_band_multiples(&intervals, opts.break_bands)?;
    }
    header.add_program(
        &format!("IntervalListTools {}", args.join(" ")),
        "IntervalListTools",
    );
    for comment in &opts.comments {
        header.add_comment(comment);
    }
    let output = IntervalList { header, intervals };
    if let Some(content) = opts.scatter_content {
        opts.scatter_count = scatter_count_for_content(opts.mode, &output, content);
    }
    let (base_count, interval_count) = match &opts.output {
        Some(dir) if opts.scatter_count > 1 => {
            let lists = scatter(&output, opts.mode, opts.scatter_count)?;
            write_scatter(dir, &lists, opts.scatter_count, &inputs)?;
            eprintln!(
                "Wrote {} scatter subdirectories to {}.",
                lists.len(),
                dir.display()
            );
            (lists.base_count(), lists.interval_count())
        }
        Some(path) => {
            write_list(&output, path, &inputs)?;
            (output.base_count(), output.len() as i64)
        }
        None => (output.base_count(), output.len() as i64),
    };
    let value = match opts.output_value.as_str() {
        "BASES" => Some(base_count),
        "INTERVALS" => Some(interval_count),
        _ => None,
    };
    if let Some(value) = value {
        match &opts.count_output {
            Some(path) => {
                std::fs::write(path, format!("{value}\n")).map_err(|source| Error::Output {
                    path: path.clone(),
                    source,
                })?
            }
            None => println!("{value}"),
        }
    }
    Ok(())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Err(error) = run(&args) {
        eprintln!("interval_ops: {error}");
        std::process::exit(1);
    }
}

//! The tool registry: Picard tools by name, with Picard-style arguments.
//!
//! A tool is a [`Tool`] constant in its own module: the Picard name, the
//! options it takes and a `run` function from [`Invocation`] to
//! [`ToolOutput`]. Adding it to [`TOOLS`] is the only edit outside that
//! module; the command line (`bammetrics <ToolName> I=... O=...`) and the
//! Python binding then both serve it, because both go through
//! [`Tool::execute`]. A tool that reads the BAM once also gives a
//! `collector`, which `CollectMultipleMetrics` runs beside the other
//! programs' collectors on one read of the file.
//!
//! Arguments are the forms Picard 3.5.0 accepts: `KEY=VALUE`, `--KEY VALUE`,
//! `--KEY=VALUE` and `-K VALUE`, with long or short names. An option no
//! entry names is an error, because dropping an option Picard applies would
//! change the numbers silently.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use noodles::sam;

use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::pass::{self, Accumulator, Collector, Limits, default_threads};

/// What an option does once recognised.
#[derive(Clone, Copy)]
pub enum Kind {
    /// Read by the tool through [`Invocation`].
    Value,
    /// The metrics file path ([`Invocation::output`]).
    Output,
    /// Never changes the metrics: accepted, dropped, noted on stderr.
    Ignored(&'static str),
    /// A [`COMMON`] option this tool does not take (Picard's parser
    /// rejects it as unrecognized).
    Rejected,
}

/// One option: its names (the long Picard name first) and kind.
#[derive(Clone, Copy)]
pub struct Opt {
    pub names: &'static [&'static str],
    pub kind: Kind,
}

impl Opt {
    pub const fn value(names: &'static [&'static str]) -> Self {
        Self {
            names,
            kind: Kind::Value,
        }
    }

    /// The option naming the metrics file. A tool that declares one takes
    /// no `OUTPUT` (Picard's `CollectDuplicateMetrics` writes to `M`).
    pub const fn output(names: &'static [&'static str]) -> Self {
        Self {
            names,
            kind: Kind::Output,
        }
    }

    pub const fn ignored(names: &'static [&'static str], why: &'static str) -> Self {
        Self {
            names,
            kind: Kind::Ignored(why),
        }
    }

    /// A common option Picard rejects for this tool.
    pub const fn rejected(names: &'static [&'static str]) -> Self {
        Self {
            names,
            kind: Kind::Rejected,
        }
    }
}

/// Options every Picard tool takes from `CommandLineProgram` and
/// `SinglePassSamProgram`. `ASSUME_SORTED=false` fails a file whose header
/// is not coordinate sorted ([`Invocation::single_pass`]); the others change
/// no metric.
const COMMON: &[Opt] = &[
    Opt::value(&["INPUT", "I"]),
    Opt::output(&["OUTPUT", "O"]),
    Opt::ignored(&["VALIDATION_STRINGENCY"], "records are not validated"),
    Opt::ignored(&["TMP_DIR"], "no temporary files are written"),
    Opt::ignored(&["QUIET"], "progress is not logged"),
    Opt::ignored(&["VERBOSITY"], "progress is not logged"),
    Opt::ignored(&["COMPRESSION_LEVEL"], "no BAM is written"),
    Opt::ignored(&["MAX_RECORDS_IN_RAM"], "records are streamed"),
    Opt::ignored(&["CREATE_INDEX"], "no BAM is written"),
    Opt::ignored(&["CREATE_MD5_FILE"], "no BAM is written"),
    Opt::ignored(&["GA4GH_CLIENT_SECRETS"], "no GA4GH access"),
    Opt::ignored(&["USE_JDK_DEFLATER", "USE_JDK_INFLATER"], "no JDK codec"),
    Opt::value(&["ASSUME_SORTED", "AS"]),
    Opt::ignored(
        &["REFERENCE_SEQUENCE", "R"],
        "this metric does not read the reference",
    ),
];

/// Builds a tool's one-pass [`Collector`] from its parsed options and the
/// BAM header (Picard's `SinglePassSamProgram.setup`).
pub type CollectorFn = fn(&Invocation, &sam::Header) -> Result<Box<dyn Collector>>;

/// A registered Picard tool.
pub struct Tool {
    /// The Picard class name, e.g. `CollectInsertSizeMetrics`.
    pub name: &'static str,
    /// Options beyond [`COMMON`]; these shadow a common option of the same name.
    pub options: &'static [Opt],
    /// Compute the metrics: the `OUTPUT` file (`None` where Picard writes
    /// none and logs a warning, which `run` prints) and any further files.
    pub run: fn(&Invocation) -> Result<ToolOutput>,
    /// For a tool that reads the BAM once: its collector, so
    /// `CollectMultipleMetrics` can run it on a shared read. `run` is then
    /// [`Invocation::run_collector`] on it.
    pub collector: Option<CollectorFn>,
}

/// Every registered tool. One line per tool.
pub const TOOLS: &[&Tool] = &[
    &crate::insert_size::TOOL,
    &crate::alignment_summary::TOOL,
    &crate::quality_yield::TOOL,
    &crate::quality_score_distribution::TOOL,
    &crate::mean_quality_by_cycle::TOOL,
    &crate::base_distribution_by_cycle::TOOL,
    &crate::duplicate_metrics::TOOL,
    &crate::accumulate_quality_yield::TOOL,
    &crate::gc_bias::TOOL,
    &crate::multiple_metrics::TOOL,
    &crate::library_complexity::TOOL,
    &crate::wgs_metrics::TOOL,
    &crate::wgs_metrics::RAW_TOOL,
    &crate::wgs_metrics::NON_ZERO_TOOL,
    &crate::bed_to_interval_list::TOOL,
    &crate::hs_metrics::TOOL,
    &crate::hs_metrics::PCR_TOOL,
];

/// The tool named `name`.
pub fn find(name: &str) -> Option<&'static Tool> {
    TOOLS.iter().copied().find(|t| t.name == name)
}

/// A parsed command line for one tool.
#[derive(Debug, Default)]
pub struct Invocation {
    pub input: PathBuf,
    pub output: Option<PathBuf>,
    /// BGZF decompression workers.
    pub threads: usize,
    /// Values by long option name, one entry per occurrence.
    values: HashMap<&'static str, Vec<String>>,
    /// `StringHeader` command line: the tool and each option as given.
    pub command: String,
    /// One line per ignored option, for standard error.
    pub notes: Vec<String>,
}

impl Invocation {
    /// The last value given for `name` (a long option name).
    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name)?.last().map(String::as_str)
    }

    /// Every value given for a repeatable option, in order.
    pub fn all(&self, name: &str) -> &[String] {
        self.values.get(name).map_or(&[], Vec::as_slice)
    }

    /// The parsed last value, or `default` when absent.
    pub fn parse<T: FromStr>(&self, name: &str, default: T) -> Result<T> {
        self.parse_opt(name).map(|v| v.unwrap_or(default))
    }

    /// The parsed last value, `None` when absent.
    pub fn parse_opt<T: FromStr>(&self, name: &str) -> Result<Option<T>> {
        self.get(name)
            .map(|v| {
                v.parse()
                    .map_err(|_| Error::Usage(format!("{name}={v}: not a valid value")))
            })
            .transpose()
    }

    /// A Picard collection option: each value given adds to `default`, and
    /// `null` empties the list so far.
    pub fn list<T: PartialEq>(
        &self,
        name: &str,
        default: Vec<T>,
        parse: impl Fn(&str) -> Result<T>,
    ) -> Result<Vec<T>> {
        picard_list(default, self.all(name), parse)
    }

    /// `STOP_AFTER` (0 when the tool has no such option) and `ASSUME_SORTED`
    /// (default true) as given.
    pub fn limits(&self) -> Result<Limits> {
        let stop_after: i64 = self.parse("STOP_AFTER", 0)?;
        Ok(Limits {
            stop_after: u64::try_from(stop_after).unwrap_or(0),
            assume_sorted: self.flag("ASSUME_SORTED", true)?,
        })
    }

    /// Read `INPUT` once with [`Invocation::limits`] applied, the
    /// `SinglePassSamProgram` loop.
    pub fn single_pass<A, F>(&self, build: F) -> Result<A::Output>
    where
        A: Accumulator,
        F: FnOnce(&sam::Header) -> Result<A>,
    {
        pass::run_limited(&self.input, self.threads, self.limits()?, build)
    }

    /// [`Invocation::single_pass`] over a tool's collector: the `run` of
    /// every tool that has one.
    pub fn run_collector(&self, build: CollectorFn) -> Result<ToolOutput> {
        self.single_pass(|header| build(self, header))
    }

    /// A Picard boolean (`true` or `false`, any case).
    pub fn flag(&self, name: &str, default: bool) -> Result<bool> {
        match self.get(name).map(str::to_ascii_lowercase).as_deref() {
            None => Ok(default),
            Some("true") => Ok(true),
            Some("false") => Ok(false),
            Some(v) => Err(Error::Usage(format!("{name}={v}: expected true or false"))),
        }
    }
}

impl Tool {
    fn lookup(&self, name: &str) -> Option<&'static Opt> {
        let own_output = self.options.iter().any(|o| matches!(o.kind, Kind::Output));
        self.options
            .iter()
            .chain(
                COMMON
                    .iter()
                    .filter(|o| !(own_output && matches!(o.kind, Kind::Output))),
            )
            .find(|o| o.names.contains(&name))
    }

    /// Long name of the option that names the metrics file.
    pub fn output_name(&self) -> &'static str {
        self.options
            .iter()
            .chain(COMMON)
            .find(|o| matches!(o.kind, Kind::Output))
            .map_or("OUTPUT", |o| o.names[0])
    }

    /// Parse `args` (Picard syntax) into an [`Invocation`].
    pub fn parse(&self, args: &[String]) -> Result<Invocation> {
        let mut inv = Invocation {
            threads: default_threads(),
            command: self.name.to_owned(),
            ..Invocation::default()
        };
        let mut input = None;
        let mut rest = args.iter();
        while let Some(arg) = rest.next() {
            let dashed = arg.starts_with('-');
            let body = arg.trim_start_matches('-');
            let (name, value) = match body.split_once('=') {
                Some((n, v)) => (n, v.to_owned()),
                None if dashed => {
                    let v = rest.next().ok_or_else(|| {
                        Error::Usage(format!("{}: {arg} needs a value", self.name))
                    })?;
                    (body, v.clone())
                }
                None => {
                    return Err(Error::Usage(format!(
                        "{}: expected NAME=VALUE, got {arg:?}",
                        self.name
                    )));
                }
            };
            if name == "threads" {
                inv.threads = value
                    .parse()
                    .map_err(|_| Error::Usage(format!("threads={value}: not a number")))?;
                continue;
            }
            let opt = self.lookup(name).ok_or_else(|| {
                Error::Usage(format!(
                    "{}: option {name} is not implemented (it may change the metrics)",
                    self.name
                ))
            })?;
            if matches!(opt.kind, Kind::Rejected) {
                return Err(Error::Usage(format!(
                    "{}: Unrecognized option: {name}",
                    self.name
                )));
            }
            inv.command.push_str(&format!(" {}={value}", opt.names[0]));
            match opt.kind {
                Kind::Ignored(why) => inv.notes.push(format!("ignoring {name}={value}: {why}")),
                Kind::Rejected => unreachable!("rejected above"),
                Kind::Output => inv.output = Some(PathBuf::from(&value)),
                Kind::Value => {
                    if opt.names[0] == "INPUT" {
                        input = Some(PathBuf::from(&value));
                    }
                    inv.values.entry(opt.names[0]).or_default().push(value);
                }
            }
        }
        inv.input =
            input.ok_or_else(|| Error::Usage(format!("{}: INPUT is required", self.name)))?;
        Ok(inv)
    }

    /// Parse and compute. The invocation is returned for its output path
    /// and header.
    pub fn execute(&self, args: &[String]) -> Result<(Invocation, ToolOutput)> {
        let inv = self.parse(args)?;
        let output = (self.run)(&inv)?;
        Ok((inv, output))
    }

    /// The command-line entry: compute and write the Picard metrics file to
    /// `OUTPUT`, or standard output without one, and every further file to
    /// its own path.
    pub fn run_cli(&self, args: &[String]) -> Result<()> {
        let (inv, output) = self.execute(args)?;
        for note in &inv.notes {
            eprintln!("bammetrics: {note}");
        }
        write_files(&inv, &output)
    }
}

/// Write a tool's files: `OUTPUT` (standard output without a path) and the
/// extra files, each with the invocation's `StringHeader` lines.
pub fn write_files(inv: &Invocation, output: &ToolOutput) -> Result<()> {
    if let Some(main) = &output.main {
        let path = inv.output.as_deref();
        write_to(path, |out| main.write(main.headers(&inv.command), out))?;
    }
    for (path, file) in &output.extra {
        write_to(Some(path), |out| {
            file.write(file.headers(&inv.command), out)
        })?;
    }
    Ok(())
}

/// The values of a Picard collection option applied to its default: each
/// value is added, `null` empties the list so far. A value already in the
/// list is added again; the options that use this only test membership.
pub fn picard_list<T, S: AsRef<str>>(
    default: Vec<T>,
    values: &[S],
    parse: impl Fn(&str) -> Result<T>,
) -> Result<Vec<T>> {
    let mut list = default;
    for value in values {
        if value.as_ref() == "null" {
            list.clear();
        } else {
            list.push(parse(value.as_ref())?);
        }
    }
    Ok(list)
}

/// Run `write` against `path`, or standard output when `path` is `None`.
pub fn write_to(
    path: Option<&Path>,
    write: impl FnOnce(&mut BufWriter<Box<dyn Write>>) -> io::Result<()>,
) -> Result<()> {
    let output_error = |source| Error::Output {
        path: path.map_or_else(|| PathBuf::from("-"), Path::to_path_buf),
        source,
    };
    let sink: Box<dyn Write> = match path {
        Some(p) => Box::new(File::create(p).map_err(output_error)?),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);
    write(&mut out).map_err(output_error)?;
    out.flush().map_err(output_error)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn tool() -> &'static Tool {
        find("CollectInsertSizeMetrics").expect("registered")
    }

    #[test]
    fn picard_and_dashed_forms_agree() {
        let forms = [
            strings(&["I=x.bam", "O=out.txt", "M=0.5"]),
            strings(&["--INPUT", "x.bam", "--OUTPUT=out.txt", "-M", "0.5"]),
            strings(&["--I=x.bam", "-O", "out.txt", "--MINIMUM_PCT", "0.5"]),
        ];
        for args in forms {
            let inv = tool().parse(&args).expect("parses");
            assert_eq!(inv.input, Path::new("x.bam"));
            assert_eq!(inv.output.as_deref(), Some(Path::new("out.txt")));
            assert_eq!(inv.get("MINIMUM_PCT"), Some("0.5"));
            assert_eq!(
                inv.command,
                "CollectInsertSizeMetrics INPUT=x.bam OUTPUT=out.txt MINIMUM_PCT=0.5"
            );
        }
    }

    #[test]
    fn repeated_options_keep_every_value() {
        let inv = tool()
            .parse(&strings(&[
                "I=x.bam",
                "LEVEL=SAMPLE",
                "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
            ]))
            .expect("parses");
        assert_eq!(
            inv.all("METRIC_ACCUMULATION_LEVEL"),
            ["SAMPLE", "READ_GROUP"]
        );
    }

    #[test]
    fn plumbing_is_noted_and_unknown_options_are_errors() {
        let inv = tool()
            .parse(&strings(&["I=x.bam", "QUIET=true", "R=ref.fa"]))
            .expect("parses");
        assert_eq!(inv.notes.len(), 2);
        let err = tool()
            .parse(&strings(&["I=x.bam", "STOP_AFTER=10"]))
            .expect_err("unknown");
        assert!(err.to_string().contains("STOP_AFTER"));
        let err = tool().parse(&strings(&["bam"])).expect_err("bare word");
        assert!(err.to_string().contains("NAME=VALUE"));
        let err = tool().parse(&strings(&["O=o.txt"])).expect_err("no input");
        assert!(err.to_string().contains("INPUT"));
    }

    #[test]
    fn typed_getters_name_the_bad_option() {
        let inv = tool()
            .parse(&strings(&["I=x.bam", "INCLUDE_DUPLICATES=maybe", "W=abc"]))
            .expect("parses");
        assert!(inv.flag("INCLUDE_DUPLICATES", false).is_err());
        let err = inv
            .parse_opt::<i64>("HISTOGRAM_WIDTH")
            .expect_err("not a number");
        assert!(err.to_string().contains("HISTOGRAM_WIDTH"));
        assert!(inv.flag("MISSING", true).expect("default"));
    }

    #[test]
    fn a_tool_with_its_own_output_option_takes_no_output() {
        let dup = find("CollectDuplicateMetrics").expect("registered");
        assert_eq!(dup.output_name(), "METRICS_FILE");
        assert_eq!(tool().output_name(), "OUTPUT");
        let inv = dup
            .parse(&strings(&["I=x.bam", "M=d.txt", "AS=false"]))
            .expect("parses");
        assert_eq!(inv.output.as_deref(), Some(Path::new("d.txt")));
        assert!(!inv.flag("ASSUME_SORTED", true).expect("bool"));
        let err = dup
            .parse(&strings(&["I=x.bam", "O=d.txt"]))
            .expect_err("no OUTPUT");
        assert!(err.to_string().contains("option O "));
    }
}

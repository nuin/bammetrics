//! The tool registry: Picard tools by name, with Picard-style arguments.
//!
//! A tool is a [`Tool`] constant in its own module: the Picard name, the
//! options it takes and a `run` function from [`Invocation`] to
//! [`MetricsOutput`]. Adding it to [`TOOLS`] is the only edit outside that
//! module; the command line (`bammetrics <ToolName> I=... O=...`) and the
//! Python binding then both serve it, because both go through
//! [`Tool::execute`].
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

use crate::error::{Error, Result};
use crate::metric::MetricsOutput;
use crate::metrics_file::headers;
use crate::pass::default_threads;

/// What an option does once recognised.
#[derive(Clone, Copy)]
pub enum Kind {
    /// Read by the tool through [`Invocation`].
    Value,
    /// Never changes the metrics: accepted, dropped, noted on stderr.
    Ignored(&'static str),
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

    pub const fn ignored(names: &'static [&'static str], why: &'static str) -> Self {
        Self {
            names,
            kind: Kind::Ignored(why),
        }
    }
}

/// Options every Picard tool takes from `CommandLineProgram` and
/// `SinglePassSamProgram`; none changes a metric.
const COMMON: &[Opt] = &[
    Opt::value(&["INPUT", "I"]),
    Opt::value(&["OUTPUT", "O"]),
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
    Opt::ignored(
        &["ASSUME_SORTED", "AS"],
        "one pass over the file in any order",
    ),
    Opt::ignored(
        &["REFERENCE_SEQUENCE", "R"],
        "this metric does not read the reference",
    ),
];

/// A registered Picard tool.
pub struct Tool {
    /// The Picard class name, e.g. `CollectInsertSizeMetrics`.
    pub name: &'static str,
    /// Options beyond [`COMMON`]; these shadow a common option of the same name.
    pub options: &'static [Opt],
    /// Compute the metrics. `None` means Picard writes no file (it logs a
    /// warning instead, which `run` prints).
    pub run: fn(&Invocation) -> Result<Option<MetricsOutput>>,
}

/// Every registered tool. One line per tool.
pub const TOOLS: &[&Tool] = &[&crate::insert_size::TOOL, &crate::alignment_summary::TOOL];

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
        self.options
            .iter()
            .chain(COMMON)
            .find(|o| o.names.contains(&name))
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
            inv.command.push_str(&format!(" {}={value}", opt.names[0]));
            match opt.kind {
                Kind::Ignored(why) => inv.notes.push(format!("ignoring {name}={value}: {why}")),
                Kind::Value => {
                    if opt.names[0] == "INPUT" {
                        input = Some(PathBuf::from(&value));
                    } else if opt.names[0] == "OUTPUT" {
                        inv.output = Some(PathBuf::from(&value));
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
    pub fn execute(&self, args: &[String]) -> Result<(Invocation, Option<MetricsOutput>)> {
        let inv = self.parse(args)?;
        let output = (self.run)(&inv)?;
        Ok((inv, output))
    }

    /// The command-line entry: compute and write the Picard metrics file to
    /// `OUTPUT`, or standard output without one.
    pub fn run_cli(&self, args: &[String]) -> Result<()> {
        let (inv, output) = self.execute(args)?;
        for note in &inv.notes {
            eprintln!("bammetrics: {note}");
        }
        let Some(output) = output else {
            return Ok(());
        };
        let path = inv.output.as_deref();
        write_to(path, |out| output.write(headers(inv.command.clone()), out))
    }
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
}

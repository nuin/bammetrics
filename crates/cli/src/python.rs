//! The option tables of the Python `qc`, `check` and `report` (docs/05
//! section 5). Each keyword argument is a clap argument of [`QcArgs`],
//! [`CheckArgs`] or [`ReportArgs`], named by its long flag in snake case
//! (`--fail-on` is `fail_on`); [`parse`] runs the command's own parser on
//! the argv the keywords spell, so a Python call and a command line are
//! checked by one table. [`python_stubs`] writes the `bammetrics.pyi`
//! entries from the same arguments.

use std::any::TypeId;
use std::path::{Path, PathBuf};

use bammetrics_core::multiple_metrics::file_extensions;
use clap::{ArgAction, Args, Command, FromArgMatches};

use crate::check::CheckArgs;
use crate::error::VerbError;
use crate::manifest::{Output, StepRecord};
use crate::qc::QcArgs;
use crate::report::ReportArgs;
use crate::verbs;

/// How a keyword turns into arguments.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Shape {
    /// `True` gives `--flag`, `False` nothing.
    Flag,
    /// One value: `--long=VALUE`.
    Value,
    /// Pairs, repeatable: `--hs BAITS TARGETS`.
    Pairs,
    /// One value or several: `--long=A,B`.
    List,
    /// The positional argument, one or several values.
    Positional { many: bool },
}

/// One keyword argument of a Python verb function.
#[derive(Clone, Debug)]
pub struct Keyword {
    /// The Python name: the long flag in snake case.
    pub name: String,
    /// The long flag, without `--`; empty for the positional argument.
    pub long: String,
    pub shape: Shape,
    /// The `.pyi` annotation.
    pub annotation: String,
    /// The `.pyi` default, as Python source; empty for a required one.
    pub default: String,
}

impl Keyword {
    /// True when the values are paths (`_Path` in the annotation): a str or
    /// an `os.PathLike`, never another object turned into text.
    pub fn takes_path(&self) -> bool {
        self.annotation.contains("_Path")
    }
}

/// A Python verb function: the clap verb it runs and what it adds or
/// leaves out.
pub struct Verb {
    pub name: &'static str,
    command: fn() -> Command,
    /// Flags with no meaning in Python (output format, progress line).
    skipped: &'static [&'static str],
    /// Python-only keywords: name, annotation, default.
    extra: &'static [(&'static str, &'static str, &'static str)],
    returns: &'static str,
    doc: &'static str,
}

fn qc_command() -> Command {
    QcArgs::augment_args(Command::new("qc"))
}

fn check_command() -> Command {
    CheckArgs::augment_args(Command::new("check"))
}

fn report_command() -> Command {
    ReportArgs::augment_args(Command::new("report"))
}

/// The Python verb functions, in `.pyi` order.
pub const VERBS: [Verb; 3] = [
    Verb {
        name: "qc",
        command: qc_command,
        skipped: &["format", "quiet"],
        extra: &[("report", "bool", "True")],
        returns: "Run[Verdict]",
        doc: "bammetrics qc: the collectors the preset needs on one read of the BAM, the verdict \
              and the page, the files the command writes. report=False leaves out \
              PREFIX.report.html. A failed verdict is a value; VerdictError only when fail_on is \
              given and reached. An incomplete verdict returns the Run (its exit_code is 11). A \
              failed collector raises InputError (its .run holds the rest).",
    },
    Verb {
        name: "check",
        command: check_command,
        skipped: &["format", "output"],
        extra: &[],
        returns: "Verdict | list[Verdict]",
        doc: "bammetrics check: verdicts from directories of metrics files, verdict.json files or \
              qc prefixes. One path with one sample gives a Verdict, else a list. VerdictError \
              when a verdict reaches fail_on (fail_on=None never raises). An incomplete verdict \
              is returned (its exit_code(fail_on) is 11).",
    },
    Verb {
        name: "report",
        command: report_command,
        skipped: &[],
        extra: &[],
        returns: "list[str]",
        doc: "bammetrics report: SAMPLE.verdict.json and SAMPLE.report.html per sample and \
              index.html, under output (default: the one directory given). Returns the paths \
              written.",
    },
];

/// The Python verb named `name`.
pub fn verb(name: &str) -> Option<&'static Verb> {
    VERBS.iter().find(|v| v.name == name)
}

/// `fail-on` -> `fail_on`.
fn snake(long: &str) -> String {
    long.replace('-', "_")
}

/// The `.pyi` type of one value of `arg`.
fn value_type(arg: &clap::Arg) -> String {
    let values = arg.get_possible_values();
    if !values.is_empty() {
        let names: Vec<String> = values
            .iter()
            .map(|v| format!("\"{}\"", v.get_name()))
            .collect();
        return format!("Literal[{}]", names.join(", "));
    }
    let id = arg.get_value_parser().type_id();
    if id == TypeId::of::<PathBuf>() {
        "_Path".to_owned()
    } else if id == TypeId::of::<usize>() || id == TypeId::of::<u64>() {
        "int".to_owned()
    } else {
        "str".to_owned()
    }
}

impl Verb {
    /// The clap command, built, so every argument has its parser.
    fn built(&self) -> Command {
        let mut command = (self.command)().no_binary_name(true);
        command.build();
        command
    }

    /// The keyword arguments: the positional argument first, then every
    /// flag in clap order, then the Python-only ones.
    pub fn keywords(&self) -> Vec<Keyword> {
        let command = self.built();
        let mut keywords = Vec::new();
        for arg in command.get_arguments() {
            let id = arg.get_id().as_str();
            if id == "help" || id == "version" {
                continue;
            }
            let many = arg.get_num_args().is_some_and(|n| n.max_values() > 1)
                || matches!(arg.get_action(), ArgAction::Append);
            if arg.is_positional() {
                let one = value_type(arg);
                let annotation = if many {
                    format!("{one} | Sequence[{one}]")
                } else {
                    one
                };
                keywords.push(Keyword {
                    name: snake(id),
                    long: String::new(),
                    shape: Shape::Positional { many },
                    annotation,
                    default: String::new(),
                });
                continue;
            }
            let Some(long) = arg.get_long() else {
                continue;
            };
            if self.skipped.contains(&long) {
                continue;
            }
            let (shape, annotation, default) = if matches!(arg.get_action(), ArgAction::SetTrue) {
                (Shape::Flag, "bool".to_owned(), "False".to_owned())
            } else if arg.get_num_args().is_some_and(|n| n.min_values() == 2) {
                (
                    Shape::Pairs,
                    "tuple[_Path, _Path] | Sequence[tuple[_Path, _Path]] | None".to_owned(),
                    "None".to_owned(),
                )
            } else if arg.get_value_delimiter().is_some() {
                let ty = value_type(arg);
                (
                    Shape::List,
                    format!("{ty} | Sequence[{ty}] | None"),
                    "None".to_owned(),
                )
            } else {
                let ty = value_type(arg);
                // A text default prints as the command's; a number the
                // command computes (threads) stays None.
                let default = match arg.get_default_values() {
                    [one] if ty != "int" => format!("\"{}\"", one.to_string_lossy()),
                    _ => "None".to_owned(),
                };
                (Shape::Value, format!("{ty} | None"), default)
            };
            keywords.push(Keyword {
                name: snake(long),
                long: long.to_owned(),
                shape,
                annotation,
                default,
            });
        }
        for (name, annotation, default) in self.extra {
            keywords.push(Keyword {
                name: (*name).to_owned(),
                long: String::new(),
                shape: Shape::Flag,
                annotation: (*annotation).to_owned(),
                default: (*default).to_owned(),
            });
        }
        keywords
    }

    /// The command's arguments from `argv` (flags, then `--` and the
    /// positional values), through the command's parser. A clap error is
    /// the usage error the command prints, `bammetrics qc: ...`.
    pub fn parse<T: FromArgMatches>(&self, argv: &[String]) -> Result<T, VerbError> {
        let matches = self
            .built()
            .try_get_matches_from(argv)
            .map_err(|e| crate::clap_error(self.name, &e))?;
        T::from_arg_matches(&matches).map_err(|e| crate::clap_error(self.name, &e))
    }

    /// The `.pyi` entry: the positional argument, `*`, every keyword.
    pub fn stub(&self) -> String {
        let mut text = format!("def {}(\n", self.name);
        let mut star = false;
        for keyword in self.keywords() {
            if matches!(keyword.shape, Shape::Positional { .. }) {
                text.push_str(&format!("    {}: {},\n", keyword.name, keyword.annotation));
                continue;
            }
            if !star {
                text.push_str("    *,\n");
                star = true;
            }
            text.push_str(&format!(
                "    {}: {} = {},\n",
                keyword.name, keyword.annotation, keyword.default
            ));
        }
        text.push_str(&format!(") -> {}:\n", self.returns));
        let mut doc = bammetrics_core::registry::wrap(self.doc, "    ", 88);
        doc[0] = format!("    \"\"\"{}", doc[0].trim_start());
        let last = doc.len() - 1;
        doc[last].push_str("\"\"\"");
        for line in doc {
            text.push_str(&line);
            text.push('\n');
        }
        text
    }
}

/// First line of the generated verb block in `bammetrics.pyi`.
pub const STUB_BEGIN: &str = "# BEGIN verbs: generated from crates/cli/src/python.rs (UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs)";
/// Last line of the block.
pub const STUB_END: &str = "# END verbs";

/// The names the Python `collect(metrics=[...])` takes: the verb of every
/// tool `qc` can run, in [`crate::qc::collector_tools`] order, then `hs` as
/// `bammetrics collect --metric hs` spells it.
pub fn collector_names() -> Vec<&'static str> {
    let mut names: Vec<&str> = crate::qc::collector_tools()
        .iter()
        .filter_map(|t| verbs::for_class(t).map(|v| v.name))
        .collect();
    names.push("hs");
    names
}

/// The generated verb block of `bammetrics.pyi`.
pub fn python_stubs() -> String {
    let mut text = format!(
        "{STUB_BEGIN}\n# Keywords are the command's long flags in snake case; None leaves an option \
         at the command's default.\n"
    );
    text.push_str(
        "\n# The metric names with a dash: collect(metrics=[...]) given one runs the collectors \
         as qc does.\n_CollectorVerb = Literal[\n",
    );
    for name in collector_names().into_iter().filter(|n| n.contains('-')) {
        text.push_str(&format!("    \"{name}\",\n"));
    }
    text.push_str("]\n");
    for verb in &VERBS {
        text.push('\n');
        text.push_str(&verb.stub());
    }
    text.push('\n');
    text.push_str(STUB_END);
    text.push('\n');
    text
}

/// The key of a metrics file in `Run.metrics`: the verb of the tool that
/// wrote it (`insert-size`), then `.BAIT_SET` for one of several hs runs,
/// then `.summary` for a tool's further file (`gc-bias.summary`).
pub fn metric_key(record: &StepRecord, output: &Output) -> String {
    let mut key = verbs::for_class(&record.tool).map_or_else(
        || bammetrics_core::registry::snake_case(&record.tool),
        |v| v.name.to_owned(),
    );
    if let Some((_, set)) = record.id.split_once('.') {
        key.push('.');
        key.push_str(set);
    }
    let name = output.path.to_string_lossy();
    let further = file_extensions(&record.tool)
        .unwrap_or_default()
        .iter()
        .skip(1)
        .find(|(_, ext)| name.ends_with(ext));
    if let Some((option, _)) = further {
        let option = option.to_ascii_lowercase();
        key.push('.');
        key.push_str(option.strip_suffix("_output").unwrap_or(&option));
    }
    key
}

/// `PREFIX.verdict.json` for a path that names no file or directory but
/// is the prefix of a `qc` run (`check("out/s")`); the path itself
/// otherwise.
pub fn check_input(path: &Path) -> PathBuf {
    if path.exists() {
        return path.to_path_buf();
    }
    let verdict = crate::qc::verdict_path(path);
    if verdict.is_file() {
        verdict
    } else {
        path.to_path_buf()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::manifest::Status;

    /// `bammetrics.pyi` holds [`python_stubs`] between its markers.
    /// `UPDATE_PYI=1` rewrites the block.
    #[test]
    fn python_stubs() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../bammetrics.pyi");
        let text = std::fs::read_to_string(&path).expect("bammetrics.pyi");
        let generated = super::python_stubs();
        let (begin, end) = match (text.find(STUB_BEGIN), text.find(STUB_END)) {
            (Some(b), Some(e)) => (b, e + STUB_END.len() + 1),
            _ => (text.len(), text.len()),
        };
        if std::env::var_os("UPDATE_PYI").is_some() {
            let updated = format!("{}{generated}{}", &text[..begin], &text[end..]);
            std::fs::write(&path, updated).expect("write bammetrics.pyi");
            return;
        }
        assert_eq!(
            text.get(begin..end).unwrap_or_default(),
            generated,
            "bammetrics.pyi is stale: UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs"
        );
    }

    #[test]
    fn qc_keywords_are_the_long_flags_in_snake_case() {
        let names: Vec<String> = verb("qc")
            .expect("qc")
            .keywords()
            .into_iter()
            .map(|k| k.name)
            .collect();
        for expected in [
            "bam",
            "output",
            "reference",
            "preset",
            "hs",
            "panel",
            "ref_flat",
            "fail_on",
            "resume",
            "force",
            "threads",
            "sample",
            "all_plots",
            "date",
            "report",
        ] {
            assert!(names.iter().any(|n| n == expected), "{expected}: {names:?}");
        }
        assert!(!names.iter().any(|n| n == "format" || n == "quiet"));
    }

    #[test]
    fn arguments_parse_through_the_command_and_errors_read_as_its_own() {
        let qc = verb("qc").expect("qc");
        let args: QcArgs = qc
            .parse(&[
                "--fail-on=warn".to_owned(),
                "--".to_owned(),
                "x.bam".to_owned(),
            ])
            .expect("parses");
        assert_eq!(args.bam, PathBuf::from("x.bam"));
        assert!(args.write_page);
        let Err(err) = qc.parse::<QcArgs>(&["--fail-on=maybe".to_owned(), "x.bam".to_owned()])
        else {
            panic!("--fail-on=maybe parsed");
        };
        assert_eq!(err.code, 2);
        assert!(
            err.to_string()
                .starts_with("bammetrics qc: invalid value \"maybe\" (from --fail-on)"),
            "{err}"
        );
    }

    #[test]
    fn metric_keys_name_the_verb_and_the_further_file() {
        let record = |tool: &str, id: &str| StepRecord {
            id: id.to_owned(),
            stage: "pass".to_owned(),
            tool: tool.to_owned(),
            key: crate::manifest::Key {
                argv_hash: String::new(),
                tool_version: String::new(),
                preset_hash: String::new(),
                input: PathBuf::new(),
                fingerprint: crate::manifest::Fingerprint {
                    size: 0,
                    mtime: String::new(),
                    edges: String::new(),
                },
                files: Vec::new(),
            },
            argv: Vec::new(),
            outputs: Vec::new(),
            status: Status::Ok,
            error: None,
        };
        let output = |path: &str| Output {
            path: PathBuf::from(path),
            tool: String::new(),
            kind: "metrics".to_owned(),
        };
        let gc = record("CollectGcBiasMetrics", "CollectGcBiasMetrics");
        assert_eq!(
            metric_key(&gc, &output("o/s.gc_bias.detail_metrics")),
            "gc-bias"
        );
        assert_eq!(
            metric_key(&gc, &output("o/s.gc_bias.summary_metrics")),
            "gc-bias.summary"
        );
        let hs = record("CollectHsMetrics", "CollectHsMetrics.exome");
        assert_eq!(
            metric_key(&hs, &output("o/s.exome.hs_metrics")),
            "hs-metrics.exome"
        );
        let dup = record("CollectDuplicateMetrics", "CollectDuplicateMetrics");
        assert_eq!(
            metric_key(&dup, &output("o/s.duplicate_metrics")),
            "duplicate-metrics"
        );
    }
}

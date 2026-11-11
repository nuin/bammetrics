//! The plan file of `bammetrics chain --plan FILE` (docs/chaining.md): the
//! source, the chain's options and the steps of one chain as TOML, read
//! through [`bammetrics_core::toml::parse_tables`] into the [`Spec`] the
//! command line makes.
//!
//! ```toml
//! [chain]
//! source = "in.bam"
//! output-prefix = "out/NA12878"
//!
//! [[step]]
//! verb = "mark-duplicates"
//!
//! [[step]]
//! verb = "insert-size"
//! minimum-pct = 0.1
//! level = ["SAMPLE", "LIBRARY"]
//!
//! [[step]]
//! verb = "write"
//! index = true
//! ```
//!
//! A key of `[chain]` is the long name of a chain option. A key of a
//! `[[step]]` is `verb`, or the long name of an option of that verb (its
//! Picard name works too). A repeatable option takes a list. Every key is
//! checked here, with the line it is on; an unknown one is an error with
//! the closest name.

use std::collections::HashSet;
use std::path::Path;

use bammetrics_core::registry::{Ty, did_you_mean};
use bammetrics_core::toml::{Entry, Scalar, Table, parse_tables};

use super::args::{find_step, step_parser_flags, unknown_step, write_flags};
use super::{Spec, StepSpec, VERB, WRITE};
use crate::error::VerbError;
use crate::grammar::{self, Flag};

/// The keys of `[chain]`.
const CHAIN_KEYS: [&str; 10] = [
    "source",
    "output-prefix",
    "reference",
    "threads",
    "tmp-dir",
    "format",
    "also",
    "sample",
    "quiet",
    "verbose",
];

/// A plan file and the line an error is on.
struct At<'a> {
    file: &'a str,
}

impl At<'_> {
    fn error(&self, line: usize, what: impl AsRef<str>) -> VerbError {
        VerbError::usage(VERB, format!("{}:{line}: {}", self.file, what.as_ref()))
    }
}

/// The text of a scalar that is one value of an option.
fn scalar_text(value: &Scalar) -> Option<String> {
    match value {
        Scalar::Str(text) => Some(text.clone()),
        Scalar::Int(n) => Some(n.to_string()),
        Scalar::Float(x) => Some(x.to_string()),
        Scalar::Bool(b) => Some(b.to_string()),
        Scalar::List(_) => None,
    }
}

/// True when `value` is what an option of type `ty` takes.
fn fits(ty: Ty, value: &Scalar) -> bool {
    match (ty, value) {
        (Ty::Bool, Scalar::Bool(_))
        | (Ty::Int, Scalar::Int(_))
        | (Ty::Float, Scalar::Int(_) | Scalar::Float(_))
        | (Ty::Text | Ty::Path, Scalar::Str(_)) => true,
        // A text option takes a number as its text (`program-record-id = 7`).
        (Ty::Text, Scalar::Int(_) | Scalar::Float(_)) => true,
        _ => false,
    }
}

fn wanted(ty: Ty) -> &'static str {
    match ty {
        Ty::Bool => "true or false",
        Ty::Int => "an integer",
        Ty::Float => "a number",
        Ty::Text | Ty::Path => "a string",
    }
}

/// A string, or a list of strings joined by commas: `format` and `also`.
fn name_list(at: &At, entry: &Entry) -> Result<String, VerbError> {
    let names: Option<Vec<&str>> = match &entry.value {
        Scalar::Str(text) => Some(vec![text.as_str()]),
        Scalar::List(items) => items
            .iter()
            .map(|item| match item {
                Scalar::Str(text) => Some(text.as_str()),
                _ => None,
            })
            .collect(),
        _ => None,
    };
    names.map(|names| names.join(",")).ok_or_else(|| {
        at.error(
            entry.line,
            format!(
                "{}: expected a string or a list of strings, got {}",
                entry.key,
                entry.value.kind()
            ),
        )
    })
}

fn string(at: &At, entry: &Entry) -> Result<String, VerbError> {
    match &entry.value {
        Scalar::Str(text) => Ok(text.clone()),
        other => Err(at.error(
            entry.line,
            format!("{}: expected a string, got {}", entry.key, other.kind()),
        )),
    }
}

fn boolean(at: &At, entry: &Entry) -> Result<bool, VerbError> {
    match &entry.value {
        Scalar::Bool(value) => Ok(*value),
        other => Err(at.error(
            entry.line,
            format!(
                "{}: expected true or false, got {}",
                entry.key,
                other.kind()
            ),
        )),
    }
}

/// A key given twice in one table is an error, as in TOML.
fn check_no_repeat<'t>(
    at: &At,
    seen: &mut HashSet<&'t str>,
    entry: &'t Entry,
) -> Result<(), VerbError> {
    if seen.insert(entry.key.as_str()) {
        return Ok(());
    }
    Err(at.error(entry.line, format!("{} is given twice", entry.key)))
}

/// The `[chain]` table into `spec`.
fn chain_table(at: &At, table: &Table, spec: &mut Spec) -> Result<(), VerbError> {
    let mut seen = HashSet::new();
    for entry in &table.entries {
        check_no_repeat(at, &mut seen, entry)?;
        match entry.key.as_str() {
            "source" => spec.source = Some(string(at, entry)?),
            "output-prefix" => spec.prefix = Some(string(at, entry)?),
            "reference" => spec.reference = Some(string(at, entry)?),
            "tmp-dir" => spec.tmp_dir = Some(string(at, entry)?),
            "sample" => spec.sample = Some(string(at, entry)?),
            "format" => spec.formats = Some(name_list(at, entry)?),
            "also" => spec.also = Some(name_list(at, entry)?),
            "quiet" => spec.quiet = boolean(at, entry)?,
            "verbose" => spec.verbose = boolean(at, entry)?,
            "threads" => {
                let threads = match &entry.value {
                    Scalar::Int(n) => usize::try_from(*n).ok(),
                    _ => None,
                };
                spec.threads = Some(threads.ok_or_else(|| {
                    at.error(entry.line, "threads: expected a whole number, 0 or more")
                })?);
            }
            other => {
                let mut error = at.error(entry.line, format!("[chain]: unknown key {other}"));
                if let Some(best) = did_you_mean(other, CHAIN_KEYS) {
                    error = error.hint(format!("did you mean {best}?"));
                }
                return Err(error.hint(format!("keys of [chain]: {}", CHAIN_KEYS.join(", "))));
            }
        }
    }
    Ok(())
}

/// One option of a step as arguments: `--long=VALUE`, once per item of a
/// list.
fn option_args(at: &At, label: &str, flag: &Flag, entry: &Entry) -> Result<Vec<String>, VerbError> {
    let bad = |what: String| at.error(entry.line, format!("{label}: {}: {what}", entry.key));
    let items: &[Scalar] = match &entry.value {
        Scalar::List(_) if !flag.repeated => {
            return Err(bad("takes one value, got a list".to_owned()));
        }
        Scalar::List(items) => items,
        one => std::slice::from_ref(one),
    };
    items
        .iter()
        .map(|item| {
            let text = scalar_text(item).filter(|_| fits(flag.ty, item));
            let text = text
                .ok_or_else(|| bad(format!("expected {}, got {}", wanted(flag.ty), item.kind())))?;
            Ok(format!("--{}={text}", flag.long))
        })
        .collect()
}

/// One `[[step]]` table as a step: `number` counts from 1.
fn step_table(at: &At, table: &Table, number: usize) -> Result<StepSpec, VerbError> {
    let verb_entry = table
        .entries
        .iter()
        .find(|entry| entry.key == "verb")
        .ok_or_else(|| {
            at.error(
                table.line,
                format!("step {number}: missing key verb (a tool verb or write)"),
            )
        })?;
    let name = string(at, verb_entry)?;
    let flags = if name == WRITE {
        write_flags()
    } else {
        let verb = find_step(&name).ok_or_else(|| {
            let mut error = unknown_step(&name, None);
            error.what = format!("{}:{}: {}", at.file, verb_entry.line, error.what);
            error
        })?;
        step_parser_flags(verb)
    };
    let label = format!("step {number} ({name})");
    let mut args = Vec::new();
    let mut seen = HashSet::new();
    for entry in &table.entries {
        check_no_repeat(at, &mut seen, entry)?;
        if entry.key == "verb" {
            continue;
        }
        let Some(flag) = grammar::lookup(&flags, &entry.key) else {
            let mut error = at.error(entry.line, format!("{label}: unknown key {}", entry.key));
            if let Some(best) = grammar::closest(&flags, &entry.key) {
                error = error.hint(format!("did you mean {best}?"));
            }
            return Err(error.hint(format!(
                "every option: bammetrics chain BAM + {name} --help"
            )));
        };
        args.extend(option_args(at, &label, flag, entry)?);
    }
    Ok(StepSpec {
        verb: name,
        args,
        origin: Some(format!("{}:{}", at.file, table.line)),
    })
}

/// The plan in `text`; `file` names it in errors.
pub fn parse(text: &str, file: &str) -> Result<Spec, VerbError> {
    let at = At { file };
    let tables = parse_tables(text).map_err(|e| at.error(e.line, e.message))?;
    let mut spec = Spec::default();
    let mut chain_seen = false;
    for table in &tables {
        match (table.name.as_str(), table.array) {
            ("chain", false) if chain_seen => {
                return Err(at.error(table.line, "[chain] is given twice"));
            }
            ("chain", false) => {
                chain_seen = true;
                chain_table(&at, table, &mut spec)?;
            }
            ("step", true) => {
                let number = spec.steps.len() + 1;
                spec.steps.push(step_table(&at, table, number)?);
            }
            ("chain", true) => {
                return Err(at.error(table.line, "[[chain]]: one chain per file, write [chain]"));
            }
            ("step", false) => {
                return Err(at.error(table.line, "[step]: steps are a list, write [[step]]"));
            }
            (other, _) => {
                let mut error = at.error(table.line, format!("unknown table {other}"));
                if let Some(best) = did_you_mean(other, ["chain", "step"]) {
                    error = error.hint(format!("did you mean {best}?"));
                }
                return Err(
                    error.hint("a plan file holds one [chain] table and one [[step]] per step")
                );
            }
        }
    }
    Ok(spec)
}

/// The plan file at `path`.
pub fn read(path: &Path) -> Result<Spec, VerbError> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| VerbError::input(VERB, format!("{}: {e}", path.display())).from("--plan"))?;
    parse(&text, &path.display().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    const PLAN: &str = "# the chain\n[chain]\nsource = \"in.bam\"\noutput-prefix = \"out/s\"\n\
        reference = \"ref.fa\"\nthreads = 4\nformat = [\"picard\", \"json\"]\nalso = \"csv\"\n\
        quiet = true\n\n[[step]]\nverb = \"mark-duplicates\"\nremove-duplicates = true\n\
        comment = [\"a b\", \"c\"]\nas = \"dedup\"\n\n[[step]]\nverb = \"insert-size\"\n\
        MINIMUM_PCT = 0.1\nlevel = [\"SAMPLE\"]\ndeviations = 5\n\n[[step]]\n\
        verb = \"write\"\noutput = \"out/s.bam\"\nindex = true\ncompression-level = 9\n";

    #[test]
    fn a_plan_is_the_spec_of_its_command_line() {
        let spec = parse(PLAN, "plan.toml").expect("parses");
        assert_eq!(spec.source.as_deref(), Some("in.bam"));
        assert_eq!(spec.prefix.as_deref(), Some("out/s"));
        assert_eq!(spec.reference.as_deref(), Some("ref.fa"));
        assert_eq!(spec.threads, Some(4));
        assert_eq!(spec.formats.as_deref(), Some("picard,json"));
        assert_eq!(spec.also.as_deref(), Some("csv"));
        assert!(spec.quiet && !spec.verbose);
        let steps: Vec<(&str, Vec<&str>, Option<&str>)> = spec
            .steps
            .iter()
            .map(|step| {
                (
                    step.verb.as_str(),
                    step.args.iter().map(String::as_str).collect(),
                    step.origin.as_deref(),
                )
            })
            .collect();
        assert_eq!(
            steps,
            [
                (
                    "mark-duplicates",
                    vec![
                        "--remove-duplicates=true",
                        "--comment=a b",
                        "--comment=c",
                        "--as=dedup"
                    ],
                    Some("plan.toml:11")
                ),
                (
                    "insert-size",
                    vec![
                        "--minimum-pct=0.1",
                        "--metric-accumulation-level=SAMPLE",
                        "--deviations=5"
                    ],
                    Some("plan.toml:17")
                ),
                (
                    "write",
                    vec![
                        "--output=out/s.bam",
                        "--index=true",
                        "--compression-level=9"
                    ],
                    Some("plan.toml:23")
                ),
            ]
        );
    }

    #[test]
    fn unknown_names_are_errors_with_the_closest_name() {
        let cases = [
            (
                "[chain]\nsorce = \"in.bam\"\n",
                "bammetrics chain: p.toml:2: [chain]: unknown key sorce\n  did you mean source?",
            ),
            (
                "[chian]\n",
                "bammetrics chain: p.toml:1: unknown table chian\n  did you mean chain?",
            ),
            (
                "[[step]]\nverb = \"inser-size\"\n",
                "bammetrics chain: p.toml:2: unknown step \"inser-size\"\n  did you mean \
                 insert-size?",
            ),
            (
                "[[step]]\nverb = \"insert-size\"\nminmum-pct = 0.1\n",
                "bammetrics chain: p.toml:3: step 1 (insert-size): unknown key minmum-pct\n  \
                 did you mean minimum-pct (Picard MINIMUM_PCT)?",
            ),
            (
                "[[step]]\nverb = \"write\"\nindx = true\n",
                "bammetrics chain: p.toml:3: step 1 (write): unknown key indx\n  did you mean \
                 index?",
            ),
        ];
        for (text, start) in cases {
            let error = parse(text, "p.toml").expect_err(text);
            assert_eq!(error.code, 2);
            assert!(error.to_string().starts_with(start), "{error}");
        }
    }

    #[test]
    fn values_are_checked_against_the_option() {
        let cases = [
            (
                "[chain]\nsource = 3\n",
                "p.toml:2: source: expected a string, got an integer",
            ),
            (
                "[chain]\nthreads = \"4\"\n",
                "p.toml:2: threads: expected a whole number",
            ),
            (
                "[chain]\nquiet = \"yes\"\n",
                "p.toml:2: quiet: expected true or false",
            ),
            (
                "[chain]\nsource = \"a\"\nsource = \"b\"\n",
                "p.toml:3: source is given twice",
            ),
            ("[chain]\n[chain]\n", "p.toml:2: [chain] is given twice"),
            (
                "[step]\n",
                "p.toml:1: [step]: steps are a list, write [[step]]",
            ),
            (
                "[[step]]\nindex = true\n",
                "p.toml:1: step 1: missing key verb",
            ),
            (
                "[[step]]\nverb = \"insert-size\"\nminimum-pct = \"x\"\n",
                "p.toml:3: step 1 (insert-size): minimum-pct: expected a number, got a string",
            ),
            (
                "[[step]]\nverb = \"insert-size\"\nminimum-pct = [0.1, 0.2]\n",
                "p.toml:3: step 1 (insert-size): minimum-pct: takes one value, got a list",
            ),
            (
                "[[step]]\nverb = \"write\"\nindex = 1\n",
                "p.toml:3: step 1 (write): index: expected true or false, got an integer",
            ),
            (
                "[[step]]\nverb = unquoted\n",
                "p.toml:2: verb: expected a quoted string",
            ),
        ];
        for (text, needle) in cases {
            let error = parse(text, "p.toml").expect_err(text).to_string();
            assert!(error.contains(needle), "{text:?}: {error}");
        }
    }

    #[test]
    fn a_missing_file_is_an_input_error() {
        let error = read(Path::new("tests/data/no-such-plan.toml")).expect_err("absent");
        assert_eq!(error.code, 1);
        assert!(error.to_string().ends_with("(from --plan)"), "{error}");
    }
}

//! The Python face of the chain (docs/python.md): the keyword tables of
//! `bammetrics.chain` and `bammetrics.step`, the key of a file in
//! `Run.metrics`, and the `bammetrics.pyi` entries.
//!
//! Every keyword is a flag of the command: [`chain_keywords`] reads
//! [`global_flags`], [`step_table`] reads the tool's Python signature
//! ([`python_keywords`], the list the generated registry stub prints) and
//! finds each name in [`step_flags`]. The Python module spells the argv
//! those flags take and hands it to [`super::parse`], so one parser checks
//! a command line, a plan file and a Python call.

use std::path::Path;

use bammetrics_core::registry::{Ty, did_you_mean, python_keywords, snake_case, wrap};

use super::args::{
    ALSO, AS, find_step, global_flags, is_step, is_transform, not_a_step, step_flags, write_flags,
};
use super::prepare::{Chain, Planned, Role};
use super::{WRITE, usage};
use crate::error::VerbError;
use crate::grammar::{Flag, Key};
use crate::verbs::{Files, VERBS, Verb};

/// How many values a keyword takes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Many {
    /// One value: `--long=VALUE`.
    One,
    /// One value or a list, each its own occurrence: `--long=A --long=B`.
    Repeated,
    /// One value or a list, joined: `--long=A,B`.
    Joined,
}

/// One keyword argument of `chain` or `step`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Keyword {
    /// The Python name (`minimum_pct`).
    pub name: String,
    /// The long flag it spells, without `--` (`minimum-pct`).
    pub long: String,
    pub ty: Ty,
    pub many: Many,
}

impl Keyword {
    fn of(name: &str, flag: &Flag) -> Self {
        let many = if flag.repeated {
            Many::Repeated
        } else if flag.key == Key::Format || flag.key == ALSO {
            Many::Joined
        } else {
            Many::One
        };
        Self {
            name: name.to_owned(),
            long: flag.long.clone(),
            ty: flag.ty,
            many,
        }
    }

    /// The `.pyi` annotation of a given value.
    fn annotation(&self) -> String {
        let one = match self.ty {
            Ty::Bool => "bool",
            Ty::Int => "int",
            Ty::Float => "float",
            Ty::Text => "str",
            Ty::Path => "_Path",
        };
        match self.many {
            Many::One => one.to_owned(),
            Many::Repeated | Many::Joined => format!("{one} | Sequence[{one}]"),
        }
    }
}

/// The Python name of a chain option: the long flag in snake case, and
/// `output` for `--output-prefix`, as `qc` names its prefix.
fn chain_name(flag: &Flag) -> String {
    if flag.key == Key::Output {
        "output".to_owned()
    } else {
        flag.long.replace('-', "_")
    }
}

/// The chain options with no Python keyword: the source is positional,
/// and nothing is printed.
const NOT_KEYWORDS: [Key; 4] = [Key::Input, Key::Quiet, Key::Verbose, Key::Help];

/// The positional `source` of `chain`.
pub fn source_keyword() -> Keyword {
    let flags = global_flags();
    let flag = flags
        .iter()
        .find(|flag| flag.key == Key::Input)
        .expect("the chain has a source");
    Keyword::of("source", flag)
}

/// The keyword-only arguments of `chain`, in [`global_flags`] order.
pub fn chain_keywords() -> Vec<Keyword> {
    global_flags()
        .iter()
        .filter(|flag| !NOT_KEYWORDS.contains(&flag.key))
        .map(|flag| Keyword::of(&chain_name(flag), flag))
        .collect()
}

/// The error for a keyword `known` does not hold, under `label`.
fn unknown_keyword(label: &str, name: &str, known: &[Keyword]) -> VerbError {
    let names: Vec<&str> = known.iter().map(|k| k.name.as_str()).collect();
    let mut error = usage(format!("{label}unknown keyword {name:?}"));
    if let Some(best) = did_you_mean(name, names.iter().copied()) {
        error = error.hint(format!("did you mean {best}?"));
    }
    error.hint(format!("keywords: {}", names.join(", ")))
}

/// The keyword `name` of `chain`, or the usage error with the closest
/// name.
pub fn chain_keyword(name: &str) -> Result<Keyword, VerbError> {
    let known = chain_keywords();
    match known.iter().find(|k| k.name == name) {
        Some(keyword) => Ok(keyword.clone()),
        None => Err(unknown_keyword("", name, &known)),
    }
}

/// What `step(verb, ...)` takes for one verb.
#[derive(Clone, Debug)]
pub struct StepTable {
    /// The verb as the command spells it: `insert-size`, or `write`.
    pub verb: &'static str,
    pub role: Role,
    pub keywords: Vec<Keyword>,
    /// A keyword of the tool's function that a step refuses, with the
    /// reason: the input every step gets from the chain, the BAM a
    /// transform leaves to `write`.
    refused: Vec<(&'static str, String)>,
}

impl StepTable {
    /// The keyword `name`, or the usage error with the closest name.
    pub fn keyword(&self, name: &str) -> Result<&Keyword, VerbError> {
        if let Some(keyword) = self.keywords.iter().find(|k| k.name == name) {
            return Ok(keyword);
        }
        let label = format!("step {}: ", self.verb);
        match self.refused.iter().find(|(refused, _)| *refused == name) {
            Some((_, why)) => Err(usage(format!(
                "{label}{name} is not a keyword of a step: {why}"
            ))),
            None => Err(unknown_keyword(&label, name, &self.keywords)),
        }
    }
}

/// The verbs `step` takes, as Python spells them: every tool verb that
/// can be a step in snake case, in `bammetrics tools` order, then `write`.
pub fn step_names() -> Vec<String> {
    VERBS
        .iter()
        .filter(|verb| is_step(verb))
        .map(|verb| verb.name.replace('-', "_"))
        .chain([WRITE.to_owned()])
        .collect()
}

/// The verb `name` refers to: the verb in snake or kebab case
/// (`insert_size`, `insert-size`), its Picard class, or the name of its
/// Python function (`collect_insert_size_metrics`).
fn find_verb(name: &str) -> Option<&'static Verb> {
    let function = name.replace('-', "_");
    find_step(&name.replace('_', "-"))
        .or_else(|| find_step(name))
        .or_else(|| VERBS.iter().find(|v| snake_case(v.class) == function))
}

/// The error for a name that is no verb: the closest step, in the
/// spelling given (snake case unless the name holds a dash).
fn unknown_step(name: &str) -> VerbError {
    let names = step_names();
    let wanted = name.replace('-', "_");
    let mut error = usage(format!("unknown step {name:?}"));
    if let Some(best) = did_you_mean(&wanted, names.iter().map(String::as_str)) {
        let best = if name.contains('-') {
            best.replace('_', "-")
        } else {
            best.to_owned()
        };
        error = error.hint(format!("did you mean {best}?"));
    }
    error.hint(format!("steps: {}", names.join(", ")))
}

/// The verb flag `--panel FILE` of `hs-metrics`: baits and targets in one
/// file.
const PANEL: Key = Key::Own("panel");

/// The keywords of a tool step: those of the tool's Python function that
/// a step can take, then `panel` on `hs_metrics`, `format`, and `as_` on a
/// transform.
fn tool_keywords(verb: &Verb) -> (Vec<Keyword>, Vec<(&'static str, String)>) {
    let flags = step_flags(verb);
    let tool = verb.tool();
    let mut refused = vec![("input", "the chain's source feeds every step".to_owned())];
    let mut keywords = Vec::new();
    for keyword in python_keywords(tool) {
        // A transform that writes a BAM as a tool writes none in a chain.
        if verb.files == Files::Bam && keyword.picard == "OUTPUT" {
            refused.push((
                "output",
                format!(
                    "{} writes no BAM in a chain, add step(\"write\", output=...) after it",
                    verb.name
                ),
            ));
            continue;
        }
        let key = if keyword.output && verb.files != Files::Bam {
            Key::Output
        } else {
            Key::Picard(keyword.picard)
        };
        if let Some(flag) = flags.iter().find(|flag| flag.key == key) {
            keywords.push(Keyword::of(&keyword.name, flag));
        }
    }
    // `--panel` of `hs-metrics` is no row of the tool's table, so its
    // function has no such keyword; the verb's missing-intervals error
    // names it, and `qc` takes `panel=` too.
    let own = [(PANEL, "panel"), (Key::Format, "format"), (AS, "as_")];
    for (key, name) in own {
        if let Some(flag) = flags.iter().find(|flag| flag.key == key) {
            keywords.push(Keyword::of(name, flag));
        }
    }
    (keywords, refused)
}

/// What `step(name, ...)` takes, or the usage error for a name that is no
/// step.
pub fn step_table(name: &str) -> Result<StepTable, VerbError> {
    if name == WRITE {
        let keywords = write_flags()
            .iter()
            .filter(|flag| flag.key != Key::Help)
            .map(|flag| Keyword::of(&flag.long.replace('-', "_"), flag))
            .collect();
        return Ok(StepTable {
            verb: WRITE,
            role: Role::Write,
            keywords,
            refused: Vec::new(),
        });
    }
    let verb = find_verb(name).ok_or_else(|| unknown_step(name))?;
    if !is_step(verb) {
        return Err(not_a_step(verb, &format!("step {}", verb.name))
            .hint(format!("steps: {}", step_names().join(", "))));
    }
    let role = if is_transform(verb) {
        Role::Transform
    } else {
        Role::Tap
    };
    let (keywords, refused) = tool_keywords(verb);
    Ok(StepTable {
        verb: verb.name,
        role,
        keywords,
        refused,
    })
}

/// The part of a metrics file's name that tells it from the step's other
/// files: the file option for a tool's further file (`summary` for
/// `SUMMARY_OUTPUT`), or what follows the prefix for a tool that names its
/// files when it finishes (`pre_adapter_detail_metrics`).
fn further_file(step: &Planned, path: &Path) -> Option<String> {
    let name = path.to_string_lossy();
    if let Some(prefix) = &step.late_prefix {
        let rest = name.strip_prefix(prefix.to_string_lossy().as_ref())?;
        let rest = rest.strip_prefix('.').unwrap_or(rest);
        return (!rest.is_empty()).then(|| rest.to_owned());
    }
    let tool = step.tool?;
    crate::python::further_file(tool.name, path)
}

/// The key of a metrics file in the `Run.metrics` of a chain: the step's
/// verb (`insert-size`), with `STREAM.` in front for a tap that does not
/// read the last stream (`in.insert-size` in front of a transform, as its
/// file is `PREFIX.in.insert_size_metrics`), and `.FILE` after it for a
/// step's further file (`gc-bias.summary`).
pub fn metric_key(chain: &Chain, step: &Planned, path: &Path) -> String {
    let mut key = String::new();
    if step.role != Role::Transform && step.reads != chain.last_stream() {
        key.push_str(&step.reads);
        key.push('.');
    }
    key.push_str(step.verb);
    if let Some(further) = further_file(step, path) {
        key.push('.');
        key.push_str(&further);
    }
    key
}

/// The docstring of `step`.
pub const STEP_DOC: &str = "One step of chain(): a tool verb in snake case (insert_size), its \
     Picard class, or write. The keywords are those of the tool's function (minimum_pct=0.1) \
     without its input, plus output (the step's own file, or prefix for a tool that writes \
     several), format, panel on hs_metrics (baits and targets in one file), and as_ on a \
     transform (the name of the stream it leaves). write takes output, index, md5, \
     compression_level and output_format (\"bam\" or \"sam\"). An unknown verb or keyword raises UsageError here, with the closest \
     name; so does a float that is NaN or infinite, and a list of texts with an item that \
     is not a str.";

/// The docstring of `chain`.
pub const CHAIN_DOC: &str = "bammetrics chain: several steps on one BAM in one process, the \
     records handed from step to step in memory. steps holds verb names (insert_size) or \
     step(...) objects, in order: mark_duplicates is a transform, a one-pass tool is a tap on \
     the records of its position, write writes them as a BAM. plan names a TOML plan file in \
     place of steps; a keyword given replaces the file's. Returns a Run: files, metrics (keyed \
     by verb, STREAM.verb in front of a transform) and stages (one Stage per step). \
     dry_run=True returns the plan text and runs nothing. A failed step raises InputError \
     after the other steps wrote their files (its .run holds the Run).";

/// A docstring as `.pyi` lines under a `def`.
fn doc_lines(doc: &str) -> String {
    let mut lines = wrap(doc, "    ", 88);
    lines[0] = format!("    \"\"\"{}", lines[0].trim_start());
    let last = lines.len() - 1;
    lines[last].push_str("\"\"\"");
    lines.join("\n") + "\n"
}

/// One `chain` overload: `dry_run` is the last keyword, and its value
/// picks the return type.
fn chain_overload(dry_run: &str, returns: &str, doc: Option<&str>) -> String {
    let source = source_keyword();
    let mut text = format!(
        "@overload\ndef chain(\n    {}: {} | None = None,\n    steps: Step | str | \
         Sequence[Step | str] | None = None,\n    *,\n",
        source.name,
        source.annotation()
    );
    for keyword in chain_keywords() {
        if keyword.name == "dry_run" {
            text.push_str(&format!("    dry_run: {dry_run},\n"));
        } else {
            text.push_str(&format!(
                "    {}: {} | None = None,\n",
                keyword.name,
                keyword.annotation()
            ));
        }
    }
    match doc {
        Some(doc) => text.push_str(&format!(") -> {returns}:\n{}", doc_lines(doc))),
        None => text.push_str(&format!(") -> {returns}: ...\n")),
    }
    text
}

/// The `bammetrics.pyi` entries of `step` and `chain`, inside the
/// generated verb block.
pub fn stubs() -> String {
    let mut text = "# The steps of chain(): each tool verb that can be a step, in snake case, \
                    and write.\n_StepVerb = Literal[\n"
        .to_owned();
    for name in step_names() {
        text.push_str(&format!("    \"{name}\",\n"));
    }
    text.push_str("]\n\n");
    text.push_str("def step(verb: _StepVerb | str, /, **options: Any) -> Step:\n");
    text.push_str(&doc_lines(STEP_DOC));
    text.push('\n');
    text.push_str(&chain_overload(
        "Literal[False] = False",
        "Run[None]",
        Some(CHAIN_DOC),
    ));
    text.push_str(&chain_overload("Literal[True]", "str", None));
    text.push_str(&chain_overload("bool", "Run[None] | str", None));
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names(keywords: &[Keyword]) -> Vec<&str> {
        keywords.iter().map(|k| k.name.as_str()).collect()
    }

    #[test]
    fn chain_keywords_are_the_chain_options_in_snake_case() {
        assert_eq!(
            names(&chain_keywords()),
            [
                "output",
                "reference",
                "threads",
                "tmp_dir",
                "format",
                "also",
                "sample",
                "plan",
                "dry_run"
            ]
        );
        let keyword = |name: &str| chain_keyword(name).expect(name);
        assert_eq!(keyword("output").long, "output-prefix");
        assert_eq!(keyword("format").many, Many::Joined);
        assert_eq!(keyword("also").many, Many::Joined);
        assert_eq!(keyword("threads").ty, Ty::Int);
        assert_eq!(keyword("dry_run").ty, Ty::Bool);
        assert_eq!(source_keyword().long, "input");
        let error = chain_keyword("outptu").expect_err("unknown").to_string();
        assert!(
            error.starts_with(
                "bammetrics chain: unknown keyword \"outptu\"\n  did you mean output?"
            ),
            "{error}"
        );
    }

    #[test]
    fn a_step_takes_the_keywords_of_its_tools_function() {
        let table = step_table("insert_size").expect("a step");
        assert_eq!((table.verb, table.role), ("insert-size", Role::Tap));
        assert_eq!(
            names(&table.keywords),
            [
                "deviations",
                "histogram_width",
                "minimum_histogram_width",
                "minimum_pct",
                "include_duplicates",
                "metric_accumulation_level",
                "output",
                "assume_sorted",
                "reference_sequence",
                "format"
            ]
        );
        let level = table.keyword("metric_accumulation_level").expect("level");
        assert_eq!(
            (level.long.as_str(), level.many),
            ("metric-accumulation-level", Many::Repeated)
        );
        let marked = step_table("mark_duplicates").expect("a step");
        assert_eq!(marked.role, Role::Transform);
        let keywords = names(&marked.keywords);
        assert!(keywords.contains(&"metrics_file") && keywords.contains(&"as_"));
        assert!(!keywords.contains(&"output") && !keywords.contains(&"input"));
        let write = step_table("write").expect("a step");
        assert_eq!(
            names(&write.keywords),
            [
                "output",
                "index",
                "md5",
                "compression_level",
                "output_format"
            ]
        );
        // `--panel`, which the error for missing intervals names, is a
        // keyword of the one verb that has the flag.
        let hs = step_table("hs_metrics").expect("a step");
        let panel = hs.keyword("panel").expect("panel");
        assert_eq!(
            (panel.long.as_str(), panel.ty, panel.many),
            ("panel", Ty::Path, Many::One)
        );
        assert!(!names(&table.keywords).contains(&"panel"));
    }

    /// Every keyword of the tool's Python function is a keyword of its
    /// step, apart from the BAM of a transform.
    #[test]
    fn no_keyword_of_a_step_verb_is_lost() {
        for verb in VERBS.iter().filter(|verb| is_step(verb)) {
            let table = step_table(verb.name).expect(verb.name);
            for keyword in python_keywords(verb.tool()) {
                let bam = verb.files == Files::Bam && keyword.picard == "OUTPUT";
                assert_eq!(
                    table.keywords.iter().any(|k| k.name == keyword.name),
                    !bam,
                    "{}: {}",
                    verb.name,
                    keyword.name
                );
            }
        }
    }

    #[test]
    fn a_step_is_named_by_its_verb_its_class_or_its_function() {
        for name in [
            "insert_size",
            "insert-size",
            "CollectInsertSizeMetrics",
            "collect_insert_size_metrics",
        ] {
            assert_eq!(step_table(name).expect(name).verb, "insert-size");
        }
        let cases = [
            (
                "inser_size",
                "bammetrics chain: unknown step \"inser_size\"\n  did you mean insert_size?\n  \
                 steps: insert_size, ",
            ),
            (
                "inser-size",
                "bammetrics chain: unknown step \"inser-size\"\n  did you mean insert-size?",
            ),
            (
                "wgs",
                "bammetrics chain: step wgs: CollectWgsMetrics has no stage and does not read a \
                 BAM in one pass: it cannot be a step\n  wave W14 of docs/07",
            ),
            (
                "fifo_buffer",
                "bammetrics chain: step fifo-buffer: FifoBuffer buffers the pipe between two \
                 commands and a chain has no pipe: it cannot be a step\n  leave the step out, \
                 or put the buffer on a pipe outside the chain: bammetrics chain BAM + write - \
                 | bammetrics fifo-buffer | ...\n  steps: insert_size, ",
            ),
        ];
        for (name, start) in cases {
            let error = step_table(name).expect_err(name);
            assert_eq!(error.code, 2);
            assert!(error.to_string().starts_with(start), "{error}");
        }
        assert!(step_names().ends_with(&[WRITE.to_owned()]));
    }

    #[test]
    fn an_unknown_keyword_gets_the_closest_name_or_the_reason() {
        let table = step_table("insert_size").expect("a step");
        let error = |table: &StepTable, name: &str| {
            let error = table.keyword(name).expect_err(name);
            assert_eq!(error.code, 2);
            error.to_string()
        };
        assert!(error(&table, "minmum_pct").starts_with(
            "bammetrics chain: step insert-size: unknown keyword \"minmum_pct\"\n  did you \
                 mean minimum_pct?"
        ));
        assert_eq!(
            error(&table, "input"),
            "bammetrics chain: step insert-size: input is not a keyword of a step: the chain's \
             source feeds every step"
        );
        let marked = step_table("mark_duplicates").expect("a step");
        assert!(error(&marked, "output").contains("add step(\"write\", output=...) after it"));
    }

    #[test]
    fn the_stub_lists_the_steps_and_the_chain_overloads() {
        let text = stubs();
        assert!(text.contains("    \"mark_duplicates\",\n"));
        assert!(text.contains("    \"write\",\n]\n"));
        assert!(text.contains("def step(verb: _StepVerb | str, /, **options: Any) -> Step:\n"));
        assert!(text.contains("    dry_run: Literal[False] = False,\n) -> Run[None]:\n"));
        assert!(text.contains("    dry_run: Literal[True],\n) -> str: ...\n"));
        assert!(text.contains("    dry_run: bool,\n) -> Run[None] | str: ...\n"));
        assert!(text.contains("    format: str | Sequence[str] | None = None,\n"));
    }
}

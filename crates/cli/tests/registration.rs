//! The registration rule of the full parity program plan: a tool in `TOOLS`
//! is complete only when every place that lists tools lists it. For each
//! registered tool the test asks for
//!
//! 1. a verb in `VERBS`;
//! 2. an argument list in `crates/core/tests/glossary.rs`, the test that
//!    checks every column the tool writes against the glossary, or a line in
//!    that file's `EXEMPT` table. [`glossary_tools`] reads the keys of
//!    `args()` and the first field of each `EXEMPT` line, so a name in a
//!    comment or in an argument value does not count;
//! 3. a concordance module that names it: `tests/concordance/test_<verb>.py`
//!    or `test_<picard_name_in_snake_case>.py`, or the module [`CONCORDANCE`]
//!    gives it, or a line in [`NO_PICARD_ORACLE`] with the reason;
//! 4. its name in `docs/cli.md`;
//! 5. a row of the README tool table that starts with its name.
//!
//! Every gap of every tool is reported in one failure, each with the edit
//! that closes it. A tool added to `TOOLS` alone fails with five lines.

use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_cli::verbs;
use bammetrics_core::registry::{self, TOOLS};

/// Tools whose concordance module is named after a group of tools: the
/// tool, then the file under `tests/concordance`. A tool with a module
/// named after itself needs no line.
const CONCORDANCE: &[(&str, &str)] = &[
    ("CollectInsertSizeMetrics", "test_foundation.py"),
    ("CollectQualityYieldMetrics", "test_read_level.py"),
    ("QualityScoreDistribution", "test_read_level.py"),
    ("MeanQualityByCycle", "test_read_level.py"),
    ("CollectBaseDistributionByCycle", "test_read_level.py"),
    ("CollectDuplicateMetrics", "test_read_level.py"),
    ("AccumulateQualityYieldMetrics", "test_read_level.py"),
    ("CollectWgsMetrics", "test_wgs_metrics.py"),
    ("CollectRawWgsMetrics", "test_wgs_metrics.py"),
    (
        "CollectWgsMetricsWithNonZeroCoverage",
        "test_wgs_metrics.py",
    ),
    ("CollectTargetedPcrMetrics", "test_hs_metrics.py"),
    ("CollectRnaSeqMetrics", "test_rna_seq_metrics.py"),
    ("CollectVariantCallingMetrics", "test_variant_tools.py"),
    ("AccumulateVariantCallingMetrics", "test_variant_tools.py"),
    ("GenotypeConcordance", "test_variant_tools.py"),
    (
        "ConvertSequencingArtifactToOxoG",
        "test_sequencing_artifact.py",
    ),
    ("CollectOxoGMetrics", "test_oxog_metrics.py"),
    (
        "CollectJumpingLibraryMetrics",
        "test_jumping_library_metrics.py",
    ),
];

/// Tools Picard 3.5.0 has no program for, each with the test that stands in
/// for the concordance module.
const NO_PICARD_ORACLE: &[(&str, &str)] = &[
    (
        "CollectCoverageGaps",
        "not a Picard tool (docs/06); crates/core/tests/coverage_gaps.rs checks every gap \
         base against the PER_BASE_COVERAGE file of CollectHsMetrics, a concordant tool",
    ),
    (
        "CollectSexCheckMetrics",
        "not a Picard tool (docs/06); the unit tests of crates/core/src/sex_check.rs check \
         the thresholds and the rows, tests/python/test_clinical_checks.py the fixture row",
    ),
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The texts and the directory the rule reads. `glossary` holds the tools
/// the glossary test runs or exempts.
struct Sources {
    glossary: BTreeSet<String>,
    cli_doc: String,
    readme: String,
    concordance: PathBuf,
}

impl Sources {
    fn read() -> Self {
        let text = |path: &str| {
            fs::read_to_string(root().join(path)).unwrap_or_else(|e| panic!("{path}: {e}"))
        };
        Self {
            glossary: glossary_tools(&text("crates/core/tests/glossary.rs")),
            cli_doc: text("docs/cli.md"),
            readme: text("README.md"),
            concordance: root().join("tests/concordance"),
        }
    }
}

/// `source` from the first `start` to the first `end` after it, empty
/// without `start`.
fn section<'a>(source: &'a str, start: &str, end: &str) -> &'a str {
    let Some(at) = source.find(start) else {
        return "";
    };
    let rest = &source[at..];
    &rest[..rest.find(end).unwrap_or(rest.len())]
}

/// The string literals of `code` that open a tuple or a function call: each
/// follows `(` and is followed by `,`, as in `("CollectX", vec![])` and
/// `args.insert("CollectX", ...)`. A macro's first argument
/// (`format!("S={}", path)`) is not one. Lines that start with `//` are
/// skipped.
fn leading_literals(code: &str) -> Vec<String> {
    let code: Vec<&str> = code
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .collect();
    let code = code.join("\n");
    let mut found = Vec::new();
    let mut rest = code.as_str();
    while let Some(open) = rest.find('"') {
        let opens = rest[..open]
            .trim_end()
            .strip_suffix('(')
            .is_some_and(|before| !before.ends_with('!'));
        let body = &rest[open + 1..];
        let mut escaped = false;
        let close = body.char_indices().find(|&(_, c)| {
            let ends = c == '"' && !escaped;
            escaped = c == '\\' && !escaped;
            ends
        });
        let Some((close, _)) = close else {
            break;
        };
        rest = &body[close + 1..];
        if opens && rest.trim_start().starts_with(',') {
            found.push(body[..close].to_owned());
        }
    }
    found
}

/// The tools the glossary test (`crates/core/tests/glossary.rs`) covers:
/// the keys of its `args()` map and the tools of its `EXEMPT` table.
fn glossary_tools(source: &str) -> BTreeSet<String> {
    let args = section(source, "fn args(", "\n}\n");
    let exempt = section(source, "const EXEMPT", "];");
    leading_literals(args)
        .into_iter()
        .chain(leading_literals(exempt))
        .collect()
}

fn is_word(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

/// Whether `text` holds `name` as a whole word, so `CollectWgsMetrics` is
/// not found inside `CollectWgsMetricsWithNonZeroCoverage`.
fn names(text: &str, name: &str) -> bool {
    text.match_indices(name).any(|(at, _)| {
        let before = text[..at].chars().next_back();
        let after = text[at + name.len()..].chars().next();
        !before.is_some_and(is_word) && !after.is_some_and(is_word)
    })
}

/// Whether a table row of the README starts with `name`.
fn has_readme_row(readme: &str, name: &str) -> bool {
    readme.lines().any(|line| {
        line.strip_prefix("| ")
            .and_then(|rest| rest.strip_prefix(name))
            .is_some_and(|rest| !rest.chars().next().is_some_and(is_word))
    })
}

/// The concordance module files a tool may be named in: the one
/// [`CONCORDANCE`] gives it, else the two named after it.
fn modules(name: &str, verb: Option<&str>) -> Vec<String> {
    if let Some((_, module)) = CONCORDANCE.iter().find(|(tool, _)| *tool == name) {
        return vec![(*module).to_owned()];
    }
    let mut files = vec![format!("test_{}.py", registry::snake_case(name))];
    if let Some(verb) = verb {
        files.insert(0, format!("test_{}.py", verb.replace('-', "_")));
    }
    files
}

/// What is missing for the tool `name`, one line per gap.
fn gaps(name: &str, sources: &Sources) -> Vec<String> {
    let mut gaps = Vec::new();
    let verb = verbs::for_class(name).map(|v| v.name);
    if verb.is_none() {
        gaps.push(format!(
            "{name}: no verb; add verb(\"...\", \"{name}\", Files::Single) to VERBS in \
             crates/cli/src/verbs.rs"
        ));
    }
    if !sources.glossary.contains(name) {
        gaps.push(format!(
            "{name}: no glossary argument list; add one to args() in \
             crates/core/tests/glossary.rs, or a line with the reason to its EXEMPT table"
        ));
    }
    let oracle = NO_PICARD_ORACLE.iter().any(|(tool, _)| *tool == name);
    let modules = modules(name, verb);
    let covered = modules.iter().any(|file| {
        fs::read_to_string(sources.concordance.join(file)).is_ok_and(|text| names(&text, name))
    });
    if !oracle && !covered {
        gaps.push(format!(
            "{name}: no concordance module names it; write tests/concordance/{}, or add a \
             line to CONCORDANCE or NO_PICARD_ORACLE in crates/cli/tests/registration.rs",
            modules[0]
        ));
    }
    if !names(&sources.cli_doc, name) {
        gaps.push(format!("{name}: docs/cli.md does not mention it"));
    }
    if !has_readme_row(&sources.readme, name) {
        gaps.push(format!(
            "{name}: no row of the README tool table starts with it"
        ));
    }
    gaps
}

#[test]
fn every_registered_tool_is_registered_everywhere() {
    let sources = Sources::read();
    let gaps: Vec<String> = TOOLS
        .iter()
        .flat_map(|tool| gaps(tool.name, &sources))
        .collect();
    assert!(
        gaps.is_empty(),
        "{} registration gaps:\n{}",
        gaps.len(),
        gaps.join("\n")
    );
    assert!(TOOLS.len() >= 32, "{} tools registered", TOOLS.len());
}

#[test]
fn the_exception_tables_name_registered_tools_once_with_a_reason() {
    for (tool, module) in CONCORDANCE {
        assert!(registry::find(tool).is_some(), "CONCORDANCE: {tool}");
        assert!(
            root().join("tests/concordance").join(module).is_file(),
            "CONCORDANCE: {tool}: no file {module}"
        );
        assert!(
            !NO_PICARD_ORACLE.iter().any(|(name, _)| name == tool),
            "{tool} is in both tables"
        );
    }
    for (tool, reason) in NO_PICARD_ORACLE {
        assert!(registry::find(tool).is_some(), "NO_PICARD_ORACLE: {tool}");
        assert!(reason.len() >= 40, "NO_PICARD_ORACLE: {tool}: say why");
    }
    for table in [CONCORDANCE, NO_PICARD_ORACLE] {
        for (i, (tool, _)) in table.iter().enumerate() {
            assert!(
                !table[..i].iter().any(|(name, _)| name == tool),
                "{tool} is listed twice"
            );
        }
    }
}

#[test]
fn a_tool_registered_nowhere_else_has_five_gaps() {
    let found = gaps("CollectNewThingMetrics", &Sources::read());
    let expected = [
        "no verb",
        "no glossary argument list",
        "no concordance module names it; write tests/concordance/test_collect_new_thing_metrics.py",
        "docs/cli.md does not mention it",
        "no row of the README tool table starts with it",
    ];
    assert_eq!(found.len(), expected.len(), "{found:?}");
    for (gap, start) in found.iter().zip(expected) {
        assert!(
            gap.starts_with(&format!("CollectNewThingMetrics: {start}")),
            "{gap}"
        );
    }
}

#[test]
fn each_surface_is_checked_on_its_own() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("registration");
    fs::create_dir_all(&dir).expect("work directory");
    // `wgs` is the verb of CollectWgsMetrics, whose module CONCORDANCE names.
    fs::write(
        dir.join("test_wgs_metrics.py"),
        "TOOL = 'CollectWgsMetricsWithNonZeroCoverage'\n",
    )
    .expect("module");
    let sources = Sources {
        glossary: glossary_tools("fn args() {\n(\"CollectWgsMetrics\", vec![]),\n}\n"),
        cli_doc: "`wgs` runs CollectWgsMetricsWithNonZeroCoverage.".to_owned(),
        readme: "| CollectWgsMetricsWithNonZeroCoverage | | |\n\
                 The table above lists CollectWgsMetrics.\n"
            .to_owned(),
        concordance: dir.clone(),
    };
    let found = gaps("CollectWgsMetrics", &sources);
    assert_eq!(found.len(), 3, "{found:?}");
    assert!(found[0].contains("no concordance module names it"));
    assert!(found[1].contains("docs/cli.md"));
    assert!(found[2].contains("README"));

    fs::write(
        dir.join("test_wgs_metrics.py"),
        "TOOLS = ['CollectWgsMetrics']\n",
    )
    .expect("module");
    let sources = Sources {
        cli_doc: "`wgs` runs CollectWgsMetrics.".to_owned(),
        readme: "| CollectWgsMetrics (`INTERVALS`) | | `collect_wgs_metrics()` |\n".to_owned(),
        ..sources
    };
    assert_eq!(gaps("CollectWgsMetrics", &sources), Vec::<String>::new());
    fs::remove_dir_all(&dir).expect("remove the work directory");
}

#[test]
fn the_glossary_surface_is_an_args_key_or_an_exempt_line() {
    let source = r#"
//! A tool added to `TOOLS` needs an entry: ("InADocComment", vec![]).
const CALLS: &str = "a \"quoted\" word";
const NOTE: (&str, &str) = ("OutsideBoth", "a tuple that is neither table");
const EXEMPT: &[(&str, &str)] = &[
    // ("ExemptInAComment", "a comment line"),
    ("ExemptTool", "no committed fixture drives it"),
    (
        "ExemptOnItsOwnLine",
        "the reason, which names \"ReasonTool\", is not a tool",
    ),
];
fn args(dir: &Path) -> BTreeMap<&'static str, Vec<String>> {
    for (tool, extra) in [
        ("OneLine", vec![]),
        (
            "SplitTuple",
            vec![reference.clone(), "PROGRAM=AsAValue".to_owned()],
        ),
        ("WithAMacro", vec![format!("MacroArgument={}", dir.display())]),
    ] {
        args.insert(tool, list);
    }
    // The fixture has no @RG line, which ("InAComment", ...) refuses.
    args.insert("Inserted", vec![bam.clone()]);
    args.insert(
        "InsertedSplit",
        vec![format!("I={DATA}/a.bam"), out("AsAnArgument")],
    );
    args
}
fn files() {
    covered.insert("AfterArgs", vec![]);
}
"#;
    let found: Vec<String> = glossary_tools(source).into_iter().collect();
    assert_eq!(
        found,
        [
            "ExemptOnItsOwnLine",
            "ExemptTool",
            "Inserted",
            "InsertedSplit",
            "OneLine",
            "SplitTuple",
            "WithAMacro"
        ]
    );
    assert!(glossary_tools("no tables here").is_empty());

    let real = Sources::read().glossary;
    for tool in TOOLS {
        assert!(real.contains(tool.name), "{}", tool.name);
    }
    for name in &real {
        assert!(registry::find(name).is_some(), "{name} is not a tool");
    }
}

#[test]
fn whole_words_only() {
    assert!(names("TOOL = \"MarkDuplicates\"", "MarkDuplicates"));
    assert!(names("MarkDuplicates", "MarkDuplicates"));
    assert!(!names("MarkDuplicatesWithMateCigar", "MarkDuplicates"));
    assert!(!names("SimpleMarkDuplicates", "MarkDuplicates"));
    assert!(!names("mark_duplicates", "MarkDuplicates"));
    assert!(has_readme_row(
        "| MarkDuplicates (writes) | | |",
        "MarkDuplicates"
    ));
    assert!(!has_readme_row(
        "| MarkDuplicatesWithMateCigar | | |",
        "MarkDuplicates"
    ));
    assert!(!has_readme_row(
        "| run | MarkDuplicates | 1 s |",
        "MarkDuplicates"
    ));
}

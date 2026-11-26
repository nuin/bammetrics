//! The frozen JSON schemas of `docs/schema` (docs/formats.md, "Frozen
//! schemas"): every object a writer emits is valid under the committed
//! schema of its version, and every object schema closes its members, so
//! a writer that emits a member the schema lacks fails here. The schema's
//! title and `schema` constant are the writer's version.
//!
//! Writers: `qc` (verdict.json, the manifest of schema 1, the events of
//! `--format jsonl`), `check --format jsonl`, and `chain` (the manifest of
//! schema 2: a run that keeps a stream, one that ends in qc, a resumed one,
//! one from a pipe). `metrics-1.json` is checked by `crates/formats/tests`.

use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use bammetrics_formats::json::Json;

#[path = "../../formats/tests/validator/mod.rs"]
mod validator;
use validator::validate;

const INSERT_SIZE: &str = "tests/data/insert_size.bam";
const MARKDUP: &str = "tests/data/markdup.bam";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-schemas-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn schema(name: &str) -> Json {
    let text = fs::read_to_string(root().join("docs/schema").join(name)).expect("read the schema");
    Json::parse(&text).expect("the schema is JSON")
}

fn json(path: &Path) -> Json {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    Json::parse(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn valid(schema: &Json, value: &Json, what: &str) {
    validate(schema, value, schema, "$").unwrap_or_else(|e| panic!("{what}: {e}"));
}

/// Every schema that lists `properties` says what other members may be:
/// `additionalProperties` is `false`, or a schema for a row's cells.
fn closed(schema: &Json, path: &str, open: &mut Vec<String>) {
    match schema {
        Json::Object(members) => {
            let has = |key: &str| members.iter().any(|(k, _)| k == key);
            if has("properties") && !has("additionalProperties") {
                open.push(path.to_owned());
            }
            for (key, member) in members {
                closed(member, &format!("{path}.{key}"), open);
            }
        }
        Json::Array(items) => {
            for (i, item) in items.iter().enumerate() {
                closed(item, &format!("{path}[{i}]"), open);
            }
        }
        _ => {}
    }
}

/// The committed schemas, each with the version its writer writes.
fn versions() -> [(&'static str, &'static str); 6] {
    [
        ("metrics-1.json", "bammetrics.metrics/1"),
        ("verdict-1.json", bammetrics_verdict::verdict::SCHEMA),
        ("manifest-1.json", bammetrics_cli::manifest::SCHEMA),
        ("manifest-2.json", bammetrics_cli::manifest::v2::SCHEMA),
        ("events-1.json", "bammetrics events, version 1"),
        ("tools-1.json", "bammetrics.tools/1"),
    ]
}

#[test]
fn each_schema_is_closed_and_names_its_writers_version() {
    for (name, version) in versions() {
        let schema = schema(name);
        assert_eq!(
            schema.get("title").and_then(Json::as_str),
            Some(version),
            "{name}"
        );
        let constant = schema
            .get("properties")
            .and_then(|p| p.get("schema"))
            .and_then(|s| s.get("const"))
            .and_then(Json::as_str);
        if name != "events-1.json" {
            assert_eq!(constant, Some(version), "{name}");
        }
        let mut open = Vec::new();
        closed(&schema, "$", &mut open);
        assert!(open.is_empty(), "{name}: open objects at {open:?}");
        let description = schema
            .get("description")
            .and_then(Json::as_str)
            .unwrap_or("");
        assert!(description.contains("Frozen"), "{name}");
    }
}

/// A member the schema lacks fails validation: the freeze holds.
#[test]
fn a_member_the_schema_lacks_fails() {
    let dir = fresh_dir("member");
    let out = dir.join("s");
    let output = bammetrics(&["qc", INSERT_SIZE, "-o", out.to_str().expect("utf-8")]);
    assert!(output.status.success());
    let verdict = fs::read_to_string(dir.join("s.verdict.json")).expect("verdict");
    let added = verdict.replacen("{\n", "{\n  \"new_member\": 1,\n", 1);
    let schema = schema("verdict-1.json");
    let value = Json::parse(&added).expect("JSON");
    assert!(validate(&schema, &value, &schema, "$").is_err());
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn qc_writes_its_verdict_manifest_and_events_under_their_schemas() {
    let dir = fresh_dir("qc");
    let out = dir.join("s");
    let prefix = out.to_str().expect("utf-8");
    let events = schema("events-1.json");
    let mut lines = 0;
    for resume in [false, true] {
        let mut args = vec!["qc", INSERT_SIZE, "-o", prefix, "--format", "jsonl"];
        if resume {
            args.push("--resume");
        }
        let output = bammetrics(&args);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        for line in String::from_utf8_lossy(&output.stdout).lines() {
            valid(&events, &Json::parse(line).expect("a JSON line"), line);
            lines += 1;
        }
    }
    assert!(lines > 10, "{lines} events");
    valid(
        &schema("verdict-1.json"),
        &json(&dir.join("s.verdict.json")),
        "verdict",
    );
    valid(
        &schema("manifest-1.json"),
        &json(&dir.join("s.manifest.json")),
        "manifest",
    );
    // check's lines, and the two lines of a run that stops on an error.
    let check = bammetrics(&["check", dir.to_str().expect("utf-8"), "--format", "jsonl"]);
    let failed = bammetrics(&["qc", "tests/data/nope.bam", "--format", "jsonl"]);
    assert_eq!(failed.status.code(), Some(1));
    let both = [check.stdout, failed.stdout].concat();
    let text = String::from_utf8_lossy(&both);
    assert!(text.contains("\"event\":\"verdict\"") && text.contains("\"step\":null"));
    for line in text.lines() {
        valid(&events, &Json::parse(line).expect("a JSON line"), line);
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn chain_writes_its_manifest_under_schema_2() {
    let dir = fresh_dir("chain");
    let manifest = schema("manifest-2.json");
    let prefix = dir.join("s");
    let prefix = prefix.to_str().expect("utf-8");
    let bam = dir.join("out.bam");
    let bam = bam.to_str().expect("utf-8");
    let chain = [
        "chain",
        MARKDUP,
        "-o",
        prefix,
        "--resume",
        "--format",
        "picard,json",
        "+",
        "mark-duplicates",
        "--keep",
        "+",
        "insert-size",
        "+",
        "write",
        bam,
        "--md5",
        "+",
        "qc",
        "--also",
        "markdown",
    ];
    for _ in 0..2 {
        let output = bammetrics(&chain);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        valid(&manifest, &json(&dir.join("s.manifest.json")), "chain");
        valid(
            &schema("verdict-1.json"),
            &json(&dir.join("s.verdict.json")),
            "qc step",
        );
    }
    // A piped source: no fingerprint, a kept stream all the same.
    let piped = dir.join("p");
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["chain", "-", "-o", piped.to_str().expect("utf-8")])
        .args([
            "+",
            "sort-sam",
            "-SO",
            "queryname",
            "--keep",
            "+",
            "quality-yield",
        ])
        .current_dir(root())
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn");
    let bytes = fs::read(root().join(MARKDUP)).expect("the source");
    child
        .stdin
        .take()
        .expect("stdin")
        .write_all(&bytes)
        .expect("feed the pipe");
    let output = child.wait_with_output().expect("wait");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let written = json(&dir.join("p.manifest.json"));
    valid(&manifest, &written, "piped");
    assert_eq!(
        written.get("source").and_then(|s| s.get("size")),
        Some(&Json::Null)
    );
    let _ = fs::remove_dir_all(&dir);
}

fn tools_listing() -> Json {
    let out = bammetrics(&["tools", "--format", "json"]);
    assert_eq!(out.status.code(), Some(0));
    Json::parse(std::str::from_utf8(&out.stdout).expect("utf-8")).expect("JSON")
}

fn array<'a>(value: &'a Json, key: &str) -> &'a [Json] {
    value.get(key).and_then(Json::as_array).expect(key)
}

fn text<'a>(value: &'a Json, key: &str) -> &'a str {
    value.get(key).and_then(Json::as_str).expect(key)
}

fn first_name(option: &Json) -> &str {
    array(option, "names")[0].as_str().expect("name")
}

#[test]
fn tools_listing_is_valid_under_tools_1() {
    let doc = tools_listing();
    valid(&schema("tools-1.json"), &doc, "tools --format json");
}

#[test]
fn tools_listing_equals_the_registry() {
    use bammetrics_core::registry::{Kind, Ty};
    let doc = tools_listing();
    let tools = array(&doc, "tools");
    assert_eq!(tools.len(), 77);
    assert_eq!(tools.len(), bammetrics_cli::verbs::VERBS.len());
    for (entry, verb) in tools.iter().zip(bammetrics_cli::verbs::VERBS) {
        assert_eq!(text(entry, "verb"), verb.name);
        assert_eq!(text(entry, "class"), verb.class);
        let tool = verb.tool();
        let rows = tool.table_with_rejected();
        // Every Tool::table() row and every rejected row, once each.
        let table = tool.table();
        assert!(
            table
                .iter()
                .all(|t| rows.iter().any(|r| r.names == t.names)),
            "{}",
            verb.name
        );
        for row in &rows {
            assert_eq!(rows.iter().filter(|r| r.names == row.names).count(), 1);
        }
        let options = array(entry, "options");
        assert_eq!(options.len(), rows.len(), "{}", verb.name);
        for (option, row) in options.iter().zip(rows) {
            let names: Vec<&str> = array(option, "names")
                .iter()
                .map(|n| n.as_str().expect("name"))
                .collect();
            assert_eq!(
                names[..row.names.len()],
                *row.names,
                "{}: the Picard names come first",
                verb.name
            );
            let kind = match row.kind {
                Kind::Value => "value",
                Kind::Output => "output",
                Kind::Ignored(_) => "ignored",
                Kind::Rejected => "rejected",
            };
            assert_eq!(text(option, "kind"), kind, "{} {names:?}", verb.name);
            let ty = match row.ty {
                Ty::Bool => "bool",
                Ty::Int => "int",
                Ty::Float => "float",
                Ty::Text => "text",
                Ty::Path => "path",
            };
            assert_eq!(text(option, "type"), ty, "{} {names:?}", verb.name);
            assert_eq!(option.get("repeated"), Some(&Json::Bool(row.repeated)));
            assert_eq!(option.get("required"), Some(&Json::Bool(row.required)));
        }
    }
}

#[test]
fn tools_listing_keywords_equal_the_python_stub() {
    use bammetrics_core::registry::python_keywords;
    let doc = tools_listing();
    for (entry, verb) in array(&doc, "tools")
        .iter()
        .zip(bammetrics_cli::verbs::VERBS)
    {
        let options = array(entry, "options");
        let keywords = python_keywords(verb.tool());
        for keyword in keywords.iter().filter(|k| !k.output) {
            let found = options
                .iter()
                .find(|o| first_name(o) == keyword.picard)
                .unwrap_or_else(|| panic!("{} lacks {}", verb.name, keyword.picard));
            assert_eq!(first_name(found).to_ascii_lowercase(), keyword.name);
        }
        for option in options.iter().filter(|o| text(o, "kind") == "value") {
            let lower = first_name(option).to_ascii_lowercase();
            let positional = lower == "input" || option.get("required") == Some(&Json::Bool(true));
            assert!(
                positional || keywords.iter().any(|k| k.name == lower),
                "{}: {lower} is no keyword",
                verb.name
            );
        }
    }
}

#[test]
fn listing_marks_aliases_ignored_and_rejected_options() {
    let doc = tools_listing();
    let tools = array(&doc, "tools");
    let tool = |verb: &str| tools.iter().find(|t| text(t, "verb") == verb).expect(verb);
    let option = |verb: &str, name: &str| {
        array(tool(verb), "options")
            .iter()
            .find(|o| first_name(o) == name)
            .unwrap_or_else(|| panic!("{verb} {name}"))
    };
    let histogram = option("insert-size", "HISTOGRAM_FILE");
    assert_eq!(text(histogram, "kind"), "ignored");
    assert_eq!(array(histogram, "names").len(), 2);
    assert_eq!(text(option("wgs", "ASSUME_SORTED"), "kind"), "rejected");
    assert_eq!(text(tool("insert-size"), "chain_role"), "tap");
    assert_eq!(text(tool("mark-duplicates"), "chain_role"), "transform");
    assert_eq!(text(tool("split-sam-by-library"), "chain_role"), "split");
    assert_eq!(text(tool("build-bam-index"), "chain_role"), "none");
    assert_eq!(
        tool("jumping-library").get("reads_stdin"),
        Some(&Json::Bool(false))
    );
    let codes: Vec<&str> = array(&doc, "exit_codes")
        .iter()
        .map(|c| match c.get("code") {
            Some(Json::Number(n)) => n.as_str(),
            _ => panic!("code"),
        })
        .collect();
    assert_eq!(
        codes,
        ["0", "1", "2", "10", "11", "100", "129", "130", "143"]
    );
}

/// Every flag and short of one tool entry, options and verb options.
fn spellings(tool: &Json) -> Vec<String> {
    let mut all = Vec::new();
    for key in ["options", "verb_options"] {
        for option in array(tool, key) {
            for member in ["flag", "short"] {
                if let Some(spelling) = option.get(member).and_then(Json::as_str) {
                    all.push(spelling.to_owned());
                }
            }
        }
    }
    all
}

#[test]
fn no_two_options_of_a_tool_share_a_flag_or_short() {
    let doc = tools_listing();
    for tool in array(&doc, "tools") {
        let mut seen = std::collections::HashSet::new();
        for spelling in spellings(tool) {
            assert!(
                seen.insert(spelling.clone()),
                "{}: {spelling} is listed twice",
                text(tool, "verb")
            );
        }
    }
}

#[test]
fn rejected_options_have_no_flag_and_verb_options_are_listed() {
    let doc = tools_listing();
    let tools = array(&doc, "tools");
    let tool = |verb: &str| tools.iter().find(|t| text(t, "verb") == verb).expect(verb);
    for t in tools {
        for option in array(t, "options") {
            if text(option, "kind") == "rejected" {
                assert_eq!(option.get("flag"), Some(&Json::Null), "{}", text(t, "verb"));
                assert_eq!(
                    option.get("short"),
                    Some(&Json::Null),
                    "{}",
                    text(t, "verb")
                );
            }
        }
    }
    let flags = |verb: &str| -> Vec<&str> {
        array(tool(verb), "verb_options")
            .iter()
            .map(|o| text(o, "flag"))
            .collect()
    };
    for verb in ["view-sam", "bam-index-stats", "split-vcfs"] {
        assert!(flags(verb).contains(&"--output"), "{verb}");
    }
    assert!(flags("insert-size").contains(&"--threads"));
    // `-q` works on the verb; the `QUIET` row it does not set has no flag.
    assert!(flags("insert-size").contains(&"--quiet"));
    let quiet = array(tool("insert-size"), "options")
        .iter()
        .find(|o| first_name(o) == "QUIET")
        .expect("row");
    assert_eq!(quiet.get("flag"), Some(&Json::Null));
    let extensions: Vec<&str> = array(tool("multiple-metrics"), "extensions")
        .iter()
        .filter_map(|e| e.as_str())
        .collect();
    assert!(extensions.contains(&".bait_bias_detail_metrics"));
    let width = array(tool("insert-size"), "options")
        .iter()
        .find(|o| first_name(o) == "MINIMUM_HISTOGRAM_WIDTH")
        .expect("row");
    assert!(
        array(width, "names")
            .iter()
            .any(|n| n.as_str() == Some("MIN_HISTOGRAM_WIDTH"))
    );
}

#[test]
fn recipes_member_lists_the_compiled_in_recipes() {
    let doc = tools_listing();
    let names: Vec<&str> = array(&doc, "recipes")
        .iter()
        .map(|r| text(r, "name"))
        .collect();
    assert_eq!(names, ["prep", "qc-only", "downsample"]);
    for recipe in array(&doc, "recipes") {
        assert!(!text(recipe, "description").is_empty());
    }
}

/// Every file a golden record of a prefix verb writes is the `-o` prefix
/// plus a listed suffix (and, for a `--format` list, a format extension).
#[test]
fn prefix_verbs_write_the_listed_suffixes() {
    let doc = tools_listing();
    let prefix: Vec<(&str, Vec<&str>)> = array(&doc, "tools")
        .iter()
        .filter(|t| text(t, "files") == "prefix")
        .map(|t| {
            let suffixes = array(t, "extensions")
                .iter()
                .map(|e| e.as_str().expect("suffix"))
                .collect();
            (text(t, "verb"), suffixes)
        })
        .collect();
    let mut checked = std::collections::HashSet::new();
    for entry in fs::read_dir(root().join("tests/golden")).expect("golden") {
        let path = entry.expect("entry").path();
        let record = fs::read_to_string(&path).expect("record");
        let mut lines = record.lines();
        if lines.next() != Some("## argv") {
            continue;
        }
        let argv: Vec<&str> = lines
            .by_ref()
            .take_while(|l| !l.starts_with("## "))
            .collect();
        let Some((_, suffixes)) = prefix.iter().find(|(verb, _)| Some(verb) == argv.first()) else {
            continue;
        };
        let stem = argv
            .iter()
            .position(|a| *a == "-o")
            .and_then(|at| argv.get(at + 1))
            .and_then(|value| value.strip_prefix("OUT/"))
            .filter(|stem| !stem.is_empty() && !stem.ends_with('/'));
        let Some(stem) = stem else { continue };
        for name in record.lines().filter_map(|l| l.strip_prefix("## file ")) {
            let rest = name
                .strip_prefix(stem)
                .unwrap_or_else(|| panic!("{}: {name} lacks the prefix {stem}", path.display()));
            assert!(
                suffixes
                    .iter()
                    .any(|s| rest == *s || rest.starts_with(&format!("{s}."))),
                "{}: {name} is {stem} plus none of {suffixes:?}",
                path.display()
            );
            checked.insert(argv[0].to_owned());
        }
    }
    // Golden records cover these prefix verbs; the rest have none that writes.
    for verb in [
        "gc-bias",
        "multiple-metrics",
        "sequencing-artifact",
        "convert-sequencing-artifact-to-oxo-g",
        "rrbs",
        "create-baf-regress-file",
        "create-verify-i-d-intensity-contamination-file",
        "arrays-variant-calling",
        "split-vcfs",
    ] {
        assert!(checked.contains(verb), "no golden record checks {verb}");
    }
}

fn parsed(out: &Output) -> Json {
    Json::parse(std::str::from_utf8(&out.stdout).expect("utf-8")).expect("JSON")
}

#[test]
fn version_json_equals_the_listing_member() {
    let member = tools_listing().get("bammetrics").expect("member").clone();
    for args in [["--version", "--json"], ["--json", "--version"]] {
        let out = bammetrics(&args);
        assert_eq!(out.status.code(), Some(0), "{args:?}");
        let doc = parsed(&out);
        assert_eq!(doc, member, "{args:?}");
        assert_eq!(text(&doc, "version"), env!("CARGO_PKG_VERSION"));
    }
}

#[test]
fn plain_version_is_unchanged() {
    let out = bammetrics(&["--version"]);
    assert_eq!(out.status.code(), Some(0));
    assert_eq!(
        String::from_utf8(out.stdout).expect("utf-8"),
        format!("bammetrics {}\n", env!("CARGO_PKG_VERSION"))
    );
}

#[test]
fn json_flag_alone_or_with_more_is_a_usage_error() {
    for args in [
        &["--json"][..],
        &["--version", "--json", "extra"],
        &["--json", "--version", "extra"],
        &["--json", "tools"],
    ] {
        let out = bammetrics(args);
        assert_eq!(out.status.code(), Some(2), "{args:?}");
        assert!(out.stdout.is_empty(), "{args:?}");
        assert!(!out.stderr.is_empty(), "{args:?}");
    }
}

#[test]
fn short_version_json_equals_the_long_form_and_help_json_takes_one_topic() {
    let long = bammetrics(&["--version", "--json"]);
    for args in [&["-V", "--json"][..], &["--json", "-V"]] {
        let out = bammetrics(args);
        assert_eq!(out.status.code(), Some(0), "{args:?}");
        assert_eq!(out.stdout, long.stdout, "{args:?}");
    }
    let out = bammetrics(&["help", "insert-size", "--format", "json", "extra"]);
    assert_eq!(out.status.code(), Some(2));
    assert!(out.stdout.is_empty());
}

#[test]
fn tools_and_help_declare_format_and_the_listing_omits_help() {
    let doc = tools_listing();
    let verbs = array(&doc, "verbs");
    for name in ["tools", "help"] {
        let verb = verbs.iter().find(|v| text(v, "name") == name).expect(name);
        let flags: Vec<&str> = array(verb, "options")
            .iter()
            .map(|o| text(o, "flag"))
            .collect();
        assert_eq!(flags, ["--format"], "{name}");
    }
    for name in ["collect", "chain"] {
        let verb = verbs.iter().find(|v| text(v, "name") == name).expect(name);
        assert!(
            array(verb, "options")
                .iter()
                .all(|o| text(o, "flag") != "--help"),
            "{name}"
        );
    }
    let out = bammetrics(&["help", "tools"]);
    assert!(String::from_utf8_lossy(&out.stdout).contains("--format"));
}

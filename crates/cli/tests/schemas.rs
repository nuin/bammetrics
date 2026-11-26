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
    use bammetrics_core::registry::{COMMON, Kind, Ty};
    let doc = tools_listing();
    let tools = array(&doc, "tools");
    assert_eq!(tools.len(), 77);
    assert_eq!(tools.len(), bammetrics_cli::verbs::VERBS.len());
    for (entry, verb) in tools.iter().zip(bammetrics_cli::verbs::VERBS) {
        assert_eq!(text(entry, "verb"), verb.name);
        assert_eq!(text(entry, "class"), verb.class);
        let tool = verb.tool();
        let own_output = tool.options.iter().any(|o| o.kind == Kind::Output);
        let common = COMMON.iter().filter(|c| {
            !(own_output && c.kind == Kind::Output)
                && !tool.options.iter().any(|o| o.names.contains(&c.picard()))
        });
        let rows: Vec<_> = tool.options.iter().chain(common).collect();
        let options = array(entry, "options");
        assert_eq!(options.len(), rows.len(), "{}", verb.name);
        for (option, row) in options.iter().zip(rows) {
            let names: Vec<&str> = array(option, "names")
                .iter()
                .map(|n| n.as_str().expect("name"))
                .collect();
            assert_eq!(names, row.names, "{}", verb.name);
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
    assert_eq!(codes, ["0", "1", "2", "10", "11", "100", "130"]);
}

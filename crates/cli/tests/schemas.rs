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
fn versions() -> [(&'static str, &'static str); 5] {
    [
        ("metrics-1.json", "bammetrics.metrics/1"),
        ("verdict-1.json", bammetrics_verdict::verdict::SCHEMA),
        ("manifest-1.json", bammetrics_cli::manifest::SCHEMA),
        ("manifest-2.json", bammetrics_cli::manifest::v2::SCHEMA),
        ("events-1.json", "bammetrics events, version 1"),
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

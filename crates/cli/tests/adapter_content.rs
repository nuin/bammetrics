//! `adapter-content` end to end: the lab preset of docs/06 section 7.6
//! judges the file through `check` and drives `qc`, which then runs the
//! tool as a source. The byte-exact forms are in the golden suite; the
//! planted truth and the FastQC comparison are
//! `tests/python/test_adapter_content.py`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "tests/data/adapter_content.bam";

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

fn scratch(test: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-adapter-{test}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch");
    dir
}

/// The TOML block between the preset markers of docs/06.
fn documented_preset() -> String {
    let doc = fs::read_to_string(root().join("docs/06_clinical_checks.md")).expect("doc");
    let start = doc
        .find("<!-- adapter preset start -->")
        .expect("start marker");
    let end = doc.find("<!-- adapter preset end -->").expect("end marker");
    let body = doc[start..end]
        .split_once("```toml\n")
        .expect("toml fence")
        .1;
    body.split_once("```").expect("closing fence").0.to_owned()
}

fn utf8(path: &Path) -> &str {
    path.to_str().expect("utf-8")
}

#[test]
fn the_documented_preset_fails_on_the_dimers() {
    let dir = scratch("check");
    let run = dir.join("run");
    fs::create_dir_all(&run).expect("run dir");
    let metrics = run.join("s1.adapter_content_metrics");
    let out = bammetrics(&["adapter-content", BAM, "-o", utf8(&metrics)]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let preset = dir.join("adapter.toml");
    fs::write(&preset, documented_preset()).expect("preset");
    let out = bammetrics(&["check", utf8(&run), "--preset", utf8(&preset)]);
    let text = String::from_utf8_lossy(&out.stdout);
    // The fixture's unpaired reads are half Universal dimers: 50 percent.
    assert_eq!(out.status.code(), Some(10), "{text}");
    assert!(text.contains("preset adapter-content v1"), "{text}");
    assert!(
        text.contains("adapter.max-pct  50  warn > 5, fail > 10"),
        "{text}"
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn qc_runs_the_tool_for_a_preset_that_names_the_class() {
    let dir = scratch("qc");
    let preset = dir.join("adapter.toml");
    fs::write(&preset, documented_preset()).expect("preset");
    let prefix = dir.join("s1");
    let out = bammetrics(&[
        "qc",
        BAM,
        "-o",
        utf8(&prefix),
        "--preset",
        utf8(&preset),
        "--fail-on",
        "fail",
        "-q",
    ]);
    assert_eq!(
        out.status.code(),
        Some(10),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let metrics = fs::read_to_string(dir.join("s1.adapter_content_metrics")).expect("metrics");
    let direct = bammetrics(&["adapter-content", BAM]);
    let body = |text: &str| text.split_once("\n\n").map(|(_, b)| b.to_owned());
    assert_eq!(
        body(&metrics),
        body(&String::from_utf8_lossy(&direct.stdout)),
        "qc writes what the verb writes"
    );
    let _ = fs::remove_dir_all(&dir);
}

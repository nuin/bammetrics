//! `report` and `check` end to end on `tests/data/picard`, a directory
//! Picard 3.5.0 wrote (CollectMultipleMetrics, CollectHsMetrics and
//! CollectDuplicateMetrics on the fixture BAM, plus the per-base coverage
//! file as a decoy). The exit codes, the files written, the verdict.json
//! round trip, a two-sample tree, and a file that looks like metrics but
//! does not parse. The byte-exact text and jsonl forms of `check` are in
//! the golden suite.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_verdict::{Band, Verdict};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

const PICARD: &str = "tests/data/picard";

fn bammetrics(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn code(output: &Output) -> i32 {
    output.status.code().expect("an exit code")
}

fn stdout(output: &Output) -> String {
    String::from_utf8_lossy(&output.stdout).into_owned()
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-verbs-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    dir
}

#[test]
fn report_on_a_picard_directory_writes_three_files_and_exits_zero() {
    let out = fresh_dir("report");
    let run = bammetrics(&[
        "report",
        PICARD,
        "--preset",
        "clinical-panel",
        "-o",
        out.to_str().expect("utf-8"),
        "--date",
        "2026-01-02 03:04:05 UTC",
    ]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    let text = stdout(&run);
    assert!(text.starts_with("1 sample, preset clinical-panel v1 (3d0bcd), bammetrics "));
    assert!(text.contains("FAIL        NA12878  5 fail, 8 pass, 1 not measured\n"));
    assert!(
        text.contains("  FAIL  Mean target coverage is 0.377778, below the fail limit of 100\n")
    );
    assert!(text.contains("  +2 more, see NA12878.report.html\n"));
    assert!(text.contains("wrote 3 files to "));
    assert!(text.ends_with("bammetrics check gates on it\n"));

    let mut names: Vec<String> = fs::read_dir(&out)
        .expect("output directory")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    assert_eq!(
        names,
        ["NA12878.report.html", "NA12878.verdict.json", "index.html"]
    );

    let verdict = Verdict::from_json(
        &fs::read_to_string(out.join("NA12878.verdict.json")).expect("verdict.json"),
    )
    .expect("parses");
    assert_eq!(verdict.sample, "NA12878");
    assert_eq!(verdict.status, Band::Fail);
    assert_eq!(verdict.counts.fail, 5);
    assert_eq!(
        verdict.files.len(),
        10,
        "every Picard metrics file, the per-base file left out"
    );
    assert!(
        verdict
            .files
            .iter()
            .all(|f| !f.path.ends_with(".per_base_coverage"))
    );
    assert!(
        verdict
            .files
            .iter()
            .any(|f| f.path.ends_with("NA12878.quality_by_cycle_metrics"))
    );

    let page = fs::read_to_string(out.join("NA12878.report.html")).expect("page");
    assert!(page.contains("<title>NA12878: FAIL (bammetrics)</title>"));
    assert!(page.contains("2026-01-02 03:04:05 UTC"));
    assert!(page.contains("id=\"plot-coverage\""));
    let index = fs::read_to_string(out.join("index.html")).expect("index");
    assert!(index.contains("<title>1 sample: FAIL (bammetrics)</title>"));
    assert!(index.contains("<a href=\"NA12878.report.html\">NA12878</a>"));
    assert!(index.contains("2026-01-02 03:04:05 UTC"));

    // check reads the verdict.json back and gates on it.
    let json_path = out.join("NA12878.verdict.json");
    let gate = bammetrics(&["check", json_path.to_str().expect("utf-8")]);
    assert_eq!(code(&gate), 10);
    assert!(stdout(&gate).ends_with("exit 10: a sample reaches --fail-on fail\n"));

    // A second run overwrites in place.
    let again = bammetrics(&["report", PICARD, "-o", out.to_str().expect("utf-8")]);
    assert_eq!(code(&again), 0);
    assert!(stdout(&again).contains("PASS        NA12878  6 pass\n"));
    assert!(!out.join("NA12878.report.html.tmp").exists());
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn report_defaults_to_the_directory_given_and_makes_the_output_directory() {
    let base = fresh_dir("tree");
    let a = base.join("runs/A");
    let b = base.join("runs/B");
    fs::create_dir_all(&a).expect("mkdir");
    fs::create_dir_all(&b).expect("mkdir");
    // Files named by tool, not by sample, so each directory names its sample.
    for (name, copy) in [
        ("NA12878.insert_size_metrics", "insert_size.txt"),
        ("NA12878.quality_yield_metrics", "quality_yield.txt"),
    ] {
        let text =
            fs::read_to_string(root().join(PICARD).join("NA12878").join(name)).expect("fixture");
        fs::write(a.join(copy), &text).expect("copy");
        // B lacks the insert size file, so it is incomplete under generic.
        if copy != "insert_size.txt" {
            fs::write(b.join(copy), &text).expect("copy");
        }
    }
    // A file that looks like metrics and does not parse is skipped.
    fs::write(
        a.join("broken.txt"),
        "## htsjdk.samtools.metrics.StringHeader\n# x\n\n## METRICS CLASS\tX\n",
    )
    .expect("write");
    let runs = base.join("runs");
    let run = bammetrics(&["report", runs.to_str().expect("utf-8")]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert!(stderr(&run).contains("bammetrics: skipped "));
    assert!(stderr(&run).contains("broken.txt"));
    let text = stdout(&run);
    assert!(text.starts_with("2 samples, preset generic v1"));
    let index = fs::read_to_string(runs.join("index.html")).expect("index in DIR");
    assert!(runs.join("A.report.html").exists());
    assert!(runs.join("B.verdict.json").exists());
    let b_at = index
        .find("<a href=\"B.report.html\">B</a>")
        .expect("B row");
    let a_at = index
        .find("<a href=\"A.report.html\">A</a>")
        .expect("A row");
    assert!(a_at < b_at, "both incomplete with no fail, so by name");

    let nested = base.join("site/deeper");
    let run = bammetrics(&[
        "report",
        runs.to_str().expect("utf-8"),
        "-o",
        nested.to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert!(nested.join("index.html").exists());

    let gate = bammetrics(&["check", runs.to_str().expect("utf-8")]);
    assert_eq!(code(&gate), 11, "both samples are incomplete under generic");
    let _ = fs::remove_dir_all(&base);
}

#[test]
fn check_formats_and_exit_codes() {
    let json = bammetrics(&[
        "check",
        PICARD,
        "--format",
        "json",
        "--preset",
        "clinical-panel",
    ]);
    assert_eq!(code(&json), 10);
    let text = stdout(&json);
    assert!(text.starts_with("[{\n  \"schema\": \"bammetrics.verdict/1\""));
    assert!(text.trim_end().ends_with("}]"));

    let jsonl = bammetrics(&["check", PICARD, "--format", "jsonl"]);
    assert_eq!(code(&jsonl), 0);
    let lines: Vec<String> = stdout(&jsonl).lines().map(str::to_owned).collect();
    assert_eq!(lines.len(), 7, "6 checks and 1 verdict");
    assert!(lines[0].starts_with(
        "{\"event\":\"check\",\"sample\":\"NA12878\",\"rule\":\"alignment.pf-reads-aligned\""
    ));
    assert!(
        lines[6].starts_with("{\"event\":\"verdict\",\"sample\":\"NA12878\",\"status\":\"pass\"")
    );

    let warn = bammetrics(&["check", PICARD, "--fail-on", "warn"]);
    assert_eq!(code(&warn), 0, "the fixture passes generic outright");
    let wgs = bammetrics(&["check", PICARD, "--preset", "wgs"]);
    assert_eq!(code(&wgs), 11);
    assert!(stdout(&wgs).ends_with("exit 11: a sample is incomplete\n"));
    let usage = bammetrics(&["report", PICARD, PICARD]);
    assert_eq!(code(&usage), 2);
    let missing = bammetrics(&["check", "tests/data/nope"]);
    assert_eq!(code(&missing), 1);
    assert!(stderr(&missing).contains("tests/data/nope"));
    let no_preset = bammetrics(&["check", PICARD, "--preset", "nope"]);
    assert_eq!(code(&no_preset), 1);
}

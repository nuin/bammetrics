//! `IntervalListTools` beside the rest of the command line: `qc --panel`
//! takes the lists it writes (docs/07 W12 acceptance), and a run that fails
//! leaves no scatter directory. Concordance with Picard is
//! tests/concordance/test_interval_list_tools.py; the command forms are in
//! the golden suite.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "tests/data/insert_size.bam";
const TARGETS: &str = "tests/data/targets.interval_list";
const LIST_A: &str = "tests/data/intervals_a.interval_list";
const REF: &str = "tests/data/ref_100k.fa";

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
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-interval-list-tools-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

fn ok(output: &Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// The hs metrics `qc --panel PANEL` writes, without its `StringHeader`
/// lines (they name the panel).
fn hs_metrics(dir: &Path, panel: &Path) -> String {
    let prefix = dir.join("s");
    let run = bammetrics(&[
        "qc",
        BAM,
        "-o",
        &prefix.display().to_string(),
        "--panel",
        &panel.display().to_string(),
        "-R",
        REF,
        "--preset",
        "clinical-panel",
        "--date",
        "2026-01-02 03:04:05 UTC",
    ]);
    ok(&run);
    let text = fs::read_to_string(dir.join("s.hs_metrics")).expect("hs metrics");
    text.lines()
        .filter(|line| !line.starts_with("# "))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn qc_panel_takes_the_lists_it_writes() {
    let dir = fresh_dir("panel");
    // The file stem names the bait set, so the union keeps the stem.
    let union = dir.join("union/targets.interval_list");
    ok(&bammetrics(&[
        "interval-list-tools",
        TARGETS,
        "--action",
        "UNION",
        "--comment",
        "panel",
        "-o",
        &union.display().to_string(),
    ]));
    let written = fs::read_to_string(&union).expect("union");
    assert!(written.contains("@PG\tID:1\tCL:IntervalListTools INPUT=["));
    assert!(written.contains("@CO\tpanel\n"));
    let original = hs_metrics(&dir.join("original"), &root().join(TARGETS));
    assert_eq!(hs_metrics(&dir.join("union"), &union), original);

    // One list of a scatter is a panel too.
    let scatter = dir.join("scatter");
    ok(&bammetrics(&[
        "IntervalListTools",
        &format!("I={TARGETS}"),
        &format!("I={LIST_A}"),
        "SCATTER_COUNT=2",
        &format!("O={}", scatter.display()),
    ]));
    let first = scatter.join("temp_0001_of_2/scattered.interval_list");
    assert!(!hs_metrics(&dir.join("first"), &first).is_empty());
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_failed_scatter_leaves_no_directory() {
    let dir = fresh_dir("failed");
    let scatter = dir.join("scatter");
    let run = bammetrics(&[
        "IntervalListTools",
        &format!("I={LIST_A}"),
        "SCATTER_COUNT=3",
        &format!("O={}", scatter.display()),
        "OUTPUT_VALUE=BASES",
        &format!("COUNT_OUTPUT={}", dir.join("missing/n.txt").display()),
    ]);
    assert_eq!(run.status.code(), Some(1));
    assert!(!scatter.exists(), "the scatter directory was left");
    assert!(!dir.join("missing").exists());
    // The verb makes the directory of --count-output, as it makes -o's.
    let made = bammetrics(&[
        "interval-list-tools",
        LIST_A,
        "--scatter-count",
        "3",
        "-o",
        &scatter.display().to_string(),
        "--output-value",
        "BASES",
        "--count-output",
        &dir.join("missing/n.txt").display().to_string(),
    ]);
    ok(&made);
    assert_eq!(
        fs::read_to_string(dir.join("missing/n.txt")).expect("count"),
        "1555\n"
    );
    assert_eq!(fs::read_dir(&scatter).expect("scatter").count(), 3);
    let _ = fs::remove_dir_all(&dir);
}

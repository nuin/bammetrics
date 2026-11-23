//! Wave W14 (docs/chaining.md sections 22 to 24): the chain's manifest,
//! `--resume`, `--keep` and `qc` as the last step.
//!
//! The gate of docs/07: in a five-step chain, change one option of step 3
//! and only steps 3, 4 and 5 run; the taps in front of step 3 keep their
//! files untouched, and every file equals what a fresh run of the changed
//! chain writes. A kept stream is the BAM `write` writes at that position.
//! The qc step writes the verdict and the report the `qc` verb writes for
//! the same taps on the same records.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::SystemTime;

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;

const MARKDUP: &str = "tests/data/markdup.bam";
const INSERT_SIZE: &str = "tests/data/insert_size.bam";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn run(args: &[String]) -> Output {
    let output = bammetrics(args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-resume-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

/// A metrics file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.split_inclusive('\n')
        .skip_while(|line| line.starts_with('#'))
        .collect()
}

/// The header lines of a BAM with the `CL` field of each `@PG` line
/// removed, and every record as its bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<RawRecord>) {
    pass::with_record_stream(path, None, 1, |header, reader, _| {
        let text = bammetrics_core::bam_writer::header_bytes(header);
        let lines = String::from_utf8_lossy(&text)
            .lines()
            .map(|line| {
                line.split('\t')
                    .filter(|field| !(line.starts_with("@PG") && field.starts_with("CL:")))
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect();
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader)? {
            records.push(record.clone());
        }
        Ok((lines, records))
    })
    .expect("a BAM")
}

/// A report with each metrics file's command line taken out: the chain's
/// taps name the chain's source and `CHAIN=`, the verb's name its BAM.
fn without_command_lines(text: &str) -> String {
    text.lines()
        .filter(|line| !line.contains(" INPUT="))
        .map(|line| {
            let mut kept = String::new();
            let mut rest = line;
            let open = "<span class=\"muted mono small\">";
            while let Some(start) = rest.find(open) {
                kept.push_str(&rest[..start + open.len()]);
                rest = &rest[start + open.len()..];
                rest = &rest[rest.find('<').unwrap_or(rest.len())..];
            }
            kept.push_str(rest);
            kept
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn modified(path: &Path) -> SystemTime {
    fs::metadata(path)
        .and_then(|m| m.modified())
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The five steps of the gate, step 3 with `records` as its record bound.
fn gate_chain(prefix: &Path, records: &str, keep: bool, extra: &[&str]) -> Vec<String> {
    let mut args: Vec<String> = ["chain", MARKDUP, "-o"]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
    args.push(prefix.display().to_string());
    args.extend(extra.iter().map(|a| (*a).to_owned()));
    let mut steps = vec!["+", "quality-yield", "+", "mark-duplicates"];
    if keep {
        steps.push("--keep");
    }
    steps.extend([
        "--program-group-command-line",
        "pinned",
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "--max-records-in-ram",
        records,
        "+",
        "insert-size",
        "+",
        "write",
    ]);
    args.extend(steps.iter().map(|a| (*a).to_owned()));
    args
}

/// The files of the gate chain under `prefix`, by step.
fn gate_files(prefix: &Path) -> [(usize, PathBuf); 5] {
    let file = |suffix: &str| PathBuf::from(format!("{}{suffix}", prefix.display()));
    [
        (1, file(".in.quality_yield_metrics")),
        (2, file(".marked.duplication_metrics")),
        (2, file(".marked.bam")),
        (4, file(".insert_size_metrics")),
        (5, file(".sorted.bam")),
    ]
}

/// The gate: change one option of step 3 and only steps 3, 4 and 5 run.
/// Steps 1 and 2 keep their files untouched, the read starts at the kept
/// stream of step 2, and every file equals a fresh run of the new chain.
#[test]
fn an_option_of_step_3_reruns_steps_3_4_and_5_only() {
    let dir = fresh_dir("gate");
    let prefix = dir.join("s");
    run(&gate_chain(&prefix, "5", true, &["--resume"]));
    let before: Vec<(usize, SystemTime)> = gate_files(&prefix)
        .iter()
        .map(|(step, path)| (*step, modified(path)))
        .collect();
    // The same chain skips every step and reads nothing.
    let again = run(&gate_chain(&prefix, "5", true, &["--resume"]));
    let stderr = String::from_utf8_lossy(&again.stderr);
    assert!(
        stderr.contains("5 steps finished, 5 of them skipped by --resume"),
        "{stderr}"
    );
    for ((_, path), (_, time)) in gate_files(&prefix).iter().zip(&before) {
        assert_eq!(modified(path), *time, "{}", path.display());
    }
    // A changed option of step 3.
    std::thread::sleep(std::time::Duration::from_millis(20));
    let changed = run(&gate_chain(&prefix, "7", true, &["--resume"]));
    let stderr = String::from_utf8_lossy(&changed.stderr);
    assert!(
        stderr.contains("step 3 (sort-sam): 36 records in"),
        "{stderr}"
    );
    assert!(stderr.contains("2 of them skipped by --resume"), "{stderr}");
    for ((step, path), (_, time)) in gate_files(&prefix).iter().zip(&before) {
        let rewritten = modified(path) != *time;
        assert_eq!(rewritten, *step >= 3, "step {step}: {}", path.display());
    }
    // A fresh run of the changed chain writes the same files.
    let fresh = dir.join("fresh/s");
    run(&gate_chain(&fresh, "7", true, &[]));
    for ((_, ours), (_, theirs)) in gate_files(&prefix).iter().zip(gate_files(&fresh)) {
        if ours.extension().is_some_and(|e| e == "bam") {
            assert_eq!(bam(ours), bam(&theirs), "{}", ours.display());
        } else {
            assert_eq!(body(ours), body(&theirs), "{}", ours.display());
        }
    }
    let manifest = fs::read_to_string(dir.join("s.manifest.json")).expect("the manifest");
    assert!(manifest.contains("\"schema\": \"bammetrics.manifest/2\""));
    assert_eq!(manifest.matches("\"skipped\": true").count(), 2);
    let _ = fs::remove_dir_all(&dir);
}

/// Without `--keep` the read starts at the source: the transform in front
/// of the changed step runs again to feed it, and the tap on the source
/// stays skipped.
#[test]
fn without_a_kept_stream_the_transforms_in_front_run_to_feed() {
    let dir = fresh_dir("nokeep");
    let prefix = dir.join("s");
    run(&gate_chain(&prefix, "5", false, &["--resume"]));
    let tap = PathBuf::from(format!("{}.in.quality_yield_metrics", prefix.display()));
    let marked = PathBuf::from(format!("{}.marked.duplication_metrics", prefix.display()));
    let (tap_time, marked_time) = (modified(&tap), modified(&marked));
    std::thread::sleep(std::time::Duration::from_millis(20));
    let plan = run(&gate_chain(&prefix, "7", false, &["--resume", "--dry-run"]));
    let text = String::from_utf8_lossy(&plan.stdout);
    assert!(
        text.contains("resume     step 2 mark-duplicates: runs: its stream feeds step 3"),
        "{text}"
    );
    let changed = run(&gate_chain(&prefix, "7", false, &["--resume"]));
    let stderr = String::from_utf8_lossy(&changed.stderr);
    assert!(
        stderr.contains("step 1 (quality-yield): skipped"),
        "{stderr}"
    );
    assert_eq!(modified(&tap), tap_time);
    assert_ne!(modified(&marked), marked_time);
    let _ = fs::remove_dir_all(&dir);
}

/// The kept stream of a transform is the BAM `write` writes right after
/// it, byte for byte, with its index for a coordinate-sorted stream.
#[test]
fn a_kept_stream_is_the_bam_write_writes_there() {
    let dir = fresh_dir("keep");
    let prefix = dir.join("s");
    let written = dir.join("written.bam");
    let args: Vec<String> = [
        "chain",
        MARKDUP,
        "-o",
        &prefix.display().to_string(),
        "+",
        "mark-duplicates",
        "--keep",
        "+",
        "write",
        &written.display().to_string(),
        "--index",
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    run(&args);
    let kept = dir.join("s.marked.bam");
    assert_eq!(
        fs::read(&kept).expect("kept"),
        fs::read(&written).expect("written")
    );
    assert_eq!(
        fs::read(dir.join("s.marked.bai")).expect("kept index"),
        fs::read(dir.join("written.bai")).expect("written index")
    );
    let manifest = fs::read_to_string(dir.join("s.manifest.json")).expect("the manifest");
    assert!(manifest.contains("\"stream\": \"marked\""), "{manifest}");
    let _ = fs::remove_dir_all(&dir);
}

/// The qc step after a transform writes the verdict and the report the
/// `qc` verb writes for the BAM the chain wrote at that position.
#[test]
fn the_qc_step_equals_the_qc_verb_on_the_same_records() {
    let dir = fresh_dir("qc");
    let date = "2026-01-02 03:04:05 UTC";
    let marked = dir.join("marked.bam");
    let chain: Vec<String> = [
        "chain",
        MARKDUP,
        "-o",
        &dir.join("a/s").display().to_string(),
        "+",
        "mark-duplicates",
        "+",
        "write",
        &marked.display().to_string(),
        "+",
        "qc",
        "--date",
        date,
        "--also",
        "junit,markdown",
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    run(&chain);
    let qc: Vec<String> = [
        "qc",
        &marked.display().to_string(),
        "-o",
        &dir.join("b/s").display().to_string(),
        "--date",
        date,
        "--also",
        "junit,markdown",
    ]
    .iter()
    .map(|a| (*a).to_owned())
    .collect();
    run(&qc);
    let (a, b) = (
        dir.join("a").display().to_string(),
        dir.join("b").display().to_string(),
    );
    for name in [
        "s.verdict.json",
        "s.report.html",
        "s.junit.xml",
        "s.report.md",
    ] {
        let ours = fs::read_to_string(dir.join("a").join(name)).expect("chain file");
        let theirs = fs::read_to_string(dir.join("b").join(name)).expect("qc file");
        assert_eq!(
            without_command_lines(&ours.replace(&a, &b)),
            without_command_lines(&theirs),
            "{name}"
        );
    }
    for name in [
        "s.alignment_summary_metrics",
        "s.insert_size_metrics",
        "s.duplicate_metrics",
    ] {
        assert_eq!(
            body(&dir.join("a").join(name)),
            body(&dir.join("b").join(name)),
            "{name}"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

/// A source that changed reruns every step; a manifest that is no
/// manifest is a note, never a failure; a file the run before did not
/// write, or that changed since, is never removed.
#[test]
fn what_resume_trusts_and_what_it_leaves_alone() {
    let dir = fresh_dir("trust");
    let source = dir.join("in.bam");
    fs::copy(root().join(INSERT_SIZE), &source).expect("copy the source");
    let prefix = dir.join("s");
    let args = |extra: &[&str]| -> Vec<String> {
        let mut args: Vec<String> = vec![
            "chain".to_owned(),
            source.display().to_string(),
            "-o".to_owned(),
            prefix.display().to_string(),
            "--resume".to_owned(),
        ];
        args.extend(extra.iter().map(|a| (*a).to_owned()));
        args
    };
    run(&args(&["+", "insert-size", "+", "quality-yield"]));
    // The source is touched: every step runs.
    let bytes = fs::read(&source).expect("read the source");
    std::thread::sleep(std::time::Duration::from_millis(20));
    fs::write(&source, &bytes).expect("rewrite the source");
    let rerun = run(&args(&["+", "insert-size", "+", "quality-yield"]));
    let stderr = String::from_utf8_lossy(&rerun.stderr);
    assert!(stderr.contains("changed: every step runs"), "{stderr}");
    // A damaged manifest: every step runs, no failure.
    let manifest = dir.join("s.manifest.json");
    fs::write(&manifest, "[".repeat(10_000)).expect("damage the manifest");
    let damaged = run(&args(&["+", "insert-size", "+", "quality-yield"]));
    let stderr = String::from_utf8_lossy(&damaged.stderr);
    assert!(stderr.contains("not JSON"), "{stderr}");
    // quality-yield leaves the chain; its file was edited since, so it is
    // no longer the file the manifest lists, and it stays.
    let quality = dir.join("s.quality_yield_metrics");
    fs::write(&quality, "edited by hand\n").expect("edit");
    run(&args(&["+", "insert-size"]));
    assert_eq!(
        fs::read_to_string(&quality).expect("kept"),
        "edited by hand\n"
    );
    let _ = fs::remove_dir_all(&dir);
}

/// A step that failed in the run before reruns; one that finished stays.
#[test]
fn a_step_that_failed_reruns() {
    let dir = fresh_dir("failed");
    let prefix = dir.join("s");
    let args = |level: &str| -> Vec<String> {
        [
            "chain",
            "tests/data/duplicate_marking.bam",
            "-o",
            &prefix.display().to_string(),
            "--resume",
            "+",
            "quality-yield",
            level,
            "+",
            "insert-size",
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect()
    };
    let first = bammetrics(&args("--no-assume-sorted"));
    assert_eq!(first.status.code(), Some(1));
    let again = bammetrics(&args("--no-assume-sorted"));
    let stderr = String::from_utf8_lossy(&again.stderr);
    assert!(stderr.contains("step 2 (insert-size): skipped"), "{stderr}");
    assert!(
        stderr.contains("step 1 (quality-yield): failed alone"),
        "{stderr}"
    );
    let _ = fs::remove_dir_all(&dir);
}

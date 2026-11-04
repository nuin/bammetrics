//! `qc` and `plan` end to end on the committed fixture BAM: the files are
//! the Picard face's bytes, `--resume` after a killed run and after one
//! lost file or an interval file edited in place, per-collector isolation,
//! the plan writes nothing, the exit codes with and without `--fail-on`,
//! errors as jsonl events, a closed stdout, and the verdict and page equal
//! what `report` writes for the same files. The byte-exact stdout of each form
//! is in the golden suite.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_verdict::json::Json;

const BAM: &str = "tests/data/insert_size.bam";
const PANEL: &str = "tests/data/targets.interval_list";
const REF: &str = "tests/data/ref_100k.fa";
const DATE: &str = "2026-01-02 03:04:05 UTC";

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
    let dir = std::env::temp_dir().join(format!("bammetrics-qc-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    dir
}

/// `qc BAM -o DIR/s --date DATE` plus `extra`.
fn qc(dir: &Path, extra: &[&str]) -> Output {
    let prefix = dir.join("s");
    let mut args = vec![
        "qc",
        BAM,
        "-o",
        prefix.to_str().expect("utf-8"),
        "--date",
        DATE,
    ];
    args.extend_from_slice(extra);
    bammetrics(&args)
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("directory")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

fn manifest(dir: &Path) -> Json {
    Json::parse(&fs::read_to_string(dir.join("s.manifest.json")).expect("manifest"))
        .expect("manifest parses")
}

fn steps(dir: &Path) -> Vec<Json> {
    manifest(dir)
        .get("steps")
        .and_then(Json::as_array)
        .expect("steps")
        .to_vec()
}

fn text<'a>(json: &'a Json, key: &str) -> &'a str {
    json.get(key).and_then(Json::as_str).expect(key)
}

/// The steps a `--format jsonl` run wrote metrics files for.
fn written_steps(output: &Output) -> Vec<String> {
    stdout(output)
        .lines()
        .map(|l| Json::parse(l).expect("a JSON line"))
        .filter(|e| {
            text(e, "event") == "file"
                && text(e, "kind") == "metrics"
                && e.get("written") == Some(&Json::Bool(true))
        })
        .map(|e| text(&e, "step").to_owned())
        .collect()
}

/// The text without its `# Started on:` line.
fn without_clock(text: &str) -> String {
    text.lines()
        .filter(|l| !l.starts_with("# Started on: "))
        .collect::<Vec<_>>()
        .join("\n")
}

const GENERIC_FILES: [&str; 7] = [
    "s.alignment_summary_metrics",
    "s.duplicate_metrics",
    "s.insert_size_metrics",
    "s.manifest.json",
    "s.quality_yield_metrics",
    "s.report.html",
    "s.verdict.json",
];

#[test]
fn every_metrics_file_is_the_picard_face_output_for_its_recorded_argv() {
    let dir = fresh_dir("picard-face");
    let run = qc(
        &dir,
        &["--panel", PANEL, "-R", REF, "--preset", "clinical-panel"],
    );
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    let all = steps(&dir);
    let collectors: Vec<&Json> = all.iter().filter(|s| text(s, "stage") == "pass").collect();
    assert_eq!(collectors.len(), 4);
    for step in collectors {
        let tool = text(step, "tool");
        let ours: Vec<String> = step
            .get("outputs")
            .and_then(Json::as_array)
            .expect("outputs")
            .iter()
            .map(|o| fs::read_to_string(text(o, "path")).expect("qc file"))
            .collect();
        let face = dir.join("face");
        let _ = fs::remove_dir_all(&face);
        fs::create_dir_all(&face).expect("face dir");
        let argv: Vec<String> = step
            .get("argv")
            .and_then(Json::as_array)
            .expect("argv")
            .iter()
            .map(|a| a.as_str().expect("text").to_owned())
            .collect();
        // Rerun the recorded command with the files moved to `face`.
        let prefix = dir.join("s").display().to_string();
        let moved: Vec<String> = argv
            .iter()
            .map(|a| a.replace(&prefix, &face.join("s").display().to_string()))
            .collect();
        let mut args: Vec<&str> = vec![tool];
        args.extend(moved.iter().map(String::as_str));
        let picard = bammetrics(&args);
        assert_eq!(code(&picard), 0, "{tool}: {}", stderr(&picard));
        assert_eq!(ours.len(), names(&face).len(), "{tool}");
        for (ours, name) in ours.iter().zip(names(&face)) {
            let theirs = fs::read_to_string(face.join(&name)).expect("face file");
            assert_eq!(
                without_clock(ours),
                without_clock(&theirs)
                    .replace(&face.display().to_string(), &dir.display().to_string()),
                "{tool} {name}"
            );
        }
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_killed_run_resumes_from_nothing_and_a_lost_file_reruns_its_collector_only() {
    let dir = fresh_dir("resume");
    let first = qc(&dir, &[]);
    assert_eq!(code(&first), 0, "{}", stderr(&first));
    assert_eq!(names(&dir), GENERIC_FILES);

    // Without --resume or --force an earlier run is not replaced.
    let again = qc(&dir, &[]);
    assert_eq!(code(&again), 2);
    assert!(stderr(&again).contains("--resume skips its finished steps"));

    // Every collector matches: no BAM read, every metrics file kept.
    let kept = qc(&dir, &["--resume", "--format", "jsonl"]);
    assert_eq!(code(&kept), 0, "{}", stderr(&kept));
    assert!(written_steps(&kept).is_empty(), "{}", stdout(&kept));

    // One lost file reruns that collector alone.
    fs::remove_file(dir.join("s.insert_size_metrics")).expect("remove");
    let one = qc(&dir, &["--resume", "--format", "jsonl"]);
    assert_eq!(code(&one), 0, "{}", stderr(&one));
    assert_eq!(written_steps(&one), ["CollectInsertSizeMetrics"]);
    assert_eq!(names(&dir), GENERIC_FILES);

    // A killed run: files gone and the manifest cut short. Every step reruns.
    fs::remove_file(dir.join("s.quality_yield_metrics")).expect("remove");
    fs::remove_file(dir.join("s.report.html")).expect("remove");
    let manifest_text = fs::read_to_string(dir.join("s.manifest.json")).expect("manifest");
    fs::write(
        dir.join("s.manifest.json"),
        &manifest_text[..manifest_text.len() / 3],
    )
    .expect("truncate");
    let killed = qc(&dir, &["--resume", "--format", "jsonl"]);
    assert_eq!(code(&killed), 0, "{}", stderr(&killed));
    assert!(stderr(&killed).contains("is not a whole manifest; running every step"));
    assert_eq!(written_steps(&killed).len(), 4);
    assert_eq!(names(&dir), GENERIC_FILES);
    assert!(steps(&dir).iter().all(|s| text(s, "status") == "ok"));

    // --force reruns everything over an earlier run.
    let forced = qc(&dir, &["--force", "--format", "jsonl"]);
    assert_eq!(code(&forced), 0, "{}", stderr(&forced));
    assert_eq!(written_steps(&forced).len(), 4);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_changed_preset_or_input_reruns_the_collectors() {
    let dir = fresh_dir("keys");
    assert_eq!(code(&qc(&dir, &[])), 0);
    let lab = qc(
        &dir,
        &[
            "--resume",
            "--format",
            "jsonl",
            "--preset",
            "tests/data/presets/lab.toml",
        ],
    );
    assert_eq!(code(&lab), 0, "{}", stderr(&lab));
    assert_eq!(written_steps(&lab), ["CollectInsertSizeMetrics"]);
    // The lab preset reads insert size and hs, and hs has no --panel: the
    // manifest lists one collector.
    let collectors = steps(&dir)
        .iter()
        .filter(|s| text(s, "stage") == "pass")
        .count();
    assert_eq!(collectors, 1);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_failing_collector_fails_alone_and_the_run_exits_1() {
    let dir = fresh_dir("isolation");
    fs::create_dir_all(&dir).expect("dir");
    let bad = dir.join("bad.bed");
    fs::write(&bad, "chrNope\t10\t100\n").expect("bed");
    let run = qc(
        &dir,
        &[
            "--panel",
            bad.to_str().expect("utf-8"),
            "--preset",
            "clinical-panel",
        ],
    );
    assert_eq!(code(&run), 1, "{}", stdout(&run));
    assert!(stderr(&run).contains("CollectHsMetrics failed: "));
    let summary = stdout(&run);
    assert!(
        summary.contains("ERROR         CollectHsMetrics"),
        "{summary}"
    );
    assert!(summary.contains("(9 rules not measured)"), "{summary}");
    assert!(
        summary.ends_with("exit 1: 1 collector(s) failed\n"),
        "{summary}"
    );
    assert!(summary.lines().count() <= 12, "{summary}");
    for name in [
        "s.alignment_summary_metrics",
        "s.insert_size_metrics",
        "s.quality_yield_metrics",
        "s.verdict.json",
        "s.report.html",
    ] {
        assert!(dir.join(name).is_file(), "{name}");
    }
    assert!(!dir.join("s.hs_metrics").exists());
    let hs = steps(&dir)
        .into_iter()
        .find(|s| text(s, "id") == "CollectHsMetrics")
        .expect("hs step");
    assert_eq!(text(&hs, "status"), "failed");
    assert!(text(&hs, "error").contains("chrNope"));
    // The failed step is not resumed: it runs again, alone.
    let again = qc(
        &dir,
        &[
            "--panel",
            bad.to_str().expect("utf-8"),
            "--preset",
            "clinical-panel",
            "--resume",
            "--format",
            "jsonl",
        ],
    );
    assert_eq!(code(&again), 1);
    assert!(stdout(&again).contains("\"event\":\"error\",\"step\":\"CollectHsMetrics\""));
    assert!(written_steps(&again).is_empty());
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn plan_runs_nothing_and_names_what_resume_would_skip() {
    let dir = fresh_dir("plan");
    let prefix = dir.join("s");
    let prefix = prefix.to_str().expect("utf-8");
    let plan = bammetrics(&["plan", BAM, "-o", prefix, "-R", REF, "--preset", "wgs"]);
    assert_eq!(code(&plan), 0, "{}", stderr(&plan));
    assert!(!dir.exists(), "plan wrote nothing");
    let text = stdout(&plan);
    assert!(
        text.contains("pass      4 collectors, 1 shared BAM read\n"),
        "{text}"
    );
    assert!(
        text.contains("  reads     2 BAM reads: 1 for the pass, 1 own read\n"),
        "{text}"
    );
    assert!(text.contains("own read  run CollectWgsMetrics"), "{text}");
    assert!(
        text.contains("not in the single pass, it walks reference loci"),
        "{text}"
    );
    assert!(text.ends_with("nothing ran: this is a plan\n"));

    assert_eq!(code(&qc(&dir, &["-R", REF, "--preset", "wgs"])), 0);
    fs::remove_file(dir.join("s.duplicate_metrics")).expect("remove");
    let plan = bammetrics(&[
        "plan", BAM, "-o", prefix, "-R", REF, "--preset", "wgs", "--resume",
    ]);
    let text = stdout(&plan);
    assert!(text.contains("--resume skips 4 of 5 collectors"), "{text}");
    assert!(
        text.contains("pass      4 collectors, 3 skipped by --resume, 1 shared BAM read\n"),
        "{text}"
    );
    assert!(
        text.contains("  reads     1 BAM read: 1 for the pass, 0 own reads\n"),
        "{text}"
    );
    assert!(text.contains("  run   CollectDuplicateMetrics"), "{text}");
    assert!(text.contains("  skip  CollectInsertSizeMetrics"), "{text}");
    let lines = bammetrics(&[
        "plan", BAM, "-o", prefix, "-R", REF, "--preset", "wgs", "--format", "jsonl",
    ]);
    let events: Vec<Json> = stdout(&lines)
        .lines()
        .map(|l| Json::parse(l).expect("JSON line"))
        .collect();
    let wgs = events
        .iter()
        .find(|e| e.get("id").and_then(Json::as_str) == Some("CollectWgsMetrics"))
        .expect("wgs step");
    assert_eq!(wgs.get("read").and_then(Json::as_str), Some("own"));
    assert_eq!(wgs.get("resumable"), Some(&Json::Bool(true)));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn exit_codes_follow_fail_on_and_qc_alone_never_gates() {
    let cases: [(&[&str], i32); 8] = [
        (&[], 0),
        (&["--fail-on", "warn"], 0),
        (&["--preset", "clinical-panel"], 0),
        (&["--preset", "clinical-panel", "--fail-on", "fail"], 11),
        (&["--preset", "clinical-panel", "--panel", PANEL], 0),
        (
            &[
                "--preset",
                "clinical-panel",
                "--panel",
                PANEL,
                "--fail-on",
                "fail",
            ],
            10,
        ),
        (&["--preset", "nope"], 1),
        (&["--resume", "--force"], 2),
    ];
    for (i, (args, expected)) in cases.iter().enumerate() {
        let dir = fresh_dir(&format!("exit-{i}"));
        let run = qc(&dir, args);
        assert_eq!(code(&run), *expected, "{args:?}: {}", stderr(&run));
        let _ = fs::remove_dir_all(&dir);
    }
    let missing = bammetrics(&["qc", "tests/data/nope.bam", "-o", "/nonexistent/x"]);
    assert_eq!(code(&missing), 1);
    assert!(stderr(&missing).contains("tests/data/nope.bam not found (from BAM)"));
}

#[test]
fn the_verdict_and_page_equal_what_report_writes_for_the_same_files() {
    let dir = fresh_dir("same-as-report");
    let run = qc(
        &dir,
        &[
            "--sample",
            "s",
            "--panel",
            PANEL,
            "--preset",
            "clinical-panel",
        ],
    );
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    let out = fresh_dir("same-as-report-out");
    let report = bammetrics(&[
        "report",
        dir.to_str().expect("utf-8"),
        "--preset",
        "clinical-panel",
        "-o",
        out.to_str().expect("utf-8"),
        "--date",
        DATE,
    ]);
    assert_eq!(code(&report), 0, "{}", stderr(&report));
    for name in ["s.verdict.json", "s.report.html"] {
        assert_eq!(
            fs::read_to_string(dir.join(name)).expect("qc"),
            fs::read_to_string(out.join(name)).expect("report"),
            "{name}"
        );
    }
    let _ = fs::remove_dir_all(&dir);
    let _ = fs::remove_dir_all(&out);
}

#[test]
fn without_o_the_prefix_is_the_bam_name_in_the_current_directory() {
    let dir = fresh_dir("default-prefix");
    fs::create_dir_all(&dir).expect("dir");
    let bam = root().join(BAM);
    let run = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["qc", bam.to_str().expect("utf-8"), "--date", DATE])
        .current_dir(&dir)
        .output()
        .expect("run");
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert!(dir.join("insert_size.insert_size_metrics").is_file());
    assert!(
        stdout(&run).contains("to ./insert_size.* (no -o: named after the BAM)"),
        "{}",
        stdout(&run)
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn an_interval_file_edited_in_place_reruns_hs_under_resume() {
    let dir = fresh_dir("aux-edit");
    fs::create_dir_all(&dir).expect("dir");
    let panel = dir.join("p.interval_list");
    fs::copy(root().join(PANEL), &panel).expect("copy panel");
    let panel = panel.to_str().expect("utf-8");
    let args = [
        "--panel",
        panel,
        "--preset",
        "clinical-panel",
        "--format",
        "jsonl",
    ];
    assert_eq!(code(&qc(&dir, &args)), 0);
    let hs = steps(&dir)
        .into_iter()
        .find(|s| text(s, "id") == "CollectHsMetrics")
        .expect("hs step");
    let recorded: Vec<&str> = hs
        .get("files")
        .and_then(Json::as_array)
        .expect("files")
        .iter()
        .map(|f| text(f, "path"))
        .collect();
    assert_eq!(recorded, [panel], "the panel is fingerprinted once");

    let mut resume = args.to_vec();
    resume.push("--resume");
    assert!(written_steps(&qc(&dir, &resume)).is_empty());
    let mut text_of = fs::read_to_string(panel).expect("panel");
    let last = text_of.lines().last().expect("an interval").to_owned();
    text_of.push_str(&last);
    text_of.push('\n');
    fs::write(panel, text_of).expect("edit panel");
    assert_eq!(written_steps(&qc(&dir, &resume)), ["CollectHsMetrics"]);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_collector_the_new_options_drop_loses_its_old_file() {
    let dir = fresh_dir("dropped");
    assert_eq!(code(&qc(&dir, &["--panel", PANEL])), 0);
    assert!(dir.join("s.hs_metrics").is_file());
    let forced = qc(&dir, &["--force"]);
    assert_eq!(code(&forced), 0, "{}", stderr(&forced));
    assert!(
        stdout(&forced).contains("1 stale removed"),
        "{}",
        stdout(&forced)
    );
    assert_eq!(names(&dir), GENERIC_FILES);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn an_output_directory_gets_the_bam_name_as_the_file_part() {
    for (i, trailing) in ["/", ""].into_iter().enumerate() {
        let dir = fresh_dir(&format!("dir-output-{i}"));
        if trailing.is_empty() {
            fs::create_dir_all(&dir).expect("existing directory");
        }
        let output = format!("{}{trailing}", dir.display());
        let run = bammetrics(&["qc", BAM, "-o", &output, "--date", DATE]);
        assert_eq!(code(&run), 0, "{}", stderr(&run));
        assert!(
            stdout(&run).contains("(-o is a directory: named after the BAM)"),
            "{}",
            stdout(&run)
        );
        let expected: Vec<String> = GENERIC_FILES
            .iter()
            .map(|name| name.replacen("s.", "insert_size.", 1))
            .collect();
        assert_eq!(names(&dir), expected);
        let _ = fs::remove_dir_all(&dir);
    }
}

#[test]
fn a_prepare_error_is_a_jsonl_error_and_done_event() {
    for verb in ["qc", "plan"] {
        let run = bammetrics(&[verb, "tests/data/nope.bam", "--format", "jsonl"]);
        assert_eq!(code(&run), 1);
        assert_eq!(
            stderr(&run),
            format!("bammetrics {verb}: tests/data/nope.bam not found (from BAM)\n")
        );
        let events: Vec<Json> = stdout(&run)
            .lines()
            .map(|l| Json::parse(l).expect("a JSON line"))
            .collect();
        assert_eq!(events.len(), 2, "{}", stdout(&run));
        assert_eq!(text(&events[0], "event"), "error");
        assert_eq!(text(&events[1], "event"), "done");
        assert_eq!(events[1].get("exit"), Some(&Json::Number(1.0)));
    }
}

#[test]
fn a_closed_stdout_does_not_stop_the_run() {
    let dir = fresh_dir("closed-stdout");
    let prefix = dir.join("s");
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "qc",
            BAM,
            "-o",
            prefix.to_str().expect("utf-8"),
            "--format",
            "jsonl",
        ])
        .current_dir(root())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect("spawn");
    drop(child.stdout.take());
    let run = child.wait_with_output().expect("wait");
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert!(!stderr(&run).contains("panicked"), "{}", stderr(&run));
    assert_eq!(names(&dir), GENERIC_FILES);
    let _ = fs::remove_dir_all(&dir);
}

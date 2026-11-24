//! `downsample-sam` and `position-based-downsample-sam` (wave W13): the
//! `downsample-sam` step of a chain equals the staged tools at the first, a
//! middle and the last position beside `mark-duplicates` and `sort-sam`,
//! for each strategy; its `--metrics-file` in a chain equals the tool's;
//! `position-based-downsample-sam` reads its option in three spellings as
//! the Picard face does, refuses standard input and is refused as a step.
//!
//! Staged means each tool from file to file, then the taps on the last
//! file; chained means one `chain` command. The records of the written BAM
//! as bytes, its header lines (the `CL` of `@PG` set aside) and each
//! metrics file after its `StringHeader` lines must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const T1: &str = "tests/corpus/out/t1.bam";

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

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_owned()).collect()
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-downsample-{}-{name}",
        std::process::id()
    ));
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
    pass::with_record_stream_at(path, None, 1, Stringency::Silent, |header, reader, _| {
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

/// One step of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `downsample-sam` with this strategy at probability 0.6 and accuracy
    /// 0.05 (chunks of 20 names).
    Downsample(&'static str),
    Sort(&'static str),
    MarkDuplicates,
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
    match step {
        Step::Downsample(strategy) => {
            let strategy = format!("STRATEGY={strategy}");
            let metrics = dir.join(format!("staged.{n}.quality_yield_metrics"));
            let chained = dir.join(format!("chain.{n}.quality_yield_metrics"));
            (
                strings(&[
                    "DownsampleSam",
                    &strategy,
                    "P=0.6",
                    "A=0.05",
                    &format!("M={}", metrics.display()),
                ]),
                strings(&[
                    "downsample-sam",
                    &strategy,
                    "--probability",
                    "0.6",
                    "--accuracy",
                    "0.05",
                    "--metrics-file",
                    &chained.display().to_string(),
                    "--as",
                    &format!("d{n}"),
                ]),
            )
        }
        Step::Sort(order) => {
            let option = format!("SORT_ORDER={order}");
            (
                strings(&["SortSam", &option]),
                strings(&["sort-sam", &option]),
            )
        }
        Step::MarkDuplicates => {
            let metrics = dir.join(format!("staged.{n}.duplication_metrics"));
            (
                strings(&["MarkDuplicates", &format!("M={}", metrics.display())]),
                strings(&["mark-duplicates", "--as", &format!("m{n}")]),
            )
        }
    }
}

/// Run `steps` staged and chained on `source` and compare.
fn chained_equals_staged(name: &str, source: &str, steps: &[Step]) {
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus)");
        return;
    }
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        run(&staged);
        input = out;
    }
    run(&[
        "CollectQualityYieldMetrics".to_owned(),
        format!("I={}", input.display()),
        format!("O={}", dir.join("staged.quality_yield_metrics").display()),
    ]);
    let written = dir.join("chained.bam");
    chain.extend(strings(&["+", "quality-yield", "+", "write"]));
    chain.push(written.display().to_string());
    run(&chain);
    assert_eq!(
        body(&dir.join("chain.quality_yield_metrics")),
        body(&dir.join("staged.quality_yield_metrics")),
        "{name}: the tap"
    );
    for (n, step) in steps.iter().enumerate() {
        let (staged, chained) = match step {
            Step::Downsample(_) => (
                format!("staged.{n}.quality_yield_metrics"),
                format!("chain.{n}.quality_yield_metrics"),
            ),
            Step::MarkDuplicates => (
                format!("staged.{n}.duplication_metrics"),
                format!("chain.m{n}.duplication_metrics"),
            ),
            Step::Sort(_) => continue,
        };
        assert_eq!(
            body(&dir.join(chained)),
            body(&dir.join(staged)),
            "{name}: step {n}"
        );
    }
    let (chain_header, chain_records) = bam(&written);
    let (staged_header, staged_records) = bam(&input);
    assert_eq!(chain_header, staged_header, "{name}: header");
    assert_eq!(chain_records.len(), staged_records.len(), "{name}");
    assert!(chain_records == staged_records, "{name}: records differ");
    assert!(!chain_records.is_empty(), "{name}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{Downsample, MarkDuplicates as Md, Sort};

/// Each strategy first, in the middle and last, beside `mark-duplicates`
/// and `sort-sam`.
#[test]
fn a_downsample_step_at_any_position_equals_the_staged_tools() {
    for strategy in ["ConstantMemory", "HighAccuracy", "Chained"] {
        let cases: [(&str, Vec<Step>); 4] = [
            ("first", vec![Downsample(strategy), Md]),
            ("middle", vec![Sort("queryname"), Downsample(strategy), Md]),
            ("last", vec![Md, Sort("queryname"), Downsample(strategy)]),
            ("alone", vec![Downsample(strategy)]),
        ];
        for (position, steps) in cases {
            chained_equals_staged(&format!("{strategy}-{position}"), MARKDUP, &steps);
        }
    }
}

#[test]
fn the_corpus_downsamples_as_the_staged_tools() {
    for strategy in ["ConstantMemory", "HighAccuracy", "Chained"] {
        chained_equals_staged(&format!("t1-{strategy}"), T1, &[Downsample(strategy), Md]);
    }
}

/// `--fraction` in its three spellings gives the records the Picard face
/// writes.
#[test]
fn position_based_reads_three_spellings_as_the_picard_face() {
    let dir = fresh_dir("spellings");
    let picard = dir.join("picard.bam");
    run(&strings(&[
        "PositionBasedDownsampleSam",
        &format!("I={MARKDUP}"),
        &format!("O={}", picard.display()),
        "FRACTION=0.4",
    ]));
    let (picard_header, picard_records) = bam(&picard);
    assert!(!picard_records.is_empty());
    for (i, form) in [["--fraction", "0.4"], ["--FRACTION", "0.4"], ["-F", "0.4"]]
        .iter()
        .enumerate()
    {
        let out = dir.join(format!("verb{i}.bam"));
        let mut args = strings(&["position-based-downsample-sam", MARKDUP]);
        args.extend(strings(form));
        args.extend(["-o".to_owned(), out.display().to_string()]);
        run(&args);
        let (header, records) = bam(&out);
        assert_eq!(header, picard_header, "{form:?}");
        assert!(records == picard_records, "{form:?}");
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn position_based_refuses_standard_input() {
    let dir = fresh_dir("stdin");
    let out = dir.join("out.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "position-based-downsample-sam",
            "-",
            "--fraction",
            "0.5",
            "-o",
            &out.display().to_string(),
        ])
        .current_dir(root())
        .stdin(Stdio::from(
            fs::File::open(root().join(MARKDUP)).expect("fixture"),
        ))
        .output()
        .expect("run bammetrics");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("reads its input twice"), "{stderr}");
    assert!(!out.exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// `position-based-downsample-sam` reads its input twice: in a chain it
/// is a materialized step (wave W14) that runs on a file of its stream,
/// here behind `mark-duplicates`, and its BAM equals the staged tool's.
#[test]
fn position_based_is_a_materialized_step() {
    let dir = fresh_dir("step");
    let path = |name: &str| dir.join(name).display().to_string();
    run(&strings(&[
        "chain",
        MARKDUP,
        "-o",
        &path("c"),
        "+",
        "mark-duplicates",
        "+",
        "position-based-downsample-sam",
        "--fraction",
        "0.5",
        "+",
        "write",
        &path("chained.bam"),
    ]));
    run(&strings(&[
        "chain",
        MARKDUP,
        "-o",
        &path("s"),
        "+",
        "mark-duplicates",
        "+",
        "write",
        &path("marked.bam"),
    ]));
    run(&strings(&[
        "position-based-downsample-sam",
        &path("marked.bam"),
        "-o",
        &path("staged.bam"),
        "--fraction",
        "0.5",
    ]));
    assert_eq!(bam(&dir.join("chained.bam")), bam(&dir.join("staged.bam")));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A step's options are checked before the first read, as the tool checks
/// them.
#[test]
fn a_downsample_step_refuses_what_the_tool_refuses() {
    let dir = fresh_dir("refuse");
    for (args, says) in [
        (vec!["--probability", "1.5"], "PROBABILITY"),
        (vec!["--strategy", "Bogus"], "STRATEGY"),
        (vec!["--strategy", "Chained", "--accuracy", "1"], "accuracy"),
        (vec!["--compression-level", "1"], "COMPRESSION_LEVEL"),
    ] {
        let mut argv = strings(&[
            "chain",
            MARKDUP,
            "-o",
            &dir.join("c").display().to_string(),
            "+",
            "downsample-sam",
        ]);
        argv.extend(strings(&args));
        argv.extend(strings(&["+", "write"]));
        argv.push(dir.join("out.bam").display().to_string());
        let output = bammetrics(&argv);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{args:?}: {stderr}");
        assert!(stderr.contains(says), "{args:?}: {stderr}");
        assert!(!dir.join("out.bam").exists());
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// `OUTPUT` is renamed into place before the registry writes
/// `METRICS_FILE`: a metrics directory the run cannot write fails the run
/// before the read, so no BAM is left.
#[test]
#[cfg(unix)]
fn an_unwritable_metrics_directory_leaves_no_bam() {
    use std::os::unix::fs::PermissionsExt;
    let dir = fresh_dir("metrics-dir");
    let locked = dir.join("ro");
    fs::create_dir(&locked).expect("dir");
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o555)).expect("chmod");
    // A user who writes anyway (root) cannot see the refusal.
    if fs::write(locked.join("probe"), b"").is_ok() {
        fs::set_permissions(&locked, fs::Permissions::from_mode(0o755)).expect("chmod");
        fs::remove_dir_all(&dir).expect("cleanup");
        return;
    }
    let path = |name: &str| dir.join(name).display().to_string();
    for (argv, out) in [
        (
            vec![
                "DownsampleSam".to_owned(),
                format!("I={MARKDUP}"),
                format!("O={}", path("d.bam")),
                format!("M={}", path("ro/d.txt")),
                "P=0.5".to_owned(),
            ],
            "d.bam",
        ),
        (
            vec![
                "MarkDuplicates".to_owned(),
                format!("I={MARKDUP}"),
                format!("O={}", path("m.bam")),
                format!("M={}", path("ro/m.txt")),
            ],
            "m.bam",
        ),
    ] {
        let output = bammetrics(&argv);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(1), "{argv:?}: {stderr}");
        assert!(stderr.contains("ro/"), "{argv:?}: {stderr}");
        assert!(!dir.join(out).exists(), "{argv:?} left {out}");
    }
    assert_eq!(fs::read_dir(&locked).expect("list").count(), 0);
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o755)).expect("chmod");
    fs::remove_dir_all(&dir).expect("cleanup");
}

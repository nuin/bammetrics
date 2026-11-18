//! `add-oa-tag` in a chain (docs/chaining.md, "Record-editing transforms"):
//! chained equals staged at the first, a middle and the last position,
//! beside `mark-duplicates` and `sort-sam`, with and without
//! `--interval-list`, on BAM and SAM text sources.
//!
//! Staged runs each transform as its own tool from file to file, then the
//! taps on the last file; chained runs one `chain` command. The header
//! lines of the written BAM (the `CL` of `@PG` set aside), each record as
//! bytes, each metrics file after its `StringHeader` lines and each
//! `MarkDuplicates` metrics file must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const MARKDUP_SAM: &str = "tests/data/markdup.sam";
const TARGETS: &str = "tests/data/targets.interval_list";
const T1: &str = "tests/corpus/out/t1.bam";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";

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
    let dir = std::env::temp_dir().join(format!("bammetrics-oa-{}-{name}", std::process::id()));
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

/// The header lines of a BAM without the `CL` field of `@PG`, and every
/// record as its bytes.
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `add-oa-tag`, with `--interval-list` when true.
    Oa(bool),
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
    match step {
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
        Step::Oa(intervals) => {
            let mut staged = strings(&["AddOATag"]);
            let mut chained = strings(&["add-oa-tag", "--as", &format!("o{n}")]);
            if intervals {
                staged.push(format!("INTERVAL_LIST={TARGETS}"));
                chained.extend(strings(&["-L", TARGETS]));
            }
            (staged, chained)
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
    for (class, extension) in [
        ("CollectQualityYieldMetrics", ".quality_yield_metrics"),
        ("CollectInsertSizeMetrics", ".insert_size_metrics"),
    ] {
        run(&[
            class.to_owned(),
            format!("I={}", input.display()),
            format!("O={}", dir.join(format!("staged{extension}")).display()),
        ]);
    }
    let written = dir.join("chained.bam");
    chain.extend(strings(&[
        "+",
        "quality-yield",
        "+",
        "insert-size",
        "+",
        "write",
    ]));
    chain.push(written.display().to_string());
    run(&chain);
    for extension in [".quality_yield_metrics", ".insert_size_metrics"] {
        let staged = dir.join(format!("staged{extension}"));
        let chained = dir.join(format!("chain{extension}"));
        if staged.exists() || chained.exists() {
            assert_eq!(body(&chained), body(&staged), "{name}: {extension}");
        }
    }
    for (n, step) in steps.iter().enumerate() {
        if matches!(step, Step::MarkDuplicates) {
            assert_eq!(
                body(&dir.join(format!("chain.m{n}.duplication_metrics"))),
                body(&dir.join(format!("staged.{n}.duplication_metrics"))),
                "{name}: mark-duplicates {n}"
            );
        }
    }
    let (chain_header, chain_records) = bam(&written);
    let (staged_header, staged_records) = bam(&input);
    assert_eq!(chain_header, staged_header, "{name}: header");
    assert_eq!(chain_records.len(), staged_records.len(), "{name}");
    assert!(chain_records == staged_records, "{name}: records differ");
    assert!(!chain_records.is_empty(), "{name}");
    let tagged = chain_records
        .iter()
        .filter(|record| record.tag(*b"OA").expect("tags").is_some())
        .count();
    assert!(tagged > 0, "{name}: no record has an OA tag");
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{MarkDuplicates as Md, Oa, Sort};

/// The step first, in the middle and last, beside `mark-duplicates` and
/// `sort-sam`, with and without intervals, and twice in one chain (each
/// run appends one more entry).
#[test]
fn the_step_at_any_position_equals_the_staged_tool() {
    let cases: [(&str, &[Step]); 7] = [
        ("first", &[Oa(false), Md]),
        ("middle", &[Sort("queryname"), Oa(true), Md]),
        ("last", &[Md, Sort("queryname"), Oa(false)]),
        ("alone", &[Oa(true)]),
        ("intervals-first", &[Oa(true), Sort("queryname")]),
        ("twice", &[Oa(false), Md, Oa(true)]),
        (
            "between-sorts",
            &[Sort("queryname"), Oa(false), Sort("coordinate")],
        ),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, MARKDUP, steps);
    }
}

#[test]
fn a_sam_source_tags_as_the_staged_tool() {
    chained_equals_staged("sam", MARKDUP_SAM, &[Oa(true), Md]);
}

#[test]
fn the_corpus_tags_as_the_staged_tool() {
    for source in [T1, FOUNDATION] {
        let name = Path::new(source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("corpus");
        chained_equals_staged(
            &format!("{name}-oa"),
            source,
            &[Md, Oa(false), Sort("queryname")],
        );
    }
}

fn chain(args: &[&str]) -> (Option<i32>, String) {
    let mut argv = vec!["chain".to_owned()];
    argv.extend(strings(args));
    let output = bammetrics(&argv);
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// The options of the BAM file belong to the `write` step: a usage error
/// before the first read, with and without `--dry-run`.
#[test]
fn the_step_refuses_what_belongs_to_the_write_step() {
    for (options, message) in [
        (
            &["--create-md5-file"][..],
            "CREATE_MD5_FILE belongs to the write step in a chain",
        ),
        (
            &["-o", "x.bam"][..],
            "names a BAM, and a transform writes none in a chain",
        ),
    ] {
        for dry_run in [false, true] {
            let mut args = vec![MARKDUP];
            if dry_run {
                args.push("--dry-run");
            }
            args.extend(["+", "add-oa-tag"]);
            args.extend(options);
            args.extend(["+", "write", "/dev/null"]);
            let (code, stderr) = chain(&args);
            assert_eq!(code, Some(2), "{options:?} {dry_run}: {stderr}");
            assert!(stderr.contains(message), "{message}\n{stderr}");
        }
    }
}

/// The step counts what it tagged: every record without intervals, the
/// records over them with.
#[test]
fn the_run_reports_the_records_tagged() {
    let dir = fresh_dir("counts");
    let prefix = dir.join("s").display().to_string();
    let (code, stderr) = chain(&[
        MARKDUP,
        "-o",
        &prefix,
        "+",
        "add-oa-tag",
        "+",
        "add-oa-tag",
        "-L",
        TARGETS,
        "+",
        "write",
        "/dev/null",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert!(
        stderr.contains("step 1 (add-oa-tag): 36 records in, 36 out, 36 changed\n"),
        "{stderr}"
    );
    let line = stderr
        .lines()
        .find(|line| line.contains("step 2 (add-oa-tag): 36 records in, 36 out, "))
        .unwrap_or_else(|| panic!("{stderr}"));
    assert!(
        !line.ends_with(" 36 changed") && !line.ends_with(" 0 changed"),
        "{line}"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

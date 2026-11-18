//! `clean-sam` in a chain (docs/chaining.md section 16): chained equals
//! staged with the step at the first, a middle and the last position,
//! beside `mark-duplicates` and `sort-sam`, on `tests/data/clean_sam.bam`
//! (reads off the end of chr1), its SAM text, corpus t1 and the foundation
//! corpus.
//!
//! Staged means each transform as its own tool from file to file, then the
//! taps on the last file. Chained means one `chain` command. The header
//! lines of the written BAM (the `CL` of `@PG` set aside), its records as
//! bytes and each metrics file after its `StringHeader` lines must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const CLEAN: &str = "tests/data/clean_sam.bam";
const CLEAN_SAM: &str = "tests/data/clean_sam.sam";
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
    let dir = std::env::temp_dir().join(format!("bammetrics-clean-{}-{name}", std::process::id()));
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    Clean,
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
        Step::Clean => (
            strings(&["CleanSam"]),
            strings(&["clean-sam", "--as", &format!("c{n}")]),
        ),
    }
}

/// Run `steps` staged and chained on `source` and compare. A SAM text
/// source is read at LENIENT by the chain: `CleanSam` reads it so.
fn chained_equals_staged(name: &str, source: &str, steps: &[Step]) {
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus)");
        return;
    }
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "--validation-stringency", "LENIENT", "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        staged.push("VALIDATION_STRINGENCY=LENIENT".to_owned());
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
        "{name}: quality yield"
    );
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
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{Clean, MarkDuplicates as Md, Sort};

/// The step first, in the middle and last, beside `mark-duplicates` and
/// `sort-sam`, alone and twice.
#[test]
fn clean_sam_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &[Step]); 6] = [
        ("first", &[Clean, Md]),
        ("middle", &[Sort("queryname"), Clean, Md]),
        ("last", &[Md, Sort("queryname"), Clean]),
        ("alone", &[Clean]),
        ("twice", &[Clean, Md, Clean]),
        ("sorted-after", &[Clean, Sort("coordinate"), Md]),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, CLEAN, steps);
    }
}

#[test]
fn a_sam_source_cleans_as_the_staged_tool() {
    chained_equals_staged("sam", CLEAN_SAM, &[Clean, Md]);
    chained_equals_staged("sam-last", CLEAN_SAM, &[Md, Clean]);
}

#[test]
fn the_corpus_cleans_as_the_staged_tools_clean_it() {
    for source in [T1, FOUNDATION] {
        let name = Path::new(source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("corpus");
        chained_equals_staged(&format!("{name}-first"), source, &[Clean, Md]);
        chained_equals_staged(&format!("{name}-last"), source, &[Md, Clean]);
    }
}

/// The step counts every record in and out and the records it changed: on
/// the fixture, eight of ten (five clipped, one `MC`, two MAPQ).
#[test]
fn the_run_reports_the_step_with_its_counts() {
    let dir = fresh_dir("counts");
    let output = run(&strings(&[
        "chain",
        CLEAN,
        "-o",
        &dir.join("s").display().to_string(),
        "+",
        "clean-sam",
        "+",
        "write",
        "/dev/null",
    ]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("step 1 (clean-sam): 10 records in, 10 out, 8 changed\n"),
        "{stderr}"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A read past the end of its reference fails the step and the chain
/// writes nothing.
#[test]
fn a_read_past_the_end_fails_the_step() {
    let dir = fresh_dir("past");
    let out = dir.join("out.bam");
    let output = bammetrics(&strings(&[
        "chain",
        "tests/data/clean_sam_past.sam",
        "--validation-stringency",
        "SILENT",
        "+",
        "clean-sam",
        "+",
        "write",
        &out.display().to_string(),
    ]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("CIGAR 69S"), "{stderr}");
    assert!(!out.exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

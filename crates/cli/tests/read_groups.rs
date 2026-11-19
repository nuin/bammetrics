//! `add-or-replace-read-groups` and
//! `revert-original-base-qualities-and-add-mate-cigar` in a chain (wave
//! W10): chained equals staged with each step first, in the middle and
//! last, beside `mark-duplicates` and `sort-sam`, streaming and sorting, on
//! `tests/data/revert.bam` (pairs without `MC`, `OQ`, an unmapped mate),
//! `tests/data/markdup.bam`, corpus t1 and the foundation corpus.
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

const REVERT: &str = "tests/data/revert.bam";
const MARKDUP: &str = "tests/data/markdup.bam";
const T1: &str = "tests/corpus/out/t1.bam";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";

const RG: [&str; 4] = ["RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1"];
const RG_FLAGS: [&str; 8] = [
    "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1",
];

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
        "bammetrics-read-groups-{}-{name}",
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `add-or-replace-read-groups`, with a sort order or none.
    ReadGroups(Option<&'static str>),
    /// `revert-original-base-qualities-and-add-mate-cigar`, with a sort
    /// order or none.
    Revert(Option<&'static str>),
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
    let order = |o: Option<&str>| o.map(|o| format!("SORT_ORDER={o}"));
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
        Step::ReadGroups(o) => {
            let mut staged = strings(&["AddOrReplaceReadGroups"]);
            staged.extend(strings(&RG));
            let mut chained = strings(&["add-or-replace-read-groups", "--as", &format!("g{n}")]);
            chained.extend(strings(&RG_FLAGS));
            staged.extend(order(o));
            chained.extend(order(o));
            (staged, chained)
        }
        Step::Revert(o) => {
            let mut staged = strings(&[
                "RevertOriginalBaseQualitiesAndAddMateCigar",
                "MAX_RECORDS_TO_EXAMINE=0",
                "CREATE_INDEX=false",
                "CREATE_MD5_FILE=false",
            ]);
            let mut chained = strings(&[
                "revert-original-base-qualities-and-add-mate-cigar",
                "--as",
                &format!("r{n}"),
            ]);
            staged.extend(order(o));
            chained.extend(order(o));
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
    let mut chain = strings(&["chain", source, "--validation-stringency", "SILENT", "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        staged.push("VALIDATION_STRINGENCY=SILENT".to_owned());
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

use Step::{MarkDuplicates as Md, ReadGroups as Rg, Revert, Sort};

#[test]
fn read_groups_at_any_position_equal_the_staged_tools() {
    let cases: [(&str, &[Step]); 6] = [
        ("rg-first", &[Rg(None), Md]),
        ("rg-middle", &[Sort("queryname"), Rg(None), Md]),
        ("rg-last", &[Md, Sort("queryname"), Rg(None)]),
        ("rg-sorts", &[Rg(Some("queryname")), Md]),
        (
            "rg-sorts-back",
            &[Sort("queryname"), Rg(Some("coordinate")), Md],
        ),
        ("rg-unsorted", &[Rg(Some("unsorted"))]),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, MARKDUP, steps);
    }
}

#[test]
fn revert_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &[Step]); 5] = [
        ("revert-first", &[Revert(None), Md]),
        ("revert-middle", &[Sort("queryname"), Revert(None), Md]),
        ("revert-last", &[Rg(None), Revert(Some("queryname"))]),
        ("revert-unsorted", &[Revert(Some("unsorted"))]),
        (
            "revert-twice",
            &[
                Revert(None),
                Rg(Some("queryname")),
                Revert(Some("coordinate")),
            ],
        ),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, REVERT, steps);
    }
}

#[test]
fn the_corpus_equals_the_staged_tools() {
    chained_equals_staged("t1-rg", T1, &[Rg(Some("queryname")), Md]);
    chained_equals_staged("t1-revert", T1, &[Revert(None), Md]);
    chained_equals_staged("foundation-rg", FOUNDATION, &[Md, Rg(None)]);
}

/// t1 and the foundation corpus were sorted by samtools, whose order of
/// records at one position differs from htsjdk's: a `SORT_ORDER` that is
/// the input's streams them unchanged, as the tool's presorted writer does.
#[test]
fn the_inputs_own_order_keeps_samtools_ties() {
    chained_equals_staged("t1-rg-same", T1, &[Rg(Some("coordinate")), Md]);
    chained_equals_staged(
        "foundation-rg-same",
        FOUNDATION,
        &[Rg(Some("coordinate")), Rg(Some("queryname"))],
    );
}

/// The sorting steps count every record in and out.
#[test]
fn the_run_reports_each_step_with_its_counts() {
    let dir = fresh_dir("counts");
    let mut args = strings(&["chain", REVERT, "-o", &dir.join("s").display().to_string()]);
    args.extend(strings(&[
        "+",
        "add-or-replace-read-groups",
        "--sort-order",
        "queryname",
    ]));
    args.extend(strings(&RG_FLAGS));
    args.extend(strings(&[
        "+",
        "revert-original-base-qualities-and-add-mate-cigar",
        "+",
        "write",
        "/dev/null",
    ]));
    let output = run(&args);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("step 1 (add-or-replace-read-groups): 13 records in, 13 out, 13 changed\n"),
        "{stderr}"
    );
    assert!(
        stderr.contains(
            "step 2 (revert-original-base-qualities-and-add-mate-cigar): 13 records in, 13 out"
        ),
        "{stderr}"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A name with two primary first reads fails the reverting step and the
/// chain writes nothing.
#[test]
fn two_first_reads_fail_the_step() {
    let dir = fresh_dir("two-first");
    let out = dir.join("out.bam");
    let output = bammetrics(&strings(&[
        "chain",
        MARKDUP,
        "+",
        "revert-original-base-qualities-and-add-mate-cigar",
        "+",
        "write",
        &out.display().to_string(),
    ]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("first of the pair"), "{stderr}");
    assert!(!out.exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

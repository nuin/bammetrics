//! `bammetrics chain` end to end (docs/chaining.md).
//!
//! Chained equals staged: one process running `mark-duplicates`, the
//! one-pass tools on its records and `write` gives the metrics files and
//! the BAM that `MarkDuplicates` to a file followed by each tool on that
//! file gives. On the committed fixture always; on corpus t1 and the
//! foundation corpus when `tests/corpus/out` holds them (`make corpus`,
//! `make concordance`).
//!
//! Also: a chain without a transform equals `qc`; the plan file equals the
//! command line; a tap's `STOP_AFTER` in the middle of a chain; what a
//! failing step does to the others; every usage error before the first
//! read.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;

const FIXTURE: (&str, &str) = ("tests/data/markdup.bam", "tests/data/ref_100k.fa");
const T1: (&str, &str) = ("tests/corpus/out/t1.bam", "tests/corpus/out/reference.fa");
const FOUNDATION: (&str, &str) = (
    "tests/corpus/out/foundation/corpus.bam",
    "tests/corpus/out/foundation/ref.fa",
);

/// The taps: the verb, its Picard class, the file options with the
/// extension `qc` gives each, and whether it reads the reference.
type Tap = (
    &'static str,
    &'static str,
    &'static [(&'static str, &'static str)],
    bool,
);

const TAPS: [Tap; 8] = [
    (
        "alignment-summary",
        "CollectAlignmentSummaryMetrics",
        &[("O", ".alignment_summary_metrics")],
        true,
    ),
    (
        "insert-size",
        "CollectInsertSizeMetrics",
        &[("O", ".insert_size_metrics")],
        false,
    ),
    (
        "quality-yield",
        "CollectQualityYieldMetrics",
        &[("O", ".quality_yield_metrics")],
        false,
    ),
    (
        "duplicate-metrics",
        "CollectDuplicateMetrics",
        &[("M", ".duplicate_metrics")],
        false,
    ),
    (
        "quality-score-distribution",
        "QualityScoreDistribution",
        &[("O", ".quality_distribution_metrics")],
        false,
    ),
    (
        "mean-quality-by-cycle",
        "MeanQualityByCycle",
        &[("O", ".quality_by_cycle_metrics")],
        false,
    ),
    (
        "base-distribution-by-cycle",
        "CollectBaseDistributionByCycle",
        &[("O", ".base_distribution_by_cycle_metrics")],
        false,
    ),
    (
        "gc-bias",
        "CollectGcBiasMetrics",
        &[
            ("O", ".gc_bias.detail_metrics"),
            ("S", ".gc_bias.summary_metrics"),
        ],
        true,
    ),
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

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_owned()).collect()
}

/// Run and require exit 0.
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
    let dir = std::env::temp_dir().join(format!("bammetrics-chain-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("directory")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

fn arg(key: &str, path: &Path) -> String {
    format!("{key}={}", path.display())
}

fn text(path: &Path) -> String {
    fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// A metrics file without its `StringHeader` lines (the command line and
/// the start time).
fn body(path: &Path) -> String {
    let text = text(path);
    let mut lines = text.split_inclusive('\n').peekable();
    while lines.peek().is_some_and(|line| line.starts_with('#')) {
        lines.next();
    }
    lines.collect()
}

/// A file without the line that carries the clock.
fn without_clock(path: &Path) -> String {
    text(path)
        .lines()
        .filter(|line| !line.contains("Started on: "))
        .collect::<Vec<_>>()
        .join("\n")
}

/// The header lines of a BAM with the `CL` field of each `@PG` line
/// removed, and every record as its bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<RawRecord>) {
    let path = root().join(path);
    pass::with_record_stream(&path, None, 1, |header, reader, _| {
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

/// `MarkDuplicates` to a file and each tap on that file, against one
/// chain under the prefix `chain`. `options` go to `MarkDuplicates` on
/// both sides, in Picard spelling.
fn chained_equals_staged(name: &str, input: (&str, &str), options: &[&str]) {
    let (source, reference) = input;
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus, make concordance)");
        return;
    }
    let dir = fresh_dir(name);
    let options = strings(options);
    let staged_bam = dir.join("staged.marked.bam");
    let mut staged = vec![
        "MarkDuplicates".to_owned(),
        format!("I={source}"),
        arg("O", &staged_bam),
        arg("M", &dir.join("staged.marked.duplication_metrics")),
    ];
    staged.extend(options.iter().cloned());
    run(&staged);
    let prefix = dir.join("chain").display().to_string();
    let mut chain = strings(&["chain", source, "-R", reference, "-o", &prefix]);
    chain.extend(strings(&["+", "mark-duplicates"]));
    chain.extend(options.iter().cloned());
    for (verb, class, files, reads_reference) in TAPS {
        let mut alone = vec![class.to_owned(), arg("I", &staged_bam)];
        for (option, extension) in files {
            alone.push(arg(option, &dir.join(format!("staged{extension}"))));
        }
        if reads_reference {
            alone.push(format!("R={reference}"));
        }
        run(&alone);
        chain.extend(strings(&["+", verb]));
    }
    chain.extend(strings(&["+", "write"]));
    run(&chain);

    let files = names(&dir);
    let staged_files: Vec<&String> = files
        .iter()
        .filter(|f| f.starts_with("staged.") && f.ends_with("_metrics"))
        .collect();
    assert_eq!(staged_files.len(), TAPS.len() + 2, "{files:?}");
    for file in staged_files {
        let chained = file.replacen("staged.", "chain.", 1);
        assert_eq!(
            body(&dir.join(&chained)),
            body(&dir.join(file)),
            "{name}: {chained}"
        );
    }
    let (chain_header, chain_records) = bam(&dir.join("chain.marked.bam"));
    let (staged_header, staged_records) = bam(&staged_bam);
    assert_eq!(chain_header, staged_header, "{name}: header");
    assert_eq!(chain_records.len(), staged_records.len(), "{name}");
    assert!(chain_records == staged_records, "{name}: records differ");
    assert!(!chain_records.is_empty());
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn chained_equals_staged_on_the_fixture() {
    chained_equals_staged("fixture", FIXTURE, &[]);
    chained_equals_staged(
        "fixture-options",
        FIXTURE,
        &[
            "REMOVE_SEQUENCING_DUPLICATES=true",
            "TAGGING_POLICY=All",
            "TAG_DUPLICATE_SET_MEMBERS=true",
            "COMMENT=chained",
        ],
    );
    chained_equals_staged(
        "fixture-queryname",
        ("tests/data/duplicate_marking.bam", FIXTURE.1),
        &[],
    );
}

#[test]
fn chained_equals_staged_on_corpus_t1() {
    chained_equals_staged("t1", T1, &[]);
    chained_equals_staged("t1-remove", T1, &["REMOVE_DUPLICATES=true"]);
}

#[test]
fn chained_equals_staged_on_the_foundation_corpus() {
    chained_equals_staged("foundation", FOUNDATION, &[]);
    chained_equals_staged(
        "foundation-spill",
        FOUNDATION,
        &["MAX_RECORDS_IN_RAM=500", "TAGGING_POLICY=OpticalOnly"],
    );
}

/// With the `@PG` command line pinned the two BAMs, their indexes and
/// their MD5 files are the same bytes. The chain is spelled as Picard
/// spells it: class names and `KEY=VALUE`.
#[test]
fn the_written_bam_is_the_tools_bam_byte_for_byte() {
    let dir = fresh_dir("bytes");
    let pinned = "PROGRAM_GROUP_COMMAND_LINE=pinned";
    run(&[
        "MarkDuplicates".to_owned(),
        format!("I={}", FIXTURE.0),
        arg("O", &dir.join("staged.bam")),
        arg("M", &dir.join("staged.txt")),
        pinned.to_owned(),
        "CREATE_INDEX=true".to_owned(),
        "CREATE_MD5_FILE=true".to_owned(),
        "COMPRESSION_LEVEL=9".to_owned(),
    ]);
    run(&[
        "chain".to_owned(),
        format!("I={}", FIXTURE.0),
        "+".to_owned(),
        "MarkDuplicates".to_owned(),
        pinned.to_owned(),
        arg("M", &dir.join("chain.txt")),
        "+".to_owned(),
        "write".to_owned(),
        arg("O", &dir.join("chain.bam")),
        "CREATE_INDEX=true".to_owned(),
        "CREATE_MD5_FILE=true".to_owned(),
        "COMPRESSION_LEVEL=9".to_owned(),
    ]);
    assert_eq!(
        names(&dir),
        [
            "chain.bai",
            "chain.bam",
            "chain.bam.md5",
            "chain.txt",
            "staged.bai",
            "staged.bam",
            "staged.bam.md5",
            "staged.txt"
        ]
    );
    for file in ["bam", "bai", "bam.md5"] {
        assert_eq!(
            fs::read(dir.join(format!("chain.{file}"))).expect("chained"),
            fs::read(dir.join(format!("staged.{file}"))).expect("staged"),
            "{file}"
        );
    }
    assert_eq!(body(&dir.join("chain.txt")), body(&dir.join("staged.txt")));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A chain without a transform is the pass of `qc`: the same four metrics
/// files under the same prefix, byte for byte apart from the line with the
/// clock. Checked on a BAM and on a CRAM with its reference.
#[test]
fn a_chain_without_a_transform_equals_qc() {
    let files = [
        "s.quality_yield_metrics",
        "s.alignment_summary_metrics",
        "s.insert_size_metrics",
        "s.duplicate_metrics",
    ];
    for (name, source, reference) in [
        ("qc-bam", "tests/data/insert_size.bam", None),
        (
            "qc-bam-reference",
            "tests/data/markdup.bam",
            Some(FIXTURE.1),
        ),
        ("qc-cram", "tests/data/insert_size.cram", Some(FIXTURE.1)),
    ] {
        let dir = fresh_dir(name);
        let prefix = dir.join("run").join("s");
        let mut qc = strings(&["qc", source, "-o"]);
        qc.push(prefix.display().to_string());
        let mut chain = strings(&["chain", source, "-o"]);
        chain.push(prefix.display().to_string());
        for args in [&mut qc, &mut chain] {
            args.extend(reference.iter().flat_map(|r| strings(&["-R", r])));
        }
        run(&qc);
        let kept = dir.join("qc");
        fs::rename(dir.join("run"), &kept).expect("keep the qc files");
        chain.extend(strings(&[
            "+",
            "quality-yield",
            "+",
            "alignment-summary",
            "+",
            "insert-size",
            "+",
            "duplicate-metrics",
        ]));
        run(&chain);
        assert_eq!(names(&dir.join("run")), {
            let mut sorted = files.map(str::to_owned);
            sorted.sort();
            sorted
        });
        for file in files {
            let chained = dir.join("run").join(file);
            assert_eq!(
                without_clock(&chained),
                without_clock(&kept.join(file)),
                "{name}: {file}"
            );
            assert!(!body(&chained).is_empty(), "{name}: {file}");
        }
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}

/// The plan file and the command line of one chain write the same files,
/// byte for byte apart from the clock, and print the same plan.
#[test]
fn the_plan_file_equals_the_command_line() {
    let cases: [(&str, &[&str], &[&str]); 2] = [
        (
            "tests/data/chain_markdup.toml",
            &[FIXTURE.0, "-R", FIXTURE.1],
            &[
                "+",
                "mark-duplicates",
                "+",
                "alignment-summary",
                "+",
                "insert-size",
                "+",
                "quality-yield",
                "+",
                "duplicate-metrics",
                "+",
                "write",
                "--index",
            ],
        ),
        (
            "tests/data/chain_taps.toml",
            &[
                "tests/data/insert_size.bam",
                "--format",
                "picard,json",
                "--sample",
                "S1",
            ],
            &[
                "+",
                "quality-yield",
                "--no-use-original-qualities",
                "+",
                "alignment-summary",
                "--level",
                "SAMPLE",
                "--level",
                "LIBRARY",
                "--stop-after",
                "8",
                "+",
                "insert-size",
                "--minimum-pct",
                "0.1",
                "+",
                "duplicate-metrics",
            ],
        ),
    ];
    for (index, (plan, globals, steps)) in cases.into_iter().enumerate() {
        let dir = fresh_dir(&format!("plan-file-{index}"));
        let prefix = dir.join("run").join("s").display().to_string();
        let from_file = strings(&["chain", "--plan", plan, "-o", &prefix]);
        let mut from_line = strings(&["chain"]);
        from_line.extend(strings(globals));
        from_line.extend(strings(&["-o", &prefix]));
        from_line.extend(strings(steps));

        let dry = |args: &[String]| {
            let mut args = args.to_vec();
            args.insert(1, "--dry-run".to_owned());
            String::from_utf8(run(&args).stdout).expect("text")
        };
        let plan_text = dry(&from_file);
        assert_eq!(plan_text, dry(&from_line), "{plan}");
        assert!(plan_text.ends_with("nothing ran: this is a plan\n"));
        assert!(!dir.join("run").exists(), "a dry run writes nothing");
        // `plan` prints the same text for either form.
        let mut planned = from_file.clone();
        planned[0] = "plan".to_owned();
        assert_eq!(String::from_utf8_lossy(&run(&planned).stdout), plan_text);

        run(&from_file);
        let kept = dir.join("file");
        fs::rename(dir.join("run"), &kept).expect("keep the plan-file run");
        run(&from_line);
        let written = names(&dir.join("run"));
        assert_eq!(written, names(&kept), "{plan}");
        assert!(written.len() >= 6, "{written:?}");
        for file in written {
            let (a, b) = (dir.join("run").join(&file), kept.join(&file));
            if file.ends_with(".bam") || file.ends_with(".bai") {
                assert_eq!(fs::read(&a).expect("a"), fs::read(&b).expect("b"), "{file}");
            } else {
                assert_eq!(without_clock(&a), without_clock(&b), "{plan}: {file}");
            }
        }
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}

/// Names under the prefix (docs/07 section 2.2): a tap on the source in
/// front of a transform is `PREFIX.in.*`, the transform's metrics and the
/// kept BAM carry its label, a tap on the last stream has the name `qc`
/// gives it. `--as` renames the stream, a step's `-o` its file, and
/// `--format` and `--also` add siblings beside each metrics file.
#[test]
fn files_are_named_by_the_prefix_the_stream_and_the_step() {
    let dir = fresh_dir("names");
    let prefix = dir.join("s").display().to_string();
    let own = dir.join("own.txt").display().to_string();
    let output = run(&strings(&[
        "chain",
        FIXTURE.0,
        "-o",
        &prefix,
        "--format",
        "picard,json",
        "--also",
        "csv,markdown",
        "--sample",
        "S1",
        "+",
        "quality-yield",
        "+",
        "mark-duplicates",
        "--as",
        "dedup",
        "--remove-duplicates",
        "+",
        "insert-size",
        "+",
        "quality-yield",
        "-o",
        &own,
        "--format",
        "picard,tsv",
        "+",
        "write",
        "--md5",
    ]));
    assert_eq!(
        names(&dir),
        [
            "own.txt",
            "own.txt.csv",
            "own.txt.tsv",
            "s.dedup.bam",
            "s.dedup.bam.md5",
            "s.dedup.duplication_metrics",
            "s.dedup.duplication_metrics.csv",
            "s.dedup.duplication_metrics.json",
            "s.in.quality_yield_metrics",
            "s.in.quality_yield_metrics.csv",
            "s.in.quality_yield_metrics.json",
            "s.insert_size_metrics",
            "s.insert_size_metrics.csv",
            "s.insert_size_metrics.json",
            "s.report.md",
        ]
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.ends_with("bammetrics chain: 5 steps finished, 15 files written\n"),
        "{stderr}"
    );
    let json = text(&dir.join("s.insert_size_metrics.json"));
    assert!(json.contains("\"sample\": \"S1\""), "{json}");
    assert!(json.contains("\"tool\": \"CollectInsertSizeMetrics\""));
    let report = text(&dir.join("s.report.md"));
    assert!(report.contains("QualityYieldMetrics") && report.contains("DuplicationMetrics"));
    // Duplicates are removed: fewer records after the transform.
    assert_ne!(
        body(&dir.join("s.in.quality_yield_metrics")),
        body(&dir.join("own.txt"))
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Without `-o` the prefix is `./` plus the source name, as `qc` and
/// `collect` name their files, and standard error says so.
#[test]
fn the_default_prefix_is_the_source_name() {
    let dir = fresh_dir("default-prefix");
    let source = root().join("tests/data/insert_size.bam");
    let output = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["chain", &source.display().to_string(), "+", "quality-yield"])
        .current_dir(&dir)
        .output()
        .expect("run bammetrics");
    assert!(output.status.success());
    assert_eq!(names(&dir), ["insert_size.quality_yield_metrics"]);
    assert_eq!(
        String::from_utf8_lossy(&output.stderr),
        "bammetrics chain: no -o, files go to ./insert_size.*\n\
         bammetrics chain: 1 step finished, 1 file written\n"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// `write -` puts the BAM on standard output and `-o -` a tap's Picard
/// file: the pipe to another tool keeps working. Nothing else is printed
/// there.
#[test]
fn standard_output_takes_the_bam_or_one_metrics_file() {
    let dir = fresh_dir("stdout");
    let prefix = dir.join("s").display().to_string();
    let piped = run(&strings(&[
        "chain",
        FIXTURE.0,
        "-o",
        &prefix,
        "+",
        "mark-duplicates",
        "PROGRAM_GROUP_COMMAND_LINE=pinned",
        "+",
        "write",
        "-",
    ]));
    run(&strings(&[
        "MarkDuplicates",
        &format!("I={}", FIXTURE.0),
        &arg("O", &dir.join("tool.bam")),
        &arg("M", &dir.join("tool.txt")),
        "PROGRAM_GROUP_COMMAND_LINE=pinned",
    ]));
    assert_eq!(
        piped.stdout,
        fs::read(dir.join("tool.bam")).expect("the tool's BAM")
    );
    let metrics = run(&strings(&[
        "chain",
        FIXTURE.0,
        "-q",
        "+",
        "quality-yield",
        "-o",
        "-",
    ]));
    let printed = String::from_utf8(metrics.stdout).expect("text");
    assert!(printed.starts_with("## htsjdk.samtools.metrics.StringHeader\n"));
    assert!(printed.contains("QualityYieldMetrics"));
    assert!(metrics.stderr.is_empty());
    assert_eq!(
        names(&dir),
        ["s.marked.duplication_metrics", "tool.bam", "tool.txt"]
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A tap's `STOP_AFTER` counts the records its position emits and never
/// shortens the stream for the steps after it (docs/chaining.md).
#[test]
fn stop_after_in_the_middle_of_a_chain_belongs_to_its_tap() {
    let dir = fresh_dir("stop-after");
    let source = FIXTURE.0;
    let quality = "CollectQualityYieldMetrics";
    let prefix = dir.join("chain").display().to_string();
    run(&strings(&[
        "chain",
        source,
        "-o",
        &prefix,
        "+",
        "quality-yield",
        "--stop-after",
        "5",
        "+",
        "mark-duplicates",
        "--remove-duplicates",
        "+",
        "quality-yield",
        "STOP_AFTER=7",
        "+",
        "alignment-summary",
        "+",
        "write",
    ]));
    let staged_bam = dir.join("staged.bam");
    run(&[
        "MarkDuplicates".to_owned(),
        format!("I={source}"),
        arg("O", &staged_bam),
        arg("M", &dir.join("staged.dup.txt")),
        "REMOVE_DUPLICATES=true".to_owned(),
    ]);
    let staged: [(&str, String, &[&str], &str); 3] = [
        (
            quality,
            format!("I={source}"),
            &["STOP_AFTER=5"],
            "chain.in.quality_yield_metrics",
        ),
        (
            quality,
            arg("I", &staged_bam),
            &["STOP_AFTER=7"],
            "chain.quality_yield_metrics",
        ),
        (
            "CollectAlignmentSummaryMetrics",
            arg("I", &staged_bam),
            &[],
            "chain.alignment_summary_metrics",
        ),
    ];
    for (index, (tool, input, options, chained)) in staged.into_iter().enumerate() {
        let alone = dir.join(format!("staged.{index}.txt"));
        let mut args = vec![tool.to_owned(), input, arg("O", &alone)];
        args.extend(strings(options));
        run(&args);
        assert_eq!(body(&dir.join(chained)), body(&alone), "{chained}");
    }
    // The tap on the source writes the command line the tool alone writes.
    let on_source = dir.join("chain.in.quality_yield_metrics");
    assert_eq!(
        text(&on_source).lines().nth(1).expect("the command line"),
        format!(
            "# {quality} INPUT={source} OUTPUT={} STOP_AFTER=5",
            on_source.display()
        )
    );
    let behind = text(&dir.join("chain.quality_yield_metrics"));
    assert!(
        behind
            .lines()
            .nth(1)
            .expect("the command line")
            .ends_with(" STOP_AFTER=7 CHAIN=MarkDuplicates"),
        "{behind}"
    );
    // Seven records after the removal are not the first seven of the source.
    assert_ne!(
        body(&dir.join("chain.quality_yield_metrics")),
        body(&on_source)
    );
    let (_, chained) = bam(&dir.join("chain.marked.bam"));
    let (_, alone) = bam(&staged_bam);
    assert_eq!(chained.len(), alone.len());
    assert!(chained.len() < 36 && chained.len() > 7);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Run `chain` with `args`; the exit code and standard error.
fn chain(args: &[&str]) -> (Option<i32>, String) {
    let mut argv = vec!["chain".to_owned()];
    argv.extend(strings(args));
    let output = bammetrics(&argv);
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

#[test]
fn a_failing_tap_fails_alone_and_a_failing_transform_stops_its_segment() {
    // The BAM is queryname sorted: a tap that asks for a coordinate-sorted
    // stream fails when it is built, on the source and behind the
    // transform. Every other step writes its files.
    let dir = fresh_dir("failures");
    let prefix = dir.join("s").display().to_string();
    let (code, stderr) = chain(&[
        "tests/data/duplicate_marking.bam",
        "-o",
        &prefix,
        "+",
        "quality-yield",
        "--no-assume-sorted",
        "+",
        "mark-duplicates",
        "+",
        "quality-yield",
        "ASSUME_SORTED=false",
        "+",
        "insert-size",
        "+",
        "write",
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    let lines: Vec<&str> = stderr.lines().collect();
    assert_eq!(lines.len(), 4, "{stderr}");
    assert!(lines[0].starts_with(
        "bammetrics chain: step 1 (quality-yield): failed alone: \
         tests/data/duplicate_marking.bam should be coordinate sorted"
    ));
    assert!(lines[1].starts_with("bammetrics chain: step 2 (mark-duplicates): 13 records in"));
    assert!(lines[2].starts_with("bammetrics chain: step 3 (quality-yield): failed alone: "));
    assert_eq!(
        lines[3],
        "bammetrics chain: 3 of 5 steps finished, 2 failed; 3 files written"
    );
    assert_eq!(
        names(&dir),
        [
            "s.insert_size_metrics",
            "s.marked.bam",
            "s.marked.duplication_metrics"
        ]
    );
    fs::remove_dir_all(&dir).expect("cleanup");

    // insert_size.bam has no read group: MarkDuplicates fails on its
    // second record, as Picard does. Nothing of its segment or after it
    // is written, the tap of the source included.
    let dir = fresh_dir("failures-transform");
    let prefix = dir.join("s").display().to_string();
    let (code, stderr) = chain(&[
        "tests/data/insert_size.bam",
        "-o",
        &prefix,
        "+",
        "quality-yield",
        "+",
        "mark-duplicates",
        "+",
        "insert-size",
        "+",
        "write",
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    let cause =
        "the return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null (read r0)";
    let lines: Vec<&str> = stderr.lines().collect();
    assert_eq!(lines.len(), 5, "{stderr}");
    assert!(lines[0].starts_with(
        "bammetrics chain: step 1 (quality-yield): stopped: step 2 (mark-duplicates) failed: "
    ));
    assert!(lines[1].starts_with(
        "bammetrics chain: step 2 (mark-duplicates): failed and stopped its segment: Cannot \
         invoke"
    ));
    assert!(lines[2].starts_with(
        "bammetrics chain: step 3 (insert-size): not run: step 2 (mark-duplicates) failed: "
    ));
    assert!(lines[3].starts_with(
        "bammetrics chain: step 4 (write): not run: step 2 (mark-duplicates) failed: "
    ));
    assert!(
        lines[..4].iter().all(|line| line.ends_with(cause)),
        "{stderr}"
    );
    assert_eq!(
        lines[4],
        "bammetrics chain: 0 of 4 steps finished, 1 failed, 1 stopped, 2 not run; 0 files written"
    );
    assert!(names(&dir).is_empty(), "{:?}", names(&dir));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A pipe cannot be read twice: it is copied to `--tmp-dir` once, as the
/// tool copies it, and the copy is removed.
#[test]
fn a_pipe_source_is_copied_once_for_the_second_read() {
    use std::io::Write as _;
    use std::process::Stdio;

    let dir = fresh_dir("pipe");
    let pinned = "PROGRAM_GROUP_COMMAND_LINE=pinned";
    run(&[
        "MarkDuplicates".to_owned(),
        format!("I={}", FIXTURE.0),
        arg("O", &dir.join("file.bam")),
        arg("M", &dir.join("file.txt")),
        pinned.to_owned(),
    ]);
    let tmp = dir.join("tmp");
    fs::create_dir_all(&tmp).expect("temporary directory");
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "chain".to_owned(),
            "-".to_owned(),
            "--tmp-dir".to_owned(),
            tmp.display().to_string(),
            "-o".to_owned(),
            dir.join("pipe").display().to_string(),
            "+".to_owned(),
            "mark-duplicates".to_owned(),
            pinned.to_owned(),
            "+".to_owned(),
            "quality-yield".to_owned(),
            "+".to_owned(),
            "write".to_owned(),
        ])
        .current_dir(root())
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run bammetrics");
    let bytes = fs::read(root().join(FIXTURE.0)).expect("fixture");
    child
        .stdin
        .take()
        .expect("a pipe")
        .write_all(&bytes)
        .expect("write the BAM");
    let output = child.wait_with_output().expect("exit");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read(dir.join("pipe.marked.bam")).expect("chained"),
        fs::read(dir.join("file.bam")).expect("staged")
    );
    assert_eq!(
        body(&dir.join("pipe.marked.duplication_metrics")),
        body(&dir.join("file.txt"))
    );
    // The tap behind the transform names standard input as the verbs do.
    let line = text(&dir.join("pipe.quality_yield_metrics"));
    assert!(
        line.lines()
            .nth(1)
            .expect("the command line")
            .starts_with("# CollectQualityYieldMetrics INPUT=- OUTPUT="),
        "{line}"
    );
    assert!(names(&tmp).is_empty());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A CRAM source is decoded once per segment with the chain's `-R`: the
/// tap and the written BAM hold what the tool reads from the CRAM.
#[test]
fn a_cram_source_feeds_taps_and_write() {
    let dir = fresh_dir("cram");
    let cram = "tests/data/insert_size.cram";
    let reference = format!("R={}", FIXTURE.1);
    let tool = "CollectInsertSizeMetrics";
    let prefix = dir.join("chain").display().to_string();
    run(&strings(&[
        "chain",
        cram,
        "-R",
        FIXTURE.1,
        "-o",
        &prefix,
        "+",
        "insert-size",
        "+",
        "write",
    ]));
    run(&[
        tool.to_owned(),
        format!("I={cram}"),
        reference,
        arg("O", &dir.join("cram.txt")),
    ]);
    run(&[
        tool.to_owned(),
        arg("I", &dir.join("chain.in.bam")),
        arg("O", &dir.join("bam.txt")),
    ]);
    let chained = body(&dir.join("chain.insert_size_metrics"));
    assert_eq!(chained, body(&dir.join("cram.txt")));
    assert_eq!(chained, body(&dir.join("bam.txt")));
    assert_eq!(bam(&dir.join("chain.in.bam")).1.len(), 10);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Two steps that name one file are refused before the first read, so no
/// step replaces the file of another (the metrics are written after the
/// run, which would put a tap's text over the renamed BAM).
#[test]
fn two_steps_cannot_write_one_file() {
    let dir = fresh_dir("one-file");
    let source = FIXTURE.0;
    let file = |name: &str| dir.join(name).display().to_string();
    let prefix = file("s");
    let x = file("x.out");
    let bam = file("x.bam");
    let bai = file("x.bai");
    let md5 = file("x.bam.md5");
    let gc = file("gc");
    let cases: [(&[&str], &str); 9] = [
        (
            &[
                "+",
                "mark-duplicates",
                "+",
                "quality-yield",
                "-o",
                &x,
                "+",
                "write",
                &x,
            ],
            "step 2 (quality-yield) and step 3 (write) both write",
        ),
        (
            &["+", "quality-yield", "-o", &x, "+", "insert-size", "-o", &x],
            "step 1 (quality-yield) and step 2 (insert-size) both write",
        ),
        // Two taps of one tool at one position, both named by the prefix.
        (
            &["+", "quality-yield", "+", "quality-yield"],
            "step 1 (quality-yield) and step 2 (quality-yield) both write",
        ),
        (
            &[
                "+",
                "mark-duplicates",
                "--metrics-file",
                &x,
                "+",
                "quality-yield",
                "-o",
                &x,
            ],
            "step 1 (mark-duplicates) and step 2 (quality-yield) both write",
        ),
        (
            &["+", "write", &bam, "+", "write", &bam],
            "step 1 (write) and step 2 (write) both write",
        ),
        // The index and the MD5 of `write` are its files too.
        (
            &[
                "+",
                "write",
                &bam,
                "--index",
                "+",
                "quality-yield",
                "-o",
                &bai,
            ],
            "step 1 (write) and step 2 (quality-yield) both write",
        ),
        (
            &[
                "+",
                "quality-yield",
                "-o",
                &md5,
                "+",
                "write",
                &bam,
                "--md5",
            ],
            "step 1 (quality-yield) and step 2 (write) both write",
        ),
        // The second file of a tool (`SUMMARY_OUTPUT`).
        (
            &[
                "-R",
                FIXTURE.1,
                "+",
                "gc-bias",
                "-o",
                &gc,
                "--summary-output",
                &x,
                "+",
                "quality-yield",
                "-o",
                &x,
            ],
            "step 1 (gc-bias) and step 2 (quality-yield) both write",
        ),
        // A further format beside one file is the file of another step.
        (
            &[
                "+",
                "quality-yield",
                "-o",
                &x,
                "--format",
                "picard,json",
                "+",
                "insert-size",
                "-o",
                &file("x.out.json"),
            ],
            "step 1 (quality-yield) and step 2 (insert-size) both write",
        ),
    ];
    for (steps, text) in cases {
        let mut args = vec![source, "-o", &prefix];
        args.extend(steps);
        let (code, stderr) = chain(&args);
        assert_eq!(code, Some(2), "{args:?}: {stderr}");
        assert!(stderr.starts_with("bammetrics chain: "), "{stderr}");
        assert!(stderr.contains(text), "{args:?}: {stderr}");
        assert!(
            stderr.contains("every output of a chain needs its own file"),
            "{stderr}"
        );
        assert!(names(&dir).is_empty(), "{args:?}: {:?}", names(&dir));
    }
    // One step with its own -o: two taps of one tool at one position run.
    let (code, stderr) = chain(&[
        source,
        "-o",
        &prefix,
        "+",
        "quality-yield",
        "+",
        "quality-yield",
        "--no-use-original-qualities",
        "-o",
        &x,
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(names(&dir), ["s.quality_yield_metrics", "x.out"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A tool whose `OUTPUT` is a prefix names its files when it finishes. One
/// that another step writes fails that tool alone, and the file keeps the
/// other step's content.
#[test]
fn a_file_under_a_prefix_never_replaces_the_file_of_another_step() {
    let dir = fresh_dir("prefix");
    let taken = dir.join("art.pre_adapter_detail_metrics");
    let (code, stderr) = chain(&[
        "tests/data/artifact.bam",
        "-R",
        FIXTURE.1,
        "-o",
        &dir.join("art").display().to_string(),
        "+",
        "sequencing-artifact",
        "+",
        "quality-yield",
        "-o",
        &taken.display().to_string(),
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "bammetrics chain: step 1 (sequencing-artifact): failed alone: step 2 \
             (quality-yield) and step 1 (sequencing-artifact) both write {}",
            taken.display()
        )),
        "{stderr}"
    );
    assert_eq!(names(&dir), ["art.pre_adapter_detail_metrics"]);
    assert!(body(&taken).contains("QualityYieldMetrics"));
    fs::remove_dir_all(&dir).expect("cleanup");

    // With no clash the five files land under the prefix, with a sibling
    // per further format.
    let dir = fresh_dir("prefix-ok");
    let (code, stderr) = chain(&[
        "tests/data/artifact.bam",
        "-R",
        FIXTURE.1,
        "-o",
        &dir.join("art").display().to_string(),
        "--format",
        "picard,csv",
        "+",
        "sequencing-artifact",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    let written = names(&dir);
    assert_eq!(written.len(), 10, "{written:?}");
    assert!(written.contains(&"art.error_summary_metrics.csv".to_owned()));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// What `write` cannot do is said: an index of a stream that is not
/// coordinate sorted is not written, with a note.
#[test]
fn write_says_when_it_writes_no_index() {
    let dir = fresh_dir("no-index");
    let (code, stderr) = chain(&[
        "tests/data/duplicate_marking.bam",
        "+",
        "write",
        &dir.join("q.bam").display().to_string(),
        "--index",
        "--md5",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "bammetrics chain: step 1 (write): CREATE_INDEX: no index written, the sort order is \
         queryname and a .bai needs coordinate\n\
         bammetrics chain: 1 step finished, 2 files written\n"
    );
    assert_eq!(names(&dir), ["q.bam", "q.bam.md5"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Missing directories above an output are created, as every verb creates
/// them. A file that cannot be written is found before the first read and
/// exits 1 under the step's label, as the Picard face exits for its
/// `OUTPUT`.
#[test]
fn an_output_that_cannot_be_written_exits_1_before_anything_runs() {
    let dir = fresh_dir("unwritable");
    let source = FIXTURE.0;
    let blocker = dir.join("blocker");
    fs::write(&blocker, "a file, not a directory").expect("blocker");
    let under = |name: &str| blocker.join(name).display().to_string();
    let good = dir.join("good.txt").display().to_string();
    let cases: [(&[&str], &str); 3] = [
        (
            &[
                "+",
                "quality-yield",
                "-o",
                &good,
                "+",
                "write",
                &under("x.bam"),
            ],
            "bammetrics chain: step 2 (write): ",
        ),
        (
            &[
                "+",
                "mark-duplicates",
                "--metrics-file",
                &good,
                "+",
                "quality-yield",
                "-o",
                &under("x.txt"),
            ],
            "bammetrics chain: step 2 (quality-yield): ",
        ),
        (
            &["+", "mark-duplicates", "--metrics-file", &under("dup.txt")],
            "bammetrics chain: step 1 (mark-duplicates): ",
        ),
    ];
    for (steps, start) in cases {
        let mut args = vec![source];
        args.extend(steps);
        let (code, stderr) = chain(&args);
        assert_eq!(code, Some(1), "{args:?}: {stderr}");
        assert!(stderr.starts_with(start), "{stderr}");
        assert_eq!(stderr.lines().count(), 1, "{stderr}");
        assert_eq!(names(&dir), ["blocker"], "{args:?}");
    }
    let deep = dir.join("made").join("deeper").join("s");
    let (code, stderr) = chain(&[
        source,
        "-o",
        &deep.display().to_string(),
        "+",
        "quality-yield",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        names(&dir.join("made").join("deeper")),
        ["s.quality_yield_metrics"]
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Every usage error exits 2 before the first read, names the step by its
/// number and verb, and writes nothing.
#[test]
fn usage_errors_exit_2_before_anything_runs() {
    let dir = fresh_dir("usage");
    let file = |name: &str| dir.join(name).display().to_string();
    let prefix = file("s");
    let sam = file("out.sam");
    let cram = file("out.cram");
    let fastq = file("out.fastq.gz");
    let source = FIXTURE.0;
    let cases: [(&[&str], &str); 24] = [
        (&[source, "+", "write", &sam], "step 1 (write): "),
        (
            &[source, "+", "mark-duplicates", "+", "write", &cram],
            "out.cram: bammetrics writes BAM only: CRAM output is wave W15",
        ),
        (
            &[source, "+", "write", &fastq],
            "out.fastq.gz: bammetrics writes BAM only: FASTQ output is wave W22",
        ),
        // Refused from the command line: the first mark-duplicates reads
        // nothing for a chain that cannot finish.
        (
            &[
                source,
                "+",
                "mark-duplicates",
                "+",
                "quality-yield",
                "+",
                "mark-duplicates",
            ],
            "step 3 (mark-duplicates): it reads its input twice and cannot follow step 1 \
             (mark-duplicates)",
        ),
        (&[source], "no step: a chain needs at least one + STEP"),
        (
            &[source, "+", "mark-duplicates", "+", "write", source],
            "step 2 (write): tests/data/markdup.bam is an input; writing it would destroy it",
        ),
        (&["+", "write", &sam], "missing BAM"),
        (
            &[source, "+", "NoSuchTool"],
            "step 1: unknown step \"NoSuchTool\"",
        ),
        (
            &[source, "+", "inser-size"],
            "step 1: unknown step \"inser-size\"\n  did you mean insert-size?",
        ),
        (
            &[source, "+", "wgs"],
            "step 1 (wgs): CollectWgsMetrics has no stage and does not read a BAM in one pass",
        ),
        (
            &[
                source,
                "+",
                "quality-yield",
                "+",
                "insert-size",
                "--minmum-pct",
                "0.1",
            ],
            "step 2 (insert-size): unknown option --minmum-pct\n  did you mean --minimum-pct \
             (Picard MINIMUM_PCT)?",
        ),
        (
            &[source, "+", "insert-size", "NO_SUCH=1"],
            "step 1 (insert-size): unknown option NO_SUCH=",
        ),
        (
            &[source, "+", "insert-size", "--minimum-pct", "abc"],
            "step 1 (insert-size): MINIMUM_PCT=abc: not a valid value (from --minimum-pct)",
        ),
        (
            &[source, "+", "quality-yield", "--stop-after", "many"],
            "step 1 (quality-yield): STOP_AFTER=many: not a valid value (from --stop-after)",
        ),
        // The options of a transform are checked here, not when it is built.
        (
            &[source, "+", "mark-duplicates", "--read-name-regex", "("],
            "step 1 (mark-duplicates): READ_NAME_REGEX",
        ),
        (
            &[source, "+", "mark-duplicates", "--create-index"],
            "CREATE_INDEX belongs to the write step",
        ),
        (
            &[source, "+", "mark-duplicates", "-o", &file("x.bam")],
            "step 1 (mark-duplicates): -o names a BAM, and a transform writes none in a chain",
        ),
        (
            &[source, "+", "insert-size", "tests/data/insert_size.bam"],
            "step 1 (insert-size): a step takes no input",
        ),
        (
            &[source, "+", "insert-size", "--threads", "2"],
            "step 1 (insert-size): --threads belongs to the chain",
        ),
        (
            &[source, "+", "hs-metrics"],
            "step 1 (hs-metrics): missing --target-intervals (Picard: TARGET_INTERVALS)",
        ),
        (
            &[source, "--format", "json", "+", "insert-size"],
            "--format json starts with json: a chain names Picard files",
        ),
        (
            &[source, "--also", "junit", "+", "insert-size"],
            "--also junit prints a verdict, and a chain has none",
        ),
        (
            &[source, "+", "quality-yield", "-o", "-", "+", "write", "-"],
            "step 1 (quality-yield) and step 2 (write) both write standard output",
        ),
        (
            &[source, "+", "mark-duplicates", "--as", "in"],
            "step 1 (mark-duplicates): the stream name in is taken",
        ),
    ];
    for (args, text) in cases {
        let mut argv = vec!["-o", prefix.as_str()];
        argv.extend(args);
        let (code, stderr) = chain(&argv);
        assert_eq!(code, Some(2), "{argv:?}: {stderr}");
        assert!(stderr.starts_with("bammetrics chain: "), "{stderr}");
        assert!(stderr.contains(text), "{argv:?}: {stderr}");
    }
    assert!(names(&dir).is_empty());
    // An input that is not there exits 1, as on every verb.
    for args in [
        &["tests/data/nope.bam", "+", "quality-yield"][..],
        &[source, "-R", "tests/data/nope.fa", "+", "quality-yield"],
        &["--plan", "tests/data/nope.toml"],
    ] {
        let (code, stderr) = chain(args);
        assert_eq!(code, Some(1), "{args:?}: {stderr}");
        assert!(
            stderr.starts_with("bammetrics chain: tests/data/nope."),
            "{stderr}"
        );
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// The hidden `__chain` entry of 0.21.0 is gone, and the verb is listed
/// where the other verbs are.
#[test]
fn the_verb_is_listed_and_the_test_entry_is_gone() {
    let output = bammetrics(&strings(&["__chain", "I=tests/data/markdup.bam"]));
    assert_eq!(output.status.code(), Some(2));
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .starts_with("bammetrics: unknown verb or Picard class \"__chain\"")
    );
    let tools = String::from_utf8(run(&strings(&["tools"])).stdout).expect("text");
    assert!(tools.contains("\n  chain        Several steps on one BAM in one process"));
    let help = String::from_utf8(run(&strings(&["help", "chain"])).stdout).expect("text");
    assert!(help.starts_with("bammetrics chain: several steps on one BAM in one process."));
    assert_eq!(
        String::from_utf8(run(&strings(&["chain", "--help"])).stdout).expect("text"),
        help
    );
    let typo = bammetrics(&strings(&["chian", "x.bam"]));
    assert!(String::from_utf8_lossy(&typo.stderr).contains("did you mean chain?"));
    for shell in ["bash", "zsh", "fish"] {
        let script =
            String::from_utf8(run(&strings(&["completions", shell])).stdout).expect("a script");
        assert!(script.contains("chain"), "{shell}");
        assert!(script.contains("dry-run"), "{shell}");
    }
}

/// The library entry a caller that is not the command line uses (the
/// Python module): a `Spec` in, `prepare`, then `execute`, which returns
/// one outcome per step with its status, its files and, for a transform,
/// its record counters. Nothing is printed.
#[test]
fn the_library_entry_returns_one_outcome_per_step() {
    use bammetrics_cli::chain::{Role, Spec, StepSpec, StepStatus, execute, prepare};

    let dir = fresh_dir("library");
    let source = root().join(FIXTURE.0).display().to_string();
    let spec = Spec {
        source: Some(source.clone()),
        prefix: Some(dir.join("s").display().to_string()),
        threads: Some(2),
        steps: vec![
            StepSpec::new("quality-yield", &[]),
            StepSpec::new("mark-duplicates", &["--remove-duplicates"]),
            StepSpec::new("insert-size", &["--minimum-pct=0.1"]),
            StepSpec::new("write", &["--index"]),
        ],
        ..Spec::default()
    };
    let chain = prepare(&spec).expect("a chain");
    assert_eq!(chain.segments(), 2);
    assert!(chain.replayable() && chain.prefix_used);
    // Each step: its number, verb, role, the stream it reads, the stream
    // it leaves and its segment.
    type Step<'a> = (usize, &'a str, Role, &'a str, Option<&'a str>, usize);
    let planned: Vec<Step> = chain
        .steps
        .iter()
        .map(|s| {
            let leaves = s.leaves.as_deref();
            (
                s.number,
                s.verb,
                s.role,
                s.reads.as_str(),
                leaves,
                s.segment,
            )
        })
        .collect();
    assert_eq!(
        planned,
        [
            (1, "quality-yield", Role::Tap, "in", None, 1),
            (
                2,
                "mark-duplicates",
                Role::Transform,
                "in",
                Some("marked"),
                1
            ),
            (3, "insert-size", Role::Tap, "marked", None, 2),
            (4, "write", Role::Write, "marked", None, 2),
        ]
    );
    assert_eq!(chain.steps[2].label, "step 3 (insert-size)");
    let insert_size = dir.join("s.insert_size_metrics");
    assert_eq!(
        chain.steps[2].command,
        Some(format!(
            "CollectInsertSizeMetrics INPUT={source} OUTPUT={} MINIMUM_PCT=0.1 \
             CHAIN=MarkDuplicates",
            insert_size.display()
        ))
    );
    let kinds: Vec<&str> = chain.steps[3].files.iter().map(|f| f.kind).collect();
    assert_eq!(kinds, ["bam", "index"]);
    assert!(names(&dir).is_empty(), "prepare writes nothing");

    let outcome = execute(chain).expect("the chain runs");
    assert_eq!(outcome.code, 0);
    assert_eq!(outcome.steps.len(), 4);
    assert!(outcome.steps.iter().all(|s| s.status == StepStatus::Ok));
    let counters = outcome.steps[1].counters.expect("a transform counts");
    assert_eq!(counters.records_in, 36);
    assert!(counters.records_out < 36 && counters.records_out > 7);
    assert!(outcome.steps[0].counters.is_none());
    let written: Vec<(&str, String)> = outcome
        .steps
        .iter()
        .flat_map(|s| &s.files)
        .map(|f| {
            let name = f.path.file_name().expect("a name").to_string_lossy();
            (f.kind, name.into_owned())
        })
        .collect();
    assert_eq!(
        written,
        [
            ("metrics", "s.in.quality_yield_metrics".to_owned()),
            ("metrics", "s.marked.duplication_metrics".to_owned()),
            ("metrics", "s.insert_size_metrics".to_owned()),
            ("bam", "s.marked.bam".to_owned()),
            ("index", "s.marked.bai".to_owned()),
        ]
    );
    assert_eq!(names(&dir).len(), 5);

    // A step that fails is an outcome, not an error: the exit code is 1.
    let no_read_group = root().join("tests/data/insert_size.bam");
    let failing = Spec {
        source: Some(no_read_group.display().to_string()),
        prefix: Some(dir.join("t").display().to_string()),
        steps: vec![
            StepSpec::new("quality-yield", &[]),
            StepSpec::new("mark-duplicates", &[]),
            StepSpec::new("write", &[]),
        ],
        ..Spec::default()
    };
    let outcome = execute(prepare(&failing).expect("a chain")).expect("it runs");
    assert_eq!(outcome.code, 1);
    let statuses: Vec<StepStatus> = outcome.steps.iter().map(|s| s.status).collect();
    assert_eq!(
        statuses,
        [StepStatus::Stopped, StepStatus::Failed, StepStatus::NotRun]
    );
    let cause = outcome.steps[1].error.as_deref().unwrap_or_default();
    assert!(cause.starts_with("Cannot invoke"), "{cause}");
    // A usage error is found by prepare, with the step in the message.
    let bad = Spec {
        steps: vec![StepSpec::new("insert-size", &["--minmum-pct", "1"])],
        ..failing
    };
    let Err(error) = prepare(&bad) else {
        panic!("an unknown option");
    };
    assert_eq!(error.code, 2);
    assert!(
        error
            .what
            .starts_with("step 1 (insert-size): unknown option --minmum-pct")
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

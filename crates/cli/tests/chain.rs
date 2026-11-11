//! `__chain` end to end (docs/chaining.md): chained equals staged. One
//! process running `MarkDuplicates`, the one-pass tools on its records and
//! `write` gives the metrics files and the BAM that `MarkDuplicates` to a
//! file followed by each tool on that file gives. On the committed fixture
//! always; on corpus t1 and the foundation corpus when `tests/corpus/out`
//! holds them (`make corpus`, `make concordance`). Also: a tap's
//! `STOP_AFTER` in the middle of a chain, and what a failing step does to
//! the others.

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

/// The taps: a one-pass tool, the option that names its file, and whether
/// it reads the reference.
const TAPS: [(&str, &str, bool); 8] = [
    ("CollectAlignmentSummaryMetrics", "O", true),
    ("CollectInsertSizeMetrics", "O", false),
    ("CollectQualityYieldMetrics", "O", false),
    ("CollectDuplicateMetrics", "M", false),
    ("QualityScoreDistribution", "O", false),
    ("MeanQualityByCycle", "O", false),
    ("CollectBaseDistributionByCycle", "O", false),
    ("CollectGcBiasMetrics", "O", true),
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

/// Run and require exit 0.
fn run(args: &[String]) {
    let output = bammetrics(args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
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

/// A metrics file without its `StringHeader` lines (the command line and
/// the start time).
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.split_inclusive('\n').peekable();
    while lines.peek().is_some_and(|line| line.starts_with('#')) {
        lines.next();
    }
    lines.collect()
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

/// The Picard arguments of tap `tool` writing under `dir` with file
/// prefix `side`.
fn tap_args(
    tool: &str,
    key: &str,
    reads_reference: bool,
    reference: &str,
    dir: &Path,
    side: &str,
) -> Vec<String> {
    let mut args = vec![arg(key, &dir.join(format!("{side}.{tool}.txt")))];
    if reads_reference {
        args.push(format!("R={reference}"));
    }
    if tool == "CollectGcBiasMetrics" {
        args.push(arg("S", &dir.join(format!("{side}.{tool}.summary.txt"))));
    }
    args
}

/// `MarkDuplicates` to a file and each tap on that file, against one
/// chain. `options` go to `MarkDuplicates` on both sides.
fn chained_equals_staged(name: &str, input: (&str, &str), options: &[&str]) {
    let (source, reference) = input;
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus, make concordance)");
        return;
    }
    let dir = fresh_dir(name);
    let options: Vec<String> = options.iter().map(|o| (*o).to_owned()).collect();
    let staged_bam = dir.join("staged.bam");
    let mut staged = vec![
        "MarkDuplicates".to_owned(),
        format!("I={source}"),
        arg("O", &staged_bam),
        arg("M", &dir.join("staged.MarkDuplicates.txt")),
    ];
    staged.extend(options.iter().cloned());
    run(&staged);
    let mut chain = vec![
        "__chain".to_owned(),
        format!("I={source}"),
        "STEP=MarkDuplicates".to_owned(),
        arg("M", &dir.join("chain.MarkDuplicates.txt")),
    ];
    chain.extend(options.iter().cloned());
    for (tool, key, reads_reference) in TAPS {
        let mut alone = vec![tool.to_owned(), arg("I", &staged_bam)];
        alone.extend(tap_args(
            tool,
            key,
            reads_reference,
            reference,
            &dir,
            "staged",
        ));
        run(&alone);
        chain.push(format!("STEP={tool}"));
        chain.extend(tap_args(
            tool,
            key,
            reads_reference,
            reference,
            &dir,
            "chain",
        ));
    }
    let chain_bam = dir.join("chain.bam");
    chain.extend(["STEP=write".to_owned(), arg("O", &chain_bam)]);
    run(&chain);

    let files = names(&dir);
    let staged_files: Vec<&String> = files
        .iter()
        .filter(|f| f.starts_with("staged.") && f.ends_with(".txt"))
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
    let (chain_header, chain_records) = bam(&chain_bam);
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
/// their MD5 files are the same bytes.
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
        "__chain".to_owned(),
        format!("I={}", FIXTURE.0),
        "STEP=MarkDuplicates".to_owned(),
        pinned.to_owned(),
        "STEP=write".to_owned(),
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
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A tap's `STOP_AFTER` counts the records its position emits and never
/// shortens the stream for the steps after it (docs/chaining.md).
#[test]
fn stop_after_in_the_middle_of_a_chain_belongs_to_its_tap() {
    let dir = fresh_dir("stop-after");
    let source = FIXTURE.0;
    let quality = "CollectQualityYieldMetrics";
    run(&[
        "__chain".to_owned(),
        format!("I={source}"),
        format!("STEP={quality}"),
        arg("O", &dir.join("chain.source.txt")),
        "STOP_AFTER=5".to_owned(),
        "STEP=MarkDuplicates".to_owned(),
        "REMOVE_DUPLICATES=true".to_owned(),
        format!("STEP={quality}"),
        arg("O", &dir.join("chain.marked.txt")),
        "STOP_AFTER=7".to_owned(),
        "STEP=CollectAlignmentSummaryMetrics".to_owned(),
        arg("O", &dir.join("chain.all.txt")),
        "STEP=write".to_owned(),
        arg("O", &dir.join("chain.bam")),
    ]);
    let staged_bam = dir.join("staged.bam");
    run(&[
        "MarkDuplicates".to_owned(),
        format!("I={source}"),
        arg("O", &staged_bam),
        arg("M", &dir.join("staged.dup.txt")),
        "REMOVE_DUPLICATES=true".to_owned(),
    ]);
    let staged: [(&str, String, &[&str]); 3] = [
        (quality, format!("I={source}"), &["STOP_AFTER=5"]),
        (quality, arg("I", &staged_bam), &["STOP_AFTER=7"]),
        ("CollectAlignmentSummaryMetrics", arg("I", &staged_bam), &[]),
    ];
    for ((tool, input, options), name) in staged.into_iter().zip(["source", "marked", "all"]) {
        let mut args = vec![
            tool.to_owned(),
            input,
            arg("O", &dir.join(format!("staged.{name}.txt"))),
        ];
        args.extend(options.iter().map(|o| (*o).to_owned()));
        run(&args);
        assert_eq!(
            body(&dir.join(format!("chain.{name}.txt"))),
            body(&dir.join(format!("staged.{name}.txt"))),
            "{name}"
        );
    }
    // The tap on the source writes the command line the tool alone writes.
    let line = fs::read_to_string(dir.join("chain.source.txt")).expect("file");
    let command = line.lines().nth(1).expect("the command line");
    assert_eq!(
        command,
        format!(
            "# {quality} INPUT={source} OUTPUT={} STOP_AFTER=5",
            dir.join("chain.source.txt").display()
        )
    );
    let line = fs::read_to_string(dir.join("chain.marked.txt")).expect("file");
    assert!(
        line.lines()
            .nth(1)
            .expect("the command line")
            .ends_with(" STOP_AFTER=7 CHAIN=MarkDuplicates"),
        "{line}"
    );
    // Seven records after the removal are not the first seven of the source.
    assert_ne!(
        body(&dir.join("chain.marked.txt")),
        body(&dir.join("chain.source.txt"))
    );
    let (_, chained) = bam(&dir.join("chain.bam"));
    let (_, alone) = bam(&staged_bam);
    assert_eq!(chained.len(), alone.len());
    assert!(chained.len() < 36 && chained.len() > 7);
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_failing_tap_fails_alone_and_a_failing_stage_stops_its_segment() {
    let dir = fresh_dir("failures");
    // CollectGcBiasMetrics has no reference: it fails, the others write.
    let output = bammetrics(&[
        "__chain".to_owned(),
        format!("I={}", FIXTURE.0),
        "STEP=CollectGcBiasMetrics".to_owned(),
        arg("O", &dir.join("gc.txt")),
        arg("S", &dir.join("gc.summary.txt")),
        "STEP=MarkDuplicates".to_owned(),
        "STEP=CollectGcBiasMetrics".to_owned(),
        arg("O", &dir.join("gc2.txt")),
        arg("S", &dir.join("gc2.summary.txt")),
        "STEP=CollectInsertSizeMetrics".to_owned(),
        arg("O", &dir.join("insert.txt")),
        "STEP=write".to_owned(),
        arg("O", &dir.join("marked.bam")),
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    // Each failure names its step once: a name given twice is numbered.
    for step in ["CollectGcBiasMetrics#1", "CollectGcBiasMetrics#2"] {
        let line = format!("bammetrics __chain: {step}: REFERENCE_SEQUENCE is required");
        assert!(stderr.lines().any(|l| l.starts_with(&line)), "{stderr}");
    }
    assert_eq!(names(&dir), ["insert.txt", "marked.bam"]);
    fs::remove_dir_all(&dir).expect("cleanup");

    // insert_size.bam has no read group: MarkDuplicates fails on its
    // second record, as Picard does. Nothing of its segment or after it
    // is written, the tap of the source included.
    let dir = fresh_dir("failures-stage");
    let output = bammetrics(&[
        "__chain".to_owned(),
        "I=tests/data/insert_size.bam".to_owned(),
        "STEP=CollectQualityYieldMetrics".to_owned(),
        arg("O", &dir.join("source.txt")),
        "STEP=MarkDuplicates".to_owned(),
        arg("M", &dir.join("dup.txt")),
        "STEP=CollectInsertSizeMetrics".to_owned(),
        arg("O", &dir.join("insert.txt")),
        "STEP=write".to_owned(),
        arg("O", &dir.join("marked.bam")),
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    let cause =
        "the return value of \"htsjdk.samtools.SAMRecord.getReadGroup()\" is null (read r0)";
    let lines: Vec<&str> = stderr.lines().collect();
    assert_eq!(lines.len(), 4, "{stderr}");
    assert!(lines[0].starts_with(
        "bammetrics __chain: CollectQualityYieldMetrics: stopped: MarkDuplicates failed: "
    ));
    assert!(lines[1].starts_with("bammetrics __chain: MarkDuplicates: Cannot invoke"));
    assert!(lines[2].starts_with(
        "bammetrics __chain: CollectInsertSizeMetrics: not run: MarkDuplicates failed: "
    ));
    assert!(lines[3].starts_with("bammetrics __chain: write: not run: MarkDuplicates failed: "));
    assert!(lines.iter().all(|line| line.ends_with(cause)), "{stderr}");
    assert!(names(&dir).is_empty(), "{:?}", names(&dir));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A pipe cannot be read twice: it is copied to `TMP_DIR` once, as the
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
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "__chain".to_owned(),
            "I=/dev/stdin".to_owned(),
            arg("TMP_DIR", &dir.join("tmp")),
            "STEP=MarkDuplicates".to_owned(),
            pinned.to_owned(),
            arg("M", &dir.join("pipe.txt")),
            "STEP=CollectQualityYieldMetrics".to_owned(),
            arg("O", &dir.join("pipe.quality.txt")),
            "STEP=write".to_owned(),
            arg("O", &dir.join("pipe.bam")),
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
        fs::read(dir.join("pipe.bam")).expect("chained"),
        fs::read(dir.join("file.bam")).expect("staged")
    );
    assert_eq!(body(&dir.join("pipe.txt")), body(&dir.join("file.txt")));
    assert!(names(&dir.join("tmp")).is_empty());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A CRAM source is decoded once per segment with the chain's `R=`: the
/// tap and the written BAM hold what the tool reads from the CRAM.
#[test]
fn a_cram_source_feeds_taps_and_write() {
    let dir = fresh_dir("cram");
    let cram = "tests/data/insert_size.cram";
    let reference = format!("R={}", FIXTURE.1);
    let tool = "CollectInsertSizeMetrics";
    run(&[
        "__chain".to_owned(),
        format!("I={cram}"),
        reference.clone(),
        format!("STEP={tool}"),
        arg("O", &dir.join("chain.txt")),
        "STEP=write".to_owned(),
        arg("O", &dir.join("chain.bam")),
    ]);
    run(&[
        tool.to_owned(),
        format!("I={cram}"),
        reference,
        arg("O", &dir.join("cram.txt")),
    ]);
    run(&[
        tool.to_owned(),
        arg("I", &dir.join("chain.bam")),
        arg("O", &dir.join("bam.txt")),
    ]);
    let chained = body(&dir.join("chain.txt"));
    assert_eq!(chained, body(&dir.join("cram.txt")));
    assert_eq!(chained, body(&dir.join("bam.txt")));
    assert_eq!(bam(&dir.join("chain.bam")).1.len(), 10);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Run `__chain` with `args`; the exit code and standard error.
fn chain(args: &[&str]) -> (Option<i32>, String) {
    let mut argv = vec!["__chain".to_owned()];
    argv.extend(args.iter().map(|a| (*a).to_owned()));
    let output = bammetrics(&argv);
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// Two steps that name one file are refused before the first read, so no
/// step replaces the file of another (the metrics are written after the
/// run, which would put a tap's text over the renamed BAM).
#[test]
fn two_steps_cannot_write_one_file() {
    let dir = fresh_dir("one-file");
    let input = format!("I={}", FIXTURE.0);
    let quality = "STEP=CollectQualityYieldMetrics";
    let insert = "STEP=CollectInsertSizeMetrics";
    let gc = "STEP=CollectGcBiasMetrics";
    let reference = format!("R={}", FIXTURE.1);
    let file = |name: &str, key: &str| arg(key, &dir.join(name));
    let x = file("x.out", "O");
    let bam = file("x.bam", "O");
    let cases: [(&[&str], &str); 8] = [
        (
            &["STEP=MarkDuplicates", quality, &x, "STEP=write", &x],
            "CollectQualityYieldMetrics and write both write",
        ),
        (
            &[quality, &x, insert, &x],
            "CollectQualityYieldMetrics and CollectInsertSizeMetrics both write",
        ),
        (
            &[quality, &x, quality, &x],
            "CollectQualityYieldMetrics#1 and CollectQualityYieldMetrics#2 both write",
        ),
        (
            &["STEP=MarkDuplicates", &file("x.out", "M"), quality, &x],
            "MarkDuplicates and CollectQualityYieldMetrics both write",
        ),
        (
            &["STEP=write", &bam, "STEP=write", &bam],
            "write#1 and write#2 both write",
        ),
        // The index and the MD5 of `write` are its files too.
        (
            &[
                "STEP=write",
                &bam,
                "CREATE_INDEX=true",
                quality,
                &file("x.bai", "O"),
            ],
            "write and CollectQualityYieldMetrics both write",
        ),
        (
            &[
                quality,
                &file("x.bam.md5", "O"),
                "STEP=write",
                &bam,
                "CREATE_MD5_FILE=true",
            ],
            "CollectQualityYieldMetrics and write both write",
        ),
        // The second file of a tool (`SUMMARY_OUTPUT`).
        (
            &[
                gc,
                &file("gc.txt", "O"),
                &file("x.out", "S"),
                &reference,
                quality,
                &x,
            ],
            "CollectGcBiasMetrics and CollectQualityYieldMetrics both write",
        ),
    ];
    for (steps, text) in cases {
        let mut args = vec![input.as_str()];
        args.extend(steps);
        let (code, stderr) = chain(&args);
        assert_eq!(code, Some(2), "{args:?}: {stderr}");
        assert!(stderr.starts_with("bammetrics __chain: "), "{stderr}");
        assert!(stderr.contains(text), "{args:?}: {stderr}");
        assert!(
            stderr.contains("every output of a chain needs its own file"),
            "{stderr}"
        );
        assert!(names(&dir).is_empty(), "{args:?}: {:?}", names(&dir));
    }
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
        "I=tests/data/artifact.bam",
        &format!("R={}", FIXTURE.1),
        "STEP=CollectSequencingArtifactMetrics",
        &arg("O", &dir.join("art")),
        "STEP=CollectQualityYieldMetrics",
        &arg("O", &taken),
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!(
            "bammetrics __chain: CollectSequencingArtifactMetrics: {} is a file of \
             CollectQualityYieldMetrics",
            taken.display()
        )),
        "{stderr}"
    );
    assert_eq!(names(&dir), ["art.pre_adapter_detail_metrics"]);
    assert!(body(&taken).contains("QualityYieldMetrics"));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// What `write` cannot do is said: an index of a stream that is not
/// coordinate sorted is not written, with a note.
#[test]
fn write_says_when_it_writes_no_index() {
    let dir = fresh_dir("no-index");
    let (code, stderr) = chain(&[
        "I=tests/data/duplicate_marking.bam",
        "STEP=write",
        &arg("O", &dir.join("q.bam")),
        "CREATE_INDEX=true",
        "CREATE_MD5_FILE=true",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        stderr,
        "bammetrics __chain: write: CREATE_INDEX: no index written, the sort order is queryname \
         and a .bai needs coordinate\n"
    );
    assert_eq!(names(&dir), ["q.bam", "q.bam.md5"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A file that cannot be written is found before the first read and exits
/// 1 under the step's label, as the Picard face exits for its `OUTPUT`.
#[test]
fn an_output_that_cannot_be_written_exits_1_before_anything_runs() {
    let dir = fresh_dir("unwritable");
    let input = format!("I={}", FIXTURE.0);
    let missing = dir.join("missing");
    let dup = arg("M", &dir.join("dup.txt"));
    let cases: [(&[&str], &str); 3] = [
        (
            &["STEP=write", &arg("O", &missing.join("x.bam"))],
            "bammetrics __chain: write: ",
        ),
        (
            &[
                "STEP=MarkDuplicates",
                &dup,
                "STEP=CollectQualityYieldMetrics",
                &arg("O", &missing.join("x.txt")),
            ],
            "bammetrics __chain: CollectQualityYieldMetrics: ",
        ),
        (
            &["STEP=MarkDuplicates", &arg("M", &missing.join("dup.txt"))],
            "bammetrics __chain: MarkDuplicates: ",
        ),
    ];
    for (steps, start) in cases {
        let mut args = vec![input.as_str()];
        args.extend(steps);
        let (code, stderr) = chain(&args);
        assert_eq!(code, Some(1), "{args:?}: {stderr}");
        assert!(stderr.starts_with(start), "{stderr}");
        assert!(stderr.contains("No such file or directory"), "{stderr}");
        assert_eq!(stderr.lines().count(), 1, "{stderr}");
        assert!(names(&dir).is_empty(), "{args:?}: {:?}", names(&dir));
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn usage_errors_exit_2_before_anything_runs() {
    let dir = fresh_dir("usage");
    let out = arg("O", &dir.join("out.txt"));
    let sam = arg("O", &dir.join("out.sam"));
    let cram = arg("O", &dir.join("out.cram"));
    let input = format!("I={}", FIXTURE.0);
    let onto_input = format!("O={}", FIXTURE.0);
    let cases: [(&[&str], &str); 11] = [
        (
            &[&input, "STEP=write", &sam],
            "out.sam: bammetrics writes BAM only",
        ),
        (
            &[&input, "STEP=MarkDuplicates", "STEP=write", &cram],
            "out.cram: bammetrics writes BAM only",
        ),
        // Refused from the command line: the first MarkDuplicates reads
        // nothing for a chain that cannot finish.
        (
            &[
                &input,
                "STEP=MarkDuplicates",
                "STEP=CollectQualityYieldMetrics",
                &out,
                "STEP=MarkDuplicates",
            ],
            "MarkDuplicates#2 reads its input twice and cannot follow MarkDuplicates#1",
        ),
        (&[&input], "no STEP= given"),
        (
            &[&input, "STEP=MarkDuplicates", "STEP=write", &onto_input],
            "is an input; writing it would destroy it",
        ),
        (&["STEP=write", &out], "INPUT is required"),
        (&[&input, "STEP=NoSuchTool"], "unknown step \"NoSuchTool\""),
        (
            &[&input, "STEP=CollectInsertSizeMetrics"],
            "a tap needs OUTPUT",
        ),
        (
            &[&input, "STEP=CollectWgsMetrics", &out],
            "does not read its input in one pass",
        ),
        (&[&input, "STEP=write"], "write: OUTPUT is required"),
        (
            &[&input, "STEP=CollectInsertSizeMetrics", &out, "NO_SUCH=1"],
            "option NO_SUCH is not implemented",
        ),
    ];
    for (args, text) in cases {
        let mut argv = vec!["__chain".to_owned()];
        argv.extend(args.iter().map(|a| (*a).to_owned()));
        let output = bammetrics(&argv);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{argv:?}: {stderr}");
        assert!(stderr.starts_with("bammetrics __chain: "), "{stderr}");
        assert!(stderr.contains(text), "{argv:?}: {stderr}");
    }
    assert!(names(&dir).is_empty());
    // The entry is hidden: no listing names it.
    for listing in ["tools", "help", "--help"] {
        let output = bammetrics(&[listing.to_owned()]);
        assert!(!String::from_utf8_lossy(&output.stdout).contains("__chain"));
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

//! The format surface end to end (docs/formats.md): `convert` through
//! standard input and output gives every Picard fixture back byte for
//! byte; a `--format` list writes its siblings with one start time; `qc
//! --also` lists its files in the manifest and a later run without the
//! option removes them; `check -o` writes each format of its list; and
//! MultiQC, when installed, reads `PREFIX.verdict_mqc.json`. The
//! byte-exact output of each form is in the golden suite.

use std::fs;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

use bammetrics_formats::json::Json;

const BAM: &str = "tests/data/insert_size.bam";
const PICARD: &str = "tests/data/picard";
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

/// `bammetrics ARGS` with `input` on standard input.
fn piped(args: &[&str], input: &[u8]) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run bammetrics");
    child
        .stdin
        .take()
        .expect("a piped stdin")
        .write_all(input)
        .expect("write stdin");
    child.wait_with_output().expect("wait for bammetrics")
}

fn code(output: &Output) -> i32 {
    output.status.code().expect("an exit code")
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-formats-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temporary directory");
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

fn picard_files() -> Vec<PathBuf> {
    let mut files: Vec<PathBuf> = fs::read_dir(root().join(PICARD).join("NA12878"))
        .expect("the Picard fixture directory")
        .map(|e| e.expect("entry").path())
        .collect();
    files.sort();
    files
}

/// `convert FILE | convert -` is `FILE`, for every metrics file Picard
/// 3.5.0 wrote for the fixture BAM.
#[test]
fn convert_twice_gives_every_picard_file_back() {
    let mut checked = 0;
    for path in picard_files() {
        let original = fs::read(&path).expect("read fixture");
        let text = String::from_utf8_lossy(&original);
        if !text.contains("## METRICS CLASS\t") && !text.contains("## HISTOGRAM\t") {
            continue;
        }
        let json = bammetrics(&["convert", path.to_str().expect("utf-8")]);
        assert_eq!(code(&json), 0, "{}: {}", path.display(), stderr(&json));
        assert!(
            json.stdout
                .starts_with(b"{\n  \"schema\": \"bammetrics.metrics/1\"")
        );
        let back = piped(&["convert", "-"], &json.stdout);
        assert_eq!(code(&back), 0, "{}", stderr(&back));
        assert_eq!(back.stdout, original, "{}", path.display());
        // Picard text in and out is the same file too.
        let same = piped(&["convert", "-", "--format", "picard"], &original);
        assert_eq!(same.stdout, original, "{}", path.display());
        checked += 1;
    }
    assert!(checked >= 10, "{checked} files");
}

/// Every format of a list is written, the first to `-o`, and the JSON
/// sibling converts back to the Picard file beside it.
#[test]
fn a_format_list_writes_siblings_that_agree() {
    let dir = fresh_dir("list");
    let out = dir.join("insert_size.txt");
    let run = bammetrics(&[
        "insert-size",
        BAM,
        "--format",
        "picard,json,tsv,csv,jsonl,markdown",
        "-o",
        out.to_str().expect("utf-8"),
        "--verbose",
    ]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert_eq!(
        names(&dir),
        [
            "insert_size.txt",
            "insert_size.txt.csv",
            "insert_size.txt.json",
            "insert_size.txt.jsonl",
            "insert_size.txt.md",
            "insert_size.txt.tsv",
        ]
    );
    for name in names(&dir) {
        assert!(stderr(&run).contains(&format!("wrote {}", dir.join(&name).display())));
    }
    let json = dir.join("insert_size.txt.json");
    let back = bammetrics(&["convert", json.to_str().expect("utf-8")]);
    assert_eq!(back.stdout, fs::read(&out).expect("the Picard file"));
    let csv = fs::read_to_string(dir.join("insert_size.txt.csv")).expect("csv");
    assert!(csv.starts_with("MEDIAN_INSERT_SIZE,MODE_INSERT_SIZE,"));
    assert!(csv.contains(",8.16496580927726,5,FR,"), "{csv}");
    let tsv = fs::read_to_string(dir.join("insert_size.txt.tsv")).expect("tsv");
    assert!(tsv.contains("\t8.164966\t5\tFR\t"), "{tsv}");
    let jsonl = fs::read_to_string(dir.join("insert_size.txt.jsonl")).expect("jsonl");
    assert_eq!(jsonl.lines().count(), 4);
    assert_eq!(
        jsonl
            .lines()
            .filter(|l| l.starts_with("{\"histogram\":true,"))
            .count(),
        3
    );
    let _ = fs::remove_dir_all(&dir);
}

fn manifest_outputs(dir: &Path, step: &str) -> Vec<(String, String)> {
    let text = fs::read_to_string(dir.join("s.manifest.json")).expect("manifest");
    let manifest = Json::parse(&text).expect("manifest parses");
    manifest
        .get("steps")
        .and_then(Json::as_array)
        .expect("steps")
        .iter()
        .filter(|s| s.get("id").and_then(Json::as_str) == Some(step))
        .flat_map(|s| s.get("outputs").and_then(Json::as_array).expect("outputs"))
        .map(|o| {
            let text = |key: &str| o.get(key).and_then(Json::as_str).expect(key).to_owned();
            let name = Path::new(&text("path"))
                .file_name()
                .expect("a file name")
                .to_string_lossy()
                .into_owned();
            (name, text("kind"))
        })
        .collect()
}

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

/// `--also` files are listed under the manifest step `also`, rewritten by
/// `--resume`, and removed by a run that no longer asks for them.
#[test]
fn qc_also_lists_its_files_and_a_later_run_removes_them() {
    let dir = fresh_dir("also");
    let first = qc(
        &dir,
        &["--also", "json,csv", "--also", "markdown,junit,multiqc"],
    );
    assert_eq!(code(&first), 0, "{}", stderr(&first));
    let metrics = [
        "s.alignment_summary_metrics",
        "s.duplicate_metrics",
        "s.insert_size_metrics",
        "s.quality_yield_metrics",
    ];
    let mut expected: Vec<(String, String)> = Vec::new();
    for kind in ["json", "csv"] {
        for name in metrics {
            expected.push((format!("{name}.{kind}"), kind.to_owned()));
        }
    }
    for (name, kind) in [
        ("s.report.md", "markdown"),
        ("s.junit.xml", "junit"),
        ("s.verdict_mqc.json", "multiqc"),
    ] {
        expected.push((name.to_owned(), kind.to_owned()));
    }
    assert_eq!(manifest_outputs(&dir, "also"), expected);
    for (name, _) in &expected {
        assert!(dir.join(name).is_file(), "{name}");
    }
    // Each JSON sibling is its Picard file.
    for name in metrics {
        let json = dir.join(format!("{name}.json"));
        let back = bammetrics(&["convert", json.to_str().expect("utf-8")]);
        assert_eq!(
            back.stdout,
            fs::read(dir.join(name)).expect("metrics"),
            "{name}"
        );
    }
    let markdown = fs::read_to_string(dir.join("s.report.md")).expect("markdown");
    assert!(markdown.starts_with("# insert_size: PASS\n"));
    assert!(markdown.contains("\n## picard.analysis.InsertSizeMetrics\n"));
    assert!(markdown.chars().count() <= 65_536);

    // --resume keeps the collectors' files and writes the also files again.
    fs::remove_file(dir.join("s.junit.xml")).expect("remove one");
    let resumed = qc(&dir, &["--resume", "--also", "junit"]);
    assert_eq!(code(&resumed), 0, "{}", stderr(&resumed));
    assert_eq!(
        manifest_outputs(&dir, "also"),
        [("s.junit.xml".to_owned(), "junit".to_owned())]
    );
    assert!(dir.join("s.junit.xml").is_file());
    assert!(!dir.join("s.report.md").exists());
    assert!(!dir.join("s.insert_size_metrics.json").exists());

    // A file no run of this prefix wrote stays, whatever its name.
    fs::write(dir.join("s.insert_size_metrics.csv"), "mine\n").expect("write");
    let plain = qc(&dir, &["--force"]);
    assert_eq!(code(&plain), 0, "{}", stderr(&plain));
    assert!(manifest_outputs(&dir, "also").is_empty());
    assert!(!dir.join("s.junit.xml").exists());
    assert_eq!(
        fs::read_to_string(dir.join("s.insert_size_metrics.csv")).expect("kept"),
        "mine\n"
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn check_writes_every_format_of_its_list_beside_the_output() {
    let base = fresh_dir("check");
    // A missing output directory is created.
    let dir = base.join("ci");
    let out = dir.join("check.xml");
    let run = bammetrics(&[
        "check",
        PICARD,
        "--format",
        "junit,text,md,multiqc",
        "-o",
        out.to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    assert!(run.stdout.is_empty());
    assert_eq!(
        names(&dir),
        [
            "check.xml",
            "check.xml.md",
            "check.xml.txt",
            "check.xml.verdict_mqc.json"
        ]
    );
    let junit = fs::read_to_string(&out).expect("junit");
    assert!(junit.starts_with("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<testsuites "));
    let printed = bammetrics(&["check", PICARD]);
    assert_eq!(
        fs::read(dir.join("check.xml.txt")).expect("text"),
        printed.stdout
    );
    // The gate decides the exit code, whatever the format.
    let gated = bammetrics(&[
        "check",
        PICARD,
        "--preset",
        "clinical-panel",
        "--format",
        "junit",
    ]);
    assert_eq!(code(&gated), 10);
    let xml = String::from_utf8_lossy(&gated.stdout).into_owned();
    assert!(xml.contains("<failure type=\"fail\""), "{xml}");
    let _ = fs::remove_dir_all(&base);
}

/// `tests/data/no_reads.bam` is the header of the fixture BAM and no read
/// (`samtools view -b -H --no-PG tests/data/insert_size.bam`).
/// `CollectDuplicateMetrics` writes a file with headers only for it, and
/// `qc` exits 0. `--also` exits 0 too, in every format, and each metrics
/// file has its sibling.
#[test]
fn qc_also_on_a_bam_with_no_read_exits_as_qc_does() {
    let dir = fresh_dir("no-reads");
    let prefix = dir.join("s");
    let run = |extra: &[&str]| {
        let mut args = vec![
            "qc",
            "tests/data/no_reads.bam",
            "-o",
            prefix.to_str().expect("utf-8"),
            "--date",
            DATE,
            "--force",
        ];
        args.extend_from_slice(extra);
        bammetrics(&args)
    };
    let plain = run(&[]);
    assert_eq!(code(&plain), 0, "{}", stderr(&plain));
    let duplicate = dir.join("s.duplicate_metrics");
    let headers_only = fs::read_to_string(&duplicate).expect("duplicate metrics");
    assert!(!headers_only.contains("## METRICS CLASS"), "{headers_only}");
    for format in [
        "json", "jsonl", "tsv", "csv", "markdown", "junit", "multiqc",
    ] {
        let also = run(&["--also", format]);
        assert_eq!(code(&also), 0, "--also {format}: {}", stderr(&also));
    }
    let all = run(&["--also", "json,jsonl,tsv,csv,markdown,junit,multiqc"]);
    assert_eq!(code(&all), 0, "{}", stderr(&all));
    let metrics: Vec<String> = manifest_outputs(&dir, "also")
        .into_iter()
        .filter(|(_, kind)| kind == "json")
        .map(|(name, _)| name)
        .collect();
    assert!(metrics.contains(&"s.duplicate_metrics.json".to_owned()));
    for name in &metrics {
        let json = dir.join(name);
        let back = bammetrics(&["convert", json.to_str().expect("utf-8")]);
        assert_eq!(code(&back), 0, "{name}: {}", stderr(&back));
        let picard = dir.join(name.trim_end_matches(".json"));
        assert_eq!(back.stdout, fs::read(picard).expect("metrics"), "{name}");
    }
    for extension in ["jsonl", "tsv", "csv"] {
        let sibling = dir.join(format!("s.duplicate_metrics.{extension}"));
        assert_eq!(fs::read(sibling).expect(extension), b"", "{extension}");
    }
    let _ = fs::remove_dir_all(&dir);
}

/// A JSON cell with a tab or a line end has no Picard text form: `convert
/// --format picard` exits 1, names the cell and writes no file. The forms
/// with an escape take it.
#[test]
fn convert_refuses_picard_text_for_a_cell_it_cannot_hold() {
    let dir = fresh_dir("unwritable");
    let input = dir.join("in.json");
    fs::write(
        &input,
        r#"{"schema": "bammetrics.metrics/1", "headers": ["x"],
            "class": "picard.sam.DuplicationMetrics",
            "columns": ["LIBRARY", "READ_PAIRS_EXAMINED"], "types": ["string", "int"],
            "metrics": [{"LIBRARY": "ok", "READ_PAIRS_EXAMINED": 1},
                        {"LIBRARY": "a\tb\nc", "READ_PAIRS_EXAMINED": 2}],
            "non_finite": [], "histogram": null}"#,
    )
    .expect("write input");
    let input = input.to_str().expect("utf-8");
    let out = dir.join("out.tsv");
    let refused = bammetrics(&[
        "convert",
        input,
        "--format",
        "tsv,picard",
        "-o",
        out.to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&refused), 1);
    assert!(refused.stdout.is_empty());
    assert!(
        stderr(&refused).contains(
            "Picard text cannot hold the tab in metrics row 1, column LIBRARY: \"a\\tb\\nc\""
        ),
        "{}",
        stderr(&refused)
    );
    assert_eq!(names(&dir), ["in.json"]);
    let tsv = bammetrics(&["convert", input, "--format", "tsv"]);
    assert_eq!(code(&tsv), 0, "{}", stderr(&tsv));
    assert_eq!(
        String::from_utf8_lossy(&tsv.stdout),
        "LIBRARY\tREAD_PAIRS_EXAMINED\nok\t1\na\\tb\\nc\t2\n"
    );
    let _ = fs::remove_dir_all(&dir);
}

/// A byte order mark before the JSON is dropped; text that starts as a
/// JSON array is reported as JSON.
#[test]
fn convert_reads_json_after_a_byte_order_mark() {
    let picard = root()
        .join(PICARD)
        .join("NA12878/NA12878.insert_size_metrics");
    let original = fs::read(&picard).expect("read fixture");
    let json = bammetrics(&["convert", picard.to_str().expect("utf-8")]);
    let mut marked = b"\xef\xbb\xbf".to_vec();
    marked.extend_from_slice(&json.stdout);
    let back = piped(&["convert", "-"], &marked);
    assert_eq!(code(&back), 0, "{}", stderr(&back));
    assert_eq!(back.stdout, original);
    let array = piped(&["convert", "-"], "[".repeat(200_000).as_bytes());
    assert_eq!(code(&array), 1);
    assert!(
        stderr(&array).starts_with("bammetrics convert: standard input: json: "),
        "{}",
        stderr(&array)
    );
}

/// In a list whose first format is not `picard`, the Picard sibling is the
/// file `-o` of that name alone writes: its header names it as `OUTPUT`.
#[test]
fn the_picard_sibling_names_itself_as_output() {
    let dir = fresh_dir("sibling");
    let json = dir.join("a.json");
    let listed = bammetrics(&[
        "insert-size",
        BAM,
        "--format",
        "json,picard",
        "-o",
        json.to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&listed), 0, "{}", stderr(&listed));
    assert_eq!(names(&dir), ["a.json", "a.json.txt"]);
    let sibling = dir.join("a.json.txt");
    let command = |path: &Path| {
        fs::read_to_string(path)
            .expect("a Picard file")
            .lines()
            .nth(1)
            .expect("the command line")
            .to_owned()
    };
    let expected = format!("OUTPUT={}", sibling.display());
    assert!(
        command(&sibling).ends_with(&expected),
        "{}",
        command(&sibling)
    );
    let alone = dir.join("alone.txt");
    let plain = bammetrics(&["insert-size", BAM, "-o", alone.to_str().expect("utf-8")]);
    assert_eq!(code(&plain), 0, "{}", stderr(&plain));
    assert_eq!(
        command(&alone),
        command(&sibling).replace("a.json.txt", "alone.txt")
    );
    // The JSON beside it converts to the sibling byte for byte.
    let back = bammetrics(&["convert", json.to_str().expect("utf-8")]);
    assert_eq!(back.stdout, fs::read(&sibling).expect("the sibling"));
    let _ = fs::remove_dir_all(&dir);
}

/// `check -o` refuses a metrics file or a verdict it reads, though the
/// input argument is its directory; a new file in that directory is fine.
#[test]
fn check_refuses_to_overwrite_a_file_it_reads() {
    let dir = fresh_dir("check-input");
    let source = root().join(PICARD).join("NA12878");
    for entry in fs::read_dir(&source).expect("fixture directory") {
        let path = entry.expect("entry").path();
        fs::copy(&path, dir.join(path.file_name().expect("a name"))).expect("copy");
    }
    let dir_text = dir.to_str().expect("utf-8");
    let metrics = dir.join("NA12878.insert_size_metrics");
    let before = fs::read(&metrics).expect("metrics");
    for format in ["text", "junit,text"] {
        let run = bammetrics(&[
            "check",
            dir_text,
            "--format",
            format,
            "-o",
            metrics.to_str().expect("utf-8"),
        ]);
        assert_eq!(code(&run), 2, "{format}: {}", stderr(&run));
        assert!(
            stderr(&run).contains("is an input; writing it would destroy it"),
            "{}",
            stderr(&run)
        );
        assert_eq!(fs::read(&metrics).expect("metrics"), before);
    }
    // A sibling of the list that is an input is refused too.
    let hs = dir.join("NA12878.hs_metrics");
    fs::copy(&hs, dir.join("check.txt")).expect("an input named as the sibling");
    let sibling = bammetrics(&[
        "check",
        dir_text,
        "--format",
        "junit,text",
        "-o",
        dir.join("check").to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&sibling), 2, "{}", stderr(&sibling));
    fs::remove_file(dir.join("check.txt")).expect("remove");
    let fine = bammetrics(&[
        "check",
        dir_text,
        "-o",
        dir.join("check.txt").to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&fine), 0, "{}", stderr(&fine));

    // A directory of verdicts: each *.verdict.json is an input.
    let verdicts = fresh_dir("check-verdicts");
    let report = bammetrics(&["report", dir_text, "-o", verdicts.to_str().expect("utf-8")]);
    assert_eq!(code(&report), 0, "{}", stderr(&report));
    let verdict = verdicts.join("NA12878.verdict.json");
    let before = fs::read(&verdict).expect("verdict");
    let run = bammetrics(&[
        "check",
        verdicts.to_str().expect("utf-8"),
        "--format",
        "json",
        "-o",
        verdict.to_str().expect("utf-8"),
    ]);
    assert_eq!(code(&run), 2, "{}", stderr(&run));
    assert_eq!(fs::read(&verdict).expect("verdict"), before);
    let _ = fs::remove_dir_all(&dir);
    let _ = fs::remove_dir_all(&verdicts);
}

/// A 0600 output keeps its mode when it is replaced, under any umask: the
/// main file and the siblings of `convert`, of a tool verb and of `check`,
/// and every file `qc --also` rewrites.
#[cfg(unix)]
#[test]
fn a_replaced_output_keeps_its_mode() {
    use std::os::unix::fs::PermissionsExt;
    let dir = fresh_dir("mode");
    let private = |path: &Path| {
        if !path.exists() {
            fs::write(path, "old").expect("an earlier file");
        }
        fs::set_permissions(path, fs::Permissions::from_mode(0o600)).expect("chmod");
    };
    let mode = |path: &Path| fs::metadata(path).expect("metadata").permissions().mode() & 0o777;
    let picard = root()
        .join(PICARD)
        .join("NA12878/NA12878.insert_size_metrics");

    let converted = dir.join("conv.json");
    let verb = dir.join("verb.json");
    let listed = dir.join("list.txt");
    let checked = dir.join("check.txt");
    let files = [
        converted.clone(),
        dir.join("conv.json.csv"),
        verb.clone(),
        listed.clone(),
        dir.join("list.txt.json"),
        checked.clone(),
        dir.join("check.txt.junit.xml"),
    ];
    for file in &files {
        private(file);
    }
    let text = |path: &Path| path.to_str().expect("utf-8").to_owned();
    let runs = [
        vec![
            "convert".to_owned(),
            text(&picard),
            "--format".to_owned(),
            "json,csv".to_owned(),
            "-o".to_owned(),
            text(&converted),
        ],
        vec![
            "insert-size".to_owned(),
            BAM.to_owned(),
            "--format".to_owned(),
            "json".to_owned(),
            "-o".to_owned(),
            text(&verb),
        ],
        vec![
            "insert-size".to_owned(),
            BAM.to_owned(),
            "--format".to_owned(),
            "picard,json".to_owned(),
            "-o".to_owned(),
            text(&listed),
        ],
        vec![
            "check".to_owned(),
            PICARD.to_owned(),
            "--format".to_owned(),
            "text,junit".to_owned(),
            "-o".to_owned(),
            text(&checked),
        ],
    ];
    for argv in &runs {
        let args: Vec<&str> = argv.iter().map(String::as_str).collect();
        let run = bammetrics(&args);
        assert_eq!(code(&run), 0, "{argv:?}: {}", stderr(&run));
    }
    for file in &files {
        assert_ne!(
            fs::read(file).expect("written"),
            b"old",
            "{}",
            file.display()
        );
        assert_eq!(mode(file), 0o600, "{}", file.display());
    }

    let run_dir = dir.join("qc");
    fs::create_dir(&run_dir).expect("run directory");
    let first = qc(&run_dir, &["--also", "json,junit"]);
    assert_eq!(code(&first), 0, "{}", stderr(&first));
    let written = names(&run_dir);
    for name in &written {
        private(&run_dir.join(name));
    }
    let again = qc(&run_dir, &["--also", "json,junit", "--force"]);
    assert_eq!(code(&again), 0, "{}", stderr(&again));
    assert_eq!(names(&run_dir), written, "no temporary file stays");
    for name in &written {
        assert_eq!(mode(&run_dir.join(name)), 0o600, "{name}");
    }
    let _ = fs::remove_dir_all(&dir);
}

/// The `multiqc` on PATH, or the one `MULTIQC` names.
fn multiqc() -> Option<PathBuf> {
    let program =
        std::env::var_os("MULTIQC").map_or_else(|| PathBuf::from("multiqc"), PathBuf::from);
    Command::new(&program)
        .arg("--version")
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|_| program)
}

/// MultiQC reads the custom-content file and the Picard files of a `qc`
/// run: its data directory holds the verdict table with the sample's row.
#[test]
fn multiqc_reads_the_verdict_file() {
    let Some(multiqc) = multiqc() else {
        eprintln!("skipped: multiqc is not installed (set MULTIQC=/path/to/multiqc)");
        return;
    };
    let dir = fresh_dir("multiqc");
    let run = qc(&dir, &["--also", "multiqc", "--sample", "NA12878"]);
    assert_eq!(code(&run), 0, "{}", stderr(&run));
    let report = dir.join("mqc");
    let output = Command::new(multiqc)
        .arg(&dir)
        .arg("--outdir")
        .arg(&report)
        .args(["--force", "--no-ansi", "--quiet"])
        .output()
        .expect("run multiqc");
    assert!(output.status.success(), "{}", stderr(&output));
    let data = report.join("multiqc_data");
    let table = names(&data)
        .into_iter()
        .find(|n| n.contains("bammetrics_verdict") && n.ends_with(".txt"))
        .unwrap_or_else(|| panic!("no bammetrics_verdict table in {:?}", names(&data)));
    let text = fs::read_to_string(data.join(table)).expect("the table");
    let mut lines = text.lines();
    let header = lines.next().expect("a header line");
    assert!(header.contains("verdict"), "{header}");
    assert!(header.contains("insert.median"), "{header}");
    let row = lines.next().expect("a sample row");
    assert!(row.starts_with("NA12878\t"), "{row}");
    assert!(row.contains("pass"), "{row}");
    assert!(row.contains("310"), "{row}");
    let _ = fs::remove_dir_all(&dir);
}

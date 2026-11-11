//! `BuildBamIndex` and `build-bam-index` end to end, without Picard: the
//! index of each coordinate-sorted fixture BAM is the file Picard 3.5.0
//! wrote for it (`tests/data/bai`, `tests/data/make_bai.py`), under both
//! faces and every spelling; the default name is Picard's; a refused run
//! writes nothing; and the `.bai` that `MarkDuplicates CREATE_INDEX=true`
//! and `chain ... + write --index` put beside a BAM is the index of that
//! BAM's bytes.

use std::fs::{self, DirBuilder, OpenOptions};
use std::io::Write;
use std::os::unix::fs::{DirBuilderExt, OpenOptionsExt, PermissionsExt};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

/// The fixture BAMs whose header says `SO:coordinate`.
const SORTED: [&str; 6] = [
    "artifact",
    "insert_size",
    "markdup",
    "no_reads",
    "orientation",
    "oxog",
];
const QUERYNAME: &str = "tests/data/duplicate_marking.bam";
const CRAM: &str = "tests/data/insert_size.cram";
const MARKDUP: &str = "tests/data/markdup.bam";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn command(args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bammetrics"));
    command.args(args).current_dir(root());
    command
}

fn bammetrics(args: &[&str]) -> Output {
    command(args).output().expect("run bammetrics")
}

fn ok(args: &[&str]) -> Output {
    let run = bammetrics(args);
    assert_eq!(run.status.code(), Some(0), "{args:?}: {run:?}");
    run
}

/// A directory of this test's own under the temp dir: mode 0700, created
/// without recursion so an existing name fails.
fn fresh_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-build-index-{}-{name}-{nanos}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("a new temp dir");
    dir
}

/// `bytes` to a new file at `path` with mode 0600; an existing path fails.
fn write_new(path: &Path, bytes: &[u8]) {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .mode(0o600)
        .open(path)
        .and_then(|mut file| file.write_all(bytes))
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
}

fn fixture(name: &str) -> String {
    format!("tests/data/{name}.bam")
}

fn picard_index(name: &str) -> Vec<u8> {
    fs::read(root().join(format!("tests/data/bai/{name}.bai"))).expect("Picard's index")
}

fn text(path: &Path) -> &str {
    path.to_str().expect("utf-8")
}

/// The names under `dir`, sorted.
fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("list")
        .map(|entry| {
            entry
                .expect("entry")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    names.sort();
    names
}

#[test]
fn every_fixture_index_is_the_file_picard_wrote() {
    let dir = fresh_dir("fixtures");
    for name in SORTED {
        let expected = picard_index(name);
        let bam = fixture(name);
        let picard_face = dir.join(format!("{name}.picard.bai"));
        let run = ok(&[
            "BuildBamIndex",
            &format!("I={bam}"),
            &format!("O={}", picard_face.display()),
        ]);
        assert_eq!(run.stdout, b"", "{name}: the Picard face prints nothing");
        assert_eq!(fs::read(&picard_face).expect("index"), expected, "{name}");

        let verb = dir.join(format!("{name}.verb.bai"));
        ok(&["build-bam-index", &bam, "-o", text(&verb)]);
        assert_eq!(fs::read(&verb).expect("index"), expected, "{name}");

        let printed = ok(&["build-bam-index", &bam, "-o", "-"]);
        assert_eq!(printed.stdout, expected, "{name}: -o - prints the index");
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn every_spelling_writes_the_same_index() {
    let dir = fresh_dir("spellings");
    let expected = picard_index("markdup");
    let out = dir.join("x.bai");
    let o = format!("O={}", out.display());
    let i = format!("I={MARKDUP}");
    let forms: [Vec<&str>; 8] = [
        vec!["build-bam-index", MARKDUP, "-o", text(&out)],
        vec!["build-bam-index", "-i", MARKDUP, "--output", text(&out)],
        vec![
            "build-bam-index",
            "--INPUT",
            MARKDUP,
            "--OUTPUT",
            text(&out),
        ],
        vec!["build-bam-index", &i, &o, "VALIDATION_STRINGENCY=SILENT"],
        vec![
            "build-bam-index",
            MARKDUP,
            "-O",
            text(&out),
            "--validation-stringency",
            "SILENT",
            "--create-index",
            "--threads",
            "1",
        ],
        vec!["BuildBamIndex", &i, &o],
        vec!["BuildBamIndex", "--INPUT", MARKDUP, "-O", text(&out)],
        vec![
            "BuildBamIndex",
            "-I",
            MARKDUP,
            &o,
            "CREATE_MD5_FILE=true",
            "CREATE_INDEX=true",
            "R=tests/data/ref_100k.fa",
            "TMP_DIR=/nowhere",
            "COMPRESSION_LEVEL=1",
        ],
    ];
    for form in forms {
        ok(&form);
        assert_eq!(fs::read(&out).expect("index"), expected, "{form:?}");
        // CREATE_MD5_FILE writes no digest of an index, as in Picard.
        assert_eq!(listing(&dir), ["x.bai"], "{form:?}");
        fs::remove_file(&out).expect("remove the index");
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn the_default_name_is_beside_the_bam_as_picard_names_it() {
    let dir = fresh_dir("default");
    let bytes = fs::read(root().join(MARKDUP)).expect("fixture BAM");
    let expected = picard_index("markdup");
    for (bam, index) in [
        ("x.bam", "x.bai"),
        ("a.b.bam", "a.b.bai"),
        ("plain", "plain.bai"),
    ] {
        for face in ["BuildBamIndex", "build-bam-index"] {
            let path = dir.join(bam);
            write_new(&path, &bytes);
            let input = match face {
                "BuildBamIndex" => format!("I={}", path.display()),
                _ => text(&path).to_owned(),
            };
            let run = ok(&[face, &input]);
            assert_eq!(run.stdout, b"", "{face} {bam}");
            let mut names = vec![bam.to_owned(), index.to_owned()];
            names.sort();
            assert_eq!(listing(&dir), names, "{face} {bam}");
            assert_eq!(fs::read(dir.join(index)).expect("index"), expected);
            fs::remove_file(dir.join(index)).expect("remove the index");
            fs::remove_file(&path).expect("remove the BAM");
        }
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn standard_input_is_indexed_and_needs_an_output_name() {
    let dir = fresh_dir("stdin");
    let out = dir.join("x.bai");
    let piped = |args: &[&str]| {
        let mut child = command(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("spawn");
        let bytes = fs::read(root().join(MARKDUP)).expect("fixture BAM");
        child
            .stdin
            .take()
            .expect("stdin")
            .write_all(&bytes)
            .expect("write the BAM");
        child.wait_with_output().expect("wait")
    };
    let run = piped(&["build-bam-index", "-", "-o", text(&out)]);
    assert_eq!(run.status.code(), Some(0), "{run:?}");
    assert_eq!(fs::read(&out).expect("index"), picard_index("markdup"));
    fs::remove_file(&out).expect("remove the index");

    let run = piped(&["build-bam-index", "-"]);
    assert_eq!(run.status.code(), Some(2), "{run:?}");
    assert_eq!(
        String::from_utf8_lossy(&run.stderr),
        "bammetrics build-bam-index: missing -o FILE: standard input has no name to put the \
         BAM index beside\n  -o - writes the BAM index to standard output\n"
    );
    let run = piped(&["BuildBamIndex", "I=/dev/stdin"]);
    assert_eq!(run.status.code(), Some(1), "{run:?}");
    assert!(
        String::from_utf8_lossy(&run.stderr).contains("give OUTPUT"),
        "{run:?}"
    );
    assert!(listing(&dir).is_empty());
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn a_refused_run_writes_nothing_and_keeps_an_earlier_index() {
    let dir = fresh_dir("refused");
    let out = dir.join("x.bai");
    let o = format!("O={}", out.display());
    write_new(&out, b"an earlier index");
    let cases: [(&[&str], i32, &str); 8] = [
        (
            &["BuildBamIndex", &format!("I={QUERYNAME}"), &o],
            1,
            "Input bam file must be sorted by coordinate (its header says the sort order is \
             queryname)",
        ),
        (
            &["build-bam-index", QUERYNAME, "-o", text(&out)],
            1,
            "Input bam file must be sorted by coordinate",
        ),
        (
            &["build-bam-index", CRAM, "-o", text(&out)],
            2,
            "a .bai indexes a BAM and this is a CRAM",
        ),
        (
            &[
                "build-bam-index",
                "tests/data/ref_100k.fa",
                "-o",
                text(&out),
            ],
            1,
            "not a BAM or CRAM file",
        ),
        (
            &[
                "build-bam-index",
                MARKDUP,
                "-o",
                text(&out),
                "--format",
                "json",
            ],
            2,
            "--format json has nothing to write",
        ),
        (
            &[
                "BuildBamIndex",
                &format!("I={MARKDUP}"),
                &o,
                "ASSUME_SORTED=true",
            ],
            1,
            "Unrecognized option: ASSUME_SORTED",
        ),
        (
            &["BuildBamIndex", &format!("I={MARKDUP}"), &o, "STOP_AFTER=5"],
            1,
            "option STOP_AFTER is not implemented",
        ),
        (
            &[
                "build-bam-index",
                "tests/data/missing.bam",
                "-o",
                text(&out),
            ],
            1,
            "tests/data/missing.bam not found",
        ),
    ];
    for (args, code, message) in cases {
        let run = bammetrics(args);
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert_eq!(run.status.code(), Some(code), "{args:?}: {stderr}");
        assert!(stderr.contains(message), "{args:?}: {stderr}");
        assert_eq!(run.stdout, b"", "{args:?}");
        assert_eq!(listing(&dir), ["x.bai"], "{args:?}");
        assert_eq!(fs::read(&out).expect("index"), b"an earlier index");
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn the_output_is_never_the_input_and_keeps_its_mode() {
    let dir = fresh_dir("files");
    let bam = dir.join("x.bam");
    let bytes = fs::read(root().join(MARKDUP)).expect("fixture BAM");
    write_new(&bam, &bytes);
    let link = dir.join("link.bai");
    fs::hard_link(&bam, &link).expect("hard link");
    for output in [&bam, &link] {
        for args in [
            vec![
                "BuildBamIndex".to_owned(),
                format!("I={}", bam.display()),
                format!("O={}", output.display()),
            ],
            vec![
                "build-bam-index".to_owned(),
                text(&bam).to_owned(),
                "-o".to_owned(),
                text(output).to_owned(),
            ],
        ] {
            let args: Vec<&str> = args.iter().map(String::as_str).collect();
            let run = bammetrics(&args);
            let stderr = String::from_utf8_lossy(&run.stderr);
            assert_ne!(run.status.code(), Some(0), "{args:?}");
            assert!(stderr.contains("would destroy it"), "{args:?}: {stderr}");
            assert_eq!(fs::read(&bam).expect("the BAM"), bytes, "{args:?}");
        }
    }
    fs::remove_file(&link).expect("remove the link");

    // An index replaced by rename keeps the mode the earlier file had.
    let index = dir.join("x.bai");
    write_new(&index, b"old");
    fs::set_permissions(&index, fs::Permissions::from_mode(0o640)).expect("chmod");
    ok(&["build-bam-index", text(&bam)]);
    assert_eq!(fs::read(&index).expect("index"), picard_index("markdup"));
    let mode = fs::metadata(&index).expect("stat").permissions().mode() & 0o777;
    assert_eq!(mode, 0o640);
    assert_eq!(listing(&dir), ["x.bai", "x.bam"]);

    // A directory that does not exist is an error before the BAM is read.
    let missing = dir.join("nowhere/x.bai");
    let run = bammetrics(&[
        "BuildBamIndex",
        &format!("I={}", bam.display()),
        &format!("O={}", missing.display()),
    ]);
    assert_eq!(run.status.code(), Some(1), "{run:?}");
    assert_eq!(listing(&dir), ["x.bai", "x.bam"]);
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

/// The hex MD5 of `bytes`.
fn md5(bytes: &[u8]) -> String {
    let mut md5 = bammetrics_core::md5::Md5::new();
    md5.update(bytes);
    bammetrics_core::md5::hex(&md5.finish())
}

#[test]
fn the_index_beside_a_written_bam_is_the_index_of_its_bytes() {
    let dir = fresh_dir("written");
    let marked = dir.join("marked.bam");
    ok(&[
        "MarkDuplicates",
        &format!("I={MARKDUP}"),
        &format!("O={}", marked.display()),
        &format!("M={}", dir.join("metrics.txt").display()),
        "CREATE_INDEX=true",
        "CREATE_MD5_FILE=true",
    ]);
    let chained = dir.join("chained.bam");
    ok(&[
        "chain",
        MARKDUP,
        "+",
        "write",
        text(&chained),
        "--index",
        "--md5",
        "--compression-level",
        "1",
    ]);
    for bam in [&marked, &chained] {
        let beside = bam.with_extension("bai");
        let again = dir.join("again.bai");
        ok(&["build-bam-index", text(bam), "-o", text(&again)]);
        let index = fs::read(&beside).expect("the index beside the BAM");
        assert_eq!(
            index,
            fs::read(&again).expect("the index built from the BAM")
        );
        assert_eq!(&index[..8], b"BAI\x01\x01\0\0\0");
        let digest = fs::read_to_string(format!("{}.md5", bam.display())).expect("the MD5");
        assert_eq!(digest, md5(&fs::read(bam).expect("the BAM")));
        fs::remove_file(&again).expect("remove the index");
    }
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

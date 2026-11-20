//! The split tools (docs/chaining.md, "Split steps"): a split as the last
//! step of a chain equals the tool run on the file the chain would have
//! written, file for file; the files of a split that fails are never in
//! place; the files of two steps never collide; the verb makes its
//! directory.
//!
//! Staged means each transform as its own tool from file to file, then
//! the split tool on the last file. Chained means one `chain` command.
//! Every file of the two directories must have the same name, the same
//! header lines (the `CL` of `@PG` set aside) and the same records as
//! bytes.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const OXOG: &str = "tests/data/oxog.bam";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";

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

fn run(args: &[&str]) -> Output {
    let output = bammetrics(args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-split-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn text(path: &Path) -> String {
    path.to_str().expect("UTF-8 path").to_owned()
}

/// The names in `dir`, sorted.
fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .map(|entries| {
            entries
                .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

/// The header lines of a BAM with the `CL` field of each `@PG` line
/// removed, and every record as its bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<Vec<u8>>) {
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
            records.push(record.as_bytes().to_vec());
        }
        Ok((lines, records))
    })
    .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Two directories of BAMs hold the same files.
fn same_files(chained: &Path, staged: &Path) -> usize {
    let listed = names(chained);
    assert_eq!(listed, names(staged), "the file names");
    for name in &listed {
        if name.ends_with(".bam") {
            assert_eq!(bam(&chained.join(name)), bam(&staged.join(name)), "{name}");
        }
    }
    listed.len()
}

#[test]
fn a_library_split_after_mark_duplicates_equals_the_staged_run() {
    for source in [MARKDUP, OXOG, FOUNDATION] {
        if !root().join(source).exists() {
            continue;
        }
        let dir = fresh_dir(&format!("library-{}", source.len()));
        let (chained, staged) = (dir.join("chained"), dir.join("staged"));
        let marked = dir.join("marked.bam");
        run(&[
            "chain",
            source,
            "-o",
            &text(&dir.join("c")),
            "+",
            "mark-duplicates",
            "+",
            "split-sam-by-library",
            "-o",
            &text(&chained),
            "--create-index",
            "--create-md5-file",
        ]);
        run(&[
            "mark-duplicates",
            source,
            "-o",
            &text(&marked),
            "--metrics-file",
            &text(&dir.join("m.txt")),
        ]);
        run(&[
            "split-sam-by-library",
            &text(&marked),
            "-o",
            &text(&staged),
            "--create-index",
            "--create-md5-file",
        ]);
        let files = same_files(&chained, &staged);
        assert!(files >= 6, "{source}: {files} files");
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}

#[test]
fn a_count_split_in_a_chain_equals_the_staged_run() {
    let dir = fresh_dir("count");
    let (chained, staged) = (dir.join("chained"), dir.join("staged"));
    run(&[
        "chain",
        MARKDUP,
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "+",
        "split-sam-by-number-of-reads",
        "-o",
        &text(&chained),
        "--n-reads",
        "10",
        "--total-reads",
        "36",
        "--out-prefix",
        "part",
    ]);
    let sorted = dir.join("sorted.bam");
    run(&[
        "sort-sam",
        MARKDUP,
        "-o",
        &text(&sorted),
        "-SO",
        "queryname",
    ]);
    run(&[
        "split-sam-by-number-of-reads",
        &text(&sorted),
        "-o",
        &text(&staged),
        "--n-reads",
        "10",
        "--out-prefix",
        "part",
    ]);
    assert!(same_files(&chained, &staged) >= 3);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A split beside a tap and `write` on one stream: each gets every record.
#[test]
fn a_split_takes_the_stream_beside_a_tap_and_a_writer() {
    let dir = fresh_dir("beside");
    let libs = dir.join("libs");
    run(&[
        "chain",
        MARKDUP,
        "-o",
        &text(&dir.join("s")),
        "+",
        "insert-size",
        "+",
        "split-sam-by-library",
        "-o",
        &text(&libs),
        "+",
        "write",
        &text(&dir.join("all.bam")),
    ]);
    let (_, all) = bam(&dir.join("all.bam"));
    let split: usize = ["lib1.bam", "lib2.bam"]
        .iter()
        .map(|name| bam(&libs.join(name)).1.len())
        .sum();
    assert_eq!(split, all.len());
    assert!(dir.join("s.insert_size_metrics").exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A record out of the header's order fails the split at that record:
/// no file of the split is in place and no temporary file is left.
#[test]
fn a_split_that_fails_leaves_no_file() {
    let dir = fresh_dir("fails");
    let libs = dir.join("libs");
    // The duplicate-marking fixture says SO:queryname; sorted by
    // coordinate its records break that order in every library.
    let output = bammetrics(&[
        "chain",
        MARKDUP,
        "+",
        "replace-sam-header",
        "--header",
        "tests/data/duplicate_marking.bam",
        "+",
        "split-sam-by-library",
        "-o",
        &text(&libs),
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(names(&libs).is_empty(), "{:?}", names(&libs));
    let output = bammetrics(&[
        "SplitSamByLibrary",
        "I=tests/data/insert_size.bam",
        &format!("O={}", text(&dir)),
    ]);
    // A header without libraries is Picard's exit 2 before any file.
    assert_eq!(output.status.code(), Some(2));
    // The failed chain removed the directory it made.
    assert!(names(&dir).is_empty(), "{:?}", names(&dir));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A count above the 100,000-file bound is taken when the records give
/// fewer files, as Picard takes it: one record per file at most.
#[test]
fn a_count_above_the_file_bound_writes_the_files_the_records_give() {
    let dir = fresh_dir("huge-count");
    run(&[
        "SplitSamByNumberOfReads",
        &format!("I={MARKDUP}"),
        &format!("O={}", text(&dir)),
        "SPLIT_TO_N_FILES=2147483647",
    ]);
    let files = names(&dir);
    assert!((2..=36).contains(&files.len()), "{files:?}");
    assert_eq!(files[0], "shard_0001.bam");
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// SIGINT while a split of 20,000 files renames them into place: the
/// files renamed so far stay registered for removal until the last
/// rename, so the run leaves no file, renamed or temporary.
#[cfg(unix)]
#[test]
fn a_signal_during_the_renames_leaves_no_file() {
    use std::fmt::Write as _;
    use std::os::unix::process::ExitStatusExt;
    let dir = fresh_dir("signal");
    let sam = dir.join("unmapped.sam");
    let mut lines = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n");
    for i in 0..20_000 {
        writeln!(lines, "r{i:05}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII").expect("line");
    }
    fs::write(&sam, lines).expect("SAM");
    let out = dir.join("out");
    fs::create_dir(&out).expect("output directory");
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "SplitSamByNumberOfReads",
            &format!("I={}", text(&sam)),
            &format!("O={}", text(&out)),
            "SPLIT_TO_N_FILES=20000",
        ])
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("run bammetrics");
    let renamed = |dir: &Path| {
        fs::read_dir(dir).is_ok_and(|mut entries| {
            entries.any(|e| e.is_ok_and(|e| !e.file_name().to_string_lossy().starts_with('.')))
        })
    };
    while child.try_wait().expect("wait").is_none() && !renamed(&out) {
        std::thread::sleep(std::time::Duration::from_millis(2));
    }
    let signalled = Command::new("kill")
        .args(["-INT", &child.id().to_string()])
        .status()
        .expect("kill runs");
    assert!(signalled.success());
    let status = child.wait().expect("wait");
    assert_eq!(status.signal(), Some(2), "the run ended before the signal");
    assert!(names(&out).is_empty(), "{} files left", names(&out).len());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A `write` whose file is one a split makes once the stream is read:
/// the split finds it when it makes the file and fails.
#[test]
fn a_split_claims_its_files_against_the_other_steps() {
    let dir = fresh_dir("claims");
    let libs = dir.join("libs");
    let output = bammetrics(&[
        "chain",
        MARKDUP,
        "+",
        "write",
        &text(&libs.join("lib1.bam")),
        "+",
        "split-sam-by-library",
        "-o",
        &text(&libs),
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("both write"), "{stderr}");
    assert!(!libs.join("lib2.bam").exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn the_verb_makes_its_directory_and_refuses_a_cram_per_library() {
    let dir = fresh_dir("verb");
    let libs = dir.join("new").join("libs");
    run(&["split-sam-by-library", OXOG, "-o", &text(&libs)]);
    assert_eq!(names(&libs), ["L1.bam", "L2.bam"]);
    let refused = dir.join("cram");
    let output = bammetrics(&[
        "split-sam-by-library",
        "tests/data/insert_size.cram",
        "-R",
        "tests/data/ref_100k.fa",
        "-o",
        &text(&refused),
    ]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("W15"));
    assert!(!refused.exists(), "the directory it made is removed");
    fs::remove_dir_all(&dir).expect("cleanup");
}

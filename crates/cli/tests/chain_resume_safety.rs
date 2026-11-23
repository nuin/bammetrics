//! What the review of wave W14 added to `chain --resume` (docs/chaining.md
//! section 27, "File safety"): a respelled source or prefix resumes, a
//! source swapped during a run reruns the next one, one run per prefix, the
//! files a resumed run removes and the ones it leaves, and the documented
//! examples run as written.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const MARKDUP: &str = "tests/data/markdup.bam";
const INSERT_SIZE: &str = "tests/data/insert_size.bam";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics_in(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(dir)
        .output()
        .expect("run bammetrics")
}

fn run_in(dir: &Path, args: &[&str]) -> String {
    let output = bammetrics_in(dir, args);
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{args:?}: {stderr}");
    stderr
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-resume-safety-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir.canonicalize().expect("scratch dir")
}

/// A metrics file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.split_inclusive('\n')
        .skip_while(|line| line.starts_with('#'))
        .collect()
}

fn all_skipped(stderr: &str, steps: usize) -> bool {
    stderr.contains(&format!(
        "{steps} steps finished, {steps} of them skipped by --resume"
    ))
}

/// `./m.bam`, an absolute path, a link and `../m.bam` from a subdirectory
/// name the source of the run before: every step is skipped.
#[test]
fn a_respelled_source_or_prefix_resumes() {
    let dir = fresh_dir("respelled");
    fs::copy(root().join(MARKDUP), dir.join("m.bam")).expect("copy");
    #[cfg(unix)]
    std::os::unix::fs::symlink("m.bam", dir.join("link.bam")).expect("link");
    fs::create_dir_all(dir.join("sub")).expect("sub");
    let steps = |source: &str, prefix: &str| -> Vec<String> {
        [
            "chain",
            source,
            "-o",
            prefix,
            "--resume",
            "+",
            "quality-yield",
            "+",
            "mark-duplicates",
            "--keep",
            "+",
            "insert-size",
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect()
    };
    let call = |at: &Path, source: &str, prefix: &str| {
        let args = steps(source, prefix);
        run_in(at, &args.iter().map(String::as_str).collect::<Vec<_>>())
    };
    call(&dir, "m.bam", "out/s");
    let absolute = dir.join("m.bam").display().to_string();
    let prefix = dir.join("out/s").display().to_string();
    let mut forms = vec![
        (dir.clone(), "./m.bam".to_owned(), "out/s".to_owned()),
        (dir.clone(), absolute, prefix),
        (
            dir.join("sub"),
            "../m.bam".to_owned(),
            "../out/s".to_owned(),
        ),
        (dir.clone(), "m.bam".to_owned(), "./out/../out/s".to_owned()),
    ];
    if cfg!(unix) {
        forms.push((dir.clone(), "link.bam".to_owned(), "out/s".to_owned()));
    }
    for (at, source, prefix) in forms {
        let stderr = call(&at, &source, &prefix);
        assert!(all_skipped(&stderr, 3), "{source} -o {prefix}: {stderr}");
    }
    // Another file under the name: the note names both reads.
    fs::copy(root().join(INSERT_SIZE), dir.join("other.bam")).expect("copy");
    let args = [
        "chain",
        "other.bam",
        "-o",
        "out/s",
        "--resume",
        "+",
        "quality-yield",
    ];
    let stderr = run_in(&dir, &args);
    assert!(
        stderr.contains("; this run reads other.bam: every step runs"),
        "{stderr}"
    );
    let _ = fs::remove_dir_all(&dir);
}

/// The source is fingerprinted before the read: when it is replaced while
/// the run reads it, the next `--resume` reruns and its files equal a fresh
/// run on the file now there, whatever the moment of the swap.
#[test]
fn a_source_swapped_during_a_run_reruns_the_next_one() {
    let dir = fresh_dir("swap");
    let source = dir.join("src.bam");
    fs::copy(root().join(MARKDUP), &source).expect("copy");
    let args = |prefix: &str| -> Vec<String> {
        [
            "chain",
            "src.bam",
            "-o",
            prefix,
            "--resume",
            "+",
            "quality-yield",
            "+",
            "insert-size",
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect()
    };
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args("o/s"))
        .current_dir(&dir)
        .stderr(std::process::Stdio::null())
        .spawn()
        .expect("spawn");
    fs::copy(root().join(INSERT_SIZE), dir.join("other.bam")).expect("copy");
    fs::rename(dir.join("other.bam"), &source).expect("swap");
    assert!(child.wait().expect("wait").success());
    let resumed = args("o/s");
    run_in(
        &dir,
        &resumed.iter().map(String::as_str).collect::<Vec<_>>(),
    );
    let fresh = args("f/s");
    run_in(&dir, &fresh.iter().map(String::as_str).collect::<Vec<_>>());
    for name in ["s.quality_yield_metrics", "s.insert_size_metrics"] {
        assert_eq!(
            body(&dir.join("o").join(name)),
            body(&dir.join("f").join(name)),
            "{name}"
        );
    }
    let _ = fs::remove_dir_all(&dir);
}

/// While one run holds the prefix's lock, a second exits 2 before it reads.
#[test]
fn a_second_run_on_one_prefix_is_refused() {
    let dir = fresh_dir("lock");
    fs::create_dir_all(dir.join("c")).expect("dir");
    let lock = fs::File::create(dir.join("c/s.manifest.json.lock")).expect("lock file");
    lock.try_lock().expect("held by the test");
    let source = root().join(MARKDUP).display().to_string();
    let args = [
        "chain",
        &source,
        "-o",
        "c/s",
        "--resume",
        "+",
        "insert-size",
    ];
    let output = bammetrics_in(&dir, &args);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(
        stderr.contains("another chain is writing c/s.manifest.json"),
        "{stderr}"
    );
    assert!(!dir.join("c/s.insert_size_metrics").exists());
    drop(lock);
    // A lock file nobody holds is no lock; the run removes it at its end.
    run_in(&dir, &args);
    assert!(!dir.join("c/s.manifest.json.lock").exists());
    let _ = fs::remove_dir_all(&dir);
}

/// A dropped `write` keeps its BAM; a dropped tap's file under the prefix
/// goes; the `.bai` of a kept stream no longer coordinate sorted goes.
#[test]
fn what_a_resumed_run_removes_and_what_it_leaves() {
    let dir = fresh_dir("stale");
    let source = root().join(MARKDUP).display().to_string();
    let first = [
        "chain",
        &source,
        "-o",
        "out/s",
        "--resume",
        "+",
        "sort-sam",
        "-SO",
        "coordinate",
        "--keep",
        "+",
        "quality-yield",
        "+",
        "write",
        "deliver/final.bam",
        "--index",
    ];
    run_in(&dir, &first);
    assert!(dir.join("out/s.sorted.bai").is_file());
    let second = [
        "chain",
        &source,
        "-o",
        "out/s",
        "--resume",
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "--keep",
    ];
    let stderr = run_in(&dir, &second);
    for gone in ["out/s.sorted.bai", "out/s.quality_yield_metrics"] {
        assert!(!dir.join(gone).exists(), "{gone}: {stderr}");
        assert!(
            stderr.contains(&format!("--resume: removed {gone}:")),
            "{stderr}"
        );
    }
    for kept in ["deliver/final.bam", "deliver/final.bai"] {
        assert!(dir.join(kept).is_file(), "{kept}");
        assert!(
            stderr.contains(&format!("--resume: left in place: {kept} (a write")),
            "{stderr}"
        );
    }
    assert!(dir.join("out/s.sorted.bam").is_file());
    let _ = fs::remove_dir_all(&dir);
}

/// A record grafted into the manifest that names a file outside the
/// prefix's files is never removed, and the note names it.
#[test]
fn a_hand_edited_manifest_removes_nothing_outside_the_prefix() {
    let dir = fresh_dir("graft");
    let source = root().join(MARKDUP).display().to_string();
    run_in(
        &dir,
        &[
            "chain",
            &source,
            "-o",
            "v/x",
            "--resume",
            "+",
            "insert-size",
        ],
    );
    let args = [
        "chain",
        &source,
        "-o",
        "o/s",
        "--resume",
        "+",
        "mark-duplicates",
        "--keep",
        "+",
        "insert-size",
    ];
    run_in(&dir, &args);
    // The record of v/x.insert_size_metrics, under a path through `..`.
    let theirs = fs::read_to_string(dir.join("v/x.manifest.json")).expect("manifest");
    let start = theirs.find("\"outputs\": [").expect("outputs") + "\"outputs\": [".len();
    let end = start + theirs[start..].find(']').expect("end");
    let record = theirs[start..end].trim().replace(
        "\"v/x.insert_size_metrics\"",
        "\"o/../v/x.insert_size_metrics\"",
    );
    let ours = fs::read_to_string(dir.join("o/s.manifest.json")).expect("manifest");
    let at = ours.rfind("\"outputs\": [").expect("outputs") + "\"outputs\": [".len();
    let grafted = format!("{}\n{record},{}", &ours[..at], &ours[at..]);
    fs::write(dir.join("o/s.manifest.json"), grafted).expect("graft");
    let changed = [
        "chain",
        &source,
        "-o",
        "o/s",
        "--resume",
        "+",
        "mark-duplicates",
        "--keep",
        "+",
        "insert-size",
        "--minimum-pct",
        "0.2",
    ];
    let stderr = run_in(&dir, &changed);
    assert!(dir.join("v/x.insert_size_metrics").is_file(), "{stderr}");
    assert!(
        stderr.contains("left in place: o/../v/x.insert_size_metrics (outside the files"),
        "{stderr}"
    );
    let _ = fs::remove_dir_all(&dir);
}

/// A temporary file a killed run left beside the prefix is named on the
/// next `--resume`, and left in place.
#[test]
fn a_killed_run_s_temporary_file_is_named() {
    let dir = fresh_dir("orphan");
    fs::create_dir_all(dir.join("k")).expect("dir");
    // No process has an id this large on Linux or macOS.
    let orphan = "k/.s.marked.bam.999999999.0.0123456789abcdef.tmp";
    fs::write(dir.join(orphan), "x").expect("orphan");
    let source = root().join(MARKDUP).display().to_string();
    let stderr = run_in(
        &dir,
        &[
            "chain",
            &source,
            "-o",
            "k/s",
            "--resume",
            "+",
            "insert-size",
        ],
    );
    assert!(
        stderr.contains(&format!(
            "--resume: left in place: {orphan} (a temporary file of a run that was killed"
        )),
        "{stderr}"
    );
    assert!(dir.join(orphan).is_file());
    let _ = fs::remove_dir_all(&dir);
}

/// A directory that takes no new file fails before the read, naming it.
#[cfg(unix)]
#[test]
fn an_unwritable_prefix_directory_fails_before_the_read() {
    use std::os::unix::fs::PermissionsExt;
    let dir = fresh_dir("readonly");
    fs::create_dir_all(dir.join("ro")).expect("dir");
    fs::set_permissions(dir.join("ro"), fs::Permissions::from_mode(0o555)).expect("chmod");
    let source = root().join(MARKDUP).display().to_string();
    let output = bammetrics_in(
        &dir,
        &[
            "chain",
            &source,
            "-o",
            "ro/s",
            "--resume",
            "+",
            "insert-size",
        ],
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    fs::set_permissions(dir.join("ro"), fs::Permissions::from_mode(0o755)).expect("chmod");
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains("the manifest ro/s.manifest.json: cannot write in ro:"),
        "{stderr}"
    );
    assert!(!stderr.contains(".tmp"), "{stderr}");
    let _ = fs::remove_dir_all(&dir);
}

/// The chain examples of docs/cli.md and the README that use `--resume`,
/// `--keep` and `qc`, as written there, on the fixture BAM.
#[test]
fn the_documented_resume_examples_run() {
    let cli = fs::read_to_string(root().join("docs/cli.md")).expect("docs/cli.md");
    let readme = fs::read_to_string(root().join("README.md")).expect("README.md");
    let examples: [(&str, &str, &[&str]); 2] = [
        (
            &cli,
            "bammetrics chain in.bam -o out/NA12878 --resume + mark-duplicates --keep \\\n    + \
             insert-size + write out/NA12878.bam + qc --preset generic --fail-on fail",
            &[
                "-o",
                "out/NA12878",
                "--resume",
                "+",
                "mark-duplicates",
                "--keep",
                "+",
                "insert-size",
                "+",
                "write",
                "out/NA12878.bam",
                "+",
                "qc",
                "--preset",
                "generic",
                "--fail-on",
                "fail",
            ],
        ),
        (
            &readme,
            "bammetrics chain sample.bam -o out/NA12878 --resume \\\n    + mark-duplicates --keep \
             + sort-sam --sort-order queryname + insert-size \\\n    + qc --preset generic \
             --fail-on fail",
            &[
                "-o",
                "out/NA12878",
                "--resume",
                "+",
                "mark-duplicates",
                "--keep",
                "+",
                "sort-sam",
                "--sort-order",
                "queryname",
                "+",
                "insert-size",
                "+",
                "qc",
                "--preset",
                "generic",
                "--fail-on",
                "fail",
            ],
        ),
    ];
    let source = root().join(MARKDUP).display().to_string();
    for (doc, text, args) in examples {
        assert!(doc.contains(text), "the docs no longer show {text}");
        let dir = fresh_dir("docs");
        let mut argv = vec!["chain", source.as_str()];
        argv.extend_from_slice(args);
        for _ in 0..2 {
            let output = bammetrics_in(&dir, &argv);
            let stderr = String::from_utf8_lossy(&output.stderr);
            // 0, or 10 and 11 of the qc gate: never a usage error.
            assert!(
                matches!(output.status.code(), Some(0 | 10 | 11)),
                "{argv:?}: {stderr}"
            );
            assert!(dir.join("out/NA12878.verdict.json").is_file());
        }
        let _ = fs::remove_dir_all(&dir);
    }
}

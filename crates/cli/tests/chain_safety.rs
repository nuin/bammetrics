//! What `bammetrics chain` refuses before its first read, from the 0.22.0
//! release review (docs/chaining.md section 2.5):
//!
//! * an output that is a symbolic link another user owns, with the refusal
//!   and the exit code of the tool verb;
//! * an output spelled through a directory that does not exist and a `..`
//!   behind it, which would name an input once the directory is made;
//! * a plan file over 1 MiB;
//! * a thread count over 256, on every command that takes one;
//! * a text value only the tool can judge (`--level NOPE`).
//!
//! And what a failed run leaves: no directory it made and left empty.
//!
//! "Before the first read" is shown with a source cut after 100 bytes: a
//! run that read it would fail with a read error in place of the refusal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "tests/data/insert_size.bam";
const MARKDUP: &str = "tests/data/markdup.bam";
/// Queryname sorted: a tap with `--no-assume-sorted` fails when it is built.
const QUERYNAME: &str = "tests/data/duplicate_marking.bam";

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
        .stdin(std::process::Stdio::null())
        .output()
        .expect("run bammetrics")
}

/// The exit code and standard error of `bammetrics ARGS`.
fn run(args: &[&str]) -> (Option<i32>, String) {
    let output = bammetrics(args);
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-chain-safety-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    // `/var` is a link on macOS: the messages print the path as given.
    dir.canonicalize().expect("scratch dir")
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("directory")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

fn text(path: &Path) -> String {
    path.display().to_string()
}

/// The first 100 bytes of the fixture BAM under `dir`: its magic is whole,
/// its first record is not.
fn cut_bam(dir: &Path) -> String {
    let bytes = fs::read(root().join(BAM)).expect("fixture");
    let path = dir.join("cut.bam");
    fs::write(&path, &bytes[..100]).expect("cut");
    text(&path)
}

/// A symbolic link on this machine that another user owns (`/etc` on
/// macOS, `/bin` on a merged-usr Linux).
#[cfg(unix)]
fn foreign_link() -> Option<&'static str> {
    use std::os::unix::fs::MetadataExt;
    // This user owns the scratch directory it has just made.
    let scratch = fresh_dir("uid");
    let me = fs::metadata(&scratch).map(|m| m.uid());
    fs::remove_dir(&scratch).expect("cleanup");
    let me = me.ok()?;
    ["/etc", "/tmp", "/var", "/bin", "/sbin", "/lib"]
        .into_iter()
        .find(|path| {
            fs::symlink_metadata(path).is_ok_and(|m| m.file_type().is_symlink() && m.uid() != me)
        })
}

/// A chain output that is a link another user owns gets the refusal the
/// tool verb gives, and its exit code, before the first read. The test
/// needs such a link on the machine and passes without one.
#[cfg(unix)]
#[test]
fn a_link_of_another_user_is_refused_as_the_verb_refuses_it() {
    let Some(link) = foreign_link() else {
        eprintln!("skipped: no symbolic link of another user on this machine");
        return;
    };
    let dir = fresh_dir("link");
    let refusal =
        format!("{link}: {link} is a symbolic link owned by another user; not written in place");
    let (verb_code, verb) = run(&["insert-size", BAM, "-o", link]);
    assert_eq!(
        verb,
        format!("bammetrics insert-size: {refusal} (from -o)\n")
    );
    let cut = cut_bam(&dir);
    let good = text(&dir.join("made/good.txt"));
    let cases: [(&[&str], &str); 4] = [
        (&["+", "insert-size", "-o", link], "step 1 (insert-size)"),
        (&["+", "write", link], "step 1 (write)"),
        (
            &["+", "mark-duplicates", "--metrics-file", link],
            "step 1 (mark-duplicates)",
        ),
        // The directory of the first step is made, then removed again.
        (
            &[
                "+",
                "quality-yield",
                "-o",
                &good,
                "+",
                "insert-size",
                "-o",
                link,
            ],
            "step 2 (insert-size)",
        ),
    ];
    for (steps, label) in cases {
        let mut args = vec!["chain", cut.as_str()];
        args.extend(steps);
        let (code, stderr) = run(&args);
        assert_eq!(stderr, format!("bammetrics chain: {label}: {refusal}\n"));
        assert_eq!(code, verb_code, "{args:?}");
        assert_eq!(code, Some(1));
        assert_eq!(names(&dir), ["cut.bam"], "{args:?}");
    }
    // The Picard face fails with the same words.
    let (code, stderr) = run(&[
        "CollectQualityYieldMetrics",
        &format!("I={BAM}"),
        &format!("O={link}"),
    ]);
    assert_eq!(code, Some(1));
    assert_eq!(stderr, format!("bammetrics: {refusal}\n"));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// `missing/../in.bam` names no file while `missing` does not exist, and
/// names the input once the run has made it. It is a usage error before
/// the first read, on a chain and on a tool verb: no directory is made and
/// the input keeps its bytes.
#[test]
fn an_output_through_a_missing_directory_and_back_is_refused() {
    let dir = fresh_dir("unresolved");
    let source = dir.join("in.bam");
    let bytes = fs::read(root().join(BAM)).expect("fixture");
    fs::write(&source, &bytes).expect("copy");
    let source = text(&source);
    let onto_input = text(&dir.join("missing/../in.bam"));
    let prefix = text(&dir.join("missing/../s"));
    let through = text(&dir.join("missing/../a.txt"));
    let plain = text(&dir.join("a.txt"));
    let behind = text(&dir.join("made/deeper/../../a.txt"));
    let missing = text(&dir.join("missing"));
    let hint = "\n  name the output without the ..: a run makes the directories that are missing\n";
    let wrong = |path: &str, missing: &str| {
        format!("{path}: {missing} does not exist, so the .. behind it names no directory")
    };
    let chains: [(&[&str], String); 5] = [
        (
            &["+", "write", &onto_input],
            format!("step 1 (write): {}", wrong(&onto_input, &missing)),
        ),
        (
            &["+", "insert-size", "-o", &onto_input],
            format!("step 1 (insert-size): {}", wrong(&onto_input, &missing)),
        ),
        (
            &["-o", &prefix, "--also", "markdown", "+", "quality-yield"],
            format!(
                "step 1 (quality-yield): {}",
                wrong(&format!("{prefix}.quality_yield_metrics"), &missing)
            ),
        ),
        (
            &[
                "+",
                "quality-yield",
                "-o",
                &plain,
                "+",
                "insert-size",
                "-o",
                &through,
            ],
            format!("step 2 (insert-size): {}", wrong(&through, &missing)),
        ),
        (
            &[
                "+",
                "sequencing-artifact",
                "-R",
                "tests/data/ref_100k.fa",
                "-o",
                &behind,
            ],
            format!(
                "step 1 (sequencing-artifact): {}",
                wrong(&behind, &text(&dir.join("made")))
            ),
        ),
    ];
    for (steps, what) in chains {
        let mut args = vec!["chain", source.as_str()];
        args.extend(steps);
        let (code, stderr) = run(&args);
        assert_eq!(
            stderr,
            format!("bammetrics chain: {what}{hint}"),
            "{args:?}"
        );
        assert_eq!(code, Some(2), "{args:?}");
        assert_eq!(names(&dir), ["in.bam"], "{args:?}");
    }
    let (code, stderr) = run(&["insert-size", &source, "-o", &onto_input]);
    assert_eq!(
        stderr,
        format!(
            "bammetrics insert-size: {} (from -o){hint}",
            wrong(&onto_input, &missing)
        )
    );
    assert_eq!(code, Some(2));
    assert_eq!(names(&dir), ["in.bam"]);
    assert_eq!(fs::read(dir.join("in.bam")).expect("input"), bytes);
    // A `..` behind a directory that exists resolves, and is checked as
    // any spelling of the file it names.
    fs::create_dir(dir.join("there")).expect("dir");
    let resolved = text(&dir.join("there/../in.bam"));
    let (code, stderr) = run(&["chain", &source, "+", "write", &resolved]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.contains("is an input; writing it would destroy it"),
        "{stderr}"
    );
    let fine = text(&dir.join("there/../qy.txt"));
    let (code, stderr) = run(&["chain", &source, "+", "quality-yield", "-o", &fine]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(names(&dir), ["in.bam", "qy.txt", "there"]);
    assert_eq!(fs::read(dir.join("in.bam")).expect("input"), bytes);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Two spellings of one file under a directory that does not exist yet are
/// one file: the deepest directory that exists is resolved.
#[cfg(unix)]
#[test]
fn one_new_file_through_a_link_and_without_it_is_claimed_once() {
    let dir = fresh_dir("alias");
    std::os::unix::fs::symlink(".", dir.join("self")).expect("symlink");
    let linked = text(&dir.join("self/new/a.txt"));
    let plain = text(&dir.join("new/a.txt"));
    let (code, stderr) = run(&[
        "chain",
        BAM,
        "+",
        "insert-size",
        "-o",
        &linked,
        "+",
        "quality-yield",
        "-o",
        &plain,
    ]);
    assert_eq!(code, Some(2), "{stderr}");
    assert!(
        stderr.starts_with(&format!(
            "bammetrics chain: step 1 (insert-size) and step 2 (quality-yield) both write {plain}: \
             every output of a chain needs its own file"
        )),
        "{stderr}"
    );
    assert_eq!(names(&dir), ["self"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A run that fails removes the directories it made that are still empty.
/// A directory with a file in it stays, and so does one that was there.
#[test]
fn a_failed_run_leaves_no_directory_it_made_empty() {
    let dir = fresh_dir("dirs");
    let blocker = dir.join("blocker");
    fs::write(&blocker, "a file, not a directory").expect("blocker");
    fs::create_dir(dir.join("old")).expect("old");
    let kept = ["blocker", "old"];
    // An output that cannot be written: found after the directories of the
    // steps in front were made, before the first read.
    let (code, stderr) = run(&[
        "chain",
        &cut_bam(&dir),
        "-o",
        &text(&dir.join("made/deeper/s")),
        "--also",
        "markdown",
        "+",
        "quality-yield",
        "+",
        "write",
        &text(&blocker.join("x.bam")),
    ]);
    fs::remove_file(dir.join("cut.bam")).expect("cut");
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.starts_with("bammetrics chain: step 2 (write): "),
        "{stderr}"
    );
    assert_eq!(names(&dir), kept);
    // Every step fails during the run: a transform on a BAM without a read
    // group stops its segment.
    let (code, stderr) = run(&[
        "chain",
        BAM,
        "-o",
        &text(&dir.join("gone/deeper/s")),
        "+",
        "mark-duplicates",
        "+",
        "insert-size",
        "+",
        "write",
        &text(&dir.join("other/x.bam")),
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stderr.contains("0 of 3 steps finished"), "{stderr}");
    assert_eq!(names(&dir), kept);
    // A directory that was there stays, empty as it was.
    let (code, stderr) = run(&[
        "chain",
        QUERYNAME,
        "-o",
        &text(&dir.join("old/s")),
        "+",
        "quality-yield",
        "--no-assume-sorted",
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stderr.contains("failed alone"), "{stderr}");
    assert_eq!(names(&dir), kept);
    assert!(names(&dir.join("old")).is_empty());
    // One step fails alone and one finishes: the directory with a file
    // stays, the empty one goes.
    let (code, stderr) = run(&[
        "chain",
        QUERYNAME,
        "+",
        "quality-yield",
        "--no-assume-sorted",
        "-o",
        &text(&dir.join("empty/qy.txt")),
        "+",
        "mean-quality-by-cycle",
        "-o",
        &text(&dir.join("full/deeper/mq.txt")),
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(stderr.contains("1 of 2 steps finished"), "{stderr}");
    assert_eq!(names(&dir), ["blocker", "full", "old"]);
    assert_eq!(names(&dir.join("full/deeper")), ["mq.txt"]);
    // A run that finishes keeps what it made.
    let (code, stderr) = run(&[
        "chain",
        MARKDUP,
        "-o",
        &text(&dir.join("done/deeper/s")),
        "+",
        "mark-duplicates",
        "+",
        "quality-yield",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(
        names(&dir.join("done/deeper")),
        ["s.marked.duplication_metrics", "s.quality_yield_metrics"]
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A plan file holds at most 1 MiB: a larger file is a usage error with
/// the cap in it, and at most one byte past the cap is read.
#[test]
fn a_plan_file_over_one_mebibyte_is_refused() {
    let dir = fresh_dir("plan");
    let plan = format!(
        "[chain]\nsource = \"{BAM}\"\n\n[[step]]\nverb = \"quality-yield\"\noutput = \"{}\"\n# ",
        text(&dir.join("qy.txt"))
    );
    let padded = |bytes: usize| format!("{plan}{}\n", "x".repeat(bytes - plan.len() - 1));
    let path = dir.join("plan.toml");
    let name = text(&path);
    fs::write(&path, padded((1 << 20) + 1)).expect("write");
    for verb in ["chain", "plan"] {
        let (code, stderr) = run(&[verb, "--plan", &name]);
        assert_eq!(
            stderr,
            format!(
                "bammetrics chain: {name}: not a plan file: it is larger than 1 MiB (from --plan)\n  \
                 a plan file is TOML text (docs/chaining.md)\n"
            )
        );
        assert_eq!(code, Some(2));
    }
    assert_eq!(names(&dir), ["plan.toml"]);
    fs::write(&path, padded(1 << 20)).expect("write");
    let (code, stderr) = run(&["chain", "--plan", &name]);
    assert_eq!(code, Some(0), "{stderr}");
    assert_eq!(names(&dir), ["plan.toml", "qy.txt"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A thread count over 256 is a usage error where it is parsed, on every
/// command that takes one. Before, 100000000 was a panic in the BGZF
/// reader (exit 101) after two minutes of starting threads.
#[test]
fn a_thread_count_over_the_maximum_is_a_usage_error_everywhere() {
    let dir = fresh_dir("threads");
    // Under directories a run would make: a usage error makes none.
    let out = text(&dir.join("made/out.txt"));
    let prefix = text(&dir.join("made/deeper/s"));
    let plan = dir.join("plan.toml");
    for threads in ["257", "100000000"] {
        fs::write(
            &plan,
            format!("[chain]\nsource = \"{BAM}\"\nthreads = {threads}\n\n[[step]]\nverb = \"quality-yield\"\n"),
        )
        .expect("plan");
        let over = format!("{threads} is more than the maximum of 256 threads");
        let clap = format!("invalid value \"{threads}\": more than the maximum of 256 threads");
        let picard_threads = format!("threads={threads}");
        let cases: [(&[&str], String, i32); 9] = [
            (
                &["insert-size", BAM, "-o", &out, "--threads", threads],
                format!("bammetrics insert-size: {over} (from --threads)\n"),
                2,
            ),
            (
                &[
                    "mark-duplicates",
                    MARKDUP,
                    "-o",
                    &out,
                    "--metrics-file",
                    &prefix,
                    "-@",
                    threads,
                ],
                format!("bammetrics mark-duplicates: {over} (from -@)\n"),
                2,
            ),
            (
                &[
                    "wgs",
                    BAM,
                    "-R",
                    "tests/data/ref_100k.fa",
                    "-o",
                    &out,
                    "--threads",
                    threads,
                ],
                format!("bammetrics wgs: {over} (from --threads)\n"),
                2,
            ),
            (
                &[
                    "chain",
                    BAM,
                    "--threads",
                    threads,
                    "-o",
                    &prefix,
                    "+",
                    "quality-yield",
                ],
                format!("bammetrics chain: {over} (from --threads)\n"),
                2,
            ),
            (
                &["chain", "--plan", &text(&plan), "-o", &prefix],
                format!("bammetrics chain: {}:3: threads: {over}\n", text(&plan)),
                2,
            ),
            (
                &["qc", BAM, "-o", &prefix, "--threads", threads],
                format!(
                    "bammetrics qc: {clap} (from --threads)\n  every option: bammetrics help qc\n"
                ),
                2,
            ),
            (
                &["plan", BAM, "-o", &prefix, "--threads", threads],
                format!(
                    "bammetrics plan: {clap} (from --threads)\n  every option: bammetrics help plan\n"
                ),
                2,
            ),
            (
                &["collect", BAM, "-o", &prefix, "--threads", threads],
                format!("bammetrics collect: {over} (from --threads)\n"),
                2,
            ),
            // The Picard face exits 1 for every error, as Picard does.
            (
                &[
                    "CollectInsertSizeMetrics",
                    "I=tests/data/insert_size.bam",
                    &picard_threads,
                ],
                format!("bammetrics: threads={threads}: more than the maximum of 256 threads\n"),
                1,
            ),
        ];
        for (args, stderr, code) in cases {
            let (got_code, got) = run(args);
            assert_eq!(got, stderr, "{args:?}");
            assert_eq!(got_code, Some(code), "{args:?}");
        }
        assert_eq!(names(&dir), ["plan.toml"], "threads {threads}");
    }
    // The maximum itself runs.
    let (code, stderr) = run(&[
        "chain",
        BAM,
        "-@",
        "256",
        "-o",
        &prefix,
        "+",
        "quality-yield",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    let (code, stderr) = run(&["insert-size", BAM, "-o", &out, "--threads", "256"]);
    assert_eq!(code, Some(0), "{stderr}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A text value only the tool can judge is a usage error before the first
/// read, with the words and the exit code of the tool verb: on the source,
/// behind a transform, and in the plan.
#[test]
fn a_value_only_the_tool_can_judge_exits_2_before_the_first_read() {
    let dir = fresh_dir("values");
    let cut = cut_bam(&dir);
    let prefix = text(&dir.join("made/s"));
    let cases: [(&str, &[&str], &str); 4] = [
        (
            "insert-size",
            &["--level", "NOPE"],
            "METRIC_ACCUMULATION_LEVEL \"NOPE\": expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP",
        ),
        (
            "alignment-summary",
            &["--expected-pair-orientations", "UP"],
            "pair orientation \"UP\": expected FR, RF or TANDEM",
        ),
        (
            "gc-bias",
            &["-R", "tests/data/ref_100k.fa", "LEVEL=NOPE"],
            "METRIC_ACCUMULATION_LEVEL \"NOPE\": expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP",
        ),
        (
            "sex-check",
            &["--expected-sex", "NOPE"],
            "EXPECTED_SEX=NOPE: expected MALE or FEMALE",
        ),
    ];
    for (verb, options, what) in cases {
        let mut args = vec![verb, BAM];
        args.extend(options);
        let (code, of_verb) = run(&args);
        assert_eq!(code, Some(2), "{args:?}: {of_verb}");
        // A verb that writes several files says first where they go, and
        // one names the flag the value came from.
        let line = format!("bammetrics {verb}: {what}");
        assert!(of_verb.lines().any(|l| l.starts_with(&line)), "{of_verb}");
        for front in [&[][..], &["+", "mark-duplicates"]] {
            for command in ["chain", "plan"] {
                let mut args = vec![command, cut.as_str(), "-o", prefix.as_str()];
                args.extend(front);
                args.extend(["+", verb]);
                args.extend(options);
                let (code, stderr) = run(&args);
                let number = front.len() / 2 + 1;
                let line = format!("bammetrics chain: step {number} ({verb}): {what}");
                assert!(stderr.starts_with(&line), "{args:?}: {stderr}");
                assert_eq!(stderr.lines().count(), 1, "{args:?}: {stderr}");
                assert_eq!(code, Some(2), "{args:?}");
            }
        }
    }
    assert_eq!(names(&dir), ["cut.bam"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

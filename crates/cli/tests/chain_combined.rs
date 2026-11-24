//! The pieces of 0.32.0 together (docs/chaining.md section 32): the simple
//! grammar (steps without `+`, `-o` naming the BAM, the config file, the
//! recipes) with the manifest (`--resume`, `--keep`, `qc` as the last step)
//! and the materialized steps.
//!
//! Each case is checked against the long form it stands for, or against a
//! fresh run of the changed chain: the same files, apart from the path of
//! the output directory, the clock line of a metrics file and the manifest.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::SystemTime;

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;

const IS: &str = "tests/data/insert_size.bam";
const MD: &str = "tests/data/markdup.bam";
const OXOG: &str = "tests/data/oxog.bam";
const REFERENCE: &str = "tests/data/ref_100k.fa";
const DATE: &str = "2026-01-01 00:00:00 UTC";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

/// A fresh directory under the system temporary directory.
fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-chain-combined-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("directory");
    dir
}

/// A config directory (`XDG_CONFIG_HOME`), with `text` as its config file.
fn config_dir(name: &str, text: Option<&str>) -> PathBuf {
    let dir = fresh_dir(name);
    fs::create_dir_all(dir.join("bammetrics/recipes")).expect("config");
    if let Some(text) = text {
        fs::write(dir.join("bammetrics/config.toml"), text).expect("config file");
    }
    dir
}

/// The binary from the repository root with `config` as its config
/// directory and no `BAMMETRICS_REFERENCE`.
fn run(args: &[&str], config: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .env("XDG_CONFIG_HOME", config)
        .env_remove("BAMMETRICS_REFERENCE")
        .output()
        .expect("run bammetrics")
}

fn stderr_of(output: &Output, code: i32) -> String {
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(code), "{stderr}");
    stderr
}

fn ok(output: &Output) -> String {
    stderr_of(output, 0)
}

fn usage_error(output: &Output) -> String {
    stderr_of(output, 2)
}

/// The files of `dir`, sorted by name.
fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("read dir")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

fn modified(path: &Path) -> SystemTime {
    fs::metadata(path)
        .and_then(|m| m.modified())
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// The header lines of a BAM with the `CL` field of each `@PG` line
/// removed, and every record as its bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<RawRecord>) {
    pass::with_record_stream(path, None, 1, |header, reader, _| {
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
    .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// A text file with the path of its directory as `DIR`, without the lines
/// that hold a command line or the clock (a metrics file's `#` lines, the
/// command line a report page prints).
fn text(path: &Path, dir: &Path) -> String {
    let bytes = fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    String::from_utf8_lossy(&bytes)
        .replace(&dir.display().to_string(), "DIR")
        .lines()
        .filter(|line| !line.starts_with('#') && !line.contains("INPUT="))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Every file of `a` and `b` but the manifest: the same names, each BAM the
/// same header lines and records, each index the same bytes, each text the
/// same lines.
fn same_files(a: &Path, b: &Path) {
    assert_eq!(names(a), names(b));
    for name in names(a) {
        let (left, right) = (a.join(&name), b.join(&name));
        if name.ends_with(".manifest.json") || left.is_dir() {
            continue;
        }
        if name.ends_with(".bam") {
            assert!(bam(&left) == bam(&right), "{name} differs");
        } else if name.ends_with(".bai") {
            // An index holds the offsets of its BAM: the same bytes when
            // the BAMs are (a `@PG` line with a path of the directory
            // moves every offset).
            let of = |dir: &Path| fs::read(dir.join(&name).with_extension("bam")).expect("bam");
            if of(a) == of(b) {
                assert!(
                    fs::read(&left).expect("a") == fs::read(&right).expect("b"),
                    "{name} differs"
                );
            }
        } else {
            assert_eq!(text(&left, a), text(&right, b), "{name}");
        }
    }
}

/// The arguments of a command built as owned text.
fn refs(args: &[String]) -> Vec<&str> {
    args.iter().map(String::as_str).collect()
}

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_owned()).collect()
}

/// The owner's example with `--resume`: `chain BAM reorder-sam insert-size
/// -o r.bam --resume`, the reference from the config file. The second run
/// reads nothing. A changed option of `insert-size` reruns that step alone,
/// reading `r.bam` in place of the source, and leaves `r.bam` untouched.
#[test]
fn the_owners_example_resumes() {
    let reference = root().join(REFERENCE);
    let config = config_dir(
        "owner-config",
        Some(&format!("reference = \"{}\"\n", reference.display())),
    );
    let dir = fresh_dir("owner");
    let out = dir.join("r.bam").display().to_string();
    let command = |step: &[&str]| {
        let mut args = vec!["chain", IS, "reorder-sam", "insert-size"];
        args.extend(step);
        args.extend(["-o", &out, "--resume"]);
        strings(&args)
    };
    let first = ok(&run(&refs(&command(&[])), &config));
    assert!(first.contains("--resume: no manifest at"), "{first}");
    assert_eq!(
        names(&dir),
        ["r.bai", "r.bam", "r.insert_size_metrics", "r.manifest.json"]
    );
    let manifest = fs::read_to_string(dir.join("r.manifest.json")).expect("manifest");
    assert!(manifest.contains("\"schema\": \"bammetrics.manifest/2\""));
    let files = ["r.bam", "r.bai", "r.insert_size_metrics"];
    let before: Vec<SystemTime> = files.iter().map(|f| modified(&dir.join(f))).collect();
    // The same command: every step is skipped and nothing is read.
    let again = ok(&run(&refs(&command(&[])), &config));
    assert!(
        again.contains("3 steps finished, 3 of them skipped by --resume"),
        "{again}"
    );
    for (file, time) in files.iter().zip(&before) {
        assert_eq!(modified(&dir.join(file)), *time, "{file}");
    }
    // An option of insert-size: the plan says where the read starts.
    let mut plan = command(&["--minimum-pct", "0.2"]);
    plan.push("--dry-run".to_owned());
    let planned = run(&refs(&plan), &config);
    ok(&planned);
    let planned = String::from_utf8_lossy(&planned.stdout).into_owned();
    assert!(
        planned.contains("step 1 reorder-sam: skipped; the BAM a write step wrote of its stream"),
        "{planned}"
    );
    assert!(
        planned.contains("step 2 insert-size: runs: its options changed"),
        "{planned}"
    );
    std::thread::sleep(std::time::Duration::from_millis(20));
    let changed = ok(&run(&refs(&command(&["--minimum-pct", "0.2"])), &config));
    for skipped in ["step 1 (reorder-sam): skipped", "step 3 (write): skipped"] {
        assert!(changed.contains(skipped), "{changed}");
    }
    assert!(
        changed.contains("2 of them skipped by --resume"),
        "{changed}"
    );
    for (file, time) in files.iter().zip(&before) {
        let rewritten = modified(&dir.join(file)) != *time;
        assert_eq!(rewritten, *file == "r.insert_size_metrics", "{file}");
    }
    // A fresh run of the changed chain in the long form writes the same files.
    let fresh = fresh_dir("owner-fresh");
    ok(&run(
        &[
            "chain",
            IS,
            "-R",
            REFERENCE,
            "-o",
            &fresh.join("r").display().to_string(),
            "--resume",
            "+",
            "reorder-sam",
            "--sequence-dictionary",
            &format!("{REFERENCE}.fai"),
            "+",
            "insert-size",
            "--minimum-pct",
            "0.2",
            "+",
            "write",
            &fresh.join("r.bam").display().to_string(),
            "--index",
        ],
        &config,
    ));
    same_files(&dir, &fresh);
    // r.bam changed by another hand: it starts no read. The write step runs
    // again and reorder-sam feeds it from the source.
    fs::write(dir.join("r.bam"), fs::read(root().join(MD)).expect("a BAM")).expect("replace");
    let mut plan = command(&["--minimum-pct", "0.3"]);
    plan.push("--dry-run".to_owned());
    let planned = run(&refs(&plan), &config);
    ok(&planned);
    let planned = String::from_utf8_lossy(&planned.stdout).into_owned();
    assert!(
        planned.contains("step 1 reorder-sam: runs: its stream feeds step 2"),
        "{planned}"
    );
    assert!(planned.contains("r.bam is missing or changed"), "{planned}");
    ok(&run(&refs(&command(&["--minimum-pct", "0.2"])), &config));
    same_files(&dir, &fresh);
}

/// A recipe resumed after an option change: `prep BAM -o p.bam --resume`,
/// then the same with `insert-size.minimum-pct=0.2`. Seven of the eight
/// steps are skipped, `insert-size` reads `p.bam`, and the files equal a
/// fresh run of the changed recipe. An option of a transform reruns every
/// step behind it.
#[test]
fn a_recipe_resumes_after_an_option_change() {
    let config = config_dir("recipe-config", None);
    let dir = fresh_dir("recipe");
    let out = dir.join("p.bam").display().to_string();
    ok(&run(&["prep", OXOG, "-o", &out, "--resume"], &config));
    let written = [
        "p.alignment_summary_metrics",
        "p.bai",
        "p.bam",
        "p.duplicate_metrics",
        "p.insert_size_metrics",
        "p.manifest.json",
        "p.marked.duplication_metrics",
        "p.quality_yield_metrics",
    ];
    assert_eq!(names(&dir), written);
    let before: Vec<SystemTime> = written.iter().map(|f| modified(&dir.join(f))).collect();
    std::thread::sleep(std::time::Duration::from_millis(20));
    let change = "insert-size.minimum-pct=0.2";
    let changed = ok(&run(
        &["prep", OXOG, "-o", &out, "--resume", change],
        &config,
    ));
    assert!(
        changed.contains("8 steps finished, 7 of them skipped by --resume"),
        "{changed}"
    );
    for (step, verb) in [
        (1, "fix-mate-information"),
        (2, "mark-duplicates"),
        (3, "sort-sam"),
        (4, "quality-yield"),
        (5, "alignment-summary"),
        (7, "duplicate-metrics"),
        (8, "write"),
    ] {
        assert!(
            changed.contains(&format!("step {step} ({verb}): skipped")),
            "{verb}: {changed}"
        );
    }
    for (file, time) in written.iter().zip(&before) {
        let rewritten = modified(&dir.join(file)) != *time;
        let expected = matches!(*file, "p.insert_size_metrics" | "p.manifest.json");
        assert_eq!(rewritten, expected, "{file}");
    }
    let fresh = fresh_dir("recipe-fresh");
    let fresh_out = fresh.join("p.bam").display().to_string();
    ok(&run(
        &["prep", OXOG, "-o", &fresh_out, "--resume", change],
        &config,
    ));
    same_files(&dir, &fresh);
    // An option of a transform: no step behind it is skipped.
    let dedup = "mark-duplicates.remove-duplicates=true";
    let rerun = ok(&run(
        &["prep", OXOG, "-o", &out, "--resume", change, dedup],
        &config,
    ));
    assert!(!rerun.contains("skipped"), "{rerun}");
    assert!(rerun.contains("8 steps finished"), "{rerun}");
    ok(&run(
        &["prep", OXOG, "-o", &fresh_out, "--resume", change, dedup],
        &config,
    ));
    same_files(&dir, &fresh);
}

/// `qc` as the last step of the prep recipe: `prep BAM -o p.bam qc` writes
/// what the same chain by hand writes, and the verdict `qc` gives for the
/// BAM the recipe wrote. Any other step beside a recipe is refused, and so
/// is a second qc step.
#[test]
fn qc_is_the_last_step_of_the_prep_recipe() {
    let config = config_dir("prep-qc-config", None);
    let (recipe, hand) = (fresh_dir("prep-qc"), fresh_dir("prep-qc-hand"));
    let out = recipe.join("p.bam").display().to_string();
    let ran = ok(&run(
        &[
            "prep", OXOG, "-o", &out, "qc", "--preset", "generic", "--date", DATE,
        ],
        &config,
    ));
    assert!(ran.contains("step 9 (qc): verdict pass"), "{ran}");
    assert_eq!(
        names(&recipe),
        [
            "p.alignment_summary_metrics",
            "p.bai",
            "p.bam",
            "p.duplicate_metrics",
            "p.insert_size_metrics",
            "p.manifest.json",
            "p.marked.duplication_metrics",
            "p.quality_yield_metrics",
            "p.report.html",
            "p.verdict.json",
        ]
    );
    ok(&run(
        &[
            "chain",
            OXOG,
            "-o",
            &hand.join("p").display().to_string(),
            "+",
            "fix-mate-information",
            "+",
            "mark-duplicates",
            "+",
            "sort-sam",
            "--sort-order",
            "coordinate",
            "+",
            "quality-yield",
            "+",
            "alignment-summary",
            "+",
            "insert-size",
            "+",
            "duplicate-metrics",
            "+",
            "write",
            &hand.join("p.bam").display().to_string(),
            "--index",
            "+",
            "qc",
            "--preset",
            "generic",
            "--date",
            DATE,
        ],
        &config,
    ));
    same_files(&recipe, &hand);
    // The qc verb on the BAM the recipe wrote: the same verdict.
    let verb = fresh_dir("prep-qc-verb");
    ok(&run(
        &[
            "qc",
            &out,
            "-o",
            &verb.join("p").display().to_string(),
            "--preset",
            "generic",
            "--date",
            DATE,
        ],
        &config,
    ));
    assert_eq!(
        text(&recipe.join("p.verdict.json"), &recipe),
        text(&verb.join("p.verdict.json"), &verb)
    );
    // --fail-on gives the chain's exit code; the verdict passes here.
    let gate = fresh_dir("prep-qc-gate");
    let gate_out = gate.join("p.bam").display().to_string();
    ok(&run(
        &["prep", OXOG, "-o", &gate_out, "qc", "--fail-on", "warn"],
        &config,
    ));
    // Any other step beside a recipe, and a second qc step, are refused.
    let other = usage_error(&run(
        &["prep", OXOG, "-o", &gate_out, "insert-size"],
        &config,
    ));
    assert!(
        other.contains("--recipe prep gives the steps: a step on the command line cannot be added"),
        "{other}"
    );
    assert!(
        other.contains("a last qc step may follow a recipe"),
        "{other}"
    );
    fs::write(
        config.join("bammetrics/recipes/checked.toml"),
        "[recipe]\ndescription = \"marked and judged\"\noutput = \"bam\"\n\n[[step]]\n\
         verb = \"mark-duplicates\"\n\n[[step]]\nverb = \"qc\"\npreset = \"generic\"\n",
    )
    .expect("recipe");
    let twice = usage_error(&run(&["checked", MD, "-o", &gate_out, "qc"], &config));
    assert!(
        twice.contains("recipe checked ends in a qc step: a second one cannot follow it"),
        "{twice}"
    );
}

/// A user recipe that ends in qc: `qc.option=value` sets an option of its
/// qc step, `-o x.bam` puts the write in front of it, and the recipe equals
/// its chain in the short form.
#[test]
fn a_user_recipe_ends_in_qc() {
    let config = config_dir("user-qc-config", None);
    fs::write(
        config.join("bammetrics/recipes/checked.toml"),
        "[recipe]\ndescription = \"marked and judged\"\noutput = \"bam\"\n\n[[step]]\n\
         verb = \"mark-duplicates\"\n\n[[step]]\nverb = \"qc\"\npreset = \"generic\"\n",
    )
    .expect("recipe");
    let (recipe, hand) = (fresh_dir("user-qc"), fresh_dir("user-qc-hand"));
    let date = format!("qc.date={DATE}");
    ok(&run(
        &[
            "checked",
            MD,
            "-o",
            &recipe.join("m.bam").display().to_string(),
            &date,
        ],
        &config,
    ));
    ok(&run(
        &[
            "chain",
            MD,
            "mark-duplicates",
            "qc",
            "--preset",
            "generic",
            "--date",
            DATE,
            "-o",
            &hand.join("m.bam").display().to_string(),
        ],
        &config,
    ));
    same_files(&recipe, &hand);
    let files = names(&recipe);
    for name in ["m.bam", "m.bai", "m.verdict.json", "m.report.html"] {
        assert!(files.iter().any(|f| f == name), "{name}: {files:?}");
    }
    // The plan: write is step 2, qc step 3, the taps qc adds follow.
    let plan = run(
        &[
            "checked",
            MD,
            "-o",
            &recipe.join("m.bam").display().to_string(),
            "--dry-run",
        ],
        &config,
    );
    ok(&plan);
    let plan = String::from_utf8_lossy(&plan.stdout).into_owned();
    let at = |what: &str| plan.find(what).unwrap_or_else(|| panic!("{what}: {plan}"));
    assert!(at("step 2  write") < at("step 3  qc"));
    assert!(at("step 3  qc") < at("step 4  quality-yield"));
    // A step that would write the recipe file or the config file is refused.
    let recipe_file = config.join("bammetrics/recipes/checked.toml");
    fs::write(
        config.join("bammetrics/recipes/over.toml"),
        format!(
            "[recipe]\ndescription = \"writes itself\"\n\n[[step]]\nverb = \"insert-size\"\n\
             output = \"{}\"\n",
            config.join("bammetrics/recipes/over.toml").display()
        ),
    )
    .expect("recipe");
    let refused = usage_error(&run(
        &["over", IS, "-o", &hand.join("x").display().to_string()],
        &config,
    ));
    assert!(
        refused.contains("over.toml is an input; writing it would destroy it"),
        "{refused}"
    );
    assert!(recipe_file.is_file());
    fs::write(config.join("bammetrics/config.toml"), "threads = 2\n").expect("config");
    let config_file = config.join("bammetrics/config.toml").display().to_string();
    let refused = usage_error(&run(
        &[
            "chain",
            IS,
            "insert-size",
            "-o",
            &config_file,
            "write",
            "x.bam",
        ],
        &config,
    ));
    assert!(refused.contains("is an input"), "{refused}");
    assert_eq!(
        fs::read_to_string(config.join("bammetrics/config.toml")).expect("config"),
        "threads = 2\n"
    );
}

/// The short form with a materialized tool that writes a BAM: a file after
/// `merge-sam-files` is its further input, `-o m.bam` names the result,
/// `--keep` keeps the merged stream and `--resume` starts there.
#[test]
fn a_merge_in_the_short_form_keeps_and_resumes() {
    let config = config_dir("merge-config", None);
    let (short, long) = (fresh_dir("merge"), fresh_dir("merge-long"));
    let out = short.join("m.bam").display().to_string();
    let command = |pct: &'static str| {
        vec![
            "chain",
            MD,
            "add-comments-to-bam",
            "--comment",
            "c",
            "merge-sam-files",
            OXOG,
            "--keep",
            "insert-size",
            "--minimum-pct",
            pct,
            "-o",
            &out,
            "--resume",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<String>>()
    };
    ok(&run(&refs(&command("0.05")), &config));
    assert_eq!(
        names(&short),
        [
            "m.bai",
            "m.bam",
            "m.insert_size_metrics",
            "m.manifest.json",
            "m.merged.bai",
            "m.merged.bam",
        ]
    );
    // The kept stream is the BAM -o names: the same records.
    assert!(bam(&short.join("m.bam")) == bam(&short.join("m.merged.bam")));
    let before = modified(&short.join("m.merged.bam"));
    std::thread::sleep(std::time::Duration::from_millis(20));
    let changed = ok(&run(&refs(&command("0.3")), &config));
    for skipped in [
        "step 1 (add-comments-to-bam): skipped",
        "step 2 (merge-sam-files): skipped",
        "step 4 (write): skipped",
    ] {
        assert!(changed.contains(skipped), "{changed}");
    }
    assert_eq!(modified(&short.join("m.merged.bam")), before);
    // The long form of the changed chain.
    ok(&run(
        &[
            "chain",
            MD,
            "-o",
            &long.join("m").display().to_string(),
            "--resume",
            "+",
            "add-comments-to-bam",
            "--comment",
            "c",
            "+",
            "merge-sam-files",
            "--input",
            OXOG,
            "--keep",
            "+",
            "insert-size",
            "--minimum-pct",
            "0.3",
            "+",
            "write",
            &long.join("m.bam").display().to_string(),
            "--index",
        ],
        &config,
    ));
    same_files(&short, &long);
    // A merge alone changes the records: -o m.bam writes them.
    let alone = fresh_dir("merge-alone");
    ok(&run(
        &[
            "chain",
            MD,
            "merge-sam-files",
            OXOG,
            "-o",
            &alone.join("a.bam").display().to_string(),
        ],
        &config,
    ));
    assert_eq!(names(&alone), ["a.bai", "a.bam"]);
    assert!(bam(&alone.join("a.bam")).1.len() == bam(&short.join("m.bam")).1.len());
}

/// The qc step and the tools with no collector. `qc --preset wgs` behind a
/// transform needs `CollectWgsMetrics` on the marked records: it reads the
/// BAM `-o` names, the kept stream, or a copy of the stream (a materialized
/// step), and a `wgs` step the chain already runs is not added twice. Each
/// way gives the metrics `wgs` gives for the BAM the chain wrote.
#[test]
fn the_qc_step_measures_wgs_behind_a_transform() {
    let config = config_dir("wgs-config", None);
    let dir = fresh_dir("wgs");
    let path = |name: &str| dir.join(name).display().to_string();
    // `steps` follow mark-duplicates, `last` follows -o.
    let chain = |steps: &[&str], out: &str, last: &[&str]| {
        let mut all = vec!["chain", MD, "-R", REFERENCE, "mark-duplicates"];
        all.extend(steps);
        all.extend(["qc", "--preset", "wgs", "--date", DATE, "-o", out]);
        all.extend(last);
        let output = run(&all, &config);
        (
            ok(&output),
            String::from_utf8_lossy(&output.stdout).into_owned(),
        )
    };
    // -o names a BAM: the tool reads it after the chain.
    let (_, plan) = chain(&[], &path("a.bam"), &["--dry-run"]);
    assert!(
        plan.contains("wgs                own read, after the chain"),
        "{plan}"
    );
    chain(&[], &path("a.bam"), &[]);
    // --keep: the tool reads the kept stream.
    let (_, plan) = chain(&["--keep"], &path("b"), &["--dry-run"]);
    assert!(plan.contains("own read, after the chain"), "{plan}");
    chain(&["--keep"], &path("b"), &[]);
    // No file holds the stream: a materialized step of the read.
    let (_, plan) = chain(&[], &path("c"), &["--dry-run"]);
    assert!(plan.contains("materialize on marked, for qc"), "{plan}");
    assert!(
        plan.contains("materialize step 7 wgs: writes its stream to a file"),
        "{plan}"
    );
    chain(&[], &path("c"), &[]);
    // A wgs step of the chain's own: qc adds none.
    let (_, plan) = chain(&["wgs"], &path("d"), &["--dry-run"]);
    assert_eq!(plan.matches(" wgs ").count(), 1, "{plan}");
    assert!(
        plan.contains("step 2  wgs                materialize on marked"),
        "{plan}"
    );
    chain(&["wgs"], &path("d"), &[]);
    // The staged run: wgs on the BAM the first chain wrote.
    ok(&run(
        &[
            "wgs",
            &path("a.bam"),
            "-R",
            REFERENCE,
            "-o",
            &path("staged.wgs_metrics"),
        ],
        &config,
    ));
    let staged = text(&dir.join("staged.wgs_metrics"), &dir);
    for prefix in ["a", "b", "c", "d"] {
        let name = format!("{prefix}.wgs_metrics");
        assert_eq!(text(&dir.join(&name), &dir), staged, "{name}");
        let verdict =
            fs::read_to_string(dir.join(format!("{prefix}.verdict.json"))).expect("verdict");
        assert!(verdict.contains("WgsMetrics"), "{prefix}");
        assert_eq!(
            verdict.replace(&format!("/{prefix}."), "/P."),
            fs::read_to_string(dir.join("a.verdict.json"))
                .expect("verdict")
                .replace("/a.", "/P."),
            "{prefix}"
        );
    }
}

/// A materialized tap under `--resume`: its record is kept when nothing
/// changed, and a changed option reruns it alone, on the kept stream.
#[test]
fn a_materialized_step_resumes() {
    let config = config_dir("mat-config", None);
    let dir = fresh_dir("mat");
    let prefix = dir.join("s").display().to_string();
    let command = |cap: &'static str| {
        vec![
            "chain",
            MD,
            "-R",
            REFERENCE,
            "mark-duplicates",
            "--keep",
            "library-complexity",
            "wgs",
            "--coverage-cap",
            cap,
            "-o",
            &prefix,
            "--resume",
        ]
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<String>>()
    };
    ok(&run(&refs(&command("10")), &config));
    let again = ok(&run(&refs(&command("10")), &config));
    assert!(
        again.contains("3 steps finished, 3 of them skipped by --resume"),
        "{again}"
    );
    let before = modified(&dir.join("s.library_complexity_metrics"));
    std::thread::sleep(std::time::Duration::from_millis(20));
    let changed = ok(&run(&refs(&command("20")), &config));
    for skipped in [
        "step 1 (mark-duplicates): skipped",
        "step 2 (library-complexity): skipped",
    ] {
        assert!(changed.contains(skipped), "{changed}");
    }
    assert_eq!(modified(&dir.join("s.library_complexity_metrics")), before);
    let fresh = fresh_dir("mat-fresh");
    let fresh_prefix = fresh.join("s").display().to_string();
    let mut fresh_command = command("20");
    let at = fresh_command
        .iter()
        .position(|a| *a == prefix)
        .expect("prefix");
    fresh_command[at] = fresh_prefix;
    ok(&run(&refs(&fresh_command), &config));
    same_files(&dir, &fresh);
}

/// The errors of the short form with the new steps, each before the first
/// read: a step after qc, `--keep` on a tap, and `--resume` on a pipe.
#[test]
fn the_short_form_refuses_what_the_long_form_refuses() {
    let config = config_dir("errors-config", None);
    let dir = fresh_dir("errors");
    let out = dir.join("r.bam").display().to_string();
    let after = usage_error(&run(
        &[
            "chain",
            MD,
            "mark-duplicates",
            "qc",
            "insert-size",
            "-o",
            &out,
        ],
        &config,
    ));
    assert!(
        after.contains("step 2 (qc) is the last step: no step may follow it"),
        "{after}"
    );
    let keep = usage_error(&run(
        &[
            "chain",
            MD,
            "insert-size",
            "--keep",
            "-o",
            &dir.join("r").display().to_string(),
        ],
        &config,
    ));
    assert!(
        keep.contains("--keep writes the stream a transform leaves, and insert-size leaves none"),
        "{keep}"
    );
    let pipe = usage_error(&run(
        &["chain", "-", "mark-duplicates", "-o", &out, "--resume"],
        &config,
    ));
    assert!(
        pipe.contains("--resume compares the source's fingerprint, and standard input has none"),
        "{pipe}"
    );
    let write = usage_error(&run(
        &["chain", MD, "mark-duplicates", "qc", "--index", "-o", &out],
        &config,
    ));
    assert!(
        write.contains("step 2 (qc): --index is an option of write"),
        "{write}"
    );
    assert!(write.contains("write OUT.bam --index qc"), "{write}");
    assert!(names(&dir).is_empty(), "{:?}", names(&dir));
}

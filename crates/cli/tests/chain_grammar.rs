//! The simpler chain grammar (docs/chaining.md, sections 22 to 26): steps
//! without `+`, `-o` naming the result BAM, the user config file and
//! `BAMMETRICS_REFERENCE`, the dictionary defaults and the recipes. Each
//! new form is checked against the old form it stands for: the same files,
//! byte for byte, apart from the path of the output directory in the
//! command lines the metrics files carry.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

/// A fresh directory under the system temporary directory.
fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-chain-grammar-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("directory");
    dir
}

/// The binary from the repository root, with `config` as its config
/// directory (`XDG_CONFIG_HOME`) and no `BAMMETRICS_REFERENCE` unless
/// `env_reference` names one.
fn run_with(args: &[&str], config: &Path, env_reference: Option<&str>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bammetrics"));
    command
        .args(args)
        .current_dir(root())
        .env("XDG_CONFIG_HOME", config)
        .env_remove("BAMMETRICS_REFERENCE");
    if let Some(reference) = env_reference {
        command.env("BAMMETRICS_REFERENCE", reference);
    }
    command.output().expect("run bammetrics")
}

fn ok(output: &Output) -> String {
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(0), "{stderr}");
    stderr
}

fn usage_error(output: &Output) -> String {
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    stderr
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

/// The bytes of each file of `a` and `b`, the directory paths and the
/// clock line of a metrics file aside.
fn same_files(a: &Path, b: &Path) {
    assert_eq!(names(a), names(b));
    for name in names(a) {
        let (left, right) = (
            fs::read(a.join(&name)).expect("a"),
            fs::read(b.join(&name)).expect("b"),
        );
        if name.ends_with(".bam") || name.ends_with(".bai") {
            assert!(left == right, "{name} differs");
            continue;
        }
        let clean = |bytes: &[u8], dir: &Path| {
            String::from_utf8_lossy(bytes)
                .replace(&dir.display().to_string(), "DIR")
                .lines()
                .filter(|line| !line.starts_with("# Started on"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(clean(&left, a), clean(&right, b), "{name}");
    }
}

fn config_dir(name: &str, text: Option<&str>) -> PathBuf {
    let dir = fresh_dir(name);
    if let Some(text) = text {
        fs::create_dir_all(dir.join("bammetrics")).expect("config");
        fs::write(dir.join("bammetrics/config.toml"), text).expect("config file");
    }
    dir
}

const SOURCE: &str = "tests/data/insert_size.bam";
const REFERENCE: &str = "tests/data/ref_100k.fa";

/// The owner's example: `chain BAM reorder-sam insert-size -o r.bam` with
/// the reference from the config file, against `chain BAM -R FASTA -o r +
/// reorder-sam --sequence-dictionary FAI + insert-size + write r.bam
/// --index`.
#[test]
fn the_short_form_writes_what_the_long_form_writes() {
    let reference = root().join(REFERENCE);
    let config = config_dir(
        "short-config",
        Some(&format!("reference = \"{}\"\n", reference.display())),
    );
    let (new, old) = (fresh_dir("short-new"), fresh_dir("short-old"));
    let short = new.join("r.bam");
    let stderr = ok(&run_with(
        &[
            "chain",
            SOURCE,
            "reorder-sam",
            "insert-size",
            "-o",
            short.to_str().expect("utf-8"),
        ],
        &config,
        None,
    ));
    assert!(stderr.contains("SEQUENCE_DICTIONARY: the dictionary is the reference's index"));
    let prefix = old.join("r");
    let fai = format!("{REFERENCE}.fai");
    ok(&run_with(
        &[
            "chain",
            SOURCE,
            "-R",
            REFERENCE,
            "-o",
            prefix.to_str().expect("utf-8"),
            "+",
            "reorder-sam",
            "--sequence-dictionary",
            &fai,
            "+",
            "insert-size",
            "+",
            "write",
            old.join("r.bam").to_str().expect("utf-8"),
            "--index",
        ],
        &config,
        None,
    ));
    assert_eq!(names(&new), ["r.bai", "r.bam", "r.insert_size_metrics"]);
    same_files(&new, &old);
}

/// Steps without `+` chain as the `+` form does: mark-duplicates, two
/// taps and write, records handed on in memory.
#[test]
fn steps_without_plus_equal_steps_with_plus() {
    let config = config_dir("bare-config", None);
    let (bare, plus) = (fresh_dir("bare"), fresh_dir("plus"));
    let out = |dir: &Path, name: &str| dir.join(name).display().to_string();
    ok(&run_with(
        &[
            "chain",
            "tests/data/markdup.bam",
            "mark-duplicates",
            "--program-group-command-line",
            "pinned",
            "--as",
            "insert-size",
            "insert-size",
            "--minimum-pct",
            "0.1",
            "quality-yield",
            "write",
            &out(&bare, "m.bam"),
            "--index",
            "-o",
            &out(&bare, "s"),
        ],
        &config,
        None,
    ));
    ok(&run_with(
        &[
            "chain",
            "tests/data/markdup.bam",
            "-o",
            &out(&plus, "s"),
            "+",
            "mark-duplicates",
            "--program-group-command-line",
            "pinned",
            "--as",
            "insert-size",
            "+",
            "insert-size",
            "--minimum-pct",
            "0.1",
            "+",
            "quality-yield",
            "+",
            "write",
            &out(&plus, "m.bam"),
            "--index",
        ],
        &config,
        None,
    ));
    same_files(&bare, &plus);
    assert!(names(&bare).contains(&"s.insert-size.duplication_metrics".to_owned()));
}

/// `BAMMETRICS_REFERENCE` replaces the file's reference, `-R` both, and
/// `--no-config` leaves the file unread.
#[test]
fn the_reference_comes_from_r_then_the_variable_then_the_file() {
    let config = config_dir(
        "order-config",
        Some("reference = \"/nowhere/file.fa\"\nthreads = 2\n"),
    );
    let printed = ok(&run_with(&["config"], &config, None));
    assert_eq!(printed, "");
    let text = String::from_utf8(run_with(&["config"], &config, None).stdout).expect("text");
    assert!(
        text.contains("reference    /nowhere/file.fa (from "),
        "{text}"
    );
    assert!(text.contains("config.toml:1)"), "{text}");
    let env = root().join(REFERENCE).display().to_string();
    let text = String::from_utf8(run_with(&["config"], &config, Some(&env)).stdout).expect("text");
    assert!(
        text.contains(&format!("reference    {env} (from BAMMETRICS_REFERENCE)")),
        "{text}"
    );
    let text = String::from_utf8(run_with(&["config", "--no-config"], &config, None).stdout)
        .expect("text");
    assert!(text.contains("(not read (--no-config))") && text.contains("reference    none"));
    // The file's reference is missing: an input error that names the file.
    let dir = fresh_dir("order-run");
    let r = dir.join("r.bam").display().to_string();
    let output = run_with(&["chain", SOURCE, "reorder-sam", "-o", &r], &config, None);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains("/nowhere/file.fa not found (the reference of "),
        "{stderr}"
    );
    // The variable or -R wins over the file.
    ok(&run_with(
        &["chain", SOURCE, "reorder-sam", "-o", &r],
        &config,
        Some(&env),
    ));
    ok(&run_with(
        &["chain", SOURCE, "reorder-sam", "-R", REFERENCE, "-o", &r],
        &config,
        None,
    ));
    // --no-config: no reference, so reorder-sam has no dictionary.
    let stderr = usage_error(&run_with(
        &["chain", SOURCE, "reorder-sam", "-o", &r, "--no-config"],
        &config,
        None,
    ));
    assert!(stderr.contains("missing --sequence-dictionary"), "{stderr}");
}

/// A tool verb that reads a reference takes the config's; the Picard face
/// never does.
#[test]
fn a_tool_verb_takes_the_config_reference_and_the_picard_face_does_not() {
    let reference = root().join(REFERENCE);
    let config = config_dir(
        "verb-config",
        Some(&format!("reference = \"{}\"\n", reference.display())),
    );
    let dir = fresh_dir("verb");
    let verb = dir.join("verb.txt");
    ok(&run_with(
        &[
            "alignment-summary",
            SOURCE,
            "-o",
            verb.to_str().expect("utf-8"),
        ],
        &config,
        None,
    ));
    let text = fs::read_to_string(&verb).expect("metrics");
    assert!(
        text.contains(&format!("REFERENCE_SEQUENCE={}", reference.display())),
        "{text}"
    );
    let picard = dir.join("picard.txt");
    ok(&run_with(
        &[
            "CollectAlignmentSummaryMetrics",
            &format!("I={SOURCE}"),
            &format!("O={}", picard.display()),
        ],
        &config,
        None,
    ));
    let text = fs::read_to_string(&picard).expect("metrics");
    assert!(!text.contains("REFERENCE_SEQUENCE="), "{text}");
    let unread = dir.join("unread.txt");
    ok(&run_with(
        &[
            "alignment-summary",
            SOURCE,
            "-o",
            unread.to_str().expect("utf-8"),
            "--no-config",
        ],
        &config,
        None,
    ));
    let text = fs::read_to_string(&unread).expect("metrics");
    assert!(!text.contains("REFERENCE_SEQUENCE="), "{text}");
}

/// The recipe qc-only writes the four files of the same taps chained by
/// hand; a user recipe replaces a shipped one of its name.
#[test]
fn a_recipe_is_its_chain() {
    let config = config_dir("recipe-config", None);
    let (recipe, hand) = (fresh_dir("recipe"), fresh_dir("hand"));
    let prefix = |dir: &Path| dir.join("s").display().to_string();
    ok(&run_with(
        &["qc-only", SOURCE, "-o", &prefix(&recipe)],
        &config,
        None,
    ));
    ok(&run_with(
        &[
            "chain",
            SOURCE,
            "quality-yield",
            "alignment-summary",
            "insert-size",
            "duplicate-metrics",
            "-o",
            &prefix(&hand),
        ],
        &config,
        None,
    ));
    same_files(&recipe, &hand);
    // A user recipe of the same name runs in its place.
    let recipes = config.join("bammetrics/recipes");
    fs::create_dir_all(&recipes).expect("recipes");
    fs::write(
        recipes.join("qc-only.toml"),
        "[recipe]\ndescription = \"mine\"\n\n[[step]]\nverb = \"insert-size\"\n",
    )
    .expect("recipe");
    let mine = fresh_dir("mine");
    ok(&run_with(
        &["qc-only", SOURCE, "-o", &prefix(&mine)],
        &config,
        None,
    ));
    assert_eq!(names(&mine), ["s.insert_size_metrics"]);
    let listed = String::from_utf8(run_with(&["recipes"], &config, None).stdout).expect("text");
    assert!(listed.contains("qc-only     mine (yours)"), "{listed}");
    // A user recipe over the cap is refused, never read whole.
    fs::write(recipes.join("big.toml"), vec![b'#'; 70_000]).expect("big");
    let stderr = usage_error(&run_with(
        &["big", SOURCE, "-o", &prefix(&mine)],
        &config,
        None,
    ));
    assert!(
        stderr.contains("not a recipe: larger than 65536 bytes"),
        "{stderr}"
    );
}

/// The downsample recipe asks for its probability, then writes the BAM
/// the same chain by hand writes.
#[test]
fn the_downsample_recipe_takes_its_probability_from_the_command_line() {
    let config = config_dir("down-config", None);
    let (recipe, hand) = (fresh_dir("down-recipe"), fresh_dir("down-hand"));
    let out = |dir: &Path| dir.join("d.bam").display().to_string();
    let stderr = usage_error(&run_with(
        &["downsample", "tests/data/oxog.bam", "-o", &out(&recipe)],
        &config,
        None,
    ));
    assert!(stderr.contains("recipe downsample needs downsample-sam.probability=VALUE"));
    ok(&run_with(
        &[
            "downsample",
            "tests/data/oxog.bam",
            "-o",
            &out(&recipe),
            "downsample-sam.probability=0.5",
        ],
        &config,
        None,
    ));
    ok(&run_with(
        &[
            "chain",
            "tests/data/oxog.bam",
            "downsample-sam",
            "--strategy",
            "ConstantMemory",
            "--probability",
            "0.5",
            "-o",
            &out(&hand),
        ],
        &config,
        None,
    ));
    same_files(&recipe, &hand);
    assert_eq!(names(&recipe), ["d.bai", "d.bam"]);
}

//! `interval-list-to-bed` and `vcf-to-interval-list` (wave W12) without a
//! JVM: `qc --panel` takes what each writes, a replaced file keeps its mode,
//! a failed run leaves nothing (the directories it would have made
//! included), and an output that is the input is refused. The Picard bytes
//! are in tests/golden and tests/concordance.

use std::fs::{self, DirBuilder};
use std::os::unix::fs::{DirBuilderExt, PermissionsExt};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const BAM: &str = "tests/data/insert_size.bam";
const TARGETS: &str = "tests/data/targets.interval_list";

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

fn ok(args: &[&str]) -> Output {
    let run = bammetrics(args);
    assert_eq!(run.status.code(), Some(0), "{args:?}: {run:?}");
    run
}

fn fresh_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-interval-converters-{}-{name}-{nanos}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("a new temp dir");
    dir
}

fn text(path: &Path) -> String {
    path.to_str().expect("utf-8 path").to_owned()
}

/// The rows of a metrics file without the `BAIT_SET` column, which is the
/// panel file's name.
fn hs_rows(path: &Path) -> Vec<String> {
    let body = fs::read_to_string(path).expect("hs metrics");
    let lines: Vec<&str> = body
        .lines()
        .skip_while(|l| !l.starts_with("BAIT_SET"))
        .take_while(|l| !l.is_empty())
        .collect();
    lines
        .iter()
        .map(|l| l.split('\t').skip(1).collect::<Vec<_>>().join("\t"))
        .collect()
}

fn qc_hs(panel: &Path, dir: &Path, name: &str) -> Vec<String> {
    let prefix = dir.join(name);
    ok(&["qc", BAM, "-o", &text(&prefix), "--panel", &text(panel)]);
    hs_rows(&dir.join(format!("{name}.hs_metrics")))
}

#[test]
fn qc_panel_takes_a_bed_from_interval_list_to_bed() {
    let dir = fresh_dir("bed");
    let bed = dir.join("targets.bed");
    ok(&["interval-list-to-bed", TARGETS, "-o", &text(&bed)]);
    assert_eq!(
        qc_hs(&bed, &dir, "from_bed"),
        qc_hs(&root().join(TARGETS), &dir, "from_list")
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn qc_panel_takes_a_list_from_vcf_to_interval_list() {
    let dir = fresh_dir("vcf");
    let vcf = dir.join("sites.vcf");
    fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         chr1\t1001\tt1\tA\t<DEL>\t.\t.\tEND=1200\nchr1\t2001\tt2\tA\t<DEL>\t.\t.\tEND=2400\n\
         chr1\t4001\tt3\tA\t<DEL>\t.\t.\tEND=4300\nchr1\t5001\tt4\tA\t<DEL>\t.\t.\tEND=5600\n\
         chr1\t8001\tt5\tA\t<DEL>\t.\t.\tEND=8300\n",
    )
    .expect("vcf");
    let list = dir.join("sites.interval_list");
    ok(&["vcf-to-interval-list", &text(&vcf), "-o", &text(&list)]);
    assert_eq!(
        qc_hs(&list, &dir, "from_vcf"),
        qc_hs(&root().join(TARGETS), &dir, "from_list")
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_replaced_file_keeps_its_mode_and_a_failed_run_writes_nothing() {
    let dir = fresh_dir("files");
    let bed = dir.join("kept.bed");
    fs::write(&bed, "old\n").expect("old");
    fs::set_permissions(&bed, fs::Permissions::from_mode(0o640)).expect("mode");
    ok(&["interval-list-to-bed", TARGETS, "-o", &text(&bed)]);
    let mode = fs::metadata(&bed).expect("bed").permissions().mode() & 0o777;
    assert_eq!(mode, 0o640);
    assert!(
        fs::read_to_string(&bed)
            .expect("bed")
            .starts_with("chr1\t1000\t1200\tt1\t500\t+\n")
    );

    let missing = dir.join("new/sub/out.interval_list");
    let run = bammetrics(&[
        "vcf-to-interval-list",
        "tests/data/intervals_bad.vcf",
        "-o",
        &text(&missing),
    ]);
    assert_eq!(run.status.code(), Some(1), "{run:?}");
    assert!(
        !dir.join("new").exists(),
        "a failed run left its directories"
    );

    let run = bammetrics(&["interval-list-to-bed", TARGETS, "-o", TARGETS]);
    assert_eq!(run.status.code(), Some(2), "{run:?}");
    assert!(String::from_utf8_lossy(&run.stderr).contains("is an input"));
    let _ = fs::remove_dir_all(&dir);
}

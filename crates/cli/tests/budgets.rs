//! The size and speed budgets of docs/05 section 4.8, measured on the
//! committed fixtures. Each budget has a target; the test fails when a
//! measured value reaches twice its target (decision 11 of docs/05 section
//! 8) and prints `over target` when it passes the target alone.
//!
//! Ignored by `cargo test`, because the times only mean something on a
//! release build and the binary growth needs the two probe examples built:
//!
//! ```text
//! make budgets
//! ```
//!
//! - Render times are the median of [`RUNS`] renders after [`WARMUP`]
//!   unmeasured ones: the failing fixture page with every plot (the
//!   largest page the fixtures give), and an index of 100 samples.
//! - Page sizes are the largest of the three fixture pages (passing,
//!   failing, failing with every plot); the index is 100 samples with the
//!   14 rules of `clinical-panel`.
//! - Binary growth is the size of the stripped `size_report` example less
//!   the stripped `size_base` example (`crates/report/examples`): the code
//!   and data the verdict and report crates add to a release binary that
//!   already reads metrics files. Both are built with the release profile
//!   of the workspace, on the machine that runs the test, so the number
//!   holds on any CI runner. The unstripped difference is printed too.
//! - Added BAM passes: `qc` reports one BAM read for the collectors, the
//!   verdict and the page together, and `report` and `check` run on the
//!   `qc` output after the BAM is deleted.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{Duration, Instant};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_report::{Input, Options, SampleSummary, Verdict, page_name, render, render_index};
use bammetrics_verdict::{Preset, evaluate};

const DATE: &str = "2026-01-02 03:04:05 UTC";
const RUNS: usize = 31;
const WARMUP: usize = 3;
const KB: f64 = 1024.0;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

/// The metrics file `name` recorded in `tests/golden/<record>.txt`.
fn golden_input(record: &str, name: &str) -> Input {
    let path = root().join("tests/golden").join(format!("{record}.txt"));
    let text = fs::read_to_string(&path).expect("golden record");
    let marker = format!("## file {name}\n");
    let start = text.find(&marker).expect("file section") + marker.len();
    let rest = &text[start..];
    let end = rest.find("\n## file ").map_or(rest.len(), |i| i + 1);
    Input {
        path: format!("fixture.{name}"),
        doc: MetricsDoc::parse(&rest[..end]).expect("metrics file"),
    }
}

/// The seven metrics files bammetrics writes for `tests/data/insert_size.bam`,
/// as the report tests use them.
fn fixture_set() -> Vec<Input> {
    vec![
        golden_input("picard_alignment_summary", "alignment_summary.txt"),
        golden_input("picard_insert_size", "insert_size.txt"),
        golden_input("picard_quality_yield", "quality_yield.txt"),
        golden_input("picard_duplicate_metrics", "duplicate_metrics.txt"),
        golden_input("picard_hs_metrics", "hs_metrics.txt"),
        golden_input("picard_gc_bias", "gc_bias.txt"),
        golden_input("picard_mean_quality_by_cycle", "mqbc.txt"),
    ]
}

fn verdict(preset: &str, inputs: &[Input]) -> Verdict {
    let preset = Preset::builtin(preset)
        .expect("a shipped preset")
        .expect("parses");
    evaluate(&preset, inputs, Some("S1"))
}

/// The median time of `RUNS` calls of `f`, after `WARMUP` calls.
fn median(mut f: impl FnMut() -> usize) -> Duration {
    for _ in 0..WARMUP {
        assert!(f() > 0);
    }
    let mut times: Vec<Duration> = (0..RUNS)
        .map(|_| {
            let start = Instant::now();
            assert!(f() > 0);
            start.elapsed()
        })
        .collect();
    times.sort();
    times[RUNS / 2]
}

/// One measured budget: its name, target, measured value and unit.
struct Budget {
    name: &'static str,
    target: f64,
    measured: f64,
    unit: &'static str,
}

impl Budget {
    /// At twice the target or more; for a target of 0, any value above 0.
    fn fails(&self) -> bool {
        if self.target == 0.0 {
            self.measured > 0.0
        } else {
            self.measured >= 2.0 * self.target
        }
    }

    fn status(&self) -> &'static str {
        if self.fails() {
            "FAIL (2x target)"
        } else if self.measured >= self.target && self.target > 0.0 {
            "over target"
        } else {
            "ok"
        }
    }
}

/// The directory of the release examples, beside the directory of this
/// test binary (`target/release/deps`).
fn examples_dir() -> PathBuf {
    let exe = std::env::current_exe().expect("test binary path");
    exe.parent()
        .and_then(Path::parent)
        .expect("target profile directory")
        .join("examples")
}

/// The size of `path` after `strip`, from a copy in `dir`.
fn stripped_size(path: &Path, dir: &Path) -> u64 {
    let out = dir.join(path.file_name().expect("file name"));
    let status = Command::new("strip")
        .arg("-o")
        .arg(&out)
        .arg(path)
        .status()
        .expect("run strip");
    assert!(status.success(), "strip {}", path.display());
    fs::metadata(&out).expect("stripped probe").len()
}

fn size(path: &Path) -> u64 {
    fs::metadata(path)
        .unwrap_or_else(|e| {
            panic!(
                "{}: {e}; build the probes first (make budgets)",
                path.display()
            )
        })
        .len()
}

fn bammetrics(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

/// BAM reads the report adds: `qc` on a copy of the fixture BAM must say
/// one BAM read, and `report` and `check` must run on its output once the
/// copy is gone. Returns the added reads (0) for the table.
fn added_bam_passes(dir: &Path) -> f64 {
    let bam = dir.join("copy.bam");
    fs::copy(root().join("tests/data/insert_size.bam"), &bam).expect("copy the fixture BAM");
    let out = dir.join("run");
    let prefix = out.join("S1");
    let qc = bammetrics(&[
        "qc",
        bam.to_str().expect("utf-8"),
        "-R",
        "tests/data/ref_100k.fa",
        "--panel",
        "tests/data/targets.interval_list",
        "--preset",
        "clinical-panel",
        "-o",
        prefix.to_str().expect("utf-8"),
        "--date",
        DATE,
        "-q",
    ]);
    let summary = String::from_utf8_lossy(&qc.stdout).into_owned();
    assert_eq!(qc.status.code(), Some(0), "{summary}");
    let header = summary.lines().next().unwrap_or_default();
    assert!(out.join("S1.report.html").is_file());
    fs::remove_file(&bam).expect("delete the BAM copy");
    let report = bammetrics(&[
        "report",
        out.to_str().expect("utf-8"),
        "--preset",
        "clinical-panel",
        "--date",
        DATE,
    ]);
    assert_eq!(
        report.status.code(),
        Some(0),
        "{}",
        String::from_utf8_lossy(&report.stderr)
    );
    let check = bammetrics(&[
        "check",
        out.to_str().expect("utf-8"),
        "--preset",
        "clinical-panel",
    ]);
    // 10: the fixture fails clinical-panel. 1 would be an unreadable input.
    assert_eq!(check.status.code(), Some(10));
    let reads: f64 = header
        .find(" BAM read")
        .and_then(|end| header[..end].rsplit(' ').next())
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("no BAM read count in {header:?}"));
    reads - 1.0
}

#[test]
fn a_budget_fails_at_twice_its_target() {
    let at = |target: f64, measured: f64| Budget {
        name: "t",
        target,
        measured,
        unit: "ms",
    };
    assert_eq!(at(20.0, 19.9).status(), "ok");
    assert_eq!(at(20.0, 20.0).status(), "over target");
    assert_eq!(at(20.0, 39.9).status(), "over target");
    assert!(!at(20.0, 39.9).fails());
    assert!(at(20.0, 40.0).fails());
    assert!(!at(0.0, 0.0).fails());
    assert!(at(0.0, 1.0).fails());
}

#[test]
#[ignore = "release build and probe examples: run make budgets"]
fn budgets_of_docs_05_section_4_8() {
    let dir = std::env::temp_dir().join(format!("bammetrics-budgets-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temporary directory");

    let inputs = fixture_set();
    let pass = verdict("generic", &inputs);
    let fail = verdict("clinical-panel", &inputs);
    let plain = Options {
        all_plots: false,
        date: DATE,
    };
    let all = Options {
        all_plots: true,
        date: DATE,
    };
    let pages = [
        render(&pass, &inputs, &plain).len(),
        render(&fail, &inputs, &plain).len(),
        render(&fail, &inputs, &all).len(),
    ];
    let page_time = median(|| render(&fail, &inputs, &all).len());

    let samples: Vec<SampleSummary> = (0..100)
        .map(|i| {
            let mut v = fail.clone();
            v.sample = format!("NA{:05}", 12878 + i);
            let href = page_name(&v.sample);
            SampleSummary::new(v, href)
        })
        .collect();
    let index = render_index(&samples, &plain).len();
    let index_time = median(|| render_index(&samples, &plain).len());

    let examples = examples_dir();
    let (base, report) = (examples.join("size_base"), examples.join("size_report"));
    let raw_growth = size(&report) as f64 - size(&base) as f64;
    let growth = stripped_size(&report, &dir) as f64 - stripped_size(&base, &dir) as f64;

    let passes = added_bam_passes(&dir);
    let _ = fs::remove_dir_all(&dir);

    let budgets = [
        Budget {
            name: "render, sample page",
            target: 20.0,
            measured: page_time.as_secs_f64() * 1000.0,
            unit: "ms",
        },
        Budget {
            name: "render, index of 100 samples",
            target: 100.0,
            measured: index_time.as_secs_f64() * 1000.0,
            unit: "ms",
        },
        Budget {
            name: "sample page size",
            target: 100.0,
            measured: *pages.iter().max().expect("three pages") as f64 / KB,
            unit: "KB",
        },
        Budget {
            name: "index size, 100 samples",
            target: 200.0,
            measured: index as f64 / KB,
            unit: "KB",
        },
        Budget {
            name: "binary growth",
            target: 300.0,
            measured: growth / KB,
            unit: "KB",
        },
        Budget {
            name: "added passes over the BAM",
            target: 0.0,
            measured: passes,
            unit: "",
        },
    ];
    println!("budget                         target    limit (2x)  measured    status");
    for b in &budgets {
        println!(
            "{:<30} {:>6} {:<2} {:>6} {:<2}  {:>9.3} {:<2}  {}",
            b.name,
            b.target,
            b.unit,
            2.0 * b.target,
            b.unit,
            b.measured,
            b.unit,
            b.status()
        );
    }
    println!(
        "pages: {} bytes passing, {} failing, {} failing with every plot; index {} bytes",
        pages[0], pages[1], pages[2], index
    );
    println!(
        "binary growth: {:.0} bytes stripped, {:.0} bytes unstripped ({} less {})",
        growth,
        raw_growth,
        report.display(),
        base.display()
    );
    let failed: Vec<&str> = budgets
        .iter()
        .filter(|b| b.fails())
        .map(|b| b.name)
        .collect();
    assert!(failed.is_empty(), "at 2x target or more: {failed:?}");
}

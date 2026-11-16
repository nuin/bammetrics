//! Write the PDF report of the committed fixtures into a directory, the
//! PDF twin of `crates/report/examples/fixture_pages.rs`:
//!
//! ```text
//! cargo run -p bammetrics-pdf --features pdf --example fixture_pdfs -- OUT_DIR [--letter] [--font PATH]
//! ```
//!
//! Writes `pass.pdf` (preset `generic`), `fail.pdf` (`clinical-panel`) and
//! `fail-all-plots.pdf` from the metrics files recorded in `tests/golden`,
//! and prints each path with its page count and size.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_pdf::{Options, Paper, write};
use bammetrics_report::Input;
use bammetrics_verdict::{Preset, evaluate};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

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

fn main() {
    let mut args = std::env::args().skip(1);
    let out = PathBuf::from(
        args.next()
            .expect("usage: fixture_pdfs OUT_DIR [--letter] [--font PATH]"),
    );
    let mut paper = Paper::A4;
    let mut font = None;
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--letter" => paper = Paper::Letter,
            "--font" => font = args.next().map(PathBuf::from),
            other => panic!("unknown argument {other}"),
        }
    }
    fs::create_dir_all(&out).expect("output directory");
    let inputs = vec![
        golden_input("picard_alignment_summary", "alignment_summary.txt"),
        golden_input("picard_insert_size", "insert_size.txt"),
        golden_input("picard_quality_yield", "quality_yield.txt"),
        golden_input("picard_duplicate_metrics", "duplicate_metrics.txt"),
        golden_input("picard_hs_metrics", "hs_metrics.txt"),
        golden_input("picard_gc_bias", "gc_bias.txt"),
        golden_input("picard_mean_quality_by_cycle", "mqbc.txt"),
    ];
    let preset = |name: &str| Preset::builtin(name).expect("shipped").expect("parses");
    let files = [
        ("pass.pdf", "generic", false),
        ("fail.pdf", "clinical-panel", false),
        ("fail-all-plots.pdf", "clinical-panel", true),
    ];
    for (file, preset_name, all_plots) in files {
        let verdict = evaluate(&preset(preset_name), &inputs, Some("S1"));
        let options = Options {
            date: "2026-01-02 03:04:05 UTC".to_owned(),
            font: font.clone(),
            all_plots,
            paper,
        };
        let path = out.join(file);
        if let Err(e) = write(&path, &verdict, &inputs, &options) {
            eprintln!("{}: {e}", path.display());
            std::process::exit(1);
        }
        let bytes = fs::read(&path).expect("written PDF");
        let pages = bytes
            .windows(11)
            .filter(|w| w.starts_with(b"/Type/Page") && w[10] != b's')
            .count();
        println!("{}: {} bytes, {pages} pages", path.display(), bytes.len());
    }
}

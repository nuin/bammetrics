//! Write the sample page for the committed fixtures into a directory, for
//! looking at it in a browser without a directory of metrics files
//!
//! ```text
//! cargo run -p bammetrics-report --example fixture_pages -- OUT_DIR
//! ```
//!
//! Writes `pass.html` (preset `generic`), `fail.html` (`clinical-panel`) and
//! `fail-all-plots.html` from the metrics files recorded in `tests/golden`.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_report::{Input, Options, write};
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
    let out = PathBuf::from(
        std::env::args()
            .nth(1)
            .expect("usage: fixture_pages OUT_DIR"),
    );
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
    let date = "2026-01-02 03:04:05 UTC";
    let pages = [
        ("pass.html", "generic", false),
        ("fail.html", "clinical-panel", false),
        ("fail-all-plots.html", "clinical-panel", true),
    ];
    for (file, preset_name, all_plots) in pages {
        let verdict = evaluate(&preset(preset_name), &inputs, Some("S1"));
        let path = out.join(file);
        write(&path, &verdict, &inputs, &Options { all_plots, date }).expect("write page");
        println!("{}", path.display());
    }
}

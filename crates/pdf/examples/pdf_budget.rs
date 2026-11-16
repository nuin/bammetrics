//! The PDF rows of the size and speed budgets (docs/05 section 4.8), for
//! `make budgets`: renders the three fixture samples of `fixture_pdfs` and
//! prints two lines the budget test reads,
//!
//! ```text
//! pdf_bytes PASS FAIL FAIL_ALL_PLOTS
//! pdf_render_ms MEDIAN
//! ```
//!
//! the sizes of the three PDFs and the median time of 31 renders (after 3
//! unmeasured) of the failing sample with every plot, on A4.
//!
//! ```text
//! cargo build --release -p bammetrics-pdf --features pdf --example pdf_budget
//! ```

use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_pdf::{Options, render};
use bammetrics_report::Input;
use bammetrics_verdict::{Preset, evaluate};

const RUNS: usize = 31;
const WARMUP: usize = 3;

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
    let pass = evaluate(&preset("generic"), &inputs, Some("S1"));
    let fail = evaluate(&preset("clinical-panel"), &inputs, Some("S1"));
    let options = |all_plots| Options {
        date: "2026-01-02 03:04:05 UTC".to_owned(),
        all_plots,
        ..Options::default()
    };
    let bytes = |verdict, all_plots| {
        render(verdict, &inputs, &options(all_plots))
            .expect("renders")
            .len()
    };
    println!(
        "pdf_bytes {} {} {}",
        bytes(&pass, false),
        bytes(&fail, false),
        bytes(&fail, true)
    );
    for _ in 0..WARMUP {
        bytes(&fail, true);
    }
    let mut times: Vec<f64> = (0..RUNS)
        .map(|_| {
            let start = Instant::now();
            bytes(&fail, true);
            start.elapsed().as_secs_f64() * 1000.0
        })
        .collect();
    times.sort_by(f64::total_cmp);
    println!("pdf_render_ms {:.3}", times[RUNS / 2]);
}

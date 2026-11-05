//! Binary-growth probe, with the report: `size_base` plus every entry point
//! the command calls in the verdict and report crates (each shipped preset,
//! `evaluate`, `verdict.json` out and back, the sample page, the index).
//! See `size_base.rs` for what the size difference measures.

use std::path::Path;

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_report::{Input, Options, SampleSummary, Tally, page_name, render, render_index};
use bammetrics_verdict::{Preset, Verdict, evaluate, names};

fn main() {
    let mut inputs = Vec::new();
    for path in std::env::args().skip(1) {
        match MetricsDoc::read(Path::new(&path)) {
            Ok(doc) => {
                println!("{path}: {} classes", doc.classes.len());
                inputs.push(Input { path, doc });
            }
            Err(e) => eprintln!("{path}: {e}"),
        }
    }
    let options = Options {
        all_plots: true,
        date: "",
    };
    let mut samples = Vec::new();
    for name in names() {
        let preset = match Preset::resolve(name) {
            Ok(preset) => preset,
            Err(e) => {
                eprintln!("{name}: {e}");
                continue;
            }
        };
        let verdict = evaluate(&preset, &inputs, None);
        let json = verdict.to_json();
        let verdict = Verdict::from_json(&json).unwrap_or(verdict);
        let page = render(&verdict, &inputs, &options);
        println!(
            "{name}: {} bytes of JSON, {} bytes of page",
            json.len(),
            page.len()
        );
        let href = page_name(&verdict.sample);
        samples.push(SampleSummary::new(verdict, href));
    }
    let index = render_index(&samples, &options);
    println!(
        "{}: {} bytes of index",
        Tally::of(&samples).line(),
        index.len()
    );
}

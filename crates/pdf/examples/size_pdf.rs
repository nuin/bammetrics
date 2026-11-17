//! Binary-growth probe of the PDF report: `size_report` (crates/report)
//! plus the PDF of every page it renders, so the stripped size of this
//! release binary less the stripped `size_report` is the code and data the
//! `pdf` feature adds (krilla, krilla-svg, the layout, the fonts). `make
//! budgets` builds both and the budget test compares them.
//!
//! ```text
//! cargo build --release -p bammetrics-pdf --features pdf --example size_pdf
//! ```

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
    let pdf_options = bammetrics_pdf::Options {
        all_plots: true,
        ..bammetrics_pdf::Options::default()
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
        match bammetrics_pdf::render(&verdict, &inputs, &pdf_options) {
            Ok(pdf) => println!(
                "{name}: {} bytes of page, {} bytes of PDF",
                page.len(),
                pdf.len()
            ),
            Err(e) => eprintln!("{name}: {e}"),
        }
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

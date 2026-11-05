//! Binary-growth probe, the baseline: read each metrics file named on the
//! command line and print its class count. `size_report` does the same and
//! then evaluates every shipped preset and renders the pages, so the size
//! difference of the two release binaries is the code and data the verdict
//! and report crates add to a binary that already reads metrics files (the
//! binary growth budget of docs/05 section 4.8; `make budgets` measures it).
//!
//! ```text
//! cargo build --release -p bammetrics-report --example size_base --example size_report
//! ```

use std::path::Path;

use bammetrics_core::metrics_file::MetricsDoc;

fn main() {
    for path in std::env::args().skip(1) {
        match MetricsDoc::read(Path::new(&path)) {
            Ok(doc) => println!("{path}: {} classes", doc.classes.len()),
            Err(e) => eprintln!("{path}: {e}"),
        }
    }
}

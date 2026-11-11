//! Binary-growth probe, the baseline: read each metrics file named on the
//! command line, print its class count, the registry tool its first
//! header names and how many of its columns the glossary types.
//! `size_formats` does the same and then writes the file in every format,
//! so the size difference of the two release binaries is the code this
//! crate adds to a binary that already reads metrics files and holds the
//! tool registry and the glossary, as `bammetrics` does (the formats
//! growth budget; `make budgets` measures it). The registry and glossary
//! lookups are in the baseline because this crate makes the same calls:
//! without them the difference would count every tool of the core and the
//! 113 KB glossary table. The upper-case name and the letter count are
//! there for the Unicode tables of `core`, which `check` and the report
//! page already link.
//!
//! ```text
//! cargo build --release -p bammetrics-formats --example size_formats_base --example size_formats
//! ```

use std::path::Path;

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_core::{glossary, registry};

fn main() {
    for path in std::env::args().skip(1) {
        match MetricsDoc::read(Path::new(&path)) {
            Ok(doc) => {
                let tool = doc
                    .headers
                    .first()
                    .and_then(|header| header.split_whitespace().next())
                    .and_then(registry::find)
                    .map_or("no registry tool", |tool| tool.name);
                let typed = doc
                    .classes
                    .iter()
                    .flat_map(|block| {
                        block
                            .columns
                            .iter()
                            .map(|column| glossary::lookup(&block.class, column))
                    })
                    .chain(doc.histogram.iter().flat_map(|histogram| {
                        let key = histogram.columns.first();
                        histogram.columns.iter().map(move |label| {
                            glossary::histogram_lookup(tool, key.unwrap_or(label), label)
                        })
                    }))
                    .flatten()
                    .count();
                println!(
                    "{path}: {} classes, {}, {} letters and digits, {typed} columns in the glossary",
                    doc.classes.len(),
                    tool.to_uppercase(),
                    tool.chars().filter(|c| c.is_alphanumeric()).count()
                );
            }
            Err(e) => eprintln!("{path}: {e}"),
        }
    }
}

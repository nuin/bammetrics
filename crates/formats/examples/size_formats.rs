//! Binary-growth probe, with the formats: `size_formats_base` plus every entry
//! point the command calls in this crate (both readers, each format of a
//! metrics file, the JSON read back, and the Markdown, JUnit and MultiQC
//! forms of a verdict). See `size_formats_base.rs` for what the size difference
//! measures.

use std::path::Path;

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_formats::{
    CheckLine, Format, Outcome, VerdictSummary, junit, markdown, multiqc, read, render,
};

/// A verdict with one check per outcome, so every branch of the three
/// verdict writers is linked.
fn verdict(sample: &str) -> VerdictSummary {
    let outcomes = [
        Outcome::Pass,
        Outcome::Warn,
        Outcome::Fail,
        Outcome::NotMeasured,
    ];
    VerdictSummary {
        sample: sample.to_owned(),
        status: "fail".to_owned(),
        preset: "generic v1 (000000)".to_owned(),
        checks: outcomes
            .into_iter()
            .map(|outcome| CheckLine {
                rule: "dup.rate".to_owned(),
                label: "Duplication rate".to_owned(),
                metric: "DuplicationMetrics PERCENT_DUPLICATION".to_owned(),
                outcome,
                value: Some(0.5),
                shown: "0.5".to_owned(),
                warn: Some("> 0.2".to_owned()),
                fail: Some("> 0.4".to_owned()),
                reason: Some("above the bound".to_owned()),
            })
            .collect(),
    }
}

fn main() {
    let mut docs = Vec::new();
    for path in std::env::args().skip(1) {
        match MetricsDoc::read(Path::new(&path)) {
            Ok(doc) => println!("{path}: {} classes", doc.classes.len()),
            Err(e) => eprintln!("{path}: {e}"),
        }
        let doc = match std::fs::read_to_string(&path)
            .map_err(|e| e.to_string())
            .and_then(|text| read(&text).map_err(|e| e.to_string()))
        {
            Ok(doc) => doc,
            Err(e) => {
                eprintln!("{path}: {e}");
                continue;
            }
        };
        for format in Format::ALL {
            match render(format, &doc, None) {
                Ok(text) => println!("{path}: {} bytes of {}", text.len(), format.name()),
                Err(e) => eprintln!("{path}: {e}"),
            }
        }
        if let Ok(json) = render(Format::Json, &doc, None) {
            match read(&json) {
                Ok(back) => println!("{path}: {} tables read back", back.tables.len()),
                Err(e) => eprintln!("{path}: {e}"),
            }
        }
        docs.push(doc);
    }
    let verdicts = [verdict("S1")];
    let all: Vec<&_> = docs.iter().collect();
    println!(
        "{} bytes of Markdown, {} of JUnit, {} of MultiQC",
        markdown::write(&all, &verdicts).len(),
        junit::write(&verdicts, true).len(),
        multiqc::write(&verdicts).len()
    );
}

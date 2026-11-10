//! Picard text to JSON to Picard text is byte for byte the first file, on
//! every file Picard 3.5.0 wrote under `tests/data/picard` and on every
//! metrics file a golden record holds (a written file or standard output).

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_formats::{Format, Ty, read, render};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn files_under(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).expect("list directory") {
        let path = entry.expect("directory entry").path();
        if path.is_dir() {
            files_under(&path, out);
        } else {
            out.push(path);
        }
    }
}

/// True for the text of a Picard metrics file: it starts as one (a
/// `StringHeader` line, or the blank line of a file without headers) and
/// has a metrics class or a histogram.
fn is_metrics_file(text: &str) -> bool {
    (text.starts_with("## htsjdk.samtools.metrics.StringHeader\n") || text.starts_with('\n'))
        && text
            .lines()
            .any(|l| l.starts_with("## METRICS CLASS\t") || l.starts_with("## HISTOGRAM\t"))
}

/// The sections of a golden record that hold a file's bytes: standard
/// output and each `## file NAME`.
fn golden_sections(record: &str) -> Vec<(String, String)> {
    let marker = |line: &str| {
        matches!(line, "## argv" | "## exit" | "## stdout" | "## stderr")
            || line.starts_with("## file ")
    };
    let mut sections: Vec<(String, String)> = Vec::new();
    for line in record.split_inclusive('\n') {
        let bare = line.trim_end_matches('\n');
        if marker(bare) {
            sections.push((bare.trim_start_matches("## ").to_owned(), String::new()));
        } else if let Some((_, text)) = sections.last_mut() {
            text.push_str(line);
        }
    }
    sections
        .into_iter()
        .filter(|(name, _)| name == "stdout" || name.starts_with("file "))
        .collect()
}

/// `text` through JSON and back; the error names what differs first.
fn round_trip(text: &str) -> Result<(), String> {
    let doc = read(text).map_err(|e| format!("does not read: {e}"))?;
    let json = render(Format::Json, &doc, None);
    let back = read(&json).map_err(|e| format!("its JSON does not read: {e}"))?;
    let again = render(Format::Picard, &back, None);
    if again == text {
        return Ok(());
    }
    let line = text
        .lines()
        .zip(again.lines())
        .position(|(a, b)| a != b)
        .map_or_else(
            || "the line count".to_owned(),
            |i| {
                format!(
                    "line {}: {:?} became {:?}",
                    i + 1,
                    text.lines().nth(i).unwrap_or_default(),
                    again.lines().nth(i).unwrap_or_default()
                )
            },
        );
    Err(format!("differs at {line}"))
}

#[test]
fn every_picard_fixture_round_trips_through_json() {
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    files.sort();
    let mut checked = 0;
    let mut failures = Vec::new();
    for path in files {
        let text = fs::read_to_string(&path).expect("read fixture");
        if !is_metrics_file(&text) {
            continue;
        }
        checked += 1;
        if let Err(e) = round_trip(&text) {
            failures.push(format!("{}: {e}", path.display()));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!("{checked} Picard files round trip");
    assert!(
        checked >= 8,
        "only {checked} Picard files under tests/data/picard"
    );
}

#[test]
fn every_golden_metrics_file_round_trips_through_json() {
    let mut records = Vec::new();
    files_under(&root().join("tests/golden"), &mut records);
    records.sort();
    let mut checked = 0;
    let mut failures = Vec::new();
    for path in records {
        let record = fs::read_to_string(&path).expect("read golden record");
        for (name, text) in golden_sections(&record) {
            if !is_metrics_file(&text) {
                continue;
            }
            checked += 1;
            if let Err(e) = round_trip(&text) {
                failures.push(format!("{} ({name}): {e}", path.display()));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!("{checked} golden metrics files round trip");
    assert!(
        checked >= 150,
        "only {checked} metrics files in tests/golden"
    );
}

/// The reason for glossary types (docs/07 section 3.1): cells whose text
/// reads as another type than the column has. MEDIAN_INSERT_SIZE prints
/// 310 and is a float.
#[test]
fn glossary_types_differ_from_the_cell_text_on_the_picard_fixture() {
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    let (mut cells, mut retyped) = (0, 0);
    for path in files {
        let text = fs::read_to_string(&path).expect("read fixture");
        if !is_metrics_file(&text) {
            continue;
        }
        let doc = read(&text).expect("reads");
        let raw = bammetrics_core::metrics_file::MetricsDoc::parse(&text).expect("parses");
        for (table, block) in doc.tables.iter().zip(&raw.classes) {
            for (typed, inferred) in table.rows.iter().zip(&block.rows) {
                for ((_, ty), cell) in table.columns.iter().zip(inferred) {
                    use bammetrics_core::metric::Value;
                    let seen = match cell {
                        Value::Int(_) => Ty::Int,
                        Value::Float(_) => Ty::Float,
                        Value::Str(_) | Value::Bool(_) => Ty::Str,
                        Value::Null => continue,
                    };
                    cells += 1;
                    retyped += usize::from(seen != *ty);
                }
                assert_eq!(typed.len(), inferred.len());
            }
        }
    }
    eprintln!("{retyped} of {cells} non-empty cells read as another type than the glossary gives");
    assert!(cells > 100, "{cells} cells");
    assert!(retyped > 20, "{retyped} of {cells} cells change type");
}

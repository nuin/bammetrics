//! Picard text to JSON to Picard text is byte for byte the first file, on
//! every file Picard 3.5.0 wrote under `tests/data/picard` and on every
//! metrics file a golden record holds (a written file or standard output).

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_formats::json::Json;
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

/// A member some object of `json` holds twice. A JSON parser that keeps
/// one value per member (Python's, jq's) would drop the others.
fn repeated_member(json: &Json) -> Option<String> {
    match json {
        Json::Array(items) => items.iter().find_map(repeated_member),
        Json::Object(members) => {
            let mut seen = HashSet::new();
            members
                .iter()
                .find(|(key, _)| !seen.insert(key.as_str()))
                .map(|(key, _)| key.clone())
                .or_else(|| members.iter().find_map(|(_, v)| repeated_member(v)))
        }
        _ => None,
    }
}

/// `text` through JSON and back; the error names what differs first. No
/// object of the JSON or of the JSON Lines form repeats a member.
fn round_trip(text: &str) -> Result<(), String> {
    let doc = read(text).map_err(|e| format!("does not read: {e}"))?;
    let json = render(Format::Json, &doc, None).map_err(|e| format!("no JSON: {e}"))?;
    let value = Json::parse(&json).map_err(|e| format!("its JSON does not parse: {e}"))?;
    if let Some(key) = repeated_member(&value) {
        return Err(format!("its JSON repeats the member {key:?} in one object"));
    }
    let lines = render(Format::Jsonl, &doc, None).map_err(|e| format!("no JSON Lines: {e}"))?;
    for line in lines.lines() {
        let value = Json::parse(line).map_err(|e| format!("a JSON line does not parse: {e}"))?;
        if let Some(key) = repeated_member(&value) {
            return Err(format!("a JSON line repeats the member {key:?}"));
        }
    }
    let back = read(&json).map_err(|e| format!("its JSON does not read: {e}"))?;
    let again = render(Format::Picard, &back, None).map_err(|e| format!("no Picard text: {e}"))?;
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

/// The histogram of `CollectHsMetrics` with two units repeats its three
/// columns. Every series is in the JSON under a member of its own, and the
/// Picard text comes back with the repeated names.
#[test]
fn a_histogram_that_repeats_its_columns_round_trips() {
    let text = "## htsjdk.samtools.metrics.StringHeader\n# CollectHsMetrics INPUT=x.bam\n\n\
                ## METRICS CLASS\tpicard.analysis.directed.HsMetrics\nBAIT_SET\tSAMPLE\n\
                panel\t\npanel\tS1\n\n## HISTOGRAM\tjava.lang.Integer\n\
                coverage_or_base_quality\thigh_quality_coverage_count\tall_bases_coverage_count\t\
                coverage_or_base_quality\thigh_quality_coverage_count\tall_bases_coverage_count\n\
                0\t3\t4\t0\t5\t?\n1\t1\t2\t1\t6\t7\n\n";
    round_trip(text).expect("round trips");
    let doc = read(text).expect("reads");
    let json = Json::parse(&render(Format::Json, &doc, None).expect("JSON")).expect("parses");
    let histogram = json.get("histogram").expect("histogram");
    let row = histogram
        .get("rows")
        .and_then(Json::as_array)
        .expect("rows")[0]
        .as_object()
        .expect("a row object");
    let members: Vec<&str> = row.iter().map(|(key, _)| key.as_str()).collect();
    assert_eq!(
        members,
        [
            "coverage_or_base_quality",
            "high_quality_coverage_count",
            "all_bases_coverage_count",
            "coverage_or_base_quality.1",
            "high_quality_coverage_count.1",
            "all_bases_coverage_count.1",
        ]
    );
    assert_eq!(
        histogram.get("columns"),
        Some(&Json::strings(members.iter().copied()))
    );
    assert_eq!(
        row[4],
        (
            "high_quality_coverage_count.1".to_owned(),
            Json::Number("5".to_owned())
        )
    );
    let lines = render(Format::Jsonl, &doc, None).expect("JSON Lines");
    assert!(
        lines.contains("\"high_quality_coverage_count.1\":6"),
        "{lines}"
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

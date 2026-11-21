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

/// `text` through Parquet and back (docs/07 section 3.1, wave W21): the
/// files `render_files` gives, written under `dir` as `name.parquet` and
/// its companions, read back with `read_path`. The document is the one
/// Picard text gave, column types included, and its Picard text is the
/// first file byte for byte.
#[cfg(feature = "parquet")]
fn parquet_round_trip(text: &str, dir: &Path, name: &str) -> Result<(), String> {
    use bammetrics_formats::{parquet, read_path, render_files};
    let doc = read(text).map_err(|e| format!("does not read: {e}"))?;
    let main = dir.join(format!("{name}.parquet"));
    let files =
        render_files(Format::Parquet, &doc, None).map_err(|e| format!("no Parquet: {e}"))?;
    for file in files {
        let path = file
            .part
            .map_or(main.clone(), |p| parquet::companion(&main, p));
        fs::write(&path, &file.bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    let back = read_path(&main).map_err(|e| format!("its Parquet does not read: {e}"))?;
    // Debug text compares NaN cells as equal, which `==` does not.
    if format!("{back:?}") != format!("{doc:?}") {
        return Err("the document read back differs".to_owned());
    }
    let again = render(Format::Picard, &back, None).map_err(|e| format!("no Picard text: {e}"))?;
    if again != text {
        return Err("its Picard text differs".to_owned());
    }
    Ok(())
}

#[cfg(feature = "parquet")]
fn temporary(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temporary directory");
    dir
}

#[cfg(feature = "parquet")]
#[test]
fn every_picard_fixture_round_trips_through_parquet() {
    let dir = temporary("parquet-picard");
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    files.sort();
    let (mut checked, mut failures) = (0, Vec::new());
    for path in files {
        let text = fs::read_to_string(&path).expect("read fixture");
        if !is_metrics_file(&text) {
            continue;
        }
        checked += 1;
        if let Err(e) = parquet_round_trip(&text, &dir, &checked.to_string()) {
            failures.push(format!("{}: {e}", path.display()));
        }
    }
    let _ = fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!("{checked} Picard files round trip through Parquet");
    assert!(checked >= 8, "only {checked} Picard files");
}

#[cfg(feature = "parquet")]
#[test]
fn every_golden_metrics_file_round_trips_through_parquet() {
    let dir = temporary("parquet-golden");
    let mut records = Vec::new();
    files_under(&root().join("tests/golden"), &mut records);
    records.sort();
    let (mut checked, mut failures) = (0, Vec::new());
    for path in records {
        let record = fs::read_to_string(&path).expect("read golden record");
        for (name, text) in golden_sections(&record) {
            if !is_metrics_file(&text) {
                continue;
            }
            checked += 1;
            if let Err(e) = parquet_round_trip(&text, &dir, &checked.to_string()) {
                failures.push(format!("{} ({name}): {e}", path.display()));
            }
        }
    }
    let _ = fs::remove_dir_all(&dir);
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!("{checked} golden metrics files round trip through Parquet");
    assert!(checked >= 150, "only {checked} metrics files");
}

/// The repeated histogram columns of `CollectHsMetrics` with two units
/// keep their Picard names through `picard_columns`.
#[cfg(feature = "parquet")]
#[test]
fn a_repeated_histogram_column_round_trips_through_parquet() {
    let text = "## htsjdk.samtools.metrics.StringHeader\n# CollectHsMetrics INPUT=x.bam\n\n\
                ## METRICS CLASS\tpicard.analysis.directed.HsMetrics\nBAIT_SET\tSAMPLE\n\
                panel\t\npanel\tS1\n\n## HISTOGRAM\tjava.lang.Integer\n\
                coverage_or_base_quality\thigh_quality_coverage_count\t\
                coverage_or_base_quality\thigh_quality_coverage_count\n\
                0\t3\t0\t?\n1\t1\t1\t7\n\n";
    let dir = temporary("parquet-repeat");
    let result = parquet_round_trip(text, &dir, "hs");
    let _ = fs::remove_dir_all(&dir);
    result.expect("round trips");
}

/// What an outside reader must see in one Parquet file: its rows, and per
/// column its name, its glossary type, its nulls and its NaNs.
#[cfg(feature = "parquet")]
struct Expected {
    rows: usize,
    columns: Vec<(String, Ty, usize, usize)>,
}

#[cfg(feature = "parquet")]
fn expected_table(table: &bammetrics_formats::Table) -> Expected {
    use bammetrics_core::metric::Value;
    let columns = table
        .keys()
        .into_iter()
        .zip(&table.columns)
        .enumerate()
        .map(|(j, (key, (_, ty)))| {
            let cells = || table.rows.iter().map(move |row| &row[j]);
            let nulls = cells().filter(|c| matches!(c, Value::Null)).count();
            let nans = cells()
                .filter(|c| matches!(c, Value::Float(v) if v.is_nan()))
                .count();
            (key, *ty, nulls, nans)
        })
        .collect();
    Expected {
        rows: table.rows.len(),
        columns,
    }
}

#[cfg(feature = "parquet")]
fn expected_long(doc: &bammetrics_formats::Document) -> Expected {
    use bammetrics_core::metric::Value;
    let cells: Vec<&Value> = doc
        .tables
        .iter()
        .chain(&doc.histogram)
        .flat_map(|t| t.rows.iter().flatten())
        .collect();
    let rows = cells.len();
    let empty = cells.iter().filter(|c| matches!(c, Value::Null)).count();
    let sample_nulls = if doc.sample.is_some() { 0 } else { rows };
    let names = ["sample", "class", "row_key", "column", "value"];
    let nulls = [sample_nulls, 0, 0, 0, empty];
    Expected {
        rows,
        columns: names
            .iter()
            .zip(nulls)
            .map(|(name, n)| ((*name).to_owned(), Ty::Str, n, 0))
            .collect(),
    }
}

/// Every Picard fixture and golden metrics file as Parquet under `dir`,
/// with what each file must read as.
#[cfg(feature = "parquet")]
fn parquet_corpus(dir: &Path) -> std::collections::BTreeMap<String, Expected> {
    use bammetrics_formats::{parquet, render_files};
    let mut texts = Vec::new();
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    files_under(&root().join("tests/golden"), &mut files);
    files.sort();
    for path in files {
        let text = fs::read_to_string(&path).expect("read a fixture");
        if path.starts_with(root().join("tests/golden")) {
            texts.extend(golden_sections(&text).into_iter().map(|(_, t)| t));
        } else {
            texts.push(text);
        }
    }
    let mut expected = std::collections::BTreeMap::new();
    for (i, text) in texts.iter().filter(|t| is_metrics_file(t)).enumerate() {
        let mut doc = read(text).expect("reads");
        // Every other file carries a sample, so the long file's sample
        // column is seen both null and filled.
        if i % 2 == 1 {
            doc.sample = Some(format!("S{i}"));
        }
        let main = dir.join(format!("f{i}.parquet"));
        for file in render_files(Format::Parquet, &doc, None).expect("writes") {
            let path = file
                .part
                .map_or(main.clone(), |p| parquet::companion(&main, p));
            fs::write(&path, &file.bytes).expect("write a Parquet file");
            let table = match file.part {
                None => doc.main_table(),
                Some("histogram") => doc.histogram.as_ref(),
                Some(_) => None,
            };
            let want = table.map_or_else(|| expected_long(&doc), expected_table);
            let name = path.file_name().expect("a name").to_string_lossy();
            expected.insert(name.into_owned(), want);
        }
    }
    expected
}

/// DuckDB's name of a glossary type.
#[cfg(feature = "parquet")]
fn duckdb_type(ty: Ty) -> &'static str {
    match ty {
        Ty::Int => "BIGINT",
        Ty::Float => "DOUBLE",
        Ty::Str => "VARCHAR",
    }
}

/// DuckDB reads every Parquet file of every fixture with the glossary's
/// types, every row and every null (docs/07 section 5, wave W21). Skipped
/// when `duckdb` is not on PATH.
#[cfg(feature = "parquet")]
#[test]
fn duckdb_reads_every_parquet_file_with_the_glossary_types() {
    use std::io::Write as _;
    use std::process::{Command, Stdio};
    if Command::new("duckdb").arg("-version").output().is_err() {
        eprintln!("skipped: duckdb is not on PATH");
        return;
    }
    let dir = temporary("parquet-duckdb");
    let expected = parquet_corpus(&dir);
    let mut sql = String::new();
    for (name, want) in &expected {
        let path = dir.join(name);
        let path = path.display();
        sql.push_str(&format!(
            "SELECT '{name}', column_name, column_type FROM (DESCRIBE SELECT * FROM read_parquet('{path}'));\n"
        ));
        let nulls: Vec<String> = want
            .columns
            .iter()
            .map(|(key, ..)| format!("count(*) - count(\"{}\")", key.replace('"', "\"\"")))
            .collect();
        sql.push_str(&format!(
            "SELECT '{name}', 'rows', concat_ws('|', count(*), {}) FROM read_parquet('{path}');\n",
            nulls.join(", ")
        ));
    }
    let mut child = Command::new("duckdb")
        .args(["-csv", "-noheader"])
        .env_remove("AI_AGENT")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("start duckdb");
    // The queries go in from a thread: duckdb answers while it reads, and
    // a full pipe each way would stop both.
    let mut stdin = child.stdin.take().expect("stdin");
    let feeder = std::thread::spawn(move || stdin.write_all(sql.as_bytes()));
    let out = child.wait_with_output().expect("duckdb runs");
    feeder
        .join()
        .expect("the feeder")
        .expect("send the queries");
    let _ = fs::remove_dir_all(&dir);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let mut seen: std::collections::BTreeMap<String, Vec<String>> = Default::default();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let (name, rest) = line.split_once(',').expect("a CSV line");
        seen.entry(name.to_owned())
            .or_default()
            .push(rest.to_owned());
    }
    let mut failures = Vec::new();
    for (name, want) in &expected {
        let mut lines: Vec<String> = want
            .columns
            .iter()
            .map(|(key, ty, ..)| {
                // DuckDB's CSV quotes a field with a comma or a quote.
                let quoted = if key.contains([',', '"', '\'']) {
                    format!("\"{}\"", key.replace('"', "\"\""))
                } else {
                    key.clone()
                };
                format!("{quoted},{}", duckdb_type(*ty))
            })
            .collect();
        let counts: Vec<String> = std::iter::once(want.rows)
            .chain(want.columns.iter().map(|c| c.2))
            .map(|n| n.to_string())
            .collect();
        lines.push(format!("rows,{}", counts.join("|")));
        if seen.get(name) != Some(&lines) {
            failures.push(format!("{name}: {:?} != {lines:?}", seen.get(name)));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!("duckdb read {} Parquet files", expected.len());
}

/// polars and pyarrow read every Parquet file of every fixture with the
/// glossary's types, every row, null and NaN. The interpreter is
/// `BAMMETRICS_PARQUET_PYTHON` (default `python3`); a reader it lacks is
/// skipped: `python3 -m venv target/pyreaders && target/pyreaders/bin/pip
/// install polars pyarrow` gives both.
#[cfg(feature = "parquet")]
#[test]
fn polars_and_pyarrow_read_every_parquet_file_with_the_glossary_types() {
    use bammetrics_formats::json::Json;
    use std::process::Command;
    let python =
        std::env::var("BAMMETRICS_PARQUET_PYTHON").unwrap_or_else(|_| "python3".to_owned());
    let dir = temporary("parquet-python");
    let expected = parquet_corpus(&dir);
    let script = root().join("crates/formats/tests/parquet_readers.py");
    let out = Command::new(&python)
        .arg("-I")
        .arg(&script)
        .arg(&dir)
        .output();
    let _ = fs::remove_dir_all(&dir);
    let Ok(out) = out else {
        eprintln!("skipped: no {python}");
        return;
    };
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let mut failures = Vec::new();
    let mut read = std::collections::BTreeMap::<String, usize>::new();
    for line in String::from_utf8_lossy(&out.stdout).lines() {
        let json = Json::parse(line).expect("a JSON line");
        let reader = json.get("reader").and_then(Json::as_str).expect("a reader");
        if json.get("missing").is_some() {
            eprintln!("skipped: {python} has no {reader}");
            continue;
        }
        *read.entry(reader.to_owned()).or_default() += 1;
        let name = json.get("file").and_then(Json::as_str).expect("a file");
        let want = &expected[name];
        let type_name = |ty: Ty| match (reader, ty) {
            ("polars", Ty::Int) => "Int64",
            ("polars", Ty::Float) => "Float64",
            ("polars", Ty::Str) => "String",
            (_, Ty::Int) => "int64",
            (_, Ty::Float) => "double",
            (_, Ty::Str) => "string",
        };
        let number = |n: usize| Json::Number(n.to_string());
        let columns: Vec<Json> = want
            .columns
            .iter()
            .map(|(key, ty, nulls, nans)| {
                Json::Array(vec![
                    Json::text(key),
                    Json::text(type_name(*ty)),
                    number(*nulls),
                    if *ty == Ty::Float {
                        number(*nans)
                    } else {
                        Json::Null
                    },
                ])
            })
            .collect();
        if json.get("rows") != Some(&number(want.rows))
            || json.get("columns") != Some(&Json::Array(columns.clone()))
        {
            failures.push(format!("{reader} {name}: {line}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    for (reader, files) in &read {
        assert_eq!(*files, expected.len(), "{reader} read every file");
        eprintln!("{reader} read {files} Parquet files");
    }
}

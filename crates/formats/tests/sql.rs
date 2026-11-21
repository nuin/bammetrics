//! The SQL dump read back by outside readers: `sqlite3` loads the dump of
//! every file Picard 3.5.0 wrote under `tests/data/picard` and of every
//! metrics file a golden record holds, and `SELECT` gives back every cell
//! (NaN and infinity as NULL), every column type and every header line.
//! DuckDB does the same on the same files. Each reader is skipped when it
//! is not on PATH.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use bammetrics_core::metric::Value;
use bammetrics_formats::json::Json;
use bammetrics_formats::{Document, Format, Table, Ty, read, render, sql};

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

/// True for the text of a Picard metrics file (as `tests/round_trip.rs`
/// tells one).
fn is_metrics_file(text: &str) -> bool {
    (text.starts_with("## htsjdk.samtools.metrics.StringHeader\n") || text.starts_with('\n'))
        && text
            .lines()
            .any(|l| l.starts_with("## METRICS CLASS\t") || l.starts_with("## HISTOGRAM\t"))
}

/// Standard output and each `## file NAME` section of a golden record.
fn golden_sections(record: &str) -> Vec<String> {
    let marker = |line: &str| {
        matches!(line, "## argv" | "## exit" | "## stdout" | "## stderr")
            || line.starts_with("## file ")
    };
    let mut sections: Vec<(String, String)> = Vec::new();
    for line in record.split_inclusive('\n') {
        let bare = line.trim_end_matches('\n');
        if marker(bare) {
            sections.push((bare.to_owned(), String::new()));
        } else if let Some((_, text)) = sections.last_mut() {
            text.push_str(line);
        }
    }
    sections
        .into_iter()
        .filter(|(name, _)| name == "## stdout" || name.starts_with("## file "))
        .map(|(_, text)| text)
        .collect()
}

/// Every metrics file of the fixtures: (name, text).
fn fixtures() -> Vec<(String, String)> {
    let mut found = Vec::new();
    let mut files = Vec::new();
    files_under(&root().join("tests/data/picard"), &mut files);
    files.sort();
    for path in files {
        let text = fs::read_to_string(&path).expect("read fixture");
        if is_metrics_file(&text) {
            found.push((path.display().to_string(), text));
        }
    }
    let mut records = Vec::new();
    files_under(&root().join("tests/golden"), &mut records);
    records.sort();
    for path in records {
        let record = fs::read_to_string(&path).expect("read golden record");
        for (i, text) in golden_sections(&record).into_iter().enumerate() {
            if is_metrics_file(&text) {
                found.push((format!("{} ({i})", path.display()), text));
            }
        }
    }
    found
}

/// An outside reader: its command and its arguments before the script.
struct Reader {
    program: &'static str,
    args: &'static [&'static str],
    /// The type names its `pragma_table_info` gives for int, float, string.
    types: [&'static str; 3],
}

const SQLITE: Reader = Reader {
    program: "sqlite3",
    args: &[":memory:"],
    types: ["BIGINT", "DOUBLE", "TEXT"],
};

const DUCKDB: Reader = Reader {
    program: "duckdb",
    args: &[":memory:", "-f"],
    types: ["BIGINT", "DOUBLE", "VARCHAR"],
};

impl Reader {
    fn present(&self) -> bool {
        let found = Command::new(self.program)
            .arg("-version")
            .output()
            .is_ok_and(|o| o.status.success());
        if !found {
            eprintln!("{} is not on PATH: skipped", self.program);
        }
        found
    }

    /// The JSON arrays the reader prints for `script`, one per result set
    /// (an empty set prints `[]` or nothing).
    fn run(&self, script: &Path) -> Result<Vec<Json>, String> {
        // DuckDB changes its output when the environment names an agent
        // (AI_AGENT, CLAUDECODE): the reader gets PATH and HOME only.
        let mut command = Command::new(self.program);
        command.args(self.args).env_clear();
        for kept in ["PATH", "HOME"] {
            if let Some(value) = std::env::var_os(kept) {
                command.env(kept, value);
            }
        }
        if self.program == "sqlite3" {
            command.arg(format!(".read {}", script.display()));
        } else {
            command.arg(script);
        }
        let output = command.output().map_err(|e| e.to_string())?;
        let stderr = String::from_utf8_lossy(&output.stderr);
        if !output.status.success() || !stderr.trim().is_empty() {
            return Err(format!("{} failed: {stderr}", self.program));
        }
        let stdout = String::from_utf8(output.stdout).map_err(|e| e.to_string())?;
        let mut arrays = Vec::new();
        let mut pending = String::new();
        for line in stdout.lines() {
            pending.push_str(line);
            pending.push('\n');
            if line.ends_with(']') {
                arrays.push(Json::parse(&pending).map_err(|e| format!("{e}: {pending}"))?);
                pending.clear();
            }
        }
        if !pending.trim().is_empty() {
            return Err(format!("output left over: {pending}"));
        }
        Ok(arrays)
    }
}

fn marker(json: &Json, key: &str) -> Option<String> {
    let rows = json.as_array()?;
    let [row] = rows else { return None };
    row.get(key)?.as_str().map(str::to_owned)
}

/// What a reader gave for one table.
#[derive(Default)]
struct Read {
    rows: Vec<Json>,
    types: Vec<Json>,
}

/// The result sets of the script [`script`] writes, by table name.
fn by_table(arrays: Vec<Json>) -> Result<Vec<(String, Read)>, String> {
    let mut tables: Vec<(String, Read)> = Vec::new();
    let mut into_types = false;
    for array in arrays {
        if let Some(name) = marker(&array, "bammetrics_rows") {
            tables.push((name, Read::default()));
            into_types = false;
        } else if marker(&array, "bammetrics_types").is_some() {
            into_types = true;
        } else {
            let (_, read) = tables.last_mut().ok_or("a result before any marker")?;
            let items = array
                .as_array()
                .ok_or("a result that is no array")?
                .to_vec();
            if into_types {
                read.types = items;
            } else {
                read.rows = items;
            }
        }
    }
    Ok(tables)
}

/// `sql` followed by one query per table: a marker, its rows, a second
/// marker and its declared column types.
fn script(sql: &str, names: &[String]) -> String {
    let mut text = sql.to_owned();
    text.push_str(".mode json\n");
    for name in names {
        let literal = name.replace('\'', "''");
        let identifier = name.replace('"', "\"\"");
        text.push_str(&format!(
            "SELECT '{literal}' AS bammetrics_rows;\nSELECT * FROM \"{identifier}\";\n\
             SELECT '{literal}' AS bammetrics_types;\n\
             SELECT name, type FROM pragma_table_info('{literal}');\n"
        ));
    }
    text
}

/// The cell a reader gives back for `value`: equal for an int and a string,
/// the same `f64` for a float, printed with a point or an exponent, and
/// null for an empty cell, NaN and infinity.
fn same_cell(value: &Value, read: &Json) -> bool {
    match (value, read) {
        (Value::Null, Json::Null) => true,
        (Value::Float(v), Json::Null) => !v.is_finite(),
        (Value::Int(v), Json::Number(text)) => text.parse::<i64>() == Ok(*v),
        (Value::Float(v), Json::Number(text)) => {
            text.contains(['.', 'e', 'E']) && text.parse::<f64>() == Ok(*v)
        }
        (Value::Str(s), Json::Str(text)) => s == text,
        (Value::Bool(b), Json::Str(text)) => text == if *b { "Y" } else { "N" },
        _ => false,
    }
}

fn check_table(
    reader: &Reader,
    name: &str,
    table: &Table,
    read: &Read,
    keys: &[String],
) -> Result<(), String> {
    let types: Vec<String> = table
        .columns
        .iter()
        .map(|(_, ty)| {
            let i = match ty {
                Ty::Int => 0,
                Ty::Float => 1,
                Ty::Str => 2,
            };
            reader.types[i].to_owned()
        })
        .collect();
    let read_types: Vec<String> = read
        .types
        .iter()
        .map(|t| {
            t.get("type")
                .and_then(Json::as_str)
                .unwrap_or("")
                .to_owned()
        })
        .collect();
    if read_types != types {
        return Err(format!("{name}: types {read_types:?}, wanted {types:?}"));
    }
    if read.rows.len() != table.rows.len() {
        return Err(format!(
            "{name}: {} rows, wanted {}",
            read.rows.len(),
            table.rows.len()
        ));
    }
    for (i, (row, got)) in table.rows.iter().zip(&read.rows).enumerate() {
        let members = got.as_object().ok_or("a row that is no object")?;
        if members.len() != keys.len() {
            return Err(format!("{name} row {i}: {} cells", members.len()));
        }
        for (j, ((key, cell), want)) in members.iter().zip(row).enumerate() {
            if key != &keys[j] || !same_cell(want, cell) {
                return Err(format!(
                    "{name} row {i} column {key}: {cell:?}, wanted {want:?}"
                ));
            }
        }
    }
    Ok(())
}

/// `doc` through `reader`: every table, the headers and the catalogue.
fn load(reader: &Reader, doc: &Document, dir: &Path, index: usize) -> Result<(), String> {
    let sql = render(Format::Sql, doc, None).map_err(|e| e.to_string())?;
    let names = sql::table_names(doc);
    let mut all = vec!["headers".to_owned(), "tables".to_owned()];
    all.extend(names.iter().cloned());
    let path = dir.join(format!("{index}.sql"));
    fs::write(&path, script(&sql, &all)).map_err(|e| e.to_string())?;
    let tables = by_table(reader.run(&path)?)?;
    fs::remove_file(&path).map_err(|e| e.to_string())?;
    let read = |name: &str| {
        tables
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, r)| r)
            .ok_or_else(|| format!("no result for {name}"))
    };
    let headers: Vec<Option<&str>> = read("headers")?
        .rows
        .iter()
        .map(|r| r.get("text").and_then(Json::as_str))
        .collect();
    let want: Vec<Option<&str>> = doc.headers.iter().map(|h| Some(h.as_str())).collect();
    if headers != want {
        return Err(format!("headers {headers:?}, wanted {want:?}"));
    }
    let listed: Vec<&str> = read("tables")?
        .rows
        .iter()
        .filter_map(|r| r.get("name").and_then(Json::as_str))
        .collect();
    if listed != names {
        return Err(format!("tables lists {listed:?}, wanted {names:?}"));
    }
    let placed = doc
        .tables
        .iter()
        .chain(&doc.histogram)
        .filter(|t| !t.columns.is_empty());
    for (name, table) in names.iter().zip(placed) {
        let read = read(name)?;
        let keys: Vec<String> = read
            .types
            .iter()
            .map(|t| {
                t.get("name")
                    .and_then(Json::as_str)
                    .unwrap_or("")
                    .to_owned()
            })
            .collect();
        check_table(reader, name, table, read, &keys)?;
    }
    Ok(())
}

fn every_fixture_loads(reader: &Reader) {
    if !reader.present() {
        return;
    }
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("sql-{}", reader.program));
    fs::create_dir_all(&dir).expect("scratch directory");
    let (mut checked, mut cells) = (0, 0);
    let mut failures = Vec::new();
    for (index, (name, text)) in fixtures().into_iter().enumerate() {
        let doc = read(&text).expect("reads");
        cells += doc
            .tables
            .iter()
            .chain(&doc.histogram)
            .map(|t| t.rows.len() * t.columns.len())
            .sum::<usize>();
        checked += 1;
        if let Err(e) = load(reader, &doc, &dir, index) {
            failures.push(format!("{name}: {e}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    eprintln!(
        "{} loads {checked} metrics files, {cells} cells",
        reader.program
    );
    assert!(checked >= 160, "only {checked} metrics files");
}

#[test]
fn sqlite3_loads_every_fixture_and_gives_back_every_cell() {
    every_fixture_loads(&SQLITE);
}

#[test]
fn duckdb_loads_every_fixture_and_gives_back_every_cell() {
    every_fixture_loads(&DUCKDB);
}

/// Strings with a quote, a carriage return, a NUL and a line that reads
/// as a shell command; a class and columns named like the catalogue.
#[test]
fn sqlite3_keeps_every_character() {
    if !SQLITE.present() {
        return;
    }
    let cells = [
        "it's",
        "a\r\nb",
        "\r",
        "x\0y",
        "line\n.quit\nGO\n/",
        "\"q\"; DROP TABLE x; --",
        "é ü 🧬",
    ];
    let doc = Document {
        headers: vec!["a'\r\nb".to_owned()],
        tables: vec![Table {
            class: "x.Headers".to_owned(),
            columns: vec![
                ("S".to_owned(), Ty::Str),
                ("s".to_owned(), Ty::Str),
                ("a\"b".to_owned(), Ty::Str),
            ],
            rows: cells
                .iter()
                .map(|c| {
                    vec![
                        Value::Str((*c).to_owned()),
                        Value::Null,
                        Value::Str(format!("{c}!")),
                    ]
                })
                .collect(),
        }],
        ..Document::default()
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sql-characters");
    fs::create_dir_all(&dir).expect("scratch directory");
    let sql = render(Format::Sql, &doc, None).expect("renders");
    let path = dir.join("chars.sql");
    let query = "SELECT hex(text) AS h FROM \"headers\";\n\
                 SELECT hex(\"S\") AS h, \"s.1\" AS n, hex(\"a\"\"b\") AS q FROM \"Headers_2\";\n";
    fs::write(&path, format!("{sql}.mode json\n{query}")).expect("write");
    let arrays = SQLITE.run(&path).expect("loads");
    let hex = |text: &str| -> String { text.bytes().map(|b| format!("{b:02X}")).collect() };
    assert_eq!(
        arrays[0].as_array().expect("headers")[0].get("h"),
        Some(&Json::text(&hex("a'\r\nb")))
    );
    let rows = arrays[1].as_array().expect("rows");
    assert_eq!(rows.len(), cells.len());
    for (row, cell) in rows.iter().zip(cells) {
        assert_eq!(row.get("h"), Some(&Json::text(&hex(cell))), "{cell:?}");
        assert_eq!(row.get("n"), Some(&Json::Null));
        assert_eq!(row.get("q"), Some(&Json::text(&hex(&format!("{cell}!")))));
    }
}

/// Classes named `sqlite_master` and `SQLite_x` load as `_sqlite_master`
/// and `_SQLite_x`: SQLite refuses a table whose name starts with `sqlite_`.
#[test]
fn sqlite3_loads_a_class_named_like_its_own_tables() {
    if !SQLITE.present() {
        return;
    }
    let class = |name: &str, n: i64| Table {
        class: format!("picard.analysis.{name}"),
        columns: vec![("N".to_owned(), Ty::Int)],
        rows: vec![vec![Value::Int(n)]],
    };
    let doc = Document {
        tables: vec![class("sqlite_master", 1), class("SQLite_x", 2)],
        ..Document::default()
    };
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("sql-reserved");
    fs::create_dir_all(&dir).expect("scratch directory");
    let sql = render(Format::Sql, &doc, None).expect("renders");
    let path = dir.join("reserved.sql");
    let query = "SELECT t.name AS name, \
                 (SELECT count(*) FROM sqlite_schema s WHERE s.name = t.name) AS made \
                 FROM \"tables\" t;\n\
                 SELECT \"N\" AS n FROM \"_sqlite_master\" \
                 UNION ALL SELECT \"N\" FROM \"_SQLite_x\";\n";
    fs::write(&path, format!("{sql}.mode json\n{query}")).expect("write");
    let arrays = SQLITE.run(&path).expect("loads");
    let pairs: Vec<_> = arrays[0]
        .as_array()
        .expect("tables")
        .iter()
        .map(|row| (row.get("name").cloned(), row.get("made").cloned()))
        .collect();
    let number = |text: &str| Some(Json::parse(text).expect("a number"));
    assert_eq!(
        pairs,
        [
            (Some(Json::text("_sqlite_master")), number("1")),
            (Some(Json::text("_SQLite_x")), number("1")),
        ]
    );
    let values: Vec<_> = arrays[1]
        .as_array()
        .expect("rows")
        .iter()
        .map(|row| row.get("n").cloned())
        .collect();
    assert_eq!(values, [number("1"), number("2")]);
}

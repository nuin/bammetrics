//! `docs/schema/metrics-1.json` describes what the JSON writer writes:
//! every Picard fixture, as JSON, is valid under it. The validator below
//! knows the keywords the schema uses and fails on any other.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_formats::json::Json;
use bammetrics_formats::{Format, read, render};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn type_matches(name: &str, value: &Json) -> bool {
    match (name, value) {
        ("null", Json::Null)
        | ("string", Json::Str(_))
        | ("number", Json::Number(_))
        | ("array", Json::Array(_))
        | ("object", Json::Object(_))
        | ("boolean", Json::Bool(_)) => true,
        ("integer", Json::Number(text)) => text.parse::<i64>().is_ok(),
        _ => false,
    }
}

/// The first way `value` breaks `schema`, with the path to it.
fn validate(schema: &Json, value: &Json, root: &Json, path: &str) -> Result<(), String> {
    let members = schema.as_object().ok_or("a schema is not an object")?;
    for (keyword, rule) in members {
        match keyword.as_str() {
            "$schema" | "$id" | "$defs" | "title" | "description" => {}
            "$ref" => {
                let name = rule
                    .as_str()
                    .and_then(|r| r.strip_prefix("#/$defs/"))
                    .ok_or("a $ref outside #/$defs")?;
                let target = root
                    .get("$defs")
                    .and_then(|d| d.get(name))
                    .ok_or_else(|| format!("no $defs/{name}"))?;
                validate(target, value, root, path)?;
            }
            "const" => {
                if rule != value {
                    return Err(format!("{path}: not the constant"));
                }
            }
            "enum" => {
                if !rule.as_array().is_some_and(|items| items.contains(value)) {
                    return Err(format!("{path}: not in the enum"));
                }
            }
            "type" => {
                let names: Vec<&str> = match rule {
                    Json::Str(name) => vec![name],
                    Json::Array(names) => names.iter().filter_map(Json::as_str).collect(),
                    _ => return Err("a type that is no name".to_owned()),
                };
                if !names.iter().any(|n| type_matches(n, value)) {
                    return Err(format!("{path}: not of type {}", names.join(" or ")));
                }
            }
            "required" => {
                for key in rule.as_array().unwrap_or_default() {
                    let key = key.as_str().unwrap_or_default();
                    if value.get(key).is_none() {
                        return Err(format!("{path}: no {key}"));
                    }
                }
            }
            "properties" => {
                for (key, sub) in rule.as_object().unwrap_or_default() {
                    if let Some(member) = value.get(key) {
                        validate(sub, member, root, &format!("{path}.{key}"))?;
                    }
                }
            }
            "additionalProperties" => {
                let known: Vec<&str> = schema
                    .get("properties")
                    .and_then(Json::as_object)
                    .unwrap_or_default()
                    .iter()
                    .map(|(k, _)| k.as_str())
                    .collect();
                for (key, member) in value.as_object().unwrap_or_default() {
                    if known.contains(&key.as_str()) {
                        continue;
                    }
                    match rule {
                        Json::Bool(false) => return Err(format!("{path}: {key} is not allowed")),
                        Json::Bool(true) => {}
                        sub => validate(sub, member, root, &format!("{path}.{key}"))?,
                    }
                }
            }
            "items" => {
                for (i, item) in value.as_array().unwrap_or_default().iter().enumerate() {
                    validate(rule, item, root, &format!("{path}[{i}]"))?;
                }
            }
            "oneOf" => {
                let passing = rule
                    .as_array()
                    .unwrap_or_default()
                    .iter()
                    .filter(|sub| validate(sub, value, root, path).is_ok())
                    .count();
                if passing != 1 {
                    return Err(format!("{path}: {passing} of the oneOf schemas match"));
                }
            }
            "minimum" => {
                let number = |j: &Json| match j {
                    Json::Number(text) => text.parse::<f64>().ok(),
                    _ => None,
                };
                if let (Some(least), Some(n)) = (number(rule), number(value))
                    && n < least
                {
                    return Err(format!("{path}: below the minimum"));
                }
            }
            other => return Err(format!("the validator does not know {other:?}")),
        }
    }
    Ok(())
}

fn schema() -> Json {
    let text =
        fs::read_to_string(root().join("docs/schema/metrics-1.json")).expect("read the schema");
    Json::parse(&text).expect("the schema is JSON")
}

fn picard_files() -> Vec<PathBuf> {
    fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
        for entry in fs::read_dir(dir).expect("list directory") {
            let path = entry.expect("directory entry").path();
            if path.is_dir() {
                walk(&path, out);
            } else {
                out.push(path);
            }
        }
    }
    let mut files = Vec::new();
    walk(&root().join("tests/data/picard"), &mut files);
    files.sort();
    files
}

#[test]
fn the_writer_output_is_valid_under_the_committed_schema() {
    let schema = schema();
    let mut checked = 0;
    for path in picard_files() {
        let text = fs::read_to_string(&path).expect("read fixture");
        let Ok(mut doc) = read(&text) else { continue };
        doc.sample = Some("S1".to_owned());
        let json = Json::parse(&render(Format::Json, &doc, None)).expect("valid JSON");
        validate(&schema, &json, &schema, "$")
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        checked += 1;
    }
    assert!(checked >= 8, "{checked} files");
}

#[test]
fn the_schema_refuses_what_the_reader_refuses() {
    let schema = schema();
    let valid = r#"{"schema": "bammetrics.metrics/1", "headers": [], "class": null, "columns": ["A"],
        "types": ["float"], "metrics": [{"A": null}],
        "non_finite": [{"row": 0, "column": "A", "value": "NaN"}], "histogram": null}"#;
    let json = Json::parse(valid).expect("JSON");
    validate(&schema, &json, &schema, "$").expect("valid");
    assert!(read(valid).is_ok());
    for (from, to) in [
        ("bammetrics.metrics/1", "bammetrics.metrics/2"),
        ("\"float\"", "\"double\""),
        ("\"NaN\"", "\"nan\""),
        ("\"headers\": [], ", ""),
        ("{\"A\": null}", "{\"A\": [1]}"),
    ] {
        let text = valid.replace(from, to);
        let json = Json::parse(&text).expect("JSON");
        assert!(validate(&schema, &json, &schema, "$").is_err(), "{to}");
    }
}

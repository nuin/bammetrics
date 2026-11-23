//! A JSON Schema validator for the keywords the committed schemas of
//! `docs/schema` use; any other keyword is an error. Shared by the schema
//! tests of `crates/formats` and `crates/cli` through `#[path]`.

use bammetrics_formats::json::Json;

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
pub fn validate(schema: &Json, value: &Json, root: &Json, path: &str) -> Result<(), String> {
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

//! docs/verdict.md lists every rule of the shipped presets, one table per
//! preset between `<!-- preset NAME -->` and `<!-- end preset NAME -->`.
//! Each table is rendered here from the compiled preset, so the doc cannot
//! drift from `presets/*.toml`: a changed bound, rule, note or hash fails
//! this test until the doc is rewritten with
//!
//! ```text
//! UPDATE_DOCS=1 cargo test -p bammetrics-verdict --test doc_presets
//! ```

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_verdict::{Preset, names};

fn doc_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/verdict.md")
}

/// The doc block of `preset`: one summary line and one table row per rule.
fn block(preset: &Preset) -> String {
    let mut out = format!(
        "Version {}, {} rules, note `{}`, hash `{}` (the report prints `{}`).\n\n",
        preset.version,
        preset.rules.len(),
        preset.note,
        preset.hash,
        &preset.hash[..6]
    );
    out.push_str("| Rule | Label | Metric | Rows | Warn | Fail | Plot |\n");
    out.push_str("|---|---|---|---|---|---|---|\n");
    let code = |text: Option<String>| text.map_or_else(String::new, |t| format!("`{t}`"));
    for rule in &preset.rules {
        let mut rows = Vec::new();
        if rule.level.name() != "ALL_READS" {
            rows.push(format!("`{}`", rule.level.name()));
        }
        if let Some(select) = rule.select_text() {
            rows.push(format!("`{select}`"));
        }
        out.push_str(&format!(
            "| `{}` | {} | `{}.{}` | {} | {} | {} | {} |\n",
            rule.id,
            rule.label,
            rule.class,
            rule.metric_text(),
            rows.join(", "),
            code(rule.warn.as_ref().map(ToString::to_string)),
            code(rule.fail.as_ref().map(ToString::to_string)),
            code(rule.plot.clone()),
        ));
    }
    out
}

/// `doc` with the block of `name` replaced by `body`, or `None` when the
/// markers are missing.
fn replace_block(doc: &str, name: &str, body: &str) -> Option<String> {
    let open = format!("<!-- preset {name} -->\n");
    let close = format!("<!-- end preset {name} -->");
    let start = doc.find(&open)? + open.len();
    let end = start + doc[start..].find(&close)?;
    Some(format!("{}{body}{}", &doc[..start], &doc[end..]))
}

#[test]
fn verdict_md_lists_every_rule_of_every_shipped_preset() {
    let path = doc_path();
    let original = fs::read_to_string(&path).expect("docs/verdict.md");
    let mut doc = original.clone();
    for name in names() {
        let preset = Preset::builtin(name)
            .expect("shipped")
            .expect("shipped preset parses");
        assert_eq!(preset.note, expected_note(name), "{name}");
        doc = replace_block(&doc, name, &block(&preset))
            .unwrap_or_else(|| panic!("docs/verdict.md has no <!-- preset {name} --> block"));
    }
    if std::env::var_os("UPDATE_DOCS").is_some() {
        fs::write(&path, &doc).expect("write docs/verdict.md");
        return;
    }
    if doc != original {
        let line = original
            .lines()
            .zip(doc.lines())
            .position(|(a, b)| a != b)
            .map_or(0, |i| i + 1);
        panic!(
            "docs/verdict.md differs from the compiled presets at line {line}; \
             run UPDATE_DOCS=1 cargo test -p bammetrics-verdict --test doc_presets and review the diff"
        );
    }
}

/// The `Shipped presets` summary table of docs/verdict.md is hand-written:
/// its `Rules` column must equal each preset's rule count, and the counts
/// in its `Reads` column (`alignment summary (3)`, one rule when no count
/// is given) must add up to it.
#[test]
fn verdict_md_summary_counts_match_the_compiled_presets() {
    let doc = fs::read_to_string(doc_path()).expect("docs/verdict.md");
    for name in names() {
        let preset = Preset::builtin(name)
            .expect("shipped")
            .expect("shipped preset parses");
        let prefix = format!("| `{name}` | ");
        let row = doc
            .lines()
            .find(|line| line.starts_with(&prefix))
            .unwrap_or_else(|| panic!("docs/verdict.md has no summary row for {name}"));
        let cells: Vec<&str> = row.split(" | ").collect();
        let rules: usize = cells[1].trim().parse().expect("a rule count");
        assert_eq!(rules, preset.rules.len(), "{name}: Rules column");
        let read: usize = cells[2]
            .split(", ")
            .map(|item| {
                item.rsplit_once(" (")
                    .and_then(|(_, n)| n.strip_suffix(')'))
                    .map_or(1, |n| n.parse().expect("a count"))
            })
            .sum();
        assert_eq!(read, preset.rules.len(), "{name}: Reads column");
    }
}

/// The shipped note: `clinical-panel` was reviewed (decision 2 of docs/05
/// section 8) and its version 3 gap and sex rules are starting values
/// (docs/06); the other presets still carry starting values.
fn expected_note(name: &str) -> &'static str {
    if name == "clinical-panel" {
        "reviewed defaults; gap and sex rules are starting values"
    } else {
        "starting values"
    }
}

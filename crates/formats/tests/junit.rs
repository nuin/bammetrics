//! The JUnit XML holds to the shape CI systems read: `testsuites` over
//! `testsuite` over `testcase`, the counts on each suite equal to the
//! cases under it, every element closed. `xmllint`, when it is on PATH,
//! parses the same text.

use std::collections::BTreeMap;
use std::io::Write as _;
use std::process::{Command, Stdio};

use bammetrics_formats::{CheckLine, Outcome, VerdictSummary, junit};

fn check(rule: &str, outcome: Outcome) -> CheckLine {
    CheckLine {
        rule: rule.to_owned(),
        label: format!("{rule} in words"),
        metric: "HsMetrics PCT_TARGET_BASES_20X".to_owned(),
        outcome,
        value: Some(0.5),
        shown: "0.5".to_owned(),
        warn: Some("< 0.9".to_owned()),
        fail: Some("< 0.8".to_owned()),
        reason: Some(format!("{rule}: 0.5 < 0.8 & \"quoted\" <tag>\nsecond line")),
    }
}

fn verdicts() -> Vec<VerdictSummary> {
    let sample = |name: &str, outcomes: &[Outcome]| VerdictSummary {
        sample: name.to_owned(),
        status: "fail".to_owned(),
        preset: "generic v3 (0123ab)".to_owned(),
        checks: outcomes
            .iter()
            .enumerate()
            .map(|(i, o)| check(&format!("rule.{i}"), *o))
            .collect(),
    };
    vec![
        sample(
            "NA12878 & <co>",
            &[
                Outcome::Fail,
                Outcome::Warn,
                Outcome::Pass,
                Outcome::NotMeasured,
            ],
        ),
        sample("S2", &[Outcome::Pass, Outcome::Pass]),
        sample("empty", &[]),
    ]
}

/// One tag: its name, its attributes, and whether it opens, closes or both.
struct Tag {
    name: String,
    attributes: BTreeMap<String, String>,
    opens: bool,
    closes: bool,
}

/// The tags of `xml` in order. Panics on text that is not the subset of
/// XML the writer produces: quoted attributes, no comments, no CDATA.
fn tags(xml: &str) -> Vec<Tag> {
    let body = xml
        .strip_prefix("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
        .expect("the XML declaration");
    let mut tags = Vec::new();
    let mut rest = body;
    while let Some(start) = rest.find('<') {
        assert!(!rest[..start].contains('>'), "a bare > in text");
        let end = rest[start..].find('>').expect("a closed tag") + start;
        let inner = &rest[start + 1..end];
        rest = &rest[end + 1..];
        let closes_only = inner.starts_with('/');
        let self_closing = inner.ends_with('/');
        let inner = inner.trim_start_matches('/').trim_end_matches('/');
        let (name, mut attribute_text) = inner.split_once(' ').unwrap_or((inner, ""));
        let mut attributes = BTreeMap::new();
        while let Some((key, after)) = attribute_text.split_once("=\"") {
            let (value, after) = after.split_once('"').expect("a closed attribute value");
            assert!(
                !value.contains('<') && !value.contains('\n'),
                "raw text in {key}"
            );
            assert!(
                attributes
                    .insert(key.trim().to_owned(), value.to_owned())
                    .is_none()
            );
            attribute_text = after;
        }
        assert!(
            attribute_text.trim().is_empty(),
            "text after the attributes"
        );
        tags.push(Tag {
            name: name.to_owned(),
            attributes,
            opens: !closes_only,
            closes: closes_only || self_closing,
        });
    }
    tags
}

/// The children each element may have, and the attributes it must have.
const SHAPE: &[(&str, &[&str], &[&str])] = &[
    (
        "testsuites",
        &["testsuite"],
        &["name", "tests", "failures", "errors", "skipped"],
    ),
    (
        "testsuite",
        &["properties", "testcase"],
        &["name", "tests", "failures", "errors", "skipped"],
    ),
    ("properties", &["property"], &[]),
    ("property", &[], &["name", "value"]),
    (
        "testcase",
        &["failure", "skipped", "system-out"],
        &["classname", "name"],
    ),
    ("failure", &[], &["type", "message"]),
    ("skipped", &[], &["message"]),
    ("system-out", &[], &[]),
];

fn check_shape(xml: &str) {
    let mut stack: Vec<String> = Vec::new();
    // Per open suite: (tests, failures, skipped) declared and counted.
    let mut declared = Vec::new();
    let mut counted = Vec::new();
    let mut roots = 0;
    for tag in tags(xml) {
        if tag.opens {
            let (_, _, required) = SHAPE
                .iter()
                .find(|(name, _, _)| *name == tag.name)
                .unwrap_or_else(|| panic!("unknown element {}", tag.name));
            for key in *required {
                assert!(
                    tag.attributes.contains_key(*key),
                    "{} has no {key}",
                    tag.name
                );
            }
            match stack.last() {
                None => {
                    assert_eq!(tag.name, "testsuites");
                    roots += 1;
                }
                Some(parent) => {
                    let (_, children, _) = SHAPE
                        .iter()
                        .find(|(name, _, _)| name == parent)
                        .expect("a known parent");
                    assert!(
                        children.contains(&tag.name.as_str()),
                        "{} in {parent}",
                        tag.name
                    );
                }
            }
            let number = |key: &str| -> usize { tag.attributes[key].parse().expect("a count") };
            match tag.name.as_str() {
                "testsuites" | "testsuite" => {
                    assert_eq!(number("errors"), 0);
                    declared.push((number("tests"), number("failures"), number("skipped")));
                    counted.push((0, 0, 0));
                }
                "testcase" => counted.iter_mut().for_each(|c| c.0 += 1),
                "failure" => counted.iter_mut().for_each(|c| c.1 += 1),
                "skipped" => counted.iter_mut().for_each(|c| c.2 += 1),
                _ => {}
            }
            stack.push(tag.name.clone());
        }
        if tag.closes {
            assert_eq!(
                stack.pop().as_deref(),
                Some(tag.name.as_str()),
                "closing tag"
            );
            if matches!(tag.name.as_str(), "testsuites" | "testsuite") {
                assert_eq!(declared.pop(), counted.pop(), "counts of {}", tag.name);
            }
        }
    }
    assert!(stack.is_empty(), "open elements: {stack:?}");
    assert_eq!(roots, 1);
}

fn xmllint_accepts(xml: &str) -> Option<bool> {
    let mut child = Command::new("xmllint")
        .args(["--noout", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .ok()?;
    child
        .stdin
        .take()
        .expect("a piped stdin")
        .write_all(xml.as_bytes())
        .ok()?;
    Some(child.wait().ok()?.success())
}

#[test]
fn the_document_has_the_junit_shape() {
    for warn_fails in [false, true] {
        let xml = junit::write(&verdicts(), warn_fails);
        check_shape(&xml);
        match xmllint_accepts(&xml) {
            Some(ok) => assert!(ok, "xmllint refuses the document"),
            None => eprintln!("skipped: xmllint is not on PATH"),
        }
    }
    check_shape(&junit::write(&[], false));
}

#[test]
fn the_shape_check_fails_on_a_planted_error() {
    let xml = junit::write(&verdicts(), false);
    let wrong_count = xml.replacen("tests=\"4\"", "tests=\"5\"", 1);
    assert!(std::panic::catch_unwind(|| check_shape(&wrong_count)).is_err());
    let unclosed = xml.replacen("  </testsuite>\n", "", 1);
    assert!(std::panic::catch_unwind(|| check_shape(&unclosed)).is_err());
}

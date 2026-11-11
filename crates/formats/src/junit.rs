//! JUnit XML for `check`: one `<testsuite>` per sample, one `<testcase>`
//! per rule, the form CI systems draw as a test report.
//!
//! A failed rule is a `<failure>`, a rule that could not be measured a
//! `<skipped>`. A warning passes with its reason in `<system-out>`, or is
//! a `<failure type="warn">` when the gate is `--fail-on warn`.

use crate::verdict::{CheckLine, Outcome, VerdictSummary};

/// `text` for XML content or a quoted attribute. A character XML 1.0 does
/// not allow is written as U+FFFD; a line break or tab in an attribute as
/// a character reference, so a parser does not fold it into a space.
fn escape(text: &str, attribute: bool) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            '\n' | '\r' | '\t' if attribute => out.push_str(&format!("&#{};", c as u32)),
            '\n' | '\r' | '\t' => out.push(c),
            c if (c as u32) < 0x20 || matches!(c, '\u{FFFE}' | '\u{FFFF}') => out.push('\u{FFFD}'),
            c => out.push(c),
        }
    }
    out
}

fn attribute(text: &str) -> String {
    escape(text, true)
}

/// `rule: value, warn BOUND, fail BOUND`, the line `check` prints.
fn detail(check: &CheckLine) -> String {
    let mut text = format!("{}: {} = {}", check.rule, check.metric, check.shown);
    for (name, bound) in [("warn", &check.warn), ("fail", &check.fail)] {
        if let Some(bound) = bound {
            text.push_str(&format!(", {name} {bound}"));
        }
    }
    text
}

fn testcase(out: &mut String, sample: &str, check: &CheckLine, warn_fails: bool) {
    out.push_str(&format!(
        "    <testcase classname=\"bammetrics.{}\" name=\"{}\"",
        attribute(sample),
        attribute(&check.rule)
    ));
    let reason = check.reason.as_deref().unwrap_or(&check.label);
    let body = escape(&detail(check), false);
    match check.outcome {
        Outcome::Pass => out.push_str("/>\n"),
        Outcome::NotMeasured => out.push_str(&format!(
            ">\n      <skipped message=\"{}\"/>\n    </testcase>\n",
            attribute(reason)
        )),
        Outcome::Warn if !warn_fails => out.push_str(&format!(
            ">\n      <system-out>WARN {}: {body}</system-out>\n    </testcase>\n",
            escape(reason, false)
        )),
        Outcome::Warn | Outcome::Fail => {
            let kind = if check.outcome == Outcome::Fail {
                "fail"
            } else {
                "warn"
            };
            out.push_str(&format!(
                ">\n      <failure type=\"{kind}\" message=\"{}\">{body}</failure>\n    </testcase>\n",
                attribute(reason)
            ));
        }
    }
}

/// `verdicts` as a JUnit XML document. `warn_fails` is the gate `--fail-on
/// warn`: a warning counts as a failure.
pub fn write(verdicts: &[VerdictSummary], warn_fails: bool) -> String {
    let failures = |v: &VerdictSummary| {
        v.count(Outcome::Fail)
            + if warn_fails {
                v.count(Outcome::Warn)
            } else {
                0
            }
    };
    let total = |pick: &dyn Fn(&VerdictSummary) -> usize| verdicts.iter().map(pick).sum::<usize>();
    let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(&format!(
        "<testsuites name=\"bammetrics check\" tests=\"{}\" failures=\"{}\" errors=\"0\" skipped=\"{}\">\n",
        total(&|v| v.checks.len()),
        total(&failures),
        total(&|v| v.count(Outcome::NotMeasured)),
    ));
    for verdict in verdicts {
        out.push_str(&format!(
            "  <testsuite name=\"{}\" tests=\"{}\" failures=\"{}\" errors=\"0\" skipped=\"{}\">\n",
            attribute(&verdict.sample),
            verdict.checks.len(),
            failures(verdict),
            verdict.count(Outcome::NotMeasured),
        ));
        out.push_str("    <properties>\n");
        for (name, value) in [("preset", &verdict.preset), ("verdict", &verdict.status)] {
            out.push_str(&format!(
                "      <property name=\"{name}\" value=\"{}\"/>\n",
                attribute(value)
            ));
        }
        out.push_str("    </properties>\n");
        for check in &verdict.checks {
            testcase(&mut out, &verdict.sample, check, warn_fails);
        }
        out.push_str("  </testsuite>\n");
    }
    out.push_str("</testsuites>\n");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verdict::tests::verdict;

    #[test]
    fn one_testcase_per_rule() {
        let text = write(&[verdict()], false);
        assert!(text.starts_with(
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <testsuites name=\"bammetrics check\" tests=\"4\" failures=\"1\" errors=\"0\" skipped=\"1\">\n  \
             <testsuite name=\"NA|12878\" tests=\"4\" failures=\"1\" errors=\"0\" skipped=\"1\">\n"
        ));
        assert_eq!(text.matches("<testcase ").count(), 4);
        assert!(text.contains(
            "<testcase classname=\"bammetrics.NA|12878\" name=\"dup.rate\">\n      \
             <failure type=\"fail\" message=\"dup.rate &lt;is&gt; &quot;off&quot; &amp; more\">\
             dup.rate: DuplicationMetrics PERCENT_DUPLICATION = 0.31, warn &gt; 0.2, fail &gt; 0.3</failure>"
        ));
        assert!(text.contains("<system-out>WARN q30 &lt;is&gt;"));
        assert!(text.contains("name=\"insert.median\"/>\n"));
        assert!(text.contains("<skipped message=\"wgs.depth &lt;is&gt;"));
        assert!(text.ends_with("  </testsuite>\n</testsuites>\n"));
    }

    #[test]
    fn a_warning_fails_under_the_warn_gate() {
        let text = write(&[verdict(), verdict()], true);
        assert!(text.contains("<testsuites name=\"bammetrics check\" tests=\"8\" failures=\"4\""));
        assert!(text.contains("<failure type=\"warn\""));
        assert!(!text.contains("<system-out>"));
    }

    #[test]
    fn text_xml_cannot_hold_is_replaced() {
        assert_eq!(escape("a\u{1}b\nc'", true), "a\u{FFFD}b&#10;c&apos;");
        assert_eq!(escape("a\tb\nc", false), "a\tb\nc");
        assert_eq!(
            write(&[], false),
            "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
             <testsuites name=\"bammetrics check\" tests=\"0\" failures=\"0\" errors=\"0\" skipped=\"0\">\n\
             </testsuites>\n"
        );
    }
}

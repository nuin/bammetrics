//! What the Markdown, JUnit and MultiQC writers need of a verdict.
//!
//! This crate depends on the core only, so the caller copies these fields
//! out of its `bammetrics_verdict::Verdict` (`crates/cli/src/format.rs`).

/// How one rule ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    Pass,
    Warn,
    Fail,
    /// The rule's file, row, column or value is absent.
    NotMeasured,
}

impl Outcome {
    /// `PASS`, `WARN`, `FAIL` or `NOT MEASURED`, as `check` prints it.
    pub fn word(self) -> &'static str {
        match self {
            Outcome::Pass => "PASS",
            Outcome::Warn => "WARN",
            Outcome::Fail => "FAIL",
            Outcome::NotMeasured => "NOT MEASURED",
        }
    }
}

/// One rule's outcome on a sample.
#[derive(Clone, Debug, PartialEq)]
pub struct CheckLine {
    /// The rule id, `dup.rate`.
    pub rule: String,
    /// The rule in words.
    pub label: String,
    /// The metrics class and what is judged in it: `HsMetrics
    /// PCT_TARGET_BASES_20X`.
    pub metric: String,
    pub outcome: Outcome,
    /// The number judged; `None` when not measured or judged as text.
    pub value: Option<f64>,
    /// The value as `check` prints it: the number, the cell text, `empty`,
    /// or `?` when not measured.
    pub shown: String,
    /// The bounds as the preset writes them: `< 0.98`.
    pub warn: Option<String>,
    pub fail: Option<String>,
    /// One line in plain words for a warn, fail or not-measured check.
    pub reason: Option<String>,
}

/// The verdict for one sample.
#[derive(Clone, Debug, PartialEq)]
pub struct VerdictSummary {
    pub sample: String,
    /// `pass`, `warn`, `fail` or `incomplete`.
    pub status: String,
    /// The preset as `check` prints it: `generic v3 (0123ab)`.
    pub preset: String,
    pub checks: Vec<CheckLine>,
}

impl VerdictSummary {
    /// How many checks ended in `outcome`.
    pub fn count(&self, outcome: Outcome) -> usize {
        self.checks.iter().filter(|c| c.outcome == outcome).count()
    }

    /// `1 fail, 0 warn, 12 pass, 2 not measured`.
    pub fn counts_line(&self) -> String {
        format!(
            "{} fail, {} warn, {} pass, {} not measured",
            self.count(Outcome::Fail),
            self.count(Outcome::Warn),
            self.count(Outcome::Pass),
            self.count(Outcome::NotMeasured)
        )
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    pub(crate) fn check(rule: &str, outcome: Outcome, value: Option<f64>) -> CheckLine {
        CheckLine {
            rule: rule.to_owned(),
            label: format!("label of {rule}"),
            metric: "DuplicationMetrics PERCENT_DUPLICATION".to_owned(),
            outcome,
            value,
            shown: value.map_or_else(|| "?".to_owned(), |v| v.to_string()),
            warn: Some("> 0.2".to_owned()),
            fail: Some("> 0.3".to_owned()),
            reason: (outcome != Outcome::Pass).then(|| format!("{rule} <is> \"off\" & more")),
        }
    }

    pub(crate) fn verdict() -> VerdictSummary {
        VerdictSummary {
            sample: "NA|12878".to_owned(),
            status: "fail".to_owned(),
            preset: "generic v3 (0123ab)".to_owned(),
            checks: vec![
                check("dup.rate", Outcome::Fail, Some(0.31)),
                check("q30", Outcome::Warn, Some(0.84)),
                check("insert.median", Outcome::Pass, Some(310.0)),
                check("wgs.depth", Outcome::NotMeasured, None),
            ],
        }
    }

    #[test]
    fn counts_follow_the_checks() {
        assert_eq!(
            verdict().counts_line(),
            "1 fail, 1 warn, 1 pass, 1 not measured"
        );
    }
}

//! Presets: a named, versioned list of [`Rule`]s read from a TOML file
//! (docs/verdict.md has the file format).
//!
//! A rule names a metrics class and a column, a row scope (an accumulation
//! level, and a `COLUMN=VALUE` selector such as `CATEGORY=PAIR`), one
//! comparator, and a warn and a fail threshold. Both thresholds are written
//! as bounds (`warn = "< 0.98"`, `fail = "< 0.95"`); the parser requires
//! the same comparator in both and a fail bound at least as strict as the
//! warn bound.

use std::fmt;
use std::fs;
use std::path::Path;

use bammetrics_core::levels::Level;
use bammetrics_core::metrics_file::format_double;

use crate::error::{Error, Result};
use crate::table::{Scalar, Table, parse_tables};

/// The shipped presets, by name, as `include_str!` text.
const BUILTIN: &[(&str, &str)] = &[
    ("generic", include_str!("../presets/generic.toml")),
    (
        "clinical-panel",
        include_str!("../presets/clinical-panel.toml"),
    ),
    ("wgs", include_str!("../presets/wgs.toml")),
    ("rna-seq", include_str!("../presets/rna-seq.toml")),
];

/// Names of the shipped presets, in the order they are listed.
pub fn names() -> Vec<&'static str> {
    BUILTIN.iter().map(|(name, _)| *name).collect()
}

/// The direction a value breaks a threshold in.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Comparator {
    /// `<`: breaks when the value is below the threshold.
    Lt,
    /// `<=`: breaks at or below.
    Le,
    /// `>`: breaks above.
    Gt,
    /// `>=`: breaks at or above.
    Ge,
    /// `outside A..B`: breaks below A or above B.
    Outside,
}

impl Comparator {
    pub fn symbol(self) -> &'static str {
        match self {
            Comparator::Lt => "<",
            Comparator::Le => "<=",
            Comparator::Gt => ">",
            Comparator::Ge => ">=",
            Comparator::Outside => "outside",
        }
    }
}

/// One threshold with its comparator, as written: `< 0.98`, `> 0.3`,
/// `outside 180..450`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Bound {
    pub comparator: Comparator,
    /// The threshold; the lower end for `outside`.
    pub low: f64,
    /// The upper end for `outside`; equal to `low` otherwise.
    pub high: f64,
}

impl Bound {
    /// `< 0.98`, `<= 1`, `> 0.3`, `>= 2`, `outside 180..450`.
    pub fn parse(text: &str) -> Option<Self> {
        let (op, rest) = text.trim().split_once(char::is_whitespace)?;
        let rest = rest.trim();
        let comparator = match op {
            "<" => Comparator::Lt,
            "<=" => Comparator::Le,
            ">" => Comparator::Gt,
            ">=" => Comparator::Ge,
            "outside" => Comparator::Outside,
            _ => return None,
        };
        let number = |s: &str| s.trim().parse::<f64>().ok().filter(|v| v.is_finite());
        let (low, high) = if comparator == Comparator::Outside {
            let (a, b) = rest.split_once("..")?;
            let (a, b) = (number(a)?, number(b)?);
            if a > b {
                return None;
            }
            (a, b)
        } else {
            let t = number(rest)?;
            (t, t)
        };
        Some(Bound {
            comparator,
            low,
            high,
        })
    }

    /// True when `value` breaks this bound.
    pub fn breaks(&self, value: f64) -> bool {
        match self.comparator {
            Comparator::Lt => value < self.low,
            Comparator::Le => value <= self.low,
            Comparator::Gt => value > self.low,
            Comparator::Ge => value >= self.low,
            Comparator::Outside => value < self.low || value > self.high,
        }
    }

    /// True when every value that breaks `other` also breaks `self`: a
    /// warn bound `self` encloses the fail bound `other` beside it.
    fn encloses(&self, other: &Bound) -> bool {
        match self.comparator {
            Comparator::Lt | Comparator::Le => other.low <= self.low,
            Comparator::Gt | Comparator::Ge => other.low >= self.low,
            Comparator::Outside => other.low <= self.low && other.high >= self.high,
        }
    }

    /// The words of a reason line: `below the fail limit of 0.95`,
    /// `outside the warn range 180 to 450`.
    pub fn describe(&self, status: &str) -> String {
        match self.comparator {
            Comparator::Lt => format!("below the {status} limit of {}", format_double(self.low)),
            Comparator::Le => format!(
                "at or below the {status} limit of {}",
                format_double(self.low)
            ),
            Comparator::Gt => format!("above the {status} limit of {}", format_double(self.low)),
            Comparator::Ge => format!(
                "at or above the {status} limit of {}",
                format_double(self.low)
            ),
            Comparator::Outside => format!(
                "outside the {status} range {} to {}",
                format_double(self.low),
                format_double(self.high)
            ),
        }
    }
}

impl fmt::Display for Bound {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.comparator {
            Comparator::Outside => write!(
                f,
                "outside {}..{}",
                format_double(self.low),
                format_double(self.high)
            ),
            op => write!(f, "{} {}", op.symbol(), format_double(self.low)),
        }
    }
}

/// One rule of a preset.
#[derive(Clone, Debug, PartialEq)]
pub struct Rule {
    /// Unique in the preset: `hs.target-20x`.
    pub id: String,
    /// Plain words for the reason line and the report: `Target bases at 20x`.
    pub label: String,
    /// The metrics class, by its short name (`HsMetrics`) or full Java name.
    pub class: String,
    /// The column read.
    pub metric: String,
    /// A second column the metric is divided by (`PF_Q30_BASES / PF_BASES`).
    pub divisor: Option<String>,
    /// The accumulation level of the rows read; `ALL_READS` by default.
    pub level: Level,
    /// `(COLUMN, VALUE)`: only rows whose `COLUMN` prints as `VALUE`.
    pub select: Option<(String, String)>,
    pub warn: Option<Bound>,
    pub fail: Option<Bound>,
    /// The plot family the report draws under a warn or fail line.
    pub plot: Option<String>,
}

impl Rule {
    /// The comparator both bounds share.
    pub fn comparator(&self) -> Comparator {
        self.fail
            .or(self.warn)
            .map_or(Comparator::Lt, |b| b.comparator)
    }

    /// `PCT_TARGET_BASES_20X`, or `PF_Q30_BASES / PF_BASES`.
    pub fn metric_text(&self) -> String {
        match &self.divisor {
            Some(d) => format!("{} / {d}", self.metric),
            None => self.metric.clone(),
        }
    }

    /// `CATEGORY=PAIR`.
    pub fn select_text(&self) -> Option<String> {
        self.select.as_ref().map(|(c, v)| format!("{c}={v}"))
    }
}

/// A parsed preset.
#[derive(Clone, Debug, PartialEq)]
pub struct Preset {
    pub name: String,
    pub version: i64,
    /// Free text from the file; the shipped presets say `starting values`.
    pub note: String,
    /// `built-in` or the path the file was read from.
    pub source: String,
    /// [`content_hash`] of the file text.
    pub hash: String,
    pub rules: Vec<Rule>,
}

impl Preset {
    /// A shipped preset by name.
    pub fn builtin(name: &str) -> Option<Result<Self>> {
        BUILTIN
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, text)| Self::parse(text, "built-in"))
    }

    /// The preset file at `path`.
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        Self::parse(&text, &path.display().to_string())
    }

    /// A shipped preset when `spec` is one of [`names`], else the file at
    /// that path.
    pub fn resolve(spec: &str) -> Result<Self> {
        if let Some(preset) = Self::builtin(spec) {
            return preset;
        }
        let path = Path::new(spec);
        if path.is_file() {
            return Self::load(path);
        }
        Err(Error::UnknownPreset {
            spec: spec.to_owned(),
        })
    }

    /// Parse preset text. `source` names it in errors.
    pub fn parse(text: &str, source: &str) -> Result<Self> {
        let tables = parse_tables(text, source)?;
        let mut preset = Preset {
            name: String::new(),
            version: 0,
            note: String::new(),
            source: source.to_owned(),
            hash: content_hash(text),
            rules: Vec::new(),
        };
        let mut seen_header = false;
        for table in &tables {
            let bad = |line: usize, message: String| Error::Preset {
                file: source.to_owned(),
                line: Some(line),
                message,
            };
            match (table.name.as_str(), table.array) {
                ("preset", false) => {
                    if seen_header {
                        return Err(bad(table.line, "a second [preset] table".to_owned()));
                    }
                    seen_header = true;
                    read_header(table, &mut preset, source)?;
                }
                ("rule", true) => {
                    let rule = read_rule(table, source)?;
                    if preset.rules.iter().any(|r| r.id == rule.id) {
                        return Err(bad(
                            table.line,
                            format!("rule id {:?} appears twice", rule.id),
                        ));
                    }
                    preset.rules.push(rule);
                }
                (name, array) => {
                    let shape = if array { "[[{name}]]" } else { "[{name}]" };
                    return Err(bad(
                        table.line,
                        format!(
                            "unknown table {}: expected [preset] or [[rule]]",
                            shape.replace("{name}", name)
                        ),
                    ));
                }
            }
        }
        if !seen_header {
            return Err(Error::Preset {
                file: source.to_owned(),
                line: None,
                message: "no [preset] table".to_owned(),
            });
        }
        if preset.rules.is_empty() {
            return Err(Error::Preset {
                file: source.to_owned(),
                line: None,
                message: "no [[rule]]".to_owned(),
            });
        }
        Ok(preset)
    }
}

/// FNV-1a 64 over the text with every `\r` removed, as 16 hex digits. The
/// same file hashes the same on every platform and checkout line ending.
pub fn content_hash(text: &str) -> String {
    bytes_hash(text.bytes().filter(|&b| b != b'\r'))
}

/// FNV-1a 64 over `bytes`, as 16 hex digits.
pub fn bytes_hash(bytes: impl IntoIterator<Item = u8>) -> String {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut hash = OFFSET;
    for byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(PRIME);
    }
    format!("{hash:016x}")
}

/// A table's entries with typed, line-numbered access and a check that
/// every key was read.
struct Entries<'a> {
    table: &'a Table,
    source: &'a str,
    used: Vec<bool>,
}

impl<'a> Entries<'a> {
    fn new(table: &'a Table, source: &'a str) -> Self {
        Entries {
            table,
            source,
            used: vec![false; table.entries.len()],
        }
    }

    fn error(&self, line: usize, message: String) -> Error {
        Error::Preset {
            file: self.source.to_owned(),
            line: Some(line),
            message,
        }
    }

    fn scalar(&mut self, key: &str) -> Option<(&'a Scalar, usize)> {
        let index = self.table.entries.iter().position(|e| e.key == key)?;
        self.used[index] = true;
        let entry = &self.table.entries[index];
        Some((&entry.value, entry.line))
    }

    fn string(&mut self, key: &str) -> Result<Option<(String, usize)>> {
        match self.scalar(key) {
            None => Ok(None),
            Some((Scalar::Str(s), line)) => Ok(Some((s.clone(), line))),
            Some((other, line)) => Err(self.error(
                line,
                format!("{key}: expected a string, got {}", other.kind()),
            )),
        }
    }

    fn required_string(&mut self, key: &str) -> Result<(String, usize)> {
        let found = self.string(key)?;
        found.ok_or_else(|| self.error(self.table.line, format!("missing key {key}")))
    }

    fn integer(&mut self, key: &str) -> Result<Option<i64>> {
        match self.scalar(key) {
            None => Ok(None),
            Some((Scalar::Int(v), _)) => Ok(Some(*v)),
            Some((other, line)) => Err(self.error(
                line,
                format!("{key}: expected an integer, got {}", other.kind()),
            )),
        }
    }

    /// An error for the first key no reader asked for.
    fn finish(self) -> Result<()> {
        match self
            .table
            .entries
            .iter()
            .zip(&self.used)
            .find(|(_, used)| !**used)
        {
            Some((entry, _)) => Err(self.error(entry.line, format!("unknown key {}", entry.key))),
            None => Ok(()),
        }
    }
}

fn read_header(table: &Table, preset: &mut Preset, source: &str) -> Result<()> {
    let mut entries = Entries::new(table, source);
    preset.name = entries.required_string("name")?.0;
    preset.version = entries
        .integer("version")?
        .ok_or_else(|| entries.error(table.line, "missing key version".to_owned()))?;
    preset.note = entries.string("note")?.map(|(s, _)| s).unwrap_or_default();
    entries.finish()
}

fn read_rule(table: &Table, source: &str) -> Result<Rule> {
    let mut entries = Entries::new(table, source);
    let (id, _) = entries.required_string("id")?;
    let label = entries
        .string("label")?
        .map_or_else(|| id.clone(), |(s, _)| s);
    let (class, _) = entries.required_string("class")?;
    let (metric_text, metric_line) = entries.required_string("metric")?;
    let (metric, divisor) = match metric_text.split_once('/') {
        Some((a, b)) => (a.trim().to_owned(), Some(b.trim().to_owned())),
        None => (metric_text.trim().to_owned(), None),
    };
    if metric.is_empty() || divisor.as_deref() == Some("") {
        return Err(entries.error(
            metric_line,
            format!("metric {metric_text:?}: expected COLUMN or COLUMN / COLUMN"),
        ));
    }
    let level = match entries.string("level")? {
        None => Level::AllReads,
        Some((text, line)) => Level::parse(&text).map_err(|_| {
            entries.error(
                line,
                format!("level {text:?}: expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP"),
            )
        })?,
    };
    let select = match entries.string("select")? {
        None => None,
        Some((text, line)) => Some(
            text.split_once('=')
                .map(|(c, v)| (c.trim().to_owned(), v.trim().to_owned()))
                .filter(|(c, _)| !c.is_empty())
                .ok_or_else(|| {
                    entries.error(line, format!("select {text:?}: expected COLUMN=VALUE"))
                })?,
        ),
    };
    let mut bound = |key: &str| -> Result<Option<Bound>> {
        match entries.string(key)? {
            None => Ok(None),
            Some((text, line)) => Bound::parse(&text).map(Some).ok_or_else(|| {
                entries.error(
                    line,
                    format!("{key} {text:?}: expected < N, <= N, > N, >= N or outside A..B"),
                )
            }),
        }
    };
    let warn = bound("warn")?;
    let fail = bound("fail")?;
    if warn.is_none() && fail.is_none() {
        return Err(entries.error(table.line, format!("rule {id}: needs warn, fail or both")));
    }
    if let (Some(w), Some(f)) = (warn, fail) {
        if w.comparator != f.comparator {
            return Err(entries.error(
                table.line,
                format!(
                    "rule {id}: warn uses {} and fail uses {}; one comparator per rule",
                    w.comparator.symbol(),
                    f.comparator.symbol()
                ),
            ));
        }
        if !w.encloses(&f) {
            return Err(entries.error(
                table.line,
                format!("rule {id}: fail bound {f} is looser than warn bound {w}"),
            ));
        }
    }
    let plot = entries.string("plot")?.map(|(s, _)| s);
    entries.finish()?;
    Ok(Rule {
        id,
        label,
        class,
        metric,
        divisor,
        level,
        select,
        warn,
        fail,
        plot,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINIMAL: &str = "[preset]\nname = \"t\"\nversion = 1\n\n[[rule]]\nid = \"a\"\nclass = \"HsMetrics\"\nmetric = \"MEAN_TARGET_COVERAGE\"\nwarn = \"< 200\"\nfail = \"< 100\"\n";

    #[test]
    fn bounds_parse_break_and_print() {
        struct Case {
            text: &'static str,
            comparator: Comparator,
            low: f64,
            high: f64,
            checks: &'static [(f64, bool)],
        }
        let cases = [
            Case {
                text: "< 0.98",
                comparator: Comparator::Lt,
                low: 0.98,
                high: 0.98,
                checks: &[(0.97, true), (0.98, false), (0.99, false)],
            },
            Case {
                text: "<= 1",
                comparator: Comparator::Le,
                low: 1.0,
                high: 1.0,
                checks: &[(1.0, true), (1.01, false)],
            },
            Case {
                text: "> 0.3",
                comparator: Comparator::Gt,
                low: 0.3,
                high: 0.3,
                checks: &[(0.3, false), (0.31, true)],
            },
            Case {
                text: ">=  2",
                comparator: Comparator::Ge,
                low: 2.0,
                high: 2.0,
                checks: &[(2.0, true), (1.99, false)],
            },
            Case {
                text: "outside 180..450",
                comparator: Comparator::Outside,
                low: 180.0,
                high: 450.0,
                checks: &[(179.0, true), (180.0, false), (450.0, false), (451.0, true)],
            },
        ];
        for case in cases {
            let bound = Bound::parse(case.text).expect(case.text);
            assert_eq!(
                (bound.comparator, bound.low, bound.high),
                (case.comparator, case.low, case.high),
                "{}",
                case.text
            );
            for &(value, breaks) in case.checks {
                assert_eq!(bound.breaks(value), breaks, "{} at {value}", case.text);
            }
        }
        assert_eq!(
            Bound::parse("outside 180..450")
                .expect("parses")
                .to_string(),
            "outside 180..450"
        );
        assert_eq!(Bound::parse(">=  2").expect("parses").to_string(), ">= 2");
        assert_eq!(
            Bound::parse("< 0.98").expect("parses").describe("fail"),
            "below the fail limit of 0.98"
        );
        assert_eq!(
            Bound::parse("<= 1").expect("parses").describe("warn"),
            "at or below the warn limit of 1"
        );
        assert_eq!(
            Bound::parse("> 5").expect("parses").describe("fail"),
            "above the fail limit of 5"
        );
        assert_eq!(
            Bound::parse(">= 5").expect("parses").describe("fail"),
            "at or above the fail limit of 5"
        );
        assert_eq!(
            Bound::parse("outside 1..2")
                .expect("parses")
                .describe("warn"),
            "outside the warn range 1 to 2"
        );
        for bad in [
            "",
            "0.98",
            "< ",
            "<0.98",
            "== 1",
            "< x",
            "outside 450..180",
            "outside 1",
            "< inf",
            "< 1 2",
        ] {
            assert_eq!(Bound::parse(bad), None, "{bad:?}");
        }
    }

    #[test]
    fn minimal_preset_parses() {
        let preset = Preset::parse(MINIMAL, "t.toml").expect("parses");
        assert_eq!(preset.name, "t");
        assert_eq!(preset.version, 1);
        assert_eq!(preset.note, "");
        assert_eq!(preset.source, "t.toml");
        assert_eq!(preset.hash.len(), 16);
        let [rule] = preset.rules.as_slice() else {
            panic!("one rule");
        };
        assert_eq!(rule.label, "a");
        assert_eq!(rule.level, Level::AllReads);
        assert_eq!(rule.select, None);
        assert_eq!(rule.divisor, None);
        assert_eq!(rule.comparator(), Comparator::Lt);
        assert_eq!(rule.metric_text(), "MEAN_TARGET_COVERAGE");
    }

    #[test]
    fn every_rule_key_parses() {
        let text = "[preset]\nname = \"t\"\nversion = 2\nnote = \"starting values\"\n[[rule]]\nid = \"yield.q30\"\nlabel = \"Q30 fraction\"\nclass = \"QualityYieldMetrics\"\nmetric = \"PF_Q30_BASES / PF_BASES\"\nlevel = \"SAMPLE\"\nselect = \"CATEGORY=PAIR\"\nwarn = \"< 0.85\"\nplot = \"cycle-quality\"\n";
        let preset = Preset::parse(text, "t").expect("parses");
        assert_eq!(preset.note, "starting values");
        let rule = &preset.rules[0];
        assert_eq!(rule.label, "Q30 fraction");
        assert_eq!(rule.metric, "PF_Q30_BASES");
        assert_eq!(rule.divisor.as_deref(), Some("PF_BASES"));
        assert_eq!(rule.metric_text(), "PF_Q30_BASES / PF_BASES");
        assert_eq!(rule.level, Level::Sample);
        assert_eq!(rule.select_text().as_deref(), Some("CATEGORY=PAIR"));
        assert_eq!(rule.fail, None);
        assert_eq!(rule.plot.as_deref(), Some("cycle-quality"));
    }

    /// A preset with the default header and `body` as its one rule's lines.
    fn with_rule(body: &str) -> String {
        format!("[preset]\nname = \"t\"\nversion = 1\n[[rule]]\n{body}")
    }

    #[test]
    fn preset_errors_name_the_line() {
        let rule = "[[rule]]\nid = \"a\"\nclass = \"X\"\nmetric = \"M\"\nfail = \"< 1\"\n";
        let cases: &[(String, &str)] = &[
            (rule.to_owned(), "no [preset] table"),
            (
                "[preset]\nname = \"t\"\nversion = 1\n".to_owned(),
                "no [[rule]]",
            ),
            (format!("[preset]\nversion = 1\n{rule}"), "missing key name"),
            (
                format!("[preset]\nname = \"t\"\n{rule}"),
                "missing key version",
            ),
            (
                format!("[preset]\nname = \"t\"\nversion = \"1\"\n{rule}"),
                "version: expected an integer",
            ),
            (
                format!("[preset]\nname = \"t\"\nversion = 1\nextra = 1\n{rule}"),
                "unknown key extra",
            ),
            (
                format!(
                    "[preset]\nname = \"t\"\nversion = 1\n[preset]\nname = \"u\"\nversion = 1\n{rule}"
                ),
                "second [preset]",
            ),
            (
                format!("[preset]\nname = \"t\"\nversion = 1\n[other]\n{rule}"),
                "unknown table [other]",
            ),
            (
                format!("[preset]\nname = \"t\"\nversion = 1\n[[rules]]\n{rule}"),
                "unknown table [[rules]]",
            ),
            (
                format!("[preset]\nname = \"t\"\nversion = 1\n{rule}{rule}"),
                "appears twice",
            ),
            (
                with_rule("class = \"X\"\nmetric = \"M\"\nfail = \"< 1\"\n"),
                "missing key id",
            ),
            (
                with_rule("id = \"a\"\nclass = \"X\"\nmetric = \"M\"\n"),
                "needs warn, fail or both",
            ),
            (
                with_rule("id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nfail = \"below 1\"\n"),
                "fail \"below 1\": expected",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nwarn = \"< 1\"\nfail = \"> 1\"\n",
                ),
                "one comparator per rule",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nwarn = \"< 1\"\nfail = \"< 2\"\n",
                ),
                "fail bound < 2 is looser than warn bound < 1",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nwarn = \"> 2\"\nfail = \"> 1\"\n",
                ),
                "looser",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nwarn = \"outside 1..4\"\nfail = \"outside 2..3\"\n",
                ),
                "looser",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nlevel = \"ALL\"\nfail = \"< 1\"\n",
                ),
                "level \"ALL\"",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nselect = \"PAIR\"\nfail = \"< 1\"\n",
                ),
                "select \"PAIR\"",
            ),
            (
                with_rule("id = \"a\"\nclass = \"X\"\nmetric = \"A /\"\nfail = \"< 1\"\n"),
                "metric \"A /\"",
            ),
            (
                with_rule(
                    "id = \"a\"\nclass = \"X\"\nmetric = \"M\"\nfail = \"< 1\"\nwhy = \"x\"\n",
                ),
                "unknown key why",
            ),
            (
                with_rule("id = 1\nclass = \"X\"\nmetric = \"M\"\nfail = \"< 1\"\n"),
                "id: expected a string, got an integer",
            ),
        ];
        for (text, needle) in cases {
            let error = Preset::parse(text, "p.toml").expect_err(text).to_string();
            assert!(error.contains(needle), "{text:?}\n  got: {error}");
            assert!(error.starts_with("p.toml:"), "{error}");
        }
    }

    #[test]
    fn hash_ignores_carriage_returns_and_is_pinned() {
        assert_eq!(content_hash(""), "cbf29ce484222325");
        assert_eq!(content_hash("a"), "af63dc4c8601ec8c");
        assert_eq!(
            content_hash(MINIMAL),
            content_hash(&MINIMAL.replace('\n', "\r\n"))
        );
        assert_ne!(
            content_hash(MINIMAL),
            content_hash(&MINIMAL.replace("200", "201"))
        );
    }

    #[test]
    fn resolve_finds_builtins_and_files() {
        assert_eq!(names(), ["generic", "clinical-panel", "wgs", "rna-seq"]);
        let generic = Preset::resolve("generic").expect("generic");
        assert_eq!(generic.source, "built-in");
        let error = Preset::resolve("no-such-preset")
            .expect_err("unknown")
            .to_string();
        assert!(
            error.contains("no preset \"no-such-preset\"") && error.contains("clinical-panel"),
            "{error}"
        );
        let dir = std::env::temp_dir().join(format!("bammetrics-verdict-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("lab.toml");
        std::fs::write(&path, MINIMAL).expect("write");
        let loaded = Preset::resolve(&path.display().to_string()).expect("file preset");
        assert_eq!(loaded.source, path.display().to_string());
        assert_eq!(loaded.hash, content_hash(MINIMAL));
        std::fs::remove_dir_all(&dir).expect("cleanup");
    }
}

//! `ConvertSequencingArtifactToOxoG`, mirroring Picard 3.5.0: the OxoG
//! rows of `CollectOxoGMetrics` (`CpcgMetrics`) from the pre-adapter and
//! bait-bias detail files of `CollectSequencingArtifactMetrics`.
//!
//! The detail files are `INPUT_BASE.pre_adapter_detail_metrics` and
//! `INPUT_BASE.bait_bias_detail_metrics` unless `PRE_ADAPTER_IN` and
//! `BAIT_BIAS_IN` name them. The output is `OXOG_OUT`, else `OUTPUT_BASE`
//! (default `INPUT_BASE`) plus `.oxog_metrics`.
//!
//! Rows: one per library and per context with a `C` at its centre, both
//! in Java `HashSet` order of the values as they first appear in the
//! pre-adapter file. `SAMPLE_ALIAS` is the first pre-adapter row's. The
//! oxidation columns come from the `G>T` pre-adapter row of the reverse
//! complement context: `TOTAL_BASES` is its four counts, `REF_OXO_BASES`
//! and `ALT_OXO_BASES` its pro counts, the `NONOXO` columns its con counts,
//! `OXIDATION_ERROR_RATE` is `max(ALT_OXO - ALT_NONOXO, 1) / TOTAL_BASES`.
//! The reference-bias columns come from the `C>A` bait-bias row of the
//! context: `C_REF` its forward counts, `G_REF` its reverse counts, and each
//! error rate the difference of the two alt fractions, at least 1e-10, NaN
//! (`?`) when a fraction has no bases. `TOTAL_SITES` is 0. A missing row, a
//! library in one file only, an empty pre-adapter file or a base outside
//! `ACGT` fails the run, as Picard fails.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::java::{self, hash_map_order};
use crate::metric::{MetricRow, MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};
use crate::sequencing_artifact::{BAIT_BIAS_DETAILS_EXT, PRE_ADAPTER_DETAILS_EXT};

/// The extension after `OUTPUT_BASE`.
pub const OXOG_METRICS_EXT: &str = ".oxog_metrics";

const PRE_ADAPTER_CLASS: &str =
    "picard.analysis.artifacts.SequencingArtifactMetrics$PreAdapterDetailMetrics";
const BAIT_BIAS_CLASS: &str =
    "picard.analysis.artifacts.SequencingArtifactMetrics$BaitBiasDetailMetrics";

/// `ConvertSequencingArtifactToOxoG` in the tool registry.
pub const TOOL: Tool = Tool {
    name: "ConvertSequencingArtifactToOxoG",
    summary: "one CpcgMetrics row per library and C-centred context, read from the \
              CollectSequencingArtifactMetrics detail files of INPUT_BASE; the file is \
              OUTPUT_BASE.oxog_metrics. Python returns one row per file: file and rows.",
    options: &[
        Opt::output(&["OUTPUT_BASE", "O"]).with_help(
            "The output prefix: the file is OUTPUT_BASE.oxog_metrics (default INPUT_BASE).",
        ),
        Opt::path(
            &["INPUT_BASE", "I"],
            "The prefix of a CollectSequencingArtifactMetrics run.",
        ),
        Opt::path(
            &["PRE_ADAPTER_IN"],
            "The pre-adapter detail file (default INPUT_BASE.pre_adapter_detail_metrics).",
        ),
        Opt::path(
            &["BAIT_BIAS_IN"],
            "The bait-bias detail file (default INPUT_BASE.bait_bias_detail_metrics).",
        ),
        Opt::path(
            &["OXOG_OUT"],
            "The output file (default OUTPUT_BASE.oxog_metrics).",
        )
        .written(),
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::rejected(&["INPUT"]),
        Opt::rejected(&["OUTPUT"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

/// `CpcgMetrics` of `CollectOxoGMetrics`.
#[derive(Clone, Debug, PartialEq)]
pub struct Row {
    pub sample_alias: String,
    pub library: String,
    pub context: String,
    pub total_sites: i64,
    pub total_bases: i64,
    pub ref_nonoxo_bases: i64,
    pub ref_oxo_bases: i64,
    pub ref_total_bases: i64,
    pub alt_nonoxo_bases: i64,
    pub alt_oxo_bases: i64,
    pub oxidation_error_rate: f64,
    pub oxidation_q: f64,
    pub c_ref_ref_bases: i64,
    pub g_ref_ref_bases: i64,
    pub c_ref_alt_bases: i64,
    pub g_ref_alt_bases: i64,
    pub c_ref_oxo_error_rate: f64,
    pub c_ref_oxo_q: f64,
    pub g_ref_oxo_error_rate: f64,
    pub g_ref_oxo_q: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.CollectOxoGMetrics$CpcgMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "SAMPLE_ALIAS",
        "LIBRARY",
        "CONTEXT",
        "TOTAL_SITES",
        "TOTAL_BASES",
        "REF_NONOXO_BASES",
        "REF_OXO_BASES",
        "REF_TOTAL_BASES",
        "ALT_NONOXO_BASES",
        "ALT_OXO_BASES",
        "OXIDATION_ERROR_RATE",
        "OXIDATION_Q",
        "C_REF_REF_BASES",
        "G_REF_REF_BASES",
        "C_REF_ALT_BASES",
        "G_REF_ALT_BASES",
        "C_REF_OXO_ERROR_RATE",
        "C_REF_OXO_Q",
        "G_REF_OXO_ERROR_RATE",
        "G_REF_OXO_Q",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.sample_alias.as_str().into(),
            self.library.as_str().into(),
            self.context.as_str().into(),
            self.total_sites.into(),
            self.total_bases.into(),
            self.ref_nonoxo_bases.into(),
            self.ref_oxo_bases.into(),
            self.ref_total_bases.into(),
            self.alt_nonoxo_bases.into(),
            self.alt_oxo_bases.into(),
            self.oxidation_error_rate.into(),
            self.oxidation_q.into(),
            self.c_ref_ref_bases.into(),
            self.g_ref_ref_bases.into(),
            self.c_ref_alt_bases.into(),
            self.g_ref_alt_bases.into(),
            self.c_ref_oxo_error_rate.into(),
            self.c_ref_oxo_q.into(),
            self.g_ref_oxo_error_rate.into(),
            self.g_ref_oxo_q.into(),
        ]
    }
}

/// The file paths after `customCommandLineValidation`.
struct Paths {
    pre_adapter: PathBuf,
    bait_bias: PathBuf,
    output: PathBuf,
}

/// `base` plus `extension`, as Java's `new File(base + extension)`.
fn with_extension(base: &str, extension: &str) -> PathBuf {
    PathBuf::from(format!("{base}{extension}"))
}

impl Paths {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let input_base = inv.get("INPUT_BASE");
        let output_base = inv
            .output
            .as_deref()
            .map(|p| p.to_string_lossy().into_owned())
            .or_else(|| input_base.map(str::to_owned));
        let mut errors = Vec::new();
        let mut path = |given: Option<&str>, base: Option<&str>, extension: &str, error: &str| {
            given.map(PathBuf::from).or_else(|| {
                if base.is_none() {
                    errors.push(error.to_owned());
                }
                base.map(|b| with_extension(b, extension))
            })
        };
        let pre_adapter = path(
            inv.get("PRE_ADAPTER_IN"),
            input_base,
            PRE_ADAPTER_DETAILS_EXT,
            "Must specify either INPUT_BASE or PRE_ADAPTER_IN",
        );
        let bait_bias = path(
            inv.get("BAIT_BIAS_IN"),
            input_base,
            BAIT_BIAS_DETAILS_EXT,
            "Must specify either INPUT_BASE or BAIT_BIAS_IN",
        );
        let output = path(
            inv.get("OXOG_OUT"),
            output_base.as_deref(),
            OXOG_METRICS_EXT,
            "Must specify either OUTPUT_BASE or OXOG_OUT",
        );
        match (pre_adapter, bait_bias, output) {
            (Some(pre_adapter), Some(bait_bias), Some(output)) => Ok(Self {
                pre_adapter,
                bait_bias,
                output,
            }),
            _ => Err(Error::Usage(errors.join("\n"))),
        }
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let paths = Paths::from_invocation(inv)?;
    let pre_adapter = read_beans(&paths.pre_adapter, &PRE_ADAPTER)?;
    let bait_bias = read_beans(&paths.bait_bias, &BAIT_BIAS)?;
    // The two detail files may come from `INPUT_BASE`, which names neither.
    crate::registry::check_writable(
        [paths.output.as_path()],
        &[paths.pre_adapter.clone(), paths.bait_bias.clone()],
    )?;
    let rows = convert(&pre_adapter, &bait_bias)?;
    Ok(ToolOutput {
        main: None,
        extra: vec![(paths.output, MetricsOutput::from_rows(&rows, None))],
    })
}

/// The type of one field of a metrics class, as `FormatUtil.parseObject`
/// reads it.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Field {
    Text,
    Char,
    Long,
    Double,
}

/// A Picard metrics class: its name and its public fields.
struct Class {
    name: &'static str,
    fields: &'static [(&'static str, Field)],
}

const PRE_ADAPTER: Class = Class {
    name: PRE_ADAPTER_CLASS,
    fields: &[
        ("SAMPLE_ALIAS", Field::Text),
        ("LIBRARY", Field::Text),
        ("REF_BASE", Field::Char),
        ("ALT_BASE", Field::Char),
        ("CONTEXT", Field::Text),
        ("PRO_REF_BASES", Field::Long),
        ("PRO_ALT_BASES", Field::Long),
        ("CON_REF_BASES", Field::Long),
        ("CON_ALT_BASES", Field::Long),
        ("ERROR_RATE", Field::Double),
        ("QSCORE", Field::Double),
    ],
};

const BAIT_BIAS: Class = Class {
    name: BAIT_BIAS_CLASS,
    fields: &[
        ("SAMPLE_ALIAS", Field::Text),
        ("LIBRARY", Field::Text),
        ("REF_BASE", Field::Char),
        ("ALT_BASE", Field::Char),
        ("CONTEXT", Field::Text),
        ("FWD_CXT_REF_BASES", Field::Long),
        ("FWD_CXT_ALT_BASES", Field::Long),
        ("REV_CXT_REF_BASES", Field::Long),
        ("REV_CXT_ALT_BASES", Field::Long),
        ("FWD_ERROR_RATE", Field::Double),
        ("REV_ERROR_RATE", Field::Double),
        ("ERROR_RATE", Field::Double),
        ("QSCORE", Field::Double),
    ],
};

/// One detail row as the converter reads it: the text fields (null when
/// the column is absent or the cell empty) and the counts (0 when the
/// column is absent).
#[derive(Clone, Debug, Default)]
struct Bean {
    sample_alias: Option<String>,
    library: Option<String>,
    ref_base: u8,
    alt_base: u8,
    context: Option<String>,
    counts: [i64; 4],
}

/// `MetricsFile.readBeans`: the rows of the file's metrics class, which
/// must be `class` (Picard's cast fails otherwise once a row is read).
/// Every cell is parsed by its field's type; an empty cell of a number or
/// character field, or a value that does not parse, fails as htsjdk does.
fn read_beans(path: &Path, class: &Class) -> Result<Vec<Bean>> {
    let text = fs::read_to_string(path).map_err(|source| Error::Input {
        path: path.to_path_buf(),
        source,
    })?;
    let fail = |message: String| Error::Metric(format!("{}: {message}", path.display()));
    let mut lines = text.lines();
    let mut found = None;
    for line in lines.by_ref() {
        let line = java::trim(line);
        if line.is_empty() || line.starts_with("# ") {
            continue;
        }
        if line.starts_with("## METRICS CLASS") {
            found = Some(line.split('\t').nth(1).unwrap_or_default().to_owned());
            break;
        }
        if line.starts_with("## HISTOGRAM") {
            break;
        }
        if !line.starts_with("## ") {
            return Err(fail(format!(
                "Illegal state. Found following string in metrics file header: {line}"
            )));
        }
    }
    let Some(found) = found else {
        return Ok(Vec::new());
    };
    let columns: Vec<&str> = lines.next().unwrap_or_default().split('\t').collect();
    let mut types = Vec::with_capacity(columns.len());
    for column in &columns {
        let field = class
            .fields
            .iter()
            .find(|(name, _)| name == column)
            .map(|(_, t)| *t);
        match field {
            Some(t) => types.push(t),
            None if found == class.name => {
                return Err(fail(format!(
                    "Could not get field with name {column} from class {found}"
                )));
            }
            None => types.push(Field::Text),
        }
    }
    let mut beans = Vec::new();
    for line in lines {
        if java::trim(line).is_empty() {
            break;
        }
        let values: Vec<&str> = line.split('\t').collect();
        if values.len() < columns.len() {
            return Err(fail(format!(
                "a row has {} fields for {} columns",
                values.len(),
                columns.len()
            )));
        }
        let mut bean = Bean::default();
        for ((column, field), value) in columns.iter().zip(&types).zip(&values) {
            let empty = value.is_empty();
            match field {
                Field::Text => {
                    let text = (!empty).then(|| (*value).to_owned());
                    match *column {
                        "SAMPLE_ALIAS" => bean.sample_alias = text,
                        "LIBRARY" => bean.library = text,
                        "CONTEXT" => bean.context = text,
                        _ => {}
                    }
                }
                _ if empty => {
                    return Err(fail(format!("{column} is empty in a row")));
                }
                Field::Char => {
                    let mut chars = value.chars();
                    let (Some(c), None) = (chars.next(), chars.next()) else {
                        return Err(fail(format!(
                            "Cannot parse string into char because length != 1 : {value}"
                        )));
                    };
                    let byte = u8::try_from(c).unwrap_or(0);
                    match *column {
                        "REF_BASE" => bean.ref_base = byte,
                        _ => bean.alt_base = byte,
                    }
                }
                Field::Long => {
                    let n = java::parse_integer(value, 64).map_err(fail)?;
                    let slot = class.fields[5..9]
                        .iter()
                        .position(|(name, _)| name == column);
                    if let Some(slot) = slot {
                        bean.counts[slot] = n;
                    }
                }
                Field::Double => {
                    if *value != "?" && *value != "-?" {
                        java::parse_double(value).map_err(fail)?;
                    }
                }
            }
        }
        beans.push(bean);
    }
    if found != class.name && !beans.is_empty() {
        return Err(fail(format!(
            "class {found} cannot be cast to class {}",
            class.name
        )));
    }
    Ok(beans)
}

/// `Transition.transitionOf(ref, alt)` is `C>A` or `G>T`; a base outside
/// `ACGT` is Picard's `IllegalArgumentException`.
fn is_oxog(bean: &Bean) -> Result<bool> {
    let acgt = |b: u8| matches!(b, b'A' | b'C' | b'G' | b'T');
    if !acgt(bean.ref_base) || !acgt(bean.alt_base) {
        return Err(Error::Metric(format!(
            "Base params should be one of {{A, C, T, G}} but ref={} and call={}",
            char::from(bean.ref_base),
            char::from(bean.alt_base)
        )));
    }
    Ok(matches!(
        (bean.ref_base, bean.alt_base),
        (b'C', b'A') | (b'G', b'T')
    ))
}

/// The key of a text field: Java's `HashMap` takes a null key too.
fn key(text: &Option<String>) -> String {
    text.clone().unwrap_or_else(|| "null".to_owned())
}

/// `SequenceUtil.reverseComplement`.
fn reverse_complement(context: &str) -> String {
    context
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            'a' => 't',
            'c' => 'g',
            'g' => 'c',
            't' => 'a',
            other => other,
        })
        .collect()
}

/// `Math.max` of two doubles: NaN when either is NaN.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `doWork` from the two lists to the rows.
fn convert(pre_adapter: &[Bean], bait_bias: &[Bean]) -> Result<Vec<Row>> {
    let missing =
        |what: String| Error::Metric(format!("{what} (Picard fails with a NullPointerException)"));
    let first = pre_adapter.first().ok_or_else(|| {
        Error::Metric("the pre-adapter detail file has no rows (Picard fails with an IndexOutOfBoundsException)".to_owned())
    })?;
    let sample_alias = first.sample_alias.clone().unwrap_or_default();
    let mut libraries: Vec<String> = Vec::new();
    let mut contexts: Vec<String> = Vec::new();
    for bean in pre_adapter {
        let library = key(&bean.library);
        if !libraries.contains(&library) {
            libraries.push(library);
        }
        if bean.ref_base == b'C' {
            let context = key(&bean.context);
            if !contexts.contains(&context) {
                contexts.push(context);
            }
        }
    }
    let mut pre_by_library: HashMap<String, HashMap<String, &Bean>> = libraries
        .iter()
        .map(|l| (l.clone(), HashMap::new()))
        .collect();
    let mut bait_by_library: HashMap<String, HashMap<String, &Bean>> = libraries
        .iter()
        .map(|l| (l.clone(), HashMap::new()))
        .collect();
    for bean in pre_adapter {
        if is_oxog(bean)? {
            pre_by_library
                .get_mut(&key(&bean.library))
                .expect("every pre-adapter library has a map")
                .insert(key(&bean.context), bean);
        }
    }
    for bean in bait_bias {
        if is_oxog(bean)? {
            let library = key(&bean.library);
            bait_by_library
                .get_mut(&library)
                .ok_or_else(|| missing(format!("library {library} is in the bait-bias file only")))?
                .insert(key(&bean.context), bean);
        }
    }
    let library_keys: Vec<&str> = libraries.iter().map(String::as_str).collect();
    let context_keys: Vec<&str> = contexts.iter().map(String::as_str).collect();
    let context_order = hash_map_order(&context_keys);
    let mut rows = Vec::new();
    for l in hash_map_order(&library_keys) {
        let library = &libraries[l];
        for &c in &context_order {
            let context = &contexts[c];
            let pre = pre_by_library[library]
                .get(&reverse_complement(context))
                .ok_or_else(|| {
                    missing(format!(
                        "no G>T pre-adapter row for library {library}, context {}",
                        reverse_complement(context)
                    ))
                })?;
            let bait = bait_by_library[library].get(context).ok_or_else(|| {
                missing(format!(
                    "no C>A bait-bias row for library {library}, context {context}"
                ))
            })?;
            rows.push(row(&sample_alias, library, context, pre, bait));
        }
    }
    Ok(rows)
}

/// One `CpcgMetrics` row from its two detail rows.
fn row(sample_alias: &str, library: &str, context: &str, pre: &Bean, bait: &Bean) -> Row {
    let [pro_ref, pro_alt, con_ref, con_alt] = pre.counts;
    let [c_ref_ref, c_ref_alt, g_ref_ref, g_ref_alt] = bait.counts;
    let total_bases = pro_ref + pro_alt + con_ref + con_alt;
    let oxidation_error_rate = (pro_alt - con_alt).max(1) as f64 / total_bases as f64;
    let c_ref_rate = c_ref_alt as f64 / (c_ref_alt + c_ref_ref) as f64;
    let g_ref_rate = g_ref_alt as f64 / (g_ref_alt + g_ref_ref) as f64;
    let c_ref_oxo_error_rate = java_max(c_ref_rate - g_ref_rate, 1e-10);
    let g_ref_oxo_error_rate = java_max(g_ref_rate - c_ref_rate, 1e-10);
    Row {
        sample_alias: sample_alias.to_owned(),
        library: library.to_owned(),
        context: context.to_owned(),
        total_sites: 0,
        total_bases,
        ref_nonoxo_bases: con_ref,
        ref_oxo_bases: pro_ref,
        ref_total_bases: pro_ref + con_ref,
        alt_nonoxo_bases: con_alt,
        alt_oxo_bases: pro_alt,
        oxidation_error_rate,
        oxidation_q: -10.0 * oxidation_error_rate.log10(),
        c_ref_ref_bases: c_ref_ref,
        g_ref_ref_bases: g_ref_ref,
        c_ref_alt_bases: c_ref_alt,
        g_ref_alt_bases: g_ref_alt,
        c_ref_oxo_error_rate,
        c_ref_oxo_q: -10.0 * c_ref_oxo_error_rate.log10(),
        g_ref_oxo_error_rate,
        g_ref_oxo_q: -10.0 * g_ref_oxo_error_rate.log10(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bean(library: &str, r: u8, a: u8, context: &str, counts: [i64; 4]) -> Bean {
        Bean {
            sample_alias: Some("S".to_owned()),
            library: Some(library.to_owned()),
            ref_base: r,
            alt_base: a,
            context: Some(context.to_owned()),
            counts,
        }
    }

    #[test]
    fn a_row_takes_the_reverse_complement_pre_adapter_counts() {
        let pre = [
            bean("L", b'C', b'A', "ACA", [1, 1, 1, 1]),
            bean("L", b'G', b'T', "TGT", [100, 7, 90, 2]),
        ];
        let bait = [bean("L", b'C', b'A', "ACA", [50, 5, 40, 0])];
        let rows = convert(&pre, &bait).expect("converts");
        assert_eq!(rows.len(), 1);
        let r = &rows[0];
        assert_eq!((r.context.as_str(), r.total_bases), ("ACA", 199));
        assert_eq!(
            (r.ref_oxo_bases, r.alt_oxo_bases, r.ref_nonoxo_bases),
            (100, 7, 90)
        );
        assert_eq!(r.oxidation_error_rate, 5.0 / 199.0);
        assert_eq!((r.c_ref_ref_bases, r.c_ref_alt_bases), (50, 5));
        assert_eq!((r.g_ref_ref_bases, r.g_ref_alt_bases), (40, 0));
        assert_eq!(r.g_ref_oxo_error_rate, 1e-10);
        assert_eq!(r.g_ref_oxo_q, 100.0);
    }

    #[test]
    fn empty_counts_give_nan_and_infinity() {
        let pre = [
            bean("L", b'C', b'A', "C", [0; 4]),
            bean("L", b'G', b'T', "G", [0; 4]),
        ];
        let bait = [bean("L", b'C', b'A', "C", [0; 4])];
        let r = &convert(&pre, &bait).expect("converts")[0];
        assert!(r.oxidation_error_rate.is_infinite());
        assert!(r.c_ref_oxo_error_rate.is_nan() && r.c_ref_oxo_q.is_nan());
    }

    #[test]
    fn missing_rows_fail() {
        assert!(convert(&[], &[]).is_err());
        let pre = [bean("L", b'C', b'A', "ACA", [1; 4])];
        assert!(convert(&pre, &[]).is_err(), "no TGT row");
        let pre = [bean("L", b'N', b'A', "ANA", [1; 4])];
        assert!(convert(&pre, &[]).is_err(), "N is no base");
    }
}

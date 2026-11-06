//! The Picard 3.5.0 metrics classes as htsjdk `MetricsFile.read` sees them:
//! `picard_fields.tsv`, written by `tests/data/PicardFields.java` (`make
//! picard-fields`) from `picard.jar` by reflection. Each class lists its
//! public fields in `Class.getFields()` order with the Java type and the
//! value a new instance holds. `CompareMetrics` reads a file with these
//! types, so `0.5` and `0.50` in a `double` column are one value and `5.0`
//! in a `long` column is an error, as in Picard.
//!
//! [`load_class`] is htsjdk's `MetricsFile.loadClass`: the name as given,
//! then the class's simple name (and, for a nested class, its inner name)
//! under each of 17 legacy and current packages, so a file naming
//! `net.sf.picard.analysis.InsertSizeMetrics` reads as
//! `picard.analysis.InsertSizeMetrics`.

use std::collections::HashMap;
use std::sync::OnceLock;

/// The generated table.
pub const TABLE: &str = include_str!("../picard_fields.tsv");

/// The packages `MetricsFile.loadClass` tries, in its order.
const PACKAGES: &[&str] = &[
    "edu.mit.broad.picard.genotype.concordance",
    "edu.mit.broad.picard.genotype.fingerprint",
    "edu.mit.broad.picard.ic",
    "edu.mit.broad.picard.illumina",
    "edu.mit.broad.picard.jumping",
    "edu.mit.broad.picard.quality",
    "edu.mit.broad.picard.samplevalidation",
    "htsjdk.samtools.analysis",
    "htsjdk.samtools.analysis.directed",
    "htsjdk.samtools.sam",
    "htsjdk.samtools.metrics",
    "picard.sam",
    "picard.metrics",
    "picard.illumina",
    "picard.analysis",
    "picard.analysis.directed",
    "picard.vcf",
];

/// The `java.lang` classes a `## HISTOGRAM` line names in practice.
pub const KEY_CLASSES: &[&str] = &[
    "java.lang.Short",
    "java.lang.Integer",
    "java.lang.Long",
    "java.lang.Float",
    "java.lang.Double",
    "java.lang.Boolean",
    "java.lang.Byte",
    "java.lang.Character",
    "java.lang.String",
];

/// The Java type of a field, as `FormatUtil.parseObject` reads it.
#[derive(Clone, Debug, PartialEq)]
pub enum FieldType {
    /// `long`, `int`, `short`, `byte` or their boxes: the width in bits.
    Integer(u32),
    Double,
    Float,
    Boolean,
    Char,
    /// `String`, and `Iso8601Date` (compared as its text here).
    Text,
    /// An enum: its canonical name (for the error) and its constants.
    Enum {
        class: &'static str,
        constants: Vec<&'static str>,
    },
    /// A type `FormatUtil` cannot parse (`Object`, a nested metrics class).
    Other(&'static str),
}

/// One public field.
#[derive(Clone, Debug)]
pub struct Field {
    pub name: &'static str,
    pub ty: FieldType,
    /// A primitive field cannot take an empty cell (htsjdk sets null).
    pub primitive: bool,
    /// `String.valueOf` of the field on a new instance, `None` for null.
    pub default: Option<String>,
}

/// One metrics class.
#[derive(Clone, Debug)]
pub struct MetricClass {
    pub name: &'static str,
    /// False for an abstract class: a row of it cannot be read.
    pub instantiable: bool,
    pub fields: Vec<Field>,
}

impl MetricClass {
    /// Index of the field `name` (`Class.getField`).
    pub fn field(&self, name: &str) -> Option<usize> {
        self.fields.iter().position(|f| f.name == name)
    }
}

struct Classes {
    metrics: HashMap<&'static str, MetricClass>,
    headers: Vec<&'static str>,
}

fn classes() -> &'static Classes {
    static CLASSES: OnceLock<Classes> = OnceLock::new();
    CLASSES.get_or_init(|| parse(TABLE))
}

fn parse(table: &'static str) -> Classes {
    let mut metrics: HashMap<&'static str, MetricClass> = HashMap::new();
    let mut headers = Vec::new();
    let mut current: Option<MetricClass> = None;
    for line in table.lines().filter(|l| !l.starts_with('#')) {
        let cells: Vec<&'static str> = line.split('\t').collect();
        match cells.as_slice() {
            ["M", name, instantiable] => {
                if let Some(done) = current.take() {
                    metrics.insert(done.name, done);
                }
                current = Some(MetricClass {
                    name,
                    instantiable: *instantiable == "Y",
                    fields: Vec::new(),
                });
            }
            ["F", name, ty, _static, default] => {
                if let Some(class) = current.as_mut() {
                    let (ty, primitive) = field_type(ty);
                    let default = (*default != "\\N").then(|| unescape(default));
                    class.fields.push(Field {
                        name,
                        ty,
                        primitive,
                        default,
                    });
                }
            }
            ["H", name] => headers.push(*name),
            _ => {}
        }
    }
    if let Some(done) = current.take() {
        metrics.insert(done.name, done);
    }
    Classes { metrics, headers }
}

/// A table type name as a [`FieldType`] and whether it is primitive.
fn field_type(name: &'static str) -> (FieldType, bool) {
    if let Some(rest) = name.strip_prefix("enum:") {
        let (class, constants) = rest.split_once(':').unwrap_or((rest, ""));
        let constants = constants.split(',').filter(|c| !c.is_empty()).collect();
        return (FieldType::Enum { class, constants }, false);
    }
    match name {
        "long" => (FieldType::Integer(64), true),
        "int" => (FieldType::Integer(32), true),
        "short" => (FieldType::Integer(16), true),
        "byte" => (FieldType::Integer(8), true),
        "double" => (FieldType::Double, true),
        "float" => (FieldType::Float, true),
        "boolean" => (FieldType::Boolean, true),
        "char" => (FieldType::Char, true),
        "java.lang.Long" => (FieldType::Integer(64), false),
        "java.lang.Integer" => (FieldType::Integer(32), false),
        "java.lang.Short" => (FieldType::Integer(16), false),
        "java.lang.Byte" => (FieldType::Integer(8), false),
        "java.lang.Double" => (FieldType::Double, false),
        "java.lang.Float" => (FieldType::Float, false),
        "java.lang.Boolean" => (FieldType::Boolean, false),
        "java.lang.Character" => (FieldType::Char, false),
        "java.lang.String" | "htsjdk.samtools.util.Iso8601Date" => (FieldType::Text, false),
        other => (FieldType::Other(other), false),
    }
}

/// The table's escapes: `\\` and `\uXXXX`.
fn unescape(text: &str) -> String {
    let mut out = String::new();
    let mut chars = text.chars();
    while let Some(c) = chars.next() {
        if c != '\\' {
            out.push(c);
            continue;
        }
        match chars.next() {
            Some('u') => {
                let code: String = chars.by_ref().take(4).collect();
                if let Some(c) = u32::from_str_radix(&code, 16).ok().and_then(char::from_u32) {
                    out.push(c);
                }
            }
            Some(other) => out.push(other),
            None => out.push('\\'),
        }
    }
    out
}

/// The metrics class named exactly `name` (`Class.getName`, `$` before a
/// nested class).
pub fn metric_class(name: &str) -> Option<&'static MetricClass> {
    classes().metrics.get(name)
}

/// Whether `name` is a `Header` class htsjdk can instantiate
/// (`StringHeader`, `VersionHeader`).
pub fn is_header_class(name: &str) -> bool {
    classes().headers.contains(&name)
}

/// Any class the table or [`KEY_CLASSES`] knows: what `Class.forName`
/// finds here.
fn known(name: &str) -> Option<&'static str> {
    let classes = classes();
    classes
        .metrics
        .get_key_value(name)
        .map(|(k, _)| *k)
        .or_else(|| classes.headers.iter().copied().find(|h| *h == name))
        .or_else(|| KEY_CLASSES.iter().copied().find(|k| *k == name))
}

/// What [`load_class`] found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Lookup {
    /// The class's name.
    Found(&'static str),
    /// `ClassNotFoundException`.
    Missing,
    /// Another exception and its message: a name without a `.` fails on
    /// `substring(-1)`.
    Thrown(String),
}

/// `MetricsFile.loadClass(name, true)`.
pub fn load_class(name: &str) -> Lookup {
    if let Some(found) = known(name) {
        return Lookup::Found(found);
    }
    let Some(dot) = name.rfind('.') else {
        let length = name.encode_utf16().count();
        return Lookup::Thrown(format!(
            "Range [-1, {length}) out of bounds for length {length}"
        ));
    };
    let simple = &name[dot..];
    let inner = name.rfind('$').map(|i| &name[i + 1..]);
    for package in PACKAGES {
        if let Some(found) = known(&format!("{package}{simple}")) {
            return Lookup::Found(found);
        }
        if let Some(found) = inner.and_then(|inner| known(&format!("{package}.{inner}"))) {
            return Lookup::Found(found);
        }
    }
    Lookup::Missing
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::metric::MetricRow;

    #[test]
    fn insert_size_fields_are_in_get_fields_order() {
        let class = metric_class("picard.analysis.InsertSizeMetrics").expect("in the table");
        assert!(class.instantiable);
        let names: Vec<&str> = class.fields.iter().map(|f| f.name).collect();
        assert_eq!(
            names[..3],
            [
                "MEDIAN_INSERT_SIZE",
                "MODE_INSERT_SIZE",
                "MEDIAN_ABSOLUTE_DEVIATION"
            ]
        );
        assert_eq!(
            names[names.len() - 3..],
            ["SAMPLE", "LIBRARY", "READ_GROUP"]
        );
        let field = |name: &str| &class.fields[class.field(name).expect("a field")];
        assert_eq!(field("READ_PAIRS").ty, FieldType::Integer(64));
        assert!(field("READ_PAIRS").primitive);
        assert_eq!(field("MIN_INSERT_SIZE").ty, FieldType::Integer(32));
        assert_eq!(field("MEDIAN_INSERT_SIZE").default.as_deref(), Some("0.0"));
        assert_eq!(field("SAMPLE").ty, FieldType::Text);
        assert_eq!(field("SAMPLE").default, None);
        let FieldType::Enum { class, constants } = &field("PAIR_ORIENTATION").ty else {
            panic!("PAIR_ORIENTATION is an enum");
        };
        assert_eq!(*class, "htsjdk.samtools.SamPairUtil.PairOrientation");
        assert_eq!(constants, &["FR", "RF", "TANDEM"]);
    }

    #[test]
    fn every_class_bammetrics_writes_is_in_the_table() {
        for class in [
            crate::insert_size::Row::CLASS,
            crate::alignment_summary::Row::CLASS,
            crate::quality_yield::CLASS,
            crate::gc_bias::DetailRow::CLASS,
            crate::gc_bias::SummaryRow::CLASS,
            crate::duplicate_metrics::Row::CLASS,
            crate::base_distribution_by_cycle::Row::CLASS,
            "picard.analysis.directed.HsMetrics",
            "picard.analysis.directed.TargetedPcrMetrics",
            "picard.analysis.WgsMetrics",
            "picard.analysis.RnaSeqMetrics",
            "picard.vcf.CollectVariantCallingMetrics$VariantCallingDetailMetrics",
            "picard.vcf.GenotypeConcordanceSummaryMetrics",
        ] {
            assert!(metric_class(class).is_some(), "{class}");
        }
        assert!(is_header_class("htsjdk.samtools.metrics.StringHeader"));
        assert!(
            !metric_class("picard.analysis.MergeableMetricBase").is_some_and(|c| c.instantiable)
        );
    }

    #[test]
    fn legacy_names_resolve_as_htsjdk_resolves_them() {
        assert_eq!(
            load_class("net.sf.picard.analysis.AlignmentSummaryMetrics"),
            Lookup::Found("picard.analysis.AlignmentSummaryMetrics")
        );
        assert_eq!(
            load_class("old.CollectWgsMetrics$WgsMetrics"),
            Lookup::Found("picard.analysis.WgsMetrics")
        );
        assert_eq!(load_class("picard.analysis.NoSuchMetrics"), Lookup::Missing);
        assert_eq!(
            load_class("Integer"),
            Lookup::Thrown("Range [-1, 7) out of bounds for length 7".to_owned())
        );
        assert_eq!(
            load_class("java.lang.Integer"),
            Lookup::Found("java.lang.Integer")
        );
    }

    #[test]
    fn char_defaults_unescape() {
        let class =
            metric_class("picard.analysis.artifacts.ErrorSummaryMetrics").expect("in the table");
        let field = &class.fields[class.field("REF_BASE").expect("a field")];
        assert_eq!(field.ty, FieldType::Char);
        assert_eq!(field.default.as_deref(), Some("\0"));
    }
}

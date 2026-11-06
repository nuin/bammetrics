//! htsjdk `MetricsFile.read` for `CompareMetrics`: the class on the
//! `## METRICS CLASS` line from [`picard_classes`], each cell parsed by its
//! field's Java type (`FormatUtil.parseObject`), the histogram bins by the
//! key class, and every failure as the message of the exception htsjdk
//! throws.

use crate::java::{self, StringHashSet};
use crate::metric::Value;
use crate::metrics_file::format_double;
use crate::picard_classes::{self, Field, FieldType, Lookup, MetricClass};

const METRICS_LINE: &str = "## METRICS CLASS\t";
const HISTOGRAM_LINE: &str = "## HISTOGRAM\t";
const MAJOR_HEADER: &str = "## ";
const MINOR_HEADER: &str = "# ";
/// The NPE a missing column line raises (`readLine().split`).
const NO_COLUMN_LINE: &str = "Cannot invoke \"String.split(String)\" because the return value of \
                              \"java.io.BufferedReader.readLine()\" is null";

/// `String.split` with limit 0: trailing empty strings dropped, and a text
/// without the separator is itself.
pub(super) fn java_split(text: &str, separator: char) -> Vec<&str> {
    if !text.contains(separator) {
        return vec![text];
    }
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

/// A field value as `Field.get` boxes it.
#[derive(Clone, Debug)]
pub(super) enum Cell {
    Null,
    /// `Long`, `Integer`, `Short`.
    Int(i64),
    Double(f64),
    Float(f32),
    Bool(bool),
    /// `String`, an enum constant's name, a `Character`.
    Text(String),
}

/// [`Cell::key`].
#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq)]
pub(super) enum CellKey<'a> {
    Null,
    Int(i64),
    Double(u64),
    Float(u32),
    Bool(bool),
    Text(&'a str),
}

impl Cell {
    /// `FormatUtil.format(Object)`, what `MetricBase.equals` compares.
    pub(super) fn format(&self) -> String {
        match self {
            Cell::Null => String::new(),
            Cell::Int(v) => v.to_string(),
            Cell::Double(v) => format_double(*v),
            Cell::Float(v) => format_double(f64::from(*v)),
            Cell::Bool(v) => if *v { "Y" } else { "N" }.to_owned(),
            Cell::Text(v) => v.clone(),
        }
    }

    /// `String.valueOf`, as the difference lines print a value.
    pub(super) fn java_string(&self) -> String {
        match self {
            Cell::Null => "null".to_owned(),
            Cell::Int(v) => v.to_string(),
            Cell::Double(v) => java::double_to_string(*v),
            Cell::Float(v) => java::float_to_string(*v),
            Cell::Bool(v) => v.to_string(),
            Cell::Text(v) => v.clone(),
        }
    }

    /// `Number.doubleValue`, for a `Number`.
    pub(super) fn number(&self) -> Option<f64> {
        match self {
            Cell::Int(v) => Some(*v as f64),
            Cell::Double(v) => Some(*v),
            Cell::Float(v) => Some(f64::from(*v)),
            _ => None,
        }
    }

    /// `Object.equals` between values of one field (`Double.equals` takes
    /// NaN as equal to NaN and 0.0 as unequal to -0.0).
    pub(super) fn same(&self, other: &Cell) -> bool {
        match (self, other) {
            (Cell::Null, Cell::Null) => true,
            (Cell::Int(a), Cell::Int(b)) => a == b,
            (Cell::Double(a), Cell::Double(b)) => same_double(*a, *b),
            (Cell::Float(a), Cell::Float(b)) => {
                (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
            }
            (Cell::Bool(a), Cell::Bool(b)) => a == b,
            (Cell::Text(a), Cell::Text(b)) => a == b,
            _ => false,
        }
    }

    /// A hashable stand-in: two cells are [`Cell::same`] exactly when
    /// their keys are equal (every NaN is one key, `-0.0` and `0.0` two).
    pub(super) fn key(&self) -> CellKey<'_> {
        match self {
            Cell::Null => CellKey::Null,
            Cell::Int(v) => CellKey::Int(*v),
            Cell::Double(v) if v.is_nan() => CellKey::Double(f64::NAN.to_bits()),
            Cell::Double(v) => CellKey::Double(v.to_bits()),
            Cell::Float(v) if v.is_nan() => CellKey::Float(f32::NAN.to_bits()),
            Cell::Float(v) => CellKey::Float(v.to_bits()),
            Cell::Bool(v) => CellKey::Bool(*v),
            Cell::Text(v) => CellKey::Text(v),
        }
    }

    pub(super) fn value(&self) -> Value {
        match self {
            Cell::Null => Value::Null,
            Cell::Int(v) => Value::Int(*v),
            Cell::Double(v) => Value::Float(*v),
            Cell::Float(v) => Value::Float(f64::from(*v)),
            Cell::Bool(v) => Value::Bool(*v),
            Cell::Text(v) => Value::Str(v.clone()),
        }
    }
}

/// `Double.compare(a, b) == 0`.
pub(super) fn same_double(a: f64, b: f64) -> bool {
    (a.is_nan() && b.is_nan()) || a.to_bits() == b.to_bits()
}

/// `FormatUtil.parseDouble`: `?` and `-?` are NaN.
fn format_util_double(text: &str) -> std::result::Result<f64, String> {
    match text {
        "?" | "-?" => Ok(f64::NAN),
        _ => java::parse_double(text),
    }
}

/// `FormatUtil.parseBoolean`: true when the first character is `Y` or `y`.
fn parse_boolean(text: &str) -> bool {
    text.chars()
        .next()
        .is_some_and(|c| c.eq_ignore_ascii_case(&'Y'))
}

/// `FormatUtil.parseObject(text, field type)` then `Field.set`.
fn parse_cell(text: &str, field: &Field, class: &str) -> std::result::Result<Cell, String> {
    let setting = || {
        format!(
            "Error setting field {} on class of type {class}",
            field.name
        )
    };
    if text.is_empty() {
        return if field.primitive {
            Err(setting())
        } else {
            Ok(Cell::Null)
        };
    }
    match &field.ty {
        FieldType::Integer(16) => java::parse_short(text).map(Cell::Int),
        // `parseObject` gives an `Integer`, which a byte field rejects.
        FieldType::Integer(8) => java::parse_integer(text, 32).and(Err(setting())),
        FieldType::Integer(bits) => java::parse_integer(text, *bits).map(Cell::Int),
        FieldType::Double => format_util_double(text).map(Cell::Double),
        FieldType::Float => match text {
            "?" | "-?" => Ok(Cell::Float(f32::NAN)),
            _ => java::parse_float(text).map(Cell::Float),
        },
        FieldType::Boolean => Ok(Cell::Bool(parse_boolean(text))),
        FieldType::Char if text.encode_utf16().count() == 1 => Ok(Cell::Text(text.to_owned())),
        FieldType::Char => Err(format!(
            "Cannot parse string into char because length != 1 : {text}"
        )),
        FieldType::Text => Ok(Cell::Text(text.to_owned())),
        FieldType::Enum { class, constants } if !constants.contains(&text) => {
            Err(format!("No enum constant {class}.{text}"))
        }
        FieldType::Enum { .. } => Ok(Cell::Text(text.to_owned())),
        FieldType::Other(name) => Err(format!("Don't know how to convert a String to a {name}")),
    }
}

/// The value of `field` on a new instance: what a column the file lacks
/// holds.
fn default_cell(field: &Field) -> Cell {
    let Some(text) = &field.default else {
        return Cell::Null;
    };
    match field.ty {
        FieldType::Integer(_) => text.parse().map_or(Cell::Null, Cell::Int),
        FieldType::Double => java::parse_double(text).map_or(Cell::Null, Cell::Double),
        FieldType::Float => java::parse_float(text).map_or(Cell::Null, Cell::Float),
        FieldType::Boolean => Cell::Bool(text == "true"),
        _ => Cell::Text(text.clone()),
    }
}

/// One `## HISTOGRAM` column: htsjdk `Histogram` with its bins in file
/// order (equality does not depend on the order).
#[derive(Clone, Debug)]
pub(super) struct Histogram {
    bin_label: String,
    value_label: String,
    /// The class of the keys as `parseObject` returns them.
    key_class: &'static str,
    bins: Vec<(Cell, f64)>,
}

impl Histogram {
    fn increment(&mut self, key: &Cell, value: f64) {
        match self.bins.iter_mut().find(|(k, _)| k.same(key)) {
            Some(bin) => bin.1 += value,
            None => self.bins.push((key.clone(), 0.0 + value)),
        }
    }

    /// `Histogram.equals`: the labels, then the `TreeMap`s. Maps of two key
    /// classes are equal only when both are empty (`AbstractMap.equals`
    /// catches the `ClassCastException`).
    pub(super) fn equals(&self, other: &Histogram) -> bool {
        if self.bin_label != other.bin_label
            || self.value_label != other.value_label
            || self.bins.len() != other.bins.len()
        {
            return false;
        }
        if self.bins.is_empty() {
            return true;
        }
        self.key_class == other.key_class
            && self.bins.iter().all(|(key, value)| {
                other
                    .bins
                    .iter()
                    .any(|(k, v)| k.same(key) && same_double(*v, *value))
            })
    }
}

/// A metrics file as `MetricsFile.read` leaves it.
#[derive(Clone, Debug)]
pub(super) struct Doc {
    pub(super) class: Option<&'static MetricClass>,
    /// `columnLabels`, a `HashSet` filled in column order.
    pub(super) columns: StringHashSet,
    /// One cell per field of `class`, in `getFields` order.
    pub(super) rows: Vec<Vec<Cell>>,
    pub(super) histograms: Vec<Histogram>,
}

/// `BufferedReader.readLine`: lines end at `\n`, `\r` or `\r\n`.
pub(super) fn java_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut rest = text;
    while !rest.is_empty() {
        match rest.find(['\n', '\r']) {
            Some(end) => {
                lines.push(&rest[..end]);
                let skip = if rest[end..].starts_with("\r\n") {
                    2
                } else {
                    1
                };
                rest = &rest[end + skip..];
            }
            None => {
                lines.push(rest);
                rest = "";
            }
        }
    }
    lines
}

/// The class a `## METRICS CLASS` line names.
fn metric_class(name: &str) -> std::result::Result<&'static MetricClass, String> {
    let missing = || format!("Could not locate class with name {name}");
    match picard_classes::load_class(name) {
        Lookup::Found(found) => picard_classes::metric_class(found).ok_or_else(missing),
        Lookup::Missing => Err(missing()),
        Lookup::Thrown(message) => Err(message),
    }
}

/// `MetricsFile.read`. The error is the message of the exception htsjdk
/// throws, which `CompareMetrics` prints.
pub(super) fn read_doc(text: &str) -> std::result::Result<Doc, String> {
    let mut lines = java_lines(text).into_iter();
    let mut doc = Doc {
        class: None,
        columns: StringHashSet::new(),
        rows: Vec::new(),
        histograms: Vec::new(),
    };
    let mut line = read_headers(&mut lines)?;
    if let Some(start) = line.filter(|l| l.starts_with(METRICS_LINE)) {
        line = read_rows(start, &mut lines, &mut doc)?;
    }
    while let Some(l) = line {
        if java::trim(l).starts_with(MAJOR_HEADER) {
            break;
        }
        line = lines.next();
    }
    if let Some(start) = line
        .map(java::trim)
        .filter(|l| l.starts_with(HISTOGRAM_LINE))
    {
        read_histogram(start, &mut lines, &mut doc)?;
    }
    Ok(doc)
}

/// The header section: `## <Header class>` and `# <text>` pairs, blank
/// lines skipped, up to the trimmed line that starts the metrics or the
/// histogram (`None` at the end of the file).
fn read_headers<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
) -> std::result::Result<Option<&'a str>, String> {
    let mut open = false;
    for raw in lines.by_ref() {
        let line = java::trim(raw);
        if line.is_empty() {
            continue;
        }
        if line.starts_with(METRICS_LINE) || line.starts_with(HISTOGRAM_LINE) {
            return Ok(Some(line));
        }
        if let Some(class) = line.strip_prefix(MAJOR_HEADER) {
            if open {
                return Err("Consecutive header class lines encountered.".to_owned());
            }
            let class = java::trim(class);
            match picard_classes::load_class(class) {
                Lookup::Found(found) if picard_classes::is_header_class(found) => open = true,
                _ => {
                    return Err(format!(
                        "Error load and/or instantiating an instance of {class}"
                    ));
                }
            }
        } else if line.starts_with(MINOR_HEADER) {
            if !open {
                return Err(format!("Header class must precede header value:{line}"));
            }
            open = false;
        } else {
            return Err(format!(
                "Illegal state. Found following string in metrics file header: {line}"
            ));
        }
    }
    Ok(None)
}

/// The metrics block from its class line: the column line, then rows to a
/// blank line. Returns the line that ended it.
fn read_rows<'a>(
    start: &str,
    lines: &mut impl Iterator<Item = &'a str>,
    doc: &mut Doc,
) -> std::result::Result<Option<&'a str>, String> {
    let class = metric_class(java_split(start, '\t')[1])?;
    let header = lines.next().ok_or_else(|| NO_COLUMN_LINE.to_owned())?;
    let columns = java_split(header, '\t');
    for column in &columns {
        doc.columns.insert(column);
    }
    let slots = columns
        .iter()
        .map(|c| {
            class.field(c).ok_or_else(|| {
                format!(
                    "Could not get field with name {c} from class {}",
                    class.name
                )
            })
        })
        .collect::<std::result::Result<Vec<usize>, String>>()?;
    doc.class = Some(class);
    let defaults: Vec<Cell> = class.fields.iter().map(default_cell).collect();
    for raw in lines.by_ref() {
        if java::trim(raw).is_empty() {
            return Ok(Some(raw));
        }
        let values: Vec<&str> = raw.split('\t').collect();
        if !class.instantiable {
            return Err(format!("Error instantiating a {}", class.name));
        }
        let mut row = defaults.clone();
        for (i, &slot) in slots.iter().enumerate() {
            let text = values
                .get(i)
                .ok_or_else(|| format!("Index {i} out of bounds for length {}", values.len()))?;
            row[slot] = parse_cell(text, &class.fields[slot], class.name)?;
        }
        doc.rows.push(row);
    }
    Ok(None)
}

/// The histogram section from its trimmed `## HISTOGRAM` line: one
/// `Histogram` per label after the key's, rows to an empty line.
fn read_histogram<'a>(
    start: &str,
    lines: &mut impl Iterator<Item = &'a str>,
    doc: &mut Doc,
) -> std::result::Result<(), String> {
    let name = java::trim(java_split(start, '\t')[1]);
    let key_class = match picard_classes::load_class(name) {
        Lookup::Found(found) => found,
        Lookup::Missing => return Err(format!("Could not load class with name {name}")),
        Lookup::Thrown(message) => return Err(message),
    };
    let header = lines.next().ok_or_else(|| NO_COLUMN_LINE.to_owned())?;
    let labels = java_split(header, '\t');
    for label in labels.iter().skip(1) {
        doc.histograms.push(Histogram {
            bin_label: labels[0].to_owned(),
            value_label: (*label).to_owned(),
            key_class: key_runtime_class(key_class),
            bins: Vec::new(),
        });
    }
    for raw in lines {
        if raw.is_empty() {
            break;
        }
        let cells = java_split(java::trim(raw), '\t');
        let key = parse_key(cells[0], key_class)?;
        let count = doc.histograms.len();
        for (i, text) in cells.iter().enumerate().skip(1) {
            let value = format_util_double(text)?;
            let histogram = doc
                .histograms
                .get_mut(i - 1)
                .ok_or_else(|| format!("Index {} out of bounds for length {count}", i - 1))?;
            histogram.increment(&key, value);
        }
    }
    Ok(())
}

/// The class `parseObject` returns for a key of `class`: a `Byte` key is
/// parsed as an `Integer`.
fn key_runtime_class(class: &'static str) -> &'static str {
    if class == "java.lang.Byte" {
        "java.lang.Integer"
    } else {
        class
    }
}

/// A histogram key parsed by its class (`FormatUtil.parseObject`).
fn parse_key(text: &str, class: &str) -> std::result::Result<Cell, String> {
    match class {
        "java.lang.Short" => java::parse_short(text).map(Cell::Int),
        "java.lang.Integer" | "java.lang.Byte" => java::parse_integer(text, 32).map(Cell::Int),
        "java.lang.Long" => java::parse_integer(text, 64).map(Cell::Int),
        "java.lang.Double" => format_util_double(text).map(Cell::Double),
        "java.lang.Float" => match text {
            "?" | "-?" => Ok(Cell::Float(f32::NAN)),
            _ => java::parse_float(text).map(Cell::Float),
        },
        "java.lang.Boolean" => Ok(Cell::Bool(parse_boolean(text))),
        "java.lang.Character" if text.encode_utf16().count() != 1 => Err(format!(
            "Cannot parse string into char because length != 1 : {text}"
        )),
        "java.lang.Character" | "java.lang.String" => Ok(Cell::Text(text.to_owned())),
        other => Err(format!("Don't know how to convert a String to a {other}")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ASM: &str = "## htsjdk.samtools.metrics.StringHeader\n# CollectAlignmentSummaryMetrics\n\n\
        ## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\n\
        CATEGORY\tTOTAL_READS\tPCT_PF_READS\tMEAN_READ_LENGTH\tSAMPLE\n\
        FIRST_OF_PAIR\t5\t1\t100\t\n\
        PAIR\t10\t0.5\t100\t\n\n\
        ## HISTOGRAM\tjava.lang.Integer\nREAD_LENGTH\tPAIRED_TOTAL_LENGTH_COUNT\n100\t10\n\n";

    #[test]
    fn typed_cells_follow_the_field_types() {
        let doc = read_doc(ASM).expect("reads");
        let class = doc.class.expect("a class");
        assert_eq!(class.name, "picard.analysis.AlignmentSummaryMetrics");
        let cell =
            |row: usize, name: &str| doc.rows[row][class.field(name).expect("field")].clone();
        assert!(matches!(cell(0, "TOTAL_READS"), Cell::Int(5)));
        assert!(matches!(cell(1, "PCT_PF_READS"), Cell::Double(v) if v == 0.5));
        assert!(matches!(cell(0, "SAMPLE"), Cell::Null));
        // A column the file lacks holds the class default.
        assert!(matches!(cell(0, "PF_READS"), Cell::Int(0)));
        assert_eq!(doc.histograms.len(), 1);
        assert_eq!(doc.columns.len(), 5);
    }

    #[test]
    fn htsjdk_errors_are_kept_word_for_word() {
        let bad = |from: &str, to: &str| read_doc(&ASM.replacen(from, to, 1)).expect_err("fails");
        assert_eq!(bad("\t5\t", "\t5.0\t"), "For input string: \"5.0\"");
        assert_eq!(
            bad("\t5\t", "\t\t"),
            "Error setting field TOTAL_READS on class of type picard.analysis.AlignmentSummaryMetrics"
        );
        assert_eq!(
            bad("FIRST_OF_PAIR\t", "BOGUS\t"),
            "No enum constant picard.analysis.AlignmentSummaryMetrics.Category.BOGUS"
        );
        assert_eq!(
            bad("\tSAMPLE", "\tNOPE"),
            "Could not get field with name NOPE from class picard.analysis.AlignmentSummaryMetrics"
        );
        assert_eq!(
            bad("100\t\n", "100\n"),
            "Index 4 out of bounds for length 4"
        );
        assert_eq!(
            bad(
                "picard.analysis.AlignmentSummaryMetrics",
                "picard.analysis.NoSuch"
            ),
            "Could not locate class with name picard.analysis.NoSuch"
        );
        assert_eq!(
            bad("# CollectAlignmentSummaryMetrics\n", "# x\nstray\n"),
            "Illegal state. Found following string in metrics file header: stray"
        );
        // Cells past the last column are ignored.
        assert!(read_doc(&ASM.replacen("100\t\n", "100\t\textra\n", 1)).is_ok());
    }

    #[test]
    fn lines_end_like_buffered_reader() {
        assert_eq!(java_lines("a\nb\r\nc\rd"), ["a", "b", "c", "d"]);
        assert_eq!(java_lines("a\n\n"), ["a", ""]);
        assert_eq!(java_split("a\t\t", '\t'), ["a"]);
        assert_eq!(java_split("", '\t'), [""]);
        assert!(java_split("\t\t", '\t').is_empty());
        assert_eq!(java_split("X:0.5:", ':'), ["X", "0.5"]);
    }
}

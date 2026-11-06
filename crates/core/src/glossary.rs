//! Column definitions for the report: `glossary.tsv`, one row per column the
//! implemented tools write (metrics class, column, type, definition in our
//! own words). The file is compiled into the binary. `tests/glossary.rs`
//! runs every tool in the registry and fails when a column or histogram
//! label it writes has no row, or a row names a column its class does not
//! write.

/// The glossary text: `#` comment lines, then [`HEADER`], then one row per
/// column.
pub const GLOSSARY: &str = include_str!("../glossary.tsv");

const HEADER: &str = "class\tcolumn\ttype\tdefinition";

/// The lines that are rows: no comment, no header, not empty.
fn rows() -> impl Iterator<Item = &'static str> {
    GLOSSARY
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty() && *line != HEADER)
}

/// One glossary row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Entry {
    /// Picard metrics class, or the program name for a histogram column.
    pub class: &'static str,
    pub column: &'static str,
    /// `int`, `float` or `string`.
    pub kind: &'static str,
    pub definition: &'static str,
}

/// Every row of the glossary, in file order. A malformed line is skipped;
/// the test in this module fails on one.
pub fn entries() -> impl Iterator<Item = Entry> {
    rows().filter_map(|line| {
        let mut fields = line.splitn(4, '\t');
        Some(Entry {
            class: fields.next()?,
            column: fields.next()?,
            kind: fields.next()?,
            definition: fields.next()?,
        })
    })
}

/// The row for `column` of `class`, by the full class name or the name
/// after its last `.` or `$`.
pub fn lookup(class: &str, column: &str) -> Option<Entry> {
    entries().find(|e| {
        e.column == column && (e.class == class || e.class.rsplit(['.', '$']).next() == Some(class))
    })
}

/// The row for histogram column `label` of a file `program` wrote, where
/// `key` is the histogram's key label: the label itself, then `*.<suffix>`
/// for a unit-prefixed label, then `*`. `CollectMultipleMetrics` writes each
/// program's histogram under its own name, so when `program` has no row for
/// `key` the program whose histogram has that key is used.
pub fn histogram_lookup(program: &str, key: &str, label: &str) -> Option<Entry> {
    let program = if lookup(program, key).is_some() {
        program
    } else {
        entries()
            .find(|e| e.column == key && is_program(e.class))?
            .class
    };
    lookup(program, label)
        .or_else(|| {
            let suffix = label.split_once('.')?.1;
            lookup(program, &format!("*.{suffix}"))
        })
        .or_else(|| lookup(program, "*"))
}

/// A histogram row's first field: a program name, not a dotted class.
fn is_program(class: &str) -> bool {
    !class.contains('.')
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;

    #[test]
    fn every_line_is_a_row() {
        assert_eq!(entries().count(), rows().count());
        assert_eq!(GLOSSARY.lines().filter(|l| *l == HEADER).count(), 1);
        let mut seen = HashSet::new();
        for e in entries() {
            assert!(["int", "float", "string"].contains(&e.kind), "{e:?}");
            assert!(e.definition.len() > 10, "{e:?}");
            assert!(seen.insert((e.class, e.column)), "duplicate {e:?}");
        }
    }

    #[test]
    fn no_definition_has_an_en_or_em_dash() {
        for e in entries() {
            assert!(!e.definition.contains(['\u{2013}', '\u{2014}']), "{e:?}");
        }
    }

    #[test]
    fn lookup_takes_the_simple_class_name() {
        let entry = lookup("QualityYieldMetrics", "READ_LENGTH").expect("row");
        assert_eq!(entry.kind, "int");
        assert_eq!(
            lookup("picard.analysis.InsertSizeMetrics", "READ_PAIRS").map(|e| e.kind),
            Some("int")
        );
        assert_eq!(lookup("InsertSizeMetrics", "NO_SUCH"), None);
    }

    #[test]
    fn histogram_lookup_tries_the_label_then_the_suffix_then_star() {
        let column = |program, key, label| histogram_lookup(program, key, label).map(|e| e.column);
        assert_eq!(
            column("CollectInsertSizeMetrics", "insert_size", "insert_size"),
            Some("insert_size")
        );
        assert_eq!(
            column("CollectInsertSizeMetrics", "insert_size", "lib1.fr_count"),
            Some("*.fr_count")
        );
        assert_eq!(
            column(
                "EstimateLibraryComplexity",
                "duplication_group_count",
                "lib1"
            ),
            Some("*")
        );
        assert_eq!(
            column("CollectInsertSizeMetrics", "insert_size", "NO_SUCH"),
            None
        );
    }

    #[test]
    fn multiple_metrics_histograms_take_the_program_of_their_key() {
        for (key, label, program) in [
            (
                "insert_size",
                "All_Reads.fr_count",
                "CollectInsertSizeMetrics",
            ),
            (
                "READ_LENGTH",
                "PAIRED_TOTAL_LENGTH_COUNT",
                "CollectAlignmentSummaryMetrics",
            ),
            ("QUALITY", "COUNT_OF_Q", "QualityScoreDistribution"),
            ("CYCLE", "MEAN_QUALITY", "MeanQualityByCycle"),
        ] {
            let entry = histogram_lookup("CollectMultipleMetrics", key, label).expect(label);
            assert_eq!(entry.class, program);
        }
        assert_eq!(
            histogram_lookup("CollectMultipleMetrics", "NO_SUCH", "x"),
            None
        );
    }
}

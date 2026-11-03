//! Column definitions for the report: `glossary.tsv`, one row per column the
//! implemented tools write (metrics class, column, type, definition in our
//! own words). The file is compiled into the binary; a test checks that
//! every column a tool emits has a row.

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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;

    use super::*;
    use crate::metric::MetricRow;

    /// Every metrics class a tool writes, with its columns.
    fn emitted_classes() -> Vec<(&'static str, &'static [&'static str])> {
        vec![
            (
                crate::insert_size::Row::CLASS,
                crate::insert_size::Row::COLUMNS,
            ),
            (
                crate::alignment_summary::Row::CLASS,
                crate::alignment_summary::Row::COLUMNS,
            ),
            (
                crate::quality_yield::Row::CLASS,
                crate::quality_yield::Row::COLUMNS,
            ),
            (
                "picard.analysis.directed.HsMetrics",
                &crate::hs_metrics::COLUMNS,
            ),
            (
                crate::gc_bias::DetailRow::CLASS,
                crate::gc_bias::DetailRow::COLUMNS,
            ),
            (
                crate::gc_bias::SummaryRow::CLASS,
                crate::gc_bias::SummaryRow::COLUMNS,
            ),
            (
                crate::duplicate_metrics::Row::CLASS,
                crate::duplicate_metrics::Row::COLUMNS,
            ),
            (
                crate::base_distribution_by_cycle::Row::CLASS,
                crate::base_distribution_by_cycle::Row::COLUMNS,
            ),
        ]
    }

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
    fn every_emitted_column_has_a_definition() {
        let mut missing = Vec::new();
        for (class, columns) in emitted_classes() {
            for column in columns {
                if lookup(class, column).is_none() {
                    missing.push(format!("{class}\t{column}"));
                }
            }
        }
        assert!(
            missing.is_empty(),
            "no glossary row for:\n{}",
            missing.join("\n")
        );
    }

    #[test]
    fn no_row_names_a_column_the_class_does_not_emit() {
        let classes = emitted_classes();
        for e in entries() {
            if let Some((_, columns)) = classes.iter().find(|(class, _)| *class == e.class) {
                assert!(columns.contains(&e.column), "stale glossary row {e:?}");
            }
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
}

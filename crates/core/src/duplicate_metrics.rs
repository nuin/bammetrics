//! `CollectDuplicateMetrics`, mirroring Picard 3.5.0
//! (`CollectDuplicateMetrics`, `AbstractMarkDuplicatesCommandLineProgram`,
//! `LibraryIdGenerator` and `DuplicationMetrics`).
//!
//! The tool reads the duplicate flags already in the file. Rules:
//! * One row per library, in Java `TreeMap` (sorted name) order. Every `@RG`
//!   library gets a row before any record is read, so a library without
//!   records prints zeros. A group without `LB`, a record without `RG` or one
//!   naming a group missing from the header counts as `Unknown Library`.
//! * Per record: unmapped reads count in `UNMAPPED_READS`, then secondary
//!   and supplementary ones in `SECONDARY_OR_SUPPLEMENTARY_RDS`, then reads
//!   unpaired or with an unmapped mate in `UNPAIRED_READS_EXAMINED`, the rest
//!   in `READ_PAIRS_EXAMINED`. A duplicate flag adds to
//!   `UNPAIRED_READ_DUPLICATES` or `READ_PAIR_DUPLICATES` on the same rule,
//!   for mapped primary records only.
//! * Pair counts are per read and halved at the end (integer division).
//!   `READ_PAIR_OPTICAL_DUPLICATES` is always 0 because no optical
//!   duplicate is detected, so the library size is the one MarkDuplicates
//!   would give without optical duplicates.
//! * With exactly one library and a library size, the file carries the ROI
//!   histogram: `java.lang.Double` keys 1.0 to 100.0, `BIN` by
//!   `CoverageMult`.
//! * The output option is `METRICS_FILE` (`M`); Picard takes no `OUTPUT`.

use std::collections::{BTreeMap, HashMap};

use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::record::data::field::{Tag, Value as FieldValue};

use crate::duplication::{estimate_library_size, estimate_roi};
use crate::error::Result;
use crate::histogram::Histogram;
use crate::levels::read_groups;
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator, Collector};
use crate::registry::{Invocation, Opt, Tool};

/// `LibraryIdGenerator.UNKNOWN_LIBRARY`.
pub const UNKNOWN_LIBRARY: &str = "Unknown Library";

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Row {
    pub library: String,
    pub unpaired_reads_examined: i64,
    pub read_pairs_examined: i64,
    pub secondary_or_supplementary_rds: i64,
    pub unmapped_reads: i64,
    pub unpaired_read_duplicates: i64,
    pub read_pair_duplicates: i64,
    pub read_pair_optical_duplicates: i64,
    pub percent_duplication: f64,
    pub estimated_library_size: Option<i64>,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.sam.DuplicationMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "LIBRARY",
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS",
        "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES",
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.library.clone().into(),
            self.unpaired_reads_examined.into(),
            self.read_pairs_examined.into(),
            self.secondary_or_supplementary_rds.into(),
            self.unmapped_reads.into(),
            self.unpaired_read_duplicates.into(),
            self.read_pair_duplicates.into(),
            self.read_pair_optical_duplicates.into(),
            self.percent_duplication.into(),
            self.estimated_library_size.into(),
        ]
    }
}

impl Row {
    /// `finalizeAndWriteMetrics` and `calculateDerivedFields`.
    fn finish(mut self) -> Self {
        self.read_pairs_examined /= 2;
        self.read_pair_duplicates /= 2;
        self.estimated_library_size = estimate_library_size(
            self.read_pairs_examined - self.read_pair_optical_duplicates,
            self.read_pairs_examined - self.read_pair_duplicates,
        );
        let examined = self.unpaired_reads_examined + self.read_pairs_examined * 2;
        self.percent_duplication = if examined == 0 {
            0.0
        } else {
            (self.unpaired_read_duplicates + self.read_pair_duplicates * 2) as f64 / examined as f64
        };
        self
    }

    /// `calculateRoiHistogram`: `None` without a library size.
    fn roi_histogram(&self) -> Option<Histogram> {
        let size = self.estimated_library_size?;
        let unique = self.read_pairs_examined - self.read_pair_duplicates;
        let mut h = Histogram::new();
        for x in 1..=100 {
            h.increment_by(
                x,
                estimate_roi(size, x as f64, self.read_pairs_examined, unique),
            );
        }
        Some(h)
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectDuplicateMetrics",
    summary: "one row per library, from the duplicate flags already in the file; the metrics \
              file goes to METRICS_FILE.",
    options: &[
        Opt::output(&["METRICS_FILE", "M"]),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(_: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
    Ok(pass::boxed(DuplicateMetrics::new(header), |rows| {
        Ok(output(&rows).into())
    }))
}

/// The metrics file: the rows and, for exactly one library with a size,
/// the ROI histogram.
fn output(rows: &[Row]) -> MetricsOutput {
    let histogram = match rows {
        [only] => only.roi_histogram().map(|histogram| HistogramSection {
            key_label: "BIN",
            key_class: KeyClass::Double,
            series: vec![HistogramSeries {
                label: "CoverageMult".to_owned(),
                histogram,
            }],
        }),
        _ => None,
    };
    MetricsOutput::from_rows(rows, histogram)
}

/// The pass accumulator.
pub struct DuplicateMetrics {
    /// Library name of each `@RG` ID.
    libraries: HashMap<Vec<u8>, String>,
    rows: BTreeMap<String, Row>,
}

impl DuplicateMetrics {
    pub fn new(header: &sam::Header) -> Self {
        let mut libraries = HashMap::new();
        let mut rows = BTreeMap::new();
        for group in read_groups(header) {
            let library = group.library.unwrap_or_else(|| UNKNOWN_LIBRARY.to_owned());
            rows.entry(library.clone()).or_insert_with(|| Row {
                library: library.clone(),
                ..Row::default()
            });
            libraries.insert(group.id.into_bytes(), library);
        }
        Self { libraries, rows }
    }

    /// `LibraryIdGenerator.getLibraryName`: the `LB` of the record's group,
    /// `Unknown Library` without one.
    fn library(&self, record: &bam::Record) -> &str {
        match record.data().get(&Tag::READ_GROUP) {
            Some(Ok(FieldValue::String(id))) => self
                .libraries
                .get(AsRef::<[u8]>::as_ref(id))
                .map_or(UNKNOWN_LIBRARY, String::as_str),
            _ => UNKNOWN_LIBRARY,
        }
    }
}

impl Accumulator for DuplicateMetrics {
    type Output = Vec<Row>;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let library = self.library(record).to_owned();
        let row = self.rows.entry(library).or_insert_with_key(|name| Row {
            library: name.clone(),
            ..Row::default()
        });
        let flags = record.flags();
        let primary = !flags.is_secondary() && !flags.is_supplementary();
        let unpaired = !flags.is_segmented() || flags.is_mate_unmapped();
        if flags.is_unmapped() {
            row.unmapped_reads += 1;
        } else if !primary {
            row.secondary_or_supplementary_rds += 1;
        } else if unpaired {
            row.unpaired_reads_examined += 1;
        } else {
            row.read_pairs_examined += 1;
        }
        if flags.is_duplicate() && primary && !flags.is_unmapped() {
            if unpaired {
                row.unpaired_read_duplicates += 1;
            } else {
                row.read_pair_duplicates += 1;
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<Vec<Row>> {
        Ok(self.rows.into_values().map(Row::finish).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pairs_are_halved_and_percent_counts_reads() {
        let row = Row {
            unpaired_reads_examined: 4,
            read_pairs_examined: 21,
            unpaired_read_duplicates: 1,
            read_pair_duplicates: 5,
            ..Row::default()
        }
        .finish();
        assert_eq!((row.read_pairs_examined, row.read_pair_duplicates), (10, 2));
        assert_eq!(row.percent_duplication, 5.0 / 24.0);
        assert!(row.estimated_library_size.is_some());
    }

    #[test]
    fn no_reads_no_size_and_zero_percent() {
        let row = Row::default().finish();
        assert_eq!(row.percent_duplication, 0.0);
        assert_eq!(row.estimated_library_size, None);
        assert!(row.roi_histogram().is_none());
    }
}

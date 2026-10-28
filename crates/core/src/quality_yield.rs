//! `CollectQualityYieldMetrics`, mirroring Picard 3.5.0
//! `CollectQualityYieldMetrics.QualityYieldMetricsCollector`.
//!
//! One row for the whole file. Every record counts, mapped or not, except
//! secondary and supplementary alignments unless the options ask for them.
//! Bases are counted from the record length, qualities from whichever
//! quality array is in use, so the two can differ (see `quality_array`).

use std::io::{self, Write};
use std::path::Path;

use noodles::bam;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};
use crate::metrics_file::MetricsFile;
use crate::pass::{self, Accumulator};
use crate::record::record_name;

/// `QualityYieldMetrics` columns in Picard 3.5.0 file order (the order the
/// fields are declared in, which `MetricsFile` reads by reflection).
pub const COLUMNS: [&str; 11] = [
    "TOTAL_READS",
    "PF_READS",
    "READ_LENGTH",
    "TOTAL_BASES",
    "PF_BASES",
    "Q20_BASES",
    "PF_Q20_BASES",
    "Q30_BASES",
    "PF_Q30_BASES",
    "Q20_EQUIVALENT_YIELD",
    "PF_Q20_EQUIVALENT_YIELD",
];

/// Quality at or above which a base counts as Q20 and Q30. A Q30 base is
/// also a Q20 base.
const Q20: i8 = 20;
const Q30: i8 = 30;
/// `Q20_EQUIVALENT_YIELD` is the sum of base qualities divided by this.
const Q20_EQUIVALENT_DIVISOR: i64 = 20;
/// Offset of a Phred quality in the SAM text of an `OQ` tag.
const PHRED_OFFSET: u8 = 33;
/// BAM stores a missing `QUAL` (`*`) as 0xFF bytes.
const MISSING_QUALITY: u8 = 0xFF;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Options {
    /// `USE_ORIGINAL_QUALITIES`: read the `OQ` tag instead of `QUAL` when
    /// the record has a non-empty one.
    pub use_original_qualities: bool,
    /// `INCLUDE_SECONDARY_ALIGNMENTS`.
    pub include_secondary_alignments: bool,
    /// `INCLUDE_SUPPLEMENTAL_ALIGNMENTS`.
    pub include_supplemental_alignments: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            use_original_qualities: true,
            include_secondary_alignments: false,
            include_supplemental_alignments: false,
        }
    }
}

/// The `QualityYieldMetrics` row. Java `long` fields are `i64`; `READ_LENGTH`
/// is a Java `int`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Row {
    pub total_reads: i64,
    pub pf_reads: i64,
    pub read_length: i32,
    pub total_bases: i64,
    pub pf_bases: i64,
    pub q20_bases: i64,
    pub pf_q20_bases: i64,
    pub q30_bases: i64,
    pub pf_q30_bases: i64,
    pub q20_equivalent_yield: i64,
    pub pf_q20_equivalent_yield: i64,
}

impl Row {
    fn picard_fields(&self) -> Vec<String> {
        vec![
            self.total_reads.to_string(),
            self.pf_reads.to_string(),
            self.read_length.to_string(),
            self.total_bases.to_string(),
            self.pf_bases.to_string(),
            self.q20_bases.to_string(),
            self.pf_q20_bases.to_string(),
            self.q30_bases.to_string(),
            self.pf_q30_bases.to_string(),
            self.q20_equivalent_yield.to_string(),
            self.pf_q20_equivalent_yield.to_string(),
        ]
    }
}

/// Write the row as a Picard `QualityYieldMetrics` file (no histogram).
pub fn write_metrics(row: &Row, headers: Vec<String>, out: &mut impl Write) -> io::Result<()> {
    MetricsFile {
        headers,
        class: "picard.analysis.CollectQualityYieldMetrics$QualityYieldMetrics",
        columns: &COLUMNS,
        rows: vec![row.picard_fields()],
        histogram: None,
    }
    .write(out)
}

/// Running totals. The two `*_q20_equivalent_yield` fields hold quality sums
/// until `finish`, which divides them as `QualityYieldMetricsCollector.finish` does.
#[derive(Default)]
struct Collector {
    row: Row,
}

impl Collector {
    /// `acceptRecord` after the secondary and supplementary checks. `length`
    /// is `getReadLength`; `qualities` may be longer or shorter than it.
    fn add(&mut self, length: i64, pass_filter: bool, qualities: &[i8]) {
        let row = &mut self.row;
        row.total_reads += 1;
        row.total_bases += length;
        if pass_filter {
            row.pf_reads += 1;
            row.pf_bases += length;
        }
        for &quality in qualities {
            let value = i64::from(quality);
            let (q20, q30) = (i64::from(quality >= Q20), i64::from(quality >= Q30));
            row.q20_equivalent_yield += value;
            row.q20_bases += q20;
            row.q30_bases += q30;
            if pass_filter {
                row.pf_q20_equivalent_yield += value;
                row.pf_q20_bases += q20;
                row.pf_q30_bases += q30;
            }
        }
    }

    /// `finish` and `calculateDerivedFields`: integer divisions, truncating
    /// like Java's `long` division. `READ_LENGTH` is 0 for an empty file.
    fn finish(self) -> Row {
        let mut row = self.row;
        row.q20_equivalent_yield /= Q20_EQUIVALENT_DIVISOR;
        row.pf_q20_equivalent_yield /= Q20_EQUIVALENT_DIVISOR;
        row.read_length = row
            .total_bases
            .checked_div(row.total_reads)
            .map_or(0, |mean| mean as i32);
        row
    }
}

/// The qualities Picard loops over for a record: the `OQ` tag when
/// original qualities are wanted and the tag is a non-empty string
/// (`SAMRecord.getOriginalBaseQualities`), else `QUAL`. Its length is not
/// checked against the read length. A missing `QUAL` yields no qualities
/// (htsjdk `NULL_QUALS`), so the read's bases count but no quality does.
fn quality_array(
    record: &bam::Record,
    use_original_qualities: bool,
    scratch: &mut Vec<i8>,
) -> Result<()> {
    scratch.clear();
    if use_original_qualities {
        match record.data().get(&Tag::ORIGINAL_QUALITY_SCORES) {
            None => {}
            Some(Err(e)) => return Err(e.into()),
            // htsjdk casts the attribute to String, which fails for other types.
            Some(Ok(Value::String(oq))) => {
                if !oq.is_empty() {
                    scratch.extend(oq.iter().map(|c| c.wrapping_sub(PHRED_OFFSET) as i8));
                    return Ok(());
                }
            }
            Some(Ok(_)) => {
                return Err(Error::MalformedRecord {
                    name: record_name(record),
                    what: "OQ tag is not a string",
                });
            }
        }
    }
    let quality = record.quality_scores();
    let bytes = quality.as_bytes();
    if bytes.first() != Some(&MISSING_QUALITY) {
        scratch.extend(bytes.iter().map(|&q| q as i8));
    }
    Ok(())
}

/// `CollectQualityYieldMetrics` as a pass [`Accumulator`].
pub struct QualityYield {
    opts: Options,
    collector: Collector,
    qualities: Vec<i8>,
}

impl QualityYield {
    pub fn new(opts: &Options) -> Self {
        Self {
            opts: *opts,
            collector: Collector::default(),
            qualities: Vec::new(),
        }
    }
}

impl Accumulator for QualityYield {
    type Output = Row;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if (flags.is_secondary() && !self.opts.include_secondary_alignments)
            || (flags.is_supplementary() && !self.opts.include_supplemental_alignments)
        {
            return Ok(());
        }
        quality_array(
            record,
            self.opts.use_original_qualities,
            &mut self.qualities,
        )?;
        self.collector.add(
            record.sequence().len() as i64,
            !flags.is_qc_fail(),
            &self.qualities,
        );
        Ok(())
    }

    fn finish(self) -> Result<Row> {
        Ok(self.collector.finish())
    }
}

/// `StringHeader` command line in Picard's shape.
pub fn command_header(opts: &Options, bam: &Path, output: Option<&Path>) -> String {
    let mut command = format!("CollectQualityYieldMetrics INPUT={}", bam.display());
    if let Some(output) = output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    command.push_str(&format!(
        " USE_ORIGINAL_QUALITIES={} INCLUDE_SECONDARY_ALIGNMENTS={} INCLUDE_SUPPLEMENTAL_ALIGNMENTS={}",
        opts.use_original_qualities,
        opts.include_secondary_alignments,
        opts.include_supplemental_alignments
    ));
    command
}

/// Compute the metrics of one BAM with the default BGZF worker count.
pub fn run_over_bam(path: impl AsRef<Path>, opts: &Options) -> Result<Row> {
    pass::run(path.as_ref(), pass::default_threads(), |_| {
        Ok(QualityYield::new(opts))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_and_equivalent_yield() {
        let mut c = Collector::default();
        // Q19 is neither, Q20 is Q20 only, Q30 is both; sum 19+20+30+93 = 162.
        c.add(4, true, &[19, 20, 30, 93]);
        let row = c.finish();
        assert_eq!((row.q20_bases, row.q30_bases), (3, 2));
        assert_eq!(row.q20_equivalent_yield, 162 / 20);
        assert_eq!(row.pf_q20_equivalent_yield, 8);
    }

    #[test]
    fn qc_fail_reads_count_in_totals_only() {
        let mut c = Collector::default();
        c.add(3, false, &[40, 40, 40]);
        c.add(2, true, &[10, 10]);
        let row = c.finish();
        assert_eq!((row.total_reads, row.pf_reads), (2, 1));
        assert_eq!((row.total_bases, row.pf_bases), (5, 2));
        assert_eq!((row.q20_bases, row.pf_q20_bases), (3, 0));
        assert_eq!(row.pf_q20_equivalent_yield, 1);
    }

    #[test]
    fn read_length_truncates_and_is_zero_when_empty() {
        let mut c = Collector::default();
        c.add(100, true, &[]);
        c.add(75, true, &[]);
        assert_eq!(c.finish().read_length, 87);
        assert_eq!(Collector::default().finish(), Row::default());
    }

    #[test]
    fn quality_array_length_is_independent_of_read_length() {
        let mut c = Collector::default();
        c.add(10, true, &[30; 25]);
        let row = c.finish();
        assert_eq!((row.total_bases, row.q30_bases), (10, 25));
    }

    #[test]
    fn row_prints_in_column_order() {
        let row = Row {
            total_reads: 2,
            read_length: 7,
            ..Row::default()
        };
        let mut out = Vec::new();
        write_metrics(&row, vec!["h".to_owned()], &mut out).expect("writes");
        let text = String::from_utf8(out).expect("utf8");
        assert!(text.contains(&COLUMNS.join("\t")));
        assert!(text.contains("\n2\t0\t7\t0\t0\t0\t0\t0\t0\t0\t0\n"));
    }
}

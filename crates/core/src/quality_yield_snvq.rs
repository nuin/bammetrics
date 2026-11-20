//! `CollectQualityYieldMetricsSNVQ`, mirroring Picard 3.5.0
//! `CollectQualityYieldMetricsSNVQ.QualityYieldMetricsCollector`.
//!
//! One row for the whole file, a `SinglePassSamProgram`. A record counts
//! as in `CollectQualityYieldMetrics`: every record, mapped or not, except
//! secondary and supplementary alignments unless the options ask for them
//! ([`crate::quality_yield::skipped`]).
//!
//! Base qualities come from `QUAL`, or from the Phred+33 text of
//! `ALTERNATE_QUALITY_ATTRIBUTE` when it is set (no `QUAL` fallback: a
//! record without the tag fails, as Picard's `NullPointerException`
//! does). Their count must equal the base count. A quality of 40 or more
//! counts in Q20, Q30 and Q40, one of 30 to 39 in Q20 and Q30, one of 20
//! to 29 in Q20.
//!
//! SNV qualities are the Phred+33 texts of the four tags `qa`, `qc`, `qg`
//! and `qt`, one character per base in the order the bases are stored (a
//! reverse-strand record is not reversed). Every record read needs all
//! four (a missing one fails the run, Picard's `NullPointerException`). At
//! each base, each of `A`, `C`, `G`, `T` the base is not adds one SNVQ
//! observation from its tag: three at an `A`, `C`, `G` or `T`, four at an
//! `N` or any other base. Lower-case bases do not reach the tool: htsjdk
//! and the BAM encoding both store them upper case. The thresholds cascade
//! as for the bases. A tag is read only where a base asks for it, so a tag
//! shorter than the read fails only when a missing position is read, and a
//! longer one is never checked.
//!
//! The texts are the bytes Java's `String.getBytes` gives in UTF-8: a
//! character from `0x80` up takes two bytes, and every byte outside `!` to
//! `~` fails when it is read (`SAMUtils.fastqToPhred`).
//!
//! `INCLUDE_BQ_HISTOGRAM=true` adds one `## HISTOGRAM` section keyed by
//! `KEY` with four series, each over its own keys and printed on the union:
//! `BQ_COUNT` (bases by quality), `SNVQ_COUNT` (observations by SNV
//! quality), `READ_INDEX_MEAN_BQ` and `READ_INDEX_MEAN_SNVQ` (the mean at
//! each read index, from 0). The means are `SeriesStats.getMean`: the sum
//! as a `double` over the count.

use std::borrow::Cow;

use noodles::sam::alignment::record::data::field::{Tag, Value as TagValue};
use noodles::{bam, sam};

use crate::bam_raw::HtsjdkTags;
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::{
    HistogramSection, HistogramSeries, MetricRow, MetricsOutput, ToolOutput, Value,
};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator};
use crate::quality::base_qualities;
use crate::quality_yield::{read_length, skipped};
use crate::record::record_name;
use crate::registry::{Invocation, Opt, Tool};

/// The class name Picard writes after `## METRICS CLASS`.
pub const CLASS: &str = "picard.analysis.CollectQualityYieldMetricsSNVQ$QualityYieldMetrics";

/// The columns in Picard 3.5.0 file order (declaration order).
pub const COLUMNS: [&str; 31] = [
    "TOTAL_READS",
    "PF_READS",
    "READ_LENGTH",
    "TOTAL_BASES",
    "PF_BASES",
    "Q20_BASES",
    "PF_Q20_BASES",
    "Q30_BASES",
    "PF_Q30_BASES",
    "Q40_BASES",
    "PF_Q40_BASES",
    "PCT_Q20_BASES",
    "PCT_Q30_BASES",
    "PCT_Q40_BASES",
    "PCT_PF_Q20_BASES",
    "PCT_PF_Q30_BASES",
    "PCT_PF_Q40_BASES",
    "TOTAL_SNVQ",
    "PF_SNVQ",
    "Q20_SNVQ",
    "PF_Q20_SNVQ",
    "Q30_SNVQ",
    "PF_Q30_SNVQ",
    "Q40_SNVQ",
    "PF_Q40_SNVQ",
    "PCT_Q20_SNVQ",
    "PCT_Q30_SNVQ",
    "PCT_Q40_SNVQ",
    "PCT_PF_Q20_SNVQ",
    "PCT_PF_Q30_SNVQ",
    "PCT_PF_Q40_SNVQ",
];

/// The histogram labels, in the order Picard adds the histograms.
pub const HISTOGRAM_KEY: &str = "KEY";
pub const SERIES: [&str; 4] = [
    "BQ_COUNT",
    "SNVQ_COUNT",
    "READ_INDEX_MEAN_BQ",
    "READ_INDEX_MEAN_SNVQ",
];

/// `CollectQualityYieldMetricsSNVQ.SNVQ_BASE_ORDER`: the bases whose tags
/// hold the SNV qualities, `q` and the base in lower case.
const BASES: [u8; 4] = *b"ACGT";
const TAGS: [Tag; 4] = [
    Tag::new(b'q', b'a'),
    Tag::new(b'q', b'c'),
    Tag::new(b'q', b'g'),
    Tag::new(b'q', b't'),
];

/// `SAMUtils.fastqToPhred`: the offset and the accepted characters.
const PHRED_OFFSET: u8 = 33;
const PHRED_LAST: u8 = 126;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Options {
    /// `ALTERNATE_QUALITY_ATTRIBUTE`: the tag read for base qualities in
    /// place of `QUAL`, as given (htsjdk checks its length per record).
    pub alternate_quality_attribute: Option<String>,
    pub include_secondary_alignments: bool,
    pub include_supplemental_alignments: bool,
    pub include_bq_histogram: bool,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        Ok(Self {
            alternate_quality_attribute: inv
                .get("ALTERNATE_QUALITY_ATTRIBUTE")
                .filter(|value| *value != "null")
                .map(str::to_owned),
            include_secondary_alignments: inv.flag("INCLUDE_SECONDARY_ALIGNMENTS", false)?,
            include_supplemental_alignments: inv.flag("INCLUDE_SUPPLEMENTAL_ALIGNMENTS", false)?,
            include_bq_histogram: inv.flag("INCLUDE_BQ_HISTOGRAM", false)?,
        })
    }
}

/// Counts of one kind (bases or SNV observations): every value and those
/// at or above 20, 30 and 40, overall and in PF reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    pub total: i64,
    pub pf: i64,
    pub q20: i64,
    pub pf_q20: i64,
    pub q30: i64,
    pub pf_q30: i64,
    pub q40: i64,
    pub pf_q40: i64,
}

impl Tally {
    /// One value: Q40 counts as Q30 and Q20, Q30 as Q20.
    fn add(&mut self, quality: i32, pass_filter: bool) {
        let (q20, q30, q40) = (
            i64::from(quality >= 20),
            i64::from(quality >= 30),
            i64::from(quality >= 40),
        );
        self.q20 += q20;
        self.q30 += q30;
        self.q40 += q40;
        if pass_filter {
            self.pf_q20 += q20;
            self.pf_q30 += q30;
            self.pf_q40 += q40;
        }
    }

    /// `PCT_Q20`, `PCT_Q30`, `PCT_Q40` and their PF forms; 0 over nothing.
    fn fractions(&self) -> [f64; 6] {
        let over = |count: i64, total: i64| {
            if total == 0 {
                0.0
            } else {
                count as f64 / total as f64
            }
        };
        [
            over(self.q20, self.total),
            over(self.q30, self.total),
            over(self.q40, self.total),
            over(self.pf_q20, self.pf),
            over(self.pf_q30, self.pf),
            over(self.pf_q40, self.pf),
        ]
    }
}

/// The `QualityYieldMetrics` row of the SNVQ tool.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Row {
    pub total_reads: i64,
    pub pf_reads: i64,
    /// `TOTAL_BASES` and `PF_BASES` are `bases.total` and `bases.pf`.
    pub bases: Tally,
    pub snvq: Tally,
}

impl MetricRow for Row {
    const CLASS: &'static str = CLASS;
    const COLUMNS: &'static [&'static str] = &COLUMNS;

    fn values(&self) -> Vec<Value> {
        let counts = |t: &Tally| [t.q20, t.pf_q20, t.q30, t.pf_q30, t.q40, t.pf_q40];
        let mut values: Vec<Value> = vec![
            self.total_reads.into(),
            self.pf_reads.into(),
            i64::from(read_length(self.bases.total, self.total_reads)).into(),
            self.bases.total.into(),
            self.bases.pf.into(),
        ];
        values.extend(counts(&self.bases).map(Value::from));
        values.extend(self.bases.fractions().map(Value::from));
        values.push(self.snvq.total.into());
        values.push(self.snvq.pf.into());
        values.extend(counts(&self.snvq).map(Value::from));
        values.extend(self.snvq.fractions().map(Value::from));
        values
    }
}

/// The running sum and count of one read index (`SeriesStats`).
#[derive(Clone, Copy, Debug, Default)]
struct Mean {
    sum: i64,
    count: u64,
}

/// A quality histogram and the per-index means (`INCLUDE_BQ_HISTOGRAM`).
#[derive(Debug, Default)]
struct ByIndex {
    counts: Histogram,
    means: Vec<Mean>,
}

impl ByIndex {
    fn add(&mut self, index: usize, quality: i32) {
        self.counts.increment(i64::from(quality));
        if self.means.len() <= index {
            self.means.resize(index + 1, Mean::default());
        }
        let mean = &mut self.means[index];
        mean.sum += i64::from(quality);
        mean.count += 1;
    }

    /// `READ_INDEX_MEAN_*`: one bin per index, the mean as its value.
    fn mean_histogram(&self) -> Histogram {
        let mut histogram = Histogram::new();
        for (index, mean) in self.means.iter().enumerate() {
            let value = if mean.count == 0 {
                f64::NAN
            } else {
                mean.sum as f64 / mean.count as f64
            };
            histogram.increment_by(index as i64, value);
        }
        histogram
    }
}

/// What the pass returns: the row and, with the option, the histograms.
#[derive(Debug, Default)]
pub struct Report {
    pub row: Row,
    bq: Option<ByIndex>,
    snvq: Option<ByIndex>,
}

impl Report {
    /// The `## HISTOGRAM` section, `None` without the option.
    pub fn histogram(&self) -> Option<HistogramSection> {
        let (bq, snvq) = (self.bq.as_ref()?, self.snvq.as_ref()?);
        let series = [
            bq.counts.clone(),
            snvq.counts.clone(),
            bq.mean_histogram(),
            snvq.mean_histogram(),
        ];
        Some(HistogramSection {
            key_label: HISTOGRAM_KEY,
            key_class: KeyClass::Integer,
            series: SERIES
                .iter()
                .zip(series)
                .map(|(label, histogram)| HistogramSeries {
                    label: (*label).to_owned(),
                    histogram,
                })
                .collect(),
        })
    }

    pub fn into_output(self) -> MetricsOutput {
        let histogram = self.histogram();
        MetricsOutput::from_rows(&[self.row], histogram)
    }
}

/// The bytes of a `Z` tag as Java's `String.getBytes` gives them: htsjdk
/// holds one character per byte, and UTF-8 writes a character from `0x80`
/// up as two bytes.
fn java_bytes(text: &[u8]) -> Cow<'_, [u8]> {
    if text.is_ascii() {
        return Cow::Borrowed(text);
    }
    let mut out = Vec::with_capacity(text.len() * 2);
    for &byte in text {
        if byte.is_ascii() {
            out.push(byte);
        } else {
            out.push(0xC0 | (byte >> 6));
            out.push(0x80 | (byte & 0x3F));
        }
    }
    Cow::Owned(out)
}

/// `SAMUtils.fastqToPhred(char)`.
fn phred(byte: u8) -> Option<i32> {
    (PHRED_OFFSET..=PHRED_LAST)
        .contains(&byte)
        .then(|| i32::from(byte - PHRED_OFFSET))
}

/// The two characters of a tag name, or `None` when htsjdk's
/// `SAMTag.makeBinaryTag` would refuse it.
fn tag_of(name: &str) -> Option<Tag> {
    let bytes = name.as_bytes();
    (bytes.len() == 2).then(|| Tag::new(bytes[0], bytes[1]))
}

/// `CollectQualityYieldMetricsSNVQ` as a pass [`Accumulator`].
pub struct QualityYieldSnvq {
    opts: Options,
    alternate: Option<Option<Tag>>,
    report: Report,
    qualities: Vec<i8>,
    bases: Vec<u8>,
    tags: HtsjdkTags,
}

impl QualityYieldSnvq {
    pub fn new(opts: &Options) -> Self {
        let histogram = || opts.include_bq_histogram.then(ByIndex::default);
        Self {
            alternate: opts.alternate_quality_attribute.as_deref().map(tag_of),
            opts: opts.clone(),
            report: Report {
                bq: histogram(),
                snvq: histogram(),
                ..Report::default()
            },
            qualities: Vec::new(),
            bases: Vec::new(),
            tags: HtsjdkTags::default(),
        }
    }

    /// The text of a `Z` tag `getStringAttribute` returns, failing as
    /// Picard fails on a record without it or with another type. The field
    /// is the one htsjdk's tag list keeps ([`HtsjdkTags`], loaded for
    /// `record`): the last of a repeated tag, or none past a lost append.
    fn string_tag<'r>(tags: &HtsjdkTags, record: &'r bam::Record, tag: Tag) -> Result<&'r [u8]> {
        let name = || String::from_utf8_lossy(tag.as_ref()).into_owned();
        let value = match tags.position(*tag.as_ref()) {
            Some(index) => match record.data().iter().nth(index) {
                Some(field) => Some(field?.1),
                None => None,
            },
            None => None,
        };
        match value {
            Some(TagValue::String(text)) => Ok(text.as_ref()),
            Some(_) => Err(Error::Metric(format!(
                "record {}: Value for tag {} is not a String",
                record_name(record),
                name()
            ))),
            None => Err(Error::Metric(format!(
                "record {}: no {} tag; CollectQualityYieldMetricsSNVQ needs it on every record \
                 it counts (Picard fails with a NullPointerException)",
                record_name(record),
                name()
            ))),
        }
    }

    /// The base qualities: `QUAL`, or the alternate tag decoded.
    fn read_qualities(&mut self, record: &bam::Record) -> Result<()> {
        let Some(alternate) = &self.alternate else {
            base_qualities(record, &mut self.qualities);
            return Ok(());
        };
        let given = self
            .opts
            .alternate_quality_attribute
            .as_deref()
            .unwrap_or("");
        let tag = alternate.ok_or_else(|| {
            Error::Metric(format!("String tag does not have length() == 2: {given}"))
        })?;
        let text = java_bytes(Self::string_tag(&self.tags, record, tag)?);
        self.qualities.clear();
        for &byte in text.iter() {
            let quality = phred(byte).ok_or_else(|| {
                Error::Metric(format!(
                    "record {}: Invalid fastq character in the {given} tag",
                    record_name(record)
                ))
            })?;
            self.qualities.push(quality as i8);
        }
        Ok(())
    }
}

impl Accumulator for QualityYieldSnvq {
    type Output = Report;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if skipped(
            flags,
            self.opts.include_secondary_alignments,
            self.opts.include_supplemental_alignments,
        ) {
            return Ok(());
        }
        let length = record.sequence().len() as i64;
        let pass_filter = !flags.is_qc_fail();
        let row = &mut self.report.row;
        row.total_reads += 1;
        row.bases.total += length;
        if pass_filter {
            row.pf_reads += 1;
            row.bases.pf += length;
        }
        self.tags.try_load(
            record
                .data()
                .iter()
                .map(|field| field.map(|(tag, _)| *tag.as_ref())),
        )?;
        self.read_qualities(record)?;
        self.bases.clear();
        self.bases.extend(record.sequence().iter());
        if self.qualities.len() != self.bases.len() {
            return Err(Error::Metric(format!(
                "record {}: quality string length does not match bases string",
                record_name(record)
            )));
        }
        let mut tags: [Cow<'_, [u8]>; 4] = Default::default();
        for (text, tag) in tags.iter_mut().zip(TAGS) {
            *text = java_bytes(Self::string_tag(&self.tags, record, tag)?);
        }
        let report = &mut self.report;
        for (index, (&quality, &base)) in self.qualities.iter().zip(&self.bases).enumerate() {
            let quality = i32::from(quality);
            report.row.bases.add(quality, pass_filter);
            for (text, &other) in tags.iter().zip(&BASES) {
                if base == other {
                    continue;
                }
                let snvq = text.get(index).copied().and_then(phred).ok_or_else(|| {
                    Error::Metric(format!(
                        "record {}: no SNV quality at read index {index} of the q{} tag",
                        record_name(record),
                        char::from(other.to_ascii_lowercase())
                    ))
                })?;
                report.row.snvq.total += 1;
                if pass_filter {
                    report.row.snvq.pf += 1;
                }
                report.row.snvq.add(snvq, pass_filter);
                if let Some(by_index) = report.snvq.as_mut() {
                    by_index.add(index, snvq);
                }
            }
            if let Some(by_index) = report.bq.as_mut() {
                by_index.add(index, quality);
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<Report> {
        Ok(self.report)
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectQualityYieldMetricsSNVQ",
    summary: "one row of read, base and Q20/Q30/Q40 totals, and the same for the SNV \
              qualities of the qa, qc, qg and qt tags (three per A, C, G or T base, four per \
              other base); INCLUDE_BQ_HISTOGRAM adds the quality histograms and the mean by \
              read index.",
    options: &[
        Opt::text(
            &["ALTERNATE_QUALITY_ATTRIBUTE", "AQA"],
            "Read base qualities from this Phred+33 tag in place of QUAL; every record needs \
             it (default null: QUAL).",
        ),
        Opt::bool(
            &["INCLUDE_BQ_HISTOGRAM"],
            "Add the BQ and SNVQ histograms and the means by read index (default false).",
        ),
        Opt::bool(
            &["INCLUDE_SECONDARY_ALIGNMENTS"],
            "Count secondary alignments (default false).",
        ),
        Opt::bool(
            &["INCLUDE_SUPPLEMENTAL_ALIGNMENTS"],
            "Count supplementary alignments (default false).",
        ),
        crate::registry::STOP_AFTER_OPTION,
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, _: &sam::Header) -> Result<Box<dyn pass::Collector>> {
    let opts = Options::from_invocation(inv)?;
    Ok(pass::boxed(QualityYieldSnvq::new(&opts), |report| {
        Ok(report.into_output().into())
    }))
}

#[cfg(test)]
mod tests;

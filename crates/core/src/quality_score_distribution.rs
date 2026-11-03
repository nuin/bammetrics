//! `QualityScoreDistribution`, mirroring Picard 3.5.0.
//!
//! The file has no metrics class, only a `java.lang.Byte` histogram:
//! `QUALITY` by `COUNT_OF_Q` and, when any counted base has an `OQ` value,
//! `COUNT_OF_OQ`. Rules from `acceptRead` and `finish`:
//! * Secondary and supplementary records never count. `PF_READS_ONLY` drops
//!   QC-failed reads, `ALIGNED_READS_ONLY` unmapped ones.
//! * Each base offset below the `QUAL` length counts once, unless the base is
//!   a no-call (`N`, `n`, `.`) and `INCLUDE_NO_CALLS` is false. The same
//!   offset of `OQ` counts too when the record carries one, so `OQ` past the
//!   `QUAL` length is ignored and a shorter `OQ` fails (Picard throws an
//!   index error).
//! * Only qualities with a non-zero count print.

use noodles::{bam, sam};

use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput, ToolOutput};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator, Collector};
use crate::quality::{base_qualities, original_qualities};
use crate::record::record_name;
use crate::registry::{Invocation, Opt, Tool};
use crate::sequence::is_no_call;

/// Picard counts into `long[128]` indexed by the Java byte.
const QUALITY_SLOTS: usize = 128;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Options {
    pub aligned_reads_only: bool,
    pub pf_reads_only: bool,
    pub include_no_calls: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        Ok(Self {
            aligned_reads_only: inv.flag("ALIGNED_READS_ONLY", false)?,
            pf_reads_only: inv.flag("PF_READS_ONLY", false)?,
            include_no_calls: inv.flag("INCLUDE_NO_CALLS", false)?,
        })
    }
}

pub const TOOL: Tool = Tool {
    name: "QualityScoreDistribution",
    options: &[
        Opt::ignored(&["CHART_OUTPUT", "CHART"], "no chart is drawn"),
        Opt::value(&["ALIGNED_READS_ONLY"]),
        Opt::value(&["PF_READS_ONLY", "PF"]),
        Opt::value(&["INCLUDE_NO_CALLS"]),
        Opt::value(&["STOP_AFTER"]),
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, _: &sam::Header) -> Result<Box<dyn Collector>> {
    let opts = Options::from_invocation(inv)?;
    Ok(pass::boxed(
        QualityScoreDistribution::new(opts),
        |section| Ok(MetricsOutput::histogram_only(section).into()),
    ))
}

/// The pass accumulator.
pub struct QualityScoreDistribution {
    opts: Options,
    q_counts: [i64; QUALITY_SLOTS],
    oq_counts: [i64; QUALITY_SLOTS],
    quals: Vec<i8>,
    oq: Vec<i8>,
}

impl QualityScoreDistribution {
    pub fn new(opts: Options) -> Self {
        Self {
            opts,
            q_counts: [0; QUALITY_SLOTS],
            oq_counts: [0; QUALITY_SLOTS],
            quals: Vec::new(),
            oq: Vec::new(),
        }
    }
}

/// The `long[128]` slot of a Java byte quality; a negative one is the index
/// error Picard throws.
fn slot(record: &bam::Record, quality: i8) -> Result<usize> {
    usize::try_from(quality).map_err(|_| Error::MalformedRecord {
        name: record_name(record),
        what: "negative base quality",
    })
}

impl Accumulator for QualityScoreDistribution {
    type Output = HistogramSection;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if (self.opts.pf_reads_only && flags.is_qc_fail())
            || (self.opts.aligned_reads_only && flags.is_unmapped())
            || flags.is_secondary()
            || flags.is_supplementary()
        {
            return Ok(());
        }
        base_qualities(record, &mut self.quals);
        let has_oq = original_qualities(record, &mut self.oq)?;
        let bases = record.sequence();
        for (i, &quality) in self.quals.iter().enumerate() {
            let base = bases.get(i).ok_or_else(|| Error::MalformedRecord {
                name: record_name(record),
                what: "fewer bases than qualities",
            })?;
            if !self.opts.include_no_calls && is_no_call(base) {
                continue;
            }
            self.q_counts[slot(record, quality)?] += 1;
            if has_oq {
                let original = *self.oq.get(i).ok_or_else(|| Error::MalformedRecord {
                    name: record_name(record),
                    what: "OQ is shorter than the qualities",
                })?;
                self.oq_counts[slot(record, original)?] += 1;
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<HistogramSection> {
        let histogram = |counts: &[i64; QUALITY_SLOTS]| {
            let mut h = Histogram::new();
            for (quality, &count) in counts.iter().enumerate() {
                if count > 0 {
                    h.increment_by(quality as i64, count as f64);
                }
            }
            h
        };
        let mut series = vec![HistogramSeries {
            label: "COUNT_OF_Q".to_owned(),
            histogram: histogram(&self.q_counts),
        }];
        let oq = histogram(&self.oq_counts);
        if !oq.is_empty() {
            series.push(HistogramSeries {
                label: "COUNT_OF_OQ".to_owned(),
                histogram: oq,
            });
        }
        Ok(HistogramSection {
            key_label: "QUALITY",
            key_class: KeyClass::Byte,
            series,
        })
    }
}

//! `MeanQualityByCycle`, mirroring Picard 3.5.0 (`MeanQualityByCycle` and
//! `HistogramGenerator`).
//!
//! The file has no metrics class, only a `java.lang.Integer` histogram:
//! `CYCLE` by `MEAN_QUALITY` and, when any record has an `OQ` tag,
//! `MEAN_ORIGINAL_QUALITY` (0 where a cycle has no `OQ` mean). Rules:
//! * Secondary and supplementary records never count. `PF_READS_ONLY` drops
//!   QC-failed reads, `ALIGNED_READS_ONLY` unmapped ones.
//! * `QUAL` and `OQ` are accumulated separately, each over its own length;
//!   records without `OQ` add nothing to the `OQ` side.
//! * Cycles follow [`crate::cycle`]: reverse reads are flipped, unpaired
//!   reads are first end, and second-end cycles print after the last
//!   first-end cycle whose quality sum is positive. A first-end cycle whose
//!   qualities are all 0 does not print and moves the second end down.

use noodles::{bam, sam};

use crate::cycle::{ByCycle, is_second_end, machine_cycle};
use crate::error::Result;
use crate::histogram::Histogram;
use crate::metric::{HistogramSection, HistogramSeries, MetricsOutput, ToolOutput};
use crate::metrics_file::KeyClass;
use crate::pass::{self, Accumulator, Collector};
use crate::quality::{base_qualities, original_qualities};
use crate::registry::{Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "MeanQualityByCycle",
    options: &[
        Opt::ignored(&["CHART_OUTPUT", "CHART"], "no chart is drawn"),
        Opt::value(&["ALIGNED_READS_ONLY"]),
        Opt::value(&["PF_READS_ONLY"]),
        Opt::value(&["STOP_AFTER"]),
    ],
    run: run_tool,
    collector: Some(collector),
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    inv.run_collector(collector)
}

fn collector(inv: &Invocation, _: &sam::Header) -> Result<Box<dyn Collector>> {
    let aligned_reads_only = inv.flag("ALIGNED_READS_ONLY", false)?;
    let pf_reads_only = inv.flag("PF_READS_ONLY", false)?;
    Ok(pass::boxed(
        MeanQualityByCycle::new(aligned_reads_only, pf_reads_only),
        |section| Ok(MetricsOutput::histogram_only(section).into()),
    ))
}

/// Quality sum and base count of one cycle.
#[derive(Clone, Copy, Debug, Default)]
struct Cell {
    total: f64,
    count: i64,
}

/// One `HistogramGenerator`: `QUAL` or `OQ`.
#[derive(Default)]
struct Generator {
    cycles: ByCycle<Cell>,
}

impl Generator {
    fn add(&mut self, quals: &[i8], reverse: bool, second: bool) {
        self.cycles.reserve(quals.len());
        for (i, &q) in quals.iter().enumerate() {
            let cell = self
                .cycles
                .cell(second, machine_cycle(reverse, quals.len(), i));
            cell.total += f64::from(q);
            cell.count += 1;
        }
    }

    /// `getMeanQualityHistogram`.
    fn mean_histogram(&self) -> Histogram {
        let mut h = Histogram::new();
        for (_, cycle, cell) in self.cycles.printed(|c| c.total > 0.0, |c| c.count > 0) {
            h.increment_by(cycle as i64, cell.total / cell.count as f64);
        }
        h
    }
}

/// The pass accumulator.
pub struct MeanQualityByCycle {
    aligned_reads_only: bool,
    pf_reads_only: bool,
    q: Generator,
    oq: Generator,
    quals: Vec<i8>,
}

impl MeanQualityByCycle {
    pub fn new(aligned_reads_only: bool, pf_reads_only: bool) -> Self {
        Self {
            aligned_reads_only,
            pf_reads_only,
            q: Generator::default(),
            oq: Generator::default(),
            quals: Vec::new(),
        }
    }
}

impl Accumulator for MeanQualityByCycle {
    type Output = HistogramSection;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        if (self.pf_reads_only && flags.is_qc_fail())
            || (self.aligned_reads_only && flags.is_unmapped())
            || flags.is_secondary()
            || flags.is_supplementary()
        {
            return Ok(());
        }
        let (reverse, second) = (flags.is_reverse_complemented(), is_second_end(flags));
        base_qualities(record, &mut self.quals);
        self.q.add(&self.quals, reverse, second);
        if original_qualities(record, &mut self.quals)? {
            self.oq.add(&self.quals, reverse, second);
        }
        Ok(())
    }

    fn finish(self) -> Result<HistogramSection> {
        let mut series = vec![HistogramSeries {
            label: "MEAN_QUALITY".to_owned(),
            histogram: self.q.mean_histogram(),
        }];
        if !self.oq.cycles.is_empty() {
            series.push(HistogramSeries {
                label: "MEAN_ORIGINAL_QUALITY".to_owned(),
                histogram: self.oq.mean_histogram(),
            });
        }
        Ok(HistogramSection {
            key_label: "CYCLE",
            key_class: KeyClass::Integer,
            series,
        })
    }
}

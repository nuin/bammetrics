//! `CollectBaseDistributionByCycle`, mirroring Picard 3.5.0.
//!
//! One `BaseDistributionByCycleMetrics` row per read end and cycle with at
//! least one base. Rules:
//! * Secondary and supplementary records never count. `PF_READS_ONLY` drops
//!   QC-failed reads, `ALIGNED_READS_ONLY` unmapped ones.
//! * `A`, `C`, `G` and `T` (either case) count as themselves; every other
//!   symbol, `=` and IUPAC codes included, counts as N.
//! * `PCT_x` is `100.0 * count / bases at the cycle`, so 0 to 100.
//! * Cycles follow [`crate::cycle`]; `READ_END` 2 cycles print after the last
//!   first-end cycle with a base.
//! * No row at all (an empty file) leaves out the `## METRICS CLASS` block.

use noodles::bam;

use crate::cycle::{ByCycle, is_second_end, machine_cycle};
use crate::error::Result;
use crate::metric::{MetricRow, MetricsOutput, Value};
use crate::pass::Accumulator;
use crate::registry::{Invocation, Opt, Tool};

pub struct Row {
    pub read_end: i64,
    pub cycle: i64,
    pub pct_a: f64,
    pub pct_c: f64,
    pub pct_g: f64,
    pub pct_t: f64,
    pub pct_n: f64,
}

impl MetricRow for Row {
    const CLASS: &'static str = "picard.analysis.BaseDistributionByCycleMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "READ_END", "CYCLE", "PCT_A", "PCT_C", "PCT_G", "PCT_T", "PCT_N",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.read_end.into(),
            self.cycle.into(),
            self.pct_a.into(),
            self.pct_c.into(),
            self.pct_g.into(),
            self.pct_t.into(),
            self.pct_n.into(),
        ]
    }
}

pub const TOOL: Tool = Tool {
    name: "CollectBaseDistributionByCycle",
    options: &[
        Opt::ignored(&["CHART_OUTPUT", "CHART"], "no chart is drawn"),
        Opt::value(&["ALIGNED_READS_ONLY"]),
        Opt::value(&["PF_READS_ONLY"]),
        Opt::value(&["STOP_AFTER"]),
    ],
    run: run_tool,
};

fn run_tool(inv: &Invocation) -> Result<Option<MetricsOutput>> {
    let aligned_reads_only = inv.flag("ALIGNED_READS_ONLY", false)?;
    let pf_reads_only = inv.flag("PF_READS_ONLY", false)?;
    let rows = inv.single_pass(|_| {
        Ok(BaseDistributionByCycle::new(
            aligned_reads_only,
            pf_reads_only,
        ))
    })?;
    Ok(Some(MetricsOutput::from_rows(&rows, None)))
}

/// Counts of A, C, G, T and N at one cycle.
type Counts = [i64; 5];

/// Index of a base in [`Counts`].
fn base_index(base: u8) -> usize {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// The pass accumulator.
pub struct BaseDistributionByCycle {
    aligned_reads_only: bool,
    pf_reads_only: bool,
    cycles: ByCycle<Counts>,
}

impl BaseDistributionByCycle {
    pub fn new(aligned_reads_only: bool, pf_reads_only: bool) -> Self {
        Self {
            aligned_reads_only,
            pf_reads_only,
            cycles: ByCycle::default(),
        }
    }
}

impl Accumulator for BaseDistributionByCycle {
    type Output = Vec<Row>;

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
        let sequence = record.sequence();
        let length = sequence.len();
        self.cycles.reserve(length);
        for (i, base) in sequence.iter().enumerate() {
            self.cycles.cell(second, machine_cycle(reverse, length, i))[base_index(base)] += 1;
        }
        Ok(())
    }

    fn finish(self) -> Result<Vec<Row>> {
        let any = |c: &Counts| c.iter().sum::<i64>() != 0;
        Ok(self
            .cycles
            .printed(any, any)
            .into_iter()
            .map(|(end, cycle, counts)| {
                let total = counts.iter().sum::<i64>() as f64;
                let pct = |i: usize| 100.0 * counts[i] as f64 / total;
                Row {
                    read_end: end.number(),
                    cycle: cycle as i64,
                    pct_a: pct(0),
                    pct_c: pct(1),
                    pct_g: pct(2),
                    pct_t: pct(3),
                    pct_n: pct(4),
                }
            })
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn only_acgt_are_bases() {
        assert_eq!(base_index(b'a'), 0);
        assert_eq!(base_index(b'T'), 3);
        for other in *b"N=R." {
            assert_eq!(base_index(other), 4);
        }
    }
}

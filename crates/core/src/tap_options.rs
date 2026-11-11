//! The options of a one-pass tool, parsed without its BAM.
//!
//! A tool's `collector` parses its options and then builds its accumulator
//! from the header of the stream. A chain builds a tap when the stream of
//! its position starts, which behind a transform is after the first read.
//! [`check`] runs the first half alone, so a value only the tool can judge
//! (`METRIC_ACCUMULATION_LEVEL=NOPE`, `STRAND_SPECIFICITY=UP`) is a usage
//! error before the first read. What a tool judges against the header or
//! against a file it reads (a contig the header lacks, an interval list
//! that does not parse) is still found when the tap is built.
//!
//! Each row calls the parser the tool's own `collector` calls, so the two
//! cannot disagree. No row opens a file. A row also runs what the tool
//! checks without a header after parsing: `FLOW_MODE=true`
//! (`quality_yield::checked_options`), and for
//! `CollectSequencingArtifactMetrics` the stand-alone validation and the
//! context count (`CONTEXT_SIZE=99`, `MINIMUM_INSERT_SIZE=-1`), which the
//! verb runs in `run_tool` and the collector alone does not.

use crate::error::Result;
use crate::registry::{Invocation, Tool};
use crate::{
    alignment_summary, base_distribution_by_cycle, coverage_gaps, duplicate_metrics, gc_bias,
    hs_metrics, insert_size, mean_quality_by_cycle, quality_score_distribution, quality_yield,
    rna_seq_metrics, sequencing_artifact, sex_check,
};

type Check = fn(&Invocation) -> Result<()>;

/// The two flags of the by-cycle tools.
fn read_filters(inv: &Invocation) -> Result<()> {
    inv.flag("ALIGNED_READS_ONLY", false)?;
    inv.flag("PF_READS_ONLY", false)?;
    Ok(())
}

/// One row per registry tool with a `collector`: its name and its parser.
const CHECKS: &[(&str, Check)] = &[
    (alignment_summary::TOOL.name, |inv| {
        alignment_summary::Options::from_invocation(inv).map(drop)
    }),
    (base_distribution_by_cycle::TOOL.name, read_filters),
    (coverage_gaps::TOOL.name, |inv| {
        coverage_gaps::Options::from_invocation(inv).map(drop)
    }),
    (duplicate_metrics::TOOL.name, |_| Ok(())),
    (gc_bias::TOOL.name, |inv| {
        gc_bias::Options::from_invocation(inv).map(drop)
    }),
    (hs_metrics::TOOL.name, |inv| {
        hs_metrics::Options::from_invocation(inv, hs_metrics::Variant::Hs).map(drop)
    }),
    (insert_size::TOOL.name, |inv| {
        insert_size::Options::from_invocation(inv).map(drop)
    }),
    (mean_quality_by_cycle::TOOL.name, read_filters),
    (quality_score_distribution::TOOL.name, |inv| {
        quality_score_distribution::Options::from_invocation(inv).map(drop)
    }),
    (quality_yield::TOOL.name, |inv| {
        quality_yield::checked_options(inv).map(drop)
    }),
    (rna_seq_metrics::TOOL.name, |inv| {
        rna_seq_metrics::Options::from_invocation(inv).map(drop)
    }),
    (sequencing_artifact::TOOL.name, |inv| {
        sequencing_artifact::Options::from_invocation(inv)?.check()
    }),
    (sex_check::TOOL.name, |inv| {
        sex_check::Options::from_invocation(inv).map(drop)
    }),
];

/// Parse the options of `tool` as its collector would, and build nothing.
/// A tool without a row passes: its collector parses when it is built.
pub fn check(tool: &Tool, inv: &Invocation) -> Result<()> {
    CHECKS
        .iter()
        .find(|(name, _)| *name == tool.name)
        .map_or(Ok(()), |(_, check)| check(inv))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error::Error;
    use crate::registry::TOOLS;

    /// A new one-pass tool gets a row, or a chain finds its bad values
    /// only when the tap is built.
    #[test]
    fn every_one_pass_tool_has_a_row() {
        let mut one_pass: Vec<&str> = TOOLS
            .iter()
            .filter(|tool| tool.collector.is_some())
            .map(|tool| tool.name)
            .collect();
        one_pass.sort_unstable();
        let mut rows: Vec<&str> = CHECKS.iter().map(|(name, _)| *name).collect();
        rows.sort_unstable();
        assert_eq!(rows, one_pass);
    }

    fn parsed(tool: &Tool, args: &[&str]) -> Invocation {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        tool.parse(&args).expect("the registry takes the text")
    }

    #[test]
    fn a_value_only_the_tool_can_judge_is_a_usage_error_without_a_bam() {
        const ARTIFACT: [&str; 3] = ["I=no.bam", "O=no", "R=no.fa"];
        let artifact = |option: &'static str| [&ARTIFACT[..], &[option]].concat();
        let (too_many, negative) = (artifact("CONTEXT_SIZE=99"), artifact("CONTEXT_SIZE=-1"));
        let (insert, order) = (
            artifact("MINIMUM_INSERT_SIZE=-1"),
            artifact("MAXIMUM_INSERT_SIZE=10"),
        );
        let cases: [(&Tool, &[&str], &str); 11] = [
            (
                &sequencing_artifact::TOOL,
                &too_many,
                "CONTEXT_SIZE=99: too many contexts to count",
            ),
            (
                &sequencing_artifact::TOOL,
                &negative,
                "CONTEXT_SIZE cannot be negative",
            ),
            (
                &sequencing_artifact::TOOL,
                &insert,
                "MINIMUM_INSERT_SIZE cannot be negative",
            ),
            (
                &sequencing_artifact::TOOL,
                &order,
                "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE",
            ),
            (
                &quality_yield::TOOL,
                &["I=no.bam", "FLOW_MODE=true"],
                "FLOW_MODE is obsolete",
            ),
            (
                &gc_bias::TOOL,
                &["I=no.bam", "R=no.fa", "S=no.txt", "SCAN_WINDOW_SIZE=0"],
                "SCAN_WINDOW_SIZE=0: a window holds at least 1 base",
            ),
            (
                &gc_bias::TOOL,
                &["I=no.bam", "R=no.fa", "S=no.txt", "WINDOW_SIZE=-5"],
                "SCAN_WINDOW_SIZE=-5: a window holds at least 1 base",
            ),
            (
                &insert_size::TOOL,
                &["I=no.bam", "METRIC_ACCUMULATION_LEVEL=NOPE"],
                "METRIC_ACCUMULATION_LEVEL \"NOPE\"",
            ),
            (
                &alignment_summary::TOOL,
                &["I=no.bam", "EXPECTED_PAIR_ORIENTATIONS=UP"],
                "UP",
            ),
            (
                &rna_seq_metrics::TOOL,
                &["I=no.bam", "REF_FLAT=no.txt", "STRAND_SPECIFICITY=UP"],
                "UP",
            ),
            (
                &sex_check::TOOL,
                &["I=no.bam", "EXPECTED_SEX=NOPE"],
                "EXPECTED_SEX=NOPE",
            ),
        ];
        for (tool, args, text) in cases {
            let error = check(tool, &parsed(tool, args)).expect_err(tool.name);
            assert!(matches!(error, Error::Usage(_)), "{}: {error}", tool.name);
            assert!(error.to_string().contains(text), "{}: {error}", tool.name);
        }
    }

    #[test]
    fn good_options_pass_and_no_file_is_opened() {
        let inv = parsed(
            &insert_size::TOOL,
            &["I=no.bam", "METRIC_ACCUMULATION_LEVEL=SAMPLE"],
        );
        check(&insert_size::TOOL, &inv).expect("parses");
        let inv = parsed(
            &rna_seq_metrics::TOOL,
            &[
                "I=no.bam",
                "REF_FLAT=no.txt",
                "STRAND_SPECIFICITY=NONE",
                "RIBOSOMAL_INTERVALS=no.interval_list",
            ],
        );
        check(&rna_seq_metrics::TOOL, &inv).expect("no file is opened");
        let inv = parsed(&duplicate_metrics::TOOL, &["I=no.bam"]);
        check(&duplicate_metrics::TOOL, &inv).expect("no option to parse");
        // The largest context size counted, and the smallest window.
        let inv = parsed(
            &sequencing_artifact::TOOL,
            &["I=no.bam", "O=no", "R=no.fa", "CONTEXT_SIZE=5"],
        );
        check(&sequencing_artifact::TOOL, &inv).expect("4^11 contexts are counted");
        let inv = parsed(
            &gc_bias::TOOL,
            &["I=no.bam", "R=no.fa", "S=no.txt", "SCAN_WINDOW_SIZE=1"],
        );
        check(&gc_bias::TOOL, &inv).expect("a window of 1 base");
    }
}

//! Picard-style command lines (`KEY=VALUE` arguments) mapped onto the
//! native flags, so a pipeline can swap `picard CollectInsertSizeMetrics ...`
//! for `bammetrics CollectInsertSizeMetrics ...`.
//!
//! Picard's CLI parser (barclay) takes `LONG=value` and `SHORT=value`
//! (`I=`, `O=`, `M=`, `W=`, `MW=`). Every option falls in one of three
//! groups: it maps to a native flag, it never changes the metrics and is
//! dropped with a note, or it changes the metrics and is not implemented,
//! which is an error naming it. A name outside all three is an error too,
//! because dropping an option Picard would have applied gives wrong numbers.

use crate::{Error, Result};

/// What an option does once recognised.
enum Action {
    /// `NAME=value` becomes `--flag=value`.
    Value(&'static str),
    /// `NAME=true|false` becomes `--flag` or nothing.
    Bool(&'static str),
    /// The only value this implementation can honour; others are errors.
    Only(&'static str),
    /// Never changes the metrics; accepted and dropped with this reason.
    Ignore(&'static str),
}

/// Option spellings (long name first) and their action.
type Entry = (&'static [&'static str], Action);

/// Options every Picard tool takes from `CommandLineProgram` and
/// `PicardCommandLine`: runtime plumbing that cannot change a metric.
const COMMON: &[Entry] = &[
    (
        &["VALIDATION_STRINGENCY"],
        Action::Ignore("records are not validated"),
    ),
    (
        &["TMP_DIR"],
        Action::Ignore("no temporary files are written"),
    ),
    (&["QUIET"], Action::Ignore("progress is not logged")),
    (&["VERBOSITY"], Action::Ignore("progress is not logged")),
    (&["COMPRESSION_LEVEL"], Action::Ignore("no BAM is written")),
    (
        &["MAX_RECORDS_IN_RAM"],
        Action::Ignore("records are streamed"),
    ),
    (&["CREATE_INDEX"], Action::Ignore("no BAM is written")),
    (&["CREATE_MD5_FILE"], Action::Ignore("no BAM is written")),
    (&["GA4GH_CLIENT_SECRETS"], Action::Ignore("no GA4GH access")),
    (
        &["USE_JDK_DEFLATER", "USE_JDK_INFLATER"],
        Action::Ignore("no JDK codec"),
    ),
    (
        &["REFERENCE_SEQUENCE", "R"],
        Action::Ignore("this metric does not read the reference"),
    ),
];

/// Options of `CollectQualityYieldMetrics` (`SinglePassSamProgram` and
/// `CollectQualityYieldMetrics` in Picard 3.5.0). The metric does not read
/// the reference, so `R` stays ignored through `COMMON`.
const QUALITY_YIELD: &[Entry] = &[
    (&["INPUT", "I"], Action::Value("--input")),
    (&["OUTPUT", "O"], Action::Value("--output")),
    (
        &["USE_ORIGINAL_QUALITIES", "OQ"],
        Action::Value("--use-original-qualities"),
    ),
    (
        &["INCLUDE_SECONDARY_ALIGNMENTS"],
        Action::Bool("--include-secondary-alignments"),
    ),
    (
        &["INCLUDE_SUPPLEMENTAL_ALIGNMENTS"],
        Action::Bool("--include-supplemental-alignments"),
    ),
    (&["FLOW_MODE"], Action::Only("false")),
    (
        &["ASSUME_SORTED", "AS"],
        Action::Ignore("one pass over the file in any order"),
    ),
];

/// Options of `CollectHsMetrics` (`CollectTargetedMetrics` and
/// `CollectHsMetrics` in Picard 3.5.0). `R` changes `AT_DROPOUT` and
/// `GC_DROPOUT`, so this entry shadows the `COMMON` one. The two interval
/// options repeat, one file per occurrence. `PER_BASE_COVERAGE` is written;
/// the other extra outputs (`PER_TARGET_COVERAGE`,
/// `THEORETICAL_SENSITIVITY_OUTPUT`) are not, so they are errors.
const HS_METRICS: &[Entry] = &[
    (&["INPUT", "I"], Action::Value("--input")),
    (&["OUTPUT", "O"], Action::Value("--output")),
    (&["BAIT_INTERVALS", "BI"], Action::Value("--bait-intervals")),
    (
        &["TARGET_INTERVALS", "TI"],
        Action::Value("--target-intervals"),
    ),
    (&["BAIT_SET_NAME", "N"], Action::Value("--bait-set-name")),
    (&["REFERENCE_SEQUENCE", "R"], Action::Value("--reference")),
    (
        &["MINIMUM_MAPPING_QUALITY", "MQ"],
        Action::Value("--minimum-mapping-quality"),
    ),
    (
        &["MINIMUM_BASE_QUALITY", "Q"],
        Action::Value("--minimum-base-quality"),
    ),
    (
        &["CLIP_OVERLAPPING_READS"],
        Action::Value("--clip-overlapping-reads"),
    ),
    (&["INCLUDE_INDELS"], Action::Bool("--include-indels")),
    (&["COVERAGE_CAP", "covMax"], Action::Value("--coverage-cap")),
    (&["NEAR_DISTANCE"], Action::Value("--near-distance")),
    (&["SAMPLE_SIZE"], Action::Value("--sample-size")),
    (&["PER_BASE_COVERAGE"], Action::Value("--per-base-coverage")),
    (
        &["ALLELE_FRACTION"],
        Action::Ignore("only THEORETICAL_SENSITIVITY_OUTPUT reads it"),
    ),
    (
        &["METRIC_ACCUMULATION_LEVEL", "LEVEL"],
        Action::Only("ALL_READS"),
    ),
];

struct Tool {
    picard_name: &'static str,
    subcommand: &'static str,
    options: &'static [Entry],
}

const TOOLS: &[Tool] = &[
    Tool {
        picard_name: "CollectQualityYieldMetrics",
        subcommand: "quality-yield",
        options: QUALITY_YIELD,
    },
    Tool {
        picard_name: "CollectHsMetrics",
        subcommand: "hs-metrics",
        options: HS_METRICS,
    },
];

/// A Picard command line rewritten for the native parser.
#[derive(Debug, PartialEq)]
pub struct Translation {
    /// Native subcommand name.
    pub subcommand: &'static str,
    /// Native flags, one `--flag=value` or `--flag` per element.
    pub args: Vec<String>,
    /// One line per dropped option, for standard error.
    pub notes: Vec<String>,
}

/// Rewrite `tool_name KEY=VALUE ...`. `Ok(None)` when `tool_name` is not a
/// Picard tool this crate implements.
pub fn translate(tool_name: &str, args: &[String]) -> Result<Option<Translation>> {
    let Some(tool) = TOOLS.iter().find(|t| t.picard_name == tool_name) else {
        return Ok(None);
    };
    let mut out = Translation {
        subcommand: tool.subcommand,
        args: Vec::new(),
        notes: Vec::new(),
    };
    for arg in args {
        let (name, value) = arg.split_once('=').ok_or_else(|| {
            Error::Usage(format!("{tool_name}: expected NAME=VALUE, got {arg:?}"))
        })?;
        let action = tool
            .options
            .iter()
            .chain(COMMON)
            .find(|(names, _)| names.contains(&name))
            .map(|(_, action)| action)
            .ok_or_else(|| {
                Error::Usage(format!(
                    "{tool_name}: option {name} is not implemented (it may change the metrics)"
                ))
            })?;
        match action {
            Action::Value(flag) => out.args.push(format!("{flag}={value}")),
            Action::Bool(flag) => {
                if parse_bool(name, value)? {
                    out.args.push((*flag).to_owned());
                }
            }
            Action::Only(allowed) => {
                if value != *allowed {
                    return Err(Error::Usage(format!(
                        "{tool_name}: {name}={value} is not implemented, only {allowed}"
                    )));
                }
            }
            Action::Ignore(why) => out.notes.push(format!("ignoring {name}={value}: {why}")),
        }
    }
    Ok(Some(out))
}

fn parse_bool(name: &str, value: &str) -> Result<bool> {
    match value.to_ascii_lowercase().as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(Error::Usage(format!(
            "{name}={value}: expected true or false"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn usage(result: Result<Option<Translation>>) -> String {
        match result {
            Err(Error::Usage(message)) => message,
            other => panic!("expected a usage error, got {other:?}"),
        }
    }

    #[test]
    fn accumulation_level_other_than_all_reads_is_an_error() {
        let ok = translate(
            "CollectHsMetrics",
            &strings(&["METRIC_ACCUMULATION_LEVEL=ALL_READS"]),
        );
        assert!(ok.is_ok());
        let message = usage(translate(
            "CollectHsMetrics",
            &strings(&["METRIC_ACCUMULATION_LEVEL=LIBRARY"]),
        ));
        assert!(message.contains("METRIC_ACCUMULATION_LEVEL"));
        assert!(message.contains("ALL_READS"));
    }

    #[test]
    fn unknown_option_names_itself() {
        let message = usage(translate(
            "CollectQualityYieldMetrics",
            &strings(&["STOP_AFTER=10"]),
        ));
        assert!(message.contains("STOP_AFTER"));
        let message = usage(translate("CollectQualityYieldMetrics", &strings(&["bam"])));
        assert!(message.contains("NAME=VALUE"));
    }

    #[test]
    fn registry_tools_are_not_translated() {
        assert_eq!(
            translate("CollectAlignmentSummaryMetrics", &[]).expect("no error"),
            None
        );
    }

    #[test]
    fn quality_yield_options() {
        let t = translate(
            "CollectQualityYieldMetrics",
            &strings(&[
                "I=x.bam",
                "O=out.txt",
                "OQ=false",
                "INCLUDE_SECONDARY_ALIGNMENTS=true",
                "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=false",
                "R=ref.fa",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(t.subcommand, "quality-yield");
        assert_eq!(
            t.args,
            [
                "--input=x.bam",
                "--output=out.txt",
                "--use-original-qualities=false",
                "--include-secondary-alignments",
            ]
        );
        assert_eq!(t.notes.len(), 1);
        let message = usage(translate(
            "CollectQualityYieldMetrics",
            &strings(&["FLOW_MODE=true"]),
        ));
        assert!(message.contains("FLOW_MODE"));
    }

    #[test]
    fn hs_metrics_options() {
        let t = translate(
            "CollectHsMetrics",
            &strings(&[
                "I=x.bam",
                "O=out.txt",
                "R=ref.fa",
                "BAIT_INTERVALS=b.interval_list",
                "TI=t1.interval_list",
                "TARGET_INTERVALS=t2.interval_list",
                "MQ=10",
                "Q=5",
                "CLIP_OVERLAPPING_READS=false",
                "INCLUDE_INDELS=true",
                "covMax=50",
                "ALLELE_FRACTION=0.1",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(t.subcommand, "hs-metrics");
        assert_eq!(
            t.args,
            [
                "--input=x.bam",
                "--output=out.txt",
                "--reference=ref.fa",
                "--bait-intervals=b.interval_list",
                "--target-intervals=t1.interval_list",
                "--target-intervals=t2.interval_list",
                "--minimum-mapping-quality=10",
                "--minimum-base-quality=5",
                "--clip-overlapping-reads=false",
                "--include-indels",
                "--coverage-cap=50",
            ]
        );
        assert_eq!(t.notes.len(), 1);
        let message = usage(translate(
            "CollectHsMetrics",
            &strings(&["PER_TARGET_COVERAGE=x.txt"]),
        ));
        assert!(message.contains("PER_TARGET_COVERAGE"));
    }

    #[test]
    fn clinical_coverage_command_translates() {
        // The pipeline's coverage step, short BI= next to long OUTPUT=.
        let t = translate(
            "CollectHsMetrics",
            &strings(&[
                "BI=panel.interval_list",
                "I=x.bam",
                "PER_BASE_COVERAGE=x.nucl.partial.out",
                "MINIMUM_MAPPING_QUALITY=0",
                "MINIMUM_BASE_QUALITY=0",
                "TARGET_INTERVALS=panel.interval_list",
                "OUTPUT=x.summary.partial.out",
                "R=ref.fa",
                "QUIET=true",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(
            t.args,
            [
                "--bait-intervals=panel.interval_list",
                "--input=x.bam",
                "--per-base-coverage=x.nucl.partial.out",
                "--minimum-mapping-quality=0",
                "--minimum-base-quality=0",
                "--target-intervals=panel.interval_list",
                "--output=x.summary.partial.out",
                "--reference=ref.fa",
            ]
        );
        assert_eq!(t.notes, ["ignoring QUIET=true: progress is not logged"]);
    }

    #[test]
    fn other_tools_are_not_ours() {
        assert_eq!(translate("insert-size", &[]).expect("no error"), None);
        assert_eq!(translate("MarkDuplicates", &[]).expect("no error"), None);
    }
}

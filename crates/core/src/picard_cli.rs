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

/// Options of `CollectInsertSizeMetrics` (`SinglePassSamProgram` and
/// `CollectInsertSizeMetrics` in Picard 3.5.0).
const INSERT_SIZE: &[Entry] = &[
    (&["INPUT", "I"], Action::Value("--input")),
    (&["OUTPUT", "O"], Action::Value("--output")),
    (
        &["HISTOGRAM_FILE", "H"],
        Action::Ignore("no histogram plot is drawn"),
    ),
    (&["DEVIATIONS"], Action::Value("--deviations")),
    (
        &["HISTOGRAM_WIDTH", "W"],
        Action::Value("--histogram-width"),
    ),
    (
        &["MINIMUM_HISTOGRAM_WIDTH", "MW"],
        Action::Value("--min-histogram-width"),
    ),
    (&["MINIMUM_PCT", "M"], Action::Value("--minimum-pct")),
    (
        &["INCLUDE_DUPLICATES"],
        Action::Bool("--include-duplicates"),
    ),
    (
        &["METRIC_ACCUMULATION_LEVEL", "LEVEL"],
        Action::Only("ALL_READS"),
    ),
    (
        &["ASSUME_SORTED", "AS"],
        Action::Ignore("one pass over the file in any order"),
    ),
];

/// Options of `CollectAlignmentSummaryMetrics` (`SinglePassSamProgram` and
/// `CollectAlignmentSummaryMetrics` in Picard 3.5.0). `R` is read here, so
/// this entry shadows the `COMMON` one. The two list options append to
/// Picard's default and `null` clears it; the native flags do the same.
const ALIGNMENT_SUMMARY: &[Entry] = &[
    (&["INPUT", "I"], Action::Value("--input")),
    (&["OUTPUT", "O"], Action::Value("--output")),
    (&["REFERENCE_SEQUENCE", "R"], Action::Value("--reference")),
    (
        &["HISTOGRAM_FILE", "H"],
        Action::Ignore("no histogram plot is drawn"),
    ),
    (&["MAX_INSERT_SIZE"], Action::Value("--max-insert-size")),
    (
        &["EXPECTED_PAIR_ORIENTATIONS"],
        Action::Value("--expected-pair-orientation"),
    ),
    (&["ADAPTER_SEQUENCE"], Action::Value("--adapter-sequence")),
    (
        &["COLLECT_ALIGNMENT_INFORMATION"],
        Action::Value("--collect-alignment-information"),
    ),
    (&["IS_BISULFITE_SEQUENCED", "BS"], Action::Only("false")),
    (
        &["METRIC_ACCUMULATION_LEVEL", "LEVEL"],
        Action::Only("ALL_READS"),
    ),
    (
        &["ASSUME_SORTED", "AS"],
        Action::Ignore("the header sort order is not checked"),
    ),
];

struct Tool {
    picard_name: &'static str,
    subcommand: &'static str,
    options: &'static [Entry],
}

const TOOLS: &[Tool] = &[
    Tool {
        picard_name: "CollectInsertSizeMetrics",
        subcommand: "insert-size",
        options: INSERT_SIZE,
    },
    Tool {
        picard_name: "CollectAlignmentSummaryMetrics",
        subcommand: "alignment-summary",
        options: ALIGNMENT_SUMMARY,
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
    fn maps_short_and_long_names() {
        let t = translate(
            "CollectInsertSizeMetrics",
            &strings(&[
                "I=x.bam",
                "O=out.txt",
                "H=hist.pdf",
                "M=0.5",
                "W=300",
                "MW=100",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(t.subcommand, "insert-size");
        assert_eq!(
            t.args,
            [
                "--input=x.bam",
                "--output=out.txt",
                "--minimum-pct=0.5",
                "--histogram-width=300",
                "--min-histogram-width=100"
            ]
        );
        assert_eq!(t.notes, ["ignoring H=hist.pdf: no histogram plot is drawn"]);
    }

    #[test]
    fn booleans_become_flags() {
        let on = translate(
            "CollectInsertSizeMetrics",
            &strings(&["INCLUDE_DUPLICATES=true"]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(on.args, ["--include-duplicates"]);
        let off = translate(
            "CollectInsertSizeMetrics",
            &strings(&["INCLUDE_DUPLICATES=FALSE"]),
        )
        .expect("translates")
        .expect("known tool");
        assert!(off.args.is_empty());
        let message = usage(translate(
            "CollectInsertSizeMetrics",
            &strings(&["INCLUDE_DUPLICATES=maybe"]),
        ));
        assert!(message.contains("INCLUDE_DUPLICATES"));
    }

    #[test]
    fn plumbing_options_are_dropped_with_a_note() {
        let t = translate(
            "CollectInsertSizeMetrics",
            &strings(&[
                "VALIDATION_STRINGENCY=SILENT",
                "TMP_DIR=/tmp",
                "QUIET=true",
                "R=ref.fa",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert!(t.args.is_empty());
        assert_eq!(t.notes.len(), 4);
        assert!(t.notes[0].starts_with("ignoring VALIDATION_STRINGENCY=SILENT"));
    }

    #[test]
    fn accumulation_level_other_than_all_reads_is_an_error() {
        let ok = translate(
            "CollectInsertSizeMetrics",
            &strings(&["METRIC_ACCUMULATION_LEVEL=ALL_READS"]),
        );
        assert!(ok.is_ok());
        let message = usage(translate(
            "CollectInsertSizeMetrics",
            &strings(&["METRIC_ACCUMULATION_LEVEL=LIBRARY"]),
        ));
        assert!(message.contains("METRIC_ACCUMULATION_LEVEL"));
        assert!(message.contains("ALL_READS"));
    }

    #[test]
    fn unknown_option_names_itself() {
        let message = usage(translate(
            "CollectInsertSizeMetrics",
            &strings(&["STOP_AFTER=10"]),
        ));
        assert!(message.contains("STOP_AFTER"));
        let message = usage(translate("CollectInsertSizeMetrics", &strings(&["bam"])));
        assert!(message.contains("NAME=VALUE"));
    }

    #[test]
    fn alignment_summary_reads_the_reference() {
        let t = translate(
            "CollectAlignmentSummaryMetrics",
            &strings(&[
                "I=x.bam",
                "O=out.txt",
                "R=ref.fa",
                "EXPECTED_PAIR_ORIENTATIONS=null",
                "EXPECTED_PAIR_ORIENTATIONS=RF",
                "MAX_INSERT_SIZE=500",
                "IS_BISULFITE_SEQUENCED=false",
            ]),
        )
        .expect("translates")
        .expect("known tool");
        assert_eq!(t.subcommand, "alignment-summary");
        assert_eq!(
            t.args,
            [
                "--input=x.bam",
                "--output=out.txt",
                "--reference=ref.fa",
                "--expected-pair-orientation=null",
                "--expected-pair-orientation=RF",
                "--max-insert-size=500",
            ]
        );
        assert!(t.notes.is_empty());
        let message = usage(translate(
            "CollectAlignmentSummaryMetrics",
            &strings(&["IS_BISULFITE_SEQUENCED=true"]),
        ));
        assert!(message.contains("IS_BISULFITE_SEQUENCED"));
    }

    #[test]
    fn other_tools_are_not_ours() {
        assert_eq!(translate("insert-size", &[]).expect("no error"), None);
        assert_eq!(translate("MarkDuplicates", &[]).expect("no error"), None);
    }
}

//! `CollectMultipleMetrics`, mirroring Picard 3.5.0: several
//! `SinglePassSamProgram` tools on one read of the BAM, each writing its
//! own file named `OUTPUT` plus a fixed extension (plus `FILE_EXTENSION`).
//!
//! Rules copied from `CollectMultipleMetrics.doWork`:
//! * `PROGRAM` is a list option: the default five stay unless `null` is
//!   given first, each name is kept once, in the order given.
//! * `METRIC_ACCUMULATION_LEVEL` goes to every program that has the option;
//!   `REFERENCE_SEQUENCE` to every program that reads one, and
//!   `CollectGcBiasMetrics` fails without it.
//! * `EXTRA_ARGUMENT=<PROGRAM>::<KEY=VALUE>` adds an option to one program.
//!   `INPUT`, `REFERENCE_SEQUENCE`, `ASSUME_SORTED` and `STOP_AFTER` in it
//!   are dropped without a word; the run's own values apply. An extra
//!   argument for a program not in `PROGRAM`, or for an unknown program, is
//!   an error, as is an option the program does not take.
//! * `STOP_AFTER` and `ASSUME_SORTED` apply to the shared read.
//! * `INTERVALS`, `DB_SNP`, `REF_FLAT`, `IGNORE_SEQUENCE` and
//!   `INCLUDE_UNPAIRED` serve only `RnaSeqMetrics` and
//!   `CollectSequencingArtifactMetrics`, which are not implemented: naming
//!   either program is an error, so the five options are accepted and
//!   dropped.
//! * Every program's `StringHeader` lines carry this tool's command line.

use std::ffi::OsString;
use std::path::PathBuf;

use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::pass::{self, Collector};
use crate::registry::{Invocation, Opt, Tool, picard_list};

/// One entry of Picard's `Program` enum: the tool it runs and the file
/// option plus extension of each file it writes.
struct Program {
    name: &'static str,
    tool: &'static Tool,
    files: &'static [(&'static str, &'static str)],
    levels: bool,
    reference: Reference,
}

/// Whether the program reads `REFERENCE_SEQUENCE`, and must have one.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Reference {
    Unused,
    Optional,
    Required,
}

const PROGRAMS: &[Program] = &[
    Program {
        name: "CollectAlignmentSummaryMetrics",
        tool: &crate::alignment_summary::TOOL,
        files: &[("OUTPUT", ".alignment_summary_metrics")],
        levels: true,
        reference: Reference::Optional,
    },
    Program {
        name: "CollectInsertSizeMetrics",
        tool: &crate::insert_size::TOOL,
        files: &[("OUTPUT", ".insert_size_metrics")],
        levels: true,
        reference: Reference::Unused,
    },
    Program {
        name: "QualityScoreDistribution",
        tool: &crate::quality_score_distribution::TOOL,
        files: &[("OUTPUT", ".quality_distribution_metrics")],
        levels: false,
        reference: Reference::Unused,
    },
    Program {
        name: "MeanQualityByCycle",
        tool: &crate::mean_quality_by_cycle::TOOL,
        files: &[("OUTPUT", ".quality_by_cycle_metrics")],
        levels: false,
        reference: Reference::Unused,
    },
    Program {
        name: "CollectBaseDistributionByCycle",
        tool: &crate::base_distribution_by_cycle::TOOL,
        files: &[("OUTPUT", ".base_distribution_by_cycle_metrics")],
        levels: false,
        reference: Reference::Unused,
    },
    Program {
        name: "CollectGcBiasMetrics",
        tool: &crate::gc_bias::TOOL,
        files: &[
            ("OUTPUT", ".gc_bias.detail_metrics"),
            ("SUMMARY_OUTPUT", ".gc_bias.summary_metrics"),
        ],
        levels: true,
        reference: Reference::Required,
    },
    Program {
        name: "CollectQualityYieldMetrics",
        tool: &crate::quality_yield::TOOL,
        files: &[("OUTPUT", ".quality_yield_metrics")],
        levels: false,
        reference: Reference::Unused,
    },
];

/// Picard's default `PROGRAM` list.
const DEFAULT_PROGRAMS: [&str; 5] = [
    "CollectAlignmentSummaryMetrics",
    "CollectBaseDistributionByCycle",
    "CollectInsertSizeMetrics",
    "MeanQualityByCycle",
    "QualityScoreDistribution",
];

/// Files of the tools `qc` runs that have no program here: the file option
/// and the extension after the prefix. `.rna_metrics` is the extension of
/// Picard's `RnaSeqMetrics` program; `.hs_metrics` is the one `collect`
/// writes; Picard has no prefix convention for the other two.
const OTHER_FILES: &[(&str, &[(&str, &str)])] = &[
    ("CollectHsMetrics", &[("OUTPUT", ".hs_metrics")]),
    (
        "CollectDuplicateMetrics",
        &[("METRICS_FILE", ".duplicate_metrics")],
    ),
    ("CollectWgsMetrics", &[("OUTPUT", ".wgs_metrics")]),
    ("CollectRnaSeqMetrics", &[("OUTPUT", ".rna_metrics")]),
    (
        "CollectCoverageGaps",
        &[
            ("OUTPUT", ".coverage_gaps"),
            ("GAPS_FILE", ".coverage_gaps.bed"),
        ],
    ),
    (
        "CollectSexCheckMetrics",
        &[("OUTPUT", ".sex_check_metrics")],
    ),
];

/// The files `tool` writes under an output prefix, as (file option,
/// extension) pairs: a program's own from the `CollectMultipleMetrics`
/// table, so `qc` and the Picard face write the same names, else from
/// [`OTHER_FILES`]. `None` for a tool with no prefix naming.
pub fn file_extensions(tool: &str) -> Option<&'static [(&'static str, &'static str)]> {
    PROGRAMS
        .iter()
        .find(|p| p.name == tool)
        .map(|p| p.files)
        .or_else(|| {
            OTHER_FILES
                .iter()
                .find(|(name, _)| *name == tool)
                .map(|(_, files)| *files)
        })
}

/// Programs Picard 3.5.0 lists that bammetrics does not have.
const NOT_IMPLEMENTED: [&str; 2] = ["RnaSeqMetrics", "CollectSequencingArtifactMetrics"];

/// Options `EXTRA_ARGUMENT` cannot set per program: Picard parses them
/// into the program and then reads the run's own values.
const RUN_LEVEL: [&str; 7] = [
    "INPUT",
    "I",
    "REFERENCE_SEQUENCE",
    "R",
    "ASSUME_SORTED",
    "AS",
    "STOP_AFTER",
];

pub const TOOL: Tool = Tool {
    name: "CollectMultipleMetrics",
    summary: "every PROGRAM on one read of the BAM, each written to OUTPUT plus its Picard \
              extension (plus FILE_EXTENSION); Python returns one row per file: file and rows.",
    options: &[
        Opt::output(&["OUTPUT", "O"])
            .with_help("The output prefix: each program writes PREFIX plus its extension.")
            .required(),
        crate::registry::REFERENCE_OPTION,
        crate::registry::STOP_AFTER_OPTION,
        crate::registry::LEVEL_OPTION,
        Opt::text(
            &["FILE_EXTENSION", "EXT"],
            "Appended to every file name after the program's extension.",
        ),
        Opt::text(
            &["PROGRAM"],
            "A program to run; each adds to Picard's default five, null clears.",
        )
        .repeated(),
        Opt::text(
            &["EXTRA_ARGUMENT"],
            "PROGRAM::KEY=VALUE, an option for one program.",
        )
        .repeated(),
        Opt::ignored(&["INTERVALS"], "only RnaSeqMetrics reads it"),
        Opt::ignored(
            &["DB_SNP"],
            "only CollectSequencingArtifactMetrics reads it",
        ),
        Opt::ignored(&["REF_FLAT"], "only RnaSeqMetrics reads it"),
        Opt::ignored(&["IGNORE_SEQUENCE"], "only RnaSeqMetrics reads it"),
        Opt::ignored_bool(
            &["INCLUDE_UNPAIRED", "UNPAIRED"],
            "only CollectSequencingArtifactMetrics reads it",
        ),
        Opt::ignored_bool(&["DO_NOT_CREATE_PLOTS"], "no plot is drawn"),
    ],
    run: run_tool,
    collector: None,
};

fn program(name: &str) -> Result<&'static Program> {
    if let Some(p) = PROGRAMS.iter().find(|p| p.name == name) {
        return Ok(p);
    }
    if NOT_IMPLEMENTED.contains(&name) {
        return Err(Error::Usage(format!(
            "PROGRAM {name} is not implemented in bammetrics"
        )));
    }
    Err(Error::Usage(format!(
        "Couldn't find program with value {name}"
    )))
}

/// `processAdditionalArguments`: `<PROGRAM>::<argument>`, where the
/// argument may carry a value after a space. The program is the text
/// before the last `::`. Picard 3.5.0 hands the argument to the program's
/// legacy parser, which takes `KEY=VALUE` only and rejects a dashed
/// `--KEY VALUE` as an invalid argument; so does this.
fn split_extra(text: &str) -> Result<(&'static Program, Vec<String>)> {
    let (name, rest) = text.rsplit_once("::").ok_or_else(|| {
        Error::Usage(format!(
            "couldn't understand EXTRA_ARGUMENT {text}: it doesn't conform to the form '<PROGRAM>::<ARGUMENT_AND_VALUE>'"
        ))
    })?;
    if rest.is_empty() {
        return Err(Error::Usage(format!(
            "couldn't understand EXTRA_ARGUMENT {text}: no argument after ::"
        )));
    }
    let program = program(name)?;
    if rest.starts_with('-') {
        return Err(Error::Usage(format!(
            "Invalid argument '{rest}' in EXTRA_ARGUMENT for {name}: use KEY=VALUE"
        )));
    }
    let mut parts = Vec::new();
    match rest.split_once(' ') {
        Some((argument, value)) if !value.trim_start_matches(' ').is_empty() => {
            parts.push(argument.to_owned());
            parts.push(value.trim_start_matches(' ').to_owned());
        }
        _ => parts.push(rest.to_owned()),
    }
    Ok((program, parts))
}

/// Drop the run-level options from one program's extra arguments.
fn without_run_level(args: Vec<String>) -> Vec<String> {
    args.into_iter()
        .filter(|arg| {
            let name = arg.split_once('=').map_or(arg.as_str(), |(n, _)| n);
            !RUN_LEVEL.contains(&name)
        })
        .collect()
}

/// The long name of the option `arg` (`NAME=VALUE`) sets, when the tool
/// takes it once. An extra argument replaces the value the program table
/// gives it (`CollectAlignmentSummaryMetrics::OUTPUT=asm.txt`).
fn single(tool: &Tool, arg: &str) -> Option<&'static str> {
    let (name, _) = arg.split_once('=')?;
    tool.option(name)
        .filter(|o| !o.repeated)
        .map(|o| o.picard())
}

/// `OUTPUT` plus `extension` plus `FILE_EXTENSION`.
fn file_name(base: &OsString, extension: &str, file_extension: &str) -> PathBuf {
    let mut name = base.clone();
    name.push(extension);
    name.push(file_extension);
    PathBuf::from(name)
}

/// The programs to run, each with the command line Picard would give it.
fn plan(inv: &Invocation) -> Result<Vec<(&'static Program, Invocation)>> {
    let names = picard_list(DEFAULT_PROGRAMS.to_vec(), inv.all("PROGRAM"), |s| {
        program(s).map(|p| p.name)
    })?;
    let mut programs: Vec<&'static Program> = Vec::new();
    for name in names {
        let p = program(name)?;
        if !programs.iter().any(|q| q.name == p.name) {
            programs.push(p);
        }
    }
    if programs.is_empty() {
        return Err(Error::Usage(
            "No programs specified with PROGRAM".to_owned(),
        ));
    }
    let mut extras: Vec<(&'static str, Vec<String>)> = Vec::new();
    for text in inv.all("EXTRA_ARGUMENT") {
        let (p, parts) = split_extra(text)?;
        match extras.iter_mut().find(|(name, _)| *name == p.name) {
            Some((_, list)) => list.extend(parts),
            None => extras.push((p.name, parts)),
        }
    }
    for (name, _) in &extras {
        if !programs.iter().any(|p| p.name == *name) {
            return Err(Error::Usage(format!(
                "EXTRA_ARGUMENT values were provided, but corresponding PROGRAM wasn't requested: {name}"
            )));
        }
    }
    let base = inv
        .output
        .as_ref()
        .ok_or_else(|| Error::Usage("CollectMultipleMetrics: OUTPUT is required".to_owned()))?
        .as_os_str()
        .to_owned();
    let file_extension = inv.get("FILE_EXTENSION").unwrap_or("");
    let reference = inv.get("REFERENCE_SEQUENCE");
    let levels = inv.all("METRIC_ACCUMULATION_LEVEL");
    let mut planned = Vec::new();
    for p in programs {
        if p.reference == Reference::Required && reference.is_none() {
            return Err(Error::Usage(format!(
                "The {} program needs a REF Sequence, please set REFERENCE_SEQUENCE in the command line",
                p.name
            )));
        }
        let mut args = vec![format!("INPUT={}", inv.input.display())];
        if let (Reference::Optional | Reference::Required, Some(r)) = (p.reference, reference) {
            args.push(format!("REFERENCE_SEQUENCE={r}"));
        }
        for (option, extension) in p.files {
            args.push(format!(
                "{option}={}",
                file_name(&base, extension, file_extension).display()
            ));
        }
        if p.levels {
            args.extend(
                levels
                    .iter()
                    .map(|l| format!("METRIC_ACCUMULATION_LEVEL={l}")),
            );
        } else if !levels.is_empty() {
            eprintln!(
                "bammetrics: the {} program does not support a metric accumulation level, but METRIC_ACCUMULATION_LEVEL was overridden in the command line; {} will be run against the entire input",
                p.name, p.name
            );
        }
        if let Some((_, extra)) = extras.iter().find(|(name, _)| *name == p.name) {
            let extra = without_run_level(extra.clone());
            let named: Vec<&str> = extra.iter().filter_map(|a| single(p.tool, a)).collect();
            args.retain(|a| single(p.tool, a).is_none_or(|n| !named.contains(&n)));
            args.extend(extra);
        }
        let mut program_inv = p
            .tool
            .parse(&args)
            .map_err(|e| Error::Usage(format!("Failed to parse arguments for {}: {e}", p.name)))?;
        program_inv.threads = inv.threads;
        planned.push((p, program_inv));
    }
    Ok(planned)
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let planned = plan(inv)?;
    let outputs = pass::run_many(
        &inv.input,
        inv.reference(),
        inv.threads,
        inv.limits()?,
        |header| {
            planned
                .iter()
                .map(|(p, program_inv)| {
                    let build = p
                        .tool
                        .collector
                        .ok_or_else(|| Error::Usage(format!("{} has no collector", p.name)))?;
                    build(program_inv, header)
                })
                .collect::<Result<Vec<Box<dyn Collector>>>>()
        },
    )?;
    let mut all = ToolOutput::NONE;
    for ((_, program_inv), output) in planned.iter().zip(outputs) {
        if let (Some(main), Some(path)) = (output.main, &program_inv.output) {
            all.extra.push((path.clone(), main));
        }
        all.extra.extend(output.extra);
    }
    Ok(all)
}

/// The programs of a parsed command line, by name, for the Python binding
/// and tests.
pub fn program_names(inv: &Invocation) -> Result<Vec<&'static str>> {
    Ok(plan(inv)?.into_iter().map(|(p, _)| p.name).collect())
}

/// Every `SinglePassSamProgram` tool `CollectMultipleMetrics` can run.
pub fn program_tools() -> impl Iterator<Item = &'static Tool> {
    PROGRAMS.iter().map(|p| p.tool)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn parse(args: &[&str]) -> Invocation {
        crate::registry::find("CollectMultipleMetrics")
            .expect("registered")
            .parse(&strings(args))
            .expect("parses")
    }

    #[test]
    fn program_is_a_list_with_the_default_five() {
        let inv = parse(&["I=x.bam", "O=out"]);
        assert_eq!(program_names(&inv).expect("plan"), DEFAULT_PROGRAMS);
        let inv = parse(&[
            "I=x.bam",
            "O=out",
            "R=ref.fa",
            "PROGRAM=CollectGcBiasMetrics",
            "PROGRAM=CollectQualityYieldMetrics",
            "PROGRAM=CollectGcBiasMetrics",
        ]);
        let names = program_names(&inv).expect("plan");
        assert_eq!(names.len(), 7);
        assert_eq!(names[5], "CollectGcBiasMetrics");
        let inv = parse(&[
            "I=x.bam",
            "O=out",
            "PROGRAM=null",
            "PROGRAM=MeanQualityByCycle",
        ]);
        assert_eq!(program_names(&inv).expect("plan"), ["MeanQualityByCycle"]);
        let inv = parse(&["I=x.bam", "O=out", "PROGRAM=null"]);
        assert!(plan(&inv).is_err());
        let inv = parse(&["I=x.bam", "O=out", "PROGRAM=RnaSeqMetrics"]);
        assert!(plan(&inv).is_err());
        let inv = parse(&["I=x.bam", "O=out", "PROGRAM=CollectGcBiasMetrics"]);
        assert!(plan(&inv).is_err());
    }

    #[test]
    fn files_are_named_by_extension() {
        let inv = parse(&[
            "I=x.bam",
            "O=dir/s1",
            "R=ref.fa",
            "EXT=.txt",
            "PROGRAM=CollectGcBiasMetrics",
            "LEVEL=LIBRARY",
        ]);
        let planned = plan(&inv).expect("plan");
        let (gc, gc_inv) = planned.last().expect("gc bias");
        assert_eq!(gc.name, "CollectGcBiasMetrics");
        assert_eq!(
            gc_inv.output.as_deref(),
            Some(std::path::Path::new("dir/s1.gc_bias.detail_metrics.txt"))
        );
        assert_eq!(
            gc_inv.get("SUMMARY_OUTPUT"),
            Some("dir/s1.gc_bias.summary_metrics.txt")
        );
        assert_eq!(gc_inv.get("REFERENCE_SEQUENCE"), Some("ref.fa"));
        assert_eq!(gc_inv.all("METRIC_ACCUMULATION_LEVEL"), ["LIBRARY"]);
        let (_, qsd) = &planned[4];
        assert!(qsd.get("METRIC_ACCUMULATION_LEVEL").is_none());
        assert!(qsd.get("REFERENCE_SEQUENCE").is_none());
        assert_eq!(
            qsd.output.as_deref(),
            Some(std::path::Path::new(
                "dir/s1.quality_distribution_metrics.txt"
            ))
        );
    }

    #[test]
    fn extra_arguments_reach_one_program() {
        let inv = parse(&[
            "I=x.bam",
            "O=out",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::MINIMUM_PCT=0.1",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::INPUT=other.bam",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::STOP_AFTER=3",
            "EXTRA_ARGUMENT=CollectAlignmentSummaryMetrics::OUTPUT=asm.txt",
        ]);
        let planned = plan(&inv).expect("plan");
        let (_, insert) = planned
            .iter()
            .find(|(p, _)| p.name == "CollectInsertSizeMetrics")
            .expect("insert size");
        assert_eq!(insert.get("HISTOGRAM_WIDTH"), Some("200"));
        assert_eq!(insert.get("MINIMUM_PCT"), Some("0.1"));
        assert_eq!(insert.input, std::path::Path::new("x.bam"));
        assert!(insert.get("STOP_AFTER").is_none());
        let (_, asm) = &planned[0];
        assert_eq!(asm.output.as_deref(), Some(std::path::Path::new("asm.txt")));
        for bad in [
            "EXTRA_ARGUMENT=CollectGcBiasMetrics::SCAN_WINDOW_SIZE=50",
            "EXTRA_ARGUMENT=Bogus::X=1",
            "EXTRA_ARGUMENT=nothing",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::BOGUS=1",
            "EXTRA_ARGUMENT=CollectInsertSizeMetrics::--HISTOGRAM_WIDTH 200",
        ] {
            let inv = parse(&["I=x.bam", "O=out", bad]);
            assert!(plan(&inv).is_err(), "{bad}");
        }
    }

    /// The name table pinned: a change here renames files on both faces
    /// (`CollectMultipleMetrics O=PREFIX`, `collect`, `qc` and the verbs).
    #[test]
    fn both_faces_name_files_from_one_table() {
        let ext = |tool| file_extensions(tool).expect("named")[0].1;
        for (tool, extension) in [
            (
                "CollectAlignmentSummaryMetrics",
                ".alignment_summary_metrics",
            ),
            ("CollectInsertSizeMetrics", ".insert_size_metrics"),
            ("QualityScoreDistribution", ".quality_distribution_metrics"),
            ("MeanQualityByCycle", ".quality_by_cycle_metrics"),
            (
                "CollectBaseDistributionByCycle",
                ".base_distribution_by_cycle_metrics",
            ),
            ("CollectQualityYieldMetrics", ".quality_yield_metrics"),
            ("CollectHsMetrics", ".hs_metrics"),
            ("CollectDuplicateMetrics", ".duplicate_metrics"),
            ("CollectWgsMetrics", ".wgs_metrics"),
            ("CollectRnaSeqMetrics", ".rna_metrics"),
        ] {
            assert_eq!(ext(tool), extension, "{tool}");
            assert_eq!(crate::collect::extension(tool), extension, "{tool}");
        }
        assert_eq!(
            file_extensions("CollectGcBiasMetrics").expect("named"),
            [
                ("OUTPUT", ".gc_bias.detail_metrics"),
                ("SUMMARY_OUTPUT", ".gc_bias.summary_metrics")
            ]
        );
        assert_eq!(
            file_extensions("CollectDuplicateMetrics").expect("named"),
            [("METRICS_FILE", ".duplicate_metrics")]
        );
        assert!(file_extensions("EstimateLibraryComplexity").is_none());
        for (tool, files) in OTHER_FILES {
            let registered = crate::registry::find(tool).expect("registered");
            assert_eq!(files[0].0, registered.output_name(), "{tool}");
            assert!(PROGRAMS.iter().all(|p| p.name != *tool), "{tool}");
        }
    }

    #[test]
    fn extra_split_keeps_the_last_double_colon() {
        let (p, parts) = split_extra("CollectInsertSizeMetrics::W=200 x").expect("splits");
        assert_eq!(p.name, "CollectInsertSizeMetrics");
        assert_eq!(parts, ["W=200", "x"]);
        let (_, parts) = split_extra("CollectInsertSizeMetrics::W=200 ").expect("splits");
        assert_eq!(parts, ["W=200 "]);
        assert!(split_extra("CollectInsertSizeMetrics::").is_err());
        assert!(split_extra("CollectInsertSizeMetrics::--W 200").is_err());
    }
}

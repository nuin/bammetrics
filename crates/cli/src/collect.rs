//! `bammetrics collect`: several metrics on one read of the BAM, one
//! Picard file each under `-o PREFIX` (default `./<bam stem>`). The flags
//! and defaults of the 0.12 subcommand, read by the verb grammar.

use std::path::{Path, PathBuf};

use bammetrics_core::Error;
use bammetrics_core::alignment_summary;
use bammetrics_core::collect::{self, Outputs, Request};
use bammetrics_core::hs_metrics;
use bammetrics_core::insert_size;
use bammetrics_core::pass;
use bammetrics_core::quality_yield;
use bammetrics_core::registry::{Ty, did_you_mean};

use crate::error::VerbError;
use crate::grammar::{self, Flag, Key, Parsed};
use crate::verbs::{option_lines, shared_flags, stem};

const VERB: &str = "collect";

/// A metric `collect` computes, by its `--metric` name.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Metric {
    QualityYield,
    AlignmentSummary,
    InsertSize,
    Hs,
}

const METRICS: [(&str, Metric); 4] = [
    ("quality-yield", Metric::QualityYield),
    ("alignment-summary", Metric::AlignmentSummary),
    ("insert-size", Metric::InsertSize),
    ("hs", Metric::Hs),
];

pub fn flags() -> Vec<Flag> {
    let mut hs = Flag::new(
        Key::Own("hs"),
        "hs",
        &[],
        Ty::Path,
        "One CollectHsMetrics run: a bait and a target interval file (interval list or BED).",
    )
    .placeholder("BAITS TARGETS")
    .repeated();
    hs.values = 2;
    let mut flags = vec![
        Flag::new(Key::Input, "input", &["i", "I"], Ty::Path, "The BAM to read; - reads standard input.")
            .names(&["INPUT", "I"])
            .placeholder("BAM"),
        Flag::new(
            Key::Output,
            "output",
            &["o", "O"],
            Ty::Path,
            "Output prefix (default ./ plus the BAM name without .bam). Files are PREFIX plus \
             .quality_yield_metrics, .alignment_summary_metrics, .insert_size_metrics (the names \
             of CollectMultipleMetrics) and .hs_metrics, or .BAIT_SET.hs_metrics with several --hs.",
        )
        .names(&["OUTPUT", "O"])
        .placeholder("PREFIX"),
        Flag::new(
            Key::Own("metric"),
            "metric",
            &["m"],
            Ty::Text,
            "Metrics, comma separated: quality-yield, alignment-summary, insert-size, hs \
             (default: the first three, and hs with --hs).",
        )
        .names(&["METRICS"])
        .placeholder("NAMES")
        .repeated(),
        Flag::new(
            Key::Own("reference"),
            "reference",
            &["R"],
            Ty::Path,
            "Reference FASTA with a .fai index, for alignment-summary and hs.",
        )
        .names(&["REFERENCE_SEQUENCE"])
        .placeholder("FASTA"),
        hs,
        Flag::new(
            Key::Own("hs-per-base"),
            "hs-per-base",
            &[],
            Ty::Bool,
            "Also write each --hs run's per-base coverage, named as its HS file with \
             .per_base_coverage in place of .hs_metrics.",
        ),
    ];
    flags.extend(shared_flags());
    flags
}

pub fn help() -> String {
    format!(
        "bammetrics collect: several metrics from one read of the BAM, one Picard file each.\n\n\
         Usage: bammetrics collect [BAM] [OPTIONS]\n\nOptions:\n{}",
        option_lines(&flags())
    )
}

/// The metrics `--metric` names, Picard's defaults for each, `-R` given
/// to the two that read a reference.
fn request(parsed: &Parsed) -> Result<Request, VerbError> {
    let mut wanted = Vec::new();
    for given in parsed.all(Key::Own("metric")) {
        for name in given.values[0].split(',') {
            let metric = METRICS
                .iter()
                .find(|(n, _)| *n == name)
                .map(|(_, m)| *m)
                .ok_or_else(|| {
                    let mut error = VerbError::usage(VERB, format!("no metric {name:?}"))
                        .from(given.spelled.clone());
                    if let Some(best) = did_you_mean(name, METRICS.iter().map(|(n, _)| *n)) {
                        error = error.hint(format!("did you mean {best}?"));
                    }
                    error.hint("metrics: quality-yield, alignment-summary, insert-size, hs")
                })?;
            wanted.push(metric);
        }
    }
    let hs: Vec<&Vec<String>> = parsed.all(Key::Own("hs")).map(|g| &g.values).collect();
    if wanted.is_empty() {
        wanted = vec![
            Metric::QualityYield,
            Metric::AlignmentSummary,
            Metric::InsertSize,
        ];
        if !hs.is_empty() {
            wanted.push(Metric::Hs);
        }
    }
    if wanted.contains(&Metric::Hs) == hs.is_empty() {
        return Err(VerbError::usage(
            VERB,
            "--metric hs and --hs BAITS TARGETS go together",
        ));
    }
    let reference = parsed
        .last(Key::Own("reference"))
        .map(|g| PathBuf::from(&g.values[0]));
    let hs_metrics = hs
        .iter()
        .map(|pair| hs_metrics::Options {
            bait_intervals: vec![PathBuf::from(&pair[0])],
            target_intervals: vec![PathBuf::from(&pair[1])],
            reference: reference.clone(),
            ..hs_metrics::Options::default()
        })
        .collect();
    Ok(Request {
        quality_yield: wanted
            .contains(&Metric::QualityYield)
            .then(quality_yield::Options::default),
        alignment_summary: wanted.contains(&Metric::AlignmentSummary).then(|| {
            alignment_summary::Options {
                reference: reference.clone(),
                ..alignment_summary::Options::default()
            }
        }),
        insert_size: wanted
            .contains(&Metric::InsertSize)
            .then(insert_size::Options::default),
        hs_metrics,
    })
}

pub fn run(args: &[String]) -> Result<u8, VerbError> {
    let parsed = grammar::parse(VERB, &flags(), args)?;
    if parsed.flag(Key::Help) {
        print!("{}", help());
        return Ok(0);
    }
    let inputs: Vec<&String> = parsed
        .positional
        .iter()
        .chain(parsed.all(Key::Input).flat_map(|g| &g.values))
        .collect();
    let bam = match inputs.as_slice() {
        [bam] => PathBuf::from(bam.as_str()),
        [] => {
            return Err(
                VerbError::usage(VERB, "missing BAM").hint("bammetrics collect BAM -o PREFIX")
            );
        }
        [_, second, ..] => {
            return Err(
                VerbError::usage(VERB, format!("one BAM expected, got {}", inputs.len()))
                    .from(second.as_str()),
            );
        }
    };
    let prefix = match parsed.last(Key::Output) {
        Some(given) => {
            let path = PathBuf::from(&given.values[0]);
            if given.values[0].ends_with('/') || path.is_dir() {
                path.join(stem(&bam.to_string_lossy()))
            } else {
                path
            }
        }
        None => {
            let path = Path::new(".").join(stem(&bam.to_string_lossy()));
            if !parsed.flag(Key::Quiet) {
                eprintln!(
                    "bammetrics collect: no -o, files go to {}.*",
                    path.display()
                );
            }
            path
        }
    };
    let threads = match parsed.last(Key::Threads) {
        None => pass::default_threads(),
        Some(given) => given.values[0].parse().map_err(|_| {
            VerbError::usage(VERB, format!("{:?} is not a number", given.values[0]))
                .from(given.spelled.clone())
        })?,
    };
    let source = if bam == Path::new("-") {
        PathBuf::from("/dev/stdin")
    } else if bam.exists() {
        bam.clone()
    } else {
        return Err(VerbError::input(VERB, format!("{} not found", bam.display())).from("BAM"));
    };
    let mut request = request(&parsed)?;
    let core = |e: Error| VerbError::core(VERB, &e);
    if let Some(parent) = prefix.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)
            .map_err(|e| VerbError::input(VERB, format!("{}: {e}", parent.display())))?;
    }
    let outputs = Outputs::new(&prefix, &request).map_err(core)?;
    if parsed.flag(Key::Own("hs-per-base")) {
        for (opts, path) in request.hs_metrics.iter_mut().zip(&outputs.hs_metrics) {
            opts.per_base_coverage = Some(collect::per_base_coverage_path(path));
        }
    }
    let results = collect::run_over_bam(&source, &request, threads).map_err(core)?;
    collect::write_files(&request, &results, &bam, &outputs).map_err(core)?;
    if let (Some(opts), Some(Ok(rows))) = (&request.insert_size, &results.insert_size)
        && rows.is_empty()
    {
        eprintln!("{}", insert_size::no_rows_warning(opts.minimum_pct));
    }
    let errors = results.errors();
    for (tool, e) in &errors {
        eprintln!("bammetrics: {tool}: {e}");
    }
    if !errors.is_empty() {
        return Err(VerbError::input(
            VERB,
            format!("{} metric(s) failed", errors.len()),
        ));
    }
    Ok(0)
}

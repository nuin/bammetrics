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
use crate::verbs::{option_lines, refuse_overwriting_inputs, shared_flags, stem};

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
        Flag::new(Key::Input, "input", &["i", "I"], Ty::Path, "The BAM or CRAM to read; - reads standard input.")
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
    // 0.12 clap: `--hs-per-base` requires `--hs`.
    if let Some(given) = parsed.last(Key::Own("hs-per-base"))
        && hs.is_empty()
    {
        return Err(
            VerbError::usage(VERB, "--hs-per-base needs --hs BAITS TARGETS")
                .from(given.spelled.clone()),
        );
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
        reference,
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
            let mut error = VerbError::usage(VERB, "missing BAM");
            if let Some(hint) = parsed.bool_word_hint() {
                error = error.hint(hint);
            }
            return Err(error.hint("bammetrics collect BAM -o PREFIX"));
        }
        [_, second, ..] => {
            let mut error =
                VerbError::usage(VERB, format!("one BAM expected, got {}", inputs.len()))
                    .from(second.as_str());
            if let Some(hint) = (0..parsed.positional.len().min(2))
                .rev()
                .find_map(|i| parsed.bool_value_hint(i))
            {
                error = error.hint(hint);
            }
            return Err(error);
        }
    };
    let prefix = match parsed.last(Key::Output) {
        Some(given) if given.values[0] == "-" => {
            return Err(VerbError::usage(
                VERB,
                "standard output takes one file; this verb writes several",
            )
            .from(given.spelled.clone())
            .hint("give -o PREFIX"));
        }
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
        Some(given) => given.values[0]
            .parse()
            .map_err(|_| format!("{:?} is not a number", given.values[0]))
            .and_then(|threads: usize| {
                pass::check_threads(threads).map_err(|e| format!("{threads} is {e}"))
            })
            .map_err(|what| VerbError::usage(VERB, what).from(given.spelled.clone()))?,
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
    let outputs = Outputs::new(&prefix, &request).map_err(core)?;
    let mut read = vec![bam.clone()];
    read.extend(
        parsed
            .last(Key::Own("reference"))
            .map(|g| PathBuf::from(&g.values[0])),
    );
    read.extend(
        parsed
            .all(Key::Own("hs"))
            .flat_map(|g| g.values.iter().map(PathBuf::from)),
    );
    let written: Vec<PathBuf> = [
        (request.quality_yield.is_some(), &outputs.quality_yield),
        (
            request.alignment_summary.is_some(),
            &outputs.alignment_summary,
        ),
        (request.insert_size.is_some(), &outputs.insert_size),
    ]
    .into_iter()
    .filter(|(wanted, _)| *wanted)
    .map(|(_, path)| path.clone())
    .chain(outputs.hs_metrics.iter().cloned())
    .collect();
    refuse_overwriting_inputs(VERB, &read, &written).map_err(|e| {
        match parsed.last(Key::Output) {
            Some(given) => e.from(given.spelled.clone()),
            None => e,
        }
    })?;
    // docs/05 section 2.3: missing output directories are created.
    if let Some(parent) = prefix.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)
            .map_err(|e| VerbError::input(VERB, format!("{}: {e}", parent.display())))?;
    }
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

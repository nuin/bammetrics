//! `qc --also FORMAT[,FORMAT]` (docs/formats.md): further formats beside
//! the Picard files of a run. `json`, `jsonl`, `tsv` and `csv` are written
//! beside every metrics file as `FILE.json`; `markdown`, `junit` and
//! `multiqc` are one file per run under the prefix. [`targets`] names the
//! files (for `plan` too), [`texts`] renders them from the metrics files
//! on disk, and `engine::write_also` writes them and records the manifest
//! step `also`.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_formats::{Document, Format, junit, markdown, multiqc};
use bammetrics_verdict::Verdict;
use clap::ValueEnum;

use crate::check::FailOn;
use crate::error::{Error, Result};
use crate::format::{sibling, summary};
use crate::qc::{manifest_path, verdict_path, with_suffix};

/// A format `qc --also` writes beside the Picard files.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Also {
    /// FILE.json beside every metrics file (schema bammetrics.metrics/1).
    Json,
    /// FILE.jsonl beside every metrics file: one object per row.
    Jsonl,
    /// FILE.tsv beside every metrics file: the metrics table.
    Tsv,
    /// FILE.csv beside every metrics file: the metrics table, RFC 4180.
    Csv,
    /// PREFIX.report.md: the verdict, then every metrics table, at most
    /// 65,536 characters.
    #[value(alias = "md")]
    Markdown,
    /// PREFIX.junit.xml: a testcase per rule.
    Junit,
    /// PREFIX.verdict_mqc.json: MultiQC custom content.
    Multiqc,
}

impl Also {
    /// The name `--also` takes, and the manifest kind of the files.
    pub fn name(self) -> &'static str {
        match self {
            Also::Json => "json",
            Also::Jsonl => "jsonl",
            Also::Tsv => "tsv",
            Also::Csv => "csv",
            Also::Markdown => "markdown",
            Also::Junit => "junit",
            Also::Multiqc => "multiqc",
        }
    }

    /// The format written beside each metrics file; `None` for a format
    /// that is one file per run.
    fn per_file(self) -> Option<Format> {
        match self {
            Also::Json => Some(Format::Json),
            Also::Jsonl => Some(Format::Jsonl),
            Also::Tsv => Some(Format::Tsv),
            Also::Csv => Some(Format::Csv),
            Also::Markdown | Also::Junit | Also::Multiqc => None,
        }
    }

    /// The file of a one-file-per-run format: the prefix plus this.
    fn run_suffix(self) -> &'static str {
        match self {
            Also::Markdown => ".report.md",
            Also::Junit => ".junit.xml",
            Also::Multiqc => ".verdict_mqc.json",
            Also::Json | Also::Jsonl | Also::Tsv | Also::Csv => "",
        }
    }
}

/// True when `path` is a file the format named `kind` writes for
/// `prefix`: the prefix plus the format's suffix, or a path with the
/// format's extension that is not the run's verdict or manifest. An
/// earlier manifest is a file on disk, so before `qc` removes a stale
/// `--also` file it checks the path against the names it computes itself.
pub fn is_target(prefix: &Path, path: &Path, kind: &str) -> bool {
    let Some(also) = Also::value_variants().iter().find(|a| a.name() == kind) else {
        return false;
    };
    match also.per_file() {
        Some(format) => {
            path.extension().and_then(|e| e.to_str()) == Some(format.extension())
                && path != verdict_path(prefix)
                && path != manifest_path(prefix)
        }
        None => path == with_suffix(prefix, also.run_suffix()),
    }
}

/// One file `--also` writes.
pub struct Target {
    pub path: PathBuf,
    pub also: Also,
    /// The metrics file a per-file format is written beside, as an index
    /// into the metrics paths; `None` for a one-file-per-run format.
    source: Option<usize>,
}

/// The files the formats of `also` write for the metrics files `metrics`
/// under `prefix`, in write order: format by format as given, a format
/// given twice once.
pub fn targets(prefix: &Path, metrics: &[&Path], also: &[Also]) -> Vec<Target> {
    let mut targets = Vec::new();
    for (i, format) in also.iter().enumerate() {
        if also[..i].contains(format) {
            continue;
        }
        match format.per_file() {
            Some(per_file) => targets.extend(metrics.iter().enumerate().map(|(j, path)| Target {
                path: sibling(path, per_file.extension()),
                also: *format,
                source: Some(j),
            })),
            None => targets.push(Target {
                path: with_suffix(prefix, format.run_suffix()),
                also: *format,
                source: None,
            }),
        }
    }
    targets
}

/// Each target with its text. The metrics files are read back from disk,
/// so a file `--resume` kept is exported like one this run wrote; the JSON
/// forms carry `sample`. `warn_fails` is the gate `--fail-on warn`, under
/// which JUnit counts a warning as a failure.
pub fn texts(
    prefix: &Path,
    metrics: &[&Path],
    also: &[Also],
    sample: &str,
    verdict: &Verdict,
    fail_on: Option<FailOn>,
) -> Result<Vec<(Target, String)>> {
    let targets = targets(prefix, metrics, also);
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    let mut docs: Vec<Document> = Vec::with_capacity(metrics.len());
    for path in metrics {
        let text = fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
        let mut doc = bammetrics_formats::read(&text)
            .map_err(|e| Error::Input(format!("{}: {e}", path.display())))?;
        doc.sample = Some(sample.to_owned());
        docs.push(doc);
    }
    let verdicts = [summary(verdict)];
    Ok(targets
        .into_iter()
        .map(|target| {
            let text = match (target.also.per_file(), target.source, target.also) {
                (Some(format), Some(j), _) => bammetrics_formats::render(format, &docs[j], None),
                (_, _, Also::Junit) => junit::write(&verdicts, fail_on == Some(FailOn::Warn)),
                (_, _, Also::Multiqc) => multiqc::write(&verdicts),
                _ => markdown::write(&docs.iter().collect::<Vec<_>>(), &verdicts),
            };
            (target, text)
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn targets_follow_the_formats_as_given() {
        let metrics = [
            Path::new("out/s.insert_size_metrics"),
            Path::new("out/s.hs_metrics"),
        ];
        let also = [Also::Markdown, Also::Csv, Also::Markdown, Also::Multiqc];
        let paths: Vec<String> = targets(Path::new("out/s"), &metrics, &also)
            .iter()
            .map(|t| format!("{} {}", t.also.name(), t.path.display()))
            .collect();
        assert_eq!(
            paths,
            [
                "markdown out/s.report.md",
                "csv out/s.insert_size_metrics.csv",
                "csv out/s.hs_metrics.csv",
                "multiqc out/s.verdict_mqc.json",
            ]
        );
        assert!(targets(Path::new("out/s"), &metrics, &[]).is_empty());
    }

    #[test]
    fn a_stale_file_is_removed_only_under_a_name_also_writes() {
        let prefix = Path::new("out/s");
        let metrics = [Path::new("out/s.insert_size_metrics")];
        for target in targets(prefix, &metrics, Also::value_variants()) {
            assert!(
                is_target(prefix, &target.path, target.also.name()),
                "{}",
                target.path.display()
            );
        }
        // The BAM, a metrics file and another run's page are never one.
        for (path, kind) in [
            ("out/s.bam", "json"),
            ("out/s.insert_size_metrics", "csv"),
            ("out/s.insert_size_metrics.json", "metrics"),
            ("out/t.junit.xml", "junit"),
            ("out/s.report.html", "markdown"),
            ("out/s.verdict.json", "multiqc"),
            ("out/s.verdict.json", "json"),
            ("out/s.manifest.json", "json"),
        ] {
            assert!(
                !is_target(prefix, Path::new(path), kind),
                "{path} as {kind}"
            );
        }
    }
}

//! `qc --also FORMAT[,FORMAT]` (docs/formats.md): further formats beside
//! the Picard files of a run. `json`, `jsonl`, `tsv`, `csv`, `parquet`,
//! `xlsx`, `docx` and `sql` are written beside every metrics file as
//! `FILE.json` (`parquet` with its companions); `markdown`, `junit` and
//! `multiqc` are one file per run under the prefix. [`targets`] names the
//! files (for `plan` too), [`texts`] renders them from the metrics files
//! on disk, and `engine::write_also` writes them and records the manifest
//! step `also`.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_formats::{Document, File, Format, junit, markdown, multiqc, parquet};
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
    /// FILE.parquet beside every metrics file, with FILE.histogram.parquet
    /// (when the file has a histogram) and FILE.long.parquet.
    Parquet,
    /// FILE.xlsx beside every metrics file: a sheet per table.
    Xlsx,
    /// FILE.docx beside every metrics file: the headers and every table.
    Docx,
    /// FILE.sql beside every metrics file: an SQL dump sqlite3 loads.
    Sql,
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
            Also::Parquet => "parquet",
            Also::Xlsx => "xlsx",
            Also::Docx => "docx",
            Also::Sql => "sql",
            Also::Markdown => "markdown",
            Also::Junit => "junit",
            Also::Multiqc => "multiqc",
        }
    }

    /// The format written beside each metrics file; `None` for a format
    /// that is one file per run.
    pub(crate) fn per_file(self) -> Option<Format> {
        match self {
            Also::Json => Some(Format::Json),
            Also::Jsonl => Some(Format::Jsonl),
            Also::Tsv => Some(Format::Tsv),
            Also::Csv => Some(Format::Csv),
            Also::Parquet => Some(Format::Parquet),
            Also::Xlsx => Some(Format::Xlsx),
            Also::Docx => Some(Format::Docx),
            Also::Sql => Some(Format::Sql),
            Also::Markdown | Also::Junit | Also::Multiqc => None,
        }
    }

    /// The file of a one-file-per-run format: the prefix plus this.
    pub(crate) fn run_suffix(self) -> &'static str {
        match self {
            Also::Markdown => ".report.md",
            Also::Junit => ".junit.xml",
            Also::Multiqc => ".verdict_mqc.json",
            Also::Json
            | Also::Jsonl
            | Also::Tsv
            | Also::Csv
            | Also::Parquet
            | Also::Xlsx
            | Also::Docx
            | Also::Sql => "",
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
    /// The companion this file is (`parquet`'s `histogram` and `long`);
    /// `None` for a format's main file.
    part: Option<&'static str>,
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
            Some(per_file) => {
                for (j, path) in metrics.iter().enumerate() {
                    let main = sibling(path, per_file.extension());
                    let parts = std::iter::once(None)
                        .chain(per_file.companions().iter().map(|part| Some(*part)));
                    targets.extend(parts.map(|part| Target {
                        path:
                            part.map_or_else(
                                || main.clone(),
                                |part| parquet::companion(&main, part),
                            ),
                        also: *format,
                        source: Some(j),
                        part,
                    }));
                }
            }
            None => targets.push(Target {
                path: with_suffix(prefix, format.run_suffix()),
                also: *format,
                source: None,
                part: None,
            }),
        }
    }
    targets
}

/// The metrics files of the run as documents, each with `sample`. A file
/// with headers only is a document with no table: `CollectDuplicateMetrics`
/// writes one for a BAM with no read, and the run that wrote it exits 0.
fn documents(metrics: &[&Path], sample: &str) -> Result<Vec<Document>> {
    metrics
        .iter()
        .map(|path| {
            let text = fs::read_to_string(path).map_err(|e| Error::io(path, e))?;
            let mut doc = bammetrics_formats::picard::read_or_empty(&text)
                .map_err(|e| Error::Input(format!("{}: {e}", path.display())))?;
            doc.sample = Some(sample.to_owned());
            Ok(doc)
        })
        .collect()
}

/// Each target with its bytes. The metrics files are read back from disk
/// when a format prints them (`junit` and `multiqc` print the verdict
/// only), so a file `--resume` kept is exported like one this run wrote;
/// the JSON and Parquet forms carry `sample`. A metrics file with no table
/// still gets its sibling: the JSON document with the headers, an empty
/// `tsv`, `csv` or `jsonl` file, a Parquet file of part `empty`. A
/// histogram companion is left out for a metrics file with no histogram.
/// `warn_fails` is the gate `--fail-on warn`, under which JUnit counts a
/// warning as a failure.
pub fn texts(
    prefix: &Path,
    metrics: &[&Path],
    also: &[Also],
    sample: &str,
    verdict: &Verdict,
    fail_on: Option<FailOn>,
) -> Result<Vec<(Target, Vec<u8>)>> {
    let targets = targets(prefix, metrics, also);
    let reads_metrics = |t: &Target| t.source.is_some() || t.also == Also::Markdown;
    let docs = if targets.iter().any(reads_metrics) {
        documents(metrics, sample)?
    } else {
        Vec::new()
    };
    let verdicts = [summary(verdict)];
    // The files of one format for one metrics file, rendered once for the
    // main target and its companions, which follow it.
    let mut rendered: Option<((Format, usize), Vec<File>)> = None;
    let mut out = Vec::with_capacity(targets.len());
    for target in targets {
        let bytes = match (target.also.per_file(), target.source, target.also) {
            (Some(format), Some(j), _) => {
                if rendered.as_ref().is_none_or(|(key, _)| *key != (format, j)) {
                    let files = bammetrics_formats::render_files(format, &docs[j], None)
                        .map_err(|e| Error::Input(format!("{}: {e}", metrics[j].display())))?;
                    rendered = Some(((format, j), files));
                }
                let files = rendered.as_mut().map(|(_, files)| files);
                let file = files.and_then(|files| {
                    let at = files.iter().position(|f| f.part == target.part)?;
                    Some(files.swap_remove(at))
                });
                match file {
                    Some(file) => file.bytes,
                    None => continue,
                }
            }
            (_, _, Also::Junit) => junit::write(&verdicts, fail_on == Some(FailOn::Warn)).into(),
            (_, _, Also::Multiqc) => multiqc::write(&verdicts).into(),
            _ => markdown::write(&docs.iter().collect::<Vec<_>>(), &verdicts).into(),
        };
        out.push((target, bytes));
    }
    Ok(out)
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

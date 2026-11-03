//! Metrics file discovery for `report` and `check`: every file under the
//! directories given that holds a `## METRICS CLASS` block (or a
//! `## HISTOGRAM` section under a Picard header, as the by-cycle tools
//! write), at any depth and under any name, grouped into samples by
//! [`SAMPLE_RULE`]. No BAM is read.

use std::collections::BTreeMap;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_verdict::Input;

use crate::error::{Error, Result};

/// How files become samples, as `--help` prints it.
pub const SAMPLE_RULE: &str = "Samples: the metrics files in one directory are grouped by the file \
name up to the first dot, so NA12878.insert_size_metrics and NA12878.hs_metrics form the sample \
NA12878 (the prefix Picard CollectMultipleMetrics writes). When no metrics class repeats across \
the groups of a directory, the directory holds one sample, named by the prefix when two or more \
files share it, else by the directory: runs/NA12878/insert_size.txt and runs/NA12878/hs.txt form \
the sample NA12878. Two samples with one name is an error.";

/// One sample: its name and the metrics files it was assembled from.
#[derive(Debug)]
pub struct Sample {
    pub name: String,
    /// The directory the files sit in.
    pub dir: PathBuf,
    pub inputs: Vec<Input>,
}

/// A file that looked like a metrics file and did not parse.
#[derive(Debug)]
pub struct Skipped {
    pub path: PathBuf,
    pub reason: String,
}

/// How many bytes the sniff reads: the Picard header lines and the first
/// `## METRICS CLASS` or `## HISTOGRAM` line sit well inside them.
const SNIFF: usize = 16 * 1024;

/// True when the head of a file is a Picard metrics file: its first line
/// is a `StringHeader` (or `VersionHeader`) section line, a
/// `## METRICS CLASS` line or a `## HISTOGRAM` line, it has no NUL byte,
/// and it holds a `## METRICS CLASS` block or a `## HISTOGRAM` section.
pub fn looks_like_metrics(head: &[u8]) -> bool {
    if head.contains(&0) {
        return false;
    }
    let text = String::from_utf8_lossy(head);
    let section =
        |line: &str| line.starts_with("## METRICS CLASS") || line.starts_with("## HISTOGRAM");
    let first = text.lines().next().unwrap_or("").trim_end();
    let header = first.starts_with("## ") && first.ends_with("Header");
    (header || section(first)) && text.lines().any(section)
}

fn sniff(path: &Path) -> Result<bool> {
    let mut file = fs::File::open(path).map_err(|e| Error::io(path, e))?;
    let mut head = Vec::with_capacity(SNIFF);
    file.by_ref()
        .take(SNIFF as u64)
        .read_to_end(&mut head)
        .map_err(|e| Error::io(path, e))?;
    Ok(looks_like_metrics(&head))
}

/// Every metrics file under `dir`, sorted by path. Names starting with a
/// dot are skipped.
pub fn find(dir: &Path) -> Result<Vec<PathBuf>> {
    let mut found = Vec::new();
    walk(dir, &mut found)?;
    found.sort();
    Ok(found)
}

fn walk(dir: &Path, found: &mut Vec<PathBuf>) -> Result<()> {
    let entries = fs::read_dir(dir).map_err(|e| Error::io(dir, e))?;
    for entry in entries {
        let entry = entry.map_err(|e| Error::io(dir, e))?;
        let path = entry.path();
        let hidden = path
            .file_name()
            .and_then(|n| n.to_str())
            .is_some_and(|n| n.starts_with('.'));
        if hidden {
            continue;
        }
        let kind = fs::metadata(&path).map_err(|e| Error::io(&path, e))?;
        if kind.is_dir() {
            walk(&path, found)?;
        } else if kind.is_file() && sniff(&path)? {
            found.push(path);
        }
    }
    Ok(())
}

/// The file name up to its first dot.
fn prefix(path: &Path) -> String {
    path.file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("")
        .split('.')
        .next()
        .unwrap_or("")
        .to_owned()
}

/// The last component of `dir`; `.` and `..` are resolved first.
fn dir_name(dir: &Path) -> String {
    let name = |p: &Path| p.file_name().map(|n| n.to_string_lossy().into_owned());
    name(dir)
        .filter(|n| n != "." && n != "..")
        .or_else(|| dir.canonicalize().ok().and_then(|p| name(&p)))
        .unwrap_or_else(|| dir.display().to_string())
}

/// The samples in `dirs`: every metrics file found, parsed and grouped by
/// [`SAMPLE_RULE`], in name order, plus the files that did not parse.
pub fn samples(dirs: &[PathBuf]) -> Result<(Vec<Sample>, Vec<Skipped>)> {
    let mut by_dir: BTreeMap<PathBuf, Vec<Input>> = BTreeMap::new();
    let mut skipped = Vec::new();
    for dir in dirs {
        if !dir.is_dir() {
            return Err(Error::Input(format!("{}: not a directory", dir.display())));
        }
        for path in find(dir)? {
            match MetricsDoc::read(&path) {
                Ok(doc) => by_dir
                    .entry(path.parent().unwrap_or(dir).to_path_buf())
                    .or_default()
                    .push(Input {
                        path: path.display().to_string(),
                        doc,
                    }),
                Err(e) => skipped.push(Skipped {
                    path,
                    reason: e.to_string(),
                }),
            }
        }
    }
    let mut samples: Vec<Sample> = Vec::new();
    for (dir, inputs) in by_dir {
        for (name, inputs) in group(&dir, inputs) {
            if let Some(other) = samples.iter().find(|s| s.name == name) {
                return Err(Error::Input(format!(
                    "sample {name:?} in two directories, {} and {}; pass them to separate runs or rename one",
                    other.dir.display(),
                    dir.display()
                )));
            }
            samples.push(Sample {
                name,
                dir: dir.clone(),
                inputs,
            });
        }
    }
    samples.sort_by(|a, b| a.name.cmp(&b.name));
    Ok((samples, skipped))
}

/// [`SAMPLE_RULE`] on the files of one directory.
fn group(dir: &Path, inputs: Vec<Input>) -> Vec<(String, Vec<Input>)> {
    let mut by_prefix: BTreeMap<String, Vec<Input>> = BTreeMap::new();
    for input in inputs {
        by_prefix
            .entry(prefix(Path::new(&input.path)))
            .or_default()
            .push(input);
    }
    let classes = |inputs: &[Input]| -> Vec<String> {
        inputs
            .iter()
            .flat_map(|i| i.doc.classes.iter().map(|b| b.class.clone()))
            .collect()
    };
    let groups: Vec<&Vec<Input>> = by_prefix.values().collect();
    let repeats = groups.iter().enumerate().any(|(i, a)| {
        let mine = classes(a);
        groups[i + 1..]
            .iter()
            .any(|b| classes(b).iter().any(|c| mine.contains(c)))
    });
    if repeats {
        return by_prefix.into_iter().collect();
    }
    let name = match by_prefix.len() {
        1 if by_prefix.values().next().is_some_and(|g| g.len() >= 2) => {
            by_prefix.keys().next().cloned().unwrap_or_default()
        }
        _ => dir_name(dir),
    };
    let all: Vec<Input> = by_prefix.into_values().flatten().collect();
    vec![(name, all)]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sniff_accepts_picard_shapes_and_rejects_the_rest() {
        assert!(looks_like_metrics(
            b"## htsjdk.samtools.metrics.StringHeader\n# CollectInsertSizeMetrics\n\n## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n"
        ));
        assert!(looks_like_metrics(
            b"## htsjdk.samtools.metrics.StringHeader\n# MeanQualityByCycle\n\n## HISTOGRAM\tjava.lang.Integer\nCYCLE\tMEAN_QUALITY\n"
        ));
        assert!(looks_like_metrics(b"## METRICS CLASS\tX\nA\n1\n"));
        assert!(!looks_like_metrics(b"chrom\tpos\ttarget\tcoverage\n"));
        assert!(!looks_like_metrics(
            b"{\"schema\": \"bammetrics.verdict/1\"}"
        ));
        assert!(!looks_like_metrics(
            b"## htsjdk\n# x\n\0\n## METRICS CLASS\tX\n"
        ));
        assert!(!looks_like_metrics(b"## a comment\nno sections\n"));
        // A golden record quotes metrics files; its first line is not a header.
        assert!(!looks_like_metrics(
            b"## argv\nx\n## file a.txt\n## METRICS CLASS\tX\n"
        ));
    }

    fn input(name: &str, class: &str) -> Input {
        let text = format!("## METRICS CLASS\t{class}\nSAMPLE\n\n");
        Input {
            path: format!("d/{name}"),
            doc: MetricsDoc::parse(&text).expect("parses"),
        }
    }

    fn names(groups: &[(String, Vec<Input>)]) -> Vec<(&str, usize)> {
        groups.iter().map(|(n, i)| (n.as_str(), i.len())).collect()
    }

    #[test]
    fn multiple_metrics_prefixes_are_samples() {
        let groups = group(
            Path::new("d"),
            vec![
                input("NA1.insert_size_metrics", "InsertSizeMetrics"),
                input("NA1.hs_metrics", "HsMetrics"),
                input("NA2.insert_size_metrics", "InsertSizeMetrics"),
                input("NA2.hs_metrics", "HsMetrics"),
            ],
        );
        assert_eq!(names(&groups), [("NA1", 2), ("NA2", 2)]);
    }

    #[test]
    fn a_directory_without_repeated_classes_is_one_sample() {
        let groups = group(
            Path::new("runs/NA12878"),
            vec![
                input("insert_size.txt", "InsertSizeMetrics"),
                input("hs.txt", "HsMetrics"),
            ],
        );
        assert_eq!(names(&groups), [("NA12878", 2)]);
    }

    #[test]
    fn one_shared_prefix_names_the_sample() {
        let groups = group(
            Path::new("out"),
            vec![
                input("s.insert_size_metrics", "InsertSizeMetrics"),
                input("s.hs_metrics", "HsMetrics"),
            ],
        );
        assert_eq!(names(&groups), [("s", 2)]);
        let one = group(Path::new("out"), vec![input("s.hs_metrics", "HsMetrics")]);
        assert_eq!(names(&one), [("out", 1)]);
    }
}

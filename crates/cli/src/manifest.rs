//! `PREFIX.manifest.json` (schema `bammetrics.manifest/1`): one record per
//! `qc` step with the fields `--resume` compares (docs/05 section 3): the
//! argv hash, the tool version, the preset hash, and the input's size,
//! mtime and a hash of its first and last 64 KiB; then the step's outputs
//! with tool and kind, and its status.

use std::fs::{self, File};
use std::io::{Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

use bammetrics_verdict::bytes_hash;
use bammetrics_verdict::json::Json;

use crate::error::{Error, Result};

pub const SCHEMA: &str = "bammetrics.manifest/1";

/// Bytes hashed at each end of an input.
const EDGE: u64 = 64 * 1024;

/// What identifies an input file without reading all of it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fingerprint {
    pub size: u64,
    /// Seconds and nanoseconds since the epoch, `1730000000.123456789`.
    pub mtime: String,
    /// FNV-1a 64 of the first 64 KiB followed by the last 64 KiB.
    pub edges: String,
}

impl Fingerprint {
    pub fn of(path: &Path) -> Result<Self> {
        let io = |e| Error::io(path, e);
        let meta = fs::metadata(path).map_err(io)?;
        let mtime = meta
            .modified()
            .ok()
            .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
            .map_or_else(String::new, |d| {
                format!("{}.{:09}", d.as_secs(), d.subsec_nanos())
            });
        let size = meta.len();
        let mut file = File::open(path).map_err(io)?;
        let mut head = Vec::new();
        (&mut file).take(EDGE).read_to_end(&mut head).map_err(io)?;
        let mut tail = Vec::new();
        file.seek(SeekFrom::Start(size.saturating_sub(EDGE)))
            .map_err(io)?;
        file.take(EDGE).read_to_end(&mut tail).map_err(io)?;
        Ok(Self {
            size,
            mtime,
            edges: bytes_hash(head.into_iter().chain(tail)),
        })
    }
}

/// The fields `--resume` compares: a step is skipped when they all match
/// the manifest's record, the record says `ok` and its outputs exist.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Key {
    pub argv_hash: String,
    pub tool_version: String,
    pub preset_hash: String,
    pub input: PathBuf,
    pub fingerprint: Fingerprint,
}

/// One file a step wrote.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Output {
    pub path: PathBuf,
    /// The Picard tool, or `bammetrics` for the verdict and the page.
    pub tool: String,
    /// `metrics`, `verdict` or `report`.
    pub kind: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Ok,
    Failed,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Ok => "ok",
            Status::Failed => "failed",
        }
    }
}

/// One step of a run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepRecord {
    /// Unique in the run: the tool name, or `CollectHsMetrics.<BAIT_SET>`
    /// with several hs runs, or `verdict` and `report`.
    pub id: String,
    /// `pass`, `own-read`, `post` or `report`.
    pub stage: String,
    pub tool: String,
    /// The tool's Picard arguments, as `bammetrics <tool> ARGV...` takes them.
    pub argv: Vec<String>,
    pub key: Key,
    pub outputs: Vec<Output>,
    pub status: Status,
    pub error: Option<String>,
}

impl StepRecord {
    /// True when `self`, the record a new run would write before running,
    /// can stand for `old`: same key, `old` finished, and every file it
    /// wrote is still there.
    pub fn resumes(&self, old: &StepRecord) -> bool {
        self.id == old.id
            && self.key == old.key
            && old.status == Status::Ok
            && old.outputs.iter().all(|o| o.path.is_file())
    }

    fn to_json(&self) -> Json {
        let text = |s: &str| Json::text(s);
        let path = |p: &Path| Json::text(&p.display().to_string());
        let f = &self.key.fingerprint;
        Json::Object(vec![
            ("id".to_owned(), text(&self.id)),
            ("stage".to_owned(), text(&self.stage)),
            ("tool".to_owned(), text(&self.tool)),
            (
                "argv".to_owned(),
                Json::Array(self.argv.iter().map(|a| text(a)).collect()),
            ),
            ("argv_hash".to_owned(), text(&self.key.argv_hash)),
            ("tool_version".to_owned(), text(&self.key.tool_version)),
            ("preset_hash".to_owned(), text(&self.key.preset_hash)),
            (
                "input".to_owned(),
                Json::Object(vec![
                    ("path".to_owned(), path(&self.key.input)),
                    ("size".to_owned(), Json::Number(f.size as f64)),
                    ("mtime".to_owned(), text(&f.mtime)),
                    ("edges_hash".to_owned(), text(&f.edges)),
                ]),
            ),
            (
                "outputs".to_owned(),
                Json::Array(
                    self.outputs
                        .iter()
                        .map(|o| {
                            Json::Object(vec![
                                ("path".to_owned(), path(&o.path)),
                                ("tool".to_owned(), text(&o.tool)),
                                ("kind".to_owned(), text(&o.kind)),
                            ])
                        })
                        .collect(),
                ),
            ),
            ("status".to_owned(), text(self.status.name())),
            ("error".to_owned(), Json::option_text(self.error.as_deref())),
        ])
    }

    fn from_json(json: &Json) -> Option<Self> {
        let text = |j: &Json, key: &str| j.get(key)?.as_str().map(str::to_owned);
        let input = json.get("input")?;
        let outputs = json
            .get("outputs")?
            .as_array()?
            .iter()
            .map(|o| {
                Some(Output {
                    path: PathBuf::from(text(o, "path")?),
                    tool: text(o, "tool")?,
                    kind: text(o, "kind")?,
                })
            })
            .collect::<Option<Vec<_>>>()?;
        let status = match json.get("status")?.as_str()? {
            "ok" => Status::Ok,
            "failed" => Status::Failed,
            _ => return None,
        };
        Some(Self {
            id: text(json, "id")?,
            stage: text(json, "stage")?,
            tool: text(json, "tool")?,
            argv: json
                .get("argv")?
                .as_array()?
                .iter()
                .map(|a| a.as_str().map(str::to_owned))
                .collect::<Option<Vec<_>>>()?,
            key: Key {
                argv_hash: text(json, "argv_hash")?,
                tool_version: text(json, "tool_version")?,
                preset_hash: text(json, "preset_hash")?,
                input: PathBuf::from(text(input, "path")?),
                fingerprint: Fingerprint {
                    size: input.get("size")?.as_f64()? as u64,
                    mtime: text(input, "mtime")?,
                    edges: text(input, "edges_hash")?,
                },
            },
            outputs,
            status,
            error: json.get("error")?.as_str().map(str::to_owned),
        })
    }
}

/// The run: who, what, and every step.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Manifest {
    pub tool_version: String,
    pub sample: String,
    pub bam: PathBuf,
    pub prefix: PathBuf,
    pub preset_name: String,
    pub preset_hash: String,
    pub steps: Vec<StepRecord>,
}

impl Manifest {
    pub fn to_json(&self) -> String {
        let path = |p: &Path| Json::text(&p.display().to_string());
        Json::Object(vec![
            ("schema".to_owned(), Json::text(SCHEMA)),
            ("tool_version".to_owned(), Json::text(&self.tool_version)),
            ("sample".to_owned(), Json::text(&self.sample)),
            ("bam".to_owned(), path(&self.bam)),
            ("prefix".to_owned(), path(&self.prefix)),
            (
                "preset".to_owned(),
                Json::Object(vec![
                    ("name".to_owned(), Json::text(&self.preset_name)),
                    ("hash".to_owned(), Json::text(&self.preset_hash)),
                ]),
            ),
            (
                "steps".to_owned(),
                Json::Array(self.steps.iter().map(StepRecord::to_json).collect()),
            ),
        ])
        .render()
    }

    /// Parse manifest text; `None` for anything that is not a whole
    /// `bammetrics.manifest/1` document.
    pub fn parse(text: &str) -> Option<Self> {
        let json = Json::parse(text).ok()?;
        if json.get("schema")?.as_str()? != SCHEMA {
            return None;
        }
        let text = |key: &str| json.get(key)?.as_str().map(str::to_owned);
        let preset = json.get("preset")?;
        Some(Self {
            tool_version: text("tool_version")?,
            sample: text("sample")?,
            bam: PathBuf::from(text("bam")?),
            prefix: PathBuf::from(text("prefix")?),
            preset_name: preset.get("name")?.as_str()?.to_owned(),
            preset_hash: preset.get("hash")?.as_str()?.to_owned(),
            steps: json
                .get("steps")?
                .as_array()?
                .iter()
                .map(StepRecord::from_json)
                .collect::<Option<Vec<_>>>()?,
        })
    }

    /// The manifest at `path`: `Ok(None)` when there is none, `Err` with
    /// the reason when it cannot be used (a killed run can leave it cut
    /// short).
    pub fn read(path: &Path) -> std::result::Result<Option<Self>, String> {
        match fs::read_to_string(path) {
            Ok(text) => Self::parse(&text)
                .map(Some)
                .ok_or_else(|| format!("{} is not a whole manifest", path.display())),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(format!("{}: {e}", path.display())),
        }
    }

    pub fn step(&self, id: &str) -> Option<&StepRecord> {
        self.steps.iter().find(|s| s.id == id)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record() -> StepRecord {
        StepRecord {
            id: "CollectInsertSizeMetrics".to_owned(),
            stage: "pass".to_owned(),
            tool: "CollectInsertSizeMetrics".to_owned(),
            argv: vec!["INPUT=a \"b\".bam".to_owned(), "OUTPUT=out/s".to_owned()],
            key: Key {
                argv_hash: "0123456789abcdef".to_owned(),
                tool_version: "0.12.0".to_owned(),
                preset_hash: "fedcba9876543210".to_owned(),
                input: PathBuf::from("a \"b\".bam"),
                fingerprint: Fingerprint {
                    size: 123_456_789_012,
                    mtime: "1730000000.000000001".to_owned(),
                    edges: "00000000000000ff".to_owned(),
                },
            },
            outputs: vec![Output {
                path: PathBuf::from("out/s.insert_size_metrics"),
                tool: "CollectInsertSizeMetrics".to_owned(),
                kind: "metrics".to_owned(),
            }],
            status: Status::Ok,
            error: None,
        }
    }

    fn manifest() -> Manifest {
        let mut failed = record();
        failed.id = "CollectHsMetrics".to_owned();
        failed.status = Status::Failed;
        failed.error = Some("no such file".to_owned());
        failed.outputs.clear();
        Manifest {
            tool_version: "0.12.0".to_owned(),
            sample: "NA12878".to_owned(),
            bam: PathBuf::from("a.bam"),
            prefix: PathBuf::from("out/s"),
            preset_name: "generic".to_owned(),
            preset_hash: "6dd89e0000000000".to_owned(),
            steps: vec![record(), failed],
        }
    }

    #[test]
    fn a_manifest_reads_back_as_written() {
        let m = manifest();
        assert_eq!(Manifest::parse(&m.to_json()), Some(m));
    }

    #[test]
    fn a_cut_short_or_foreign_manifest_does_not_parse() {
        let text = manifest().to_json();
        assert!(Manifest::parse(&text[..text.len() / 2]).is_none());
        assert!(Manifest::parse("").is_none());
        assert!(Manifest::parse(&text.replace(SCHEMA, "bammetrics.verdict/1")).is_none());
    }

    #[test]
    fn resume_needs_the_same_key_a_finished_step_and_its_files() {
        let dir = std::env::temp_dir().join(format!("bammetrics-manifest-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("dir");
        let file = dir.join("s.insert_size_metrics");
        fs::write(&file, "x").expect("file");
        let mut old = record();
        old.outputs[0].path = file.clone();
        let new = old.clone();
        assert!(new.resumes(&old));
        let mut moved = new.clone();
        moved.key.fingerprint.mtime = "1730000001.000000000".to_owned();
        assert!(!moved.resumes(&old));
        let mut other_version = new.clone();
        other_version.key.tool_version = "0.12.1".to_owned();
        assert!(!other_version.resumes(&old));
        let mut failed = old.clone();
        failed.status = Status::Failed;
        assert!(!new.resumes(&failed));
        fs::remove_file(&file).expect("remove");
        assert!(!new.resumes(&old));
        let _ = fs::remove_dir_all(&dir);
    }

    #[test]
    fn the_fingerprint_hashes_both_ends() {
        let dir = std::env::temp_dir().join(format!("bammetrics-print-{}", std::process::id()));
        fs::create_dir_all(&dir).expect("dir");
        let path = dir.join("big.bin");
        let mut bytes = vec![7u8; 200 * 1024];
        fs::write(&path, &bytes).expect("write");
        let first = Fingerprint::of(&path).expect("print");
        assert_eq!(first.size, 200 * 1024);
        bytes[100 * 1024] = 8;
        fs::write(&path, &bytes).expect("write middle");
        assert_eq!(Fingerprint::of(&path).expect("print").edges, first.edges);
        let last = bytes.len() - 1;
        bytes[last] = 8;
        fs::write(&path, &bytes).expect("write tail");
        assert_ne!(Fingerprint::of(&path).expect("print").edges, first.edges);
        let _ = fs::remove_dir_all(&dir);
    }
}

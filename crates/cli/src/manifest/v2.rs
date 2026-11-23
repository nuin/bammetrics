//! `PREFIX.manifest.json` of a chain (schema `bammetrics.manifest/2`,
//! `docs/schema/manifest-2.json`): the source's fingerprint, then one
//! record per step with what `chain --resume` compares (docs/chaining.md
//! section 27): the verb, the canonical argv and its hash, the files the
//! step reads and writes with their fingerprints, its status and timing,
//! and the ids of the transforms in front of it. The kept products follow.
//!
//! The text is read through [`bammetrics_formats::json`], which refuses
//! nesting past 64 levels, and at most [`CAP`] bytes are read, so a
//! damaged or hostile file is an error and never a crash.

use std::path::{Path, PathBuf};

use bammetrics_formats::json::Json;

use super::Fingerprint;

pub const SCHEMA: &str = "bammetrics.manifest/2";

/// The most bytes of a manifest read back: 16 MiB. A chain of 20 steps
/// on the fixture writes about 20 KB.
pub const CAP: u64 = 16 << 20;

/// A file a step reads or writes, with its fingerprint when it is a
/// regular file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FileRecord {
    pub path: PathBuf,
    /// What an output is (`metrics`, `bam`, `index`, `kept`, `verdict`);
    /// `None` for an input.
    pub kind: Option<String>,
    /// `None` for a pipe, a device or a file that is not there.
    pub print: Option<Fingerprint>,
}

impl FileRecord {
    /// An input, fingerprinted now.
    pub fn input(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
            kind: None,
            print: fingerprint(path),
        }
    }

    /// An output of `kind`, fingerprinted now.
    pub fn output(path: &Path, kind: &str) -> Self {
        Self {
            path: path.to_path_buf(),
            kind: Some(kind.to_owned()),
            print: fingerprint(path),
        }
    }

    /// True when the file on disk has the recorded fingerprint. A record
    /// without one never matches: a pipe cannot be compared.
    pub fn on_disk(&self) -> bool {
        self.print
            .as_ref()
            .is_some_and(|print| fingerprint(&self.path).as_ref() == Some(print))
    }

    /// [`FileRecord::on_disk`] for a record of a run in working directory
    /// `dir`: a relative path is read under it.
    pub fn on_disk_in(&self, dir: &Path) -> bool {
        self.print
            .as_ref()
            .is_some_and(|print| fingerprint(&dir.join(&self.path)).as_ref() == Some(print))
    }

    fn to_json(&self) -> Json {
        let mut members = vec![("path".to_owned(), Json::text(&path_text(&self.path)))];
        if let Some(kind) = &self.kind {
            members.push(("kind".to_owned(), Json::text(kind)));
        }
        let print = self.print.as_ref();
        members.extend([
            (
                "size".to_owned(),
                print.map_or(Json::Null, |p| Json::Number(p.size.to_string())),
            ),
            (
                "mtime".to_owned(),
                Json::option_text(print.map(|p| p.mtime.as_str())),
            ),
            (
                "edges_hash".to_owned(),
                Json::option_text(print.map(|p| p.edges.as_str())),
            ),
        ]);
        Json::Object(members)
    }

    fn from_json(json: &Json, output: bool) -> Option<Self> {
        let kind = match json.get("kind") {
            Some(kind) => Some(kind.as_str()?.to_owned()),
            None if output => return None,
            None => None,
        };
        let print = match (
            json.get("size")?,
            json.get("mtime")?,
            json.get("edges_hash")?,
        ) {
            (Json::Null, Json::Null, Json::Null) => None,
            (Json::Number(size), Json::Str(mtime), Json::Str(edges)) => Some(Fingerprint {
                size: size.parse().ok()?,
                mtime: mtime.clone(),
                edges: edges.clone(),
            }),
            _ => return None,
        };
        Some(Self {
            path: PathBuf::from(json.get("path")?.as_str()?),
            kind,
            print,
        })
    }
}

/// The fingerprint of a regular file; `None` for anything else.
pub fn fingerprint(path: &Path) -> Option<Fingerprint> {
    std::fs::metadata(path)
        .ok()
        .filter(|meta| meta.is_file() && !path.starts_with("/dev"))
        .and_then(|_| Fingerprint::of(path).ok())
}

fn path_text(path: &Path) -> String {
    path.display().to_string()
}

/// `1.234`: milliseconds as seconds.
fn seconds(millis: u64) -> Json {
    Json::Number(format!("{}.{:03}", millis / 1000, millis % 1000))
}

fn millis(json: &Json) -> Option<u64> {
    match json {
        Json::Number(text) => {
            let value: f64 = text.parse().ok()?;
            (value.is_finite() && value >= 0.0).then(|| (value * 1000.0).round() as u64)
        }
        _ => None,
    }
}

/// One step of a chain run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StepRecord {
    /// Unique in the chain: a transform's stream (`marked`), else
    /// `VERB@STREAM` (`insert-size@marked`), with `#2` for the second of
    /// one verb on one stream; `qc` for the qc step.
    pub id: String,
    /// Its place in the chain, from 1.
    pub number: usize,
    pub verb: String,
    /// The Picard tool; `None` for `write` and `qc`.
    pub tool: Option<String>,
    /// `transform`, `tap`, `write`, `split` or `qc`.
    pub role: String,
    /// The stream it reads.
    pub reads: String,
    /// The stream a transform leaves.
    pub leaves: Option<String>,
    /// The ids of the transforms in front of it, first first.
    pub upstream: Vec<String>,
    /// What it runs, canonical: the Picard argv of a tool step, the
    /// options of `write` and `qc` as `KEY=VALUE`.
    pub argv: Vec<String>,
    /// FNV-1a 64 of the tool version, the argv, the formats, the sample
    /// and the source's validation stringency.
    pub options_hash: String,
    /// Every file it reads by name: the reference with its `.fai`,
    /// dictionaries, interval and read lists.
    pub inputs: Vec<FileRecord>,
    pub outputs: Vec<FileRecord>,
    /// `ok`, `failed`, `stopped` or `not run`.
    pub status: String,
    pub error: Option<String>,
    /// True when this run kept the record of the run before.
    pub skipped: bool,
    /// When the run that ran it started (`2026-01-02 03:04:05 UTC`).
    pub started: String,
    /// The wall time of that run: the steps of a chain share its reads.
    pub millis: u64,
}

impl StepRecord {
    fn to_json(&self) -> Json {
        let files =
            |files: &[FileRecord]| Json::Array(files.iter().map(FileRecord::to_json).collect());
        Json::Object(vec![
            ("id".to_owned(), Json::text(&self.id)),
            ("number".to_owned(), Json::Number(self.number.to_string())),
            ("verb".to_owned(), Json::text(&self.verb)),
            ("tool".to_owned(), Json::option_text(self.tool.as_deref())),
            ("role".to_owned(), Json::text(&self.role)),
            ("reads".to_owned(), Json::text(&self.reads)),
            (
                "leaves".to_owned(),
                Json::option_text(self.leaves.as_deref()),
            ),
            (
                "upstream".to_owned(),
                Json::strings(self.upstream.iter().map(String::as_str)),
            ),
            (
                "argv".to_owned(),
                Json::strings(self.argv.iter().map(String::as_str)),
            ),
            ("options_hash".to_owned(), Json::text(&self.options_hash)),
            ("inputs".to_owned(), files(&self.inputs)),
            ("outputs".to_owned(), files(&self.outputs)),
            ("status".to_owned(), Json::text(&self.status)),
            ("error".to_owned(), Json::option_text(self.error.as_deref())),
            ("skipped".to_owned(), Json::Bool(self.skipped)),
            ("started".to_owned(), Json::text(&self.started)),
            ("seconds".to_owned(), seconds(self.millis)),
        ])
    }

    fn from_json(json: &Json) -> Option<Self> {
        let text = |key: &str| json.get(key)?.as_str().map(str::to_owned);
        let optional = |key: &str| match json.get(key)? {
            Json::Null => Some(None),
            Json::Str(text) => Some(Some(text.clone())),
            _ => None,
        };
        let strings = |key: &str| {
            json.get(key)?
                .as_array()?
                .iter()
                .map(|item| item.as_str().map(str::to_owned))
                .collect::<Option<Vec<_>>>()
        };
        let files = |key: &str, output: bool| {
            json.get(key)?
                .as_array()?
                .iter()
                .map(|item| FileRecord::from_json(item, output))
                .collect::<Option<Vec<_>>>()
        };
        let number = match json.get("number")? {
            Json::Number(text) => text.parse().ok()?,
            _ => return None,
        };
        let skipped = match json.get("skipped")? {
            Json::Bool(b) => *b,
            _ => return None,
        };
        Some(Self {
            id: text("id")?,
            number,
            verb: text("verb")?,
            tool: optional("tool")?,
            role: text("role")?,
            reads: text("reads")?,
            leaves: optional("leaves")?,
            upstream: strings("upstream")?,
            argv: strings("argv")?,
            options_hash: text("options_hash")?,
            inputs: files("inputs", false)?,
            outputs: files("outputs", true)?,
            status: text("status")?,
            error: optional("error")?,
            skipped,
            started: text("started")?,
            millis: millis(json.get("seconds")?)?,
        })
    }
}

/// A stream a transform kept as a BAM (`--keep`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Kept {
    /// The id of the transform.
    pub step: String,
    pub stream: String,
    pub file: FileRecord,
}

impl Kept {
    fn to_json(&self) -> Json {
        let mut members = vec![
            ("step".to_owned(), Json::text(&self.step)),
            ("stream".to_owned(), Json::text(&self.stream)),
        ];
        if let Json::Object(file) = self.file.to_json() {
            members.extend(file.into_iter().filter(|(key, _)| key != "kind"));
        }
        Json::Object(members)
    }

    fn from_json(json: &Json) -> Option<Self> {
        Some(Self {
            step: json.get("step")?.as_str()?.to_owned(),
            stream: json.get("stream")?.as_str()?.to_owned(),
            file: FileRecord {
                kind: Some("kept".to_owned()),
                ..FileRecord::from_json(json, false)?
            },
        })
    }
}

/// A chain run.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainManifest {
    pub tool_version: String,
    pub sample: String,
    /// The working directory of the run, resolved: every relative path of
    /// the manifest is relative to it.
    pub directory: PathBuf,
    pub source: FileRecord,
    pub reference: Option<FileRecord>,
    pub prefix: PathBuf,
    /// `STRICT`, `LENIENT` or `SILENT`: how a SAM text source was read.
    pub stringency: String,
    pub started: String,
    pub millis: u64,
    /// The exit code of the run.
    pub exit: u8,
    pub steps: Vec<StepRecord>,
    pub kept: Vec<Kept>,
}

impl ChainManifest {
    pub fn to_json(&self) -> String {
        Json::Object(vec![
            ("schema".to_owned(), Json::text(SCHEMA)),
            ("tool_version".to_owned(), Json::text(&self.tool_version)),
            ("verb".to_owned(), Json::text(crate::chain::VERB)),
            ("sample".to_owned(), Json::text(&self.sample)),
            (
                "directory".to_owned(),
                Json::text(&path_text(&self.directory)),
            ),
            ("source".to_owned(), self.source.to_json()),
            (
                "reference".to_owned(),
                self.reference
                    .as_ref()
                    .map_or(Json::Null, FileRecord::to_json),
            ),
            ("prefix".to_owned(), Json::text(&path_text(&self.prefix))),
            ("stringency".to_owned(), Json::text(&self.stringency)),
            ("started".to_owned(), Json::text(&self.started)),
            ("seconds".to_owned(), seconds(self.millis)),
            ("exit".to_owned(), Json::Number(self.exit.to_string())),
            (
                "steps".to_owned(),
                Json::Array(self.steps.iter().map(StepRecord::to_json).collect()),
            ),
            (
                "kept".to_owned(),
                Json::Array(self.kept.iter().map(Kept::to_json).collect()),
            ),
        ])
        .render()
    }

    /// The manifest in `json`; `None` for anything that is not a whole
    /// `bammetrics.manifest/2` document.
    fn from_json(json: &Json) -> Option<Self> {
        if json.get("schema")?.as_str()? != SCHEMA {
            return None;
        }
        let text = |key: &str| json.get(key)?.as_str().map(str::to_owned);
        let reference = match json.get("reference")? {
            Json::Null => None,
            other => Some(FileRecord::from_json(other, false)?),
        };
        let exit = match json.get("exit")? {
            Json::Number(text) => text.parse().ok()?,
            _ => return None,
        };
        Some(Self {
            tool_version: text("tool_version")?,
            sample: text("sample")?,
            directory: PathBuf::from(text("directory")?),
            source: FileRecord::from_json(json.get("source")?, false)?,
            reference,
            prefix: PathBuf::from(text("prefix")?),
            stringency: text("stringency")?,
            started: text("started")?,
            millis: millis(json.get("seconds")?)?,
            exit,
            steps: json
                .get("steps")?
                .as_array()?
                .iter()
                .map(StepRecord::from_json)
                .collect::<Option<Vec<_>>>()?,
            kept: json
                .get("kept")?
                .as_array()?
                .iter()
                .map(Kept::from_json)
                .collect::<Option<Vec<_>>>()?,
        })
    }

    /// Parse manifest text; `None` for anything that is not a whole
    /// `bammetrics.manifest/2` document.
    pub fn parse(text: &str) -> Option<Self> {
        Self::from_json(&Json::parse(text).ok()?)
    }

    pub fn step(&self, id: &str) -> Option<&StepRecord> {
        self.steps.iter().find(|step| step.id == id)
    }
}

/// A manifest of either schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnyManifest {
    /// `bammetrics.manifest/1`: a `qc` run.
    Qc(Box<super::Manifest>),
    /// `bammetrics.manifest/2`: a chain run.
    Chain(Box<ChainManifest>),
}

/// The manifest in `text`, of either schema; `Err` with the reason when it
/// is none.
pub fn parse_any(text: &str) -> Result<AnyManifest, String> {
    // The reader with a depth limit decides first: the reader of schema 1
    // has none.
    let json = Json::parse(text).map_err(|e| format!("not JSON: {e}"))?;
    match json.get("schema").and_then(Json::as_str) {
        Some(SCHEMA) => ChainManifest::from_json(&json)
            .map(|m| AnyManifest::Chain(Box::new(m)))
            .ok_or_else(|| format!("not a whole {SCHEMA} document")),
        Some(super::SCHEMA) => super::Manifest::parse(text)
            .map(|m| AnyManifest::Qc(Box::new(m)))
            .ok_or_else(|| format!("not a whole {} document", super::SCHEMA)),
        Some(other) => Err(format!("schema {other:?} is no manifest schema")),
        None => Err("no schema member".to_owned()),
    }
}

/// The manifest at `path`, of either schema: `Ok(None)` when there is
/// none, `Err` with the reason when it cannot be used. At most [`CAP`]
/// bytes are read.
pub fn read_any(path: &Path) -> Result<Option<AnyManifest>, String> {
    match crate::chain::read_capped(path, CAP) {
        Ok(text) => parse_any(&text)
            .map(Some)
            .map_err(|why| format!("{}: {why}", path.display())),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
        Err(e) => Err(format!("{}: {e}", path.display())),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn print(n: u64) -> Option<Fingerprint> {
        Some(Fingerprint {
            size: n,
            mtime: format!("1730000000.00000000{n}"),
            edges: format!("{n:016x}"),
        })
    }

    fn manifest() -> ChainManifest {
        let step = StepRecord {
            id: "insert-size@marked".to_owned(),
            number: 2,
            verb: "insert-size".to_owned(),
            tool: Some("CollectInsertSizeMetrics".to_owned()),
            role: "tap".to_owned(),
            reads: "marked".to_owned(),
            leaves: None,
            upstream: vec!["marked".to_owned()],
            argv: vec![
                "INPUT=a \"b\".bam".to_owned(),
                "OUTPUT=out/s.txt".to_owned(),
            ],
            options_hash: "0123456789abcdef".to_owned(),
            inputs: vec![FileRecord {
                path: PathBuf::from("ref.fa"),
                kind: None,
                print: print(3),
            }],
            outputs: vec![FileRecord {
                path: PathBuf::from("out/s.insert_size_metrics"),
                kind: Some("metrics".to_owned()),
                print: print(4),
            }],
            status: "ok".to_owned(),
            error: None,
            skipped: false,
            started: "2026-11-22 10:00:00 UTC".to_owned(),
            millis: 1234,
        };
        let mut failed = step.clone();
        failed.id = "write@marked".to_owned();
        failed.tool = None;
        failed.status = "failed".to_owned();
        failed.error = Some("no space".to_owned());
        failed.outputs.clear();
        failed.skipped = true;
        ChainManifest {
            tool_version: "0.31.0".to_owned(),
            sample: "NA12878".to_owned(),
            directory: PathBuf::from("/work"),
            source: FileRecord {
                path: PathBuf::from("in.bam"),
                kind: None,
                print: print(1),
            },
            reference: None,
            prefix: PathBuf::from("out/s"),
            stringency: "STRICT".to_owned(),
            started: "2026-11-22 10:00:00 UTC".to_owned(),
            millis: 50,
            exit: 1,
            steps: vec![step, failed],
            kept: vec![Kept {
                step: "marked".to_owned(),
                stream: "marked".to_owned(),
                file: FileRecord {
                    path: PathBuf::from("out/s.marked.bam"),
                    kind: Some("kept".to_owned()),
                    print: print(5),
                },
            }],
        }
    }

    #[test]
    fn a_manifest_reads_back_as_written() {
        let m = manifest();
        assert_eq!(ChainManifest::parse(&m.to_json()), Some(m.clone()));
        let mut piped = m;
        piped.source.print = None;
        assert_eq!(ChainManifest::parse(&piped.to_json()), Some(piped));
    }

    #[test]
    fn either_schema_is_read_and_anything_else_is_refused() {
        let text = manifest().to_json();
        assert!(matches!(parse_any(&text), Ok(AnyManifest::Chain(_))));
        assert!(parse_any(&text[..text.len() / 2]).is_err());
        assert!(parse_any("{}").is_err());
        let foreign = text.replace(SCHEMA, "bammetrics.verdict/1");
        assert!(
            parse_any(&foreign)
                .unwrap_err()
                .contains("no manifest schema")
        );
        // Nesting past the reader's limit is an error, not a crash.
        let deep = format!("{}{}", "[".repeat(100_000), "]".repeat(100_000));
        assert!(parse_any(&deep).is_err());
        // A schema 1 document of the qc verb reads as one.
        let qc = super::super::Manifest {
            tool_version: "0.31.0".to_owned(),
            sample: "s".to_owned(),
            bam: PathBuf::from("a.bam"),
            prefix: PathBuf::from("out/s"),
            preset_name: "generic".to_owned(),
            preset_hash: "6dd89e0000000000".to_owned(),
            steps: Vec::new(),
        };
        assert!(matches!(parse_any(&qc.to_json()), Ok(AnyManifest::Qc(_))));
    }

    #[test]
    fn a_record_on_disk_matches_its_file_only() {
        let dir = std::env::temp_dir().join(format!("bammetrics-v2-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let path = dir.join("s.txt");
        std::fs::write(&path, "x").expect("write");
        let record = FileRecord::output(&path, "metrics");
        assert!(record.on_disk());
        std::fs::write(&path, "xy").expect("rewrite");
        assert!(!record.on_disk());
        std::fs::remove_file(&path).expect("remove");
        assert!(!FileRecord::output(&path, "metrics").on_disk());
        assert!(FileRecord::input(Path::new("/dev/stdin")).print.is_none());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_file_over_the_cap_is_refused() {
        let dir = std::env::temp_dir().join(format!("bammetrics-v2cap-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let path = dir.join("big.manifest.json");
        std::fs::write(&path, vec![b' '; usize::try_from(CAP).expect("fits") + 1]).expect("write");
        assert!(read_any(&path).unwrap_err().contains("larger than 16 MiB"));
        assert_eq!(read_any(&dir.join("none.json")), Ok(None));
        let _ = std::fs::remove_dir_all(&dir);
    }
}

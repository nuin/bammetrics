//! `doctor`: a pre-flight check of one BAM, CRAM or SAM file before a long
//! run (docs/doctor.md, docs/07 W07). It reads the header, the last bytes
//! of the file and the first records (10,000 by default), never the whole
//! file, and gives one [`Check`] per question with a status and a reason.
//!
//! Nothing here writes a file. A file that cannot be opened is an error;
//! everything a file that opens can be wrong about is a check.

mod checks;
pub mod input;
pub mod records;
pub mod reference;

use std::collections::HashSet;
use std::fmt;
use std::path::{Path, PathBuf};

use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::genome_build::Build;
use crate::sam_text::Stringency;

pub use checks::NAMES;

/// The shortest magic of the kinds read: `BAM\1` inflated, `CRAM`, and
/// the gzip magic and flags. A shorter file is none of them.
const MAGIC_BYTES: u64 = 4;

/// Records read by default.
pub const DEFAULT_RECORDS: u64 = 10_000;

/// The status of one check.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Status {
    /// The check does not apply (a CRAM check on a BAM), or what it needs
    /// could not be read.
    Skip,
    Pass,
    Warn,
    Fail,
}

impl Status {
    pub fn name(self) -> &'static str {
        match self {
            Status::Skip => "skip",
            Status::Pass => "pass",
            Status::Warn => "warn",
            Status::Fail => "fail",
        }
    }
}

impl fmt::Display for Status {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// One check: its name (one of [`NAMES`]), its status, what was found and
/// why that status.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Check {
    pub name: &'static str,
    pub status: Status,
    pub value: String,
    pub reason: String,
}

impl Check {
    fn new(name: &'static str, status: Status, value: impl Into<String>) -> Self {
        Self {
            name,
            status,
            value: value.into(),
            reason: String::new(),
        }
    }

    fn because(mut self, reason: impl Into<String>) -> Self {
        self.reason = reason.into();
        self
    }
}

/// The sex a sample is expected to be.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sex {
    Male,
    Female,
}

impl Sex {
    /// `male`, `m`, `xy`, `female`, `f`, `xx`, any case; `None` otherwise.
    pub fn parse(text: &str) -> Option<Self> {
        match text.trim().to_ascii_lowercase().as_str() {
            "male" | "m" | "xy" => Some(Sex::Male),
            "female" | "f" | "xx" => Some(Sex::Female),
            _ => None,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Sex::Male => "male",
            Sex::Female => "female",
        }
    }
}

/// The sample sheet row(s) for the sample asked for (docs/07 D13: one
/// sample per call). The command reads the sheet; this is what it found.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SheetRow {
    /// The sheet, as named on the command line.
    pub sheet: PathBuf,
    /// The sample asked for.
    pub sample: String,
    /// Rows that name the sample: 1 when the sheet is right.
    pub rows: usize,
    /// The `bam` column, a path relative to the sheet's directory or
    /// absolute; empty or absent when the sheet has no such column.
    pub bam: Option<String>,
    pub sex: Option<String>,
    pub build: Option<String>,
    pub library: Option<String>,
}

/// What `doctor` is asked.
#[derive(Clone, Debug)]
pub struct Options {
    /// Records read from the start of the file; 0 reads none.
    pub records: u64,
    /// `-R`: the reference FASTA (with its `.fai`).
    pub reference: Option<PathBuf>,
    /// Compute the MD5 of every `@SQ` contig with an `M5` from the FASTA
    /// in place of the `.dict` beside it.
    pub check_md5: bool,
    pub expect_build: Option<Build>,
    pub expect_sex: Option<Sex>,
    pub sheet: Option<SheetRow>,
    /// Threads of the CRAM decoder and the MD5 workers.
    pub threads: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            records: DEFAULT_RECORDS,
            reference: None,
            check_md5: false,
            expect_build: None,
            expect_sex: None,
            sheet: None,
            threads: 1,
        }
    }
}

/// Every check of one file, in [`NAMES`] order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Report {
    pub path: PathBuf,
    /// `BAM`, `CRAM 3.0`, `SAM`, `SAM (gzip)`; empty when unknown.
    pub kind: String,
    pub checks: Vec<Check>,
}

impl Report {
    /// The worst status of any check; `pass` when every check is skipped.
    pub fn status(&self) -> Status {
        self.checks
            .iter()
            .map(|c| c.status)
            .max()
            .filter(|s| *s != Status::Skip)
            .unwrap_or(Status::Pass)
    }

    /// Checks of `status`.
    pub fn count(&self, status: Status) -> usize {
        self.checks.iter().filter(|c| c.status == status).count()
    }

    pub fn get(&self, name: &str) -> Option<&Check> {
        self.checks.iter().find(|c| c.name == name)
    }
}

/// What the file is, its header and the first records, before any check is
/// judged.
pub(crate) struct Facts {
    pub kind: Option<Kind>,
    pub kind_error: Option<String>,
    /// The size of a file shorter than any magic.
    pub short: Option<u64>,
    pub cram_version: Option<(u8, u8)>,
    pub eof: Option<Eof>,
    pub header: Option<std::result::Result<input::RawHeader, String>>,
    pub lines: Vec<input::Line>,
    pub stats: Option<records::Stats>,
    /// Why no record was read, when none could be.
    pub records_error: Option<String>,
    /// Why no record was read on purpose.
    pub records_skipped: Option<String>,
}

/// How the file ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Eof {
    Bgzf(input::BgzfEnd),
    Cram(bool),
    /// Plain SAM text: whether it ends with a line end.
    Text(bool),
    /// gzip that is not BGZF: no marker to look for.
    Gzip,
}

/// Run every check on `path`.
pub fn run(path: &Path, options: &Options) -> Result<Report> {
    let opened = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let metadata = std::fs::metadata(path).map_err(opened)?;
    if !metadata.is_file() {
        return Err(opened(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "not a regular file (doctor reads the end of the file, which a pipe does not have)",
        )));
    }
    let facts = read_facts(path, options).map_err(opened)?;
    let kind = facts
        .kind
        .map_or_else(String::new, |k| kind_name(k, &facts));
    Ok(Report {
        path: path.to_path_buf(),
        kind,
        checks: checks::judge(path, &facts, options),
    })
}

fn kind_name(kind: Kind, facts: &Facts) -> String {
    match kind {
        Kind::Bam => "BAM".to_owned(),
        Kind::Cram => match facts.cram_version {
            Some((major, minor)) => format!("CRAM {major}.{minor}"),
            None => "CRAM".to_owned(),
        },
        Kind::Sam => "SAM".to_owned(),
        Kind::GzipSam if matches!(facts.eof, Some(Eof::Bgzf(_))) => "SAM (bgzip)".to_owned(),
        Kind::GzipSam => "SAM (gzip)".to_owned(),
    }
}

/// Read what the checks judge. An I/O error opening the file is the only
/// error; a file that opens and holds the wrong bytes gives facts that
/// say so.
fn read_facts(path: &Path, options: &Options) -> std::io::Result<Facts> {
    let mut facts = Facts {
        kind: None,
        kind_error: None,
        short: None,
        cram_version: None,
        eof: None,
        header: None,
        lines: Vec::new(),
        stats: None,
        records_error: None,
        records_skipped: None,
    };
    let size = std::fs::metadata(path)?.len();
    if size < MAGIC_BYTES {
        facts.short = Some(size);
        facts.kind_error = Some(format!(
            "the file holds {}, under the {MAGIC_BYTES}-byte magic of BAM, CRAM or SAM",
            match size {
                1 => "1 byte".to_owned(),
                n => format!("{n} bytes"),
            }
        ));
        return Ok(facts);
    }
    let kind = match input::kind(path)? {
        Ok(kind) => kind,
        Err(e) => {
            facts.kind_error = Some(e);
            return Ok(facts);
        }
    };
    facts.kind = Some(kind);
    facts.eof = Some(match kind {
        Kind::Bam => Eof::Bgzf(input::bgzf_end(path)?),
        Kind::Cram => Eof::Cram(input::cram_has_eof(path)?),
        Kind::Sam => Eof::Text(input::text_ends_with_newline(path)?),
        Kind::GzipSam if input::is_bgzf(path)? => Eof::Bgzf(input::bgzf_end(path)?),
        Kind::GzipSam => Eof::Gzip,
    });
    let header = match kind {
        Kind::Bam => input::bam_header(path).map(|(header, stream)| (header, Some(stream))),
        Kind::Cram => input::cram_header(path).map(|header| (header, None)),
        Kind::Sam => input::sam_header(path, false).map(|header| (header, None)),
        Kind::GzipSam => input::sam_header(path, true).map(|header| (header, None)),
    };
    let (header, stream) = match header {
        Ok(found) => found,
        Err(e) => {
            facts.header = Some(Err(if e.kind() == std::io::ErrorKind::UnexpectedEof {
                "the file ends inside the header".to_owned()
            } else {
                e.to_string()
            }));
            return Ok(facts);
        }
    };
    facts.cram_version = header.cram_version;
    facts.lines = input::lines(&header.text);
    let read_groups: HashSet<String> = facts
        .lines
        .iter()
        .filter(|l| l.kind == "RG")
        .filter_map(|l| l.get("ID").map(str::to_owned))
        .collect();
    let references = checks::reference_names(&facts.lines, &header.binary);
    facts.header = Some(Ok(header));
    if options.records == 0 {
        return Ok(facts);
    }
    // A CRAM decoded against other contigs fails inside the decoder; the
    // `reference` check names the contig instead.
    if kind == Kind::Cram
        && let Some(fasta) = &options.reference
        && let Some(why) = checks::fasta_disagrees(fasta, &facts.lines)
    {
        facts.records_skipped = Some(why);
        return Ok(facts);
    }
    let mut tally = records::Tally::new(&read_groups, &references);
    match stream {
        Some(mut stream) => {
            tally.read_from(&mut stream, options.records);
            facts.stats = Some(tally.stats);
        }
        None => {
            let read = crate::pass::with_record_stream_at(
                path,
                options.reference.as_deref(),
                options.threads.max(1),
                Stringency::Silent,
                |_, stream, _| {
                    tally.read_from(stream, options.records);
                    Ok(())
                },
            );
            match read {
                Ok(()) => facts.stats = Some(tally.stats),
                Err(e) => facts.records_error = Some(e.to_string()),
            }
        }
    }
    Ok(facts)
}

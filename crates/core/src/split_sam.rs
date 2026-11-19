//! The split tools of docs/07 wave W10, `SplitSamByLibrary` and
//! `SplitSamByNumberOfReads`, and the fan-out sink of a chain that runs
//! either on the stream at its position: one stream in, several files out.
//!
//! [`Fanout`] writes the files. Each is an [`AtomicFile`]: written under a
//! temporary name in the output directory, renamed into place by
//! [`Fanout::finish`] only after the last record, with its `.bai` and
//! `.md5` when asked. A run that fails, or a fan-out dropped before its
//! finish, leaves none of its files. Picard holds one open writer per
//! library and may run out of file descriptors; a fan-out holds at most
//! [`MAX_OPEN_FILES`] files open at once. A file closed to make room gets
//! its staged records written as a short BGZF block, and is opened again
//! (append, the same device and inode checked) for its next record: BGZF
//! blocks concatenate, so the file holds the same records. A split writes
//! at most [`MAX_FILES`] files.
//!
//! Records are written as htsjdk's presorted `SAMFileWriter` writes them:
//! the header with `@HD VN:1.6` (and `SO:unsorted` when it has no order),
//! each record with the `bin` htsjdk's codec computes, each file checking
//! the header's order ([`Presorted`]).

use std::collections::HashMap;
use std::io;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::header::record::value::map::read_group::tag::LIBRARY;

use crate::bam_raw::{RawRecord, TagValue};
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::mark_duplicates::SortOrder;
use crate::pass;
use crate::record_edit::{self, Mode, Restore};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text::Stringency;
use crate::stage::{BuildSink, Frame, Products, Sink, Step};

mod fanout;
#[cfg(test)]
mod tests;

pub use fanout::{Claims, Fanout, Format, MAX_FILES, MAX_OPEN_FILES, WriterOptions};

/// `IOUtil.makeFileNameSafe` on a header value: each byte is the Java
/// character of its code (htsjdk reads header text one byte per
/// character), the string is trimmed of characters up to the space, and
/// each whitespace character and each of ``!"#$%&'()*/:;<=>?@[]\^`{|}~``
/// becomes `_`. The name is the result encoded in UTF-8, as Java names a
/// file: the byte 0xE9 gives `é`.
pub fn file_name_safe(value: &[u8]) -> String {
    let start = value.iter().position(|b| *b > b' ').unwrap_or(value.len());
    let end = value
        .iter()
        .rposition(|b| *b > b' ')
        .map_or(start, |i| i + 1);
    value[start..end.max(start)]
        .iter()
        .map(|&byte| match byte {
            b' ' | b'\t' | b'\n' | 0x0B | 0x0C | b'\r' => '_',
            b'!' | b'"' | b'#' | b'$' | b'%' | b'&' | b'\'' | b'(' | b')' | b'*' | b'/' | b':'
            | b';' | b'<' | b'=' | b'>' | b'?' | b'@' | b'[' | b']' | b'\\' | b'^' | b'`'
            | b'{' | b'|' | b'}' | b'~' => '_',
            other => char::from(other),
        })
        .collect()
}

/// Picard's message for a header without a library, before the input.
pub const NO_LIBRARIES: &str = "No individual libraries are specified in the header of ";

/// The `SplitSamByNumberOfReads` options that size the files.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CountOptions {
    /// `SPLIT_TO_N_FILES`, 0 when not given.
    pub files: i32,
    /// `SPLIT_TO_N_READS`, 0 when not given.
    pub reads: i32,
    /// `TOTAL_READS_IN_INPUT`, 0 when not given.
    pub total: i64,
    /// `OUT_PREFIX`.
    pub prefix: String,
}

impl CountOptions {
    /// The options of `inv`, checked as Picard's command line checks them.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let given = |name| inv.get(name).is_some();
        if given("SPLIT_TO_N_FILES") && given("SPLIT_TO_N_READS") {
            return Err(Error::Usage(
                "SPLIT_TO_N_FILES and SPLIT_TO_N_READS cannot be used together".to_owned(),
            ));
        }
        let int = |name| -> Result<i32> {
            match inv.get(name) {
                None | Some("null") => Ok(0),
                Some(text) => crate::java::parse_integer(text, 32)
                    .map(|n| n as i32)
                    .map_err(|e| Error::Usage(format!("{name}={text}: not a valid value: {e}"))),
            }
        };
        let options = Self {
            files: int("SPLIT_TO_N_FILES")?,
            reads: int("SPLIT_TO_N_READS")?,
            total: inv.parse("TOTAL_READS_IN_INPUT", 0i64)?,
            prefix: inv.get("OUT_PREFIX").unwrap_or("shard").to_owned(),
        };
        if options.total < 0 {
            return Err(Error::Usage(format!(
                "Cannot set TOTAL_READS_IN_INPUT to a number less than 1, found {}.",
                options.total
            )));
        }
        if options.files <= 1 && options.reads <= 1 {
            return Err(Error::Usage(format!(
                "One of SPLIT_TO_N_FILES or SPLIT_TO_N_READS must be greater than 0. Found \
                 SPLIT_TO_N_FILES is {} and SPLIT_TO_N_READS is {}.",
                options.files, options.reads
            )));
        }
        Ok(options)
    }

    /// Picard's arithmetic: the number of files, then the records per file,
    /// each `(int) Math.ceil` of a double division.
    pub fn sizes(&self, total: i64) -> (i32, i32) {
        #[allow(clippy::cast_precision_loss, clippy::cast_possible_truncation)]
        let ceil = |a: f64, b: f64| (a / b).ceil() as i32;
        let files = if self.files != 0 {
            self.files
        } else {
            ceil(total as f64, f64::from(self.reads))
        };
        #[allow(clippy::cast_precision_loss)]
        let per_file = ceil(total as f64, f64::from(files));
        (files, per_file)
    }
}

/// How records are sent to files.
enum Route {
    Library(ByLibrary),
    Count(ByCount),
}

/// `SplitSamByLibrary`: a file per library, made before the first record;
/// `unknown` for a record whose read group is missing, not in the header,
/// or without `LB`, made when the first such record comes.
struct ByLibrary {
    /// Read group ID to its file; `None` for a group without `LB`.
    groups: HashMap<Vec<u8>, Option<usize>>,
    unknown: Option<usize>,
    unknown_header: sam::Header,
    dir: PathBuf,
    extension: String,
    format: Format,
}

/// `SplitSamByNumberOfReads`: records in file order, a new file once
/// `per_file` records are written and the name changes.
struct ByCount {
    per_file: i32,
    written: i32,
    index: i32,
    last_name: Vec<u8>,
    current: usize,
    prefix: String,
    dir: PathBuf,
    extension: String,
    format: Format,
    header: sam::Header,
}

impl ByCount {
    fn path(&self) -> PathBuf {
        self.dir.join(format!(
            "{}_{:04}.{}",
            self.prefix, self.index, self.extension
        ))
    }
}

/// What a split sends where.
#[derive(Clone, Debug)]
pub enum SplitKind {
    Library,
    Count(CountOptions),
}

/// One split: its route and its files.
pub struct Splitter {
    route: Route,
    fanout: Fanout,
    record: RawRecord,
}

/// The header each file of a split gets: htsjdk's writer sets `VN:1.6`
/// and `SO:unsorted` when there is no order.
fn written_header(header: sam::Header) -> sam::Header {
    record_edit::with_version(header, true)
}

/// `base` (a header without read groups) with the read groups of `header`
/// at `positions`, in header order.
fn with_groups(base: &sam::Header, header: &sam::Header, positions: &[usize]) -> sam::Header {
    let mut out = base.clone();
    let groups = out.read_groups_mut();
    for &at in positions {
        if let Some((id, map)) = header.read_groups().get_index(at) {
            groups.insert(id.clone(), map.clone());
        }
    }
    out
}

/// What the files of a split are named and hold.
pub struct Layout<'a> {
    pub dir: &'a Path,
    /// `bam`, `sam` or `cram`: the input's kind.
    pub extension: &'a str,
    pub format: Format,
    /// The input, for the error of a header without libraries.
    pub input: &'a str,
    /// The total of a count split, known before the first record.
    pub total: i64,
}

impl Splitter {
    /// The split of a stream with `header` into `layout`'s files.
    pub fn new(
        kind: &SplitKind,
        header: &sam::Header,
        layout: &Layout,
        fanout: Fanout,
    ) -> Result<Self> {
        let mut splitter = Self {
            route: Route::Library(ByLibrary {
                groups: HashMap::new(),
                unknown: None,
                unknown_header: sam::Header::default(),
                dir: layout.dir.to_path_buf(),
                extension: layout.extension.to_owned(),
                format: layout.format,
            }),
            fanout,
            record: RawRecord::default(),
        };
        match kind {
            SplitKind::Library => splitter.libraries(header, layout)?,
            SplitKind::Count(options) => splitter.count(options, header, layout)?,
        }
        Ok(splitter)
    }

    fn libraries(&mut self, header: &sam::Header, layout: &Layout) -> Result<()> {
        // Libraries in header order, each with the positions of its read
        // groups; one hash lookup per read group.
        let mut libraries: Vec<(Vec<u8>, Vec<usize>)> = Vec::new();
        let mut library_at: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut without_library = Vec::new();
        for (at, (_, map)) in header.read_groups().iter().enumerate() {
            let Some(library) = map.other_fields().get(&LIBRARY) else {
                without_library.push(at);
                continue;
            };
            let next = libraries.len();
            let index = *library_at.entry(library.to_vec()).or_insert(next);
            if index == next {
                libraries.push((library.to_vec(), Vec::new()));
            }
            libraries[index].1.push(at);
        }
        let mut groups: HashMap<Vec<u8>, Option<usize>> = HashMap::new();
        for &at in &without_library {
            if let Some((id, _)) = header.read_groups().get_index(at) {
                groups.insert(id.to_vec(), None);
            }
        }
        if libraries.is_empty() {
            return Err(Error::Usage(format!("{NO_LIBRARIES}{}", layout.input)));
        }
        let mut base = header.clone();
        base.read_groups_mut().clear();
        for (library, positions) in &libraries {
            let only = with_groups(&base, header, positions);
            let name = format!("{}.{}", file_name_safe(library), layout.extension);
            let what = format!("library {:?}", String::from_utf8_lossy(library));
            let index = self.fanout.create(
                &layout.dir.join(name),
                &what,
                &written_header(only),
                layout.format,
            )?;
            for &at in positions {
                if let Some((id, _)) = header.read_groups().get_index(at) {
                    groups.insert(id.to_vec(), Some(index));
                }
            }
        }
        self.route = Route::Library(ByLibrary {
            groups,
            unknown: None,
            unknown_header: written_header(with_groups(&base, header, &without_library)),
            dir: layout.dir.to_path_buf(),
            extension: layout.extension.to_owned(),
            format: layout.format,
        });
        Ok(())
    }

    fn count(
        &mut self,
        options: &CountOptions,
        header: &sam::Header,
        layout: &Layout,
    ) -> Result<()> {
        // A count above MAX_FILES is taken: the files written are at most
        // the records, and `Fanout::create` refuses the one past the bound.
        let (_, per_file) = options.sizes(layout.total);
        let mut route = ByCount {
            per_file,
            written: 0,
            index: 1,
            last_name: Vec::new(),
            current: 0,
            prefix: options.prefix.clone(),
            dir: layout.dir.to_path_buf(),
            extension: layout.extension.to_owned(),
            format: layout.format,
            header: written_header(header.clone()),
        };
        route.current = self.fanout.create(
            &route.path(),
            &format!("file {}", route.index),
            &route.header,
            route.format,
        )?;
        self.route = Route::Count(route);
        Ok(())
    }

    /// Send one record to its file, with the `bin` htsjdk's writer gives it.
    pub fn push(&mut self, record: &RawRecord) -> Result<()> {
        self.record.clone_from(record);
        self.record.set_codec_bin();
        let record = &self.record;
        let fanout = &mut self.fanout;
        let index = match &mut self.route {
            Route::Library(route) => {
                let group = match record.tag(*b"RG").map_err(Error::Record)? {
                    None => None,
                    Some(TagValue::String(id)) => route.groups.get(&id).copied().flatten(),
                    Some(_) => {
                        return Err(Error::Record(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "record {}: its RG tag is not a string",
                                String::from_utf8_lossy(record.name())
                            ),
                        )));
                    }
                };
                match (group, route.unknown) {
                    (Some(index), _) | (None, Some(index)) => index,
                    (None, None) => {
                        let name = format!("unknown.{}", route.extension);
                        let index = fanout.create(
                            &route.dir.join(name),
                            "the records without a library",
                            &route.unknown_header,
                            route.format,
                        )?;
                        route.unknown = Some(index);
                        index
                    }
                }
            }
            Route::Count(route) => {
                if route.written >= route.per_file && route.last_name != record.name() {
                    fanout.end(route.current)?;
                    route.index = route.index.wrapping_add(1);
                    route.current = fanout.create(
                        &route.path(),
                        &format!("file {}", route.index),
                        &route.header,
                        route.format,
                    )?;
                    route.written = 0;
                }
                route.last_name.clear();
                route.last_name.extend_from_slice(record.name());
                route.written = route.written.wrapping_add(1);
                route.current
            }
        };
        fanout.write(index, record)
    }

    /// A count split's files so far, records per file and records in the
    /// last file: the numbers of Picard's warning when
    /// `TOTAL_READS_IN_INPUT` is not the count read.
    pub fn count_state(&self) -> Option<(i32, i32, i32)> {
        match &self.route {
            Route::Count(route) => Some((route.index, route.per_file, route.written)),
            Route::Library(_) => None,
        }
    }

    /// End every file and rename all of them into place.
    pub fn finish(self) -> Result<Vec<PathBuf>> {
        self.fanout.finish()
    }
}

/// `OUTPUT` of a split: an existing directory, as htsjdk's
/// `IOUtil.assertDirectoryIsWritable` asks.
pub fn check_directory(dir: &Path) -> Result<()> {
    let refuse = |message: &str| Error::Output {
        path: dir.to_path_buf(),
        source: io::Error::other(message.to_owned()),
    };
    match std::fs::metadata(dir) {
        Err(_) => Err(refuse("Directory does not exist")),
        Ok(meta) if !meta.is_dir() => Err(refuse("not a directory")),
        Ok(_) => Ok(()),
    }
}

/// The directory `SplitSamByLibrary` writes into without `OUTPUT`, as
/// Picard's default: the current directory, whatever the input.
pub fn current_directory(_input: &Path) -> PathBuf {
    PathBuf::from(".")
}

/// `OUTPUT` of a split tool, `default` when not given.
pub fn output_dir(inv: &Invocation, default: Option<&str>) -> Result<PathBuf> {
    inv.get("OUTPUT")
        .or(default)
        .map(PathBuf::from)
        .ok_or_else(|| Error::Usage("OUTPUT is required".to_owned()))
}

/// The options of a split run on its own.
pub struct Run<'a> {
    pub inv: &'a Invocation,
    pub kind: SplitKind,
    pub dir: PathBuf,
    /// The format of each kind of input, by extension: `SplitSamByLibrary`
    /// writes the input's kind (a CRAM is refused, docs/07 wave W15);
    /// `SplitSamByNumberOfReads` writes SAM text for SAM, else a BAM.
    pub format: fn(&str) -> Result<Format>,
}

/// A split tool on its own: read `INPUT` once (twice for a count split
/// without `TOTAL_READS_IN_INPUT`), write the files. Returns their paths.
pub fn run(run: &Run) -> Result<Vec<PathBuf>> {
    let inv = run.inv;
    let input = &inv.input;
    let reference = inv.reference();
    let stringency = inv.stringency()?;
    let options = WriterOptions::from_invocation(inv)?;
    let mut total = match &run.kind {
        SplitKind::Count(count) => count.total,
        SplitKind::Library => 0,
    };
    let stream = !std::fs::metadata(input).is_ok_and(|m| m.is_file()) || input.starts_with("/dev");
    if matches!(run.kind, SplitKind::Count(_)) && total == 0 && stream {
        return Err(Error::Input {
            path: input.clone(),
            source: io::Error::other(
                "INPUT is not a regular file. If TOTAL_READS_IN_INPUT is not supplied, INPUT \
                 cannot be a stream.",
            ),
        });
    }
    check_directory(&run.dir)?;
    if let SplitKind::Count(_) = &run.kind
        && total == 0
    {
        total = count_records(input, reference, inv.threads, stringency)?;
    }
    let mut reads = vec![input.clone()];
    if let Some(reference) = reference {
        reads.extend(crate::reference::files_read(reference));
    }
    pass::with_record_stream_kind(
        input,
        reference,
        inv.threads,
        stringency,
        |header, reader, kind| {
            record_edit::check_sort_order(header, stringency, input)?;
            let extension = match kind {
                Kind::Bam => "bam",
                Kind::Cram => "cram",
                Kind::Sam | Kind::GzipSam => "sam",
            };
            let format = (run.format)(extension)?;
            if matches!(run.kind, SplitKind::Count(_)) {
                warn_count(inv, header);
            }
            let cram = kind == Kind::Cram;
            let mut restore = Restore::new(cram, Mode::Rewrite, reference, header)?;
            let fanout = Fanout::new(options, header, reads.clone(), None)?;
            let absolute = std::path::absolute(input).unwrap_or_else(|_| input.clone());
            let layout = Layout {
                dir: &run.dir,
                extension,
                format,
                input: &absolute.display().to_string(),
                total,
            };
            let mut splitter = Splitter::new(&run.kind, header, &layout, fanout)?;
            let mut record = RawRecord::default();
            let mut seen: i64 = 0;
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(input, cram, e))?
            {
                restore.apply(&mut record)?;
                splitter.push(&record)?;
                seen += 1;
            }
            let counted = splitter.count_state();
            let files = splitter.finish()?;
            if let (Some((generated, per_file, last)), false) = (counted, inv.quiet)
                && seen != total
            {
                eprintln!(
                    "bammetrics: The totalReads ({total}) provided does not match the reads found \
                 in the input file ({seen}). Files may not be split evenly or number of files \
                 may not match what was requested. There were {generated} files generated each \
                 with around {per_file} reads except the last file which contained {last} \
                 reads."
                );
            }
            Ok(files)
        },
    )
}

/// The warnings Picard's `SplitSamByNumberOfReads` logs about its input.
fn warn_count(inv: &Invocation, header: &sam::Header) {
    if inv.quiet {
        return;
    }
    if SortOrder::of(header) == SortOrder::Coordinate {
        eprintln!(
            "bammetrics: Splitting a coordinate sorted bam may result in invalid bams that do \
             not always contain each read's mate in the same bam."
        );
    }
    let version = header.header().map(|hd| hd.version());
    let shown = version.map_or_else(
        || "1.6".to_owned(),
        |v| format!("{}.{}", v.major(), v.minor()),
    );
    if shown != "1.6" {
        eprintln!(
            "bammetrics: Input file's version is {shown}, but the current SAM format version \
             is 1.6. Outputs will be written with current version."
        );
    }
}

/// The records of `input`: the first read of a count split without
/// `TOTAL_READS_IN_INPUT`.
fn count_records(
    input: &Path,
    reference: Option<&Path>,
    threads: usize,
    stringency: Stringency,
) -> Result<i64> {
    pass::with_record_stream_kind(input, reference, threads, stringency, |_, reader, kind| {
        let mut record = RawRecord::default();
        let mut count = 0i64;
        while record
            .read_from(reader)
            .map_err(|e| pass::read_error(input, kind == Kind::Cram, e))?
        {
            count += 1;
        }
        Ok(count)
    })
}

/// The options of a split as a step of a chain.
pub struct ChainSplit {
    pub kind: SplitKind,
    pub dir: PathBuf,
    pub writer: WriterOptions,
    /// Every file the chain reads.
    pub reads: Vec<PathBuf>,
    /// The step's label and the chain's claims.
    pub claims: Option<(String, Claims)>,
}

/// The options of split tool `tool` as a chain step, checked before the
/// first read: a count split needs `TOTAL_READS_IN_INPUT`, since a chain
/// reads its stream once (Picard refuses a stream without it).
pub fn chain_options(tool: &Tool, inv: &Invocation) -> Result<(SplitKind, PathBuf)> {
    match tool.name {
        "SplitSamByLibrary" => Ok((SplitKind::Library, output_dir(inv, Some("."))?)),
        "SplitSamByNumberOfReads" => {
            let count = CountOptions::from_invocation(inv)?;
            if count.total == 0 {
                return Err(Error::Usage(
                    "TOTAL_READS_IN_INPUT is required in a chain: the stream is read once, and \
                     Picard needs it for a stream"
                        .to_owned(),
                ));
            }
            Ok((SplitKind::Count(count), output_dir(inv, None)?))
        }
        other => Err(Error::Usage(format!("{other} is no split"))),
    }
}

/// True for a tool that is a split: in a chain, a sink that writes a file
/// per library or per run of records.
pub fn is_split(tool: &Tool) -> bool {
    matches!(tool.name, "SplitSamByLibrary" | "SplitSamByNumberOfReads")
}

/// The split as a chain step: a sink on the records of its position. Its
/// files are BAMs named `.bam`, whatever the chain's source.
pub fn chain_step(label: &str, split: ChainSplit) -> Step {
    let build: BuildSink = Box::new(move |frame: &Frame| {
        let ChainSplit {
            kind,
            dir,
            writer,
            mut reads,
            claims,
        } = split;
        check_directory(&dir)?;
        reads.extend(frame.reads.iter().cloned());
        let total = match &kind {
            SplitKind::Count(count) => count.total,
            SplitKind::Library => 0,
        };
        let fanout = Fanout::new(writer, &frame.header, reads, claims)?;
        let layout = Layout {
            dir: &dir,
            extension: "bam",
            format: Format::Bam,
            input: &frame.label,
            total,
        };
        let splitter = Splitter::new(&kind, &frame.header, &layout, fanout)?;
        Ok(Box::new(SplitSink { splitter }) as Box<dyn Sink>)
    });
    Step::sink(label, build)
}

/// A split on the thread of its chain step.
struct SplitSink {
    splitter: Splitter,
}

impl Sink for SplitSink {
    fn accept_batch(&mut self, batch: &[RawRecord]) -> Result<()> {
        batch
            .iter()
            .try_for_each(|record| self.splitter.push(record))
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        Ok(Products {
            files: self.splitter.finish()?,
            ..Products::default()
        })
    }
}

/// `COMPRESSION_LEVEL` of the files a split writes.
pub const COMPRESSION_OPTION: Opt = Opt::int(
    &["COMPRESSION_LEVEL"],
    "BGZF compression level of each BAM written, 0 to 9 (default 5).",
);
/// `CREATE_INDEX` of the files a split writes.
pub const INDEX_OPTION: Opt = Opt::bool(
    &["CREATE_INDEX"],
    "Write a .bai beside each coordinate-sorted BAM written (default false).",
);
/// `CREATE_MD5_FILE` of the files a split writes.
pub const MD5_OPTION: Opt = Opt::bool(
    &["CREATE_MD5_FILE"],
    "Write FILE.md5 beside each file written (default false).",
);
/// `MAX_RECORDS_IN_RAM`, which Picard takes and a split never reads.
pub const RECORDS_IN_RAM_OPTION: Opt = Opt::ignored_int(
    &["MAX_RECORDS_IN_RAM"],
    "records are written as they stream",
);

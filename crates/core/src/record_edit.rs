//! The shape of a Picard tool that edits a SAM, BAM or CRAM header and its
//! records one by one, then writes them (docs/07 waves W09 and W10:
//! `AddCommentsToBam`, `ReplaceSamHeader`, `AddOATag`, `CleanSam`,
//! `SetNmMdAndUqTags`, `AddOrReplaceReadGroups`, ...).
//!
//! A tool supplies a [`RecordEdit`]: the output header from the input
//! header, then each record edited in place, kept or dropped, and records
//! to add after the last one. It is built once per run from a [`Setup`]:
//! the input header, the reference, whether the records came from a CRAM,
//! and the [`Mode`]. This module runs it two ways:
//!
//! * [`run`], the tool on its own: one read of `INPUT` (BAM, CRAM with
//!   `REFERENCE_SEQUENCE`, SAM text at `VALIDATION_STRINGENCY`), then a BAM
//!   through [`BamWriter`] (`COMPRESSION_LEVEL`, `CREATE_INDEX`,
//!   `CREATE_MD5_FILE`), or SAM text through [`SamWriter`] for an `OUTPUT`
//!   ending in `.sam`. A `.cram` `OUTPUT` is refused with the wave that adds
//!   it. Every file is written under a temporary name and renamed into place
//!   when complete ([`crate::cleanup::AtomicFile`]), so a failed run leaves
//!   none. An `OUTPUT`, `.bai` or `.md5` that is an input (by device and
//!   inode, or by resolved path) is refused before the first read.
//! * [`stage`], the same edit as a [`Transform`] of a chain: each batch is
//!   edited where it stands, the stream goes on in the same read, and the
//!   output header is handed down through [`Transform::header`]. The
//!   options of the BAM file belong to the `write` step ([`check_stage`]).
//!
//! Two ways of writing, as Picard's tools write ([`Mode`]):
//! * [`Mode::BlockCopy`]: htsjdk's `BamFileIoUtils.reheaderBamFile`
//!   (`AddCommentsToBam`, `ReplaceSamHeader` on a `.bam`). The input must be
//!   a BAM; its record bytes go out unchanged, in the order read, with no
//!   order check; the index is `<OUTPUT>.bai` (`out.bam.bai`). Picard writes
//!   BAM bytes whatever the `OUTPUT` name, so a `.sam` name is refused.
//! * [`Mode::Rewrite`]: a `SamReader` into a presorted `SAMFileWriter`.
//!   Records of a CRAM get `MD` and `NM` as htsjdk decodes them, every
//!   record gets the `bin` htsjdk's `BAMRecordCodec` computes, the writer
//!   checks the header's order, and the index is `out.bai`.
//!
//! No `@PG` line is added: none of the six W09 tools adds one in Picard
//! 3.5.0 (none calls `getPGRecord`). A tool that adds one does it in
//! [`RecordEdit::header`].

use std::io;
use std::path::{Path, PathBuf};

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::header::Version;

use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, BamWriter, Companions};
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::mark_duplicates::{ContigCache, restore_md_and_nm, same_file};
use crate::pass;
use crate::reference::{Reference, files_read};
use crate::registry::{Invocation, Opt};
use crate::sam_text::{self, SamWriter, Stringency};
use crate::spill::TempDir;
use crate::stage::{Counters, Frame, Stage, Transform, WriteFormat, later_kind};

mod reorder;
pub use reorder::{
    DEFAULT_RECORDS_IN_RAM, Gate, Regroup, Reorder, WriterOrder, presorted_stage, sort_limits,
    sorting_stage, stage_limits,
};

/// What [`RecordEdit::record`] did to one record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// Handed on as it came.
    Same,
    /// Handed on with a change.
    Changed,
    /// Left out of the output.
    Dropped,
}

/// The edit of one record-editing tool, built once per run ([`Setup`]).
pub trait RecordEdit {
    /// The output header, from the header of the records read.
    fn header(&self, input: &sam::Header) -> Result<sam::Header>;

    /// Edit `record` in place; the records come in file order.
    fn record(&mut self, record: &mut RawRecord) -> Result<Outcome>;

    /// After the last record: records to write after it.
    fn finish(&mut self, _out: &mut Vec<RawRecord>) -> Result<()> {
        Ok(())
    }

    /// True for a tool whose stand-alone run is a [`Mode::Rewrite`]: as a
    /// chain step it gives each record, and each record of
    /// [`RecordEdit::finish`], the `bin` htsjdk's writer computes, so the
    /// stream after it holds the records the tool writes.
    fn rewrites(&self) -> bool {
        false
    }
}

/// How the records are read and written; Picard's tools use one of two
/// htsjdk paths, and a chain step a third.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    /// `BamFileIoUtils.reheaderBamFile`: a new header, then the record
    /// bytes of a BAM as they are.
    BlockCopy,
    /// A `SamReader` into a presorted `SAMFileWriter`.
    Rewrite,
    /// A step of a chain ([`stage`]): the records of the stream at its
    /// position, as BAM records.
    Chain,
}

/// What a [`RecordEdit`] is built from.
#[derive(Clone, Copy, Debug)]
pub struct Setup<'a> {
    /// The header of the records read.
    pub header: &'a sam::Header,
    /// `REFERENCE_SEQUENCE`.
    pub reference: Option<&'a Path>,
    /// The records were decoded from a CRAM.
    pub cram: bool,
    pub mode: Mode,
}

/// `OUTPUT` of a record-editing tool.
pub const OUTPUT_OPTION: Opt = Opt::path(
    &["OUTPUT", "O"],
    "The BAM to write, or SAM text for a name ending in .sam.",
)
.written()
.required();
/// `COMPRESSION_LEVEL` of the BAM written.
pub const COMPRESSION_OPTION: Opt = Opt::int(
    &["COMPRESSION_LEVEL"],
    "BGZF compression level of OUTPUT, 0 to 9 (default 5).",
);
/// `CREATE_INDEX` of the BAM written.
pub const INDEX_OPTION: Opt = Opt::bool(
    &["CREATE_INDEX"],
    "Write a .bai beside a coordinate-sorted BAM OUTPUT (default false).",
);
/// `CREATE_MD5_FILE` of the file written.
pub const MD5_OPTION: Opt = Opt::bool(&["CREATE_MD5_FILE"], "Write OUTPUT.md5 (default false).");
/// `MAX_RECORDS_IN_RAM`, which Picard takes and a record stream never reads.
pub const RECORDS_IN_RAM_OPTION: Opt =
    Opt::ignored_int(&["MAX_RECORDS_IN_RAM"], "records are edited as they stream");

/// The options of the file a stand-alone run writes, which a chain step
/// leaves to its `write` step.
pub const WRITER_OPTIONS: [&str; 4] = [
    "OUTPUT",
    "COMPRESSION_LEVEL",
    "CREATE_INDEX",
    "CREATE_MD5_FILE",
];

/// The options of a stand-alone run.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: PathBuf,
    pub output: PathBuf,
    pub reference: Option<PathBuf>,
    /// Files beside `INPUT` the tool reads (`HEADER`): no output may be one.
    pub reads: Vec<PathBuf>,
    pub compression_level: u8,
    pub create_index: bool,
    pub create_md5: bool,
    pub threads: usize,
    /// How SAM text input is checked.
    pub stringency: Stringency,
    pub mode: Mode,
    /// Sorts between the edit and the file ([`Reorder`]); none by default.
    pub reorder: Reorder,
}

impl Options {
    /// The common options of `inv`, in [`Mode::Rewrite`]. A tool sets
    /// `mode`, `reads` and `stringency` afterwards when it differs.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let output = inv
            .get("OUTPUT")
            .map(PathBuf::from)
            .ok_or_else(|| Error::Usage("OUTPUT is required".to_owned()))?;
        let level: i64 = inv.parse("COMPRESSION_LEVEL", 5)?;
        let compression_level = u8::try_from(level)
            .ok()
            .filter(|l| *l <= 9)
            .ok_or_else(|| Error::Usage(format!("COMPRESSION_LEVEL={level}: expected 0 to 9")))?;
        Ok(Self {
            input: inv.input.clone(),
            output,
            reference: inv.reference().map(Path::to_path_buf),
            reads: Vec::new(),
            compression_level,
            create_index: inv.flag("CREATE_INDEX", false)?,
            create_md5: inv.flag("CREATE_MD5_FILE", false)?,
            threads: inv.threads,
            stringency: inv.stringency()?,
            mode: Mode::Rewrite,
            reorder: Reorder::default(),
        })
    }

    /// The format written: BAM, or SAM text for a `.sam` name. A block
    /// copy writes a BAM.
    fn format(&self) -> WriteFormat {
        match self.mode {
            Mode::BlockCopy => WriteFormat::Bam,
            _ => WriteFormat::of(&self.output),
        }
    }

    /// Where `CREATE_INDEX` writes: `out.bam.bai` for a block copy, as
    /// htsjdk's `StreamInflatingIndexingOutputStream`, else the name
    /// [`bam_writer::index_path`] gives (`out.bai`, `out.dat.bai`).
    pub fn index_path(&self) -> PathBuf {
        match self.mode {
            Mode::BlockCopy => {
                let mut name = self.output.as_os_str().to_owned();
                name.push(".bai");
                PathBuf::from(name)
            }
            _ => bam_writer::index_path(&self.output),
        }
    }

    /// Every file the run may write.
    pub fn written(&self) -> Vec<PathBuf> {
        let mut files = vec![self.output.clone()];
        if self.create_index && self.format() == WriteFormat::Bam {
            files.push(self.index_path());
        }
        if self.create_md5 {
            files.push(bam_writer::md5_path(&self.output));
        }
        files
    }
}

/// Refuse an `OUTPUT` the run cannot write, or one that would destroy a
/// file it reads, before anything is read.
fn check_output(opts: &Options) -> Result<()> {
    if let Some((kind, wave)) = later_kind(&opts.output) {
        return Err(Error::Usage(format!(
            "OUTPUT={}: bammetrics writes BAM and SAM only: {kind} output is wave {wave} of \
             docs/07_full_parity_program.md",
            opts.output.display()
        )));
    }
    if opts.mode == Mode::BlockCopy && WriteFormat::of(&opts.output) == WriteFormat::Sam {
        return Err(Error::Usage(format!(
            "OUTPUT={}: the records of a BAM are copied as BAM blocks, and Picard writes BAM \
             bytes under this name: name the OUTPUT .bam",
            opts.output.display()
        )));
    }
    let written = opts.written();
    let named = std::iter::once(&opts.input)
        .chain(&opts.reference)
        .chain(&opts.reads);
    for input in named {
        if !input.exists() {
            return Err(Error::Input {
                path: input.clone(),
                source: io::ErrorKind::NotFound.into(),
            });
        }
    }
    let read = std::iter::once(opts.input.clone())
        .chain(opts.reference.iter().flat_map(|r| files_read(r)))
        .chain(opts.reads.iter().cloned());
    for input in &read.collect::<Vec<_>>() {
        if let Some(path) = written.iter().find(|path| same_file(input, path)) {
            return Err(Error::Usage(format!(
                "OUTPUT={}: it is the input {}; writing it would destroy it",
                path.display(),
                input.display()
            )));
        }
    }
    Ok(())
}

/// Where OUTPUT goes, told apart as htsjdk's `SAMFileWriterFactory` tells
/// it: by whether the path is a regular file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Target {
    File,
    /// A pipe, a terminal or another device: htsjdk writes no index and no
    /// MD5 beside it and says so in a warning.
    Stream,
    /// A path under `/dev` that is a regular file (`/dev/stdout` redirected
    /// to a file): htsjdk opens `/dev/stdout.bai` or `.md5` and fails.
    Redirected,
}

impl Target {
    fn of(path: &Path) -> Self {
        match std::fs::metadata(path) {
            Ok(meta) if !meta.is_file() => Self::Stream,
            Ok(_) if path.starts_with("/dev") => Self::Redirected,
            _ => Self::File,
        }
    }
}

/// The index and MD5 a run writes beside OUTPUT. A stream gets neither,
/// with htsjdk's warning for each one asked, as Picard writes the BAM to a
/// pipe and exits 0. Standard output redirected to a file is a usage error
/// before the first record is written: Picard fails opening
/// `/dev/stdout.bai` or `/dev/stdout.md5`. The index is asked only of a
/// coordinate-sorted BAM, as [`bam_writer::finish_with_index`] writes it.
fn companions(opts: &Options, header: &sam::Header) -> Result<Companions> {
    use crate::header_merge::SortOrder;
    let index = opts.create_index
        && opts.format() == WriteFormat::Bam
        && SortOrder::of(header) == SortOrder::Coordinate;
    let asked: Vec<&str> = [(index, "index"), (opts.create_md5, "MD5 file")]
        .into_iter()
        .filter_map(|(asked, what)| asked.then_some(what))
        .collect();
    match Target::of(&opts.output) {
        Target::Stream => {
            for what in asked {
                eprintln!(
                    "bammetrics: Cannot create {what} for BAM because output file is not a \
                     regular file: {}",
                    opts.output.display()
                );
            }
            Ok(Companions {
                index: false,
                md5: false,
            })
        }
        Target::Redirected if !asked.is_empty() => Err(Error::Usage(format!(
            "OUTPUT={}: a regular file under /dev, where no {} can be made beside it; set \
             CREATE_INDEX and CREATE_MD5_FILE false (--no-create-index, --no-create-md5-file)",
            opts.output.display(),
            asked.join(" or ")
        ))),
        Target::File | Target::Redirected => Ok(Companions {
            index,
            md5: opts.create_md5,
        }),
    }
}

/// htsjdk `SAMRecord.setHeaderStrict`, which a writer of decoded records
/// calls on each one: a reference or mate reference index that is neither
/// -1 nor one of the header's `references` sequences fails, whatever
/// `VALIDATION_STRINGENCY`. A block copy writes the bytes and checks none.
fn check_reference_ids(record: &RawRecord, references: usize) -> Result<()> {
    for index in [record.reference_id(), record.mate_reference_id()] {
        if index != -1
            && usize::try_from(index)
                .ok()
                .is_none_or(|at| at >= references)
        {
            return Err(Error::Metric(format!(
                "record {}: Reference name for '{index}' not found in sequence dictionary.",
                String::from_utf8_lossy(record.name())
            )));
        }
    }
    Ok(())
}

/// The error for a block copy of a file that is no BAM.
fn not_a_bam(path: &Path, kind: Kind) -> Error {
    let what = match kind {
        Kind::Cram => "a CRAM",
        Kind::Sam | Kind::GzipSam => "SAM text",
        Kind::Bam => "a BAM",
    };
    Error::Input {
        path: path.to_path_buf(),
        source: io::Error::new(
            io::ErrorKind::InvalidData,
            format!("{what}: this tool copies the records of a BAM, as Picard does"),
        ),
    }
}

/// `MD` and `NM` of CRAM records, rebuilt as htsjdk decodes them.
struct Restore(Option<ContigCache>);

impl Restore {
    fn new(cram: bool, mode: Mode, reference: Option<&Path>, header: &sam::Header) -> Result<Self> {
        Ok(Self(match (cram, mode, reference) {
            (true, Mode::Rewrite | Mode::Chain, Some(fasta)) => {
                Some(ContigCache::new(Reference::open(fasta, header)?))
            }
            _ => None,
        }))
    }

    fn apply(&mut self, record: &mut RawRecord) -> Result<()> {
        let Some(contigs) = self.0.as_mut() else {
            return Ok(());
        };
        let mut tags = record.tags().map_err(Error::Record)?;
        restore_md_and_nm(record, &mut tags, contigs)?;
        record.set_tags(&tags).map_err(Error::Record)
    }
}

/// Where the records go: a BAM or SAM text, with the companions it gets.
enum Out {
    Bam(BamWriter, Companions),
    Sam(SamWriter, Companions),
}

impl Out {
    fn create(opts: &Options, header: &sam::Header) -> Result<Self> {
        let companions = companions(opts, header)?;
        Ok(match opts.format() {
            WriteFormat::Sam => {
                let text = sam_text::header_text(header, false);
                Self::Sam(SamWriter::create(&opts.output, header, &text)?, companions)
            }
            WriteFormat::Bam => {
                let writer =
                    BamWriter::create(&opts.output, opts.compression_level, opts.threads, header)?;
                let writer = match opts.mode {
                    Mode::BlockCopy => writer.unchecked(),
                    _ => writer,
                };
                Self::Bam(writer, companions)
            }
        })
    }

    /// Write `record`. A rewrite gives it the `bin` htsjdk's writer computes;
    /// a block copy keeps its bytes.
    fn write(&mut self, record: &mut RawRecord, mode: Mode) -> Result<()> {
        if mode == Mode::Rewrite {
            record.set_codec_bin();
        }
        match self {
            Self::Bam(writer, _) => writer.write(record),
            Self::Sam(writer, _) => writer.write(record),
        }
    }

    /// Close the file, write its companions, rename everything into place.
    fn finish(self, header: &sam::Header, opts: &Options) -> Result<()> {
        match self {
            Self::Bam(writer, companions) => {
                bam_writer::finish_with_index(writer, header, companions, &opts.index_path())
                    .map(|_| ())
            }
            Self::Sam(writer, companions) => {
                let output = writer.close()?;
                let target = output.path().to_path_buf();
                let md5 = companions
                    .md5
                    .then(|| {
                        bam_writer::write_md5(output.writing(), &bam_writer::md5_path(&target))
                    })
                    .transpose()?;
                for file in std::iter::once(output).chain(md5) {
                    let path = file.path().to_path_buf();
                    file.commit()
                        .map_err(|source| Error::Output { path, source })?;
                }
                Ok(())
            }
        }
    }
}

/// The `SO` values htsjdk's `SAMFileHeader.SortOrder` names.
const SORT_ORDERS: [&[u8]; 5] = [
    b"unsorted",
    b"queryname",
    b"coordinate",
    b"duplicate",
    b"unknown",
];

/// `SAMTextHeaderCodec` on the `SO` of the header read: a value
/// `SortOrder` does not name fails at `STRICT` and is reported at
/// `LENIENT`; [`with_version`] writes it `unknown`. SAM text gets the same
/// check as it is read, so its header holds `unknown` by now.
fn check_sort_order(header: &sam::Header, stringency: Stringency, path: &Path) -> Result<()> {
    use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
    let Some(order) = header
        .header()
        .and_then(|hd| hd.other_fields().get(&SORT_ORDER))
        .filter(|order| !SORT_ORDERS.contains(&order.as_slice()))
    else {
        return Ok(());
    };
    let shown = String::from_utf8_lossy(order);
    match stringency {
        Stringency::Strict => Err(Error::Input {
            path: path.to_path_buf(),
            source: io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Found non-conforming header SO tag: {shown}, exiting because \
                     VALIDATION_STRINGENCY=STRICT"
                ),
            ),
        }),
        Stringency::Lenient => {
            eprintln!("Found non-conforming header SO tag: {shown}. Treating as 'unknown'.");
            Ok(())
        }
        Stringency::Silent => Ok(()),
    }
}

/// The header as htsjdk holds it: a `SAMFileHeader` starts with `@HD
/// VN:1.6`, so a file without an `@HD` line is written with that one.
/// `rewrite`: htsjdk's `SAMFileWriter` also sets `VN:1.6` over the
/// header's own version (`VN:1.4` is written `VN:1.6`) and `SO:unsorted`
/// after the other fields of an `@HD` without `SO`; a block copy keeps the
/// line as it is. An `SO` value htsjdk's `SAMFileHeader.SortOrder` does not
/// name becomes `unknown`, as `SAMTextHeaderCodec` reads it below `STRICT`.
fn with_version(mut header: sam::Header, rewrite: bool) -> sam::Header {
    use noodles::sam::header::record::value::map::header::tag::SORT_ORDER;
    let current = Version::new(1, 6);
    let hd = header
        .header_mut()
        .get_or_insert_with(|| Map::<sam::header::record::value::map::Header>::new(current));
    if let Some(order) = hd.other_fields_mut().get_mut(&SORT_ORDER)
        && !SORT_ORDERS.contains(&order.as_slice())
    {
        *order = "unknown".into();
    }
    if rewrite {
        *hd.version_mut() = current;
        hd.other_fields_mut()
            .entry(SORT_ORDER)
            .or_insert_with(|| "unsorted".into());
    }
    header
}

/// The bytes htsjdk writes for a Java string from the command line (a
/// `COMMENT`, a read group field): the low byte of each UTF-16 unit, so
/// `é` is the byte 0xE9.
pub fn java_text(text: &str) -> Vec<u8> {
    text.encode_utf16().map(|unit| unit as u8).collect()
}

/// Count one record's outcome; true when it is written.
fn tally(counters: &mut Counters, outcome: Outcome) -> bool {
    counters.records_in += 1;
    if outcome == Outcome::Dropped {
        return false;
    }
    counters.records_out += 1;
    if outcome == Outcome::Changed {
        counters.changed += 1;
    }
    true
}

/// The tool on its own: read `opts.input`, edit the header and each record
/// with the edit `build` makes, write `opts.output` (and its index and MD5
/// when asked) under temporary names, renamed when every file is complete.
/// Returns what the edit did.
pub fn run(
    opts: &Options,
    build: impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>>,
) -> Result<Counters> {
    run_gated(opts, build, None).map(Option::unwrap_or_default)
}

/// [`run`] with a [`Gate`] that may skip the run: `None` when it did, and
/// then no file is written. An edit that fails before the gate decides is
/// reported only when the run goes on, as Picard decides from its own read
/// of the first records before it edits any.
pub fn run_gated(
    opts: &Options,
    build: impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>>,
    mut gate: Option<&mut dyn Gate>,
) -> Result<Option<Counters>> {
    check_output(opts)?;
    let temp = TempDir::new(opts.reorder.tmp_dir.as_deref(), "record-edit");
    pass::with_record_stream_kind(
        &opts.input,
        opts.reference.as_deref(),
        opts.threads,
        opts.stringency,
        |header, reader, kind| {
            if opts.mode == Mode::BlockCopy && kind != Kind::Bam {
                return Err(not_a_bam(&opts.input, kind));
            }
            check_sort_order(header, opts.stringency, &opts.input)?;
            let cram = kind == Kind::Cram;
            let mut restore = Restore::new(cram, opts.mode, opts.reference.as_deref(), header)?;
            let mut edit = build(&Setup {
                header,
                reference: opts.reference.as_deref(),
                cram,
                mode: opts.mode,
            })?;
            let out_header = with_version(edit.header(header)?, opts.mode == Mode::Rewrite);
            // A run a gate may skip makes its file once the gate lets it go
            // on, so a skip writes nothing, not even to standard output.
            let out = match gate {
                None => Some(Out::create(opts, &out_header)?),
                Some(_) => None,
            };
            let input_order = crate::header_merge::SortOrder::of(header);
            let mut out = reorder::Pipeline::new(opts, input_order, &out_header, &temp, out);
            let mut counters = Counters::default();
            let mut record = RawRecord::default();
            let references = header.reference_sequences().len();
            // An edit's error while the gate looks on: Picard reports it
            // only when it does not skip.
            let mut held: Option<Error> = None;
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&opts.input, cram, e))?
            {
                if let Some(look) = gate.as_deref_mut() {
                    match look.look(&record)? {
                        Some(true) => return Ok(None),
                        Some(false) => {
                            gate = None;
                            if let Some(error) = held.take() {
                                return Err(error);
                            }
                            out.open()?;
                        }
                        None => {}
                    }
                }
                if held.is_some() {
                    continue;
                }
                let mut step = || -> Result<()> {
                    if opts.mode == Mode::Rewrite {
                        check_reference_ids(&record, references)?;
                    }
                    restore.apply(&mut record)?;
                    if tally(&mut counters, edit.record(&mut record)?) {
                        out.push(&mut record)?;
                    }
                    Ok(())
                };
                match step() {
                    Err(error) if gate.is_some() => held = Some(error),
                    other => other?,
                }
            }
            if let Some(look) = gate
                && look.at_end()
            {
                return Ok(None);
            }
            if let Some(error) = held {
                return Err(error);
            }
            let mut tail = Vec::new();
            edit.finish(&mut tail)?;
            for record in &mut tail {
                counters.records_out += 1;
                out.push(record)?;
            }
            out.finish()?.finish(&out_header, opts)?;
            Ok(Some(counters))
        },
    )
}

/// A usage error for an option of the file a stand-alone run writes: in a
/// chain the `write` step writes the records.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    match WRITER_OPTIONS.iter().find(|name| inv.get(name).is_some()) {
        Some(name) => Err(Error::Usage(format!(
            "{name} belongs to the write step in a chain; this step hands its records to the \
             next step"
        ))),
        None => Ok(()),
    }
}

/// The edit `build` makes as a step of a chain: a [`Transform`] on the
/// stream at its position, which is handed down with the edited header.
/// Records of a CRAM get `MD` and `NM` when the chain has a reference, as
/// the `SortSam` step gives them.
pub fn stage<'t>(
    inv: &Invocation,
    frame: &Frame,
    build: impl FnOnce(&Setup) -> Result<Box<dyn RecordEdit>>,
) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let reference = inv.reference();
    let restore = Restore::new(frame.cram, Mode::Chain, reference, &frame.header)?;
    let edit = build(&Setup {
        header: &frame.header,
        reference,
        cram: frame.cram,
        mode: Mode::Chain,
    })?;
    Ok(Stage::Transform(Box::new(Editing {
        edit,
        restore,
        counters: Counters::default(),
        references: frame.header.reference_sequences().len(),
    })))
}

/// The transform of [`stage`].
struct Editing {
    edit: Box<dyn RecordEdit>,
    restore: Restore,
    counters: Counters,
    /// The sequences of the input dictionary, for [`check_reference_ids`].
    references: usize,
}

impl Transform for Editing {
    fn header(&self, up: &sam::Header) -> Result<sam::Header> {
        let rewrites = self.edit.rewrites();
        self.edit
            .header(up)
            .map(|header| with_version(header, rewrites))
    }

    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        let rewrites = self.edit.rewrites();
        let mut kept = 0;
        for at in 0..batch.len() {
            let record = &mut batch[at];
            if rewrites {
                check_reference_ids(record, self.references)?;
            }
            self.restore.apply(record)?;
            let outcome = self.edit.record(record)?;
            if tally(&mut self.counters, outcome) {
                if rewrites {
                    record.set_codec_bin();
                }
                batch.swap(kept, at);
                kept += 1;
            }
        }
        batch.truncate(kept);
        Ok(())
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        self.edit.finish(out)?;
        if self.edit.rewrites() {
            out.iter_mut().for_each(RawRecord::set_codec_bin);
        }
        self.counters.records_out += out.len() as u64;
        Ok(())
    }

    fn counters(&self) -> Counters {
        self.counters
    }
}

#[cfg(test)]
mod tests;

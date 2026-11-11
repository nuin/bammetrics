//! `MarkDuplicates`, mirroring Picard 3.5.0 (`MarkDuplicates`,
//! `AbstractMarkDuplicatesCommandLineProgram`, `ReadEndsForMarkDuplicates`,
//! `LibraryIdGenerator`, `OpticalDuplicateFinder`, `DuplicationMetrics`)
//! and the htsjdk 5.0.0 pieces it runs on.
//!
//! Three steps, as in Picard:
//! 1. [`build`]: one read of the inputs makes a `ReadEnds` per mapped
//!    primary record (fragments) and per pair whose two mates are mapped
//!    (pairs). Each holds the library, the unclipped 5' position and
//!    strand of each end, the score, the tile, x and y of the read name and
//!    the index of each record in the file. Both lists go through a sort
//!    that spills to `TMP_DIR` ([`crate::spill`]), so memory holds at most
//!    `MAX_RECORDS_IN_RAM` entries per list.
//! 2. [`mark`]: the sorted pairs, then the sorted fragments, are cut into
//!    duplicate sets; in each set the first entry of highest score stays,
//!    every other record index goes to the duplicate list, optical
//!    duplicates to their own list, and set sizes to the histograms.
//! 3. [`write`]: a second read copies every record to `OUTPUT` with the
//!    duplicate flag set or cleared from those lists, the `DT`, `DS`, `DI`,
//!    `MI` and `PG` tags as asked, and counts the `DuplicationMetrics`.
//!
//! The rules that decide which record is flagged are in the submodules.

mod build;
mod inputs;
mod mark;
mod mates;
mod md_nm;
mod read_ends;
mod write;

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};

use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::optical::{
    DEFAULT_MAX_DUPLICATE_SET_SIZE, DEFAULT_OPTICAL_DUPLICATE_DISTANCE, OpticalDuplicateFinder,
    ReadNameParser,
};
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_text::Stringency;
use crate::spill::TempDir;
use crate::stage::{Barrier, Frame, Next, Sealed, Stage};

pub use build::unclipped;
pub use inputs::SortOrder;
pub use md_nm::{ContigCache, restore as restore_md_and_nm};

/// `DuplicateScoringStrategy.ScoringStrategy`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scoring {
    SumOfBaseQualities,
    TotalMappedReferenceLength,
    Random,
}

/// `MarkDuplicates.DuplicateTaggingPolicy`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TaggingPolicy {
    DontTag,
    OpticalOnly,
    All,
}

/// Entries held in memory per sorted list before a run is spilled.
pub const DEFAULT_RECORDS_IN_RAM: usize = 4_000_000;

/// Every option, parsed.
#[derive(Clone, Debug)]
pub struct Options {
    pub inputs: Vec<PathBuf>,
    pub output: PathBuf,
    pub metrics: PathBuf,
    pub reference: Option<PathBuf>,
    pub threads: usize,
    /// `VALIDATION_STRINGENCY`: how SAM text inputs are checked.
    pub stringency: Stringency,
    pub remove_duplicates: bool,
    pub remove_sequencing_duplicates: bool,
    /// `ASSUME_SORT_ORDER`, or `coordinate` for `ASSUME_SORTED=true`.
    pub assume_sort_order: Option<SortOrder>,
    pub scoring: Scoring,
    /// `None` for `PROGRAM_RECORD_ID=null`: no `@PG` record, no `PG` tag.
    pub program_record_id: Option<String>,
    pub program_group_version: Option<String>,
    pub program_group_command_line: Option<String>,
    /// `None` for `PROGRAM_GROUP_NAME=null`: the `@PG` record has no `PN`.
    pub program_group_name: Option<String>,
    pub add_pg_tag_to_reads: bool,
    pub comments: Vec<String>,
    pub barcode_tag: Option<String>,
    pub read_one_barcode_tag: Option<String>,
    pub read_two_barcode_tag: Option<String>,
    pub tag_duplicate_set_members: bool,
    pub tagging_policy: TaggingPolicy,
    pub clear_dt: bool,
    pub duplex_umi: bool,
    pub molecular_identifier_tag: Option<String>,
    pub finder: OpticalDuplicateFinder,
    pub compression_level: u8,
    pub create_index: bool,
    pub create_md5_file: bool,
    pub tmp_dir: Option<PathBuf>,
    pub records_in_ram: usize,
    /// The `@PG CL` default: the command line as bammetrics read it.
    pub command: String,
}

impl Options {
    pub fn use_barcodes(&self) -> bool {
        self.barcode_tag.is_some()
            || self.read_one_barcode_tag.is_some()
            || self.read_two_barcode_tag.is_some()
    }

    /// Optical duplicates get an index list of their own.
    pub fn index_optical(&self) -> bool {
        self.remove_sequencing_duplicates || self.tagging_policy != TaggingPolicy::DontTag
    }

    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        Self::parse(inv, false)
    }

    /// The options of the tool as a chain step ([`stage`]): it writes no
    /// BAM itself, so `OUTPUT` and the options of that file are refused,
    /// and `METRICS_FILE` is optional.
    pub fn for_stage(inv: &Invocation) -> Result<Self> {
        Self::parse(inv, true)
    }

    fn parse(inv: &Invocation, staged: bool) -> Result<Self> {
        let text = |name: &str| inv.get(name).filter(|v| *v != "null").map(str::to_owned);
        let given = |name: &str| inv.get(name).is_some();
        if given("ASSUME_SORTED") && given("ASSUME_SORT_ORDER") {
            return Err(Error::Usage(
                "ASSUME_SORTED and ASSUME_SORT_ORDER cannot be used together".to_owned(),
            ));
        }
        let assume_sort_order = match inv.get("ASSUME_SORT_ORDER") {
            None | Some("null") => inv
                .flag("ASSUME_SORTED", false)?
                .then_some(SortOrder::Coordinate),
            Some(value) => Some(SortOrder::parse(value).ok_or_else(|| {
                Error::Usage(format!(
                    "ASSUME_SORT_ORDER={value}: expected unsorted, queryname, coordinate, \
                     duplicate or unknown"
                ))
            })?),
        };
        let scoring = match inv.get("DUPLICATE_SCORING_STRATEGY") {
            None | Some("SUM_OF_BASE_QUALITIES") => Scoring::SumOfBaseQualities,
            Some("TOTAL_MAPPED_REFERENCE_LENGTH") => Scoring::TotalMappedReferenceLength,
            Some("RANDOM") => Scoring::Random,
            Some(value) => {
                return Err(Error::Usage(format!(
                    "DUPLICATE_SCORING_STRATEGY={value}: expected SUM_OF_BASE_QUALITIES, \
                     TOTAL_MAPPED_REFERENCE_LENGTH or RANDOM"
                )));
            }
        };
        let tagging_policy = match inv.get("TAGGING_POLICY") {
            None | Some("DontTag") => TaggingPolicy::DontTag,
            Some("OpticalOnly") => TaggingPolicy::OpticalOnly,
            Some("All") => TaggingPolicy::All,
            Some(value) => {
                return Err(Error::Usage(format!(
                    "TAGGING_POLICY={value}: expected DontTag, OpticalOnly or All"
                )));
            }
        };
        let compression_level: i64 = inv.parse("COMPRESSION_LEVEL", 5)?;
        let compression_level = u8::try_from(compression_level)
            .ok()
            .filter(|l| *l <= 9)
            .ok_or_else(|| {
                Error::Usage(format!(
                    "COMPRESSION_LEVEL={compression_level}: expected 0 to 9"
                ))
            })?;
        let records_in_ram: i64 = inv.parse("MAX_RECORDS_IN_RAM", DEFAULT_RECORDS_IN_RAM as i64)?;
        let records_in_ram = usize::try_from(records_in_ram)
            .ok()
            .filter(|n| *n >= 1)
            .ok_or_else(|| {
                Error::Usage(format!(
                    "MAX_RECORDS_IN_RAM={records_in_ram}: expected 1 or more"
                ))
            })?;
        refuse_flow_mode(inv)?;
        let (output, metrics) = if staged {
            stage_files(inv)?
        } else {
            let required = |path: Option<PathBuf>, name: &str| {
                path.ok_or_else(|| Error::Usage(format!("MarkDuplicates: {name} is required")))
            };
            (
                required(inv.get("OUTPUT").map(PathBuf::from), "OUTPUT")?,
                required(inv.output.clone(), "METRICS_FILE")?,
            )
        };
        Ok(Self {
            inputs: inv.all("INPUT").iter().map(PathBuf::from).collect(),
            output,
            metrics,
            reference: inv.reference().map(Path::to_path_buf),
            threads: inv.threads,
            remove_duplicates: inv.flag("REMOVE_DUPLICATES", false)?,
            remove_sequencing_duplicates: inv.flag("REMOVE_SEQUENCING_DUPLICATES", false)?,
            assume_sort_order,
            scoring,
            program_record_id: match inv.get("PROGRAM_RECORD_ID") {
                None => Some("MarkDuplicates".to_owned()),
                Some("null") => None,
                Some(id) => Some(id.to_owned()),
            },
            program_group_version: text("PROGRAM_GROUP_VERSION"),
            program_group_command_line: text("PROGRAM_GROUP_COMMAND_LINE"),
            program_group_name: match inv.get("PROGRAM_GROUP_NAME") {
                None => Some("MarkDuplicates".to_owned()),
                Some("null") => None,
                Some(name) => Some(name.to_owned()),
            },
            add_pg_tag_to_reads: inv.flag("ADD_PG_TAG_TO_READS", true)?,
            comments: crate::registry::picard_list(Vec::new(), inv.all("COMMENT"), |c| {
                Ok(c.to_owned())
            })?,
            barcode_tag: text("BARCODE_TAG"),
            read_one_barcode_tag: text("READ_ONE_BARCODE_TAG"),
            read_two_barcode_tag: text("READ_TWO_BARCODE_TAG"),
            tag_duplicate_set_members: inv.flag("TAG_DUPLICATE_SET_MEMBERS", false)?,
            tagging_policy,
            clear_dt: inv.flag("CLEAR_DT", true)?,
            duplex_umi: inv.flag("DUPLEX_UMI", false)?,
            molecular_identifier_tag: text("MOLECULAR_IDENTIFIER_TAG"),
            finder: OpticalDuplicateFinder {
                parser: ReadNameParser::new(inv.get("READ_NAME_REGEX"))?,
                pixel_distance: inv.parse(
                    "OPTICAL_DUPLICATE_PIXEL_DISTANCE",
                    DEFAULT_OPTICAL_DUPLICATE_DISTANCE,
                )?,
                max_set_size: inv.parse(
                    "MAX_OPTICAL_DUPLICATE_SET_SIZE",
                    DEFAULT_MAX_DUPLICATE_SET_SIZE,
                )?,
            },
            compression_level,
            create_index: inv.flag("CREATE_INDEX", false)?,
            create_md5_file: inv.flag("CREATE_MD5_FILE", false)?,
            stringency: inv.stringency()?,
            tmp_dir: inv.all("TMP_DIR").first().map(PathBuf::from),
            records_in_ram,
            command: inv.command.clone(),
        })
    }
}

/// The options of the BAM file, which a chain step does not write.
const WRITER_OPTIONS: [&str; 4] = [
    "OUTPUT",
    "COMPRESSION_LEVEL",
    "CREATE_INDEX",
    "CREATE_MD5_FILE",
];

/// `OUTPUT` and `METRICS_FILE` of a chain step: no BAM, and the metrics
/// path when one was given (the caller writes the file).
fn stage_files(inv: &Invocation) -> Result<(PathBuf, PathBuf)> {
    if let Some(name) = WRITER_OPTIONS.iter().find(|name| inv.get(name).is_some()) {
        return Err(Error::Usage(format!(
            "MarkDuplicates in a chain: {name} belongs to the write step; this step hands its \
             records to the next step"
        )));
    }
    if inv.all("INPUT").len() > 1 {
        return Err(Error::Usage(
            "MarkDuplicates in a chain reads the chain's one source: a second INPUT is not merged"
                .to_owned(),
        ));
    }
    Ok((PathBuf::new(), inv.output.clone().unwrap_or_default()))
}

/// The flow-based options with their Picard defaults. They are accepted at
/// the default only: flow-based duplicate marking is not implemented.
const FLOW_DEFAULTS: [(&str, &str); 9] = [
    ("FLOW_MODE", "false"),
    ("FLOW_DUP_STRATEGY", "FLOW_QUALITY_SUM_STRATEGY"),
    ("FLOW_EFFECTIVE_QUALITY_THRESHOLD", "15"),
    ("FLOW_Q_IS_KNOWN_END", "false"),
    ("FLOW_SKIP_FIRST_N_FLOWS", "0"),
    ("FLOW_UNPAIRED_END_UNCERTAINTY", "0"),
    ("FLOW_UNPAIRED_START_UNCERTAINTY", "0"),
    ("FLOW_USE_END_IN_UNPAIRED_READS", "false"),
    ("FLOW_USE_UNPAIRED_CLIPPED_END", "false"),
];

fn refuse_flow_mode(inv: &Invocation) -> Result<()> {
    for (name, default) in FLOW_DEFAULTS {
        let Some(value) = inv.get(name) else { continue };
        let at_default = match default {
            "false" => !inv.flag(name, false)?,
            "FLOW_QUALITY_SUM_STRATEGY" => value == default,
            _ => inv.parse::<i64>(name, 0)?.to_string() == default,
        };
        if !at_default {
            return Err(Error::Usage(format!(
                "{name}={value}: flow-based duplicate marking is not implemented; only the \
                 default {default} is accepted"
            )));
        }
    }
    Ok(())
}

pub const TOOL: Tool = Tool {
    name: "MarkDuplicates",
    summary: "one row per library (DuplicationMetrics) plus the set_size histograms; OUTPUT \
              gets every record with the duplicate flag set or cleared, METRICS_FILE the \
              metrics.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A BAM or CRAM to read; repeat to merge several.",
        )
        .repeated()
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The BAM to write: every record, duplicates flagged (0x400).",
        )
        .written()
        .required(),
        Opt::output(&["METRICS_FILE", "M"])
            .with_help("The duplication metrics file.")
            .required(),
        Opt::bool(
            &["REMOVE_DUPLICATES"],
            "Leave duplicates out of OUTPUT instead of flagging them (default false).",
        ),
        Opt::bool(
            &["REMOVE_SEQUENCING_DUPLICATES"],
            "Leave optical duplicates out of OUTPUT (default false).",
        ),
        Opt::bool(
            &["ASSUME_SORTED", "AS"],
            "true reads the input as coordinate sorted whatever its header says (default false).",
        ),
        Opt::text(
            &["ASSUME_SORT_ORDER", "ASO"],
            "Read the input in this order (coordinate or queryname) whatever its header says.",
        ),
        Opt::text(
            &["DUPLICATE_SCORING_STRATEGY", "DS"],
            "SUM_OF_BASE_QUALITIES (default), TOTAL_MAPPED_REFERENCE_LENGTH or RANDOM.",
        ),
        Opt::text(
            &["PROGRAM_RECORD_ID", "PG"],
            "ID of the @PG record added (default MarkDuplicates); null adds none.",
        ),
        Opt::text(
            &["PROGRAM_GROUP_VERSION", "PG_VERSION"],
            "VN of the @PG record (default: the bammetrics version).",
        ),
        Opt::text(
            &["PROGRAM_GROUP_COMMAND_LINE", "PG_COMMAND"],
            "CL of the @PG record (default: the command line).",
        ),
        Opt::text(
            &["PROGRAM_GROUP_NAME", "PG_NAME"],
            "PN of the @PG record (default MarkDuplicates).",
        ),
        Opt::bool(
            &["ADD_PG_TAG_TO_READS"],
            "Set each record's PG tag to the new @PG ID (default true).",
        ),
        Opt::text(&["COMMENT", "CO"], "A @CO line for the OUTPUT header.").repeated(),
        Opt::text(
            &["BARCODE_TAG"],
            "Tag of the molecular barcode (UMI), or null.",
        ),
        Opt::text(
            &["READ_ONE_BARCODE_TAG"],
            "Tag of the first read's barcode, or null.",
        ),
        Opt::text(
            &["READ_TWO_BARCODE_TAG"],
            "Tag of the second read's barcode, or null.",
        ),
        Opt::bool(
            &["TAG_DUPLICATE_SET_MEMBERS"],
            "Add DS (set size) and DI (representative index) to duplicate set members (default \
             false).",
        ),
        Opt::text(
            &["TAGGING_POLICY"],
            "DT tag on duplicates: DontTag (default), OpticalOnly (SQ) or All (SQ and LB).",
        ),
        Opt::bool(
            &["CLEAR_DT"],
            "Remove DT tags found in the input (default true).",
        ),
        Opt::bool(
            &["DUPLEX_UMI"],
            "Read BARCODE_TAG as a duplex UMI X-Y, swapped for the bottom strand (default false).",
        ),
        Opt::text(
            &["MOLECULAR_IDENTIFIER_TAG"],
            "Tag for contig:position/UMI molecule identifiers; needs BARCODE_TAG.",
        ),
        Opt::text(
            &["READ_NAME_REGEX"],
            "Regex with tile, x and y groups; null for no optical duplicates (default: \
             5- and 7-field Illumina names). Rust regex syntax: no look-around, no \\Q...\\E.",
        ),
        Opt::int(
            &["OPTICAL_DUPLICATE_PIXEL_DISTANCE"],
            "Largest pixel distance of optical duplicates (default 100).",
        ),
        Opt::int(
            &["MAX_OPTICAL_DUPLICATE_SET_SIZE"],
            "Larger duplicate sets are not checked for optical duplicates (default 300000).",
        ),
        Opt::int(
            &["COMPRESSION_LEVEL"],
            "BGZF compression level of OUTPUT, 0 to 9 (default 5).",
        ),
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside a coordinate-sorted OUTPUT (default false).",
        ),
        Opt::bool(&["CREATE_MD5_FILE"], "Write OUTPUT.md5 (default false)."),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled sort runs (default: the system temporary directory).",
        )
        .repeated(),
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Entries each sorted list holds in memory before it spills to TMP_DIR, 1 or more \
             (default 4000000).",
        ),
        Opt::ignored_int(
            &["MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP", "MAX_SEQS"],
            "obsolete in Picard",
        ),
        Opt::ignored_int(
            &["MAX_FILE_HANDLES_FOR_READ_ENDS_MAP", "MAX_FILE_HANDLES"],
            "waiting mates spill to one file per reference, opened one at a time",
        ),
        Opt::ignored_float(
            &["SORTING_COLLECTION_SIZE_RATIO"],
            "MAX_RECORDS_IN_RAM sets the sort buffers",
        ),
        Opt::bool(
            &["FLOW_MODE"],
            "Flow-based reads are not supported: only false, the default, is accepted.",
        ),
        Opt::text(
            &["FLOW_DUP_STRATEGY"],
            "Only FLOW_QUALITY_SUM_STRATEGY, the default, is accepted.",
        ),
        Opt::int(
            &["FLOW_EFFECTIVE_QUALITY_THRESHOLD"],
            "Only 15, the default, is accepted.",
        ),
        Opt::bool(
            &["FLOW_Q_IS_KNOWN_END"],
            "Only false, the default, is accepted.",
        ),
        Opt::int(
            &["FLOW_SKIP_FIRST_N_FLOWS"],
            "Only 0, the default, is accepted.",
        ),
        Opt::int(
            &["FLOW_UNPAIRED_END_UNCERTAINTY"],
            "Only 0, the default, is accepted.",
        ),
        Opt::int(
            &["FLOW_UNPAIRED_START_UNCERTAINTY"],
            "Only 0, the default, is accepted.",
        ),
        Opt::bool(
            &["FLOW_USE_END_IN_UNPAIRED_READS"],
            "Only false, the default, is accepted.",
        ),
        Opt::bool(
            &["FLOW_USE_UNPAIRED_CLIPPED_END"],
            "Only false, the default, is accepted.",
        ),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    run(&opts).map(ToolOutput::from)
}

/// Mark duplicates as `opts` say: write `OUTPUT` (and its index and MD5
/// when asked) and return the metrics file for `METRICS_FILE`. `OUTPUT`,
/// its index and its MD5 are written under temporary names and renamed
/// into place at the end, so a failed run leaves no partial file.
pub fn run(opts: &Options) -> Result<crate::metric::MetricsOutput> {
    for path in [&opts.output, &opts.metrics] {
        check_writable(path)?;
    }
    if has_extension(&opts.output, "cram") || has_extension(&opts.output, "sam") {
        return Err(Error::Usage(format!(
            "OUTPUT={}: bammetrics writes BAM only (SAM and CRAM output are not supported)",
            opts.output.display()
        )));
    }
    for input in opts.inputs.iter().chain(&opts.reference) {
        if !input.exists() {
            return Err(Error::Input {
                path: input.clone(),
                source: io::ErrorKind::NotFound.into(),
            });
        }
    }
    refuse_overlaps(opts)?;
    let temp = TempDir::new(opts.tmp_dir.as_deref());
    let inputs = spool_streams(&opts.inputs, &temp)?;
    // An input copied from a pipe fails under the name the user gave.
    let named = |mut error: Error| {
        for (copy, given) in inputs.iter().zip(&opts.inputs) {
            error = error.naming_input(copy, given);
        }
        error
    };
    let built = build::run(opts, &inputs, &temp).map_err(named)?;
    let marks = mark::run(opts, built, &temp)?;
    write::run(opts, &inputs, marks).map_err(named)
}

/// `MarkDuplicates` as a step of a chain ([`crate::stage`]): a barrier.
/// The records of the step before go through the first read
/// ([`build::Builder::push`]); sealing it finds the duplicate sets and
/// leaves the second read ([`write::Copier`]) as the transform the source
/// is replayed through, so the next step sees every record with its
/// duplicate flag, tags and `PG` decided, under the header with the new
/// `@PG` line. Spilled runs go to `temp`.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    let opts = Options::for_stage(inv)?;
    let mut header = frame.header.clone();
    if let Some(order) = opts.assume_sort_order {
        inputs::set_header_tag(&mut header, inputs::HeaderTag::SortOrder, order.name());
    }
    let builder = build::Builder::new(&opts, &header, temp)?;
    Ok(Stage::Barrier(Box::new(Marking {
        opts,
        header,
        cram: frame.cram,
        builder,
        temp,
    })))
}

/// The barrier of [`stage`]: the first read, fed in batches.
struct Marking<'t> {
    opts: Options,
    /// The header of the records pushed, `ASSUME_SORT_ORDER` applied.
    header: noodles::sam::Header,
    cram: bool,
    builder: build::Builder<'t>,
    temp: &'t TempDir,
}

impl<'t> Barrier<'t> for Marking<'t> {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()> {
        batch
            .iter()
            .try_for_each(|record| self.builder.push(record))
    }

    fn seal(self: Box<Self>) -> Result<Sealed<'t>> {
        let Marking {
            opts,
            header,
            cram,
            builder,
            temp,
        } = *self;
        let marks = mark::run(&opts, builder.finish(), temp)?;
        let (output, chained) = write::output_header(&header, &opts, &marks.pg_ids);
        let copier = write::Copier::new(&opts, &header, chained, cram, marks)?;
        Ok(Sealed {
            header: output,
            next: Next::Replay(Box::new(copier)),
        })
    }
}

/// The files a run writes, with the option that names each.
fn written_files(opts: &Options) -> Vec<(&'static str, PathBuf)> {
    let mut files = vec![
        ("OUTPUT", opts.output.clone()),
        ("METRICS_FILE", opts.metrics.clone()),
    ];
    if opts.create_index {
        files.push(("CREATE_INDEX", opts.output.with_extension("bai")));
    }
    if opts.create_md5_file {
        files.push(("CREATE_MD5_FILE", crate::bam_writer::md5_path(&opts.output)));
    }
    files
}

/// Refuse a written file that is an input (writing it would destroy the
/// input before the second read) or that another written file names too.
/// Paths compare after symbolic links and `..` are resolved; on Unix a
/// hard link to an input is found by device and inode. Device files are
/// not destroyed by writing and are not compared.
fn refuse_overlaps(opts: &Options) -> Result<()> {
    let device = |path: &Path| same_file_key(path).starts_with("/dev");
    let read: Vec<&PathBuf> = opts
        .inputs
        .iter()
        .chain(&opts.reference)
        .filter(|p| !device(p))
        .collect();
    let written = written_files(opts);
    for (i, (name, path)) in written.iter().enumerate() {
        if device(path) {
            continue;
        }
        if let Some(input) = read.iter().find(|input| same_file(input, path)) {
            return Err(Error::Usage(format!(
                "{name}={}: it is the input {}; writing it would destroy it",
                path.display(),
                input.display()
            )));
        }
        if let Some((other, _)) = written[..i].iter().find(|(_, o)| same_file(o, path)) {
            return Err(Error::Usage(format!(
                "{name}={}: {other} names the same file; one file cannot hold both",
                path.display()
            )));
        }
    }
    Ok(())
}

/// `path` absolute with its links resolved, or with its directory resolved
/// when the file does not exist yet.
pub(crate) fn same_file_key(path: &Path) -> PathBuf {
    if let Ok(resolved) = path.canonicalize() {
        return resolved;
    }
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    match (absolute.parent(), absolute.file_name()) {
        (Some(dir), Some(name)) => dir
            .canonicalize()
            .map_or_else(|_| absolute.clone(), |d| d.join(name)),
        _ => absolute,
    }
}

pub(crate) fn same_file(a: &Path, b: &Path) -> bool {
    #[cfg(unix)]
    {
        use std::os::unix::fs::MetadataExt;
        if let (Ok(x), Ok(y)) = (fs::metadata(a), fs::metadata(b))
            && (x.dev(), x.ino()) == (y.dev(), y.ino())
        {
            return true;
        }
    }
    same_file_key(a) == same_file_key(b)
}

fn has_extension(path: &Path, extension: &str) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case(extension))
}

/// `IOUtil.assertFileIsWritable`: an existing file must not be a
/// directory, and a new one needs an existing directory.
fn check_writable(path: &Path) -> Result<()> {
    let error = |message: &str| Error::Output {
        path: path.to_path_buf(),
        source: io::Error::other(message.to_owned()),
    };
    if path.is_dir() {
        return Err(error("is a directory"));
    }
    if path.exists() || path.starts_with("/dev") {
        return Ok(());
    }
    let parent = path.parent().filter(|p| !p.as_os_str().is_empty());
    match parent {
        Some(dir) if !dir.is_dir() => Err(error("its directory does not exist")),
        _ => Ok(()),
    }
}

/// The inputs as files that can be read twice: a pipe or standard input
/// is copied into the temporary directory first.
pub(crate) fn spool_streams(inputs: &[PathBuf], temp: &TempDir) -> Result<Vec<PathBuf>> {
    inputs
        .iter()
        .map(|input| {
            let regular = fs::metadata(input).is_ok_and(|m| m.is_file());
            if regular && !input.starts_with("/dev") {
                return Ok(input.clone());
            }
            let read_error = |source| Error::Input {
                path: input.clone(),
                source,
            };
            let mut source = File::open(input).map_err(read_error)?;
            let (copy, mut sink) = temp.create_file("input")?;
            copy_all(&mut source, &mut sink).map_err(read_error)?;
            Ok(copy)
        })
        .collect()
}

fn copy_all(source: &mut impl Read, sink: &mut File) -> io::Result<()> {
    io::copy(source, sink).map(|_| ())
}

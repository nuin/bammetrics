//! `SortSam`, mirroring Picard 3.5.0 (`picard.sam.SortSam`): the records of
//! a BAM, CRAM or SAM in `queryname`, `coordinate` or `duplicate` order, as
//! a BAM, or as SAM text for an `OUTPUT` ending in `.sam`.
//!
//! The output header is the input's with `SO` set to the order and `@HD
//! VN:1.6`, the version htsjdk's writers give it. `GO`, `SS` and every
//! other line stay as they are, and no `@PG` line is added (Picard 3.5.0
//! adds none). The records go through [`crate::sam_sort`]: htsjdk's
//! comparator for the order, Java's `TimSort` per run, runs of at most
//! `MAX_RECORDS_IN_RAM` records spilled to `TMP_DIR` and merged with ties to
//! the earlier run, so records with equal keys keep their input order as
//! in Picard. A CRAM input gets `MD` and `NM` rebuilt as htsjdk decodes
//! them. `CREATE_INDEX` writes a `.bai` beside a coordinate-sorted BAM;
//! beside SAM text it writes nothing, as in Picard. `CREATE_MD5_FILE`
//! writes `OUTPUT.md5` for either.
//!
//! In a chain ([`stage`]) the tool is a barrier. It takes every record of
//! its position, then hands on the sorted records ([`Next::Records`]) under
//! the sorted header. Its buffer holds at most `--sort-memory` bytes
//! ([`Invocation::sort_memory`]) and, when the step gives it,
//! `MAX_RECORDS_IN_RAM` records. `SORT_ORDER` is checked before the first
//! read ([`Options::for_stage`]).

use std::io;
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, BamWriter, Companions};
use crate::error::{Error, Result};
use crate::mark_duplicates::{ContigCache, restore_md_and_nm, same_file};
use crate::metric::ToolOutput;
use crate::pass;
use crate::reference::Reference;
use crate::registry::{Invocation, Opt, Tool};
use crate::sam_order::Order;
use crate::sam_sort::{self, Limits, RecordSorter};
use crate::sam_text::{self, SamWriter};
use crate::spill::{SortStats, TempDir};
use crate::stage::{
    Barrier, Counters, Frame, Next, Products, Sealed, Stage, WriteFormat, later_kind,
};

/// Picard's `MAX_RECORDS_IN_RAM` default.
pub const DEFAULT_RECORDS_IN_RAM: usize = 500_000;

pub const TOOL: Tool = Tool {
    name: "SortSam",
    summary: "no rows; OUTPUT gets every record in SORT_ORDER (queryname, coordinate or \
              duplicate) under the input header with SO set, as a BAM, or as SAM text for a \
              name ending in .sam.",
    options: &[
        Opt::path(
            &["OUTPUT", "O"],
            "The sorted BAM, or SAM text for a name ending in .sam.",
        )
        .written()
        .required(),
        Opt::text(
            &["SORT_ORDER", "SO"],
            "queryname, coordinate or duplicate (duplicate needs the MC tag on paired records \
             with a mapped mate).",
        )
        .required(),
        Opt::int(
            &["COMPRESSION_LEVEL"],
            "BGZF compression level of OUTPUT, 0 to 9 (default 5).",
        ),
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside a coordinate-sorted BAM OUTPUT (default false).",
        ),
        Opt::bool(&["CREATE_MD5_FILE"], "Write OUTPUT.md5 (default false)."),
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more (default \
             500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs (default: the system temporary directory).",
        )
        .repeated(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The options of one run.
#[derive(Clone, Debug)]
pub struct Options {
    pub input: PathBuf,
    pub output: PathBuf,
    pub order: Order,
    pub reference: Option<PathBuf>,
    pub compression_level: u8,
    pub create_index: bool,
    pub create_md5: bool,
    pub limits: Limits,
    pub tmp_dir: Option<PathBuf>,
    pub stringency: sam_text::Stringency,
}

/// The options of the BAM file, which a chain step does not write.
const WRITER_OPTIONS: [&str; 4] = [
    "OUTPUT",
    "COMPRESSION_LEVEL",
    "CREATE_INDEX",
    "CREATE_MD5_FILE",
];

/// `SORT_ORDER`, with Picard's message for a value `SortSam` does not take.
fn sort_order(inv: &Invocation) -> Result<Order> {
    let value = inv
        .get("SORT_ORDER")
        .ok_or_else(|| Error::Usage("SortSam: SORT_ORDER is required".to_owned()))?;
    match Order::parse(value) {
        Some(Order::Unsorted) | None => Err(Error::Usage(format!(
            "'{value}' is not a valid value for SortOrder. Possible values: queryname, \
             coordinate, duplicate"
        ))),
        Some(order) => Ok(order),
    }
}

/// `MAX_RECORDS_IN_RAM`, `None` when not given.
fn records_in_ram(inv: &Invocation) -> Result<Option<usize>> {
    let Some(n) = inv.parse_opt::<i64>("MAX_RECORDS_IN_RAM")? else {
        return Ok(None);
    };
    usize::try_from(n)
        .ok()
        .filter(|n| *n >= 1)
        .map(Some)
        .ok_or_else(|| Error::Usage(format!("MAX_RECORDS_IN_RAM={n}: expected 1 or more")))
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let order = sort_order(inv)?;
        let output = inv
            .get("OUTPUT")
            .map(PathBuf::from)
            .ok_or_else(|| Error::Usage("SortSam: OUTPUT is required".to_owned()))?;
        let compression_level: i64 = inv.parse("COMPRESSION_LEVEL", 5)?;
        let compression_level = u8::try_from(compression_level)
            .ok()
            .filter(|l| *l <= 9)
            .ok_or_else(|| {
                Error::Usage(format!(
                    "COMPRESSION_LEVEL={compression_level}: expected 0 to 9"
                ))
            })?;
        let records = records_in_ram(inv)?.unwrap_or(DEFAULT_RECORDS_IN_RAM);
        Ok(Self {
            input: inv.input.clone(),
            output,
            order,
            reference: inv.reference().map(Path::to_path_buf),
            compression_level,
            create_index: inv.flag("CREATE_INDEX", false)?,
            create_md5: inv.flag("CREATE_MD5_FILE", false)?,
            limits: Limits {
                records,
                bytes: inv.sort_memory.unwrap_or(usize::MAX),
                threads: inv.threads,
            },
            tmp_dir: inv.all("TMP_DIR").first().map(PathBuf::from),
            stringency: inv.stringency()?,
        })
    }

    /// The options of the tool as a chain step: no BAM of its own, so
    /// `OUTPUT` and its options are refused; the buffer holds at most
    /// `--sort-memory` bytes and, when given, `MAX_RECORDS_IN_RAM` records.
    pub fn for_stage(inv: &Invocation) -> Result<(Order, Limits)> {
        if let Some(name) = WRITER_OPTIONS.iter().find(|name| inv.get(name).is_some()) {
            return Err(Error::Usage(format!(
                "SortSam in a chain: {name} belongs to the write step; this step hands its \
                 records to the next step"
            )));
        }
        let order = sort_order(inv)?;
        Ok((
            order,
            Limits {
                records: records_in_ram(inv)?.unwrap_or(usize::MAX),
                bytes: inv.sort_memory.unwrap_or(usize::MAX),
                threads: inv.threads,
            },
        ))
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    run(&opts)?;
    Ok(ToolOutput::NONE)
}

/// Where the sorted records go: a BAM or SAM text.
enum Out {
    Bam(BamWriter),
    Sam(SamWriter),
}

impl Out {
    fn write(&mut self, record: &RawRecord) -> Result<()> {
        match self {
            Self::Bam(writer) => writer.write(record),
            Self::Sam(writer) => writer.write(record),
        }
    }
}

/// Refuse an `OUTPUT` that would destroy an input, or that `SortSam`
/// cannot write, before anything is read.
fn check_output(opts: &Options) -> Result<()> {
    if let Some((kind, wave)) = later_kind(&opts.output) {
        return Err(Error::Usage(format!(
            "OUTPUT={}: bammetrics writes BAM and SAM only: {kind} output is wave {wave} of \
             docs/07_full_parity_program.md",
            opts.output.display()
        )));
    }
    let written = [
        Some(opts.output.clone()),
        opts.create_index.then(|| opts.output.with_extension("bai")),
        opts.create_md5.then(|| bam_writer::md5_path(&opts.output)),
    ];
    let read = std::iter::once(&opts.input).chain(&opts.reference);
    for input in read {
        if !input.exists() {
            return Err(Error::Input {
                path: input.clone(),
                source: io::ErrorKind::NotFound.into(),
            });
        }
        if let Some(path) = written.iter().flatten().find(|path| same_file(input, path)) {
            return Err(Error::Usage(format!(
                "OUTPUT={}: it is the input {}; writing it would destroy it",
                path.display(),
                input.display()
            )));
        }
    }
    Ok(())
}

/// Sort `opts.input` into `opts.output` (and its index and MD5 when asked),
/// each under a temporary name until every file is complete. Returns what
/// the sort did.
pub fn run(opts: &Options) -> Result<SortStats> {
    check_output(opts)?;
    let temp = TempDir::new(opts.tmp_dir.as_deref());
    pass::with_record_stream_at(
        &opts.input,
        opts.reference.as_deref(),
        opts.limits.threads,
        opts.stringency,
        |header, reader, cram| {
            let mut contigs = match (cram, &opts.reference) {
                (true, Some(fasta)) => Some(ContigCache::new(Reference::open(fasta, header)?)),
                _ => None,
            };
            let mut sorter = RecordSorter::new(opts.order, header, &temp, opts.limits);
            let mut record = RawRecord::default();
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&opts.input, cram, e))?
            {
                if let Some(contigs) = contigs.as_mut() {
                    restore(&mut record, contigs)?;
                }
                sorter.push(std::mem::take(&mut record))?;
            }
            let (mut sorted, stats) = sorter.finish()?;
            let out_header = sam_sort::sorted_header(header, opts.order);
            let mut out = match WriteFormat::of(&opts.output) {
                WriteFormat::Sam => {
                    let text = sam_text::header_text(&out_header, false);
                    Out::Sam(SamWriter::create(&opts.output, &out_header, &text)?)
                }
                WriteFormat::Bam => Out::Bam(BamWriter::create(
                    &opts.output,
                    opts.compression_level,
                    opts.limits.threads,
                    &out_header,
                )?),
            };
            while let Some(record) = sorted.next_record()? {
                out.write(&record)?;
            }
            finish(out, &out_header, opts)?;
            Ok(stats)
        },
    )
}

/// `MD` and `NM` of a CRAM record, rebuilt as htsjdk decodes them; the tags
/// come back in htsjdk's order.
fn restore(record: &mut RawRecord, contigs: &mut ContigCache) -> Result<()> {
    let mut tags = record.tags().map_err(Error::Record)?;
    restore_md_and_nm(record, &mut tags, contigs)?;
    record.set_tags(&tags).map_err(Error::Record)
}

/// Close the output, write its companions, rename everything into place.
fn finish(out: Out, header: &sam::Header, opts: &Options) -> Result<()> {
    match out {
        Out::Bam(writer) => {
            let companions = Companions {
                index: opts.create_index,
                md5: opts.create_md5,
            };
            bam_writer::finish(writer, header, companions).map(|_| ())
        }
        Out::Sam(writer) => {
            let output = writer.close()?;
            let target = output.path().to_path_buf();
            let md5 = opts
                .create_md5
                .then(|| bam_writer::write_md5(output.writing(), &bam_writer::md5_path(&target)))
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

/// `SortSam` as a step of a chain: a barrier that holds the records of its
/// position and hands them on sorted, under the sorted header. Runs spill
/// to `temp`.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    let (order, limits) = Options::for_stage(inv)?;
    let contigs = match (frame.cram, inv.reference()) {
        (true, Some(fasta)) => Some(ContigCache::new(Reference::open(fasta, &frame.header)?)),
        _ => None,
    };
    Ok(Stage::Barrier(Box::new(Sorting {
        sorter: RecordSorter::new(order, &frame.header, temp, limits),
        header: sam_sort::sorted_header(&frame.header, order),
        contigs,
        records: 0,
    })))
}

/// The barrier of [`stage`].
struct Sorting<'t> {
    sorter: RecordSorter<'t>,
    header: sam::Header,
    /// The reference that rebuilds `MD` and `NM` of CRAM records.
    contigs: Option<ContigCache>,
    records: u64,
}

impl<'t> Barrier<'t> for Sorting<'t> {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()> {
        for record in batch {
            let mut record = record.clone();
            if let Some(contigs) = self.contigs.as_mut() {
                restore(&mut record, contigs)?;
            }
            self.sorter.push(record)?;
            self.records += 1;
        }
        Ok(())
    }

    fn seal(self: Box<Self>) -> Result<Sealed<'t>> {
        let Sorting {
            sorter,
            header,
            records,
            ..
        } = *self;
        let (sorted, stats) = sorter.finish()?;
        let mut notes = Vec::new();
        if stats.spilled_bytes > 0 {
            notes.push(format!(
                "spilled {} bytes in sorted runs, {} merge level{}",
                stats.spilled_bytes,
                stats.merge_levels,
                if stats.merge_levels == 1 { "" } else { "s" }
            ));
        }
        Ok(Sealed {
            header,
            next: Next::Records(Box::new(sorted)),
            products: Products {
                counters: Some(Counters {
                    records_in: records,
                    records_out: records,
                    changed: 0,
                }),
                notes,
                ..Products::default()
            },
        })
    }
}

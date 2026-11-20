//! `GatherBamFiles`, mirroring Picard 3.5.0 `picard.sam.GatherBamFiles`:
//! the records of every `INPUT`, in the order given, under the header of
//! the first.
//!
//! `INPUT` is unrolled as htsjdk's `IOUtil.unrollFiles(INPUT, ".bam",
//! ".sam", ".cram")` unrolls it: a path whose file name ends in one of the
//! three is an input; any other path is a list of paths, one per line,
//! each line trimmed of Java whitespace, empty lines skipped, a path on a
//! line resolved from the current directory, and a list may name lists.
//!
//! When every input name ends in `.bam`, Picard block copies
//! (`BamFileIoUtils.gatherWithBlockCopying`, [`crate::bgzf_copy`]): the
//! first file from its first byte, each later file from the block of its
//! first record, each without its end-of-file block, then one end-of-file
//! block. Nothing is decoded but the headers. With `COMPRESSION_LEVEL=0`,
//! or with one input, the output is Picard's file byte for byte; at levels
//! 1 to 9 the one block per later input that htsjdk deflates again holds
//! zlib-rs bytes, which inflate to Picard's. `CREATE_INDEX` writes
//! `<basename>.bai` beside the output (`out.bai` for `out.bam` and for
//! `out.dat`) when the first header says `SO:coordinate` or `SO:unsorted`
//! or has no `SO`, as htsjdk's `BAMIndexer` takes them
//! ([`bam_writer::finish_file`]); `CREATE_MD5_FILE` writes `<OUTPUT>.md5`.
//!
//! Any other name sends Picard down `gatherNormally`: the header of the
//! first input, then every record of every input read one by one (SAM text
//! at `VALIDATION_STRINGENCY`, a CRAM with `REFERENCE_SEQUENCE`) into a
//! presorted `SAMFileWriter`, which takes each record's reference by name
//! (`setHeaderStrict`) and checks the header's order: a BAM, or SAM text
//! for an `OUTPUT` ending in `.sam`, as [`crate::record_edit`] writes them.
//!
//! What bammetrics does differently (docs/07 decision rows):
//!
//! - A block copy of an input with records whose binary sequence
//!   dictionary (names and lengths) is not the first input's is refused.
//!   Picard writes the file, and each of those records keeps a reference
//!   index of another dictionary.
//! - A block copy into an `OUTPUT` named `.sam` or `.cram` is refused:
//!   Picard writes BAM bytes under that name.
//! - An `OUTPUT` that is an input, a list file, the reference or its
//!   `.fai` is refused before the first read; Picard empties it. A failed
//!   run leaves no output; Picard leaves a part of the file.
//! - The copied blocks are walked header by header: a byte that starts no
//!   BGZF block fails the run, where Picard copies it.
//!
//! No `@PG` line is added; Picard adds none.

use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, Companions};
use crate::bgzf_copy::{self, Source};
use crate::cleanup::AtomicFile;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::mark_duplicates::same_file;
use crate::metric::ToolOutput;
use crate::pass;
use crate::record_edit::{self, Mode, Out, Restore};
use crate::reference::files_read;
use crate::registry::{self, Invocation, Opt, Tool};
use crate::stage::{WriteFormat, later_kind};

pub const TOOL: Tool = Tool {
    name: "GatherBamFiles",
    summary: "no rows; OUTPUT gets the records of every INPUT in the order given under the \
              header of the first: BAM inputs block copied as Picard copies them, byte for \
              byte at COMPRESSION_LEVEL=0; with a SAM or CRAM input every record read and \
              written again.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A BAM, SAM or CRAM to gather, or a text file of such paths, one per line; repeat \
             in the order to gather.",
        )
        .repeated()
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The BAM to write, or SAM text for a name ending in .sam when an input is SAM or \
             CRAM.",
        )
        .written()
        .required(),
        record_edit::COMPRESSION_OPTION,
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside a BAM OUTPUT, named OUTPUT less its extension plus .bai for \
             a block copy, which is indexed under SO:coordinate, SO:unsorted or no @HD; \
             written records under SO:coordinate (default false).",
        ),
        record_edit::MD5_OPTION,
        record_edit::RECORDS_IN_RAM_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The input names `IOUtil.unrollFiles` keeps as they are.
const EXTENSIONS: [&str; 3] = [".bam", ".sam", ".cram"];
/// Inputs after the lists are read, held and waiting together: a list
/// that names itself would never end.
pub const MAX_INPUTS: usize = 100_000;
/// List files read in one run.
pub const MAX_LISTS: usize = 10_000;
/// The bytes of one list file.
pub const MAX_LIST_BYTES: u64 = 64 << 20;
/// The bytes of every list file of one run together.
pub const MAX_TOTAL_LIST_BYTES: u64 = 256 << 20;

/// `INPUT` unrolled as `IOUtil.unrollPaths` unrolls it: a stack of paths,
/// each list replaced by its lines. Returns the inputs in the order given,
/// then the list files read.
pub fn unroll(given: &[String]) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    unroll_within(given, MAX_TOTAL_LIST_BYTES)
}

/// [`unroll`], reading at most `total` bytes of list files in all.
fn unroll_within(given: &[String], total: u64) -> Result<(Vec<PathBuf>, Vec<PathBuf>)> {
    let mut left = total;
    let mut stack: Vec<PathBuf> = given.iter().map(PathBuf::from).collect();
    let mut inputs = Vec::new();
    let mut lists = Vec::new();
    while let Some(path) = stack.pop() {
        let name = path
            .file_name()
            .map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_default();
        if EXTENSIONS.iter().any(|extension| name.ends_with(extension)) {
            inputs.push(path);
            continue;
        }
        if lists.len() == MAX_LISTS {
            return Err(Error::Usage(format!(
                "INPUT: more than {MAX_LISTS} list files to read; does a list name itself?"
            )));
        }
        let text = list_text(&path, left, total)?;
        left -= text.len() as u64;
        for line in list_lines(&text) {
            if stack.len() + inputs.len() == MAX_INPUTS {
                return Err(Error::Usage(format!(
                    "INPUT: more than {MAX_INPUTS} inputs; does a list name itself?"
                )));
            }
            stack.push(PathBuf::from(line));
        }
        lists.push(path);
    }
    inputs.reverse();
    Ok((inputs, lists))
}

/// The paths of a list's text: `Files.lines` (lines ended by `\n`, `\r` or
/// `\r\n`), each `String.trim`med, empty lines dropped. Lazy, so the caller
/// stops at [`MAX_INPUTS`] and no line beyond it is held.
fn list_lines(text: &str) -> impl Iterator<Item = &str> {
    text.split('\n')
        .flat_map(|line| line.split('\r'))
        .map(|line| line.trim_matches(|c: char| c <= ' '))
        .filter(|line| !line.is_empty())
}

/// The text of the list file `path`, at most [`MAX_LIST_BYTES`] of UTF-8
/// (`Files.lines` fails on other bytes) and at most `left`, what remains
/// of the run's `total` for list files.
fn list_text(path: &Path, left: u64, total: u64) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    let mut bytes = Vec::new();
    BufReader::new(file)
        .take(MAX_LIST_BYTES.min(left) + 1)
        .read_to_end(&mut bytes)
        .map_err(input)?;
    let read = bytes.len() as u64;
    if read > MAX_LIST_BYTES {
        return Err(input(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("a list of inputs over {MAX_LIST_BYTES} bytes"),
        )));
    }
    if read > left {
        return Err(input(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("the lists of inputs of this run hold over {total} bytes together"),
        )));
    }
    String::from_utf8(bytes).map_err(|_| {
        input(io::Error::new(
            io::ErrorKind::InvalidData,
            "a list of inputs that is not UTF-8 text (Files.lines fails on it)",
        ))
    })
}

/// `IOUtil.assertFileIsReadable`: each input exists, is no directory and
/// opens.
fn check_readable(inputs: &[PathBuf]) -> Result<()> {
    for path in inputs {
        let input = |source| Error::Input {
            path: path.clone(),
            source,
        };
        if std::fs::metadata(path).map_err(input)?.is_dir() {
            return Err(input(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Cannot read file because it is a directory",
            )));
        }
        File::open(path).map_err(input)?;
    }
    Ok(())
}

/// htsjdk's `IOUtil.basename`: the file name without its last extension.
fn basename(output: &Path) -> String {
    let name = output
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_default();
    match name.rfind('.') {
        Some(dot) if dot > 0 => name[..dot].to_owned(),
        _ => name,
    }
}

/// Where a block copy's `CREATE_INDEX` writes: `<basename>.bai` in the
/// output's directory.
pub fn block_copy_index(output: &Path) -> PathBuf {
    output.with_file_name(format!("{}.bai", basename(output)))
}

/// The options of one gather.
struct Gather {
    inputs: Vec<PathBuf>,
    /// The files read by name besides the inputs: list files, the
    /// reference and its `.fai`.
    others: Vec<PathBuf>,
    opts: record_edit::Options,
}

impl Gather {
    fn block_copy(&self) -> bool {
        self.inputs.iter().all(|path| {
            path.file_name()
                .is_some_and(|n| n.to_string_lossy().ends_with(".bam"))
        })
    }

    fn index_path(&self) -> PathBuf {
        if self.block_copy() {
            block_copy_index(&self.opts.output)
        } else {
            self.opts.index_path()
        }
    }

    /// Every file the run may write.
    fn written(&self) -> Vec<PathBuf> {
        let mut files = vec![self.opts.output.clone()];
        if self.opts.create_index {
            files.push(self.index_path());
        }
        if self.opts.create_md5 {
            files.push(bam_writer::md5_path(&self.opts.output));
        }
        files
    }

    /// Refuse an `OUTPUT` the run cannot write, or one that would destroy
    /// a file it reads, before anything is read.
    fn check_output(&self) -> Result<()> {
        let output = &self.opts.output;
        if self.block_copy() && matches!(WriteFormat::of(output), WriteFormat::Sam)
            || self.block_copy() && later_kind(output).is_some()
        {
            return Err(Error::Usage(format!(
                "OUTPUT={}: every input is a BAM, whose blocks are copied, and Picard writes BAM \
                 bytes under this name: name the OUTPUT .bam",
                output.display()
            )));
        }
        if let Some((kind, wave)) = later_kind(output) {
            return Err(Error::Usage(format!(
                "OUTPUT={}: bammetrics writes BAM and SAM only: {kind} output is wave {wave} of \
                 docs/07_full_parity_program.md",
                output.display()
            )));
        }
        let written = self.written();
        registry::check_writable(written.iter().map(PathBuf::as_path), &[])?;
        for input in self.inputs.iter().chain(&self.others) {
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
}

/// `gatherWithBlockCopying`: every input checked first (its end, its
/// header, its dictionary against the first's), then the copy.
fn gather_blocks(gather: &Gather) -> Result<()> {
    let opts = &gather.opts;
    let first = Source::open(&gather.inputs[0])?;
    for path in &gather.inputs[1..] {
        let source = Source::open(path)?;
        if source.start().has_records && source.start().references != first.start().references {
            return Err(Error::Input {
                path: path.clone(),
                source: io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "its sequence dictionary is not the one of the first INPUT {}: a block \
                         copy keeps each record's reference index, so its records would name \
                         other references",
                        first.path().display()
                    ),
                ),
            });
        }
    }
    // The sort order decides the index (htsjdk's `BAMIndexer` reads the
    // output header, the first input's); without one the header is not
    // parsed, as htsjdk's block copy never parses the first file's.
    let header = if opts.create_index {
        pass::read_header_at(first.path(), crate::sam_text::Stringency::Silent)?
    } else {
        sam::Header::default()
    };
    let output = |source| Error::Output {
        path: opts.output.clone(),
        source,
    };
    let file = AtomicFile::new(&opts.output);
    let mut out = bgzf_copy::writer(file.create().map_err(output)?);
    let naming = |e| match e {
        Error::Record(source) => output(source),
        other => other,
    };
    first
        .copy(&mut out, false, true, opts.compression_level)
        .map_err(naming)?;
    drop(first);
    for path in &gather.inputs[1..] {
        Source::open(path)?
            .copy(&mut out, true, true, opts.compression_level)
            .map_err(naming)?;
    }
    out.write_all(&bgzf_copy::EOF_BLOCK).map_err(output)?;
    drop(out.into_inner().map_err(|e| output(e.into_error()))?);
    let companions = Companions {
        index: opts.create_index,
        md5: opts.create_md5,
    };
    bam_writer::finish_file(
        file,
        &header,
        companions,
        &gather.index_path(),
        opts.threads,
    )?;
    Ok(())
}

/// The index in `to` of each reference of `from`, by name; `None` for a
/// name `to` lacks.
fn reference_map(from: &sam::Header, to: &sam::Header) -> Vec<Option<i32>> {
    from.reference_sequences()
        .keys()
        .map(|name| {
            to.reference_sequences()
                .get_index_of(name)
                .and_then(|i| i32::try_from(i).ok())
        })
        .collect()
}

/// `SAMRecord.setHeaderStrict`: the reference and mate reference of
/// `record` found again by name in the output header.
fn take_header(record: &mut RawRecord, map: &[Option<i32>], from: &sam::Header) -> Result<()> {
    let rename = |id: i32| -> Result<i32> {
        let Ok(at) = usize::try_from(id) else {
            return Ok(-1);
        };
        map.get(at).copied().flatten().ok_or_else(|| {
            let name = from
                .reference_sequences()
                .get_index(at)
                .map_or_else(|| "?".to_owned(), |(name, _)| name.to_string());
            Error::Metric(format!(
                "Reference index for '{name}' not found in sequence dictionary."
            ))
        })
    };
    let (reference, mate) = (record.reference_id(), record.mate_reference_id());
    let renamed = (rename(reference)?, rename(mate)?);
    if renamed != (reference, mate) {
        record.set_reference_ids(renamed.0, renamed.1);
    }
    Ok(())
}

/// `gatherNormally`: the first input's header, then every record of every
/// input into one presorted writer.
fn gather_records(gather: &Gather) -> Result<()> {
    let opts = &gather.opts;
    let first = &gather.inputs[0];
    let header = pass::read_header_at(first, opts.stringency)?;
    record_edit::check_sort_order(&header, opts.stringency, first)?;
    let out_header = record_edit::with_version(header, true);
    let mut out = Out::create(opts, &out_header)?;
    for path in &gather.inputs {
        pass::with_record_stream_kind(
            path,
            opts.reference.as_deref(),
            opts.threads,
            opts.stringency,
            |header, reader, kind| {
                let cram = kind == Kind::Cram;
                let mut restore =
                    Restore::new(cram, Mode::Rewrite, opts.reference.as_deref(), header)?;
                let map = reference_map(header, &out_header);
                let mut record = RawRecord::default();
                while record
                    .read_from(reader)
                    .map_err(|e| pass::read_error(path, cram, e))?
                {
                    restore.apply(&mut record)?;
                    take_header(&mut record, &map, header)?;
                    out.write(&mut record, Mode::Rewrite)?;
                }
                Ok(())
            },
        )?;
    }
    out.finish(&out_header, opts)
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let mut opts = record_edit::Options::from_invocation(inv)?;
    let (inputs, lists) = unroll(inv.all("INPUT"))?;
    if inputs.is_empty() {
        return Err(Error::Usage(
            "INPUT: the lists name no input to gather".to_owned(),
        ));
    }
    check_readable(&inputs)?;
    let others = lists
        .into_iter()
        .chain(opts.reference.iter().flat_map(|r| files_read(r)))
        .collect();
    opts.input = inputs[0].clone();
    let gather = Gather {
        inputs,
        others,
        opts,
    };
    gather.check_output()?;
    if gather.block_copy() {
        gather_blocks(&gather)?;
    } else {
        gather_records(&gather)?;
    }
    Ok(ToolOutput::NONE)
}

#[cfg(test)]
mod tests;

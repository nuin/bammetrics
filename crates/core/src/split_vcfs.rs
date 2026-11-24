//! Picard 3.5.0 `SplitVcfs`: the records of a VCF read as `VCFCodec`
//! reads them, each `INDEL` record written to `INDEL_OUTPUT` and each `SNP`
//! record to `SNP_OUTPUT`, both under the input's header, as htsjdk's
//! `VCFWriter` writes them ([`crate::vcf_writer`]).
//!
//! - The type is `VariantContext.getType` (`isIndel` is asked first). A
//!   record of any other type (`MNP`, `MIXED`, `SYMBOLIC`,
//!   `NO_VARIATION`) fails the run under `STRICT` (the default) with
//!   Picard's `Found a record with type MNP`, and is dropped without it.
//!   Picard logs `Found N records that didn't match SNP or INDEL` at
//!   `DEBUG`; so does this tool, under `VERBOSITY=DEBUG`.
//! - `CREATE_INDEX` (default true) needs a dictionary: the `##contig`
//!   lines, or `SEQUENCE_DICTIONARY` (a `.dict`, SAM, BAM or CRAM header,
//!   or an interval list), which only the index reads. Without one the run
//!   fails before any record, as Picard's does, even when no output takes
//!   an index.
//! - A `.vcf.gz` output (any name htsjdk writes as BGZF) gets its `.tbi`.
//!   A plain `.vcf` gets no index: Picard writes tribble's `.idx` there,
//!   which this writer does not (docs/07 decision W16-3). `CREATE_INDEX`
//!   left at its default then writes the VCF alone with a note on standard
//!   error; `CREATE_INDEX=true` given on the command line is refused.
//! - Both outputs and their companions are written under temporary names
//!   and renamed into place as one set after the last record: a failed run
//!   leaves none of them (Picard leaves both files written so far).

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::stage::Products;
use crate::stage::variant::{BuildVariantSink, VariantSink, VcfLine};
use crate::vcf::{self, Encoder, VariantType, VcfLines};
use crate::vcf_writer::{VcfFormat, VcfOptions, VcfWriter};

pub const TOOL: Tool = Tool {
    name: "SplitVcfs",
    summary: "writes the SNP records of a VCF to SNP_OUTPUT and its indel records to \
              INDEL_OUTPUT, both under the input header; STRICT fails on any other type. \
              Returns one row: the records written to each file and those dropped.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The VCF to read (plain, gzip or BGZF by its first bytes).",
        )
        .required(),
        Opt::path(
            &["SNP_OUTPUT"],
            "The VCF of the SNP records: block-compressed with a .tbi for a name ending in .gz, \
             else plain text.",
        )
        .written()
        .required(),
        Opt::path(
            &["INDEL_OUTPUT"],
            "The VCF of the indel records: block-compressed with a .tbi for a name ending in \
             .gz, else plain text.",
        )
        .written()
        .required(),
        Opt::path(
            &["SEQUENCE_DICTIONARY", "D"],
            "The dictionary of the index in place of the VCF's ##contig lines: a .dict, a SAM, \
             BAM or CRAM header, or an interval list.",
        ),
        Opt::bool(
            &["STRICT"],
            "Fail on a record that is neither a SNP nor an indel (default true); false drops it.",
        ),
        INDEX_OPTION,
        MD5_OPTION,
        COMPRESSION_OPTION,
        VERBOSITY_OPTION,
        Opt::rejected(&["OUTPUT", "O"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        registry::REFERENCE_NOT_DECODED,
        registry::STRINGENCY_NO_RECORDS,
    ],
    run: run_tool,
    collector: None,
};

/// `CREATE_INDEX` of a Picard VCF writer, default true.
pub(crate) const INDEX_OPTION: Opt = Opt::bool(
    &["CREATE_INDEX"],
    "Write a .tbi beside each .vcf.gz output (default true; a plain .vcf gets no index, and \
     true given for one is refused).",
);
/// `CREATE_MD5_FILE` of a Picard VCF writer.
pub(crate) const MD5_OPTION: Opt = Opt::bool(
    &["CREATE_MD5_FILE"],
    "Write OUTPUT.md5 beside each VCF (default false).",
);
/// `COMPRESSION_LEVEL` of a Picard VCF writer.
pub(crate) const COMPRESSION_OPTION: Opt = Opt::int(
    &["COMPRESSION_LEVEL"],
    "BGZF compression level of a .vcf.gz output, 0 to 9 (default 5).",
);
const VERBOSITY_OPTION: Opt = Opt::text(
    &["VERBOSITY"],
    "ERROR, WARNING, INFO (default) or DEBUG: DEBUG prints Picard's count of the records \
     dropped.",
);

/// The columns of the one row the tool returns.
pub const COLUMNS: [&str; 3] = ["snp_records", "indel_records", "other_records"];

/// The writer options of a Picard VCF tool: `CREATE_INDEX` (default true),
/// `CREATE_MD5_FILE`, `COMPRESSION_LEVEL`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WriterOptions {
    pub create_index: bool,
    /// `CREATE_INDEX` was given, not left at its default.
    pub index_given: bool,
    pub md5: bool,
    pub compression_level: u8,
    pub threads: usize,
}

impl WriterOptions {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let level: i64 = inv.parse("COMPRESSION_LEVEL", 5)?;
        let compression_level = u8::try_from(level)
            .ok()
            .filter(|level| *level <= 9)
            .ok_or_else(|| Error::Usage(format!("COMPRESSION_LEVEL={level}: expected 0 to 9")))?;
        Ok(Self {
            create_index: inv.flag("CREATE_INDEX", true)?,
            index_given: inv.get("CREATE_INDEX").is_some(),
            md5: inv.flag("CREATE_MD5_FILE", false)?,
            compression_level,
            threads: inv.threads.max(1),
        })
    }

    /// The options of the VCF written to `path` (named `option` in
    /// messages): a `.tbi` for a block-compressed name when the index is
    /// asked for. A plain name with `CREATE_INDEX` given true is refused
    /// (decision W16-3); with the default it gets no index and a note.
    pub fn for_output(&self, option: &str, path: &Path) -> Result<(VcfOptions, Option<String>)> {
        let plain = VcfFormat::of(path)? == VcfFormat::Plain;
        if plain && self.create_index && self.index_given {
            return Err(Error::Usage(format!(
                "{option}={}: CREATE_INDEX=true asks for tribble's .idx beside a plain VCF, \
                 which bammetrics does not write (docs/07 decision W16-3): name the output \
                 .vcf.gz for its .tbi, or give CREATE_INDEX=false",
                path.display()
            )));
        }
        let note = (plain && self.create_index).then(|| {
            format!(
                "{option}={}: no .idx is written beside a plain VCF (docs/07 decision W16-3); \
                 Picard writes one. A .vcf.gz name gets its .tbi",
                path.display()
            )
        });
        let options = VcfOptions {
            compression_level: self.compression_level,
            threads: self.threads,
            index: self.create_index && !plain,
            md5: self.md5,
        };
        Ok((options, note))
    }

    /// The files written for `path`.
    pub fn files(&self, path: &Path) -> Result<Vec<PathBuf>> {
        Ok(self.for_output("OUTPUT", path)?.0.files(path))
    }
}

/// A BCF input is refused before it is opened: the reader takes VCF text.
pub(crate) fn refuse_bcf(input: &Path) -> Result<()> {
    if input
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bcf"))
    {
        return Err(Error::Usage(format!(
            "{}: a BCF is not read; convert it with bcftools view -Ov first",
            input.display()
        )));
    }
    Ok(())
}

/// The `A sequence dictionary must be available` refusal of a Picard VCF
/// tool asked for an index with no dictionary.
pub(crate) fn no_dictionary(input: &Path) -> Error {
    Error::Usage(format!(
        "{}: A sequence dictionary must be available (either through the input file or by \
         setting it explicitly) when creating indexed output.",
        input.display()
    ))
}

/// Every option the tool takes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub input: PathBuf,
    pub snp_output: PathBuf,
    pub indel_output: PathBuf,
    pub sequence_dictionary: Option<PathBuf>,
    pub strict: bool,
    pub writer: WriterOptions,
    /// `VERBOSITY=DEBUG`: the count of dropped records is printed.
    pub debug: bool,
    pub quiet: bool,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let path = |name: &str| inv.get(name).map(PathBuf::from);
        let required =
            |name: &str| path(name).ok_or_else(|| Error::Usage(format!("{name} is required")));
        let debug = match inv.get("VERBOSITY") {
            None | Some("ERROR" | "WARNING" | "INFO") => false,
            Some("DEBUG") => true,
            Some(other) => {
                return Err(Error::Usage(format!(
                    "VERBOSITY={other}: one of ERROR, WARNING, INFO, DEBUG"
                )));
            }
        };
        Ok(Self {
            input: inv.input.clone(),
            snp_output: required("SNP_OUTPUT")?,
            indel_output: required("INDEL_OUTPUT")?,
            sequence_dictionary: path("SEQUENCE_DICTIONARY"),
            strict: inv.flag("STRICT", true)?,
            writer: WriterOptions::from_invocation(inv)?,
            debug,
            quiet: inv.quiet,
        })
    }
}

/// What a run wrote: the records of each output and those dropped.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub snps: u64,
    pub indels: u64,
    pub others: u64,
}

/// Where a record goes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Route {
    Snp,
    Indel,
    /// Neither: dropped, or a failure under `STRICT`.
    Other(VariantType),
}

/// The output of a record of type `kind`: `isIndel` first, then `isSNP`.
pub fn route(kind: VariantType) -> Route {
    match kind {
        VariantType::Indel => Route::Indel,
        VariantType::Snp => Route::Snp,
        other => Route::Other(other),
    }
}

/// The failure of a record that is neither a SNP nor an indel under
/// `STRICT`, with Picard's text.
pub fn strict_failure(kind: VariantType) -> Error {
    Error::Metric(format!("Found a record with type {}", kind.as_str()))
}

/// Picard's `DEBUG` line for `others` dropped records.
pub fn debug_line(others: u64) -> String {
    format!("DEBUG\tSplitVcfs\tFound {others} records that didn't match SNP or INDEL")
}

/// The two writers of one split, created before the first record as
/// Picard creates them.
pub struct Writers {
    snp: VcfWriter,
    indel: VcfWriter,
    strict: bool,
    counts: Counts,
}

impl Writers {
    /// Check the outputs against each other and the inputs, and start both
    /// under `header`. `notes` gets the note of each plain output that
    /// takes no index.
    pub fn create(
        snp: &Path,
        indel: &Path,
        header: &Arc<vcf::VcfHeader>,
        writer: &WriterOptions,
        dictionary: Option<vcf::Dictionary>,
        strict: bool,
        notes: &mut Vec<String>,
    ) -> Result<Self> {
        let (snp_options, snp_note) = writer.for_output("SNP_OUTPUT", snp)?;
        let (indel_options, indel_note) = writer.for_output("INDEL_OUTPUT", indel)?;
        refuse_shared_files(&snp_options.files(snp), &indel_options.files(indel))?;
        if writer.create_index && dictionary.is_none() && header.dictionary().is_none() {
            return Err(no_dictionary(snp));
        }
        notes.extend(snp_note);
        notes.extend(indel_note);
        let snp =
            VcfWriter::create_indexed_by(snp, Arc::clone(header), snp_options, dictionary.clone())?;
        let indel =
            VcfWriter::create_indexed_by(indel, Arc::clone(header), indel_options, dictionary)?;
        Ok(Self {
            snp,
            indel,
            strict,
            counts: Counts::default(),
        })
    }

    /// Write `line`, whose type `encoder` reads, to its output.
    pub fn write(&mut self, encoder: &Encoder, line: &str) -> Result<()> {
        let read = encoder.read_type(line).map_err(|message| {
            Error::Metric(format!(
                "record {}: {message}",
                self.counts.snps + self.counts.indels + self.counts.others + 1
            ))
        })?;
        let Some(kind) = read else {
            return Ok(());
        };
        match route(kind) {
            Route::Indel => {
                self.counts.indels += 1;
                self.indel.write_line(line)
            }
            Route::Snp => {
                self.counts.snps += 1;
                self.snp.write_line(line)
            }
            Route::Other(kind) if self.strict => Err(strict_failure(kind)),
            Route::Other(_) => {
                self.counts.others += 1;
                Ok(())
            }
        }
    }

    pub fn counts(&self) -> Counts {
        self.counts
    }

    /// Close both, then rename every file into place as one set. The
    /// paths written, each output's VCF first.
    pub fn finish(self) -> Result<Vec<PathBuf>> {
        let snp = self.snp.stage()?;
        let indel = self.indel.stage()?;
        let mut written = snp.written;
        written.extend(indel.written);
        let mut files = snp.files;
        files.extend(indel.files);
        crate::split_sam::commit_all(files)?;
        Ok(written)
    }
}

/// Two outputs on one file: the SNP and indel files and their companions
/// must differ by name and by device and inode.
fn refuse_shared_files(snp: &[PathBuf], indel: &[PathBuf]) -> Result<()> {
    for a in snp {
        if let Some(b) = indel
            .iter()
            .find(|b| *b == a || crate::mark_duplicates::same_file(a, b))
        {
            return Err(Error::Usage(format!(
                "SNP_OUTPUT and INDEL_OUTPUT both write {} ({}): give each its own file",
                a.display(),
                b.display()
            )));
        }
    }
    Ok(())
}

/// The dictionary `SEQUENCE_DICTIONARY` names: any file
/// [`vcf::read_dictionary`] reads but a VCF, which Picard's SAM reader
/// cannot open.
fn read_sequence_dictionary(path: &Path) -> Result<vcf::Dictionary> {
    if vcf::is_vcf(path) {
        return Err(Error::Usage(format!(
            "SEQUENCE_DICTIONARY={}: a VCF is no dictionary here; give a .dict, a SAM, BAM or \
             CRAM header, or an interval list",
            path.display()
        )));
    }
    vcf::read_dictionary(path)
}

/// Run the tool: every record of `opts.input` to its output.
pub fn split(opts: &Options) -> Result<(Counts, Vec<PathBuf>)> {
    refuse_bcf(&opts.input)?;
    let mut inputs = vec![opts.input.clone()];
    inputs.extend(opts.sequence_dictionary.iter().cloned());
    for input in &inputs {
        if !input.exists() {
            return Err(Error::Input {
                path: input.clone(),
                source: std::io::ErrorKind::NotFound.into(),
            });
        }
    }
    let mut targets = opts.writer.files(&opts.snp_output)?;
    targets.extend(opts.writer.files(&opts.indel_output)?);
    registry::check_writable(targets.iter().map(PathBuf::as_path), &inputs)?;
    let mut lines = VcfLines::open(&opts.input)?;
    let header = lines.shared_header();
    let dictionary = opts
        .sequence_dictionary
        .as_deref()
        .map(read_sequence_dictionary)
        .transpose()?;
    let mut notes = Vec::new();
    let mut writers = Writers::create(
        &opts.snp_output,
        &opts.indel_output,
        &header,
        &opts.writer,
        dictionary,
        opts.strict,
        &mut notes,
    )?;
    if !opts.quiet {
        for note in &notes {
            eprintln!("bammetrics: {note}");
        }
    }
    let encoder = Encoder::new(&header);
    let mut line = String::new();
    let place = |e: Error| match e {
        Error::Metric(message) => Error::Metric(format!("{}: {message}", opts.input.display())),
        other => other,
    };
    while lines.next_line(&mut line)? {
        writers.write(&encoder, &line).map_err(place)?;
    }
    let counts = writers.counts();
    if opts.debug && counts.others > 0 {
        eprintln!("{}", debug_line(counts.others));
    }
    let written = writers.finish()?;
    Ok((counts, written))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let (counts, _) = split(&opts)?;
    let rows = vec![vec![
        Value::from(counts.snps as i64),
        Value::from(counts.indels as i64),
        Value::from(counts.others as i64),
    ]];
    Ok(MetricsOutput::text(&COLUMNS, rows, String::new()).into())
}

/// The tool as a split of a variant stream (docs/chaining.md section 34):
/// the records of its position through [`Writers`], both files renamed
/// into place as one set when the stream ends. The chain names
/// `SNP_OUTPUT` and `INDEL_OUTPUT` when the step does not.
pub fn variant_split(inv: &Invocation) -> Result<BuildVariantSink> {
    let opts = Options::from_invocation(inv)?;
    refuse_shared_files(
        &opts.writer.files(&opts.snp_output)?,
        &opts.writer.files(&opts.indel_output)?,
    )?;
    Ok(Box::new(move |header: &Arc<vcf::VcfHeader>| {
        let dictionary = opts
            .sequence_dictionary
            .as_deref()
            .map(read_sequence_dictionary)
            .transpose()?;
        let mut notes = Vec::new();
        let writers = Writers::create(
            &opts.snp_output,
            &opts.indel_output,
            header,
            &opts.writer,
            dictionary,
            opts.strict,
            &mut notes,
        )?;
        Ok(Box::new(SplitSink {
            encoder: Encoder::new(header),
            writers,
            notes,
        }) as Box<dyn VariantSink>)
    }))
}

/// The sink of [`variant_split`].
struct SplitSink {
    encoder: Encoder,
    writers: Writers,
    notes: Vec<String>,
}

impl VariantSink for SplitSink {
    fn accept_batch(&mut self, _header: &vcf::VcfHeader, batch: &[VcfLine]) -> Result<()> {
        batch
            .iter()
            .try_for_each(|line| self.writers.write(&self.encoder, &line.text))
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        let counts = self.writers.counts();
        let mut notes = self.notes;
        if counts.others > 0 {
            notes.push(format!(
                "{} records neither SNP nor indel were dropped (STRICT false)",
                counts.others
            ));
        }
        let files = self.writers.finish()?;
        Ok(Products {
            files,
            notes,
            ..Products::default()
        })
    }
}

/// The files a split of `snp` and `indel` writes: each VCF, then its
/// `.tbi` and `.md5` when asked for.
pub fn output_files(writer: &WriterOptions, snp: &Path, indel: &Path) -> Result<Vec<PathBuf>> {
    let mut files = writer.files(snp)?;
    files.extend(writer.files(indel)?);
    Ok(files)
}

#[cfg(test)]
mod tests;

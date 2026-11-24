//! Picard 3.5.0 `MakeSitesOnlyVcf`: every record of a VCF with its
//! genotypes dropped, or cut to the samples of `SAMPLE`
//! ([`crate::vcf::SampleSubset`]), the site columns and `INFO` as read
//! (`AN`, `AF` and the like stay as the genotypes made them), written as
//! htsjdk's `VCFWriter` writes them ([`crate::vcf_writer`]).
//!
//! - The header keeps every line of the input; its sample columns are the
//!   names of `SAMPLE` in `TreeSet` order, none by default.
//! - `CREATE_INDEX` (default true) needs the input's `##contig` lines,
//!   else the run fails before the first record, as Picard's does. A
//!   `.vcf.gz` `OUTPUT` gets its `.tbi`; a plain `.vcf` gets no index
//!   (docs/07 decision W16-3; [`crate::split_vcfs::WriterOptions`]).
//! - `OUTPUT` and its companions are renamed into place after the last
//!   record; a failed run leaves none of them.

use std::path::PathBuf;
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::split_vcfs::{self, WriterOptions};
use crate::stage::Counters;
use crate::stage::variant::{BuildVariantTransform, VariantTransform, VcfLine};
use crate::vcf::{SampleSubset, VcfHeader, VcfLines};
use crate::vcf_writer::VcfWriter;

pub const TOOL: Tool = Tool {
    name: "MakeSitesOnlyVcf",
    summary: "writes OUTPUT: every record of the VCF without its genotypes, or with those of \
              the SAMPLE names only, the site columns and INFO as read. Returns one row: the \
              records written.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The VCF to read (plain, gzip or BGZF by its first bytes).",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The VCF to write: block-compressed with a .tbi for a name ending in .gz, else \
             plain text.",
        )
        .written()
        .required(),
        Opt::text(
            &["SAMPLE", "S"],
            "A sample whose genotypes stay, in the header in sorted order; a name the VCF lacks \
             gets no-calls. Repeat for each; null clears. None by default.",
        )
        .repeated(),
        split_vcfs::INDEX_OPTION,
        split_vcfs::MD5_OPTION,
        split_vcfs::COMPRESSION_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        registry::REFERENCE_NOT_DECODED,
        registry::STRINGENCY_NO_RECORDS,
    ],
    run: run_tool,
    collector: None,
};

/// The columns of the one row the tool returns.
pub const COLUMNS: [&str; 1] = ["records"];

/// Every option the tool takes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub input: PathBuf,
    pub output: PathBuf,
    /// `SAMPLE` after `null` handling, as given.
    pub samples: Vec<String>,
    pub writer: WriterOptions,
    pub quiet: bool,
}

/// `SAMPLE`: each value added, `null` clearing the list so far.
pub fn samples_of(inv: &Invocation) -> Result<Vec<String>> {
    registry::picard_list(Vec::new(), inv.all("SAMPLE"), |name| Ok(name.to_owned()))
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let output = inv
            .get("OUTPUT")
            .map(PathBuf::from)
            .ok_or_else(|| Error::Usage("OUTPUT is required".to_owned()))?;
        Ok(Self {
            input: inv.input.clone(),
            output,
            samples: samples_of(inv)?,
            writer: WriterOptions::from_invocation(inv)?,
            quiet: inv.quiet,
        })
    }
}

/// Run the tool: the records written.
pub fn make_sites_only(opts: &Options) -> Result<u64> {
    split_vcfs::refuse_bcf(&opts.input)?;
    if !opts.input.exists() {
        return Err(Error::Input {
            path: opts.input.clone(),
            source: std::io::ErrorKind::NotFound.into(),
        });
    }
    let (options, note) = opts.writer.for_output("OUTPUT", &opts.output)?;
    let files = options.files(&opts.output);
    registry::check_writable(
        files.iter().map(PathBuf::as_path),
        std::slice::from_ref(&opts.input),
    )?;
    let mut lines = VcfLines::open(&opts.input)?;
    if opts.writer.create_index && lines.header().dictionary().is_none() {
        return Err(split_vcfs::no_dictionary(&opts.input));
    }
    if let Some(note) = note.filter(|_| !opts.quiet) {
        eprintln!("bammetrics: {note}");
    }
    let subset = SampleSubset::new(lines.header(), &opts.samples);
    let mut writer = VcfWriter::create_shared(&opts.output, subset.shared_header(), options)?;
    let mut line = String::new();
    let mut edited = String::new();
    let mut read = 0u64;
    while lines.next_line(&mut line)? {
        read += 1;
        let kept = subset.edit(&line, &mut edited).map_err(|message| {
            Error::Metric(format!(
                "{}: record {read}: {message}",
                opts.input.display()
            ))
        })?;
        if kept {
            writer.write_line(&edited)?;
        }
    }
    let records = writer.records();
    writer.finish()?;
    Ok(records)
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let records = make_sites_only(&opts)?;
    let rows = vec![vec![Value::from(records as i64)]];
    Ok(MetricsOutput::text(&COLUMNS, rows, String::new()).into())
}

/// The options a transform of a chain leaves to its `write` step.
const WRITE_STEP_OPTIONS: [&str; 3] = ["OUTPUT", "CREATE_INDEX", "COMPRESSION_LEVEL"];

/// The tool as a transform of a variant stream (docs/chaining.md section
/// 34): each record cut as the tool cuts it, the header's samples those of
/// `SAMPLE`. The VCF is the chain's `write` step's.
pub fn variant_transform(inv: &Invocation) -> Result<BuildVariantTransform> {
    if let Some(name) = WRITE_STEP_OPTIONS.iter().find(|n| inv.get(n).is_some()) {
        return Err(Error::Usage(format!(
            "{name} belongs to the VCF a write step writes: a transform of a chain writes \
             none (-o r.vcf.gz or write PATH)"
        )));
    }
    let samples = samples_of(inv)?;
    Ok(Box::new(move |header: &Arc<VcfHeader>| {
        Ok(Box::new(SitesOnly {
            subset: SampleSubset::new(header, &samples),
            counters: Counters::default(),
            edited: String::new(),
        }) as Box<dyn VariantTransform>)
    }))
}

/// The transform of [`variant_transform`].
struct SitesOnly {
    subset: SampleSubset,
    counters: Counters,
    edited: String,
}

impl VariantTransform for SitesOnly {
    fn header(&self, _up: &Arc<VcfHeader>) -> Result<Arc<VcfHeader>> {
        Ok(self.subset.shared_header())
    }

    fn apply(&mut self, batch: &mut Vec<VcfLine>) -> Result<()> {
        let mut failure = None;
        let counters = &mut self.counters;
        batch.retain_mut(|line| {
            if failure.is_some() {
                return false;
            }
            counters.records_in += 1;
            match self.subset.edit(&line.text, &mut self.edited) {
                Ok(true) => {
                    std::mem::swap(&mut line.text, &mut self.edited);
                    counters.records_out += 1;
                    counters.changed += 1;
                    true
                }
                Ok(false) => false,
                Err(message) => {
                    failure = Some(Error::Metric(format!(
                        "record {}: {message}",
                        counters.records_in
                    )));
                    false
                }
            }
        });
        failure.map_or(Ok(()), Err)
    }

    fn counters(&self) -> Counters {
        self.counters
    }
}

#[cfg(test)]
mod tests;

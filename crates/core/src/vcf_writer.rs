//! The VCF writer (wave W16): records written as htsjdk 5.0.0's `VCFWriter`
//! writes them for Picard 3.5.0, the header through
//! [`crate::vcf::VcfHeader`] and each record line through
//! [`crate::vcf::Encoder`], one line in memory at a time.
//!
//! - The kind follows the name as `VariantContextWriterBuilder` decides it:
//!   `.bcf` is BCF (refused: no BCF writer); a name ending in `.gz`,
//!   `.gzip`, `.bgz` or `.bgzf` in any case is block-compressed text
//!   (BGZF on worker threads, its end-of-file block last); any other name
//!   is plain text.
//! - The index of a block-compressed VCF is the `.tbi` beside it
//!   (`out.vcf.gz.tbi`), built as `TabixIndexCreator` builds it
//!   ([`tabix`]) from the finished file read back, so each record's offset
//!   is where the file holds it. The records must be sorted by contig and
//!   start, and the header needs `##contig` lines, as Picard's VCF writers
//!   ask (`A sequence dictionary must be available in the input file when
//!   creating indexed output`). Picard also writes a tribble `.idx` beside a plain `.vcf`
//!   (`DynamicIndexCreator`); this writer refuses that index (docs/07
//!   decision W16-3) and writes the `.vcf` alone.
//! - `md5` writes `OUTPUT.md5`, the hex MD5 of the file.
//! - The text is written as ISO-8859-1, as htsjdk's `VCFWriter` writes it
//!   (`VCFEncoder.VCF_CHARSET`): a character up to U+00FF as its one byte,
//!   any other as `?`. The reader decodes UTF-8, a malformed byte becoming
//!   U+FFFD, so it too is written `?`.
//! - Every file is written under a temporary name in its directory and
//!   renamed into place when the last one is whole
//!   ([`crate::cleanup::AtomicFile`]), the index and MD5 first and the VCF
//!   last, as one set ([`crate::split_sam::commit_all`]): a failed run or
//!   rename leaves none of them.

use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles::bgzf;

use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::java::iso_8859_1;
use crate::vcf::{Encoder, VcfHeader};

mod tabix;
#[cfg(test)]
mod tests;

/// How a VCF is written, from its name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VcfFormat {
    /// Text.
    Plain,
    /// BGZF-compressed text.
    Bgzf,
}

impl VcfFormat {
    /// The kind of `path` (`determineOutputTypeFromFile`): a BCF name is a
    /// usage error.
    pub fn of(path: &Path) -> Result<Self> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        if name.ends_with(".bcf") {
            return Err(Error::Usage(format!(
                "{}: a BCF is not written; name a .vcf or .vcf.gz",
                path.display()
            )));
        }
        let lower = name.to_ascii_lowercase();
        Ok(
            if [".gz", ".gzip", ".bgz", ".bgzf"]
                .iter()
                .any(|e| lower.ends_with(e))
            {
                Self::Bgzf
            } else {
                Self::Plain
            },
        )
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Plain => "vcf",
            Self::Bgzf => "vcf.gz",
        }
    }
}

/// `Tribble.tabixIndexFile`: the `.tbi` beside a block-compressed VCF.
pub fn index_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_owned();
    name.push(".tbi");
    PathBuf::from(name)
}

/// Where `CREATE_MD5_FILE` writes: the VCF's path with `.md5` appended.
pub fn md5_path(output: &Path) -> PathBuf {
    crate::bam_writer::md5_path(output)
}

/// What [`VcfWriter`] writes and beside it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct VcfOptions {
    /// BGZF compression level, 0 to 9 (Picard's `COMPRESSION_LEVEL`
    /// default 5).
    pub compression_level: u8,
    /// BGZF workers.
    pub threads: usize,
    /// The `.tbi` of a block-compressed VCF.
    pub index: bool,
    pub md5: bool,
}

impl Default for VcfOptions {
    fn default() -> Self {
        Self {
            compression_level: 5,
            threads: 1,
            index: false,
            md5: false,
        }
    }
}

impl VcfOptions {
    /// A usage error for options `path` cannot take: an index of a plain
    /// VCF (decision W16-3).
    pub fn check(&self, path: &Path) -> Result<()> {
        if self.compression_level > 9 {
            return Err(Error::Usage(format!(
                "compression level {} is not 0 to 9",
                self.compression_level
            )));
        }
        if self.index && VcfFormat::of(path)? == VcfFormat::Plain {
            return Err(Error::Usage(plain_index_refusal(path)));
        }
        Ok(())
    }

    /// The files written for `path`: the VCF, then its index and MD5 when
    /// asked for.
    pub fn files(&self, path: &Path) -> Vec<PathBuf> {
        let mut files = vec![path.to_path_buf()];
        if self.index {
            files.push(index_path(path));
        }
        if self.md5 {
            files.push(md5_path(path));
        }
        files
    }
}

/// The refusal of an index beside a plain VCF.
pub fn plain_index_refusal(path: &Path) -> String {
    format!(
        "{}: the index of a plain VCF is tribble's .idx, which this writer does not write \
         (decision W16-3); name the output .vcf.gz for its .tbi, or leave the index out",
        path.display()
    )
}

/// The open output.
enum Output {
    Plain(BufWriter<File>),
    Bgzf(bgzf::io::MultithreadedWriter<BufWriter<File>>),
}

impl Output {
    fn writer(&mut self) -> &mut dyn Write {
        match self {
            Self::Plain(out) => out,
            Self::Bgzf(out) => out,
        }
    }

    fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut out) => out.flush(),
            Self::Bgzf(mut out) => out.finish()?.flush(),
        }
    }
}

/// A VCF being written under its temporary name.
pub struct VcfWriter {
    path: PathBuf,
    file: AtomicFile,
    out: Output,
    format: VcfFormat,
    header: Arc<VcfHeader>,
    encoder: Encoder,
    options: VcfOptions,
    /// The dictionary of the `.tbi` when it is not the header's
    /// (`SplitVcfs SEQUENCE_DICTIONARY`).
    dictionary: Option<crate::vcf::Dictionary>,
    line: String,
    records: u64,
}

/// A closed VCF and its companions under their temporary names: the paths
/// they take, the VCF's first, and the files to rename into place.
pub struct Staged {
    pub written: Vec<PathBuf>,
    pub files: Vec<AtomicFile>,
}

impl VcfWriter {
    /// Start `path` under its temporary name and write `header` as htsjdk
    /// writes it.
    pub fn create(path: &Path, header: &VcfHeader, options: VcfOptions) -> Result<Self> {
        Self::create_shared(path, Arc::new(header.clone()), options)
    }

    /// [`VcfWriter::create`] holding the caller's copy of the header.
    pub fn create_shared(path: &Path, header: Arc<VcfHeader>, options: VcfOptions) -> Result<Self> {
        Self::create_indexed_by(path, header, options, None)
    }

    /// [`VcfWriter::create_shared`] whose `.tbi` takes its contig lengths
    /// from `dictionary` when one is given (Picard's
    /// `setReferenceDictionary`), else from the header's `##contig` lines.
    pub fn create_indexed_by(
        path: &Path,
        header: Arc<VcfHeader>,
        options: VcfOptions,
        dictionary: Option<crate::vcf::Dictionary>,
    ) -> Result<Self> {
        options.check(path)?;
        if options.index && dictionary.is_none() && header.dictionary().is_none() {
            // Picard's VCF writers refuse it before the first record.
            return Err(Error::Usage(format!(
                "{}: A sequence dictionary must be available in the input file when creating \
                 indexed output (the header has no ##contig line)",
                path.display()
            )));
        }
        let format = VcfFormat::of(path)?;
        let output = |source| Error::Output {
            path: path.to_path_buf(),
            source,
        };
        let file = AtomicFile::new(path);
        let handle = BufWriter::with_capacity(1 << 20, file.create().map_err(output)?);
        let mut out = match format {
            VcfFormat::Plain => Output::Plain(handle),
            VcfFormat::Bgzf => {
                let level = bgzf::io::writer::CompressionLevel::try_from(options.compression_level)
                    .map_err(|e| output(io::Error::new(io::ErrorKind::InvalidInput, e)))?;
                let workers = NonZero::new(options.threads).unwrap_or(NonZero::<usize>::MIN);
                Output::Bgzf(
                    bgzf::io::multithreaded_writer::Builder::default()
                        .set_compression_level(level)
                        .set_worker_count(workers)
                        .build_from_writer(handle),
                )
            }
        };
        let writer = out.writer();
        header
            .each_line(&mut |piece| writer.write_all(&iso_8859_1(piece)))
            .map_err(output)?;
        Ok(Self {
            path: path.to_path_buf(),
            file,
            out,
            format,
            encoder: Encoder::new(&header),
            header,
            options,
            dictionary,
            line: String::new(),
            records: 0,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The records written so far.
    pub fn records(&self) -> u64 {
        self.records
    }

    /// Write one record line (no line end) as the encoder writes it. A
    /// line that starts with `#` is skipped, as the codec skips it.
    pub fn write_line(&mut self, line: &str) -> Result<()> {
        self.line.clear();
        let encoded = self
            .encoder
            .encode(&self.header, line, &mut self.line)
            .map_err(|message| {
                Error::Metric(format!(
                    "{}: record {}: {message}",
                    self.path.display(),
                    self.records + 1
                ))
            })?;
        if encoded.is_none() {
            return Ok(());
        }
        self.line.push('\n');
        self.records += 1;
        self.out
            .writer()
            .write_all(&iso_8859_1(&self.line))
            .map_err(|source| Error::Output {
                path: self.path.clone(),
                source,
            })
    }

    /// Close the VCF, write its index and MD5 from the finished file, then
    /// rename every file into place. The VCF's path first.
    pub fn finish(self) -> Result<Vec<PathBuf>> {
        let Staged { written, files } = self.stage()?;
        crate::split_sam::commit_all(files)?;
        Ok(written)
    }

    /// Close the VCF and write its index and MD5 from the finished file,
    /// every file left under its temporary name for the caller to rename
    /// with others as one set ([`crate::split_sam::commit_all`]).
    pub fn stage(self) -> Result<Staged> {
        let output = |path: &Path| {
            let path = path.to_path_buf();
            move |source| Error::Output {
                path: path.clone(),
                source,
            }
        };
        self.out.finish().map_err(output(&self.path))?;
        let index = if self.options.index && self.format == VcfFormat::Bgzf {
            Some(write_tabix(
                self.file.writing(),
                &self.path,
                self.dictionary.as_ref().or(self.header.dictionary()),
                &index_path(&self.path),
            )?)
        } else {
            None
        };
        let md5 = if self.options.md5 {
            Some(crate::bam_writer::write_md5(
                self.file.writing(),
                &md5_path(&self.path),
            )?)
        } else {
            None
        };
        // The companions first and the VCF last, as one set: a failed
        // rename leaves none of them and puts back each file it replaced.
        let mut staged: Vec<AtomicFile> = index.into_iter().chain(md5).collect();
        let mut written = vec![self.path.clone()];
        written.extend(staged.iter().map(|file| file.path().to_path_buf()));
        staged.push(self.file);
        Ok(Staged {
            written,
            files: staged,
        })
    }
}

/// The `.tbi` of the block-compressed VCF at `vcf` (renamed to `shown`
/// later), into a file to be renamed to `target`. Each record line is read
/// back for `CHROM`, `POS` and its end ([`tabix::record_span`]); the writer
/// checked the rest.
fn write_tabix(
    vcf: &Path,
    shown: &Path,
    dictionary: Option<&crate::vcf::Dictionary>,
    target: &Path,
) -> Result<AtomicFile> {
    let read = |source| Error::Input {
        path: shown.to_path_buf(),
        source,
    };
    let failed = |message: String| Error::Metric(format!("{}: {message}", shown.display()));
    let mut reader = bgzf::io::Reader::new(io::BufReader::new(File::open(vcf).map_err(read)?));
    let mut builder = tabix::TabixBuilder::new(dictionary);
    let mut line = Vec::new();
    let mut end = 0u64;
    loop {
        if reader.fill_buf().map_err(read)?.is_empty() {
            break;
        }
        let start = u64::from(reader.virtual_position());
        line.clear();
        reader.read_until(b'\n', &mut line).map_err(read)?;
        end = u64::from(reader.virtual_position());
        if line.first() == Some(&b'#') {
            continue;
        }
        if line.last() == Some(&b'\n') {
            line.pop();
        }
        let (contig, first, last) = tabix::record_span(&line).map_err(failed)?;
        builder.add(contig, first, last, start).map_err(failed)?;
    }
    let bytes = builder.finish(end).map_err(failed)?;
    let file = AtomicFile::new(target);
    let written = |source| Error::Output {
        path: target.to_path_buf(),
        source,
    };
    let out = BufWriter::new(file.create().map_err(written)?);
    tabix::write_compressed(out, &bytes).map_err(written)?;
    Ok(file)
}

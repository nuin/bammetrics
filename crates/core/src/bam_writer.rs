//! The BAM file a tool or a chain writes: header and [`RawRecord`]s into
//! BGZF blocks compressed on worker threads, under a temporary name
//! ([`AtomicFile`]) until [`finish`] renames it into place with its index
//! and MD5. `MarkDuplicates` and the `write` step of a chain
//! ([`crate::stage`]) share it, so both write the same bytes for the same
//! header and records.
//!
//! The header is written byte for byte, as htsjdk writes it; noodles' own
//! writer refuses bytes outside printable ASCII, which a valid BAM may hold
//! (`@RG SM:s\xc3\xa9`, a tab in `@CO`).

use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::bgzf;
use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::mark_duplicates::SortOrder;

/// The fields of one header record as `(tag, value)` byte pairs.
fn fields<'a, K, V>(
    other: impl IntoIterator<Item = (&'a K, &'a V)>,
) -> impl Iterator<Item = (&'a [u8], &'a [u8])>
where
    K: AsRef<[u8; 2]> + 'a,
    V: AsRef<[u8]> + 'a,
{
    other
        .into_iter()
        .map(|(t, v)| (t.as_ref().as_slice(), v.as_ref()))
}

/// One header line: `kind`, then a tab and `TAG:value` per field.
fn header_line<'a>(
    text: &mut Vec<u8>,
    kind: &[u8],
    fields: impl IntoIterator<Item = (&'a [u8], &'a [u8])>,
) {
    text.extend_from_slice(kind);
    for (tag, value) in fields {
        text.push(b'\t');
        text.extend_from_slice(tag);
        text.push(b':');
        text.extend_from_slice(value);
    }
    text.push(b'\n');
}

/// The SAM header text: `@HD`, `@SQ`, `@RG`, `@PG` and `@CO` lines in
/// that order, each record's fields in the order noodles keeps them, the
/// order noodles' own writer uses. Values are written byte for byte.
pub fn header_bytes(header: &sam::Header) -> Vec<u8> {
    let mut text = Vec::new();
    if let Some(map) = header.header() {
        let version = map.version();
        let version = format!("{}.{}", version.major(), version.minor());
        let vn = std::iter::once((b"VN".as_slice(), version.as_bytes()));
        header_line(&mut text, b"@HD", vn.chain(fields(map.other_fields())));
    }
    for (name, map) in header.reference_sequences() {
        let length = map.length().to_string();
        let sn_ln = [
            (b"SN".as_slice(), name.as_ref()),
            (b"LN".as_slice(), length.as_bytes()),
        ];
        header_line(
            &mut text,
            b"@SQ",
            sn_ln.into_iter().chain(fields(map.other_fields())),
        );
    }
    for (id, map) in header.read_groups() {
        let id = std::iter::once((b"ID".as_slice(), id.as_ref()));
        header_line(&mut text, b"@RG", id.chain(fields(map.other_fields())));
    }
    for (id, map) in header.programs().as_ref() {
        let id = std::iter::once((b"ID".as_slice(), id.as_ref()));
        header_line(&mut text, b"@PG", id.chain(fields(map.other_fields())));
    }
    for comment in header.comments() {
        text.extend_from_slice(b"@CO\t");
        text.extend_from_slice(comment);
        text.push(b'\n');
    }
    text
}

/// The BAM header: magic, the text, the reference list.
fn write_bam_header(out: &mut impl Write, header: &sam::Header) -> io::Result<()> {
    let invalid = |e| io::Error::new(io::ErrorKind::InvalidInput, e);
    let text = header_bytes(header);
    out.write_all(b"BAM\x01")?;
    out.write_all(&i32::try_from(text.len()).map_err(invalid)?.to_le_bytes())?;
    out.write_all(&text)?;
    let references = header.reference_sequences();
    out.write_all(
        &i32::try_from(references.len())
            .map_err(invalid)?
            .to_le_bytes(),
    )?;
    for (name, map) in references {
        let length = u32::try_from(name.len() + 1).map_err(invalid)?;
        out.write_all(&length.to_le_bytes())?;
        out.write_all(name)?;
        out.write_all(&[0])?;
        let length = i32::try_from(usize::from(map.length())).map_err(invalid)?;
        out.write_all(&length.to_le_bytes())?;
    }
    Ok(())
}

/// The BGZF output, compressed on worker threads. Records must come in the
/// header's order, as htsjdk's presorted `SAMFileWriterImpl` asserts:
/// coordinate (reference, then start) or query name; any other order is
/// not checked.
pub struct BamWriter {
    path: PathBuf,
    /// The output under its temporary name.
    file: AtomicFile,
    inner: noodles::bam::io::Writer<bgzf::io::MultithreadedWriter<BufWriter<File>>>,
    order: SortOrder,
    /// The BGZF workers, which also read the file back for its index.
    threads: usize,
    /// Reference and start of the last record written, and its name.
    previous: Option<(i32, i32)>,
    previous_name: Vec<u8>,
}

impl BamWriter {
    /// Start `path` under its temporary name and write `header`.
    /// `compression_level` is 0 to 9; `threads` BGZF workers compress.
    pub fn create(
        path: &Path,
        compression_level: u8,
        threads: usize,
        header: &sam::Header,
    ) -> Result<Self> {
        let path = path.to_path_buf();
        let error = |source| Error::Output {
            path: path.clone(),
            source,
        };
        let atomic = AtomicFile::new(&path);
        let file = atomic.create().map_err(error)?;
        let level = bgzf::io::writer::CompressionLevel::try_from(compression_level)
            .map_err(|e| error(io::Error::new(io::ErrorKind::InvalidInput, e)))?;
        let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        let bgzf = bgzf::io::multithreaded_writer::Builder::default()
            .set_compression_level(level)
            .set_worker_count(workers)
            .build_from_writer(BufWriter::with_capacity(1 << 20, file));
        let mut inner = noodles::bam::io::Writer::from(bgzf);
        write_bam_header(inner.get_mut(), header).map_err(error)?;
        Ok(Self {
            path,
            file: atomic,
            inner,
            order: SortOrder::of(header),
            threads,
            previous: None,
            previous_name: Vec::new(),
        })
    }

    /// `assertPresorted`.
    fn check_order(&mut self, record: &RawRecord) -> Result<()> {
        if let Some((reference, start)) = self.previous {
            let out_of_order = match self.order {
                SortOrder::Coordinate => {
                    let next = record.reference_id();
                    let by_reference = match (reference, next) {
                        (-1, -1) => std::cmp::Ordering::Equal,
                        (-1, _) => std::cmp::Ordering::Greater,
                        (_, -1) => std::cmp::Ordering::Less,
                        _ => reference.cmp(&next),
                    };
                    by_reference
                        .then(start.cmp(&record.alignment_start()))
                        .is_gt()
                }
                _ => self.previous_name.as_slice() > record.name(),
            };
            if out_of_order {
                return Err(Error::Metric(format!(
                    "Alignments added out of order in SAMFileWriterImpl.addAlignment for {}. \
                     Sort order is {}. Offending records are at [{}] and [{}]",
                    self.path.display(),
                    self.order.name(),
                    String::from_utf8_lossy(&self.previous_name),
                    String::from_utf8_lossy(record.name())
                )));
            }
        }
        self.previous = Some((record.reference_id(), record.alignment_start()));
        self.previous_name.clear();
        self.previous_name.extend_from_slice(record.name());
        Ok(())
    }

    pub fn write(&mut self, record: &RawRecord) -> Result<()> {
        if matches!(self.order, SortOrder::Coordinate | SortOrder::Queryname) {
            self.check_order(record)?;
        }
        record
            .write_to(self.inner.get_mut())
            .map_err(|source| Error::Output {
                path: self.path.clone(),
                source,
            })
    }

    /// Write the last block and the EOF marker. The file keeps its
    /// temporary name until [`finish`] or [`AtomicFile::commit`].
    pub fn close(mut self) -> Result<AtomicFile> {
        let error = |source| Error::Output {
            path: self.path.clone(),
            source,
        };
        let mut file = self.inner.get_mut().finish().map_err(error)?;
        file.flush().map_err(error)?;
        Ok(self.file)
    }
}

/// What [`finish`] writes beside the BAM.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Companions {
    /// A `.bai` in place of the BAM's extension, for a coordinate-sorted
    /// header only.
    pub index: bool,
    /// `<BAM>.md5`: the hex MD5 of the file.
    pub md5: bool,
}

/// Close `writer`, write the index and the MD5 asked for from the finished
/// file, then rename all of them into place. Returns the final paths, the
/// BAM first. `header` is the one the writer was created with.
pub fn finish(
    writer: BamWriter,
    header: &sam::Header,
    companions: Companions,
) -> Result<Vec<PathBuf>> {
    let threads = writer.threads;
    let output = writer.close()?;
    let target = output.path().to_path_buf();
    let coordinate = SortOrder::of(header) == SortOrder::Coordinate;
    let index = if companions.index && coordinate {
        Some(write_index(
            output.writing(),
            &target.with_extension("bai"),
            threads,
        )?)
    } else {
        None
    };
    let md5 = if companions.md5 {
        Some(write_md5(output.writing(), &md5_path(&target))?)
    } else {
        None
    };
    let mut written = vec![target.clone()];
    output.commit().map_err(|source| Error::Output {
        path: target,
        source,
    })?;
    for file in index.into_iter().chain(md5) {
        let path = file.path().to_path_buf();
        file.commit().map_err(|source| Error::Output {
            path: path.clone(),
            source,
        })?;
        written.push(path);
    }
    Ok(written)
}

/// The BAI of the BAM at `bam`, into a file to be renamed to `target`:
/// the bytes htsjdk's `BAMIndexer` writes for that file ([`crate::bai`]).
/// The finished BAM is read once more, its blocks inflated by `threads`
/// workers, so the offsets are the ones of the file as written.
fn write_index(bam: &Path, target: &Path, threads: usize) -> Result<AtomicFile> {
    let output = |source| Error::Output {
        path: target.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(target);
    let out = BufWriter::with_capacity(1 << 16, file.create().map_err(output)?);
    let (mut out, _) = crate::bai::index_bam(bam, threads, out, target)?;
    out.flush().map_err(output)?;
    Ok(file)
}

/// Where `CREATE_MD5_FILE` writes: the BAM's path with `.md5` appended.
pub fn md5_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_owned();
    name.push(".md5");
    PathBuf::from(name)
}

/// The hex MD5 of the BAM at `path`, as htsjdk writes it, into a file to
/// be renamed to `target`.
fn write_md5(path: &Path, target: &Path) -> Result<AtomicFile> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let mut file = File::open(path).map_err(input)?;
    let mut md5 = crate::md5::Md5::new();
    let mut buffer = vec![0u8; 1 << 20];
    loop {
        let n = file.read(&mut buffer).map_err(input)?;
        if n == 0 {
            break;
        }
        md5.update(&buffer[..n]);
    }
    let file = AtomicFile::new(target);
    file.create()
        .and_then(|mut out| out.write_all(crate::md5::hex(&md5.finish()).as_bytes()))
        .map_err(|source| Error::Output {
            path: target.to_path_buf(),
            source,
        })?;
    Ok(file)
}

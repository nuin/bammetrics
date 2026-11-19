//! [`Fanout`]: the files of one split, written from one stream with at
//! most [`MAX_OPEN_FILES`] open at once, renamed into place together.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use noodles::bgzf;
use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, Presorted};
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::mark_duplicates::{SortOrder, same_file};
use crate::registry::Invocation;
use crate::sam_text;
use crate::stage::Outputs;

/// Files a fan-out holds open at once. Picard holds one per library.
pub const MAX_OPEN_FILES: usize = 64;
/// Files one split may write. Each costs a few hundred bytes until the
/// stream ends, when all are renamed into place.
pub const MAX_FILES: usize = 100_000;
/// The buffer of an open SAM text file.
const SAM_BUFFER: usize = 1 << 16;

/// What the files of a split hold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Bam,
    /// SAM text, as htsjdk's `SAMTextWriter` writes it.
    Sam,
}

/// The options of every file a fan-out writes.
#[derive(Clone, Copy, Debug)]
pub struct WriterOptions {
    /// BGZF level, 0 to 9.
    pub compression_level: u8,
    /// A `.bai` beside each coordinate-sorted BAM.
    pub create_index: bool,
    /// `<file>.md5` beside each file.
    pub create_md5: bool,
    /// Workers that read a finished BAM back for its index.
    pub threads: usize,
}

impl WriterOptions {
    /// `COMPRESSION_LEVEL`, `CREATE_INDEX` and `CREATE_MD5_FILE` of `inv`.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let level: i64 = inv.parse("COMPRESSION_LEVEL", 5)?;
        let compression_level = u8::try_from(level)
            .ok()
            .filter(|l| *l <= 9)
            .ok_or_else(|| Error::Usage(format!("COMPRESSION_LEVEL={level}: expected 0 to 9")))?;
        Ok(Self {
            compression_level,
            create_index: inv.flag("CREATE_INDEX", false)?,
            create_md5: inv.flag("CREATE_MD5_FILE", false)?,
            threads: inv.threads,
        })
    }
}

/// The files of a chain's steps, shared by the split sinks of one chain:
/// a split claims each file when it makes it, so no two steps write one
/// file whose name is known only once the stream is read.
pub type Claims = Arc<Mutex<Outputs>>;

/// An open file of a fan-out.
enum Open {
    Bam(bgzf::io::Writer<File>),
    Sam(BufWriter<File>),
}

/// One file of a fan-out.
struct Output {
    /// Dropped before `file`, so a writer never writes to a removed name.
    open: Option<Open>,
    file: AtomicFile,
    /// Device and inode of the file as created, asked again on each open.
    identity: (u64, u64),
    format: Format,
    coordinate: bool,
    order: Presorted,
    /// The end-of-file block is written.
    ended: bool,
    last_used: u64,
}

impl Output {
    fn target(&self) -> &Path {
        self.file.path()
    }

    fn error(&self, source: io::Error) -> Error {
        Error::Output {
            path: self.target().to_path_buf(),
            source,
        }
    }
}

/// Device and inode of an open file.
fn identity(file: &File) -> io::Result<(u64, u64)> {
    use std::os::unix::fs::MetadataExt;
    let meta = file.metadata()?;
    Ok((meta.dev(), meta.ino()))
}

/// Several files written from one stream ([`MAX_OPEN_FILES`] open at
/// once), renamed into place together by [`Fanout::finish`].
pub struct Fanout {
    outputs: Vec<Output>,
    /// Indices of the open outputs.
    open: Vec<usize>,
    clock: u64,
    options: WriterOptions,
    level: bgzf::io::writer::CompressionLevel,
    /// Each file name made so far and what made it, for the error of two
    /// libraries with one file name.
    names: HashMap<PathBuf, String>,
    /// Files the run reads: no output may be one.
    reads: Vec<PathBuf>,
    /// The step label and the chain's claims, in a chain.
    claims: Option<(String, Claims)>,
    references: Vec<Vec<u8>>,
    line: Vec<u8>,
}

impl Fanout {
    pub fn new(
        options: WriterOptions,
        header: &sam::Header,
        reads: Vec<PathBuf>,
        claims: Option<(String, Claims)>,
    ) -> Result<Self> {
        let level = bgzf::io::writer::CompressionLevel::try_from(options.compression_level)
            .map_err(|e| {
                Error::Usage(format!(
                    "COMPRESSION_LEVEL={}: {e}",
                    options.compression_level
                ))
            })?;
        Ok(Self {
            outputs: Vec::new(),
            open: Vec::new(),
            clock: 0,
            options,
            level,
            names: HashMap::new(),
            reads,
            claims,
            references: header
                .reference_sequences()
                .keys()
                .map(|name| name.to_vec())
                .collect(),
            line: Vec::new(),
        })
    }

    /// The files beside `target` the finish writes.
    fn companions(&self, target: &Path, format: Format, coordinate: bool) -> Vec<PathBuf> {
        let mut files = Vec::new();
        if self.options.create_index && format == Format::Bam && coordinate {
            files.push(bam_writer::index_path(target));
        }
        if self.options.create_md5 {
            files.push(bam_writer::md5_path(target));
        }
        files
    }

    /// Refuse `path` when another file of this split, a file the run
    /// reads, or a file of another step of the chain is the same file.
    fn claim(&mut self, path: &Path, what: &str) -> Result<()> {
        if let Some(other) = self.names.get(path) {
            return Err(Error::Metric(format!(
                "{}: both {other} and {what} are written to this file; Picard writes one over \
                 the other",
                path.display()
            )));
        }
        if let Some(read) = self.reads.iter().find(|read| same_file(read, path)) {
            return Err(Error::Usage(format!(
                "{}: it is the input {}; writing it would destroy it",
                path.display(),
                read.display()
            )));
        }
        if let Some((label, claims)) = &self.claims {
            let mut claims = claims.lock().unwrap_or_else(|e| e.into_inner());
            match claims.owner(path) {
                Some(owner) if owner != label => {
                    return Err(Error::Metric(format!(
                        "{owner} and {label} both write {}: every output of a chain needs its \
                         own file",
                        path.display()
                    )));
                }
                Some(_) => {}
                None => claims.claim(label, path)?,
            }
        }
        self.names.insert(path.to_path_buf(), what.to_owned());
        Ok(())
    }

    /// Start `target` under its temporary name with `header`; `what` names
    /// it in errors (the library). Returns its index.
    pub fn create(
        &mut self,
        target: &Path,
        what: &str,
        header: &sam::Header,
        format: Format,
    ) -> Result<usize> {
        if self.outputs.len() >= MAX_FILES {
            return Err(Error::Metric(format!(
                "{}: a split writes at most {MAX_FILES} files",
                target.display()
            )));
        }
        let coordinate = SortOrder::of(header) == SortOrder::Coordinate;
        for path in
            std::iter::once(target.to_path_buf()).chain(self.companions(target, format, coordinate))
        {
            self.claim(&path, what)?;
        }
        let error = |source| Error::Output {
            path: target.to_path_buf(),
            source,
        };
        let file = AtomicFile::new(target);
        let handle = file.create().map_err(error)?;
        let identity = identity(&handle).map_err(error)?;
        let open = match format {
            Format::Bam => {
                let mut writer = bgzf::io::writer::Builder::default()
                    .set_compression_level(self.level)
                    .build_from_writer(handle);
                bam_writer::write_bam_header(&mut writer, header).map_err(error)?;
                Open::Bam(writer)
            }
            Format::Sam => {
                let mut writer = BufWriter::with_capacity(SAM_BUFFER, handle);
                writer
                    .write_all(&sam_text::header_text(header, false))
                    .map_err(error)?;
                Open::Sam(writer)
            }
        };
        let index = self.outputs.len();
        self.outputs.push(Output {
            open: None,
            file,
            identity,
            format,
            coordinate,
            order: Presorted::of(header),
            ended: false,
            last_used: 0,
        });
        self.make_room()?;
        self.outputs[index].open = Some(open);
        self.open.push(index);
        self.touch(index);
        Ok(index)
    }

    /// The files open now: at most [`MAX_OPEN_FILES`].
    pub fn open_files(&self) -> usize {
        self.open.len()
    }

    fn touch(&mut self, index: usize) {
        self.clock += 1;
        self.outputs[index].last_used = self.clock;
    }

    /// Close the least recently written open file when [`MAX_OPEN_FILES`]
    /// are open.
    fn make_room(&mut self) -> Result<()> {
        if self.open.len() < MAX_OPEN_FILES {
            return Ok(());
        }
        let (at, _) = self
            .open
            .iter()
            .enumerate()
            .min_by_key(|(_, index)| self.outputs[**index].last_used)
            .expect("the open set is full");
        let index = self.open.swap_remove(at);
        self.suspend(index)
    }

    /// Close an open file without ending it: its staged records become a
    /// block, the descriptor is released.
    fn suspend(&mut self, index: usize) -> Result<()> {
        let output = &mut self.outputs[index];
        match output.open.take() {
            Some(Open::Bam(mut writer)) => {
                writer.flush().map_err(|e| output.error(e))?;
                drop(writer.into_inner());
            }
            Some(Open::Sam(mut writer)) => writer.flush().map_err(|e| output.error(e))?,
            None => {}
        }
        Ok(())
    }

    /// Open a closed file again, appending, after checking it is the file
    /// this run created.
    fn reopen(&mut self, index: usize) -> Result<()> {
        if self.outputs[index].open.is_some() {
            return Ok(());
        }
        self.make_room()?;
        let level = self.level;
        let output = &mut self.outputs[index];
        let file = File::options()
            .append(true)
            .open(output.file.writing())
            .map_err(|e| output.error(e))?;
        if identity(&file).map_err(|e| output.error(e))? != output.identity {
            return Err(output.error(io::Error::other(
                "the file was replaced while it was written",
            )));
        }
        output.open = Some(match output.format {
            Format::Bam => Open::Bam(
                bgzf::io::writer::Builder::default()
                    .set_compression_level(level)
                    .build_from_writer(file),
            ),
            Format::Sam => Open::Sam(BufWriter::with_capacity(SAM_BUFFER, file)),
        });
        self.open.push(index);
        Ok(())
    }

    /// Write `record` to file `index`.
    pub fn write(&mut self, index: usize, record: &RawRecord) -> Result<()> {
        let target = self.outputs[index].target().to_path_buf();
        self.outputs[index].order.check(record, &target)?;
        self.reopen(index)?;
        self.touch(index);
        let Self {
            outputs,
            references,
            line,
            ..
        } = self;
        let output = &mut outputs[index];
        let written = match output.open.as_mut().expect("reopened") {
            Open::Bam(writer) => record.write_to(writer),
            Open::Sam(writer) => {
                line.clear();
                sam_text::write_record(line, record, references).map_err(Error::Record)?;
                writer.write_all(line)
            }
        };
        written.map_err(|e| output.error(e))
    }

    /// Write the end of file `index` and close it. It is renamed into
    /// place by [`Fanout::finish`].
    pub fn end(&mut self, index: usize) -> Result<()> {
        if self.outputs[index].ended {
            return Ok(());
        }
        if self.outputs[index].format == Format::Bam {
            self.reopen(index)?;
        }
        self.open.retain(|open| *open != index);
        let output = &mut self.outputs[index];
        match output.open.take() {
            Some(Open::Bam(writer)) => drop(writer.finish().map_err(|e| output.error(e))?),
            Some(Open::Sam(mut writer)) => writer.flush().map_err(|e| output.error(e))?,
            None => {}
        }
        output.ended = true;
        Ok(())
    }

    /// End every file, write the index and MD5 asked for, then rename all
    /// of them into place. On any failure no file of the split is left:
    /// those renamed already are removed, the others never leave their
    /// temporary names. Returns the paths, each file then its companions.
    pub fn finish(mut self) -> Result<Vec<PathBuf>> {
        for index in 0..self.outputs.len() {
            self.end(index)?;
        }
        let mut staged: Vec<AtomicFile> = Vec::with_capacity(self.outputs.len());
        for output in std::mem::take(&mut self.outputs) {
            let target = output.target().to_path_buf();
            let temp = output.file.writing().to_path_buf();
            let index =
                (self.options.create_index && output.format == Format::Bam && output.coordinate)
                    .then(|| {
                        bam_writer::write_index(
                            &temp,
                            &target,
                            &bam_writer::index_path(&target),
                            self.options.threads,
                        )
                    })
                    .transpose()?;
            let md5 = self
                .options
                .create_md5
                .then(|| bam_writer::write_md5(&temp, &bam_writer::md5_path(&target)))
                .transpose()?;
            staged.push(output.file);
            staged.extend(index);
            staged.extend(md5);
        }
        commit_all(staged)
    }
}

/// Rename every file into place. When one fails, or when two of them
/// turn out to be one file (names that differ in case only, on a file
/// system that ignores case), every file renamed is removed again.
fn commit_all(staged: Vec<AtomicFile>) -> Result<Vec<PathBuf>> {
    let mut done: Vec<PathBuf> = Vec::with_capacity(staged.len());
    let undo = |done: &[PathBuf]| {
        for path in done {
            let _ = std::fs::remove_file(path);
        }
    };
    let mut files = staged.into_iter();
    while let Some(file) = files.next() {
        let path = file.path().to_path_buf();
        if let Err(source) = file.commit() {
            undo(&done);
            // The files not yet renamed are dropped: their temporary
            // names are removed.
            drop(files);
            return Err(Error::Output { path, source });
        }
        done.push(path);
    }
    let mut seen: HashMap<(u64, u64), &Path> = HashMap::with_capacity(done.len());
    for path in &done {
        let Ok(meta) = std::fs::metadata(path) else {
            continue;
        };
        use std::os::unix::fs::MetadataExt;
        if let Some(other) = seen.insert((meta.dev(), meta.ino()), path) {
            let message = format!(
                "{} and {} are one file on this file system (their names differ in case \
                 only); Picard writes one over the other",
                other.display(),
                path.display()
            );
            undo(&done);
            return Err(Error::Metric(message));
        }
    }
    Ok(done)
}

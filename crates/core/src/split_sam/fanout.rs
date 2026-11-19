//! [`Fanout`]: the files of one split, written from one stream with at
//! most [`MAX_OPEN_FILES`] open at once, renamed into place together.
//!
//! Each BAM file stages its bytes and cuts them into BGZF blocks of
//! [`BLOCK`] bytes; a full block of any file goes to the workers the files
//! share, and each file writes its frames in the order it sent its blocks.
//! A file closed to make room first sends its staged bytes as a short
//! block and writes every frame it has pending.

use std::collections::{HashMap, VecDeque};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, PoisonError, mpsc};
use std::thread;

use noodles::bgzf;
use noodles::bgzf::io::writer::CompressionLevel;
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
    /// Workers that compress the blocks of every file and read a finished
    /// BAM back for its index.
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

/// The uncompressed bytes of one BGZF block: at most what one block holds
/// with level 0's overhead (65,505).
const BLOCK: usize = 65_280;
/// The BGZF end-of-file block (SAM specification 4.1.2).
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// One block compressed to its BGZF frame by noodles' writer: `data` fits
/// one block, so the writer writes one frame and no end-of-file block.
fn frame(data: &[u8], level: CompressionLevel) -> io::Result<Vec<u8>> {
    let mut writer = bgzf::io::writer::Builder::default()
        .set_compression_level(level)
        .build_from_writer(Vec::with_capacity(data.len() / 2 + 64));
    writer.write_all(data)?;
    writer.flush()?;
    Ok(writer.into_inner())
}

/// A block on its way to a worker, and where its frame goes back.
type Job = (
    Vec<u8>,
    CompressionLevel,
    mpsc::SyncSender<io::Result<Vec<u8>>>,
);

/// The compression workers every file of a fan-out shares: each full block
/// of any file is one job, and each file writes its frames in the order it
/// sent its blocks.
struct Pool {
    jobs: Option<mpsc::Sender<Job>>,
    workers: Vec<thread::JoinHandle<()>>,
}

impl Pool {
    fn new(threads: usize) -> Self {
        let (jobs, queue) = mpsc::channel::<Job>();
        let queue = Arc::new(Mutex::new(queue));
        let workers = (0..threads)
            .map(|_| {
                let queue = Arc::clone(&queue);
                thread::spawn(move || {
                    loop {
                        let job = queue.lock().unwrap_or_else(PoisonError::into_inner).recv();
                        let Ok((data, level, done)) = job else {
                            break;
                        };
                        // The file that sent it may have failed and gone.
                        let _ = done.send(frame(&data, level));
                    }
                })
            })
            .collect();
        Self {
            jobs: Some(jobs),
            workers,
        }
    }

    fn send(&self, data: Vec<u8>, level: CompressionLevel) -> io::Result<Pending> {
        let (done, frame) = mpsc::sync_channel(1);
        self.jobs
            .as_ref()
            .and_then(|jobs| jobs.send((data, level, done)).ok())
            .ok_or_else(|| io::Error::other("the compression workers stopped"))?;
        Ok(frame)
    }
}

impl Drop for Pool {
    fn drop(&mut self) {
        drop(self.jobs.take());
        for worker in self.workers.drain(..) {
            let _ = worker.join();
        }
    }
}

/// A frame a worker is compressing.
type Pending = mpsc::Receiver<io::Result<Vec<u8>>>;

/// An open file of a fan-out.
enum Open {
    /// BGZF frames are written as they come back from the workers.
    Bam(File),
    Sam(BufWriter<File>),
}

/// One file of a fan-out.
struct Output {
    /// Dropped before `file`, so nothing is written to a removed name.
    open: Option<Open>,
    file: AtomicFile,
    /// Device and inode of the file as created, asked again on each open.
    identity: (u64, u64),
    format: Format,
    coordinate: bool,
    order: Presorted,
    /// BAM bytes not yet in a block: at most [`BLOCK`].
    staged: Vec<u8>,
    /// Blocks sent to the workers, oldest first; only an open file has any.
    pending: VecDeque<Pending>,
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

    /// Write `bytes` to the open file.
    fn put(&mut self, bytes: &[u8]) -> Result<()> {
        let written = match self.open.as_mut() {
            Some(Open::Bam(file)) => file.write_all(bytes),
            Some(Open::Sam(writer)) => writer.write_all(bytes),
            None => Err(io::Error::other("the file is not open")),
        };
        written.map_err(|e| self.error(e))
    }

    /// Write the frames the workers have finished, oldest first; with
    /// `all`, wait for every one. A file with more than `limit` frames
    /// pending waits for its oldest.
    fn drain(&mut self, all: bool, limit: usize) -> Result<()> {
        while let Some(front) = self.pending.front() {
            let ready = if all || self.pending.len() > limit {
                front.recv().ok()
            } else {
                match front.try_recv() {
                    Ok(frame) => Some(frame),
                    Err(mpsc::TryRecvError::Empty) => return Ok(()),
                    Err(mpsc::TryRecvError::Disconnected) => None,
                }
            };
            let frame = ready
                .unwrap_or_else(|| Err(io::Error::other("a compression worker stopped")))
                .map_err(|e| self.error(e))?;
            self.pending.pop_front();
            self.put(&frame)?;
        }
        Ok(())
    }
}

/// Device and inode of an open file.
fn identity(file: &File) -> io::Result<(u64, u64)> {
    use std::os::unix::fs::MetadataExt;
    let meta = file.metadata()?;
    Ok((meta.dev(), meta.ino()))
}

/// Several files written from one stream ([`MAX_OPEN_FILES`] open at
/// once), renamed into place together by [`Fanout::finish`]. Full blocks
/// of every file are compressed on `threads` shared workers.
pub struct Fanout {
    outputs: Vec<Output>,
    /// Indices of the open outputs.
    open: Vec<usize>,
    clock: u64,
    options: WriterOptions,
    level: CompressionLevel,
    /// `None` with one thread: blocks are compressed where they are cut.
    pool: Option<Pool>,
    /// Frames one file may have pending before it waits for its oldest.
    limit: usize,
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
        let level = CompressionLevel::try_from(options.compression_level).map_err(|e| {
            Error::Usage(format!(
                "COMPRESSION_LEVEL={}: {e}",
                options.compression_level
            ))
        })?;
        let threads = options.threads.max(1);
        Ok(Self {
            outputs: Vec::new(),
            open: Vec::new(),
            clock: 0,
            options,
            level,
            pool: (threads > 1).then(|| Pool::new(threads)),
            limit: threads + 1,
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
            let mut claims = claims.lock().unwrap_or_else(PoisonError::into_inner);
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
        self.make_room()?;
        let error = |source| Error::Output {
            path: target.to_path_buf(),
            source,
        };
        let file = AtomicFile::new(target);
        let handle = file.create().map_err(error)?;
        let identity = identity(&handle).map_err(error)?;
        let mut staged = Vec::new();
        let open = match format {
            Format::Bam => {
                bam_writer::write_bam_header(&mut staged, header).map_err(error)?;
                Open::Bam(handle)
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
            open: Some(open),
            file,
            identity,
            format,
            coordinate,
            order: Presorted::of(header),
            staged,
            pending: VecDeque::new(),
            ended: false,
            last_used: 0,
        });
        self.open.push(index);
        self.touch(index);
        self.send_full_blocks(index)?;
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

    /// Compress `block` of file `index`: on a worker, or here with one
    /// thread. The file is open.
    fn send(&mut self, index: usize, block: Vec<u8>) -> Result<()> {
        let output = &mut self.outputs[index];
        match &self.pool {
            Some(pool) => {
                let pending = pool.send(block, self.level).map_err(|e| output.error(e))?;
                output.pending.push_back(pending);
                output.drain(false, self.limit)
            }
            None => {
                let frame = frame(&block, self.level).map_err(|e| output.error(e))?;
                output.put(&frame)
            }
        }
    }

    /// Send every full block of the staged bytes of file `index`.
    fn send_full_blocks(&mut self, index: usize) -> Result<()> {
        while self.outputs[index].staged.len() >= BLOCK {
            let staged = &mut self.outputs[index].staged;
            let rest = staged.split_off(BLOCK);
            let block = std::mem::replace(staged, rest);
            self.send(index, block)?;
        }
        Ok(())
    }

    /// Write out everything file `index` holds: the staged bytes as a
    /// block, then every pending frame.
    fn settle(&mut self, index: usize) -> Result<()> {
        let staged = std::mem::take(&mut self.outputs[index].staged);
        if !staged.is_empty() {
            self.send(index, staged)?;
        }
        let limit = self.limit;
        self.outputs[index].drain(true, limit)
    }

    /// Close an open file without ending it: its staged bytes become a
    /// block, every frame is written, the descriptor is released.
    fn suspend(&mut self, index: usize) -> Result<()> {
        self.settle(index)?;
        let output = &mut self.outputs[index];
        if let Some(Open::Sam(mut writer)) = output.open.take() {
            writer.flush().map_err(|e| output.error(e))?;
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
            Format::Bam => Open::Bam(file),
            Format::Sam => Open::Sam(BufWriter::with_capacity(SAM_BUFFER, file)),
        });
        self.open.push(index);
        Ok(())
    }

    /// Write `record` to file `index`.
    pub fn write(&mut self, index: usize, record: &RawRecord) -> Result<()> {
        let output = &mut self.outputs[index];
        output.order.check(record, output.file.path())?;
        self.reopen(index)?;
        self.touch(index);
        let Self {
            outputs,
            references,
            line,
            ..
        } = self;
        let output = &mut outputs[index];
        match output.format {
            Format::Bam => {
                record
                    .write_to(&mut output.staged)
                    .map_err(|e| output.error(e))?;
                if output.staged.len() >= BLOCK {
                    self.send_full_blocks(index)?;
                }
                Ok(())
            }
            Format::Sam => {
                line.clear();
                sam_text::write_record(line, record, references).map_err(Error::Record)?;
                output.put(line)
            }
        }
    }

    /// Write the end of file `index` and close it. It is renamed into
    /// place by [`Fanout::finish`].
    pub fn end(&mut self, index: usize) -> Result<()> {
        if self.outputs[index].ended {
            return Ok(());
        }
        self.reopen(index)?;
        self.settle(index)?;
        self.open.retain(|open| *open != index);
        let output = &mut self.outputs[index];
        if output.format == Format::Bam {
            output.put(&BGZF_EOF)?;
        }
        match output.open.take() {
            Some(Open::Sam(mut writer)) => writer.flush().map_err(|e| output.error(e))?,
            Some(Open::Bam(file)) => drop(file),
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

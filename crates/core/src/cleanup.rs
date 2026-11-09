//! Files and directories this process is writing and must not leave
//! behind: spilled sort runs ([`crate::spill::TempDir`]) and the temporary
//! file of an output written in place by rename ([`AtomicFile`]).
//!
//! Each path is registered by the code that computed it and unregistered
//! when that code removes or renames it. [`remove_registered`] removes what
//! is still registered; the `bammetrics` binary calls it from its SIGINT,
//! SIGTERM and SIGHUP handler, so an interrupted run leaves no spill file
//! and no partial output. Nothing else is ever removed.

use std::fs::{self, File};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};

static REGISTERED: Mutex<Vec<PathBuf>> = Mutex::new(Vec::new());

fn registered() -> MutexGuard<'static, Vec<PathBuf>> {
    REGISTERED
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Remember `path` for [`remove_registered`].
pub fn register(path: &Path) {
    registered().push(path.to_path_buf());
}

/// Forget `path` (it was removed or renamed into place).
pub fn unregister(path: &Path) {
    let mut paths = registered();
    if let Some(i) = paths.iter().rposition(|p| p == path) {
        paths.remove(i);
    }
}

/// Remove every registered path, newest first, so the files of a
/// directory go before the directory. The returned guard keeps the
/// registry locked: a signal handler holds it until the process exits,
/// so no thread registers a new file in between.
pub fn remove_registered() -> MutexGuard<'static, Vec<PathBuf>> {
    let mut paths = registered();
    remove_newest_first(paths.drain(..));
    paths
}

/// Remove each file, or each empty directory, last path first.
fn remove_newest_first(paths: impl DoubleEndedIterator<Item = PathBuf>) {
    for path in paths.rev() {
        if fs::remove_file(&path).is_err() {
            let _ = fs::remove_dir(&path);
        }
    }
}

/// An output file written under a temporary name in its own directory and
/// renamed over `path` by [`AtomicFile::commit`]. Dropped without a commit,
/// the temporary file is removed and `path` is left as it was. A path that
/// exists and is not a regular file (`/dev/stdout`, a named pipe) or is a
/// symbolic link is written in place.
pub struct AtomicFile {
    path: PathBuf,
    /// The temporary name, `None` when writing in place.
    temp: Option<PathBuf>,
}

impl AtomicFile {
    pub fn new(path: &Path) -> Self {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let in_place = fs::symlink_metadata(path).is_ok_and(|m| !m.file_type().is_file())
            || path.starts_with("/dev");
        let temp = (!in_place).then(|| {
            let name = path
                .file_name()
                .map_or_else(|| "output".into(), |n| n.to_string_lossy().into_owned());
            let temp = path.with_file_name(format!(
                ".{name}.{}.{}.tmp",
                std::process::id(),
                COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            register(&temp);
            temp
        });
        Self {
            path: path.to_path_buf(),
            temp,
        }
    }

    /// The path written to until the commit.
    pub fn writing(&self) -> &Path {
        self.temp.as_deref().unwrap_or(&self.path)
    }

    /// The final path.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Create the file to write: the temporary one, or `path` in place.
    pub fn create(&self) -> io::Result<File> {
        match &self.temp {
            Some(temp) => File::options().write(true).create_new(true).open(temp),
            None => File::create(&self.path),
        }
    }

    /// Rename the temporary file over `path`.
    pub fn commit(mut self) -> io::Result<()> {
        let Some(temp) = self.temp.take() else {
            return Ok(());
        };
        let renamed = fs::rename(&temp, &self.path);
        if renamed.is_err() {
            let _ = fs::remove_file(&temp);
        }
        unregister(&temp);
        renamed
    }
}

impl Drop for AtomicFile {
    fn drop(&mut self) {
        if let Some(temp) = self.temp.take() {
            let _ = fs::remove_file(&temp);
            unregister(&temp);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Write `path` through an [`AtomicFile`]; `write` fills it.
    fn write_atomically(
        path: &Path,
        write: impl FnOnce(&mut File) -> io::Result<()>,
    ) -> io::Result<()> {
        let file = AtomicFile::new(path);
        write(&mut file.create()?)?;
        file.commit()
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "bammetrics-cleanup-test.{}.{name}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).expect("scratch dir");
        dir
    }

    #[test]
    fn an_uncommitted_file_leaves_the_old_one() {
        let dir = scratch("uncommitted");
        let path = dir.join("out.txt");
        fs::write(&path, "old").expect("old");
        let failed = write_atomically(&path, |out| {
            out.write_all(b"new")?;
            Err(io::Error::other("stop"))
        });
        assert!(failed.is_err());
        assert_eq!(fs::read_to_string(&path).expect("read"), "old");
        assert_eq!(fs::read_dir(&dir).expect("dir").count(), 1);
        write_atomically(&path, |out| out.write_all(b"new")).expect("write");
        assert_eq!(fs::read_to_string(&path).expect("read"), "new");
        assert_eq!(fs::read_dir(&dir).expect("dir").count(), 1);
        fs::remove_file(&path).expect("remove");
        fs::remove_dir(&dir).expect("remove dir");
    }

    #[test]
    fn registered_paths_are_removed_files_first() {
        let dir = scratch("registered");
        let inner = dir.join("inner");
        fs::create_dir(&inner).expect("inner");
        let file = inner.join("run.tmp");
        fs::write(&file, "x").expect("file");
        let kept = dir.join("kept.txt");
        fs::write(&kept, "x").expect("kept");
        remove_newest_first(vec![inner.clone(), file].into_iter());
        assert!(!inner.exists());
        assert!(kept.exists());
        fs::remove_file(&kept).expect("remove");
        fs::remove_dir(&dir).expect("remove dir");
    }
}

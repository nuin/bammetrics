//! One run per manifest: a chain that writes `PREFIX.manifest.json` holds
//! an advisory lock (`flock`) on `PREFIX.manifest.json.lock` from before
//! its first read to its end. A second run on the same prefix is refused
//! with exit 2 before it reads anything, so two runs never write one kept
//! BAM or one manifest (docs/chaining.md section 27).
//!
//! The lock file is removed when the run ends, and on SIGINT, SIGTERM and
//! SIGHUP through the cleanup registry. A file a killed run left is no
//! lock: `flock` ends with its process, and the next run takes the file.

use std::fs::{File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};

use bammetrics_core::cleanup;

use super::{VERB, usage};
use crate::error::VerbError;
use crate::qc::with_suffix;

/// The lock of one prefix, held until dropped.
#[derive(Debug)]
pub struct PrefixLock {
    path: PathBuf,
    /// Open while held: closing it releases the lock.
    file: Option<File>,
}

/// The device and inode of `path` itself, a link not followed.
#[cfg(unix)]
fn id(meta: &std::fs::Metadata) -> (u64, u64) {
    use std::os::unix::fs::MetadataExt;
    (meta.dev(), meta.ino())
}

#[cfg(not(unix))]
fn id(meta: &std::fs::Metadata) -> (u64, u64) {
    (meta.len(), 0)
}

impl PrefixLock {
    /// The lock file of `manifest`.
    pub fn path_of(manifest: &Path) -> PathBuf {
        with_suffix(manifest, ".lock")
    }

    /// Take the lock of `manifest`, or refuse: another run holds it.
    pub fn take(manifest: &Path) -> Result<Self, VerbError> {
        let path = Self::path_of(manifest);
        let failed = |e: io::Error| VerbError::input(VERB, format!("{}: {e}", path.display()));
        // A run that took the file just before another removed it holds a
        // lock on a name that is gone: take the new file then.
        for _ in 0..8 {
            if std::fs::symlink_metadata(&path).is_ok_and(|m| !m.file_type().is_file()) {
                return Err(usage(format!(
                    "{} is not a regular file: it cannot be the lock of {}",
                    path.display(),
                    manifest.display()
                )));
            }
            let file = OpenOptions::new()
                .read(true)
                .write(true)
                .create(true)
                .truncate(false)
                .open(&path)
                .map_err(failed)?;
            match file.try_lock() {
                Ok(()) => {}
                Err(std::fs::TryLockError::WouldBlock) => {
                    return Err(usage(format!(
                        "another chain is writing {} (it holds {}): one run per output prefix",
                        manifest.display(),
                        path.display()
                    ))
                    .hint("wait for it to end, or give this run another -o PREFIX"));
                }
                Err(std::fs::TryLockError::Error(e)) => return Err(failed(e)),
            }
            let held = file.metadata().map_err(failed)?;
            let named = std::fs::symlink_metadata(&path);
            if named.is_ok_and(|m| id(&m) == id(&held)) {
                // Registered once held: a signal before then leaves a file
                // that is no lock, never the lock of another run.
                cleanup::register(&path);
                return Ok(Self {
                    path,
                    file: Some(file),
                });
            }
        }
        Err(failed(io::Error::other(
            "the lock file was replaced 8 times while this run took it",
        )))
    }
}

impl Drop for PrefixLock {
    /// Remove the file while the lock is held, then release it: a run
    /// waiting on the old file finds its name gone and takes a new one.
    fn drop(&mut self) {
        if let Some(file) = self.file.take() {
            let ours = file.metadata().ok().map(|m| id(&m));
            let named = std::fs::symlink_metadata(&self.path).ok().map(|m| id(&m));
            if ours.is_some() && ours == named {
                let _ = std::fs::remove_file(&self.path);
            }
            cleanup::unregister(&self.path);
            drop(file);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_second_run_on_one_prefix_is_refused_and_the_file_goes_at_the_end() {
        let dir = std::env::temp_dir().join(format!("bammetrics-lock-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let manifest = dir.join("s.manifest.json");
        let first = PrefixLock::take(&manifest).expect("first run takes it");
        let lock = PrefixLock::path_of(&manifest);
        assert!(lock.is_file());
        let second = PrefixLock::take(&manifest).expect_err("second run");
        assert_eq!(second.code, 2);
        assert!(second.what.contains("another chain is writing"), "{second}");
        // The refused run leaves the holder's file alone.
        assert!(lock.is_file());
        drop(first);
        assert!(!lock.exists());
        let third = PrefixLock::take(&manifest).expect("taken again");
        drop(third);
        // A file a killed run left is no lock.
        std::fs::write(&lock, "").expect("stale file");
        drop(PrefixLock::take(&manifest).expect("a stale file is taken"));
        assert!(!lock.exists());
        let _ = std::fs::remove_dir_all(&dir);
    }
}

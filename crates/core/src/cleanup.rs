//! Files and directories this process is writing and must not leave
//! behind: spilled sort runs ([`crate::spill::TempDir`]) and the temporary
//! file of an output written in place by rename ([`AtomicFile`]).
//!
//! Each path is registered by the code that made it and unregistered
//! when that code removes or renames it. [`remove_registered`] removes what
//! is still registered; the `bammetrics` binary calls it from its SIGINT,
//! SIGTERM and SIGHUP handler, so an interrupted run leaves no spill file
//! and no partial output. Nothing else is ever removed.
//!
//! Three checks an output gets before the first read: [`refusal`] (the
//! owner check of [`AtomicFile::create`], asked ahead of the write),
//! [`unresolved_parent`] (a `..` behind a directory that does not exist)
//! and [`MadeDirs`] (the directories a run made above its outputs, removed
//! again when the run fails and they are still empty).

use std::fs::{self, File};
use std::io;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};
use std::time::{SystemTime, UNIX_EPOCH};

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
/// the temporary file is removed and `path` is left as it was.
///
/// Three kinds of path are written in place: a path under `/dev` with no
/// `..` component (`/dev/stdout`, `/dev/fd/3`; only root makes files
/// there), an existing file this user owns that is neither regular nor a
/// symbolic link (a named pipe, a socket, a device), and a symbolic link
/// this user owns, which the write follows. Outside `/dev`, a symbolic
/// link or a non-regular file another user owns is refused by
/// [`AtomicFile::create`]: that user chose where the bytes would go.
///
/// The temporary name carries a random part and is opened with
/// `create_new`. A path already there was made by someone else, so the
/// open fails and neither the drop, the commit nor a signal removes it.
pub struct AtomicFile {
    path: PathBuf,
    /// The temporary name, `None` when writing in place.
    temp: Option<PathBuf>,
    /// True once [`AtomicFile::create`] has made the temporary file.
    created: AtomicBool,
    /// Why `path` is not written, told by [`AtomicFile::create`].
    refused: Option<String>,
}

/// How a path is written.
enum Target {
    Rename,
    InPlace,
    Refused(String),
}

/// A number another process cannot predict: the standard library's hash
/// key, which the system seeds at random, over the time, the process id
/// and `counter`.
fn nonce(counter: u64) -> u64 {
    use std::hash::{BuildHasher, Hasher};
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let mut hasher = std::collections::hash_map::RandomState::new().build_hasher();
    hasher.write_u128(nanos);
    hasher.write_u32(std::process::id());
    hasher.write_u64(counter);
    hasher.finish()
}

/// The effective user id of this process: the owner of a file it creates
/// (in the temporary directory, else in `beside`) and removes again. The
/// standard library has no `geteuid`.
#[cfg(unix)]
fn effective_uid(beside: &Path) -> io::Result<u32> {
    use std::os::unix::fs::MetadataExt;
    static UID: std::sync::OnceLock<u32> = std::sync::OnceLock::new();
    if let Some(uid) = UID.get() {
        return Ok(*uid);
    }
    let probe = |dir: &Path| {
        let path = dir.join(format!(
            ".bammetrics-uid.{}.{:016x}.tmp",
            std::process::id(),
            nonce(0)
        ));
        let file = File::options().write(true).create_new(true).open(&path)?;
        let uid = file.metadata().map(|m| m.uid());
        let _ = fs::remove_file(&path);
        uid
    };
    let uid = probe(&std::env::temp_dir()).or_else(|_| probe(beside))?;
    Ok(*UID.get_or_init(|| uid))
}

#[cfg(unix)]
fn target(path: &Path) -> Target {
    use std::os::unix::fs::MetadataExt;
    let leaves_dev = path
        .components()
        .any(|c| matches!(c, std::path::Component::ParentDir));
    if path.starts_with("/dev") && !leaves_dev {
        return Target::InPlace;
    }
    let Ok(metadata) = fs::symlink_metadata(path) else {
        return Target::Rename;
    };
    let kind = metadata.file_type();
    if kind.is_file() {
        return Target::Rename;
    }
    let what = if kind.is_symlink() {
        "a symbolic link"
    } else {
        "not a regular file"
    };
    let parent = match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => parent,
        _ => Path::new("."),
    };
    match effective_uid(parent) {
        Ok(uid) if uid == metadata.uid() => Target::InPlace,
        Ok(_) => Target::Refused(format!(
            "{} is {what} owned by another user; not written in place",
            path.display()
        )),
        Err(e) => Target::Refused(format!(
            "{} is {what} and its owner cannot be compared with this user: {e}",
            path.display()
        )),
    }
}

#[cfg(not(unix))]
fn target(path: &Path) -> Target {
    let in_place = fs::symlink_metadata(path).is_ok_and(|m| !m.file_type().is_file())
        || path.starts_with("/dev");
    if in_place {
        Target::InPlace
    } else {
        Target::Rename
    }
}

/// The error [`AtomicFile::create`] gives for `path` when it refuses it (a
/// symbolic link or a non-regular file another user owns, outside `/dev`),
/// so a caller can fail before it reads its input. Nothing is opened.
pub fn refusal(path: &Path) -> Option<io::Error> {
    match target(path) {
        Target::Refused(reason) => Some(io::Error::new(io::ErrorKind::PermissionDenied, reason)),
        Target::Rename | Target::InPlace => None,
    }
}

/// The first directory of `path` that does not exist, when a `..` comes
/// after it. Such a path names no file yet: once the directory is made,
/// the `..` leads back out of it, to a file no check has seen (an input,
/// or the output of another step). One `stat` per component.
pub fn unresolved_parent(path: &Path) -> Option<PathBuf> {
    let mut prefix = PathBuf::new();
    let mut missing: Option<PathBuf> = None;
    for component in path.components() {
        if missing.is_some() {
            if matches!(component, std::path::Component::ParentDir) {
                return missing;
            }
            continue;
        }
        prefix.push(component);
        if fs::metadata(&prefix).is_err() {
            missing = Some(prefix.clone());
        }
    }
    None
}

/// The directories a run made above its outputs. A run that fails calls
/// [`MadeDirs::remove_empty`], which removes those still empty, deepest
/// first; a directory that was there before, or that holds a file, stays.
/// Each made directory is registered, so a signal removes it under the
/// same rule ([`remove_registered`] removes a directory only when empty).
#[derive(Debug, Default)]
pub struct MadeDirs {
    made: Vec<PathBuf>,
}

impl MadeDirs {
    /// Make the directories above `path` that do not exist, top first.
    pub fn create_parent(&mut self, path: &Path) -> io::Result<()> {
        let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) else {
            return Ok(());
        };
        let missing: Vec<&Path> = parent
            .ancestors()
            .take_while(|dir| !dir.as_os_str().is_empty() && fs::metadata(dir).is_err())
            .collect();
        for dir in missing.into_iter().rev() {
            match fs::create_dir(dir) {
                Ok(()) => {
                    register(dir);
                    self.made.push(dir.to_path_buf());
                }
                // Another process made it first, or `a/..` once `a` is there.
                Err(e) if e.kind() == io::ErrorKind::AlreadyExists && dir.is_dir() => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }

    /// The directories made so far, in the order they were made.
    pub fn made(&self) -> &[PathBuf] {
        &self.made
    }

    /// Remove every made directory that is empty, deepest first, and
    /// forget them all.
    pub fn remove_empty(&mut self) {
        for dir in self.made.drain(..).rev() {
            let _ = fs::remove_dir(&dir);
            unregister(&dir);
        }
    }

    /// Keep every made directory: the run finished.
    pub fn keep(&mut self) {
        for dir in self.made.drain(..) {
            unregister(&dir);
        }
    }
}

impl Drop for MadeDirs {
    /// Dropped without [`MadeDirs::keep`], the run did not finish.
    fn drop(&mut self) {
        self.remove_empty();
    }
}

impl AtomicFile {
    pub fn new(path: &Path) -> Self {
        static COUNTER: AtomicU64 = AtomicU64::new(0);
        let target = target(path);
        let temp = matches!(target, Target::Rename).then(|| {
            let name = path
                .file_name()
                .map_or_else(|| "output".into(), |n| n.to_string_lossy().into_owned());
            let counter = COUNTER.fetch_add(1, Ordering::Relaxed);
            path.with_file_name(format!(
                ".{name}.{}.{counter}.{:016x}.tmp",
                std::process::id(),
                nonce(counter)
            ))
        });
        Self {
            path: path.to_path_buf(),
            temp,
            created: AtomicBool::new(false),
            refused: match target {
                Target::Refused(reason) => Some(reason),
                Target::Rename | Target::InPlace => None,
            },
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
    /// When `path` already exists, the temporary file starts out readable by
    /// this user only and then takes the permissions `path` has, so the
    /// rename keeps the mode the user gave the file. A new `path` gets the
    /// umask default. The temporary file is registered for removal on a
    /// signal once this process has made it.
    pub fn create(&self) -> io::Result<File> {
        if let Some(reason) = &self.refused {
            return Err(io::Error::new(
                io::ErrorKind::PermissionDenied,
                reason.clone(),
            ));
        }
        let Some(temp) = &self.temp else {
            return File::create(&self.path);
        };
        let mut options = File::options();
        options.write(true).create_new(true);
        let existing = fs::metadata(&self.path).ok().map(|m| m.permissions());
        #[cfg(unix)]
        if existing.is_some() {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        let file = options.open(temp)?;
        self.created.store(true, Ordering::Release);
        register(temp);
        if let Some(permissions) = existing {
            file.set_permissions(permissions)?;
        }
        Ok(file)
    }

    /// The temporary path when this process made the file there, taken out
    /// so the drop does nothing more.
    fn take_created(&mut self) -> Option<PathBuf> {
        let temp = self.temp.take()?;
        self.created.load(Ordering::Acquire).then_some(temp)
    }

    /// Rename the temporary file over `path`. Without a file from
    /// [`AtomicFile::create`] nothing is renamed: a file at the temporary
    /// path is not this process's.
    pub fn commit(mut self) -> io::Result<()> {
        if self.temp.is_none() {
            return Ok(());
        }
        let Some(temp) = self.take_created() else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "no temporary file was created",
            ));
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
        if let Some(temp) = self.take_created() {
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

    #[cfg(unix)]
    #[test]
    fn a_replaced_file_keeps_its_mode() {
        use std::os::unix::fs::PermissionsExt;
        let dir = scratch("mode");
        let path = dir.join("out.txt");
        fs::write(&path, "old").expect("old");
        fs::set_permissions(&path, fs::Permissions::from_mode(0o640)).expect("chmod");
        write_atomically(&path, |out| out.write_all(b"new")).expect("write");
        let mode = fs::metadata(&path).expect("metadata").permissions().mode() & 0o777;
        assert_eq!(mode, 0o640);
        assert_eq!(fs::read_to_string(&path).expect("read"), "new");
        fs::remove_file(&path).expect("remove");
        fs::remove_dir(&dir).expect("remove dir");
    }

    /// A link this user made is written through: the target takes the
    /// text and the link stays a link.
    #[cfg(unix)]
    #[test]
    fn a_link_owned_by_this_user_is_written_through() {
        let dir = scratch("own-link");
        let real = dir.join("real.txt");
        fs::write(&real, "old").expect("old");
        let link = dir.join("link.txt");
        std::os::unix::fs::symlink(&real, &link).expect("symlink");
        let file = AtomicFile::new(&link);
        assert_eq!(file.writing(), link);
        drop(file);
        write_atomically(&link, |out| out.write_all(b"new")).expect("write");
        assert_eq!(fs::read_to_string(&real).expect("read"), "new");
        let kind = fs::symlink_metadata(&link).expect("link").file_type();
        assert!(kind.is_symlink());
        assert_eq!(fs::read_dir(&dir).expect("dir").count(), 2);
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// A link another user owns is refused with the path in the error.
    /// The test needs such a link on the machine (`/tmp` and `/etc` on
    /// macOS, `/bin` on a merged-usr Linux) and passes without one.
    #[cfg(unix)]
    #[test]
    fn a_link_owned_by_another_user_is_refused() {
        use std::os::unix::fs::MetadataExt;
        let me = effective_uid(&std::env::temp_dir()).expect("uid");
        let foreign = ["/tmp", "/etc", "/var", "/bin", "/sbin", "/lib"]
            .into_iter()
            .map(Path::new)
            .find(|p| {
                fs::symlink_metadata(p).is_ok_and(|m| m.file_type().is_symlink() && m.uid() != me)
            });
        let Some(link) = foreign else {
            return;
        };
        let file = AtomicFile::new(link);
        let error = file.create().expect_err("refused");
        assert_eq!(error.kind(), io::ErrorKind::PermissionDenied);
        let text = error.to_string();
        assert!(text.contains(&link.display().to_string()), "{text}");
        assert!(text.contains("owned by another user"), "{text}");
    }

    /// `/dev/stdout` and its kin are in place with no owner check; a
    /// path that leaves `/dev` through `..` gets the rules of any path.
    #[cfg(unix)]
    #[test]
    fn the_dev_shortcut_needs_a_path_without_parent_components() {
        for path in ["/dev/stdout", "/dev/stderr", "/dev/fd/1", "/dev/null"] {
            assert!(matches!(target(Path::new(path)), Target::InPlace), "{path}");
        }
        let dir = scratch("dev-parent");
        let outside = PathBuf::from(format!("/dev/..{}/x", dir.display()));
        assert!(outside.starts_with("/dev"));
        assert!(matches!(target(&outside), Target::Rename));
        let file = AtomicFile::new(&outside);
        assert_ne!(file.writing(), outside);
        drop(file);
        fs::write(dir.join("x"), "old").expect("old");
        assert!(matches!(target(&outside), Target::Rename));
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// A named pipe this user made is written in place. Only the
    /// classification is checked: opening a pipe with no reader blocks.
    #[cfg(unix)]
    #[test]
    fn a_pipe_owned_by_this_user_is_in_place() {
        let dir = scratch("fifo");
        let pipe = dir.join("pipe");
        let made = std::process::Command::new("mkfifo")
            .arg(&pipe)
            .status()
            .expect("mkfifo runs");
        assert!(made.success());
        assert!(matches!(target(&pipe), Target::InPlace));
        assert_eq!(AtomicFile::new(&pipe).writing(), pipe);
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// A non-regular file another user owns is refused outside `/dev`.
    /// The test needs one on the machine (a root-owned socket under
    /// `/var/run` or `/run`) and passes without one. Nothing is opened.
    #[cfg(unix)]
    #[test]
    fn a_non_regular_file_owned_by_another_user_is_refused() {
        use std::os::unix::fs::{FileTypeExt, MetadataExt};
        let me = effective_uid(&std::env::temp_dir()).expect("uid");
        let foreign = ["/var/run", "/run"]
            .into_iter()
            .filter_map(|dir| fs::read_dir(dir).ok())
            .flatten()
            .filter_map(|entry| entry.ok().map(|e| e.path()))
            .find(|p| {
                fs::symlink_metadata(p).is_ok_and(|m| {
                    let kind = m.file_type();
                    (kind.is_socket() || kind.is_fifo()) && m.uid() != me
                })
            });
        let Some(path) = foreign else {
            return;
        };
        let error = AtomicFile::new(&path).create().expect_err("refused");
        assert_eq!(error.kind(), io::ErrorKind::PermissionDenied);
        let text = error.to_string();
        assert!(text.contains(&path.display().to_string()), "{text}");
        assert!(text.contains("owned by another user"), "{text}");
    }

    #[cfg(unix)]
    #[test]
    fn the_effective_uid_is_the_owner_of_a_new_file() {
        use std::os::unix::fs::MetadataExt;
        let dir = scratch("uid");
        let path = dir.join("mine");
        fs::write(&path, "x").expect("file");
        let owner = fs::metadata(&path).expect("metadata").uid();
        assert_eq!(effective_uid(&dir).expect("uid"), owner);
        assert_eq!(fs::read_dir(&dir).expect("dir").count(), 1);
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// A file already at the temporary path is not this process's: the
    /// create fails, and the drop and the commit leave it where it is.
    #[test]
    fn a_file_planted_at_the_temporary_path_is_left_in_place() {
        let dir = scratch("planted");
        let path = dir.join("out.txt");
        for commit in [false, true] {
            let file = AtomicFile::new(&path);
            let temp = file.writing().to_path_buf();
            assert_ne!(temp, path);
            fs::write(&temp, "planted").expect("plant");
            let error = file.create().expect_err("create_new refuses it");
            assert_eq!(error.kind(), io::ErrorKind::AlreadyExists);
            assert!(!registered().contains(&temp));
            if commit {
                assert!(file.commit().is_err());
            } else {
                drop(file);
            }
            assert_eq!(fs::read_to_string(&temp).expect("planted"), "planted");
            assert!(!path.exists());
            fs::remove_file(&temp).expect("remove");
        }
        fs::remove_dir(&dir).expect("remove dir");
    }

    /// The temporary file is registered while it exists, and two names
    /// for one path differ in their random part.
    #[test]
    fn the_temporary_name_is_registered_once_created_and_has_a_random_part() {
        let dir = scratch("names");
        let path = dir.join("out.txt");
        let first = AtomicFile::new(&path);
        let second = AtomicFile::new(&path);
        let random = |file: &AtomicFile| {
            let name = file.writing().to_string_lossy().into_owned();
            let parts: Vec<&str> = name.rsplit('.').collect();
            assert_eq!(parts[0], "tmp");
            assert_eq!(parts[1].len(), 16);
            parts[1].to_owned()
        };
        assert_ne!(random(&first), random(&second));
        let temp = first.writing().to_path_buf();
        assert!(!registered().contains(&temp));
        drop(first.create().expect("create"));
        assert!(registered().contains(&temp));
        drop(first);
        assert!(!registered().contains(&temp));
        assert!(!temp.exists());
        drop(second);
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

    /// The refusal asked ahead of the write is the one the write gives.
    #[cfg(unix)]
    #[test]
    fn the_refusal_is_known_before_the_write() {
        use std::os::unix::fs::MetadataExt;
        let dir = scratch("refusal");
        let plain = dir.join("out.txt");
        assert!(refusal(&plain).is_none());
        fs::write(&plain, "old").expect("old");
        assert!(refusal(&plain).is_none());
        assert!(refusal(Path::new("/dev/stdout")).is_none());
        let me = effective_uid(&dir).expect("uid");
        let foreign = ["/tmp", "/etc", "/var", "/bin", "/sbin", "/lib"]
            .into_iter()
            .map(Path::new)
            .find(|p| {
                fs::symlink_metadata(p).is_ok_and(|m| m.file_type().is_symlink() && m.uid() != me)
            });
        if let Some(link) = foreign {
            let asked = refusal(link).expect("refused").to_string();
            let written = AtomicFile::new(link).create().expect_err("refused");
            assert_eq!(asked, written.to_string());
            assert!(asked.ends_with("owned by another user; not written in place"));
        }
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// A `..` behind a directory that does not exist is found; one behind
    /// a directory that exists, or a missing directory without one, is not.
    #[test]
    fn a_parent_component_behind_a_missing_directory_is_found() {
        let dir = scratch("unresolved");
        fs::create_dir(dir.join("there")).expect("dir");
        let found = |tail: &str| unresolved_parent(&dir.join(tail));
        assert_eq!(found("missing/../out.txt"), Some(dir.join("missing")));
        assert_eq!(
            found("there/missing/x/../../out.txt"),
            Some(dir.join("there/missing"))
        );
        assert_eq!(found("there/../out.txt"), None);
        assert_eq!(found("there/../missing/out.txt"), None);
        assert_eq!(found("missing/deeper/out.txt"), None);
        assert_eq!(found("out.txt"), None);
        assert_eq!(unresolved_parent(Path::new("../out.txt")), None);
        assert_eq!(unresolved_parent(Path::new("/dev/stdout")), None);
        // A link that leads nowhere is a directory that does not exist.
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(dir.join("nowhere"), dir.join("link")).expect("symlink");
            assert_eq!(found("link/../out.txt"), Some(dir.join("link")));
        }
        fs::remove_dir_all(&dir).expect("remove dir");
    }

    /// Only the directories the run made are removed, and only when empty.
    #[test]
    fn made_directories_are_removed_when_empty_and_never_others() {
        let dir = scratch("made");
        fs::create_dir(dir.join("old")).expect("old");
        let mut dirs = MadeDirs::default();
        dirs.create_parent(&dir.join("old/out.txt")).expect("old");
        dirs.create_parent(&dir.join("old/a/b/out.txt"))
            .expect("a/b");
        dirs.create_parent(&dir.join("c/out.txt")).expect("c");
        dirs.create_parent(&dir.join("d/e/out.txt")).expect("d/e");
        dirs.create_parent(Path::new("out.txt")).expect("no parent");
        let made = ["old/a", "old/a/b", "c", "d", "d/e"].map(|d| dir.join(d));
        assert_eq!(dirs.made(), made);
        assert!(made.iter().all(|d| registered().contains(d)));
        fs::write(dir.join("d/e/out.txt"), "x").expect("file");
        dirs.remove_empty();
        assert!(dir.join("old").is_dir());
        assert!(!dir.join("old/a").exists());
        assert!(!dir.join("c").exists());
        assert!(dir.join("d/e/out.txt").is_file());
        assert!(made.iter().all(|d| !registered().contains(d)));
        // Kept directories stay, empty or not; a drop removes the empty ones.
        let mut kept = MadeDirs::default();
        kept.create_parent(&dir.join("kept/out.txt")).expect("kept");
        kept.keep();
        drop(kept);
        assert!(dir.join("kept").is_dir());
        let mut dropped = MadeDirs::default();
        dropped
            .create_parent(&dir.join("dropped/out.txt"))
            .expect("dropped");
        drop(dropped);
        assert!(!dir.join("dropped").exists());
        assert!(!registered().iter().any(|p| p.starts_with(&dir)));
        fs::remove_dir_all(&dir).expect("remove dir");
    }
}

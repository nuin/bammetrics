//! `IntervalList.write`: the header as `SAMTextHeaderCodec.encode` writes
//! it, then one line per interval (`IntervalListWriter`).

use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use super::IntervalList;
use crate::cleanup::{AtomicFile, MadeDirs};
use crate::error::{Error, Result};

impl IntervalList {
    /// The file `IntervalList.write` writes: the header lines, then
    /// `contig`, `start`, `end`, `+` or `-` and the name (`.` without one),
    /// tab-separated, each line ended by `\n`.
    pub fn write_to(&self, out: &mut impl Write) -> io::Result<()> {
        out.write_all(self.header.encode().as_bytes())?;
        for interval in &self.intervals {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                self.contig_name(interval),
                interval.start,
                interval.end,
                if interval.negative_strand { '-' } else { '+' },
                interval.name.as_deref().unwrap_or(".")
            )?;
        }
        Ok(())
    }

    /// [`IntervalList::write_to`] as text.
    pub fn to_text(&self) -> String {
        let mut bytes = Vec::new();
        self.write_to(&mut bytes)
            .expect("writing to memory cannot fail");
        String::from_utf8(bytes).expect("the header and names are text")
    }
}

/// True for a name htsjdk's `IOUtil.openFileForWriting` compresses.
pub(super) fn gzip_name(path: &Path) -> bool {
    path.file_name().is_some_and(|name| {
        let name = name.to_string_lossy();
        [".gz", ".gzip", ".bfq"].iter().any(|e| name.ends_with(e))
    })
}

/// Refuse `target` when it is one of `inputs` (by name or by device and
/// inode), or a link or non-regular file another user owns.
pub(super) fn check_target(target: &Path, inputs: &[PathBuf]) -> Result<()> {
    if let Some(refused) = crate::cleanup::refusal(target) {
        return Err(Error::Output {
            path: target.to_path_buf(),
            source: refused,
        });
    }
    if let Some(missing) = crate::cleanup::unresolved_parent(target) {
        return Err(Error::Usage(format!(
            "{}: the directory {} does not exist and a `..` follows it",
            target.display(),
            missing.display()
        )));
    }
    if inputs
        .iter()
        .any(|input| crate::mark_duplicates::same_file(input, target))
    {
        return Err(Error::Usage(format!(
            "{} is an input; writing it would destroy it",
            target.display()
        )));
    }
    Ok(())
}

/// Write `list` into `file`: plain text, or BGZF (a gzip file htsjdk
/// reads back) under a name [`gzip_name`] accepts.
pub(super) fn stage(list: &IntervalList, file: &AtomicFile) -> Result<()> {
    let target = file.path().to_path_buf();
    let output = |source| Error::Output {
        path: target.clone(),
        source,
    };
    let handle = file.create().map_err(output)?;
    if gzip_name(&target) {
        let mut writer = noodles::bgzf::io::Writer::new(BufWriter::new(handle));
        list.write_to(&mut writer).map_err(output)?;
        writer
            .finish()
            .and_then(|mut inner| inner.flush())
            .map_err(output)?;
    } else {
        let mut writer = BufWriter::new(handle);
        list.write_to(&mut writer).map_err(output)?;
        writer.flush().map_err(output)?;
    }
    Ok(())
}

/// `IntervalList.write(path)` with the file rules of the crate: no input
/// is replaced, the file is written under a temporary name and renamed
/// into place (keeping the mode of a file it replaces), and a failed write
/// leaves `path` and its directory as they were.
pub fn write_list(list: &IntervalList, path: &Path, inputs: &[PathBuf]) -> Result<()> {
    check_target(path, inputs)?;
    let mut dirs = MadeDirs::default();
    dirs.create_parent(path).map_err(|source| Error::Output {
        path: path.to_path_buf(),
        source,
    })?;
    let file = AtomicFile::new(path);
    stage(list, &file)?;
    file.commit().map_err(|source| Error::Output {
        path: path.to_path_buf(),
        source,
    })?;
    dirs.keep();
    Ok(())
}

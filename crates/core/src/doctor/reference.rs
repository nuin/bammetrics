//! The reference FASTA given with `-R`, against the `@SQ` lines: names and
//! lengths from its `.fai`, MD5s from the `.dict` beside it or computed.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use noodles::core::Region;
use noodles::fasta;

use super::input;
use crate::md5::{Md5, hex};

/// The most bytes of a `.fai` or `.dict` read: 64 MiB.
pub const MAX_INDEX_BYTES: u64 = 64 << 20;

/// Bytes of FASTA read at a time while an MD5 is computed.
const CHUNK: usize = 1 << 20;

/// One `.fai` line: name, length, offset, bases per line, bytes per line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FaiEntry {
    pub name: String,
    pub length: u64,
    pub offset: u64,
    pub line_bases: u64,
    pub line_width: u64,
}

/// A whole small file, at most [`MAX_INDEX_BYTES`].
fn read_capped(path: &Path) -> io::Result<String> {
    let mut bytes = Vec::new();
    File::open(path)?
        .take(MAX_INDEX_BYTES + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_INDEX_BYTES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("over the {MAX_INDEX_BYTES}-byte cap"),
        ));
    }
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

fn invalid(what: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, what)
}

/// The `.fai` of `fasta` (`FASTA.fai`).
pub fn fai_path(fasta: &Path) -> PathBuf {
    let mut name = fasta.as_os_str().to_owned();
    name.push(".fai");
    PathBuf::from(name)
}

/// The entries of the `.fai` at `path`.
pub fn read_fai(path: &Path) -> io::Result<Vec<FaiEntry>> {
    let text = read_capped(path)?;
    text.lines()
        .enumerate()
        .filter(|(_, line)| !line.is_empty())
        .map(|(i, line)| {
            let fields: Vec<&str> = line.split('\t').collect();
            let number = |k: usize| -> io::Result<u64> {
                fields
                    .get(k)
                    .and_then(|f| f.trim().parse().ok())
                    .ok_or_else(|| invalid(format!("line {}: field {} is no number", i + 1, k + 1)))
            };
            Ok(FaiEntry {
                name: fields[0].to_owned(),
                length: number(1)?,
                offset: number(2)?,
                line_bases: number(3)?,
                line_width: number(4)?,
            })
        })
        .collect()
}

/// The `.dict` beside `fasta`: `ref.dict` for `ref.fa` (also `.fasta`,
/// `.fna`, each optionally `.gz`), else `ref.fa.dict`; the first that
/// exists.
pub fn dict_path(fasta: &Path) -> Option<PathBuf> {
    let name = fasta.file_name()?.to_str()?;
    let mut stem = name.strip_suffix(".gz").unwrap_or(name);
    for extension in [".fasta", ".fa", ".fna"] {
        if let Some(s) = stem.strip_suffix(extension) {
            stem = s;
            break;
        }
    }
    [format!("{stem}.dict"), format!("{name}.dict")]
        .into_iter()
        .map(|candidate| fasta.with_file_name(candidate))
        .find(|p| p.is_file())
}

/// The `M5` of each `@SQ` line of the `.dict` at `path`, by name.
pub fn read_dict_md5(path: &Path) -> io::Result<HashMap<String, String>> {
    let text = read_capped(path)?;
    Ok(input::lines(&text)
        .into_iter()
        .filter(|line| line.kind == "SQ")
        .filter_map(|line| Some((line.get("SN")?.to_owned(), line.get("M5")?.to_owned())))
        .collect())
}

/// The MD5 of the bases the `@SQ` `M5` field covers: upper case, every
/// byte outside `!` to `~` dropped (SAM specification 1.3).
fn digest_bases(md5: &mut Md5, bytes: &[u8]) {
    let kept: Vec<u8> = bytes
        .iter()
        .filter(|b| (b'!'..=b'~').contains(*b))
        .map(u8::to_ascii_uppercase)
        .collect();
    md5.update(&kept);
}

/// The MD5 of contig `entry` of the plain FASTA `file`, read 1 MiB at a
/// time.
fn plain_md5(file: &mut File, entry: &FaiEntry) -> io::Result<String> {
    if entry.line_bases == 0 || entry.length == 0 {
        return Ok(hex(&Md5::new().finish()));
    }
    // Every full line before the last with its line end, then the last
    // line's bases.
    let full_lines = (entry.length - 1) / entry.line_bases;
    let span = full_lines
        .saturating_mul(entry.line_width)
        .saturating_add(entry.length - full_lines * entry.line_bases);
    file.seek(SeekFrom::Start(entry.offset))?;
    let mut md5 = Md5::new();
    let mut left = span;
    let mut buffer = vec![0u8; CHUNK];
    while left > 0 {
        let n = CHUNK.min(usize::try_from(left).unwrap_or(CHUNK));
        file.read_exact(&mut buffer[..n])?;
        digest_bases(&mut md5, &buffer[..n]);
        left -= n as u64;
    }
    Ok(hex(&md5.finish()))
}

/// True for a FASTA read through the bgzip indexed reader.
fn compressed(fasta: &Path) -> bool {
    matches!(
        fasta.extension().and_then(|e| e.to_str()),
        Some("gz" | "bgz")
    )
}

/// A FASTA as one MD5 worker reads it.
enum Source {
    Plain(File),
    Indexed(fasta::io::IndexedReader<fasta::io::BufReader<File>>),
}

/// The MD5s of `names` in `fasta` (each in `entries`), on `threads`
/// workers. A bgzip FASTA is read through its `.gzi` one contig at a time.
pub fn compute_md5(
    fasta: &Path,
    entries: &[&FaiEntry],
    threads: usize,
) -> io::Result<HashMap<String, String>> {
    let next = AtomicUsize::new(0);
    let found = Mutex::new(HashMap::new());
    let failed: Mutex<Option<io::Error>> = Mutex::new(None);
    let work = || -> io::Result<()> {
        let mut source = if compressed(fasta) {
            Source::Indexed(
                fasta::io::indexed_reader::Builder::default()
                    .build_from_path(fasta)
                    .map_err(|e| invalid(format!("cannot open with its .fai and .gzi: {e}")))?,
            )
        } else {
            Source::Plain(File::open(fasta)?)
        };
        loop {
            let i = next.fetch_add(1, Ordering::Relaxed);
            let Some(entry) = entries.get(i) else {
                return Ok(());
            };
            let digest = match &mut source {
                Source::Plain(file) => plain_md5(file, entry)?,
                Source::Indexed(reader) => {
                    let record = reader.query(&Region::new(entry.name.as_str(), ..))?;
                    let mut md5 = Md5::new();
                    digest_bases(&mut md5, record.sequence().as_ref());
                    hex(&md5.finish())
                }
            };
            found
                .lock()
                .unwrap_or_else(|p| p.into_inner())
                .insert(entry.name.clone(), digest);
        }
    };
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| {
                if let Err(e) = work() {
                    // A failed worker stops the others at their next contig.
                    next.store(entries.len(), Ordering::Relaxed);
                    failed
                        .lock()
                        .unwrap_or_else(|p| p.into_inner())
                        .get_or_insert(e);
                }
            });
        }
    });
    if let Some(e) = failed.into_inner().unwrap_or_else(|p| p.into_inner()) {
        return Err(e);
    }
    Ok(found.into_inner().unwrap_or_else(|p| p.into_inner()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_md5_drops_line_ends_and_upper_cases() {
        let mut a = Md5::new();
        digest_bases(&mut a, b"acgt\nACGT\r\n");
        let mut b = Md5::new();
        b.update(b"ACGTACGT");
        assert_eq!(hex(&a.finish()), hex(&b.finish()));
    }

    #[test]
    fn the_dict_beside_a_fasta_is_found_by_its_stem() {
        let dir =
            std::env::temp_dir().join(format!("bammetrics-doctor-dict.{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let fasta = dir.join("ref.fa.gz");
        assert_eq!(dict_path(&fasta), None);
        std::fs::write(dir.join("ref.dict"), "@SQ\tSN:c\tLN:1\tM5:x\n").expect("dict");
        assert_eq!(dict_path(&fasta), Some(dir.join("ref.dict")));
        let md5 = read_dict_md5(&dir.join("ref.dict")).expect("read");
        assert_eq!(md5.get("c").map(String::as_str), Some("x"));
        std::fs::remove_dir_all(&dir).expect("clean");
    }
}

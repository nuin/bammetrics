//! The reference FASTA given with `-R`, against the `@SQ` lines: names and
//! lengths from its `.fai`, MD5s from the `.dict` beside it or computed.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use noodles::bgzf;

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

/// One `@SQ` line of a `.dict`: its `LN` (when a number) and its `M5`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DictEntry {
    pub length: Option<u64>,
    pub m5: Option<String>,
}

/// The `@SQ` lines of the `.dict` at `path`, by name.
pub fn read_dict(path: &Path) -> io::Result<HashMap<String, DictEntry>> {
    let text = read_capped(path)?;
    Ok(input::lines(&text)
        .into_iter()
        .filter(|line| line.kind == "SQ")
        .filter_map(|line| {
            let entry = DictEntry {
                length: line.get("LN").and_then(|v| v.parse().ok()),
                m5: line.get("M5").map(str::to_owned),
            };
            Some((line.get("SN")?.to_owned(), entry))
        })
        .collect())
}

/// The MD5 of the bases the `@SQ` `M5` field covers: upper case, every
/// byte outside `!` to `~` dropped (SAM specification 1.3). The bases are
/// compacted in place, so no second buffer is held.
fn digest_bases(md5: &mut Md5, bytes: &mut [u8]) {
    let mut kept = 0;
    for i in 0..bytes.len() {
        let b = bytes[i];
        if (b'!'..=b'~').contains(&b) {
            bytes[kept] = b.to_ascii_uppercase();
            kept += 1;
        }
    }
    md5.update(&bytes[..kept]);
}

/// The MD5 of contig `entry` of a FASTA read through `reader`, whose
/// positions are the uncompressed bytes the `.fai` offsets count (a plain
/// file, or a bgzip file through its `.gzi`), read 1 MiB at a time.
fn contig_md5<R: Read + Seek>(reader: &mut R, entry: &FaiEntry) -> io::Result<String> {
    if entry.line_bases == 0 || entry.length == 0 {
        return Ok(hex(&Md5::new().finish()));
    }
    // Every full line before the last with its line end, then the last
    // line's bases.
    let full_lines = (entry.length - 1) / entry.line_bases;
    let span = full_lines
        .saturating_mul(entry.line_width)
        .saturating_add(entry.length - full_lines * entry.line_bases);
    reader.seek(SeekFrom::Start(entry.offset))?;
    let mut md5 = Md5::new();
    let mut left = span;
    let mut buffer = vec![0u8; CHUNK];
    while left > 0 {
        let n = CHUNK.min(usize::try_from(left).unwrap_or(CHUNK));
        reader.read_exact(&mut buffer[..n])?;
        digest_bases(&mut md5, &mut buffer[..n]);
        left -= n as u64;
    }
    Ok(hex(&md5.finish()))
}

/// True for a FASTA read as bgzip through its `.gzi`.
fn compressed(fasta: &Path) -> bool {
    matches!(
        fasta.extension().and_then(|e| e.to_str()),
        Some("gz" | "bgz")
    )
}

/// The `.gzi` of `fasta` (`FASTA.gzi`).
fn gzi_path(fasta: &Path) -> PathBuf {
    let mut name = fasta.as_os_str().to_owned();
    name.push(".gzi");
    PathBuf::from(name)
}

/// The `.gzi` of a bgzip FASTA, at most [`MAX_INDEX_BYTES`].
pub(crate) fn read_gzi(fasta: &Path) -> io::Result<bgzf::gzi::Index> {
    let path = gzi_path(fasta);
    let size = std::fs::metadata(&path)
        .map_err(|e| invalid(format!("{}: {e} (bgzip -r writes it)", path.display())))?
        .len();
    if size > MAX_INDEX_BYTES {
        return Err(invalid(format!(
            "{}: over the {MAX_INDEX_BYTES}-byte cap",
            path.display()
        )));
    }
    bgzf::gzi::fs::read(&path).map_err(|e| invalid(format!("{}: {e}", path.display())))
}

/// A FASTA as one MD5 worker reads it: plain, or bgzip with the `.gzi`
/// that maps an uncompressed offset to its block. Both are read 1 MiB at a
/// time, so a worker holds 1 MiB of bases and one BGZF block whatever the
/// length of the contig.
enum Source {
    Plain(File),
    Bgzf(Box<bgzf::io::IndexedReader<BufReader<File>>>),
}

impl Source {
    fn open(fasta: &Path, gzi: Option<&bgzf::gzi::Index>) -> io::Result<Self> {
        let file = File::open(fasta)?;
        Ok(match gzi {
            Some(index) => Source::Bgzf(Box::new(bgzf::io::IndexedReader::new(
                BufReader::new(file),
                index.clone(),
            ))),
            None => Source::Plain(file),
        })
    }

    fn md5(&mut self, entry: &FaiEntry) -> io::Result<String> {
        match self {
            Source::Plain(file) => contig_md5(file, entry),
            Source::Bgzf(reader) => contig_md5(reader.as_mut(), entry),
        }
    }
}

/// The MD5s of `names` in `fasta` (each in `entries`), on `threads`
/// workers, one contig per worker at a time.
pub fn compute_md5(
    fasta: &Path,
    entries: &[&FaiEntry],
    threads: usize,
) -> io::Result<HashMap<String, String>> {
    let gzi = if compressed(fasta) {
        Some(read_gzi(fasta)?)
    } else {
        None
    };
    let next = AtomicUsize::new(0);
    let found = Mutex::new(HashMap::new());
    let failed: Mutex<Option<io::Error>> = Mutex::new(None);
    let work = || -> io::Result<()> {
        let mut source = Source::open(fasta, gzi.as_ref())?;
        loop {
            let i = next.fetch_add(1, Ordering::Relaxed);
            let Some(entry) = entries.get(i) else {
                return Ok(());
            };
            let digest = source.md5(entry)?;
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
        digest_bases(&mut a, &mut b"acgt\nACGT\r\n".to_owned());
        let mut b = Md5::new();
        b.update(b"ACGTACGT");
        assert_eq!(hex(&a.finish()), hex(&b.finish()));
    }

    /// The `.gzi` entries of BGZF bytes: (compressed, uncompressed) offset
    /// of every block after the first.
    fn gzi_entries(bytes: &[u8]) -> Vec<(u64, u64)> {
        let (mut at, mut inflated, mut entries) = (0usize, 0u64, Vec::new());
        while at < bytes.len() {
            if at > 0 {
                entries.push((at as u64, inflated));
            }
            let size = usize::from(u16::from_le_bytes([bytes[at + 16], bytes[at + 17]])) + 1;
            let isize = &bytes[at + size - 4..at + size];
            inflated += u64::from(u32::from_le_bytes(isize.try_into().expect("4 bytes")));
            at += size;
        }
        entries
    }

    #[test]
    fn bgzip_and_plain_fastas_give_the_same_md5s() {
        use std::io::Write as _;
        let dir =
            std::env::temp_dir().join(format!("bammetrics-doctor-md5.{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        // Three contigs over about 250 KB: several BGZF blocks, contigs
        // that start and end inside a block, a last line shorter than 60.
        let mut text = String::new();
        let mut fai = String::new();
        let mut wanted = Vec::new();
        for (name, length) in [("a", 120_007u64), ("b", 59), ("c", 130_000)] {
            text.push_str(&format!(">{name}\n"));
            let offset = text.len() as u64;
            let bases: String = (0..length)
                .map(|i| ['a', 'C', 'g', 'T', 'N'][(i * 7 % 5) as usize])
                .collect();
            for line in bases.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(line).expect("ascii"));
                text.push('\n');
            }
            fai.push_str(&format!("{name}\t{length}\t{offset}\t60\t61\n"));
            let mut md5 = Md5::new();
            md5.update(bases.to_ascii_uppercase().as_bytes());
            wanted.push((name.to_owned(), hex(&md5.finish())));
        }
        let plain = dir.join("ref.fa");
        std::fs::write(&plain, &text).expect("write");
        std::fs::write(fai_path(&plain), &fai).expect("write");
        let gz = dir.join("ref.fa.gz");
        let mut writer = bgzf::io::Writer::new(File::create(&gz).expect("create"));
        writer.write_all(text.as_bytes()).expect("write");
        writer.finish().expect("finish");
        let bytes = std::fs::read(&gz).expect("read");
        bgzf::gzi::fs::write(gzi_path(&gz), &bgzf::gzi::Index::from(gzi_entries(&bytes)))
            .expect("write the .gzi");
        let entries = read_fai(&fai_path(&plain)).expect("fai");
        let refs: Vec<&FaiEntry> = entries.iter().collect();
        for (fasta, threads) in [(&plain, 1), (&gz, 1), (&gz, 3)] {
            let found = compute_md5(fasta, &refs, threads).expect("md5");
            for (name, md5) in &wanted {
                assert_eq!(found.get(name), Some(md5), "{} {name}", fasta.display());
            }
        }
        std::fs::remove_file(gzi_path(&gz)).expect("remove");
        let missing = compute_md5(&gz, &refs, 1).expect_err("no .gzi");
        assert!(missing.to_string().contains("ref.fa.gz.gzi"), "{missing}");
        std::fs::remove_dir_all(&dir).expect("clean");
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
        let dict = read_dict(&dir.join("ref.dict")).expect("read");
        let entry = dict.get("c").expect("c");
        assert_eq!((entry.length, entry.m5.as_deref()), (Some(1), Some("x")));
        std::fs::remove_dir_all(&dir).expect("clean");
    }
}

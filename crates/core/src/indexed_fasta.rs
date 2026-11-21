//! A FASTA as htsjdk 5.0.0's `ReferenceSequenceFileFactory.getReferenceSequenceFile`
//! opens it for `ScatterIntervalsByNs` and `NonNFastaSize`: the `.fai`
//! beside it gives each contig's offset and line geometry, the `.dict`
//! beside it the sequence dictionary, and a contig's bases are the bytes
//! `IndexedFastaSequenceFile.getSubsequenceAt` takes at those offsets.
//!
//! What htsjdk does, and this module with it:
//!
//! - The name ends in one of [`FASTA_EXTENSIONS`], else `File is not a
//!   supported reference file type`.
//! - The `.dict` is the FASTA's name with its extension replaced by
//!   `.dict`, else the name plus `.dict`. It is decoded by
//!   `SAMTextHeaderCodec` at `SILENT` ([`Header::decode`]); one with no
//!   `@SQ` line counts as none.
//! - The `.fai` is `FASTA.fai`; a BGZF FASTA needs `FASTA.gzi` too. Without
//!   them the FASTA is not indexed ([`Opened::Unindexed`]).
//! - Each `.fai` line holds `(.+)\t+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)`, found
//!   by a `Scanner` ([`read_index`]): a name (cut at its first whitespace),
//!   a length, an offset, bases per line and bytes per line; a name given
//!   twice fails.
//! - With a dictionary, its `@SQ` lines must name the index's contigs in
//!   the index's order with the same lengths
//!   (`sanityCheckDictionaryAgainstIndex`).
//! - Base `i` (from 0) of a contig is the byte at `offset + (i / bases) *
//!   width + i % bases`, whatever that byte is. Bases per line of 0 fails
//!   (`/ by zero`), also for a contig of no base.
//!
//! Where this module differs (KNOWN_DIFFERENCES.md):
//!
//! - A FASTA that ends before a contig's last base fails. htsjdk returns
//!   zero bytes for the missing bases, which count as bases.
//! - A gzip FASTA that is not BGZF, or a BGZF FASTA whose name does not end
//!   in `.gz`, is refused. htsjdk reads the compressed bytes as bases.
//! - A FASTA or `.dict` that is not a regular file (a named pipe, a
//!   socket, a device) fails before it is opened.
//!
//! Memory: one [`CHUNK`] buffer whatever the contig length, and the `.fai`
//! and `.dict` read whole, each under [`MAX_INDEX_BYTES`].

use std::collections::HashSet;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use noodles::bgzf;

use crate::create_sequence_dictionary::{FASTA_EXTENSIONS, default_dictionary};
use crate::error::{Error, Result};
use crate::interval_ops::Header;

/// Bytes of FASTA read at a time.
pub const CHUNK: usize = 1 << 20;
/// The most bytes of a `.fai` or `.dict` read: 64 MiB.
pub const MAX_INDEX_BYTES: u64 = crate::doctor::reference::MAX_INDEX_BYTES;

/// One `.fai` line, with htsjdk's types: `long` length and offset, `int`
/// bases and bytes per line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Entry {
    pub name: String,
    pub length: i64,
    pub offset: i64,
    pub line_bases: i32,
    pub line_width: i32,
}

/// The FASTA's bytes: plain, or BGZF read at uncompressed offsets through
/// its `.gzi`.
enum Source {
    Plain(File),
    Bgzf(Box<bgzf::io::IndexedReader<BufReader<File>>>),
}

impl Source {
    fn seek_to(&mut self, offset: u64) -> io::Result<()> {
        match self {
            Source::Plain(file) => file.seek(SeekFrom::Start(offset)).map(drop),
            Source::Bgzf(reader) => reader.seek(SeekFrom::Start(offset)).map(drop),
        }
    }

    fn read_exact(&mut self, buffer: &mut [u8]) -> io::Result<()> {
        match self {
            Source::Plain(file) => file.read_exact(buffer),
            Source::Bgzf(reader) => reader.read_exact(buffer),
        }
    }
}

/// An indexed FASTA: its `.fai` entries, its dictionary when it has one,
/// and the file open for the walk.
pub struct IndexedFasta {
    path: PathBuf,
    entries: Vec<Entry>,
    dictionary: Option<Header>,
    source: Source,
    buffer: Vec<u8>,
}

/// What `getReferenceSequenceFile` returns: an indexed FASTA, or, without
/// a `.fai` (or a `.gzi` for BGZF), a FASTA whose only use is its
/// dictionary.
pub enum Opened {
    Indexed(IndexedFasta),
    Unindexed { dictionary: Option<Header> },
}

/// The `.fai` of `fasta`.
pub fn fai_path(fasta: &Path) -> PathBuf {
    with_suffix(fasta, ".fai")
}

fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

/// `AbstractFastaSequenceFile.findSequenceDictionary`: `NAME.dict` with
/// the FASTA extension replaced, else `FASTA.dict`; `None` when neither
/// exists. A name without a FASTA extension fails.
pub fn dictionary_path(fasta: &Path) -> Result<Option<PathBuf>> {
    let replaced = default_dictionary(fasta)?;
    if replaced.exists() {
        return Ok(Some(replaced));
    }
    let appended = with_suffix(fasta, ".dict");
    Ok(appended.exists().then_some(appended))
}

/// Every file [`open`] may read for `fasta`: the FASTA, its `.fai`, its
/// `.gzi` and its `.dict`. No output may be one of them.
pub fn files_read(fasta: &Path) -> Vec<PathBuf> {
    let mut files = vec![
        fasta.to_path_buf(),
        fai_path(fasta),
        with_suffix(fasta, ".gzi"),
    ];
    files.extend(default_dictionary(fasta).ok());
    files.push(with_suffix(fasta, ".dict"));
    files
}

fn reference_error(path: &Path, message: String) -> Error {
    Error::Reference {
        path: path.to_path_buf(),
        message,
    }
}

/// Refuse a path that is not a regular file (after links), before it is
/// opened: opening a named pipe with no writer blocks, and the walk seeks.
fn refuse_irregular(path: &Path) -> Result<()> {
    let metadata = std::fs::metadata(path).map_err(|source| Error::Input {
        path: path.to_path_buf(),
        source,
    })?;
    if metadata.is_dir() {
        return Err(Error::Input {
            path: path.to_path_buf(),
            source: io::Error::other("Cannot read file because it is a directory"),
        });
    }
    if !metadata.is_file() {
        return Err(reference_error(
            path,
            "not a regular file (a pipe, socket or device); the reader seeks in it".to_owned(),
        ));
    }
    Ok(())
}

/// A whole small file, at most [`MAX_INDEX_BYTES`]; bytes that are not
/// UTF-8 become U+FFFD.
fn read_capped(path: &Path) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    refuse_irregular(path)?;
    let mut bytes = Vec::new();
    File::open(path)
        .map_err(input)?
        .take(MAX_INDEX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(input)?;
    if bytes.len() as u64 > MAX_INDEX_BYTES {
        return Err(reference_error(
            path,
            format!("over the {MAX_INDEX_BYTES}-byte cap"),
        ));
    }
    Ok(String::from_utf8_lossy(&bytes).into_owned())
}

/// `loadDictionary`: the leading `@` lines of the `.dict`, decoded at
/// `SILENT`.
fn read_dictionary(path: &Path) -> Result<Header> {
    let text = read_capped(path)?;
    let lines: Vec<_> = crate::interval_ops::java_lines(text.as_bytes())
        .take_while(|l| l.starts_with('@'))
        .collect();
    Header::decode(lines.iter().map(|l| l.as_ref()))
        .map_err(|message| reference_error(path, message))
}

/// `Character.isWhitespace`: the delimiters `Scanner.hasNext` skips.
fn java_whitespace(c: char) -> bool {
    matches!(c, '\t'..='\r' | '\u{1c}'..='\u{1f}')
        || (c.is_whitespace() && !matches!(c, '\u{85}' | '\u{a0}' | '\u{2007}' | '\u{202f}'))
}

/// The line at the head of `text` and where the next one starts, as
/// `Scanner` splits lines: at `\r\n`, `\n`, `\r`, U+2028, U+2029 or
/// U+0085. `None` when no separator ends the line.
fn scanner_line(text: &str) -> (&str, Option<usize>) {
    let Some(at) = text.find(['\n', '\r', '\u{2028}', '\u{2029}', '\u{85}']) else {
        return (text, None);
    };
    let separator = if text[at..].starts_with("\r\n") {
        2
    } else {
        text[at..].chars().next().map_or(1, char::len_utf8)
    };
    (&text[..at], Some(at + separator))
}

/// The pattern `parseIndexFile` finds in each line,
/// `(.+)\t+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)`, with Java's ASCII `\d` and `\s`.
static INDEX_LINE: std::sync::LazyLock<regex::Regex> = std::sync::LazyLock::new(|| {
    let space = r"[\t\n\x0B\x0C\r ]+";
    regex::Regex::new(&format!(
        r"(?s)(.+)\t+([0-9]+){space}([0-9]+){space}([0-9]+){space}([0-9]+)"
    ))
    .expect("a valid pattern")
});

/// `SAMSequenceRecord.truncateSequenceName`: the name up to its first
/// space, tab, line feed, vertical tab, form feed or carriage return.
fn truncate_sequence_name(name: &str) -> &str {
    name.split(['\t', '\n', '\x0B', '\x0C', '\r', ' '])
        .next()
        .unwrap_or_default()
}

/// `FastaSequenceIndex.parseIndexFile`, a `Scanner` over the text. While
/// a token is left, the pattern is found anywhere in the current line and
/// the rest of the line is skipped. A line without it, a blank line before
/// a later entry included, is invalid. A last entry whose match runs to the
/// end of a file with no final line end fails (`nextLine` finds no line).
pub fn read_index(path: &Path) -> Result<Vec<Entry>> {
    let text = read_capped(path)?;
    let mut entries = Vec::new();
    let mut names = HashSet::new();
    let mut rest = text.as_str();
    while !rest.chars().all(java_whitespace) {
        let (line, next) = scanner_line(rest);
        let invalid = || reference_error(path, format!("Found invalid line in index file:{line}"));
        let found = INDEX_LINE.captures(line).ok_or_else(invalid)?;
        if next.is_none() && found.get(0).is_some_and(|m| m.end() == line.len()) {
            return Err(reference_error(
                path,
                "No line found: the last line has no line end (htsjdk's Scanner.nextLine fails)"
                    .to_owned(),
            ));
        }
        let group = |i: usize| found.get(i).map_or("", |m| m.as_str());
        let number = |text: &str, bits| {
            crate::java::parse_integer(text, bits).map_err(|e| reference_error(path, e))
        };
        let name = truncate_sequence_name(group(1));
        if !names.insert(name.to_owned()) {
            return Err(reference_error(
                path,
                format!("Contig '{name}' already exists in fasta index."),
            ));
        }
        entries.push(Entry {
            name: name.to_owned(),
            length: number(group(2), 64)?,
            offset: number(group(3), 64)?,
            line_bases: number(group(4), 32)? as i32,
            line_width: number(group(5), 32)? as i32,
        });
        rest = next.map_or("", |at| &rest[at..]);
    }
    Ok(entries)
}

/// True for the first block of a BGZF file: gzip with the `BC` extra
/// field (`IOUtil.isBlockCompressed`). `None` for a file that is not gzip.
fn sniff_gzip(path: &Path) -> Result<Option<bool>> {
    let mut head = Vec::with_capacity(18);
    File::open(path)
        .and_then(|f| f.take(18).read_to_end(&mut head))
        .map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
    if !head.starts_with(&[0x1f, 0x8b]) {
        return Ok(None);
    }
    let bgzf = head.len() == 18
        && head[2] == 8
        && head[3] & 4 != 0
        && u16::from_le_bytes([head[10], head[11]]) >= 6
        && head[12] == b'B'
        && head[13] == b'C';
    Ok(Some(bgzf))
}

/// `sanityCheckDictionaryAgainstIndex`.
fn check_dictionary(fasta: &Path, dictionary: &Header, entries: &[Entry]) -> Result<()> {
    if dictionary.sequences.len() != entries.len() {
        return Err(reference_error(
            fasta,
            "Sequence dictionary and index contain different numbers of contigs".to_owned(),
        ));
    }
    for (sequence, entry) in dictionary.sequences.iter().zip(entries) {
        if sequence.name != entry.name {
            return Err(reference_error(
                fasta,
                format!(
                    "Mismatch between sequence dictionary fasta index for {}, sequence '{}' != '{}'.",
                    fasta.display(),
                    sequence.name,
                    entry.name
                ),
            ));
        }
        if i64::from(sequence.length) != entry.length {
            return Err(reference_error(
                fasta,
                format!(
                    "Index length does not match dictionary length for contig: {}",
                    sequence.name
                ),
            ));
        }
    }
    Ok(())
}

/// `ReferenceSequenceFileFactory.getReferenceSequenceFile(path, true, true)`.
pub fn open(path: &Path) -> Result<Opened> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    refuse_irregular(path)?;
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    if !FASTA_EXTENSIONS.iter().any(|e| name.ends_with(e)) {
        return Err(reference_error(
            path,
            "File is not a supported reference file type".to_owned(),
        ));
    }
    // htsjdk keeps a dictionary with no `@SQ` line as no dictionary.
    let dictionary = dictionary_path(path)?
        .map(|dict| read_dictionary(&dict))
        .transpose()?
        .filter(|header| !header.sequences.is_empty());
    let bgzf = match sniff_gzip(path)? {
        None => false,
        Some(true) if name.ends_with(".gz") => true,
        Some(_) => {
            return Err(reference_error(
                path,
                "a gzip FASTA is read through its .fai only when it is BGZF (bgzip) and named \
                 .gz; htsjdk would read the compressed bytes as bases"
                    .to_owned(),
            ));
        }
    };
    let fai = fai_path(path);
    let gzi = with_suffix(path, ".gzi");
    if !fai.is_file() || (bgzf && !gzi.is_file()) {
        return Ok(Opened::Unindexed { dictionary });
    }
    let entries = read_index(&fai)?;
    if let Some(dictionary) = &dictionary {
        check_dictionary(path, dictionary, &entries)?;
    }
    let file = File::open(path).map_err(input)?;
    let source = if bgzf {
        let index = crate::doctor::reference::read_gzi(path)
            .map_err(|e| reference_error(&gzi, e.to_string()))?;
        Source::Bgzf(Box::new(bgzf::io::IndexedReader::new(
            BufReader::new(file),
            index,
        )))
    } else {
        Source::Plain(file)
    };
    Ok(Opened::Indexed(IndexedFasta {
        path: path.to_path_buf(),
        entries,
        dictionary,
        source,
        buffer: Vec::new(),
    }))
}

impl IndexedFasta {
    /// The `.fai` entries in their order.
    pub fn entries(&self) -> &[Entry] {
        &self.entries
    }

    /// The `.dict`, when there is one.
    pub fn dictionary(&self) -> Option<&Header> {
        self.dictionary.as_ref()
    }

    /// Hand every base of contig `index` to `take`, in order, a line's run
    /// of bases at a time, reading [`CHUNK`] bytes at a time: the bases
    /// `getSequence` returns, before `toUpperCase`.
    pub fn walk(&mut self, index: usize, mut take: impl FnMut(&[u8])) -> Result<()> {
        let entry = self.entries[index].clone();
        let malformed = |what: &str| {
            reference_error(
                &self.path,
                format!("contig {} in the .fai: {what}", entry.name),
            )
        };
        if entry.line_bases == 0 {
            return Err(malformed("0 bases per line (htsjdk fails: / by zero)"));
        }
        if entry.line_bases < 0 || entry.line_width < entry.line_bases || entry.offset < 0 {
            return Err(malformed(
                "a negative number, or fewer bytes than bases per line",
            ));
        }
        if entry.length > i64::from(i32::MAX) {
            return Err(malformed("longer than a Java array holds"));
        }
        if entry.length <= 0 {
            return Ok(());
        }
        let (length, bases, width) = (
            entry.length as u64,
            entry.line_bases as u64,
            entry.line_width as u64,
        );
        let full_lines = (length - 1) / bases;
        let span = full_lines * width + (length - full_lines * bases);
        let Self {
            path,
            source,
            buffer,
            ..
        } = self;
        let failed = |e: io::Error| {
            if e.kind() == io::ErrorKind::UnexpectedEof {
                reference_error(
                    path,
                    format!(
                        "the FASTA ends before the {length} bases of contig {} its .fai gives",
                        entry.name
                    ),
                )
            } else {
                Error::Input {
                    path: path.clone(),
                    source: e,
                }
            }
        };
        source.seek_to(entry.offset as u64).map_err(failed)?;
        buffer.resize(CHUNK.min(span as usize), 0);
        let mut column = 0u64;
        let mut left = span;
        while left > 0 {
            let n = buffer.len().min(left as usize);
            source.read_exact(&mut buffer[..n]).map_err(failed)?;
            let mut at = 0;
            while at < n {
                let room = (n - at) as u64;
                if column < bases {
                    let run = (bases - column).min(room) as usize;
                    take(&buffer[at..at + run]);
                    at += run;
                    column += run as u64;
                } else {
                    let skip = (width - column).min(room) as usize;
                    at += skip;
                    column += skip as u64;
                }
                if column == width {
                    column = 0;
                }
            }
            left -= n as u64;
        }
        Ok(())
    }
}

/// FASTA fixtures for the unit tests of the tools on this module.
#[cfg(test)]
pub(crate) mod testing {
    use std::path::{Path, PathBuf};

    /// A fresh directory under the system temporary directory.
    pub fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "bammetrics-indexed-fasta-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch directory");
        dir
    }

    /// `dir/name` with `contigs` wrapped at `width` bases and lines ended
    /// by `end`, its `.fai` as samtools writes it, and `NAME.dict` with
    /// `SN` and `LN` per contig.
    pub fn write_fasta(
        dir: &Path,
        name: &str,
        contigs: &[(&str, &[u8])],
        width: usize,
        end: &str,
    ) -> PathBuf {
        let mut fasta = Vec::new();
        let mut fai = String::new();
        let mut dict = String::from("@HD\tVN:1.6\n");
        for (contig, bases) in contigs {
            fasta.extend_from_slice(format!(">{contig} description\n").as_bytes());
            let offset = fasta.len();
            for line in bases.chunks(width) {
                fasta.extend_from_slice(line);
                fasta.extend_from_slice(end.as_bytes());
            }
            fai.push_str(&format!(
                "{contig}\t{}\t{offset}\t{width}\t{}\n",
                bases.len(),
                width + end.len()
            ));
            dict.push_str(&format!("@SQ\tSN:{contig}\tLN:{}\n", bases.len()));
        }
        let path = dir.join(name);
        std::fs::write(&path, fasta).expect("FASTA");
        std::fs::write(super::fai_path(&path), fai).expect(".fai");
        let dict_path = super::default_dictionary(&path).expect("a FASTA name");
        std::fs::write(dict_path, dict).expect(".dict");
        path
    }
}

#[cfg(test)]
mod tests;

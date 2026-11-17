//! What `doctor` reads of the file itself: its kind, its last bytes and its
//! header text, line by line, with no parser that refuses a header.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;

use noodles::bgzf;

use crate::bam_raw::read_append;
use crate::cram::{self, Kind};
use crate::gzip::GzipReader;

/// The most header text `doctor` reads: 64 MiB. A GRCh38 header with alt,
/// decoy and HLA contigs is about 200 KiB.
pub const MAX_HEADER_BYTES: usize = 64 << 20;

/// The BGZF end-of-file block (SAM specification 4.1.2).
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// The CRAM 3.x end-of-file container (CRAM specification 9).
pub const CRAM_EOF: [u8; 38] = [
    0x0f, 0x00, 0x00, 0x00, 0xff, 0xff, 0xff, 0xff, 0x0f, 0xe0, 0x45, 0x4f, 0x46, 0x00, 0x00, 0x00,
    0x00, 0x01, 0x00, 0x05, 0xbd, 0xd9, 0x4f, 0x00, 0x01, 0x00, 0x06, 0x06, 0x01, 0x00, 0x01, 0x00,
    0x01, 0x00, 0xee, 0x63, 0x01, 0x4b,
];

/// The largest BGZF block: 64 KiB.
const MAX_BLOCK: u64 = 1 << 16;

/// The file kind, or the reason it has none.
pub(crate) fn kind(path: &Path) -> io::Result<std::result::Result<Kind, String>> {
    let file = File::open(path)?;
    Ok(cram::sniff(file).map(|s| s.kind).map_err(|e| e.to_string()))
}

/// The last `n` bytes of the file (fewer for a shorter file).
fn tail(path: &Path, n: u64) -> io::Result<Vec<u8>> {
    let mut file = File::open(path)?;
    let size = file.metadata()?.len();
    let start = size.saturating_sub(n);
    file.seek(SeekFrom::Start(start))?;
    let mut bytes = Vec::new();
    file.take(n).read_to_end(&mut bytes)?;
    Ok(bytes)
}

/// True when the first bytes of the file are a BGZF block header (a gzip
/// member with the `BC` extra field): BAM and bgzip text.
pub fn is_bgzf(path: &Path) -> io::Result<bool> {
    let mut head = Vec::new();
    File::open(path)?.take(18).read_to_end(&mut head)?;
    Ok(block_size_at(&head, 0).is_some())
}

/// The size of the BGZF block whose header starts at `at`, when one does.
fn block_size_at(bytes: &[u8], at: usize) -> Option<usize> {
    let head = bytes.get(at..at + 18)?;
    let gzip = head[..4] == [0x1f, 0x8b, 0x08, 0x04];
    let extra = head[10..12] == [6, 0] && head[12..16] == [b'B', b'C', 2, 0];
    (gzip && extra).then(|| usize::from(u16::from_le_bytes([head[16], head[17]])) + 1)
}

/// True when `bytes` start a BGZF block, or are a cut start of one (a
/// prefix of the 4-byte gzip magic and flags every block begins with).
fn starts_block(bytes: &[u8]) -> bool {
    if bytes.len() >= 18 {
        return block_size_at(bytes, 0).is_some();
    }
    BGZF_EOF[..4].starts_with(&bytes[..bytes.len().min(4)])
}

/// How a BGZF file ends.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BgzfEnd {
    /// The 28-byte end-of-file block.
    Marker,
    /// A whole block, no marker: the writer stopped, or the file was cut at
    /// a block boundary.
    WholeBlock,
    /// The last block is cut.
    Cut,
    /// The end-of-file block, then this many bytes that start no block.
    Trailing(u64),
}

/// How the BGZF file at `path` ends, from its last 64 KiB and 28 bytes: a
/// block that starts in that window and ends at the last byte is whole.
/// An end-of-file block followed by bytes that start no BGZF block is a
/// file with bytes after its end; one followed by a block header is a
/// concatenation, judged by its last block.
pub fn bgzf_end(path: &Path) -> io::Result<BgzfEnd> {
    let bytes = tail(path, MAX_BLOCK + BGZF_EOF.len() as u64)?;
    if bytes.ends_with(&BGZF_EOF) {
        return Ok(BgzfEnd::Marker);
    }
    if let Some(at) = bytes
        .windows(BGZF_EOF.len())
        .rposition(|w| w == BGZF_EOF)
        .map(|at| at + BGZF_EOF.len())
        && !starts_block(&bytes[at..])
    {
        return Ok(BgzfEnd::Trailing((bytes.len() - at) as u64));
    }
    let whole = (0..bytes.len()).any(|at| block_size_at(&bytes, at) == Some(bytes.len() - at));
    Ok(if whole {
        BgzfEnd::WholeBlock
    } else {
        BgzfEnd::Cut
    })
}

/// True when the CRAM at `path` ends with the end-of-file container.
pub fn cram_has_eof(path: &Path) -> io::Result<bool> {
    Ok(tail(path, CRAM_EOF.len() as u64)? == CRAM_EOF)
}

/// True when the plain text at `path` ends with a line end (or is empty).
pub fn text_ends_with_newline(path: &Path) -> io::Result<bool> {
    let bytes = tail(path, 1)?;
    Ok(bytes.is_empty() || bytes == b"\n")
}

/// The header text and what else the file's container says of it.
pub struct RawHeader {
    pub text: String,
    /// A BAM's binary reference list, name and length.
    pub binary: Vec<(String, u64)>,
    /// A CRAM's major and minor version.
    pub cram_version: Option<(u8, u8)>,
}

/// `size` header bytes read from `reader`, or an error past the cap.
fn capped_text<R: Read + ?Sized>(reader: &mut R, size: usize) -> io::Result<String> {
    if size > MAX_HEADER_BYTES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("the header text claims {size} bytes, over the {MAX_HEADER_BYTES}-byte cap"),
        ));
    }
    let mut text = Vec::new();
    read_append(reader, &mut text, size)?;
    // The text ends at its first NUL (htsjdk and samtools pad it).
    if let Some(nul) = text.iter().position(|&b| b == 0) {
        text.truncate(nul);
    }
    Ok(String::from_utf8_lossy(&text).into_owned())
}

fn read_u32<R: Read + ?Sized>(reader: &mut R) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

/// The inflated stream of a BAM.
pub type BamStream = bgzf::io::Reader<BufReader<File>>;

/// The header of the BAM at `path` and its stream, left at the first
/// record.
pub fn bam_header(path: &Path) -> io::Result<(RawHeader, BamStream)> {
    let mut reader = bgzf::io::Reader::new(BufReader::new(File::open(path)?));
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if magic != *b"BAM\x01" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "the first block does not start with the BAM magic",
        ));
    }
    let l_text = read_u32(&mut reader)? as usize;
    let text = capped_text(&mut reader, l_text)?;
    let count = read_u32(&mut reader)?;
    let mut binary = Vec::new();
    for _ in 0..count {
        let l_name = read_u32(&mut reader)? as usize;
        let mut name = Vec::new();
        read_append(&mut reader, &mut name, l_name)?;
        if name.pop() != Some(0) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "a reference name of the binary list does not end in NUL",
            ));
        }
        let length = u64::from(read_u32(&mut reader)?);
        binary.push((String::from_utf8_lossy(&name).into_owned(), length));
    }
    let header = RawHeader {
        text,
        binary,
        cram_version: None,
    };
    Ok((header, reader))
}

/// The header of the CRAM at `path`.
pub fn cram_header(path: &Path) -> io::Result<RawHeader> {
    let mut reader = noodles::cram::io::Reader::new(BufReader::new(File::open(path)?));
    let mut header = reader.header_reader();
    header.read_magic_number()?;
    let version = header.read_format_version()?;
    header.read_file_id()?;
    let mut container = header.container_reader()?;
    let mut raw = container.raw_sam_header_reader()?;
    let mut text = Vec::new();
    (&mut raw)
        .take(MAX_HEADER_BYTES as u64 + 1)
        .read_to_end(&mut text)?;
    if text.len() > MAX_HEADER_BYTES {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("the header text is over the {MAX_HEADER_BYTES}-byte cap"),
        ));
    }
    Ok(RawHeader {
        text: String::from_utf8_lossy(&text).into_owned(),
        binary: Vec::new(),
        cram_version: Some((version.major(), version.minor())),
    })
}

/// The header lines of the SAM text at `path` (plain or gzip): every line
/// before the first that does not start with `@`.
pub fn sam_header(path: &Path, gzip: bool) -> io::Result<RawHeader> {
    let file = File::open(path)?;
    let mut reader: Box<dyn BufRead> = if gzip {
        Box::new(BufReader::new(GzipReader::new(BufReader::new(file))))
    } else {
        Box::new(BufReader::new(file))
    };
    let mut text = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        let n = (&mut reader)
            .take(MAX_HEADER_BYTES as u64 + 1)
            .read_until(b'\n', &mut line)?;
        if n == 0 || line.first() != Some(&b'@') {
            break;
        }
        text.extend_from_slice(&line);
        if text.len() > MAX_HEADER_BYTES {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("the header text is over the {MAX_HEADER_BYTES}-byte cap"),
            ));
        }
    }
    Ok(RawHeader {
        text: String::from_utf8_lossy(&text).into_owned(),
        binary: Vec::new(),
        cram_version: None,
    })
}

/// One header line: its record type (`HD`, `SQ`, `RG`, `PG`, `CO`) and its
/// `TAG:value` fields in order.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Line {
    pub kind: String,
    pub fields: Vec<(String, String)>,
}

impl Line {
    /// The value of `tag`, the first when the line repeats it.
    pub fn get(&self, tag: &str) -> Option<&str> {
        self.fields
            .iter()
            .find(|(t, _)| t == tag)
            .map(|(_, v)| v.as_str())
    }
}

/// The lines of a header text. A `@CO` line keeps its text as one field.
pub fn lines(text: &str) -> Vec<Line> {
    text.lines()
        .map(|line| line.strip_suffix('\r').unwrap_or(line))
        .filter(|line| !line.is_empty())
        .map(|line| {
            let mut parts = line.split('\t');
            let kind = parts
                .next()
                .unwrap_or_default()
                .trim_start_matches('@')
                .to_owned();
            let fields = if kind == "CO" {
                Vec::new()
            } else {
                parts
                    .map(|field| match field.split_once(':') {
                        Some((tag, value)) => (tag.to_owned(), value.to_owned()),
                        None => (field.to_owned(), String::new()),
                    })
                    .collect()
            };
            Line { kind, fields }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_lines_split_into_tags() {
        let text = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5\r\n@CO\tfree: text\n\n";
        let lines = lines(text);
        assert_eq!(lines.len(), 3);
        assert_eq!(lines[0].get("SO"), Some("coordinate"));
        assert_eq!(lines[1].get("LN"), Some("5"));
        assert_eq!(lines[2].kind, "CO");
        assert!(lines[2].fields.is_empty());
    }

    #[test]
    fn bytes_after_the_end_of_file_block_are_told_from_a_cut_block() {
        let dir =
            std::env::temp_dir().join(format!("bammetrics-doctor-end.{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let path = dir.join("x.bam");
        let end = |bytes: &[u8]| {
            std::fs::write(&path, bytes).expect("write");
            bgzf_end(&path).expect("read")
        };
        let eof = BGZF_EOF.as_slice();
        assert_eq!(end(&[eof, eof].concat()), BgzfEnd::Marker);
        assert_eq!(end(&[eof, b"junk"].concat()), BgzfEnd::Trailing(4));
        assert_eq!(end(&[eof, b"\n"].concat()), BgzfEnd::Trailing(1));
        // A second file cut inside its first block header, or after it.
        assert_eq!(end(&[eof, &eof[..3]].concat()), BgzfEnd::Cut);
        assert_eq!(end(&[eof, &eof[..20]].concat()), BgzfEnd::Cut);
        assert_eq!(end(&eof[..27]), BgzfEnd::Cut);
        std::fs::remove_dir_all(&dir).expect("clean");
    }

    #[test]
    fn a_block_header_is_found_by_its_bc_field() {
        assert_eq!(block_size_at(&BGZF_EOF, 0), Some(28));
        assert_eq!(block_size_at(&BGZF_EOF, 1), None);
        assert_eq!(block_size_at(&BGZF_EOF[..17], 0), None);
    }
}

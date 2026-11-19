//! `CalculateReadGroupChecksum`, mirroring Picard 3.5.0
//! `picard.sam.CalculateReadGroupChecksum` with htsjdk 5.0.0
//! `SAMUtils.calculateReadGroupRecordChecksum`.
//!
//! The checksum is the MD5 of the header's read groups, written as 32
//! lower-case hex digits with no line end to `OUTPUT`, or without it to
//! `INPUT` plus `.read_group_md5` beside the input. htsjdk computes it so:
//!
//! 1. The `@RG` records as `SAMTextHeaderCodec` keeps them, in header
//!    order, a repeated `ID` kept twice. Each record holds its fields but
//!    `ID` (the record's identifier, held apart), each key once with its
//!    first value; a field without a colon or a value, a key that is not
//!    two characters, and a line without `ID` are dropped. A `DT` value is
//!    read as a date and held as `Iso8601Date.toString` writes it, in the
//!    JVM's time zone ([`crate::java_date`]); one that is no date stays as
//!    written.
//! 2. A stable sort by `PU`, `LB`, `DT`, `SM`, `CN`, `PL`, `DS` (and `ID`,
//!    which is never among the fields), each compared as Java compares
//!    strings, a record without the field first.
//! 3. For each record, its fields in key order (a `TreeMap`): the key,
//!    then the value, each as `String.getBytes("UTF-8")`, into one MD5.
//!    A header without `@RG` gives the MD5 of nothing.
//! 4. `BigInteger(1, digest).toString(16)`, zeros put back in front.
//!
//! The text is the one htsjdk reads: a BAM's header text one character per
//! byte (`é` stored as `C3 A9` is two characters, four bytes in UTF-8),
//! SAM text and a CRAM's header as UTF-8. Every header line is checked at
//! `VALIDATION_STRINGENCY` as `SAMTextHeaderCodec` checks it
//! ([`crate::sam_text::Reader`]): at `STRICT` an `@RG` without `SM`, two
//! values of one key or a `DT` that is no date fails the run.
//!
//! Where bammetrics differs (docs/07 decision rows): a `DT` written in the
//! JVM locale's medium date format (`Dec 31, 1969, 5:00:00 p.m.`), which
//! htsjdk reads with `DateFormat` before it tries `Date.parse`, is no date
//! here; an `OUTPUT` that is the input is refused; a failed run writes
//! nothing.
//!
//! The Python function returns one row: the checksum.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};

use crate::bam_raw::read_append;
use crate::cleanup::AtomicFile;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::gzip::GzipReader;
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::sam_text::Stringency;

pub const TOOL: Tool = Tool {
    name: "CalculateReadGroupChecksum",
    summary: "writes the MD5 of the header's read groups as Picard computes it (32 hex \
              digits, no line end) to OUTPUT, or beside INPUT as INPUT.read_group_md5. \
              Returns one row: the checksum.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The SAM, BAM or CRAM whose @RG lines are read.",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The file to write (default: beside INPUT, its name plus .read_group_md5).",
        )
        .written(),
        Opt::text(
            &["VALIDATION_STRINGENCY"],
            "How the header lines are checked: STRICT (default), LENIENT or SILENT.",
        ),
        Opt::ignored(
            &["REFERENCE_SEQUENCE", "R"],
            "the header is read, and a CRAM's header needs no reference",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The column of the one row the tool returns.
pub const COLUMNS: [&str; 1] = ["checksum"];
/// What Picard appends to the input name without `OUTPUT`.
pub const EXTENSION: &str = ".read_group_md5";
/// The header text read at most, as `doctor` reads it.
pub const MAX_HEADER_BYTES: usize = crate::doctor::input::MAX_HEADER_BYTES;

/// The checksum file beside `input`: its path plus `.read_group_md5`.
pub fn default_output(input: &Path) -> PathBuf {
    let mut name = input.as_os_str().to_owned();
    name.push(EXTENSION);
    PathBuf::from(name)
}

/// `OUTPUT`, or the default name beside `INPUT`. A device or a pipe has no
/// name to put the file beside.
fn output_path(inv: &Invocation) -> Result<PathBuf> {
    if let Some(output) = inv.get("OUTPUT") {
        return Ok(PathBuf::from(output));
    }
    let regular = std::fs::metadata(&inv.input).is_ok_and(|m| m.is_file());
    if inv.input.starts_with("/dev") || (inv.input.exists() && !regular) {
        return Err(Error::Usage(format!(
            "INPUT={}: a pipe or a device has no name to put the checksum beside; give OUTPUT",
            inv.input.display()
        )));
    }
    Ok(default_output(&inv.input))
}

fn too_long() -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("the header text is over the {MAX_HEADER_BYTES}-byte cap"),
    )
}

/// The header text of `path` as bytes, and whether htsjdk reads it as
/// UTF-8 (SAM text, a CRAM) or one character per byte (a BAM).
fn header_bytes(path: &Path) -> io::Result<(Vec<u8>, bool)> {
    let sniffed = crate::cram::sniff(File::open(path)?)?;
    match sniffed.kind {
        Kind::Bam => {
            let mut reader = noodles::bgzf::io::Reader::new(sniffed.reader);
            let mut word = [0u8; 4];
            reader.read_exact(&mut word)?;
            if word != *b"BAM\x01" {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Invalid BAM file header",
                ));
            }
            reader.read_exact(&mut word)?;
            let length = usize::try_from(i32::from_le_bytes(word)).map_err(|_| too_long())?;
            if length > MAX_HEADER_BYTES {
                return Err(too_long());
            }
            let mut text = Vec::new();
            read_append(&mut reader, &mut text, length)?;
            Ok((text, false))
        }
        Kind::Cram => {
            let mut reader = noodles::cram::io::Reader::new(BufReader::new(sniffed.reader));
            let mut header = reader.header_reader();
            header.read_magic_number()?;
            header.read_format_version()?;
            header.read_file_id()?;
            let mut container = header.container_reader()?;
            let mut raw = container.raw_sam_header_reader()?;
            let mut text = Vec::new();
            (&mut raw)
                .take(MAX_HEADER_BYTES as u64 + 1)
                .read_to_end(&mut text)?;
            if text.len() > MAX_HEADER_BYTES {
                return Err(too_long());
            }
            Ok((text, true))
        }
        Kind::Sam | Kind::GzipSam => {
            let mut reader: Box<dyn BufRead> = if sniffed.kind == Kind::GzipSam {
                Box::new(BufReader::new(GzipReader::new(BufReader::new(
                    sniffed.reader,
                ))))
            } else {
                Box::new(BufReader::new(sniffed.reader))
            };
            let mut text = Vec::new();
            while reader.fill_buf()?.first() == Some(&b'@') {
                read_line(&mut reader, &mut text)?;
                if text.len() > MAX_HEADER_BYTES {
                    return Err(too_long());
                }
            }
            Ok((text, true))
        }
    }
}

/// One line of `reader` onto `text` with its end (`\n`, `\r` or `\r\n`),
/// as `BufferedReader.readLine` ends it; at most one more step than the
/// cap is held before the caller checks it.
fn read_line(reader: &mut dyn BufRead, text: &mut Vec<u8>) -> io::Result<()> {
    loop {
        let buffer = reader.fill_buf()?;
        if buffer.is_empty() {
            return Ok(());
        }
        let Some(end) = buffer.iter().position(|&b| b == b'\n' || b == b'\r') else {
            text.extend_from_slice(buffer);
            let taken = buffer.len();
            reader.consume(taken);
            if text.len() > MAX_HEADER_BYTES {
                return Err(too_long());
            }
            continue;
        };
        let cr = buffer[end] == b'\r';
        text.extend_from_slice(&buffer[..=end]);
        reader.consume(end + 1);
        if cr && reader.fill_buf()?.first() == Some(&b'\n') {
            text.push(b'\n');
            reader.consume(1);
        }
        return Ok(());
    }
}

/// The header text as the Java string htsjdk holds: UTF-16 units.
fn java_string(bytes: &[u8], utf8: bool) -> Vec<u16> {
    if utf8 {
        String::from_utf8_lossy(bytes).encode_utf16().collect()
    } else {
        bytes.iter().map(|&b| u16::from(b)).collect()
    }
}

/// `String.getBytes("UTF-8")`: a lone surrogate becomes `?`.
fn utf8_bytes(text: &[u16]) -> Vec<u8> {
    char::decode_utf16(text.iter().copied())
        .map(|c| c.unwrap_or('?'))
        .collect::<String>()
        .into_bytes()
}

fn units(text: &str) -> Vec<u16> {
    text.encode_utf16().collect()
}

/// One `@RG` record as `SAMReadGroupRecord` holds it.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ReadGroup {
    pub id: Vec<u16>,
    /// The fields but `ID`, each key once, in line order.
    pub fields: Vec<(Vec<u16>, Vec<u16>)>,
}

impl ReadGroup {
    fn get(&self, key: &str) -> Option<&[u16]> {
        let key = units(key);
        self.fields
            .iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.as_slice())
    }
}

/// The header lines htsjdk reads: lines ended by `\n`, `\r` or `\r\n`, up
/// to the first that does not start with `@`.
fn header_lines(text: &[u16]) -> Vec<&[u16]> {
    let (cr, lf, at) = (u16::from(b'\r'), u16::from(b'\n'), u16::from(b'@'));
    let mut lines = Vec::new();
    let mut start = 0;
    let mut i = 0;
    while start < text.len() {
        if text[start] != at {
            break;
        }
        while i < text.len() && text[i] != cr && text[i] != lf {
            i += 1;
        }
        lines.push(&text[start..i]);
        if i < text.len() && text[i] == cr && text.get(i + 1) == Some(&lf) {
            i += 1;
        }
        i += 1;
        start = i;
    }
    lines
}

/// The read groups htsjdk keeps from `text`, in header order. A `DT` that
/// is no date is reported at `stringency` (an error at `STRICT`).
pub fn read_groups(text: &[u16], stringency: Stringency) -> Result<Vec<ReadGroup>> {
    let rg = units("@RG");
    let (tab, colon) = (u16::from(b'\t'), u16::from(b':'));
    let mut groups = Vec::new();
    for line in header_lines(text) {
        let mut fields: Vec<&[u16]> = line.split(|&u| u == tab).collect();
        if fields.first() != Some(&rg.as_slice()) {
            continue;
        }
        // `StringUtil.split` drops an empty last token.
        if fields.last().is_some_and(|f| f.is_empty()) {
            fields.pop();
        }
        let mut kept: Vec<(Vec<u16>, Vec<u16>)> = Vec::new();
        for field in &fields[1..] {
            let Some(at) = field.iter().position(|&u| u == colon) else {
                continue;
            };
            let (key, value) = (&field[..at], &field[at + 1..]);
            if value.is_empty() {
                continue;
            }
            // A second value of a key is dropped (an error htsjdk reports),
            // the same value again kept once.
            if kept.iter().any(|(k, _)| k.as_slice() == key) {
                continue;
            }
            if key.len() != 2 {
                continue;
            }
            kept.push((key.to_vec(), value.to_vec()));
        }
        let id_key = units("ID");
        let Some(at) = kept.iter().position(|(k, _)| *k == id_key) else {
            continue;
        };
        let (_, id) = kept.remove(at);
        let mut group = ReadGroup { id, fields: kept };
        if let Some(value) = group.get("DT").map(<[u16]>::to_vec) {
            let text = String::from_utf16_lossy(&value);
            let decoded = crate::java_date::decode_date(&text).map_err(Error::Metric)?;
            match decoded.and_then(crate::java_date::iso8601) {
                Some(date) => {
                    let date = units(&date);
                    let dt = units("DT");
                    if let Some(slot) = group.fields.iter_mut().find(|(k, _)| *k == dt) {
                        slot.1 = date;
                    }
                }
                None => report_date(&text, line, stringency)?,
            }
        }
        groups.push(group);
    }
    Ok(groups)
}

/// `reportErrorParsingLine` for a `DT` that is no date.
fn report_date(value: &str, line: &[u16], stringency: Stringency) -> Result<()> {
    let message = format!(
        "Error parsing SAM header. DT tag value '{value}' is not parseable as a date. Line:\n{}",
        String::from_utf16_lossy(line)
    );
    match stringency {
        Stringency::Strict => Err(Error::Metric(message)),
        Stringency::Lenient => {
            eprintln!("Ignoring SAM validation error: {message}");
            Ok(())
        }
        Stringency::Silent => Ok(()),
    }
}

/// `SAMHeaderRecordComparator` over `PU`, `LB`, `DT`, `SM`, `CN`, `PL`,
/// `DS` and `ID`: a missing field first, then Java's `String.compareTo`
/// (UTF-16 units).
fn compare(a: &ReadGroup, b: &ReadGroup) -> std::cmp::Ordering {
    for key in ["PU", "LB", "DT", "SM", "CN", "PL", "DS", "ID"] {
        let order = match (a.get(key), b.get(key)) {
            (None, None) => continue,
            (None, Some(_)) => std::cmp::Ordering::Less,
            (Some(_), None) => std::cmp::Ordering::Greater,
            (Some(x), Some(y)) => x.cmp(y),
        };
        if order.is_ne() {
            return order;
        }
    }
    std::cmp::Ordering::Equal
}

/// `calculateReadGroupRecordChecksum` of `groups`: 32 lower-case hex digits.
pub fn checksum(groups: &[ReadGroup]) -> String {
    let mut sorted: Vec<&ReadGroup> = groups.iter().collect();
    sorted.sort_by(|a, b| compare(a, b));
    let mut md5 = crate::md5::Md5::new();
    for group in sorted {
        let mut fields: Vec<&(Vec<u16>, Vec<u16>)> = group.fields.iter().collect();
        fields.sort_by(|a, b| a.0.cmp(&b.0));
        for (key, value) in fields {
            md5.update(&utf8_bytes(key));
            md5.update(&utf8_bytes(value));
        }
    }
    crate::md5::hex(&md5.finish())
}

/// The header of `path` checked at `stringency` as htsjdk checks it, then
/// its checksum. Two `@RG` lines may share an `ID` (`setReadGroups` keeps
/// both); two `@SQ` lines may not share a name (`SAMSequenceDictionary`
/// refuses it at any stringency).
pub fn checksum_of(path: &Path, stringency: Stringency) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let (bytes, utf8) = header_bytes(path).map_err(input)?;
    let reader =
        crate::sam_text::Reader::new(Cursor::new(&bytes), stringency, Some(path)).map_err(input)?;
    let mut names = std::collections::HashSet::new();
    for line in reader.header_text().split(|&b| b == b'\n') {
        let name = line.strip_prefix(b"@SQ\t").and_then(|fields| {
            fields
                .split(|&b| b == b'\t')
                .find_map(|f| f.strip_prefix(b"SN:"))
        });
        if let Some(name) = name.filter(|name| !names.insert(*name)) {
            return Err(input(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Cannot add sequence that already exists in SAMSequenceDictionary: {}",
                    String::from_utf8_lossy(name)
                ),
            )));
        }
    }
    let groups = read_groups(&java_string(&bytes, utf8), stringency)?;
    Ok(checksum(&groups))
}

/// Write `text` to `output` under a temporary name, renamed into place.
fn write_checksum(output: &Path, text: &str) -> Result<()> {
    let error = |source| Error::Output {
        path: output.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(output);
    let mut out = file.create().map_err(error)?;
    out.write_all(text.as_bytes()).map_err(error)?;
    drop(out);
    file.commit().map_err(error)
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let output = output_path(inv)?;
    registry::check_writable([output.as_path()], &[inv.input.clone()])?;
    if crate::mark_duplicates::same_file(&inv.input, &output) {
        return Err(Error::Usage(format!(
            "OUTPUT={}: it is the input {}; writing it would destroy it",
            output.display(),
            inv.input.display()
        )));
    }
    let text = checksum_of(&inv.input, inv.stringency()?)?;
    write_checksum(&output, &text)?;
    let rows = vec![vec![Value::Str(text)]];
    Ok(MetricsOutput::text(&COLUMNS, rows, String::new()).into())
}

#[cfg(test)]
mod tests;

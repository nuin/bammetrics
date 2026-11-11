//! `SAMTextWriter` and `SAMTextHeaderCodec.encode`: records and the header
//! as SAM text.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles::sam;

use super::invalid;
use crate::bam_raw::{RawRecord, TagValue};
use crate::bam_writer;
use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::java;

/// `SAMFileHeader.CURRENT_VERSION`.
const CURRENT_VERSION: &str = "1.6";
/// `SAMUtils.MAX_PHRED_SCORE`.
const MAX_PHRED_SCORE: u8 = 93;
const CIGAR_OPS: &[u8; 9] = b"MIDNSHP=X";

/// The header as `SAMTextHeaderCodec.encode` writes it: an `@HD` line
/// first (`VN` then the other fields; `VN:1.6` alone for a header without
/// one), a non-conforming `SO` as `unknown`, then the `@SQ`, `@RG`, `@PG`
/// and `@CO` lines. `keep_version` false writes `VN:1.6` whatever the
/// header says, as a `SAMFileWriter` does; true keeps the header's
/// version, as `ViewSam` does. An `@RG DT` value stays as it is, where
/// htsjdk rewrites it in the JVM's time zone (KNOWN_DIFFERENCES.md).
pub fn header_text(header: &sam::Header, keep_version: bool) -> Vec<u8> {
    let body = bam_writer::header_bytes(header);
    let (first, rest) = match body.iter().position(|&b| b == b'\n') {
        Some(end) if body.starts_with(b"@HD") => (&body[..end], &body[end + 1..]),
        _ => (b"".as_slice(), body.as_slice()),
    };
    let mut text = b"@HD".to_vec();
    let mut fields = first.split(|&b| b == b'\t').skip(1);
    let version = fields.next().filter(|_| keep_version);
    text.push(b'\t');
    text.extend_from_slice(version.unwrap_or(format!("VN:{CURRENT_VERSION}").as_bytes()));
    for field in fields {
        text.push(b'\t');
        match field.strip_prefix(b"SO:") {
            Some(order)
                if !matches!(
                    order,
                    b"unsorted" | b"queryname" | b"coordinate" | b"duplicate" | b"unknown"
                ) =>
            {
                text.extend_from_slice(b"SO:unknown");
            }
            _ => text.extend_from_slice(field),
        }
    }
    text.push(b'\n');
    text.extend_from_slice(rest);
    text
}

fn push_number(out: &mut Vec<u8>, value: impl std::fmt::Display) {
    // Writing to a `Vec` cannot fail.
    let _ = write!(out, "{value}");
}

/// The elements of a `B` array as `TextTagCodec.encodeArrayValue` prints
/// them: a comma before each.
fn push_array(out: &mut Vec<u8>, subtype: u8, bytes: &[u8]) -> io::Result<()> {
    macro_rules! elements {
        ($ty:ty, $width:literal) => {
            for chunk in bytes.as_chunks::<$width>().0 {
                out.push(b',');
                push_number(out, <$ty>::from_le_bytes(*chunk));
            }
        };
    }
    match subtype {
        b'c' => elements!(i8, 1),
        b'C' => elements!(u8, 1),
        b's' => elements!(i16, 2),
        b'S' => elements!(u16, 2),
        b'i' => elements!(i32, 4),
        b'I' => elements!(u32, 4),
        b'f' => {
            for chunk in bytes.as_chunks::<4>().0 {
                out.push(b',');
                out.extend_from_slice(java::float_to_string(f32::from_le_bytes(*chunk)).as_bytes());
            }
        }
        _ => return Err(invalid("unknown array subtype".to_owned())),
    }
    Ok(())
}

/// `TextTagCodec.encode`: `TAG:TYPE:VALUE` with every integer as `i`.
fn push_tag(out: &mut Vec<u8>, name: [u8; 2], value: &TagValue) -> io::Result<()> {
    out.extend_from_slice(&name);
    out.push(b':');
    match value {
        TagValue::Char(c) => out.extend_from_slice(&[b'A', b':', *c]),
        TagValue::Int(v) => {
            out.extend_from_slice(b"i:");
            push_number(out, v);
        }
        TagValue::Float(f) => {
            out.extend_from_slice(b"f:");
            out.extend_from_slice(java::float_to_string(*f).as_bytes());
        }
        TagValue::String(text) => {
            out.extend_from_slice(b"Z:");
            out.extend_from_slice(text);
        }
        TagValue::Hex(bytes) => {
            out.extend_from_slice(b"B:c");
            push_array(out, b'c', bytes)?;
        }
        TagValue::Array(subtype, _, bytes) => {
            out.extend_from_slice(b"B:");
            out.push(*subtype);
            push_array(out, *subtype, bytes)?;
        }
    }
    Ok(())
}

/// The name of reference `id` (`*` for -1).
fn reference_name(references: &[Vec<u8>], id: i32) -> io::Result<&[u8]> {
    if id == -1 {
        return Ok(b"*");
    }
    usize::try_from(id)
        .ok()
        .and_then(|index| references.get(index))
        .map(Vec::as_slice)
        .ok_or_else(|| invalid(format!("reference index {id} is not in the header")))
}

/// `SAMTextWriter.writeAlignment`: append the SAM line of `record` and its
/// `\n` to `out`. `references` holds the header's sequence names in order.
pub fn write_record(
    out: &mut Vec<u8>,
    record: &RawRecord,
    references: &[Vec<u8>],
) -> io::Result<()> {
    out.extend_from_slice(record.name());
    out.push(b'\t');
    push_number(out, record.flags());
    out.push(b'\t');
    out.extend_from_slice(reference_name(references, record.reference_id())?);
    out.push(b'\t');
    push_number(out, record.alignment_start());
    out.push(b'\t');
    push_number(out, record.mapping_quality());
    out.push(b'\t');
    let mut operators = 0;
    for (length, op) in record.cigar() {
        push_number(out, length);
        out.push(
            *CIGAR_OPS
                .get(usize::from(op))
                .ok_or_else(|| invalid(format!("unknown CIGAR operator code {op}")))?,
        );
        operators += 1;
    }
    if operators == 0 {
        out.push(b'*');
    }
    out.push(b'\t');
    let mate = record.mate_reference_id();
    if mate == record.reference_id() && mate != -1 {
        out.push(b'=');
    } else {
        out.extend_from_slice(reference_name(references, mate)?);
    }
    out.push(b'\t');
    push_number(out, record.mate_alignment_start());
    out.push(b'\t');
    push_number(out, record.template_length());
    out.push(b'\t');
    if record.sequence_length() == 0 {
        out.extend_from_slice(b"*\t*");
    } else {
        out.extend_from_slice(&record.bases());
        out.push(b'\t');
        let qualities = record.qualities();
        if qualities[0] == 0xFF {
            out.push(b'*');
        } else {
            for &quality in qualities {
                if quality > MAX_PHRED_SCORE {
                    return Err(invalid(format!("Cannot encode phred score: {quality}")));
                }
                out.push(quality + 33);
            }
        }
    }
    for (name, value) in record.tags()?.iter() {
        out.push(b'\t');
        push_tag(out, *name, value)?;
    }
    out.push(b'\n');
    Ok(())
}

/// The sequence names of `header` in order, for [`write_record`].
pub(crate) fn reference_names(header: &sam::Header) -> Vec<Vec<u8>> {
    header
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect()
}

/// A SAM text file being written: under a temporary name
/// ([`AtomicFile`]) until [`SamWriter::finish`] renames it into place. A
/// path under `/dev` (`/dev/stdout`) is written in place.
pub struct SamWriter {
    path: PathBuf,
    file: AtomicFile,
    out: BufWriter<File>,
    references: Vec<Vec<u8>>,
    line: Vec<u8>,
}

impl SamWriter {
    /// Start `path` and write `header_text`, the header lines as they go
    /// into the file ([`header_text`], or the text a BAM of the same
    /// records would hold). `header` names the references of the records.
    pub fn create(path: &Path, header: &sam::Header, header_text: &[u8]) -> Result<Self> {
        let error = |source| Error::Output {
            path: path.to_path_buf(),
            source,
        };
        let file = AtomicFile::new(path);
        let mut out = BufWriter::with_capacity(1 << 20, file.create().map_err(error)?);
        out.write_all(header_text).map_err(error)?;
        Ok(Self {
            path: path.to_path_buf(),
            file,
            out,
            references: reference_names(header),
            line: Vec::new(),
        })
    }

    fn error(&self, source: io::Error) -> Error {
        Error::Output {
            path: self.path.clone(),
            source,
        }
    }

    pub fn write(&mut self, record: &RawRecord) -> Result<()> {
        self.line.clear();
        write_record(&mut self.line, record, &self.references).map_err(Error::Record)?;
        self.out
            .write_all(&self.line)
            .map_err(|source| self.error(source))
    }

    /// Flush the last bytes. The file keeps its temporary name until
    /// [`AtomicFile::commit`].
    pub fn close(mut self) -> Result<AtomicFile> {
        self.out.flush().map_err(|source| self.error(source))?;
        Ok(self.file)
    }

    /// Flush and rename the file into place. Returns its path.
    pub fn finish(self) -> Result<PathBuf> {
        let path = self.path.clone();
        self.close()?.commit().map_err(|source| Error::Output {
            path: path.clone(),
            source,
        })?;
        Ok(path)
    }
}

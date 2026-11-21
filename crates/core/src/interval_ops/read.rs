//! Reading an [`IntervalList`]: a Picard interval list (`IntervalList.fromPath`
//! with `IntervalListCodec.decode`), a VCF (`VCFFileReader.toIntervalList`)
//! or a BED (`BedToIntervalList` before its sort).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use super::header::{Header, Sequence, SortOrder};
use super::{Interval, IntervalList};
use crate::bed_to_interval_list;
use crate::error::{Error, Result};
use crate::vcf::VcfReader;

/// The bytes of an interval list read at most, inflated when it is gzip:
/// the text is held whole while it is parsed.
pub const MAX_LIST_BYTES: u64 = 1 << 30;
/// The longest line of a VCF read for its intervals.
pub const MAX_LINE_BYTES: usize = crate::vcf::MAX_LINE_BYTES;

/// True for a name htsjdk's `IOUtil.openFileForReading` inflates:
/// `.gz`, `.gzip` or `.bfq`, case as written.
fn has_gzip_extension(path: &Path) -> bool {
    path.file_name().is_some_and(|name| {
        let name = name.to_string_lossy();
        [".gz", ".gzip", ".bfq"].iter().any(|e| name.ends_with(e))
    })
}

/// The text of `path`, at most [`MAX_LIST_BYTES`] after inflation. Bytes
/// that are not UTF-8 become U+FFFD, as Java's decoder makes them.
fn read_text(path: &Path, cap: u64) -> Result<String> {
    let input = |source| Error::Input {
        path: path.to_path_buf(),
        source,
    };
    let file = File::open(path).map_err(input)?;
    let mut bytes = Vec::new();
    if has_gzip_extension(path) {
        crate::gzip::GzipReader::new(BufReader::new(file))
            .take(cap + 1)
            .read_to_end(&mut bytes)
            .map_err(input)?;
    } else {
        file.take(cap + 1).read_to_end(&mut bytes).map_err(input)?;
    }
    if bytes.len() as u64 > cap {
        return Err(Error::Intervals {
            path: path.to_path_buf(),
            message: format!("the interval list holds more than {cap} bytes, the most read"),
        });
    }
    Ok(match String::from_utf8(bytes) {
        Ok(text) => text,
        Err(e) => String::from_utf8_lossy(e.as_bytes()).into_owned(),
    })
}

/// The lines of `text` as `BufferedReader.readLine` ends them: at `\n`,
/// `\r` or `\r\n`.
pub(crate) fn java_lines(text: &str) -> impl Iterator<Item = &str> {
    let mut rest = text;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let end = rest.find(['\n', '\r']).unwrap_or(rest.len());
        let line = &rest[..end];
        let skip = match rest.as_bytes().get(end) {
            Some(b'\r') if rest.as_bytes().get(end + 1) == Some(&b'\n') => 2,
            Some(_) => 1,
            None => 0,
        };
        rest = &rest[end + skip..];
        Some(line)
    })
}

/// `String.split("\t")`: trailing empty fields dropped.
fn java_split_tab(line: &str) -> Vec<&str> {
    let mut fields: Vec<&str> = line.split('\t').collect();
    while fields.last().is_some_and(|f| f.is_empty()) {
        fields.pop();
    }
    fields
}

fn parse_int(text: &str) -> std::result::Result<i32, String> {
    crate::java::parse_integer(text, 32).map(|v| v as i32)
}

impl IntervalList {
    /// `IntervalList.fromPath`: read through gzip when the name ends in
    /// `.gz`, `.gzip` or `.bfq`.
    pub fn read(path: &Path) -> Result<IntervalList> {
        let text = read_text(path, MAX_LIST_BYTES)?;
        IntervalList::parse(&text).map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message,
        })
    }

    /// `IntervalList.fromReader`: the leading `@` lines are the header (at
    /// least one), then every other line is a record, a later `@` line and
    /// a blank line skipped. A record on a contig the header lacks is
    /// skipped, as htsjdk does with a warning.
    pub fn parse(text: &str) -> std::result::Result<IntervalList, String> {
        let mut lines = java_lines(text).peekable();
        let mut header_lines = Vec::new();
        while let Some(line) = lines.next_if(|l| l.starts_with('@')) {
            header_lines.push(line);
        }
        if header_lines.is_empty() {
            return Err("Interval list file must contain header. ".to_owned());
        }
        let header = Header::decode(header_lines)?;
        let index = header.name_index();
        let mut intervals = Vec::new();
        for line in lines {
            if let Some(interval) = decode_record(line, &header, &index)? {
                intervals.push(interval);
            }
        }
        Ok(IntervalList { header, intervals })
    }

    /// `VCFFileReader.toIntervalList(path, includeFiltered)`: the
    /// `##contig` lines as the dictionary (`SN`, `LN` or 0, `AS` from
    /// `assembly`), then per record not filtered (or every record with
    /// `include_filtered`) `POS` to `INFO/END` (else the reference
    /// allele's end) on `+`, named by `ID`, or `interval-N` counting the
    /// records without one.
    pub fn from_vcf(path: &Path, include_filtered: bool) -> Result<IntervalList> {
        let mut reader = VcfReader::open_capped(path, MAX_LINE_BYTES, crate::vcf::MAX_DBSNP_LINES)?;
        let intervals_error = |message: String| Error::Intervals {
            path: path.to_path_buf(),
            message,
        };
        let contigs = reader.header().contigs();
        if contigs.is_empty() {
            return Err(intervals_error(
                "the VCF has no ##contig lines, so no sequence dictionary (Picard fails with a \
                 NullPointerException)"
                    .to_owned(),
            ));
        }
        let mut sequences = Vec::with_capacity(contigs.len());
        for (name, contig) in contigs {
            let length = match contig.length() {
                Some(length) => i32::try_from(length).map_err(|_| {
                    intervals_error(format!("For input string: \"{length}\" (##contig {name})"))
                })?,
                None => 0,
            };
            let mut sequence = Sequence::new(name, length);
            if let Some((_, assembly)) = contig
                .other_fields()
                .iter()
                .find(|(key, _)| key.as_ref() == "assembly")
            {
                sequence
                    .attributes
                    .push(("AS".to_owned(), assembly.clone()));
            }
            sequences.push(sequence);
        }
        let mut list = IntervalList::new(Header::new(sequences));
        let index = list.header.name_index();
        let mut unnamed = 0u64;
        while let Some(record) = reader.next_site()? {
            if record.filtered && !include_filtered {
                continue;
            }
            let Some(&contig) = index.get(record.contig.as_str()) else {
                return Err(intervals_error(format!(
                    "Cannot add interval {}:{}-{}\t+\t{}, contig not in header",
                    record.contig,
                    record.start,
                    record.end,
                    record.id.as_deref().unwrap_or(".")
                )));
            };
            let name = match record.id {
                Some(id) => id,
                None => {
                    unnamed += 1;
                    format!("interval-{unnamed}")
                }
            };
            list.intervals.push(Interval {
                contig,
                start: record.start,
                end: record.end,
                negative_strand: false,
                name: Some(name),
            });
        }
        Ok(list)
    }

    /// The intervals of a BED as `BedToIntervalList` reads them before its
    /// sort (`opts.sort` and `opts.unique` are not applied here), with the
    /// dictionary of `opts.sequence_dictionary`; the header says
    /// `SO:coordinate`, as `BedToIntervalList` writes it.
    pub fn from_bed(path: &Path, opts: &bed_to_interval_list::Options) -> Result<IntervalList> {
        let dictionary = bed_to_interval_list::read_dictionary(&opts.sequence_dictionary)?;
        let text = read_text(path, MAX_LIST_BYTES)?;
        let bed = bed_to_interval_list::read_bed(&text, &dictionary, opts)?;
        let mut sequences = Vec::with_capacity(dictionary.len());
        for entry in &dictionary {
            let length = i32::try_from(entry.length).map_err(|_| Error::Intervals {
                path: opts.sequence_dictionary.clone(),
                message: format!("For input string: \"{}\"", entry.length),
            })?;
            let mut sequence = Sequence::new(&entry.name, length);
            sequence.attributes = entry
                .other_fields
                .iter()
                .map(|f| match f.split_once(':') {
                    Some((k, v)) => (k.to_owned(), v.to_owned()),
                    None => (f.clone(), String::new()),
                })
                .collect();
            sequences.push(sequence);
        }
        let mut header = Header::new(sequences);
        header.set_sort_order(SortOrder::Coordinate);
        let mut list = IntervalList::new(header);
        // `read_bed` checked every coordinate against the dictionary, whose
        // lengths fit an `int`.
        list.intervals = bed
            .into_iter()
            .map(|i| Interval {
                contig: i.contig,
                start: i.start as i32,
                end: i.end as i32,
                negative_strand: i.negative_strand,
                name: i.name,
            })
            .collect();
        Ok(list)
    }
}

/// `IntervalListCodec.decode`: `None` for a skipped line.
fn decode_record(
    line: &str,
    header: &Header,
    index: &HashMap<String, usize>,
) -> std::result::Result<Option<Interval>, String> {
    if line.starts_with('@') || crate::java::trim(line).is_empty() {
        return Ok(None);
    }
    let fields = java_split_tab(line);
    let [contig, start, end, strand, name] = fields[..] else {
        return Err(format!(
            "Invalid interval record contains {} fields: {line}",
            fields.len()
        ));
    };
    let start = parse_int(start)?;
    let end = parse_int(end)?;
    if start < 1 {
        return Err(format!(
            "Coordinate less than 1: start value of {start} is less than 1 and thus illegal"
        ));
    }
    if start > end.wrapping_add(1) {
        return Err(format!(
            "Start value of {start} is greater than end + 1 for end of value: {end}. I'm afraid I cannot let you do that."
        ));
    }
    let negative_strand = match strand {
        "+" => false,
        "-" => true,
        _ => return Err(format!("Invalid strand field: {strand}")),
    };
    let Some(&index) = index.get(contig) else {
        return Ok(None);
    };
    let length = header.sequences[index].length;
    if length > 0 && length < end {
        return Err(format!(
            "interval with end: {end} extends beyond end of sequence with length: {length}"
        ));
    }
    Ok(Some(Interval {
        contig: index,
        start,
        end,
        negative_strand,
        name: Some(name.to_owned()),
    }))
}

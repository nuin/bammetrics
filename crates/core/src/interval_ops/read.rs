//! Reading an [`IntervalList`]: a Picard interval list (`IntervalList.fromPath`
//! with `IntervalListCodec.decode`), a VCF (`VCFFileReader.toIntervalList`)
//! or a BED (`BedToIntervalList` before its sort).

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use super::header::{Header, Sequence, SortOrder};
use super::{Interval, IntervalList};
use crate::bed_to_interval_list;
use crate::create_sequence_dictionary::{NAME_PATTERN, legal_name};
use crate::error::{Error, Result, excerpt};
use crate::vcf::SiteReader;

/// The bytes of an interval list read at most, inflated when it is gzip:
/// the text is held whole while it is parsed.
pub const MAX_LIST_BYTES: u64 = 1 << 30;
/// The longest line of a VCF read for its intervals.
pub const MAX_LINE_BYTES: usize = crate::vcf::MAX_LINE_BYTES;
/// The intervals one list or VCF holds at most. Each costs 76 bytes (name
/// `.`) to 102 bytes (`target_N`) of peak memory in `IntervalListToBed`
/// and about 127 bytes in an `IntervalListTools` operation, so a list at
/// the cap needs 1.5 to 2.5 GB; htsjdk has no bound.
pub const MAX_INTERVALS: usize = 20_000_000;
/// The header lines of one list (`@` lines) or VCF (`#` lines) at most:
/// about 450 bytes each in memory for an `@SQ` line.
pub const MAX_HEADER_LINES: usize = 1_000_000;
/// The bytes of an error message kept: a message quotes the line or
/// field it failed on, which may be the whole file.
const MESSAGE_BYTES: usize = 300;

/// True for a name htsjdk's `IOUtil.openFileForReading` inflates:
/// `.gz`, `.gzip` or `.bfq`, case as written.
fn has_gzip_extension(path: &Path) -> bool {
    path.file_name().is_some_and(|name| {
        let name = name.to_string_lossy();
        [".gz", ".gzip", ".bfq"].iter().any(|e| name.ends_with(e))
    })
}

/// The bytes of `path`, at most `cap` after inflation.
fn read_bytes(path: &Path, cap: u64) -> Result<Vec<u8>> {
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
    Ok(bytes)
}

/// The text of `path`, at most `cap` bytes after inflation. Bytes that are
/// not UTF-8 become U+FFFD, as Java's decoder makes them.
fn read_text(path: &Path, cap: u64) -> Result<String> {
    Ok(match String::from_utf8(read_bytes(path, cap)?) {
        Ok(text) => text,
        Err(e) => crate::java::decode_utf8(e.as_bytes()).into_owned(),
    })
}

/// The lines of `text` as `BufferedReader.readLine` ends them: at `\n`,
/// `\r` or `\r\n`. Each line is decoded alone, bytes that are not UTF-8
/// as U+FFFD: no line end is part of a UTF-8 sequence, so the lines are
/// those of the decoded text, and a text that is not UTF-8 is never copied
/// whole.
pub(crate) fn java_lines(text: &[u8]) -> impl Iterator<Item = Cow<'_, str>> {
    let mut rest = text;
    std::iter::from_fn(move || {
        if rest.is_empty() {
            return None;
        }
        let end = rest
            .iter()
            .position(|&b| b == b'\n' || b == b'\r')
            .unwrap_or(rest.len());
        let line = &rest[..end];
        let skip = match rest.get(end) {
            Some(b'\r') if rest.get(end + 1) == Some(&b'\n') => 2,
            Some(_) => 1,
            None => 0,
        };
        rest = &rest[end + skip..];
        Some(String::from_utf8_lossy(line))
    })
}

/// `message` cut to its first [`MESSAGE_BYTES`] bytes, with the length
/// it had.
pub(super) fn cut(mut message: String) -> String {
    if message.len() > MESSAGE_BYTES {
        let length = message.len();
        let mut end = MESSAGE_BYTES;
        while !message.is_char_boundary(end) {
            end -= 1;
        }
        message.truncate(end);
        message.push_str(&format!("... ({length} bytes)"));
    }
    message
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
    crate::java::parse_integer(text, 32)
        .map(|v| v as i32)
        .map_err(|message| excerpt(&message).into_owned())
}

impl IntervalList {
    /// `IntervalList.fromPath`: read through gzip when the name ends in
    /// `.gz`, `.gzip` or `.bfq`.
    pub fn read(path: &Path) -> Result<IntervalList> {
        let bytes = read_bytes(path, MAX_LIST_BYTES)?;
        IntervalList::parse_bytes(&bytes).map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message,
        })
    }

    /// `IntervalList.fromReader`: the leading `@` lines are the header (at
    /// least one), then every other line is a record, a later `@` line and
    /// a blank line skipped. A record on a contig the header lacks is
    /// skipped, as htsjdk does with a warning. At most
    /// [`MAX_HEADER_LINES`] header lines and [`MAX_INTERVALS`]
    /// intervals are read.
    pub fn parse(text: &str) -> std::result::Result<IntervalList, String> {
        Self::parse_bytes(text.as_bytes())
    }

    /// [`IntervalList::parse`] of bytes, each line decoded as UTF-8 alone.
    /// An error message is cut to its first 300 bytes.
    pub fn parse_bytes(text: &[u8]) -> std::result::Result<IntervalList, String> {
        Self::parse_limited(text, &Limits::DEFAULT, true).map_err(cut)
    }

    /// The intervals of `path` as tribble's `AbstractFeatureReader` reads
    /// them with `IntervalListCodec` (`IntervalListToBed`): the rules of
    /// [`IntervalList::read`], but a file without a header line is read
    /// with an empty dictionary, so every record is decoded, checked and
    /// then skipped.
    pub fn read_features(path: &Path) -> Result<IntervalList> {
        let bytes = read_bytes(path, MAX_LIST_BYTES)?;
        Self::parse_limited(&bytes, &Limits::DEFAULT, false).map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message: cut(message),
        })
    }

    pub(super) fn parse_limited(
        text: &[u8],
        limits: &Limits,
        require_header: bool,
    ) -> std::result::Result<IntervalList, String> {
        let mut lines = java_lines(text).peekable();
        let mut header_lines = Vec::new();
        while let Some(line) = lines.next_if(|l| l.starts_with('@')) {
            if header_lines.len() == limits.header_lines {
                return Err(format!(
                    "the header holds more than {} lines, the most read",
                    limits.header_lines
                ));
            }
            header_lines.push(line);
        }
        if require_header && header_lines.is_empty() {
            return Err("Interval list file must contain header. ".to_owned());
        }
        // `IntervalListCodec.readActualHeader` peeks the first character of
        // the line after the header, which an empty line does not have.
        if !require_header && lines.peek().is_some_and(|l| l.is_empty()) {
            return Err(
                "Unable to parse header with error: Index 0 out of bounds for length 0".to_owned(),
            );
        }
        let header = Header::decode_limited(
            header_lines.iter().map(|l| l.as_ref()),
            limits.header_fields,
        )?;
        drop(header_lines);
        let index = header.name_index();
        let mut intervals = Vec::new();
        for line in lines {
            if let Some(interval) = decode_record(&line, &header, &index)? {
                if intervals.len() == limits.intervals {
                    return Err(too_many_intervals(limits.intervals));
                }
                intervals.push(interval);
            }
        }
        Ok(IntervalList { header, intervals })
    }

    /// `VCFFileReader.toIntervalList(path, includeFiltered)`: the
    /// intervals of [`VcfIntervals`] added to a list under its dictionary;
    /// a record on a contig the dictionary lacks is an error
    /// (`IntervalList.add`). At most [`MAX_HEADER_LINES`] header lines and
    /// [`MAX_INTERVALS`] intervals are read.
    pub fn from_vcf(path: &Path, include_filtered: bool) -> Result<IntervalList> {
        let (header, mut sites) = VcfIntervals::open(path, include_filtered)?;
        let mut list = IntervalList::new(header);
        let index = list.header.name_index();
        while let Some(site) = sites.next_interval()? {
            let Some(&contig) = index.get(site.contig.as_str()) else {
                return Err(Error::Intervals {
                    path: path.to_path_buf(),
                    message: format!(
                        "Cannot add interval {}:{}-{}\t+\t{}, contig not in header",
                        site.contig, site.start, site.end, site.name
                    ),
                });
            };
            if list.intervals.len() == MAX_INTERVALS {
                return Err(Error::Intervals {
                    path: path.to_path_buf(),
                    message: too_many_intervals(MAX_INTERVALS),
                });
            }
            list.intervals.push(Interval {
                contig,
                start: site.start,
                end: site.end,
                negative_strand: false,
                name: Some(site.name),
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

fn too_many_intervals(limit: usize) -> String {
    format!("the list holds more than {limit} intervals, the most read")
}

/// The counts [`IntervalList::parse_bytes`] reads at most.
pub(super) struct Limits {
    pub header_lines: usize,
    pub header_fields: usize,
    pub intervals: usize,
}

impl Limits {
    pub(super) const DEFAULT: Limits = Limits {
        header_lines: MAX_HEADER_LINES,
        header_fields: super::header::MAX_HEADER_FIELDS,
        intervals: MAX_INTERVALS,
    };
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
            "Invalid interval record contains {} fields: {}",
            fields.len(),
            excerpt(line)
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

/// One interval of [`VcfIntervals`], its contig by name: `VcfToIntervalList`
/// writes a record on a contig the dictionary lacks as it is.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VcfInterval {
    pub contig: String,
    pub start: i32,
    pub end: i32,
    pub name: String,
}

/// `VCFFileReader.toIntervals(reader, includeFiltered)` over a plain or
/// gzip VCF: per record not filtered (or every record with
/// `include_filtered`), `POS` to `INFO/END` (else the reference allele's
/// end) on `+`, named by `ID`, or `interval-N` counting the records
/// without one (an `AtomicInteger`, so it wraps as a Java `int`). Each
/// record is checked as htsjdk's codec checks its site
/// ([`SiteReader::next_site`](crate::vcf::SiteReader::next_site)); lines are read one at a time under
/// [`MAX_LINE_BYTES`].
pub struct VcfIntervals {
    reader: SiteReader,
    include_filtered: bool,
    unnamed: i32,
}

impl VcfIntervals {
    /// Open `path` and make its dictionary as `VCFHeader.getSequenceDictionary`
    /// makes it from the `##contig` lines: `SN`, `LN` from `length`
    /// (`Integer.parseInt`) or 0 without one, the name checked as
    /// `SAMSequenceRecord` checks it, `AS` from `assembly`. A VCF without
    /// `##contig` lines has none, and Picard fails there with a
    /// `NullPointerException`.
    pub fn open(path: &Path, include_filtered: bool) -> Result<(Header, Self)> {
        let reader = SiteReader::open(
            path,
            MAX_LINE_BYTES,
            crate::vcf::MAX_DBSNP_LINES,
            MAX_HEADER_LINES as u64,
        )?;
        let intervals_error = |message: String| Error::Intervals {
            path: path.to_path_buf(),
            message,
        };
        let contigs = &reader.header().contigs;
        if contigs.is_empty() {
            return Err(intervals_error(
                "the VCF has no ##contig lines, so no sequence dictionary (Picard fails with a \
                 NullPointerException)"
                    .to_owned(),
            ));
        }
        let mut sequences = Vec::with_capacity(contigs.len());
        for contig in contigs {
            let length = match &contig.length {
                Some(text) => parse_int(text).map_err(|message| {
                    intervals_error(format!("{message} (##contig {})", contig.id))
                })?,
                None => 0,
            };
            if !legal_name(&contig.id) {
                return Err(intervals_error(format!(
                    "Sequence name '{}' doesn't match regex: '{NAME_PATTERN}'",
                    contig.id
                )));
            }
            let mut sequence = Sequence::new(&contig.id, length);
            if let Some(assembly) = &contig.assembly {
                sequence
                    .attributes
                    .push(("AS".to_owned(), assembly.clone()));
            }
            sequences.push(sequence);
        }
        let sites = Self {
            reader,
            include_filtered,
            unnamed: 0,
        };
        Ok((Header::new(sequences), sites))
    }

    /// The next interval, `None` at the end of the file.
    pub fn next_interval(&mut self) -> Result<Option<VcfInterval>> {
        while let Some(site) = self.reader.next_site()? {
            if site.filtered && !self.include_filtered {
                continue;
            }
            let name = match site.id {
                Some(id) => id,
                None => {
                    self.unnamed = self.unnamed.wrapping_add(1);
                    format!("interval-{}", self.unnamed)
                }
            };
            return Ok(Some(VcfInterval {
                contig: site.contig,
                start: site.start,
                end: site.end,
                name,
            }));
        }
        Ok(None)
    }
}

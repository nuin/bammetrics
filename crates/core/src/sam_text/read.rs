//! `SAMTextReader`: the header lines, then one record per line, and the
//! same records as a BAM record stream for [`crate::pass`].

use std::collections::HashMap;
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};

use noodles::sam;

use super::{LineParser, Stringency, invalid, java_bytes};
use crate::bam_raw::RawRecord;

/// The sort orders htsjdk's `SAMFileHeader.SortOrder` names.
const SORT_ORDERS: [&[u8]; 5] = [
    b"unsorted",
    b"queryname",
    b"coordinate",
    b"duplicate",
    b"unknown",
];
/// `SAMFileHeader.GroupOrder`.
const GROUP_ORDERS: [&[u8]; 3] = [b"none", b"query", b"reference"];

/// The bytes of one line held at most: a line is read whole into memory.
/// A Java `String`, which htsjdk reads a line into, holds 2^31 characters.
const MAX_LINE_BYTES: usize = 1 << 30;

/// Lines as `BufferedReader.readLine` ends them: at `\n`, `\r` or `\r\n`.
struct LineReader<R> {
    inner: R,
    /// Lines read so far: the number of the last one, from 1.
    number: u64,
    /// The last line ended at `\r`: a `\n` next belongs to that line end.
    after_cr: bool,
}

impl<R: BufRead> LineReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            number: 0,
            after_cr: false,
        }
    }

    /// The first byte of the next line; `None` at the end.
    fn peek(&mut self) -> io::Result<Option<u8>> {
        if self.after_cr {
            if self.inner.fill_buf()?.first() == Some(&b'\n') {
                self.inner.consume(1);
            }
            self.after_cr = false;
        }
        Ok(self.inner.fill_buf()?.first().copied())
    }

    /// Read the next line into `line` without its end. False at the end.
    fn read_line(&mut self, line: &mut Vec<u8>) -> io::Result<bool> {
        line.clear();
        if self.peek()?.is_none() {
            return Ok(false);
        }
        loop {
            let buffer = self.inner.fill_buf()?;
            if buffer.is_empty() {
                break;
            }
            match buffer.iter().position(|&b| b == b'\n' || b == b'\r') {
                Some(end) => {
                    line.extend_from_slice(&buffer[..end]);
                    self.after_cr = buffer[end] == b'\r';
                    self.inner.consume(end + 1);
                    break;
                }
                None => {
                    line.extend_from_slice(buffer);
                    let taken = buffer.len();
                    self.inner.consume(taken);
                }
            }
            if line.len() > MAX_LINE_BYTES {
                return Err(invalid(format!(
                    "SAM text line {} is over {MAX_LINE_BYTES} bytes",
                    self.number + 1
                )));
            }
        }
        self.number += 1;
        Ok(true)
    }
}

/// One header line under `SAMTextHeaderCodec`'s checks.
struct HeaderLine<'a> {
    stringency: Stringency,
    file: Option<&'a Path>,
    text: &'a [u8],
    number: u64,
}

impl HeaderLine<'_> {
    /// `reportErrorParsingLine`: an error under `STRICT`, a warning under
    /// `LENIENT`, nothing under `SILENT`.
    fn report(&self, reason: &str) -> io::Result<()> {
        let reason = format!(
            "Error parsing SAM header. {reason}. Line:\n{}",
            String::from_utf8_lossy(self.text)
        );
        match self.stringency {
            Stringency::Strict => {
                let file = self
                    .file
                    .map_or_else(String::new, |f| format!("File {}", f.display()));
                Err(invalid(format!(
                    "{reason}; {file}; Line number {}",
                    self.number
                )))
            }
            Stringency::Lenient => {
                eprintln!("Ignoring SAM validation error: {reason}");
                Ok(())
            }
            Stringency::Silent => Ok(()),
        }
    }

    /// The line as htsjdk keeps it: `None` for a line it drops, else the
    /// line without the fields it drops and with a non-conforming `SO`
    /// read as `unknown`.
    fn clean(&self) -> io::Result<Option<Vec<u8>>> {
        let mut fields = self.text.split(|&b| b == b'\t');
        let kind = fields.next().unwrap_or_default();
        let required: &[&[u8]] = match kind {
            b"@CO" => return Ok(Some(self.text.to_vec())),
            b"@HD" => &[b"VN"],
            b"@SQ" => &[b"SN", b"LN"],
            b"@RG" | b"@PG" => &[b"ID"],
            _ => {
                self.report("Unrecognized header record type")?;
                return Ok(None);
            }
        };
        let kind = String::from_utf8_lossy(kind);
        let mut fields: Vec<&[u8]> = fields.collect();
        // `StringUtil.split` drops an empty last token.
        if fields.last().is_some_and(|field| field.is_empty()) {
            fields.pop();
        }
        let mut kept: Vec<(&[u8], &[u8])> = Vec::new();
        let mut index: HashMap<&[u8], usize> = HashMap::new();
        for field in fields {
            let colon = field.iter().position(|&b| b == b':');
            let Some((key, value)) = colon
                .map(|at| (&field[..at], &field[at + 1..]))
                .filter(|(_, value)| !value.is_empty())
            else {
                self.report(&format!("Problem parsing {kind} key:value pair"))?;
                continue;
            };
            let show = |text: &[u8]| String::from_utf8_lossy(text).into_owned();
            let earlier = index.get(key).map(|&at| kept[at].1);
            if let Some(earlier) = earlier.filter(|earlier| *earlier != value) {
                self.report(&format!(
                    "Problem parsing {kind} key:value pair {0}:{1} clashes with {0}:{2}",
                    show(key),
                    show(value),
                    show(earlier)
                ))?;
                continue;
            }
            if key.len() != 2 {
                self.report(&format!(
                    "Problem parsing {kind} key:value pair {}:<value> key is not two characters",
                    show(key)
                ))?;
                continue;
            }
            let value = if key == b"SO" && !SORT_ORDERS.contains(&value) {
                match self.stringency {
                    Stringency::Strict => {
                        return Err(invalid(format!(
                            "Found non-conforming header SO tag: {}, exiting because \
                             VALIDATION_STRINGENCY=STRICT",
                            show(value)
                        )));
                    }
                    Stringency::Lenient => eprintln!(
                        "Found non-conforming header SO tag: {}. Treating as 'unknown'.",
                        show(value)
                    ),
                    Stringency::Silent => {}
                }
                b"unknown".as_slice()
            } else {
                value
            };
            if earlier.is_none() {
                index.insert(key, kept.len());
                kept.push((key, value));
            }
        }
        let value_of = |tag: &[u8]| index.get(tag).map(|&at| kept[at].1);
        for tag in required {
            if value_of(tag).is_none() {
                self.report(&format!(
                    "{kind} line missing {} tag",
                    String::from_utf8_lossy(tag)
                ))?;
                return Ok(None);
            }
        }
        match kind.as_ref() {
            "@HD" => {
                if let Some(order) = value_of(b"GO").filter(|v| !GROUP_ORDERS.contains(v)) {
                    self.report(&format!(
                        "@HD line has non-conforming GO tag value: {}.",
                        String::from_utf8_lossy(order)
                    ))?;
                }
            }
            "@RG" => {
                if value_of(b"SM").is_none() {
                    self.report("@RG line missing SM tag")?;
                }
                if let Some(size) = value_of(b"PI")
                    && std::str::from_utf8(size)
                        .ok()
                        .and_then(|size| size.parse::<i32>().ok())
                        .is_none()
                {
                    self.report(&format!(
                        "PI is not numeric: {}",
                        String::from_utf8_lossy(size)
                    ))?;
                }
            }
            _ => {}
        }
        let mut line = kind.as_bytes().to_vec();
        for (key, value) in kept {
            line.push(b'\t');
            line.extend_from_slice(key);
            line.push(b':');
            line.extend_from_slice(value);
        }
        Ok(Some(line))
    }
}

/// A SAM text file: [`Reader::new`] reads the header, [`Reader::records`]
/// gives the records.
pub struct Reader<R> {
    lines: LineReader<R>,
    header: Vec<u8>,
    stringency: Stringency,
    file: Option<PathBuf>,
}

impl<R: BufRead> Reader<R> {
    /// Read the header of `inner`: every leading line that starts with
    /// `@`, checked as `SAMTextHeaderCodec.decode` checks it under
    /// `stringency`. `file` is the name error messages give.
    pub fn new(inner: R, stringency: Stringency, file: Option<&Path>) -> io::Result<Self> {
        let mut lines = LineReader::new(inner);
        let mut header = Vec::new();
        let mut line = Vec::new();
        while lines.peek()? == Some(b'@') {
            lines.read_line(&mut line)?;
            let text = java_bytes(&line);
            let checked = HeaderLine {
                stringency,
                file,
                text: &text,
                number: lines.number,
            };
            if let Some(clean) = checked.clean()? {
                header.extend_from_slice(&clean);
                header.push(b'\n');
            }
        }
        Ok(Self {
            lines,
            header,
            stringency,
            file: file.map(Path::to_path_buf),
        })
    }

    /// The header lines htsjdk keeps, each ended by `\n`.
    pub fn header_text(&self) -> &[u8] {
        &self.header
    }

    /// The header as noodles holds it. htsjdk takes the `@HD` line
    /// wherever it stands; noodles wants it first, so it is given first.
    pub fn header(&self) -> io::Result<sam::Header> {
        let mut parser = sam::header::Parser::default();
        let lines = || {
            self.header
                .split(|&b| b == b'\n')
                .filter(|line| !line.is_empty())
        };
        let is_hd = |line: &&[u8]| line.starts_with(b"@HD");
        for line in lines()
            .filter(is_hd)
            .chain(lines().filter(|line| !is_hd(line)))
        {
            parser.parse_partial(line).map_err(|e| {
                let mut causes = e.to_string();
                let mut source = std::error::Error::source(&e);
                while let Some(cause) = source {
                    causes.push_str(&format!(": {cause}"));
                    source = cause.source();
                }
                invalid(format!(
                    "the SAM header line {:?} cannot be read: {causes}",
                    String::from_utf8_lossy(line)
                ))
            })?;
        }
        Ok(parser.finish())
    }

    /// The records, parsed against `header` (the value of
    /// [`Reader::header`]).
    pub fn records(self, header: &sam::Header) -> Records<R> {
        Records {
            parser: LineParser::new(header, self.stringency, self.file.as_deref()),
            lines: self.lines,
            line: Vec::new(),
        }
    }
}

/// The record lines of a SAM text file, in file order.
pub struct Records<R> {
    lines: LineReader<R>,
    parser: LineParser,
    line: Vec<u8>,
}

impl<R: BufRead> Records<R> {
    /// The next record; `None` at the end of the file.
    pub fn next_record(&mut self) -> io::Result<Option<RawRecord>> {
        if !self.lines.read_line(&mut self.line)? {
            return Ok(None);
        }
        self.parser.parse(&self.line, self.lines.number).map(Some)
    }
}

/// Bytes of BAM records a [`RecordStream`] holds before it parses more.
const STREAM_BUFFER: usize = 1 << 16;

/// The records of a SAM text file as the stream a BAM holds after its
/// header: `block_size` and the record, one after another. The loop that
/// reads a BAM then reads a SAM.
pub struct RecordStream<R> {
    records: Records<R>,
    buffer: Vec<u8>,
    at: usize,
    done: bool,
}

impl<R: BufRead> RecordStream<R> {
    pub fn new(records: Records<R>) -> Self {
        Self {
            records,
            buffer: Vec::with_capacity(STREAM_BUFFER),
            at: 0,
            done: false,
        }
    }

    fn refill(&mut self) -> io::Result<()> {
        self.buffer.clear();
        self.at = 0;
        while !self.done && self.buffer.len() < STREAM_BUFFER {
            match self.records.next_record()? {
                Some(record) => record.write_to(&mut self.buffer)?,
                None => self.done = true,
            }
        }
        Ok(())
    }
}

impl<R: BufRead> Read for RecordStream<R> {
    fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
        if self.at == self.buffer.len() {
            self.refill()?;
        }
        let n = out.len().min(self.buffer.len() - self.at);
        out[..n].copy_from_slice(&self.buffer[self.at..self.at + n]);
        self.at += n;
        Ok(n)
    }
}

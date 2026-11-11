//! `SAMLineParser`, `TextTagCodec.decode`, `TextCigarCodec.decode` and
//! `SAMRecord.isValid`: one SAM text line into a [`RawRecord`].

use std::borrow::Cow;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::Path;
use std::sync::Arc;

use noodles::sam;

use super::{Stringency, invalid, java_bytes};
use crate::bam_raw::{Fields, RawRecord, TagValue, Tags};
use crate::java;

/// htsjdk's `mFields` array: a line is cut into at most this many fields.
const MAX_FIELDS: usize = 10_000;
const REQUIRED_FIELDS: usize = 11;

const PAIRED: i32 = 0x1;
const PROPER_PAIR: i32 = 0x2;
const UNMAPPED: i32 = 0x4;
const MATE_UNMAPPED: i32 = 0x8;
const MATE_REVERSE: i32 = 0x20;
const FIRST: i32 = 0x40;
const SECOND: i32 = 0x80;
const SECONDARY: i32 = 0x100;
const SUPPLEMENTARY: i32 = 0x800;

/// `TextCigarCodec.decode`: `(length, op)` pairs with op the BAM code
/// (`MIDNSHP=X`), empty for `*`. The error is htsjdk's message.
pub fn parse_cigar(text: &[u8]) -> Result<Vec<(u32, u8)>, String> {
    let mut cigar = Vec::new();
    parse_cigar_into(text, &mut cigar)?;
    Ok(cigar)
}

fn parse_cigar_into(text: &[u8], cigar: &mut Vec<(u32, u8)>) -> Result<(), String> {
    cigar.clear();
    if text == b"*" {
        return Ok(());
    }
    let malformed = || format!("Malformed CIGAR string: {}", String::from_utf8_lossy(text));
    let mut at = 0;
    while at < text.len() {
        if !text[at].is_ascii_digit() {
            return Err(malformed());
        }
        let mut length: u32 = 0;
        while at < text.len() && text[at].is_ascii_digit() {
            length = length
                .wrapping_mul(10)
                .wrapping_add(u32::from(text[at] - b'0'));
            at += 1;
        }
        let op = match text.get(at) {
            Some(b'M') => 0,
            Some(b'I') => 1,
            Some(b'D') => 2,
            Some(b'N') => 3,
            Some(b'S') => 4,
            Some(b'H') => 5,
            Some(b'P') => 6,
            Some(b'=') => 7,
            Some(b'X') => 8,
            Some(other) => return Err(format!("Unrecognized CigarOperator: {other}")),
            None => return Err(malformed()),
        };
        cigar.push((length, op));
        at += 1;
    }
    Ok(())
}

/// The read bases a CIGAR covers: `M`, `I`, `S`, `=` and `X`.
fn read_length(cigar: &[(u32, u8)]) -> i64 {
    cigar
        .iter()
        .filter(|(_, op)| matches!(op, 0 | 1 | 4 | 7 | 8))
        .map(|(len, _)| i64::from(*len))
        .sum()
}

/// The `ADJACENT_INDEL_IN_CIGAR` errors reported for one CIGAR at most.
/// htsjdk reports one per pair of operators, a number that grows with the
/// square of the CIGAR's length.
const MAX_ADJACENT_INDEL_ERRORS: usize = 100;

/// `Cigar.isValid`: each error as its `SAMValidationError` type (with the
/// severity) and message.
fn cigar_errors(cigar: &[(u32, u8)]) -> Vec<(&'static str, &'static str)> {
    const INVALID: &str = "ERROR::INVALID_CIGAR";
    let real = |op: u8| matches!(op, 0 | 1 | 2 | 3 | 7 | 8);
    // htsjdk reports an `I` or `D` once for each later operator of the
    // same kind before the next `M`, `N`, `=`, `X` or `P`. Counted from the
    // right, so a long CIGAR is walked once.
    let mut later_same = vec![0u32; cigar.len()];
    let (mut insertions, mut deletions) = (0u32, 0u32);
    for (i, &(_, op)) in cigar.iter().enumerate().rev() {
        match op {
            1 => {
                later_same[i] = insertions;
                insertions += 1;
            }
            2 => {
                later_same[i] = deletions;
                deletions += 1;
            }
            0 | 3 | 6 | 7 | 8 => (insertions, deletions) = (0, 0),
            _ => {}
        }
    }
    let mut adjacent = 0;
    let mut errors = Vec::new();
    let mut seen_real = false;
    let last = cigar.len().wrapping_sub(1);
    for (i, &(length, op)) in cigar.iter().enumerate() {
        if length == 0 {
            errors.push((INVALID, "CIGAR element with zero length"));
        }
        match op {
            5 if i != 0 && i != last => {
                errors.push((
                    INVALID,
                    "Hard clipping operator not at start or end of CIGAR",
                ));
            }
            4 if i != 0 && i != last => {
                let inside =
                    "Soft clipping CIGAR operator can only be inside of hard clipping operator";
                if i == 1 {
                    let both_ends = cigar.len() == 3 && cigar[2].1 == 5;
                    if !both_ends && cigar[0].1 != 5 {
                        errors.push((INVALID, inside));
                    }
                } else if i == cigar.len() - 2 {
                    if cigar[last].1 != 5 {
                        errors.push((INVALID, inside));
                    }
                } else {
                    errors.push((
                        INVALID,
                        "Soft clipping CIGAR operator can at start or end of read, or be inside of hard clipping operator",
                    ));
                }
            }
            _ if real(op) => {
                seen_real = true;
                // Zero for an operator that is no `I` and no `D`.
                let message = if op == 1 {
                    "No M or N operator between pair of I operators in CIGAR"
                } else {
                    "No M or N operator between pair of D operators in CIGAR"
                };
                let room = MAX_ADJACENT_INDEL_ERRORS.saturating_sub(adjacent);
                let more = (later_same[i] as usize).min(room);
                adjacent += more;
                errors.extend(std::iter::repeat_n(
                    ("WARNING::ADJACENT_INDEL_IN_CIGAR", message),
                    more,
                ));
            }
            6 => {
                if i == 0 {
                    errors.push((INVALID, "Padding operator not valid at start of CIGAR"));
                } else if i == last {
                    errors.push((INVALID, "Padding operator not valid at end of CIGAR"));
                } else if !real(cigar[i - 1].1) || !real(cigar[i + 1].1) {
                    errors.push((
                        INVALID,
                        "Padding operator not between real operators in CIGAR",
                    ));
                }
            }
            _ => {}
        }
    }
    if !seen_real {
        errors.push((INVALID, "No real operator (M|I|D|N) in CIGAR"));
    }
    errors
}

/// A reference field as `SAMRecord` holds it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Reference {
    /// `*`.
    None,
    /// A name of the sequence dictionary: its index.
    Known(i32),
    /// A name the dictionary lacks.
    Unknown,
}

impl Reference {
    /// `getReferenceIndex`.
    fn index(self) -> i32 {
        match self {
            Reference::Known(index) => index,
            _ => -1,
        }
    }

    /// `hasReferenceName`.
    fn named(self) -> bool {
        self != Reference::None
    }
}

/// Why a tag field was not read.
enum TagError {
    /// A `SAMFormatException`: reported by the stringency.
    Format(String),
    /// Any other exception: fatal at every stringency.
    Fatal(String),
}

/// `Long.parseLong`.
fn parse_long(text: &[u8]) -> Option<i64> {
    std::str::from_utf8(text).ok()?.parse().ok()
}

/// `Integer.parseInt`.
fn parse_int(text: &[u8]) -> Option<i32> {
    std::str::from_utf8(text).ok()?.parse().ok()
}

/// `Float.parseFloat`.
fn parse_float(text: &[u8]) -> Option<f32> {
    java::parse_float(std::str::from_utf8(text).ok()?).ok()
}

/// `StringUtil.splitConcatenateExcessTokens(text, [_; 2], delimiter)`: the
/// text before the first delimiter and the rest, `None` when the rest is
/// empty or there is no delimiter.
fn split_two(text: &[u8], delimiter: u8) -> (&[u8], Option<&[u8]>) {
    match text.iter().position(|&b| b == delimiter) {
        Some(at) => {
            let rest = &text[at + 1..];
            (&text[..at], (!rest.is_empty()).then_some(rest))
        }
        None => (text, None),
    }
}

/// `TextTagCodec.covertStringArrayToObject`.
fn parse_array(text: &[u8]) -> Result<TagValue, String> {
    let (kind, values) = split_two(text, b',');
    let &[subtype] = kind else {
        return Err(format!(
            "Unrecognized element type for array tag value: {}",
            String::from_utf8_lossy(kind)
        ));
    };
    // Java's `String.split(",")` drops trailing empty strings.
    let mut items: Vec<&[u8]> = values.map_or_else(Vec::new, |v| v.split(|&b| b == b',').collect());
    while items.last().is_some_and(|item| item.is_empty()) {
        items.pop();
    }
    let count = u32::try_from(items.len()).map_err(|_| "Array tag too long".to_owned())?;
    let mut bytes = Vec::new();
    if subtype == b'f' {
        for item in items {
            let value = parse_float(item).ok_or_else(|| {
                "Array tag of type f should have single-precision floating point value".to_owned()
            })?;
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        return Ok(TagValue::Array(b'f', count, bytes));
    }
    let (min, max, width) = match subtype {
        b'c' => (i64::from(i8::MIN), i64::from(i8::MAX), 1),
        b'C' => (0, i64::from(u8::MAX), 1),
        b's' => (i64::from(i16::MIN), i64::from(i16::MAX), 2),
        b'S' => (0, i64::from(u16::MAX), 2),
        b'i' => (i64::from(i32::MIN), i64::from(i32::MAX), 4),
        b'I' => (0, i64::from(u32::MAX), 4),
        _ => {
            return Err(format!(
                "Unrecognized array tag element type: {}",
                char::from(subtype)
            ));
        }
    };
    for item in items {
        let value = parse_long(item).ok_or_else(|| {
            format!(
                "Array tag of type {} should have integral value",
                char::from(subtype)
            )
        })?;
        if value < min || value > max {
            return Err(format!(
                "Value for element of array tag of type {} is out of allowed range: {value}",
                char::from(subtype)
            ));
        }
        bytes.extend_from_slice(&value.to_le_bytes()[..width]);
    }
    Ok(TagValue::Array(subtype, count, bytes))
}

/// `StringUtil.hexStringToBytes`.
fn parse_hex(text: &[u8]) -> Option<Vec<u8>> {
    if !text.len().is_multiple_of(2) {
        return None;
    }
    let digit = |c: u8| char::from(c).to_digit(16).map(|d| d as u8);
    text.as_chunks::<2>()
        .0
        .iter()
        .map(|pair| Some(digit(pair[0])? << 4 | digit(pair[1])?))
        .collect()
}

/// `TextTagCodec.decode` and `SAMTag.makeBinaryTag`.
fn parse_tag(field: &[u8]) -> Result<([u8; 2], TagValue), TagError> {
    let format = |message: &str| TagError::Format(message.to_owned());
    let (key, rest) = split_two(field, b':');
    let Some(rest) = rest else {
        return Err(TagError::Format(format!(
            "Not enough fields in tag '{}'",
            String::from_utf8_lossy(field)
        )));
    };
    let (kind, value) = split_two(rest, b':');
    let value = value.unwrap_or_default();
    let decoded = match kind {
        b"Z" => TagValue::String(value.to_vec()),
        b"A" => match value {
            &[c] => TagValue::Char(c),
            _ => return Err(format("Tag of type A should have a single-character value")),
        },
        b"i" => {
            let number = parse_long(value)
                .ok_or_else(|| format("Tag of type i should have signed decimal value"))?;
            if number < i64::from(i32::MIN) || number > i64::from(u32::MAX) {
                return Err(TagError::Format(format!(
                    "Integer is out of range for both a 32-bit signed and unsigned integer: {}",
                    String::from_utf8_lossy(value)
                )));
            }
            TagValue::Int(number)
        }
        b"f" => TagValue::Float(parse_float(value).ok_or_else(|| {
            format("Tag of type f should have single-precision floating point value")
        })?),
        b"H" => TagValue::Hex(parse_hex(value).ok_or_else(|| {
            format("Tag of type H should have valid hex string with even number of digits")
        })?),
        b"B" => parse_array(value).map_err(TagError::Format)?,
        other => {
            return Err(TagError::Format(format!(
                "Unrecognized tag type: {}",
                String::from_utf8_lossy(other)
            )));
        }
    };
    let &[first, second] = key else {
        return Err(TagError::Fatal(format!(
            "String tag does not have length() == 2: {}",
            String::from_utf8_lossy(key)
        )));
    };
    Ok(([first, second], decoded))
}

/// `SamFlagField.parseDefault`: decimal, `0x` hexadecimal, `0` octal, or
/// the letter form. `Err(None)` is a `NumberFormatException`.
fn parse_flag(text: &[u8]) -> Result<i32, Option<String>> {
    let radix = |digits: &[u8], radix: u32| {
        std::str::from_utf8(digits)
            .ok()
            .and_then(|digits| i32::from_str_radix(digits, radix).ok())
            .ok_or(None)
    };
    match text {
        [] => Err(Some(
            "Could not determine flag field type; saw an empty flag field".to_owned(),
        )),
        [b'0', b'x' | b'X', digits @ ..] => radix(digits, 16),
        [b'0', _, ..] => radix(text, 8),
        [first, ..] if first.is_ascii_digit() => radix(text, 10),
        _ => {
            let mut flags = 0;
            for &c in text {
                flags |= match c {
                    b'p' => PAIRED,
                    b'P' => PROPER_PAIR,
                    b'u' => UNMAPPED,
                    b'U' => MATE_UNMAPPED,
                    b'r' => 0x10,
                    b'R' => MATE_REVERSE,
                    b'1' => FIRST,
                    b'2' => SECOND,
                    b's' => SECONDARY,
                    b'x' => 0x200,
                    b'd' => 0x400,
                    b'S' => SUPPLEMENTARY,
                    b'f' | b'F' | b'm' | b'M' | b'_' => 0,
                    other => {
                        return Err(Some(format!(
                            "Unknown flag character '{}' in flag '{}'",
                            char::from(other),
                            String::from_utf8_lossy(text)
                        )));
                    }
                };
            }
            Ok(flags)
        }
    }
}

/// `SAMSequenceRecord.truncateSequenceName`: the name up to its first
/// whitespace character.
fn truncate_name(name: &[u8]) -> &[u8] {
    let end = name
        .iter()
        .position(|b| matches!(b, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r'))
        .unwrap_or(name.len());
    &name[..end]
}

fn is_valid_base(base: u8) -> bool {
    matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C'
            | b'M'
            | b'G'
            | b'R'
            | b'S'
            | b'V'
            | b'T'
            | b'W'
            | b'Y'
            | b'H'
            | b'K'
            | b'D'
            | b'B'
            | b'N'
            | b'.'
            | b'='
    )
}

/// One line being parsed: the text and number for an error message.
struct Line<'a> {
    stringency: Stringency,
    file_message: &'a str,
    text: &'a [u8],
    number: u64,
}

impl Line<'_> {
    /// `makeErrorString`.
    fn message(&self, reason: &str) -> String {
        let number = if self.number == 0 {
            "unknown".to_owned()
        } else {
            self.number.to_string()
        };
        format!(
            "Error parsing text SAM file. {reason}; {}Line {number}\nLine: {}",
            self.file_message,
            String::from_utf8_lossy(self.text)
        )
    }

    /// `reportFatalErrorParsingLine`.
    fn fatal(&self, reason: &str) -> io::Error {
        invalid(self.message(reason))
    }

    /// `reportErrorParsingLine`: an error under `STRICT`, a warning under
    /// `LENIENT`, nothing under `SILENT`.
    fn report(&self, reason: &str) -> io::Result<()> {
        match self.stringency {
            Stringency::Strict => Err(self.fatal(reason)),
            Stringency::Lenient => {
                eprintln!("Ignoring SAM validation error due to lenient parsing:");
                eprintln!("{}", self.message(reason));
                Ok(())
            }
            Stringency::Silent => Ok(()),
        }
    }

    /// `parseInt`.
    fn int(&self, text: &[u8], column: &str) -> io::Result<i32> {
        parse_int(text).ok_or_else(|| self.fatal(&format!("Non-numeric value in {column} column")))
    }
}

/// The fields of a line as `SAMRecord` holds them before `isValid`.
struct Parsed<'a> {
    name: Cow<'a, [u8]>,
    flags: i32,
    reference: Reference,
    reference_name: Cow<'a, [u8]>,
    start: i32,
    mapq: i32,
    mate: Reference,
    mate_name: Cow<'a, [u8]>,
    mate_start: i32,
    template_length: i32,
    /// `None` for `*`.
    bases: Option<Cow<'a, [u8]>>,
    /// Phred scores; `None` for `*`.
    qualities: Option<Vec<u8>>,
    tags: Tags,
}

impl Parsed<'_> {
    fn has(&self, flag: i32) -> bool {
        self.flags & flag != 0
    }

    fn read_len(&self) -> usize {
        self.bases.as_ref().map_or(0, |bases| bases.len())
    }
}

/// `SAMLineParser` for one file: the sequence dictionary and read groups
/// of its header, and the stringency.
pub struct LineParser {
    stringency: Stringency,
    /// `File <path>; `, empty without a file.
    file_message: Arc<str>,
    references: HashMap<Vec<u8>, i32>,
    lengths: Vec<i64>,
    read_groups: HashSet<Vec<u8>>,
    /// The CIGAR of the line being parsed, once [`LineParser::decode_cigar`]
    /// has run for it.
    cigar: Vec<(u32, u8)>,
    cigar_decoded: bool,
}

impl LineParser {
    /// `file` is named in error messages, as htsjdk names the `File`.
    pub fn new(header: &sam::Header, stringency: Stringency, file: Option<&Path>) -> Self {
        let references = header
            .reference_sequences()
            .keys()
            .enumerate()
            .map(|(index, name)| (java_bytes(name).into_owned(), index as i32))
            .collect();
        let lengths = header
            .reference_sequences()
            .values()
            .map(|map| usize::from(map.length()) as i64)
            .collect();
        let read_groups = header
            .read_groups()
            .keys()
            .map(|id| java_bytes(id).into_owned())
            .collect();
        Self {
            stringency,
            file_message: file
                .map_or_else(String::new, |f| format!("File {}; ", f.display()))
                .into(),
            references,
            lengths,
            read_groups,
            cigar: Vec::new(),
            cigar_decoded: false,
        }
    }

    /// `TextCigarCodec.decode` of the line's CIGAR, once per line.
    fn decode_cigar(&mut self, text: &[u8]) -> io::Result<()> {
        if !self.cigar_decoded {
            parse_cigar_into(&java_bytes(text), &mut self.cigar).map_err(invalid)?;
            self.cigar_decoded = true;
        }
        Ok(())
    }

    /// `resolveIndexFromName`.
    fn resolve(&self, name: &[u8]) -> Reference {
        if name == b"*" {
            return Reference::None;
        }
        self.references
            .get(name)
            .map_or(Reference::Unknown, |&index| Reference::Known(index))
    }

    /// `validateReferenceName`.
    fn validate_reference(&self, line: &Line<'_>, name: &[u8], field: &str) -> io::Result<()> {
        if name == b"=" {
            if field == "MRNM" {
                return Ok(());
            }
            line.report(&format!("= is not a valid value for {field} field."))?;
        }
        if !self.references.is_empty() && !self.references.contains_key(name) {
            line.report(&format!(
                "{field} '{}' not found in any SQ record",
                String::from_utf8_lossy(name)
            ))?;
        }
        Ok(())
    }

    /// `parseLine`: the record of `line` (without its line end), or the
    /// error htsjdk raises for it. `line_number` is 1-based; 0 for unknown.
    pub fn parse(&mut self, line: &[u8], line_number: u64) -> io::Result<RawRecord> {
        let file_message = Arc::clone(&self.file_message);
        let at = Line {
            stringency: self.stringency,
            file_message: &file_message,
            text: line,
            number: line_number,
        };
        let mut fields: Vec<&[u8]> = line
            .splitn(MAX_FIELDS + 1, |&b| b == b'\t')
            .take(MAX_FIELDS)
            .collect();
        // `StringUtil.split` drops an empty last token.
        if fields.len() > 1 && fields.last().is_some_and(|field| field.is_empty()) {
            fields.pop();
        }
        if fields.len() < REQUIRED_FIELDS {
            return Err(at.fatal("Not enough fields"));
        }
        if fields.len() == MAX_FIELDS {
            at.report("Too many fields in SAM text record.")?;
        }
        for (i, field) in fields.iter().enumerate() {
            if field.is_empty() {
                at.report(&format!("Empty field at position {i} (zero-based)"))?;
            }
        }
        self.cigar_decoded = false;
        let parsed = self.fields(&at, &fields)?;
        if self.stringency != Stringency::Silent {
            self.validate(&at, &parsed, fields[5])?;
        }
        self.encode(&at, &parsed, fields[5])
    }

    /// The eleven mandatory fields and the tags, with the checks
    /// `parseLine` makes on the way.
    fn fields<'a>(&self, at: &Line<'_>, fields: &[&'a [u8]]) -> io::Result<Parsed<'a>> {
        let flags = parse_flag(fields[1]).map_err(|e| match e {
            None => at.fatal("Non-numeric value in FLAG column"),
            Some(message) => at.fatal(&format!("Error in FLAG column: {message}")),
        })?;
        let unmapped = flags & UNMAPPED != 0;
        let paired = flags & PAIRED != 0;

        let mut reference = Reference::None;
        let mut reference_name = Cow::Borrowed(b"*".as_slice());
        if fields[2] != b"*" {
            let name = match java_bytes(fields[2]) {
                Cow::Borrowed(name) => Cow::Borrowed(truncate_name(name)),
                Cow::Owned(name) => Cow::Owned(truncate_name(&name).to_vec()),
            };
            self.validate_reference(at, &name, "RNAME")?;
            reference = self.resolve(&name);
            reference_name = name;
        } else if !unmapped {
            at.report("RNAME is not specified but flags indicate mapped")?;
        }

        let start = at.int(fields[3], "POS")?;
        let mapq = at.int(fields[4], "MAPQ")?;
        let cigar = fields[5];
        if reference_name.as_ref() != b"*" {
            if start == 0 {
                at.report("POS must be non-zero if RNAME is specified")?;
            }
            if !unmapped && cigar == b"*" {
                at.report("CIGAR must not be '*' if RNAME is specified")?;
            }
        } else {
            if start != 0 {
                at.report("POS must be zero if RNAME is not specified")?;
            }
            if mapq != 0 {
                at.report("MAPQ must be zero if RNAME is not specified")?;
            }
            if cigar != b"*" {
                at.report("CIGAR must be '*' if RNAME is not specified")?;
            }
        }

        let mut mate = Reference::None;
        let mut mate_name = Cow::Borrowed(b"*".as_slice());
        if fields[6] == b"*" {
            if paired && flags & MATE_UNMAPPED == 0 {
                at.report("MRNM not specified but flags indicate mate mapped")?;
            }
        } else {
            if !paired {
                at.report("MRNM specified but flags indicate unpaired")?;
            }
            if fields[6] == b"=" {
                mate = reference;
                mate_name = reference_name.clone();
            } else {
                let name = match java_bytes(fields[6]) {
                    Cow::Borrowed(name) => Cow::Borrowed(truncate_name(name)),
                    Cow::Owned(name) => Cow::Owned(truncate_name(&name).to_vec()),
                };
                self.validate_reference(at, &name, "MRNM")?;
                mate = self.resolve(&name);
                mate_name = name;
            }
        }

        let mate_start = at.int(fields[7], "MPOS")?;
        let template_length = at.int(fields[8], "ISIZE")?;
        if mate_name.as_ref() != b"*" {
            if mate_start == 0 {
                at.report("MPOS must be non-zero if MRNM is specified")?;
            }
        } else {
            if mate_start != 0 {
                at.report("MPOS must be zero if MRNM is not specified")?;
            }
            if template_length != 0 {
                at.report("ISIZE must be zero if MRNM is not specified")?;
            }
        }

        let bases = if fields[9] == b"*" {
            None
        } else {
            let bases = java_bytes(fields[9]);
            if !bases.iter().copied().all(is_valid_base) {
                at.report("Invalid character in read bases")?;
            }
            Some(bases)
        };
        let qualities = if fields[10] == b"*" {
            None
        } else {
            if bases.is_none() {
                at.report("QUAL should not be specified if SEQ is not specified")?;
            }
            // `getReadString` of a record without bases is `*`: length 1.
            let read_string_len = bases.as_ref().map_or(1, |b| b.len());
            let quals = java_bytes(fields[10]);
            if read_string_len != quals.len() {
                at.report("length(QUAL) != length(SEQ)")?;
            }
            // `fastqToPhred` tests the character, so one outside ASCII
            // fails whatever its low byte is.
            if let Some(bad) = String::from_utf8_lossy(fields[10])
                .chars()
                .find(|c| !('!'..='~').contains(c))
            {
                return Err(invalid(format!("Invalid fastq character: {bad}")));
            }
            Some(quals.iter().map(|q| q - 33).collect())
        };

        let mut tags = Vec::with_capacity(fields.len() - REQUIRED_FIELDS);
        for field in &fields[REQUIRED_FIELDS..] {
            match parse_tag(&java_bytes(field)) {
                Ok(tag) => tags.push(tag),
                Err(TagError::Format(message)) => at.report(&message)?,
                Err(TagError::Fatal(message)) => return Err(invalid(message)),
            }
        }
        let tags = Tags::from_unsorted(tags);

        Ok(Parsed {
            name: java_bytes(fields[0]),
            flags,
            reference,
            reference_name,
            start,
            mapq,
            mate,
            mate_name,
            mate_start,
            template_length,
            bases,
            qualities,
            tags,
        })
    }

    /// `isValidReferenceIndexAndPosition`.
    fn reference_errors(
        &self,
        reference: Reference,
        name: &[u8],
        start: i32,
        mate: bool,
        errors: &mut Vec<String>,
    ) {
        let prefix = if mate { "Mate " } else { "" };
        if !reference.named() {
            if start != 0 {
                errors.push(format!(
                    "{prefix}Alignment start should be 0 because reference name = *."
                ));
            }
            return;
        }
        if start == 0 {
            errors.push(format!(
                "{prefix}Alignment start should != 0 because reference name != *."
            ));
        }
        if self.references.is_empty() {
            return;
        }
        match reference {
            Reference::Known(index) => {
                let length = self.lengths[index as usize];
                if i64::from(start) > length {
                    errors.push(format!(
                        "{prefix}Alignment start ({start}) must be <= reference sequence length \
                         ({length}) on reference {}",
                        String::from_utf8_lossy(name)
                    ));
                }
            }
            _ => errors.push(format!(
                "{prefix}Reference sequence not found in sequence dictionary."
            )),
        }
    }

    /// `SAMRecord.getCigar` of a mapped record under `STRICT` or
    /// `LENIENT`: decode, then `validateCigar`, whose errors do not go
    /// through the line parser's report.
    fn validate_cigar(&mut self, parsed: &Parsed<'_>, text: &[u8]) -> io::Result<()> {
        self.decode_cigar(text)?;
        let mut errors: Vec<(&str, &str)> = cigar_errors(&self.cigar);
        if self.cigar.is_empty() {
            errors.clear();
        }
        if let Reference::Known(index) = parsed.reference {
            let length = self.lengths[index as usize];
            let mut position = i64::from(parsed.start);
            for &(len, op) in &self.cigar {
                let len = i64::from(len);
                match op {
                    0 | 7 | 8 => {
                        if position + len - 1 > length {
                            errors.push((
                                "ERROR::CIGAR_MAPS_OFF_REFERENCE",
                                "Read CIGAR M operator maps off end of reference",
                            ));
                            break;
                        }
                        position += len;
                    }
                    2 | 3 => position += len,
                    _ => {}
                }
            }
        }
        for (kind, message) in errors {
            let text = format!(
                "{kind}:Read name {}, {message}",
                String::from_utf8_lossy(&parsed.name)
            );
            if self.stringency == Stringency::Strict {
                return Err(invalid(format!("SAM validation error: {text}")));
            }
            eprintln!("Ignoring SAM validation error: {text}");
        }
        Ok(())
    }

    /// `SAMRecord.isValid`, each error reported in order.
    fn validate(&mut self, at: &Line<'_>, parsed: &Parsed<'_>, cigar: &[u8]) -> io::Result<()> {
        let mut errors = Vec::new();
        let unmapped = parsed.has(UNMAPPED);
        if !parsed.has(PAIRED) {
            for (flag, what) in [
                (PROPER_PAIR, "Proper pair"),
                (MATE_UNMAPPED, "Mate unmapped"),
                (MATE_REVERSE, "Mate negative strand"),
                (FIRST, "First of pair"),
                (SECOND, "Second of pair"),
            ] {
                if parsed.has(flag) {
                    errors.push(format!("{what} flag should not be set for unpaired read."));
                }
            }
            if parsed.mate.index() != -1 {
                errors.push("MRNM should not be set for unpaired read.".to_owned());
            }
            if parsed.mate_name.as_ref() != b"*" {
                errors.push(format!(
                    "Unpaired read mate reference is {} not * for unpaired read",
                    String::from_utf8_lossy(&parsed.mate_name)
                ));
            }
        } else {
            self.reference_errors(
                parsed.mate,
                &parsed.mate_name,
                parsed.mate_start,
                true,
                &mut errors,
            );
            if !parsed.mate.named() && !parsed.has(MATE_UNMAPPED) {
                errors.push("Mapped mate should have mate reference name".to_owned());
            }
            if !parsed.has(FIRST) && !parsed.has(SECOND) {
                errors.push(
                    "Paired read should be marked as first of pair or second of pair.".to_owned(),
                );
            }
        }
        if parsed.template_length == i32::MIN {
            errors.push("Insert size out of range".to_owned());
        }
        if unmapped {
            if parsed.has(SECONDARY) {
                errors.push(
                    "Secondary alignment flag should not be set for unmapped read.".to_owned(),
                );
            }
            if parsed.has(SUPPLEMENTARY) {
                errors.push(
                    "Supplementary alignment flag should not be set for unmapped read.".to_owned(),
                );
            }
            if parsed.mapq != 0 {
                errors.push("MAPQ should be 0 for unmapped read.".to_owned());
            }
        } else {
            if parsed.mapq >= 256 {
                errors.push("MAPQ should be < 256.".to_owned());
            }
            self.validate_cigar(parsed, cigar)?;
            if self.cigar.is_empty() {
                errors.push("CIGAR should have > zero elements for mapped read.".to_owned());
            }
            if self.references.is_empty() {
                errors.push("Empty sequence dictionary.".to_owned());
            }
            if !parsed.reference.named() {
                errors.push("Mapped read should have valid reference name".to_owned());
            }
        }
        match parsed.tags.get(*b"RG") {
            Some(TagValue::String(id)) => {
                if !self.read_groups.contains(id) {
                    errors.push(format!(
                        "RG ID on SAMRecord not found in header: {}",
                        String::from_utf8_lossy(id)
                    ));
                }
            }
            Some(_) => {
                return Err(invalid(
                    "the RG tag is not a string (java.lang.ClassCastException)".to_owned(),
                ));
            }
            None => {}
        }
        self.reference_errors(
            parsed.reference,
            &parsed.reference_name,
            parsed.start,
            false,
            &mut errors,
        );
        let read_len = parsed.read_len();
        let quals_len = parsed.qualities.as_ref().map_or(0, Vec::len);
        if parsed.qualities.is_some() && read_len != quals_len {
            errors.push("Read length does not match quals length".to_owned());
        }
        if parsed.mate_name.as_ref() == b"*" && parsed.mate_start != 0 {
            errors.push(format!(
                "The unaligned mate start position is {}, should be 0",
                parsed.start
            ));
        }
        if read_len != 0 {
            self.decode_cigar(cigar)?;
            let covered = read_length(&self.cigar);
            if covered != 0 && covered != read_len as i64 {
                errors.push(format!(
                    "CIGAR covers {covered} bases but the sequence is {read_len} read bases "
                ));
            }
        }
        if quals_len != 0 && read_len != quals_len {
            errors.push(format!(
                "Read length is  {read_len} bases but have {quals_len} qualities "
            ));
        }
        errors.iter().try_for_each(|error| at.report(error))
    }

    /// The BAM bytes of the record.
    fn encode(
        &mut self,
        at: &Line<'_>,
        parsed: &Parsed<'_>,
        cigar: &[u8],
    ) -> io::Result<RawRecord> {
        self.decode_cigar(cigar)?;
        let bases: &[u8] = parsed.bases.as_deref().unwrap_or_default();
        if parsed
            .qualities
            .as_ref()
            .is_some_and(|quals| quals.len() != bases.len())
        {
            return Err(at.fatal("length(QUAL) != length(SEQ), which a BAM record cannot hold"));
        }
        RawRecord::from_fields(&Fields {
            name: &parsed.name,
            flags: parsed.flags as u16,
            reference_id: parsed.reference.index(),
            alignment_start: parsed.start,
            mapping_quality: parsed.mapq as u8,
            cigar: &self.cigar,
            mate_reference_id: parsed.mate.index(),
            mate_alignment_start: parsed.mate_start,
            template_length: parsed.template_length,
            bases,
            qualities: parsed.qualities.as_deref(),
            tags: &parsed.tags,
        })
        .map_err(|e| at.fatal(&e.to_string()))
    }
}

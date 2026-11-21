//! The SAM header of an interval list as htsjdk 5.0.0 holds it: decoded by
//! `SAMTextHeaderCodec.decode` at `SILENT` (the codec `IntervalList.fromReader`
//! makes) and encoded by `SAMTextHeaderCodec.encode`.
//!
//! Decoding: a line htsjdk drops (an unknown record type, a required tag
//! missing) is dropped by [`crate::sam_text::clean_header_line`]. An `@SQ`
//! name is cut at its first whitespace and must match `SAMSequenceRecord`'s
//! pattern; `LN` is read by `Integer.parseInt`. `@HD` attributes keep the
//! order they were first set in, behind `VN`. An `@RG` `DT` is rewritten
//! by `Iso8601Date` in the process's zone. A second `@SQ` with a name, or
//! a second `@RG` or `@PG` with an ID, is htsjdk's error.
//!
//! Encoding: `@HD` with `VN:1.6` first, then the `@SQ`, `@RG`, `@PG` and
//! `@CO` lines, each attribute as `KEY:VALUE` in its order.

use std::collections::{HashMap, HashSet};

use crate::create_sequence_dictionary::{NAME_PATTERN, legal_name};
use crate::sam_text::{Stringency, clean_header_line};

/// `SAMFileHeader.CURRENT_VERSION`: the `VN` every written `@HD` carries.
const VERSION: &str = "1.6";

/// `SAMFileHeader.SortOrder`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unsorted => "unsorted",
            Self::Queryname => "queryname",
            Self::Coordinate => "coordinate",
            Self::Duplicate => "duplicate",
            Self::Unknown => "unknown",
        }
    }

    fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "unsorted" => Self::Unsorted,
            "queryname" => Self::Queryname,
            "coordinate" => Self::Coordinate,
            "duplicate" => Self::Duplicate,
            "unknown" => Self::Unknown,
            _ => return None,
        })
    }
}

/// One `@SQ` line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    /// `LN`, a Java `int`; 0 is htsjdk's unknown length.
    pub length: i32,
    /// The other tags in their order.
    pub attributes: Vec<(String, String)>,
}

impl Sequence {
    pub fn new(name: &str, length: i32) -> Self {
        Self {
            name: name.to_owned(),
            length,
            attributes: Vec::new(),
        }
    }

    fn attribute(&self, key: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// `SAMSequenceRecord.isSameSequence` for records at the same index:
    /// lengths equal unless one is 0, `M5` values equal as hexadecimal
    /// numbers when both are set, names equal. An `AN` alternative name is
    /// not followed (docs/07, interval layer decisions).
    fn same_as(&self, other: &Sequence) -> Result<bool, String> {
        if self.length != 0 && other.length != 0 && self.length != other.length {
            return Ok(false);
        }
        if let (Some(a), Some(b)) = (self.attribute("M5"), other.attribute("M5")) {
            let (a, b) = (hex_value(a)?, hex_value(b)?);
            if a != b {
                return Ok(false);
            }
        }
        Ok(self.name == other.name)
    }

    /// The `@SQ` line: `SN`, `LN`, the other tags.
    fn line(&self) -> String {
        let mut line = format!("@SQ\tSN:{}\tLN:{}", self.name, self.length);
        push_attributes(&mut line, &self.attributes);
        line
    }
}

/// `new BigInteger(text, 16)`: the digits without leading zeros, so two
/// `M5` values compare as numbers.
fn hex_value(text: &str) -> Result<String, String> {
    let digits = text.strip_prefix('+').unwrap_or(text);
    let (sign, digits) = match digits.strip_prefix('-') {
        Some(rest) if !text.starts_with('+') => ("-", rest),
        _ => ("", digits),
    };
    if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_hexdigit()) {
        return Err(format!("For input string: \"{text}\""));
    }
    let trimmed = digits.trim_start_matches('0').to_ascii_lowercase();
    Ok(if trimmed.is_empty() {
        "0".to_owned()
    } else {
        format!("{sign}{trimmed}")
    })
}

/// One `@RG` or `@PG` line: its `ID` and the other tags in their order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HeaderRecord {
    pub id: String,
    pub attributes: Vec<(String, String)>,
}

impl HeaderRecord {
    fn line(&self, kind: &str) -> String {
        let mut line = format!("@{kind}\tID:{}", self.id);
        push_attributes(&mut line, &self.attributes);
        line
    }
}

/// The header of an interval list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Header {
    /// `@HD` attributes, `VN` first.
    hd: Vec<(String, String)>,
    pub sequences: Vec<Sequence>,
    pub read_groups: Vec<HeaderRecord>,
    pub programs: Vec<HeaderRecord>,
    /// Each comment line whole, `@CO` included.
    pub comments: Vec<String>,
}

impl Default for Header {
    fn default() -> Self {
        Self::new(Vec::new())
    }
}

impl Header {
    /// `new SAMFileHeader(dictionary)`: `@HD VN:1.6` and the `@SQ` lines.
    pub fn new(sequences: Vec<Sequence>) -> Self {
        Self {
            hd: vec![("VN".to_owned(), VERSION.to_owned())],
            sequences,
            read_groups: Vec::new(),
            programs: Vec::new(),
            comments: Vec::new(),
        }
    }

    /// `SAMFileHeader.getSortOrder`: unsorted without `SO`.
    pub fn sort_order(&self) -> SortOrder {
        self.hd_attribute("SO").map_or(SortOrder::Unsorted, |so| {
            SortOrder::parse(so).unwrap_or(SortOrder::Unknown)
        })
    }

    /// `SAMFileHeader.setSortOrder`.
    pub fn set_sort_order(&mut self, order: SortOrder) {
        self.set_hd_attribute("SO", order.as_str());
    }

    pub fn hd_attribute(&self, key: &str) -> Option<&str> {
        self.hd
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    }

    /// `setAttribute`: an attribute already set keeps its place.
    fn set_hd_attribute(&mut self, key: &str, value: &str) {
        match self.hd.iter_mut().find(|(k, _)| k == key) {
            Some(slot) => slot.1 = value.to_owned(),
            None => self.hd.push((key.to_owned(), value.to_owned())),
        }
    }

    /// The index of the `@SQ` line named `name`.
    pub fn sequence_index(&self, name: &str) -> Option<usize> {
        self.sequences.iter().position(|s| s.name == name)
    }

    /// Each `@SQ` name with its index, for a lookup per record.
    pub fn name_index(&self) -> HashMap<String, usize> {
        self.sequences
            .iter()
            .enumerate()
            .map(|(at, s)| (s.name.clone(), at))
            .collect()
    }

    /// `SAMFileHeader.addComment`: `@CO\t` in front unless the text starts
    /// with it.
    pub fn add_comment(&mut self, text: &str) {
        if text.starts_with("@CO\t") {
            self.comments.push(text.to_owned());
        } else {
            self.comments.push(format!("@CO\t{text}"));
        }
    }

    /// The `@PG` line `IntervalListTools.doWork` adds: the smallest ID from
    /// 1 up that no `@PG` line has, `CL` then `PN`.
    pub fn add_program(&mut self, command_line: &str, program_name: &str) {
        let taken: HashSet<&str> = self.programs.iter().map(|p| p.id.as_str()).collect();
        let id = (1..=self.programs.len() + 1)
            .map(|n| n.to_string())
            .find(|id| !taken.contains(id.as_str()))
            .unwrap_or_default();
        self.programs.push(HeaderRecord {
            id,
            attributes: vec![
                ("CL".to_owned(), command_line.to_owned()),
                ("PN".to_owned(), program_name.to_owned()),
            ],
        });
    }

    /// `SequenceUtil.assertSequenceDictionariesEqual`: the same number of
    /// `@SQ` lines, each the same sequence as the other at its index.
    pub fn assert_same_dictionary(&self, other: &Header) -> Result<(), String> {
        let (a, b) = (&self.sequences, &other.sequences);
        if a.len() != b.len() {
            return Err(format!(
                "Sequence dictionaries are not the same size ({}, {})",
                a.len(),
                b.len()
            ));
        }
        for (index, (x, y)) in a.iter().zip(b).enumerate() {
            if !x.same_as(y)? {
                let tags = |s: &Sequence| {
                    s.attributes
                        .iter()
                        .map(|(k, v)| format!("/{k}={v}"))
                        .collect::<String>()
                };
                return Err(format!(
                    "Sequences at index {index} don't match: {index}/{}/{}{} {index}/{}/{}{}",
                    x.length,
                    x.name,
                    tags(x),
                    y.length,
                    y.name,
                    tags(y)
                ));
            }
        }
        Ok(())
    }

    /// `SAMTextHeaderCodec.decode` at `SILENT` of the header lines, each
    /// without its line end.
    pub fn decode<'a>(lines: impl IntoIterator<Item = &'a str>) -> Result<Header, String> {
        let mut header = Header::default();
        let mut read_group_ids = HashSet::new();
        let mut program_ids = HashSet::new();
        for line in lines {
            let Some(kept) = clean_header_line(line.as_bytes(), Stringency::Silent)
                .map_err(|e| e.to_string())?
            else {
                continue;
            };
            let kept = String::from_utf8_lossy(&kept).into_owned();
            let mut fields = kept.split('\t');
            let kind = fields.next().unwrap_or_default();
            if kind == "@CO" {
                header.comments.push(kept.clone());
                continue;
            }
            let mut pairs: Vec<(String, String)> = fields
                .filter_map(|f| f.split_once(':'))
                .map(|(k, v)| (k.to_owned(), v.to_owned()))
                .collect();
            let mut take = |key: &str| {
                pairs
                    .iter()
                    .position(|(k, _)| k == key)
                    .map(|at| pairs.remove(at).1)
            };
            match kind {
                "@HD" => {
                    for (key, value) in pairs {
                        header.set_hd_attribute(&key, &value);
                    }
                }
                "@SQ" => {
                    let name = truncate_name(&take("SN").unwrap_or_default()).to_owned();
                    if !legal_name(&name) {
                        return Err(format!(
                            "Sequence name '{name}' doesn't match regex: '{NAME_PATTERN}'"
                        ));
                    }
                    let length = take("LN").unwrap_or_default();
                    let length = crate::java::parse_integer(&length, 32)? as i32;
                    header.sequences.push(Sequence {
                        name,
                        length,
                        attributes: pairs,
                    });
                }
                "@RG" => {
                    let id = take("ID").unwrap_or_default();
                    for slot in pairs.iter_mut().filter(|(k, _)| k == "DT") {
                        if let Some(date) = crate::java_date::decode_date(&slot.1)?
                            .and_then(crate::java_date::iso8601)
                        {
                            slot.1 = date;
                        }
                    }
                    header.read_groups.push(HeaderRecord {
                        id,
                        attributes: pairs,
                    });
                }
                "@PG" => {
                    let id = take("ID").unwrap_or_default();
                    if !program_ids.insert(id.clone()) {
                        return Err(format!(
                            "Program record with group id {id} already exists in SAMFileHeader!"
                        ));
                    }
                    header.programs.push(HeaderRecord {
                        id,
                        attributes: pairs,
                    });
                }
                _ => {}
            }
        }
        let mut names = HashSet::new();
        for sequence in &header.sequences {
            if !names.insert(sequence.name.as_str()) {
                return Err(format!(
                    "Cannot add sequence that already exists in SAMSequenceDictionary: {}",
                    sequence.name
                ));
            }
        }
        for group in &header.read_groups {
            if !read_group_ids.insert(group.id.as_str()) {
                return Err(format!(
                    "Read group with group id {} already exists in SAMFileHeader!",
                    group.id
                ));
            }
        }
        Ok(header)
    }

    /// `SAMTextHeaderCodec.encode`: every header line, each ended by `\n`.
    pub fn encode(&self) -> String {
        let mut text = format!("@HD\tVN:{VERSION}");
        let rest: Vec<(String, String)> =
            self.hd.iter().filter(|(k, _)| k != "VN").cloned().collect();
        push_attributes(&mut text, &rest);
        text.push('\n');
        for sequence in &self.sequences {
            text.push_str(&sequence.line());
            text.push('\n');
        }
        for group in &self.read_groups {
            text.push_str(&group.line("RG"));
            text.push('\n');
        }
        for program in &self.programs {
            text.push_str(&program.line("PG"));
            text.push('\n');
        }
        for comment in &self.comments {
            text.push_str(comment);
            text.push('\n');
        }
        text
    }
}

fn push_attributes(line: &mut String, attributes: &[(String, String)]) {
    for (key, value) in attributes {
        line.push('\t');
        line.push_str(key);
        line.push(':');
        line.push_str(value);
    }
}

/// `SAMSequenceRecord.truncateSequenceName`: the name up to its first
/// space, line feed, vertical tab, form feed or carriage return.
fn truncate_name(name: &str) -> &str {
    let cut = name
        .find([' ', '\t', '\n', '\u{b}', '\u{c}', '\r'])
        .unwrap_or(name.len());
    &name[..cut]
}

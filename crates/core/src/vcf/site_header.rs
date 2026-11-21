//! The header of a VCF read for its sites as htsjdk 5.0.0's `VCFCodec`
//! reads it (`readActualHeader`, then `AbstractVCFCodec.parseHeaderFromLines`),
//! for `VcfToIntervalList`.
//!
//! - Every line from the start that begins with `##` is a header line.
//!   One that splits on `=` into exactly two parts, the first `format` or
//!   `fileformat` after a trim, names the version: `VCFv4.0` to `VCFv4.3`
//!   are read, the last such line wins, `VCFv3.2` (spelled `VCRv3.2`),
//!   `VCFv3.3` and `VCFv4.4` are refused, any other text is not a version.
//!   The first line that begins with one `#` is the column line and ends the
//!   header; a line before it that begins with no `#` is an error.
//! - The column line names the 8 fixed columns in order, then `FORMAT` and
//!   the samples (`String.split("\t")`, so a trailing tab is dropped); a
//!   repeated sample name is kept once (`LinkedHashSet`).
//! - `##INFO` and `##FORMAT` lines go through `VCF4Parser.parseLine` (the
//!   `<...>` text cut at commas outside quotes, keys and values trimmed, a
//!   later key replacing an earlier one in place) and the checks of
//!   `VCFCompoundHeaderLine`: the tag order, `Number`, `Type`, a missing
//!   `Description` allowed. `##FILTER`, `##ALT`, `##contig`, `##META`,
//!   `##SAMPLE` and, from VCF 4.3, `##PEDIGREE` need an `ID`. Any other
//!   `##` line with a `=` is kept as a key and value, one without is
//!   dropped. A repeated `ID` keeps the first line.
//!
//! The header is held whole while it is read, at most [`MAX_HEADER_BYTES`].

use std::collections::HashMap;
use std::io::BufRead;

use crate::error::excerpt;

/// The bytes of header text read at most (htsjdk has no bound).
pub const MAX_HEADER_BYTES: usize = 64 << 20;

/// A `##contig` line's `ID`, `length` and `assembly` as written.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContigLine {
    pub id: String,
    pub length: Option<String>,
    pub assembly: Option<String>,
}

/// What a site reader needs from the header.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SiteHeader {
    /// The minor version of `VCFv4.N`.
    pub minor: u8,
    /// The contig lines, the first of each `ID`, in file order.
    pub contigs: Vec<ContigLine>,
    /// The sample names, each once, in column order.
    pub samples: Vec<String>,
    /// The first `##INFO` line for `END` declares a `Flag`.
    pub end_is_flag: bool,
}

const COLUMNS: [&str; 8] = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"];

/// `TribbleException.InvalidHeader`'s message.
fn malformed(message: impl std::fmt::Display) -> String {
    format!("Your input file has a malformed header: {message}")
}

/// Read the header lines of `reader`, the column line last, leaving the
/// reader at the first record.
pub fn read(reader: &mut impl BufRead) -> Result<SiteHeader, String> {
    let mut lines = Vec::new();
    let mut held = 0usize;
    let mut minor = None;
    let mut raw = Vec::new();
    loop {
        raw.clear();
        let read = reader
            .read_until(b'\n', &mut raw)
            .map_err(|e| e.to_string())?;
        if read == 0 {
            return Err(malformed(
                "We never saw the required CHROM header line (starting with one #) for the input \
                 VCF file",
            ));
        }
        held += read;
        if held > MAX_HEADER_BYTES {
            return Err(format!(
                "the VCF header holds more than {MAX_HEADER_BYTES} bytes, the most read"
            ));
        }
        if raw.last() == Some(&b'\n') {
            raw.pop();
            if raw.last() == Some(&b'\r') {
                raw.pop();
            }
        }
        let line = crate::java::decode_utf8(&raw).into_owned();
        if line.starts_with("##") {
            if let Some(version) = version_line(&line)? {
                minor = Some(version);
            }
            lines.push(line);
        } else if line.starts_with('#') {
            let Some(minor) = minor else {
                return Err(malformed(
                    "We never saw a header line specifying VCF version",
                ));
            };
            return parse_lines(&lines, &line, minor);
        } else {
            return Err(malformed(
                "We never saw the required CHROM header line (starting with one #) for the input \
                 VCF file",
            ));
        }
    }
}

/// The minor version a `##fileformat` line names, `None` for any other
/// line; a version the codec does not read is an error.
fn version_line(line: &str) -> Result<Option<u8>, String> {
    let parts = java_split(&line[2..], '=');
    if parts.len() != 2 || !matches!(java_trim(parts[0]), "format" | "fileformat") {
        return Ok(None);
    }
    let known = match java_trim(parts[1]) {
        "VCRv3.2" | "VCFv3.3" | "VCFv4.4" => None,
        "VCFv4.0" => Some(0),
        "VCFv4.1" => Some(1),
        "VCFv4.2" => Some(2),
        "VCFv4.3" => Some(3),
        _ => {
            return Err(malformed(format!(
                "{} is not a supported version",
                parts[1]
            )));
        }
    };
    known.map(Some).ok_or_else(|| {
        malformed(format!(
            "This codec is strictly for VCFv4 and does not support {}",
            parts[1]
        ))
    })
}

fn parse_lines(lines: &[String], columns: &str, minor: u8) -> Result<SiteHeader, String> {
    let mut header = SiteHeader {
        minor,
        ..SiteHeader::default()
    };
    let mut contig_ids = std::collections::HashSet::new();
    let mut end_declared = false;
    for line in lines {
        if line.starts_with("##INFO") {
            let info = compound(after(line, 7)?, minor, false)?;
            if info.id == "END" && !end_declared {
                end_declared = true;
                header.end_is_flag = info.flag;
            }
        } else if line.starts_with("##FILTER") {
            simple(after(line, 9)?, "FILTER", Some(&["ID", "Description"]))?;
        } else if line.starts_with("##FORMAT") {
            compound(after(line, 9)?, minor, true)?;
        } else if line.starts_with("##contig") {
            let fields = simple(after(line, 9)?, "contig", None)?;
            let id = get(&fields, "ID").unwrap_or_default().to_owned();
            if contig_ids.insert(id.clone()) {
                header.contigs.push(ContigLine {
                    id,
                    length: get(&fields, "length").map(str::to_owned),
                    assembly: get(&fields, "assembly").map(str::to_owned),
                });
            }
        } else if line.starts_with("##ALT") {
            simple(after(line, 6)?, "ALT", Some(&["ID", "Description"]))?;
        } else if line.starts_with("##PEDIGREE") && minor >= 3 {
            simple(after(line, 11)?, "PEDIGREE", None)?;
        } else if line.starts_with("##META") {
            simple(after(line, 7)?, "META", None)?;
        } else if line.starts_with("##SAMPLE") {
            simple(after(line, 9)?, "SAMPLE", None)?;
        } else if let Some(equals) = line.find('=') {
            let key = &line[2..equals];
            if key.contains(['<', '>']) {
                return Err("VCFHeaderLine: key cannot contain angle brackets".to_owned());
            }
        }
    }
    header.samples = column_samples(columns)?;
    Ok(header)
}

/// `String.substring(begin)`: the text after the first `begin` characters,
/// an error when the line is shorter.
fn after(line: &str, begin: usize) -> Result<&str, String> {
    let length = line.chars().count();
    match line.char_indices().nth(begin) {
        Some((at, _)) => Ok(&line[at..]),
        None if begin == length => Ok(""),
        None => Err(format!("begin {begin}, end {length}, length {length}")),
    }
}

/// The column line's samples, each once.
fn column_samples(line: &str) -> Result<Vec<String>, String> {
    let fields = java_split(&line[1..], '\t');
    if fields.len() < COLUMNS.len() {
        return Err(malformed(format!(
            "there are not enough columns present in the header line: {}",
            excerpt(line)
        )));
    }
    for (expected, seen) in COLUMNS.iter().zip(&fields) {
        if !COLUMNS.contains(seen) {
            return Err(malformed(format!(
                "unknown column name '{seen}'; it does not match a legal column header name."
            )));
        }
        if expected != seen {
            return Err(malformed(format!(
                "we were expecting column name '{expected}' but we saw '{seen}'"
            )));
        }
    }
    let Some((format, names)) = fields[COLUMNS.len()..].split_first() else {
        return Ok(Vec::new());
    };
    if *format != "FORMAT" {
        return Err(malformed(format!(
            "we were expecting column name 'FORMAT' but we saw '{format}'"
        )));
    }
    if names.is_empty() {
        return Err(malformed(
            "The FORMAT field was provided but there is no genotype/sample data",
        ));
    }
    let mut seen = std::collections::HashSet::new();
    Ok(names
        .iter()
        .filter(|name| seen.insert(**name))
        .map(|name| (*name).to_owned())
        .collect())
}

/// An `##INFO` or `##FORMAT` line's `ID` and whether its type is `Flag`.
struct Compound {
    id: String,
    flag: bool,
}

/// `VCFCompoundHeaderLine(String, VCFHeaderVersion, lineType)`.
fn compound(line: &str, minor: u8, format: bool) -> Result<Compound, String> {
    let recommended: &[&str] = if minor >= 2 {
        &["Source", "Version"]
    } else {
        &[]
    };
    let fields = parse_fields(
        line,
        Some(&["ID", "Number", "Type", "Description"]),
        recommended,
    )?;
    let id = get(&fields, "ID");
    let Some(number) = get(&fields, "Number") else {
        return Err(format!(
            "the Number of the header line {} is missing",
            excerpt(line)
        ));
    };
    let count = match number {
        "A" | "R" | "G" | "." => None,
        text => {
            let count = crate::java::parse_integer(text, 32)?;
            if count < 0 {
                return Err(malformed(format!(
                    "Count < 0 for fixed size VCF header field {}",
                    id.unwrap_or("null")
                )));
            }
            Some(count)
        }
    };
    let ty = get(&fields, "Type");
    let flag = match ty {
        Some("Flag") => true,
        Some("Integer" | "Float" | "String" | "Character") => false,
        other => {
            return Err(format!(
                "{} is not a valid type in the VCF specification (note that types are \
                 case-sensitive)",
                other.unwrap_or("null")
            ));
        }
    };
    if flag && format {
        return Err(format!(
            "Flag is an unsupported type for this kind of field at line - {}",
            excerpt(line)
        ));
    }
    if !flag && count.is_some_and(|c| c <= 0) {
        return Err(format!(
            "Invalid count number, with fixed count the number should be 1 or higher: key={} \
             name={} type={} count={}",
            if format { "FORMAT" } else { "INFO" },
            id.unwrap_or("null"),
            ty.unwrap_or("null"),
            count.unwrap_or_default()
        ));
    }
    let Some(id) = id else {
        return Err(format!(
            "Invalid VCFCompoundHeaderLine: key={} name=null",
            if format { "FORMAT" } else { "INFO" }
        ));
    };
    check_id(id)?;
    Ok(Compound {
        id: id.to_owned(),
        flag,
    })
}

/// `VCFSimpleHeaderLine` over `VCF4Parser.parseLine`: an `ID` is needed.
fn simple(
    line: &str,
    key: &str,
    expected: Option<&[&str]>,
) -> Result<Vec<(String, String)>, String> {
    let fields = parse_fields(line, expected, &[])?;
    let Some(id) = get(&fields, "ID") else {
        return Err(format!("Invalid VCFSimpleHeaderLine: key={key} name=null"));
    };
    check_id(id)?;
    Ok(fields)
}

fn check_id(id: &str) -> Result<(), String> {
    if id.contains(['<', '>']) {
        return Err("VCFHeaderLine: ID cannot contain angle brackets".to_owned());
    }
    if id.contains('=') {
        return Err("VCFHeaderLine: ID cannot contain an equals sign".to_owned());
    }
    Ok(())
}

fn get<'a>(fields: &'a [(String, String)], key: &str) -> Option<&'a str> {
    fields
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.as_str())
}

/// `LinkedHashMap.put`: a key already present keeps its place. `index`
/// maps each key to its place in `fields`, so a line of many keys costs
/// one hash lookup per pair.
fn put(
    fields: &mut Vec<(String, String)>,
    index: &mut HashMap<String, usize>,
    key: &str,
    value: String,
) {
    match index.get(key) {
        Some(&at) => fields[at].1 = value,
        None => {
            index.insert(key.to_owned(), fields.len());
            fields.push((key.to_owned(), value));
        }
    }
}

/// `VCF4Parser.parseLine`: the `key=value` pairs of a `<...>` text in
/// order, then the tag order checked against `expected`.
fn parse_fields(
    line: &str,
    expected: Option<&[&str]>,
    recommended: &[&str],
) -> Result<Vec<(String, String)>, String> {
    let mut fields: Vec<(String, String)> = Vec::new();
    let mut index_of = HashMap::new();
    let mut builder = String::new();
    let mut key = String::new();
    let mut in_quote = false;
    let mut escape = false;
    let last = line.chars().count().saturating_sub(1);
    for (index, c) in line.chars().enumerate() {
        if c == '"' {
            if escape {
                builder.push(c);
                escape = false;
            } else {
                in_quote = !in_quote;
            }
        } else if in_quote {
            if escape {
                if c != '\\' {
                    builder.push('\\');
                }
                builder.push(c);
                escape = false;
            } else if c == '\\' {
                escape = true;
            } else {
                builder.push(c);
            }
        } else {
            escape = false;
            match c {
                '<' if index == 0 => {}
                '<' | '>' => {
                    if index == last {
                        put(
                            &mut fields,
                            &mut index_of,
                            &key,
                            java_trim(&builder).to_owned(),
                        );
                    }
                }
                '=' => {
                    key = java_trim(&builder).to_owned();
                    builder.clear();
                }
                ',' => {
                    put(
                        &mut fields,
                        &mut index_of,
                        &key,
                        java_trim(&builder).to_owned(),
                    );
                    builder.clear();
                }
                c => builder.push(c),
            }
        }
    }
    if in_quote {
        return Err(malformed(format!(
            "Unclosed quote in header line value {}",
            excerpt(line)
        )));
    }
    let Some(expected) = expected else {
        return Ok(fields);
    };
    if fields.is_empty() && !expected.is_empty() {
        return Err(malformed(format!(
            "Header with no tags is not supported when there are expected tags in line {}",
            excerpt(line)
        )));
    }
    for (index, (key, _)) in fields.iter().enumerate().take(expected.len()) {
        if expected[index] == key {
            continue;
        }
        if let Some(at) = expected.iter().position(|e| e == key) {
            return Err(malformed(format!(
                "Tag {key} in wrong order (was #{}, expected #{}) in line {}",
                index + 1,
                at + 1,
                excerpt(line)
            )));
        }
        if recommended.contains(&key.as_str()) {
            return Err(malformed(format!(
                "Recommended tag {key} must be listed after all expected tags in line {}",
                excerpt(line)
            )));
        }
        return Err(malformed(format!(
            "Unexpected tag {key} in line {}",
            excerpt(line)
        )));
    }
    Ok(fields)
}

/// `String.trim`: characters up to U+0020 removed from both ends.
fn java_trim(text: &str) -> &str {
    text.trim_matches(|c: char| c <= ' ')
}

/// `String.split` on one character: trailing empty parts dropped.
fn java_split(text: &str, separator: char) -> Vec<&str> {
    let mut parts: Vec<&str> = text.split(separator).collect();
    while parts.last().is_some_and(|p| p.is_empty()) {
        parts.pop();
    }
    parts
}

#[cfg(test)]
mod tests {
    use super::*;

    const COLUMN_LINE: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    fn header(lines: &str) -> Result<SiteHeader, String> {
        read(&mut lines.as_bytes())
    }

    fn with(meta: &str) -> Result<SiteHeader, String> {
        header(&format!("##fileformat=VCFv4.2\n{meta}{COLUMN_LINE}"))
    }

    #[test]
    fn a_repeated_key_keeps_its_first_place() {
        let fields = parse_fields("<ID=a,B=1,C=2,B=3>", None, &[]).expect("parses");
        assert_eq!(
            fields,
            [("ID", "a"), ("B", "3"), ("C", "2")].map(|(k, v)| (k.to_owned(), v.to_owned()))
        );
    }

    /// One `put` per pair costs a hash lookup: 200,000 distinct keys took
    /// minutes when each pair searched the pairs before it.
    #[test]
    fn many_distinct_keys_parse_in_linear_time() {
        const KEYS: usize = 200_000;
        let pairs: Vec<String> = (0..KEYS).map(|i| format!("K{i}={i}")).collect();
        let line = format!("<{}>", pairs.join(","));
        let started = std::time::Instant::now();
        let fields = parse_fields(&line, None, &[]).expect("parses");
        let elapsed = started.elapsed();
        assert_eq!(fields.len(), KEYS);
        assert_eq!(
            fields[KEYS - 1],
            (format!("K{}", KEYS - 1), (KEYS - 1).to_string())
        );
        assert!(elapsed.as_secs_f64() < 1.0, "{elapsed:?}");
    }

    #[test]
    fn versions_follow_the_codec() {
        assert_eq!(with("").map(|h| h.minor), Ok(2));
        assert_eq!(
            header(&format!(
                "##contig=<ID=c1>\n##fileformat=VCFv4.3\n{COLUMN_LINE}"
            ))
            .map(|h| h.minor),
            Ok(3)
        );
        assert_eq!(
            header(&format!(
                "##fileformat=VCFv4.1\n##fileformat=VCFv4.3\n{COLUMN_LINE}"
            ))
            .map(|h| h.minor),
            Ok(3)
        );
        assert_eq!(
            header(&format!("##fileformat = VCFv4.0 \n{COLUMN_LINE}")).map(|h| h.minor),
            Ok(0)
        );
        for bad in ["VCFv3.3", "VCFv4.4", "VCRv3.2"] {
            let error = header(&format!("##fileformat={bad}\n{COLUMN_LINE}")).unwrap_err();
            assert!(error.contains("does not support"), "{error}");
        }
        let error = header(&format!("##fileformat=VCFv5.0\n{COLUMN_LINE}")).unwrap_err();
        assert!(
            error.contains("VCFv5.0 is not a supported version"),
            "{error}"
        );
        assert!(
            header(COLUMN_LINE)
                .unwrap_err()
                .contains("specifying VCF version")
        );
        assert!(header("##fileformat=VCFv4.2\n").is_err());
        assert!(header("##fileformat=VCFv4.2\nc1\t1\n").is_err());
    }

    #[test]
    fn meta_lines_follow_htsjdk() {
        for ok in [
            "##INFO=<ID=X,Number=1,Type=Integer>\n",
            "##INFO=<ID=X,Number=1,Type=Integer,Description=\"a\">\n\
             ##INFO=<ID=X,Number=1,Type=Integer,Description=\"b\">\n",
            "##FILTER=<ID=q,Description=\"a\">\n##FILTER=<ID=q,Description=\"b\">\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"a\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"b\">\n",
            "##foo\n",
            "##foo=bar\n",
            "##ALT=<ID=DEL,Description=\"d\">\n",
            "##PEDIGREE=<Name_0=a>\n",
            "##INFO=<ID=X,Number=0,Type=Flag,Description=\"a\">\n",
            "##INFO=<ID=X,Number=.,Type=Float,Description=\"a, \\\"b\\\"\">\n",
        ] {
            assert!(with(ok).is_ok(), "{ok}");
        }
        for bad in [
            "##INFO=<ID=X,Number=1,Type=Bogus,Description=\"a\">\n",
            "##INFO=<ID=X,Number=zz,Type=Integer,Description=\"a\">\n",
            "##INFO=<ID=X,Number=0,Type=Integer,Description=\"a\">\n",
            "##INFO=<ID=X,Number=-2,Type=Integer,Description=\"a\">\n",
            "##INFO=<Number=1,ID=X,Type=Integer,Description=\"a\">\n",
            "##INFO=<ID=X,Number=1,Type=Integer,Foo=1>\n",
            "##INFO=<ID=X,Type=Integer>\n",
            "##FORMAT=<ID=X,Number=0,Type=Flag,Description=\"a\">\n",
            "##FILTER=<Description=\"a\">\n",
            "##contig=<>\n",
            "##contig=<length=5>\n",
            "##META=<Type=String>\n",
            "##INFO=<ID=X,Number=1,Type=Integer,Description=\"a>\n",
            "##fo<o=bar\n",
            "##INFO\n",
            "##ALT\n",
            "##contig=<ID=a=b>\n",
        ] {
            assert!(with(bad).is_err(), "{bad}");
        }
        let pedigree = header(&format!(
            "##fileformat=VCFv4.3\n##PEDIGREE=<Name_0=a>\n{COLUMN_LINE}"
        ));
        assert!(pedigree.is_err());
    }

    #[test]
    fn contig_lines_are_cut_as_vcf4_parser_cuts_them() {
        let contigs = |meta: &str| with(meta).expect("reads").contigs;
        let one = |id: &str, length: Option<&str>, assembly: Option<&str>| ContigLine {
            id: id.to_owned(),
            length: length.map(str::to_owned),
            assembly: assembly.map(str::to_owned),
        };
        assert_eq!(
            contigs("##contig=<ID=c1, length=1000>\n"),
            [one("c1", Some("1000"), None)]
        );
        assert_eq!(
            contigs("##contig=<ID=c1,length= 100>\n"),
            [one("c1", Some("100"), None)]
        );
        assert_eq!(
            contigs("##contig=<ID=c1,length=1000,length=2000>\n"),
            [one("c1", Some("2000"), None)]
        );
        assert_eq!(
            contigs("##contig=<ID=c1,length=1000,>\n"),
            [one("c1", Some(""), None)]
        );
        assert_eq!(
            contigs("##contig=<ID=c1,length=5,assembly=\"b 37\">\n##contig=<ID=c1,length=9>\n"),
            [one("c1", Some("5"), Some("b 37"))]
        );
    }

    #[test]
    fn columns_and_samples_follow_the_codec() {
        let samples = |line: &str| header(&format!("##fileformat=VCFv4.2\n{line}\n"));
        let cols = COLUMN_LINE.trim_end();
        assert_eq!(samples(&format!("{cols}\t")).map(|h| h.samples), Ok(vec![]));
        assert_eq!(
            samples(&format!("{cols}\tFORMAT\tS\tS\tT\t")).map(|h| h.samples),
            Ok(vec!["S".to_owned(), "T".to_owned()])
        );
        assert!(samples(&format!("{cols}\tFORMAT")).is_err());
        assert!(samples(&format!("{cols}\tS1")).is_err());
        assert!(samples("#CHROM\tPOS\tID\tREF\tALT").is_err());
        assert!(samples("#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").is_err());
        assert!(samples(&cols.replace('\t', " ")).is_err());
    }

    #[test]
    fn end_takes_the_first_info_line() {
        let flag = |meta: &str| with(meta).expect("reads").end_is_flag;
        assert!(flag(
            "##INFO=<ID=END,Number=0,Type=Flag,Description=\"a\">\n"
        ));
        assert!(!flag(
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"a\">\n\
             ##INFO=<ID=END,Number=0,Type=Flag,Description=\"a\">\n"
        ));
    }

    #[test]
    fn the_header_is_capped() {
        let big = format!(
            "##fileformat=VCFv4.2\n##x={}\n{COLUMN_LINE}",
            "a".repeat(MAX_HEADER_BYTES)
        );
        assert!(header(&big).unwrap_err().contains("most read"));
    }
}

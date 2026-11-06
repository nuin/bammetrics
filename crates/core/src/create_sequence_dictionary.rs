//! `CreateSequenceDictionary`, mirroring Picard 3.5.0
//! `CreateSequenceDictionary.doWork` with `SequenceDictionaryUtils`,
//! htsjdk's `FastaSequenceFile` and `SAMSequenceDictionaryCodec`.
//!
//! The FASTA is read from start to end. A name ending in `.gz` (or `.gzip`,
//! `.bgz`, `.bgzf`) is inflated member by member, so plain gzip and bgzip
//! both work; a gzip file under another name is read as text and fails at
//! its first byte, as in Picard. The file name must end in one of
//! [`FASTA_EXTENSIONS`]. Line ends are `\n`, `\r\n` or `\r`; blank lines
//! are skipped; a `>` starts a header only at the start of a line. The
//! name is the header text trimmed of characters up to `U+0020` at both
//! ends and, under `TRUNCATE_NAMES_AT_WHITESPACE` (default true), cut at
//! the first space, tab, line feed, vertical tab, form feed or carriage
//! return; it must match the SAM name pattern. A sequence line keeps every
//! byte but its trailing Java whitespace (tab to carriage return, `0x1C` to
//! `0x1F`, space), so a line of spaces adds nothing and a space inside a
//! line is a base. `LN` counts the bases kept, `M5` is the MD5 of them
//! upper-cased (`a` to `z` only; `N`, IUPAC codes, `*` and `-` stay as they
//! are), printed as 32 lower-case hex digits.
//!
//! The file is `@HD VN:1.6`, then one `@SQ` line per sequence: `SN`, `LN`,
//! `M5`, `AS` when given, `UR` (default `file://` and the FASTA's absolute
//! path as `java.nio.file.Path.toUri` encodes it), `SP` when given, `AN`
//! when `ALT_NAMES` names the sequence. Duplicate names and empty
//! sequences are written as they come (`LN:0` has the MD5 of no bytes).
//! `NUM_SEQUENCES` stops after that many lines, never fewer than one, and
//! the sequence after the last one written is still read, as Picard's
//! iterator reads one ahead. An existing `OUTPUT` is an error. On any error
//! no file is written: Picard leaves the lines written before the failing
//! sequence.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::ops::Range;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};
use crate::gzip::GzipReader;
use crate::java;
use crate::md5::{self, Md5};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

/// The columns of the rows returned to Python, one per `@SQ` line.
pub const COLUMNS: [&str; 7] = [
    "name",
    "length",
    "md5",
    "assembly",
    "uri",
    "species",
    "alternate_names",
];

/// htsjdk `FileExtensions.FASTA`: the FASTA name must end in one, case and
/// all.
pub const FASTA_EXTENSIONS: [&str; 10] = [
    ".fasta",
    ".fasta.gz",
    ".fas",
    ".fas.gz",
    ".fa",
    ".fa.gz",
    ".fna",
    ".fna.gz",
    ".txt",
    ".txt.gz",
];

/// htsjdk `FileExtensions.GZIP`: `IOUtil.openFileForReading` inflates a
/// file whose name ends in one.
const GZIP_EXTENSIONS: [&str; 4] = [".gz", ".gzip", ".bgz", ".bgzf"];

/// `SAMFileHeader.CURRENT_VERSION` in htsjdk 5.0.0.
const HEADER_VERSION: &str = "1.6";

/// `FastaSequenceFile.readSequenceName` reads a name into 4096 bytes.
const NAME_LIMIT: usize = 4096;

/// `SAMSequenceRecord.LEGAL_RNAME_PATTERN`, written out as its two
/// character classes.
const NAME_PATTERN: &str = "[0-9A-Za-z!#$%&+./:;?@^_|~-][0-9A-Za-z!#$%&*+./:;=?@^_|~-]*";

pub const TOOL: Tool = Tool {
    name: "CreateSequenceDictionary",
    summary: "one row per FASTA sequence (name, length, md5, assembly, uri, species, \
              alternate_names); OUTPUT gets the .dict.",
    options: &[
        Opt::path(
            &["REFERENCE", "R"],
            "The FASTA: .fasta, .fas, .fa, .fna or .txt, each also as .gz.",
        )
        .required(),
        Opt::text(
            &["GENOME_ASSEMBLY", "AS"],
            "The AS field of every @SQ line.",
        ),
        Opt::text(
            &["URI", "UR"],
            "The UR field of every @SQ line (default file:// and the FASTA's absolute path).",
        ),
        Opt::text(&["SPECIES", "SP"], "The SP field of every @SQ line."),
        Opt::bool(
            &["TRUNCATE_NAMES_AT_WHITESPACE"],
            "Name a sequence by the first word of its > line (default true).",
        ),
        Opt::int(
            &["NUM_SEQUENCES"],
            "Stop after this many sequences; at least one is written.",
        ),
        Opt::path(
            &["ALT_NAMES", "AN"],
            "Lines of name, tab, alternative name: the AN field.",
        ),
        Opt::bool(
            &["CREATE_MD5_FILE"],
            "Also write OUTPUT.md5, the MD5 of the dictionary (default false).",
        ),
        Opt::rejected(&["INPUT", "I"]),
        Opt::rejected(&["ASSUME_SORTED"]),
        Opt::rejected(&["REFERENCE_SEQUENCE"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

#[derive(Clone, Debug)]
pub struct Options {
    pub reference: PathBuf,
    pub output: Option<PathBuf>,
    pub genome_assembly: Option<String>,
    pub uri: Option<String>,
    pub species: Option<String>,
    pub truncate_names_at_whitespace: bool,
    pub num_sequences: i32,
    pub alt_names: Option<PathBuf>,
    pub create_md5_file: bool,
}

impl Options {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let reference = inv.get("REFERENCE").ok_or_else(|| {
            Error::Usage("CreateSequenceDictionary: REFERENCE (R) is required".to_owned())
        })?;
        Ok(Self {
            reference: reference.into(),
            output: inv.output.clone(),
            genome_assembly: inv.get("GENOME_ASSEMBLY").map(str::to_owned),
            uri: inv.get("URI").map(str::to_owned),
            species: inv.get("SPECIES").map(str::to_owned),
            truncate_names_at_whitespace: inv.flag("TRUNCATE_NAMES_AT_WHITESPACE", true)?,
            num_sequences: inv.parse("NUM_SEQUENCES", i32::MAX)?,
            alt_names: inv.get("ALT_NAMES").map(PathBuf::from),
            create_md5_file: inv.flag("CREATE_MD5_FILE", false)?,
        })
    }
}

/// One FASTA sequence as the dictionary describes it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sequence {
    pub name: String,
    pub length: u64,
    pub md5: String,
}

/// `path` made absolute against the working directory, as
/// `Path.toAbsolutePath` makes it (no `.` or `..` is resolved).
fn absolute(path: &Path) -> Result<PathBuf> {
    if path.is_absolute() {
        return Ok(path.to_path_buf());
    }
    let cwd = std::env::current_dir().map_err(|source| Error::Input {
        path: PathBuf::from("."),
        source,
    })?;
    Ok(cwd.join(path))
}

/// A byte `UnixUriUtils.toUri` writes as itself: letters, digits and
/// `-_.!~*'()`, `:@&=+$,` and `;/`.
fn uri_byte(b: u8) -> bool {
    b.is_ascii_alphanumeric() || b"-_.!~*'():@&=+$,;/".contains(&b)
}

/// `Path.toUri().toString()` of a file: `file://`, then the absolute path
/// with repeated and trailing slashes dropped and every other byte as
/// `%XX` in upper-case hex. Picard's default `UR` and the name in its
/// `already exists` message.
pub fn file_uri(path: &Path) -> Result<String> {
    let absolute = absolute(path)?;
    let mut uri = String::from("file://");
    let mut previous_slash = false;
    for &b in absolute.as_os_str().as_encoded_bytes() {
        if b == b'/' && previous_slash {
            continue;
        }
        previous_slash = b == b'/';
        if uri_byte(b) {
            uri.push(char::from(b));
        } else {
            uri.push_str(&format!("%{b:02X}"));
        }
    }
    if uri.len() > "file:///".len() && uri.ends_with('/') {
        uri.pop();
    }
    Ok(uri)
}

/// The extension of [`FASTA_EXTENSIONS`] that `name` ends in.
fn fasta_extension(name: &str) -> Option<&'static str> {
    FASTA_EXTENSIONS.iter().copied().find(|e| name.ends_with(e))
}

/// `ReferenceSequenceFileFactory.getFastaExtension` failing.
fn unsupported(path: &Path) -> Result<Error> {
    Ok(Error::Metric(format!(
        "File is not a supported reference file type: {}",
        absolute(path)?.display()
    )))
}

/// `ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence`:
/// the FASTA's extension replaced by `.dict`, beside it.
pub fn default_dictionary(reference: &Path) -> Result<PathBuf> {
    let name = reference
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let Some(extension) = fasta_extension(&name) else {
        return Err(unsupported(reference)?);
    };
    Ok(reference.with_file_name(format!("{}.dict", &name[..name.len() - extension.len()])))
}

/// The Picard face's arguments with `OUTPUT` added when none is given:
/// Picard's `customCommandLineValidation` names [`default_dictionary`] of
/// `REFERENCE`. Without `REFERENCE` they are returned as given, and the
/// run reports it.
pub fn with_default_output(args: &[String]) -> Result<Vec<String>> {
    let inv = TOOL.parse(args)?;
    let mut args = args.to_vec();
    if inv.output.is_none()
        && let Some(reference) = inv.get("REFERENCE")
    {
        let output = default_dictionary(Path::new(reference))?;
        args.push(format!("OUTPUT={}", output.display()));
    }
    Ok(args)
}

/// `Character.isWhitespace` of a byte read as Latin-1.
fn java_whitespace(b: u8) -> bool {
    matches!(b, b'\t'..=b'\r' | 0x1c..=0x1f | b' ')
}

/// `Character.isWhitespace` of a `char`.
fn java_whitespace_char(c: char) -> bool {
    u8::try_from(c).map_or_else(
        |_| c.is_whitespace() && !matches!(c, '\u{2007}' | '\u{202f}'),
        java_whitespace,
    )
}

/// `StringUtil.isBlank`: empty or Java whitespace only.
fn is_blank(text: &str) -> bool {
    text.chars().all(java_whitespace_char)
}

/// `SAMSequenceRecord.validateSequenceName`: [`NAME_PATTERN`].
fn legal_name(name: &str) -> bool {
    let first = |b: u8| b.is_ascii_alphanumeric() || b"!#$%&+./:;?@^_|~-".contains(&b);
    let rest = |b: u8| first(b) || b == b'*' || b == b'=';
    let mut bytes = name.bytes();
    bytes.next().is_some_and(first) && bytes.all(rest)
}

fn check_name(name: &str) -> Result<()> {
    if legal_name(name) {
        Ok(())
    } else {
        Err(Error::Metric(format!(
            "Sequence name '{name}' doesn't match regex: '{NAME_PATTERN}'"
        )))
    }
}

/// `FastaSequenceFile` over a plain or gzip stream: one [`Sequence`] per
/// header, bases hashed as they are read.
struct FastaReader {
    input: Box<dyn BufRead>,
    /// The absolute path, for messages (`getSource`).
    source: String,
    truncate: bool,
    /// One line of the file, without its `\n`.
    buf: Vec<u8>,
    /// Where the unread part of `buf` starts; `None` when `buf` is used up.
    start: Option<usize>,
    /// The header (after `>`) that ended the previous sequence.
    pending: Option<Vec<u8>>,
    index: usize,
    upper: Vec<u8>,
}

impl FastaReader {
    fn open(path: &Path, truncate: bool) -> Result<Self> {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        if fasta_extension(&name).is_none() {
            return Err(unsupported(path)?);
        }
        let file = File::open(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        let capacity = 1 << 16;
        let input: Box<dyn BufRead> = if GZIP_EXTENSIONS.iter().any(|e| name.ends_with(e)) {
            Box::new(BufReader::with_capacity(
                capacity,
                GzipReader::new(BufReader::with_capacity(capacity, file)),
            ))
        } else {
            Box::new(BufReader::with_capacity(capacity, file))
        };
        Ok(Self {
            input,
            source: absolute(path)?.display().to_string(),
            truncate,
            buf: Vec::new(),
            start: None,
            pending: None,
            index: 0,
            upper: Vec::new(),
        })
    }

    /// The next line of `buf`, ended by `\n` or `\r` as `FastLineReader`
    /// ends lines (so `\r\n` gives a line and an empty one); `None` at the
    /// end of the file.
    fn next_line(&mut self) -> Result<Option<Range<usize>>> {
        let start = match self.start {
            Some(start) => start,
            None => {
                self.buf.clear();
                let read = self.input.read_until(b'\n', &mut self.buf).map_err(|e| {
                    Error::Metric(format!("Error reading FASTA {}: {e}", self.source))
                })?;
                if read == 0 {
                    return Ok(None);
                }
                if self.buf.last() == Some(&b'\n') {
                    self.buf.pop();
                }
                0
            }
        };
        match self.buf[start..].iter().position(|&b| b == b'\r') {
            Some(end) => {
                self.start = Some(start + end + 1);
                Ok(Some(start..start + end))
            }
            None => {
                self.start = None;
                Ok(Some(start..self.buf.len()))
            }
        }
    }

    /// The header text after `>` of the first sequence, or `None` for a
    /// file of blank lines.
    fn first_header(&mut self) -> Result<Option<Vec<u8>>> {
        while let Some(range) = self.next_line()? {
            let line = &self.buf[range];
            let Some(&first) = line.first() else {
                continue;
            };
            if first != b'>' {
                return Err(Error::Metric(format!(
                    "Format exception reading FASTA {}.  Expected > but saw chr({}) at start of \
                     sequence with index {}",
                    self.source,
                    first.cast_signed(),
                    self.index
                )));
            }
            return Ok(Some(line[1..].to_vec()));
        }
        Ok(None)
    }

    /// `readSequenceName`: the two length checks, then trim and truncate.
    fn name(&self, header: &[u8]) -> Result<String> {
        if header.is_empty() {
            return Err(Error::Metric(format!(
                "Missing sequence name in FASTA {}",
                self.source
            )));
        }
        if header.len() > NAME_LIMIT {
            return Err(Error::Metric(format!(
                "Sequence name too long in FASTA {}",
                self.source
            )));
        }
        // `StringUtil.bytesToString` maps each byte to one char.
        let text: String = header.iter().map(|&b| char::from(b)).collect();
        let mut name = text.trim_matches(|c: char| c <= ' ').to_owned();
        if self.truncate
            && let Some(cut) = name.find([' ', '\t', '\n', '\u{b}', '\u{c}', '\r'])
        {
            name.truncate(cut);
        }
        Ok(name)
    }

    /// The next sequence, or `None` at the end of the file.
    fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        let header = match self.pending.take() {
            Some(header) => header,
            None if self.index == 0 => match self.first_header()? {
                Some(header) => header,
                None => return Ok(None),
            },
            None => return Ok(None),
        };
        let name = self.name(&header)?;
        let mut digest = Md5::new();
        let mut length = 0u64;
        while let Some(range) = self.next_line()? {
            let line = &self.buf[range];
            match line.first() {
                None => continue,
                Some(b'>') => {
                    self.pending = Some(line[1..].to_vec());
                    break;
                }
                Some(_) => {}
            }
            let kept = line
                .iter()
                .rposition(|&b| !java_whitespace(b))
                .map_or(0, |i| i + 1);
            self.upper.clear();
            self.upper
                .extend(line[..kept].iter().map(u8::to_ascii_uppercase));
            digest.update(&self.upper);
            length += kept as u64;
        }
        self.index += 1;
        Ok(Some(Sequence {
            name,
            length,
            md5: md5::hex(&digest.finish()),
        }))
    }
}

/// Every sequence `CreateSequenceDictionary` writes: at most
/// `num_sequences` (at least one), each name checked, one sequence read
/// past the last.
pub fn read_sequences(path: &Path, truncate: bool, num_sequences: i32) -> Result<Vec<Sequence>> {
    let mut reader = FastaReader::open(path, truncate)?;
    let mut sequences = Vec::new();
    let mut next = reader.next_sequence()?;
    while let Some(sequence) = next {
        check_name(&sequence.name)?;
        next = reader.next_sequence()?;
        sequences.push(sequence);
        if sequences.len() >= usize::try_from(num_sequences).unwrap_or(0) {
            break;
        }
    }
    Ok(sequences)
}

/// The lines of a text file as `BufferedReader.readLine` splits them: at
/// `\n`, `\r\n` or `\r`, with no empty line after a final line end.
fn read_lines(text: &str) -> Vec<&str> {
    let mut lines = Vec::new();
    let mut rest = text;
    while !rest.is_empty() {
        match rest.find(['\n', '\r']) {
            Some(end) => {
                lines.push(&rest[..end]);
                let skip = if rest[end..].starts_with("\r\n") {
                    2
                } else {
                    1
                };
                rest = &rest[end + skip..];
            }
            None => {
                lines.push(rest);
                rest = "";
            }
        }
    }
    lines
}

/// `loadContigAliasesMap`: each contig's alternative names, in the order
/// its `HashSet` iterates them. An error carries Picard's message.
pub fn read_alt_names(path: &Path) -> Result<HashMap<String, Vec<String>>> {
    let fail = |message: String| {
        Error::Metric(format!(
            "Can't read alias file {}: {message}",
            path.display()
        ))
    };
    let bytes = std::fs::read(path).map_err(|e| fail(e.to_string()))?;
    let text = String::from_utf8_lossy(&bytes);
    // Contigs in insertion order, each with its names in insertion order.
    let mut contigs: Vec<(String, Vec<String>)> = Vec::new();
    for line in read_lines(&text) {
        if is_blank(line) {
            continue;
        }
        let Some((contig, alt)) = line.split_once('\t') else {
            return Err(fail(format!("tabulation missing in {line}")));
        };
        if is_blank(contig) {
            return Err(fail(format!("empty contig in {line}")));
        }
        if is_blank(alt) {
            return Err(fail(format!("empty alternative name in {line}")));
        }
        if alt == contig {
            continue;
        }
        if !legal_name(alt) {
            return Err(fail(format!(
                "Illegal alternative reference sequence name in {line}"
            )));
        }
        if contigs.iter().any(|(c, _)| c == alt) {
            return Err(fail(format!(
                "alternate name {alt} previously defined as a contig in {line}"
            )));
        }
        if contigs
            .iter()
            .any(|(c, alts)| c != contig && alts.iter().any(|a| a == contig))
        {
            return Err(fail(format!(
                "contig  {contig} previously defined as an alternate name in {line}"
            )));
        }
        let index = match contigs.iter().position(|(c, _)| c == contig) {
            Some(index) => index,
            None => {
                contigs.push((contig.to_owned(), Vec::new()));
                contigs.len() - 1
            }
        };
        let alts = &mut contigs[index].1;
        if !alts.iter().any(|a| a == alt) {
            alts.push(alt.to_owned());
        }
    }
    Ok(contigs
        .into_iter()
        .map(|(contig, alts)| {
            let keys: Vec<&str> = alts.iter().map(String::as_str).collect();
            let ordered = java::hash_map_order(&keys)
                .into_iter()
                .map(|i| alts[i].clone())
                .collect();
            (contig, ordered)
        })
        .collect())
}

/// One `@SQ` line's fields beyond `SN` and `LN`.
#[derive(Clone, Debug, Default)]
pub struct Fields {
    pub genome_assembly: Option<String>,
    pub uri: String,
    pub species: Option<String>,
    pub alt_names: HashMap<String, Vec<String>>,
}

/// The dictionary text and the Python rows.
pub fn dictionary(sequences: &[Sequence], fields: &Fields) -> (String, Vec<Vec<Value>>) {
    let mut text = format!("@HD\tVN:{HEADER_VERSION}\n");
    let mut rows = Vec::new();
    for sequence in sequences {
        let alternate = fields
            .alt_names
            .get(&sequence.name)
            .map(|names| names.join(","));
        let mut line = format!(
            "@SQ\tSN:{}\tLN:{}\tM5:{}",
            sequence.name, sequence.length, sequence.md5
        );
        if let Some(assembly) = &fields.genome_assembly {
            line.push_str(&format!("\tAS:{assembly}"));
        }
        line.push_str(&format!("\tUR:{}", fields.uri));
        if let Some(species) = &fields.species {
            line.push_str(&format!("\tSP:{species}"));
        }
        if let Some(alternate) = &alternate {
            line.push_str(&format!("\tAN:{alternate}"));
        }
        text.push_str(&line);
        text.push('\n');
        rows.push(vec![
            Value::Str(sequence.name.clone()),
            Value::Int(i64::try_from(sequence.length).unwrap_or(i64::MAX)),
            Value::Str(sequence.md5.clone()),
            fields.genome_assembly.clone().into(),
            Value::Str(fields.uri.clone()),
            fields.species.clone().into(),
            alternate.into(),
        ]);
    }
    (text, rows)
}

/// `doWork`'s checks before the FASTA is opened: an existing `OUTPUT`,
/// then a missing directory for it (`IOUtil.assertFileIsWritable`).
fn check_output(output: &Path) -> Result<()> {
    if output.exists() {
        return Err(Error::Metric(format!(
            "{} already exists.  Delete this file and try again, or specify a different output \
             file.",
            file_uri(output)?
        )));
    }
    let parent = absolute(output)?.parent().map(Path::to_path_buf);
    if parent.is_some_and(|p| !p.is_dir()) {
        return Err(Error::Metric(format!(
            "Cannot write file: {}. Neither file nor parent directory exist.",
            absolute(output)?.display()
        )));
    }
    Ok(())
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    match &opts.output {
        Some(output) => check_output(output)?,
        None if opts.create_md5_file => {
            return Err(Error::Usage(
                "CreateSequenceDictionary: CREATE_MD5_FILE=true needs OUTPUT".to_owned(),
            ));
        }
        None => {}
    }
    let alt_names = match &opts.alt_names {
        Some(path) => read_alt_names(path)?,
        None => HashMap::new(),
    };
    let uri = match &opts.uri {
        Some(uri) => uri.clone(),
        None => file_uri(&opts.reference)?,
    };
    let sequences = read_sequences(
        &opts.reference,
        opts.truncate_names_at_whitespace,
        opts.num_sequences,
    )?;
    let fields = Fields {
        genome_assembly: opts.genome_assembly,
        uri,
        species: opts.species,
        alt_names,
    };
    let (text, rows) = dictionary(&sequences, &fields);
    let mut output = ToolOutput::file(MetricsOutput::text(&COLUMNS, rows, text.clone()));
    if opts.create_md5_file
        && let Some(path) = &opts.output
    {
        let mut digest = Md5::new();
        digest.update(text.as_bytes());
        let mut name = path.as_os_str().to_owned();
        name.push(".md5");
        output.extra.push((
            PathBuf::from(name),
            MetricsOutput::text(&[], Vec::new(), md5::hex(&digest.finish())),
        ));
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Files under the test's own temporary directory, removed on drop.
    struct Scratch(PathBuf);

    impl Scratch {
        fn new(test: &str) -> Self {
            let dir =
                std::env::temp_dir().join(format!("bammetrics-csd-{}-{test}", std::process::id()));
            std::fs::create_dir_all(&dir).expect("scratch directory");
            Self(dir)
        }

        fn write(&self, name: &str, bytes: &[u8]) -> PathBuf {
            let path = self.0.join(name);
            std::fs::write(&path, bytes).expect("write fixture");
            path
        }
    }

    impl Drop for Scratch {
        fn drop(&mut self) {
            let _ = std::fs::remove_dir_all(&self.0);
        }
    }

    fn sequences(text: &[u8], truncate: bool, num: i32) -> Result<Vec<Sequence>> {
        let mut digest = Md5::new();
        digest.update(text);
        let scratch = Scratch::new(&md5::hex(&digest.finish()));
        read_sequences(&scratch.write("x.fa", text), truncate, num)
    }

    fn summary(sequences: &[Sequence]) -> Vec<(String, u64, String)> {
        sequences
            .iter()
            .map(|s| (s.name.clone(), s.length, s.md5.clone()))
            .collect()
    }

    /// Digests Picard 3.5.0 wrote for the same files.
    #[test]
    fn bases_are_hashed_upper_case_with_line_ends_trimmed() {
        let iupac: &[u8] =
            b">s1\nacgtNNnnRYKMSWBDHVrykmswbdhv\nACGT*-.xX\n>s2 desc here\nAC GT\tAC\n  \nAC\n";
        assert_eq!(
            summary(&sequences(iupac, true, i32::MAX).expect("reads")),
            [
                ("s1".into(), 37, "746518b67f69d69e03f948eec75e686d".into()),
                ("s2".into(), 10, "dc410c7d30f5474a0dbea00efa96ee88".into()),
            ]
        );
        let ws: &[u8] = b">ws\n  \n AC \n\t\nA C\t \n";
        assert_eq!(
            summary(&sequences(ws, true, i32::MAX).expect("reads")),
            [("ws".into(), 6, "0a2040ad687fa7b1b504a1ed6d051e34".into())]
        );
        let crlf: &[u8] = b">c1 x\r\nACGT\r\nAC\r\n>c2\r\nGG\r\n";
        let read = sequences(crlf, true, i32::MAX).expect("reads");
        assert_eq!(
            summary(&read)[0],
            ("c1".into(), 6, "1617b7d879d437fa4c87da5875264b14".into())
        );
        assert_eq!(read[1].length, 2);
        let empty = sequences(b">e1\n>e2\nACGT\n", true, i32::MAX).expect("reads");
        assert_eq!(empty[0].length, 0);
        assert_eq!(empty[0].md5, "d41d8cd98f00b204e9800998ecf8427e");
        assert!(sequences(b"", true, i32::MAX).expect("reads").is_empty());
        assert!(
            sequences(b"\n\n", true, i32::MAX)
                .expect("reads")
                .is_empty()
        );
    }

    #[test]
    fn names_are_trimmed_and_cut_at_whitespace() {
        let names: &[u8] =
            b">  lead  desc\nA\n>tab\tdesc\nA\n>vt\x0bdesc\nA\n>trail   \nA\n>d1\nA\n>d1 again\nA\n";
        let read = sequences(names, true, i32::MAX).expect("reads");
        let names: Vec<&str> = read.iter().map(|s| s.name.as_str()).collect();
        assert_eq!(names, ["lead", "tab", "vt", "trail", "d1", "d1"]);
        let err = sequences(b">s2 desc\nA\n", false, i32::MAX).expect_err("space in a name");
        assert!(
            err.to_string()
                .contains("Sequence name 's2 desc' doesn't match regex")
        );
        let err = sequences(b">\nA\n", true, i32::MAX).expect_err("no name");
        assert!(
            err.to_string()
                .starts_with("Missing sequence name in FASTA /")
        );
        let err = sequences(b"ACGT\n>c\nA\n", true, i32::MAX).expect_err("no header");
        assert!(
            err.to_string()
                .contains("Expected > but saw chr(65) at start of sequence with index 0")
        );
        assert!(sequences(b">*x\nA\n", true, i32::MAX).is_err());
        assert!(sequences(b">x=y*\nA\n", true, i32::MAX).is_ok());
    }

    #[test]
    fn num_sequences_reads_one_past_the_last() {
        let three: &[u8] = b">a\nAC\n>b\nGG\n>c\nT\n";
        for (num, written) in [(1, 1), (0, 1), (-1, 1), (2, 2), (5, 3)] {
            let read = sequences(three, true, num).expect("reads");
            assert_eq!(read.len(), written, "{num}");
        }
        // The sequence read ahead fails on its name, not on the pattern.
        assert!(sequences(b">a\nAC\n>\nGG\n", true, 1).is_err());
        let read = sequences(b">a\nAC\n>*x\nGG\n", true, 1).expect("reads");
        assert_eq!(read.len(), 1);
    }

    #[test]
    fn gzip_and_names_outside_the_extension_list() {
        let scratch = Scratch::new("gzip");
        // `printf '>g\nACGT\n' | gzip -n`.
        let gz: &[u8] = &[
            0x1f, 0x8b, 0x08, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x03, 0xb3, 0x4b, 0xe7, 0x72,
            0x74, 0x76, 0x0f, 0xe1, 0x02, 0x00, 0x0a, 0xa3, 0x0a, 0xbd, 0x08, 0x00, 0x00, 0x00,
        ];
        let read = read_sequences(&scratch.write("g.fa.gz", gz), true, i32::MAX).expect("gzip");
        assert_eq!(read[0].md5, "f1f8f4bf413b16ad135722aa4591043e");
        let err = read_sequences(&scratch.write("g.fa", gz), true, i32::MAX).expect_err("text");
        assert!(err.to_string().contains("Expected > but saw chr(31)"));
        for name in ["g.FA", "g.fa.bgz", "g.faa"] {
            let err = read_sequences(&scratch.write(name, b">g\nA\n"), true, 1).expect_err(name);
            assert!(
                err.to_string()
                    .contains("not a supported reference file type"),
                "{name}"
            );
        }
    }

    #[test]
    fn uris_and_default_names_follow_java() {
        assert_eq!(
            file_uri(Path::new("/data//ref dir/../a+b%.fa")).expect("uri"),
            "file:///data/ref%20dir/../a+b%25.fa"
        );
        assert_eq!(
            file_uri(Path::new("/d/\u{e9}/")).expect("uri"),
            "file:///d/%C3%A9"
        );
        let cwd = std::env::current_dir().expect("cwd");
        assert_eq!(
            file_uri(Path::new("x.fa")).expect("uri"),
            file_uri(&cwd.join("x.fa")).expect("uri")
        );
        assert_eq!(
            default_dictionary(Path::new("dir/x.fa.gz")).expect("name"),
            Path::new("dir/x.dict")
        );
        assert_eq!(
            default_dictionary(Path::new("x.fasta")).expect("name"),
            Path::new("x.dict")
        );
        assert!(default_dictionary(Path::new("x.FA")).is_err());
        let args = vec!["R=dir/ref.fna".to_owned()];
        let with = with_default_output(&args).expect("args");
        assert_eq!(with.last().map(String::as_str), Some("OUTPUT=dir/ref.dict"));
        let given = vec!["R=x.fa".to_owned(), "O=y.dict".to_owned()];
        assert_eq!(with_default_output(&given).expect("args"), given);
    }

    #[test]
    fn alternative_names_keep_hash_set_order_and_picard_checks() {
        let scratch = Scratch::new("alt");
        let path = scratch.write(
            "alt.tsv",
            b"chr1\t1\nchr1\tNC_000001\nchr1\tCM000663\nchr2\t2\n\nchrM\tMT\nchrM\tchrM\nchr1\tone\nchr1\t1\n",
        );
        let names = read_alt_names(&path).expect("reads");
        // Picard 3.5.0 wrote AN:1,NC_000001,one,CM000663 for chr1.
        assert_eq!(names["chr1"], ["1", "NC_000001", "one", "CM000663"]);
        assert_eq!(names["chrM"], ["MT"]);
        let cases = [
            ("chr1 1\n", "tabulation missing in chr1 1"),
            ("\tX\n", "empty contig in \tX"),
            ("chr1\t \n", "empty alternative name in chr1\t "),
            (
                "chr1\t*x\n",
                "Illegal alternative reference sequence name in chr1\t*x",
            ),
            (
                "chr1\tchr2\nchr2\tz\n",
                "contig  chr2 previously defined as an alternate name",
            ),
            (
                "chr1\tq\nchr2\tchr1\n",
                "alternate name chr1 previously defined as a contig",
            ),
        ];
        for (text, message) in cases {
            let path = scratch.write("bad.tsv", text.as_bytes());
            let err = read_alt_names(&path).expect_err(text).to_string();
            assert!(err.contains(message), "{text:?}: {err}");
            assert!(err.starts_with("Can't read alias file "), "{err}");
        }
    }

    #[test]
    fn the_dictionary_text_orders_fields_as_picard() {
        let sequences = [Sequence {
            name: "chr1".into(),
            length: 4,
            md5: "f1f8f4bf413b16ad135722aa4591043e".into(),
        }];
        let fields = Fields {
            genome_assembly: Some("hg19".into()),
            uri: "foo:bar".into(),
            species: Some("human".into()),
            alt_names: HashMap::from([("chr1".to_owned(), vec!["1".to_owned(), "one".to_owned()])]),
        };
        let (text, rows) = dictionary(&sequences, &fields);
        assert_eq!(
            text,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\tAS:hg19\t\
             UR:foo:bar\tSP:human\tAN:1,one\n"
        );
        assert_eq!(rows[0].len(), COLUMNS.len());
        let bare = Fields {
            uri: "u".into(),
            ..Fields::default()
        };
        assert_eq!(
            dictionary(&sequences, &bare).0,
            "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\tM5:f1f8f4bf413b16ad135722aa4591043e\tUR:u\n"
        );
    }
}

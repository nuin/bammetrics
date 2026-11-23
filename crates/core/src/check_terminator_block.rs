//! `CheckTerminatorBlock`, mirroring Picard 3.5.0
//! (`picard.sam.CheckTerminatorBlock`).
//!
//! Picard asserts that `INPUT` is a readable file, asks htsjdk's
//! `BlockCompressedInputStream.checkTermination` how the file ends and
//! prints the answer's name on standard error: `HAS_TERMINATOR_BLOCK` (the
//! last 28 bytes are the BGZF end-of-file block), `HAS_HEALTHY_LAST_BLOCK`
//! (no end-of-file block, but the last block header found ends at the last
//! byte) or `DEFECTIVE` (anything else: under 28 bytes, a cut block, bytes
//! after the end-of-file block, a gzip file that is not BGZF). It exits 100
//! on `DEFECTIVE` and 0 otherwise. The check is
//! [`crate::bgzf_copy::termination`], which reads the last 64 KiB of the
//! file and inflates nothing, so the cost does not grow with the file.
//!
//! htsjdk looks for the 16 header bytes its own writer and htslib write (no
//! modification time, `XFL` 0, `OS` 255). A BGZF file from a writer that
//! sets another `OS` byte and has no end-of-file block is `DEFECTIVE` in
//! both tools.
//!
//! The Picard face ([`run_cli`]) prints the name on standard error, as
//! Picard does, and standard output stays empty. The verb and the Python
//! function return it as the tool's text and one row ([`COLUMNS`]). Picard
//! has no `OUTPUT`.
//!
//! A path that is not a regular file (a named pipe, a socket, a device) is
//! refused before it is opened: a pipe has no last bytes to seek to, and
//! opening one would wait for a writer. Picard opens it and reports a size
//! of 0 as `DEFECTIVE` (KNOWN_DIFFERENCES.md).

use std::fs::File;
use std::path::Path;

use crate::bgzf_copy::{Termination, termination};
use crate::compare_metrics::JavaFile;
use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

/// Picard's exit code for a `DEFECTIVE` file.
pub const DEFECTIVE_STATUS: u8 = 100;

/// The columns of the one row returned to Python and printed by `--format`.
pub const COLUMNS: [&str; 2] = ["file", "termination"];

pub const TOOL: Tool = Tool {
    name: "CheckTerminatorBlock",
    summary: "one row (file, termination) saying how a BGZF file ends: HAS_TERMINATOR_BLOCK, \
              HAS_HEALTHY_LAST_BLOCK or DEFECTIVE, from its last 64 KiB; the command prints the \
              name and exits 100 on DEFECTIVE.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The BGZF file whose end is checked: a BAM, a bgzip VCF, any bgzip file.",
        )
        .required(),
        Opt::rejected(&["OUTPUT", "O"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
        crate::registry::STRINGENCY_NO_RECORDS,
        crate::registry::REFERENCE_NOT_DECODED,
    ],
    run: run_tool,
    collector: None,
};

/// htsjdk's `FileTermination.name()`.
pub fn name(termination: Termination) -> &'static str {
    match termination {
        Termination::Terminator => "HAS_TERMINATOR_BLOCK",
        Termination::HealthyLastBlock => "HAS_HEALTHY_LAST_BLOCK",
        Termination::Defective => "DEFECTIVE",
    }
}

/// `IOUtil.assertFileIsReadable` with htsjdk's messages, then the refusal
/// of a path that is not a regular file.
fn open(input: &Path) -> Result<(File, u64)> {
    let file = JavaFile::new(&input.to_string_lossy());
    let path = Path::new(&file.absolute);
    let Ok(metadata) = std::fs::metadata(path) else {
        return Err(Error::Metric(format!(
            "Cannot read non-existent file: {}",
            file.uri()
        )));
    };
    if metadata.is_dir() {
        return Err(Error::Metric(format!(
            "Cannot read file because it is a directory: {}",
            file.uri()
        )));
    }
    if !metadata.is_file() {
        return Err(Error::Usage(format!(
            "INPUT {} is not a regular file: CheckTerminatorBlock reads the last bytes of a \
             file, and a pipe or a device has none to seek to",
            input.display()
        )));
    }
    let opened = File::open(path)
        .map_err(|_| Error::Metric(format!("File exists but is not readable: {}", file.uri())))?;
    Ok((opened, metadata.len()))
}

/// How the file at `input` ends.
pub fn check_file(input: &Path) -> Result<Termination> {
    let (mut file, size) = open(input)?;
    termination(&mut file, size).map_err(|source| Error::Input {
        path: input.to_path_buf(),
        source,
    })
}

/// One check: the answer, and the tool's output (the name as text, one
/// row).
#[derive(Clone, Debug)]
pub struct Check {
    pub termination: Termination,
    pub output: ToolOutput,
}

impl Check {
    /// Picard's exit code: 100 for `DEFECTIVE`, else 0.
    pub fn status(&self) -> u8 {
        match self.termination {
            Termination::Defective => DEFECTIVE_STATUS,
            Termination::Terminator | Termination::HealthyLastBlock => 0,
        }
    }
}

/// Check `INPUT` of `inv`.
pub fn check(inv: &Invocation) -> Result<Check> {
    let termination = check_file(&inv.input)?;
    let name = name(termination);
    let row = vec![
        Value::Str(inv.input.to_string_lossy().into_owned()),
        Value::Str(name.to_owned()),
    ];
    Ok(Check {
        termination,
        output: MetricsOutput::text(&COLUMNS, vec![row], format!("{name}\n")).into(),
    })
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    Ok(check(inv)?.output)
}

/// The Picard face: the name on standard error, nothing on standard
/// output, Picard's exit code.
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(args)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    let check = check(&inv)?;
    eprintln!("{}", name(check.termination));
    Ok(check.status())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bgzf_copy::EOF_BLOCK;
    use std::io::Write;

    /// A file of `bytes` in a directory of its own.
    fn file_of(tag: &str, bytes: &[u8]) -> (std::path::PathBuf, std::path::PathBuf) {
        use std::os::unix::fs::DirBuilderExt;
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos());
        let dir = std::env::temp_dir().join(format!(
            "bammetrics-terminator-{tag}-{}-{nanos}",
            std::process::id()
        ));
        std::fs::DirBuilder::new()
            .mode(0o700)
            .create(&dir)
            .expect("a new temp dir");
        let path = dir.join("x.bam");
        let mut file = File::create_new(&path).expect("a new file");
        file.write_all(bytes).expect("write");
        (dir, path)
    }

    fn termination_of(tag: &str, bytes: &[u8]) -> Termination {
        let (dir, path) = file_of(tag, bytes);
        let found = check_file(&path);
        std::fs::remove_dir_all(&dir).expect("remove the temp dir");
        found.expect("a regular file is checked")
    }

    /// A stored BGZF block holding `data`, as htsjdk frames one.
    fn block(data: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        crate::bgzf_copy::write_blocks(&mut out, data, 0).expect("write a block");
        out
    }

    #[test]
    fn the_three_outcomes() {
        let data = block(b"some bytes of a record");
        let whole = [data.as_slice(), &EOF_BLOCK].concat();
        assert_eq!(termination_of("whole", &whole), Termination::Terminator);
        assert_eq!(
            termination_of("no-eof", &data),
            Termination::HealthyLastBlock
        );
        for cut in [1, 8, 27, 29, 40] {
            let short = &whole[..whole.len() - cut];
            assert_eq!(
                termination_of("cut", short),
                Termination::Defective,
                "cut {cut}"
            );
        }
        assert_eq!(termination_of("empty", b""), Termination::Defective);
        assert_eq!(
            termination_of("short", &EOF_BLOCK[..27]),
            Termination::Defective
        );
        assert_eq!(
            termination_of("eof-alone", &EOF_BLOCK),
            Termination::Terminator
        );
        let twice = [EOF_BLOCK.as_slice(), &EOF_BLOCK].concat();
        assert_eq!(termination_of("eof-twice", &twice), Termination::Terminator);
        let trailing = [whole.as_slice(), b"xyz"].concat();
        assert_eq!(
            termination_of("trailing", &trailing),
            Termination::Defective
        );
    }

    #[test]
    fn names_and_exit_codes_are_picards() {
        assert_eq!(name(Termination::Terminator), "HAS_TERMINATOR_BLOCK");
        assert_eq!(
            name(Termination::HealthyLastBlock),
            "HAS_HEALTHY_LAST_BLOCK"
        );
        assert_eq!(name(Termination::Defective), "DEFECTIVE");
        let (dir, path) = file_of("status", &EOF_BLOCK[..10]);
        let args = [format!("INPUT={}", path.display())];
        let inv = TOOL.parse(&args).expect("parse");
        let check = check(&inv).expect("check");
        assert_eq!(check.status(), 100);
        let main = check.output.main.expect("the text");
        assert_eq!(main.text.as_deref(), Some("DEFECTIVE\n"));
        assert_eq!(main.rows.len(), 1);
        std::fs::remove_dir_all(&dir).expect("remove the temp dir");
    }

    #[test]
    fn a_missing_file_a_directory_and_a_pipe_are_errors() {
        let (dir, path) = file_of("errors", b"");
        let missing = check_file(&dir.join("none.bam")).expect_err("missing");
        assert!(
            missing
                .to_string()
                .starts_with("Cannot read non-existent file: file:///")
        );
        let directory = check_file(&dir).expect_err("a directory");
        assert!(
            directory
                .to_string()
                .starts_with("Cannot read file because it is a directory: file:///")
        );
        let device = check_file(Path::new("/dev/null")).expect_err("a device");
        assert!(matches!(device, Error::Usage(_)), "{device}");
        assert!(path.exists());
        std::fs::remove_dir_all(&dir).expect("remove the temp dir");
    }
}

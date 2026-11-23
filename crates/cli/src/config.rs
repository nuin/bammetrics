//! The user config file (docs/cli.md, "Config file"):
//! `$XDG_CONFIG_HOME/bammetrics/config.toml`, else
//! `$HOME/.config/bammetrics/config.toml`.
//!
//! ```toml
//! reference = "~/refs/hs37d5.fa"
//! threads = 8
//! tmp_dir = "/scratch"
//! sort_memory = "4G"
//! ```
//!
//! The file is read through [`bammetrics_core::toml::parse_tables`] with
//! a cap of [`FILE_CAP`] bytes; every key is checked, an unknown one with
//! the closest name. `~/` at the start of a path is the home directory,
//! and a relative path is relative to the directory of the file.
//! `BAMMETRICS_REFERENCE` replaces the file's `reference`, and `-R`
//! replaces both. `--no-config` leaves the file unread (the variable still
//! counts). The Picard face (`KEY=VALUE` after a Picard class) reads
//! neither: it stays Picard's command line.

use std::path::{Path, PathBuf};

use bammetrics_core::pass;
use bammetrics_core::registry::did_you_mean;
use bammetrics_core::toml::{Entry, Scalar, parse_tables};

use crate::error::VerbError;

/// The verb that prints the effective values.
pub const VERB: &str = "config";
/// The most bytes of a config file: 64 KiB. A file of the four keys holds
/// about 100 bytes.
pub const FILE_CAP: u64 = 64 << 10;
/// The variable that names the reference in place of the file's.
pub const REFERENCE_ENV: &str = "BAMMETRICS_REFERENCE";
/// The option that leaves the file unread.
pub const NO_CONFIG: &str = "--no-config";
/// The keys of the file.
pub const KEYS: [&str; 4] = ["reference", "threads", "tmp_dir", "sort_memory"];

/// Where a value came from.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Origin {
    /// The config file, at a line.
    File(PathBuf, usize),
    /// [`REFERENCE_ENV`].
    Env,
}

impl Origin {
    pub fn describe(&self) -> String {
        match self {
            Origin::File(path, line) => format!("{}:{line}", path.display()),
            Origin::Env => REFERENCE_ENV.to_owned(),
        }
    }
}

/// One value with its origin.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Sourced<T> {
    pub value: T,
    pub origin: Origin,
}

/// The effective values: the file's, with the variable over them.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Config {
    /// The file looked for; `None` without `XDG_CONFIG_HOME` and `HOME`.
    pub path: Option<PathBuf>,
    /// The file was there and read.
    pub read: bool,
    /// `--no-config`: the file was left unread.
    pub skipped: bool,
    pub reference: Option<Sourced<String>>,
    pub threads: Option<Sourced<usize>>,
    pub tmp_dir: Option<Sourced<String>>,
    pub sort_memory: Option<Sourced<String>>,
}

impl Config {
    /// The reference path, if any.
    pub fn reference(&self) -> Option<&str> {
        self.reference.as_ref().map(|r| r.value.as_str())
    }
}

/// The directory of the user's bammetrics files: `$XDG_CONFIG_HOME/bammetrics`
/// when the variable holds an absolute path, else `$HOME/.config/bammetrics`.
pub fn dir() -> Option<PathBuf> {
    let xdg = std::env::var_os("XDG_CONFIG_HOME")
        .map(PathBuf::from)
        .filter(|path| path.is_absolute());
    let base = xdg.or_else(|| {
        std::env::var_os("HOME")
            .filter(|home| !home.is_empty())
            .map(|home| PathBuf::from(home).join(".config"))
    })?;
    Some(base.join("bammetrics"))
}

/// The config file's path.
pub fn path() -> Option<PathBuf> {
    dir().map(|dir| dir.join("config.toml"))
}

/// `~/` at the start of `text` as the home directory, a relative path
/// under `base`.
fn expand(text: &str, base: &Path) -> String {
    if let Some(rest) = text.strip_prefix("~/")
        && let Some(home) = std::env::var_os("HOME")
    {
        return PathBuf::from(home).join(rest).display().to_string();
    }
    let path = Path::new(text);
    if path.is_absolute() || text.is_empty() {
        return text.to_owned();
    }
    base.join(path).display().to_string()
}

fn error(file: &Path, line: usize, what: impl AsRef<str>) -> VerbError {
    VerbError::usage(
        VERB,
        format!("{}:{line}: {}", file.display(), what.as_ref()),
    )
}

fn string_value(file: &Path, entry: &Entry) -> Result<String, VerbError> {
    match &entry.value {
        Scalar::Str(text) if !text.is_empty() => Ok(text.clone()),
        Scalar::Str(_) => Err(error(file, entry.line, format!("{} is empty", entry.key))),
        other => Err(error(
            file,
            entry.line,
            format!("{}: expected a string, got {}", entry.key, other.kind()),
        )),
    }
}

/// The file's text into `config`. The text holds keys without a table, so
/// the reader sees it under one `[config]` line in front: every line
/// number moves by one, and the errors move it back.
pub fn parse(text: &str, file: &Path, config: &mut Config) -> Result<(), VerbError> {
    let wrapped = format!(
        "[config]\n{}",
        text.strip_prefix('\u{feff}').unwrap_or(text)
    );
    let tables = parse_tables(&wrapped)
        .map_err(|e| error(file, e.line.saturating_sub(1).max(1), e.message))?;
    let base = file.parent().unwrap_or(Path::new("."));
    let mut seen: Vec<&str> = Vec::new();
    for table in &tables {
        if table.name != "config" || table.array || table.line != 1 {
            return Err(error(
                file,
                table.line - 1,
                format!("a config file holds keys alone, no table [{}]", table.name),
            )
            .hint(format!("keys: {}", KEYS.join(", "))));
        }
        for entry in &table.entries {
            let line = entry.line - 1;
            let at = || Origin::File(file.to_path_buf(), line);
            if seen.contains(&entry.key.as_str()) {
                return Err(error(file, line, format!("{} is given twice", entry.key)));
            }
            seen.push(&entry.key);
            let entry = Entry {
                line,
                ..entry.clone()
            };
            match entry.key.as_str() {
                "reference" => {
                    config.reference = Some(Sourced {
                        value: expand(&string_value(file, &entry)?, base),
                        origin: at(),
                    });
                }
                "tmp_dir" => {
                    config.tmp_dir = Some(Sourced {
                        value: expand(&string_value(file, &entry)?, base),
                        origin: at(),
                    });
                }
                "sort_memory" => {
                    let value = match &entry.value {
                        Scalar::Int(n) if *n > 0 => n.to_string(),
                        Scalar::Str(text) => text.clone(),
                        _ => String::new(),
                    };
                    if crate::chain::bytes_or_why(&value).is_err() {
                        return Err(error(
                            file,
                            line,
                            "sort_memory: expected bytes above 0, a number or a string with \
                             K, M, G or T (\"4G\")",
                        ));
                    }
                    config.sort_memory = Some(Sourced {
                        value,
                        origin: at(),
                    });
                }
                "threads" => {
                    let threads = match &entry.value {
                        Scalar::Int(n) => usize::try_from(*n).ok(),
                        _ => None,
                    }
                    .ok_or_else(|| error(file, line, "threads: expected a whole number"))?;
                    let threads = pass::check_threads(threads)
                        .map_err(|e| error(file, line, format!("threads: {threads} is {e}")))?;
                    config.threads = Some(Sourced {
                        value: threads,
                        origin: at(),
                    });
                }
                other => {
                    let mut e = error(file, line, format!("unknown key {other}"));
                    if let Some(best) = did_you_mean(other, KEYS) {
                        e = e.hint(format!("did you mean {best}?"));
                    }
                    return Err(e.hint(format!("keys: {}", KEYS.join(", "))));
                }
            }
        }
    }
    Ok(())
}

/// The effective values: the file unless `use_file` is false, then
/// [`REFERENCE_ENV`] over its reference. A file that is not there is no
/// error; one that cannot be read, is not UTF-8, holds more than
/// [`FILE_CAP`] bytes or does not parse is a usage error naming it.
pub fn load(use_file: bool) -> Result<Config, VerbError> {
    let mut config = Config {
        path: path(),
        skipped: !use_file,
        ..Config::default()
    };
    if let (true, Some(file)) = (use_file, config.path.clone())
        && file.exists()
    {
        let text = read_regular_capped(&file, FILE_CAP).map_err(|e| {
            VerbError::usage(VERB, format!("{}: not a config file: {e}", file.display()))
                .hint(format!("{NO_CONFIG} leaves it unread"))
        })?;
        parse(&text, &file, &mut config)?;
        config.read = true;
    }
    if let Some(value) = std::env::var_os(REFERENCE_ENV).filter(|v| !v.is_empty()) {
        config.reference = Some(Sourced {
            value: value.to_string_lossy().into_owned(),
            origin: Origin::Env,
        });
    }
    Ok(config)
}

/// The text of the user file at `path` (a config file or a recipe), at
/// most `cap` bytes. The file is opened without blocking and must be a
/// regular file: a FIFO, a device or a directory is refused before a byte
/// is read, so a pipe named `config.toml` never stops a command.
pub(crate) fn read_regular_capped(path: &Path, cap: u64) -> std::io::Result<String> {
    use std::io::{Error, ErrorKind, Read as _};
    use std::os::unix::fs::OpenOptionsExt as _;
    let file = std::fs::OpenOptions::new()
        .read(true)
        .custom_flags(libc::O_NONBLOCK)
        .open(path)?;
    if !file.metadata()?.is_file() {
        return Err(Error::new(ErrorKind::InvalidInput, "not a regular file"));
    }
    let mut bytes = Vec::new();
    file.take(cap.saturating_add(1)).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > cap {
        return Err(Error::new(
            ErrorKind::FileTooLarge,
            format!("larger than {cap} bytes"),
        ));
    }
    String::from_utf8(bytes).map_err(|e| Error::new(ErrorKind::InvalidData, e.utf8_error()))
}

fn renamed(verb: &str, mut error: VerbError) -> VerbError {
    verb.clone_into(&mut error.verb);
    error
}

/// The reference of the config file or [`REFERENCE_ENV`] for `verb`, which
/// was given none, checked by [`checked`]: for a command that reads it.
pub fn reference_for(verb: &str, use_file: bool) -> Result<Option<String>, VerbError> {
    let config = load(use_file).map_err(|e| renamed(verb, e))?;
    config
        .reference
        .map(|reference| checked(verb, &reference))
        .transpose()
}

/// The reference of the config file or [`REFERENCE_ENV`] for `verb`, which
/// was given none and may run without one (`qc`, `plan`, `collect`): a
/// path that is not a file is left out with a note on standard error, so
/// a stale config fails only the steps that need a reference.
pub fn optional_reference_for(verb: &str, use_file: bool) -> Result<Option<String>, VerbError> {
    let config = load(use_file).map_err(|e| renamed(verb, e))?;
    let Some(reference) = config.reference else {
        return Ok(None);
    };
    match checked(verb, &reference) {
        Ok(path) => Ok(Some(path)),
        Err(error) => {
            eprintln!("bammetrics {verb}: note: {}: not used", error.what);
            Ok(None)
        }
    }
}

/// The path of a reference from the config file or [`REFERENCE_ENV`]: an
/// input error (exit 1) naming where it came from when it is not a file.
pub fn checked(verb: &str, reference: &Sourced<String>) -> Result<String, VerbError> {
    if Path::new(&reference.value).is_file() {
        return Ok(reference.value.clone());
    }
    Err(VerbError::input(
        verb,
        format!(
            "{} not found (the reference of {})",
            reference.value,
            reference.origin.describe()
        ),
    )
    .hint(match reference.origin {
        Origin::File(..) => format!("-R FASTA names another; {NO_CONFIG} leaves the file unread"),
        Origin::Env => format!("-R FASTA names another; unset {REFERENCE_ENV} to use none"),
    }))
}

/// `args` without every [`NO_CONFIG`] in front of a `--`, and whether one
/// was there.
pub fn strip_no_config(args: &[String]) -> (Vec<String>, bool) {
    let end = args.iter().position(|a| a == "--").unwrap_or(args.len());
    let found = args[..end].iter().any(|a| a == NO_CONFIG);
    let kept = args
        .iter()
        .enumerate()
        .filter(|(i, a)| *i >= end || *a != NO_CONFIG)
        .map(|(_, a)| a.clone())
        .collect();
    (kept, found)
}

/// `bammetrics config [--no-config]`: the effective values, one line each
/// with where it came from, then the file.
pub fn run(args: &[String]) -> Result<u8, VerbError> {
    let (rest, no_config) = strip_no_config(args);
    if let Some(first) = rest.first() {
        if first == "--help" || first == "-h" {
            crate::events::print_stdout(&help());
            return Ok(0);
        }
        return Err(
            VerbError::usage(VERB, format!("unexpected argument {first}"))
                .from(first.clone())
                .hint("usage: bammetrics config [--no-config]"),
        );
    }
    let config = load(!no_config)?;
    crate::events::print_stdout(&render(&config));
    Ok(0)
}

fn help() -> String {
    format!(
        "bammetrics config: the values chain, qc, plan, collect and the tool verbs take from \
         the user config file and {REFERENCE_ENV}, and where each came from.\n\n\
         Usage: bammetrics config [--no-config]\n\n\
         The file: $XDG_CONFIG_HOME/bammetrics/config.toml, else \
         $HOME/.config/bammetrics/config.toml. Keys: reference, threads, tmp_dir, \
         sort_memory. {REFERENCE_ENV} replaces the file's reference; -R replaces both. \
         --no-config leaves the file unread. The Picard face (KEY=VALUE) reads neither.\n"
    )
}

/// The text `bammetrics config` prints.
pub fn render(config: &Config) -> String {
    let mut text = String::new();
    let file = config.path.as_ref().map_or_else(
        || "none (no XDG_CONFIG_HOME, no HOME)".to_owned(),
        |p| p.display().to_string(),
    );
    let state = if config.skipped {
        "not read (--no-config)"
    } else if config.read {
        "read"
    } else {
        "not there"
    };
    text.push_str(&format!("file         {file} ({state})\n"));
    fn line<T: std::fmt::Display>(
        text: &mut String,
        key: &str,
        value: Option<&Sourced<T>>,
        default: &str,
    ) {
        match value {
            Some(v) => text.push_str(&format!(
                "{key:<12} {} (from {})\n",
                v.value,
                v.origin.describe()
            )),
            None => text.push_str(&format!("{key:<12} {default}\n")),
        }
    }
    line(
        &mut text,
        "reference",
        config.reference.as_ref(),
        "none (give -R)",
    );
    line(
        &mut text,
        "threads",
        config.threads.as_ref(),
        "default (the available cores, at most 8)",
    );
    line(
        &mut text,
        "tmp_dir",
        config.tmp_dir.as_ref(),
        "default (the system temporary directory)",
    );
    line(
        &mut text,
        "sort_memory",
        config.sort_memory.as_ref(),
        "default (2G)",
    );
    let recipes = crate::chain::recipes::user_dir()
        .map_or_else(|| "none".to_owned(), |d| d.display().to_string());
    text.push_str(&format!("recipes      {recipes}\n"));
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed(text: &str) -> Result<Config, String> {
        let mut config = Config::default();
        parse(text, Path::new("/etc/b/config.toml"), &mut config)
            .map(|()| config)
            .map_err(|e| e.to_string())
    }

    #[test]
    fn the_four_keys_are_read_with_their_lines() {
        let config = parsed(
            "# mine\nreference = \"ref/hs37d5.fa\"\nthreads = 4\ntmp_dir = \"/scratch\"\nsort_memory = \"4G\"\n",
        )
        .expect("parses");
        let reference = config.reference.expect("reference");
        assert_eq!(reference.value, "/etc/b/ref/hs37d5.fa");
        assert_eq!(
            reference.origin,
            Origin::File(PathBuf::from("/etc/b/config.toml"), 2)
        );
        assert_eq!(config.threads.expect("threads").value, 4);
        assert_eq!(config.tmp_dir.expect("tmp").value, "/scratch");
        assert_eq!(config.sort_memory.expect("memory").value, "4G");
        assert_eq!(
            parsed("sort_memory = 1024\n")
                .expect("bytes")
                .sort_memory
                .expect("m")
                .value,
            "1024"
        );
    }

    #[test]
    fn malformed_files_are_errors_with_the_line() {
        let error = |text: &str| parsed(text).expect_err(text);
        assert!(
            error("refernce = \"x\"\n")
                .contains("/etc/b/config.toml:1: unknown key refernce\n  did you mean reference?")
        );
        assert!(error("\nthreads = \"8\"\n").contains(":2: threads: expected a whole number"));
        assert!(error("threads = 9999\n").contains(":1: threads: 9999 is more than the maximum"));
        assert!(error("sort_memory = \"lots\"\n").contains(":1: sort_memory: expected bytes"));
        assert!(error("reference = 3\n").contains("reference: expected a string, got an integer"));
        assert!(error("reference = \"\"\n").contains("reference is empty"));
        assert!(error("[chain]\nthreads = 2\n").contains(":1: a config file holds keys alone"));
        assert!(error("threads = 2\nthreads = 3\n").contains(":2: threads is given twice"));
        assert!(error("threads = [1, 2\n").contains(":1:"));
    }

    #[test]
    fn no_config_is_taken_out_in_front_of_a_double_dash() {
        let args: Vec<String> = ["a", "--no-config", "b", "--", "--no-config"]
            .iter()
            .map(|s| (*s).to_owned())
            .collect();
        let (kept, found) = strip_no_config(&args);
        assert!(found);
        assert_eq!(kept, ["a", "b", "--", "--no-config"]);
    }
}

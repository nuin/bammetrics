//! The error of the `report` and `check` verbs, with the exit code of
//! docs/05 section 2.4: 2 for a usage error, 1 for everything else. The
//! Picard face keeps its own [`bammetrics_core::Error`], printed as before
//! and exiting 1. The verbs of the grammar print [`VerbError`].

use std::fmt;
use std::path::PathBuf;

/// An error on a verb in one style: `bammetrics <verb>: <what> (from
/// <option>)`, then at most two indented hint lines.
#[derive(Debug)]
pub struct VerbError {
    pub verb: String,
    pub what: String,
    /// The option or argument the error comes from, as the user spelled it.
    pub from: Option<String>,
    pub hints: Vec<String>,
    pub code: u8,
}

impl VerbError {
    /// A usage error (exit 2): nothing ran.
    pub fn usage(verb: &str, what: impl Into<String>) -> Self {
        Self {
            verb: verb.to_owned(),
            what: what.into(),
            from: None,
            hints: Vec::new(),
            code: 2,
        }
    }

    /// An input error (exit 1).
    pub fn input(verb: &str, what: impl Into<String>) -> Self {
        Self {
            code: 1,
            ..Self::usage(verb, what)
        }
    }

    pub fn from(mut self, option: impl Into<String>) -> Self {
        self.from = Some(option.into());
        self
    }

    pub fn hint(mut self, hint: impl Into<String>) -> Self {
        if self.hints.len() < 2 {
            self.hints.push(hint.into());
        }
        self
    }

    /// A core error raised while a verb ran: usage errors exit 2, the
    /// rest 1.
    pub fn core(verb: &str, e: &bammetrics_core::Error) -> Self {
        match e {
            bammetrics_core::Error::Usage(message) => Self::usage(verb, message.clone()),
            other => Self::input(verb, other.to_string()),
        }
    }
}

impl fmt::Display for VerbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // An empty verb is an error about the first argument itself.
        if self.verb.is_empty() {
            write!(f, "bammetrics: {}", self.what)?;
        } else {
            write!(f, "bammetrics {}: {}", self.verb, self.what)?;
        }
        if let Some(from) = &self.from {
            write!(f, " (from {from})")?;
        }
        for hint in &self.hints {
            write!(f, "\n  {hint}")?;
        }
        Ok(())
    }
}

impl std::error::Error for VerbError {}

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Core(#[from] bammetrics_core::Error),
    #[error(transparent)]
    Verdict(#[from] bammetrics_verdict::Error),
    /// A file or directory could not be read or written.
    #[error("{path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    /// The arguments cannot be acted on; nothing ran.
    #[error("{0}")]
    Usage(String),
    /// Something about the inputs stops the verb (no metrics file, two
    /// samples with one name).
    #[error("{0}")]
    Input(String),
    /// An error of a grammar verb, printed whole.
    #[error(transparent)]
    Verb(#[from] VerbError),
}

impl Error {
    pub fn io(path: &std::path::Path, source: std::io::Error) -> Self {
        Error::Io {
            path: path.to_path_buf(),
            source,
        }
    }

    pub fn exit_code(&self) -> u8 {
        match self {
            Error::Usage(_) => 2,
            Error::Verb(e) => e.code,
            _ => 1,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

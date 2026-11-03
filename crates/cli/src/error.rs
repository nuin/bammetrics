//! The error of the `report` and `check` verbs, with the exit code of
//! docs/05 section 2.4: 2 for a usage error, 1 for everything else. The
//! Picard face and the five older subcommands keep their own
//! [`bammetrics_core::Error`], printed as before and exiting 1.

use std::path::PathBuf;

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
            _ => 1,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

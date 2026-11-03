//! Crate error type.

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// A preset file cannot be used: the syntax, a rule or the bounds.
    #[error("{file}:{}{message}", line.map(|l| format!("{l}: ")).unwrap_or_default())]
    Preset {
        /// The file name or `built-in`.
        file: String,
        line: Option<usize>,
        message: String,
    },
    /// No shipped preset has this name and no file exists at this path.
    #[error("no preset {spec:?}: not a shipped preset ({}) and not a file", crate::preset::names().join(", "))]
    UnknownPreset { spec: String },
    /// A preset or metrics file cannot be read.
    #[error("{path}: {source}")]
    Input {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },
    /// A metrics file does not parse.
    #[error(transparent)]
    Metrics(bammetrics_core::Error),
    /// A verdict.json does not parse or lacks a field.
    #[error("verdict.json: {0}")]
    Json(String),
}

impl From<bammetrics_core::Error> for Error {
    fn from(error: bammetrics_core::Error) -> Self {
        Error::Metrics(error)
    }
}

pub type Result<T> = std::result::Result<T, Error>;

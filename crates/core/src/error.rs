//! Crate error type.

use std::io;
use std::path::PathBuf;

/// Everything that can go wrong while computing metrics.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    /// Opening or reading an input file failed.
    #[error("{path}: {source}")]
    Input {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    /// Writing an output file failed.
    #[error("{path}: {source}")]
    Output {
        path: PathBuf,
        #[source]
        source: io::Error,
    },
    /// A record could not be decoded.
    #[error("bad record: {0}")]
    Record(#[from] io::Error),
    /// The reference FASTA cannot be used with this BAM.
    #[error("{path}: {message}")]
    Reference { path: PathBuf, message: String },
    /// An interval list or BED file cannot be used with this BAM.
    #[error("{path}: {message}")]
    Intervals { path: PathBuf, message: String },
    /// The data cannot produce the metric (Picard throws in the same case).
    #[error("{0}")]
    Metric(String),
    /// A command-line option cannot be honoured.
    #[error("{0}")]
    Usage(String),
    /// A record carries a flag combination that needs a field it does not have.
    #[error("record {name}: {what}")]
    MalformedRecord { name: String, what: &'static str },
}

pub type Result<T> = std::result::Result<T, Error>;

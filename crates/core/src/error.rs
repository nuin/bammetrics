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
    /// A record could not be decoded.
    #[error("bad record: {0}")]
    Record(#[from] io::Error),
    /// A record carries a flag combination that needs a field it does not have.
    #[error("record {name}: {what}")]
    MalformedRecord { name: String, what: &'static str },
}

pub type Result<T> = std::result::Result<T, Error>;

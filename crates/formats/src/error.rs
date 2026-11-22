//! The crate's error: text that is not JSON, JSON that is not a
//! `bammetrics.metrics/1` document, text that is not a metrics file, or a
//! document Picard text cannot hold.

use std::fmt;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The text is not JSON; the message names the byte.
    Json(String),
    /// The JSON is not a metrics document.
    Document(String),
    /// The text is not a Picard metrics file.
    Metrics(String),
    /// The document holds a cell Picard text has no form for.
    Picard(String),
    /// A Parquet file that cannot be written or read, or a build without
    /// the feature `parquet`.
    Parquet(String),
    /// A file that cannot be read; the message names it.
    Io(String),
    /// The document is over a limit of the format (an XLSX sheet's rows).
    Limit(String),
    /// A binary format asked of the text writer.
    Binary(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Json(message) => write!(f, "json: {message}"),
            Error::Parquet(message) => write!(f, "parquet: {message}"),
            Error::Document(message)
            | Error::Metrics(message)
            | Error::Picard(message)
            | Error::Io(message)
            | Error::Limit(message) => f.write_str(message),
            Error::Binary(name) => write!(f, "{name} is a binary format, not text"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

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
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Json(message) => write!(f, "json: {message}"),
            Error::Document(message) | Error::Metrics(message) | Error::Picard(message) => {
                f.write_str(message)
            }
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

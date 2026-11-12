//! Crate error type.

use std::io;
use std::path::{Path, PathBuf};

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

impl Error {
    /// The error with `shown` in place of `copy` when it names `copy`: as
    /// the input that failed, or in its message (a bad line of SAM text
    /// names its file). A pipe copied into the temporary directory is
    /// reported under the name the user gave, never under the copy's.
    #[must_use]
    pub fn naming_input(self, copy: &Path, shown: &Path) -> Self {
        let (from, to) = (copy.display().to_string(), shown.display().to_string());
        match self {
            Error::Input { path, source } if path == copy => Error::Input {
                path: shown.to_path_buf(),
                source,
            },
            Error::Record(source) if source.to_string().contains(&from) => Error::Record(
                io::Error::new(source.kind(), source.to_string().replace(&from, &to)),
            ),
            Error::Metric(text) if text.contains(&from) => Error::Metric(text.replace(&from, &to)),
            other => other,
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_input_error_of_a_copy_names_the_file_the_user_gave() {
        let copy = Path::new("/tmp/bammetrics-markdup.1.0/input.0.tmp");
        let failed = |path: &str| Error::Input {
            path: PathBuf::from(path),
            source: io::Error::other("not a BAM or CRAM file"),
        };
        let shown = failed("/tmp/bammetrics-markdup.1.0/input.0.tmp")
            .naming_input(copy, Path::new("/dev/stdin"));
        assert_eq!(shown.to_string(), "/dev/stdin: not a BAM or CRAM file");
        let other = failed("/data/other.bam").naming_input(copy, Path::new("/dev/stdin"));
        assert_eq!(other.to_string(), "/data/other.bam: not a BAM or CRAM file");
        let usage = Error::Usage("x".to_owned()).naming_input(copy, Path::new("/dev/stdin"));
        assert_eq!(usage.to_string(), "x");
        // A bad line of SAM text names its file in the message.
        let line = Error::Record(io::Error::new(
            io::ErrorKind::InvalidData,
            "Error parsing text SAM file. MAPQ should be < 256; File \
             /tmp/bammetrics-markdup.1.0/input.0.tmp; Line 2",
        ))
        .naming_input(copy, Path::new("/dev/stdin"));
        assert_eq!(
            line.to_string(),
            "bad record: Error parsing text SAM file. MAPQ should be < 256; File /dev/stdin; Line 2"
        );
    }
}

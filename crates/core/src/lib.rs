//! Picard-concordant BAM QC metrics.

pub mod error;
pub mod histogram;
pub mod insert_size;
pub mod metrics_file;
pub mod picard_cli;

pub use error::{Error, Result};

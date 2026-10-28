//! Picard-concordant BAM QC metrics.

pub mod error;
pub mod histogram;
pub mod insert_size;
pub mod metrics_file;

pub use error::{Error, Result};

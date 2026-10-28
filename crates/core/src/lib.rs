//! Picard-concordant BAM QC metrics.

pub mod adapter;
pub mod alignment_summary;
pub mod duplication;
pub mod error;
pub mod histogram;
pub mod hs_metrics;
pub mod insert_size;
pub mod intervals;
pub mod metrics_file;
pub mod picard_cli;
pub mod quality_yield;
pub mod record;
pub mod reference;
pub mod sensitivity;
pub mod sequence;

pub use error::{Error, Result};

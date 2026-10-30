//! Picard-concordant BAM QC metrics.

pub mod adapter;
pub mod alignment_summary;
pub mod blocks;
pub mod collect;
pub mod duplication;
pub mod error;
pub mod filters;
pub mod histogram;
pub mod hs_metrics;
pub mod insert_size;
pub mod intervals;
pub mod levels;
pub mod metric;
pub mod metrics_file;
pub mod pass;
pub mod picard_cli;
pub mod quality_yield;
pub mod record;
pub mod reference;
pub mod registry;
pub mod sensitivity;
pub mod sequence;

pub use error::{Error, Result};

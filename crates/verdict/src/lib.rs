//! Verdicts: a preset of threshold rules evaluated over one sample's Picard
//! metrics files, written as `verdict.json` (docs/05 section 4.3, docs/verdict.md).
//!
//! Rules of this crate, from docs/05 section 4.9: it reads the metrics files
//! it is given and nothing else. No other run directory, no history, no user
//! identity, no network, no comment field. A feature that needs one of those
//! belongs to the hosted dashboard.

pub mod error;
pub mod json;
pub mod preset;
pub mod table;
pub mod verdict;

pub use error::{Error, Result};
pub use preset::{Bound, Comparator, Preset, Rule, content_hash, names};
pub use verdict::{
    Band, Check, Counts, FileInfo, Input, PresetInfo, SCHEMA, Status, Verdict, evaluate,
    evaluate_files,
};

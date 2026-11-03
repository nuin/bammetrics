//! The verdict-first per-sample report (docs/05 section 4): one
//! self-contained HTML page from a [`Verdict`] and the parsed metrics files
//! it was evaluated over. No script, no network request, no dependency
//! beyond the core and verdict crates; HTML by `write!` into a `String`,
//! CSS by `include_str!`, the four plots as inline SVG.
//!
//! The rules of docs/05 section 4.9 apply as in the verdict crate: the page
//! is built from the inputs it is given and nothing else.

pub mod html;
pub mod page;
pub mod plots;
pub mod svg;
pub mod theme;

pub use bammetrics_verdict::{Input, Verdict};
pub use html::escape;
pub use page::{Options, render, write};
pub use plots::{Family, figure};

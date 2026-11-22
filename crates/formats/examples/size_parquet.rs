//! Binary-growth probe, with Parquet: `size_formats` built with the feature
//! `parquet`, so it also writes each file as Parquet and reads it back. The
//! size difference of this example over `size_formats` is the Parquet growth
//! of `make budgets`.
//!
//! ```text
//! cargo build --release -p bammetrics-formats --features parquet --example size_parquet
//! ```

#[path = "size_formats.rs"]
mod probe;

fn main() {
    probe::main();
}

//! `IntervalListTools` of Picard 3.5.0 over `bammetrics_core::interval_ops`,
//! for the concordance module `tests/concordance/test_interval_ops.py`:
//!
//!     cargo run -p bammetrics-core --example interval_ops -- \
//!         ACTION=SUBTRACT I=a.interval_list SI=b.interval_list O=out.interval_list
//!
//! Arguments are Picard's (`KEY=VALUE`, long or short names); the run is the
//! registry tool `crates/core/src/interval_list_tools.rs`, so the `@PG`
//! line holds Picard's command line for them. Exit 0, or 1 with the error.

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if let Err(error) = bammetrics_core::interval_list_tools::TOOL.run_cli(&args) {
        eprintln!("interval_ops: {error}");
        std::process::exit(1);
    }
}

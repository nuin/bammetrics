use std::path::PathBuf;
use std::process::ExitCode;

use bammetrics_core::insert_size::{self, Options};
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "bammetrics",
    version,
    about = "Picard-concordant BAM QC metrics"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// CollectInsertSizeMetrics: one line per pair orientation.
    InsertSize {
        /// Input BAM.
        bam: PathBuf,
    },
}

fn run(cli: Cli) -> bammetrics_core::Result<()> {
    match cli.command {
        Command::InsertSize { bam } => {
            for r in insert_size::run_over_bam(&bam, &Options::default())? {
                println!(
                    "ORIENT={} READ_PAIRS={} MIN={} MAX={} MEDIAN={} MODE={} MAD={} MEAN={:.4} SD={:.4} W50={}",
                    r.orientation.as_str(),
                    r.read_pairs,
                    r.min_insert_size,
                    r.max_insert_size,
                    r.median_insert_size,
                    r.mode_insert_size,
                    r.median_absolute_deviation,
                    r.mean_insert_size,
                    r.standard_deviation,
                    r.width_of_50_percent
                );
            }
            Ok(())
        }
    }
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("bammetrics: {e}");
            ExitCode::FAILURE
        }
    }
}

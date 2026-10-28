use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bammetrics_core::Error;
use bammetrics_core::insert_size::{self, Options};
use clap::{Args, Parser, Subcommand};

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
    /// CollectInsertSizeMetrics: a Picard InsertSizeMetrics file, one row per
    /// pair orientation plus the trimmed histograms.
    InsertSize(InsertSizeArgs),
}

#[derive(Args)]
struct InsertSizeArgs {
    /// Input BAM.
    bam: PathBuf,
    /// Metrics file to write (Picard OUTPUT). Standard output when absent.
    #[arg(short = 'O', long)]
    output: Option<PathBuf>,
    /// Drop orientations below this fraction of all pairs (Picard MINIMUM_PCT).
    #[arg(short = 'M', long, default_value_t = 0.05)]
    minimum_pct: f32,
    /// Trim mean and SD to median + DEVIATIONS * MAD (Picard DEVIATIONS).
    #[arg(long, default_value_t = 10.0)]
    deviations: f64,
    /// Trim to this insert size instead (Picard HISTOGRAM_WIDTH).
    #[arg(short = 'W', long)]
    histogram_width: Option<i64>,
    /// Lower bound on the computed trim width (Picard MIN_HISTOGRAM_WIDTH).
    #[arg(long)]
    min_histogram_width: Option<i64>,
    /// Count reads flagged as duplicates (Picard INCLUDE_DUPLICATES).
    #[arg(long)]
    include_duplicates: bool,
}

/// `StringHeader` lines in Picard's shape: the tool name and its arguments
/// (parsers such as MultiQC read `INPUT=` from it), then the producer.
fn insert_size_headers(args: &InsertSizeArgs) -> Vec<String> {
    let mut command = format!(
        "CollectInsertSizeMetrics INPUT={} MINIMUM_PCT={:?} DEVIATIONS={:?} INCLUDE_DUPLICATES={} METRIC_ACCUMULATION_LEVEL=[ALL_READS]",
        args.bam.display(),
        args.minimum_pct,
        args.deviations,
        args.include_duplicates
    );
    if let Some(output) = &args.output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(w) = args.histogram_width {
        command.push_str(&format!(" HISTOGRAM_WIDTH={w}"));
    }
    if let Some(w) = args.min_histogram_width {
        command.push_str(&format!(" MIN_HISTOGRAM_WIDTH={w}"));
    }
    vec![
        command,
        format!("Started by bammetrics {}", env!("CARGO_PKG_VERSION")),
    ]
}

type Output = BufWriter<Box<dyn Write>>;

/// Run `write` against `path`, or standard output when `path` is `None`.
fn write_to(
    path: Option<&Path>,
    write: impl FnOnce(&mut Output) -> io::Result<()>,
) -> bammetrics_core::Result<()> {
    let output_error = |source| Error::Output {
        path: path.map_or_else(|| PathBuf::from("-"), Path::to_path_buf),
        source,
    };
    let sink: Box<dyn Write> = match path {
        Some(p) => Box::new(File::create(p).map_err(output_error)?),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);
    write(&mut out).map_err(output_error)?;
    out.flush().map_err(output_error)
}

fn insert_size(args: InsertSizeArgs) -> bammetrics_core::Result<()> {
    let opts = Options {
        minimum_pct: args.minimum_pct,
        deviations: args.deviations,
        histogram_width: args.histogram_width,
        min_histogram_width: args.min_histogram_width,
        include_duplicates: args.include_duplicates,
    };
    let rows = insert_size::run_over_bam(&args.bam, &opts)?;
    // CollectInsertSizeMetrics.finish writes no file when every orientation
    // was dropped (or there were no pairs), and only logs a warning.
    if rows.is_empty() {
        eprintln!(
            "bammetrics: all data categories were discarded because they contained < {:?} of the total aligned paired data; no metrics written",
            args.minimum_pct
        );
        return Ok(());
    }
    let headers = insert_size_headers(&args);
    write_to(args.output.as_deref(), |out| {
        insert_size::write_metrics(&rows, headers, out)
    })
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::InsertSize(args) => insert_size(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("bammetrics: {e}");
            ExitCode::FAILURE
        }
    }
}

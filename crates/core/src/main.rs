use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use bammetrics_core::Error;
use bammetrics_core::alignment_summary::{self, DEFAULT_MAX_INSERT_SIZE};
use bammetrics_core::insert_size::{self, Options};
use bammetrics_core::metrics_file::started_on;
use bammetrics_core::picard_cli;
use bammetrics_core::quality_yield;
use clap::{ArgAction, Args, Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "bammetrics",
    version,
    about = "Picard-concordant BAM QC metrics",
    after_help = "Picard syntax also works: bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt [M=0.05 ...]\n\
                  bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt [R=ref.fa ...]\n\
                  bammetrics CollectQualityYieldMetrics I=in.bam O=out.txt [OQ=false ...]"
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
    /// CollectAlignmentSummaryMetrics: a Picard AlignmentSummaryMetrics file,
    /// one row per read category plus the read-length histograms.
    AlignmentSummary(AlignmentSummaryArgs),
    /// CollectQualityYieldMetrics: a Picard QualityYieldMetrics file, one
    /// row of read, base and Q20/Q30 totals.
    QualityYield(QualityYieldArgs),
}

#[derive(Args)]
struct InsertSizeArgs {
    /// Input BAM (or --input).
    #[arg(required_unless_present = "input", conflicts_with = "input")]
    bam: Option<PathBuf>,
    /// Input BAM (Picard INPUT).
    #[arg(short = 'I', long)]
    input: Option<PathBuf>,
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

#[derive(Args)]
struct AlignmentSummaryArgs {
    /// Input BAM (or --input).
    #[arg(required_unless_present = "input", conflicts_with = "input")]
    bam: Option<PathBuf>,
    /// Input BAM (Picard INPUT).
    #[arg(short = 'I', long)]
    input: Option<PathBuf>,
    /// Metrics file to write (Picard OUTPUT). Standard output when absent.
    #[arg(short = 'O', long)]
    output: Option<PathBuf>,
    /// Reference FASTA with a .fai index (Picard REFERENCE_SEQUENCE). Without
    /// it no base counts as a mismatch.
    #[arg(short = 'R', long)]
    reference: Option<PathBuf>,
    /// Pairs with a larger |TLEN| are chimeric (Picard MAX_INSERT_SIZE).
    #[arg(long, default_value_t = DEFAULT_MAX_INSERT_SIZE)]
    max_insert_size: i32,
    /// Adds FR, RF or TANDEM to the expected orientations, FR by default;
    /// "null" clears the list (Picard EXPECTED_PAIR_ORIENTATIONS).
    #[arg(long = "expected-pair-orientation", value_name = "ORIENTATION")]
    expected_pair_orientations: Vec<String>,
    /// Adds an adapter to Picard's default list; "null" clears the list
    /// (Picard ADAPTER_SEQUENCE).
    #[arg(long = "adapter-sequence", value_name = "SEQUENCE")]
    adapter_sequences: Vec<String>,
    /// false counts reads only (Picard COLLECT_ALIGNMENT_INFORMATION).
    #[arg(long, default_value = "true", value_parser = parse_bool, action = ArgAction::Set)]
    collect_alignment_information: bool,
}

#[derive(Args)]
struct QualityYieldArgs {
    /// Input BAM (or --input).
    #[arg(required_unless_present = "input", conflicts_with = "input")]
    bam: Option<PathBuf>,
    /// Input BAM (Picard INPUT).
    #[arg(short = 'I', long)]
    input: Option<PathBuf>,
    /// Metrics file to write (Picard OUTPUT). Standard output when absent.
    #[arg(short = 'O', long)]
    output: Option<PathBuf>,
    /// false reads QUAL even where an OQ tag exists (Picard USE_ORIGINAL_QUALITIES).
    #[arg(long, default_value = "true", value_parser = parse_bool, action = ArgAction::Set)]
    use_original_qualities: bool,
    /// Count secondary alignments (Picard INCLUDE_SECONDARY_ALIGNMENTS).
    #[arg(long)]
    include_secondary_alignments: bool,
    /// Count supplementary alignments (Picard INCLUDE_SUPPLEMENTAL_ALIGNMENTS).
    #[arg(long)]
    include_supplemental_alignments: bool,
}

/// A Picard collection argument: values add to the default, and `null`
/// empties the list so far.
fn picard_list<T: PartialEq>(
    default: Vec<T>,
    values: &[String],
    parse: impl Fn(&str) -> bammetrics_core::Result<T>,
    distinct: bool,
) -> bammetrics_core::Result<Vec<T>> {
    let mut list = default;
    for value in values {
        if value == "null" {
            list.clear();
            continue;
        }
        let item = parse(value)?;
        if !(distinct && list.contains(&item)) {
            list.push(item);
        }
    }
    Ok(list)
}

/// Booleans as Picard reads them, in any case.
fn parse_bool(value: &str) -> Result<bool, String> {
    match value.to_ascii_lowercase().as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        _ => Err(format!("expected true or false, got {value:?}")),
    }
}

impl AlignmentSummaryArgs {
    fn bam(&self) -> bammetrics_core::Result<&Path> {
        self.bam
            .as_deref()
            .or(self.input.as_deref())
            .ok_or_else(|| Error::Usage("no input BAM".to_owned()))
    }

    fn options(&self) -> bammetrics_core::Result<alignment_summary::Options> {
        let defaults = alignment_summary::Options::default();
        Ok(alignment_summary::Options {
            reference: self.reference.clone(),
            max_insert_size: self.max_insert_size,
            expected_orientations: picard_list(
                defaults.expected_orientations,
                &self.expected_pair_orientations,
                str::parse,
                true,
            )?,
            adapter_sequences: picard_list(
                defaults.adapter_sequences,
                &self.adapter_sequences,
                |s| Ok(s.to_owned()),
                false,
            )?,
            collect_alignment_information: self.collect_alignment_information,
        })
    }
}

impl InsertSizeArgs {
    /// The BAM given positionally or with `--input`; clap guarantees one.
    fn bam(&self) -> bammetrics_core::Result<&Path> {
        self.bam
            .as_deref()
            .or(self.input.as_deref())
            .ok_or_else(|| Error::Usage("no input BAM".to_owned()))
    }
}

/// `StringHeader` lines in Picard's shape: the tool name and its arguments
/// (parsers such as MultiQC read `INPUT=` from it), then `Started on: <time>`.
fn insert_size_headers(args: &InsertSizeArgs, bam: &Path) -> Vec<String> {
    let mut command = format!(
        "CollectInsertSizeMetrics INPUT={} MINIMUM_PCT={:?} DEVIATIONS={:?} INCLUDE_DUPLICATES={} METRIC_ACCUMULATION_LEVEL=[ALL_READS]",
        bam.display(),
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
    vec![command, started_on(now())]
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
    let bam = args.bam()?;
    let rows = insert_size::run_over_bam(bam, &opts)?;
    // CollectInsertSizeMetrics.finish writes no file when every orientation
    // was dropped (or there were no pairs), and only logs a warning.
    if rows.is_empty() {
        eprintln!(
            "bammetrics: all data categories were discarded because they contained < {:?} of the total aligned paired data; no metrics written",
            args.minimum_pct
        );
        return Ok(());
    }
    let headers = insert_size_headers(&args, bam);
    write_to(args.output.as_deref(), |out| {
        insert_size::write_metrics(&rows, headers, out)
    })
}

/// `StringHeader` lines in Picard's shape for `CollectAlignmentSummaryMetrics`.
fn alignment_summary_headers(
    args: &AlignmentSummaryArgs,
    bam: &Path,
    opts: &alignment_summary::Options,
) -> Vec<String> {
    let mut command = format!("CollectAlignmentSummaryMetrics INPUT={}", bam.display());
    if let Some(output) = &args.output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    if let Some(reference) = &opts.reference {
        command.push_str(&format!(" REFERENCE_SEQUENCE={}", reference.display()));
    }
    let orientations: Vec<&str> = opts
        .expected_orientations
        .iter()
        .map(|o| o.as_str())
        .collect();
    command.push_str(&format!(
        " MAX_INSERT_SIZE={} EXPECTED_PAIR_ORIENTATIONS=[{}] ADAPTER_SEQUENCE=[{}] METRIC_ACCUMULATION_LEVEL=[ALL_READS] IS_BISULFITE_SEQUENCED=false COLLECT_ALIGNMENT_INFORMATION={}",
        opts.max_insert_size,
        orientations.join(", "),
        opts.adapter_sequences.join(", "),
        opts.collect_alignment_information
    ));
    vec![command, started_on(now())]
}

/// `StringHeader` lines in Picard's shape for `CollectQualityYieldMetrics`.
fn quality_yield_headers(
    args: &QualityYieldArgs,
    bam: &Path,
    opts: &quality_yield::Options,
) -> Vec<String> {
    let mut command = format!("CollectQualityYieldMetrics INPUT={}", bam.display());
    if let Some(output) = &args.output {
        command.push_str(&format!(" OUTPUT={}", output.display()));
    }
    command.push_str(&format!(
        " USE_ORIGINAL_QUALITIES={} INCLUDE_SECONDARY_ALIGNMENTS={} INCLUDE_SUPPLEMENTAL_ALIGNMENTS={}",
        opts.use_original_qualities,
        opts.include_secondary_alignments,
        opts.include_supplemental_alignments
    ));
    vec![command, started_on(now())]
}

fn now() -> u64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_secs())
}

fn alignment_summary(args: AlignmentSummaryArgs) -> bammetrics_core::Result<()> {
    let bam = args.bam()?;
    let opts = args.options()?;
    let report = alignment_summary::run_over_bam(bam, &opts)?;
    let headers = alignment_summary_headers(&args, bam, &opts);
    write_to(args.output.as_deref(), |out| {
        alignment_summary::write_metrics(&report, headers, out)
    })
}

fn quality_yield(args: QualityYieldArgs) -> bammetrics_core::Result<()> {
    let bam = args
        .bam
        .as_deref()
        .or(args.input.as_deref())
        .ok_or_else(|| Error::Usage("no input BAM".to_owned()))?;
    let opts = quality_yield::Options {
        use_original_qualities: args.use_original_qualities,
        include_secondary_alignments: args.include_secondary_alignments,
        include_supplemental_alignments: args.include_supplemental_alignments,
    };
    let row = quality_yield::run_over_bam(bam, &opts)?;
    let headers = quality_yield_headers(&args, bam, &opts);
    write_to(args.output.as_deref(), |out| {
        quality_yield::write_metrics(&row, headers, out)
    })
}

/// The command line as clap sees it. A first argument naming a Picard tool
/// switches to Picard syntax and is rewritten to the native flags.
fn command_line() -> bammetrics_core::Result<Vec<String>> {
    let mut argv = std::env::args();
    let program = argv.next().unwrap_or_else(|| "bammetrics".to_owned());
    let rest: Vec<String> = argv.collect();
    let Some(translation) = rest
        .split_first()
        .map(|(tool, options)| picard_cli::translate(tool, options))
        .transpose()?
        .flatten()
    else {
        return Ok(std::iter::once(program).chain(rest).collect());
    };
    for note in &translation.notes {
        eprintln!("bammetrics: {note}");
    }
    Ok([program, translation.subcommand.to_owned()]
        .into_iter()
        .chain(translation.args)
        .collect())
}

fn run() -> bammetrics_core::Result<()> {
    match Cli::parse_from(command_line()?).command {
        Command::InsertSize(args) => insert_size(args),
        Command::AlignmentSummary(args) => alignment_summary(args),
        Command::QualityYield(args) => quality_yield(args),
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("bammetrics: {e}");
            ExitCode::FAILURE
        }
    }
}

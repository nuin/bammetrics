use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use bammetrics_core::Error;
use bammetrics_core::alignment_summary::{self, AlignmentSummary, DEFAULT_MAX_INSERT_SIZE};
use bammetrics_core::collect::{self, Outputs, Request};
use bammetrics_core::hs_metrics::{self, HsMetrics};
use bammetrics_core::insert_size::{self, InsertSize};
use bammetrics_core::metrics_file::headers;
use bammetrics_core::pass;
use bammetrics_core::picard_cli;
use bammetrics_core::quality_yield::{self, QualityYield};
use clap::{ArgAction, Args, Parser, Subcommand, ValueEnum};

#[derive(Parser)]
#[command(
    name = "bammetrics",
    version,
    about = "Picard-concordant BAM QC metrics",
    after_help = "Picard syntax also works: bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt [M=0.05 ...]\n\
                  bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt [R=ref.fa ...]\n\
                  bammetrics CollectQualityYieldMetrics I=in.bam O=out.txt [OQ=false ...]\n\
                  bammetrics CollectHsMetrics I=in.bam O=out.txt BAIT_INTERVALS=b.interval_list TARGET_INTERVALS=t.interval_list [R=ref.fa PER_BASE_COVERAGE=bases.txt ...]"
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
    /// CollectHsMetrics: a Picard HsMetrics file, one row of bait, target and
    /// coverage metrics plus the depth and base-quality histograms.
    HsMetrics(HsMetricsArgs),
    /// Several metrics from one pass over the BAM, one Picard file each:
    /// PREFIX.quality_yield_metrics, .alignment_summary_metrics,
    /// .insert_size_metrics and .hs_metrics.
    Collect(CollectArgs),
}

/// BGZF decompression workers. Records are processed on the main thread in
/// file order, so the metrics do not depend on this.
#[derive(Args)]
struct Threads {
    /// BGZF decompression threads (default: available cores, at most 8).
    #[arg(long, default_value_t = pass::default_threads())]
    threads: usize,
}

/// A metric `collect` can compute.
#[derive(Clone, Copy, PartialEq, Eq, ValueEnum)]
enum MetricName {
    /// CollectQualityYieldMetrics.
    QualityYield,
    /// CollectAlignmentSummaryMetrics.
    AlignmentSummary,
    /// CollectInsertSizeMetrics.
    InsertSize,
    /// CollectHsMetrics, once per --hs pair.
    Hs,
}

#[derive(Args)]
struct CollectArgs {
    /// Input BAM (or --input).
    #[arg(required_unless_present = "input", conflicts_with = "input")]
    bam: Option<PathBuf>,
    /// Input BAM.
    #[arg(short = 'I', long)]
    input: Option<PathBuf>,
    /// Output prefix. Files are PREFIX plus .quality_yield_metrics,
    /// .alignment_summary_metrics, .insert_size_metrics (the names Picard
    /// CollectMultipleMetrics uses) and .hs_metrics, or .BAIT_SET.hs_metrics
    /// with several --hs pairs.
    #[arg(short = 'O', long)]
    output: PathBuf,
    /// Metrics to compute, comma separated (default: quality-yield,
    /// alignment-summary, insert-size, and hs when --hs is given).
    #[arg(long = "metric", value_enum, value_delimiter = ',')]
    metrics: Vec<MetricName>,
    /// Reference FASTA with a .fai index, for alignment-summary and hs.
    #[arg(short = 'R', long)]
    reference: Option<PathBuf>,
    /// One CollectHsMetrics run: a bait and a target interval file
    /// (interval list or BED). Repeat for more runs.
    #[arg(long = "hs", num_args = 2, value_names = ["BAITS", "TARGETS"], action = ArgAction::Append)]
    hs: Vec<PathBuf>,
    /// Also write each --hs run's per-base target coverage (Picard
    /// PER_BASE_COVERAGE), named as its HS file with .per_base_coverage in
    /// place of .hs_metrics.
    #[arg(long, requires = "hs")]
    hs_per_base: bool,
    #[command(flatten)]
    threads: Threads,
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
    #[command(flatten)]
    threads: Threads,
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
    #[command(flatten)]
    threads: Threads,
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
    #[command(flatten)]
    threads: Threads,
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

#[derive(Args)]
struct HsMetricsArgs {
    /// Input BAM (or --input).
    #[arg(required_unless_present = "input", conflicts_with = "input")]
    bam: Option<PathBuf>,
    /// Input BAM (Picard INPUT).
    #[arg(short = 'I', long)]
    input: Option<PathBuf>,
    /// Metrics file to write (Picard OUTPUT). Standard output when absent.
    #[arg(short = 'O', long)]
    output: Option<PathBuf>,
    #[command(flatten)]
    threads: Threads,
    /// Bait intervals: an interval list, or a BED file (Picard BAIT_INTERVALS).
    /// Repeat for several files.
    #[arg(long = "bait-intervals", value_name = "FILE", required = true)]
    bait_intervals: Vec<PathBuf>,
    /// Target intervals, as for --bait-intervals (Picard TARGET_INTERVALS).
    #[arg(long = "target-intervals", value_name = "FILE", required = true)]
    target_intervals: Vec<PathBuf>,
    /// BAIT_SET column; by default the bait file names up to the first dot
    /// (Picard BAIT_SET_NAME).
    #[arg(long)]
    bait_set_name: Option<String>,
    /// Reference FASTA with a .fai index, for AT_DROPOUT and GC_DROPOUT
    /// (Picard REFERENCE_SEQUENCE).
    #[arg(short = 'R', long)]
    reference: Option<PathBuf>,
    /// Reads below this MAPQ add no coverage (Picard MINIMUM_MAPPING_QUALITY).
    #[arg(long, default_value_t = hs_metrics::DEFAULT_MINIMUM_MAPPING_QUALITY)]
    minimum_mapping_quality: i32,
    /// Bases below this quality add no coverage (Picard MINIMUM_BASE_QUALITY).
    #[arg(long, default_value_t = hs_metrics::DEFAULT_MINIMUM_BASE_QUALITY)]
    minimum_base_quality: i32,
    /// false counts both mates where they overlap (Picard CLIP_OVERLAPPING_READS).
    #[arg(long, default_value = "true", value_parser = parse_bool, action = ArgAction::Set)]
    clip_overlapping_reads: bool,
    /// Count inserted bases as on target and deleted bases as covered
    /// (Picard INCLUDE_INDELS).
    #[arg(long)]
    include_indels: bool,
    /// Depth cap of the het SNP sensitivity model (Picard COVERAGE_CAP).
    #[arg(long, default_value_t = hs_metrics::DEFAULT_COVERAGE_CAP)]
    coverage_cap: i32,
    /// Reads this close to a bait are near bait (Picard NEAR_DISTANCE).
    #[arg(long, default_value_t = hs_metrics::DEFAULT_NEAR_DISTANCE)]
    near_distance: i32,
    /// Samples of the het SNP sensitivity simulation (Picard SAMPLE_SIZE).
    #[arg(long, default_value_t = hs_metrics::DEFAULT_SAMPLE_SIZE)]
    sample_size: usize,
    /// Also write the high-quality depth of every target base to this file;
    /// gzip when it ends in .gz (Picard PER_BASE_COVERAGE).
    #[arg(long, value_name = "FILE")]
    per_base_coverage: Option<PathBuf>,
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

/// The BAM given positionally or with `--input`; clap guarantees one.
fn input_bam<'a>(
    bam: &'a Option<PathBuf>,
    input: &'a Option<PathBuf>,
) -> bammetrics_core::Result<&'a Path> {
    bam.as_deref()
        .or(input.as_deref())
        .ok_or_else(|| Error::Usage("no input BAM".to_owned()))
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

/// The warning `CollectInsertSizeMetrics.finish` logs when it writes no file.
fn no_insert_size_rows(minimum_pct: f32) {
    eprintln!(
        "bammetrics: all data categories were discarded because they contained < {minimum_pct:?} of the total aligned paired data; no metrics written"
    );
}

impl InsertSizeArgs {
    fn options(&self) -> insert_size::Options {
        insert_size::Options {
            minimum_pct: self.minimum_pct,
            deviations: self.deviations,
            histogram_width: self.histogram_width,
            min_histogram_width: self.min_histogram_width,
            include_duplicates: self.include_duplicates,
        }
    }
}

impl QualityYieldArgs {
    fn options(&self) -> quality_yield::Options {
        quality_yield::Options {
            use_original_qualities: self.use_original_qualities,
            include_secondary_alignments: self.include_secondary_alignments,
            include_supplemental_alignments: self.include_supplemental_alignments,
        }
    }
}

impl HsMetricsArgs {
    fn options(&self) -> hs_metrics::Options {
        hs_metrics::Options {
            bait_intervals: self.bait_intervals.clone(),
            target_intervals: self.target_intervals.clone(),
            bait_set_name: self.bait_set_name.clone(),
            reference: self.reference.clone(),
            minimum_mapping_quality: self.minimum_mapping_quality,
            minimum_base_quality: self.minimum_base_quality,
            clip_overlapping_reads: self.clip_overlapping_reads,
            include_indels: self.include_indels,
            coverage_cap: self.coverage_cap,
            near_distance: self.near_distance,
            sample_size: self.sample_size,
            per_base_coverage: self.per_base_coverage.clone(),
        }
    }
}

fn insert_size(args: InsertSizeArgs) -> bammetrics_core::Result<()> {
    let bam = input_bam(&args.bam, &args.input)?;
    let opts = args.options();
    let rows = pass::run(bam, args.threads.threads, |_| Ok(InsertSize::new(&opts)))?;
    // CollectInsertSizeMetrics.finish writes no file when every orientation
    // was dropped (or there were no pairs), and only logs a warning.
    if rows.is_empty() {
        no_insert_size_rows(opts.minimum_pct);
        return Ok(());
    }
    let output = args.output.as_deref();
    let headers = headers(insert_size::command_header(&opts, bam, output));
    write_to(output, |out| {
        insert_size::write_metrics(&rows, headers, out)
    })
}

fn alignment_summary(args: AlignmentSummaryArgs) -> bammetrics_core::Result<()> {
    let bam = input_bam(&args.bam, &args.input)?;
    let opts = args.options()?;
    let report = pass::run(bam, args.threads.threads, |header| {
        AlignmentSummary::new(&opts, header)
    })?;
    let output = args.output.as_deref();
    let headers = headers(alignment_summary::command_header(&opts, bam, output));
    write_to(output, |out| {
        alignment_summary::write_metrics(&report, headers, out)
    })
}

fn quality_yield(args: QualityYieldArgs) -> bammetrics_core::Result<()> {
    let bam = input_bam(&args.bam, &args.input)?;
    let opts = args.options();
    let row = pass::run(bam, args.threads.threads, |_| Ok(QualityYield::new(&opts)))?;
    let output = args.output.as_deref();
    let headers = headers(quality_yield::command_header(&opts, bam, output));
    write_to(output, |out| {
        quality_yield::write_metrics(&row, headers, out)
    })
}

fn hs_metrics(args: HsMetricsArgs) -> bammetrics_core::Result<()> {
    let bam = input_bam(&args.bam, &args.input)?;
    let opts = args.options();
    let report = pass::run(bam, args.threads.threads, |header| {
        HsMetrics::new(&opts, header)
    })?;
    let output = args.output.as_deref();
    let headers = headers(hs_metrics::command_header(&opts, bam, output));
    write_to(output, |out| {
        hs_metrics::write_metrics(&report, headers, out)
    })
}

impl CollectArgs {
    /// Picard's defaults for every metric, with `-R` given to the two that
    /// read a reference (`CollectInsertSizeMetrics` does not).
    fn request(&self) -> bammetrics_core::Result<Request> {
        let wanted = if self.metrics.is_empty() {
            let mut all = vec![
                MetricName::QualityYield,
                MetricName::AlignmentSummary,
                MetricName::InsertSize,
            ];
            if !self.hs.is_empty() {
                all.push(MetricName::Hs);
            }
            all
        } else {
            self.metrics.clone()
        };
        let has_hs = wanted.contains(&MetricName::Hs);
        if has_hs == self.hs.is_empty() {
            return Err(Error::Usage(
                "--metric hs and --hs BAITS TARGETS go together".to_owned(),
            ));
        }
        let hs_metrics = self
            .hs
            .chunks(2)
            .map(|pair| hs_metrics::Options {
                bait_intervals: pair[..1].to_vec(),
                target_intervals: pair[1..].to_vec(),
                reference: self.reference.clone(),
                ..hs_metrics::Options::default()
            })
            .collect();
        Ok(Request {
            quality_yield: wanted
                .contains(&MetricName::QualityYield)
                .then(quality_yield::Options::default),
            alignment_summary: wanted.contains(&MetricName::AlignmentSummary).then(|| {
                alignment_summary::Options {
                    reference: self.reference.clone(),
                    ..alignment_summary::Options::default()
                }
            }),
            insert_size: wanted
                .contains(&MetricName::InsertSize)
                .then(insert_size::Options::default),
            hs_metrics,
        })
    }
}

fn collect(args: CollectArgs) -> bammetrics_core::Result<()> {
    let bam = input_bam(&args.bam, &args.input)?;
    let mut request = args.request()?;
    let outputs = Outputs::new(&args.output, &request)?;
    if args.hs_per_base {
        for (opts, path) in request.hs_metrics.iter_mut().zip(&outputs.hs_metrics) {
            opts.per_base_coverage = Some(collect::per_base_coverage_path(path));
        }
    }
    let results = collect::run_over_bam(bam, &request, args.threads.threads)?;
    collect::write_files(&request, &results, bam, &outputs)?;
    if let (Some(opts), Some(Ok(rows))) = (&request.insert_size, &results.insert_size)
        && rows.is_empty()
    {
        no_insert_size_rows(opts.minimum_pct);
    }
    let errors = results.errors();
    for (tool, e) in &errors {
        eprintln!("bammetrics: {tool}: {e}");
    }
    if !errors.is_empty() {
        return Err(Error::Metric(format!("{} metric(s) failed", errors.len())));
    }
    Ok(())
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
        Command::HsMetrics(args) => hs_metrics(args),
        Command::Collect(args) => collect(args),
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

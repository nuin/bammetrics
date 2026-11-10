//! The tool registry: Picard tools by name, with Picard-style arguments.
//!
//! A tool is a [`Tool`] constant in its own module: the Picard name, the
//! options it takes and a `run` function from [`Invocation`] to
//! [`ToolOutput`]. Adding it to [`TOOLS`] is the only edit outside that
//! module; the command line (`bammetrics <ToolName> I=... O=...`) and the
//! Python binding then both serve it, because both go through
//! [`Tool::execute`]. A tool that reads the BAM once also gives a
//! `collector`, which `CollectMultipleMetrics` runs beside the other
//! programs' collectors on one read of the file.
//!
//! Arguments are the forms Picard 3.5.0 accepts: `KEY=VALUE`, `--KEY VALUE`,
//! `--KEY=VALUE` and `-K VALUE`, with long or short names. An option no
//! entry names is an error, because dropping an option Picard applies would
//! change the numbers silently.

use std::collections::HashMap;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use noodles::sam;

use crate::cram;
use crate::error::{Error, Result};
use crate::metric::ToolOutput;
use crate::pass::{self, Accumulator, Collector, Limits, default_threads};

/// What an option does once recognised.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Read by the tool through [`Invocation`].
    Value,
    /// The metrics file path ([`Invocation::output`]).
    Output,
    /// Never changes the metrics: accepted, dropped, noted on stderr.
    Ignored(&'static str),
    /// A [`COMMON`] option this tool does not take (Picard's parser
    /// rejects it as unrecognized).
    Rejected,
}

/// The type of an option's value. The verb grammar reads `Bool` options
/// as flags (`--x`, `--no-x`, `--x=false`); the Python stub types the
/// keyword argument from it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ty {
    Bool,
    Int,
    Float,
    Text,
    Path,
}

impl Ty {
    /// The Python annotation of one value.
    pub fn python(self) -> &'static str {
        match self {
            Ty::Bool => "bool",
            Ty::Int => "int",
            Ty::Float => "float",
            Ty::Text => "str",
            Ty::Path => "_Path",
        }
    }

    /// The placeholder `--help` prints after the flag.
    pub fn placeholder(self) -> &'static str {
        match self {
            Ty::Bool => "",
            Ty::Int => "INT",
            Ty::Float => "NUM",
            Ty::Text => "TEXT",
            Ty::Path => "FILE",
        }
    }
}

/// One option row: its Picard names (the long name first), kind, value
/// type and help text. The row is the only place an option is named: the
/// Picard parser, the verb grammar (`--minimum-pct`), `--help`, the Python
/// keyword argument (`minimum_pct`) and the `bammetrics.pyi` stub all come
/// from it.
#[derive(Clone, Copy, Debug)]
pub struct Opt {
    pub names: &'static [&'static str],
    pub kind: Kind,
    pub ty: Ty,
    /// A Picard collection option: every occurrence adds a value.
    pub repeated: bool,
    /// Picard fails without it.
    pub required: bool,
    /// A path option naming a file the tool writes. Every other path
    /// option names a file it reads, which [`write_files`] refuses to
    /// replace.
    pub written: bool,
    pub help: &'static str,
}

impl Opt {
    const fn new(names: &'static [&'static str], kind: Kind, ty: Ty, help: &'static str) -> Self {
        Self {
            names,
            kind,
            ty,
            repeated: false,
            required: false,
            written: false,
            help,
        }
    }

    /// The path option names a file the tool writes.
    pub const fn written(mut self) -> Self {
        self.written = true;
        self
    }

    pub const fn bool(names: &'static [&'static str], help: &'static str) -> Self {
        Self::new(names, Kind::Value, Ty::Bool, help)
    }

    pub const fn int(names: &'static [&'static str], help: &'static str) -> Self {
        Self::new(names, Kind::Value, Ty::Int, help)
    }

    pub const fn float(names: &'static [&'static str], help: &'static str) -> Self {
        Self::new(names, Kind::Value, Ty::Float, help)
    }

    pub const fn text(names: &'static [&'static str], help: &'static str) -> Self {
        Self::new(names, Kind::Value, Ty::Text, help)
    }

    pub const fn path(names: &'static [&'static str], help: &'static str) -> Self {
        Self::new(names, Kind::Value, Ty::Path, help)
    }

    /// Every occurrence adds a value (Picard's collection options).
    pub const fn repeated(mut self) -> Self {
        self.repeated = true;
        self
    }

    /// Picard fails without the option.
    pub const fn required(mut self) -> Self {
        self.required = true;
        self
    }

    /// The row with its own `--help` text.
    pub const fn with_help(mut self, help: &'static str) -> Self {
        self.help = help;
        self
    }

    /// The option naming the metrics file. A tool that declares one takes
    /// no `OUTPUT` (Picard's `CollectDuplicateMetrics` writes to `M`).
    pub const fn output(names: &'static [&'static str]) -> Self {
        Self::new(
            names,
            Kind::Output,
            Ty::Path,
            "The metrics file; standard output without one.",
        )
    }

    pub const fn ignored(names: &'static [&'static str], why: &'static str) -> Self {
        Self::new(names, Kind::Ignored(why), Ty::Text, why)
    }

    /// An ignored Picard integer: a value that is not one is an error, as
    /// in Picard.
    pub const fn ignored_int(names: &'static [&'static str], why: &'static str) -> Self {
        Self::new(names, Kind::Ignored(why), Ty::Int, why)
    }

    /// An ignored Picard number: a value that is not one is an error.
    pub const fn ignored_float(names: &'static [&'static str], why: &'static str) -> Self {
        Self::new(names, Kind::Ignored(why), Ty::Float, why)
    }

    /// An ignored Picard boolean (`CREATE_INDEX`): the verb grammar reads
    /// it as a flag, so `--create-index in.bam` keeps the BAM as the input.
    pub const fn ignored_bool(names: &'static [&'static str], why: &'static str) -> Self {
        Self::new(names, Kind::Ignored(why), Ty::Bool, why)
    }

    /// A common option Picard rejects for this tool.
    pub const fn rejected(names: &'static [&'static str]) -> Self {
        Self::new(names, Kind::Rejected, Ty::Text, "")
    }

    /// The Picard long name, `MINIMUM_PCT`.
    pub fn picard(&self) -> &'static str {
        self.names[0]
    }

    /// The verb flag without its dashes, `minimum-pct`.
    pub fn kebab(&self) -> String {
        self.picard().to_ascii_lowercase().replace('_', "-")
    }

    /// The Python keyword argument, `minimum_pct`.
    pub fn kwarg(&self) -> String {
        self.picard().to_ascii_lowercase()
    }

    /// Whether `name` is one of the option's Picard names once normalised
    /// ([`normalize`]): `minimum-pct`, `MINIMUM_PCT` and `m` all name
    /// `MINIMUM_PCT`, whose short name is `M`.
    pub fn matches(&self, name: &str) -> bool {
        let wanted = normalize(name);
        self.names.iter().any(|n| n.eq_ignore_ascii_case(&wanted))
    }
}

/// An option name as the verb grammar compares it: upper case, `-` as `_`.
pub fn normalize(name: &str) -> String {
    name.to_ascii_uppercase().replace('-', "_")
}

/// Options every Picard tool takes from `CommandLineProgram` and
/// `SinglePassSamProgram`. `ASSUME_SORTED=false` fails a file whose header
/// is not coordinate sorted ([`Invocation::single_pass`]); the others change
/// no metric.
pub const COMMON: &[Opt] = &[
    Opt::path(&["INPUT", "I"], "The BAM or CRAM to read.").required(),
    Opt::output(&["OUTPUT", "O"]),
    Opt::ignored(&["VALIDATION_STRINGENCY"], "records are not validated"),
    Opt::ignored(&["TMP_DIR"], "no temporary files are written"),
    Opt::ignored_bool(&["QUIET"], "progress is not logged"),
    Opt::ignored(&["VERBOSITY"], "progress is not logged"),
    Opt::ignored(&["COMPRESSION_LEVEL"], "no BAM is written"),
    Opt::ignored(&["MAX_RECORDS_IN_RAM"], "records are streamed"),
    Opt::ignored_bool(&["CREATE_INDEX"], "no BAM is written"),
    Opt::ignored_bool(&["CREATE_MD5_FILE"], "no BAM is written"),
    Opt::ignored(&["GA4GH_CLIENT_SECRETS"], "no GA4GH access"),
    Opt::ignored_bool(&["USE_JDK_DEFLATER", "USE_JDK_INFLATER"], "no JDK codec"),
    Opt::bool(
        &["ASSUME_SORTED", "AS"],
        "false fails a BAM whose header is not coordinate sorted (default true).",
    ),
    Opt::ignored(&["REFERENCE_SEQUENCE", "R"], REFERENCE_UNREAD),
];

/// Why `REFERENCE_SEQUENCE` is ignored by a tool whose metrics do not read
/// it. A CRAM input is decoded with it, and then no note is written.
const REFERENCE_UNREAD: &str = "this metric does not read the reference";

/// `REFERENCE_SEQUENCE` of a tool whose input is no BAM or CRAM (a BED, a
/// VCF, a metrics file): Picard takes it, nothing reads it, and Python
/// offers no keyword for it.
pub const REFERENCE_NOT_DECODED: Opt =
    Opt::ignored(&["REFERENCE_SEQUENCE", "R"], "no BAM or CRAM is read");

/// `METRIC_ACCUMULATION_LEVEL`, the same row in every tool that has it.
pub const LEVEL_OPTION: Opt = Opt::text(
    &["METRIC_ACCUMULATION_LEVEL", "LEVEL"],
    "ALL_READS, SAMPLE, LIBRARY or READ_GROUP; each adds to the default ALL_READS, null clears.",
)
.repeated();

/// `STOP_AFTER` of a `SinglePassSamProgram`.
pub const STOP_AFTER_OPTION: Opt = Opt::int(
    &["STOP_AFTER"],
    "Read only the first N records; 0 reads them all.",
);

/// `ALIGNED_READS_ONLY` of the three by-cycle and quality tools.
pub const ALIGNED_READS_ONLY_OPTION: Opt = Opt::bool(
    &["ALIGNED_READS_ONLY"],
    "Count only aligned reads (default false).",
);

/// `REFERENCE_SEQUENCE` for a tool that reads it when given.
pub const REFERENCE_OPTION: Opt = Opt::path(
    &["REFERENCE_SEQUENCE", "R"],
    "Reference FASTA with a .fai index.",
);

/// Builds a tool's one-pass [`Collector`] from its parsed options and the
/// BAM header (Picard's `SinglePassSamProgram.setup`).
pub type CollectorFn = fn(&Invocation, &sam::Header) -> Result<Box<dyn Collector>>;

/// A registered Picard tool.
pub struct Tool {
    /// The Picard class name, e.g. `CollectInsertSizeMetrics`.
    pub name: &'static str,
    /// What the tool returns, one or two sentences: the `tools` listing,
    /// `--help` and the Python docstring.
    pub summary: &'static str,
    /// Options beyond [`COMMON`]; these shadow a common option of the same name.
    pub options: &'static [Opt],
    /// Compute the metrics: the `OUTPUT` file (`None` where Picard writes
    /// none and logs a warning, which `run` prints) and any further files.
    pub run: fn(&Invocation) -> Result<ToolOutput>,
    /// For a tool that reads the BAM once: its collector, so
    /// `CollectMultipleMetrics` can run it on a shared read. `run` is then
    /// [`Invocation::run_collector`] on it.
    pub collector: Option<CollectorFn>,
}

/// Every registered tool. One line per tool.
pub const TOOLS: &[&Tool] = &[
    &crate::insert_size::TOOL,
    &crate::alignment_summary::TOOL,
    &crate::quality_yield::TOOL,
    &crate::quality_score_distribution::TOOL,
    &crate::mean_quality_by_cycle::TOOL,
    &crate::base_distribution_by_cycle::TOOL,
    &crate::duplicate_metrics::TOOL,
    &crate::accumulate_quality_yield::TOOL,
    &crate::gc_bias::TOOL,
    &crate::multiple_metrics::TOOL,
    &crate::library_complexity::TOOL,
    &crate::wgs_metrics::TOOL,
    &crate::wgs_metrics::RAW_TOOL,
    &crate::wgs_metrics::NON_ZERO_TOOL,
    &crate::bed_to_interval_list::TOOL,
    &crate::hs_metrics::TOOL,
    &crate::hs_metrics::PCR_TOOL,
    &crate::rna_seq_metrics::TOOL,
    &crate::variant_calling_metrics::TOOL,
    &crate::variant_calling_metrics::ACCUMULATE_TOOL,
    &crate::genotype_concordance::TOOL,
    &crate::compare_metrics::TOOL,
    &crate::create_sequence_dictionary::TOOL,
    &crate::mark_duplicates::TOOL,
    &crate::sequencing_artifact::TOOL,
    &crate::artifact_to_oxog::TOOL,
];

/// The tool named `name`.
pub fn find(name: &str) -> Option<&'static Tool> {
    TOOLS.iter().copied().find(|t| t.name == name)
}

/// A parsed command line for one tool.
#[derive(Debug, Default)]
pub struct Invocation {
    pub input: PathBuf,
    pub output: Option<PathBuf>,
    /// BGZF decompression workers.
    pub threads: usize,
    /// Values by long option name, one entry per occurrence.
    values: HashMap<&'static str, Vec<String>>,
    /// `StringHeader` command line: the tool and each option as given.
    pub command: String,
    /// One line per ignored option, for standard error.
    pub notes: Vec<String>,
    /// The files the tool reads by name: `INPUT`, `REFERENCE_SEQUENCE` and
    /// every path option not marked [`Opt::written`].
    pub inputs: Vec<PathBuf>,
}

impl Invocation {
    /// The last value given for `name` (a long option name).
    pub fn get(&self, name: &str) -> Option<&str> {
        self.values.get(name)?.last().map(String::as_str)
    }

    /// Every value given for a repeatable option, in order.
    pub fn all(&self, name: &str) -> &[String] {
        self.values.get(name).map_or(&[], Vec::as_slice)
    }

    /// The parsed last value, or `default` when absent.
    pub fn parse<T: FromStr>(&self, name: &str, default: T) -> Result<T> {
        self.parse_opt(name).map(|v| v.unwrap_or(default))
    }

    /// The parsed last value, `None` when absent.
    pub fn parse_opt<T: FromStr>(&self, name: &str) -> Result<Option<T>> {
        self.get(name)
            .map(|v| {
                v.parse()
                    .map_err(|_| Error::Usage(format!("{name}={v}: not a valid value")))
            })
            .transpose()
    }

    /// A Picard collection option: each value given adds to `default`, and
    /// `null` empties the list so far.
    pub fn list<T: PartialEq>(
        &self,
        name: &str,
        default: Vec<T>,
        parse: impl Fn(&str) -> Result<T>,
    ) -> Result<Vec<T>> {
        picard_list(default, self.all(name), parse)
    }

    /// `STOP_AFTER` (0 when the tool has no such option) and `ASSUME_SORTED`
    /// (default true) as given.
    pub fn limits(&self) -> Result<Limits> {
        let stop_after: i64 = self.parse("STOP_AFTER", 0)?;
        Ok(Limits {
            stop_after: u64::try_from(stop_after).unwrap_or(0),
            assume_sorted: self.flag("ASSUME_SORTED", true)?,
        })
    }

    /// Read `INPUT` once with [`Invocation::limits`] applied, the
    /// `SinglePassSamProgram` loop.
    pub fn single_pass<A, F>(&self, build: F) -> Result<A::Output>
    where
        A: Accumulator,
        F: FnOnce(&sam::Header) -> Result<A>,
    {
        pass::run_limited(
            &self.input,
            self.reference(),
            self.threads,
            self.limits()?,
            build,
        )
    }

    /// Read `path` whole, no `STOP_AFTER` and no sort check, with the
    /// invocation's workers and reference: the read of a tool that is not
    /// a `SinglePassSamProgram`.
    pub fn read_whole<A, F>(&self, path: &Path, build: F) -> Result<A::Output>
    where
        A: Accumulator,
        F: FnOnce(&sam::Header) -> Result<A>,
    {
        pass::run_limited(
            path,
            self.reference(),
            self.threads,
            Limits::default(),
            build,
        )
    }

    /// `REFERENCE_SEQUENCE` as given, also for a tool whose metrics do not
    /// read it: a CRAM input is decoded with it.
    pub fn reference(&self) -> Option<&Path> {
        self.get("REFERENCE_SEQUENCE").map(Path::new)
    }

    /// [`Invocation::single_pass`] over a tool's collector: the `run` of
    /// every tool that has one.
    pub fn run_collector(&self, build: CollectorFn) -> Result<ToolOutput> {
        self.single_pass(|header| build(self, header))
    }

    /// Read an `INPUT` of `-` from standard input (`/dev/stdin`). The
    /// `StringHeader` command line keeps `-`. Only the verb grammar calls
    /// this: the Picard face reads a file named `-`, as Picard does.
    pub fn read_dash_from_stdin(&mut self) {
        let stdin = "/dev/stdin";
        if self.input == Path::new("-") {
            self.input = PathBuf::from(stdin);
        }
        if let Some(inputs) = self.values.get_mut("INPUT") {
            for input in inputs.iter_mut().filter(|i| *i == "-") {
                *input = stdin.to_owned();
            }
        }
    }

    /// A Picard boolean (`true` or `false`, any case).
    pub fn flag(&self, name: &str, default: bool) -> Result<bool> {
        match self.get(name).map(str::to_ascii_lowercase).as_deref() {
            None => Ok(default),
            Some("true") => Ok(true),
            Some("false") => Ok(false),
            Some(v) => Err(Error::Usage(format!("{name}={v}: expected true or false"))),
        }
    }
}

impl Tool {
    fn lookup(&self, name: &str) -> Option<&'static Opt> {
        let own_output = self.options.iter().any(|o| matches!(o.kind, Kind::Output));
        self.options
            .iter()
            .chain(
                COMMON
                    .iter()
                    .filter(|o| !(own_output && matches!(o.kind, Kind::Output))),
            )
            .find(|o| o.names.contains(&name))
    }

    /// Every option the tool takes: its own rows, then the [`COMMON`] rows
    /// it neither shadows nor rejects. Rejected rows are left out.
    pub fn table(&self) -> Vec<&'static Opt> {
        let own_output = self.options.iter().any(|o| o.kind == Kind::Output);
        let own = self.options.iter().filter(|o| o.kind != Kind::Rejected);
        let common = COMMON.iter().filter(|c| {
            !(own_output && c.kind == Kind::Output)
                && !self.options.iter().any(|o| o.names.contains(&c.picard()))
        });
        own.chain(common).collect()
    }

    /// The row a verb-grammar or Python name refers to ([`Opt::matches`]).
    pub fn option(&self, name: &str) -> Option<&'static Opt> {
        self.table().into_iter().find(|o| o.matches(name))
    }

    /// The row whose name is closest to `name`, for a did-you-mean hint.
    pub fn suggest(&self, name: &str) -> Option<&'static Opt> {
        let table = self.table();
        let names = table.iter().flat_map(|o| o.names.iter().copied());
        let best = did_you_mean(&normalize(name), names)?;
        self.option(best)
    }

    /// The Python keyword of the metrics file: `output`, or the metrics
    /// row's own keyword (`metrics_file`) for a tool whose `OUTPUT` is a
    /// value of its own (`MarkDuplicates` writes its BAM there).
    pub fn metrics_kwarg(&self) -> String {
        let own_output = self
            .options
            .iter()
            .any(|o| o.kind == Kind::Value && o.picard() == "OUTPUT");
        match self.lookup(self.output_name()) {
            Some(metrics) if own_output => metrics.kwarg(),
            _ => "output".to_owned(),
        }
    }

    /// Long name of the option that names the metrics file.
    pub fn output_name(&self) -> &'static str {
        self.options
            .iter()
            .chain(COMMON)
            .find(|o| matches!(o.kind, Kind::Output))
            .map_or("OUTPUT", |o| o.names[0])
    }

    /// Parse `args` (Picard syntax) into an [`Invocation`].
    pub fn parse(&self, args: &[String]) -> Result<Invocation> {
        let mut inv = Invocation {
            threads: default_threads(),
            command: self.name.to_owned(),
            ..Invocation::default()
        };
        let mut input = None;
        let mut rest = args.iter();
        while let Some(arg) = rest.next() {
            let dashed = arg.starts_with('-');
            let body = arg.trim_start_matches('-');
            let (name, value) = match body.split_once('=') {
                Some((n, v)) => (n, v.to_owned()),
                None if dashed => {
                    let v = rest.next().ok_or_else(|| {
                        Error::Usage(format!("{}: {arg} needs a value", self.name))
                    })?;
                    (body, v.clone())
                }
                None => {
                    return Err(Error::Usage(format!(
                        "{}: expected NAME=VALUE, got {arg:?}",
                        self.name
                    )));
                }
            };
            if name == "threads" {
                inv.threads = value
                    .parse()
                    .map_err(|_| Error::Usage(format!("threads={value}: not a number")))?;
                continue;
            }
            let opt = self.lookup(name).ok_or_else(|| {
                Error::Usage(format!(
                    "{}: option {name} is not implemented (it may change the metrics)",
                    self.name
                ))
            })?;
            if matches!(opt.kind, Kind::Rejected) {
                return Err(Error::Usage(format!(
                    "{}: Unrecognized option: {name}",
                    self.name
                )));
            }
            let given_before = match opt.kind {
                Kind::Output => inv.output.is_some(),
                Kind::Value => inv.values.contains_key(opt.names[0]),
                Kind::Ignored(_) | Kind::Rejected => false,
            };
            if given_before && !opt.repeated {
                return Err(Error::Usage(format!(
                    "Option '{}' cannot be specified more than once.",
                    opt.names[0]
                )));
            }
            inv.command.push_str(&format!(" {}={value}", opt.names[0]));
            match opt.kind {
                Kind::Ignored(why) => {
                    let typed = value == "null"
                        || match opt.ty {
                            Ty::Int => value.parse::<i64>().is_ok(),
                            Ty::Float => value.parse::<f64>().is_ok(),
                            _ => true,
                        };
                    if !typed {
                        return Err(Error::Usage(format!("{name}={value}: not a valid value")));
                    }
                    inv.notes.push(format!("ignoring {name}={value}: {why}"));
                    if opt.names[0] == "REFERENCE_SEQUENCE" {
                        inv.inputs.push(PathBuf::from(&value));
                        inv.values.entry(opt.names[0]).or_default().push(value);
                    }
                }
                Kind::Rejected => unreachable!("rejected above"),
                Kind::Output => inv.output = Some(PathBuf::from(&value)),
                Kind::Value => {
                    if opt.names[0] == "INPUT" {
                        input = Some(PathBuf::from(&value));
                    }
                    if opt.ty == Ty::Path && !opt.written && value != "null" {
                        inv.inputs.push(PathBuf::from(&value));
                    }
                    inv.values.entry(opt.names[0]).or_default().push(value);
                }
            }
        }
        inv.input = match input {
            Some(path) => path,
            None if self.takes_input() => {
                return Err(Error::Usage(format!("{}: INPUT is required", self.name)));
            }
            None => PathBuf::new(),
        };
        if inv
            .all("INPUT")
            .iter()
            .any(|path| cram::may_be_cram(Path::new(path)))
        {
            // A CRAM is decoded with the reference: it is not ignored. A
            // pipe may hold one.
            inv.notes.retain(|note| !note.ends_with(REFERENCE_UNREAD));
        }
        Ok(inv)
    }

    /// False for a tool that declares `INPUT` as [`Opt::rejected`]
    /// (`GenotypeConcordance` reads `TRUTH_VCF` and `CALL_VCF`); its
    /// positional Python arguments then fill its own value options in
    /// declaration order.
    pub fn takes_input(&self) -> bool {
        self.lookup("INPUT")
            .is_some_and(|o| matches!(o.kind, Kind::Value))
    }

    /// The long names of the tool's own value options, in declaration
    /// order.
    pub fn value_options(&self) -> impl Iterator<Item = &'static str> {
        self.options
            .iter()
            .filter(|o| matches!(o.kind, Kind::Value))
            .map(|o| o.names[0])
    }

    /// Parse and compute. The invocation is returned for its output path
    /// and header.
    pub fn execute(&self, args: &[String]) -> Result<(Invocation, ToolOutput)> {
        let inv = self.parse(args)?;
        let output = (self.run)(&inv)?;
        Ok((inv, output))
    }

    /// The command-line entry: compute and write the Picard metrics file to
    /// `OUTPUT`, or standard output without one, and every further file to
    /// its own path.
    pub fn run_cli(&self, args: &[String]) -> Result<()> {
        let (inv, output) = self.execute(args)?;
        for note in &inv.notes {
            eprintln!("bammetrics: {note}");
        }
        write_files(&inv, &output)
    }
}

/// Write a tool's files: `OUTPUT` (standard output without a path) and the
/// extra files, each with the invocation's `StringHeader` lines.
pub fn write_files(inv: &Invocation, output: &ToolOutput) -> Result<()> {
    let targets = output
        .main
        .as_ref()
        .and(inv.output.as_deref())
        .into_iter()
        .chain(output.extra.iter().map(|(path, _)| path.as_path()));
    check_writable(targets, &inv.inputs)?;
    if let Some(main) = &output.main {
        let path = inv.output.as_deref();
        write_to(path, |out| main.write(main.headers(&inv.command), out))?;
    }
    for (path, file) in &output.extra {
        write_to(Some(path), |out| {
            file.write(file.headers(&inv.command), out)
        })?;
    }
    Ok(())
}

/// Fail before any file is written when one of `targets` cannot be
/// written, as Picard's `assertFileIsWritable` does for every output first:
/// a directory, a file without write permission, a missing parent
/// directory. A target that is one of `inputs` fails too, so a run never
/// replaces a file it reads.
pub fn check_writable<'a>(
    targets: impl IntoIterator<Item = &'a Path>,
    inputs: &[PathBuf],
) -> Result<()> {
    let inputs: Vec<PathBuf> = inputs
        .iter()
        .filter_map(|p| p.canonicalize().ok())
        .collect();
    for target in targets {
        let output_error = |source| Error::Output {
            path: target.to_path_buf(),
            source,
        };
        match std::fs::metadata(target) {
            Ok(meta) if meta.is_dir() || meta.is_file() => {
                if target.canonicalize().is_ok_and(|t| inputs.contains(&t)) {
                    return Err(Error::Usage(format!(
                        "{} is an input; writing it would destroy it",
                        target.display()
                    )));
                }
                // Opening for write without truncation changes nothing; a
                // directory fails with the error `File::create` would give.
                std::fs::OpenOptions::new()
                    .write(true)
                    .open(target)
                    .map_err(output_error)?;
            }
            // A device or pipe (`/dev/stdout`) is opened once, when written.
            Ok(_) => {}
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                if let Some(parent) = target.parent().filter(|p| !p.as_os_str().is_empty()) {
                    std::fs::metadata(parent).map_err(output_error)?;
                }
            }
            Err(error) => return Err(output_error(error)),
        }
    }
    Ok(())
}

/// The values of a Picard collection option applied to its default: each
/// value is added, `null` empties the list so far. A value already in the
/// list is added again; the options that use this only test membership.
pub fn picard_list<T, S: AsRef<str>>(
    default: Vec<T>,
    values: &[S],
    parse: impl Fn(&str) -> Result<T>,
) -> Result<Vec<T>> {
    let mut list = default;
    for value in values {
        if value.as_ref() == "null" {
            list.clear();
        } else {
            list.push(parse(value.as_ref())?);
        }
    }
    Ok(list)
}

/// Run `write` against `path`, or standard output when `path` is `None`.
/// A file is written under a temporary name and renamed over `path` once
/// complete ([`crate::cleanup::AtomicFile`]), so a failed write leaves the
/// old file, or none.
pub fn write_to(
    path: Option<&Path>,
    write: impl FnOnce(&mut BufWriter<Box<dyn Write>>) -> io::Result<()>,
) -> Result<()> {
    let output_error = |source| Error::Output {
        path: path.map_or_else(|| PathBuf::from("-"), Path::to_path_buf),
        source,
    };
    let file = path.map(crate::cleanup::AtomicFile::new);
    let sink: Box<dyn Write> = match &file {
        Some(file) => Box::new(file.create().map_err(output_error)?),
        None => Box::new(io::stdout().lock()),
    };
    let mut out = BufWriter::new(sink);
    write(&mut out).map_err(output_error)?;
    out.flush().map_err(output_error)?;
    drop(out);
    match file {
        Some(file) => file.commit().map_err(output_error),
        None => Ok(()),
    }
}

/// The candidate closest to `word` by edit distance, ignoring case, when it
/// is close enough to be a typo: at most 2 edits, or a third of the word.
pub fn did_you_mean<'a>(
    word: &str,
    candidates: impl IntoIterator<Item = &'a str>,
) -> Option<&'a str> {
    let word = word.to_ascii_lowercase();
    let limit = (word.chars().count() / 3).max(2);
    candidates
        .into_iter()
        .map(|c| (edit_distance(&word, &c.to_ascii_lowercase()), c))
        .filter(|(d, _)| *d <= limit)
        .min_by_key(|(d, _)| *d)
        .map(|(_, c)| c)
}

/// Levenshtein distance over chars.
fn edit_distance(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut row: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut diagonal = row[0];
        row[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let above = row[j + 1];
            row[j + 1] = (above + 1)
                .min(row[j] + 1)
                .min(diagonal + usize::from(ca != *cb));
            diagonal = above;
        }
    }
    row[b.len()]
}

/// `CollectInsertSizeMetrics` -> `collect_insert_size_metrics`, the Python
/// function name of a tool.
pub fn snake_case(name: &str) -> String {
    let mut out = String::new();
    for (i, c) in name.chars().enumerate() {
        if c.is_ascii_uppercase() {
            if i > 0 {
                out.push('_');
            }
            out.push(c.to_ascii_lowercase());
        } else {
            out.push(c);
        }
    }
    out
}

/// First line of the generated block in `bammetrics.pyi`.
pub const STUB_BEGIN: &str = "# BEGIN registry tools: generated from crates/core/src/registry.rs (UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs)";
/// Last line of the generated block.
pub const STUB_END: &str = "# END registry tools";

/// The Python annotation of one row's keyword argument.
fn python_annotation(opt: &Opt) -> String {
    let one = opt.ty.python();
    if opt.repeated {
        format!("{one} | Sequence[{one}]")
    } else {
        one.to_owned()
    }
}

/// `text` wrapped to lines of at most `width` characters, each indented.
pub fn wrap(text: &str, indent: &str, width: usize) -> Vec<String> {
    let mut lines = Vec::new();
    let mut line = String::new();
    for word in text.split_whitespace() {
        if !line.is_empty() && indent.len() + line.len() + 1 + word.len() > width {
            lines.push(format!("{indent}{line}"));
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    if !line.is_empty() {
        lines.push(format!("{indent}{line}"));
    }
    lines
}

/// The `bammetrics.pyi` entry of one tool: its positional inputs, one
/// keyword argument per value row of [`Tool::table`] (`output` for the
/// metrics file), and the summary as the docstring. `None` leaves an
/// option out, so Picard's default applies.
pub fn python_stub(tool: &Tool) -> String {
    let table = tool.table();
    let mut params = Vec::new();
    let positional: Vec<&Opt> = if tool.takes_input() {
        table
            .iter()
            .copied()
            .filter(|o| o.picard() == "INPUT")
            .collect()
    } else {
        table
            .iter()
            .copied()
            .filter(|o| o.kind == Kind::Value && o.required)
            .collect()
    };
    for opt in &positional {
        let ty = if opt.picard() == "INPUT" && opt.repeated {
            "_Paths".to_owned()
        } else {
            python_annotation(opt)
        };
        params.push(format!("{}: {ty}", opt.kwarg()));
    }
    let mut keywords = Vec::new();
    for opt in &table {
        if positional.iter().any(|p| p.picard() == opt.picard()) {
            continue;
        }
        match opt.kind {
            Kind::Output if opt.required => {
                keywords.push(format!("{}: _Path", tool.metrics_kwarg()));
            }
            Kind::Output => keywords.push(format!("{}: _Path | None = None", tool.metrics_kwarg())),
            Kind::Value if opt.required => {
                keywords.push(format!("{}: {}", opt.kwarg(), python_annotation(opt)));
            }
            Kind::Value => keywords.push(format!(
                "{}: {} | None = None",
                opt.kwarg(),
                python_annotation(opt)
            )),
            // The metrics do not read it, a CRAM input is decoded with it.
            Kind::Ignored(REFERENCE_UNREAD) if opt.picard() == "REFERENCE_SEQUENCE" => {
                keywords.push(format!("{}: _Path | None = None", opt.kwarg()));
            }
            Kind::Ignored(_) | Kind::Rejected => {}
        }
    }
    if !keywords.is_empty() {
        params.push("*".to_owned());
        params.extend(keywords);
    }
    let mut text = format!("def {}(\n", snake_case(tool.name));
    for param in params {
        text.push_str(&format!("    {param},\n"));
    }
    text.push_str(") -> list[dict[str, Any]]:\n");
    let mut doc = wrap(&format!("{}: {}", tool.name, tool.summary), "    ", 88);
    doc[0] = format!("    \"\"\"{}", doc[0].trim_start());
    let last = doc.len() - 1;
    doc[last].push_str("\"\"\"");
    for line in doc {
        text.push_str(&line);
        text.push('\n');
    }
    text
}

/// The generated block of `bammetrics.pyi`: every tool in [`TOOLS`] order,
/// between [`STUB_BEGIN`] and [`STUB_END`].
pub fn python_stubs() -> String {
    let mut text = format!(
        "{STUB_BEGIN}\n# Each returns its rows as dicts keyed by Picard's column names; output= also \
         writes the files.\n"
    );
    for tool in TOOLS {
        text.push('\n');
        text.push_str(&python_stub(tool));
    }
    text.push('\n');
    text.push_str(STUB_END);
    text.push('\n');
    text
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn tool() -> &'static Tool {
        find("CollectInsertSizeMetrics").expect("registered")
    }

    #[test]
    fn picard_and_dashed_forms_agree() {
        let forms = [
            strings(&["I=x.bam", "O=out.txt", "M=0.5"]),
            strings(&["--INPUT", "x.bam", "--OUTPUT=out.txt", "-M", "0.5"]),
            strings(&["--I=x.bam", "-O", "out.txt", "--MINIMUM_PCT", "0.5"]),
        ];
        for args in forms {
            let inv = tool().parse(&args).expect("parses");
            assert_eq!(inv.input, Path::new("x.bam"));
            assert_eq!(inv.output.as_deref(), Some(Path::new("out.txt")));
            assert_eq!(inv.get("MINIMUM_PCT"), Some("0.5"));
            assert_eq!(
                inv.command,
                "CollectInsertSizeMetrics INPUT=x.bam OUTPUT=out.txt MINIMUM_PCT=0.5"
            );
        }
    }

    #[test]
    fn repeated_options_keep_every_value() {
        let inv = tool()
            .parse(&strings(&[
                "I=x.bam",
                "LEVEL=SAMPLE",
                "METRIC_ACCUMULATION_LEVEL=READ_GROUP",
            ]))
            .expect("parses");
        assert_eq!(
            inv.all("METRIC_ACCUMULATION_LEVEL"),
            ["SAMPLE", "READ_GROUP"]
        );
    }

    /// Picard's parser: `Option 'MINIMUM_PCT' cannot be specified more than once.`
    #[test]
    fn a_single_option_given_twice_is_an_error() {
        for args in [
            ["I=x.bam", "M=0.5", "MINIMUM_PCT=0.1"],
            ["I=x.bam", "I=y.bam", "M=0.5"],
            ["I=x.bam", "O=a.txt", "OUTPUT=b.txt"],
        ] {
            let err = tool().parse(&strings(&args)).expect_err("given twice");
            assert!(
                matches!(err, Error::Usage(_)) && err.to_string().contains("more than once"),
                "{args:?}: {err}"
            );
        }
        let inv = tool()
            .parse(&strings(&["I=x.bam", "QUIET=true", "QUIET=false"]))
            .expect("an ignored option may repeat (Picard's TMP_DIR is a list)");
        assert_eq!(inv.notes.len(), 2);
    }

    #[test]
    fn plumbing_is_noted_and_unknown_options_are_errors() {
        let inv = tool()
            .parse(&strings(&["I=x.bam", "QUIET=true", "R=ref.fa"]))
            .expect("parses");
        assert_eq!(inv.notes.len(), 2);
        let err = tool()
            .parse(&strings(&["I=x.bam", "STOP_AFTER=10"]))
            .expect_err("unknown");
        assert!(err.to_string().contains("STOP_AFTER"));
        let err = tool().parse(&strings(&["bam"])).expect_err("bare word");
        assert!(err.to_string().contains("NAME=VALUE"));
        let err = tool().parse(&strings(&["O=o.txt"])).expect_err("no input");
        assert!(err.to_string().contains("INPUT"));
    }

    #[test]
    fn typed_getters_name_the_bad_option() {
        let inv = tool()
            .parse(&strings(&["I=x.bam", "INCLUDE_DUPLICATES=maybe", "W=abc"]))
            .expect("parses");
        assert!(inv.flag("INCLUDE_DUPLICATES", false).is_err());
        let err = inv
            .parse_opt::<i64>("HISTOGRAM_WIDTH")
            .expect_err("not a number");
        assert!(err.to_string().contains("HISTOGRAM_WIDTH"));
        assert!(inv.flag("MISSING", true).expect("default"));
    }

    /// `bammetrics.pyi` holds [`python_stubs`] between its markers, so the
    /// stub and the option table cannot disagree. `UPDATE_PYI=1` rewrites
    /// the block.
    #[test]
    fn python_stubs() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../bammetrics.pyi");
        let text = std::fs::read_to_string(&path).expect("bammetrics.pyi");
        let begin = text.find(STUB_BEGIN).expect("begin marker");
        let end = text.find(STUB_END).expect("end marker") + STUB_END.len() + 1;
        let generated = super::python_stubs();
        if std::env::var_os("UPDATE_PYI").is_some() {
            let updated = format!("{}{generated}{}", &text[..begin], &text[end..]);
            std::fs::write(&path, updated).expect("write bammetrics.pyi");
            return;
        }
        assert_eq!(
            &text[begin..end],
            generated,
            "bammetrics.pyi is stale: UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs"
        );
    }

    #[test]
    fn path_options_that_name_outputs_are_written() {
        for tool in TOOLS {
            for opt in tool.table() {
                let name = opt.picard();
                let looks_written = name.contains("OUT") || name.starts_with("PER_");
                if opt.kind == Kind::Value && opt.ty == Ty::Path {
                    assert_eq!(opt.written, looks_written, "{}: {name}", tool.name);
                }
            }
        }
    }

    #[test]
    fn an_output_that_is_an_input_or_a_directory_is_refused_first() {
        let dir = std::env::temp_dir().join(format!("bammetrics-writable-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("in.txt");
        std::fs::write(&input, "keep").unwrap();
        let fresh = dir.join("fresh.txt");
        let error = check_writable(
            [fresh.as_path(), input.as_path()],
            std::slice::from_ref(&input),
        );
        assert!(error.unwrap_err().to_string().contains("is an input"));
        assert_eq!(std::fs::read_to_string(&input).unwrap(), "keep");
        let error = check_writable([dir.as_path()], &[]).unwrap_err();
        assert!(matches!(error, Error::Output { .. }), "{error}");
        let error = check_writable([dir.join("no/such/file").as_path()], &[]).unwrap_err();
        assert!(matches!(error, Error::Output { .. }), "{error}");
        assert!(check_writable([fresh.as_path(), Path::new("/dev/stdout")], &[]).is_ok());
        assert!(!fresh.exists());
        std::fs::remove_file(&input).unwrap();
        std::fs::remove_dir(&dir).unwrap();
    }

    #[test]
    fn every_row_has_help_and_one_meaning_per_name() {
        for tool in TOOLS {
            assert!(!tool.summary.is_empty(), "{}: no summary", tool.name);
            let table = tool.table();
            for opt in &table {
                assert!(
                    !opt.help.is_empty(),
                    "{}: {} has no help",
                    tool.name,
                    opt.picard()
                );
                for name in opt.names {
                    let rows = table.iter().filter(|o| o.matches(name)).count();
                    assert_eq!(rows, 1, "{}: {name} names {rows} rows", tool.name);
                    // The Picard face resolves the name to the same row.
                    assert_eq!(
                        tool.lookup(name).map(Opt::picard),
                        Some(opt.picard()),
                        "{}: {name}",
                        tool.name
                    );
                }
            }
        }
    }

    #[test]
    fn spellings_name_one_row() {
        let tool = tool();
        for spelling in ["minimum-pct", "MINIMUM_PCT", "Minimum_Pct", "M", "m"] {
            assert_eq!(tool.option(spelling).map(Opt::picard), Some("MINIMUM_PCT"));
        }
        let pct = tool.option("M").expect("row");
        assert_eq!(pct.kebab(), "minimum-pct");
        assert_eq!(pct.kwarg(), "minimum_pct");
        assert_eq!(pct.ty, Ty::Float);
        assert_eq!(
            tool.option("LEVEL").map(|o| (o.picard(), o.repeated)),
            Some(("METRIC_ACCUMULATION_LEVEL", true))
        );
        assert!(
            tool.option("STOP_AFTER").is_none(),
            "not a CollectInsertSizeMetrics option"
        );
        let dup = find("CollectDuplicateMetrics").expect("registered");
        assert!(
            dup.option("o").is_none(),
            "OUTPUT is replaced by METRICS_FILE"
        );
        assert_eq!(dup.option("m").map(Opt::picard), Some("METRICS_FILE"));
    }

    #[test]
    fn did_you_mean_finds_the_typo_only() {
        let tool = tool();
        assert_eq!(
            tool.suggest("minmum-pct").map(Opt::picard),
            Some("MINIMUM_PCT")
        );
        assert_eq!(
            tool.suggest("deviation").map(Opt::picard),
            Some("DEVIATIONS")
        );
        assert_eq!(
            tool.suggest("histogram-widht").map(Opt::picard),
            Some("HISTOGRAM_WIDTH")
        );
        assert!(tool.suggest("reference-sequence-file-name").is_none());
        assert_eq!(did_you_mean("tols", ["tools", "presets"]), Some("tools"));
        assert_eq!(did_you_mean("xyz", ["tools", "presets"]), None);
    }

    #[test]
    fn stdin_dash_is_read_from_dev_stdin_with_the_header_kept() {
        let mut inv = tool().parse(&strings(&["I=-"])).expect("parses");
        inv.read_dash_from_stdin();
        assert_eq!(inv.input, Path::new("/dev/stdin"));
        assert_eq!(inv.all("INPUT"), ["/dev/stdin"]);
        assert_eq!(inv.command, "CollectInsertSizeMetrics INPUT=-");
    }

    #[test]
    fn a_tool_with_its_own_output_option_takes_no_output() {
        let dup = find("CollectDuplicateMetrics").expect("registered");
        assert_eq!(dup.output_name(), "METRICS_FILE");
        assert_eq!(tool().output_name(), "OUTPUT");
        let inv = dup
            .parse(&strings(&["I=x.bam", "M=d.txt", "AS=false"]))
            .expect("parses");
        assert_eq!(inv.output.as_deref(), Some(Path::new("d.txt")));
        assert!(!inv.flag("ASSUME_SORTED", true).expect("bool"));
        let err = dup
            .parse(&strings(&["I=x.bam", "O=d.txt"]))
            .expect_err("no OUTPUT");
        assert!(err.to_string().contains("option O "));
    }
}

//! The tool verbs (docs/05 section 2.1): one kebab-case verb per registry
//! tool, parsed by [`grammar`] and run through the registry, so a verb and
//! the Picard class it names compute the same thing.
//!
//! A verb's flags are the tool's option table rows plus the shared options
//! of section 2.2. Output naming follows section 2.3: a single-file verb
//! writes `-o FILE` (standard output without it, or with `-o -`); a
//! multi-file verb writes under `-o PREFIX`, `./<input stem>` by default,
//! with the names of [`file_extensions`]. `-` as the input reads standard
//! input.
//!
//! The four tool verbs that were 0.12 subcommands (`insert-size`,
//! `alignment-summary`, `quality-yield`, `hs-metrics`) keep their flags and
//! the `StringHeader` line they wrote, which lists every option with its
//! value.

use std::ffi::OsString;
use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::pass;
use bammetrics_core::registry::{self, Invocation, Kind, Tool, Ty};

use crate::error::VerbError;
use crate::format::{self, Context, Format};
use crate::grammar::{self, Flag, Key, Parsed};

/// How a verb names what it writes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Files {
    /// `-o FILE`, or standard output.
    Single,
    /// `-o PREFIX`, default `./<input stem>`.
    Prefix,
}

/// Builds the `StringHeader` command line of a 0.12 subcommand.
type Header = fn(&Invocation) -> bammetrics_core::Result<String>;

/// One tool verb.
pub struct Verb {
    pub name: &'static str,
    /// The Picard class it runs.
    pub class: &'static str,
    pub files: Files,
    /// What the input is, for `--help` and errors: `BAM`, `VCF`, ...
    pub input: &'static str,
    /// 0.12 subcommand spellings and the Picard name each stands for.
    pub aliases: &'static [(&'static str, &'static str)],
    /// The `StringHeader` line of a 0.12 subcommand.
    pub header: Option<Header>,
}

const fn verb(name: &'static str, class: &'static str, files: Files) -> Verb {
    Verb {
        name,
        class,
        files,
        input: "BAM",
        aliases: &[],
        header: None,
    }
}

const fn with_input(mut verb: Verb, input: &'static str) -> Verb {
    verb.input = input;
    verb
}

/// The verb table. The names of docs/05 section 2.1 first, then one verb
/// per other registry tool, named by the rule: the Picard class without
/// `Collect` and `Metrics`, in kebab case.
pub const VERBS: &[Verb] = &[
    Verb {
        aliases: &[("MIN_HISTOGRAM_WIDTH", "MINIMUM_HISTOGRAM_WIDTH")],
        header: Some(bammetrics_core::insert_size::subcommand_header),
        ..verb("insert-size", "CollectInsertSizeMetrics", Files::Single)
    },
    Verb {
        aliases: &[
            ("EXPECTED_PAIR_ORIENTATION", "EXPECTED_PAIR_ORIENTATIONS"),
            ("REFERENCE", "REFERENCE_SEQUENCE"),
        ],
        header: Some(bammetrics_core::alignment_summary::subcommand_header),
        ..verb(
            "alignment-summary",
            "CollectAlignmentSummaryMetrics",
            Files::Single,
        )
    },
    Verb {
        header: Some(bammetrics_core::quality_yield::subcommand_header),
        ..verb("quality-yield", "CollectQualityYieldMetrics", Files::Single)
    },
    Verb {
        aliases: &[
            ("BAITS", "BAIT_INTERVALS"),
            ("TARGETS", "TARGET_INTERVALS"),
            ("REFERENCE", "REFERENCE_SEQUENCE"),
        ],
        header: Some(bammetrics_core::hs_metrics::subcommand_header),
        ..verb("hs-metrics", "CollectHsMetrics", Files::Single)
    },
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb("gc-bias", "CollectGcBiasMetrics", Files::Prefix)
    },
    verb(
        "quality-score-distribution",
        "QualityScoreDistribution",
        Files::Single,
    ),
    verb("mean-quality-by-cycle", "MeanQualityByCycle", Files::Single),
    verb(
        "base-distribution-by-cycle",
        "CollectBaseDistributionByCycle",
        Files::Single,
    ),
    verb(
        "duplicate-metrics",
        "CollectDuplicateMetrics",
        Files::Single,
    ),
    verb(
        "library-complexity",
        "EstimateLibraryComplexity",
        Files::Single,
    ),
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb("multiple-metrics", "CollectMultipleMetrics", Files::Prefix)
    },
    with_input(
        verb(
            "accumulate-quality-yield",
            "AccumulateQualityYieldMetrics",
            Files::Single,
        ),
        "FILE",
    ),
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb("wgs", "CollectWgsMetrics", Files::Single)
    },
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb("raw-wgs", "CollectRawWgsMetrics", Files::Single)
    },
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb(
            "wgs-with-non-zero-coverage",
            "CollectWgsMetricsWithNonZeroCoverage",
            Files::Single,
        )
    },
    Verb {
        aliases: &[("REFERENCE", "REFERENCE_SEQUENCE")],
        ..verb("targeted-pcr", "CollectTargetedPcrMetrics", Files::Single)
    },
    verb("rna-seq", "CollectRnaSeqMetrics", Files::Single),
    with_input(
        verb(
            "variant-calling",
            "CollectVariantCallingMetrics",
            Files::Prefix,
        ),
        "VCF",
    ),
    with_input(
        verb(
            "accumulate-variant-calling",
            "AccumulateVariantCallingMetrics",
            Files::Prefix,
        ),
        "PREFIX",
    ),
    with_input(
        verb("genotype-concordance", "GenotypeConcordance", Files::Prefix),
        "VCF",
    ),
    with_input(
        verb("bed-to-interval-list", "BedToIntervalList", Files::Single),
        "BED",
    ),
];

/// The verb named `name`.
pub fn find(name: &str) -> Option<&'static Verb> {
    VERBS.iter().find(|v| v.name == name)
}

/// The verb that runs Picard class `class`.
pub fn for_class(class: &str) -> Option<&'static Verb> {
    VERBS.iter().find(|v| v.class == class)
}

/// The shared flags of section 2.2 every verb takes.
pub fn shared_flags() -> Vec<Flag> {
    vec![
        Flag::new(
            Key::Threads,
            "threads",
            &["@"],
            Ty::Int,
            &format!(
                "BGZF decompression threads (default: available cores, at most 8; {} here).",
                pass::default_threads()
            ),
        )
        .names(&["THREADS"])
        .placeholder("N"),
        Flag::new(
            Key::Quiet,
            "quiet",
            &["q"],
            Ty::Bool,
            "No notes on standard error.",
        ),
        Flag::new(
            Key::Verbose,
            "verbose",
            &["v"],
            Ty::Bool,
            "Print the Picard command line and each file written on standard error.",
        ),
        Flag::new(Key::Help, "help", &["h"], Ty::Bool, "This help."),
    ]
}

impl Verb {
    pub fn tool(&self) -> &'static Tool {
        registry::find(self.class).expect("every verb names a registered tool")
    }

    /// The verb's flags: input, output, `--format` and `--sample`, every
    /// row of the tool's option table, then the shared flags.
    pub fn flags(&self) -> Vec<Flag> {
        let tool = self.tool();
        let table = tool.table();
        let mut flags = Vec::new();
        if let Some(input) = table.iter().find(|o| o.picard() == "INPUT") {
            let mut flag = Flag::new(
                Key::Input,
                "input",
                &["i"],
                Ty::Path,
                &format!(
                    "{} The positional argument; - reads standard input.",
                    input.help
                ),
            )
            .names(input.names)
            .placeholder(self.input);
            flag.repeated = input.repeated;
            flags.push(flag);
        }
        let output_names: Vec<&'static str> = table
            .iter()
            .filter(|o| o.kind == Kind::Output)
            .flat_map(|o| o.names.iter().copied())
            .chain(["OUTPUT", "O"])
            .collect();
        let (placeholder, help) = match self.files {
            Files::Single => (
                "FILE",
                "The metrics file; standard output without one, or with -.",
            ),
            Files::Prefix => (
                "PREFIX",
                "Output prefix (default ./ plus the input name without its extension); a \
                 directory gets that name as the file part.",
            ),
        };
        flags.push(
            Flag::new(Key::Output, "output", &["o", "O"], Ty::Path, help)
                .names(&output_names)
                .placeholder(placeholder),
        );
        if self.files == Files::Single {
            flags.push(
                Flag::new(
                    Key::Format,
                    "format",
                    &[],
                    Ty::Text,
                    "picard (the metrics file), json, tsv, jsonl or human.",
                )
                .placeholder("FORMAT"),
            );
            flags.push(
                Flag::new(
                    Key::Sample,
                    "sample",
                    &[],
                    Ty::Text,
                    "Sample name for the json and jsonl forms; never written into a Picard file.",
                )
                .placeholder("NAME"),
            );
        }
        if self.name == "hs-metrics" {
            flags.push(
                Flag::new(
                    Key::Own("panel"),
                    "panel",
                    &[],
                    Ty::Path,
                    "Baits and targets in one file: --bait-intervals FILE --target-intervals FILE.",
                )
                .placeholder("FILE"),
            );
        }
        let shared = shared_flags();
        for opt in table {
            // `QUIET=true` reads as the shared `--quiet`.
            let is_shared = shared.iter().any(|f| f.long == opt.kebab());
            if opt.picard() == "INPUT" || opt.kind == Kind::Output || is_shared {
                continue;
            }
            let mut names: Vec<&'static str> = opt.names.to_vec();
            names.extend(
                self.aliases
                    .iter()
                    .filter(|(_, target)| *target == opt.picard())
                    .map(|(alias, _)| *alias),
            );
            let help = match opt.kind {
                Kind::Ignored(why) => format!("Accepted and ignored: {why}."),
                _ if opt.required => format!("{} Required.", opt.help),
                _ => opt.help.to_owned(),
            };
            let ty = if matches!(opt.kind, Kind::Ignored(_)) {
                Ty::Text
            } else {
                opt.ty
            };
            let mut flag = Flag::new(Key::Picard(opt.picard()), &opt.kebab(), &[], ty, &help)
                .names(&names)
                .placeholder(ty.placeholder());
            flag.repeated = opt.repeated;
            flag.required = opt.required;
            flags.push(flag);
        }
        flags.extend(shared);
        flags
    }

    /// `bammetrics <verb> --help`.
    pub fn help(&self) -> String {
        let tool = self.tool();
        let mut text = format!(
            "bammetrics {}: {}, {}\n\nUsage: bammetrics {} [{}] [OPTIONS]\n       \
             bammetrics {} {}=... [KEY=VALUE ...]   (Picard syntax)\n\n",
            self.name,
            tool.name,
            tool.summary,
            self.name,
            self.input,
            tool.name,
            if tool.takes_input() { "I" } else { "KEY" },
        );
        text.push_str(
            "Every option also reads as --PICARD_KEY VALUE and PICARD_KEY=VALUE. A boolean \
             takes --x, --no-x, --x=false or --x false.\n\nOptions:\n",
        );
        text.push_str(&option_lines(&self.flags()));
        text
    }

    /// Parse `args` and run the tool. The exit code is 0; errors carry
    /// theirs (2 for usage, 1 for inputs and failed tools).
    pub fn run(&self, args: &[String]) -> Result<u8, VerbError> {
        let flags = self.flags();
        let parsed = grammar::parse(self.name, &flags, args)?;
        if parsed.flag(Key::Help) {
            print!("{}", self.help());
            return Ok(0);
        }
        let tool = self.tool();
        let format = match parsed.last(Key::Format) {
            None => Format::Picard,
            Some(given) => Format::parse(&given.values[0]).ok_or_else(|| {
                let mut error =
                    VerbError::usage(self.name, format!("no format {:?}", given.values[0]))
                        .from(given.spelled.clone());
                if let Some(best) = registry::did_you_mean(&given.values[0], format::NAMES) {
                    error = error.hint(format!("did you mean {best}?"));
                }
                error.hint(format!("formats: {}", format::NAMES.join(", ")))
            })?,
        };
        let mut picard = self.picard_arguments(&parsed, &flags)?;
        let output = self.output_path(&parsed, &picard)?;
        if let Some(path) = &output {
            create_parent(self.name, path)?;
            self.name_outputs(path, &parsed, &mut picard, format);
        }
        self.check_required(&flags, &picard)?;
        self.check_inputs(&picard)?;
        let mut inv = tool
            .parse(&picard)
            .map_err(|e| VerbError::core(self.name, &e))?;
        inv.threads = self.threads(&parsed)?;
        if let Some(header) = self.header {
            inv.command = header(&inv).map_err(|e| VerbError::core(self.name, &e))?;
        }
        inv.read_dash_from_stdin();
        let quiet = parsed.flag(Key::Quiet);
        let verbose = parsed.flag(Key::Verbose);
        if verbose {
            eprintln!("bammetrics {}: {}", self.name, inv.command);
        }
        let result = (tool.run)(&inv).map_err(|e| VerbError::core(self.name, &e))?;
        if !quiet {
            for note in &inv.notes {
                eprintln!("bammetrics {}: {note}", self.name);
            }
        }
        let sample = parsed.last(Key::Sample).map(|g| g.values[0].as_str());
        let context = Context {
            tool: tool.name,
            command: &inv.command,
            sample,
        };
        let written = if self.files == Files::Single {
            write_single(&inv, &result, format, &context, output.as_deref())
        } else {
            registry::write_files(&inv, &result).map(|()| written_paths(&inv, &result))
        }
        .map_err(|e| VerbError::core(self.name, &e))?;
        if verbose {
            for path in written {
                eprintln!("bammetrics {}: wrote {}", self.name, path.display());
            }
        }
        Ok(0)
    }

    /// The Picard arguments of every option given, in order: `KEY=VALUE`
    /// with the row's long name. Positional arguments are the inputs
    /// (`INPUT`, or the tool's first required options when it has no
    /// `INPUT`).
    fn picard_arguments(&self, parsed: &Parsed, flags: &[Flag]) -> Result<Vec<String>, VerbError> {
        let tool = self.tool();
        let mut picard = Vec::new();
        let inputs: Vec<&String> = parsed
            .positional
            .iter()
            .chain(parsed.all(Key::Input).flat_map(|g| &g.values))
            .collect();
        if tool.takes_input() {
            let many = flags.iter().any(|f| f.key == Key::Input && f.repeated);
            if inputs.len() > 1 && !many {
                return Err(VerbError::usage(
                    self.name,
                    format!("one {} expected, got {}", self.input, inputs.len()),
                )
                .from(inputs[1].clone()));
            }
            picard.extend(inputs.iter().map(|i| format!("INPUT={i}")));
        } else {
            let slots: Vec<&str> = flags
                .iter()
                .filter(|f| f.required)
                .filter_map(Flag::picard)
                .collect();
            if inputs.len() > slots.len() {
                return Err(VerbError::usage(
                    self.name,
                    format!(
                        "{} positional arguments, at most {}",
                        inputs.len(),
                        slots.len()
                    ),
                )
                .hint(format!("they are {}", slots.join(", "))));
            }
            picard.extend(
                slots
                    .iter()
                    .zip(&inputs)
                    .map(|(slot, i)| format!("{slot}={i}")),
            );
        }
        for given in &parsed.given {
            match given.key {
                Key::Picard(name) => {
                    picard.extend(given.values.iter().map(|v| format!("{name}={v}")));
                }
                Key::Own("panel") => {
                    let file = &given.values[0];
                    picard.push(format!("BAIT_INTERVALS={file}"));
                    picard.push(format!("TARGET_INTERVALS={file}"));
                }
                _ => {}
            }
        }
        Ok(picard)
    }

    /// The `-o` path: `None` for standard output (single-file verbs
    /// without `-o`, or `-o -`). A prefix verb defaults to `./<stem>`, and
    /// a directory gets the stem as the file part.
    fn output_path(
        &self,
        parsed: &Parsed,
        picard: &[String],
    ) -> Result<Option<PathBuf>, VerbError> {
        let given = parsed.last(Key::Output);
        match self.files {
            Files::Single => Ok(given
                .map(|g| PathBuf::from(&g.values[0]))
                .filter(|p| p != Path::new("-"))),
            Files::Prefix => {
                if let Some(g) = given.filter(|g| g.values[0] == "-") {
                    return Err(VerbError::usage(
                        self.name,
                        "standard output takes one file; this verb writes several",
                    )
                    .from(g.spelled.clone())
                    .hint("give -o PREFIX"));
                }
                let stem = picard
                    .iter()
                    .find_map(|a| {
                        a.strip_prefix("INPUT=")
                            .or_else(|| a.strip_prefix("CALL_VCF="))
                    })
                    .map(stem);
                let path = match (given, stem) {
                    (Some(g), stem) => {
                        let path = PathBuf::from(&g.values[0]);
                        let directory = g.values[0].ends_with('/') || path.is_dir();
                        match stem {
                            Some(stem) if directory => path.join(stem),
                            _ => path,
                        }
                    }
                    (None, Some(stem)) => {
                        let path = Path::new(".").join(&stem);
                        if !parsed.flag(Key::Quiet) {
                            eprintln!(
                                "bammetrics {}: no -o, files go to {}.*",
                                self.name,
                                path.display()
                            );
                        }
                        path
                    }
                    (None, None) => {
                        return Err(VerbError::usage(self.name, "missing -o PREFIX"));
                    }
                };
                Ok(Some(path))
            }
        }
    }

    /// The tool's file options for output `path`: `OUTPUT` for a
    /// single-file verb in the Picard form (another form is written by the
    /// verb itself), the prefix as `OUTPUT` for a tool whose Picard
    /// `OUTPUT` is a prefix, or one name per file of [`file_extensions`]
    /// (an option the user gave keeps its value).
    fn name_outputs(&self, path: &Path, parsed: &Parsed, picard: &mut Vec<String>, format: Format) {
        let tool = self.tool();
        let output = tool.output_name();
        if self.files == Files::Single {
            if format == Format::Picard {
                picard.push(format!("{output}={}", path.display()));
            }
            return;
        }
        let Some(files) = file_extensions(tool.name) else {
            picard.push(format!("{output}={}", path.display()));
            return;
        };
        for (option, extension) in files {
            let given = parsed.given.iter().any(|g| g.key == Key::Picard(option));
            if *option == "OUTPUT" || !given {
                let mut name = OsString::from(path.as_os_str());
                name.push(extension);
                let option = if *option == "OUTPUT" { output } else { option };
                picard.push(format!("{option}={}", Path::new(&name).display()));
            }
        }
    }

    /// A required row with no value: an error that names its flag and its
    /// Picard name.
    fn check_required(&self, flags: &[Flag], picard: &[String]) -> Result<(), VerbError> {
        for flag in flags.iter().filter(|f| f.required) {
            let Some(name) = flag.picard() else { continue };
            if picard.iter().any(|a| a.starts_with(&format!("{name}="))) {
                continue;
            }
            let mut error = VerbError::usage(
                self.name,
                format!("missing --{} (Picard: {name})", flag.long),
            );
            if name == "BAIT_INTERVALS" || name == "TARGET_INTERVALS" {
                error = error.hint("Use --panel FILE when baits and targets are the same file.");
            }
            return Err(error.hint(format!("every option: bammetrics help {}", self.name)));
        }
        Ok(())
    }

    /// Every `INPUT` that names a file must exist (`-` is standard input).
    fn check_inputs(&self, picard: &[String]) -> Result<(), VerbError> {
        for input in picard.iter().filter_map(|a| a.strip_prefix("INPUT=")) {
            if input == "-" || self.input == "PREFIX" || Path::new(input).exists() {
                continue;
            }
            return Err(VerbError::input(self.name, format!("{input} not found")).from(self.input));
        }
        Ok(())
    }

    fn threads(&self, parsed: &Parsed) -> Result<usize, VerbError> {
        let Some(given) = parsed.last(Key::Threads) else {
            return Ok(pass::default_threads());
        };
        given.values[0].parse().map_err(|_| {
            VerbError::usage(self.name, format!("{:?} is not a number", given.values[0]))
                .from(given.spelled.clone())
        })
    }
}

/// The `Options:` lines of `--help`: spellings, Picard names, help.
pub fn option_lines(flags: &[Flag]) -> String {
    let width = flags.iter().map(|f| f.spellings().len()).max().unwrap_or(0);
    let mut text = String::new();
    for flag in flags {
        let picard = if flag.names.is_empty() || !matches!(flag.key, Key::Picard(_)) {
            String::new()
        } else {
            format!("[{}] ", flag.names.join(", "))
        };
        let repeat = if flag.repeated { " Repeatable." } else { "" };
        text.push_str(&format!(
            "  {:<width$}  {picard}{}{repeat}\n",
            flag.spellings(),
            flag.help
        ));
    }
    text
}

/// The input name without its directory and extension: the default prefix.
pub fn stem(input: &str) -> String {
    if input == "-" {
        return "stdin".to_owned();
    }
    let name = Path::new(input)
        .file_name()
        .map_or_else(|| input.to_owned(), |n| n.to_string_lossy().into_owned());
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    for extension in [".bam", ".sam", ".cram", ".vcf", ".bcf", ".bed"] {
        if let Some(stem) = name.strip_suffix(extension) {
            return stem.to_owned();
        }
    }
    name.to_owned()
}

/// Create the directories above `path` that do not exist.
fn create_parent(verb: &str, path: &Path) -> Result<(), VerbError> {
    match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => fs::create_dir_all(parent)
            .map_err(|e| VerbError::input(verb, format!("{}: {e}", parent.display()))),
        _ => Ok(()),
    }
}

/// The main file in `format` to `target` (standard output for `None`),
/// and every further file as Picard text.
fn write_single(
    inv: &Invocation,
    result: &bammetrics_core::metric::ToolOutput,
    format: Format,
    context: &Context<'_>,
    target: Option<&Path>,
) -> bammetrics_core::Result<Vec<PathBuf>> {
    let mut written = Vec::new();
    if let Some(main) = &result.main {
        let text = format::render(format, context, main);
        let output_error = |source| bammetrics_core::Error::Output {
            path: target.map_or_else(|| PathBuf::from("-"), Path::to_path_buf),
            source,
        };
        match target {
            None => io::stdout()
                .lock()
                .write_all(text.as_bytes())
                .map_err(output_error)?,
            Some(path) => {
                fs::write(path, text).map_err(output_error)?;
                written.push(path.to_path_buf());
            }
        }
    }
    for (path, file) in &result.extra {
        registry::write_to(Some(path), |out| {
            file.write(file.headers(&inv.command), out)
        })?;
        written.push(path.clone());
    }
    Ok(written)
}

fn written_paths(inv: &Invocation, result: &bammetrics_core::metric::ToolOutput) -> Vec<PathBuf> {
    let main = result.main.as_ref().and(inv.output.clone());
    main.into_iter()
        .chain(result.extra.iter().map(|(p, _)| p.clone()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `CollectInsertSizeMetrics` -> `insert-size`: the class without
    /// `Collect` and `Metrics`, in kebab case.
    fn rule_name(class: &str) -> String {
        let short = class
            .strip_prefix("Collect")
            .unwrap_or(class)
            .replace("Metrics", "");
        let mut name = String::new();
        for (i, c) in short.chars().enumerate() {
            if c.is_ascii_uppercase() && i > 0 {
                name.push('-');
            }
            name.push(c.to_ascii_lowercase());
        }
        name
    }

    /// docs/05 section 2.1: the first twelve names, then the rule.
    #[test]
    fn verb_names_follow_the_table_then_the_rule() {
        let table = [
            ("insert-size", "CollectInsertSizeMetrics"),
            ("alignment-summary", "CollectAlignmentSummaryMetrics"),
            ("quality-yield", "CollectQualityYieldMetrics"),
            ("hs-metrics", "CollectHsMetrics"),
            ("gc-bias", "CollectGcBiasMetrics"),
            ("quality-score-distribution", "QualityScoreDistribution"),
            ("mean-quality-by-cycle", "MeanQualityByCycle"),
            (
                "base-distribution-by-cycle",
                "CollectBaseDistributionByCycle",
            ),
            ("duplicate-metrics", "CollectDuplicateMetrics"),
            ("library-complexity", "EstimateLibraryComplexity"),
            ("multiple-metrics", "CollectMultipleMetrics"),
            ("accumulate-quality-yield", "AccumulateQualityYieldMetrics"),
        ];
        for (name, class) in table {
            assert_eq!(find(name).map(|v| v.class), Some(class), "{name}");
        }
        for verb in VERBS
            .iter()
            .filter(|v| !table.iter().any(|(n, _)| *n == v.name))
        {
            assert_eq!(verb.name, rule_name(verb.class), "{}", verb.class);
        }
        for tool in registry::TOOLS {
            assert!(for_class(tool.name).is_some(), "{} has no verb", tool.name);
        }
        assert_eq!(
            rule_name("CollectWgsMetricsWithNonZeroCoverage"),
            "wgs-with-non-zero-coverage"
        );
    }

    #[test]
    fn every_verb_flag_has_one_long_spelling() {
        for verb in VERBS {
            let flags = verb.flags();
            for flag in &flags {
                let same = flags.iter().filter(|f| f.long == flag.long).count();
                assert_eq!(same, 1, "{}: --{}", verb.name, flag.long);
            }
        }
    }

    #[test]
    fn stems() {
        assert_eq!(stem("dir/NA12878.bam"), "NA12878");
        assert_eq!(stem("calls.vcf.gz"), "calls");
        assert_eq!(stem("-"), "stdin");
        assert_eq!(stem("out/shard.1"), "shard.1");
    }
}

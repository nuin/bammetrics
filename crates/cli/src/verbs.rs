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

use bammetrics_core::fifo_buffer;
use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::registry::{self, Invocation, Kind, Tool, Ty, normalize};
use bammetrics_core::{
    bam_index_stats, build_bam_index, check_duplicate_marking, compare_metrics, cram, pass,
};

use crate::error::VerbError;
use crate::format::{self, Context, Format};
use crate::grammar::{self, Flag, Given, Key, Parsed};
use crate::report::write_atomic;

/// How a verb names what it writes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Files {
    /// `-o FILE`, or standard output.
    Single,
    /// `-o PREFIX`, default `./<input stem>`.
    Prefix,
    /// `-o BAM` (`-o -` for standard output) is the tool's `OUTPUT`; the
    /// metrics file has a flag of its own (`mark-duplicates --metrics`).
    Bam,
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
    /// The `--help` summary when it differs from [`Tool::summary`], which
    /// describes the rows Python returns.
    pub about: Option<&'static str>,
    /// What the main file is, for the `-o` and `--format` help.
    pub writes: &'static str,
    /// Options whose spellings name one source: two of them on one command
    /// line are a usage error (`--from-bam` and `--dict` both give
    /// `SEQUENCE_DICTIONARY`). One spelling given twice keeps the last.
    pub one_spelling: &'static [&'static str],
    /// Picard options beside `INPUT` that name a file the verb reads: a
    /// missing one is an input error from the spelling given, and `-o`
    /// may not name it (`SEQUENCE_DICTIONARY` of `bed-to-interval-list`).
    pub reads: &'static [&'static str],
    /// Picard options beside `OUTPUT` that name a file the verb writes: it
    /// may not be an input, nor the `-o` file (`OUTPUT_TABLE` of
    /// `compare-metrics`).
    pub also_writes: &'static [&'static str],
    /// Shared flags the verb accepts and ignores, with the reason `--help`
    /// prints.
    pub ignores: &'static [(Key, &'static str)],
    /// `-` reads standard input. False for a tool that reads each input
    /// twice (`jumping-library`): its `--help` says standard input is
    /// refused, and the tool exits 2 on `-`, `/dev/stdin` or a pipe.
    pub reads_stdin: bool,
    /// For a [`Files::Bam`] verb whose `-o` may be left out: the file it
    /// then writes, from the input path (`build-bam-index` puts the index
    /// beside the BAM, as Picard names it).
    pub default_output: Option<fn(&Path) -> PathBuf>,
}

const fn verb(name: &'static str, class: &'static str, files: Files) -> Verb {
    Verb {
        name,
        class,
        files,
        input: "BAM",
        aliases: &[],
        header: None,
        about: None,
        writes: "metrics file",
        one_spelling: &[],
        reads: &[],
        also_writes: &[],
        ignores: &[],
        reads_stdin: true,
        default_output: None,
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
        aliases: &[("EXPECTED_PAIR_ORIENTATION", "EXPECTED_PAIR_ORIENTATIONS")],
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
        aliases: &[("BAITS", "BAIT_INTERVALS"), ("TARGETS", "TARGET_INTERVALS")],
        header: Some(bammetrics_core::hs_metrics::subcommand_header),
        reads: &["BAIT_INTERVALS", "TARGET_INTERVALS"],
        ..verb("hs-metrics", "CollectHsMetrics", Files::Single)
    },
    verb("gc-bias", "CollectGcBiasMetrics", Files::Prefix),
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
    verb("multiple-metrics", "CollectMultipleMetrics", Files::Prefix),
    with_input(
        verb(
            "accumulate-quality-yield",
            "AccumulateQualityYieldMetrics",
            Files::Single,
        ),
        "FILE",
    ),
    verb("wgs", "CollectWgsMetrics", Files::Single),
    verb("raw-wgs", "CollectRawWgsMetrics", Files::Single),
    verb(
        "wgs-with-non-zero-coverage",
        "CollectWgsMetricsWithNonZeroCoverage",
        Files::Single,
    ),
    verb("targeted-pcr", "CollectTargetedPcrMetrics", Files::Single),
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
    Verb {
        input: "BED",
        about: Some(
            "writes an interval list: the @SQ lines of the dictionary (--from-bam BAM, --dict \
             FILE or --sd FILE), then one line per BED interval, 1-based and closed.",
        ),
        writes: "interval list",
        aliases: &[
            ("FROM_BAM", "SEQUENCE_DICTIONARY"),
            ("DICT", "SEQUENCE_DICTIONARY"),
        ],
        one_spelling: &["SEQUENCE_DICTIONARY"],
        reads: &["SEQUENCE_DICTIONARY"],
        ignores: &[
            (
                Key::Threads,
                "no BAM record is read, a BAM gives only its header",
            ),
            (Key::Sample, "an interval list has no sample"),
        ],
        ..verb("bed-to-interval-list", "BedToIntervalList", Files::Single)
    },
    Verb {
        input: "FILE",
        about: Some(
            "compares two metrics files as Picard does: the comparison (both files, equal or NOT \
             equal, one line per difference) to -o or standard output, exit 1 when they differ.",
        ),
        writes: "comparison",
        also_writes: &["OUTPUT_TABLE"],
        ignores: &[
            (Key::Threads, "no BAM is read"),
            (Key::Sample, "a comparison has no sample"),
        ],
        ..verb("compare-metrics", "CompareMetrics", Files::Single)
    },
    Verb {
        input: "FASTA",
        about: Some(
            "writes a sequence dictionary: @HD, then one @SQ line per FASTA sequence with SN, \
             LN, M5 and UR. In picard format -o FILE refuses an existing file, as Picard \
             does; the other formats replace it, as every verb does.",
        ),
        writes: "sequence dictionary",
        reads: &["REFERENCE", "ALT_NAMES"],
        ignores: &[
            (Key::Threads, "no BAM is read"),
            (Key::Sample, "a sequence dictionary has no sample"),
        ],
        ..verb(
            "create-sequence-dictionary",
            "CreateSequenceDictionary",
            Files::Single,
        )
    },
    Verb {
        about: Some(
            "marks duplicate reads: -o gets every record with the duplicate flag (0x400) set \
             or cleared and a @PG line, --metrics-file the DuplicationMetrics rows per library \
             and the set_size histograms. Several inputs are merged; the input is read twice.",
        ),
        writes: "marked BAM",
        aliases: &[("METRICS", "METRICS_FILE")],
        ..verb("mark-duplicates", "MarkDuplicates", Files::Bam)
    },
    verb(
        "sequencing-artifact",
        "CollectSequencingArtifactMetrics",
        Files::Prefix,
    ),
    Verb {
        input: "PREFIX",
        reads: &["PRE_ADAPTER_IN", "BAIT_BIAS_IN"],
        also_writes: &["OXOG_OUT"],
        ignores: &[(Key::Threads, "no BAM is read")],
        ..verb(
            "convert-sequencing-artifact-to-oxo-g",
            "ConvertSequencingArtifactToOxoG",
            Files::Prefix,
        )
    },
    verb("oxo-g", "CollectOxoGMetrics", Files::Single),
    Verb {
        aliases: &[
            ("TARGETS", "TARGET_INTERVALS"),
            ("PANEL", "TARGET_INTERVALS"),
        ],
        reads: &["TARGET_INTERVALS"],
        also_writes: &["GAPS_FILE"],
        ..verb("coverage-gaps", "CollectCoverageGaps", Files::Single)
    },
    Verb {
        aliases: &[("TARGETS", "TARGET_INTERVALS"), ("SEX", "EXPECTED_SEX")],
        reads: &["TARGET_INTERVALS"],
        ..verb("sex-check", "CollectSexCheckMetrics", Files::Single)
    },
    Verb {
        about: Some(
            "prints Picard's index statistics: one line per reference of the BAM header with \
             its length and the Aligned= and Unaligned= counts of the .bai or .csi beside the \
             BAM, then NoCoordinateCount=. A CRAM, or a BAM on standard input, is read whole \
             and its records counted.",
        ),
        writes: "index statistics",
        ..verb("bam-index-stats", "BamIndexStats", Files::Single)
    },
    Verb {
        about: Some(
            "checks that every record of a queryname carries the same duplicate flag, as \
             Picard's CheckDuplicateMarking does: the bad querynames (one line per record that \
             disagrees) to -o or standard output, the summary on standard error (dropped by \
             --verbosity as Picard's log line is), exit 1 when any record disagrees. A BAM \
             that is not queryname sorted is sorted in memory.",
        ),
        writes: "bad querynames",
        ..verb(
            "check-duplicate-marking",
            "CheckDuplicateMarking",
            Files::Single,
        )
    },
    Verb {
        reads_stdin: false,
        ..verb(
            "jumping-library",
            "CollectJumpingLibraryMetrics",
            Files::Single,
        )
    },
    Verb {
        about: Some(
            "writes the .bai index of a coordinate-sorted BAM, byte for byte the file Picard's \
             BuildBamIndex writes for the same BAM: to -o, or beside the BAM (x.bai for x.bam, \
             else the input name plus .bai). A BAM whose header does not say SO:coordinate is \
             refused, as Picard refuses it. The Picard face without O= writes that file name \
             into the current directory, as Picard 3.5.0 does.",
        ),
        writes: "BAM index",
        default_output: Some(build_bam_index::default_output),
        ..verb("build-bam-index", "BuildBamIndex", Files::Bam)
    },
    // No input and no -o: standard input to standard output (`copy_stream`).
    Verb {
        input: "",
        about: Some(
            "copies standard input to standard output through a memory buffer of \
             --buffer-size bytes, as Picard's FifoBuffer does, so a producer and a consumer \
             that work in bursts do not wait for each other. Any bytes pass unchanged. One \
             line on standard error says that bammetrics chain needs no such pipe; -q drops \
             it.",
        ),
        writes: "copied bytes",
        ignores: &[(Key::Threads, "one thread reads and one writes")],
        ..verb("fifo-buffer", "FifoBuffer", Files::Single)
    },
];

/// What `fifo-buffer` says on standard error unless `-q` is given: the
/// pipe it buffers is one `bammetrics chain` does not have.
pub const CHAIN_NOTE: &str = "bammetrics chain runs its steps in one process and removes the \
                              pipe this buffer serves (-q drops this line)";

/// Verb spellings every verb's rows take, beside a verb's own
/// [`Verb::aliases`]: `--reference` is `REFERENCE_SEQUENCE` as in 0.12.
const SHARED_ALIASES: &[(&str, &str)] = &[("REFERENCE", "REFERENCE_SEQUENCE")];

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
                "BGZF decompression threads, or CRAM container decoders (default: available \
                 cores, at most 8; {} here; 0 counts as 1).",
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
        let mut flags = self.flags_with_files();
        if self.copies_stream() {
            // Bytes to standard output: no file, no format, no sample.
            flags.retain(|f| !names_a_file(f.key));
        }
        flags
    }

    /// [`Verb::flags`] with the `-o`, `--format` and `--sample` a stream
    /// verb drops: [`Verb::file_option_error`] parses against them.
    fn flags_with_files(&self) -> Vec<Flag> {
        let tool = self.tool();
        let table = tool.table();
        let mut flags = Vec::new();
        if let Some(input) = table.iter().find(|o| o.picard() == "INPUT") {
            let stdin = if self.reads_stdin {
                "; - reads standard input."
            } else {
                ". Standard input and pipes are refused: the tool reads each input twice."
            };
            let mut flag = Flag::new(
                Key::Input,
                "input",
                &["i"],
                Ty::Path,
                &format!("{} The positional argument{stdin}", input.help),
            )
            .names(input.names)
            .placeholder(self.input);
            flag.repeated = input.repeated;
            flags.push(flag);
        }
        let output_names: Vec<&'static str> = table
            .iter()
            .filter(|o| o.kind == Kind::Output && self.files != Files::Bam)
            .flat_map(|o| o.names.iter().copied())
            .chain(["OUTPUT", "O"])
            .collect();
        let (placeholder, help) = match self.files {
            Files::Bam if self.default_output.is_some() => (
                "FILE",
                format!(
                    "The {} to write; - writes standard output. Without -o it goes beside the \
                     input, named as Picard names it.",
                    self.writes
                ),
            ),
            Files::Bam => (
                "BAM",
                "The BAM to write; - writes standard output. Required.".to_owned(),
            ),
            Files::Single => (
                "FILE",
                format!(
                    "The {}; standard output without one, or with -.",
                    self.writes
                ),
            ),
            Files::Prefix => (
                "PREFIX",
                "Output prefix (default ./ plus the input name without its extension); a \
                 directory gets that name as the file part."
                    .to_owned(),
            ),
        };
        flags.push(
            Flag::new(Key::Output, "output", &["o", "O"], Ty::Path, &help)
                .names(&output_names)
                .placeholder(placeholder),
        );
        let format_help = match self.files {
            Files::Single => format!(
                "picard (the {}), json, tsv, jsonl, human, csv or markdown. A list (picard,json) \
                 writes the first to -o and each further one beside it, as FILE.json.",
                self.writes
            ),
            Files::Prefix => "picard first: this verb writes several Picard files. A list \
                              (picard,json) adds each further format beside every metrics file, \
                              as FILE.json."
                .to_owned(),
            Files::Bam if !self.writes_metrics() => {
                format!("picard only: the {} has no other form.", self.writes)
            }
            Files::Bam => "picard first: this verb writes a BAM and a Picard metrics file. A list \
                           (picard,json) adds each further format beside the metrics file, as \
                           FILE.json."
                .to_owned(),
        };
        flags.push(
            Flag::new(Key::Format, "format", &[], Ty::Text, &format_help).placeholder("FORMAT"),
        );
        if self.files == Files::Single {
            flags.push(
                Flag::new(
                    Key::Sample,
                    "sample",
                    &[],
                    Ty::Text,
                    "Sample name for the json and jsonl forms (default the first @RG SM, else the input \
                     name); never written into a Picard file.",
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
            // A BAM verb's `-o` is its `OUTPUT`; its metrics file row stays a flag.
            let output_row = match self.files {
                Files::Bam => opt.picard() == "OUTPUT",
                _ => opt.kind == Kind::Output,
            };
            if opt.picard() == "INPUT" || output_row || is_shared {
                continue;
            }
            let mut names: Vec<&'static str> = opt.names.to_vec();
            names.extend(
                self.aliases
                    .iter()
                    .chain(SHARED_ALIASES)
                    .filter(|(_, target)| *target == opt.picard())
                    .map(|(alias, _)| *alias),
            );
            let help = match opt.kind {
                // A tool's own ignored row (`build-bam-index`) decodes nothing.
                Kind::Ignored(why)
                    if opt.picard() == "REFERENCE_SEQUENCE"
                        && opt.kind != registry::REFERENCE_NOT_DECODED.kind
                        && !tool.options.iter().any(|own| own.picard() == opt.picard()) =>
                {
                    format!("Decodes a CRAM input; for a BAM accepted and ignored: {why}.")
                }
                Kind::Ignored(why) => format!("Accepted and ignored: {why}."),
                _ if opt.required => format!("{} Required.", opt.help),
                _ => opt.help.to_owned(),
            };
            let mut flag = Flag::new(Key::Picard(opt.picard()), &opt.kebab(), &[], opt.ty, &help)
                .names(&names)
                .placeholder(opt.ty.placeholder());
            flag.repeated = opt.repeated;
            flag.required = opt.required;
            flags.push(flag);
        }
        flags.extend(shared);
        for flag in &mut flags {
            if let Some((_, why)) = self.ignores.iter().find(|(key, _)| *key == flag.key) {
                flag.help = format!("Accepted and ignored: {why}.");
            }
        }
        flags
    }

    fn ignored(&self, key: Key) -> bool {
        self.ignores.iter().any(|(k, _)| *k == key)
    }

    /// `bammetrics <verb> --help`.
    pub fn help(&self) -> String {
        let tool = self.tool();
        let (input, redirects) = if self.copies_stream() {
            (String::new(), " < IN > OUT")
        } else {
            (format!("[{}] ", self.input), "")
        };
        let mut text = format!(
            "bammetrics {}: {}, {}\n\nUsage: bammetrics {} {input}[OPTIONS]{redirects}\n       \
             bammetrics {} {}=... [KEY=VALUE ...]   (Picard syntax)\n\n",
            self.name,
            tool.name,
            self.about.unwrap_or(tool.summary),
            self.name,
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
        let parsed = grammar::parse(self.name, &flags, args)
            .map_err(|e| self.file_option_error(args).unwrap_or(e))?;
        if parsed.flag(Key::Help) {
            print!("{}", self.help());
            return Ok(0);
        }
        if self.copies_stream() {
            return self.copy_stream(&parsed, &flags);
        }
        let tool = self.tool();
        let formats = match parsed.last(Key::Format) {
            None => vec![Format::Picard],
            Some(given) => format::parse_list(self.name, &given.values[0], &given.spelled)?,
        };
        let format = formats[0];
        if self.files == Files::Bam
            && !self.writes_metrics()
            && formats != [Format::Picard]
            && let Some(given) = parsed.last(Key::Format)
        {
            return Err(VerbError::usage(
                self.name,
                format!("--format {} has nothing to write", given.values[0]),
            )
            .from(given.spelled.clone())
            .hint(format!(
                "this verb writes the {} and no metrics file",
                self.writes
            )));
        }
        if self.files != Files::Single
            && format != Format::Picard
            && let Some(given) = parsed.last(Key::Format)
        {
            let hint = match self.files {
                Files::Bam => "this verb writes a BAM to -o and its metrics to --metrics-file",
                _ => "this verb writes several Picard files under -o PREFIX",
            };
            return Err(VerbError::usage(
                self.name,
                format!("--format {} needs a single-file verb", given.values[0]),
            )
            .from(given.spelled.clone())
            .hint(hint)
            .hint(format!(
                "picard,{} writes the Picard files and each as FILE.{} beside it",
                format.name(),
                format.extension()
            )));
        }
        self.check_one_spelling(&parsed, &flags)?;
        let mut picard = self.picard_arguments(&parsed, &flags)?;
        let output = self.output_path(&parsed, &picard)?;
        let output_flag = || parsed.last(Key::Output).map(|g| g.spelled.clone());
        if let Some(path) = &output {
            self.name_outputs(path, &parsed, &mut picard, &formats);
        }
        self.check_required(&flags, &picard)?;
        self.check_inputs(&parsed, &picard)?;
        let core = |e: bammetrics_core::Error| self.core_error(&parsed, &e);
        let mut inv = tool.parse(&picard).map_err(core)?;
        let read = self.input_paths(&inv);
        // A single-file verb writes each further format beside `-o`.
        let siblings: Vec<PathBuf> = match (&output, self.files) {
            (Some(path), Files::Single) => formats[1..]
                .iter()
                .map(|f| format::sibling(path, f.extension()))
                .collect(),
            (None, Files::Single) if formats.len() > 1 => {
                let given = parsed.last(Key::Format).expect("a list was given");
                return Err(VerbError::usage(
                    self.name,
                    format!(
                        "--format {} writes {} files; standard output holds one",
                        given.values[0],
                        formats.len()
                    ),
                )
                .from(given.spelled.clone())
                .hint("give -o FILE: the further formats are written beside it"));
            }
            _ => Vec::new(),
        };
        let mut writes = self.output_paths(&inv, output.as_deref());
        writes.extend(siblings);
        for written in writes {
            refuse_overwriting_inputs(self.name, &read, std::slice::from_ref(&written))
                .map_err(|e| with_from(e, flag_for_path(&parsed, &written).or_else(output_flag)))?;
        }
        self.refuse_one_file_twice(&inv, output.as_deref())?;
        if let Some(path) = &output {
            create_parent(self.name, path).map_err(|e| with_from(e, output_flag()))?;
        }
        inv.threads = self.threads(&parsed)?;
        if let Some(header) = self.header {
            inv.command = header(&inv).map_err(core)?;
        }
        inv.read_dash_from_stdin();
        let quiet = parsed.flag(Key::Quiet);
        let verbose = parsed.flag(Key::Verbose);
        if verbose {
            eprintln!("bammetrics {}: {}", self.name, inv.command);
        }
        let (result, status) = if tool.name == compare_metrics::TOOL.name {
            let picard = formats.contains(&Format::Picard);
            self.compare(&inv, output.as_deref(), picard, quiet)
                .map_err(core)?
        } else if tool.name == check_duplicate_marking::TOOL.name {
            let check = check_duplicate_marking::check(&inv).map_err(core)?;
            if let Some(summary) = check.summary_shown().filter(|_| !quiet) {
                eprintln!("bammetrics {}: {summary}", self.name);
            }
            let status = check.status();
            (check.output, status)
        } else {
            ((tool.run)(&inv).map_err(core)?, 0)
        };
        if !quiet {
            for note in &inv.notes {
                eprintln!("bammetrics {}: {note}", self.name);
            }
        }
        let sample = match parsed.last(Key::Sample) {
            _ if self.ignored(Key::Sample) => None,
            Some(given) => Some(given.values[0].clone()),
            None if formats
                .iter()
                .any(|f| matches!(f, Format::Json | Format::Jsonl)) =>
            {
                self.default_sample(&inv)
            }
            None => None,
        };
        let context = Context {
            tool: tool.name,
            command: &inv.command,
            sample: sample.as_deref(),
        };
        let written = match self.files {
            Files::Single => write_single(&inv, &result, &formats, &context, output.as_deref()),
            Files::Prefix => registry::write_files(&inv, &result)
                .and_then(|()| write_siblings(&inv, &result, &formats, &context))
                .map(|siblings| {
                    let mut paths = written_paths(&inv, &result);
                    paths.extend(siblings);
                    paths
                }),
            Files::Bam => registry::write_files(&inv, &result)
                .and_then(|()| write_siblings(&inv, &result, &formats, &context))
                .map(|siblings| {
                    let mut paths: Vec<PathBuf> = output.iter().cloned().collect();
                    paths.extend(written_paths(&inv, &result));
                    paths.extend(siblings);
                    paths
                }),
        }
        .map_err(core)?;
        if verbose {
            for path in written {
                eprintln!("bammetrics {}: wrote {}", self.name, path.display());
            }
        }
        Ok(status)
    }

    /// `compare-metrics`: the comparison as the main file (standard output
    /// gets it too, where the Picard face prints nothing, with a final
    /// newline the file lacks), the summary line on standard error when it
    /// goes to a file, and Picard's exit code.
    fn compare(
        &self,
        inv: &Invocation,
        output: Option<&Path>,
        picard: bool,
        quiet: bool,
    ) -> bammetrics_core::Result<(bammetrics_core::metric::ToolOutput, u8)> {
        let mut comparison = compare_metrics::compare(inv)?;
        if picard && let Some(main) = comparison.output.main.as_mut() {
            let mut text = comparison.report.clone();
            if output.is_none() && !text.ends_with('\n') {
                text.push('\n');
            }
            main.text = Some(text);
        }
        if output.is_some() && !quiet {
            eprintln!("bammetrics {}: {}", self.name, comparison.summary);
        }
        let status = comparison.status();
        Ok((comparison.output, status))
    }

    /// An option of [`Verb::one_spelling`] given under two of its
    /// spellings: a usage error from the second. A verb alias is its own
    /// spelling; the row's Picard names and their kebab forms are one.
    pub(crate) fn check_one_spelling(
        &self,
        parsed: &Parsed,
        flags: &[Flag],
    ) -> Result<(), VerbError> {
        for &name in self.one_spelling {
            let spelling = |given: &Given| {
                let word = normalize(given.spelled.trim_start_matches('-').trim_end_matches('='));
                self.aliases
                    .iter()
                    .find(|(alias, target)| *target == name && *alias == word)
                    .map_or(name, |(alias, _)| *alias)
            };
            let mut given = parsed.all(Key::Picard(name));
            let Some(first) = given.next() else { continue };
            let Some(second) = given.find(|g| spelling(g) != spelling(first)) else {
                continue;
            };
            let long = flags
                .iter()
                .find(|f| f.key == Key::Picard(name))
                .map_or_else(|| name.to_owned(), |f| f.long.clone());
            let spellings: Vec<String> = self
                .aliases
                .iter()
                .filter(|(_, target)| *target == name)
                .map(|(alias, _)| format!("--{}", alias.to_ascii_lowercase().replace('_', "-")))
                .chain([format!("--{long}")])
                .collect();
            return Err(VerbError::usage(
                self.name,
                format!("{} and {} both set {name}", first.spelled, second.spelled),
            )
            .from(second.spelled.clone())
            .hint(format!("give one of {}", spellings.join(", "))));
        }
        Ok(())
    }

    /// The Picard arguments of every option given, in order: `KEY=VALUE`
    /// with the row's long name. Positional arguments are the inputs
    /// (`INPUT`, or the tool's first required options when it has no
    /// `INPUT`, or its first value option when it has neither).
    pub(crate) fn picard_arguments(
        &self,
        parsed: &Parsed,
        flags: &[Flag],
    ) -> Result<Vec<String>, VerbError> {
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
                let mut error = VerbError::usage(
                    self.name,
                    format!("one {} expected, got {}", self.input, inputs.len()),
                )
                .from(inputs[1].clone());
                if let Some(hint) = (0..parsed.positional.len().min(2))
                    .rev()
                    .find_map(|i| parsed.bool_value_hint(i))
                {
                    error = error.hint(hint);
                }
                return Err(error);
            }
            if inputs.is_empty() {
                let mut error = VerbError::usage(self.name, format!("missing {}", self.input));
                if let Some(hint) = parsed.bool_word_hint() {
                    error = error.hint(hint);
                }
                return Err(
                    error.hint(format!("bammetrics {} {} [OPTIONS]", self.name, self.input))
                );
            }
            picard.extend(inputs.iter().map(|i| format!("INPUT={i}")));
        } else {
            let mut slots: Vec<&str> = flags
                .iter()
                .filter(|f| f.required)
                .filter_map(Flag::picard)
                .collect();
            if slots.is_empty() {
                // No required option (`INPUT_BASE` of the OxoG converter is
                // optional): the first value option takes the argument.
                slots.extend(tool.value_options().next());
            }
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
                // The verb keeps the last value of an option Picard takes once;
                // the Picard face refuses a second one.
                Key::Picard(name) if !tool.option(name).is_some_and(|o| o.repeated) => {
                    let prefix = format!("{name}=");
                    picard.retain(|a| !a.starts_with(&prefix));
                    picard.extend(given.values.last().map(|v| format!("{name}={v}")));
                }
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
            Files::Bam => match given {
                Some(g) if g.values[0] == "-" => Ok(Some(PathBuf::from("/dev/stdout"))),
                Some(g) => Ok(Some(PathBuf::from(&g.values[0]))),
                None => match self.default_output {
                    Some(beside) => self.output_beside_input(picard, beside).map(Some),
                    None => Err(
                        VerbError::usage(self.name, "missing -o BAM (Picard: OUTPUT)")
                            .hint("-o - writes the BAM to standard output"),
                    ),
                },
            },
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
                            .or_else(|| a.strip_prefix("INPUT_BASE="))
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
                    // `--oxog-out` names the one file: no prefix is needed.
                    (None, None)
                        if self.also_writes.iter().any(|name| {
                            picard
                                .iter()
                                .any(|a| a.strip_prefix(name).is_some_and(|v| v.starts_with('=')))
                        }) =>
                    {
                        return Ok(None);
                    }
                    (None, None) => {
                        return Err(VerbError::usage(self.name, "missing -o PREFIX"));
                    }
                };
                Ok(Some(path))
            }
        }
    }

    /// The file a verb with [`Verb::default_output`] writes without `-o`:
    /// `beside` the one input. Standard input has no name to derive it from.
    fn output_beside_input(
        &self,
        picard: &[String],
        beside: fn(&Path) -> PathBuf,
    ) -> Result<PathBuf, VerbError> {
        match picard.iter().find_map(|a| a.strip_prefix("INPUT=")) {
            Some(input) if input != "-" => Ok(beside(Path::new(input))),
            _ => Err(VerbError::usage(
                self.name,
                format!(
                    "missing -o FILE: standard input has no name to put the {} beside",
                    self.writes
                ),
            )
            .hint(format!(
                "-o - writes the {} to standard output",
                self.writes
            ))),
        }
    }

    /// True when the tool has a metrics file option of its own
    /// (`mark-duplicates --metrics-file`); false for a verb whose only
    /// file is the one `-o` names (`build-bam-index`).
    fn writes_metrics(&self) -> bool {
        self.tool().table().iter().any(|o| o.kind == Kind::Output)
    }

    /// The tool's file options for output `path`: `OUTPUT` for a
    /// single-file verb that writes the Picard form, naming the Picard file
    /// (`path` when `picard` is the first format, its `.txt` sibling when
    /// it is a further one, so the sibling's `StringHeader` is the one
    /// `-o` of that file alone writes; a list without `picard` gets none,
    /// the verb writes those forms itself, and so does a tool without an
    /// output option, `BamIndexStats`), the prefix as `OUTPUT` for a tool
    /// whose Picard `OUTPUT` is a prefix, or one name per file of
    /// [`file_extensions`] (an option the user gave keeps its value).
    fn name_outputs(
        &self,
        path: &Path,
        parsed: &Parsed,
        picard: &mut Vec<String>,
        formats: &[Format],
    ) {
        let tool = self.tool();
        let output = tool.output_name();
        if self.files == Files::Bam {
            picard.push(format!("OUTPUT={}", path.display()));
            return;
        }
        if self.files == Files::Single {
            let has_output = tool.table().iter().any(|o| o.kind == Kind::Output);
            let picard_file = match formats.iter().position(|f| *f == Format::Picard) {
                Some(0) => Some(path.to_path_buf()),
                Some(_) => Some(format::sibling(path, Format::Picard.extension())),
                None => None,
            };
            if let Some(file) = picard_file.filter(|_| has_output) {
                picard.push(format!("{output}={}", file.display()));
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
    pub(crate) fn check_required(
        &self,
        flags: &[Flag],
        picard: &[String],
    ) -> Result<(), VerbError> {
        for flag in flags.iter().filter(|f| f.required) {
            let Some(name) = flag.picard() else { continue };
            if picard.iter().any(|a| a.starts_with(&format!("{name}="))) {
                continue;
            }
            let mut error = VerbError::usage(
                self.name,
                format!("missing --{} (Picard: {name})", flag.long),
            );
            if self.name == "hs-metrics" && (name == "BAIT_INTERVALS" || name == "TARGET_INTERVALS")
            {
                error = error.hint("Use --panel FILE when baits and targets are the same file.");
            } else if name == "TARGET_INTERVALS" {
                error = error.hint("--targets FILE and --panel FILE spell it too.");
            }
            if name == "SEQUENCE_DICTIONARY" {
                error = error.hint(
                    "--from-bam BAM copies the @SQ lines of a BAM header; --dict FILE reads a .dict.",
                );
            }
            return Err(error.hint(format!("every option: bammetrics help {}", self.name)));
        }
        Ok(())
    }

    /// Every `INPUT` that names a file must exist (`-` is standard input).
    /// The error says where the name came from: `-i` as spelled, else the
    /// positional argument by its placeholder (`BAM`, as docs/05 example 12).
    pub(crate) fn check_inputs(&self, parsed: &Parsed, picard: &[String]) -> Result<(), VerbError> {
        let stdin_reads = picard.iter().filter(|a| *a == "INPUT=-").count();
        if stdin_reads > 1 {
            return Err(VerbError::usage(
                self.name,
                format!("- names standard input {stdin_reads} times; it can be read once"),
            ));
        }
        for input in picard.iter().filter_map(|a| a.strip_prefix("INPUT=")) {
            if input == "-" || self.input == "PREFIX" || Path::new(input).exists() {
                continue;
            }
            let from = parsed
                .all(Key::Input)
                .find(|g| g.values.iter().any(|v| v == input))
                .map_or_else(|| self.input.to_owned(), |g| g.spelled.clone());
            return Err(VerbError::input(self.name, format!("{input} not found")).from(from));
        }
        for &name in self.reads {
            let prefix = format!("{name}=");
            // One option given twice keeps the last value.
            let last = picard
                .iter()
                .rev()
                .find_map(|a| a.strip_prefix(prefix.as_str()));
            if let Some(path) = last.filter(|p| !Path::new(p).exists()) {
                let from = parsed
                    .all(Key::Picard(name))
                    .filter(|g| g.values.iter().any(|v| v == path))
                    .last()
                    .map_or_else(|| name.to_owned(), |g| g.spelled.clone());
                return Err(VerbError::input(self.name, format!("{path} not found")).from(from));
            }
        }
        Ok(())
    }

    /// A core error with the option it names: `MINIMUM_PCT=abc: ...` comes
    /// from the flag that gave `MINIMUM_PCT`, a write error from `-o`.
    pub(crate) fn core_error(&self, parsed: &Parsed, e: &bammetrics_core::Error) -> VerbError {
        let from = match e {
            bammetrics_core::Error::Usage(message) => {
                message.split_once('=').and_then(|(name, _)| {
                    parsed
                        .given
                        .iter()
                        .rev()
                        .find(|g| matches!(g.key, Key::Picard(n) if n == name))
                        .map(|g| g.spelled.clone())
                })
            }
            bammetrics_core::Error::Output { path, .. } => flag_for_path(parsed, path)
                .or_else(|| parsed.last(Key::Output).map(|g| g.spelled.clone())),
            _ => None,
        };
        with_from(VerbError::core(self.name, e), from)
    }

    /// The files the run reads: every `INPUT`, the reference and its `.fai`,
    /// for a tool without `INPUT` its required options, and for
    /// `BamIndexStats` the `.bai` or `.csi` it reads beside each BAM.
    fn input_paths(&self, inv: &Invocation) -> Vec<PathBuf> {
        let tool = self.tool();
        let mut names = vec!["INPUT", "REFERENCE_SEQUENCE"];
        names.extend(self.reads);
        if !tool.takes_input() {
            names.extend(
                tool.table()
                    .iter()
                    .filter(|o| o.required && o.kind == Kind::Value)
                    .map(|o| o.picard()),
            );
        }
        let indexes = inv
            .all("REFERENCE_SEQUENCE")
            .iter()
            .map(|fasta| format!("{fasta}.fai"));
        let mut paths: Vec<PathBuf> = names
            .iter()
            .flat_map(|name| inv.all(name).iter().cloned())
            .chain(indexes)
            .map(PathBuf::from)
            .collect();
        if tool.name == bam_index_stats::TOOL.name {
            paths.extend(
                inv.all("INPUT")
                    .iter()
                    .map(Path::new)
                    .filter(|bam| !cram::is_cram(bam))
                    .filter_map(bam_index_stats::find_index),
            );
        }
        paths
    }

    /// The files the run writes under names it knows before it runs: `-o`,
    /// the metrics file and each file option of [`file_extensions`].
    fn output_paths(&self, inv: &Invocation, output: Option<&Path>) -> Vec<PathBuf> {
        let mut paths: Vec<PathBuf> = output
            .map(Path::to_path_buf)
            .into_iter()
            .chain(inv.output.clone())
            .collect();
        for (option, _) in file_extensions(self.tool().name).unwrap_or_default() {
            paths.extend(inv.all(option).iter().map(PathBuf::from));
        }
        paths.extend(self.further_outputs(inv));
        paths.extend(bammetrics_core::create_sequence_dictionary::md5_output(inv));
        paths
    }

    /// The files of [`Verb::also_writes`].
    fn further_outputs(&self, inv: &Invocation) -> Vec<PathBuf> {
        self.also_writes
            .iter()
            .flat_map(|name| inv.all(name))
            .map(PathBuf::from)
            .collect()
    }

    /// A file of [`Verb::also_writes`] that is the main file: the second
    /// write would replace the first.
    fn refuse_one_file_twice(
        &self,
        inv: &Invocation,
        output: Option<&Path>,
    ) -> Result<(), VerbError> {
        let Some(main) = output.map(Path::to_path_buf).or_else(|| inv.output.clone()) else {
            return Ok(());
        };
        let main_path = same_file_key(&main);
        let metrics = inv
            .output
            .iter()
            .filter(|_| self.files == Files::Bam)
            .cloned();
        for further in self.further_outputs(inv).into_iter().chain(metrics) {
            if same_file_key(&further) == main_path {
                return Err(VerbError::usage(
                    self.name,
                    format!(
                        "{} is the main output too; one file cannot hold both",
                        further.display()
                    ),
                ));
            }
        }
        Ok(())
    }

    /// `--sample` for json and jsonl when none is given (docs/05 section
    /// 2.2): the first `@RG` `SM` of a BAM, else the input name without its
    /// extension. Standard input has neither.
    fn default_sample(&self, inv: &Invocation) -> Option<String> {
        let input = inv.all("INPUT").first()?;
        if input == "/dev/stdin" {
            return None;
        }
        let from_header = (self.input == "BAM")
            .then(|| pass::read_header(Path::new(input)).ok())
            .flatten()
            .and_then(|header| {
                bammetrics_core::levels::read_groups(&header)
                    .into_iter()
                    .find_map(|rg| rg.sample)
            });
        Some(from_header.unwrap_or_else(|| stem(input)))
    }

    fn threads(&self, parsed: &Parsed) -> Result<usize, VerbError> {
        let given = parsed
            .last(Key::Threads)
            .filter(|_| !self.ignored(Key::Threads));
        let Some(given) = given else {
            return Ok(pass::default_threads());
        };
        given.values[0].parse().map_err(|_| {
            VerbError::usage(self.name, format!("{:?} is not a number", given.values[0]))
                .from(given.spelled.clone())
        })
    }
}

/// The flags of a file a verb writes, which a stream verb has not.
fn names_a_file(key: Key) -> bool {
    matches!(key, Key::Output | Key::Format | Key::Sample)
}

/// `fifo-buffer`: the one verb that moves bytes and writes no file.
impl Verb {
    /// True for the verb of `FifoBuffer`: standard input to standard
    /// output, so it has no input argument, no `-o` and no `--format`.
    pub(crate) fn copies_stream(&self) -> bool {
        self.class == fifo_buffer::TOOL.name
    }

    /// A usage error of a stream verb, with the usage line as its hint.
    fn stream_usage(&self, what: String) -> VerbError {
        VerbError::usage(self.name, what)
            .hint(format!("bammetrics {} [OPTIONS] < IN > OUT", self.name))
            .hint(format!("every option: bammetrics help {}", self.name))
    }

    /// The error for `-o`, `--format` or `--sample` on a stream verb, in
    /// any spelling: `args` parse once those flags are back and one of
    /// them was given. `None` for every other verb and every other error.
    fn file_option_error(&self, args: &[String]) -> Option<VerbError> {
        if !self.copies_stream() {
            return None;
        }
        let parsed = grammar::parse(self.name, &self.flags_with_files(), args).ok()?;
        let given = parsed.given.iter().find(|g| names_a_file(g.key))?;
        Some(self.stream_usage(format!(
            "{} is no option of this verb: the bytes go to standard output unchanged",
            given.spelled
        )))
    }

    /// Copy standard input to standard output through the buffer. On
    /// standard error: [`CHAIN_NOTE`] and the notes of ignored options
    /// (not with `-q`), the log lines of `--debug-frequency`, and with
    /// `-v` the Picard command line and the byte count. A closed standard
    /// output ends the copy with exit 0 and no message.
    fn copy_stream(&self, parsed: &Parsed, flags: &[Flag]) -> Result<u8, VerbError> {
        if let Some(word) = parsed.positional.first() {
            return Err(self.stream_usage(format!(
                "{word:?} is no option: the bytes come from standard input"
            )));
        }
        let core = |e: bammetrics_core::Error| self.core_error(parsed, &e);
        let picard = self.picard_arguments(parsed, flags)?;
        let inv = self.tool().parse(&picard).map_err(core)?;
        let options = fifo_buffer::Options::from_invocation(&inv).map_err(core)?;
        let verbose = parsed.flag(Key::Verbose);
        if verbose {
            eprintln!("bammetrics {}: {}", self.name, inv.command);
        }
        if !parsed.flag(Key::Quiet) {
            for note in &inv.notes {
                eprintln!("bammetrics {}: {note}", self.name);
            }
            eprintln!("bammetrics {}: {CHAIN_NOTE}", self.name);
        }
        let name = self.name;
        let log = move |line: &str| eprintln!("bammetrics {name}: {line}");
        let copied = fifo_buffer::copy_standard_streams(&options, log).map_err(core)?;
        if verbose {
            eprintln!(
                "bammetrics {}: copied {} of {} bytes read",
                self.name, copied.written, copied.read
            );
        }
        Ok(0)
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

/// `error` with `from` set when it has none.
/// The flag, as spelled, whose value is `path`, or else the last one whose
/// value is a directory above it (`--tmp-dir` for a file under it).
fn flag_for_path(parsed: &Parsed, path: &Path) -> Option<String> {
    let named = |test: &dyn Fn(&Path) -> bool| {
        parsed
            .given
            .iter()
            .rev()
            .find(|g| g.values.iter().any(|v| test(Path::new(v))))
            .map(|g| g.spelled.clone())
    };
    named(&|v| v == path).or_else(|| named(&|v| !v.as_os_str().is_empty() && path.starts_with(v)))
}

fn with_from(error: VerbError, from: Option<String>) -> VerbError {
    match from {
        Some(from) if error.from.is_none() => error.from(from),
        _ => error,
    }
}

/// An output that is one of the inputs would overwrite it: a usage error.
/// Only existing files compare equal, after resolving links and `..`.
pub fn refuse_overwriting_inputs(
    verb: &str,
    inputs: &[PathBuf],
    outputs: &[PathBuf],
) -> Result<(), VerbError> {
    let inputs: Vec<PathBuf> = inputs
        .iter()
        .filter_map(|p| p.canonicalize().ok())
        .collect();
    for output in outputs {
        let Ok(resolved) = output.canonicalize() else {
            continue;
        };
        if inputs.contains(&resolved) {
            return Err(VerbError::usage(
                verb,
                format!(
                    "{} is an input; writing it would destroy it",
                    output.display()
                ),
            ));
        }
    }
    Ok(())
}

/// `path` with its directory resolved where it exists, so `a.txt` and
/// `./a.txt` compare equal before either file is written.
fn same_file_key(path: &Path) -> PathBuf {
    let absolute = std::path::absolute(path).unwrap_or_else(|_| path.to_path_buf());
    match (absolute.parent(), absolute.file_name()) {
        (Some(dir), Some(name)) => dir
            .canonicalize()
            .map_or_else(|_| absolute.clone(), |d| d.join(name)),
        _ => absolute,
    }
}

/// Create the directories above `path` that do not exist (docs/05 section
/// 2.3: missing output directories are created).
pub(crate) fn create_parent(verb: &str, path: &Path) -> Result<(), VerbError> {
    match path.parent() {
        Some(parent) if !parent.as_os_str().is_empty() => fs::create_dir_all(parent)
            .map_err(|e| VerbError::input(verb, format!("{}: {e}", parent.display()))),
        _ => Ok(()),
    }
}

/// `text` to `path` through a temporary file, as a core output error.
fn write_text(path: &Path, text: String) -> bammetrics_core::Result<()> {
    write_atomic(path, text).map_err(|e| bammetrics_core::Error::Output {
        path: path.to_path_buf(),
        source: match e {
            crate::error::Error::Io { source, .. } => source,
            other => io::Error::other(other.to_string()),
        },
    })
}

/// The main file in the first of `formats` to `target` (standard output
/// for `None`) and in each further one beside it, and every further file
/// as Picard text.
fn write_single(
    inv: &Invocation,
    result: &bammetrics_core::metric::ToolOutput,
    formats: &[Format],
    context: &Context<'_>,
    target: Option<&Path>,
) -> bammetrics_core::Result<Vec<PathBuf>> {
    let mut written = Vec::new();
    if let Some(main) = &result.main {
        let mut texts = format::render_all(formats, context, main).into_iter();
        let text = texts.next().expect("one text per format");
        match target {
            None => io::stdout()
                .lock()
                .write_all(text.as_bytes())
                .map_err(|source| bammetrics_core::Error::Output {
                    path: PathBuf::from("-"),
                    source,
                })?,
            Some(path) => {
                write_text(path, text)?;
                written.push(path.to_path_buf());
                for (format, text) in formats[1..].iter().zip(texts) {
                    let sibling = format::sibling(path, format.extension());
                    write_text(&sibling, text)?;
                    written.push(sibling);
                }
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

/// Every metrics file `result` holds in each format after the first, as
/// `FILE.<extension>` beside the Picard file a multi-file verb wrote. A
/// file that is not a metrics file (a BED of gaps) gets none.
fn write_siblings(
    inv: &Invocation,
    result: &bammetrics_core::metric::ToolOutput,
    formats: &[Format],
    context: &Context<'_>,
) -> bammetrics_core::Result<Vec<PathBuf>> {
    let mut written = Vec::new();
    if formats.len() < 2 {
        return Ok(written);
    }
    let main = result.main.as_ref().zip(inv.output.as_ref());
    let files = main
        .into_iter()
        .chain(result.extra.iter().map(|(path, file)| (file, path)))
        .filter(|(file, _)| file.text.is_none());
    for (file, path) in files {
        let texts = format::render_all(&formats[1..], context, file);
        for (format, text) in formats[1..].iter().zip(texts) {
            let sibling = format::sibling(path, format.extension());
            write_text(&sibling, text)?;
            written.push(sibling);
        }
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
            ("compare-metrics", "CompareMetrics"),
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

//! The command line of `bammetrics chain`: the chain's options and the
//! source in front of the first `+`, then one step per `+`.

use bammetrics_core::pass;
use bammetrics_core::registry::Ty;
use bammetrics_core::stage;

use super::{STEP, Spec, StepSpec, VERB, WRITE, bare, planfile, recipes, usage};
use crate::error::VerbError;
use crate::grammar::{self, Flag, Key, Parsed};
use crate::verbs::{self, Files, VERBS, Verb, option_lines};

/// The chain's own options, by the key the parser files them under.
pub(crate) const REFERENCE: Key = Key::Own("reference");
pub(crate) const TMP_DIR: Key = Key::Own("tmp-dir");
pub(crate) const ALSO: Key = Key::Own("also");
pub(crate) const PLAN: Key = Key::Own("plan");
pub(crate) const DRY_RUN: Key = Key::Own("dry-run");
pub(crate) const SORT_MEMORY: Key = Key::Own("sort-memory");
pub(crate) const FORCE_SPACE: Key = Key::Own("force-space");
pub(crate) const STRINGENCY: Key = Key::Own("validation-stringency");
pub(crate) const RECIPE: Key = Key::Own("recipe");
pub(crate) const RESUME: Key = Key::Own("resume");
const LOWER: Key = Key::Own("lower");
const FREE_SPACE: Key = Key::Own("free-space");
/// A transform's `--as NAME`.
pub(crate) const AS: Key = Key::Own("as");
/// A transform's `--keep`: its stream as `PREFIX.<stream>.bam`.
pub(crate) const KEEP: Key = Key::Own("keep");
/// The options of `write`.
pub(crate) const INDEX: Key = Key::Own("index");
pub(crate) const MD5: Key = Key::Own("md5");
pub(crate) const COMPRESSION_LEVEL: Key = Key::Own("compression-level");
pub(crate) const OUTPUT_FORMAT: Key = Key::Own("output-format");

const USAGE: &str = "bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]";
/// The last hint of a command-line error about a step that cannot be one.
pub(crate) const STEPS_HINT: &str = "the steps: bammetrics chain --help";

/// What a command line asks for.
pub enum Asked {
    /// `--help`, on the chain or on a step: the text to print.
    Help(String),
    /// A chain to run, or to print the plan of with `--dry-run`.
    Chain { spec: Box<Spec>, dry_run: bool },
}

/// The chain's options: the flags in front of the first `+`.
pub fn global_flags() -> Vec<Flag> {
    vec![
        Flag::new(
            Key::Input,
            "input",
            &["i"],
            Ty::Path,
            "The source: a BAM, or a CRAM with -R. The positional argument; - reads standard \
             input.",
        )
        .names(&["INPUT", "I"])
        .placeholder("BAM"),
        Flag::new(
            Key::Output,
            "output-prefix",
            &["o", "O"],
            Ty::Path,
            "Output prefix of every step that names no file of its own (default ./ plus the \
             source name without its extension); a directory gets that name as the file part. \
             Names follow qc: PREFIX.insert_size_metrics.",
        )
        .names(&["OUTPUT", "O"])
        .placeholder("PREFIX"),
        Flag::new(
            REFERENCE,
            "reference",
            &[],
            Ty::Path,
            "Reference FASTA with a .fai index: decodes a CRAM source, and goes to every step \
             that reads a reference and names none.",
        )
        .names(&["REFERENCE_SEQUENCE", "R"])
        .placeholder("FASTA"),
        Flag::new(
            Key::Threads,
            "threads",
            &["@"],
            Ty::Int,
            &format!(
                "BGZF decompression threads of each read, and compression threads of write \
                 (default: available cores, at most 8; {} here).",
                pass::default_threads()
            ),
        )
        .names(&["THREADS"])
        .placeholder("N"),
        Flag::new(
            TMP_DIR,
            "tmp-dir",
            &[],
            Ty::Path,
            "Where a transform spills and a piped source is copied (default: the system \
             temporary directory).",
        )
        .names(&["TMP_DIR"])
        .placeholder("DIR"),
        Flag::new(
            Key::Format,
            "format",
            &[],
            Ty::Text,
            "picard first: a chain names Picard files. A list (picard,json) adds each further \
             format beside every metrics file, as FILE.json. A step takes its own --format.",
        )
        .placeholder("FORMAT"),
        Flag::new(
            ALSO,
            "also",
            &[],
            Ty::Text,
            "Further formats, one or a list: json, jsonl, tsv, csv, parquet, xlsx, docx or sql \
             beside every metrics file; markdown as PREFIX.report.md.",
        )
        .placeholder("FORMAT"),
        Flag::new(
            Key::Sample,
            "sample",
            &[],
            Ty::Text,
            "Sample name for the json and jsonl forms (default the first @RG SM, else the \
             source name); never written into a Picard file.",
        )
        .placeholder("NAME"),
        Flag::new(
            SORT_MEMORY,
            "sort-memory",
            &[],
            Ty::Text,
            "Bytes a sort-sam step holds before a sorted run spills to --tmp-dir: a number, or \
             one with K, M, G or T (1024 based; 2G, 512MiB). Default 2G.",
        )
        .names(&["SORT_MEMORY"])
        .placeholder("BYTES"),
        Flag::new(
            FORCE_SPACE,
            "force-space",
            &[],
            Ty::Bool,
            "Run when the estimated temporary and output bytes exceed the free space of their \
             file system; without it the chain refuses before the first read.",
        )
        .names(&["FORCE_SPACE"]),
        Flag::new(
            STRINGENCY,
            "validation-stringency",
            &[],
            Ty::Text,
            "How a SAM text source is checked: STRICT (default), LENIENT or SILENT. A BAM or \
             CRAM source is read the same at every value.",
        )
        .names(&["VALIDATION_STRINGENCY"])
        .placeholder("LEVEL"),
        Flag::new(
            RESUME,
            "resume",
            &[],
            Ty::Bool,
            "Skip each step whose record in PREFIX.manifest.json matches (verb, options, input \
             and output fingerprints) and whose upstream steps are skipped too; a later \
             segment reads a kept stream (--keep) in place of the source. Writes the manifest.",
        ),
        Flag::new(
            PLAN,
            "plan",
            &[],
            Ty::Path,
            "Read the source, the options and the steps from a TOML plan file \
             (docs/chaining.md). An option given here replaces the file's.",
        )
        .placeholder("FILE"),
        Flag::new(
            RECIPE,
            "recipe",
            &[],
            Ty::Text,
            "Run a named chain (bammetrics recipes lists them): its steps in place of the \
             command line's. STEP.option=value sets an option of one of its steps.",
        )
        .placeholder("NAME"),
        Flag::new(
            DRY_RUN,
            "dry-run",
            &[],
            Ty::Bool,
            "Print the plan and run nothing: the segments, which transform replays or spools, \
             the files of each step, the estimated temporary disk and memory.",
        ),
        Flag::new(
            Key::Quiet,
            "quiet",
            &["q"],
            Ty::Bool,
            "No notes and no summary on standard error; errors stay.",
        ),
        Flag::new(
            Key::Verbose,
            "verbose",
            &["v"],
            Ty::Bool,
            "Print each step's Picard command line and each file written on standard error.",
        ),
        Flag::new(Key::Help, "help", &["h"], Ty::Bool, "This help."),
    ]
}

/// The flags the parser takes in front of the first `+`: the documented
/// ones, `--no-lower`, which the measurements of docs/chaining.md use, and
/// `--free-space`, which the tests of the free-space refusal use.
pub(crate) fn parser_flags() -> Vec<Flag> {
    let mut flags = global_flags();
    flags.push(
        Flag::new(
            LOWER,
            "lower",
            &[],
            Ty::Bool,
            "Run a chain of taps on the read qc uses.",
        )
        .names(&["LOWER"]),
    );
    flags.push(
        Flag::new(
            FREE_SPACE,
            "free-space",
            &[],
            Ty::Int,
            "The free bytes the planner assumes in every directory.",
        )
        .names(&["FREE_SPACE"]),
    );
    flags
}

/// The flags of the `write` step.
pub fn write_flags() -> Vec<Flag> {
    vec![
        Flag::new(
            Key::Output,
            "output",
            &["o", "O"],
            Ty::Path,
            "The BAM to write, or SAM text for a name ending in .sam: the positional \
             argument; - writes standard output. Without one, PREFIX.<stream>.bam \
             (PREFIX.marked.bam after mark-duplicates).",
        )
        .names(&["OUTPUT", "O"])
        .placeholder("PATH"),
        Flag::new(
            INDEX,
            "index",
            &[],
            Ty::Bool,
            "Also write the .bai (PATH with .bai in place of .bam) of a coordinate-sorted \
             stream.",
        )
        .names(&["CREATE_INDEX"]),
        Flag::new(MD5, "md5", &[], Ty::Bool, "Also write PATH.md5.").names(&["CREATE_MD5_FILE"]),
        Flag::new(
            COMPRESSION_LEVEL,
            "compression-level",
            &[],
            Ty::Int,
            "BGZF level, 0 to 9 (default 5).",
        )
        .names(&["COMPRESSION_LEVEL"])
        .placeholder("N"),
        Flag::new(
            OUTPUT_FORMAT,
            "output-format",
            &[],
            Ty::Text,
            "bam or sam: what is written, whatever the name says. write - --output-format sam \
             prints SAM text on standard output.",
        )
        .placeholder("KIND"),
        Flag::new(Key::Help, "help", &["h"], Ty::Bool, "This help."),
    ]
}

/// True for a verb whose tool is a transform in a chain.
pub(crate) fn is_transform(verb: &Verb) -> bool {
    stage::stage_of(verb.tool()).is_some()
}

/// True for a verb that can be a step: a transform, a split or a one-pass
/// tool.
pub(crate) fn is_step(verb: &Verb) -> bool {
    is_transform(verb) || is_split(verb) || verb.tool().collector.is_some()
}

/// True for a split tool: a step that writes the records of its position
/// over several BAMs in a directory.
pub(crate) fn is_split(verb: &Verb) -> bool {
    bammetrics_core::split_sam::is_split(verb.tool())
}

/// The usage error for a verb that cannot be a step, named `place`
/// (`step 2 (wgs)`), with one hint. The caller adds the second and last
/// hint, the one that lists the steps.
///
/// `build-bam-index` names the `write --index` that does its work in a
/// chain, `view-sam` the `write - --output-format sam` that prints the
/// stream. `fifo-buffer` has its own reason. It copies the bytes of a pipe between
/// two commands, and the steps of a chain hand records on in memory, so no
/// later wave makes it a step.
pub(crate) fn not_a_step(verb: &Verb, place: &str) -> VerbError {
    if verb.copies_stream() {
        return usage(format!(
            "{place}: {} buffers the pipe between two commands and a chain has no pipe: it \
             cannot be a step",
            verb.class
        ))
        .hint(
            "leave the step out, or put the buffer on a pipe outside the chain: bammetrics \
             chain BAM + write - | bammetrics fifo-buffer | ...",
        );
    }
    if verb.class == bammetrics_core::build_bam_index::TOOL.name {
        return usage(format!(
            "{place}: {} indexes a BAM file and a chain hands on records: it cannot be a step",
            verb.class
        ))
        .hint(
            "write --index writes the .bai of the BAM a chain writes: bammetrics chain BAM + \
             ... + write out.bam --index",
        );
    }
    if verb.class == bammetrics_core::view_sam::TOOL.name {
        return usage(format!(
            "{place}: {} prints the records of a file as SAM text and a chain hands records to \
             write: it cannot be a step",
            verb.class
        ))
        .hint(
            "write - --output-format sam prints the stream as SAM text: bammetrics chain BAM + \
             ... + write - --output-format sam",
        );
    }
    if verb.class == bammetrics_core::merge_sam_files::TOOL.name {
        return usage(format!(
            "{place}: {} merges several INPUT files and a chain reads one source: it cannot be \
             a step",
            verb.class
        ))
        .hint(
            "merge the files first, then chain the merged BAM; wave W14 of docs/07 adds steps \
             that read more than the stream",
        );
    }
    if matches!(verb.input, "TEXT" | "VCF" | "LIST" | "FASTA") {
        return usage(format!(
            "{place}: {} reads a {} file and no BAM, and a chain hands on the records of a BAM: \
             it cannot be a step",
            verb.class, verb.input
        ))
        .hint(format!(
            "run it on its own: bammetrics {} --help",
            verb.name
        ));
    }
    usage(format!(
        "{place}: {} has no stage and does not read a BAM in one pass: it cannot be a step",
        verb.class
    ))
    .hint("wave W14 of docs/07 runs such a tool on a spooled copy of the stream")
}

/// The flags the parser takes for a step of `verb`: the verb's own, plus
/// `--as` on a transform. The chain-level ones stay in, so the error for
/// one given on a step says where it belongs.
pub(crate) fn step_parser_flags(verb: &Verb) -> Vec<Flag> {
    let mut flags = verb.flags();
    if is_transform(verb) {
        let default = stage::label_of(verb.tool()).unwrap_or_default();
        flags.push(
            Flag::new(
                AS,
                "as",
                &[],
                Ty::Text,
                &format!(
                    "The name of the stream this step leaves (default {default}): later files \
                     are PREFIX.NAME.<extension>."
                ),
            )
            .placeholder("NAME"),
        );
        flags.push(Flag::new(
            KEEP,
            "keep",
            &[],
            Ty::Bool,
            "Also write the stream this step leaves as PREFIX.<stream>.bam, with its .bai when \
             coordinate sorted, and record it in PREFIX.manifest.json: chain --resume reads it \
             in place of the source.",
        ));
    }
    flags
}

/// The keys of a verb's flags that belong to the chain, not to a step.
pub(crate) const CHAIN_KEYS: [Key; 6] = [
    Key::Input,
    Key::Threads,
    Key::Quiet,
    Key::Verbose,
    Key::Sample,
    // The chain reads the source; a step reads none.
    Key::Picard("VALIDATION_STRINGENCY"),
];

/// The error for a chain option given on the step at `place`.
pub(crate) fn belongs_to_chain(place: &str, spelled: &str) -> VerbError {
    usage(format!(
        "{place}: {spelled} belongs to the chain: give it before the first +"
    ))
    .from(spelled)
}

/// The chain option `name` (a flag as spelled, or a plan file's key) is
/// when no flag of the step has that name: `--also`, `--dry-run`, `--plan`
/// and `--output-prefix` everywhere, and on `write` every chain option.
pub(crate) fn chain_option(name: &str) -> Option<Flag> {
    let flags = global_flags();
    let flag = grammar::spelled(&flags, name).filter(|flag| flag.key != Key::Help)?;
    Some(flag.clone())
}

/// An unknown-option error of a step, when the option is one of the
/// chain's: the error that says where it belongs. The step's own parser
/// would only offer its closest name (`--as` for `--also`).
pub(crate) fn chain_option_on_step(error: &VerbError, place: &str) -> Option<VerbError> {
    let spelled = grammar::unknown_option(error)?;
    let flag = chain_option(spelled)?;
    let mut error = belongs_to_chain(place, spelled);
    if flag.key == Key::Output {
        error = error.hint("a step's own file is -o FILE");
    }
    Some(error)
}

/// The flags of `verb` as a step, for `--help` and for a caller that
/// builds a step's arguments: the verb's own without the chain-level ones
/// and, for a transform, without the options of the BAM `write` writes.
pub fn step_flags(verb: &Verb) -> Vec<Flag> {
    let mut flags: Vec<Flag> = step_parser_flags(verb)
        .into_iter()
        .filter(|flag| !CHAIN_KEYS.contains(&flag.key))
        .filter(|flag| !(verb.files == Files::Bam && flag.key == Key::Output && !is_split(verb)))
        .collect();
    for flag in &mut flags {
        match (flag.key, verb.files) {
            (Key::Output, _) if is_split(verb) => {
                "The directory this step writes its files into; made when missing."
                    .clone_into(&mut flag.help);
            }
            (Key::Output, Files::Single) => {
                flag.help = format!(
                    "The {} of this step; - writes standard output. Without one the chain's \
                     prefix names it.",
                    verb.writes
                );
            }
            (Key::Output, _) => {
                "An output prefix for this step alone, in place of the chain's."
                    .clone_into(&mut flag.help);
            }
            (Key::Format, _) => {
                "picard first; a list (picard,json) adds each further format beside this \
                 step's metrics files. Default: the chain's --format."
                    .clone_into(&mut flag.help);
            }
            _ => {}
        }
    }
    flags
}

/// The step verbs of one kind, for the help text.
fn step_names(wanted: fn(&Verb) -> bool) -> String {
    VERBS
        .iter()
        .filter(|verb| wanted(verb))
        .map(|verb| verb.name)
        .collect::<Vec<_>>()
        .join(", ")
}

/// `bammetrics chain --help`.
pub fn help() -> String {
    let taps = |verb: &Verb| is_step(verb) && !is_transform(verb) && !is_split(verb);
    let mut text = format!(
        "bammetrics chain: several steps on one BAM in one process. Records move from step \
         to step in memory: no pipe, no intermediate file.\n\n\
         Usage: {USAGE}\n       \
         bammetrics chain [OPTIONS] BAM STEP [STEP OPTIONS] [STEP ...] [OPTIONS]\n       \
         bammetrics chain --plan FILE [OPTIONS] [BAM]\n       \
         bammetrics chain BAM --recipe NAME [OPTIONS] [STEP.option=value ...]\n\n\
         Every + starts a step; so does a bare word that names a step verb (bammetrics \
         chain BAM mark-duplicates insert-size -o out.bam). The chain's options go in front \
         of the first step or after the last. After the last bare step, -o, -R, --threads, \
         --format and --tmp-dir are the chain's (write without a path and a split keep their \
         -o); in the middle, a step's own -o, -R and --format are its own. A step sees the \
         records the step before it emits.\n  \
         transform  changes the records later steps see: {}\n  \
         tap        reads the records of its position, changes nothing: {}\n  \
         write      write [PATH] writes the records of its position as a BAM\n  \
         split      writes the records of its position over BAMs in its -o DIR: {}\n  \
         qc         the last step: the verdict and report over the taps of its stream, \
         as bammetrics qc writes them; it adds the taps its --preset needs\n\n\
         A step takes its verb's own options in every spelling (--minimum-pct 0.1, \
         --MINIMUM_PCT 0.1, MINIMUM_PCT=0.1); bammetrics chain BAM + VERB --help lists them. \
         A value that is a lone + is written --option=+.\n\n\
         Files: a tap on the last stream writes PREFIX.<extension>, the names qc writes; a \
         tap in front of a transform writes PREFIX.<stream>.<extension> (PREFIX.in.* on the \
         source); mark-duplicates writes PREFIX.marked.duplication_metrics; write without a \
         path writes PREFIX.<stream>.bam. A step's own -o replaces the name. -o out.bam on \
         a chain with a transform and no write adds write out.bam --index as the last step \
         and names the metrics out.<extension>. --keep on a transform writes \
         PREFIX.<stream>.bam; --resume reads PREFIX.manifest.json and skips the steps it \
         proves done.\n\n\
         Defaults: -R, --threads, --tmp-dir and --sort-memory come from the user config file \
         when not given, -R also from BAMMETRICS_REFERENCE (bammetrics config prints them); \
         reorder-sam's dictionary from the reference's .dict, else its .fai.\n\n\
         Exit 0 when every step finished; 1 when a step failed (a tap fails alone, a \
         transform or write stops its segment; the summary says which); 2 for a usage \
         error, found before the first read.\n\nOptions:\n",
        step_names(is_transform),
        step_names(taps),
        step_names(is_split),
    );
    // `--no-config` is taken out before the parser sees the arguments
    // (`crate::config::strip_no_config`); it is listed with the others.
    let mut listed = global_flags();
    listed.push(Flag::new(
        Key::Own("no-config"),
        "no-config",
        &[],
        Ty::Bool,
        "Leave the user config file unread (bammetrics config prints what it gives).",
    ));
    text.push_str(&option_lines(&listed));
    text.push_str("\nOptions of write:\n");
    text.push_str(&option_lines(&write_flags()));
    text
}

/// The help of step `number`: `write`, or a tool verb as a step. A verb
/// that cannot be a step gets the refusal a run gives it.
pub fn step_help(name: &str, number: usize) -> Result<String, VerbError> {
    if name == super::QC {
        return Ok(super::qc_step::help());
    }
    if name == WRITE {
        return Ok(format!(
            "bammetrics chain ... + write [PATH] [OPTIONS]: the records of its position as a \
             BAM. SAM, CRAM and FASTQ output are later waves.\n\nOptions:\n{}",
            option_lines(&write_flags())
        ));
    }
    let verb =
        find_step(name).ok_or_else(|| unknown_step(name, Some(&format!("step {number}"))))?;
    if !is_step(verb) {
        let place = format!("step {number} ({})", verb.name);
        return Err(not_a_step(verb, &place).hint(STEPS_HINT));
    }
    let tool = verb.tool();
    let kind = if is_transform(verb) {
        "a transform: later steps see the records it emits"
    } else if is_split(verb) {
        "a split: it writes the records of its position over BAMs in its -o DIR, renamed \
         into place when the stream ends"
    } else {
        "a tap: it reads the records of its position and changes nothing"
    };
    Ok(format!(
        "bammetrics chain ... + {} [OPTIONS]: {}, {kind}.\n\n\
         Every option also reads as --PICARD_KEY VALUE and PICARD_KEY=VALUE. A boolean takes \
         --x, --no-x, --x=false or --x false. The source, --output-prefix, --threads, \
         --sample, --also, --plan, --dry-run, -q and -v belong to the chain: give them before \
         the first +.\n\nOptions:\n{}",
        verb.name,
        tool.name,
        option_lines(&step_flags(verb))
    ))
}

/// The verb a step names: by verb name or by Picard class.
pub(crate) fn find_step(name: &str) -> Option<&'static Verb> {
    verbs::find(name).or_else(|| verbs::for_class(name))
}

/// The error for a step name that is no verb, no Picard class and not
/// `write`, under `label` when the caller has one.
pub(crate) fn unknown_step(name: &str, label: Option<&str>) -> VerbError {
    let what = match label {
        Some(label) => format!("{label}: unknown step {name:?}"),
        None => format!("unknown step {name:?}"),
    };
    let candidates = VERBS
        .iter()
        .flat_map(|verb| [verb.name, verb.class])
        .chain([WRITE, super::QC]);
    let mut error = usage(what);
    if let Some(best) = bammetrics_core::registry::did_you_mean(name, candidates) {
        error = error.hint(format!("did you mean {best}?"));
    }
    error.hint("a step is a tool verb (bammetrics tools lists them) or write")
}

/// The one value of the last occurrence of `key`.
fn last_value(parsed: &Parsed, key: Key) -> Option<String> {
    parsed.last(key).map(|given| given.values[0].clone())
}

/// The chain's options as given, over `spec` (a plan file's, or empty).
fn apply_globals(parsed: &Parsed, spec: &mut Spec) -> Result<(), VerbError> {
    let inputs: Vec<&String> = parsed
        .positional
        .iter()
        .chain(parsed.all(Key::Input).flat_map(|given| &given.values))
        .collect();
    if let [_, second, ..] = inputs.as_slice() {
        return Err(usage(format!("one source expected, got {}", inputs.len()))
            .from((*second).clone())
            .hint("every step starts with +")
            .hint(format!("usage: {USAGE}")));
    }
    if let Some(source) = inputs.first() {
        spec.source = Some((*source).clone());
    }
    let texts = [
        (Key::Output, &mut spec.prefix),
        (REFERENCE, &mut spec.reference),
        (TMP_DIR, &mut spec.tmp_dir),
        (Key::Format, &mut spec.formats),
        (ALSO, &mut spec.also),
        (Key::Sample, &mut spec.sample),
        (SORT_MEMORY, &mut spec.sort_memory),
        (STRINGENCY, &mut spec.stringency),
    ];
    for (key, slot) in texts {
        if let Some(value) = last_value(parsed, key) {
            *slot = Some(value);
        }
    }
    if let Some(given) = parsed.last(Key::Threads) {
        let threads = given.values[0].parse().map_err(|_| {
            usage(format!("{:?} is not a number", given.values[0])).from(given.spelled.clone())
        })?;
        let threads = pass::check_threads(threads)
            .map_err(|e| usage(format!("{threads} is {e}")).from(given.spelled.clone()))?;
        spec.threads = Some(threads);
    }
    if let Some(given) = parsed.last(FREE_SPACE) {
        let bytes = given.values[0].parse().map_err(|_| {
            usage(format!("{:?} is not a number", given.values[0])).from(given.spelled.clone())
        })?;
        spec.free_space = Some(bytes);
    }
    let flags = [
        (Key::Quiet, &mut spec.quiet, true),
        (Key::Verbose, &mut spec.verbose, true),
        (LOWER, &mut spec.stage_executor, false),
        (FORCE_SPACE, &mut spec.force_space, true),
        (RESUME, &mut spec.resume, true),
    ];
    for (key, slot, when) in flags {
        if parsed.last(key).is_some() {
            *slot = parsed.flag(key) == when;
        }
    }
    Ok(())
}

/// The steps after the first `+`: each a verb and its own arguments.
fn steps(segments: &[&[String]]) -> Result<Vec<StepSpec>, VerbError> {
    segments
        .iter()
        .enumerate()
        .map(|(index, segment)| {
            let (verb, args) = segment.split_first().ok_or_else(|| {
                usage(format!(
                    "step {} is empty: a verb follows each +",
                    index + 1
                ))
                .hint(format!("usage: {USAGE}"))
            })?;
            Ok(StepSpec {
                verb: verb.clone(),
                args: args.to_vec(),
                origin: None,
            })
        })
        .collect()
}

/// The report options of `report` and `qc`. A chain writes its verdict
/// and report in its last step, `qc` (docs/chaining.md section 29).
const REPORT_OPTIONS: [&str; 5] = ["report", "pdf", "pdf-only", "pdf-font", "page-size"];

/// An unknown chain option that is one of [`REPORT_OPTIONS`]: the refusal
/// that names the qc step; `--keep`: the refusal that names the transform
/// step. Any other error as it is.
fn no_report_step(error: VerbError) -> VerbError {
    let Some(spelled) = grammar::unknown_option(&error) else {
        return error;
    };
    let name = spelled.trim_start_matches('-').split('=').next();
    if name.is_some_and(|name| name.eq_ignore_ascii_case("keep")) {
        return usage("--keep belongs to a transform step: give it after the verb")
            .from(spelled)
            .hint("bammetrics chain BAM -o PREFIX + mark-duplicates --keep + insert-size");
    }
    if !name.is_some_and(|name| REPORT_OPTIONS.contains(&name)) {
        return error;
    }
    usage("a chain writes its verdict and report in its last step, qc: give the option there")
        .from(spelled)
        .hint("bammetrics chain BAM + STEP ... + qc --report html,pdf")
        .hint("bammetrics report DIR --pdf writes them for the metrics files a chain wrote to DIR")
}

/// The recipe `--recipe NAME` names, with the `STEP.option=value` words
/// applied: its name, whether `-o` must name its BAM, and its chain. A
/// recipe gives the steps, so it takes no `--plan` and no step on the
/// command line but a last `qc`, which is added to its steps; a
/// `STEP.option=value` word without a recipe is refused.
fn recipe_spec(
    parsed: &Parsed,
    steps: &[StepSpec],
    overrides: &[(String, recipes::Override)],
) -> Result<Option<(String, bool, Spec)>, VerbError> {
    let Some(given) = parsed.last(RECIPE) else {
        if let Some((word, _)) = overrides.first() {
            return Err(usage(format!(
                "{word} sets an option of a recipe's step, and no recipe is named"
            ))
            .from(word.clone())
            .hint("bammetrics chain BAM --recipe NAME STEP.option=value; on a step: STEP --option VALUE"));
        }
        return Ok(None);
    };
    let named: Vec<&str> = parsed
        .all(RECIPE)
        .map(|given| given.values[0].as_str())
        .collect();
    if let [first, .., last] = named.as_slice() {
        return Err(
            usage(format!("--recipe is given twice: {first} and {last}"))
                .from(given.spelled.clone())
                .hint("a chain runs one recipe"),
        );
    }
    if let Some(plan) = parsed.last(PLAN) {
        return Err(
            usage(format!("{} and --recipe both give the steps", plan.spelled))
                .from(plan.spelled.clone()),
        );
    }
    // One step may follow a recipe: `qc`, the last step of any chain, which
    // reads the metrics files the recipe's taps write.
    let added = match steps {
        [] => None,
        [qc] if qc.verb == super::QC => Some(qc),
        _ => {
            return Err(usage(format!(
                "--recipe {} gives the steps: a step on the command line cannot be added",
                given.values[0]
            ))
            .from(given.spelled.clone())
            .hint(
                "STEP.option=value sets an option of a recipe's step: bammetrics recipes lists \
                 them",
            )
            .hint("a last qc step may follow a recipe: bammetrics NAME BAM -o OUT qc --preset P"));
        }
    };
    let (recipe, mut spec) = recipes::load(&given.values[0], overrides)?;
    if let Some(qc) = added {
        if spec.steps.last().is_some_and(|last| last.verb == super::QC) {
            return Err(usage(format!(
                "recipe {} ends in a qc step: a second one cannot follow it",
                recipe.name
            ))
            .from(super::QC)
            .hint(format!(
                "qc.option=value sets an option of the recipe's qc step: bammetrics {} --help",
                recipe.name
            )));
        }
        spec.steps.push(qc.clone());
    }
    Ok(Some((recipe.name.clone(), recipe.output_bam, spec)))
}

/// Parse the arguments after `chain`.
pub fn parse(args: &[String]) -> Result<Asked, VerbError> {
    let (args, no_config) = crate::config::strip_no_config(args);
    let (args, overrides) = recipes::take_overrides(&args)?;
    let args = bare::with_plus(&args)?;
    let mut segments = args.split(|arg| arg == STEP);
    let head = segments.next().unwrap_or_default();
    let rest: Vec<&[String]> = segments.collect();
    let parsed = grammar::parse(VERB, &parser_flags(), head).map_err(no_report_step)?;
    if parsed.flag(Key::Help) {
        return Ok(Asked::Help(help()));
    }
    let steps = steps(&rest)?;
    let wants_help = |step: &StepSpec| step.args.iter().any(|a| a == "--help" || a == "-h");
    if let Some(index) = steps.iter().position(wants_help) {
        return step_help(&steps[index].verb, index + 1).map(Asked::Help);
    }
    let recipe = recipe_spec(&parsed, &steps, &overrides)?;
    let mut spec = match parsed.last(PLAN) {
        None => Spec::default(),
        Some(given) if !steps.is_empty() => {
            return Err(usage(format!(
                "{} gives the steps: a + step on the command line cannot be added",
                given.spelled
            ))
            .from(given.spelled.clone())
            .hint("options and the source may follow --plan FILE; they replace the file's"));
        }
        Some(given) => planfile::read(std::path::Path::new(&given.values[0]))?,
    };
    if parsed.last(PLAN).is_none() {
        spec.steps = steps;
    }
    if let Some((name, output_bam, recipe)) = recipe {
        spec = recipe;
        spec.recipe_bam = output_bam.then_some(name);
    }
    apply_globals(&parsed, &mut spec)?;
    spec.no_config = no_config;
    Ok(Asked::Chain {
        spec: Box::new(spec),
        dry_run: parsed.flag(DRY_RUN),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|a| (*a).to_owned()).collect()
    }

    fn spec(args: &[&str]) -> (Spec, bool) {
        match parse(&strings(args)).expect("parses") {
            Asked::Chain { spec, dry_run } => (*spec, dry_run),
            Asked::Help(_) => panic!("help"),
        }
    }

    #[test]
    fn the_source_and_the_options_come_before_the_first_plus() {
        let (spec, dry_run) = spec(&[
            "-o",
            "out/s",
            "in.bam",
            "-R",
            "ref.fa",
            "--threads=3",
            "TMP_DIR=/scratch",
            "--format",
            "picard,json",
            "--also",
            "csv",
            "--dry-run",
            "+",
            "mark-duplicates",
            "--remove-duplicates",
            "+",
            "insert-size",
            "M=0.1",
            "+",
            "write",
            "-",
        ]);
        assert!(dry_run);
        assert_eq!(spec.source.as_deref(), Some("in.bam"));
        assert_eq!(spec.prefix.as_deref(), Some("out/s"));
        assert_eq!(spec.reference.as_deref(), Some("ref.fa"));
        assert_eq!(spec.threads, Some(3));
        assert_eq!(spec.tmp_dir.as_deref(), Some("/scratch"));
        assert_eq!(spec.formats.as_deref(), Some("picard,json"));
        assert_eq!(spec.also.as_deref(), Some("csv"));
        assert_eq!(
            spec.steps,
            [
                StepSpec::new("mark-duplicates", &["--remove-duplicates"]),
                StepSpec::new("insert-size", &["M=0.1"]),
                StepSpec::new("write", &["-"]),
            ]
        );
    }

    #[test]
    fn picard_spellings_name_the_source_and_the_executor() {
        let (spec, _) = spec(&["I=-", "R=ref.fa", "LOWER=false", "+", "quality-yield"]);
        assert_eq!(spec.source.as_deref(), Some("-"));
        assert_eq!(spec.reference.as_deref(), Some("ref.fa"));
        assert!(spec.stage_executor);
        let (spec, _) = super::tests::spec(&["in.bam", "--no-lower", "-q", "+", "write", "x.bam"]);
        assert!(spec.stage_executor && spec.quiet);
    }

    #[test]
    fn usage_errors_name_what_is_wrong() {
        let error = |args: &[&str]| match parse(&strings(args)) {
            Err(e) => e.to_string(),
            Ok(_) => panic!("{args:?} parses"),
        };
        assert!(error(&["a.bam", "b.bam", "+", "write"]).starts_with(
            "bammetrics chain: one source expected, got 2 (from b.bam)\n  every step starts \
             with +"
        ));
        assert!(
            error(&["a.bam", "+", "+", "write"])
                .starts_with("bammetrics chain: step 1 is empty: a verb follows each +")
        );
        assert!(
            error(&["a.bam", "--thread", "2", "+", "write"]).starts_with(
                "bammetrics chain: unknown option --thread\n  did you mean --threads?"
            )
        );
        assert!(
            error(&["--plan", "p.toml", "+", "write"])
                .starts_with("bammetrics chain: --plan gives the steps")
        );
        assert!(error(&["a.bam", "--threads", "x", "+", "write"]).contains("is not a number"));
        for threads in ["257", "100000000"] {
            assert_eq!(
                error(&["a.bam", "-@", threads, "+", "write"]),
                format!(
                    "bammetrics chain: {threads} is more than the maximum of 256 threads (from -@)"
                )
            );
        }
        // The help of a step that is no verb names the step by its number,
        // as a run names it.
        assert!(
            error(&["a.bam", "+", "write", "+", "inser-size", "--help"]).starts_with(
                "bammetrics chain: step 2: unknown step \"inser-size\"\n  did you mean insert-size?"
            )
        );
    }

    #[test]
    fn the_most_threads_are_taken() {
        let (spec, _) = spec(&["a.bam", "--threads", "256", "+", "write"]);
        assert_eq!(spec.threads, Some(pass::MAX_THREADS));
    }

    /// A chain option no step has, given after the first `+`, gets the
    /// error that says where it belongs in place of the step's closest
    /// name (`--as` for `--also`).
    #[test]
    fn a_chain_option_on_a_step_says_where_it_belongs() {
        let verb = find_step("insert-size").expect("a verb");
        let on_step = |flags: &[Flag], args: &[&str]| {
            let error = grammar::parse("insert-size", flags, &strings(args)).expect_err("unknown");
            chain_option_on_step(&error, "step 2 (insert-size)").map(|e| e.to_string())
        };
        let step = step_parser_flags(verb);
        for spelled in ["--also", "--dry-run", "--no-dry-run", "--plan", "ALSO="] {
            let arg = if spelled.ends_with('=') {
                format!("{spelled}json")
            } else {
                spelled.to_owned()
            };
            assert_eq!(
                on_step(&step, &[&arg, "x"]).as_deref(),
                Some(
                    format!(
                        "bammetrics chain: step 2 (insert-size): {spelled} belongs to the chain: \
                         give it before the first + (from {spelled})"
                    )
                    .as_str()
                ),
                "{spelled}"
            );
        }
        let prefix = on_step(&step, &["--output-prefix", "x"]).expect("the chain's");
        assert!(
            prefix.ends_with("(from --output-prefix)\n  a step's own file is -o FILE"),
            "{prefix}"
        );
        // `write` has none of the chain's options among its flags.
        for spelled in ["--threads", "-q", "--sample", "-@"] {
            let error = on_step(&write_flags(), &[spelled, "2"]).expect(spelled);
            assert!(error.contains(&format!("{spelled} belongs to the chain")));
        }
        // A name that is no option of the chain stays an unknown option.
        assert_eq!(on_step(&step, &["--als", "json"]), None);
        assert_eq!(on_step(&step, &["--lower"]), None);
    }

    #[test]
    fn help_is_asked_on_the_chain_or_on_a_step() {
        let text = |args: &[&str]| match parse(&strings(args)).expect("parses") {
            Asked::Help(text) => text,
            Asked::Chain { .. } => panic!("no help"),
        };
        let chain = text(&["--help"]);
        assert!(chain.starts_with("bammetrics chain: several steps"));
        assert!(chain.contains("transform  changes the records later steps see: mark-duplicates"));
        assert!(chain.contains("--dry-run"));
        assert!(!chain.contains("--lower"));
        let step = text(&["in.bam", "+", "insert-size", "--help"]);
        assert!(step.contains("+ insert-size [OPTIONS]: CollectInsertSizeMetrics, a tap"));
        assert!(step.contains("--minimum-pct"));
        assert!(!step.contains("-@, --threads"));
        let marked = text(&["in.bam", "+", "mark-duplicates", "-h"]);
        assert!(marked.contains("a transform") && marked.contains("--as NAME"));
        // No option line for the BAM; `--output-prefix` is named as the
        // chain's.
        assert!(!marked.contains("--output "));
        assert!(marked.contains("--output-prefix, --threads"));
        assert!(text(&["in.bam", "+", "write", "--help"]).contains("--compression-level"));
    }

    #[test]
    fn every_step_verb_has_one_long_spelling_per_flag() {
        for verb in VERBS.iter().filter(|verb| is_step(verb)) {
            let flags = step_parser_flags(verb);
            for flag in &flags {
                let same = flags.iter().filter(|f| f.long == flag.long).count();
                assert_eq!(same, 1, "{}: --{}", verb.name, flag.long);
            }
        }
    }
}

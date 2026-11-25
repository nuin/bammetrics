//! The steps of a chain whose source is a VCF (wave W16, docs/chaining.md
//! section 33): the stream holds VCF records, so its steps are `write`
//! (a VCF, plain or block-compressed by its name) and the taps of
//! [`stage::variant`] (`variant-calling`, `genotype-concordance`). A step
//! of the other kind is refused before the first read, naming the
//! stream's kind, and so is a variant tap on a stream of alignments.

use std::path::PathBuf;

use bammetrics_core::split_vcfs;
use bammetrics_core::stage::{self, Kind};
use bammetrics_core::vcf_writer::{VcfFormat, VcfOptions};

use super::{
    COMPRESSION_LEVEL, Context, Draft, INDEX, MD5, OUTPUT_FORMAT, Planned, PlannedFile, Position,
    Role, check_types, prefix_under, step_formats,
};
use crate::error::VerbError;
use crate::format;
use crate::grammar::Key;
use crate::qc::with_suffix;
use crate::verbs::Verb;

/// The kind of a VCF a `write` step writes.
pub const KIND: &str = "vcf";

/// True for a verb whose tool is a tap of a variant stream.
pub(crate) fn is_tap(verb: &Verb) -> bool {
    stage::variant::tap_of(verb.tool()).is_some()
}

/// The role of a verb on a variant stream: a tap, a transform
/// (`make-sites-only-vcf`) or a split (`split-vcfs`); `None` for a verb
/// that is no step of one.
pub(super) fn role_of(verb: &Verb) -> Option<Role> {
    let tool = verb.tool();
    if is_tap(verb) {
        Some(Role::Tap)
    } else if stage::variant::transform_label(tool).is_some() {
        Some(Role::Transform)
    } else if stage::variant::is_split(tool) {
        Some(Role::Split)
    } else {
        None
    }
}

/// The refusal of `--keep` on a step of a variant stream, named `place`.
pub(super) fn keep_refusal(place: &str) -> VerbError {
    super::usage(format!(
        "{place}: --keep writes the stream a transform leaves as a BAM, and this stream holds \
         VCF records"
    ))
    .from("--keep")
    .hint("write PATH.vcf.gz after the step writes its records as a VCF")
}

/// The Picard arguments of a variant step: the stream as `INPUT`, then
/// the step's own, parsed and checked as the tool parses them.
fn step_invocation(
    cx: &Context,
    draft: &mut Draft,
    verb: &'static Verb,
    named: &[String],
) -> Result<(Vec<String>, bammetrics_core::registry::Invocation), VerbError> {
    verb.check_one_spelling(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    check_types(draft)?;
    draft.parsed.positional.insert(0, cx.source.to_owned());
    let picard = verb
        .picard_arguments(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    let (input, own) = picard.split_first().expect("the source is the first input");
    let mut argv = vec![input.clone()];
    argv.extend(named.iter().cloned());
    argv.extend(own.iter().cloned());
    verb.check_inputs(&draft.parsed, &argv)
        .map_err(|e| draft.error(e))?;
    let mut inv = verb.tool().parse(&argv).map_err(|e| draft.core(verb, &e))?;
    inv.threads = cx.threads;
    Ok((argv, inv))
}

/// A transform of a variant stream (`make-sites-only-vcf`): it writes no
/// file; the stream it leaves goes on to the steps after it.
pub(super) fn transform_step(
    cx: &Context,
    mut draft: Draft,
    verb: &'static Verb,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let tool = verb.tool();
    if cx.stringency != bammetrics_core::sam_text::Stringency::Strict {
        // Picard's MakeSitesOnlyVcf writes a key the header lacks under a
        // stringency other than STRICT; the write step and split-vcfs check
        // every key as SplitVcfs does.
        return Err(draft
            .usage(format!(
                "{} writes keys the header lacks under --validation-stringency other than \
                 STRICT, and the steps after it check every key",
                tool.name
            ))
            .hint(format!(
                "run it on its own: bammetrics {} IN.vcf --validation-stringency LENIENT -o OUT.vcf.gz",
                verb.name
            )));
    }
    let (argv, inv) = step_invocation(cx, &mut draft, verb, &[])?;
    stage::variant::VariantStep::tool_transform(tool, &inv).map_err(|e| draft.core(verb, &e))?;
    let planned = Planned {
        number: draft.number,
        verb: verb.name,
        tool: Some(tool),
        role: Role::Transform,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: position.leaves,
        segment: position.segment,
        files: Vec::new(),
        late_prefix: None,
        stdout: false,
        command: Some(inv.command),
        argv,
        output: None,
        notes: inv.notes,
        formats: Vec::new(),
        write: None,
        vcf: None,
        inputs: inv.inputs,
    };
    Ok((planned, false))
}

/// A split of a variant stream (`split-vcfs`): its two VCFs are the step's
/// `--snp-output` and `--indel-output`, else `PREFIX.snp.vcf.gz` and
/// `PREFIX.indel.vcf.gz` under the step's `-o` or the chain's prefix.
pub(super) fn split_step(
    cx: &Context,
    mut draft: Draft,
    verb: &'static Verb,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let tool = verb.tool();
    let own_output = draft.parsed.last(Key::Output).map(|g| g.values[0].clone());
    let (prefix, uses_prefix) = match own_output.as_deref() {
        Some("-") => {
            return Err(draft
                .usage("standard output takes one file; this step writes several")
                .hint("give the step -o PREFIX"));
        }
        Some(path) => (prefix_under(path, cx.source).0, false),
        None => (cx.base(&position.reads), true),
    };
    let given = |key: &'static str| draft.parsed.given.iter().any(|g| g.key == Key::Picard(key));
    let mut named = Vec::new();
    let mut takes_prefix = false;
    for (option, extension) in
        bammetrics_core::multiple_metrics::file_extensions(tool.name).unwrap_or_default()
    {
        if !given(option) {
            named.push(format!(
                "{option}={}",
                with_suffix(&prefix, extension).display()
            ));
            takes_prefix = true;
        }
    }
    let (mut argv, inv) = step_invocation(cx, &mut draft, verb, &named)?;
    argv.retain(|arg| !arg.starts_with("OUTPUT=") && !arg.starts_with("O="));
    stage::variant::VariantStep::tool_split(tool, &inv).map_err(|e| draft.core(verb, &e))?;
    let options = split_vcfs::Options::from_invocation(&inv).map_err(|e| draft.core(verb, &e))?;
    let mut notes = inv.notes;
    let mut files = Vec::new();
    for (option, path) in [
        ("SNP_OUTPUT", &options.snp_output),
        ("INDEL_OUTPUT", &options.indel_output),
    ] {
        let (vcf, note) = options
            .writer
            .for_output(option, path)
            .map_err(|e| draft.core(verb, &e))?;
        notes.extend(note);
        for file in vcf.files(path) {
            let kind = if file == *path {
                KIND
            } else if file == bammetrics_core::vcf_writer::index_path(path) {
                "index"
            } else {
                "md5"
            };
            files.push(PlannedFile { path: file, kind });
        }
    }
    let planned = Planned {
        number: draft.number,
        verb: verb.name,
        tool: Some(tool),
        role: Role::Split,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: None,
        segment: position.segment,
        files,
        late_prefix: None,
        stdout: false,
        command: Some(inv.command),
        argv,
        output: None,
        notes,
        formats: Vec::new(),
        write: None,
        vcf: None,
        inputs: inv.inputs,
    };
    Ok((planned, uses_prefix && takes_prefix))
}

/// The words for what a stream of `kind` holds.
fn holds(kind: Kind) -> &'static str {
    match kind {
        Kind::Alignments => "alignments (the source is a BAM, CRAM or SAM file)",
        Kind::Variants => "VCF records (the source is a VCF)",
    }
}

/// The refusal of `verb` as a step on a stream of `kind`: a step reads
/// one kind of record, and a chain hands on one kind.
pub(super) fn refusal(verb: &Verb, place: &str, kind: Kind) -> VerbError {
    if verb.class == bammetrics_core::make_vcf_sample_name_map::TOOL.name {
        return super::not_a_step(verb, place);
    }
    if kind == Kind::Variants && !super::is_step(verb) {
        return super::usage(format!(
            "{place}: {} reads a {} file of its own: it cannot be a step; a VCF stream takes \
             write, the taps variant-calling and genotype-concordance, the transform \
             make-sites-only-vcf and the split split-vcfs",
            verb.class, verb.input
        ))
        .hint(format!(
            "run it on its own: bammetrics {} --help",
            verb.name
        ));
    }
    let reads = if role_of(verb).is_some() {
        "VCF records"
    } else {
        "alignments"
    };
    let error = super::usage(format!(
        "{place}: {} reads {reads}, and this chain's stream holds {}: it cannot be a step of \
         this chain",
        verb.class,
        holds(kind)
    ));
    match kind {
        Kind::Variants => error.hint(
            "on a VCF stream: write [PATH] writes a VCF; make-sites-only-vcf is a transform, \
             split-vcfs a split, variant-calling and \
             genotype-concordance are taps",
        ),
        Kind::Alignments => error.hint(format!(
            "give the chain a VCF source: bammetrics chain calls.vcf.gz {} ...",
            verb.name
        )),
    }
}

/// The refusal of the qc step on a variant stream.
pub(super) fn qc_refusal(draft: &Draft) -> VerbError {
    draft
        .usage(format!(
            "qc reads the metrics of the taps of a BAM, and this chain's stream holds {}",
            holds(Kind::Variants)
        ))
        .hint("bammetrics check DIR reads the metrics files a chain wrote")
}

/// The `write` step of a variant stream: a VCF at its path, or at
/// `PREFIX.<stream>.vcf.gz`, with its `.tbi` and `.md5` when asked for.
pub(super) fn write_step(
    cx: &Context,
    draft: &Draft,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let parsed = &draft.parsed;
    let paths: Vec<&String> = parsed
        .positional
        .iter()
        .chain(parsed.all(Key::Output).flat_map(|given| &given.values))
        .collect();
    if let [_, second, ..] = paths.as_slice() {
        return Err(draft
            .usage(format!("one PATH expected, got {}", paths.len()))
            .from((*second).clone()));
    }
    if let Some(given) = parsed.last(OUTPUT_FORMAT) {
        return Err(draft
            .usage(format!(
                "{} picks a BAM, SAM or CRAM, and this stream holds {}: the name picks a plain \
                 (.vcf) or block-compressed (.vcf.gz) VCF",
                given.spelled,
                holds(Kind::Variants)
            ))
            .from(given.spelled.clone()));
    }
    let (path, uses_prefix) = match paths.first().map(|path| path.as_str()) {
        Some("-") => {
            return Err(draft
                .usage("a VCF is written to a file: standard output takes no VCF in a chain")
                .from("-"));
        }
        Some(path) => (PathBuf::from(path), false),
        None => (
            with_suffix(cx.prefix, &format!(".{}.vcf.gz", position.reads)),
            true,
        ),
    };
    let mut options = VcfOptions {
        threads: cx.threads,
        index: parsed.flag(INDEX),
        md5: parsed.flag(MD5),
        ..VcfOptions::default()
    };
    if let Some(given) = parsed.last(COMPRESSION_LEVEL) {
        options.compression_level = given.values[0]
            .parse()
            .ok()
            .filter(|level| *level <= 9)
            .ok_or_else(|| {
                draft
                    .usage(format!("expected 0 to 9, got {:?}", given.values[0]))
                    .from(given.spelled.clone())
            })?;
    }
    options.check(&path).map_err(|e| {
        let error = draft.error(VerbError::core(super::VERB, &e));
        if VcfFormat::of(&path).is_ok_and(|f| f == VcfFormat::Plain) && options.index {
            error.hint(format!(
                "write {}.gz --index writes a block-compressed VCF and its .tbi",
                path.display()
            ))
        } else {
            error
        }
    })?;
    let kinds = [KIND, "index", "md5"];
    let mut files = vec![PlannedFile {
        path: path.clone(),
        kind: kinds[0],
    }];
    let all = options.files(&path);
    for extra in &all[1..] {
        let kind = if *extra == bammetrics_core::vcf_writer::index_path(&path) {
            kinds[1]
        } else {
            kinds[2]
        };
        files.push(PlannedFile {
            path: extra.clone(),
            kind,
        });
    }
    let planned = Planned {
        number: draft.number,
        verb: super::WRITE,
        tool: None,
        role: Role::Write,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: None,
        segment: position.segment,
        files,
        late_prefix: None,
        stdout: false,
        command: None,
        argv: Vec::new(),
        output: None,
        notes: Vec::new(),
        formats: Vec::new(),
        write: None,
        vcf: Some(options),
        inputs: Vec::new(),
    };
    Ok((planned, uses_prefix))
}

/// A tap of a variant stream: its Picard command line with the stream as
/// its VCF (`INPUT`, or `CALL_VCF` of `genotype-concordance`) and `OUTPUT`
/// the prefix of its files, its options checked as the tool checks them.
pub(super) fn tap_step(
    cx: &Context,
    mut draft: Draft,
    verb: &'static Verb,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let tool = verb.tool();
    let mut formats = step_formats(cx, &draft)?;
    verb.check_one_spelling(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    check_types(&draft)?;
    let call_vcf = tool.option("CALL_VCF").is_some();
    if call_vcf {
        if let Some(given) = draft
            .parsed
            .given
            .iter()
            .find(|given| given.key == Key::Picard("CALL_VCF"))
        {
            return Err(draft
                .usage(format!(
                    "{} names the call VCF, and the chain's stream is the call VCF",
                    given.spelled
                ))
                .from(given.spelled.clone()));
        }
    } else {
        draft.parsed.positional.insert(0, cx.source.to_owned());
    }
    let picard = verb
        .picard_arguments(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    let own_output = draft.parsed.last(Key::Output).map(|g| g.values[0].clone());
    let (prefix, uses_prefix) = match own_output.as_deref() {
        Some("-") => {
            return Err(draft
                .usage("standard output takes one file; this step writes several")
                .hint("give the step -o PREFIX"));
        }
        Some(path) => (prefix_under(path, cx.source).0, false),
        None => (cx.base(&position.reads), true),
    };
    let mut argv: Vec<String> = picard
        .into_iter()
        .filter(|arg| !arg.starts_with("OUTPUT=") && !arg.starts_with("O="))
        .collect();
    if call_vcf {
        // After `TRUTH_VCF`, as the tool's own command line orders them.
        let at = argv
            .iter()
            .position(|arg| arg.starts_with("TRUTH_VCF="))
            .map_or(0, |at| at + 1);
        argv.insert(at, format!("CALL_VCF={}", cx.source));
    }
    argv.push(format!("OUTPUT={}", prefix.display()));
    verb.check_required(&draft.flags, &argv)
        .and_then(|()| verb.check_inputs(&draft.parsed, &argv))
        .map_err(|e| draft.error(e))?;
    let mut inv = tool.parse(&argv).map_err(|e| draft.core(verb, &e))?;
    inv.threads = cx.threads;
    stage::variant::VariantStep::tool_tap(tool, &inv).map_err(|e| draft.core(verb, &e))?;
    if formats.is_empty() {
        formats.push(format::Format::Picard);
    }
    let mut files = Vec::new();
    for path in stage::variant::tap_files(tool, &prefix) {
        files.push(PlannedFile {
            path: path.clone(),
            kind: "metrics",
        });
        for f in &formats[1..] {
            files.extend(
                format::files_of(&format::sibling(&path, f.extension()), *f)
                    .into_iter()
                    .map(|path| PlannedFile {
                        path,
                        kind: f.name(),
                    }),
            );
        }
    }
    let planned = Planned {
        number: draft.number,
        verb: verb.name,
        tool: Some(tool),
        role: Role::Tap,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: None,
        segment: position.segment,
        files,
        // Names the keys of the tap's files (`Run.metrics` of the Python chain).
        late_prefix: Some(prefix),
        stdout: false,
        command: Some(inv.command),
        argv,
        output: None,
        notes: inv.notes,
        formats,
        write: None,
        vcf: None,
        inputs: inv.inputs,
    };
    Ok((planned, uses_prefix))
}

//! What a [`Spec`] leaves to defaults, filled in before [`super::prepare`]
//! plans it (docs/chaining.md):
//! * The user config file and `BAMMETRICS_REFERENCE` ([`crate::config`])
//!   give the reference, the threads, the temporary directory and the
//!   sort memory the command line does not.
//! * `-o r.bam` names the result: a chain with a transform and no `write`
//!   step gets `write r.bam --index` as its last step, in front of a `qc`
//!   step, which stays last, and its metrics files, its manifest and the
//!   files of its qc step are named from `r` (`r.insert_size_metrics`,
//!   `r.manifest.json`, `r.verdict.json`). `-o r.sam` writes SAM text and
//!   `-o r.cram` a CRAM against the chain's reference (`-R`, the config
//!   file or `BAMMETRICS_REFERENCE`), each with no index. A chain that
//!   changes no record refuses such a name.
//! * A step option that names a dictionary defaults to the reference's
//!   ([`DICTIONARY_DEFAULTS`]): the `.dict` beside the FASTA, else its
//!   `.fai`, with a note that says which.

use std::path::{Path, PathBuf};

use bammetrics_core::registry::Kind as OptKind;

use super::args::{find_step, is_split, is_transform, materialize_leaves, step_parser_flags};
use super::{QC, Spec, StepSpec, WRITE, usage};
use crate::config;
use crate::error::VerbError;
use crate::grammar::{self, Key};

/// The step options that default to the reference's dictionary: the
/// Picard class and its option. The options that take the reference
/// itself (`REFERENCE_SEQUENCE`, `REFERENCE`) get the chain's reference in
/// [`super::prepare`].
pub const DICTIONARY_DEFAULTS: [(&str, &str); 1] = [("ReorderSam", "SEQUENCE_DICTIONARY")];

/// The FASTA extensions htsjdk replaces with `.dict` to find a reference's
/// dictionary (`ReferenceSequenceFileFactory.getDefaultDictionaryForReferenceSequence`).
const FASTA_EXTENSIONS: [&str; 10] = [
    ".fasta.gz",
    ".fasta",
    ".fas.gz",
    ".fas",
    ".fa.gz",
    ".fa",
    ".fna.gz",
    ".fna",
    ".txt.gz",
    ".txt",
];

/// A spec with its defaults filled in, and the notes of each step (by its
/// index) that say where a default came from.
pub(super) struct Resolved {
    pub spec: Spec,
    pub notes: Vec<(usize, String)>,
}

/// The kind of file a result name asks for.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Kind {
    Bam,
    Sam,
    Cram,
    /// `.SAM`, `.Sam` or `.CRAM`: the write step, as htsjdk's writer,
    /// writes SAM text or a CRAM only for a name ending in lower-case
    /// `.sam` or `.cram`, so such a name is refused. The lower-case
    /// extension.
    Case(&'static str),
}

/// The kind of a name that is a result, as `-o` gives it: one ending in
/// `.bam`, `.sam` or `.cram` in any case, and not a directory.
/// `.sam` or `.cram` in another case than lower is [`Kind::Case`].
fn result_of(text: &str) -> Option<Kind> {
    let lower = text.to_ascii_lowercase();
    if Path::new(text).is_dir() || text.ends_with('/') {
        return None;
    }
    if lower.ends_with(".bam") {
        return Some(Kind::Bam);
    }
    if text.ends_with(".sam") {
        return Some(Kind::Sam);
    }
    if text.ends_with(".cram") {
        return Some(Kind::Cram);
    }
    if lower.ends_with(".sam") {
        return Some(Kind::Case("sam"));
    }
    if lower.ends_with(".cram") {
        return Some(Kind::Case("cram"));
    }
    None
}

/// `-o r.bam`: the implicit `write` and the prefix `r`.
fn result_output(spec: &mut Spec) -> Result<(), VerbError> {
    let Some(text) = spec.prefix.clone() else {
        return Ok(());
    };
    let Some(result) = result_of(&text) else {
        return Ok(());
    };
    if spec.steps.iter().any(|step| step.verb == WRITE) {
        return Ok(());
    }
    // A transform, or a materialized tool that writes a BAM
    // (`merge-sam-files`): either leaves records no file holds.
    let changes = spec.steps.iter().any(|step| {
        find_step(&step.verb).is_some_and(|verb| is_transform(verb) || materialize_leaves(verb))
    });
    if !changes {
        return Err(usage(format!(
            "-o {text} names the BAM a chain writes, and no step of this chain changes the records"
        ))
        .from("-o")
        .hint(
            "-o PREFIX names the metrics files: -o out/sample gives out/sample.insert_size_metrics",
        )
        .hint(format!("write {text} as a step copies the records")));
    }
    if let Kind::Case(extension) = result {
        let stem = &text[..text.len() - extension.len() - 1];
        let (kind, hint) = if extension == "sam" {
            ("SAM text", format!("-o {stem}.sam writes SAM text"))
        } else {
            ("a CRAM", format!("-o {stem}.cram writes a CRAM against -R"))
        };
        return Err(usage(format!(
            "-o {text}: the write step writes {kind} only for a name ending in .{extension} \
             (lower case), as htsjdk does, and a BAM under any other name"
        ))
        .from("-o")
        .hint(format!("{hint}; -o {stem}.bam writes a BAM")));
    }
    let mut args = vec![text.clone()];
    if result == Kind::Bam {
        args.push("--index".to_owned());
    }
    // A qc step is the last step of its chain: the write reads the stream
    // qc's taps read, so it goes in front of it.
    let at = match spec.steps.last() {
        Some(last) if last.verb == QC => spec.steps.len() - 1,
        _ => spec.steps.len(),
    };
    spec.steps.insert(
        at,
        StepSpec {
            verb: WRITE.to_owned(),
            args,
            origin: None,
        },
    );
    let extension = if result == Kind::Cram { 5 } else { 4 };
    spec.prefix = Some(text[..text.len() - extension].to_owned());
    Ok(())
}

/// The dictionary of `reference`: its `.dict` (the FASTA extension
/// replaced, as htsjdk finds it), else its `.fai`; with the note.
fn dictionary_of(reference: &str) -> Option<(PathBuf, String)> {
    let lower = reference.to_ascii_lowercase();
    let stem = FASTA_EXTENSIONS
        .iter()
        .find(|ext| lower.ends_with(*ext))
        .map_or(reference, |ext| &reference[..reference.len() - ext.len()]);
    let dict = PathBuf::from(format!("{stem}.dict"));
    if dict.is_file() {
        let note = format!("the dictionary is the reference's: {}", dict.display());
        return Some((dict, note));
    }
    let fai = PathBuf::from(format!("{reference}.fai"));
    fai.is_file().then(|| {
        let note = format!(
            "the dictionary is the reference's index {} (no {} beside the reference): names \
             and lengths only",
            fai.display(),
            dict.display()
        );
        (fai, note)
    })
}

/// A step option of [`DICTIONARY_DEFAULTS`] the step leaves out: the
/// reference's dictionary in its place.
/// A step without the option and no reference to take it from is a usage
/// error that says both ways to give it.
fn dictionary_defaults(spec: &mut Spec) -> Result<Vec<(usize, String)>, VerbError> {
    let mut notes = Vec::new();
    let chain_reference = spec.reference.clone();
    for (index, step) in spec.steps.iter_mut().enumerate() {
        let Some(verb) = find_step(&step.verb) else {
            continue;
        };
        for (class, option) in DICTIONARY_DEFAULTS {
            if verb.class != class {
                continue;
            }
            let flags = step_parser_flags(verb);
            // A step that does not parse is reported by prepare as written.
            let Ok(parsed) = grammar::parse(verb.name, &flags, &step.args) else {
                continue;
            };
            if parsed.last(Key::Picard(option)).is_some() {
                continue;
            }
            // The step's own -R comes first: in the middle of a bare
            // chain it stays with the step.
            let own = parsed
                .last(Key::Picard("REFERENCE_SEQUENCE"))
                .map(|given| given.values[0].clone());
            let reference = own.or_else(|| chain_reference.clone());
            let flag = flags
                .iter()
                .find(|flag| flag.key == Key::Picard(option))
                .expect("the option is a flag of its verb");
            let label = format!("step {} ({})", index + 1, verb.name);
            let Some(reference) = &reference else {
                return Err(usage(format!(
                    "{label}: missing --{}: no reference to take it from",
                    flag.long
                ))
                .hint(format!(
                    "give --{} FILE on the step (--from-bam BAM copies the @SQ lines of a BAM \
                     header; --dict FILE reads a .dict)",
                    flag.long
                ))
                .hint("or -R FASTA, or reference in the config file: its .dict, else its .fai"));
            };
            let Some((path, note)) = dictionary_of(reference) else {
                return Err(usage(format!(
                    "{label}: missing --{}: {reference} has no .dict and no .fai beside it",
                    flag.long
                ))
                .hint(format!("give --{} FILE on the step", flag.long))
                .hint("samtools faidx FASTA writes the .fai; bammetrics create-sequence-dictionary the .dict"));
            };
            step.args
                .push(format!("--{}={}", flag.long, path.display()));
            notes.push((index, format!("{option}: {note}")));
        }
    }
    Ok(notes)
}

/// The config file's values (unless `--no-config`) and the variable's
/// reference, under what the command line gave.
fn config_defaults(spec: &mut Spec) -> Result<(), VerbError> {
    let config = config::load(!spec.no_config).map_err(|mut e| {
        super::VERB.clone_into(&mut e.verb);
        e
    })?;
    // The chain read the file: no step may write it.
    if let (true, Some(path)) = (config.read, &config.path) {
        spec.reads.push(path.display().to_string());
    }
    if spec.reference.is_none()
        && let Some(reference) = &config.reference
    {
        // A path that is not a file fails the chain only when a step
        // reads the reference; any other chain runs without it.
        match config::checked(super::VERB, reference) {
            Ok(path) => spec.reference = Some(path),
            Err(error) if needs_reference(spec) => return Err(error),
            Err(_) => {}
        }
    }
    if spec.threads.is_none() {
        spec.threads = config.threads.map(|t| t.value);
    }
    if spec.tmp_dir.is_none() {
        spec.tmp_dir = config.tmp_dir.map(|t| t.value);
    }
    if spec.sort_memory.is_none() {
        spec.sort_memory = config.sort_memory.map(|m| m.value);
    }
    Ok(())
}

/// True when a step of `spec` reads the chain's reference: a tap whose
/// `REFERENCE_SEQUENCE` or `REFERENCE` is a value option, a transform
/// that requires one, a step of
/// [`DICTIONARY_DEFAULTS`] without its option, any step on a CRAM
/// source, which the reference decodes, and a chain that writes a CRAM.
fn needs_reference(spec: &Spec) -> bool {
    if spec
        .source
        .as_deref()
        .is_some_and(|source| bammetrics_core::cram::is_cram(Path::new(source)))
        || writes_cram(spec)
    {
        return true;
    }
    spec.steps.iter().any(|step| {
        let Some(verb) = find_step(&step.verb) else {
            return false;
        };
        let tool = verb.tool();
        // A transform or split takes `REFERENCE_SEQUENCE` to decode a CRAM
        // (the shared `registry::REFERENCE_OPTION`), unless it cannot run
        // without it (`set-nm-md-and-uq-tags`). A tap reads it as given.
        let edits = is_transform(verb) || is_split(verb);
        let reads = ["REFERENCE", "REFERENCE_SEQUENCE"].iter().any(|name| {
            tool.option(name)
                .is_some_and(|opt| opt.kind == OptKind::Value && (opt.required || !edits))
        });
        let defaults = DICTIONARY_DEFAULTS.iter().any(|(class, option)| {
            verb.class == *class
                && grammar::parse(verb.name, &step_parser_flags(verb), &step.args)
                    .is_ok_and(|parsed| parsed.last(Key::Picard(option)).is_none())
        });
        reads || defaults
    })
}

/// True when `spec` writes a CRAM: `-o r.cram` with no `write` step of
/// its own, or a `write` step whose path ends in `.cram` or whose
/// `--output-format` is `cram`. The reference encodes it.
fn writes_cram(spec: &Spec) -> bool {
    let writes = spec.steps.iter().filter(|step| step.verb == WRITE);
    if writes.clone().next().is_none() {
        return spec
            .prefix
            .as_deref()
            .and_then(result_of)
            .is_some_and(|kind| kind == Kind::Cram);
    }
    writes.into_iter().any(|step| {
        let mut args = step.args.iter();
        let (mut named, mut format) = (false, None);
        while let Some(arg) = args.next() {
            match arg.split_once('=') {
                Some(("--output-format", value)) => format = Some(value),
                _ if arg == "--output-format" => format = args.next().map(String::as_str),
                _ => named |= arg.ends_with(".cram"),
            }
        }
        format.map_or(named, |format| format == "cram")
    })
}

/// The note the run prints before the first read when a `.bam` or `.sam`
/// prefix meets a `write` step of the chain's own: the prefix names the
/// metrics files there (D162).
fn chain_notes(spec: &Spec) -> Vec<(usize, String)> {
    let mut notes = Vec::new();
    if let Some(prefix) = spec.prefix.as_deref().filter(|p| result_of(p).is_some())
        && let Some(index) = spec.steps.iter().rposition(|step| step.verb == WRITE)
    {
        notes.push((
            index,
            format!(
                "-o {prefix} is the metrics prefix: this write step names the BAM, so metrics \
                 files are named {prefix}.<extension>"
            ),
        ));
    }
    notes
}

/// `spec` with its defaults: the config, the result BAM, the dictionaries.
pub(super) fn resolve(spec: &Spec) -> Result<Resolved, VerbError> {
    let mut spec = spec.clone();
    config_defaults(&mut spec)?;
    if let Some(name) = spec.recipe_bam.clone() {
        let names_bam = spec
            .prefix
            .as_deref()
            .and_then(result_of)
            .is_some_and(|kind| kind == Kind::Bam);
        if !names_bam {
            return Err(usage(format!(
                "recipe {name} writes a BAM: -o names it (-o out.bam)"
            ))
            .hint(format!("bammetrics {name} BAM -o out.bam")));
        }
    }
    let mut notes = chain_notes(&spec);
    result_output(&mut spec)?;
    notes.extend(dictionary_defaults(&mut spec)?);
    Ok(Resolved { spec, notes })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(prefix: &str, verbs: &[&str]) -> Spec {
        Spec {
            prefix: Some(prefix.to_owned()),
            steps: verbs.iter().map(|v| StepSpec::new(v, &[])).collect(),
            no_config: true,
            ..Spec::default()
        }
    }

    #[test]
    fn a_result_name_adds_the_write_and_names_the_prefix() {
        let mut chain = spec("out/r.bam", &["reorder-sam", "insert-size"]);
        result_output(&mut chain).expect("resolves");
        assert_eq!(chain.prefix.as_deref(), Some("out/r"));
        assert_eq!(
            chain.steps.last(),
            Some(&StepSpec::new("write", &["out/r.bam", "--index"]))
        );
        let mut sam = spec("r.sam", &["mark-duplicates"]);
        result_output(&mut sam).expect("resolves");
        assert_eq!(sam.prefix.as_deref(), Some("r"));
        assert_eq!(sam.steps.last(), Some(&StepSpec::new("write", &["r.sam"])));
        let mut cram = spec("out/r.cram", &["mark-duplicates", "insert-size"]);
        result_output(&mut cram).expect("resolves");
        assert_eq!(cram.prefix.as_deref(), Some("out/r"));
        assert_eq!(
            cram.steps.last(),
            Some(&StepSpec::new("write", &["out/r.cram"]))
        );
        assert!(writes_cram(&spec("r.cram", &["mark-duplicates"])));
        assert!(writes_cram(&cram));
        let mut format = spec("r", &["mark-duplicates", "write"]);
        format.steps[1].args = vec!["x".to_owned(), "--output-format=cram".to_owned()];
        assert!(writes_cram(&format));
        format.steps[1].args = vec![
            "x.cram".to_owned(),
            "--output-format".to_owned(),
            "bam".to_owned(),
        ];
        assert!(!writes_cram(&format));
        assert!(!writes_cram(&spec("r.bam", &["mark-duplicates"])));
        let mut upper = spec("r.BAM", &["mark-duplicates"]);
        result_output(&mut upper).expect("resolves");
        assert_eq!(
            upper.steps.last(),
            Some(&StepSpec::new("write", &["r.BAM", "--index"]))
        );
        // A prefix, or a chain with its own write, is as before.
        for (prefix, verbs) in [
            ("out/s", &["mark-duplicates"][..]),
            ("out/s.bam", &["mark-duplicates", "write"]),
        ] {
            let mut same = spec(prefix, verbs);
            result_output(&mut same).expect("resolves");
            assert_eq!(same, spec(prefix, verbs));
        }
    }

    /// A qc step stays last: the write goes in front of it. A materialized
    /// tool that writes a BAM changes the records, so `-o m.bam` writes them.
    #[test]
    fn the_write_goes_in_front_of_qc_and_a_merge_counts_as_a_change() {
        let mut chain = spec("out/r.bam", &["mark-duplicates", "insert-size", "qc"]);
        result_output(&mut chain).expect("resolves");
        let verbs: Vec<&str> = chain.steps.iter().map(|s| s.verb.as_str()).collect();
        assert_eq!(verbs, ["mark-duplicates", "insert-size", "write", "qc"]);
        assert_eq!(chain.steps[2].args, ["out/r.bam", "--index"]);
        assert_eq!(chain.prefix.as_deref(), Some("out/r"));
        let mut merge = spec("m.bam", &["merge-sam-files"]);
        result_output(&mut merge).expect("a merge changes the records");
        assert_eq!(
            merge.steps.last(),
            Some(&StepSpec::new("write", &["m.bam", "--index"]))
        );
        // Taps and qc alone change nothing.
        let error = result_output(&mut spec("r.bam", &["insert-size", "qc"])).expect_err("taps");
        assert!(
            error
                .to_string()
                .contains("no step of this chain changes the records")
        );
    }

    #[test]
    fn a_result_name_without_a_transform_or_in_upper_case_is_refused() {
        let error = |prefix: &str, verbs: &[&str]| {
            result_output(&mut spec(prefix, verbs))
                .expect_err(prefix)
                .to_string()
        };
        assert!(error("r.bam", &["insert-size"]).starts_with(
            "bammetrics chain: -o r.bam names the BAM a chain writes, and no step of this \
                 chain changes the records (from -o)\n  -o PREFIX names the metrics files"
        ));
        assert!(
            error("r.CRAM", &["mark-duplicates"])
                .contains("writes a CRAM only for a name ending in .cram (lower case)")
        );
        for name in ["r.SAM", "r.Sam"] {
            assert!(
                error(name, &["mark-duplicates"])
                    .contains("writes SAM text only for a name ending in .sam (lower case)"),
                "{name}"
            );
        }
    }

    #[test]
    fn the_dictionary_is_the_references_dict_then_its_index() {
        let dir = std::env::temp_dir().join(format!("bammetrics-defaults-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let fasta = dir.join("ref.fa");
        std::fs::write(&fasta, ">c\nACGT\n").expect("fasta");
        std::fs::write(dir.join("ref.fa.fai"), "c\t4\t3\t4\t5\n").expect("fai");
        let mut chain = spec("r.bam", &["reorder-sam", "insert-size"]);
        chain.reference = Some(fasta.display().to_string());
        let notes = dictionary_defaults(&mut chain).expect("a default");
        assert_eq!(
            chain.steps[0].args,
            [format!("--sequence-dictionary={}.fai", fasta.display())]
        );
        assert!(notes[0].1.contains("names and lengths only"), "{notes:?}");
        std::fs::write(dir.join("ref.dict"), "@HD\tVN:1.6\n").expect("dict");
        let mut chain = spec("r.bam", &["reorder-sam"]);
        chain.reference = Some(fasta.display().to_string());
        let notes = dictionary_defaults(&mut chain).expect("a default");
        assert_eq!(
            chain.steps[0].args,
            [format!(
                "--sequence-dictionary={}",
                dir.join("ref.dict").display()
            )]
        );
        assert_eq!(notes[0].0, 0);
        // The step's own dictionary stays.
        let mut own = spec("r.bam", &[]);
        own.steps = vec![StepSpec::new("reorder-sam", &["--sd", "x.dict"])];
        own.reference = Some(fasta.display().to_string());
        assert!(dictionary_defaults(&mut own).expect("its own").is_empty());
        let mut none = spec("r.bam", &["reorder-sam"]);
        let error = dictionary_defaults(&mut none)
            .expect_err("no reference")
            .to_string();
        assert!(
            error.contains("step 1 (reorder-sam): missing --sequence-dictionary: no reference")
        );
        assert_eq!(own.steps[0].args, ["--sd", "x.dict"]);
        let _ = std::fs::remove_dir_all(&dir);
    }
}

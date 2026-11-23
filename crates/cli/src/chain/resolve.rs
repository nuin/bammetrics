//! What a [`Spec`] leaves to defaults, filled in before [`super::prepare`]
//! plans it (docs/chaining.md):
//! * The user config file and `BAMMETRICS_REFERENCE` ([`crate::config`])
//!   give the reference, the threads, the temporary directory and the
//!   sort memory the command line does not.
//! * `-o r.bam` names the result: a chain with a transform and no `write`
//!   step gets `write r.bam --index` as its last step, and its metrics
//!   files are named from `r` (`r.insert_size_metrics`). `-o r.sam` writes
//!   SAM text, with no index. A chain that changes no record refuses such
//!   a name.
//! * A step option that names a dictionary defaults to the reference's
//!   ([`DICTIONARY_DEFAULTS`]): the `.dict` beside the FASTA, else its
//!   `.fai`, with a note that says which.

use std::path::{Path, PathBuf};

use bammetrics_core::stage::later_kind;

use super::args::{find_step, is_transform, step_parser_flags};
use super::{Spec, StepSpec, WRITE, usage};
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
}

/// The extension of a name that is a result, as `-o` gives it: `.bam`,
/// `.sam` or `.cram` in any case. `Err` holds the kind and wave of a name
/// the writer cannot make yet.
fn result_of(text: &str) -> Option<Result<(Kind, usize), (&'static str, &'static str)>> {
    let lower = text.to_ascii_lowercase();
    if Path::new(text).is_dir() || text.ends_with('/') {
        return None;
    }
    if lower.ends_with(".bam") {
        return Some(Ok((Kind::Bam, 4)));
    }
    if lower.ends_with(".sam") {
        return Some(Ok((Kind::Sam, 4)));
    }
    if lower.ends_with(".cram") {
        return Some(Err(later_kind(Path::new(text)).unwrap_or(("CRAM", "W15"))));
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
    let changes = spec
        .steps
        .iter()
        .any(|step| find_step(&step.verb).is_some_and(is_transform));
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
    let (kind, cut) = match result {
        Ok(found) => found,
        Err((kind, wave)) => {
            return Err(usage(format!(
                "-o {text}: bammetrics writes BAM and SAM only: {kind} output is wave {wave} of \
                 docs/07_full_parity_program.md"
            ))
            .from("-o")
            .hint("-o out.bam writes a BAM"));
        }
    };
    let mut args = vec![text.clone()];
    if kind == Kind::Bam {
        args.push("--index".to_owned());
    }
    spec.steps.push(StepSpec {
        verb: WRITE.to_owned(),
        args,
        origin: None,
    });
    spec.prefix = Some(text[..text.len() - cut].to_owned());
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
    let reference = spec.reference.clone();
    for (index, step) in spec.steps.iter_mut().enumerate() {
        let Some(verb) = find_step(&step.verb) else {
            continue;
        };
        for (class, option) in DICTIONARY_DEFAULTS {
            if verb.class != class {
                continue;
            }
            let flags = step_parser_flags(verb);
            let given = grammar::parse(verb.name, &flags, &step.args)
                .map(|parsed| parsed.last(Key::Picard(option)).is_some());
            // A step that does not parse is reported by prepare as written.
            if given.unwrap_or(true) {
                continue;
            }
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
                .hint(format!("give --{} FILE on the step", flag.long))
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
    if spec.reference.is_none()
        && let Some(reference) = &config.reference
    {
        spec.reference = Some(config::checked(super::VERB, reference)?);
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

/// `spec` with its defaults: the config, the result BAM, the dictionaries.
pub(super) fn resolve(spec: &Spec) -> Result<Resolved, VerbError> {
    let mut spec = spec.clone();
    config_defaults(&mut spec)?;
    if let Some(name) = spec.recipe_bam.clone() {
        let names_bam = spec
            .prefix
            .as_deref()
            .and_then(result_of)
            .is_some_and(|r| r.is_ok());
        if !names_bam {
            return Err(usage(format!(
                "recipe {name} writes a BAM: -o names it (-o out.bam)"
            ))
            .hint(format!("bammetrics {name} BAM -o out.bam")));
        }
    }
    result_output(&mut spec)?;
    let notes = dictionary_defaults(&mut spec)?;
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
        let mut sam = spec("r.SAM", &["mark-duplicates"]);
        result_output(&mut sam).expect("resolves");
        assert_eq!(sam.prefix.as_deref(), Some("r"));
        assert_eq!(sam.steps.last(), Some(&StepSpec::new("write", &["r.SAM"])));
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

    #[test]
    fn a_result_name_without_a_transform_or_as_cram_is_refused() {
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
            error("r.cram", &["mark-duplicates"])
                .contains("CRAM output is wave W15 of docs/07_full_parity_program.md")
        );
    }

    #[test]
    fn the_dictionary_is_the_references_dict_then_its_index() {
        let dir = std::env::temp_dir().join(format!("bammetrics-resolve-{}", std::process::id()));
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

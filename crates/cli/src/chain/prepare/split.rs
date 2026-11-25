//! A split step (`split-sam-by-library`, `split-sam-by-number-of-reads`):
//! its directory from its own `-o`, its options checked as the tool checks
//! them, and one planned name that stands for the files it makes once the
//! stream is read. Two splits into one directory with one naming rule
//! claim that name twice: a usage error before the first read.

use bammetrics_core::split_sam::{self, SplitKind};

use super::{Context, Draft, Planned, PlannedFile, Position, Role, check_types};
use crate::error::VerbError;
use crate::grammar::Key;
use crate::verbs::Verb;

/// The kind of the planned name of a split's files.
pub const KIND: &str = "split";

/// The split step `draft` of `verb`.
pub(super) fn split_step(
    cx: &Context,
    mut draft: Draft,
    verb: &'static Verb,
    position: Position,
) -> Result<(Planned, bool), VerbError> {
    let tool = verb.tool();
    verb.check_one_spelling(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    check_types(&draft)?;
    if let Some(given) = draft.parsed.last(Key::Format) {
        return Err(draft
            .usage(format!(
                "{} has nothing to write: this step writes BAM files",
                given.spelled
            ))
            .from(given.spelled.clone()));
    }

    draft.parsed.positional.push(cx.source.to_owned());
    let picard = verb
        .picard_arguments(&draft.parsed, &draft.flags)
        .map_err(|e| draft.error(e))?;
    let (input, own) = picard.split_first().expect("the source is the first input");
    let mut argv = vec![input.clone()];
    match draft.parsed.last(Key::Output) {
        Some(given) if given.values[0] == "-" => {
            return Err(draft
                .usage("standard output takes one file; this step writes a directory of them")
                .from(given.spelled.clone())
                .hint("give the step -o DIR"));
        }
        Some(given) => argv.push(format!("OUTPUT={}", given.values[0])),
        None => {}
    }
    argv.extend(own.iter().cloned());
    let mut inv = tool.parse(&argv).map_err(|e| draft.core(verb, &e))?;
    inv.threads = cx.threads;
    let (kind, dir) = split_sam::chain_options(tool, &inv).map_err(|e| {
        let error = draft.core(verb, &e);
        if matches!(e, bammetrics_core::Error::Usage(ref m) if m.starts_with("OUTPUT")) {
            error.hint("give the step -o DIR")
        } else {
            error
        }
    })?;
    split_sam::WriterOptions::from_invocation(&inv).map_err(|e| draft.core(verb, &e))?;
    let pattern = match &kind {
        SplitKind::Library => "<LIBRARY>.bam".to_owned(),
        SplitKind::Count(count) => format!("{}_<N>.bam", count.prefix),
    };
    let planned = Planned {
        number: draft.number,
        verb: verb.name,
        tool: Some(tool),
        role: Role::Split,
        label: draft.label.clone(),
        reads: position.reads,
        leaves: None,
        segment: position.segment,
        files: vec![PlannedFile {
            path: dir.join(pattern),
            kind: KIND,
        }],
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

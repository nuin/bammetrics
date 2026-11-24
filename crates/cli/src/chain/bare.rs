//! Steps without `+` (docs/chaining.md, "Steps without +"): in
//! `bammetrics chain BAM mark-duplicates insert-size -o r.bam`, a bare word
//! that names a step verb starts a step. [`with_plus`] writes such a command
//! line in the `+` form, so one parser reads both forms and a command line
//! that has a `+` before every step comes back unchanged.
//!
//! The rules, in the order they apply to each word:
//! * `+` starts a step, and the word after it is the step's verb, as before.
//! * After `--` every word up to the next `+` is positional.
//! * An option takes its values with it (`--level SAMPLE`), so a value that
//!   names a verb stays a value.
//! * The first bare word in front of the first step is the source, unless
//!   `-i`, `--input` or `I=` gave it.
//! * Any other bare word that names a step verb, a Picard class, `write` or
//!   `qc` starts a step; one a letter or two away from a verb is an unknown
//!   step with the closest verb; any other bare word stays with its step
//!   (the path of `write`, a further input of `merge-sam-files`,
//!   `gather-bam-files` or `jumping-library`).
//! * In the last step, when a bare word started it, an option of the chain
//!   goes to the chain (`-o r.bam`, `-R ref.fa`, `--threads 8`, `--also
//!   csv`): one the step has no flag for, and one the step has with the
//!   chain's meaning (its own `-o`, `--format`, `--threads`, its reference
//!   and its `TMP_DIR`). A step option of another meaning stays
//!   (`downsample-sam -R 7` is its `RANDOM_SEED`). Two steps keep their
//!   `-o` in last place: `write` without a path (`write -o out.bam`) and a
//!   split step (`split-sam-by-library -o DIR`), whose `-o` names what
//!   they write. `qc` in last place keeps its own options (`--preset`,
//!   `--report`, its `--also` of verdict forms) and hands on the rest
//!   (`-o r.bam`, `--resume`).
//!
//! A command line in the new form (a bare word started a step, or no `+`
//! and no `--plan` at all) gets errors in its own terms: a source that
//! names a step, a word that is no verb where a step or the source was
//! expected, a chain option in the middle of the chain, an option of
//! `write` on the last step, `-o` given twice. The `+` form keeps the
//! errors of the parser.

use std::path::Path;

use bammetrics_core::registry::did_you_mean;

use super::args::{
    DRY_RUN, find_step, is_chain_key, is_split, parser_flags, step_parser_flags,
    takes_further_inputs, unknown_step, write_flags,
};
use super::{QC, STEP, WRITE, qc_step, usage};
use crate::error::VerbError;
use crate::grammar::{self, Flag, Key};
use crate::verbs::VERBS;

/// The step keys whose flag means what the chain's flag of the same
/// spelling means, so the chain takes the option in last place.
const SAME_MEANING: [Key; 11] = [
    Key::Input,
    Key::Output,
    Key::Threads,
    Key::Format,
    Key::Quiet,
    Key::Verbose,
    Key::Sample,
    Key::Picard("REFERENCE_SEQUENCE"),
    Key::Picard("REFERENCE"),
    Key::Picard("TMP_DIR"),
    Key::Picard("VALIDATION_STRINGENCY"),
];

/// The form of the bare command line, for its errors.
const BARE_USAGE: &str = "bammetrics chain BAM STEP [STEP OPTIONS] [STEP ...] [-o OUT]";

/// One step as written: whether a bare word started it, its verb and the
/// words after the verb.
struct Written {
    bare: bool,
    verb: Option<String>,
    args: Vec<String>,
}

/// The flags of the step named `verb`: its verb's, `write`'s, or none for
/// a word no verb has (the parser of `prepare` reports it).
fn flags_of(verb: Option<&str>) -> Vec<Flag> {
    match verb {
        Some(WRITE) => write_flags(),
        Some(QC) => qc_step::flags(),
        Some(name) => find_step(name).map(step_parser_flags).unwrap_or_default(),
        None => Vec::new(),
    }
}

/// True when the bare word `word` starts a step: a step verb, a Picard
/// class, `write` or `qc`.
fn starts_step(word: &str) -> bool {
    word == WRITE || word == QC || find_step(word).is_some()
}

/// True for a word that may name a file or set a Picard key: the parser,
/// not this module, says what is wrong with it.
fn path_like(word: &str) -> bool {
    word.contains(['/', '.', '='])
}

/// The verb names and `write`, for the closest name to a misspelled word.
/// `qc` is left out: every word of two letters is two letters from it.
fn verb_names() -> impl Iterator<Item = &'static str> {
    VERBS.iter().map(|verb| verb.name).chain([WRITE])
}

/// The verb closest to `word`: the shortest verb `word` starts (3 letters
/// or more, so `ins` is `insert-size`), else the closest by spelling.
/// `loose` adds a verb that holds `word` (`dup` in `mark-duplicates`).
fn closest_verb(word: &str, loose: bool) -> Option<&'static str> {
    let word = word.to_ascii_lowercase();
    let long_enough = word.chars().count() >= 3;
    let shortest = |found: &dyn Fn(&str) -> bool| {
        verb_names()
            .filter(|name| found(name))
            .min_by_key(|name| name.len())
    };
    if long_enough && let Some(name) = shortest(&|name| name.starts_with(&word)) {
        return Some(name);
    }
    if let Some(name) = did_you_mean(&word, verb_names()) {
        return Some(name);
    }
    if loose && long_enough {
        return shortest(&|name| name.contains(&word))
            .or_else(|| common_prefix(&word, verb_names()));
    }
    None
}

/// The name of `names` that shares the longest start with `word`, 4
/// letters or more (`mark-dups` and `mark-duplicates`).
pub(super) fn common_prefix<'a>(
    word: &str,
    names: impl Iterator<Item = &'a str>,
) -> Option<&'a str> {
    let shared = |name: &str| {
        word.bytes()
            .zip(name.bytes())
            .take_while(|(a, b)| a.eq_ignore_ascii_case(b))
            .count()
    };
    names
        .map(|name| (shared(name), name))
        .filter(|(count, _)| *count >= 4)
        .max_by_key(|(count, name)| (*count, std::cmp::Reverse(name.len())))
        .map(|(_, name)| name)
}

/// The verb a word is a misspelling of: a bare word that names no file
/// and no verb, a letter or two from a verb name. The verb it names is
/// [`closest_verb`]'s.
fn misspelled(word: &str) -> Option<&'static str> {
    if word.contains(['/', '.']) {
        return None;
    }
    did_you_mean(word, verb_names())?;
    closest_verb(word, false)
}

/// The option `args[at]` with the words it takes, from the step's flags,
/// then the chain's: the flag and the count of words.
fn option_at<'f>(
    step: &'f [Flag],
    chain: &'f [Flag],
    args: &[String],
    at: usize,
) -> Option<(&'f Flag, usize)> {
    let next = args.get(at + 1).map(String::as_str);
    grammar::option_of(step, &args[at], next).or_else(|| grammar::option_of(chain, &args[at], next))
}

/// A bare word that is neither the source nor a step's verb: where it
/// stands (0 in front of the first step, else the step's number) and the
/// word.
struct Stray {
    step: usize,
    word: String,
}

/// The command line in the `+` form: every step behind a `+`, the chain's
/// options of the last bare step moved in front of the first.
pub(super) fn with_plus(args: &[String]) -> Result<Vec<String>, VerbError> {
    let chain = parser_flags();
    let mut head: Vec<String> = Vec::new();
    let mut steps: Vec<Written> = Vec::new();
    let mut source = false;
    let mut source_word: Option<String> = None;
    let mut strays: Vec<Stray> = Vec::new();
    let mut literal = false;
    // The word before was an option no flag has: the word after it may be
    // its value, so it starts no step.
    let mut after_unknown = false;
    let mut at = 0;
    while at < args.len() {
        let arg = &args[at];
        at += 1;
        let value_of_unknown = std::mem::take(&mut after_unknown);
        if arg == STEP {
            steps.push(Written {
                bare: false,
                verb: None,
                args: Vec::new(),
            });
            literal = false;
            continue;
        }
        if let Some(step) = steps.last_mut()
            && step.verb.is_none()
            && !step.bare
        {
            step.verb = Some(arg.clone());
            continue;
        }
        let words = match steps.last_mut() {
            Some(step) => &mut step.args,
            None => &mut head,
        };
        if literal || arg == "--" {
            literal = true;
            words.push(arg.clone());
            continue;
        }
        let step_flags = flags_of(steps.last().and_then(|s| s.verb.as_deref()));
        let in_head = steps.is_empty();
        let own: &[Flag] = if in_head { &chain } else { &step_flags };
        if let Some((flag, takes)) = option_at(own, &chain, args, at - 1) {
            if in_head && flag.key == Key::Input {
                source = true;
            }
            let words = match steps.last_mut() {
                Some(step) => &mut step.args,
                None => &mut head,
            };
            words.push(arg.clone());
            let end = (at + takes).min(args.len());
            words.extend(args[at..end].iter().cloned());
            at = end;
            continue;
        }
        let bare = !arg.starts_with('-') || arg == "-";
        if in_head && !source && bare {
            source = true;
            source_word = Some(arg.clone());
            head.push(arg.clone());
            continue;
        }
        after_unknown = !bare && !arg.contains('=');
        if bare && arg != "-" && !value_of_unknown && starts_step(arg) {
            steps.push(Written {
                bare: true,
                verb: Some(arg.clone()),
                args: Vec::new(),
            });
            continue;
        }
        // A step that takes words of its own: the path of `write`, the
        // further inputs of a step that reads several (`merge-sam-files
        // b.bam`).
        let is_write = steps
            .last()
            .and_then(|s| s.verb.as_deref())
            .is_some_and(|verb| verb == WRITE || find_step(verb).is_some_and(takes_further_inputs));
        if bare
            && !is_write
            && !value_of_unknown
            && let Some(best) = misspelled(arg)
        {
            let number = steps.len() + 1;
            let mut error = unknown_step(arg, Some(&format!("step {number}")));
            error.hints.retain(|hint| !hint.starts_with("did you mean"));
            error.hints.insert(0, format!("did you mean {best}?"));
            return Err(error);
        }
        if bare && arg != "-" && !is_write && !value_of_unknown {
            strays.push(Stray {
                step: steps.len(),
                word: arg.clone(),
            });
        }
        let words = match steps.last_mut() {
            Some(step) => &mut step.args,
            None => &mut head,
        };
        words.push(arg.clone());
    }
    let has_plan = args
        .iter()
        .any(|arg| arg == "--plan" || arg.starts_with("--plan="));
    let new_form =
        steps.iter().any(|step| step.bare) || !(args.iter().any(|a| a == STEP) || has_plan);
    if new_form {
        new_form_errors(&steps, source_word.as_deref(), &strays)?;
    }
    let number = steps.len();
    if let Some(last) = steps.last_mut().filter(|step| step.bare) {
        let (kept, moved) = split_last(last.verb.as_deref(), number, &chain, &last.args)?;
        last.args = kept;
        if new_form {
            twice(&head, &moved, &chain)?;
        }
        head.extend(moved);
    }
    if let Some(last) = steps.last_mut().filter(|step| !step.bare) {
        let (kept, moved) = dry_run_after_last(last.verb.as_deref(), &chain, &last.args);
        last.args = kept;
        head.extend(moved);
    }
    let mut out = head;
    for step in steps {
        out.push(STEP.to_owned());
        out.extend(step.verb);
        out.extend(step.args);
    }
    Ok(out)
}

/// The errors of a command line in the new form, found before the parser
/// would give its own in the terms of the `+` form.
fn new_form_errors(
    steps: &[Written],
    source: Option<&str>,
    strays: &[Stray],
) -> Result<(), VerbError> {
    if let Some(word) = source.filter(|word| starts_step(word) && !Path::new(word).exists()) {
        return Err(
            usage(format!("{word} names a step: the source comes first"))
                .from(word)
                .hint(format!("usage: {BARE_USAGE}")),
        );
    }
    for (index, step) in steps.iter().enumerate() {
        if step.bare && index + 1 < steps.len() {
            middle_chain_option(step, index + 1)?;
        }
    }
    let Some(stray) = strays.first() else {
        return Ok(());
    };
    let word = stray.word.as_str();
    let suggest = |error: VerbError| match closest_verb(word, true) {
        Some(best) => error.hint(format!("did you mean {best}?")),
        None => error,
    };
    if stray.step == 0 {
        if path_like(word) || Path::new(word).exists() {
            let count = 1 + strays.iter().filter(|s| s.step == 0).count();
            return Err(usage(format!("one source expected, got {count}"))
                .from(word)
                .hint("the source comes first, then the steps")
                .hint(format!("usage: {BARE_USAGE}")));
        }
        let error = usage(format!(
            "{word:?} is no step verb, and the chain takes one source"
        ))
        .from(word);
        return Err(suggest(error).hint(
            "the source comes first, then the steps; bammetrics tools lists the step verbs",
        ));
    }
    if path_like(word) {
        // The parser's error: a step takes no input.
        return Ok(());
    }
    let step = &steps[stray.step - 1];
    let verb = step.verb.as_deref().unwrap_or_default();
    let error = usage(format!(
        "step {} ({verb}): {word:?} is no step verb, and a step takes no positional value",
        stray.step
    ))
    .from(word);
    Err(suggest(error).hint(format!(
        "an option's value follows its option (--option VALUE); bammetrics tools lists the step \
         verbs, bammetrics chain BAM + {verb} --help the options of {verb}"
    )))
}

/// The error for an option of the chain in a bare step that is not the
/// last: the chain's options go in front of the first step or after the
/// last.
fn middle_chain_option(step: &Written, number: usize) -> Result<(), VerbError> {
    let verb = step.verb.as_deref().unwrap_or_default();
    if verb == WRITE {
        return Ok(());
    }
    let flags = flags_of(Some(verb));
    let chain = super::args::global_flags();
    let mut at = 0;
    while at < step.args.len() {
        let arg = &step.args[at];
        if arg == "--" {
            break;
        }
        let next = step.args.get(at + 1).map(String::as_str);
        let (chain_option, takes) = match grammar::option_of(&flags, arg, next) {
            Some((flag, takes)) => (is_chain_key(find_step(verb), flag.key), takes),
            None => match grammar::option_of(&chain, arg, next) {
                Some((flag, takes)) => (flag.key != Key::Help, takes),
                None => (false, 0),
            },
        };
        if chain_option {
            let spelled = arg.split_once('=').map_or(arg.as_str(), |(name, _)| name);
            return Err(usage(format!(
                "step {number} ({verb}): {spelled} belongs to the chain: give it before the first \
                 step or after the last"
            ))
            .from(spelled));
        }
        at += 1 + takes;
    }
    Ok(())
}

/// The error for `-o` given in front of the first step and again after
/// the last: the parser would keep the last without a word.
fn twice(head: &[String], moved: &[String], chain: &[Flag]) -> Result<(), VerbError> {
    let outputs = |words: &[String]| -> Vec<String> {
        let mut found = Vec::new();
        let mut at = 0;
        while at < words.len() {
            let next = words.get(at + 1).map(String::as_str);
            match grammar::option_of(chain, &words[at], next) {
                Some((flag, takes)) => {
                    if flag.key == Key::Output {
                        let value = match takes {
                            0 => words[at].split_once('=').map_or("", |(_, v)| v).to_owned(),
                            _ => words.get(at + 1).cloned().unwrap_or_default(),
                        };
                        found.push(value);
                    }
                    at += 1 + takes;
                }
                None => at += 1,
            }
        }
        found
    };
    let (before, after) = (outputs(head), outputs(moved));
    if let (Some(first), Some(last)) = (before.first(), after.last()) {
        return Err(usage(format!(
            "-o is given twice: {first} in front of the first step, {last} after the last"
        ))
        .from("-o")
        .hint("give it once"));
    }
    Ok(())
}

/// True when the words of a `write` step name its file: a word that is no
/// option nor an option's value.
fn names_path(args: &[String]) -> bool {
    let flags = write_flags();
    let mut at = 0;
    while at < args.len() {
        if args[at] == "--" {
            return at + 1 < args.len();
        }
        let next = args.get(at + 1).map(String::as_str);
        match grammar::option_of(&flags, &args[at], next) {
            Some((_, takes)) => at += 1 + takes,
            None if !args[at].starts_with('-') || args[at] == "-" => return true,
            None => at += 1,
        }
    }
    false
}

/// The words of the last step `verb` of the `+` form: those that stay
/// with it and a `--dry-run` it hands to the chain, when the step has no
/// flag of that spelling. Every other chain option stays and gets the
/// parser's error.
fn dry_run_after_last(
    verb: Option<&str>,
    chain: &[Flag],
    args: &[String],
) -> (Vec<String>, Vec<String>) {
    let step = flags_of(verb);
    let (mut kept, mut moved) = (Vec::new(), Vec::new());
    let mut at = 0;
    while at < args.len() {
        if args[at] == "--" {
            kept.extend(args[at..].iter().cloned());
            break;
        }
        let next = args.get(at + 1).map(String::as_str);
        let (target, takes) = match grammar::option_of(&step, &args[at], next) {
            Some((_, takes)) => (&mut kept, takes),
            None => match grammar::option_of(chain, &args[at], next) {
                Some((flag, takes)) if flag.key == DRY_RUN => (&mut moved, takes),
                _ => (&mut kept, 0),
            },
        };
        let end = (at + 1 + takes).min(args.len());
        target.extend(args[at..end].iter().cloned());
        at = end;
    }
    (kept, moved)
}

/// The words of the last bare step `verb`, step `number`: those that stay
/// with it and the chain's options it hands to the chain. An option of
/// `write` the step has no flag for is an error that names the `write`
/// step: the short form's write takes `--index` alone.
fn split_last(
    verb: Option<&str>,
    number: usize,
    chain: &[Flag],
    args: &[String],
) -> Result<(Vec<String>, Vec<String>), VerbError> {
    let step = flags_of(verb);
    let keeps_output = match verb {
        Some(WRITE) => !names_path(args),
        Some(name) => find_step(name).is_some_and(is_split),
        None => false,
    };
    // `--input` of a step that reads several files names a further input.
    let keeps_input = verb.and_then(find_step).is_some_and(takes_further_inputs);
    let writes = write_flags();
    let (mut kept, mut moved) = (Vec::new(), Vec::new());
    let mut at = 0;
    while at < args.len() {
        if args[at] == "--" {
            kept.extend(args[at..].iter().cloned());
            break;
        }
        let next = args.get(at + 1).map(String::as_str);
        let own = grammar::option_of(&step, &args[at], next);
        let of_chain = grammar::option_of(chain, &args[at], next).filter(|(flag, _)| {
            flag.key != Key::Help
                && own.is_none_or(|(own, _)| {
                    SAME_MEANING.contains(&own.key)
                        && !(keeps_output && own.key == Key::Output)
                        && !(keeps_input && own.key == Key::Input)
                })
        });
        if own.is_none()
            && of_chain.is_none()
            && verb.is_some_and(|v| v != WRITE)
            && let Some((flag, _)) = grammar::option_of(&writes, &args[at], next)
                .filter(|(flag, _)| !matches!(flag.key, Key::Output | Key::Help))
        {
            let verb = verb.unwrap_or_default();
            let spelled = args[at]
                .split_once('=')
                .map_or(args[at].as_str(), |(n, _)| n);
            return Err(usage(format!(
                "step {number} ({verb}): {spelled} is an option of write"
            ))
            .from(spelled)
            .hint(if verb == QC {
                format!(
                    "-o OUT.bam adds write OUT.bam --index in front of qc; for --{}, name the \
                     step there: ... write OUT.bam --{} qc",
                    flag.long, flag.long
                )
            } else {
                format!(
                    "-o OUT.bam adds write OUT.bam --index; for --{}, name the step: {verb} ... \
                     write OUT.bam --{}",
                    flag.long, flag.long
                )
            }));
        }
        let (target, takes) = match (of_chain, own) {
            (Some((_, takes)), _) => (&mut moved, takes),
            (None, Some((_, takes))) => (&mut kept, takes),
            (None, None) => (&mut kept, 0),
        };
        let end = (at + 1 + takes).min(args.len());
        target.extend(args[at..end].iter().cloned());
        at = end;
    }
    Ok((kept, moved))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|a| (*a).to_owned()).collect()
    }

    fn plus(args: &[&str]) -> Vec<String> {
        with_plus(&strings(args)).expect("splits")
    }

    #[test]
    fn a_command_line_in_the_plus_form_comes_back_unchanged() {
        for args in [
            &["IS", "-o", "OUT/s", "+", "insert-size", "-o", "x.txt"][..],
            &[
                "MD",
                "+",
                "mark-duplicates",
                "+",
                "write",
                "OUT/m.bam",
                "--index",
            ],
            &[
                "I=MD",
                "O=OUT/s",
                "+",
                "MarkDuplicates",
                "M=OUT/dup.txt",
                "+",
                "write",
                "O=x",
            ],
            &["IS", "OR", "+", "insert-size"],
            &["IS", "+", "+", "write"],
            &["IS", "+", "insert-size", "OR"],
            &["--plan", "p.toml", "-o", "OUT/s"],
            &["-", "--tmp-dir", "OUT/tmp", "--dry-run", "+", "write", "-"],
        ] {
            assert_eq!(plus(args), strings(args), "{args:?}");
        }
    }

    #[test]
    fn a_bare_verb_starts_a_step_and_the_last_hands_on_the_chains_options() {
        assert_eq!(
            plus(&["in.bam", "reorder-sam", "insert-size", "-o", "r.bam"]),
            strings(&[
                "in.bam",
                "-o",
                "r.bam",
                "+",
                "reorder-sam",
                "+",
                "insert-size"
            ])
        );
        assert_eq!(
            plus(&[
                "-R",
                "ref.fa",
                "in.bam",
                "mark-duplicates",
                "--remove-duplicates",
                "insert-size",
                "--level",
                "SAMPLE",
                "-o",
                "out/s",
                "--threads",
                "4",
                "--also",
                "csv"
            ]),
            strings(&[
                "-R",
                "ref.fa",
                "in.bam",
                "-o",
                "out/s",
                "--threads",
                "4",
                "--also",
                "csv",
                "+",
                "mark-duplicates",
                "--remove-duplicates",
                "+",
                "insert-size",
                "--level",
                "SAMPLE"
            ])
        );
        // Mixed: a + step, then bare ones; write keeps its path.
        assert_eq!(
            plus(&[
                "in.bam",
                "+",
                "mark-duplicates",
                "quality-yield",
                "write",
                "m.bam",
                "--index"
            ]),
            strings(&[
                "in.bam",
                "+",
                "mark-duplicates",
                "+",
                "quality-yield",
                "+",
                "write",
                "m.bam",
                "--index"
            ])
        );
    }

    /// `qc` is a step word: in last place it keeps its own options and hands
    /// the chain's on. A file after `merge-sam-files` stays with the step,
    /// and so does its `--input`.
    #[test]
    fn qc_and_a_step_with_further_inputs_in_the_bare_form() {
        assert_eq!(
            plus(&[
                "in.bam",
                "mark-duplicates",
                "--keep",
                "qc",
                "--preset",
                "wgs",
                "--also",
                "junit",
                "-o",
                "r.bam",
                "--resume",
                "-R",
                "ref.fa"
            ]),
            strings(&[
                "in.bam",
                "-o",
                "r.bam",
                "--resume",
                "-R",
                "ref.fa",
                "+",
                "mark-duplicates",
                "--keep",
                "+",
                "qc",
                "--preset",
                "wgs",
                "--also",
                "junit"
            ])
        );
        assert_eq!(
            plus(&[
                "in.bam",
                "merge-sam-files",
                "second",
                "third.bam",
                "insert-size",
                "-o",
                "m.bam"
            ]),
            strings(&[
                "in.bam",
                "-o",
                "m.bam",
                "+",
                "merge-sam-files",
                "second",
                "third.bam",
                "+",
                "insert-size"
            ])
        );
        assert_eq!(
            plus(&[
                "in.bam",
                "merge-sam-files",
                "--input",
                "b.bam",
                "-o",
                "m.bam"
            ]),
            strings(&[
                "in.bam",
                "-o",
                "m.bam",
                "+",
                "merge-sam-files",
                "--input",
                "b.bam"
            ])
        );
        // A word of two letters is no misspelled qc.
        assert_eq!(
            plus(&["IS", "OR", "+", "insert-size"]),
            strings(&["IS", "OR", "+", "insert-size"])
        );
    }

    #[test]
    fn a_value_that_names_a_verb_stays_a_value() {
        assert_eq!(
            plus(&[
                "in.bam",
                "mark-duplicates",
                "--as",
                "insert-size",
                "quality-yield"
            ]),
            strings(&[
                "in.bam",
                "+",
                "mark-duplicates",
                "--as",
                "insert-size",
                "+",
                "quality-yield"
            ])
        );
        // A step option of another meaning stays with the last step.
        assert_eq!(
            plus(&[
                "in.bam",
                "downsample-sam",
                "-R",
                "7",
                "-P",
                "0.5",
                "-o",
                "d.bam"
            ]),
            strings(&[
                "in.bam",
                "-o",
                "d.bam",
                "+",
                "downsample-sam",
                "-R",
                "7",
                "-P",
                "0.5"
            ])
        );
        // In the middle a step keeps its own -o; after -- every word stays.
        assert_eq!(
            plus(&[
                "in.bam",
                "insert-size",
                "-o",
                "x.txt",
                "quality-yield",
                "--",
                "write"
            ]),
            strings(&[
                "in.bam",
                "+",
                "insert-size",
                "-o",
                "x.txt",
                "+",
                "quality-yield",
                "--",
                "write"
            ])
        );
    }

    #[test]
    fn a_misspelled_step_names_its_number_and_the_closest_verb() {
        let error = with_plus(&strings(&["in.bam", "insert-size", "quality-yeld"]))
            .expect_err("misspelled");
        assert!(
            error
                .to_string()
                .starts_with("bammetrics chain: step 2: unknown step \"quality-yeld\"\n  did you mean quality-yield?"),
            "{error}"
        );
        let error = with_plus(&strings(&["in.bam", "inser-size"])).expect_err("misspelled");
        assert!(
            error
                .to_string()
                .contains("step 1: unknown step \"inser-size\"")
        );
        // A path is no misspelling: a second source in the new form, the
        // parser's own error in the + form.
        assert!(
            error_of(&["in.bam", "other.bam", "insert-size"])
                .starts_with("bammetrics chain: one source expected, got 2 (from other.bam)")
        );
        assert_eq!(
            plus(&["in.bam", "other.bam", "+", "insert-size"]),
            strings(&["in.bam", "other.bam", "+", "insert-size"])
        );
        // A prefix ranks first: ins is insert-size, not wgs.
        assert!(error_of(&["in.bam", "insert-size", "ins"]).contains("did you mean insert-size?"));
    }

    fn error_of(args: &[&str]) -> String {
        with_plus(&strings(args)).expect_err("an error").to_string()
    }

    #[test]
    fn a_word_that_is_no_verb_gets_an_error_in_the_terms_of_the_bare_form() {
        let error = error_of(&["in.bam", "insert-size", "mark", "-o", "x"]);
        assert!(
            error.starts_with(
                "bammetrics chain: step 1 (insert-size): \"mark\" is no step verb, and a step \
                 takes no positional value (from mark)\n  did you mean mark-duplicates?"
            ),
            "{error}"
        );
        assert!(
            error_of(&["in.bam", "sort-sam", "coordinate"])
                .contains("\"coordinate\" is no step verb")
        );
        let error = error_of(&["in.bam", "bogus"]);
        assert!(error.contains("\"bogus\" is no step verb, and the chain takes one source"));
        assert!(!error.contains('+'), "{error}");
        assert!(
            error_of(&["insert-size", "in.bam"])
                .starts_with("bammetrics chain: insert-size names a step: the source comes first")
        );
        // The + form keeps the parser's errors.
        assert_eq!(
            plus(&["in.bam", "+", "sort-sam", "coordinate"]),
            strings(&["in.bam", "+", "sort-sam", "coordinate"])
        );
        assert_eq!(
            plus(&["IS", "OR", "+", "insert-size"]),
            strings(&["IS", "OR", "+", "insert-size"])
        );
    }

    #[test]
    fn write_and_split_keep_their_output_in_last_place() {
        assert_eq!(
            plus(&["in.bam", "mark-duplicates", "write", "-o", "x.bam"]),
            strings(&[
                "in.bam",
                "+",
                "mark-duplicates",
                "+",
                "write",
                "-o",
                "x.bam"
            ])
        );
        // With its path, write hands -o to the chain.
        assert_eq!(
            plus(&["in.bam", "mark-duplicates", "write", "m.bam", "-o", "out/s"]),
            strings(&[
                "in.bam",
                "-o",
                "out/s",
                "+",
                "mark-duplicates",
                "+",
                "write",
                "m.bam"
            ])
        );
        assert_eq!(
            plus(&[
                "in.bam",
                "split-sam-by-library",
                "-o",
                "dir",
                "--threads",
                "2"
            ]),
            strings(&[
                "in.bam",
                "--threads",
                "2",
                "+",
                "split-sam-by-library",
                "-o",
                "dir"
            ])
        );
    }

    #[test]
    fn chain_options_in_the_middle_write_options_last_and_two_outputs_are_errors() {
        assert!(
            error_of(&["in.bam", "insert-size", "-q", "alignment-summary"]).starts_with(
                "bammetrics chain: step 1 (insert-size): -q belongs to the chain: give it before \
                 the first step or after the last (from -q)"
            )
        );
        assert!(
            error_of(&["in.bam", "insert-size", "--also=csv", "alignment-summary"])
                .contains("--also belongs to the chain")
        );
        let error = error_of(&["in.bam", "mark-duplicates", "-o", "x.bam", "--md5"]);
        assert!(
            error.starts_with(
                "bammetrics chain: step 1 (mark-duplicates): --md5 is an option of write (from --md5)"
            ),
            "{error}"
        );
        assert!(error.contains("write OUT.bam --md5"));
        assert!(
            error_of(&["-o", "a", "in.bam", "insert-size", "-o", "b"])
                .starts_with("bammetrics chain: -o is given twice: a in front of the first step, b after the last")
        );
    }
}

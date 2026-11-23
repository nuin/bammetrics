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
//! * Any other bare word that names a step verb, a Picard class or `write`
//!   starts a step; one a letter or two away from a verb is an unknown
//!   step with the closest verb; any other bare word stays with its step
//!   (the path of `write`).
//! * In the last step, when a bare word started it, an option of the chain
//!   goes to the chain (`-o r.bam`, `-R ref.fa`, `--threads 8`, `--also
//!   csv`): one the step has no flag for, and one the step has with the
//!   chain's meaning (its own `-o`, `--format`, `--threads`, its reference
//!   and its `TMP_DIR`). A step option of another meaning stays
//!   (`downsample-sam -R 7` is its `RANDOM_SEED`). A step's own `-o` in last
//!   place needs the `+` form.

use bammetrics_core::registry::did_you_mean;

use super::args::{find_step, parser_flags, step_parser_flags, unknown_step, write_flags};
use super::{STEP, WRITE};
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
        Some(name) => find_step(name).map(step_parser_flags).unwrap_or_default(),
        None => Vec::new(),
    }
}

/// True when the bare word `word` starts a step.
fn starts_step(word: &str) -> bool {
    word == WRITE || find_step(word).is_some()
}

/// The verb a word is a misspelling of: a bare word that names no file
/// and no verb, close to a verb name.
fn misspelled(word: &str) -> Option<&'static str> {
    if word.contains(['/', '.']) {
        return None;
    }
    let names = VERBS.iter().map(|verb| verb.name).chain([WRITE]);
    did_you_mean(word, names)
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

/// The command line in the `+` form: every step behind a `+`, the chain's
/// options of the last bare step moved in front of the first.
pub(super) fn with_plus(args: &[String]) -> Result<Vec<String>, VerbError> {
    let chain = parser_flags();
    let mut head: Vec<String> = Vec::new();
    let mut steps: Vec<Written> = Vec::new();
    let mut source = false;
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
        let is_write = steps.last().and_then(|s| s.verb.as_deref()) == Some(WRITE);
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
        let words = match steps.last_mut() {
            Some(step) => &mut step.args,
            None => &mut head,
        };
        words.push(arg.clone());
    }
    if let Some(last) = steps.last_mut().filter(|step| step.bare) {
        let flags = flags_of(last.verb.as_deref());
        let (kept, moved) = split_last(&flags, &chain, &last.args);
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

/// The words of the last bare step: those that stay with it and the
/// chain's options it hands to the chain.
fn split_last(step: &[Flag], chain: &[Flag], args: &[String]) -> (Vec<String>, Vec<String>) {
    let (mut kept, mut moved) = (Vec::new(), Vec::new());
    let mut at = 0;
    while at < args.len() {
        if args[at] == "--" {
            kept.extend(args[at..].iter().cloned());
            break;
        }
        let next = args.get(at + 1).map(String::as_str);
        let own = grammar::option_of(step, &args[at], next);
        let of_chain = grammar::option_of(chain, &args[at], next).filter(|(flag, _)| {
            flag.key != Key::Help && own.is_none_or(|(own, _)| SAME_MEANING.contains(&own.key))
        });
        let (target, takes) = match (of_chain, own) {
            (Some((_, takes)), _) => (&mut moved, takes),
            (None, Some((_, takes))) => (&mut kept, takes),
            (None, None) => (&mut kept, 0),
        };
        let end = (at + 1 + takes).min(args.len());
        target.extend(args[at..end].iter().cloned());
        at = end;
    }
    (kept, moved)
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
        // A path is no misspelling: the parser's own error follows.
        assert_eq!(
            plus(&["in.bam", "other.bam", "insert-size"]),
            strings(&["in.bam", "other.bam", "+", "insert-size"])
        );
    }
}

//! The verb grammar of docs/05 section 2.2: one parser for every tool verb
//! and `collect`.
//!
//! A verb's flags come from the registry option table plus the shared
//! options. One option has several spellings, all read the same way:
//! `--minimum-pct 0.1`, `--minimum-pct=0.1`, `--MINIMUM_PCT 0.1`,
//! `MINIMUM_PCT=0.1`, `-M 0.1` and `M=0.1`. A name is normalised (upper
//! case, `-` as `_`) and compared with the row's Picard names, so the
//! registry parser then sees `MINIMUM_PCT=0.1`.
//!
//! A one-letter spelling takes its value attached as clap read it in 0.12:
//! `-Oout.txt`, `-M0.1`, `-@4` (and `-O=out.txt`).
//!
//! A boolean takes `--x` (true), `--no-x`, `--x=false`, `X=false`, or the
//! space form `--x false` of the 0.12 subcommands. Only a literal `true` or
//! `false` (any case) after a boolean is read as its value, so a BAM path
//! after a flag stays a positional argument. [`Parsed`] records both
//! cases, so an error about the input can say what the flag did.

use bammetrics_core::registry::{Ty, did_you_mean, normalize};

use crate::error::VerbError;

/// What a flag sets.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Key {
    /// A row of the tool's option table, by its Picard long name.
    Picard(&'static str),
    /// A verb's own option, by its long name (`collect --metric`).
    Own(&'static str),
    Input,
    Output,
    Threads,
    Sample,
    Format,
    Quiet,
    Verbose,
    Help,
}

/// One flag of a verb.
#[derive(Clone, Debug)]
pub struct Flag {
    pub key: Key,
    /// The `--long` spelling, kebab case.
    pub long: String,
    /// Further `--long` spellings `--help` lists: the kebab form of each
    /// longer Picard name (`--level` for `LEVEL`).
    pub aliases: Vec<String>,
    /// Spellings with one dash, compared as given: `o`, `@`, `q`, and each
    /// one-letter Picard name (`M`).
    pub shorts: Vec<&'static str>,
    /// Picard names and aliases, compared after normalising.
    pub names: Vec<&'static str>,
    pub ty: Ty,
    /// Values each occurrence takes: 1, or 2 for `--hs BAITS TARGETS`.
    pub values: usize,
    pub repeated: bool,
    pub required: bool,
    pub help: String,
    /// What `--help` prints after the flag (`FILE`, `BAITS TARGETS`).
    pub placeholder: String,
}

impl Flag {
    /// A flag that is not a table row.
    pub fn new(key: Key, long: &str, shorts: &'static [&'static str], ty: Ty, help: &str) -> Self {
        Self {
            key,
            long: long.to_owned(),
            aliases: Vec::new(),
            shorts: shorts.to_vec(),
            names: Vec::new(),
            ty,
            values: 1,
            repeated: false,
            required: false,
            help: help.to_owned(),
            placeholder: ty.placeholder().to_owned(),
        }
    }

    /// Add Picard names: a one-letter name is also a short spelling, a
    /// longer one a `--kebab` alias.
    pub fn names(mut self, names: &[&'static str]) -> Self {
        for &name in names {
            if !self.names.contains(&name) {
                self.names.push(name);
            }
            if name.chars().count() == 1 {
                if !self.shorts.contains(&name) {
                    self.shorts.push(name);
                }
                continue;
            }
            let kebab = name.to_ascii_lowercase().replace('_', "-");
            if kebab != self.long && !self.aliases.contains(&kebab) {
                self.aliases.push(kebab);
            }
        }
        self
    }

    pub fn placeholder(mut self, placeholder: &str) -> Self {
        self.placeholder = placeholder.to_owned();
        self
    }

    pub fn repeated(mut self) -> Self {
        self.repeated = true;
        self
    }

    /// The spellings `--help` lists: `-o, --output FILE`.
    pub fn spellings(&self) -> String {
        let mut parts: Vec<String> = self.shorts.iter().map(|s| format!("-{s}")).collect();
        parts.push(format!("--{}", self.long));
        parts.extend(self.aliases.iter().map(|a| format!("--{a}")));
        let mut text = parts.join(", ");
        if !self.placeholder.is_empty() {
            text.push(' ');
            text.push_str(&self.placeholder);
        }
        text
    }

    /// The Picard spelling, `MINIMUM_PCT=`, for a table row.
    pub fn picard(&self) -> Option<&'static str> {
        match self.key {
            Key::Picard(name) => Some(name),
            _ => None,
        }
    }

    fn is_named(&self, name: &str) -> bool {
        let wanted = normalize(name);
        normalize(&self.long) == wanted
            || self.names.iter().any(|n| n.eq_ignore_ascii_case(&wanted))
    }
}

/// One occurrence of a flag: its key, the values and the spelling given.
#[derive(Clone, Debug)]
pub struct Given {
    pub key: Key,
    pub values: Vec<String>,
    pub spelled: String,
}

/// A parsed command line.
#[derive(Debug, Default)]
pub struct Parsed {
    pub positional: Vec<String>,
    pub given: Vec<Given>,
    /// A boolean given as `--x true` or `--x false` (space form): the flag
    /// as spelled and the word it took.
    pub bool_word: Option<(String, String)>,
    /// Positional arguments right after a bare boolean flag: the index in
    /// `positional` and the flag as spelled (`--x maybe` reads `maybe` as
    /// an input).
    pub after_bool: Vec<(usize, String)>,
}

impl Parsed {
    /// The last value given for `key`.
    pub fn last(&self, key: Key) -> Option<&Given> {
        self.given.iter().rev().find(|g| g.key == key)
    }

    /// Every occurrence of `key`, in order.
    pub fn all(&self, key: Key) -> impl Iterator<Item = &Given> {
        self.given.iter().filter(move |g| g.key == key)
    }

    /// Whether a boolean flag was given as true.
    pub fn flag(&self, key: Key) -> bool {
        self.last(key).is_some_and(|g| g.values[0] == "true")
    }

    /// The hint for positional argument `index` when a bare boolean flag
    /// came right before it: `--x maybe` reads `maybe` as an input.
    pub fn bool_value_hint(&self, index: usize) -> Option<String> {
        self.after_bool
            .iter()
            .find(|(i, _)| *i == index)
            .map(|(_, flag)| format!("expected true or false after {flag}"))
    }

    /// The hint for a missing input when a boolean took a literal word
    /// that may have been meant as the input file.
    pub fn bool_word_hint(&self) -> Option<String> {
        self.bool_word.as_ref().map(|(flag, word)| {
            format!("{flag} took {word} as its value; write ./{word} for an input named {word}")
        })
    }
}

/// The flag `name` (the argument without its dashes) refers to: a
/// one-dash short spelling exactly, then the long kebab spelling, then any
/// normalised name.
fn find<'a>(flags: &'a [Flag], dashes: usize, name: &str) -> Option<&'a Flag> {
    if dashes == 1
        && let Some(flag) = flags.iter().find(|f| f.shorts.contains(&name))
    {
        return Some(flag);
    }
    flags
        .iter()
        .find(|f| f.long == name)
        .or_else(|| flags.iter().find(|f| f.is_named(name)))
}

/// A one-dash argument with its value attached, as clap read one in 0.12:
/// `-Oout.txt` is `-O out.txt` and `-O=out.txt` is `-O out.txt`, when the
/// whole word names no flag and its first letter is the short spelling of
/// a flag that takes one value.
fn attached<'a>(flags: &[Flag], dashes: usize, arg: &'a str) -> Option<(&'a str, &'a str)> {
    if dashes != 1 {
        return None;
    }
    let body = &arg[1..];
    let whole = body.split_once('=').map_or(body, |(n, _)| n);
    if find(flags, 1, whole).is_some() {
        return None;
    }
    let first = body.chars().next()?;
    let (short, rest) = body.split_at(first.len_utf8());
    let takes_value = flags
        .iter()
        .any(|f| f.shorts.contains(&short) && f.ty != Ty::Bool && f.values == 1);
    let value = rest.strip_prefix('=').unwrap_or(rest);
    (takes_value && !value.is_empty()).then_some((short, value))
}

/// `true` or `false` in any case, as Picard reads a boolean.
fn bool_text(text: &str) -> Option<&'static str> {
    if text.eq_ignore_ascii_case("true") {
        Some("true")
    } else if text.eq_ignore_ascii_case("false") {
        Some("false")
    } else {
        None
    }
}

/// A bare word that reads as `NAME=VALUE`: the part before `=` has no path
/// characters, so `a/b=c.bam` stays a positional argument.
fn option_word(arg: &str) -> Option<(&str, &str)> {
    let (name, value) = arg.split_once('=')?;
    let plain = !name.is_empty()
        && name
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    plain.then_some((name, value))
}

/// The flag a name without dashes refers to, as `--name` does: the long
/// kebab spelling, then any normalised name (a plan file's keys).
pub fn lookup<'a>(flags: &'a [Flag], name: &str) -> Option<&'a Flag> {
    find(flags, 2, name)
}

/// The long spelling closest to `name`, a name no flag has, with the
/// Picard name of its flag as ` (Picard NAME)`: the did-you-mean of an
/// unknown option.
pub fn closest(flags: &[Flag], name: &str) -> Option<String> {
    let kebab = name.to_ascii_lowercase().replace('_', "-");
    let spellings = flags
        .iter()
        .flat_map(|f| std::iter::once(&f.long).chain(&f.aliases))
        .map(String::as_str);
    let best = did_you_mean(&kebab, spellings)?;
    let flag = flags
        .iter()
        .find(|f| f.long == best || f.aliases.iter().any(|a| a == best))?;
    let picard = flag
        .picard()
        .map_or_else(String::new, |p| format!(" (Picard {p})"));
    Some(format!("{best}{picard}"))
}

/// How the message of [`unknown`] starts.
const UNKNOWN_OPTION: &str = "unknown option ";

/// The option an error of [`unknown`] is about, as it was spelled
/// (`--also`, `-q`, `ALSO=`); `None` for any other error.
pub fn unknown_option(error: &VerbError) -> Option<&str> {
    error.what.strip_prefix(UNKNOWN_OPTION)
}

/// The flag an option as spelled on a command line (`--also`, `-q`,
/// `ALSO=`, `--no-dry-run`) refers to.
pub fn spelled<'a>(flags: &'a [Flag], spelled: &str) -> Option<&'a Flag> {
    let dashes = spelled.len() - spelled.trim_start_matches('-').len();
    let name = spelled[dashes..].trim_end_matches('=');
    find(flags, dashes, name).or_else(|| {
        let positive = name.strip_prefix("no-").filter(|_| dashes == 2)?;
        find(flags, 2, positive).filter(|flag| flag.ty == Ty::Bool)
    })
}

/// The error for an option no flag has, with the closest spelling.
fn unknown(verb: &str, flags: &[Flag], spelled: &str, name: &str) -> VerbError {
    let mut error = VerbError::usage(verb, format!("{UNKNOWN_OPTION}{spelled}"));
    if let Some(best) = closest(flags, name) {
        error = error.hint(format!("did you mean --{best}?"));
    }
    error.hint(format!("every option: bammetrics help {verb}"))
}

/// Parse `args` against `flags`. Arguments that are not options are
/// positional; `--` makes every later argument positional.
pub fn parse(verb: &str, flags: &[Flag], args: &[String]) -> Result<Parsed, VerbError> {
    let mut parsed = Parsed::default();
    let mut rest = args.iter().peekable();
    // The spelling of a boolean flag given bare as the previous argument.
    let mut bare_bool: Option<String> = None;
    while let Some(arg) = rest.next() {
        let previous_bool = bare_bool.take();
        if arg == "--" {
            parsed.positional.extend(rest.by_ref().cloned());
            break;
        }
        let dashes = arg.len() - arg.trim_start_matches('-').len();
        let short_value = attached(flags, dashes, arg);
        let (name, inline) = if let Some((short, value)) = short_value {
            (short, Some(value))
        } else if arg == "-" || dashes > 2 {
            parsed.positional.push(arg.clone());
            continue;
        } else if dashes > 0 {
            match arg[dashes..].split_once('=') {
                Some((n, v)) => (n, Some(v)),
                None => (&arg[dashes..], None),
            }
        } else if let Some((n, v)) = option_word(arg) {
            (n, Some(v))
        } else {
            if let Some(flag) = previous_bool {
                parsed.after_bool.push((parsed.positional.len(), flag));
            }
            parsed.positional.push(arg.clone());
            continue;
        };
        let spelled = if short_value.is_some() {
            format!("-{name}")
        } else if dashes > 0 {
            arg.split('=').next().unwrap_or(arg).to_owned()
        } else {
            format!("{name}=")
        };
        // `--no-x` sets the boolean `x` to false.
        let negated = (dashes == 2 && inline.is_none())
            .then(|| name.strip_prefix("no-"))
            .flatten()
            .and_then(|positive| find(flags, 2, positive))
            .filter(|f| f.ty == Ty::Bool);
        if let Some(flag) = negated {
            parsed.given.push(Given {
                key: flag.key,
                values: vec!["false".to_owned()],
                spelled,
            });
            continue;
        }
        let flag = find(flags, dashes, name).ok_or_else(|| unknown(verb, flags, &spelled, name))?;
        let values = if flag.ty == Ty::Bool {
            let value = match inline {
                Some(text) => bool_text(text).ok_or_else(|| {
                    VerbError::usage(verb, format!("expected true or false, got {text:?}"))
                        .from(spelled.clone())
                })?,
                None => match rest.peek().and_then(|next| bool_text(next)) {
                    Some(value) => {
                        let word = rest.next().cloned().unwrap_or_default();
                        parsed.bool_word = Some((spelled.clone(), word));
                        value
                    }
                    None => {
                        bare_bool = Some(spelled.clone());
                        "true"
                    }
                },
            };
            vec![value.to_owned()]
        } else if let Some(text) = inline {
            if flag.values > 1 {
                return Err(VerbError::usage(
                    verb,
                    format!("{spelled} takes {} values after a space", flag.values),
                ));
            }
            vec![text.to_owned()]
        } else {
            let mut values = Vec::new();
            for _ in 0..flag.values {
                let value = rest.next().ok_or_else(|| {
                    VerbError::usage(verb, "missing value")
                        .from(spelled.clone())
                        .hint(format!("{} {}", spelled, flag.placeholder))
                })?;
                values.push(value.clone());
            }
            values
        };
        parsed.given.push(Given {
            key: flag.key,
            values,
            spelled,
        });
    }
    Ok(parsed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn flags() -> Vec<Flag> {
        vec![
            Flag::new(Key::Output, "output", &["o", "O"], Ty::Path, "out").names(&["OUTPUT", "O"]),
            Flag::new(Key::Quiet, "quiet", &["q"], Ty::Bool, "quiet"),
            Flag::new(
                Key::Picard("MINIMUM_PCT"),
                "minimum-pct",
                &[],
                Ty::Float,
                "pct",
            )
            .names(&["MINIMUM_PCT", "M"]),
            Flag::new(
                Key::Picard("CLIP_OVERLAPPING_READS"),
                "clip-overlapping-reads",
                &[],
                Ty::Bool,
                "clip",
            )
            .names(&["CLIP_OVERLAPPING_READS"]),
            Flag::new(Key::Own("hs"), "hs", &[], Ty::Path, "hs").placeholder("BAITS TARGETS"),
        ]
        .into_iter()
        .map(|mut f| {
            if f.long == "hs" {
                f.values = 2;
            }
            f
        })
        .collect()
    }

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| (*s).to_owned()).collect()
    }

    fn values(args: &[&str], key: Key) -> Vec<String> {
        let parsed = parse("v", &flags(), &strings(args)).expect("parses");
        parsed.all(key).flat_map(|g| g.values.clone()).collect()
    }

    #[test]
    fn every_spelling_of_one_option_is_one_option() {
        let pct = Key::Picard("MINIMUM_PCT");
        for args in [
            &["--minimum-pct", "0.1"][..],
            &["--minimum-pct=0.1"],
            &["--MINIMUM_PCT", "0.1"],
            &["--MINIMUM_PCT=0.1"],
            &["MINIMUM_PCT=0.1"],
            &["minimum_pct=0.1"],
            &["-M", "0.1"],
            &["M=0.1"],
        ] {
            assert_eq!(values(args, pct), ["0.1"], "{args:?}");
        }
    }

    #[test]
    fn booleans_take_four_forms_and_never_a_path() {
        let clip = Key::Picard("CLIP_OVERLAPPING_READS");
        assert_eq!(values(&["--clip-overlapping-reads"], clip), ["true"]);
        assert_eq!(values(&["--no-clip-overlapping-reads"], clip), ["false"]);
        assert_eq!(values(&["--clip-overlapping-reads=false"], clip), ["false"]);
        assert_eq!(values(&["CLIP_OVERLAPPING_READS=False"], clip), ["false"]);
        assert_eq!(
            values(&["--clip-overlapping-reads", "FALSE"], clip),
            ["false"]
        );
        let parsed = parse(
            "v",
            &flags(),
            &strings(&["--clip-overlapping-reads", "true.bam", "-q"]),
        )
        .expect("parses");
        assert_eq!(parsed.positional, ["true.bam"]);
        assert!(parsed.flag(Key::Quiet));
        let err = parse("v", &flags(), &strings(&["--clip-overlapping-reads=maybe"]))
            .expect_err("not a bool");
        assert_eq!(
            err.to_string(),
            "bammetrics v: expected true or false, got \"maybe\" (from --clip-overlapping-reads)"
        );
    }

    #[test]
    fn positionals_dashes_and_paths_with_equals() {
        let parsed = parse(
            "v",
            &flags(),
            &strings(&[
                "-",
                "dir/a=b.bam",
                "--hs",
                "b.bed",
                "t.bed",
                "--",
                "--quiet",
            ]),
        )
        .expect("parses");
        assert_eq!(parsed.positional, ["-", "dir/a=b.bam", "--quiet"]);
        assert_eq!(values(&["--hs", "b", "t"], Key::Own("hs")), ["b", "t"]);
        assert_eq!(
            values(&["-o", "x", "--output=y", "O=z"], Key::Output),
            ["x", "y", "z"]
        );
    }

    #[test]
    fn one_letter_flags_take_attached_values() {
        let pct = Key::Picard("MINIMUM_PCT");
        assert_eq!(values(&["-M0.1"], pct), ["0.1"]);
        assert_eq!(values(&["-M=0.1"], pct), ["0.1"]);
        assert_eq!(values(&["-Oout.txt"], Key::Output), ["out.txt"]);
        assert_eq!(values(&["-Oa=b.txt"], Key::Output), ["a=b.txt"]);
        let parsed = parse("v", &flags(), &strings(&["-M0.1"])).expect("parses");
        assert_eq!(parsed.given[0].spelled, "-M");
        // A boolean takes no attached value; `-qx` stays unknown.
        assert!(parse("v", &flags(), &strings(&["-qx"])).is_err());
    }

    #[test]
    fn boolean_words_are_recorded_for_the_input_errors() {
        let parsed = parse(
            "v",
            &flags(),
            &strings(&["--clip-overlapping-reads", "true"]),
        )
        .expect("parses");
        assert!(parsed.positional.is_empty());
        assert_eq!(
            parsed.bool_word_hint().as_deref(),
            Some(
                "--clip-overlapping-reads took true as its value; write ./true for an input \
                 named true"
            )
        );
        let parsed = parse(
            "v",
            &flags(),
            &strings(&["in.bam", "--clip-overlapping-reads", "maybe"]),
        )
        .expect("parses");
        assert_eq!(parsed.positional, ["in.bam", "maybe"]);
        assert_eq!(parsed.bool_value_hint(0), None);
        assert_eq!(
            parsed.bool_value_hint(1).as_deref(),
            Some("expected true or false after --clip-overlapping-reads")
        );
    }

    #[test]
    fn unknown_options_suggest_the_closest_flag() {
        let err = parse("insert-size", &flags(), &strings(&["--minmum-pct", "0.1"]))
            .expect_err("unknown");
        assert_eq!(err.code, 2);
        assert_eq!(
            err.to_string(),
            "bammetrics insert-size: unknown option --minmum-pct\n  \
             did you mean --minimum-pct (Picard MINIMUM_PCT)?\n  \
             every option: bammetrics help insert-size"
        );
        let err = parse("v", &flags(), &strings(&["MINIMUM_PCX=1"])).expect_err("unknown");
        assert!(
            err.to_string()
                .starts_with("bammetrics v: unknown option MINIMUM_PCX=")
        );
        let mut flags = flags();
        flags.push(
            Flag::new(
                Key::Picard("SEQUENCE_DICTIONARY"),
                "sequence-dictionary",
                &[],
                Ty::Path,
                "dict",
            )
            .names(&["SEQUENCE_DICTIONARY", "FROM_BAM", "DICT"]),
        );
        for (typo, alias) in [("--from-bm", "--from-bam"), ("--dic", "--dict")] {
            let err = parse("v", &flags, &strings(&[typo, "x"])).expect_err("unknown");
            assert!(
                err.to_string().contains(&format!(
                    "did you mean {alias} (Picard SEQUENCE_DICTIONARY)?"
                )),
                "{err}"
            );
        }
        let err = parse("v", &flags, &strings(&["--minimum-pct"])).expect_err("no value");
        assert_eq!(
            err.to_string(),
            "bammetrics v: missing value (from --minimum-pct)\n  --minimum-pct NUM"
        );
    }
}

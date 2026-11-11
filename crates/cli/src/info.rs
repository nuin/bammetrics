//! The info verbs of docs/05 section 2.1: `tools`, `presets`, `help <verb
//! or Picard class>` and `completions <shell>`, plus the clap command that
//! lists every verb (top-level `--help` and the completion scripts).

use bammetrics_core::registry::{Ty, did_you_mean};
use bammetrics_verdict::Preset;
use clap::{Arg, ArgAction, Command};
use clap_complete::Shell;

use crate::collect;
use crate::error::VerbError;
use crate::grammar::{Flag, Key};
use crate::verbs::{self, VERBS};

/// The verbs that are not tool verbs, with one line each.
const OTHER_VERBS: &[(&str, &str)] = &[
    (
        "collect",
        "Several metrics from one read of the BAM, one Picard file each.",
    ),
    ("tools", "List the verbs and the Picard class each runs."),
    (
        "presets",
        "List the built-in verdict presets with version and hash.",
    ),
    ("help", "Help for a verb or a Picard class."),
    (
        "completions",
        "A completion script for bash, zsh, fish, elvish or powershell.",
    ),
];

/// The verbs clap parses (`main.rs`), in `tools` order.
pub const CLAP_VERBS: [&str; 5] = ["qc", "plan", "check", "report", "convert"];

/// `text` up to its first `.`, `;` or `:` outside parentheses and before
/// a space: one line for a verb list.
fn first_clause(text: &str) -> String {
    let mut depth = 0usize;
    let mut end = text.len();
    for (i, c) in text.char_indices() {
        match c {
            '(' => depth += 1,
            ')' => depth = depth.saturating_sub(1),
            '.' | ';' | ':'
                if depth == 0 && text[i + 1..].chars().next().is_none_or(char::is_whitespace) =>
            {
                end = i;
                break;
            }
            _ => {}
        }
    }
    let clause = text[..end].split_whitespace().collect::<Vec<_>>().join(" ");
    format!("{clause}.")
}

/// A clap argument for a grammar flag, for completion scripts: every
/// short spelling (`-M`), every `--alias` (`--level`).
fn argument(flag: &Flag) -> Arg {
    let mut arg = Arg::new(flag.long.clone()).long(flag.long.clone());
    let mut shorts = flag.shorts.iter().filter_map(|s| s.chars().next());
    if let Some(first) = shorts.next() {
        arg = arg.short(first);
    }
    arg = arg.visible_short_aliases(shorts);
    arg = arg.visible_aliases(flag.aliases.clone());
    if flag.ty == Ty::Bool {
        return arg.action(ArgAction::SetTrue).help(flag.help.clone());
    }
    let arg = arg
        .num_args(flag.values)
        .action(ArgAction::Append)
        .help(flag.help.clone());
    if flag.ty == Ty::Path {
        arg.value_hint(clap::ValueHint::FilePath)
    } else {
        arg
    }
}

/// A clap subcommand for a verb's flags: the input positional and as
/// `-i`, and `--no-x` beside each boolean `--x`.
fn subcommand(name: &'static str, about: String, long_about: String, flags: &[Flag]) -> Command {
    let mut command = Command::new(name)
        .about(about)
        .long_about(long_about)
        .disable_help_flag(true);
    if flags.iter().any(|f| f.key == Key::Input) {
        command = command.arg(
            Arg::new("INPUT")
                .num_args(0..)
                .value_hint(clap::ValueHint::FilePath),
        );
    }
    for flag in flags {
        command = command.arg(argument(flag));
        let negated = format!("no-{}", flag.long);
        if flag.ty == Ty::Bool && flag.key != Key::Help && !flags.iter().any(|f| f.long == negated)
        {
            command = command.arg(
                Arg::new(negated.clone())
                    .long(negated)
                    .action(ArgAction::SetTrue)
                    .help(format!("--{} false.", flag.long)),
            );
        }
    }
    command
}

/// `base` (the clap verbs `qc`, `plan`, `check`, `report`) with every
/// grammar verb added, so top-level `--help` and the completion scripts
/// list all of them.
pub fn command(base: Command) -> Command {
    let mut command = base.disable_help_subcommand(true);
    for name in CLAP_VERBS {
        command = command.mut_subcommand(name, |sub| {
            let about = sub.get_about().map(ToString::to_string).unwrap_or_default();
            let long = sub
                .get_long_about()
                .map_or_else(|| about.clone(), ToString::to_string);
            sub.about(first_clause(&about)).long_about(long)
        });
    }
    for verb in VERBS {
        let tool = verb.tool();
        command = command.subcommand(subcommand(
            verb.name,
            format!("{}: {}", tool.name, first_clause(tool.summary)),
            format!("{}: {}", tool.name, tool.summary),
            &verb.flags(),
        ));
    }
    command = command.subcommand(subcommand(
        "collect",
        OTHER_VERBS[0].1.to_owned(),
        OTHER_VERBS[0].1.to_owned(),
        &collect::flags(),
    ));
    command
        .subcommand(Command::new("tools").about(OTHER_VERBS[1].1))
        .subcommand(Command::new("presets").about(OTHER_VERBS[2].1))
        .subcommand(
            Command::new("help")
                .about(OTHER_VERBS[3].1)
                .arg(Arg::new("TOPIC")),
        )
        .subcommand(
            Command::new("completions").about(OTHER_VERBS[4].1).arg(
                Arg::new("SHELL")
                    .required(true)
                    .value_parser(clap::value_parser!(Shell)),
            ),
        )
}

/// `bammetrics tools`: the tool verbs with their Picard class, then the
/// other verbs with what they do (`base` holds the clap verbs' text).
pub fn tools(base: &Command) -> String {
    let width = VERBS.iter().map(|v| v.name.len()).max().unwrap_or(0);
    let mut text = format!("{:<width$}  PICARD CLASS\n", "VERB");
    for verb in VERBS {
        text.push_str(&format!("{:<width$}  {}\n", verb.name, verb.class));
    }
    text.push_str("\nOther verbs:\n");
    for (name, about) in OTHER_VERBS {
        text.push_str(&format!("  {name:<12} {about}\n"));
    }
    for name in CLAP_VERBS {
        let about = base
            .find_subcommand(name)
            .and_then(Command::get_about)
            .map(|a| first_clause(&a.to_string()))
            .unwrap_or_default();
        text.push_str(&format!("  {name:<12} {about}\n"));
    }
    text.push_str("\nPicard syntax: bammetrics <PICARD CLASS> I=... O=... [KEY=VALUE ...]\n");
    text
}

/// The shipped presets, parsed, in `bammetrics presets` order.
pub fn shipped_presets() -> Result<Vec<Preset>, VerbError> {
    bammetrics_verdict::names()
        .into_iter()
        .map(|name| {
            Preset::builtin(name)
                .expect("a shipped name")
                .map_err(|e| VerbError::input("presets", e.to_string()))
        })
        .collect()
}

/// `bammetrics presets`.
pub fn presets() -> Result<String, VerbError> {
    let mut text = String::new();
    for preset in shipped_presets()? {
        text.push_str(&format!(
            "{:<16} v{}  {}  {} rules  {}\n",
            preset.name,
            preset.version,
            preset.hash.get(..6).unwrap_or(&preset.hash),
            preset.rules.len(),
            preset.note
        ));
    }
    Ok(text)
}

/// Every name `help` knows: verbs, Picard classes and the clap verbs.
fn topics() -> Vec<&'static str> {
    VERBS
        .iter()
        .flat_map(|v| [v.name, v.class])
        .chain(OTHER_VERBS.iter().map(|(n, _)| *n))
        .chain(CLAP_VERBS)
        .collect()
}

/// The error for a first argument that is no verb and no Picard class
/// (docs/05 section 6): the closest name, and where the list is.
pub fn unknown_verb(name: &str) -> VerbError {
    let mut error = VerbError::usage("", format!("unknown verb or Picard class {name:?}"));
    if let Some(best) = did_you_mean(name, topics()) {
        error = error.hint(format!("did you mean {best}?"));
    }
    error.hint("bammetrics tools lists them")
}

/// `bammetrics help [TOPIC]`: the help of a verb, of the verb for a Picard
/// class, or of the whole command.
pub fn help(base: Command, args: &[String]) -> Result<String, VerbError> {
    let Some(topic) = args.first() else {
        return Ok(command(base).render_help().to_string());
    };
    if let Some(verb) = verbs::find(topic).or_else(|| verbs::for_class(topic)) {
        return Ok(verb.help());
    }
    if topic == "collect" {
        return Ok(collect::help());
    }
    let mut full = command(base).bin_name("bammetrics");
    full.build();
    if let Some(sub) = full.find_subcommand_mut(topic) {
        return Ok(sub.render_long_help().to_string());
    }
    let mut error = VerbError::usage("help", format!("no verb or Picard class {topic:?}"));
    if let Some(best) = did_you_mean(topic, topics()) {
        error = error.hint(format!("did you mean {best}?"));
    }
    Err(error.hint("bammetrics tools lists them"))
}

/// `bammetrics completions <shell>`: the script on standard output, dropped
/// on a closed pipe.
pub fn completions(base: Command, args: &[String]) -> Result<(), VerbError> {
    let name = args.first().ok_or_else(|| {
        VerbError::usage("completions", "missing shell")
            .hint("bash, zsh, fish, elvish or powershell")
    })?;
    let shell: Shell = name.parse().map_err(|_| {
        VerbError::usage("completions", format!("no shell {name:?}"))
            .hint("bash, zsh, fish, elvish or powershell")
    })?;
    let mut full = command(base);
    let mut script = Vec::new();
    clap_complete::generate(shell, &mut full, "bammetrics", &mut script);
    crate::events::print_stdout(&String::from_utf8_lossy(&script));
    Ok(())
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    /// clap's own checks (one short and one long per argument) pass for the
    /// command that lists every verb, so `--help` and completions build.
    #[test]
    fn the_full_command_is_valid() {
        command(crate::Cli::command()).debug_assert();
    }

    #[test]
    fn first_clauses_stop_outside_parentheses() {
        assert_eq!(first_clause("a b. c"), "a b.");
        assert_eq!(first_clause("x (y; z) w; v"), "x (y; z) w.");
        assert_eq!(
            first_clause("PREFIX.<files>, more"),
            "PREFIX.<files>, more."
        );
        assert_eq!(first_clause("the page: PREFIX.x"), "the page.");
    }

    #[test]
    fn unknown_first_words_get_the_closest_verb_or_class() {
        assert_eq!(
            unknown_verb("inser-size").to_string(),
            "bammetrics: unknown verb or Picard class \"inser-size\"\n  \
             did you mean insert-size?\n  bammetrics tools lists them"
        );
        let class = unknown_verb("CollectInsertSizeMetric").to_string();
        assert!(class.contains("did you mean CollectInsertSizeMetrics?"));
    }
}

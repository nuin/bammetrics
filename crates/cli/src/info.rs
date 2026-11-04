//! The info verbs of docs/05 section 2.1: `tools`, `presets`, `help <verb
//! or Picard class>` and `completions <shell>`, plus the clap command that
//! lists every verb (top-level `--help` and the completion scripts).

use std::io;

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

/// A clap argument for a grammar flag, for completion scripts.
fn argument(flag: &Flag) -> Arg {
    let mut arg = Arg::new(flag.long.clone()).long(flag.long.clone());
    let mut shorts = flag.shorts.iter().filter_map(|s| s.chars().next());
    if let Some(first) = shorts.next() {
        arg = arg.short(first);
    }
    for alias in shorts {
        arg = arg.short_alias(alias);
    }
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

/// A clap subcommand for a verb's flags; the input is positional.
fn subcommand(name: &'static str, about: String, flags: &[Flag]) -> Command {
    let mut command = Command::new(name).about(about).disable_help_flag(true);
    if flags.iter().any(|f| f.key == Key::Input) {
        command = command.arg(
            Arg::new("INPUT")
                .num_args(0..)
                .value_hint(clap::ValueHint::FilePath),
        );
    }
    for flag in flags {
        if flag.key != Key::Input {
            command = command.arg(argument(flag));
        }
    }
    command
}

/// `base` (the clap verbs `qc`, `plan`, `check`, `report`) with every
/// grammar verb added, so top-level `--help` and the completion scripts
/// list all of them.
pub fn command(base: Command) -> Command {
    let mut command = base.disable_help_subcommand(true);
    for verb in VERBS {
        let tool = verb.tool();
        command = command.subcommand(subcommand(
            verb.name,
            format!("{}: {}", tool.name, tool.summary),
            &verb.flags(),
        ));
    }
    command = command.subcommand(subcommand(
        "collect",
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

/// `bammetrics tools`.
pub fn tools() -> String {
    let width = VERBS.iter().map(|v| v.name.len()).max().unwrap_or(0);
    let mut text = format!("{:<width$}  PICARD CLASS\n", "VERB");
    for verb in VERBS {
        text.push_str(&format!("{:<width$}  {}\n", verb.name, verb.class));
    }
    text.push_str("\nOther verbs:\n");
    for (name, about) in OTHER_VERBS {
        text.push_str(&format!("  {name:<12} {about}\n"));
    }
    for name in ["qc", "plan", "check", "report"] {
        text.push_str(&format!("  {name:<12} bammetrics help {name}\n"));
    }
    text.push_str("\nPicard syntax: bammetrics <PICARD CLASS> I=... O=... [KEY=VALUE ...]\n");
    text
}

/// `bammetrics presets`.
pub fn presets() -> Result<String, VerbError> {
    let mut text = String::new();
    for name in bammetrics_verdict::names() {
        let preset = Preset::builtin(name)
            .expect("a shipped name")
            .map_err(|e| VerbError::input("presets", e.to_string()))?;
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
        .chain(["qc", "plan", "check", "report"])
        .collect()
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
    let mut full = command(base);
    if let Some(sub) = full.find_subcommand_mut(topic) {
        return Ok(sub.render_help().to_string());
    }
    let mut error = VerbError::usage("help", format!("no verb or Picard class {topic:?}"));
    if let Some(best) = did_you_mean(topic, topics()) {
        error = error.hint(format!("did you mean {best}?"));
    }
    Err(error.hint("bammetrics tools lists them"))
}

/// `bammetrics completions <shell>`: the script on standard output.
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
    clap_complete::generate(shell, &mut full, "bammetrics", &mut io::stdout());
    Ok(())
}

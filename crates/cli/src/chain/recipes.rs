//! Recipes (docs/chaining.md, "Recipes"): named chains in the plan file
//! format, compiled in from `crates/cli/recipes/*.toml` and read from the
//! user's `recipes` directory beside the config file
//! (`~/.config/bammetrics/recipes/NAME.toml`). A user recipe with the name
//! of a shipped one replaces it.
//!
//! A recipe is a plan file with one more table and no source:
//!
//! ```toml
//! [recipe]
//! description = "one line for bammetrics recipes"
//! output = "bam"                                # optional: -o names the BAM
//! requires = ["downsample-sam.probability"]     # optional
//!
//! [[step]]
//! verb = "downsample-sam"
//! strategy = "ConstantMemory"
//! ```
//!
//! `bammetrics chain BAM --recipe NAME -o out.bam` and the shorthand
//! `bammetrics NAME BAM -o out.bam` run one. A word `STEP.option=value`
//! on the command line sets an option of the recipe's step `STEP`,
//! replacing the recipe's value of that option.

use std::path::{Path, PathBuf};

use bammetrics_core::registry::did_you_mean;
use bammetrics_core::toml::parse_tables;

use super::args::{find_step, step_parser_flags, write_flags};
use super::planfile;
use super::{Spec, VERB, WRITE, read_capped, usage};
use crate::error::VerbError;
use crate::grammar::{self, Flag};

/// The recipes compiled in: the name and the text.
pub const SHIPPED: [(&str, &str); 3] = [
    ("prep", include_str!("../../recipes/prep.toml")),
    ("qc-only", include_str!("../../recipes/qc-only.toml")),
    ("downsample", include_str!("../../recipes/downsample.toml")),
];

/// The most bytes of a user recipe file: 64 KiB. The longest shipped
/// recipe holds under 1 KB.
pub const FILE_CAP: u64 = 64 << 10;
/// The most `.toml` files of the user's directory that are read.
pub const MAX_USER_RECIPES: usize = 256;
/// The verb that lists the recipes.
pub const LIST_VERB: &str = "recipes";

/// One recipe, read.
#[derive(Clone, Debug)]
pub struct Recipe {
    pub name: String,
    pub description: String,
    /// `shipped`, or the file it was read from.
    pub origin: String,
    /// `-o` must name the BAM the recipe writes.
    pub output_bam: bool,
    /// `STEP.option` names the command line must set.
    pub requires: Vec<String>,
    /// The chain: no source, the steps and the `[chain]` options.
    pub spec: Spec,
}

/// The user's recipe directory: `recipes` beside the config file.
pub fn user_dir() -> Option<PathBuf> {
    crate::config::dir().map(|dir| dir.join("recipes"))
}

/// A recipe name: lower-case letters, digits and `-`, starting with a
/// letter.
fn plain_name(name: &str) -> bool {
    name.starts_with(|c: char| c.is_ascii_lowercase())
        && name
            .chars()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-')
}

/// The recipe in `text`; `file` names it in errors.
pub fn parse(name: &str, text: &str, file: &str) -> Result<Recipe, VerbError> {
    let at = |line: usize, what: String| VerbError::usage(VERB, format!("{file}:{line}: {what}"));
    let tables = parse_tables(text).map_err(|e| at(e.line, e.message))?;
    let (recipe, rest): (Vec<_>, Vec<_>) = tables.into_iter().partition(|t| t.name == "recipe");
    let table = match recipe.as_slice() {
        [one] if !one.array => one,
        [one] => {
            return Err(at(
                one.line,
                "[[recipe]]: one recipe per file, write [recipe]".into(),
            ));
        }
        [] => {
            return Err(VerbError::usage(VERB, format!("{file}: no [recipe] table"))
                .hint("a recipe starts with [recipe] and description = \"...\""));
        }
        [_, second, ..] => return Err(at(second.line, "[recipe] is given twice".into())),
    };
    let mut description = None;
    let mut output_bam = false;
    let mut requires = Vec::new();
    for entry in &table.entries {
        use bammetrics_core::toml::Scalar;
        match (entry.key.as_str(), &entry.value) {
            ("description", Scalar::Str(text)) => description = Some(text.clone()),
            ("output", Scalar::Str(kind)) if kind == "bam" => output_bam = true,
            ("output", _) => {
                return Err(at(entry.line, "output: expected \"bam\"".into()));
            }
            ("requires", Scalar::List(items)) => {
                for item in items {
                    let Scalar::Str(text) = item else {
                        return Err(at(
                            entry.line,
                            "requires: expected a list of strings".into(),
                        ));
                    };
                    requires.push(text.clone());
                }
            }
            ("description" | "requires", other) => {
                return Err(at(
                    entry.line,
                    format!("{}: unexpected {}", entry.key, other.kind()),
                ));
            }
            (other, _) => {
                let keys = ["description", "output", "requires"];
                let mut error = at(entry.line, format!("[recipe]: unknown key {other}"));
                if let Some(best) = did_you_mean(other, keys) {
                    error = error.hint(format!("did you mean {best}?"));
                }
                return Err(error.hint("keys of [recipe]: description, output, requires"));
            }
        }
    }
    let description =
        description.ok_or_else(|| at(table.line, "[recipe]: missing key description".into()))?;
    let spec = planfile::from_tables(&rest, file)?;
    if spec.source.is_some() {
        return Err(VerbError::usage(
            VERB,
            format!("{file}: a recipe names no source: the command line gives it"),
        ));
    }
    if spec.steps.is_empty() {
        return Err(VerbError::usage(
            VERB,
            format!("{file}: a recipe needs a [[step]]"),
        ));
    }
    for required in &requires {
        let (step, option) = split_target(required).ok_or_else(|| {
            VerbError::usage(
                VERB,
                format!("{file}: requires: {required:?} is not STEP.option"),
            )
        })?;
        option_flag(&spec, name, step, option)?;
    }
    Ok(Recipe {
        name: name.to_owned(),
        description,
        origin: file.to_owned(),
        output_bam,
        requires,
        spec,
    })
}

fn shipped(name: &str) -> Option<Recipe> {
    let (name, text) = SHIPPED.iter().find(|(n, _)| *n == name)?;
    let recipe = parse(name, text, &format!("recipe {name}"))
        .expect("the shipped recipes parse (recipes::tests)");
    Some(Recipe {
        origin: "shipped".to_owned(),
        ..recipe
    })
}

/// The user recipe file `name`, when there is one.
fn user_file(name: &str) -> Option<PathBuf> {
    let path = user_dir()?.join(format!("{name}.toml"));
    path.is_file().then_some(path)
}

fn read_user(name: &str, path: &Path) -> Result<Recipe, VerbError> {
    let file = path.display().to_string();
    let text = read_capped(path, FILE_CAP).map_err(|e| {
        VerbError::usage(VERB, format!("{file}: not a recipe: {e}"))
            .hint("a recipe is TOML text (docs/chaining.md)")
    })?;
    parse(name, &text, &file)
}

/// The recipe `name`: the user's file, else a shipped one; `None` when
/// neither is there.
pub fn find(name: &str) -> Result<Option<Recipe>, VerbError> {
    if !plain_name(name) {
        return Ok(None);
    }
    match user_file(name) {
        Some(path) => read_user(name, &path).map(Some),
        None => Ok(shipped(name)),
    }
}

/// True when `name` names a recipe: the shorthand verb.
pub fn exists(name: &str) -> bool {
    plain_name(name) && (SHIPPED.iter().any(|(n, _)| *n == name) || user_file(name).is_some())
}

/// The names of every recipe: the shipped ones, then the user's, sorted.
fn names() -> Vec<String> {
    let mut names: Vec<String> = SHIPPED.iter().map(|(n, _)| (*n).to_owned()).collect();
    if let Some(dir) = user_dir()
        && let Ok(entries) = std::fs::read_dir(dir)
    {
        let mut user: Vec<String> = entries
            .filter_map(Result::ok)
            .filter_map(|entry| {
                let name = entry
                    .file_name()
                    .to_str()?
                    .strip_suffix(".toml")?
                    .to_owned();
                plain_name(&name).then_some(name)
            })
            .take(MAX_USER_RECIPES)
            .collect();
        user.sort();
        for name in user {
            if !names.contains(&name) {
                names.push(name);
            }
        }
    }
    names
}

/// The error for a name that is no recipe, with the closest.
fn unknown(name: &str) -> VerbError {
    let names = names();
    let mut error = usage(format!("no recipe {name:?}")).from("--recipe");
    if let Some(best) = did_you_mean(name, names.iter().map(String::as_str)) {
        error = error.hint(format!("did you mean {best}?"));
    }
    error.hint("the recipes: bammetrics recipes")
}

/// `STEP.option` split at its first `.`.
fn split_target(text: &str) -> Option<(&str, &str)> {
    let (step, option) = text.split_once('.')?;
    (!step.is_empty() && !option.is_empty()).then_some((step, option))
}

/// The verb name a step word resolves to: the verb, or `write`.
fn verb_name(word: &str) -> Option<&'static str> {
    if word == WRITE {
        return Some(WRITE);
    }
    find_step(word).map(|verb| verb.name)
}

fn flags_of(verb: &str) -> Vec<Flag> {
    if verb == WRITE {
        return write_flags();
    }
    find_step(verb).map(step_parser_flags).unwrap_or_default()
}

/// The flag `option` names on the recipe's step `step`, checked: the
/// recipe has such a step and the step such an option.
fn option_flag(spec: &Spec, recipe: &str, step: &str, option: &str) -> Result<Flag, VerbError> {
    let wanted = verb_name(step);
    let has = spec
        .steps
        .iter()
        .any(|s| wanted.is_some() && verb_name(&s.verb) == wanted);
    if !has {
        let steps: Vec<&str> = spec.steps.iter().map(|s| s.verb.as_str()).collect();
        let mut error = usage(format!("recipe {recipe} has no step {step}"));
        if let Some(best) = did_you_mean(step, steps.iter().copied()) {
            error = error.hint(format!("did you mean {best}?"));
        }
        return Err(error.hint(format!("its steps: {}", steps.join(", "))));
    }
    let flags = flags_of(step);
    grammar::lookup(&flags, option).cloned().ok_or_else(|| {
        let mut error = usage(format!("{step} has no option {option}"));
        if let Some(best) = grammar::closest(&flags, option) {
            error = error.hint(format!("did you mean {step}.{best}?"));
        }
        error.hint(format!(
            "every option: bammetrics chain BAM + {step} --help"
        ))
    })
}

/// One `STEP.option=value` word of the command line.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Override {
    pub step: String,
    pub option: String,
    pub value: String,
}

impl Override {
    /// The word as `STEP.option=value`, when its step part names a step
    /// verb or `write`.
    pub fn of(word: &str) -> Option<Self> {
        let (target, value) = word.split_once('=')?;
        let (step, option) = split_target(target)?;
        let plain = option
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c == '-' || c == '_');
        (plain && verb_name(step).is_some()).then(|| Self {
            step: step.to_owned(),
            option: option.to_owned(),
            value: value.to_owned(),
        })
    }
}

/// The `STEP.option=value` words of `args` in front of a `--`, taken out.
pub(super) fn take_overrides(args: &[String]) -> (Vec<String>, Vec<(String, Override)>) {
    let end = args.iter().position(|a| a == "--").unwrap_or(args.len());
    let mut kept = Vec::new();
    let mut overrides = Vec::new();
    for (index, arg) in args.iter().enumerate() {
        match Override::of(arg).filter(|_| index < end) {
            Some(found) => overrides.push((arg.clone(), found)),
            None => kept.push(arg.clone()),
        }
    }
    (kept, overrides)
}

/// The recipe's chain with `overrides` applied and its required options
/// checked: each override replaces every value the recipe gives that
/// option on every step of its verb.
pub fn apply(recipe: &Recipe, overrides: &[(String, Override)]) -> Result<Spec, VerbError> {
    let mut spec = recipe.spec.clone();
    for (word, set) in overrides {
        let flag = option_flag(&spec, &recipe.name, &set.step, &set.option)
            .map_err(|e| e.from(word.clone()))?;
        let wanted = verb_name(&set.step);
        for step in spec
            .steps
            .iter_mut()
            .filter(|s| verb_name(&s.verb) == wanted)
        {
            let flags = flags_of(&step.verb);
            step.args.retain(|arg| {
                grammar::option_of(&flags, arg, None).is_none_or(|(f, _)| f.key != flag.key)
            });
            step.args.push(format!("--{}={}", flag.long, set.value));
        }
    }
    for required in &recipe.requires {
        let (step, option) = split_target(required).unwrap_or_default();
        let flag = option_flag(&spec, &recipe.name, step, option)?;
        let wanted = verb_name(step);
        let given = spec
            .steps
            .iter()
            .filter(|s| verb_name(&s.verb) == wanted)
            .all(|s| {
                let flags = flags_of(&s.verb);
                s.args.iter().any(|arg| {
                    grammar::option_of(&flags, arg, None).is_some_and(|(f, _)| f.key == flag.key)
                })
            });
        if !given {
            return Err(
                usage(format!("recipe {} needs {required}=VALUE", recipe.name))
                    .hint(format!("{required}: {}", flag.help))
                    .hint(format!(
                        "bammetrics {} BAM -o out.bam {required}=VALUE",
                        recipe.name
                    )),
            );
        }
    }
    Ok(spec)
}

/// The recipe `name` with `overrides`, or the error for a name that is no
/// recipe.
pub(super) fn load(
    name: &str,
    overrides: &[(String, Override)],
) -> Result<(Recipe, Spec), VerbError> {
    let recipe = find(name)?.ok_or_else(|| unknown(name))?;
    let spec = apply(&recipe, overrides)?;
    Ok((recipe, spec))
}

/// `bammetrics recipes`: one line per recipe, then where they come from.
pub fn list() -> Result<String, VerbError> {
    let recipes: Vec<Recipe> = names()
        .iter()
        .filter_map(|name| find(name).transpose())
        .collect::<Result<_, _>>()?;
    let width = recipes.iter().map(|r| r.name.len()).max().unwrap_or(0);
    let mut text = String::new();
    for recipe in &recipes {
        let from = if recipe.origin == "shipped" {
            ""
        } else {
            " (yours)"
        };
        text.push_str(&format!(
            "{:<width$}  {}{from}\n",
            recipe.name, recipe.description
        ));
        let steps: Vec<&str> = recipe.spec.steps.iter().map(|s| s.verb.as_str()).collect();
        let mut line = steps.join(" ");
        if recipe.output_bam {
            line.push_str(" -o OUT.bam");
        }
        for required in &recipe.requires {
            line.push_str(&format!(" {required}=VALUE"));
        }
        text.push_str(&format!("{:<width$}    steps: {line}\n", ""));
    }
    let dir = user_dir().map_or_else(|| "none".to_owned(), |d| d.display().to_string());
    text.push_str(&format!(
        "\nRun one: bammetrics NAME BAM -o OUT, or bammetrics chain BAM --recipe NAME -o OUT. \
         STEP.option=value sets an option of a step. Your recipes: {dir}/NAME.toml.\n"
    ));
    Ok(text)
}

/// `bammetrics recipes [--help]`.
pub fn run_list(args: &[String]) -> Result<u8, VerbError> {
    if let Some(first) = args.first() {
        if first == "--help" || first == "-h" {
            crate::events::print_stdout(
                "bammetrics recipes: the named chains bammetrics NAME BAM -o OUT runs, one \
                 line each with its steps. docs/chaining.md, \"Recipes\".\n\nUsage: \
                 bammetrics recipes\n",
            );
            return Ok(0);
        }
        return Err(
            VerbError::usage(LIST_VERB, format!("unexpected argument {first}"))
                .from(first.clone())
                .hint("usage: bammetrics recipes"),
        );
    }
    crate::events::print_stdout(&list()?);
    Ok(0)
}

/// `bammetrics NAME ARGS` for a recipe: `chain --recipe NAME ARGS`.
pub fn run_shorthand(name: &str, args: &[String]) -> Result<u8, VerbError> {
    let mut argv = vec!["--recipe".to_owned(), name.to_owned()];
    argv.extend(args.iter().cloned());
    super::run(&argv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_shipped_recipes_parse() {
        for (name, _) in SHIPPED {
            let recipe = shipped(name).expect("shipped");
            assert!(!recipe.description.is_empty(), "{name}");
            assert!(plain_name(name), "{name}");
        }
        let prep = shipped("prep").expect("prep");
        assert!(prep.output_bam);
        let verbs: Vec<&str> = prep.spec.steps.iter().map(|s| s.verb.as_str()).collect();
        assert_eq!(
            verbs,
            [
                "fix-mate-information",
                "mark-duplicates",
                "sort-sam",
                "quality-yield",
                "alignment-summary",
                "insert-size",
                "duplicate-metrics"
            ]
        );
        assert_eq!(
            shipped("downsample").expect("downsample").requires,
            ["downsample-sam.probability"]
        );
    }

    #[test]
    fn an_override_replaces_the_recipes_value_and_a_required_one_is_asked_for() {
        let recipe = shipped("downsample").expect("downsample");
        let error = apply(&recipe, &[]).expect_err("no probability").to_string();
        assert!(
            error.starts_with(
                "bammetrics chain: recipe downsample needs downsample-sam.probability=VALUE"
            ),
            "{error}"
        );
        let word = "downsample-sam.probability=0.25";
        let set = Override::of(word).expect("an override");
        let spec = apply(&recipe, &[(word.to_owned(), set)]).expect("applies");
        assert_eq!(
            spec.steps[0].args,
            ["--strategy=ConstantMemory", "--probability=0.25"]
        );
        let word = "DownsampleSam.strategy=Chained";
        let spec = apply(
            &recipe,
            &[(word.to_owned(), Override::of(word).expect("set"))],
        );
        assert!(spec.is_err(), "probability still missing");
        let prep = shipped("prep").expect("prep");
        let word = "sort-sam.sort-order=queryname";
        let spec = apply(
            &prep,
            &[(word.to_owned(), Override::of(word).expect("set"))],
        )
        .expect("applies");
        assert_eq!(spec.steps[2].args, ["--sort-order=queryname"]);
    }

    #[test]
    fn a_bad_override_or_recipe_is_an_error() {
        let prep = shipped("prep").expect("prep");
        let bad = |word: &str| {
            apply(&prep, &[(word.to_owned(), Override::of(word).expect(word))])
                .expect_err(word)
                .to_string()
        };
        assert!(bad("insert-size.minimum-pt=0.1").contains("did you mean insert-size.minimum-pct"));
        assert!(bad("gc-bias.scan-window-size=3").contains("recipe prep has no step gc-bias"));
        assert_eq!(Override::of("out/a.b=c"), None);
        assert_eq!(Override::of("nope.x=1"), None);
        let error = |text: &str| parse("x", text, "x.toml").expect_err(text).to_string();
        assert!(error("[[step]]\nverb = \"insert-size\"\n").contains("no [recipe] table"));
        assert!(
            error("[recipe]\n[[step]]\nverb = \"insert-size\"\n")
                .contains("missing key description")
        );
        assert!(error("[recipe]\ndescription = \"d\"\n").contains("a recipe needs a [[step]]"));
        assert!(
            error("[recipe]\ndescription = \"d\"\noutput = \"cram\"\n")
                .contains("output: expected")
        );
        assert!(error("[recipe]\ndescripton = \"d\"\n").contains("did you mean description?"));
        assert!(
            error("[recipe]\ndescription = \"d\"\n[chain]\nsource = \"a.bam\"\n[[step]]\nverb = \"write\"\n")
                .contains("a recipe names no source")
        );
    }
}

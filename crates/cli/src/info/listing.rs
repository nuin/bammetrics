//! The JSON listing of tools, options, exit codes and version
//! (`bammetrics.tools/1`, `docs/schema/tools-1.json`), built from the same
//! tables the verbs, `--help` and the Python stub are built from.

use bammetrics_core::multiple_metrics::file_extensions;
use bammetrics_core::registry::{Kind, Opt, Ty};
use bammetrics_formats::json::Json;
use clap::{ArgAction, Command};

use super::{CLAP_VERBS, OTHER_VERBS, about, command, first_clause};
use crate::grammar::{Flag, Key};
use crate::verbs::{Files, VERBS, Verb};
use crate::{chain, collect};

/// The schema name of the listing.
const SCHEMA: &str = "bammetrics.tools/1";

/// The schemas the binary writes, in the order of `docs/schema`.
const SCHEMAS: &[&str] = &[
    "bammetrics.metrics/1",
    "bammetrics.verdict/1",
    "bammetrics.manifest/1",
    "bammetrics.manifest/2",
    "bammetrics.events/1",
    SCHEMA,
];

/// Every exit code of the command line with its meaning (docs/cli.md).
pub const EXIT_CODES: &[(i32, &str)] = &[
    (
        0,
        "success: every verdict passes, or none reaches --fail-on",
    ),
    (
        1,
        "an input could not be read, a tool failed, the preset does not parse, or a qc collector failed",
    ),
    (2, "usage error, nothing ran"),
    (
        10,
        "a verdict reached --fail-on: fail, or warn with --fail-on warn; for doctor, a check",
    ),
    (11, "a verdict is incomplete and none reached --fail-on"),
    (100, "CheckTerminatorBlock found a DEFECTIVE file end"),
    (
        130,
        "interrupted by SIGINT, SIGTERM or SIGHUP; partial outputs removed",
    ),
];

fn member(key: &str, value: Json) -> (String, Json) {
    (key.to_owned(), value)
}

/// The `bammetrics` member: name, version, commit, features, schemas.
pub fn version_object() -> Json {
    let mut features = Vec::new();
    if cfg!(feature = "parquet") {
        features.push("parquet");
    }
    if cfg!(feature = "pdf") {
        features.push("pdf");
    }
    Json::Object(vec![
        member("name", Json::text("bammetrics")),
        member("version", Json::text(env!("CARGO_PKG_VERSION"))),
        member(
            "commit",
            Json::option_text(option_env!("BAMMETRICS_COMMIT")),
        ),
        member("features", Json::strings(features)),
        member("schemas", Json::strings(SCHEMAS.iter().copied())),
    ])
}

fn type_name(ty: Ty) -> &'static str {
    match ty {
        Ty::Bool => "bool",
        Ty::Int => "int",
        Ty::Float => "float",
        Ty::Text => "text",
        Ty::Path => "path",
    }
}

fn kind_name(kind: Kind) -> &'static str {
    match kind {
        Kind::Value => "value",
        Kind::Output => "output",
        Kind::Ignored(_) => "ignored",
        Kind::Rejected => "rejected",
    }
}

fn short_of(flag: &Flag) -> Json {
    Json::option_text(flag.shorts.first().map(|s| format!("-{s}")).as_deref())
}

fn option_entry(opt: &Opt, flags: &[Flag]) -> Json {
    let picard = opt.picard();
    let flag = flags
        .iter()
        .find(|f| f.key == Key::Picard(picard) || f.names.contains(&picard));
    let (long, short) = match flag {
        Some(f) => (Json::text(&format!("--{}", f.long)), short_of(f)),
        None if opt.kind == Kind::Rejected => (Json::Null, Json::Null),
        None => (Json::text(&format!("--{}", opt.kebab())), Json::Null),
    };
    Json::Object(vec![
        member("names", Json::strings(opt.names.iter().copied())),
        member("flag", long),
        member("short", short),
        member("type", Json::text(type_name(opt.ty))),
        member("kind", Json::text(kind_name(opt.kind))),
        member("repeated", Json::Bool(opt.repeated)),
        member("required", Json::Bool(opt.required)),
        member("written", Json::Bool(opt.written)),
        member("help", Json::text(opt.help)),
    ])
}

/// One tool verb: its class, files, chain role and every option row.
pub fn tool_entry(verb: &Verb) -> Json {
    let tool = verb.tool();
    let flags = verb.flags();
    let files = match verb.files {
        Files::Single => "single",
        Files::Prefix => "prefix",
        Files::Bam => "bam",
    };
    let extensions = file_extensions(tool.name)
        .unwrap_or(&[])
        .iter()
        .map(|(_, extension)| *extension);
    let role = chain::role_of(verb).map_or("none", chain::Role::name);
    let options = tool
        .table_with_rejected()
        .into_iter()
        .map(|opt| option_entry(opt, &flags))
        .collect();
    Json::Object(vec![
        member("verb", Json::text(verb.name)),
        member("class", Json::text(verb.class)),
        member("summary", Json::text(verb.about.unwrap_or(tool.summary))),
        member("input", Json::text(verb.input)),
        member("reads_stdin", Json::Bool(verb.reads_stdin)),
        member("files", Json::text(files)),
        member("writes", Json::text(verb.writes)),
        member("extensions", Json::strings(extensions)),
        member("chain_role", Json::text(role)),
        member("options", Json::Array(options)),
    ])
}

fn flag_option(flag: &Flag) -> Json {
    Json::Object(vec![
        member("flag", Json::text(&format!("--{}", flag.long))),
        member("short", short_of(flag)),
        member("takes_value", Json::Bool(flag.ty != Ty::Bool)),
        member("repeated", Json::Bool(flag.repeated)),
        member("required", Json::Bool(flag.required)),
        member("help", Json::text(&flag.help)),
    ])
}

fn verb_entry(name: &str, about: &str, options: Vec<Json>) -> Json {
    Json::Object(vec![
        member("name", Json::text(name)),
        member("about", Json::text(about)),
        member("options", Json::Array(options)),
    ])
}

/// The long options of a clap subcommand, as listing options.
fn clap_options(sub: &Command) -> Vec<Json> {
    sub.get_arguments()
        .filter(|a| a.get_long().is_some() && !matches!(a.get_id().as_str(), "help" | "version"))
        .map(|a| {
            let takes = !matches!(
                a.get_action(),
                ArgAction::SetTrue | ArgAction::SetFalse | ArgAction::Count
            );
            Json::Object(vec![
                member(
                    "flag",
                    Json::text(&format!("--{}", a.get_long().unwrap_or(""))),
                ),
                member(
                    "short",
                    Json::option_text(a.get_short().map(|c| format!("-{c}")).as_deref()),
                ),
                member("takes_value", Json::Bool(takes)),
                member(
                    "repeated",
                    Json::Bool(matches!(a.get_action(), ArgAction::Append)),
                ),
                member("required", Json::Bool(a.is_required_set())),
                member(
                    "help",
                    Json::text(&a.get_help().map(ToString::to_string).unwrap_or_default()),
                ),
            ])
        })
        .collect()
}

/// The `--format` option the info verbs `tools` and `help` take.
fn format_option() -> Json {
    Json::Object(vec![
        member("flag", Json::text("--format")),
        member("short", Json::Null),
        member("takes_value", Json::Bool(true)),
        member("repeated", Json::Bool(false)),
        member("required", Json::Bool(false)),
        member(
            "help",
            Json::text("The output form: text (default) or json."),
        ),
    ])
}

/// The `--no-config` flag every verb takes before its other arguments.
fn no_config_option() -> Json {
    Json::Object(vec![
        member("flag", Json::text(crate::config::NO_CONFIG)),
        member("short", Json::Null),
        member("takes_value", Json::Bool(false)),
        member("repeated", Json::Bool(false)),
        member("required", Json::Bool(false)),
        member("help", Json::text("Leave the user config file unread.")),
    ])
}

fn other_verb(name: &str) -> Vec<Json> {
    match name {
        "collect" => collect::flags().iter().map(flag_option).collect(),
        n if n == chain::VERB => chain::global_flags().iter().map(flag_option).collect(),
        crate::config::VERB => vec![no_config_option()],
        "tools" | "help" => vec![format_option()],
        _ => Vec::new(),
    }
}

/// The whole listing: version, exit codes, every tool, every other verb.
/// `base` is the clap command of the clap verbs (`Cli::command()`).
pub fn listing(base: &Command) -> Json {
    let full = command(base.clone());
    let exit_codes = EXIT_CODES
        .iter()
        .map(|(code, meaning)| {
            Json::Object(vec![
                member("code", Json::int(i64::from(*code))),
                member("meaning", Json::text(meaning)),
            ])
        })
        .collect();
    let mut verbs = Vec::new();
    for name in CLAP_VERBS {
        let sub = full.find_subcommand(name);
        let about = sub
            .and_then(Command::get_about)
            .map(|a| first_clause(&a.to_string()))
            .unwrap_or_default();
        verbs.push(verb_entry(
            name,
            &about,
            sub.map(clap_options).unwrap_or_default(),
        ));
    }
    for (name, _) in OTHER_VERBS {
        verbs.push(verb_entry(name, &about(name), other_verb(name)));
    }
    Json::Object(vec![
        member("schema", Json::text(SCHEMA)),
        member("bammetrics", version_object()),
        member("exit_codes", Json::Array(exit_codes)),
        member("tools", Json::Array(VERBS.iter().map(tool_entry).collect())),
        member("verbs", Json::Array(verbs)),
    ])
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    fn text<'a>(value: &'a Json, key: &str) -> &'a str {
        value.get(key).and_then(Json::as_str).expect(key)
    }

    fn tool<'a>(doc: &'a Json, verb: &str) -> &'a Json {
        doc.get("tools")
            .and_then(Json::as_array)
            .expect("tools")
            .iter()
            .find(|t| text(t, "verb") == verb)
            .expect(verb)
    }

    #[test]
    fn one_entry_per_verb_in_verbs_order() {
        let doc = listing(&crate::Cli::command());
        let tools = doc.get("tools").and_then(Json::as_array).expect("tools");
        let names: Vec<&str> = tools.iter().map(|t| text(t, "verb")).collect();
        let expected: Vec<&str> = VERBS.iter().map(|v| v.name).collect();
        assert_eq!(names, expected);
        assert_eq!(tools.len(), 77);
        let members: Vec<&str> = doc
            .as_object()
            .expect("object")
            .iter()
            .map(|(k, _)| k.as_str())
            .collect();
        assert_eq!(
            members,
            ["schema", "bammetrics", "exit_codes", "tools", "verbs"]
        );
    }

    #[test]
    fn roles_and_options_are_marked() {
        let doc = listing(&crate::Cli::command());
        assert_eq!(text(tool(&doc, "insert-size"), "chain_role"), "tap");
        assert_eq!(
            text(tool(&doc, "mark-duplicates"), "chain_role"),
            "transform"
        );
        assert_eq!(text(tool(&doc, "build-bam-index"), "chain_role"), "none");
        assert_eq!(
            tool(&doc, "jumping-library").get("reads_stdin"),
            Some(&Json::Bool(false))
        );
        let options = tool(&doc, "insert-size")
            .get("options")
            .and_then(Json::as_array)
            .expect("options");
        let histogram = options
            .iter()
            .find(|o| {
                o.get("names")
                    .and_then(Json::as_array)
                    .is_some_and(|n| n[0].as_str() == Some("HISTOGRAM_FILE"))
            })
            .expect("HISTOGRAM_FILE");
        assert_eq!(text(histogram, "kind"), "ignored");
    }

    #[test]
    fn verbs_list_clap_and_other_verbs() {
        let doc = listing(&crate::Cli::command());
        let verbs = doc.get("verbs").and_then(Json::as_array).expect("verbs");
        let names: Vec<&str> = verbs.iter().map(|v| text(v, "name")).collect();
        for expected in ["qc", "doctor", "chain", "collect", "tools", "completions"] {
            assert!(names.contains(&expected), "{expected}");
        }
        let qc = verbs.iter().find(|v| text(v, "name") == "qc").expect("qc");
        let flags: Vec<&str> = qc
            .get("options")
            .and_then(Json::as_array)
            .expect("options")
            .iter()
            .map(|o| text(o, "flag"))
            .collect();
        assert!(flags.contains(&"--preset"));
    }

    #[test]
    fn exit_codes_are_the_documented_seven() {
        let codes: Vec<i32> = EXIT_CODES.iter().map(|(c, _)| *c).collect();
        assert_eq!(codes, [0, 1, 2, 10, 11, 100, 130]);
    }
}

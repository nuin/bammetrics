# E1 Machine-Readable Interface Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the command line describe itself as JSON, so wrapper generators and language bindings read the tools, options, exit codes and version from the binary.

**Architecture:** One new module, `crates/cli/src/info/listing.rs`, builds a `bammetrics_formats::json::Json` value from `registry::TOOLS`, `verbs::VERBS` and the clap command that `info::command` already builds. `tools`, `help` and `--version` print it. A new frozen schema `docs/schema/tools-1.json` describes it, and `crates/cli/tests/schemas.rs` validates it with the existing validator.

**Tech Stack:** Rust 1.91, clap 4 (already a dependency), the hand-written JSON of `crates/formats/src/json.rs`. No new crate.

**Spec:** `docs/superpowers/specs/2026-11-25-embedding-design.md` section 5.1 and the E1 row of section 6.

## Global Constraints

- Version 0.34.1 (EQ4 default: E waves take patch versions; W17 stays 0.35.0).
- Schema id `bammetrics.tools/1`, file `docs/schema/tools-1.json`. Every object in it closes its members (`additionalProperties: false`), as the five existing schemas do.
- Plain `bammetrics --version` still prints `bammetrics X.Y.Z` and one line end.
- Plain `bammetrics tools` and `bammetrics help VERB` output stays byte for byte: existing golden records do not change. New golden records are added for the new forms.
- `make concordance` unchanged. No change to any tool's output.
- No new dependency. Nothing is written twice: tool options come from `Opt` rows and `Verb::flags`; options of the other verbs come from the clap command.
- Exit codes in the listing: 0, 1, 2, 10, 11, 100, 130.
- Commit messages carry no AI reference.

## Review Focus

1. A tool option with several names (`INPUT` and `I`, `METRIC_ACCUMULATION_LEVEL` and `LEVEL`): the first name is the Picard key, the rest are aliases, in `Opt::names` order. Test in Task 1.
2. Ignored and rejected options (`HISTOGRAM_FILE` of `CollectInsertSizeMetrics`, `ASSUME_SORTED` of `CollectWgsMetrics`): listed with their kind, so a wrapper generator can leave them out. Test in Task 1.
3. Help text with quotes, backslashes, tabs or non-ASCII bytes: the listing must still parse. Test in Task 1 parses the whole listing back through `bammetrics_formats::json`.
4. `--format` with a value other than `text` or `json` (`yaml`, `JSON`, empty): exit 2 with a usage line, nothing on standard output. Test in Task 2.
5. `--json` without `--version`, or `--version --json` with a further argument: exit 2. Plain `--version` unchanged. Test in Task 3.

## Designer's calls in this plan

- `default` is left out of an option entry. `Opt` has no default field; the default lives only in the help text ("(default 0.05)"). Parsing it out of help text would break on the first reworded line. A later wave may add `Opt::default`.
- The version object is the `bammetrics` member of the listing and the whole output of `--version --json`, one shape under one schema.
- `tools --format text` is accepted and equals plain `tools`.

---

### Task 1: The listing and its schema

**Files:**
- Create: `crates/cli/src/info/listing.rs` (move `crates/cli/src/info.rs` to `crates/cli/src/info/mod.rs` unchanged, then add `mod listing; pub use listing::{listing, tool_entry, version_object};`)
- Modify: `crates/cli/src/chain/prepare.rs:330-360`: split the role decision out of the step-role function into `pub(crate) fn role_of(verb: &Verb) -> Option<Role>` (variant role first, then transform, split, tap, materialize; `None` for a refused tool). The existing function calls it, so `chain` behaviour is unchanged.
- Create: `docs/schema/tools-1.json`
- Test: `crates/cli/tests/schemas.rs`, unit tests at the end of `listing.rs`

**Interfaces:**
- Produces: `pub fn listing(base: &clap::Command) -> Json`; `pub fn tool_entry(verb: &Verb) -> Json`; `pub fn version_object() -> Json`; `pub const EXIT_CODES: &[(i32, &str)]` in `listing.rs`.

Listing shape (members in this order):

```json
{
  "schema": "bammetrics.tools/1",
  "bammetrics": {"name": "bammetrics", "version": "0.34.1", "commit": "abc1234" | null,
                 "features": ["parquet", "pdf"], "schemas": ["bammetrics.metrics/1", "...", "bammetrics.tools/1"]},
  "exit_codes": [{"code": 0, "meaning": "..."}, ...],
  "tools": [{"verb": "insert-size", "class": "CollectInsertSizeMetrics", "summary": "...",
             "input": "BAM", "reads_stdin": true, "files": "single" | "prefix" | "bam",
             "writes": "...", "extensions": [".insert_size_metrics"],
             "chain_role": "transform" | "tap" | "split" | "materialize" | "none",
             "options": [{"names": ["HISTOGRAM_WIDTH", "W"], "flag": "--histogram-width", "short": "-W" | null,
                          "type": "bool" | "int" | "float" | "text" | "path",
                          "kind": "value" | "output" | "ignored" | "rejected",
                          "repeated": false, "required": false, "written": false, "help": "..."}]}],
  "verbs": [{"name": "qc", "about": "...", "options": [{"flag": "--preset", "short": null,
             "takes_value": true, "repeated": false, "required": false, "help": "..."}]}]
}
```

`reads_stdin` is `Verb::reads_stdin`: false means the tool reads its input twice and a wrapper must stage a file. `extensions` is `multiple_metrics::file_extensions(tool)`, empty when it returns none. `features` lists `parquet` and `pdf` by `cfg!(feature = ...)` in `bammetrics-cli`. `commit` comes from `option_env!("BAMMETRICS_COMMIT")`, set by a new `crates/cli/build.rs` from `git rev-parse --short HEAD` when that succeeds, else `null` (an sdist build has no git). `verbs` holds the `CLAP_VERBS` (`qc`, `report`, `check`, `plan`, `convert`, `doctor` and the rest), each with its options read from the clap subcommand of `info::command(base)`, then the `OTHER_VERBS` of `info/mod.rs` with their one line: `chain` with its options from `chain::args::global_flags()`, `collect` with its grammar flags, and `recipes`, `config`, `tools`, `presets`, `help`, `completions` with the options they take (none beyond `--format` for `tools` and `help`).

- [ ] **Step 1: Write the failing tests in `crates/cli/tests/schemas.rs`**

```rust
#[test]
fn tools_listing_is_valid_under_tools_1() {
    let out = bammetrics(&["tools", "--format", "json"]);
    assert_eq!(out.status.code(), Some(0));
    let doc = Json::parse(std::str::from_utf8(&out.stdout).unwrap()).unwrap();
    validate(&schema("tools-1.json"), &doc).unwrap();
}

#[test]
fn tools_listing_equals_the_registry() {
    // one entry per VERBS row, in VERBS order; each entry's option names equal
    // the tool's Opt rows (names in order, kind, type, repeated, required).
    // 77 tools at 0.34.1.
}

#[test]
fn tools_listing_keywords_equal_the_python_stub() {
    // for every tool, the snake_case of each value option's first name equals
    // registry::python_keywords(tool), the source of the generated .pyi block.
}

#[test]
fn listing_marks_aliases_ignored_and_rejected_options() {
    // CollectInsertSizeMetrics: names ["HISTOGRAM_FILE", "H"] has kind "ignored";
    // CollectWgsMetrics: ASSUME_SORTED has kind "rejected";
    // insert-size chain_role "tap", mark-duplicates "transform",
    // split-sam-by-library "split", build-bam-index "none", jumping-library reads_stdin false.
    // exit_codes codes == [0, 1, 2, 10, 11, 100, 130].
}
```

- [ ] **Step 2: Run them and see them fail**

Run: `cargo test -p bammetrics-cli --test schemas tools_listing`
Expected: FAIL, `tools --format json` exits 2 (`--format` unknown) or the schema file is missing.

- [ ] **Step 3: Write `docs/schema/tools-1.json`** in the style of `docs/schema/verdict-1.json`: title `bammetrics.tools/1`, `schema` constant, every object `additionalProperties: false`, enums for `files`, `chain_role`, `type`, `kind`.

- [ ] **Step 4: Implement `role_of`, `listing`, `tool_entry`, `version_object`, `EXIT_CODES`, `crates/cli/build.rs`.** The flag and short spelling of an option come from `Verb::flags`, matched to the `Opt` row by its Picard key. `EXIT_CODES` meanings copy the "Meaning" column of the docs/cli.md exit code table, plus `100: CheckTerminatorBlock found a DEFECTIVE file end` and `130: interrupted by SIGINT, SIGTERM or SIGHUP; partial outputs removed`.

- [ ] **Step 5: Run the tests**

Run: `cargo test -p bammetrics-cli --test schemas && cargo test -p bammetrics-cli --lib info && cargo test -p bammetrics-cli --test chain_positions`
Expected: PASS.

- [ ] **Step 6: Commit** `feat: tools listing as JSON under schema tools-1`

### Task 2: `tools --format` and `help VERB --format`

**Files:**
- Modify: `crates/cli/src/lib.rs:249-260` (the `tools` and `help` arms)
- Modify: `crates/cli/src/info/mod.rs` (`tools`, `help`, the `OTHER_VERBS` lines of `tools` and `help`)
- Test: `crates/cli/tests/golden.rs` (cases next to `verb_tools`, line 251), `tests/golden/*.txt`

**Interfaces:**
- Consumes: `listing`, `tool_entry` from Task 1.
- Produces: `pub fn info_format(verb: &str, args: &[String]) -> Result<(InfoFormat, Vec<String>), VerbError>` in `info/mod.rs`, `enum InfoFormat { Text, Json }`: pulls `--format VALUE` or `--format=VALUE` from any position and returns the remaining arguments.

- [ ] **Step 1: Add golden cases**

```rust
case!("verb_tools_json", "tools", "--format", "json"),
case!("verb_tools_text", "tools", "--format", "text"),
case!("verb_help_insert_size_json", "help", "insert-size", "--format", "json"),
case!("verb_help_class_json", "help", "--format=json", "CollectHsMetrics"),
case!("error_tools_format_yaml", "tools", "--format", "yaml"),
case!("error_help_format_upper", "help", "insert-size", "--format", "JSON"),
case!("error_help_qc_json", "help", "qc", "--format", "json"),
```

`help qc --format json` is a usage error in E1: a clap verb's entry lives in the `verbs` array of `tools --format json`, and the hint says so.

- [ ] **Step 2: Run** `cargo test --release -p bammetrics-cli --test golden verb_tools` and see the new cases fail (no record).

- [ ] **Step 3: Implement `info_format` and wire both arms.** `tools --format json` prints `listing(&Cli::command())` and a line end through `events::print_stdout`; `help TOPIC --format json` prints `tool_entry` of the verb `verbs::find` or `verbs::for_class` gives. A bad value: `VerbError::usage("tools", "--format takes text or json, not \"yaml\"")`, exit 2, standard output empty. The plain forms take the old code path unchanged.

- [ ] **Step 4: Write the records and read the diff**

Run: `make golden-update && git diff --stat tests/golden`
Expected: only the 7 new files added; `git diff tests/golden` shows no changed record. Then `make golden` passes.

- [ ] **Step 5: Commit** `feat: tools and help take --format json`

### Task 3: `--version --json`

**Files:**
- Modify: `crates/cli/src/lib.rs` (an arm before clap runs, next to the `tools` arm)
- Test: `crates/cli/tests/schemas.rs`, `crates/cli/tests/verbs.rs`

**Interfaces:**
- Consumes: `version_object` from Task 1.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn version_json_equals_the_listing_member() {
    // stdout of `--version --json` and of `--json --version` parse to the
    // `bammetrics` member of `tools --format json`; "version" is env!("CARGO_PKG_VERSION").
}

#[test]
fn plain_version_is_unchanged() {
    // `--version` prints exactly format!("bammetrics {}\n", env!("CARGO_PKG_VERSION")).
}

#[test]
fn json_flag_alone_or_with_more_is_a_usage_error() {
    // `--json`, `--version --json extra`: exit 2, stdout empty.
}
```

- [ ] **Step 2: Run** `cargo test -p bammetrics-cli --test schemas version` and see them fail.
- [ ] **Step 3: Implement the arm:** `args` exactly `["--version", "--json"]` in either order prints `version_object()`; `--json` in any other first-argument position exits 2.
- [ ] **Step 4: Run** `cargo test -p bammetrics-cli --test schemas --test verbs` and expect PASS.
- [ ] **Step 5: Commit** `feat: --version --json`

### Task 4: Docs and release 0.34.1

**Files:**
- Create: `docs/embedding.md`
- Modify: `docs/cli.md:224-251` (exit code table), `docs/cli.md` new section "Machine-readable listing" after "## Exit codes", `docs/formats.md` ("Frozen schemas" list gains `tools-1`), `CHANGELOG.md` (`## 0.34.1`), `README.md` (one line under the docs list), `CLAUDE.md` (one sentence after the W16 paragraphs), `Cargo.toml` and `pyproject.toml` version 0.34.1
- Regenerate: `tests/report/*.html` (they embed the version), `tests/golden/report_pdf*.txt` sizes, `bammetrics.pyi` if its header carries the version

- [ ] **Step 1: Add rows 100 and 130 to the docs/cli.md exit code table**, with the meanings `EXIT_CODES` holds. Then a unit test in `listing.rs`, `exit_codes_match_docs_cli`, reads `docs/cli.md` through `include_str!`, takes the rows between `## Exit codes` and the next `##`, and asserts their codes equal `EXIT_CODES` codes. Run it and see it pass.
- [ ] **Step 2: Write `docs/embedding.md`.** Sections: what a program reads (`tools --format json`, `help VERB --format json`, `--version --json`, the schema), the exit code table and what Nextflow and WDL do with each code (copy the spec's table in 5.1), the frozen schemas a program parses (`metrics-1`, `verdict-1`, `manifest-1`, `manifest-2`, `events-1`, `tools-1`), resume through manifest 2 (`docs/chaining.md` sections 27 to 30), the Python module (`docs/python.md`). One example: `bammetrics tools --format json | jq '.tools[] | select(.chain_role=="tap") | .verb'`. State that layers 2 to 4 (packages, wrappers, in-process library) come in E2 to E7.
- [ ] **Step 3: Bump to 0.34.1 and regenerate snapshots**

Run: `UPDATE_REPORT=1 cargo test -p bammetrics-report && make golden-update && UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs && UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs`
Expected: `git diff` outside the new files shows only version strings and the two PDF sizes.

- [ ] **Step 4: Full local gate**

Run: `cargo fmt --all --check && cargo clippy --workspace --all-targets -- -D warnings && cargo test --workspace && make golden && make concordance`
Expected: all pass; `make concordance` output equal to 0.34.0's.

- [ ] **Step 5: Commit** `docs: embedding guide, exit codes 100 and 130; release 0.34.1`

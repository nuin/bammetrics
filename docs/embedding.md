# Embedding BAMmetrics

A workflow engine or wrapper generator reads three things from the binary: a JSON listing of
every tool, an exit code table, and frozen JSON schemas for what the tools write. Layer 1
(this page) ships in 0.34.1. Layers 2 to 4 (packages, workflow wrappers, an in-process
library) come in E2 to E7.

## What a program reads

| Command | Prints |
|---|---|
| `bammetrics tools --format json` | The listing, schema `bammetrics.tools/1` ([tools-1.json](schema/tools-1.json)) |
| `bammetrics tools --format text` | Plain `tools`, a two-column text table |
| `bammetrics help VERB --format json` | One tool's entry from the listing |
| `bammetrics --version --json` | The `bammetrics` member: name, version, commit, features, schema names |

The listing has five members: `schema`, `bammetrics`, `exit_codes`, `tools` and `verbs`. Each
entry of `tools` (a registry tool) holds `verb`, `class`, `summary`, `input`, `reads_stdin`,
`files`, `writes`, `extensions`, `chain_role` and `options`. `chain_role` is `transform`,
`tap`, `split`, `materialize` or `none`. `verbs` lists the verbs with no registry row
(`chain`, `qc`, `report`, `check`, `doctor`, `convert`, `plan`, `prep`) with their options,
read from the clap command. `help qc --format json` exits 2 and points to that array.
Option defaults are in the help text only.

The features member tells a workflow whether the binary has `pdf` and `parquet`. Plain
`--version` stays `bammetrics X.Y.Z`.

Example, the verbs of every chain tap:

```bash
bammetrics tools --format json | jq '.tools[] | select(.chain_role=="tap") | .verb'
```

## Exit codes

The listing carries the table in `exit_codes`; [docs/cli.md](cli.md) holds the full text.

| Code | Meaning | Nextflow | WDL |
|---|---|---|---|
| 0 | Done | success | success |
| 1 | Input unreadable or tool failed | fail, retry only on out-of-memory codes the engine sets | fail |
| 2 | Usage error, nothing ran | fail, never retry | fail |
| 10 | `--fail-on` reached | wrappers accept it as a result unless a `fail_on` input is set | `returnCodes: [0, 10, 11]` unless gated |
| 11 | Verdict incomplete | same as 10 | same as 10 |
| 100 | `CheckTerminatorBlock` found a `DEFECTIVE` file end | fail | fail |
| 130 | Interrupted by SIGINT, SIGTERM or SIGHUP; partial outputs removed | fail, the engine killed the task | fail |

`qc` and `report` exit 0 on any verdict unless `--fail-on` is set, so a default wrapper never
fails a task on a bad sample. That suits a cohort run.

## Frozen schemas

A program may parse these documents. Each has a committed JSON Schema in `docs/schema`
([docs/formats.md](formats.md), "Frozen schemas").

| Schema | Written by |
|---|---|
| `bammetrics.metrics/1` ([metrics-1.json](schema/metrics-1.json)) | `--format json`, `--also json`, `convert` |
| `bammetrics.verdict/1` ([verdict-1.json](schema/verdict-1.json)) | `qc`, `report`, a chain's qc step |
| `bammetrics.manifest/1` ([manifest-1.json](schema/manifest-1.json)) | `qc` |
| `bammetrics.manifest/2` ([manifest-2.json](schema/manifest-2.json)) | `chain` that resumes, keeps a stream or ends in `qc` |
| events, version 1 ([events-1.json](schema/events-1.json)) | `qc --format jsonl`, `check --format jsonl` |
| `bammetrics.tools/1` ([tools-1.json](schema/tools-1.json)) | `tools --format json` |

A chain has no event stream. Per-step status is in manifest 2, one record per step.

## Resume

A chain that ends in `qc`, keeps a stream or runs with `--resume` writes `PREFIX.manifest.json`
(`bammetrics.manifest/2`). A rerun skips each step whose record matches and whose upstream
transforms were skipped. A task retried by the engine therefore redoes only what changed.
Details: [docs/chaining.md](chaining.md) sections 27 to 30.

## Python

`import bammetrics` gives the same tools as functions, `chain` and `step` for chains, and
`qc`, `check` and `report`. See [docs/python.md](python.md). `bammetrics.tools()` lists the
tools to Python callers only; a Nextflow or WDL author uses `tools --format json`.

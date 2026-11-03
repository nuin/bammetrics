# Verdicts, presets and verdict.json

Phase 3 of [docs/05_cli_and_ui_design.md](05_cli_and_ui_design.md) (sections 4.1 to 4.3).
The crate `crates/verdict` (`bammetrics-verdict`) evaluates a preset of threshold rules
over one sample's Picard metrics files and writes `verdict.json`. The `report` and
`check` verbs that call it arrive in phase 6; today the entry points are the Rust API
and `bammetrics.evaluate` in Python.

## Statuses

Each rule gives one check with one status:

| Status | When |
|---|---|
| `pass` | the value breaks neither bound |
| `warn` | the value breaks the warn bound and not the fail bound |
| `fail` | the value breaks the fail bound |
| `not_measured` | the metrics file, the row, the column or the value is absent (`?`, empty, text), or the ratio's divisor is 0. Never a pass. |

The band over all checks is `fail` when any check fails, else `incomplete` when any
check is not measured, else `warn` when any check warns, else `pass`. A warn or fail
check carries one reason line in plain words with the value and the limit:

```
Target bases at 20x is 0.91, below the fail limit of 0.95
Median insert size is 900, outside the fail range 120 to 600
Duplicate reads is 0.6 in LibC, above the fail limit of 0.5
```

A not-measured check names what is missing: `no HsMetrics file`, `HsMetrics has no
column X`, `no AlignmentSummaryMetrics row with CATEGORY=PAIR at ALL_READS`,
`FOLD_80_BASE_PENALTY is ?`, `PF_BASES is 0`.

## Preset file

A preset is a TOML file: one `[preset]` table and one `[[rule]]` per rule.

```toml
[preset]
name = "clinical-panel"
version = 1
note = "starting values"        # free text; the shipped presets say this

[[rule]]
id = "hs.target-20x"            # unique in the preset
label = "Target bases at 20x"   # plain words for the reason line; default: the id
class = "HsMetrics"             # metrics class, short name or full Java name
metric = "PCT_TARGET_BASES_20X" # a column, or a ratio: "PF_Q30_BASES / PF_BASES"
level = "ALL_READS"             # ALL_READS (default), SAMPLE, LIBRARY or READ_GROUP
select = "CATEGORY=PAIR"        # optional: only rows whose column prints as the value
warn = "< 0.98"                 # bound: < N, <= N, > N, >= N, outside A..B
fail = "< 0.95"                 # at least one of warn and fail; same comparator in both
plot = "coverage"               # optional plot family for the report (phase 4)
```

Rules the parser enforces, each with the file and line in the error: unknown keys and
unknown tables are errors; `fail` must be at least as strict as `warn` (`< 0.95` beside
`< 0.98`, `outside 120..600` beside `outside 180..450`); a bound breaks when the value
is on the named side of the threshold (`< 0.98` breaks at 0.97, holds at 0.98).

Row scope. `level` selects rows by their `SAMPLE`, `LIBRARY` and `READ_GROUP` cells:
`ALL_READS` rows leave all three empty, a `SAMPLE` row fills only `SAMPLE`, a `LIBRARY`
row fills `LIBRARY` and not `READ_GROUP`, a `READ_GROUP` row fills `READ_GROUP`. A
column the class lacks counts as empty. `CollectDuplicateMetrics` writes one row per
`LIBRARY` and no `SAMPLE` column, so a duplication rule says `level = "LIBRARY"`. When
several rows match (two libraries), every row is judged and the worst is reported with
its unit in `unit` and in the reason.

The file reader accepts the subset of TOML the presets use: comments, `[table]`,
`[[table]]`, and `key = value` with a double-quoted string, an integer, a decimal, `true`
or `false`. A preset file is valid TOML, so any TOML tool reads it. The parser is
`table::parse_tables` in the crate, the one function a switch to the `toml` crate (open
decision 4 of docs/05) would replace.

The preset hash is FNV-1a 64 of the file text with every `\r` removed, printed as 16 hex
digits. The same file hashes the same on every platform and under either line ending;
a change to any byte, a comment included, changes it. The report prints the first 6
digits; `verdict.json` holds all 16.

A lab owns its file: `bammetrics.evaluate("lab.toml", paths)` today, `--preset lab.toml`
on the verbs of phase 6. The shipped presets are compiled into the binary.

## Shipped presets

All four carry `note = "starting values"`. The numbers are not validated limits for any
assay; they need maintainer and lab sign-off before anyone calls them clinical, and the
verdict is not a release decision.

| Preset | Rules | Reads | Band on a bammetrics run today |
|---|---|---|---|
| `generic` | 6 | alignment summary (3), insert size, quality yield, duplicate metrics | complete |
| `clinical-panel` | 14 | the table in docs/05 section 4.3: hs metrics (9), alignment summary (3), insert size, quality yield | complete with `CollectHsMetrics` |
| `wgs` | 10 | `WgsMetrics` (4), alignment summary (3), insert size, quality yield, duplicate metrics | INCOMPLETE until `CollectWgsMetrics` exists; a Picard `WgsMetrics` file completes it |
| `rna-seq` | 8 | `RnaSeqMetrics` (4), alignment summary (2), quality yield, duplicate metrics | INCOMPLETE until `CollectRnaSeqMetrics` exists; a Picard `RnaSeqMetrics` file completes it |

The alignment rules select `CATEGORY=PAIR`; a single-end BAM has no such row, so they
read not measured and the band is INCOMPLETE until the lab's preset selects `UNPAIRED`.

## verdict.json

Schema `bammetrics.verdict/1`. The report (phase 4) and the hosted dashboard read this
file and nothing else from a run. Written by `Verdict::to_json`, read back by
`Verdict::from_json`; `bammetrics.evaluate` returns the same structure as dicts.

```json
{
  "schema": "bammetrics.verdict/1",
  "tool_version": "0.6.0",
  "sample": "NA12878",
  "preset": {"name": "clinical-panel", "version": 1, "hash": "3d0bcdc8430441db",
             "source": "built-in", "note": "starting values"},
  "status": "fail",
  "counts": {"pass": 8, "warn": 0, "fail": 5, "not_measured": 1},
  "files": [{"path": "s.hs_metrics", "classes": ["picard.analysis.directed.HsMetrics"]}],
  "checks": [
    {"rule": "hs.target-20x", "label": "Target bases at 20x", "class": "HsMetrics",
     "metric": "PCT_TARGET_BASES_20X", "level": "ALL_READS", "select": null, "unit": null,
     "value": 0.91, "warn": "< 0.98", "fail": "< 0.95", "status": "fail",
     "reason": "Target bases at 20x is 0.91, below the fail limit of 0.95",
     "plot": "coverage"}
  ]
}
```

| Field | Content |
|---|---|
| `schema` | `bammetrics.verdict/1`; a reader rejects any other value |
| `tool_version` | the bammetrics version that wrote the file |
| `sample` | as given, else the first non-empty `SAMPLE` cell in any class block, else the first file's name up to its first `.` |
| `preset` | `name`, `version`, `hash` (16 hex digits), `source` (`built-in` or the file path), `note` |
| `status` | `pass`, `warn`, `fail` or `incomplete` |
| `counts` | checks per status: `pass`, `warn`, `fail`, `not_measured` |
| `files` | each metrics file read, with its path and the metrics classes it holds, in the order given |
| `checks` | one per rule, in preset order: `rule`, `label`, `class`, `metric` (`A / B` for a ratio), `level`, `select` (`COLUMN=VALUE` or null), `unit` (the row's sample, library or read group, or null at `ALL_READS`), `value` (a number, or null when not measured), `warn` and `fail` (the bounds as written, or null), `status`, `reason` (null on a pass), `plot` (or null) |

Numbers are JSON numbers; a value that is NaN or infinite is written as null, and the
check is not measured. The input fingerprint (BAM size, mtime and head/tail hash) joins
in phase 7 with the manifest; `files` is the bridge until then.

## Rust API

```rust
use bammetrics_verdict::{Preset, Input, evaluate, evaluate_files};

let preset = Preset::resolve("clinical-panel")?;        // a shipped name, or a file path
let verdict = evaluate_files(&preset, &[hs, alignment, insert, yield_], Some("NA12878"))?;
std::fs::write("NA12878.verdict.json", verdict.to_json())?;
```

`evaluate(&preset, &[Input { path, doc }], sample)` takes parsed `MetricsDoc`s, which is
what `qc` will hand it from memory in phase 7. `bammetrics_verdict::names()` lists the
shipped presets; `Preset::builtin(name)`, `Preset::load(path)` and `Preset::parse(text,
source)` are the other constructors.

Rules for this crate, from docs/05 section 4.9: it reads the metrics files it is given
and nothing else. No other run directory, no history store, no user identity, no network
call, no comment field. A feature that needs one of those belongs to the hosted dashboard.

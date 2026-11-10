# Verdicts, presets and verdict.json

Phase 3 of [docs/05_cli_and_ui_design.md](05_cli_and_ui_design.md) (sections 4.1 to 4.3).
The crate `crates/verdict` (`bammetrics-verdict`) evaluates a preset of threshold rules
over one sample's Picard metrics files and writes `verdict.json`. Entry points: `qc`,
`report` and `check` with `--preset` ([docs/cli.md](cli.md)), `bammetrics presets`, the
Rust API below, and `bammetrics.evaluate`, `qc`, `check` and `report` in Python
([docs/python.md](python.md)).

## Statuses

Each rule gives one check with one status:

| Status | When |
|---|---|
| `pass` | the value breaks neither bound |
| `warn` | the value breaks the warn bound and not the fail bound |
| `fail` | the value breaks the fail bound |
| `not_measured` | the metrics file, the row, the column or the value is absent (`?`, empty, text under a numeric bound), or the ratio's divisor is 0. Never a pass. |

A rule with an `is TEXT` bound judges the cell's text instead of a number: it breaks
when the cell prints as `TEXT`, an empty cell is the empty string, and any other text
passes (docs/06 section 2.5). The match is exact and case sensitive. The shipped text
rules read `SexCheckMetrics.SEX_MATCH` (`Y`, `N` or empty) and
`SexCheckMetrics.INFERRED_SEX` (`FEMALE`, `MALE` or `UNKNOWN`); an `is` bound on a column
the glossary types as a number is a preset error. `check` prints the cell in its value
column, `empty` for an empty cell.

The band over all checks is `fail` when any check fails, else `incomplete` when any
check is not measured, else `warn` when any check warns, else `pass`. A warn or fail
check carries one reason line in plain words with the value and the limit:

```
Target bases at 20x is 0.91, below the fail limit of 0.95
Median insert size is 900, outside the fail range 120 to 600
Duplicate reads is 0.6 in LibC, above the fail limit of 0.5
Declared sex matches: N, a fail value
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
note = "starting values"        # free text; shown on the page and in verdict.json

[[rule]]
id = "hs.target-20x"            # unique in the preset
label = "Target bases at 20x"   # plain words for the reason line; default: the id
class = "HsMetrics"             # metrics class, short name or full Java name
metric = "PCT_TARGET_BASES_20X" # a column, or a ratio: "PF_Q30_BASES / PF_BASES"
level = "ALL_READS"             # ALL_READS (default), SAMPLE, LIBRARY or READ_GROUP
select = "CATEGORY=PAIR"        # optional: only rows whose column prints as the value
warn = "< 0.98"                 # bound: < N, <= N, > N, >= N, outside A..B, is TEXT
fail = "< 0.95"                 # at least one of warn and fail; same comparator in both
plot = "coverage"               # optional plot family for the report
```

Rules the parser enforces, each with the file and line in the error: unknown keys and
unknown tables are errors; `fail` must be at least as strict as `warn` (`< 0.95` beside
`< 0.98`, `outside 120..600` beside `outside 180..450`; two `is` bounds need the same
text); a bound breaks when the value is on the named side of the threshold (`< 0.98`
breaks at 0.97, holds at 0.98); an `is` bound takes one column, not a ratio, and not a
column the glossary types as `int` or `float`. A rule with an unknown key is reported by
that key before a missing bound, so `expect = "Y"` says `unknown key expect`.

Row scope. `level` selects rows by their `SAMPLE`, `LIBRARY` and `READ_GROUP` cells:
`ALL_READS` rows leave all three empty, a `SAMPLE` row fills only `SAMPLE`, a `LIBRARY`
row fills `LIBRARY` and not `READ_GROUP`, a `READ_GROUP` row fills `READ_GROUP`. A
column the class lacks counts as empty. `CollectDuplicateMetrics` writes one row per
`LIBRARY` and no `SAMPLE` column, so a duplication rule says `level = "LIBRARY"`. When
several rows match (two libraries), every row is judged and the worst is reported with
its unit in `unit` and in the reason.

Plot families: `coverage`, `insert-size`, `gc-bias` and `cycle-quality` (docs/05
section 4.4). The page draws a family under the first warn or fail check whose rule
names it, once per page.

The file reader accepts the subset of TOML the presets use: comments, `[table]`,
`[[table]]`, and `key = value` with a double-quoted string, an integer, a decimal, `true`
or `false`. A preset file is valid TOML, so any TOML tool reads it. The parser is
`table::parse_tables` in the crate, the one function a switch to the `toml` crate (open
decision 4 of docs/05) would replace.

The preset hash is FNV-1a 64 of the file text with every `\r` removed, printed as 16 hex
digits. The same file hashes the same on every platform and under either line ending;
a change to any byte, a comment included, changes it. The page, the index, the `qc`
summary and `bammetrics presets` print the first 6 digits; `verdict.json` and the
manifest hold all 16.

A lab owns its file: `--preset lab.toml` on `qc`, `report` and `check`, or
`bammetrics.evaluate("lab.toml", paths)`. The shipped presets are compiled into the
binary; to change one bound, copy the shipped file and edit it (the `--set` override of
docs/05 section 4.3 is not built).

## Shipped presets

`clinical-panel` version 3 carries `note = "reviewed defaults; gap and sex rules are
starting values"`: the maintainer reviewed its first 14 thresholds on 2026-11-06 and kept
the version 1 values (decision 2 of docs/05 section 8); the four rules of version 3
(`gaps.reportable-range`, `gaps.targets-complete`, `sex.match`, `sex.inferred`,
[docs/06](06_clinical_checks.md)) are starting values. The sex rules live in
`clinical-panel` only: `generic` and `wgs` stay at version 1, so a `report` or `check` over
Picard output alone still completes under them. `generic`, `wgs` and `rna-seq` carry
`note = "starting values"`.
No shipped number is a validated limit for any assay: a lab validates and owns its own
copy, and the verdict is not a release decision.

| Preset | Rules | Reads | What completes it |
|---|---|---|---|
| `generic` | 6 | alignment summary (3), insert size, quality yield, duplicate metrics | `qc BAM` |
| `clinical-panel` | 18 | hs metrics (9), alignment summary (3), insert size, quality yield, coverage gaps (2), sex check (2) | `qc BAM --panel FILE` (or `--hs BAITS TARGETS`); `-R` for `GC_DROPOUT`; `--sex male|female` for `sex.match` to judge a declaration |
| `wgs` | 10 | `WgsMetrics` (4), alignment summary (3), insert size, quality yield, duplicate metrics | `qc BAM -R FASTA`; `CollectWgsMetrics` reads the BAM on its own |
| `rna-seq` | 8 | `RnaSeqMetrics` (4), alignment summary (2), quality yield, duplicate metrics | `qc BAM --ref-flat FILE --strand S`; `PCT_RIBOSOMAL_BASES` prints empty without `--ribosomal-intervals`, so that rule reads not measured |

The alignment rules select `CATEGORY=PAIR`; a single-end BAM has no such row, so they
read not measured and the band is INCOMPLETE until the lab's preset selects `UNPAIRED`.
The insert size rules select `PAIR_ORIENTATION=FR`. `sex.match` (`fail = "is N"`) passes
on `Y` and on the empty cell a run without `--sex` writes; `sex.inferred`
(`warn = "is UNKNOWN"`) warns when the reads give no sex.

The four tables below are rendered from the compiled presets by
`crates/verdict/tests/doc_presets.rs`, which fails when a table differs from
`presets/*.toml`; `UPDATE_DOCS=1 cargo test -p bammetrics-verdict --test doc_presets`
rewrites them. Bounds print as the parser reads them (`< 0.80` prints `< 0.8`), the
form `verdict.json` and the page use. Rows is empty for `ALL_READS` with no selector.

### generic

<!-- preset generic -->
Version 1, 6 rules, note `starting values`, hash `6dd89ed0478fd5e4` (the report prints `6dd89e`).

| Rule | Label | Metric | Rows | Warn | Fail | Plot |
|---|---|---|---|---|---|---|
| `alignment.pf-reads-aligned` | PF reads aligned | `AlignmentSummaryMetrics.PCT_PF_READS_ALIGNED` | `CATEGORY=PAIR` | `< 0.9` | `< 0.8` |  |
| `alignment.adapter` | Adapter reads | `AlignmentSummaryMetrics.PCT_ADAPTER` | `CATEGORY=PAIR` | `> 0.01` | `> 0.05` | `insert-size` |
| `alignment.chimeras` | Chimeric pairs | `AlignmentSummaryMetrics.PCT_CHIMERAS` | `CATEGORY=PAIR` | `> 0.05` | `> 0.1` |  |
| `insert.median` | Median insert size | `InsertSizeMetrics.MEDIAN_INSERT_SIZE` | `PAIR_ORIENTATION=FR` | `outside 100..1000` | `outside 50..2000` | `insert-size` |
| `yield.q30` | Q30 fraction of PF bases | `QualityYieldMetrics.PF_Q30_BASES / PF_BASES` |  | `< 0.8` | `< 0.7` | `cycle-quality` |
| `duplication.percent` | Duplicate reads | `DuplicationMetrics.PERCENT_DUPLICATION` | `LIBRARY` | `> 0.3` | `> 0.5` |  |
<!-- end preset generic -->

### clinical-panel

The table of docs/05 section 4.3, with the rule ids and selectors.

<!-- preset clinical-panel -->
Version 3, 18 rules, note `reviewed defaults; gap and sex rules are starting values`, hash `6f0d59319750923f` (the report prints `6f0d59`).

| Rule | Label | Metric | Rows | Warn | Fail | Plot |
|---|---|---|---|---|---|---|
| `hs.mean-target-coverage` | Mean target coverage | `HsMetrics.MEAN_TARGET_COVERAGE` |  | `< 200` | `< 100` | `coverage` |
| `hs.target-20x` | Target bases at 20x | `HsMetrics.PCT_TARGET_BASES_20X` |  | `< 0.98` | `< 0.95` | `coverage` |
| `hs.target-100x` | Target bases at 100x | `HsMetrics.PCT_TARGET_BASES_100X` |  | `< 0.8` | `< 0.6` | `coverage` |
| `hs.evenness` | Fold-80 base penalty | `HsMetrics.FOLD_80_BASE_PENALTY` |  | `> 2` | `> 3` | `coverage` |
| `hs.usable-on-target` | Usable bases on target | `HsMetrics.PCT_USABLE_BASES_ON_TARGET` |  | `< 0.5` | `< 0.3` |  |
| `hs.selected-bases` | Selected bases | `HsMetrics.PCT_SELECTED_BASES` |  | `< 0.7` | `< 0.5` |  |
| `hs.excluded-duplicates` | Bases excluded as duplicates | `HsMetrics.PCT_EXC_DUPE` |  | `> 0.3` | `> 0.5` |  |
| `hs.zero-coverage-targets` | Targets with zero coverage | `HsMetrics.ZERO_CVG_TARGETS_PCT` |  | `> 0.01` | `> 0.03` | `coverage` |
| `hs.gc-dropout` | GC dropout | `HsMetrics.GC_DROPOUT` |  | `> 5` | `> 10` | `gc-bias` |
| `alignment.pf-reads-aligned` | PF reads aligned | `AlignmentSummaryMetrics.PCT_PF_READS_ALIGNED` | `CATEGORY=PAIR` | `< 0.98` | `< 0.95` |  |
| `alignment.adapter` | Adapter reads | `AlignmentSummaryMetrics.PCT_ADAPTER` | `CATEGORY=PAIR` | `> 0.005` | `> 0.02` | `insert-size` |
| `alignment.chimeras` | Chimeric pairs | `AlignmentSummaryMetrics.PCT_CHIMERAS` | `CATEGORY=PAIR` | `> 0.05` | `> 0.1` |  |
| `insert.median` | Median insert size | `InsertSizeMetrics.MEDIAN_INSERT_SIZE` | `PAIR_ORIENTATION=FR` | `outside 180..450` | `outside 120..600` | `insert-size` |
| `yield.q30` | Q30 fraction of PF bases | `QualityYieldMetrics.PF_Q30_BASES / PF_BASES` |  | `< 0.85` | `< 0.75` | `cycle-quality` |
| `gaps.reportable-range` | Target bases at MIN_DEPTH | `CoverageGapMetrics.PCT_REPORTABLE` |  | `< 0.98` | `< 0.95` | `coverage` |
| `gaps.targets-complete` | Targets without a gap | `CoverageGapMetrics.PCT_TARGETS_COMPLETE` |  | `< 0.95` | `< 0.9` | `coverage` |
| `sex.match` | Declared sex matches | `SexCheckMetrics.SEX_MATCH` |  |  | `is N` |  |
| `sex.inferred` | Inferred sex | `SexCheckMetrics.INFERRED_SEX` |  | `is UNKNOWN` |  |  |
<!-- end preset clinical-panel -->

### wgs

<!-- preset wgs -->
Version 1, 10 rules, note `starting values`, hash `8a5fc4360f357d10` (the report prints `8a5fc4`).

| Rule | Label | Metric | Rows | Warn | Fail | Plot |
|---|---|---|---|---|---|---|
| `wgs.mean-coverage` | Mean coverage | `WgsMetrics.MEAN_COVERAGE` |  | `< 30` | `< 20` | `coverage` |
| `wgs.bases-10x` | Genome bases at 10x | `WgsMetrics.PCT_10X` |  | `< 0.95` | `< 0.9` | `coverage` |
| `wgs.bases-30x` | Genome bases at 30x | `WgsMetrics.PCT_30X` |  | `< 0.8` | `< 0.6` | `coverage` |
| `wgs.evenness` | Fold-80 base penalty | `WgsMetrics.FOLD_80_BASE_PENALTY` |  | `> 1.5` | `> 2` | `coverage` |
| `alignment.pf-reads-aligned` | PF reads aligned | `AlignmentSummaryMetrics.PCT_PF_READS_ALIGNED` | `CATEGORY=PAIR` | `< 0.98` | `< 0.95` |  |
| `alignment.adapter` | Adapter reads | `AlignmentSummaryMetrics.PCT_ADAPTER` | `CATEGORY=PAIR` | `> 0.005` | `> 0.02` | `insert-size` |
| `alignment.chimeras` | Chimeric pairs | `AlignmentSummaryMetrics.PCT_CHIMERAS` | `CATEGORY=PAIR` | `> 0.02` | `> 0.05` |  |
| `insert.median` | Median insert size | `InsertSizeMetrics.MEDIAN_INSERT_SIZE` | `PAIR_ORIENTATION=FR` | `outside 250..600` | `outside 150..900` | `insert-size` |
| `yield.q30` | Q30 fraction of PF bases | `QualityYieldMetrics.PF_Q30_BASES / PF_BASES` |  | `< 0.85` | `< 0.75` | `cycle-quality` |
| `duplication.percent` | Duplicate reads | `DuplicationMetrics.PERCENT_DUPLICATION` | `LIBRARY` | `> 0.15` | `> 0.3` |  |
<!-- end preset wgs -->

### rna-seq

<!-- preset rna-seq -->
Version 1, 8 rules, note `starting values`, hash `2d1283287b120a2e` (the report prints `2d1283`).

| Rule | Label | Metric | Rows | Warn | Fail | Plot |
|---|---|---|---|---|---|---|
| `rna.ribosomal` | Ribosomal bases | `RnaSeqMetrics.PCT_RIBOSOMAL_BASES` |  | `> 0.1` | `> 0.3` |  |
| `rna.mrna-bases` | mRNA bases | `RnaSeqMetrics.PCT_MRNA_BASES` |  | `< 0.6` | `< 0.4` |  |
| `rna.intergenic` | Intergenic bases | `RnaSeqMetrics.PCT_INTERGENIC_BASES` |  | `> 0.15` | `> 0.3` |  |
| `rna.five-to-three-bias` | Median 5' to 3' bias | `RnaSeqMetrics.MEDIAN_5PRIME_TO_3PRIME_BIAS` |  | `outside 0.5..2` | `outside 0.25..4` |  |
| `alignment.pf-reads-aligned` | PF reads aligned | `AlignmentSummaryMetrics.PCT_PF_READS_ALIGNED` | `CATEGORY=PAIR` | `< 0.9` | `< 0.8` |  |
| `alignment.adapter` | Adapter reads | `AlignmentSummaryMetrics.PCT_ADAPTER` | `CATEGORY=PAIR` | `> 0.01` | `> 0.05` | `insert-size` |
| `yield.q30` | Q30 fraction of PF bases | `QualityYieldMetrics.PF_Q30_BASES / PF_BASES` |  | `< 0.85` | `< 0.75` | `cycle-quality` |
| `duplication.percent` | Duplicate reads | `DuplicationMetrics.PERCENT_DUPLICATION` | `LIBRARY` | `> 0.5` | `> 0.7` |  |
<!-- end preset rna-seq -->

## verdict.json

Schema `bammetrics.verdict/1`. The report and the hosted dashboard read this file and
nothing else from a run. Written by `Verdict::to_json`, read back by
`Verdict::from_json` (`check` reads a `verdict.json` given to it this way);
`bammetrics.evaluate` returns the same structure as dicts.

```json
{
  "schema": "bammetrics.verdict/1",
  "tool_version": "0.16.0",
  "sample": "NA12878",
  "preset": {"name": "clinical-panel", "version": 3, "hash": "6f0d59319750923f",
             "source": "built-in", "note": "reviewed defaults; gap and sex rules are starting values"},
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
| `checks` | one per rule, in preset order: `rule`, `label`, `class`, `metric` (`A / B` for a ratio), `level`, `select` (`COLUMN=VALUE` or null), `unit` (the row's sample, library or read group, or null at `ALL_READS`), `value` (a number, or null when not measured or judged as text), `text` (the cell's text, present only on a check with an `is TEXT` bound), `warn` and `fail` (the bounds as parsed, or null), `status`, `reason` (null on a pass), `plot` (or null) |

Numbers are JSON numbers; a value that is NaN or infinite is written as null, and the
check is not measured. The input fingerprint (BAM size, mtime and the hash of its first
and last 64 KiB) is in `PREFIX.manifest.json`, which `qc` writes beside the verdict;
`verdict.json` lists the metrics files in `files`.

## Rust API

```rust
use bammetrics_verdict::{Preset, Input, evaluate, evaluate_files};

let preset = Preset::resolve("clinical-panel")?;        // a shipped name, or a file path
let verdict = evaluate_files(&preset, &[hs, alignment, insert, yield_], Some("NA12878"))?;
std::fs::write("NA12878.verdict.json", verdict.to_json())?;
```

`evaluate(&preset, &[Input { path, doc }], sample)` takes parsed `MetricsDoc`s, which is
what `qc` hands it from memory. `bammetrics_verdict::names()` lists the shipped presets;
`Preset::builtin(name)`, `Preset::load(path)` and `Preset::parse(text, source)` are the
other constructors.

Rules for this crate, from docs/05 section 4.9: it reads the metrics files it is given
and nothing else. No other run directory, no history store, no user identity, no network
call, no comment field. A feature that needs one of those belongs to the hosted dashboard.

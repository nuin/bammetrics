# CLI and report design (docs/05_cli_and_ui_design.md)

Status: proposal, read against origin/main 87e8d97. Items marked NEW do not exist yet. Facts about the current CLI come from a surface audit, not from a fresh read of every source file.

## 1. Goals and non-goals

Goal: one verb-first CLI that covers every tool in the registry, chains steps by pipes and a shared output prefix, and ends in a verdict-first static HTML report. The Picard syntax stays byte-identical. The report is one self-contained file per sample plus one run index. It is written by the binary in under 20 ms (a target, to be measured), with no server, JVM, Python or JavaScript. It opens with pass, warn or fail per metric against a preset, plots only what explains a failure, and keeps every Picard column behind a collapsed table with its definition. It rejects the two models the maintainer named. MultiQC draws a wall of generic plots, gives no verdict and ships a JavaScript runtime. Qualimap writes dozens of plots, gives no verdict, needs a JVM and looks dated. Non-goals: a workflow engine (Nextflow and Snakemake own that; `chain`, docs/chaining.md, is a line of steps with taps, with no branching, loops, scatter or cluster), run history, trends, comments, users, audit trails and uploads (the paid dashboard in ROADMAP.md owns those), and any change to Picard-face output. (Until docs/07 the list also held "new metric tools beyond the bed-to-interval-list step"; D7 of docs/07 section 7 removed it.)

## 2. CLI grammar

### 2.1 Shape

```
bammetrics <verb> [INPUT] [options]       NEW kebab-case verbs
bammetrics <PicardClass> KEY=VALUE ...    unchanged drop-in
bammetrics <old-subcommand> ...           unchanged (5 subcommands)
```

Dispatch order in main.rs: (1) first argument is a Picard class in the registry or in picard_cli: legacy path, untouched. (2) first argument is a NEW verb (`chain` among them since W04 of docs/07, and `plan` when its arguments describe a chain). (3) the 5 clap subcommands. Class names are CamelCase and verbs are kebab, so they cannot collide. The old subcommands `insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics` and `collect` become the same entries in the new verb table, with every flag and default kept.

Verbs:

| Group | Verbs |
|---|---|
| Run | `qc` (preset or metric set, one pass, files, verdict, report), `collect` (same engine, no verdict, kept), `chain` (steps split at `+` in one process, records handed on in memory; W04 of docs/07, docs/chaining.md), one verb per tool |
| Prepare | `bed-to-interval-list` (P4 library piece pulled forward; 0.14.0), `mark-duplicates` (0.18.0), `fifo-buffer` (docs/07 wave W04: standard input to standard output, no input argument and no `-o`), `build-bam-index` (docs/07 wave W05: `-o` names the `.bai`, default beside the BAM, no metrics file and no `--format`) |
| Post | `accumulate-quality-yield`, `compare-metrics`, `report`, `check` |
| Info | `tools`, `presets`, `plan`, `completions <shell>`, `help <verb or Picard class>` |

Tool verbs. The name is the Picard class without `Collect` and `Metrics`, in kebab case.

| Verb | Picard class |
|---|---|
| insert-size | CollectInsertSizeMetrics |
| alignment-summary | CollectAlignmentSummaryMetrics |
| quality-yield | CollectQualityYieldMetrics |
| hs-metrics | CollectHsMetrics |
| gc-bias | CollectGcBiasMetrics (needs `-R`) |
| quality-score-distribution | QualityScoreDistribution |
| mean-quality-by-cycle | MeanQualityByCycle |
| base-distribution-by-cycle | CollectBaseDistributionByCycle |
| duplicate-metrics | CollectDuplicateMetrics |
| library-complexity | EstimateLibraryComplexity |
| multiple-metrics | CollectMultipleMetrics (expands to a metric set) |
| accumulate-quality-yield | AccumulateQualityYieldMetrics |
| compare-metrics | CompareMetrics (the rule would give `compare`; the verb keeps the class's noun) |

### 2.2 Shared options

| Option | Meaning | Picard spelling |
|---|---|---|
| INPUT or `-i` | BAM path, or `-` for stdin | I |
| `-o` | file, prefix, or `-` for stdout (`-O` kept) | O |
| `-@`, `--threads` | default min(cores, 8), help text states it | none |
| `-R` | FASTA with .fai | R, REFERENCE_SEQUENCE |
| `--sample NAME` | sample name, default first @RG SM, else file stem | none |
| `--level` | accumulation level, repeatable | METRIC_ACCUMULATION_LEVEL |
| `--stop-after N` | read limit | STOP_AFTER |
| `--format` | `picard`, `json`, `tsv`, `jsonl`, `human` | none |
| `-q`, `-v` | quiet, verbose | QUIET, VERBOSITY |

Spelling rules:

- One option has one table row in registry.rs. The row generates the kebab flag, the Picard key, `--help`, the Python kwarg and the `.pyi` stub (0.13.0; they were three hand-kept lists). This is the StepSpec idea from the Run-first design.
- `--minimum-pct`, `--MINIMUM_PCT` and `MINIMUM_PCT=` are one option. The normaliser upper-cases and swaps `-` for `_`, then calls the existing registry parser.
- Booleans take `--x`, `--x=false` and `--no-x`. The old space form `--clip-overlapping-reads false` still parses. A BAM path is never read as a boolean value.
- An unknown option exits 2 with a did-you-mean from the option table.
- Errors use one style on new verbs: `bammetrics <verb>: <what> (from <option>)`, then at most 2 hint lines.

### 2.3 Output naming

- Single-file verbs: `-o FILE`. No `-o` prints Picard text to stdout, as today.
- Multi-file verbs (`qc`, `collect`, `multiple-metrics`, `gc-bias`): `-o PREFIX`. No `-o` means `./<bam stem>`, and the summary says so. Today `collect` errors here.
- Extensions come from one table in registry.rs. They reuse the CollectMultipleMetrics extensions, so both faces write the same names: `.alignment_summary_metrics`, `.insert_size_metrics`, `.quality_yield_metrics`, `.hs_metrics`.
- Several Hs sets: `PREFIX.<name>.hs_metrics`. The name is the bait file stem (`qc --hs BAITS TARGETS`, as `collect --hs` names it); two bait files with one stem are an error. A `--hs name=baits,targets` form is not built (0.12.0). A test pins the name table so it cannot drift from the Picard face.
- Run files: `PREFIX.manifest.json`, `PREFIX.verdict.json`, `PREFIX.report.html`.
- Missing output directories are created. A `qc -o` that ends in `/` or names an existing directory gets the BAM stem as the file part. Files are written to a temporary name and renamed. The summary counts overwrites (`6 replaced`); it does not list them.
- The sample name is not written into Picard files, which stay byte-compatible. It goes in the manifest.

### 2.4 Exit codes

The table `check.rs` implements (`EXIT_CODES`), shared by `check` and `qc`:

| Code | Meaning |
|---|---|
| 0 | success: every verdict passes, or none reaches `--fail-on` (`qc` without `--fail-on`: any verdict) |
| 1 | an input could not be read (a missing BAM, FASTA, .fai or interval file), the preset does not parse, or a `qc` collector failed |
| 2 | usage error, nothing ran (bad arguments; `qc` over an earlier run without `--resume` or `--force`) |
| 10 | a verdict reached `--fail-on`: fail, or warn with `--fail-on warn` (`check` defaults to `--fail-on fail`) |
| 11 | a verdict is incomplete (a rule's file, row or column is missing) and none reached `--fail-on` |

`qc` exits 0 on a failed verdict unless `--fail-on warn|fail` is set. A Nextflow task is not killed by a bad sample. `check` always gates. The summary prints `use --fail-on fail or bammetrics check to gate`. A `qc` run with a failed collector exits 1 after it writes the other files, the verdict and the page. The Picard face keeps 0 and 1. Both faces use the same defaults for the same verdict: `qc` never gates, `check` always gates.

### 2.5 Twelve worked examples

Output is illustrative.

**1. One tool, as today.**
```
bammetrics insert-size sample.bam
```
Picard text on stdout. Byte-identical to the current subcommand.

**2. One tool, machine-readable.**
```
bammetrics insert-size sample.bam --format json | jq '.metrics[0].MEDIAN_INSERT_SIZE'
310
bammetrics insert-size sample.bam --MINIMUM_PCT 0.1 --format tsv
```
Kebab and Picard spellings mix freely.

**3. Quick look, no preset.**
```
$ bammetrics qc sample.bam -R ref.fa
S1  WARN  preset generic v1 (6dd89e)  4 collectors, 1 BAM read (41.8M records)
WARN  alignment.PCT_ADAPTER   0.012   above 0.005   adapter in 1.2% of reads
NOT MEASURED  gc-bias.AT_DROPOUT  needs -R
PASS  4 checks
wrote 4 files, bammetrics-out/S1.*
```
With no `--preset`, `qc` uses the named preset `generic` and prints its name. It never picks a clinical preset by itself. A metric it cannot compute is listed as not measured. It is not skipped silently.

**4. Clinical panel, one sample.**
```
$ bammetrics qc sample.bam --preset clinical-panel --panel panel.bed -R ref.fa -o out/NA12878
NA12878  FAIL  preset clinical-panel v2 (e62755)  4 collectors, 1 BAM read (41.8M records)
FAIL  hs.PCT_TARGET_BASES_20X  0.874  below 0.95   12.6% of panel bases are under 20x
WARN  insert.MEDIAN_INSERT_SIZE  171  below 180   fragments shorter than usual
PASS  12 checks
wrote 7 files, report out/NA12878.report.html
```
`--panel` is baits and targets. The collectors follow the preset's rules: `clinical-panel` reads quality yield, alignment summary, insert size and hs, so those 4 run in 1 read. `--panel` adds hs under a preset with no hs rule too. A rule set that reads gc bias or duplication adds those collectors. Hs per-base coverage is not requested. The summary has no elapsed time, so the same inputs print the same stdout; the progress line on stderr shows the rate.

**5. Vendor baits as a second Hs run.**
```
bammetrics qc sample.bam --preset clinical-panel --panel panel.bed \
    --hs vendor.interval_list panel.bed -R ref.fa -o out/NA12878
```
Writes `out/NA12878.vendor.hs_metrics` and `out/NA12878.panel.hs_metrics`, each named by its bait file stem. The verdict reads the one whose path sorts first (`panel`).

**6. Cohort: one `qc` per sample, then one `report`.**
```
$ for bam in bams/*.bam; do
    bammetrics qc "$bam" --preset clinical-panel --panel panel.bed -R ref.fa \
        -o runs/2026-11-02/$(basename "$bam" .bam)
  done
$ bammetrics report runs/2026-11-02/ --preset clinical-panel
index: 3 fail, 5 warn, 40 pass (worst first)
wrote runs/2026-11-02/index.html
```
Cohorts stay in the shell, make, Snakemake or Nextflow. `qc` runs one sample; `report` builds the index from any directory of per-sample outputs. There is no sample sheet runner (decision 6 in section 8).

**7. Chain: BedToIntervalList to CollectHsMetrics.**
```
bammetrics bed-to-interval-list panel.bed --from-bam sample.bam -o panel.interval_list
bammetrics hs-metrics sample.bam -R ref.fa \
    --baits panel.interval_list --targets panel.interval_list -o hs.txt
```
`--from-bam` copies @SQ lines from the BAM header. `--dict ref.dict` also works. `-o -` writes to stdout, so process substitution joins the steps. `qc --panel` does this internally and keeps `PREFIX.panel.interval_list` (0.14.0: a prepare step with the BAM header as the dictionary, recorded in the manifest). Hs reads BED directly, so the step is for Picard interop only. Giving both `--from-bam` and `--dict` is a usage error.

**8. Chain: scattered quality yield to Accumulate.**
```
parallel -j8 'bammetrics quality-yield {} -o out/{/.}.qy' ::: shards/*.bam
bammetrics accumulate-quality-yield out/*.qy -o out/sample.qy
```
Files in, file out. `AccumulateQualityYieldMetrics I=... I=... O=...` still works.

**9. Chain: MarkDuplicates to CollectDuplicateMetrics.**
```
bammetrics chain in.bam -o out/s + mark-duplicates + duplicate-metrics + write
```
One process (W04 of docs/07): `out/s.marked.duplication_metrics`, `out/s.duplicate_metrics` and `out/s.marked.bam`. `duplicate-metrics` gets the marked records in memory and reads the flags the marker set. The marker reads its input twice, so a file source is read twice and a piped one is copied to `--tmp-dir` once. The two-process form of 0.18.0 still works: `bammetrics mark-duplicates in.bam --metrics out/s.markdup_metrics -o - | tee out/s.marked.bam | bammetrics duplicate-metrics - -o out/s.dup`.

**10. One read, seven collectors, stdin.**
```
samtools view -b in.bam chr20 | bammetrics collect - -R ref.fa -o out/chr20 \
    -m quality-yield,alignment-summary,insert-size,gc-bias,mean-quality-by-cycle,quality-score-distribution,base-distribution-by-cycle
```
1 read of the input. 9 files: gc-bias writes detail, summary and chart data. `-m multiple-metrics` expands to the Picard default five.
Not built in 0.13.0: `collect -m` knows `quality-yield`, `alignment-summary`, `insert-size` and `hs`. Today `bammetrics multiple-metrics - -R ref.fa -o out/chr20 --program ...` runs those programs on one read of standard input.

**11. Report and index from files, including real Picard output.**
```
bammetrics report runs/*/ --preset clinical-panel -o site/
site/index.html  site/NA12878.html  site/NA12879.html ...
```
`report` finds metrics by `## METRICS CLASS` blocks, not by file name. It never reads a BAM.

**12. Errors and the CI gate.**
```
$ bammetrics hs-metrics sample.bam --baits panel.bed -R ref.fa
bammetrics hs-metrics: missing --targets (Picard: TARGET_INTERVALS)
  Use --panel panel.bed when baits and targets are the same file.
exit status 2

$ bammetrics qc nope.bam
bammetrics qc: nope.bam not found (from BAM)
exit status 1

$ bammetrics check out/NA12878 --preset clinical-panel --format jsonl > checks.jsonl; echo $?
11
```
`plan` is the dry run. `bammetrics plan BAM [the qc options] [--format text|jsonl]` lists the stages, each collector with its argv, outputs and rules, the BAM reads (the shared pass plus each own read), and the steps `--resume` would skip. It writes nothing. It reads the BAM header and the first and last 64 KiB of each input, as `qc` prepare does, so a file that is not a BAM exits 1 there.

## 3. Chaining model

Three layers. No workflow language. Since W04 of docs/07 `chain` joins them (docs/chaining.md): a line of steps in one process, with the records handed from step to step in memory. Layers 1 and 2 stay as they are. Layer 3 becomes a preset of `chain` in W14, when `qc` can be the last step of a chain; until then `qc` keeps its fixed graph, and a chain without a transform runs the same pass and writes the same metrics files.

1. **Pipes for BAM streams.** `-` is stdin or stdout. A BAM written to `-` uses compression level 1. `tee` saves a copy.
2. **Prefix-named files for handoff.** A step names its outputs from the prefix. The next step takes the prefix or the exact file. `PREFIX.manifest.json` lists each output with tool and kind.
3. **A fixed stage graph inside `qc`.** Order: prepare, pass, post, report.

| Stage | Content | BAM reads | As a `chain` step (stage) |
|---|---|---|---|
| prepare | bed-to-interval-list (reads BED; `qc --panel BED` since 0.14.0), mark-duplicates (writes a BAM) | 0 | mark-duplicates is a transform (a barrier: it reads its input twice); bed-to-interval-list reads no BAM and is no step |
| pass | quality-yield, alignment-summary, insert-size, any number of hs, gc-bias, by-cycle tools, duplicate-metrics | 1, shared through pass.rs | each is a tap; `write` is the tap that writes the BAM |
| post | accumulate-quality-yield, verdict | 0 | no step (W14: `qc` as the last step) |
| report | HTML and JSON | 0 | no step (W14) |

Rules:

- A collector that needs a file from an earlier step waits for it. Hs reads the `--panel` BED itself, so it does not wait for `PREFIX.panel.interval_list` (0.14.0). A collector joins the pass only when its input BAM, reference and read filter match. `plan` prints the reason when it does not.
- Per-collector error isolation (`pass::run_isolated`, 0.12.0): one failing collector marks only itself failed, the others finish. A failure in reading the BAM still fails the pass. The concordance suite gates this change.
- Resume: `qc --resume` reads `PREFIX.manifest.json`. Each step stores the argv hash, tool version, preset hash, input size, mtime, and a hash of the first and last 64 KiB of the input. A step is skipped when all match and its outputs exist. If 1 of 6 collector files is missing, the pass reruns for that collector only. `--force` reruns all. `report` always reruns because it takes milliseconds.
- Users who need loops, scatter, conditionals or clusters use the shell, make, Snakemake or Nextflow.

## 4. UI: the verdict-first per-sample report

### 4.1 What it is

`bammetrics report` and the end of `qc` write the same files: one `SAMPLE.report.html` per sample, one `index.html` per set, and `SAMPLE.verdict.json`. Each HTML file is self-contained, makes no network request and holds no script. Targets: under 20 ms to render once the metrics are in memory, under 100 KB per sample page, under 300 KB binary growth. CI measures all three with `make budgets` since 0.16.0 (section 4.8).

### 4.2 The page, top to bottom

| # | Section | Content | Data source |
|---|---|---|---|
| 1 | Header | sample, preset name, version and hash, input name, reads, tool version | manifest, quality-yield `TOTAL_READS` |
| 2 | Verdict band | one word (PASS, WARN, FAIL or INCOMPLETE) and counts, for example `1 fail, 1 warn, 12 pass, 2 not measured` | all checks |
| 3 | Failure list | one line per warn or fail: metric, value, limit, reason | rule `why` templates |
| 4 | Plots | only under the failure line they explain | see 4.4 |
| 5 | Passed checks | one collapsed line, expands to a compact list | all checks |
| 6 | Tables | one `<details>` per tool, collapsed. Each Picard column is a row with value and definition. Histograms available as collapsed tables | every metrics file |
| 7 | Provenance | command line, input size and mtime, preset hash, all thresholds | manifest, preset |

Section data sources:

- Hs table and coverage checks: CollectHsMetrics `MEAN_TARGET_COVERAGE`, `PCT_TARGET_BASES_20X`, `PCT_TARGET_BASES_100X`, `FOLD_80_BASE_PENALTY`, `PCT_USABLE_BASES_ON_TARGET`, `PCT_SELECTED_BASES`, `PCT_EXC_DUPE`, `ZERO_CVG_TARGETS_PCT`, `AT_DROPOUT`, `GC_DROPOUT`.
- Alignment table: CollectAlignmentSummaryMetrics, row `CATEGORY=PAIR` (or `UNPAIRED`), columns `PCT_PF_READS_ALIGNED`, `PCT_ADAPTER`, `PCT_CHIMERAS`.
- Insert table: CollectInsertSizeMetrics `MEDIAN_INSERT_SIZE`, `MEAN_INSERT_SIZE`, `STANDARD_DEVIATION`, plus the histogram.
- Yield table: CollectQualityYieldMetrics `PF_Q30_BASES`, `PF_BASES`, `TOTAL_READS`.
- Duplicate table: CollectDuplicateMetrics `PERCENT_DUPLICATION`.
- GC and cycle plots: CollectGcBiasMetrics detail, MeanQualityByCycle.

A column definition lives in `glossary.tsv` (tool, column, definition, unit). The text is written in our own words from Picard's column behavior, not copied from Javadoc or the Picard website. One reviewer checks it before release.

### 4.3 Thresholds and presets

A preset is a TOML file. A rule has a metric selector, an operator, a warn bound, a fail bound and a reason template:

```toml
[[rule]]
id    = "hs.target-20x"
tool  = "hs"
metric = "PCT_TARGET_BASES_20X"
warn  = "< 0.98"
fail  = "< 0.95"
why   = "{1-value:pct1} of panel bases are under 20x"
plot  = "coverage"
```

Rule language: one comparison per bound (`<`, `<=`, `>`, `>=`, `outside A..B`), an optional single ratio of 2 columns (`PF_Q30_BASES / PF_BASES`), and a selector `tool[COLUMN=VALUE].METRIC` for rows such as `CATEGORY=PAIR`. There is no `&&`, `||` or scripting, so the parser stays near 150 lines and a verdict can be read at a glance.

Status logic:

- Pass when the value is inside the warn bound. Warn when it breaks the warn bound but not the fail bound. Fail when it breaks the fail bound.
- A rule whose metrics file is absent has status `not measured`. It never counts as pass. If any rule is not measured, the band reads INCOMPLETE unless a rule already failed, in which case it reads FAIL with the not-measured count beside it.
- Presets shipped: `generic`, `clinical-panel`, `wgs`, `rna-seq`. `wgs` and `rna-seq` carry only rules for tools that exist today. Coverage rules for CollectWgsMetrics and ribosomal rules for CollectRnaSeqMetrics arrive with those tools, and until then the band reads INCOMPLETE.
- The lab owns the file. `--preset my.toml` replaces a shipped preset. `--set hs.target-20x.fail="< 0.90"` overrides one bound. The page prints the preset name, a 6-character hash of its content, and every threshold.
- Shipped numbers are not validated clinical limits. The maintainer reviewed the `clinical-panel` values (2026-11-06, version 2, note `reviewed defaults`); the other presets stay starting values. A lab validates its own copy. The verdict is not a release decision.

Default thresholds, preset `clinical-panel` version 2 (14 rules, reviewed defaults; the values of version 1):

| Rule | Source | Warn | Fail |
|---|---|---|---|
| Mean target coverage | hs `MEAN_TARGET_COVERAGE` | < 200 | < 100 |
| Target bases at 20x | hs `PCT_TARGET_BASES_20X` | < 0.98 | < 0.95 |
| Target bases at 100x | hs `PCT_TARGET_BASES_100X` | < 0.80 | < 0.60 |
| Evenness | hs `FOLD_80_BASE_PENALTY` | > 2.0 | > 3.0 |
| Usable bases on target | hs `PCT_USABLE_BASES_ON_TARGET` | < 0.50 | < 0.30 |
| Selected bases | hs `PCT_SELECTED_BASES` | < 0.70 | < 0.50 |
| Duplicate bases excluded | hs `PCT_EXC_DUPE` | > 0.30 | > 0.50 |
| Zero coverage targets | hs `ZERO_CVG_TARGETS_PCT` | > 0.01 | > 0.03 |
| GC dropout | hs `GC_DROPOUT` | > 5 | > 10 |
| Aligned reads | alignment `PCT_PF_READS_ALIGNED` | < 0.98 | < 0.95 |
| Adapter | alignment `PCT_ADAPTER` | > 0.005 | > 0.020 |
| Chimeras | alignment `PCT_CHIMERAS` | > 0.05 | > 0.10 |
| Median insert size | insert `MEDIAN_INSERT_SIZE` | outside 180..450 | outside 120..600 |
| Q30 fraction | yield `PF_Q30_BASES / PF_BASES` | < 0.85 | < 0.75 |

### 4.4 Plots

Each rule names at most 1 plot family. A plot appears only under a warn or fail line from a rule that names it. A passing sample has no plot and 1 line: `4 plots hidden, all related checks pass.` `--plots always` forces them. At most 4 plots render, which keeps print at 2 pages.

| Family | Shown when | Source | Shape |
|---|---|---|---|
| coverage | a coverage or evenness rule is warn or fail | hs `high_quality_coverage` histogram | cumulative curve, share of target at or above depth, markers at 20x and 30x |
| insert-size | insert size or adapter rule is warn or fail | insert-size histogram | bars, median line, rule limits as dashed lines |
| gc-bias | GC dropout or GC bias rule is warn or fail | gc-bias detail | normalized coverage by GC bin over a window-count strip, 1.0 line |
| cycle-quality | Q30 rule is warn or fail | mean-quality-by-cycle | 1 line per read end, Q20 and Q30 guides |

### 4.5 The run index

One row per sample, worst first. Sort: verdict rank (fail, warn, incomplete, pass), then fail count, then the largest miss as a fraction of the limit. Columns: sample (link to its page), verdict, `F/W/P` counts, then 1 cell per metric group in the preset (coverage, insert size, quality, duplicates). Each cell has a glyph and the value. A failing cell names the rule. Counts per verdict sit on top. Order is fixed at build time. A filter box would need JavaScript, which is not justified for tens of rows. The index covers only the directories passed in.

### 4.6 Visual design

- Type scale 1.2 on a 15 px base: 12.5, 15, 18, 21.6, 26 px. System font stack. Tabular numerals for every value. A monospace stack for metric names. No web font download.
- Palette: a 9-step neutral ramp, 1 ink-blue accent for data, 3 status hues used only for verdicts. Status is a glyph (check, triangle, cross, dash), a word and a color, so color is never the only signal. Text contrast is at least 4.5:1 in light and dark, checked by a unit test that computes ratios from the tokens.
- All colors are CSS custom properties on `:root`. `@media (prefers-color-scheme: dark)` swaps them. A `data-theme` attribute overrides it. Plots read the same tokens, so there is 1 SVG per plot.
- SVG: hand-built in Rust, fixed viewBox, 3 to 5 ticks per axis, direct labels instead of legends, the rule limit as a thin dashed line with its value, the failing region shaded at 8% opacity. No gridlines beyond the baseline.
- Print: `@media print` hides collapsed tables and the footer detail, drops status backgrounds, keeps the band, failure list and plots. Targets: 1 page for a pass, 2 pages with plots, on A4 and Letter.
- Layout: 1 column, text to 72 ch, tables to 1100 px. Works at 360 px and 1280 px with no horizontal scroll. No hover-only content. `make a11y` checks both pages with axe-core at 320, 360 and 1280 px, the reflow at 320 and 360 px with every `<details>` open, and the keyboard order (docs/benchmarks.md).
- No framework CSS. About 8 KB of CSS inlined per file (8,153 bytes in 0.16.0; a test fails at 8 KB). No JavaScript. A sort script (under 2 KB, inline) is justified only for indexes over about 200 samples, and is out of v1.
- Security: a `<meta>` CSP of `default-src 'none'; style-src 'unsafe-inline'`. Every string from a BAM header (sample, read group, path) is HTML-escaped, with a test using `<script>` in @RG SM.

### 4.7 How it is rendered from Rust

- New crate `bammetrics-report`, std only. HTML by `write!` into a `String`, CSS and glossary by `include_str!`, SVG by 4 functions (about 400 lines). It has no noodles or BAM dependency, so the wheel and a future hosted service can reuse it.
- New crate `bammetrics-verdict`: preset parser, rule evaluation, `verdict.json`. It needs the `toml` crate. JSON is written by hand.
- New reader in `core/src/metrics_file.rs`: parses `## METRICS CLASS` blocks by header name, multi-row classes and `## HISTOGRAM` sections. It skips unknown columns and tolerates empty histograms. `accumulate_quality_yield.rs` moves onto it. Tests run on Picard 3.5.0 output.
- Both crates take one parsed `MetricsDoc`. `report` builds it from files. `qc` builds it from rows in memory, so `qc` does not re-parse what it just wrote.
- New dependencies in total: `toml` and `clap_complete`. Each gets a review.

### 4.8 Size and speed targets

| Target | Value | Status (0.16.0, Mac Studio, release build) |
|---|---|---|
| Render, sample page | under 20 ms | 1.1 ms, the failing fixture page with every plot, median of 31 |
| Render, index of 100 samples | under 100 ms | 0.63 ms, median of 31 |
| Sample page size | under 100 KB | 88.6 KB, the largest of the three fixture pages |
| Index size, 100 samples | under 200 KB | 151.9 KB, 14 `clinical-panel` rules per sample |
| Binary growth | under 300 KB | 290.8 KB stripped (314.8 KB unstripped), the `size_report` probe less `size_base`. The probe holds the two crates and leaves out the verbs' own code (`report.rs`, `check.rs`, `discover.rs`, the `qc` engine). The whole binary grew 339.4 KB stripped when the `report` and `check` verbs first linked the crates (28647a8 to f82df56): over the target, under 2x. |
| Formats growth | under 150 KB | 129.2 KB stripped (2026-11-10), the `size_formats` probe less `size_formats_base` (`crates/formats/examples`): the format writers and readers of docs/formats.md. The whole binary grew 193.8 KB stripped with the format verbs and options. |
| Added passes over the BAM | 0 | 0: `qc` reports 1 BAM read with the report; `report` and `check` run with the BAM deleted |

`make budgets` (`crates/cli/tests/budgets.rs`) measures every row on the test fixtures and fails the build at 2x any target; CI runs it on Linux and macOS (decision 11). Each row's method and the history of the binary size are in docs/benchmarks.md.

### 4.9 Out of scope: the hosted dashboard

| Stays in the open report | Belongs to the hosted dashboard (ROADMAP.md) |
|---|---|
| verdict for 1 run | runs over time, trends per metric |
| index of the directories passed in | cohort and lab-wide views, cross-run comparison |
| provenance footer (command, preset hash) | audit trail with users and sign-off |
| thresholds in a file | shared threshold sets with change history |
| static file, no state | comments, assignments, alerts |
| none | accounts, uploads, history |

`verdict.json` (schema `bammetrics.verdict/1`: sample, preset hash, per-rule value, status, reason, input fingerprint) is the only bridge. Rules for the report crate, written into CLAUDE.md: no read of other run directories, no history store, no user identity, no network call, no comment field. One exception, from docs/06: the page reads the gaps BED beside a `CoverageGapMetrics` file (`<metrics path>.bed`, in the same run directory) for its `Coverage gaps` table, so the page bytes depend on that sidecar too. A feature that needs one of the rest goes to the hosted product.

### 4.10 Terminal output

Progress goes to stderr, only on a terminal, 10 updates per second: `reading 12.4M reads  1.9M/s  62%`. Percent is compressed bytes read over file size, so stdin shows none. `NO_COLOR` and `-q` are honored. The summary is at most 12 lines: header (sample, band, preset, collectors, BAM reads, records read; no elapsed time, so stdout is reproducible), every fail, every warn, a pass count, files written, report path. Extra failures collapse to `+N more, see report`. `--format jsonl` emits `start`, `progress`, `file`, `error`, `check` and `done` events, one per line. An error that stops the run emits an `error` event and a `done` event with its exit code.

### 4.11 Clinical checks (docs/06)

Two single-sample checks Picard does not have, specified in
[docs/06_clinical_checks.md](06_clinical_checks.md): `CollectCoverageGaps` (verb
`coverage-gaps`, class `bammetrics.qc.CoverageGapMetrics`, a BED of the target runs under
`MIN_DEPTH`, the Hs per-base depth) and `CollectSexCheckMetrics` (verb `sex-check`, class
`bammetrics.qc.SexCheckMetrics`, reads per base on X, Y and the autosomes against
`EXPECTED_SEX`, `qc --sex`). The rule language gains one text comparator, `is TEXT`, for
`sex.match` (fail `is N`) and `sex.inferred` (warn `is UNKNOWN`); `verdict.json` checks
judged as text carry a `text` field. `clinical-panel` version 3 adds the four rules
(`gaps.reportable-range`, `gaps.targets-complete`, `sex.match`, `sex.inferred`); `generic`
and `wgs` are unchanged. The page lists the 200 longest gaps under the first
warn or fail gap check, from the BED beside the metrics file. Both stay inside the rules
of section 4.9: one run's files, no history, no cohort.

## 5. Python mirror

Same verbs, same snake_case option names, same results. The generated registry functions and the older typed functions stay.

```python
import bammetrics as bm

r = bm.insert_size("s.bam", minimum_pct=0.05)       # typed result, as today
run = bm.collect("s.bam", metrics=["quality-yield", "insert-size", "gc-bias"],
                 reference="ref.fa", output="out/s", threads=4)
run.metrics["insert-size"]; run.files

run = bm.qc("s.bam", preset="clinical-panel", panel="panel.bed",
            reference="ref.fa", output="out/s", report=True)
run.verdict.status                  # "pass" | "warn" | "fail" | "incomplete"
for c in run.verdict.checks: ...    # Check(rule, status, value, limit, why)
run                                 # in a notebook: the report body via _repr_html_

m = bm.read_metrics("out/s.insert_size_metrics")    # rows and histograms
v = bm.check("out/s", preset="clinical-panel")
bm.report(["runs/a", "runs/b"], preset="clinical-panel", output="site/")
bm.tools(); bm.presets()

run = bm.chain("s.bam", [bm.step("mark_duplicates"), "insert_size", bm.step("write", index=True)],
               output="out/s")       # W04 of docs/07: the chain verb, one process
run.stages[0].records_out           # one Stage per step
```

- `qc`, `collect`, `report` and `chain` release the GIL.
- `chain` and `step` mirror the `chain` verb (docs/python.md, docs/chaining.md section 5.3). A step is a verb name or `step(verb, **options)`, checked when it is made against the keywords of the tool's function. This is a list of steps, as on the command line. It is no fluent builder.
- A failed verdict is a value. It raises only with `fail_on="fail"`.
- Errors: `UsageError(ValueError)` and `InputError(OSError)`. Existing `except ValueError` and `except OSError` keep working.
- Option names come from the same table as the CLI, and `.pyi` is generated from it. This ends `reference` against `reference_sequence` and `path` against `input` for new functions. Old functions keep their signatures and are documented as aliases.
- The notebook repr is the report body in `<div class="bm">` with scoped CSS and inline SVG. Histograms are downsampled to 200 points.
- CI builds the wheel, imports it and compares `dir(bammetrics)` with the stub. This catches the stale wheel now in `.venv` (`make wheel-check`, 0.15.0).
- Not added: a fluent builder, a `Cohort` class, pandas frames, Python callables as checks. A verdict must reproduce from `verdict.json` alone.

## 6. Compatibility and migration

- The Picard face is untouched. `registry.rs` and `picard_cli.rs` change no output byte and no exit code on a command that works today.
- Three Picard-face error paths become features: `--help` and `-h` print the option table (they error today), an unknown class name lists the tools with a did-you-mean, and `--version` works. `O=-` still creates a file named `-`, as listed in KNOWN_DIFFERENCES.md. `-` as stdout and stdin exists only on the new verbs.
- The 5 clap subcommands keep every flag, default and output byte. `-O` stays beside `-o`. `collect` keeps `--metric`, `--hs BAITS TARGETS` and `--hs-per-base`. It gains metric names and a default prefix. A command that errored before now works. No command that worked changes.
- Guard: a golden suite records argv, stdout, written files and exit code for every command in the docs and the concordance corpus. CI fails on any byte change. It lands in the first PR and gates all later ones. The existing 90-case concordance suite keeps running.
- Docs: `docs/cli.md` is rewritten around verbs, with a table from each Picard class to its verb. README gets one `qc` example.
- Defaults that differ by design: native verbs write to stdout without `-o` like the current subcommands, and Picard syntax keeps Picard's rules.

## 7. Implementation plan

Each phase is 1 PR-sized unit. Days are for 1 engineer. Total: 73 days (15 weeks). The first release with a visible report is after phase 6, at day 40.

| # | Phase | Content | Tests | Days |
|---|---|---|---|---|
| 1 | Golden suite and fixtures | record argv, stdout, files, exit code for all documented commands and the corpus. Generate the corpus under tests/corpus/out. | the suite itself, run in CI | 4 |
| 2 | Metrics reader and glossary | `metrics_file.rs` reader, `glossary.tsv`, `accumulate_quality_yield` moved to it | round trip on bammetrics output, parse of Picard 3.5.0 files, extra columns, empty histograms, several class blocks, locale decimals | 7 |
| 3 | Verdict crate and presets | rule parser, evaluation, `generic` and `clinical-panel`, `verdict.json` | parser cases, every status incl. not measured and INCOMPLETE, preset hash stable | 7 |
| 4 | Sample page and plots | HTML builders, 4 SVG plots, CSS, dark mode, print | snapshot HTML, contrast unit test, escaping test with `<script>` in @RG SM, render timing | 12 |
| 5 | Run index | worst-first table, group columns, links | sort order cases, 100-sample size and time check | 4 |
| 6 | `report` and `check` verbs | file discovery by class block, exit codes 10 and 11, release 0.7 | end to end on Picard output, exit code table | 4 |
| 7 | `qc` engine | stage graph, manifest, resume, atomic writes, `plan`, summary, progress, per-collector error isolation, jsonl | resume after a killed run, 1 missing file reruns 1 collector, golden suite unchanged, concordance unchanged | 10 |
| 8 | Verb grammar | option table in registry.rs, normaliser, bool forms, `tools`, `presets`, `help`, did-you-mean, completions, `-` for stdin and stdout, output naming table | every old argv still byte-identical, spelling matrix per tool, name-table pin vs Picard face, stdin BAM case | 10 |
| 9 | `bed-to-interval-list` | BED to interval list with @SQ from BAM or dict (0.14.0) | compare against Picard BedToIntervalList output on 3 BEDs | 2.5 |
| 10 | Python mirror | `qc`, `check`, `report`, `read_metrics`, `_repr_html_`, generated stubs, wheel check (0.15.0; docs/python.md) | stub vs `dir()` test, GIL release, error classes | 8 |
| 11 | Docs and benchmark | `cli.md` rewrite, preset docs, CI size and time budgets, accessibility pass (0.16.0; `make budgets`, docs/benchmarks.md) | benchmark fails at 2x target | 4.5 |

Order of delivery: phases 1 to 6 as 0.7 (report and check on any Picard directory, no change to collection). Phase 7 as 0.8. Phases 8 and 9 as 0.9. Phases 10 and 11 as 0.10. `mark-duplicates` is outside the 73 days because it waits on parity plan P9 (6 to 8 weeks).

Risks, each with its mitigation:

- Wrong preset defaults mislead a lab: starting-value label, hash on the page, lab-owned file.
- `wgs` and `rna-seq` stay INCOMPLETE until their collectors exist: not-measured status.
- The bool normaliser could swallow a path: only literal `true` or `false` after a known boolean option, plus the golden suite.
- Output name table drifts between faces: a pinned test.
- Scope creep into the hosted product: the table in 4.9.
- Added dependencies: `toml` and `clap_complete` only, each reviewed.

## 8. Open decisions for the maintainer

1. Does `qc` without `--preset` use `generic` (this proposal), or must it require a preset? Decided: `generic` (0.12.0).
2. Are the 14 clinical-panel starting values in 4.3 acceptable as shipped defaults, or does a lab review them first? Decided: they ship labelled "starting values" after the maintainer reviews them; the review closed on 2026-11-06: the 14 values are kept as `clinical-panel` version 2 with the note `reviewed defaults`, and a lab still validates its own copy.
3. Should `qc` exit 0 on a failed verdict by default (this proposal), or gate by default? Decided: exit 0 unless `--fail-on` is given; `check` always gates (0.12.0).
4. Is `toml` acceptable as a dependency, or should presets use a hand-written format? Decided: hand-written. `crates/verdict/src/table.rs` reads the TOML subset the presets use, with no dependency; `parse_tables` is the one function a later switch to the `toml` crate replaces (2026-11-06). The decision holds in docs/07 (D12): in W04 `parse_tables` moved to `crates/core/src/toml.rs`, where the plan file of `chain --plan` reads the same subset plus one-line arrays; the verdict crate calls it and names the preset file in the error.
5. Is the 73-day plan acceptable, or should phases 9 and 10 move after 0.8? Decided: the plan ran in order; phases 7 to 11 shipped as 0.12.0 to 0.16.0.
6. Do you want a sample sheet runner in this scope, or leave cohorts to Nextflow plus `report`? Decided: no sample sheet runner; `qc` runs one sample, cohorts loop in the shell or a workflow engine and finish with `report` (example 6). Its phase is dropped.
7. Should the old subcommands be listed as aliases in `help` and docs, or hidden after 0.9? Decided: neither. `insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics` and `collect` are verb-table entries since 0.13.0, listed in `help`, `tools` and docs/cli.md as verbs like the rest (2026-11-06).
8. Who writes and reviews the column glossary text (about 400 entries, 1 week of review)? Decided: `crates/core/glossary.tsv` (317 rows) is cross-checked against the Picard metric classes and the code that computes each column, and the maintainer approves the disagreements found (2026-11-06).
9. Is a `data-theme` override for dark mode wanted, or only `prefers-color-scheme`? Decided: both. The page follows `prefers-color-scheme`; `data-theme="light"` or `"dark"` on the html element forces one (2026-11-06).
10. Does `mark-duplicates` ship in the P9 PR train, and should the verb slot exist as an error stub from 0.9? Decided: MarkDuplicates is the next feature release after 0.17.0; no stub, so `mark-duplicates` stays an unknown verb until then (2026-11-06).
11. Do the size and speed budgets in 4.8 become hard CI failures? Decided: yes, at 2x the target. `make budgets` fails the CI job when a measured value reaches twice its target and prints `over target` when it passes the target alone (0.16.0).
12. Which name for the report file: `SAMPLE.report.html` (this proposal) or `SAMPLE.html`? Decided: `SAMPLE.report.html`, as shipped since 0.12.0 (2026-11-06).
13. Do the coverage gap and sex check rules of docs/06 join `clinical-panel`, and how are their thresholds labelled? Decided: `clinical-panel` version 3 (hash `6f0d59`) carries all four with the note `reviewed defaults; gap and sex rules are starting values`; `generic` and `wgs` are unchanged at version 1 (hashes `6dd89e` and `8a5fc4`, byte-identical to 0.18.0), so `qc --sex` under them exits 2 and no not-measured sex check appears under them. The sex thresholds are starting values measured on HG002 and NA12878 (docs/06 section 4); the 14 reviewed values are unchanged (reviewer decision, 0.19.0, 2026-11-09).

Result: 12 verbs and 1 prefix cover every chain in the audit, the Picard face stays byte-identical, and the verdict report ships first at day 40.
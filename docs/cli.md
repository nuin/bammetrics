# Command line and Python reference

Every verb, its options and the Picard-syntax forms. The short version is in the [README](../README.md).

## Build and test

Needs Rust 1.85 or newer (edition 2024).

```bash
cargo build --release            # target/release/bammetrics
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
```

## Shape

```
bammetrics <verb> [INPUT] [options]       kebab-case verbs
bammetrics <PicardClass> KEY=VALUE ...    the Picard drop-in, unchanged
```

`main.rs` dispatches in this order: a first argument that names a Picard class in the
registry runs the Picard face, untouched; then a grammar verb (a tool verb, `collect`,
`tools`, `presets`, `help`, `completions`); then the clap verbs `qc`, `plan`, `check` and
`report`. Class names are CamelCase and verbs are kebab case, so the two cannot collide.

```bash
bammetrics insert-size in.bam                               # Picard text on stdout
bammetrics insert-size in.bam -o in.insert_size_metrics --minimum-pct 0.1
bammetrics insert-size in.bam --MINIMUM_PCT 0.1 --format json | jq '.metrics[0].MEDIAN_INSERT_SIZE'
samtools view -b in.bam chr20 | bammetrics quality-yield - -o chr20.quality_yield_metrics
bammetrics gc-bias in.bam -R ref.fa -o out/NA12878          # out/NA12878.gc_bias.*_metrics
bammetrics hs-metrics in.bam --panel panel.bed -R ref.fa -o out/NA12878.hs_metrics
bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt M=0.5
```

## Tool verbs

One verb per tool in the registry (`crates/core/src/registry.rs`). The name is the Picard
class without `Collect` and `Metrics`, in kebab case; the first twelve are the names of
docs/05 section 2.1. `bammetrics tools` prints this table.

| Verb | Picard class | Files |
|---|---|---|
| `insert-size` | CollectInsertSizeMetrics | one |
| `alignment-summary` | CollectAlignmentSummaryMetrics | one |
| `quality-yield` | CollectQualityYieldMetrics | one |
| `hs-metrics` | CollectHsMetrics | one (plus the per-base, per-target and sensitivity files asked for) |
| `gc-bias` | CollectGcBiasMetrics (needs `-R`) | `PREFIX.gc_bias.detail_metrics`, `PREFIX.gc_bias.summary_metrics` |
| `quality-score-distribution` | QualityScoreDistribution | one |
| `mean-quality-by-cycle` | MeanQualityByCycle | one |
| `base-distribution-by-cycle` | CollectBaseDistributionByCycle | one |
| `duplicate-metrics` | CollectDuplicateMetrics | one |
| `library-complexity` | EstimateLibraryComplexity | one |
| `multiple-metrics` | CollectMultipleMetrics | `PREFIX.<program extension>` per program |
| `accumulate-quality-yield` | AccumulateQualityYieldMetrics | one |
| `wgs` | CollectWgsMetrics (needs `-R`) | one |
| `raw-wgs` | CollectRawWgsMetrics (needs `-R`) | one |
| `wgs-with-non-zero-coverage` | CollectWgsMetricsWithNonZeroCoverage (needs `-R`) | one |
| `targeted-pcr` | CollectTargetedPcrMetrics | one |
| `rna-seq` | CollectRnaSeqMetrics | one |
| `variant-calling` | CollectVariantCallingMetrics (input is a VCF) | `PREFIX.variant_calling_{summary,detail}_metrics` |
| `accumulate-variant-calling` | AccumulateVariantCallingMetrics (inputs are prefixes) | as `variant-calling` |
| `genotype-concordance` | GenotypeConcordance (truth and call VCF as positionals) | `PREFIX.genotype_concordance_{summary,detail,contingency}_metrics` |
| `bed-to-interval-list` | BedToIntervalList | one interval list |

A verb runs its tool through the registry, so `bammetrics gc-bias BAM -R ref.fa -o P` and
`bammetrics CollectGcBiasMetrics I=BAM R=ref.fa O=P.gc_bias.detail_metrics
S=P.gc_bias.summary_metrics` write the same metrics (the `StringHeader` command line lists
the options in the order the verb passes them). `crates/cli/tests/grammar.rs` checks that
for every tool.

## Options

Each option has one row in the registry option table: its Picard names, value type,
whether it repeats or is required, and its help text. The row gives the kebab flag, the
Picard key, the `--help` line, the Python keyword argument and the `bammetrics.pyi` stub.
`bammetrics help VERB` (or `bammetrics VERB --help`, or `bammetrics PicardClass --help`)
lists every row of a verb.

Spellings. `--minimum-pct 0.1`, `--minimum-pct=0.1`, `--MINIMUM_PCT 0.1`, `MINIMUM_PCT=0.1`,
`minimum_pct=0.1`, `-M 0.1` and `M=0.1` are one option. A name is upper-cased with `-` read
as `_`, then matched against the row's Picard names, and the registry parser gets
`MINIMUM_PCT=0.1`. A bare word holding `=` is an option only when the part before `=` has no
path character, so `dir/a=b.bam` stays an input.

Booleans take `--x` (true), `--no-x`, `--x=false`, `X=false` and the space form `--x false`
of the 0.12 subcommands. Only a literal `true` or `false` (any case) after a boolean is read
as its value: in `--include-duplicates in.bam`, `in.bam` is the input.

Shared options, on every tool verb:

| Option | Meaning | Picard spelling |
|---|---|---|
| INPUT or `-i` | input file; `-` reads standard input | `I`, `INPUT` |
| `-o` (`-O` kept) | file, prefix, or `-` for standard output | `O`, `OUTPUT` (`M` for duplicate-metrics) |
| `-@`, `--threads N` | BGZF decompression workers; default the available cores, at most 8 (`--help` prints the number) | none |
| `-R` | FASTA with a `.fai` | `R`, `REFERENCE_SEQUENCE` |
| `--level` | accumulation level, repeatable | `LEVEL`, `METRIC_ACCUMULATION_LEVEL` |
| `--stop-after N` | read limit | `STOP_AFTER` |
| `--format` | `picard`, `json`, `tsv`, `jsonl`, `human` (single-file verbs) | none |
| `--sample NAME` | sample name in the json and jsonl forms; never written into a Picard file | none |
| `-q`, `-v` | no notes on stderr; the Picard command line and each file written on stderr | `QUIET` |

`--level` and `--stop-after` exist on a verb whose tool has the Picard option; elsewhere
they are unknown options, as the Picard face rejects them. `-R` on a tool that reads no
reference is accepted and noted, as `R=` is on the Picard face. Records are decoded and counted
on one thread in file order, so the metrics do not depend on `--threads`.

Options Picard has that cannot change a metric (`VALIDATION_STRINGENCY`, `TMP_DIR`,
`VERBOSITY`, `H`, `CHART`, `R` on a tool that reads no reference, ...) are accepted,
dropped and noted on stderr. An option no row names is an error, because dropping an option
Picard applies would change the numbers silently.

Errors on a verb read `bammetrics <verb>: <what> (from <option>)`, then at most two hint
lines, and exit 2 for a usage error (nothing ran) and 1 for an input that cannot be read or
a tool that fails:

```
$ bammetrics insert-size in.bam --minmum-pct 0.1
bammetrics insert-size: unknown option --minmum-pct
  did you mean --minimum-pct (Picard MINIMUM_PCT)?
  every option: bammetrics help insert-size
$ bammetrics hs-metrics in.bam --baits panel.bed
bammetrics hs-metrics: missing --target-intervals (Picard: TARGET_INTERVALS)
  Use --panel FILE when baits and targets are the same file.
  every option: bammetrics help hs-metrics
$ bammetrics quality-yield nope.bam
bammetrics quality-yield: nope.bam not found (from BAM)
```

A usage error that clap finds on `qc`, `plan`, `check` or `report` uses the same form, with
the usage line as a hint.

## Output naming

- A single-file verb writes `-o FILE`. Without `-o`, or with `-o -`, the file goes to
  standard output.
- A multi-file verb (`gc-bias`, `multiple-metrics`, `variant-calling`,
  `accumulate-variant-calling`, `genotype-concordance`, `collect`) writes under `-o PREFIX`.
  Without `-o` the prefix is `./` plus the input name without its extension, and stderr says
  so. An `-o` that ends in `/` or names an existing directory gets that name as the file
  part. `-o -` is an error.
- The names after the prefix come from one table, `multiple_metrics::file_extensions`, the
  `CollectMultipleMetrics` program extensions (`.alignment_summary_metrics`,
  `.insert_size_metrics`, `.quality_yield_metrics`, `.gc_bias.detail_metrics`, ...) plus
  `.hs_metrics`, `.duplicate_metrics`, `.wgs_metrics` and `.rna_metrics`. `collect`, `qc`
  and the verbs read it; a test pins it, and another checks the verbs against
  `CollectMultipleMetrics O=PREFIX`.
- Missing directories above an output are created.
- `-` as the input reads the BAM from standard input on the verbs only; the header line
  reads `INPUT=-`. The Picard face reads a file named `-` and writes `O=-` to a file named
  `-`, as listed in KNOWN_DIFFERENCES.md.

## Formats

`--format picard` (the default) is the Picard metrics file. `json` is one object: `tool`,
`sample` (with `--sample`), `class`, `metrics` (one object per row, keyed by column) and
`histogram` (`columns`, `rows`, or `null`). `jsonl` is one object per row, `tsv` a header
line and one line per row, `human` an aligned table (a row per line up to 6 columns, else
`COLUMN  value` lines per row). Numbers are JSON numbers, NaN and infinity `null`; the text
forms print Picard's cell text. For `QualityScoreDistribution` and `MeanQualityByCycle`,
whose file has no metrics class, the rows are the histogram.

## The 0.12 subcommands

`insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics` and `collect` were clap
subcommands up to 0.12. They are verbs now, with every flag and default kept, and every
command that worked writes the same bytes (the golden suite holds them). The four tool verbs
also write the `StringHeader` line the subcommand wrote, which lists every option with its
value.

- `insert-size`: `-M` (MINIMUM_PCT), `--deviations`, `-W` (HISTOGRAM_WIDTH),
  `--min-histogram-width` (MINIMUM_HISTOGRAM_WIDTH) and `--include-duplicates`. Like Picard,
  it writes no file when no orientation is kept (single-end input, for example).
- `alignment-summary`: `-R` (a FASTA with a `.fai`, read one contig at a time),
  `--max-insert-size`, `--expected-pair-orientation`, `--adapter-sequence`,
  `--collect-alignment-information`, `--is-bisulfite-sequenced`, `--level` and
  `--stop-after`. The two list options add to Picard's default (`[FR]`, six Illumina
  adapters) and `null` clears it. Without `-R` every row is written and no base counts as a
  mismatch. With `-R` the BAM `@SQ` lines must match the `.fai` (names, lengths, order) and
  the records must be coordinate sorted, or it exits 1 as Picard does.
- `quality-yield`: `--use-original-qualities <bool>` (`OQ`; default true, so an `OQ` tag
  replaces `QUAL`), `--include-secondary-alignments`, `--include-supplemental-alignments`.
  It needs no index and no sort order.
- `hs-metrics`: `--bait-intervals` and `--target-intervals` (repeatable; `BI`, `TI`;
  `--baits` and `--targets` too), `--panel FILE` for both, `--bait-set-name` (`N`), `-R`
  (needed for `AT_DROPOUT` and `GC_DROPOUT`), `--minimum-mapping-quality` (`MQ`, 20),
  `--minimum-base-quality` (`Q`, 20), `--clip-overlapping-reads <bool>` (true),
  `--include-indels`, `--coverage-cap` (`covMax`, 200), `--near-distance` (250) and
  `--sample-size` (10000). An interval file is an interval list when its first line starts
  with `@` (the `@SQ` lines must equal the BAM's), otherwise BED (0-based, read against the
  BAM header). Picard 3.5.0 refuses BED files; bammetrics reads them as `BedToIntervalList`
  would convert them. The `HET_SNP_SENSITIVITY` simulation draws 10000 x 201 qualities from
  the same `java.util.Random(51)` stream as Picard, about 0.7 s. `--per-base-coverage FILE`
  writes Picard's per-base file: a `chrom pos target coverage` header, then one line per base
  of every merged target in interval order with the 1-based position, the target name (names
  of merged targets joined by `|`) and the high-quality depth; a name ending in `.gz` is
  BGZF-compressed.
- `collect`: several metrics on one read of the BAM, one Picard file each under `-o PREFIX`
  (default `./<bam stem>`; it was a usage error before 0.13.0): `PREFIX.quality_yield_metrics`,
  `.alignment_summary_metrics`, `.insert_size_metrics` and `.hs_metrics`. `--hs BAITS
  TARGETS` adds one `CollectHsMetrics` run and may be repeated; with more than one, each file
  is `PREFIX.<BAIT_SET>.hs_metrics`, BAIT_SET being the bait file name up to its first dot.
  `--hs-per-base` also writes each run's per-base coverage, named as its HS file with
  `.per_base_coverage` in place of `.hs_metrics`. `-R` goes to alignment summary and hs.
  `--metric` (or `-m`) picks a subset, comma separated (default the first three, plus hs
  with `--hs`). Every other option is Picard's default. A metric that fails (an unsorted BAM
  with `-R`, for example) writes no file and exits 1 after the others are written. In
  `collect` each metric runs on a thread of its own, fed the same record batches in file
  order. Timings are in [docs/benchmarks.md](benchmarks.md).

```bash
bammetrics collect in.bam -o out/sample -R ref.fa \
    --hs panel.bed panel.bed --hs baits.interval_list panel.bed
```

## Picard syntax

```bash
bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt H=hist.pdf M=0.5
bammetrics CollectInsertSizeMetrics --INPUT in.bam --OUTPUT out.txt \
    --METRIC_ACCUMULATION_LEVEL SAMPLE --METRIC_ACCUMULATION_LEVEL READ_GROUP
bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt R=ref.fa
bammetrics CollectHsMetrics I=in.bam O=out.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
```

The Picard face takes `KEY=VALUE`, `--KEY VALUE`, `--KEY=VALUE` and `-K VALUE` with long
or short names, writes `O=` (standard output without it) and exits 0 or 1. It is
byte-identical to 0.12, apart from `--help` alone, which now prints the verb help.
`ASSUME_SORTED=false` exits 1 on a BAM whose header is not coordinate sorted, as Picard
does. `METRIC_ACCUMULATION_LEVEL` (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`; `null`
clears the default) works where Picard has it.

## Info verbs

- `bammetrics tools`: every verb and the Picard class it runs, then the other verbs.
- `bammetrics presets`: the built-in verdict presets with version, the first 6 characters
  of the content hash, rule count and note (`clinical-panel   v1  3d0bcd  14 rules  starting
  values`).
- `bammetrics help VERB|PicardClass`: the verb's options in every spelling; `help qc` (or
  `plan`, `check`, `report`) prints clap's help.
- `bammetrics completions bash|zsh|fish|elvish|powershell`: a completion script on stdout
  (`bammetrics completions zsh > ~/.zfunc/_bammetrics`).

## Python

```python
import bammetrics as bm

rows = bm.collect_insert_size_metrics("s.bam", minimum_pct=0.1, metric_accumulation_level=["SAMPLE"])
files = bm.collect_gc_bias_metrics("s.bam", reference_sequence="ref.fa",
                                   output="out/s.gc", summary_output="out/s.gc_summary")
```

One function per registry tool, named after the Picard class in snake case. A keyword is a
row's Picard name in lower case (`minimum_pct`; a short name such as `m` works too), a list
is a repeated option, `None` leaves the option out, `output=` writes the files. The
generated block of `bammetrics.pyi` holds one typed signature per tool from the option table;
`UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs` rewrites it and the test fails
while the stub and the table disagree. The older typed functions stay:
`insert_size(path, minimum_pct=0.05, ...)`, `alignment_summary(path, reference=None, ...)`,
`quality_yield(path, ...)` (one `QualityYieldMetrics`), `hs_metrics(path, bait_intervals,
target_intervals, reference=None, ..., per_base_coverage=None)` (one `HsMetrics`; every
column an attribute in lower case, `to_dict()` gives all of them, `high_quality_coverage`
and `unfiltered_base_quality` are the two histograms) and `collect(path, metrics=None,
reference=None, hs=None, output_prefix=None, threads=None, hs_per_base=False)` (`{name:
result}`, writing the files when `output_prefix` is given).

`bammetrics.read_metrics(path)` reads a metrics file written by Picard or bammetrics:
`{"headers": [...], "classes": [{"class", "columns", "rows"}], "histogram": {"key_class",
"columns", "rows"} or None}`, every row a dict keyed by column. Cells are typed from
their text: int, float, str, or None for an empty field; `?` is `nan`, and `0,5` from a
comma-decimal locale is `0.5`. The same reader is `metrics_file::MetricsDoc` in Rust,
and `crates/core/glossary.tsv` defines every column the tools write.

`bammetrics.evaluate(preset, paths, sample=None)` evaluates a verdict preset (a shipped
name from `bammetrics.presets()`, or the path of a preset file) over one sample's
metrics files and returns `verdict.json` as dicts: band `status`, `counts`, and one
check per rule with its value, bounds, status and reason. Presets, the rule language
and the JSON schema are in [docs/verdict.md](verdict.md). The CLI writes the same
verdict through `report` and gates on it through `check`, below. `bammetrics.report()`
in Python arrives with the Python mirror (phase 10 of docs/05).

## qc and plan

```bash
bammetrics qc BAM [-o PREFIX] [-R FASTA] [--preset NAME|FILE] [--hs BAITS TARGETS]... [--panel FILE]
              [--ref-flat FILE --strand STRAND [--ribosomal-intervals FILE]] [--sample NAME]
              [--fail-on warn|fail] [--resume | --force] [--threads N] [--format text|jsonl] [-q]
              [--all-plots] [--date ISO]
bammetrics plan BAM [the same options]
```

`qc` runs the collectors the preset's rules need (`generic` without `--preset`), on one
read of the BAM, then writes the verdict and the page. `--hs` or `--panel` adds
`CollectHsMetrics` and `--ref-flat` adds `CollectRnaSeqMetrics` whether or not a rule reads
them. A rule whose collector cannot run is listed as not measured with the reason (`needs
--hs BAITS TARGETS or --panel FILE`, `needs -R FASTA`), and the verdict reads INCOMPLETE.

| Metrics class | Collector | Needs | File |
|---|---|---|---|
| QualityYieldMetrics | CollectQualityYieldMetrics | | `PREFIX.quality_yield_metrics` |
| AlignmentSummaryMetrics | CollectAlignmentSummaryMetrics | (`-R` used when given) | `PREFIX.alignment_summary_metrics` |
| InsertSizeMetrics | CollectInsertSizeMetrics | | `PREFIX.insert_size_metrics` |
| DuplicationMetrics | CollectDuplicateMetrics | | `PREFIX.duplicate_metrics` |
| HsMetrics | CollectHsMetrics, once per `--hs`/`--panel` | intervals | `PREFIX.hs_metrics`, `PREFIX.BAIT_SET.hs_metrics` with several |
| GcBiasSummaryMetrics, GcBiasDetailMetrics | CollectGcBiasMetrics | `-R` | `PREFIX.gc_bias.detail_metrics`, `PREFIX.gc_bias.summary_metrics` |
| RnaSeqMetrics | CollectRnaSeqMetrics | `--ref-flat`, `--strand` | `PREFIX.rna_metrics` |
| WgsMetrics | CollectWgsMetrics, on a read of its own | `-R` | `PREFIX.wgs_metrics` |

File names come from the `CollectMultipleMetrics` extension table
(`multiple_metrics::file_extensions`), so both faces write the same names. Each metrics
file is byte for byte what `bammetrics <Collector> INPUT=BAM OUTPUT=FILE ...` writes for
the argv the manifest records, apart from the `Started on` line. The header command line
follows that argv: `qc` passes `REFERENCE_SEQUENCE` only to the collectors that read `-R`
(alignment summary, hs, gc bias, wgs), so the insert size and duplicate headers lack the
`REFERENCE_SEQUENCE=` that `bammetrics CollectInsertSizeMetrics ... R=ref.fa` prints; the
data lines are the same. The sample name (`--sample`, else the first `@RG SM`, else the
BAM name) goes only into the verdict, the page and the manifest. No `-o` means `./<BAM name
without .bam>`, and the summary says so. An `-o` that ends in `/` or names an existing
directory gets the BAM name as the file part (`-o out/` writes `out/NA12878.*`), and the
summary says so. Missing directories are created. Every file is written to a temporary
name opened with `create_new` and renamed into place. A file a collector plans and does
not write this time (it failed, or insert size dropped every orientation and Picard writes
no file) is removed, and so is a metrics file of an earlier run's collector that the new
options drop (`--panel` given before, not now), so the prefix holds what the manifest
lists; the summary counts both as stale removed. With several `--hs` runs the verdict
reads the HsMetrics file whose path sorts first (`PREFIX.<bait stem>.hs_metrics`); every
file is written and listed in the manifest. Two bait files with one stem are an error.

Stages, in order: prepare (the BAM header, the input fingerprint, the preset, every file an
option names; nothing written), pass (every one-pass collector on one read; a collector
that fails to build, fails on a record or panics fails alone, and a BAM read error fails
the pass), post (the verdict over the files), report (`PREFIX.verdict.json` and
`PREFIX.report.html`, the bytes `bammetrics report` writes for the same files and sample).
`CollectWgsMetrics` walks reference loci and reads the BAM on its own.

`PREFIX.manifest.json` (schema `bammetrics.manifest/1`) holds per step its stage, tool,
Picard argv and the argv hash, the tool version, the preset hash, the input path, size,
mtime and a hash of its first and last 64 KiB, the same three for every further file the
argv names (`files`: the FASTA and its `.fai`, the bait and target intervals, the refFlat,
the ribosomal intervals), the outputs with tool and kind, and the status (`ok` or
`failed`, with the error). `--resume` skips a collector when every one of those fields
matches and its files exist: one lost file reruns that collector alone, in a pass of its
own, and an interval file edited in place reruns the collectors that read it. A BAM (or
FASTA) larger than 128 KiB edited in the middle with its size and mtime kept (`touch -r`,
`cp -p` over it) is not detected; use `--force`. A manifest that does not parse (a killed run) is reported on stderr and
every step runs. `--force` reruns every step. Without either, `qc` refuses a prefix that
holds an earlier run's manifest (exit 2). The verdict and the page are always rewritten.

`plan` prints the stages, the BAM reads (`2 BAM reads: 1 for the pass, 1 own read`), each
collector with its files and rule count, which collectors share the read and why one does
not, the classes no collector can write, and which steps `--resume` would skip. It reads
the BAM header and the first and last 64 KiB of the BAM and of every further input file,
and writes nothing.
`--format jsonl` gives one `plan` event, a `step` event per collector (with `argv`,
`outputs`, `rules`, `read`, `action`, `resumable`) and an `unmet` event per class.

The `qc` summary on stdout is at most 12 lines: the header (sample, band, preset, collectors,
BAM reads, records the pass read: `4 collectors, 1 BAM read (41.8M records)`), every error, fail, warn and not-measured line (collapsed into `+N more, see
PAGE`), the pass count, the files written (kept by `--resume`, replaced, stale removed), and
the gate hint or the exit reason. Progress (`reading 12.4M reads  1.9M/s  62%`) goes to
stderr only on a terminal, at most 10 times a second; `-q` turns it off. `--format jsonl`
prints `start`, `progress` (every million records), `file` (with `written` false for a
file `--resume` kept), `error`, `check` and `done` events on stdout instead. A closed
stdout (`| head -1`) drops the rest of the output; the run still writes its files and keeps
its exit code. An error that stops `qc` or `plan` prints `bammetrics qc: WHAT` (or
`bammetrics plan:`) on stderr and, with `--format jsonl`, an `error` event (`step` null) and
a `done` event with only `exit` on stdout.

`qc` exits 0 on any verdict unless `--fail-on` is given; then 10 when the verdict reaches it
and 11 when the verdict is incomplete, as `check` does. A failed collector exits 1 after the
other files, the verdict and the page are written.

## report and check

```bash
bammetrics report DIR... [--preset NAME|FILE] [-o OUT] [--all-plots] [--date ISO]
bammetrics check DIR|VERDICT.json... [--preset NAME|FILE] [--fail-on warn|fail] [--format text|json|jsonl]
```

Both work on any directory of Picard output; neither reads a BAM. A metrics file is
any file whose text holds a `## METRICS CLASS` block (or a `## HISTOGRAM` section
under a Picard header), whatever its name, at any depth under `DIR`. Files that look
like metrics and do not parse are listed on stderr and left out. Samples: the files in
one directory are grouped by the file name up to the first dot, so
`NA12878.insert_size_metrics` and `NA12878.hs_metrics` form the sample `NA12878` (the
prefix `CollectMultipleMetrics` writes). When no metrics class repeats across the
groups of a directory, the directory holds one sample, named by the prefix when two or
more files share it, else by the directory: `runs/NA12878/insert_size.txt` and
`runs/NA12878/hs.txt` form the sample `NA12878`. Two samples with one name is an error.

`report` evaluates the preset (`generic` without `--preset`; a shipped name or a TOML
file) and writes `OUT/SAMPLE.verdict.json`, `OUT/SAMPLE.report.html` and
`OUT/index.html`, `OUT` being `DIR` when one directory is given. Each file is written
whole and renamed into place. The summary on stdout is one screen: a header line, one
line per sample with its band and counts, the fail and warn reasons (three per sample,
then `+N more`), and the files written. The index lists the samples worst first (fail,
incomplete, warn, pass, then by fail count, then name) with one cell per rule.
`--all-plots` draws every plot family the files have data for; `--date` fixes the date
printed on the pages, so the same inputs give the same bytes. `report` exits 0 on any
verdict: a Nextflow task is not killed by a bad sample.

`check` is the gate. A directory is discovered and evaluated as `report` does; a
`verdict.json` file is read back as is. `--format text` prints a band line per sample
and every check under it; `json` prints a JSON array of verdicts, each as
`verdict.json`; `jsonl` prints one object per line, a `check` event per check and a
`verdict` event per sample.

| Exit | Meaning |
|---|---|
| 0 | every sample passes, or none reaches `--fail-on` |
| 1 | an input could not be read, or the preset does not parse |
| 2 | usage error; nothing ran |
| 10 | a sample reached `--fail-on`: its verdict is fail (or warn with `--fail-on warn`) |
| 11 | a sample is incomplete (a rule's file, row or column is missing) and none reached `--fail-on` |

`qc` uses the same table: 10 and 11 only with `--fail-on`, and 1 also for a failed collector. The tool verbs and `collect` exit 0, 1 or 2; the Picard face keeps 0 and 1.

`.cargo/config.toml` lets the Python cdylib link under `cargo test` on macOS.

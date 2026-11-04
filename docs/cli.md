# Command line and Python reference

Every subcommand, its options and the Picard-syntax forms. The short version is in the [README](../README.md).

## Build and test

Needs Rust 1.85 or newer (edition 2024).

```bash
cargo build --release            # target/release/bammetrics
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
target/release/bammetrics insert-size tests/data/insert_size.bam            # Picard-format file on stdout
target/release/bammetrics insert-size in.bam -O in.insert_size_metrics.txt -M 0.05
target/release/bammetrics insert-size --input in.bam --output out.txt
target/release/bammetrics alignment-summary in.bam -R ref.fa -O in.alignment_summary_metrics.txt
target/release/bammetrics quality-yield in.bam -O in.quality_yield_metrics.txt
target/release/bammetrics hs-metrics in.bam -O in.hs_metrics.txt -R ref.fa \
    --bait-intervals baits.interval_list --target-intervals targets.interval_list
# Picard syntax, a drop-in for the jar:
target/release/bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt H=hist.pdf M=0.5
target/release/bammetrics CollectInsertSizeMetrics --INPUT in.bam --OUTPUT out.txt \
    --METRIC_ACCUMULATION_LEVEL SAMPLE --METRIC_ACCUMULATION_LEVEL READ_GROUP
target/release/bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt R=ref.fa
target/release/bammetrics CollectQualityYieldMetrics I=in.bam O=out.txt
target/release/bammetrics CollectHsMetrics I=in.bam O=out.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
```

Picard options that cannot change a metric (`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`,
`VERBOSITY`, `H`, `CHART`, `R`, ...) are accepted, dropped, and noted on stderr.
`ASSUME_SORTED=false` makes a registry tool exit 1 on a BAM whose header is not
coordinate sorted, as Picard does; `STOP_AFTER` works in CollectAlignmentSummaryMetrics,
CollectQualityYieldMetrics, the three by-cycle tools and CollectDuplicateMetrics.
Any other option exits 1 naming the option. `CollectInsertSizeMetrics` and
`CollectAlignmentSummaryMetrics` take `METRIC_ACCUMULATION_LEVEL` (`ALL_READS`, `SAMPLE`,
`LIBRARY`, `READ_GROUP`; `null` clears the default); `CollectHsMetrics` accepts `ALL_READS` only. Tools in the registry
(`crates/core/src/registry.rs`) also take Picard 3's `--NAME value` and `--NAME=value` forms.

`alignment-summary` takes `-R` (REFERENCE_SEQUENCE; a FASTA with a `.fai`, read one
contig at a time), `--max-insert-size`, `--expected-pair-orientation`, `--adapter-sequence`,
`--collect-alignment-information`, `--is-bisulfite-sequenced`, `--level` and `--stop-after`. The two list options add to Picard's default
(`[FR]`, six Illumina adapters) and `null` clears it, as Picard's own parser does.
Without `-R` Picard still writes every row and counts no mismatches; so does
`alignment-summary`. With `-R` the BAM `@SQ` lines must match the `.fai` (names, lengths,
order) and the records must be coordinate sorted, or it exits 1 as Picard does.

`quality-yield` takes `--use-original-qualities <bool>` (USE_ORIGINAL_QUALITIES, `OQ=`; default
true, so an `OQ` tag replaces `QUAL`), `--include-secondary-alignments` and
`--include-supplemental-alignments`. It reads no reference, so `R=` is dropped with a note. It
needs no index and no sort order. Python: `bammetrics.quality_yield(path, ...)` returns one
`QualityYieldMetrics`.

`hs-metrics` takes `--bait-intervals` and `--target-intervals` (repeatable; Picard
`BAIT_INTERVALS`/`BI`, `TARGET_INTERVALS`/`TI`), `--bait-set-name` (`N`), `-R` (needed for
`AT_DROPOUT` and `GC_DROPOUT`), `--minimum-mapping-quality` (`MQ`, 20),
`--minimum-base-quality` (`Q`, 20), `--clip-overlapping-reads <bool>` (true),
`--include-indels`, `--coverage-cap` (`covMax`, 200), `--near-distance` (250) and
`--sample-size` (10000). An interval file is an interval list when its first line starts
with `@` (the `@SQ` lines must equal the BAM's), otherwise BED (0-based, read against the
BAM header). Picard 3.5.0 refuses BED files; bammetrics reads them as `BedToIntervalList`
would convert them. The run time has a fixed part of about 0.7 s: the `HET_SNP_SENSITIVITY`
simulation draws 10000 x 201 qualities from the same `java.util.Random(51)` stream as
Picard. `--per-base-coverage FILE` (`PER_BASE_COVERAGE`) also writes Picard's per-base file: a
`chrom pos target coverage` header, then one line per base of every merged target in
interval order, with the 1-based position, the target name (names of merged targets joined
by `|`) and the high-quality depth. A name ending in `.gz` is gzip-compressed (BGZF).
Python: `bammetrics.hs_metrics(path, bait_intervals, target_intervals,
reference=None, ..., per_base_coverage=None)` returns one `HsMetrics`; every column is an attribute in lower case,
`to_dict()` gives all of them, and `high_quality_coverage` and `unfiltered_base_quality`
are the two histograms as dicts.

`collect` computes several metrics in one pass over the BAM and writes one Picard file
per metric: `PREFIX.quality_yield_metrics`, `.alignment_summary_metrics`,
`.insert_size_metrics` (the names Picard `CollectMultipleMetrics` uses) and `.hs_metrics`.
`--hs BAITS TARGETS` adds one `CollectHsMetrics` run and may be repeated; with more than
one, each file is `PREFIX.<BAIT_SET>.hs_metrics`, where BAIT_SET is the bait file name up
to its first dot, as in Picard. `--hs-per-base` also writes each run's per-base coverage
to the same name with `.per_base_coverage` in place of `.hs_metrics`. `-R` goes to
alignment-summary and hs. `--metric` picks a
subset (comma separated; default all three, plus hs when `--hs` is given). Every other
option is Picard's default. A metric that fails (for example an unsorted BAM with `-R`)
writes no file and exits 1 after the others are written.

```bash
target/release/bammetrics collect in.bam -O out/sample -R ref.fa \
    --hs panel.bed panel.bed --hs baits.interval_list panel.bed
```

Every subcommand takes `--threads N`: the BGZF decompression workers (default: the
available cores, at most 8). Records are decoded and counted on one thread in file order,
so results do not depend on N. In `collect` each metric also runs on a thread of its own,
fed the same record batches in file order. Python: `bammetrics.collect(path,
metrics=None, reference=None, hs=None, output_prefix=None, threads=None,
hs_per_base=False)` returns
`{name: result}`, each result as the single-metric function returns it (`hs_metrics` is a
list), and writes the files when `output_prefix` is given. Timings are in
[docs/benchmarks.md](benchmarks.md).

`insert-size` takes Picard's options: `-M` (MINIMUM_PCT), `--deviations`, `-W`
(HISTOGRAM_WIDTH), `--min-histogram-width` and `--include-duplicates`. Like Picard,
it writes no file when no orientation is kept (single-end input, for example).

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
in Python arrives with the Python mirror (phase 11 of docs/05).

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

`qc` uses the same table: 10 and 11 only with `--fail-on`, and 1 also for a failed collector. The Picard face and the five older subcommands keep exit codes 0 and 1.

`.cargo/config.toml` lets the Python cdylib link under `cargo test` on macOS.

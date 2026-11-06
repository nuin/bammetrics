# Command line reference

Every verb, its options, the files it writes, its exit codes and the Picard-syntax
forms, grouped as docs/05 section 2.1 groups them: run, prepare, post and info. The short
version is in the [README](../README.md); presets and `verdict.json` are in
[docs/verdict.md](verdict.md), the Python module in [docs/python.md](python.md).

## Build and test

Needs Rust 1.88 or newer (edition 2024).

```bash
cargo build --release            # target/release/bammetrics
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
make budgets                     # the size and speed budgets of docs/05 section 4.8
```

`.cargo/config.toml` lets the Python cdylib link under `cargo test` on macOS.

## Shape

```
bammetrics <verb> [INPUT] [options]       kebab-case verbs
bammetrics <PicardClass> KEY=VALUE ...    the Picard drop-in, unchanged
```

`run` in `crates/cli/src/lib.rs` dispatches in this order: a first argument that names a Picard class in the
registry runs the Picard face, untouched; then a grammar verb (a tool verb, `collect`,
`tools`, `presets`, `help`, `completions`); then the clap verbs `qc`, `plan`, `check` and
`report`. Class names are CamelCase and verbs are kebab case, so the two cannot collide.

| Group | Verbs |
|---|---|
| Run | `qc` (the collectors a preset needs, one pass, files, verdict, page), `collect` (the same pass, no verdict), one verb per registry tool |
| Prepare | `bed-to-interval-list`; `mark-duplicates` is reserved for parity plan P9 and not built |
| Post | `accumulate-quality-yield`, `accumulate-variant-calling`, `report`, `check` |
| Info | `tools`, `presets`, `plan`, `help`, `completions` |

```bash
bammetrics insert-size in.bam                               # Picard text on stdout
bammetrics insert-size in.bam -o in.insert_size_metrics --minimum-pct 0.1
bammetrics insert-size in.bam --MINIMUM_PCT 0.1 --format json | jq '.metrics[0].MEDIAN_INSERT_SIZE'
samtools view -b in.bam chr20 | bammetrics quality-yield - -o chr20.quality_yield_metrics
bammetrics gc-bias in.bam -R ref.fa -o out/NA12878          # out/NA12878.gc_bias.*_metrics
bammetrics hs-metrics in.bam --panel panel.bed -R ref.fa -o out/NA12878.hs_metrics
bammetrics qc in.bam --preset clinical-panel --panel panel.bed -R ref.fa -o out/NA12878
bammetrics report out/ --preset clinical-panel
bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt M=0.5
```

## Shared options

Each option has one row in the registry option table: its Picard names, value type,
whether it repeats or is required, and its help text. The row gives the kebab flag, the
Picard key, the `--help` line, the Python keyword argument and the `bammetrics.pyi` stub.
`bammetrics help VERB` (or `bammetrics VERB --help`, or `bammetrics PicardClass --help`)
lists every row of a verb.

Spellings. `--minimum-pct 0.1`, `--minimum-pct=0.1`, `--MINIMUM_PCT 0.1`, `MINIMUM_PCT=0.1`,
`minimum_pct=0.1`, `-M 0.1`, `-M0.1` and `M=0.1` are one option. A name is upper-cased with
`-` read as `_`, then matched against the row's Picard names, and the registry parser gets
`MINIMUM_PCT=0.1`. A one-letter Picard name is a short flag (`-M`, `-W`, `-R`), and its value
may be attached as clap read it in 0.12 (`-Oout.txt`, `-Rref.fa`, `-@4`). A longer Picard
name is a `--kebab` alias (`--level` for `LEVEL`). A bare word holding `=` is an option only
when the part before `=` has no path character, so `dir/a=b.bam` stays an input.

Booleans take `--x` (true), `--no-x`, `--x=false`, `X=false` and the space form `--x false`
of the 0.12 subcommands. Only a literal `true` or `false` (any case) after a boolean is read
as its value: in `--include-duplicates in.bam`, `in.bam` is the input. Picard's plumbing
booleans (`--create-index`, `--create-md5-file`, ...) are flags too. An input file named
`true` or `false` right after a boolean is read as the value; write `./true`, and the
missing-input error says so. A word after a bare boolean that is neither (`--x maybe`) is
read as an input, and the error names the boolean.

On every tool verb:

| Option | Meaning | Picard spelling |
|---|---|---|
| INPUT or `-i` | input file; `-` reads standard input | `I`, `INPUT` |
| `-o` (`-O` kept) | file, prefix, or `-` for standard output | `O`, `OUTPUT` (`M` for duplicate-metrics) |
| `-@`, `--threads N` | BGZF decompression workers; default the available cores, at most 8 (`--help` prints the number); 0 counts as 1, as in 0.12 | none |
| `-R` | FASTA with a `.fai` | `R`, `REFERENCE_SEQUENCE` |
| `--level` | accumulation level, repeatable | `LEVEL`, `METRIC_ACCUMULATION_LEVEL` |
| `--stop-after N` | read limit | `STOP_AFTER` |
| `--format` | `picard`, `json`, `tsv`, `jsonl`, `human` (single-file verbs) | none |
| `--sample NAME` | sample name in the json and jsonl forms, default the first `@RG SM`, else the input name without its extension (none for standard input); never written into a Picard file | none |
| `-q`, `-v` | no notes on stderr; the Picard command line and each file written on stderr | `QUIET` |

`--level` and `--stop-after` exist on a verb whose tool has the option; elsewhere they are
unknown options, as the Picard face rejects them. `insert-size` has no `--stop-after`:
Picard's `CollectInsertSizeMetrics` takes `STOP_AFTER`, bammetrics does not implement it
yet, and `quality-yield` has no `--level`, as in Picard. `-R` on a tool that reads no
reference is accepted and noted, as `R=` is on the Picard face. Records are decoded and counted
on one thread in file order, so the metrics do not depend on `--threads`.

Options Picard has that cannot change a metric (`VALIDATION_STRINGENCY`, `TMP_DIR`,
`VERBOSITY`, `H`, `CHART`, `R` on a tool that reads no reference, ...) are accepted,
dropped and noted on stderr. An option no row names is an error, because dropping an option
Picard applies would change the numbers silently.

`qc`, `plan`, `report` and `check` are clap verbs with their own options, listed under each
verb below and by `bammetrics help qc` (or `plan`, `check`, `report`). They share `-o`,
`-R`, `--preset NAME|FILE`, `--sample`, `--threads`, `--format` and `-q` with the meaning
above where a verb has them.

## Output naming

- A single-file verb writes `-o FILE`. Without `-o`, or with `-o -`, the file goes to
  standard output.
- A multi-file verb (`gc-bias`, `multiple-metrics`, `variant-calling`,
  `accumulate-variant-calling`, `genotype-concordance`, `collect`, `qc`) writes under `-o
  PREFIX`. Without `-o` the prefix is `./` plus the input name without its extension, and
  stderr (the summary for `qc`) says so. An `-o` that ends in `/` or names an existing
  directory gets that name as the file part (`-o out/` writes `out/NA12878.*`). `-o -` is
  an error.
- The names after the prefix come from one table, `multiple_metrics::file_extensions`, the
  `CollectMultipleMetrics` program extensions (`.alignment_summary_metrics`,
  `.insert_size_metrics`, `.quality_yield_metrics`, `.gc_bias.detail_metrics`, ...) plus
  `.hs_metrics`, `.duplicate_metrics`, `.wgs_metrics` and `.rna_metrics`. `collect`, `qc`
  and the verbs read it; a test pins it, and another checks the verbs against
  `CollectMultipleMetrics O=PREFIX`.
- Several Hs runs write `PREFIX.<BAIT_SET>.hs_metrics`, `BAIT_SET` being the bait file name
  up to its first dot. Two bait files with one stem are an error.
- Run files of `qc`: `PREFIX.manifest.json`, `PREFIX.verdict.json`, `PREFIX.report.html`,
  and `PREFIX.panel.interval_list` for a `--panel` BED. `report` writes
  `OUT/SAMPLE.verdict.json`, `OUT/SAMPLE.report.html` and `OUT/index.html`.
- Missing directories above an output are created. `qc` and `report` write each file to a
  temporary name and rename it into place.
- An output that resolves to an input file (`-o in.bam`) is a usage error, and nothing is
  written.
- `-` as the input reads the BAM from standard input on the verbs only; the header line
  reads `INPUT=-`. The Picard face reads a file named `-` and writes `O=-` to a file named
  `-`, as listed in KNOWN_DIFFERENCES.md.
- The sample name is never written into a Picard file, which stays byte-compatible. It goes
  into the json forms, the verdict, the page and the manifest.

## Formats

`--format picard` (the default) is the Picard metrics file. `json` is one object: `tool`,
`sample` (`--sample`, else its default), `class`, `metrics` (one object per row, keyed by column) and
`histogram` (`columns`, `rows`, or `null`). `jsonl` is one object per row, `tsv` a header
line and one line per row, `human` an aligned table (a row per line up to 6 columns, else
`COLUMN  value` lines per row). Numbers are JSON numbers, NaN and infinity `null`; the text
forms print Picard's cell text. For `QualityScoreDistribution` and `MeanQualityByCycle`,
whose file has no metrics class, the rows are the histogram. `--format` on a multi-file
verb is a usage error that names the problem.

`qc` and `plan` take `--format text|jsonl`, `check` takes `text|json|jsonl`; each is
described under its verb.

## Exit codes

One table for every verb, `EXIT_CODES` in `check.rs`:

| Code | Meaning | Verbs |
|---|---|---|
| 0 | success: every verdict passes, or none reaches `--fail-on` (`qc` and `report` without `--fail-on`: any verdict) | all |
| 1 | an input could not be read (a missing BAM, FASTA, `.fai` or interval file), a tool failed, the preset does not parse, or a `qc` collector failed | all |
| 2 | usage error, nothing ran (bad arguments; `qc` over an earlier run without `--resume` or `--force`) | all |
| 10 | a verdict reached `--fail-on`: fail, or warn with `--fail-on warn` (`check` defaults to `--fail-on fail`) | `check`, `qc --fail-on` |
| 11 | a verdict is incomplete (a rule's file, row or column is missing) and none reached `--fail-on` | `check`, `qc --fail-on` |

`qc` exits 0 on any verdict unless `--fail-on warn|fail` is set, so a Nextflow task is not
killed by a bad sample, and its summary prints `use --fail-on fail or bammetrics check to
gate`. `report` exits 0 on any verdict. `check` always gates. A `qc` run with a failed
collector exits 1 after it writes the other files, the verdict and the page. The tool verbs
and `collect` exit 0, 1 or 2; the Picard face keeps 0 and 1.

## Errors

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
$ bammetrics inser-size in.bam
bammetrics: unknown verb or Picard class "inser-size"
  did you mean insert-size?
  bammetrics tools lists them
```

A usage error that clap finds on `qc`, `plan`, `check` or `report` uses the same form, with
the usage line as a hint. An error that stops `qc` or `plan` prints `bammetrics qc: WHAT`
(or `bammetrics plan:`) on stderr and, with `--format jsonl`, an `error` event (`step`
null) and a `done` event with only `exit` on stdout.

## Run verbs

### Tool verbs

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

`multiple-metrics --program X` adds to Picard's default five programs, as `PROGRAM=` does;
`--program null` first clears them. The accumulate verbs and `bed-to-interval-list` are
listed here because the registry runs them; they are described under post and prepare.

### The 0.12 subcommands

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

### collect

Several metrics on one read of the BAM, one Picard file each under `-o PREFIX` (default
`./<bam stem>`; it was a usage error before 0.13.0): `PREFIX.quality_yield_metrics`,
`.alignment_summary_metrics`, `.insert_size_metrics` and `.hs_metrics`. `--hs BAITS
TARGETS` adds one `CollectHsMetrics` run and may be repeated; with more than one, each file
is `PREFIX.<BAIT_SET>.hs_metrics`. `--hs-per-base` also writes each run's per-base coverage,
named as its HS file with `.per_base_coverage` in place of `.hs_metrics`; without `--hs` it
is a usage error. `-R` goes to alignment summary and hs. `-o -` is an error, as on the
other multi-file verbs. `--metric` (or `-m`) picks a subset, comma separated, of
`quality-yield`, `alignment-summary`, `insert-size` and `hs` (default the first three, plus
hs with `--hs`); another name exits 2 with the list. Every other option is Picard's
default. A metric that fails (an unsorted BAM with `-R`, for example) writes no file and
exits 1 after the others are written. In `collect` each metric runs on a thread of its own,
fed the same record batches in file order. `collect` takes `I=` and `O=` as well. Timings
are in [docs/benchmarks.md](benchmarks.md).

```bash
bammetrics collect in.bam -o out/sample -R ref.fa \
    --hs panel.bed panel.bed --hs baits.interval_list panel.bed
```

For the by-cycle tools and gc bias on one read, use `multiple-metrics` (see the chaining
examples).

### qc

```bash
bammetrics qc BAM [-o PREFIX] [-R FASTA] [--preset NAME|FILE] [--hs BAITS TARGETS]... [--panel FILE]
              [--ref-flat FILE --strand STRAND [--ribosomal-intervals FILE]] [--sample NAME]
              [--fail-on warn|fail] [--resume | --force] [--threads N] [--format text|jsonl] [-q]
              [--all-plots] [--date ISO]
```

`qc` runs the collectors the preset's rules need (`generic` without `--preset`; it never
picks a clinical preset by itself), on one read of the BAM, then writes the verdict and
the page. `--hs` or `--panel` adds `CollectHsMetrics` and `--ref-flat` adds
`CollectRnaSeqMetrics` whether or not a rule reads them. A rule whose collector cannot run
is listed as not measured with the reason (`needs --hs BAITS TARGETS or --panel FILE`,
`needs -R FASTA`), and the verdict reads INCOMPLETE.

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

Each metrics file is byte for byte what `bammetrics <Collector> INPUT=BAM OUTPUT=FILE ...`
writes for the argv the manifest records, apart from the `Started on` line. The header
command line follows that argv: `qc` passes `REFERENCE_SEQUENCE` only to the collectors that
read `-R` (alignment summary, hs, gc bias, wgs), so the insert size and duplicate headers
lack the `REFERENCE_SEQUENCE=` that `bammetrics CollectInsertSizeMetrics ... R=ref.fa`
prints; the data lines are the same. The sample name (`--sample`, else the first `@RG SM`,
else the BAM name) goes only into the verdict, the page and the manifest. No `-o` means
`./<BAM name without .bam>`, and the summary says so; a directory `-o` gets the BAM name as
the file part, and the summary says so. Every file is written to a temporary name opened
with `create_new` and renamed into place. A file a collector plans and does not write this
time (it failed, or insert size dropped every orientation and Picard writes no file) is
removed, and so is a metrics file of an earlier run's collector that the new options drop
(`--panel` given before, not now), so the prefix holds what the manifest lists; the summary
counts both as stale removed. With several `--hs` runs the verdict reads the HsMetrics file
whose path sorts first (`PREFIX.<bait stem>.hs_metrics`); every file is written and listed
in the manifest. Two bait files with one stem are an error.

Stages, in order: prepare (the BAM header, the input fingerprint, the preset, every file an
option names; then, when `--panel` names a BED, the prepare step of `bed-to-interval-list`
below; no BAM record is read), pass (every one-pass collector on one read; a collector that
fails to build, fails on a record or panics fails alone, and a BAM read error fails the
pass), post (the verdict over the files), report (`PREFIX.verdict.json` and
`PREFIX.report.html`, the bytes `bammetrics report` writes for the same files and sample).
`CollectWgsMetrics` walks reference loci and reads the BAM on its own. The verdict and the
page add no BAM read: `make budgets` checks that `qc` reports one BAM read with the report
included and that `report` and `check` run on its output once the BAM is deleted.

`PREFIX.manifest.json` (schema `bammetrics.manifest/1`) holds per step its stage, tool,
Picard argv and the argv hash, the tool version, the preset hash, the input path, size,
mtime and a hash of its first and last 64 KiB, the same three for every further file the
argv names (`files`: the FASTA and its `.fai`, the bait and target intervals, the refFlat,
the ribosomal intervals), the outputs with tool and kind, and the status (`ok` or
`failed`, with the error). `--resume` skips a collector when every one of those fields
matches and its files exist: one lost file reruns that collector alone, in a pass of its
own, and an interval file edited in place reruns the collectors that read it. A BAM (or
FASTA) larger than 128 KiB edited in the middle with its size and mtime kept (`touch -r`,
`cp -p` over it) is not detected; use `--force`. A manifest that does not parse (a killed
run) is reported on stderr and every step runs. `--force` reruns every step. Without
either, `qc` refuses a prefix that holds an earlier run's manifest (exit 2: `out/s.manifest.json
holds an earlier run: --resume skips its finished steps, --force reruns them`), and names a
manifest it cannot read with the reason (`out/s.manifest.json is not a whole manifest:
--resume or --force runs every step`). The verdict and the page are always rewritten.

The summary on stdout is at most 12 lines: the header (sample, band, preset, collectors,
BAM reads, records the pass read: `4 collectors, 1 BAM read (41.8M records)`), every error,
fail, warn and not-measured line (collapsed into `+N more, see PAGE`), the pass count, the
files written (kept by `--resume`, replaced, stale removed), and the gate hint or the exit
reason. It has no elapsed time, so the same inputs print the same stdout. Progress
(`reading 12.4M reads  1.9M/s  62%`) goes to stderr only on a terminal, at most 10 times a
second; `-q` turns it off. The output has no color. `--format jsonl` prints
`start`, `progress` (every million records), `file` (with `written` false for a file
`--resume` kept), `error`, `check` and `done` events on stdout instead. A closed stdout
(`| head -1`) drops the rest of the output; the run still writes its files and keeps its
exit code. `--all-plots` and `--date` act as on `report`.

`qc` exits 0 on any verdict unless `--fail-on` is given; then 10 when the verdict reaches it
and 11 when the verdict is incomplete, as `check` does. A failed collector exits 1 after the
other files, the verdict and the page are written.

## Prepare verbs

### bed-to-interval-list

```bash
bammetrics bed-to-interval-list panel.bed --from-bam sample.bam -o panel.interval_list
bammetrics bed-to-interval-list panel.bed --dict ref.dict -o -
bammetrics CollectHsMetrics I=sample.bam O=hs.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
```

The verb writes the interval list `BedToIntervalList I=BED SD=FILE O=FILE` writes, byte for
byte: `@HD VN:1.6 SO:coordinate`, the dictionary's `@SQ` lines with their fields in its
order, then one line per BED interval (1-based, closed, `+` or `-`, the name or `.`).
Without `-o`, or with `-o -`, the same bytes go to standard output. `-` as the BED reads
standard input.

- `--from-bam FILE`, `--dict FILE` and `--sd FILE` (`--sequence-dictionary`,
  `SEQUENCE_DICTIONARY=`, `SD=`) are spellings of one option. The file type comes from its
  name, as the Picard face reads it: a `.bam` name is read as a BAM header, a FASTA name
  (`.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`) finds the `.dict` beside it (its `.fai`
  without one), and any other name is read as text with `@SQ` lines (`.dict`, SAM, interval list).
  So `--from-bam x.dict` and `--dict x.bam` both work.
- Two of those spellings on one command line exit 2, from the second one:
  `bammetrics bed-to-interval-list: --from-bam and --dict both set SEQUENCE_DICTIONARY (from
  --dict)`. One spelling given twice keeps the last value, as every option does. A
  dictionary that does not exist exits 1 from the spelling given (`nothere.bam not found
  (from --from-bam)`), and `-o` naming the dictionary or the BED exits 2. A typo of
  `--from-bam` or `--dict` gets a `did you mean` line, as a typo of the long name does.
- A BED line the dictionary rejects (a contig it lacks, a start or end past the contig, a
  field that is not a number) exits 1 with Picard's message, on the verb and on the Picard
  face; exit 2 stays for command-line mistakes.
- `--sort` (default true), `--unique`, `--drop-missing-contigs` and
  `--keep-length-zero-intervals` are Picard's options.
- `--format json`, `jsonl`, `tsv` and `human` print the intervals as rows (`contig`,
  `start`, `end`, `strand`, `name`), with no `sample`. `--threads` and `--sample` are
  accepted and ignored: no BAM record is read, and an interval list has no sample.

`hs-metrics`, `collect` and `qc` read a BED directly, so the verb is for Picard tools and
other programs that take only interval lists.

In `qc`, `--panel BED` runs the same conversion with the BAM as the dictionary and keeps
`PREFIX.panel.interval_list`, the bytes `bammetrics bed-to-interval-list BED --from-bam
BAM` writes. A panel whose first non-blank line starts with `@` is an interval list
already, and prepare writes nothing for it. A `--panel` BED that is itself named
`PREFIX.panel.interval_list` exits 2 (`is an input; writing it would destroy it (from
--panel)`), as does any planned file that is an input. The hs collector reads the `--panel`
file itself, so its file, header line and `BAIT_SET` are the same with or without the
prepare step. The prepare step is a step of the manifest (stage `prepare`, output kind
`intervals`, the BED under `files`), so `--resume` keeps the interval list, a lost one is
written again, and a later run without the BED removes it as a stale file. It is not
counted as a collector; when it fails (a contig the BAM header lacks), the summary lists it
as an `ERROR` line and `qc` exits 1 (`exit 1: 1 collector(s) and 1 prepare step(s) failed`
when hs fails on the same BED).

### mark-duplicates

Reserved for parity plan P9 (docs/05 section 2.5, example 9). It is not built:
`bammetrics mark-duplicates` exits 2 as an unknown verb. `duplicate-metrics marked.bam`
works on a BAM marked by another tool.

## Post verbs

### accumulate-quality-yield and accumulate-variant-calling

```bash
bammetrics accumulate-quality-yield out/a.qy out/b.qy -o out/sample.qy
bammetrics AccumulateQualityYieldMetrics I=out/a.qy I=out/b.qy O=out/sample.qy
```

Files in, file out; no BAM is read. The two forms write the same metrics.
`accumulate-variant-calling` takes the prefixes `variant-calling` wrote and writes the two
files under `-o PREFIX`, inverting the printed ratios into counts as Picard does.

### report

```bash
bammetrics report DIR... [--preset NAME|FILE] [-o OUT] [--all-plots] [--date ISO]
```

`report` works on any directory of Picard output and never reads a BAM. A metrics file is
any file whose text holds a `## METRICS CLASS` block (or a `## HISTOGRAM` section under a
Picard header), whatever its name, at any depth under `DIR`. Files that look like metrics
and do not parse are listed on stderr and left out. Samples: the files in one directory are
grouped by the file name up to the first dot, so `NA12878.insert_size_metrics` and
`NA12878.hs_metrics` form the sample `NA12878` (the prefix `CollectMultipleMetrics`
writes). When no metrics class repeats across the groups of a directory, the directory
holds one sample, named by the prefix when two or more files share it, else by the
directory: `runs/NA12878/insert_size.txt` and `runs/NA12878/hs.txt` form the sample
`NA12878`. Two samples with one name is an error.

`report` evaluates the preset (`generic` without `--preset`; a shipped name or a TOML
file) and writes `OUT/SAMPLE.verdict.json`, `OUT/SAMPLE.report.html` and `OUT/index.html`,
`OUT` being `DIR` when one directory is given. Each file is written whole and renamed into
place. The summary on stdout is one screen: a header line, one line per sample with its
band and counts, the fail and warn reasons (three per sample, then `+N more`), and the
files written. The index lists the samples worst first (fail, incomplete, warn, pass, then
by fail count, then name) with one cell per rule. `--all-plots` draws every plot family the
files have data for; `--date` fixes the date printed on the pages, so the same inputs give
the same bytes. `report` exits 0 on any verdict: a Nextflow task is not killed by a bad
sample.

The page and the index are self-contained HTML with no script and no network request: a
`<meta>` CSP of `default-src 'none'; style-src 'unsafe-inline'`, the CSS inlined, plots as
inline SVG with a text alternative. They follow the light or dark system theme, print on a
light palette, and pass axe-core 4.10.2 with no violation at 320, 360 and 1280 px in light,
dark and print (`make a11y`, recorded in [docs/benchmarks.md](benchmarks.md)). Every text color pair meets WCAG AA (4.5:1) in light, dark and print, checked
by a unit test over the tokens of `report.css`. Size and render time are budgets of docs/05
section 4.8, measured by `make budgets` ([docs/benchmarks.md](benchmarks.md)).

### check

```bash
bammetrics check DIR|VERDICT.json... [--preset NAME|FILE] [--fail-on warn|fail] [--format text|json|jsonl]
```

`check` is the gate. A directory is discovered and evaluated as `report` does; a
`verdict.json` file is read back as is. `--format text` prints a band line per sample and
every check under it; `json` prints a JSON array of verdicts, each as `verdict.json`;
`jsonl` prints one object per line, a `check` event per check and a `verdict` event per
sample. It writes no file. Exit codes: the table above, with `--fail-on fail` by default.

## Info verbs

- `bammetrics tools`: every verb and the Picard class it runs, then the other verbs.
- `bammetrics presets`: the built-in verdict presets with version, the first 6 characters
  of the content hash, rule count and note (`clinical-panel   v2  e62755  14 rules  reviewed
  defaults`). Every rule is listed in [docs/verdict.md](verdict.md).
- `bammetrics help VERB|PicardClass`: the verb's options in every spelling; `help qc` (or
  `plan`, `check`, `report`) prints clap's help.
- `bammetrics completions bash|zsh|fish|elvish|powershell`: a completion script on stdout
  (`bammetrics completions zsh > ~/.zfunc/_bammetrics`).

### plan

```bash
bammetrics plan BAM [the qc options] [--format text|jsonl]
```

`plan` is the dry run of `qc`. It prints the stages, the BAM reads (`2 BAM reads: 1 for the
pass, 1 own read`), each collector with its files and rule count, which collectors share
the read and why one does not, the classes no collector can write, and which steps
`--resume` would skip. It reads the BAM header and the first and last 64 KiB of the BAM and
of every further input file, and writes nothing; a file that is not a BAM exits 1 there.
`--format jsonl` gives one `plan` event, a `step` event per collector (with `argv`,
`outputs`, `rules`, `read`, `action`, `resumable`) and an `unmet` event per class. The
prepare step of a `--panel` BED is listed under `prepare` in the text plan and as a `step`
event with `read` `header`, before the collectors.

## The Picard face

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
`ASSUME_SORTED=false` exits
1 on a BAM whose header is not coordinate sorted, as Picard does.
`METRIC_ACCUMULATION_LEVEL` (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`; `null`
clears the default) works where Picard has it. Where the output differs from Picard on
purpose is listed in KNOWN_DIFFERENCES.md.

## Chaining examples

The worked examples of docs/05 section 2.5 that work today, written against the committed
fixtures so each runs from the repository root as shown (`tests/data/insert_size.bam` has
no `@RG SM`, so the sample is the BAM name). Output is what 0.16.0 prints; long lines are
cut.

**1 and 2. One tool, as Picard text or machine-readable.**

```
$ bammetrics insert-size tests/data/insert_size.bam --format json | jq '.metrics[0].MEDIAN_INSERT_SIZE'
310
$ bammetrics insert-size tests/data/insert_size.bam --MINIMUM_PCT 0.1 --format tsv | cut -f1-4
MEDIAN_INSERT_SIZE	MODE_INSERT_SIZE	MEDIAN_ABSOLUTE_DEVIATION	MIN_INSERT_SIZE
310	310	10	300
```

**3. Quick look, no preset.**

```
$ bammetrics qc tests/data/insert_size.bam -R tests/data/ref_100k.fa -o out/quick/S1
insert_size  PASS  preset generic v1 (6dd89e)  4 collectors, 1 BAM read (10 records)
PASS          6 checks
wrote 7 files to out/quick/S1.*, report out/quick/S1.report.html
```

**4. Clinical panel, one sample.**

```
$ bammetrics qc tests/data/insert_size.bam --preset clinical-panel --panel tests/data/targets.bed \
      -R tests/data/ref_100k.fa -o out/NA12878
insert_size  FAIL  preset clinical-panel v2 (e62755)  4 collectors, 1 BAM read (10 records)
FAIL          hs.mean-target-coverage   Mean target coverage is 0.377778, below the fail limit of 100
FAIL          hs.target-20x             Target bases at 20x is 0, below the fail limit of 0.95
FAIL          hs.target-100x            Target bases at 100x is 0, below the fail limit of 0.6
FAIL          hs.zero-coverage-targets  Targets with zero coverage is 0.2, above the fail limit of 0.03
FAIL          hs.gc-dropout             GC dropout is 20.588235, above the fail limit of 10
NOT MEASURED  hs.evenness               FOLD_80_BASE_PENALTY is ?
PASS          8 checks
wrote 8 files to out/NA12878.*, report out/NA12878.report.html
use --fail-on fail or bammetrics check to gate
```

**5. Vendor baits as a second Hs run.** Two bait files with one stem are an error, so the
second run's baits get their own name.

```
$ cp tests/data/targets.interval_list out/vendor.interval_list
$ bammetrics qc tests/data/insert_size.bam --preset clinical-panel --panel tests/data/targets.bed \
      --hs out/vendor.interval_list tests/data/targets.bed -R tests/data/ref_100k.fa -o out/two/NA12878 \
      > /dev/null
$ ls out/two | grep hs_metrics
NA12878.targets.hs_metrics
NA12878.vendor.hs_metrics
```

The verdict reads the file whose path sorts first (`targets`).

**6. Cohort: one `qc` per sample, then one `report`.**

```
$ mkdir -p bams && cp tests/data/insert_size.bam bams/NA1.bam && cp tests/data/insert_size.bam bams/NA2.bam
$ for bam in bams/*.bam; do
    bammetrics qc "$bam" --preset clinical-panel --panel tests/data/targets.bed -R tests/data/ref_100k.fa \
        -o runs/2026-11-02/$(basename "$bam" .bam) -q > /dev/null
  done
$ bammetrics report runs/2026-11-02/ --preset clinical-panel
2 samples, preset clinical-panel v2 (e62755), bammetrics 0.16.0
FAIL        NA1  5 fail, 8 pass, 1 not measured
  FAIL  Mean target coverage is 0.377778, below the fail limit of 100
  ...
wrote 5 files to runs/2026-11-02/, index runs/2026-11-02/index.html
the verdict does not change the exit code; bammetrics check gates on it
```

There is no sample sheet runner (decision 6 of docs/05 section 8): cohorts loop in the
shell, make, Snakemake or Nextflow.

**7. Chain: BedToIntervalList to CollectHsMetrics.**

```bash
bammetrics bed-to-interval-list tests/data/targets.bed --from-bam tests/data/insert_size.bam \
    -o out/panel.interval_list
bammetrics hs-metrics tests/data/insert_size.bam -R tests/data/ref_100k.fa \
    --baits out/panel.interval_list --targets out/panel.interval_list -o out/hs.txt
# the same step through process substitution: with no -o, the first verb writes to stdout
bammetrics hs-metrics tests/data/insert_size.bam -R tests/data/ref_100k.fa \
    --baits <(bammetrics bed-to-interval-list tests/data/targets.bed --from-bam tests/data/insert_size.bam) \
    --targets tests/data/targets.bed -o out/hs2.txt
```

**8. Chain: scattered quality yield to Accumulate.**

```bash
mkdir -p shards && cp tests/data/insert_size.bam shards/a.bam && cp tests/data/insert_size.bam shards/b.bam
for bam in shards/*.bam; do
    bammetrics quality-yield "$bam" -o "out/$(basename "$bam" .bam).qy" &
done; wait
bammetrics accumulate-quality-yield out/*.qy -o out/sample.qy
```

`AccumulateQualityYieldMetrics I=out/a.qy I=out/b.qy O=...` writes the same metrics. GNU
`parallel` works in place of the loop.

**9. MarkDuplicates to CollectDuplicateMetrics.** Not built: `mark-duplicates` waits for
parity plan P9. `bammetrics duplicate-metrics tests/data/insert_size.bam -o out/s.dup`
works on a BAM whose duplicates another tool marked.

**10. One read of standard input, several collectors.**

```
$ samtools view -b tests/data/insert_size.bam | bammetrics multiple-metrics - -R tests/data/ref_100k.fa \
      -o out/mm --program CollectAlignmentSummaryMetrics --program CollectInsertSizeMetrics \
      --program QualityScoreDistribution
$ ls out | grep '^mm'
mm.alignment_summary_metrics
mm.base_distribution_by_cycle_metrics
mm.insert_size_metrics
mm.quality_by_cycle_metrics
mm.quality_distribution_metrics
$ samtools view -b tests/data/insert_size.bam | bammetrics collect - -R tests/data/ref_100k.fa \
      -o out/stdin -m quality-yield,alignment-summary,insert-size
```

`--program` adds to the default five, so five files come out. docs/05 writes this example
with `collect -m ...,gc-bias,mean-quality-by-cycle,...`; `collect -m` knows only
`quality-yield`, `alignment-summary`, `insert-size` and `hs` and exits 2 on another name.

**11. Report and index into another directory.**

```
$ bammetrics report runs/*/ --preset clinical-panel -o site/ > /dev/null; ls site
index.html
NA1.report.html
NA1.verdict.json
NA2.report.html
NA2.verdict.json
```

**12. Errors and the CI gate.**

```
$ bammetrics hs-metrics tests/data/insert_size.bam --baits tests/data/targets.bed -R tests/data/ref_100k.fa
bammetrics hs-metrics: missing --target-intervals (Picard: TARGET_INTERVALS)
  Use --panel FILE when baits and targets are the same file.
  every option: bammetrics help hs-metrics
$ echo $?
2
$ bammetrics qc nope.bam
bammetrics qc: nope.bam not found (from BAM)
$ echo $?
1
$ bammetrics check out/quick --preset clinical-panel --format jsonl > checks.jsonl; echo $?
11
$ bammetrics check runs/2026-11-02 --preset clinical-panel > /dev/null; echo $?
10
```

`out/quick` holds the `generic` run of example 3, so `clinical-panel` finds no HsMetrics
file there and the verdict is incomplete (11). The cohort fails its hs rules (10). The
error names `--target-intervals`, the option's long name; docs/05 writes `--targets`, an
alias.

## Python

```python
import bammetrics as bm

rows = bm.collect_insert_size_metrics("s.bam", minimum_pct=0.1, metric_accumulation_level=["SAMPLE"])
gc_rows = bm.collect_gc_bias_metrics("s.bam", reference_sequence="ref.fa",
                                     output="out/s.gc", summary_output="out/s.gc_summary")
run = bm.qc("s.bam", preset="clinical-panel", panel="panel.bed", reference="ref.fa", output="out/s")
```

One function per registry tool, named after the Picard class in snake case. A keyword is a
row's Picard name in lower case (`minimum_pct`; a short name such as `m` works too), a list
is a repeated option, `None` leaves the option out, `output=` writes the files (into a
directory that must exist) and the function returns the rows either way. The
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
check per rule with its value, bounds, status and reason. `qc`, `check`, `report`,
`tools` and `presets` run the command's code (0.15.0). Every function and type is in
[docs/python.md](python.md).

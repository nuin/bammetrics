# Command line reference

Every verb, its options, the files it writes, its exit codes and the Picard-syntax
forms, grouped as docs/05 section 2.1 groups them: run, prepare, post and info. The short
version is in the [README](../README.md); presets and `verdict.json` are in
[docs/verdict.md](verdict.md), the Python module in [docs/python.md](python.md).

## Build and test

Needs Rust 1.91 or newer (edition 2024; noodles 0.119 declares 1.91). The PDF report
(cargo feature `pdf`, [docs/pdf.md](pdf.md)) needs Rust 1.92 (krilla 0.8.2).

```bash
cargo build --release            # target/release/bammetrics
cargo build --release -p bammetrics-cli --features pdf   # the same binary with report --pdf and qc --report pdf
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
`chain`, `tools`, `presets`, `help`, `completions`, and `plan` when its arguments describe
a chain); then the clap verbs `qc`, `plan`, `check`,
`report`, `convert` and `doctor`. Class names are CamelCase and verbs are kebab case, so the two cannot collide.

| Group | Verbs |
|---|---|
| Run | `qc` (the collectors a preset needs, one pass, files, verdict, page), `collect` (the same pass, no verdict), `chain` (several steps on one BAM in one process, records handed on in memory), one verb per registry tool |
| Prepare | `bed-to-interval-list`, `mark-duplicates`, `build-bam-index`, `fifo-buffer`, `doctor` |
| Post | `accumulate-quality-yield`, `accumulate-variant-calling`, `compare-metrics`, `check-duplicate-marking`, `report`, `check`, `convert` |
| Info | `tools`, `presets`, `plan`, `help`, `completions` |

```bash
bammetrics insert-size in.bam                               # Picard text on stdout
bammetrics insert-size in.bam -o in.insert_size_metrics --minimum-pct 0.1
bammetrics insert-size in.bam --MINIMUM_PCT 0.1 --format json | jq '.metrics[0].MEDIAN_INSERT_SIZE'
samtools view -b in.bam chr20 | bammetrics quality-yield - -o chr20.quality_yield_metrics
bammetrics gc-bias in.bam -R ref.fa -o out/NA12878          # out/NA12878.gc_bias.*_metrics
bammetrics hs-metrics in.bam --panel panel.bed -R ref.fa -o out/NA12878.hs_metrics
bammetrics qc in.bam --preset clinical-panel --panel panel.bed -R ref.fa -o out/NA12878
bammetrics chain in.bam -o out/NA12878 + mark-duplicates + insert-size + write --index
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
| `-o` (`-O` kept) | file, prefix, or `-` for standard output; a verb refuses one that names an input, the `-R` FASTA or its `.fai` (exit 2), the Picard face writes it as Picard does | `O`, `OUTPUT` (`M` for duplicate-metrics) |
| `-@`, `--threads N` | BGZF decompression workers, or CRAM container decoders; default the available cores, at most 8 (`--help` prints the number); 0 counts as 1, as in 0.12 | none |
| `-R` | FASTA with a `.fai` | `R`, `REFERENCE_SEQUENCE` |
| `--level` | accumulation level, repeatable | `LEVEL`, `METRIC_ACCUMULATION_LEVEL` |
| `--stop-after N` | read limit | `STOP_AFTER` |
| `--format` | `picard`, `json`, `tsv`, `jsonl`, `human`, `csv`, `markdown`, or a list (`picard,json`): the first to `-o`, each further one beside it as `FILE.json` | none |
| `--sample NAME` | sample name in the json and jsonl forms, default the first `@RG SM`, else the input name without its extension (none for standard input); never written into a Picard file | none |
| `-q`, `-v` | no notes on stderr; the Picard command line and each file written on stderr | `QUIET` |

`--level` and `--stop-after` exist on a verb whose tool has the option; elsewhere they are
unknown options, as the Picard face rejects them. `insert-size` has no `--stop-after`:
Picard's `CollectInsertSizeMetrics` takes `STOP_AFTER`, bammetrics does not implement it
yet, and `quality-yield` has no `--level`, as in Picard. `-R` on a tool that reads no
reference is accepted and noted, as `R=` is on the Picard face. Records are decoded and counted
on one thread in file order, so the metrics do not depend on `--threads`.

CRAM. Every verb and Picard-face tool that reads a BAM reads a CRAM, and so do `collect`,
`qc` and `plan`. A file that starts with the bytes `CRAM` is a CRAM whatever its name;
anything else is read as a BAM. `-R` (`R=`) names the FASTA the CRAM was written with,
with a `.fai` beside it, and is then used on every tool: no `ignoring R=` note is written for
a CRAM, nor for an input read from a pipe, whose format is known only once it is read. Without `-R` the run fails at the first slice that needs reference
bases, as Picard does (`the records on chr1 need reference bases: a CRAM input needs
REFERENCE_SEQUENCE (R)`, exit 1); `qc` and `plan` refuse the CRAM before they start. A FASTA
whose bases do not match a slice's MD5 fails there (`the MD5 of the reference bases does not
match the CRAM slice`), as one without the CRAM's contigs does (`no contig chr1 in the
reference`). `--threads` is then the number of container decoders. CRAM 3.0 and 3.1 are
read; another version fails (`CRAM 2.1 is not supported (3.0 and 3.1 are)`), and so does a
cut file (`truncated CRAM: the file ends inside a container`). The `.crai` is never opened.

SAM text. Every verb and Picard-face tool that reads a BAM reads SAM text, plain or in
gzip members (`bgzip` or `gzip`), from a file, a pipe or `-`. The first bytes tell: a
gzip member that inflates to `BAM\1` is a BAM, any other gzip file whose text starts with
`@` or with a line of 11 tab-separated fields is SAM, and so is plain text that starts
that way. Each line is read as htsjdk's `SAMLineParser` reads it into the record a BAM
would hold, so a tool writes the same file for the SAM as for the BAM Picard makes of
it. `--validation-stringency` (`VALIDATION_STRINGENCY`) applies to SAM text: `STRICT`
(the default) exits 1 at the first line htsjdk refuses, with htsjdk's message (`Error
parsing text SAM file. MRNM specified but flags indicate unpaired; File x.sam; Line 7`);
`LENIENT` prints `Ignoring SAM validation error` on standard error and reads on;
`SILENT` checks only what cannot be read at all. A BAM or a CRAM is read the same at
every value, and the note `ignoring VALIDATION_STRINGENCY=...` is printed for those. `-` reads a CRAM from standard input too:
`samtools view -C -T ref.fa in.bam | bammetrics quality-yield - -R ref.fa`.

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
  `accumulate-variant-calling`, `genotype-concordance`, `sequencing-artifact`,
  `convert-sequencing-artifact-to-oxo-g`, `collect`, `qc`) writes under `-o PREFIX`.
  Without `-o` the prefix is `./` plus the input name without its extension, and stderr
  (the summary for `qc`) says so. An `-o` that ends in `/` or names an existing directory
  gets that name as the file part (`-o out/` writes `out/NA12878.*`). `-o -` is an error.
- The names after the prefix come from one table, `multiple_metrics::file_extensions`, the
  `CollectMultipleMetrics` program extensions (`.alignment_summary_metrics`,
  `.insert_size_metrics`, `.quality_yield_metrics`, `.gc_bias.detail_metrics`, ...) plus
  `.hs_metrics`, `.duplicate_metrics`, `.wgs_metrics` and `.rna_metrics`. `collect`, `qc`
  and the verbs read it; a test pins it, and another checks the verbs against
  `CollectMultipleMetrics O=PREFIX`. The `CollectSequencingArtifactMetrics` entry is
  `("OUTPUT", "")`: the program takes the prefix itself and adds its five extensions.
- Several Hs runs write `PREFIX.<BAIT_SET>.hs_metrics`, `BAIT_SET` being the bait file name
  up to its first dot. Two bait files with one stem are an error.
- `chain` writes under `-o PREFIX` too, with the same table: a tap on the last stream
  writes the name `qc` writes, a tap in front of a transform `PREFIX.<stream>.<extension>`
  (`PREFIX.in.insert_size_metrics` on the source), `mark-duplicates`
  `PREFIX.marked.duplication_metrics`, `write` without a path `PREFIX.marked.bam`. A
  step's own `-o` replaces the name ([docs/chaining.md](chaining.md) section 2.3).
- Run files of `qc`: `PREFIX.manifest.json`, `PREFIX.verdict.json`, `PREFIX.report.html`,
  and `PREFIX.panel.interval_list` for a `--panel` BED. `report` writes
  `OUT/SAMPLE.verdict.json`, `OUT/SAMPLE.report.html` and `OUT/index.html`.
- Missing directories above an output are created. `qc` and `report` write each file to a
  temporary name and rename it into place.
- An output that resolves to an input file (`-o in.bam`) is a usage error, and nothing is
  written.
- `-` as the input reads the BAM or CRAM from standard input on the verbs only; the header line
  reads `INPUT=-`. The Picard face reads a file named `-` and writes `O=-` to a file named
  `-`, as listed in KNOWN_DIFFERENCES.md. `jumping-library` is the one exception: it reads
  each input twice, so `-`, `/dev/stdin` and a pipe are a usage error (exit 2).
- The sample name is never written into a Picard file, which stays byte-compatible. It goes
  into the json forms, the verdict, the page and the manifest.

## Formats

[docs/formats.md](formats.md) is the reference: the typed table every format is written
from, the JSON schema, the escapes of each text form and the round trip.

`--format picard` (the default) is the Picard metrics file. `json` is one object (schema
`bammetrics.metrics/1`): `tool`, `sample` (`--sample`, else its default), `headers`, `class`,
`columns`, `types` (`int`, `float` or `string` per column, from the glossary), `metrics` (one
object per row, keyed by column), `non_finite` and `histogram` (`key_class`, `columns`,
`types`, `rows`, `non_finite`, or `null`). `jsonl` is one object per row, a histogram row
with `"histogram": true` first. `tsv` is a header line and one line per row, with tab,
newline and backslash in a string cell written `\t`, `\n`, `\\`. `csv` is RFC 4180 with
full doubles. `markdown` is a GFM table per class. `human` is an aligned table (a row per
line up to 6 columns, else `COLUMN  value` lines per row). Numbers are JSON numbers, NaN and
infinity `null` and listed in `non_finite`; `tsv` and `human` print Picard's cell text. For
`QualityScoreDistribution` and `MeanQualityByCycle`, whose file has no metrics class, the
rows of `tsv`, `csv` and `human` are the histogram.

`--format` takes a list (`picard,json,csv`). The first format goes to `-o` or standard
output, each further one beside `-o` as `FILE.json`, `FILE.csv`; a list needs `-o`. A
multi-file verb and `mark-duplicates` take `picard` first and write each further format
beside every metrics file. Another first format there is a usage error that names the
problem.

`qc` and `plan` take `--format text|jsonl`, the event stream on standard output; `qc --also`
writes file formats. `check` takes `text`, `json`, `jsonl`, `junit`, `markdown` and
`multiqc`, as a list too. Each is described under its verb. `chain` takes `--format` with
`picard` first, on the chain and on a step, and `--also json|jsonl|tsv|csv|markdown`.

## Exit codes

One table for every verb, `EXIT_CODES` in `check.rs`:

| Code | Meaning | Verbs |
|---|---|---|
| 0 | success: every verdict passes, or none reaches `--fail-on` (`qc` and `report` without `--fail-on`: any verdict) | all |
| 1 | an input could not be read (a missing BAM, FASTA, `.fai` or interval file), a tool failed, the preset does not parse, or a `qc` collector failed | all |
| 2 | usage error, nothing ran (bad arguments; `qc` over an earlier run without `--resume` or `--force`) | all |
| 10 | a verdict reached `--fail-on`: fail, or warn with `--fail-on warn` (`check` defaults to `--fail-on fail`); for `doctor`, a check | `check`, `qc --fail-on`, `doctor` |
| 11 | a verdict is incomplete (a rule's file, row or column is missing) and none reached `--fail-on` | `check`, `qc --fail-on` |

`qc` exits 0 on any verdict unless `--fail-on warn|fail` is set, so a Nextflow task is not
killed by a bad sample, and its summary prints `use --fail-on fail or bammetrics check to
gate`. `report` exits 0 on any verdict. `check` always gates. A `qc` run with a failed
collector exits 1 after it writes the other files, the verdict and the page. The tool verbs
and `collect` exit 0, 1 or 2; the Picard face keeps 0 and 1. `chain` exits 0, 1 or 2: 2
for a usage error in any step, found before the first read, and 1 when a step failed, after
the other steps wrote their files. `compare-metrics` and
`CompareMetrics` exit 1 when the two files are NOT equal, as Picard does;
`check-duplicate-marking` and `CheckDuplicateMarking` exit 1 when a record's duplicate flag
differs from the first record of its queryname.

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
class without `Collect` and `Metrics`, in kebab case; the first twelve and
`compare-metrics` are the names of docs/05 section 2.1. `bammetrics tools` prints this
table.

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
| `compare-metrics` | CompareMetrics (two metrics files as positionals) | the comparison text |
| `create-sequence-dictionary` | CreateSequenceDictionary (input is a FASTA) | one `.dict` |
| `sequencing-artifact` | CollectSequencingArtifactMetrics (needs `-R`) | `PREFIX.{pre_adapter,bait_bias}_{detail,summary}_metrics`, `PREFIX.error_summary_metrics` |
| `convert-sequencing-artifact-to-oxo-g` | ConvertSequencingArtifactToOxoG (input is a `sequencing-artifact` prefix) | `PREFIX.oxog_metrics` |
| `oxo-g` | CollectOxoGMetrics (needs `-R`) | one |
| `coverage-gaps` | CollectCoverageGaps (not a Picard tool, [docs/06](06_clinical_checks.md); `--targets` or `--panel`) | one (plus the `--gaps-file` BED) |
| `sex-check` | CollectSexCheckMetrics (not a Picard tool, docs/06; `--sex male|female`, `--targets` for a panel) | one |
| `mt-copy-number` | CollectMtCopyNumberMetrics (not a Picard tool, docs/06 section 6; `--targets` for an exome, `--exclude` for NUMTs, `--mt-contig`, `--autosomes`) | one |
| `adapter-content` | CollectAdapterContentMetrics (not a Picard tool, docs/06 section 7; FastQC's Adapter Content on a BAM; `--adapter-sequence NAME=SEQ`, `--adapter-kmer`) | one |
| `bam-index-stats` | BamIndexStats (reads the `.bai` or `.csi` beside the BAM) | Picard's text on standard output |
| `check-duplicate-marking` | CheckDuplicateMarking | the bad querynames, one line per record |
| `jumping-library` | CollectJumpingLibraryMetrics (one or more coordinate-sorted BAMs as positionals, read in turn; standard input refused) | one metrics file |
| `build-bam-index` | BuildBamIndex (a coordinate-sorted BAM; `-o` names the index, default beside the BAM) | one `.bai`, no metrics file |
| `fifo-buffer` | FifoBuffer (standard input to standard output; no input argument, no `-o`) | no file: the bytes it read, on standard output |
| `view-sam` | ViewSam (SAM, BAM or CRAM in) | SAM text on standard output, or in `-o` |
| `sort-sam` | SortSam (`--sort-order` queryname, coordinate or duplicate) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `exon-coverage` | CollectExonCoverageMetrics (not a Picard tool, [exon_coverage.md](exon_coverage.md); `--annotation` GTF, GFF3, BED12 or BED) | one (plus the `--summary-output` gene and transcript file) |
| `duplicate-by-tile` | CollectDuplicateByTileMetrics (not a Picard tool, docs/06 section 8; `--read-name-regex`, `--outlier-sd`) | one |
| `telomere` | CollectTelomereMetrics (not a Picard tool, [telomere.md](telomere.md); `-K` the repeat threshold, `--repeat`, `--count-duplicates`) | one |
| `index-depth` | CollectIndexDepthMetrics (not a Picard tool, [docs/06](06_clinical_checks.md) section 10; a BAM with its `.bai` or `.csi`, never read whole; a CRAM refused) | one |
| `add-comments-to-bam` | AddCommentsToBam (a BAM; `--comment` repeated) | `-o` BAM, no metrics file |
| `replace-sam-header` | ReplaceSamHeader (a SAM, BAM or CRAM; `--header` names the file whose header replaces the input's) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `add-o-a-tag` | AddOATag (a SAM, BAM or CRAM; `-L` an interval list) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `clean-sam` | CleanSam (a SAM, BAM or CRAM; soft-clips alignments past the end of their reference, MAPQ 0 on unmapped records) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `set-nm-md-and-uq-tags` | SetNmMdAndUqTags (a coordinate-sorted SAM, BAM or CRAM; `-R` required; `--is-bisulfite-sequence`, `--set-only-uq`) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `set-nm-and-uq-tags` | SetNmAndUqTags (Picard's deprecated name of SetNmMdAndUqTags: the same options and output) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `add-or-replace-read-groups` | AddOrReplaceReadGroups (a SAM, BAM or CRAM; `--rglb`, `--rgpl`, `--rgpu`, `--rgsm` required; `--sort-order`) | `-o` BAM, or SAM text for a name ending in `.sam` |
| `revert-original-base-qualities-and-add-mate-cigar` | RevertOriginalBaseQualitiesAndAddMateCigar (a SAM, BAM or CRAM; `--sort-order`, `--restore-original-qualities`, `--max-records-to-examine`) | `-o` BAM (with `.bai` and `.md5` by default), or SAM text for a name ending in `.sam`; nothing for a file in shape |
| `split-sam-by-library` | SplitSamByLibrary (a SAM, BAM or CRAM; a CRAM exits 2 naming wave W15) | `-o DIR`, default the current directory: one BAM, or SAM text, per library |
| `split-sam-by-number-of-reads` | SplitSamByNumberOfReads (a SAM, BAM or CRAM; `--n-files` or `--n-reads`) | `-o DIR`: `shard_0001.bam` and on |
| `gather-bam-files` | GatherBamFiles (BAMs as positionals or `-i`, in the order to gather; a text file of paths, one per line; SAM or CRAM read record by record) | `-o` BAM, no metrics file |
| `calculate-read-group-checksum` | CalculateReadGroupChecksum (a SAM, BAM or CRAM; only its header is read) | the 32 hex digits of the checksum in `-o`, default beside the input |
| `quality-yield-s-n-v-q` | CollectQualityYieldMetricsSNVQ (`qa`, `qc`, `qg` and `qt` on every counted record; `--alternate-quality-attribute`, `--include-bq-histogram`) | one |
| `umi-prevalence` | CollectUmiPrevalenceMetrics (`--barcode-tag`, `--barcode-bq`, `--minimum-mq`, `--minimum-barcode-bq`, `--filter-unpaired-reads`) | one histogram, no metrics row |
| `rrbs` | CollectRrbsMetrics (a coordinate-sorted SAM, BAM or CRAM; `-R` required) | `PREFIX.rrbs_summary_metrics`, `PREFIX.rrbs_detail_metrics` |
| `create-baf-regress-file` | CreateBafRegressMetricsFile (the text bafRegress prints) | `PREFIX.bafregress_metrics` |
| `create-verify-i-d-intensity-contamination-file` | CreateVerifyIDIntensityContaminationMetricsFile (the text VerifyIDIntensity prints) | `PREFIX.verifyidintensity_metrics` |
| `arrays-variant-calling` | CollectArraysVariantCallingMetrics (an arrays VCF with its index; `--dbsnp` required) | `PREFIX.arrays_variant_calling_{detail,summary}_metrics`, `PREFIX.arrays_control_code_summary_metrics` |
| `interval-list-tools` | IntervalListTools (interval lists as positionals or `-i`; `--second-input` for SUBTRACT, SYMDIFF and OVERLAPS) | `-o` interval list, or with `--scatter-count` above 1 the directory `-o` of `temp_NNNN_of_N/scattered.interval_list`; the `--output-value` count on standard output or in `--count-output`; no file without `-o` |
| `scatter-intervals-by-ns` | ScatterIntervalsByNs (a FASTA with its `.fai` and `.dict`; `--max-to-merge`, `--output-type`) | `-o` interval list of `ACGTmer` and `Nmer` runs, `-o -` standard output |
| `non-n-fasta-size` | NonNFastaSize (a FASTA with its `.fai` and `.dict`; `--intervals`) | the count and a line end, in `-o` or on standard output |

A verb runs its tool through the registry, so `bammetrics gc-bias BAM -R ref.fa -o P` and
`bammetrics CollectGcBiasMetrics I=BAM R=ref.fa O=P.gc_bias.detail_metrics
S=P.gc_bias.summary_metrics` write the same metrics (the `StringHeader` command line lists
the options in the order the verb passes them). `crates/cli/tests/grammar.rs` checks that
for every tool.

`oxo-g` follows the naming rule: `CollectOxoGMetrics` without `Collect` and `Metrics` is
`OxoG`, two words in kebab case. It writes one `CpcgMetrics` row per context and library.
`--context-size N` sets the bases on each side of the C (default 1, at most 4 without
`--contexts`), and each `--contexts ACA` adds one context to report:

```bash
bammetrics oxo-g sample.bam -R ref.fa -o sample.oxog_metrics
bammetrics oxo-g sample.bam -R ref.fa --contexts CCG --contexts TCG --db-snp dbsnp.vcf
```

`multiple-metrics --program X` adds to Picard's default five programs, as `PROGRAM=` does;
`--program null` first clears them. The accumulate verbs and `bed-to-interval-list` are
listed here because the registry runs them; they are described under post and prepare,
as is `create-sequence-dictionary`.

`coverage-gaps` and `sex-check` are the two clinical checks of
[docs/06](06_clinical_checks.md), with no Picard counterpart. `bammetrics coverage-gaps
BAM --panel panel.bed --gaps-file gaps.bed -o gaps.txt` writes one
`CoverageGapMetrics` row (`PCT_REPORTABLE`, `GAP_COUNT`, `LARGEST_GAP`,
`PCT_TARGETS_COMPLETE`, a histogram of gap lengths) and the BED of maximal runs of target
bases under `--min-depth` (default 20); its depth is the `CollectHsMetrics` per-base depth,
so `--min-depth 20` gives `PCT_TARGET_BASES_20X`. `bammetrics sex-check BAM --sex male`
writes one `SexCheckMetrics` row: reads per base on X, Y and the autosomes against the
`@SQ` lengths (or the target bases with `--targets`), the ratios, `INFERRED_SEX`
(`FEMALE`, `MALE`, `UNKNOWN`) and `SEX_MATCH` against the declaration.

`bammetrics mt-copy-number BAM -o mt.txt` writes one `MtCopyNumberMetrics` row per
accumulation unit: mean depth on the mitochondrial contig (`chrM`, `MT`, `M` or `chrMT`,
or `--mt-contig`) and on the autosomes from aligned bases at base and mapping quality 20,
and `MT_COPY_NUMBER`, 2 x their ratio ([docs/06](06_clinical_checks.md) section 6). With
`--targets` on an exome the value is a ratio between off-target and on-target depth, not
a copy number. A header without a mitochondrial contig is a usage error naming
`MT_CONTIG` (exit 2; the Picard face exits 1, as for every error there).

`bammetrics adapter-content BAM -o adapters.txt` writes one `AdapterContentMetrics` row
per accumulation unit, read end (`FIRST`, `SECOND`, `UNPAIRED`) and adapter: reads, reads
whose bases in sequencing order hold the adapter's first 12 bases exactly, their mean
start and `PAIR_READ_THROUGH`, the proper pairs with an insert shorter than the read
([docs/06](06_clinical_checks.md) section 7). The `## HISTOGRAM` section is FastQC's
cumulative Adapter Content curve, one column per unit, read end and adapter. The
default adapters are FastQC's six; each `--adapter-sequence NAME=SEQ` replaces them.
Unmapped records count, secondary and supplementary ones do not.

`bammetrics duplicate-by-tile marked.bam -o tiles.txt` writes one `DuplicateByTileMetrics`
row per flowcell lane and tile of the read names, one per lane, one `UNPARSED` row for
names without a tile and one `ALL` row ([docs/06](06_clinical_checks.md) section 8). The
counts are those of `CollectDuplicateMetrics`, from the duplicate flags in the file, so the
`ALL` row equals its totals. `OPTICAL_DUPLICATES` needs the `DT:Z:SQ` tags of
`mark-duplicates --tagging-policy All`; `chain in.bam + mark-duplicates --tagging-policy
All + duplicate-by-tile` marks and counts in one read. `TILE_OUTLIER` is `Y` for a tile
more than `--outlier-sd` (default 3) robust standard deviations, and at least 0.005, above
its lane's median.

`exon-coverage` reports the depth of every exon of an annotation, with no Picard
counterpart ([exon_coverage.md](exon_coverage.md)). `bammetrics exon-coverage BAM
--annotation genes.gtf.gz --gene BRCA1 --gene BRCA2 --canonical --summary-output genes.txt
-o exons.txt` writes one `ExonCoverageMetrics` row per exon and one `GeneCoverageMetrics`
row per gene and per transcript, with the `CollectHsMetrics` per-base depth, so an exon's
`MEAN_DEPTH` and `MINIMUM_DEPTH` are the `mean_coverage` and `min_coverage` of
`PER_TARGET_COVERAGE` for the same interval.

`telomere` measures telomeric read content as TelSeq does, with no Picard counterpart
([telomere.md](telomere.md)). `bammetrics telomere BAM --level READ_GROUP -o tel.txt`
writes one `TelomereMetrics` row per accumulation unit: primary reads that are not QC
failures, unmapped reads included, the reads with at least 7 copies of TTAGGG or CCCTAA
(`-K` sets the 7), the reads of 48 to 52 percent GC and `ESTIMATED_LENGTH_KB`, TelSeq's
telomere length in kb. The histogram counts reads by repeat count.

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

### bam-index-stats

```bash
bammetrics bam-index-stats sample.bam
bammetrics bam-index-stats sample.bam --format tsv       # the lines samtools idxstats prints
bammetrics bam-index-stats sample.cram -R ref.fa -o stats.txt
bammetrics BamIndexStats I=sample.bam
```

Picard's `BamIndexStats` reads the index beside the BAM and prints one line per `@SQ`
line of the header, `NAME length=<tab>LENGTH<tab>Aligned= N<tab>Unaligned= N`, then
`NoCoordinateCount= N`. The counts are the mapped and unmapped record counts `samtools
index` stores in the metadata bin of each reference and the trailing count of records
without a coordinate. Both faces print that text byte for byte; the verb writes it to `-o`.
Picard has no `OUTPUT` option, so `O=` on the Picard face exits 1 with `Unrecognized
option: O`, as in Picard.

- The index is found as htsjdk finds it: `x.bai`, `x.csi`, `x.bam.bai`, `x.bam.csi` beside
  `x.bam`, in that order, then the same names beside the file a symbolic link points to.
  A CSI index holds the same counts. Without one both tools exit 1 with `No index for bam
  file x.bam`.
- A reference the index has no metadata for prints `Aligned= 0` and `Unaligned= 0`. A header
  with more `@SQ` lines than the index has references exits 1 (Picard throws
  `ArrayIndexOutOfBoundsException` after printing the lines it could; bammetrics prints
  nothing and names both counts). An index with more references than the header prints the
  header's references, as Picard does.
- A CRAM has no `.bai`: with `-R` its records are counted (a record without a reference has
  no coordinate, one with a reference is unaligned when its unmapped flag is set, aligned
  otherwise), which gives the text Picard prints for the BAM of the same records. Picard
  3.5.0 opens the CRAM as a BAM and fails. A BAM on standard input (`-`) is counted the same
  way; `--threads` applies to those reads only.
- `--format json`, `jsonl`, `tsv` and `human` print one row per reference (`reference`,
  `length`, `aligned`, `unaligned`) and a last row `*` whose `unaligned` is the
  `NoCoordinateCount`, the row `samtools idxstats` prints; the tsv body is that command's
  output. The json forms carry `--sample` (default the first `@RG SM`).
- `ASSUME_SORTED` is rejected, as Picard rejects it. `VALIDATION_STRINGENCY`, `QUIET` and
  the other plumbing options are accepted and ignored.

### view-sam

```bash
bammetrics view-sam sample.bam | less -S
bammetrics view-sam sample.bam --header-only
bammetrics view-sam sample.cram -R ref.fa --alignment-status Unaligned --records-only
bammetrics view-sam sample.bam --interval-list targets.interval_list -o on_target.sam
bammetrics view-sam reads.sam.gz --validation-stringency LENIENT
bammetrics ViewSam I=sample.bam PF_STATUS=PF
```

Picard's `ViewSam`: the header and every record as SAM text, as htsjdk writes it. Both
faces print it on standard output as the file is read; the verb writes it to `-o`. The
text is byte for byte what Picard prints (`tests/concordance/test_view_sam.py`).

- The header is re-encoded: an `@HD` line first with `VN` in front (`@HD VN:1.6` when the
  input has none), then the `@SQ`, `@RG`, `@PG` and `@CO` lines, each with `SN` and `LN`
  or `ID` first. The version the file states is kept.
- A record's tags come in htsjdk's order (by the second character of the name, then the
  first), a repeated tag once with its last value, every integer as `i`, an `H` value as
  `B:c`, a float as Java prints it (`1.5`, `1.0E-10`).
- `--alignment-status Aligned|Unaligned|All` keeps records by the unmapped flag,
  `--pf-status PF|NonPF|All` by the QC-fail flag. The values are case-sensitive, as in
  Picard.
- `--header-only` prints the header and reads no record. `--records-only` prints no
  header. Both together exit 1 (`Cannot specify both HEADER_ONLY=true and
  RECORDS_ONLY=true.`).
- `--interval-list FILE` keeps the records that overlap an interval of the list. Picard
  has two paths and bammetrics follows the one Picard would take. With a `.bai`, `.csi`
  or `.crai` beside the input, Picard queries the index: an unmapped record placed with
  its mate matches by its start. Without one (SAM text, a pipe), Picard scans: reading
  stops after the last interval and an unmapped record never matches. bammetrics opens
  no index. It reads the whole file and applies the same test.
- Picard has no `OUTPUT` option: `O=` on the Picard face exits 1 (`Unrecognized option:
  O`). `ASSUME_SORTED` is rejected too.
- `--format json`, `jsonl`, `tsv` and `human` give one row per record (`QNAME`, `FLAG`,
  `RNAME`, `POS`, `MAPQ`, `CIGAR`, `RNEXT`, `PNEXT`, `TLEN`, `SEQ`, `QUAL`, `TAGS`) and
  no header line of the SAM. Rows are held in memory: at most 256 MiB of SAM text. The
  default form streams and has no limit.
- A closed standard output (`| head`) ends the run with exit 1 and no message, as
  Picard's `PrintStream.checkError` does.

### sort-sam

```bash
bammetrics sort-sam sample.bam -o sample.qname.bam --sort-order queryname
bammetrics sort-sam sample.qname.bam -o sample.sorted.bam --sort-order coordinate --create-index
bammetrics sort-sam reads.sam -o reads.sorted.sam --sort-order coordinate
bammetrics sort-sam big.bam -o big.sorted.bam -SO coordinate --max-records-in-ram 2000000 --tmp-dir /scratch
bammetrics SortSam I=sample.bam O=sample.sorted.bam SORT_ORDER=coordinate
```

Picard's `SortSam`: every record of a BAM, CRAM or SAM in `queryname`, `coordinate` or
`duplicate` order (the three Picard 3.5.0 takes; `unsorted` and `unknown` exit as Picard
exits). The output is a BAM, or SAM text for a name ending in `.sam`, the records Picard
writes in Picard's order (`tests/concordance/test_sort_sam.py`: decoded records for BAM,
bytes for SAM text).

- The header is the input's with `SO` set and `@HD VN:1.6`. `GO`, `SS`, every other line
  and every field stay, and no `@PG` line is added, as in Picard.
- `queryname` is htsjdk's `SAMRecordQueryNameComparator` (name by bytes, first of pair
  before second, forward before reverse, primary before secondary and supplementary, `HI`),
  `coordinate` its coordinate comparator (unplaced records last), `duplicate` its duplicate
  comparator, which needs the `MC` tag on every paired record with a mapped mate.
  Records with equal keys keep their input order.
- `--max-records-in-ram N` (default 500,000) records are sorted in memory per run; more
  spill to `--tmp-dir` as sorted runs (BGZF level 3, mode 0600, removed at the end, on a
  failure and on a signal) and are merged with ties to the earlier run.
- `--create-index` writes `out.bai` beside a coordinate-sorted BAM and nothing beside SAM
  text, as Picard; `--create-md5-file` writes `out.bam.md5` (or `out.sam.md5`);
  `--compression-level` sets the BGZF level. A CRAM input needs `-R`, and its records get
  the `MD` and `NM` tags htsjdk rebuilds.
- `-o out.cram` exits 2 naming wave W15. An `-o` that is the input (a link or a hard link
  included) exits 2 and replaces nothing.

### add-comments-to-bam

```bash
bammetrics add-comments-to-bam sample.bam -o sample.noted.bam --comment "library prep v2"
bammetrics add-comments-to-bam sample.bam -o sample.noted.bam -C first -C second --create-index
bammetrics AddCommentsToBam I=sample.bam O=sample.noted.bam COMMENT="library prep v2"
```

Picard's `AddCommentsToBam`: the input header with one `@CO` line per `--comment` after the
header's own, then every record of the BAM, byte for byte (`tests/concordance/test_add_comments_to_bam.py`).

- The input must be a BAM, as Picard's block copy needs: a name ending in `.sam` exits 1
  with `SAM files are not supported`, and a CRAM or SAM text under another name exits 1.
- A comment that starts with `@CO` and a tab keeps one prefix; a newline exits 2. Other
  characters go out as htsjdk writes them (`é` as the byte 0xE9).
- No `@PG` line is added. A header without `@HD` gets `@HD VN:1.6`, as htsjdk writes it.
- `--create-index` writes `out.bam.bai` (Picard's block copy names it so) under a header
  that says `SO:coordinate`, `SO:unsorted` or has no `@HD`, and no index under another
  order, as htsjdk's indexer does; `--create-md5-file` writes `out.bam.md5`. An `-o` ending in `.sam`
  exits 2: Picard writes BAM bytes under that name.
- In a chain the step adds the comments to the stream at its position: `bammetrics chain
  in.bam + add-comments-to-bam -C note + write out.bam`.

### replace-sam-header

```bash
bammetrics replace-sam-header sample.bam --header new_header.sam -o sample.reheadered.bam
bammetrics replace-sam-header reads.sam --header other.bam -o reads.reheadered.sam
bammetrics ReplaceSamHeader I=sample.bam HEADER=new_header.sam O=sample.reheadered.bam
```

Picard's `ReplaceSamHeader`: the header of `--header` (SAM text, a BAM or a CRAM, read at
`--validation-stringency`; its records are not read), then every record of the input
(`tests/concordance/test_replace_sam_header.py`).

- An input named `.bam` is a block copy, as in Picard: each record keeps its reference
  index and takes the name the new dictionary gives it. The new `SO` must be `unsorted` or
  the input's. `--create-index` writes `out.bam.bai`.
- Any other input (SAM text, a CRAM with `-R`, a BAM under another name) is read at
  `SILENT`: each record keeps its reference name and takes that name's index in the new
  dictionary. A name the dictionary lacks exits 1 (`Reference index for 'chr1' not found
  in sequence dictionary.`). The two `SO` values must be equal. The output may be SAM text,
  and `--create-index` writes `out.bai`.
- No `@PG` line is added; `SO`, `VN` and every other line are the ones `--header` gives.
- Two cases where Picard writes a file no reader opens exit 1 here: a block-copied record
  whose reference index the new dictionary does not have, and a block copy into a `.sam`
  name (exit 2).
- In a chain the step replaces the header of the stream at its position and every record
  keeps its reference index, as in a block copy.

### add-o-a-tag

```bash
bammetrics add-o-a-tag sample.bam -o sample.oa.bam
bammetrics add-o-a-tag sample.cram -R ref.fa -L targets.interval_list -o sample.oa.sam
bammetrics AddOATag I=sample.bam O=sample.oa.bam INTERVAL_LIST=targets.interval_list
```

Picard's `AddOATag`: each record's alignment appended to its `OA` tag, then every record
(`tests/concordance/test_add_o_a_tag.py`).

- The entry is `RNAME,POS,STRAND,CIGAR,MAPQ,NM;`, for example `chr1,1001,+,50M,60,2;`. An
  unmapped record, placed or not, gets `*,0,+,*,255,;` (`-` with the reverse flag).
- `NM` prints as Java prints the tag's value: an integer in decimal, a float as
  `Float.toString`, a `Z` or `A` value as its text. Without `NM` the field is empty
  (`chr1,1001,+,50M,60,;`).
- The entry goes after the record's own `OA` value. `OA` becomes a `Z` tag at its place in
  htsjdk's tag order. Secondary and supplementary records are tagged as any other.
- `-L` (`INTERVAL_LIST`) takes an interval list with its own `@SQ` lines; contigs match by
  name. Only a record without the unmapped flag that shares a base with an interval gets
  the tag; its span runs from `POS` over the reference bases of its CIGAR, so a CIGAR
  without one never overlaps. A BED exits 1 with `Interval list file must contain
  header.`, as in Picard. A name ending in `.gz`, `.gzip` or `.bfq` is read through gzip,
  as htsjdk reads it. In a chain the list is read before the first read of the source.
- An `NM` or `OA` value that is an array or `H` exits 1: Picard prints it as a Java object
  reference (`[B@1a2b3c`), which changes from run to run.
- The writer checks the header's order, as Picard's presorted writer does. No `@PG` line
  is added. `--create-index` writes `out.bai` for `out.bam` and `out.dat.bai` for
  `out.dat`, as htsjdk names it. CRAM records get `MD` and `NM` as htsjdk decodes them, so
  their entries carry the `NM` htsjdk computes.
- SAM text is written for a name ending in `.sam`, lower case; `OUT.SAM` is a BAM, as in
  Picard. A `.cram` `OUTPUT` is refused before the first read until wave W15: the verb
  exits 2, the Picard face exits 1. Picard writes the CRAM.
- An `@HD SO` value htsjdk does not name fails at `VALIDATION_STRINGENCY=STRICT` (`Found
  non-conforming header SO tag`), as in Picard. At `LENIENT` and `SILENT` it is written
  `SO:unknown`.
- In a chain the step tags the records of the stream at its position:
  `bammetrics chain in.bam + add-o-a-tag -L targets.interval_list + write out.bam`.

### clean-sam

```bash
bammetrics clean-sam sample.bam -o sample.clean.bam
bammetrics clean-sam reads.sam -o reads.clean.sam
bammetrics CleanSam I=sample.bam O=sample.clean.bam CREATE_INDEX=true
```

Picard's `CleanSam`: every record of the input with three fixes
(`tests/concordance/test_clean_sam.py`).

- A mapped record whose alignment ends past the end of its reference is soft-clipped at
  the end of its CIGAR, whatever its strand. Its start stays. The clip point is the read
  length less the overhang, so a deletion or insertion before it moves it, as in Picard:
  `10M5D10M` 5 bases past the end becomes `10M5D6M4S`. A trailing `H` stays after the new
  `S` (`20M5H` becomes `11M9S5H`).
- The `MC` tag of a paired record whose mate is mapped and ends past the end of its
  reference gets the mate's CIGAR clipped the same way.
- An unmapped record gets MAPQ 0.
- `NM` and `MD` are left as they were. No `@PG` line is added. A changed record has its
  tags in htsjdk's order; a record with a new MAPQ only keeps its tag bytes.
- `--validation-stringency STRICT`, the default, reads SAM text at `LENIENT`, as Picard
  does, so a line past the end of its reference is read and named on standard error.
- A read that starts past the end of its reference gets a CIGAR longer than its bases
  (`20M` at 100050 on a 100000-base contig becomes `69S`). Picard writes that record and
  samtools refuses the file; bammetrics exits 1 and writes nothing.
- `--create-index` writes `out.bai`; `-o` ending in `.cram` exits 1 naming wave W15.
- In a chain the step cleans the stream at its position: `bammetrics chain in.bam +
  clean-sam + mark-duplicates + write out.bam`. Its stream is `cleaned`.

### set-nm-md-and-uq-tags

```bash
bammetrics set-nm-md-and-uq-tags sample.bam -R ref.fa -o sample.tagged.bam
bammetrics set-nm-md-and-uq-tags sample.bam -R ref.fa -o sample.tagged.sam --set-only-uq
bammetrics SetNmMdAndUqTags I=sample.bam R=ref.fa O=sample.tagged.bam IS_BISULFITE_SEQUENCE=true
```

Picard's `SetNmMdAndUqTags`: each mapped record gets `MD`, `NM` and `UQ` computed against
`-R` as htsjdk's `SequenceUtil` computes them; every other record is written as read
(`tests/concordance/test_set_nm_md_and_uq_tags.py`). `set-nm-and-uq-tags`
(`SetNmAndUqTags`) is Picard's deprecated name for the same tool.

- The header must say `SO:coordinate`, else exit 1 with `Input must be coordinate-sorted
  for this program to run. Found: queryname`. The reference is read forward, one contig at
  a time; its `.fai` contigs must be the header's.
- A mismatch writes the reference base as the FASTA holds it (`MD:Z:3t4` on a lower-case
  base). Bases compare by IUPAC code, ignoring case: `N` equals `N`, `R` is no `A`.
- `UQ` is the sum of the qualities of the mismatched aligned bases. A record without
  qualities gets no `UQ`. An aligned base past the end of its contig exits 1, as Picard
  does; `MD` and `NM` stop at the contig end.
- `--is-bisulfite-sequence`: a read `T` on a reference `C` (forward strand) or `A` on `G`
  (reverse strand) is no mismatch for `NM` and `UQ`. `MD` is unchanged.
- `--set-only-uq` sets `UQ` alone. A record that gets a tag has its tags in htsjdk's order;
  any other keeps its bytes. No `@PG` line is added.
- In a chain the step reads the chain's `-R`: `bammetrics chain in.bam -R ref.fa +
  sort-sam -SO coordinate + set-nm-md-and-uq-tags + write out.bam`. Its stream is
  `tagged`.

### add-or-replace-read-groups

```bash
bammetrics add-or-replace-read-groups sample.bam -o sample.rg.bam --rglb lib1 --rgpl ILLUMINA --rgpu unit1 --rgsm NA12878
bammetrics add-or-replace-read-groups sample.bam -o sample.rg.bam --rgid 4 --rglb lib1 --rgpl ILLUMINA --rgpu unit1 --rgsm 20 --rgdt 2020-01-02 --sort-order queryname
bammetrics AddOrReplaceReadGroups I=sample.bam O=sample.rg.bam RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=20
```

Picard's `AddOrReplaceReadGroups`: every `@RG` line is replaced by one read group and
every record's `RG` tag is set to its ID
(`tests/concordance/test_add_or_replace_read_groups.py`).

- `RGLB`, `RGPL`, `RGPU` and `RGSM` are required; `RGID` defaults to `1`. The `@RG` line
  holds `ID`, `LB`, `PL`, `SM`, `PU`, then `CN`, `DS`, `DT`, `PI`, `PG`, `PM`, `KS`, `FO`
  when given, in that order.
- Every value must be printable ASCII or a space, one character or more (`^[ -~]+$`);
  each value that is not is named in one usage error, as in Picard.
- `RGDT` takes htsjdk's ISO 8601 forms (`2020-01-02`, `2020-01-02T03:04:05`,
  `2020-01-02T03:04:05.123+05:30`, `2020-13-45` carried into 2021-02-14) and is written in
  UTC: `DT:2020-01-02T00:00:00+0000`. Picard writes the same instant in the JVM's time
  zone. Dates before 15 October 1582 are Julian, as in Java (`1582-10-10` is written
  `1582-10-20T00:00:00+0000`); a year past 100,000,000 is refused.
- `RGID=null` writes `ID:null` and removes every record's `RG` tag, as Picard does.
- `--sort-order` sets `SO`. An order that differs from the input's is sorted as `sort-sam`
  sorts it (`--max-records-in-ram`, `--tmp-dir`); `unsorted` is written as read; `unknown`
  exits 1 on two records or more, where Picard throws a `ClassCastException`. Without a
  sort the records keep their order and a header that names a sort is checked.
- Every record is encoded again as htsjdk encodes a record whose tag changed: tags in
  htsjdk's order and widths. No `@PG` line is added.
- In a chain the step leaves `rg` and streams; with `--sort-order` queryname, coordinate or
  duplicate it is a barrier that sorts. `--sort-order unknown` is a usage error there.

### revert-original-base-qualities-and-add-mate-cigar

```bash
bammetrics revert-original-base-qualities-and-add-mate-cigar sample.bam -o sample.reverted.bam
bammetrics revert-original-base-qualities-and-add-mate-cigar sample.bam -o sample.reverted.bam --max-records-to-examine 0 --sort-order queryname
bammetrics RevertOriginalBaseQualitiesAndAddMateCigar I=sample.bam O=sample.reverted.bam CREATE_INDEX=false
```

Picard's `RevertOriginalBaseQualitiesAndAddMateCigar`
(`tests/concordance/test_revert_original_base_qualities_and_add_mate_cigar.py`).

- The first `--max-records-to-examine` records (10,000) decide whether the file needs
  work: a record with `OQ`, or a first paired record with a mapped mate and no `MC`, starts
  the run; that record with `MC`, or the end of the file, skips it. A skipped file exits 0,
  writes nothing (not even `-o`) and says `Can skip the SAM/BAM/CRAM file` on standard
  error. `0` runs any file with a record.
- Each record gets the soft clip of `clean-sam` where it ends past its reference, and
  with `--restore-original-qualities` (default true) `QUAL` from `OQ`, which is removed.
- The records are sorted by name; each pair gets the mate's reference, position, strand,
  unmapped flag, `MQ`, `MC` and `TLEN` as htsjdk's `SamPairUtil.setMateInfo` sets them; a
  paired supplementary record gets them from its mate. A name with two primary first
  records, or a primary first record without its second, exits 1.
- `--sort-order` (default the input's) is applied after the mate pass: `unsorted` writes
  the name order, `unknown` exits 1 on two records or more.
- `--create-index` and `--create-md5-file` default to true. Every record is encoded again
  with tags in htsjdk's order and widths. No `@PG` line is added.
- In a chain the step is a barrier that leaves `reverted`. It never skips its stream, so
  `--max-records-to-examine` is a usage error there.

### split-sam-by-library

```bash
bammetrics split-sam-by-library sample.bam -o by_library/
bammetrics SplitSamByLibrary I=sample.bam O=by_library CREATE_INDEX=true
bammetrics chain sample.bam + mark-duplicates + split-sam-by-library -o by_library
```

Picard's `SplitSamByLibrary`: one file per library of the header's `@RG` lines
(`tests/concordance/test_split_sam_by_library.py`).

- A library's file is `LB.bam` (`LB.sam` for SAM text input) under the header with that
  library's read groups only. The name is the `LB` value made safe as htsjdk's
  `makeFileNameSafe` makes it: trimmed, then each space and each of
  ``!"#$%&'()*/:;<=>?@[]\^`{|}~`` becomes `_` (`Lib A` gives `Lib_A.bam`, `lib/x` gives
  `lib_x.bam`). A byte above 127 in a BAM header is one character, as htsjdk reads it, so
  `café` in UTF-8 gives `cafÃ©.bam`, Picard's name. SAM text is read as UTF-8, so
  `LB:中文` gives `中文.sam`; the header inside keeps the low byte of each UTF-16 unit, as
  htsjdk writes it.
- Every library gets its file before the first record; a library without records gets a
  header-only file.
- A record without an `RG` tag, with one the header lacks, or whose read group has no `LB`
  goes to `unknown.bam`, made with the first such record, under the header with the read
  groups that have no `LB`.
- A header with no `LB` on any read group exits 2, `No individual libraries are specified
  in the header of INPUT`, as Picard does, and writes nothing.
- Picard writes a CRAM per library for a CRAM input; bammetrics exits 2 naming wave W15.
- Two libraries with one safe name (`a b` and `a_b`), a library `unknown` beside records
  without a library, and names that differ in case only on a file system that ignores case
  make Picard write one file over the other. bammetrics exits 1 and leaves no file.
- Every file is written under a temporary name and renamed into place after the last
  record, with its `.bai` (`--create-index`, a coordinate-sorted BAM) and `.md5`. A run that
  fails leaves none. At most 64 files are open at once; a file closed to make room is
  opened again for its next record. A split writes at most 100,000 files.
- `-o` names the directory (`-o -` exits 2); the verb makes it when missing, the Picard
  face needs it to exist. Without `-o` the files go to the current directory, as Picard's.
- In a chain the step is a split: it writes a BAM per library of the stream at its
  position into its own `-o DIR`. Its files are claimed when it makes them, so no other
  step may write one.

### split-sam-by-number-of-reads

```bash
bammetrics split-sam-by-number-of-reads sample.bam -o shards --n-files 8
bammetrics SplitSamByNumberOfReads I=unmapped.bam O=shards N_READS=48000000 TOTAL_READS=800000000
bammetrics chain sample.bam + split-sam-by-number-of-reads -o shards --n-files 8 --total-reads 12918
```

Picard's `SplitSamByNumberOfReads`: the records in file order over files of about one size
(`tests/concordance/test_split_sam_by_number_of_reads.py`).

- The total is `--total-reads`, else the records of a first read of the input. A pipe
  without `--total-reads` exits 1, as Picard does.
- The number of files is `--n-files`, else `ceil(total / --n-reads)`; each file takes
  `ceil(total / files)` records, then every further record with the name of the last one
  written. A query-grouped input keeps each template in one file; a coordinate-sorted one
  does not, and the run prints Picard's warning. Every record counts, secondary and
  supplementary ones too.
- Files are `OUT_PREFIX_N.EXT`: `N` from 1 with at least four digits, `EXT` the input's
  kind. SAM text input gives SAM text files; a CRAM input gives BAMs named
  `shard_0001.cram`, as Picard writes them. The first file is made before the first record,
  so an empty input gives one header-only file.
- `--n-files` and `--n-reads` together exit 2. One of them must be above 1.
- A `--total-reads` other than the count read changes the sizes and prints Picard's
  warning with the files written.
- Files are renamed into place after the last record and a run that fails leaves none, as
  for `split-sam-by-library`. A run that would make more than 100,000 files exits 1 and
  leaves none; a larger `--n-files` is taken when the records give fewer files.
- In a chain the step needs `--total-reads`: the stream is read once.

### gather-bam-files

```bash
bammetrics gather-bam-files part1.bam part2.bam part3.bam -o sample.bam
bammetrics gather-bam-files -i parts.list -o sample.bam --create-index --create-md5-file
bammetrics GatherBamFiles I=part1.bam I=part2.bam O=sample.bam COMPRESSION_LEVEL=0
```

Picard's `GatherBamFiles`: the records of every input, in the order given, under the
header of the first (`tests/concordance/test_gather_bam_files.py`).

- An input whose name does not end in `.bam`, `.sam` or `.cram` is a list: one path per
  line, each line trimmed, empty lines skipped, a relative path read from the current
  directory. A list may name lists. At most 100,000 inputs and 10,000 lists are read, each
  list at most 64 MiB and all lists of a run together at most 256 MiB.
- When every input name ends in `.bam`, the BGZF blocks are copied as bytes, as Picard
  copies them: the first file whole, each later file from the block of its first record,
  one end-of-file block at the end. Picard deflates that one block per later file again;
  bammetrics deflates it with zlib-rs, so at `--compression-level` 1 to 9 those blocks hold
  other bytes than Picard's and inflate to the same. At level 0, and with one input, the
  file is Picard's byte for byte.
- An input with records whose sequence dictionary (names and lengths) is not the first
  input's exits 1: its records would name other references. Picard writes that file.
- An input whose last bytes are no BGZF block exits 1 with `does not have a valid GZIP
  block at the end of the file`, as Picard does. A file cut at a block boundary is copied.
- With a SAM or CRAM input every record is read (a CRAM with `-R`) and written again into
  a BAM, or SAM text for an `-o` ending in `.sam`; each record takes its reference by name
  in the first header, and a coordinate-sorted header must stay in order, as Picard's
  writer checks.
- `--create-index` writes a `.bai` named as Picard names it: the output name less its
  extension plus `.bai` for a block copy (`out.bai` for `out.bam` and for `out.dat`),
  `out.bai` or `out.dat.bai` otherwise. `--create-md5-file` writes `out.bam.md5`. A block
  copy indexes the records under a first header that says `SO:coordinate`, `SO:unsorted`
  or has no `@HD`, as htsjdk's indexer does, and writes no index under another order. A
  record the indexer refuses (a reference that goes back) exits 1 with no output; Picard
  exits 0 with an empty `.bai`.
- No input, or lists that name none, exits 1. Picard writes a 28-byte file with no BAM
  header.
- An `-o` that is an input or a list exits 2. A block copy into an `-o` ending in `.sam`
  or `.cram` exits 2: Picard writes BAM bytes under that name. A failed run leaves no
  output. Not a chain step: a chain reads one source.

### calculate-read-group-checksum

```bash
bammetrics calculate-read-group-checksum sample.bam
bammetrics calculate-read-group-checksum sample.cram -o sample.rg.md5
bammetrics CalculateReadGroupChecksum I=sample.bam
```

Picard's `CalculateReadGroupChecksum`: the MD5 of the header's `@RG` records as htsjdk
computes it, 32 hex digits with no line end, in `-o` or beside the input as
`sample.bam.read_group_md5` (`tests/concordance/test_calculate_read_group_checksum.py`).

- The records are sorted by `PU`, `LB`, `DT`, `SM`, `CN`, `PL` and `DS`, a missing field
  first, the header order kept for ties. Each record gives its fields in key order, key
  then value, `ID` left out. A header without `@RG` gives `d41d8cd98f00b204e9800998ecf8427e`.
- A BAM's header text counts one character per byte, so `é` stored as UTF-8 is two
  characters; SAM text and a CRAM header are read as UTF-8.
- A `DT` date is written back as htsjdk writes it, in the time zone of the process:
  `2020-01-02` becomes `2020-01-01T17:00:00-0700` in Edmonton. The checksum then depends on
  `TZ`, as Picard's depends on the JVM's zone.
- Every header line is checked at `--validation-stringency`: at `STRICT` (the default) an
  `@RG` without `SM`, two values of one key or a `DT` that is no date exits 1.
- Only the header is read. An `-o` that is the input exits 2. Not a chain step.

### quality-yield-s-n-v-q

```bash
bammetrics quality-yield-s-n-v-q sample.bam -o sample.quality_yield_snvq_metrics
bammetrics quality-yield-s-n-v-q sample.bam --include-bq-histogram --alternate-quality-attribute XQ
bammetrics CollectQualityYieldMetricsSNVQ I=sample.bam O=out.txt INCLUDE_BQ_HISTOGRAM=true
```

Picard's `CollectQualityYieldMetricsSNVQ`: one row of read, base and Q20, Q30 and Q40
totals, and the same counts for the SNV qualities of the Ultima tags `qa`, `qc`, `qg` and
`qt` (`tests/concordance/test_quality_yield_s_n_v_q.py`).

- Each tag holds one Phred+33 character per base, in stored order. At each base the tags
  of the three other bases count, so an `A`, `C`, `G` or `T` gives 3 observations and an
  `N`, an IUPAC code or `=` gives 4. Lower-case bases are read upper case, as htsjdk reads
  them. A quality of 40 counts in Q20, Q30 and Q40.
- Every counted record needs all four tags; a record without one exits 1 naming the record
  and the tag (Picard throws a `NullPointerException`). A secondary or supplementary record
  is skipped before its tags are read unless its `--include-*` flag is set. A tag shorter
  than the read fails only where a base reads past its end; a character outside `!` to `~`
  fails where it is read.
- A tag given twice in a record is read as htsjdk 5.0.0 lists it: the last value counts.
  When the second copy replaces the field htsjdk last appended, the tags after it that sort
  later are lost, so a record `qa qa qc qg qt` has no `qc` and exits 1 in both tools.
  `umi-prevalence` reads `RX` and `BQ` the same way.
- `--alternate-quality-attribute TAG` reads the base qualities from a Phred+33 tag in place
  of `QUAL`; every counted record needs it, with one character per base.
- `--include-bq-histogram` adds one `KEY` histogram holding four series: `BQ_COUNT`,
  `SNVQ_COUNT`, and the mean base and SNV quality at each read index. The quality keys and
  the read indexes share the `KEY` column, as in Picard.
- A chain tap: `bammetrics chain in.bam -o s + quality-yield-s-n-v-q` writes
  `s.quality_yield_snvq_metrics`. Not a `CollectMultipleMetrics` program in Picard 3.5.0.

### umi-prevalence

```bash
bammetrics umi-prevalence sample.bam -o sample.umi_prevalence_metrics
bammetrics umi-prevalence sample.bam --filter-unpaired-reads false --minimum-mq 10
bammetrics CollectUmiPrevalenceMetrics I=sample.bam O=out.txt BARCODE_TAG=RX
```

Picard's `CollectUmiPrevalenceMetrics`: a histogram of duplicate sets by the number of
distinct UMIs in each, `numUmis` against `duplicateSets`, with no metrics row
(`tests/concordance/test_umi_prevalence.py`).

- Records count when mapped, at or above `--minimum-mq` (30), primary, with the
  `--barcode-tag` (RX) tag, past the barcode quality filter, and with
  `--filter-unpaired-reads` (on) paired with a mapped mate.
- The barcode quality filter is Picard's as written: a record with the `--barcode-bq` (BQ)
  tag is kept only when one quality (Phred+33, spaces ignored) is below
  `--minimum-barcode-bq` (30). A record whose barcode qualities are all good is dropped; one
  without the tag is kept (KNOWN_DIFFERENCES.md, "CollectUmiPrevalenceMetrics").
- The input need not be sorted or indexed: the records are sorted in htsjdk's duplicate
  order, `--max-records-in-ram` at a time spilled under `--tmp-dir`, and a set is the run of
  records equal to its first under `duplicateSetCompare` (library, reference, unclipped 5'
  end, orientation, mate). Each read of a pair is its own record. A fragment or a read with
  an unmapped mate joins the pairs on its strand when `--filter-unpaired-reads false` keeps
  it.
- A paired record with a mapped mate needs `MC`, as in Picard; a run that keeps one record
  never compares it and passes without one.
- A sorted run also spills at 256 MiB of keys and UMIs; Picard bounds it by records alone.
  The distinct UMIs of one set are held in memory, as in Picard, so peak memory grows with
  the largest set: 1,000,000 distinct 16-base UMIs at one position took 176 MiB.
- `--progress-step-interval 0` exits 1 once a set is found, as Picard's progress logger
  divides by it. Not a chain step: the tool reads the stream whole before its first set.

### rrbs

```bash
bammetrics rrbs sample.bam -R ref.fa -o out/sample --level READ_GROUP
bammetrics CollectRrbsMetrics I=sample.bam R=ref.fa M=out/sample DO_NOT_CREATE_PLOTS=true
bammetrics chain sample.bam -R ref.fa -o out/sample + rrbs + mark-duplicates + rrbs
```

Picard's `CollectRrbsMetrics`: bisulfite conversion at CpG and non-CpG Cs against the
reference, `PREFIX.rrbs_summary_metrics` with one row per accumulation unit and
`PREFIX.rrbs_detail_metrics` with one row per CpG site and unit
(`tests/concordance/test_rrbs.py`). `-o` (Picard `METRICS_FILE_PREFIX`, `M`) gets a `.`
unless it ends in one. No plot is drawn: `--do-not-create-plots` is accepted.

- Every record not flagged unmapped counts, secondary, supplementary, duplicate and
  QC-failed records included; `--sequence-names` limits the contigs. Reads under
  `--minimum-read-length` (5) and reads with more bisulfite mismatches than
  `Math.round(length x --max-mismatch-rate)` (0.1) are counted and set aside.
- A CpG counts when the read has a C or a T at the C, the G matched, the C at
  `--c-quality-threshold` (20) and the base after it at `--next-base-quality-threshold`
  (10). A reverse-strand read is read reverse complemented. `POSITION` is the 0-based C of
  the forward strand. Detail rows sort by contig name, then position.
- The non-CpG quality test reads the record's qualities in file order at the index within
  the alignment block, as Picard does.
- The reference contig of a record is taken by its index in the FASTA, with no dictionary
  check, as Picard's walker takes it. The `.fai` is required. A record on an earlier contig
  than the one before it exits 1.
- A header that is not `SO:coordinate` exits 1 unless `--assume-sorted`. A record without
  base qualities and a second alignment block exits 1, as Picard does. A chain step: the
  chain's `-R` is its reference.

### create-baf-regress-file, create-verify-i-d-intensity-contamination-file and arrays-variant-calling

```bash
bammetrics create-baf-regress-file bafregress.txt -o chip
bammetrics create-verify-i-d-intensity-contamination-file verifyidintensity.txt -o chip
bammetrics arrays-variant-calling chip.vcf.gz --dbsnp dbsnp.vcf.gz -o chip
bammetrics CollectArraysVariantCallingMetrics I=chip.vcf.gz DBSNP=dbsnp.vcf.gz O=chip
```

The three `picard.arrays` metrics tools of Picard 3.5.0, each writing under the prefix
`-o` (`tests/concordance/test_create_baf_regress_metrics_file.py`,
`test_create_verify_i_d_intensity_contamination_metrics_file.py`,
`test_arrays_variant_calling.py`).

- `CreateBafRegressMetricsFile` reads the text bafRegress prints: the first line must be
  `sample estimate stderr tval pval callrate Nhom` with tabs, and every later line must
  split at spaces and tabs into 7 fields, numbers as Java reads them (`NaN`, `1e-5`).
  `LOG10_PVAL` is `log10(pval)`. The file is `PREFIX.bafregress_metrics`.
- `CreateVerifyIDIntensityContaminationMetricsFile` reads the text VerifyIDIntensity
  prints: the `ID %Mix LLK LLK0` header, a rule of dashes, then one line of four numbers
  per result. The file is `PREFIX.verifyidintensity_metrics`.
- The two text inputs are read whole, at most 64 MiB, 1 MiB per line and 100,000 lines,
  decoded as UTF-8 (a malformed sequence becomes U+FFFD, as in Java). A header alone gives a file
  with no row; any other line Picard refuses exits 1 and writes nothing.
- `CollectArraysVariantCallingMetrics` reads an arrays VCF as GtcToVcf writes it. Picard
  opens it with its index, so a bgzip VCF needs its `.tbi` and a plain one its `.idx`;
  without one the run exits 1. The index is not read. It writes
  `PREFIX.arrays_variant_calling_detail_metrics` (one row per sample, in the order of a
  Java `HashMap` filled in reverse header order),
  `PREFIX.arrays_variant_calling_summary_metrics` (the sum) and
  `PREFIX.arrays_control_code_summary_metrics` (the 23 Illumina controls from their header
  lines). The header must hold `sampleAlias`, `arrayType`, `autocallGender`,
  `autocallVersion`, `autocallDate`, `imagingDate`, `extendedIlluminaManifestVersion`,
  `clusterFile`, `p95Green`, `p95Red`, `scannerName` and the control lines.
- `AUTOCALL_DATE` and `IMAGING_DATE` print in the time zone of the process, as Picard's
  print in the JVM's: set `TZ` to compare runs from two machines. `##imagingDate` is read
  in New York time with the `AM` and `PM` markers GtcToVcf writes; an `en_CA` JVM reads
  `a.m.` and `p.m.` and fails on them (`KNOWN_DIFFERENCES.md`).
- `--call-rate-pf-threshold` (default 0.98, over 0 and at most 1) sets `AUTOCALL_PF`.
  `--num-processors` is taken and changes nothing. `--sequence-dictionary` gives the
  contigs the dbSNP file is checked against. The VCF is read at most 8 MiB per line,
  100,000,000 lines and 50,000 samples. A record counts when it overlaps 1 to the length
  of its `##contig` line. A tab or line end in a header value exits 1 and writes nothing,
  as in Picard.
- None of the three reads a BAM: none is a chain step.

### qc

```bash
bammetrics qc BAM [-o PREFIX] [-R FASTA] [--preset NAME|FILE] [--hs BAITS TARGETS]... [--panel FILE]
              [--ref-flat FILE --strand STRAND [--ribosomal-intervals FILE]] [--sex male|female]
              [--sample NAME] [--fail-on warn|fail] [--resume | --force] [--threads N]
              [--format text|jsonl] [--also FORMAT[,FORMAT]] [-q] [--all-plots] [--date ISO]
              [--report html|pdf|html,pdf] [--pdf-font FONT] [--page-size a4|letter]
```

`qc` runs the collectors the preset's rules need (`generic` without `--preset`; it never
picks a clinical preset by itself), on one read of the BAM, then writes the verdict and
the page. `--hs` or `--panel` adds `CollectHsMetrics` and `--ref-flat` adds
`CollectRnaSeqMetrics` whether or not a rule reads them. `--sex` is the declared sex for
`CollectSexCheckMetrics` (`EXPECTED_SEX`), which only the `clinical-panel` preset runs
(`generic` and `wgs` have no `SexCheckMetrics` rule, so `--sex` under them exits 2); the
panel targets (`--panel`, else the first `--hs` pair's targets) reach it and
`CollectCoverageGaps`. A rule whose collector cannot run
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
| SexCheckMetrics | CollectSexCheckMetrics | (`--sex`, the panel targets when given) | `PREFIX.sex_check_metrics` |
| CoverageGapMetrics | CollectCoverageGaps | the panel targets | `PREFIX.coverage_gaps`, `PREFIX.coverage_gaps.bed` (manifest kind `gaps`, not a verdict input) |

Each metrics file is byte for byte what `bammetrics <Collector> INPUT=BAM OUTPUT=FILE ...`
writes for the argv the manifest records, apart from the `Started on` line. The header
command line follows that argv: `qc` passes `REFERENCE_SEQUENCE` only to the collectors that
read `-R` (alignment summary, hs, gc bias, wgs), so the insert size and duplicate headers
lack the `REFERENCE_SEQUENCE=` that `bammetrics CollectInsertSizeMetrics ... R=ref.fa`
prints; the data lines are the same. A CRAM input needs `-R`, and then every collector gets
`REFERENCE_SEQUENCE`, since each one decodes the CRAM with it. The sample name (`--sample`,
else the first `@RG SM`, else the BAM name) goes only into the verdict, the page and the
manifest. No `-o` means `./<BAM name without .bam or .cram>`, and the summary says so; a directory `-o` gets the BAM name as
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

`--report` picks the sample report: `html` (the default, `PREFIX.report.html`), `pdf`
(`PREFIX.report.pdf` alone) or `html,pdf`. The PDF needs a binary built with the cargo
feature `pdf`; without it `--report pdf` exits 2 before anything is read and names the build
command. The manifest lists the PDF under the step `report` with kind `pdf`, and a later run
without it removes it. `--pdf-font` and `--page-size` set the PDF and need `--report pdf`
([docs/pdf.md](pdf.md)).

`--also` writes further formats beside the Picard files, one or a list (`--also
json,markdown`). `json`, `jsonl`, `tsv` and `csv` add `FILE.json` (`.jsonl`, `.tsv`, `.csv`)
beside every metrics file. `markdown` adds `PREFIX.report.md` (the verdict, then every
metrics table, at most 65,536 characters), `junit` adds `PREFIX.junit.xml` (a testcase per
rule) and `multiqc` adds `PREFIX.verdict_mqc.json` (MultiQC custom content). The manifest
lists them under the step `also`. A later run without `--also` removes the ones the earlier
manifest lists ([docs/formats.md](formats.md)). `--also` does not change the exit code: a
BAM with no read exits 0 with it as without it, and a metrics file with headers only gets a
sibling with no row.

`qc` exits 0 on any verdict unless `--fail-on` is given; then 10 when the verdict reaches it
and 11 when the verdict is incomplete, as `check` does. A failed collector exits 1 after the
other files, the verdict and the page are written.

### chain

```bash
bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]
bammetrics chain --plan FILE [OPTIONS] [BAM]

bammetrics chain in.bam -R ref.fa -o out/NA12878 \
    + mark-duplicates \
    + alignment-summary + insert-size + quality-yield + duplicate-metrics \
    + write --index
bammetrics chain in.bam -o out/NA12878 --sort-memory 4G \
    + sort-sam --sort-order queryname + mark-duplicates \
    + sort-sam --sort-order coordinate + write out/NA12878.bam --index
samtools view -b in.bam chr20 | bammetrics chain - + mark-duplicates + write - | samtools flagstat -
bammetrics chain --plan plan.toml -o out/NA12878 --dry-run
```

Several steps on one BAM in one process. Records move from step to step in memory: no
shell pipe and no intermediate file. [docs/chaining.md](chaining.md) is the reference; this
is the short form.

- The chain's options and the source (a BAM, a CRAM with `-R`, `-` for standard input)
  come before the first `+`. Every `+` starts a step.
- A step is a tool verb, its Picard class, or `write`. A transform (`mark-duplicates`,
  `sort-sam`, the six record-editing verbs of wave W09: `add-comments-to-bam`,
  `replace-sam-header`, `add-o-a-tag`, `clean-sam`, `set-nm-md-and-uq-tags`,
  `set-nm-and-uq-tags`, and the two of wave W10: `add-or-replace-read-groups`,
  `revert-original-base-qualities-and-add-mate-cigar`) changes the records later steps
  see; any number may follow one another. A tap (the 20 one-pass tool verbs) reads the
  records of its position and changes nothing. A split (`split-sam-by-library`,
  `split-sam-by-number-of-reads`) writes them over BAMs in a directory. `write [PATH]` writes those records as a
  BAM, or as SAM text for a name ending in `.sam`; `-` is standard output, and `write -
  --output-format sam` prints SAM text there. A `.cram` or `.fastq` name is refused with
  the wave of docs/07 that adds it.
- Each step takes its verb's own options in every spelling (`--minimum-pct 0.1`,
  `--MINIMUM_PCT 0.1`, `MINIMUM_PCT=0.1`), with did-you-mean per step. `bammetrics chain
  BAM + VERB --help` lists them.
- Files go under `-o PREFIX` with the names `qc` writes (Output naming, above). A step's
  own `-o` names its file. `--format picard,json` and `--also csv,markdown` add further
  formats beside every metrics file.
- The chain's options (`-o`, `--threads`, `--sample`, `--also`, `--plan`, `--dry-run`,
  `-q`, `-v`) given after the first `+` are a usage error that says so: `step 1
  (insert-size): --also belongs to the chain: give it before the first +`. An empty `-o`
  and an empty file name of a step are usage errors too.
- `--plan FILE` reads the source, the options and the steps from a TOML file. An option
  on the command line replaces the file's. An unknown key is an error with the line and
  the closest name. Paths in the file are read relative to the working directory. The
  plan file is an input: no step may write it. docs/chaining.md section 3 lists the part
  of TOML the reader takes.
- `--dry-run` prints the plan: the segments, one line per barrier (`sort-sam` holds its
  records; `mark-duplicates` replays the source, a copy of the pipe, or a spool of its
  input behind another transform), the bytes each spills, the files of each step, the
  estimated temporary and output disk and the peak memory. `bammetrics plan` with a
  chain's arguments prints the same.
- `--sort-memory BYTES` (default 2G) bounds what a `sort-sam` step holds before it
  spills. A chain whose estimate exceeds the free space of a file system it writes to
  exits 1 before the first read; `--force-space` runs it. `--validation-stringency`
  sets how a SAM text source is checked.
- Standard error gets each transform's record counts, each failed step and one summary
  line. `-q` leaves only the errors.

Exit 0 when every step finished. Exit 2 for a usage error in any step: it is found before
the first read, names the step (`step 2 (insert-size): unknown option --minmum-pct`), and
nothing is written. Exit 1 when a step failed: a tap fails alone and every other step
writes its files; a failing transform or `write` stops its segment, and the summary says
which (`failed alone`, `failed and stopped its segment`, `stopped`, `not run`).

A chain without a transform is the pass of `qc`: the same metrics files, byte for byte
apart from the `Started on` line. It writes no verdict, page or manifest.

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
  `start`, `end`, `strand`, `name`), with no `sample`. `--threads` is checked (at most
  256) and ignored: no BAM record is read. `--sample` is accepted and ignored: an interval
  list has no sample.

`hs-metrics`, `collect` and `qc` read a BED directly, so the verb is for Picard tools and
other programs that take only interval lists.

In `qc`, `--panel BED` runs the same conversion with the BAM as the dictionary and keeps
`PREFIX.panel.interval_list`, the bytes `bammetrics bed-to-interval-list BED --from-bam
BAM` writes. A panel whose first non-blank line starts with `@` is an interval list
already, and prepare writes nothing for it: the lists `IntervalListTools`,
`ScatterIntervalsByNs` and `VcfToIntervalList` write go to `qc` and `hs-metrics --panel` as
they are. A panel named `.gz`, `.gzip` or `.bfq` is inflated first, as htsjdk inflates it,
then told by its first line in the same way; a gzip BED goes to `BedToIntervalList`, which
inflates it as Picard does. A VCF (`##fileformat=VCF` on its first line, gzip or not)
exits 2 (`a VCF is no panel`): `VcfToIntervalList` converts it. Gzip data under another
name (`.bgz`) exits 2, as does a panel named `.interval_list` whose first line is no `@`
header line: the hs reader takes either name at its word. A `--panel` BED that is itself named
`PREFIX.panel.interval_list` exits 2 (`is an input; writing it would destroy it (from
--panel)`), as does any planned file that is an input. The hs collector reads the `--panel`
file itself, so its file, header line and `BAIT_SET` are the same with or without the
prepare step. The prepare step is a step of the manifest (stage `prepare`, output kind
`intervals`, the BED under `files`), so `--resume` keeps the interval list, a lost one is
written again, and a later run without the BED removes it as a stale file. It is not
counted as a collector; when it fails (a contig the BAM header lacks), the summary lists it
as an `ERROR` line and `qc` exits 1 (`exit 1: 1 collector(s) and 1 prepare step(s) failed`
when hs fails on the same BED).

### interval-list-to-bed

```bash
bammetrics interval-list-to-bed panel.interval_list -o panel.bed
bammetrics interval-list-to-bed panel.interval_list -o - --score 0 --sort false
bammetrics IntervalListToBed I=panel.interval_list O=panel.bed.gz
```

Picard's `IntervalListToBed`: one BED line per interval of the list, byte for byte the
lines Picard writes. Each line holds the contig, the start less 1, the end, the name (`.`
stays `.`), `--score` (`SCORE=`, default 500, any integer) and the strand. `--sort`
(default true) orders the lines by contig, start, end, strand and name; `--sort false`
keeps the file order.

- The list is read as Picard's feature reader reads it: a record on a contig the header
  lacks is skipped, and a file with no `@` line has no contig, so its BED is empty. A
  record Picard rejects (not five fields, a start below 1, an end below start - 1, a strand
  other than `+` or `-`, an end past a contig length other than 0) exits 1.
- `-o FILE` is required; `-o -` writes standard output. A name ending in `.gz` is
  compressed (BGZF, which inflates to Picard's text). `-o` naming the input exits 2, and a
  failed run leaves no file.
- A list of more than 500,000 intervals goes through Picard's sort spill, which keeps the
  low byte of each name character: a name past U+00FF prints as Picard prints it. A name
  with a character whose low byte is 0 exits 1.
- An empty first line, or an empty line right after the `@` lines, exits 1, as in Picard.
  A list of more than 20,000,000 intervals exits 1.
- `--threads` is checked (at most 256) and ignored. `--format` takes `picard` only.

`bed-to-interval-list --dict LIST` converts the BED back. The intervals come back in
sorted order under `@HD VN:1.6 SO:coordinate` when each has a length above 0 and a name
without a space or a double quote (`tests/concordance/test_interval_list_to_bed.py`
checks both tools). `qc --panel` reads the BED.

### vcf-to-interval-list

```bash
bammetrics vcf-to-interval-list calls.vcf.gz -o calls.interval_list
bammetrics vcf-to-interval-list calls.vcf -o - --include-filtered --variant-id-method USE_FIRST
bammetrics VcfToIntervalList I=calls.vcf O=calls.interval_list
```

Picard's `VcfToIntervalList`, byte for byte: `@HD VN:1.6`, one `@SQ` line per `##contig`
line (`LN:0` without a length, `AS` from `assembly`), then one interval per record from
`POS` to `END` (else the end of the reference allele) on `+`, named by its `ID` or
`interval-N` for the Nth record without one. A record that overlaps or abuts the interval
before it in file order joins it: with `--variant-id-method CONCAT_ALL` (the default) the
interval spans both and joins the distinct names with `|`; with `USE_FIRST` it keeps the
first start and name. An unsorted VCF merges only neighbours, as Picard does.

- `--include-filtered` (`INCLUDE_FILTERED=`, `IF=`) keeps records whose `FILTER` is neither
  `.` nor `PASS`.
- A record on a contig without a `##contig` line, or past its length, is written as it is.
  A VCF without `##contig` lines exits 1 (Picard fails with a `NullPointerException`).
- Each record is checked as htsjdk's codec checks its site: a missing column, a `POS`,
  `QUAL` or `END` that is not a number, an empty `ID` or `INFO`, a bad `FILTER` name, an
  allele that is not one, or a repeated allele exits 1 with htsjdk's message. A failed run
  leaves no file.
- The header is read as htsjdk's codec reads it: `VCFv4.0` to `VCFv4.3`, the last
  `##fileformat` line winning; a repeated `ID`, a `##INFO` line without `Description` and
  a `##` line without `=` are read. `VCFv3.3`, `VCFv4.4`, `VCFv5.0`, a `##INFO` line with
  a bad `Type` or `Number` and an unclosed quote exit 1. A header over 64 MiB exits 1.
- The input is plain text or gzip (BGZF or plain), told by its first bytes. A `.bcf` input
  exits 2: BCF is not read.
- Each interval is written as it is made, so memory does not grow with the VCF.
- `-o FILE` is required; `-o -` writes standard output; a name ending in `.gz` is compressed.
  `--threads` is checked (at most 256) and ignored. `--format` takes `picard` only.

`qc --panel` takes the interval list as it is.

### create-sequence-dictionary

```bash
bammetrics create-sequence-dictionary ref.fa -o ref.dict
bammetrics create-sequence-dictionary ref.fa.gz --genome-assembly GRCh38 --species "Homo sapiens"
bammetrics CreateSequenceDictionary R=ref.fa        # writes ref.dict, as Picard does
```

The verb writes the bytes `CreateSequenceDictionary R=FASTA O=FILE` writes: `@HD VN:1.6`,
then one `@SQ` line per FASTA sequence with `SN`, `LN`, `M5` (the MD5 of the bases in
upper case), `AS` when given, `UR` (default `file://` and the FASTA's absolute path,
encoded as Java encodes it) and `SP` when given, then `AN` for a sequence `--alt-names`
names. Without `-o`, or with `-o -`, the same bytes go to standard output. The Picard face
without `O=` writes the FASTA name with its extension replaced by `.dict`, beside it.

- The FASTA name ends in `.fasta`, `.fas`, `.fa`, `.fna` or `.txt`, each also with `.gz`;
  any other name exits 1 (`File is not a supported reference file type`). A `.gz` file is
  read as plain gzip or bgzip, member after member. The `.fai` is never read.
- `O=FILE`, and `-o FILE` in picard format, refuse a file that exists (exit 1, `already
  exists.  Delete this file and try again`). With `--format json`, `jsonl`, `tsv` or
  `human`, `-o FILE` replaces an existing file, as in every verb. The verb exits 2 when
  `-o` or `FILE.md5` names the FASTA or the `--alt-names` file.
- A line feed in `--genome-assembly`, `--uri` or `--species` exits 2 before any file is
  written; Picard fails there too (`Supplied String contains illegal character`).
- `--truncate-names-at-whitespace` (default true) names a sequence by the first word of
  its `>` line. `--num-sequences N` stops after N sequences, at least one.
  `--create-md5-file` also writes `FILE.md5` and needs `-o`.
- `--format json`, `jsonl`, `tsv` and `human` print one row per sequence (`name`,
  `length`, `md5`, `assembly`, `uri`, `species`, `alternate_names`), with no `sample`.
  `--threads` and `--sample` are accepted and ignored.

The tools that look for a `.dict` beside a FASTA (`bed-to-interval-list --dict ref.fa`,
`SD=ref.fa`) read the file this verb writes, and Picard's CollectWgsMetrics takes it as
the dictionary of `R=ref.fa`.

### scatter-intervals-by-ns and non-n-fasta-size

```bash
bammetrics scatter-intervals-by-ns ref.fa -o ref.acgt.interval_list --output-type ACGT --max-to-merge 50
bammetrics ScatterIntervalsByNs R=ref.fa O=ref.interval_list
bammetrics non-n-fasta-size ref.fa
bammetrics NonNFastaSize I=ref.fa INTERVALS=panel.interval_list O=panel.non_n
```

Both read the FASTA as htsjdk's indexed reader does: the name ends in a FASTA extension,
the `.fai` beside it gives each contig's offset and line geometry, and the `.dict` (the
FASTA name with its extension replaced by `.dict`, else the name plus `.dict`) must list
the `.fai`'s contigs in its order with its lengths. A bgzip FASTA named `.gz` is read
through its `.gzi`. The bases are read 1 MiB at a time, whatever the contig length; the
`.fai` and `.dict` are read whole, at most 64 MiB each
(`tests/concordance/test_scatter_intervals_by_ns.py`, `test_non_n_fasta_size.py`,
fixtures from `tests/concordance/n_run_refs.py`).

- `scatter-intervals-by-ns` writes `@HD VN:1.6 SO:coordinate`, the `.dict`'s `@SQ` lines,
  then one interval per run: `Nmer` for a run of `N`, `n` or `.`, `ACGTmer` for a run of any
  other byte (IUPAC codes, `-` and `*` included). An `Nmer` of at most `--max-to-merge`
  bases (default 1) between two `ACGTmer`s joins them, from the start of the list, so a
  joined interval may join the next. `--output-type` (`ACGT`, `N` or `BOTH`, the default;
  upper case only) picks the intervals written. Each interval is written as soon as it is
  decided: at most three are held. `-o` is required; `-o -` writes standard output, and a
  name ending in `.gz` is written as BGZF. `--format` takes `picard` only. Python returns
  one row per contig: its length and the `ACGTmer` and `Nmer` intervals and bases.
- `non-n-fasta-size` counts every base but `N` and `n` (IUPAC codes and `.` count) and
  writes the number and a line end. `--intervals FILE` (an interval list) counts only the
  bases inside it; its `@SQ` line `k` stands for the FASTA's contig `k`, whatever the name,
  as Picard's mask matches them. `--format json` and the other forms print one row,
  `non_n_bases`. `--threads` is checked (at most 256) and ignored; `--sample` is accepted
  and ignored.
- Both exit 1 without a `.fai`, without a `.dict` or with one of no `@SQ` line, with a
  `.dict` whose contigs or lengths differ from the `.fai`, and on a `.fai` line with 0
  bases per line, as Picard fails. `scatter-intervals-by-ns` also exits 1 on a contig of
  no base, which Picard cannot index. A FASTA that ends before the bases its `.fai` gives,
  and a gzip FASTA that is not bgzip, exit 1 where Picard reads bytes that are no bases
  (`KNOWN_DIFFERENCES.md`). A `.fai` line is read as htsjdk's `Scanner` reads it (the
  pattern `(.+)\t+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)` found in the line, a blank line before a
  later entry and a last entry with no line end refused). A FASTA or `.dict` that is no
  regular file exits 1 before it is opened. No output may be the FASTA, its `.fai`, `.gzi`
  or `.dict`, or the `--intervals` file; a failed run writes no file, and with `-o -` prints
  nothing: the list is copied to standard output once it is whole.
- Neither reads a BAM, so neither is a step of `chain` (exit 2).

### mark-duplicates

```bash
bammetrics mark-duplicates in.bam -o marked.bam --metrics marked.markdup_metrics
bammetrics mark-duplicates a.bam b.bam -o - --metrics ab.markdup_metrics --tagging-policy All | samtools flagstat -
bammetrics MarkDuplicates I=in.bam O=marked.bam M=marked.markdup_metrics TMP_DIR=/scratch
```

The verb runs `MarkDuplicates`, Picard 3.5.0's duplicate marker (0.18.0). `-o BAM` is
`OUTPUT`, the BAM it writes; `-o -` writes it to standard output. `--metrics-file`
(`--metrics`, `-M`, `METRICS_FILE=`) is the `DuplicationMetrics` file. Both are required.
`--format` takes `picard` first; `--format picard,json` adds `METRICS_FILE.json` beside the
metrics file.

- `OUTPUT` holds every input record in input order. The duplicate flag (0x400) is set on
  each duplicate and cleared on every other record, so a flag from an earlier run does
  not survive. `--remove-duplicates` leaves duplicates out, and
  `--remove-sequencing-duplicates` leaves out optical duplicates.
- The header gains one `@PG` line per distinct `PG` value of the input records, chained
  with `PP` to that value, and each record's `PG` tag names its new line. The ID is
  `MarkDuplicates`, then `MarkDuplicates.1` and so on when the ID is taken, as Picard
  numbers them. `VN` is `bammetrics` and the version, `CL` the bammetrics command line;
  `--program-group-version` and `--program-group-command-line` replace them.
  `--program-record-id null` adds no line and changes no `PG` tag;
  `--no-add-pg-tag-to-reads` keeps the lines and leaves the tags.
- The metrics file has one row per library and the `set_size` histograms `all_sets`,
  `optical_sets` and `non_optical_sets`; with one library and a size estimate the
  `CoverageMult` column comes first and the key reads `BIN`.
- The input may be coordinate sorted or query-name sorted (by the header or by
  `--assume-sort-order`). In query-name order the secondary, supplementary and
  unmapped records of a duplicate template are flagged too. Any other order exits 1
  after the first read. Several `-i` inputs are merged as Picard merges them; their
  dictionaries must be equal.
- The input is read twice. Standard input and pipes are copied to `--tmp-dir` first.
  Read ends and index lists are sorted in runs of `--max-records-in-ram` entries (1 or
  more, default 4,000,000) spilled to `--tmp-dir` (default the system temporary
  directory). Mates waiting for a mate on a later reference are spilled there too, one
  file per reference; memory holds those waiting on the reference being read. The
  directory `bammetrics-mark-duplicates.<pid>.<16 hex digits>` is made at the first
  spill with mode 0700. Its last part is random, and a name already there is an error,
  never reused. It is removed at the end, on an error, and on SIGINT, SIGTERM or SIGHUP;
  a SIGKILL leaves it. `sort-sam` names its directory `bammetrics-sort-sam.`,
  `check-duplicate-marking` `bammetrics-check-duplicate-marking.` and `chain`
  `bammetrics-chain.`.
- `--tmp-dir` needs about 55 bytes per input record at its peak, during the first read:
  6.5 GiB for the 150.5M-record HG002 exome, 42.1 GiB for a 768.6M-record genome (about
  the size of the output BAM). Peak memory was 0.8 GB on the exome and 2.75 GB on the
  genome at the default `--max-records-in-ram`.
- `-o` and `--metrics-file` are written under temporary names in their directories and
  renamed at the end, so a failed or interrupted run leaves an earlier file whole. An
  output that is an input (by path, symbolic link or hard link), or two outputs that
  name one file, exit 2 before anything is read.
- A CRAM input needs `-R`. Its records come out with the `MD` and `NM` tags htsjdk
  rebuilds while decoding. The output is always BAM: an `-o` ending in `.sam` or `.cram`
  exits 2 (the Picard face exits 1); `bammetrics chain IN + mark-duplicates + write
  out.sam` writes SAM text. SAM text input is read, plain or in gzip members, at
  `--validation-stringency`.
- `--create-index` writes `marked.bai` beside a coordinate-sorted BAM, the index Picard
  `BuildBamIndex` writes for that BAM, byte for byte;
  `--create-md5-file` writes `marked.bam.md5`. `--compression-level` sets the BGZF level
  (default 5). `--threads` sets the BGZF decompression threads, or the CRAM container
  decoders, and the BGZF compression workers of the output.

KNOWN_DIFFERENCES.md lists what differs from Picard: the `@PG` `VN` and `CL`, the flow
options it refuses, the `READ_NAME_REGEX` syntax and the `MAX_RECORDS_IN_RAM` meaning.

### build-bam-index

```bash
bammetrics build-bam-index sample.bam                         # sample.bai beside the BAM
bammetrics build-bam-index sample.bam -o index/sample.bam.bai
bammetrics build-bam-index sample.bam -o - > sample.bai
samtools sort in.bam | tee sorted.bam | bammetrics build-bam-index - -o sorted.bai
bammetrics BuildBamIndex I=sample.bam O=sample.bai
```

The verb runs `BuildBamIndex`. It writes the `.bai` Picard 3.5.0 writes for the same BAM,
byte for byte: the bin htsjdk computes from each record's start and end, one chunk per
run of records in one BGZF block, the linear index, the metadata bin 37450 of each
reference and the count of records without a coordinate (`crates/core/src/bai.rs`).
samtools and htsjdk read the index as they read their own. `samtools index` writes
another file for the same BAM: 3,694,104 bytes against 7,186,968 on the HG002 exome.

- The header must say `SO:coordinate`. `queryname`, `unsorted`, `unknown` and a header
  without the tag exit 1 with Picard's `Input bam file must be sorted by coordinate`.
- Without `-o` the verb writes beside the BAM: `x.bai` for `x.bam`, the name plus `.bai`
  for any other name. `-o -` prints the index. Standard input (`-`) needs `-o`.
- The Picard face without `O=` writes that file name into the current directory
  (`BuildBamIndex I=dir/x.bam` writes `./x.bai`), as Picard 3.5.0 does. Picard's help says
  beside the BAM.
- The index is written under a temporary name and renamed. A refused or failed run
  leaves no file and keeps an earlier index. `-o` may not name the input.
- `--format` takes `picard` only: an index has one form. `build_bam_index()` in Python
  returns the rows `bam-index-stats` prints for the index written.
- A CRAM is refused: its index is a `.crai`. In a chain the index comes from `write
  --index`; `build-bam-index` as a step exits 2 and says so.
- `--threads` sets the BGZF decompression threads. The index is the same at any count.
  The HG002 exome BAM (10.3 GB, 150.5 M records) takes 11.8 s and 9 MiB at 8 threads;
  Picard takes 85.6 s and 487 MiB (docs/benchmarks.md).

KNOWN_DIFFERENCES.md lists what differs from Picard: records are not validated, an empty
BGZF block in the middle of a file does not end the read, and the verb's default
directory.

### fifo-buffer

```bash
bwa mem -t 8 ref.fa r1.fq r2.fq | bammetrics fifo-buffer | samtools sort -o sorted.bam -
producer | bammetrics fifo-buffer --buffer-size 1073741824 --debug-frequency 10 --name bwa | consumer
producer | bammetrics FifoBuffer BUFFER_SIZE=1073741824 DEBUG_FREQUENCY=10 NAME=bwa | consumer
```

The verb runs `FifoBuffer`, Picard 3.5.0's memory buffer between two commands of a pipe. It
copies standard input to standard output. Any bytes pass unchanged: text, BAM, an empty
input. The verb takes no input argument, no `-o`, no `--format` and no `--sample`; each
exits 2 with the usage line `bammetrics fifo-buffer [OPTIONS] < IN > OUT`.

- `--buffer-size` (`BUFFER_SIZE`) is the most bytes held, 1 to 2147483645, default
  536870912 (512 MiB). The producer waits only when that many bytes are held. The bytes
  are reserved, and written only as the buffer fills: the buffer's resident memory is at
  most twice the most bytes ever held and at most the buffer size. A 4 GiB copy to a
  reader that keeps up peaked at 7 to 11 MiB, a full default buffer at 515 MiB
  (docs/benchmarks.md).
- `--io-size` (`IO_SIZE`) is the most bytes one read or one write moves, default 65536.
  The reader and the writer each hold one chunk of that size, cut to the buffer size:
  128 KiB at the defaults. On an empty input `--io-size 1073741824 --buffer-size 10`
  peaked at 2.8 MiB and `--io-size 1073741824 --buffer-size 1073741824` at 1.76 GiB.
- `--debug-frequency N` (`DEBUG_FREQUENCY`) logs the bytes held on standard error: one line
  before the first read, then one every N seconds. 0 or less logs none. `--name`
  (`NAME`) puts a name in the line. `--verbosity WARNING` or `ERROR` drops the lines.
- The Picard face prints the line as Picard does, without the timestamp:
  `INFO<TAB>FifoBuffer<TAB>Fifo buffer bwa used 1,024 / 1,073,741,824 (0%).` The verb prints
  `bammetrics fifo-buffer: Fifo buffer bwa used 1,024 / 1,073,741,824 (0%).`
- The verb prints one more line on standard error when it starts:
  `bammetrics fifo-buffer: bammetrics chain runs its steps in one process and removes the
  pipe this buffer serves (-q drops this line)`. `-q` or `QUIET=true` drops it. The Picard
  face never prints it. `-v` adds the Picard command line and `copied N of M bytes read`.
- When the reader of standard output closes it, the copy stops and exits 0 with no
  message, and the writer of standard input sees its pipe close. Any other failed write
  exits 1 (`standard output: ...`). A failed read exits 1 after the bytes read before it
  are written.
- Exit codes: 0 after the copy. A bad value or an option Picard refuses (`INPUT`,
  `OUTPUT`, `ASSUME_SORTED`, `GA4GH_CLIENT_SECRETS`) exits 2 on the verb and 1 on the
  Picard face, before one byte is copied.

`fifo-buffer` is no step of `bammetrics chain`. A chain hands records from step to step
in memory and has no pipe to buffer, so `chain BAM + fifo-buffer` exits 2 and says so
([docs/chaining.md](chaining.md) section 1). The verb still serves a pipe into or out of
a chain: `bammetrics chain in.bam + mark-duplicates + write - | bammetrics fifo-buffer |
samtools flagstat -`.

`fifo_buffer()` in Python raises `UsageError`: a Python call has no pipe to buffer
([docs/python.md](python.md)). KNOWN_DIFFERENCES.md lists what differs from Picard: the
first log line, a size of 0, a closed or failed standard output.

### doctor

```bash
bammetrics doctor FILE [-R FASTA] [--records N] [--check-md5] [--expect-build BUILD]
    [--expect-sex male|female] [--sample-sheet FILE --sample NAME]
    [--format text|json|jsonl] [--fail-on warn|fail] [--threads N]
```

`doctor` checks one BAM, CRAM or SAM file before a long run: the kind and version, the
end-of-file marker, the header, `@SQ` count and length, the genome build (GRCh37, hg19,
GRCh38, hg38, T2T-CHM13 or unknown), the `@RG` lines, the sort order against the records,
then over the first 10,000 records the read groups, paired, mapped and duplicate fractions,
read length, quality range (Phred+64), MD and NM; with `-R` the reference's names, lengths
and MD5s against `@SQ`; an index beside the file; the sample sheet row of `--sample`. It
reads no more than that: 0.35 s on a 16 GB CRAM. Each check prints one line with `PASS`,
`WARN`, `FAIL` or `SKIP`, the value and the reason. It writes no file. Exit 10 when a check
fails (or warns, with `--fail-on warn`), 1 when the file, FASTA or sheet cannot be read.
Every check, threshold and form is in [docs/doctor.md](doctor.md).

```bash
bammetrics doctor sample.cram -R ref.fa --expect-build hg38
bammetrics doctor sample.bam --sample-sheet samples.csv --sample NA12878 --format json
```

### interval-list-tools

`IntervalListTools` of Picard 3.5.0 (wave W12, `crates/core/src/interval_list_tools.rs`).
Each input list is padded by `--padding`, the lists are combined by `--action`, then the
result is sorted, inverted (`--invert`), merged (`--unique`, `--dont-merge-abutting`) and
broken at `--break-bands-at-multiples-of`, in Picard's order. `-o` gets the list, or with
`--scatter-count` (or `--scatter-content`) above 1 a directory of
`temp_NNNN_of_N/scattered.interval_list` files, renamed into place as one set; a run that
fails leaves none. `--output-value BASES` or `INTERVALS` prints the count of what was
written to standard output, or to `--count-output`; without `-o` no list is written. The
`@PG` line holds Picard's command line for the same options, so the files equal
Picard's byte for byte. A file named exactly `.vcf`, `.vcf.gz` or `.vcf.bgz` is read as a
VCF (`--include-filtered` keeps its filtered records); every other name is an interval
list, as in Picard. A chain refuses the verb as a step: it reads no BAM.

Stated limits: a list or VCF holds at most 20,000,000 intervals, 1,000,000 header lines
and 8,000,000 header fields; an intersection (`INTERSECT`, `SUBTRACT`, `SYMDIFF`) makes
at most 10,000,000 overlapping pairs and 256 MiB of names; a scatter makes at most
100,000 lists and 32 GiB of text. A value with a line feed exits 2 before any file, as
Picard fails before it reads. A `--count-output` inside a list directory of the scatter
is refused.

```bash
bammetrics interval-list-tools panel_a.interval_list -i panel_b.interval_list --action UNION -o union.interval_list
bammetrics IntervalListTools I=targets.interval_list SI=blacklist.interval_list ACTION=SUBTRACT PADDING=50 O=kept.interval_list
bammetrics interval-list-tools genome.interval_list --scatter-count 24 -M BALANCING_WITHOUT_INTERVAL_SUBDIVISION -o shards
bammetrics interval-list-tools targets.interval_list --unique --output-value BASES
```

`qc --panel` takes the list it writes.

## Post verbs

### accumulate-quality-yield and accumulate-variant-calling

```bash
bammetrics accumulate-quality-yield out/a.qy out/b.qy -o out/sample.qy
bammetrics AccumulateQualityYieldMetrics I=out/a.qy I=out/b.qy O=out/sample.qy
```

Files in, file out; no BAM is read. The two forms write the same metrics.
`accumulate-variant-calling` takes the prefixes `variant-calling` wrote and writes the two
files under `-o PREFIX`, inverting the printed ratios into counts as Picard does.

### convert-sequencing-artifact-to-oxo-g

```bash
bammetrics sequencing-artifact sample.bam -R ref.fa -o out/sample
bammetrics convert-sequencing-artifact-to-oxo-g out/sample -o out/sample
bammetrics ConvertSequencingArtifactToOxoG I=out/sample
```

Reads `PREFIX.pre_adapter_detail_metrics` and `PREFIX.bait_bias_detail_metrics` and writes
`PREFIX.oxog_metrics`, the `CpcgMetrics` rows `CollectOxoGMetrics` writes: one per library
and C-centred context. The positional argument is `INPUT_BASE`; `--pre-adapter-in` and
`--bait-bias-in` name the two files when `FILE_EXTENSION` renamed them, `--oxog-out` the
output. With `--pre-adapter-in`, `--bait-bias-in` and `--oxog-out` the verb needs neither a
positional prefix nor `-o`. The Picard face defaults `OUTPUT_BASE` to `INPUT_BASE`; the verb,
like every prefix verb, defaults `-o` to `./` plus the prefix's file name. `--oxog-out` may
not name either detail file: the run fails and leaves both in place. A detail file cut by
`CONTEXTS_TO_PRINT` that lacks a context the conversion needs fails the run, as in Picard.
The name is the naming rule applied to the class (`oxo-g` from `OxoG`).

### compare-metrics

```bash
bammetrics compare-metrics picard.insert_size_metrics ours.insert_size_metrics
bammetrics compare-metrics a.hs_metrics b.hs_metrics -o diff.txt --marc MEAN_TARGET_COVERAGE:0.01
bammetrics CompareMetrics I=a.txt I=b.txt O=diff.txt OUTPUT_TABLE=diff.table KEY=CATEGORY
```

Two metrics files in, the comparison out; Picard 3.5.0's `CompareMetrics`, exit code
included: 0 equal, 1 NOT equal or an error. The verb writes the comparison (`Comparison of
<class> metrics between files A and B`, `Metrics are equal` or `NOT equal`, then one line
per difference) to `-o`, or to standard output. The Picard face writes it to `O=` only and
prints Picard's log lines on standard error without their timestamp: `INFO` for the
summary, `ERROR` for each difference. `OUTPUT_TABLE` (`--output-table`) gets the value
differences as a `CompareMetrics$MetricComparisonDifferences` metrics file.

- Values are compared typed, by the Java type of the class's field, from a table read out
  of `picard.jar` (`crates/core/picard_fields.tsv`, `make picard-fields`). `0.5` and `0.50`
  in a double column are one value; `5.0` in a long column is an error, as in Picard.
- Rows equal when every field prints the same 6-decimal text, so a change below the sixth
  decimal does not count unless another field of the files differs.
- Rows pair by position (`Key: 0`), or by the `--key` columns' values. A key value that is
  empty is an error in Picard and here.
- `--marc NAME:FRACTION` (`METRIC_ALLOWABLE_RELATIVE_CHANGE`) allows a relative change
  below FRACTION. A change from 0 is a relative change of `1.7976931348623157E308`; a NaN
  relative change is within any allowance.
- `--metrics-to-ignore` leaves a column out but both files must have it;
  `--metrics-not-required` lets either file lack it. `--ignore-histogram-differences`
  stops a histogram difference from failing the files.
- `--format json`, `jsonl`, `tsv` and `human` print one row per difference line: `KEY`,
  `METRIC`, `VALUE1`, `VALUE2` for a value, `DIFFERENCE` for every line. `--threads` and
  `--sample` are accepted and ignored.
- The text on standard output ends with a newline; the `-o` file is Picard's `OUTPUT`
  byte for byte, which has none.
- An `-o` or `--output-table` that names an input, the two naming one file, and `- -`
  (standard input read twice) are usage errors (exit 2). The Picard face writes over an
  input as Picard does.
- `-h` and `--help` print this help and exit 0; Picard prints its usage and exits 1.

### check-duplicate-marking

```bash
bammetrics check-duplicate-marking marked.qname.bam
bammetrics check-duplicate-marking marked.bam -o bad_querynames.txt --mode PRIMARY_MAPPED_ONLY
bammetrics CheckDuplicateMarking I=marked.qname.bam O=bad_querynames.txt MODE=PRIMARY_ONLY
```

A marked BAM or CRAM in, the querynames whose records disagree on the duplicate flag
(0x400) out; Picard 3.5.0's `CheckDuplicateMarking`, exit code included: 0 when every
queryname agrees, 1 when a record differs from the first record of its queryname or on an
error. Whether the marking is right is not checked. The verb writes one line per record
that disagrees (a name repeats when two of its records do) to `-o`, or to standard output,
and the summary (`All records' duplicate markings agree.` or `Found N records that do not
agree on their duplicate flag.`) on standard error. The Picard face writes `O=` only and
prints Picard's log lines on standard error without their timestamp: one `WARNING` per bad
record for the first 100 (`Reads with queryname NAME have different duplicate flags (at
CONTIG:START)`, `null` for an unmapped record), then `Further warnings will be
suppressed.`, then the summary at `ERROR` or `INFO`.

- `SO:queryname` is read in file order, and only adjacent records are one queryname: a
  name that comes back later starts a new group, as Picard's single current name does.
  Any other header (coordinate, unsorted, unknown, no `@HD`) logs `Input file isn't
  queryname sorted. Sorting into temp space.` and is sorted first, in Picard's
  `SAMRecordQueryNameComparator` order (name by bytes, first of pair before second,
  paired before unpaired, forward before reverse, primary before supplementary before
  secondary, `HI` ascending after records without one, file order between ties).
  bammetrics keeps the sort keys in memory and writes no temporary file.
- `--mode` (`MODE`) is `ALL` (default), `PRIMARY_ONLY` (secondary and supplementary records
  dropped), `PRIMARY_MAPPED_ONLY` (unmapped records dropped too) or
  `PRIMARY_PROPER_PAIR_ONLY` (proper pairs only; an unpaired record is an error, as
  Picard's `getProperPairFlag` throws). A dropped record does not change the current name.
- `--verbosity ERROR` drops the INFO and WARNING log lines, `WARNING` the INFO lines. The
  verb's summary on standard error follows the same rule: it is an INFO line when every
  name agrees and an ERROR line otherwise, so `--verbosity ERROR` prints it only on a bad
  marking.
- `--tmp-dir` and `--max-records-in-ram` are accepted and ignored: a file without
  `SO:queryname` is sorted in memory, so memory grows with the file.
- Read names print byte by byte, as htsjdk decodes them: a name holding UTF-8 bytes
  above 127 prints each byte as one character, as Picard prints it.
- An `-o` (or `O=`) that is the input file is refused before the read.
- `--format json`, `jsonl`, `tsv` and `human` print one row per bad record: `read_name`,
  `contig` (empty for an unmapped record) and `start` (0 without a position).
- `ASSUME_SORTED` and `STOP_AFTER` are not options of Picard's tool and exit 1 on the
  Picard face, 2 on the verb.
- Python `check_duplicate_marking(bam, mode=..., output=...)` returns the rows; an empty
  list means every queryname agrees. It raises nothing for a bad marking, only for an
  error (a bad `mode`, an unpaired record under `PRIMARY_PROPER_PAIR_ONLY`).

### report

```bash
bammetrics report DIR... [--preset NAME|FILE] [-o OUT] [--all-plots] [--date ISO]
                  [--pdf | --pdf-only] [--pdf-font FONT] [--page-size a4|letter]
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

`--pdf` adds `OUT/SAMPLE.report.pdf`, the page's content as a tagged PDF/UA-1 file with an
appendix of every metrics table. `--pdf-only` writes `OUT/SAMPLE.verdict.json` and
`OUT/SAMPLE.report.pdf` per sample and no page or index. `--page-size` is `a4` (default) or
`letter`; `--pdf-font FONT` sets the text from one TTF or OTF file for sample names outside
Latin. Both need a binary built with the cargo feature `pdf`; a binary without it exits 2
and names `cargo build --release -p bammetrics-cli --features pdf`
([docs/pdf.md](pdf.md)).

The page and the index are self-contained HTML with no script and no network request: a
`<meta>` CSP of `default-src 'none'; style-src 'unsafe-inline'`, the CSS inlined, plots as
inline SVG with a text alternative. They follow the light or dark system theme, print on a
light palette, and pass axe-core 4.10.2 with no violation at 320, 360 and 1280 px in light,
dark and print (`make a11y`, recorded in [docs/benchmarks.md](benchmarks.md)). Every text color pair meets WCAG AA (4.5:1) in light, dark and print, checked
by a unit test over the tokens of `report.css`. Size and render time are budgets of docs/05
section 4.8, measured by `make budgets` ([docs/benchmarks.md](benchmarks.md)).

### check

```bash
bammetrics check DIR|VERDICT.json... [--preset NAME|FILE] [--fail-on warn|fail]
                 [--format FORMAT[,FORMAT]] [-o FILE]
```

`check` is the gate. A directory is discovered and evaluated as `report` does; a
`verdict.json` file is read back as is. `--format text` prints a band line per sample and
every check under it; `json` prints a JSON array of verdicts, each as `verdict.json`;
`jsonl` prints one object per line, a `check` event per check and a `verdict` event per
sample. `junit` prints JUnit XML, a testsuite per sample and a testcase per rule: a fail is
a `<failure>`, a rule not measured a `<skipped>`, a warn a `<failure>` only under `--fail-on
warn`. `markdown` (`md`) prints the checks as GFM tables, at most 65,536 characters.
`multiqc` prints MultiQC custom content, a table with a row per sample; save it under a name
that ends in `_mqc.json`.

Without `-o` the output goes to standard output and no file is written. `-o FILE` writes
the first format to `FILE`. `--format` takes a list (`text,junit`): each further format is
written beside `-o` as `FILE.txt`, `FILE.json`, `FILE.jsonl`, `FILE.junit.xml`, `FILE.md`
or `FILE.verdict_mqc.json`, and a list needs `-o`. `check` refuses an `-o` or a sibling that
is a file it reads (exit 2): a metrics file under an input directory, or a `verdict.json`
given or found there. A new file in an input directory is written. Exit codes: the table
above, with `--fail-on fail` by default, whatever the format.

### convert

```bash
bammetrics convert FILE [--format FORMAT[,FORMAT]] [-o FILE] [--table NAME] [--sample NAME]
```

`convert` reads one metrics file and writes it in another format. The input is Picard text
or the JSON of `--format json` (schema `bammetrics.metrics/1`), told apart by its first
character (`{` or `[` is JSON; a byte order mark before it is dropped); `-` reads standard
input. `--format` takes `picard`, `json`, `jsonl`, `tsv`,
`csv`, `markdown` or `human`, one or a list; without it the output is `json` for Picard
text and `picard` for JSON. Picard text to JSON to Picard text is byte for byte the first
file, for a file Picard or bammetrics wrote.

```bash
bammetrics convert S1.hs_metrics --format csv > S1.hs_metrics.csv
bammetrics convert S1.insert_size_metrics --format csv --table histogram
bammetrics convert S1.insert_size_metrics | bammetrics convert - | cmp - S1.insert_size_metrics
```

- `-o FILE` writes the first format there; each further format of a list goes beside it.
  `convert` refuses to write its input (exit 2).
- `--table metrics|histogram|CLASS` picks the table `tsv`, `csv` and `human` print.
- `--sample NAME` sets the sample of the `json` and `jsonl` forms.
- A file that is neither form exits 1 with the reason. [docs/formats.md](formats.md) lists
  what a round trip does not keep.
- `--format picard` exits 1 and writes nothing for a JSON cell with a tab or a line end, a
  row that would be a blank line and a first cell that starts with `## `: Picard text has no
  escape for them. The error names the row and the column.

## Info verbs

- `bammetrics tools`: every verb and the Picard class it runs, then the other verbs.
- `bammetrics presets`: the built-in verdict presets with version, the first 6 characters
  of the content hash, rule count and note (`clinical-panel   v3  6f0d59  18 rules  reviewed
  defaults; gap and sex rules are starting values`). Every rule is listed in
  [docs/verdict.md](verdict.md).
- `bammetrics help VERB|PicardClass`: the verb's options in every spelling; `help qc` (or
  `plan`, `check`, `report`, `convert`, `doctor`) prints clap's help. `help chain` lists the chain's
  options, the step kinds and the options of `write`; `bammetrics chain BAM + VERB --help`
  lists the options of one step.
- `bammetrics completions bash|zsh|fish|elvish|powershell`: a completion script on stdout
  (`bammetrics completions zsh > ~/.zfunc/_bammetrics`). `chain` completes its own
  options; the words after a `+` are not completed.

### plan

```bash
bammetrics plan BAM [the qc options] [--format text|jsonl]
bammetrics plan BAM [the chain options] + STEP [STEP OPTIONS] [+ STEP ...]
bammetrics plan --plan FILE [the chain options]
```

With a chain's arguments (a `+` step, `--plan FILE`, or the word `chain` first) `plan`
prints the plan `chain --dry-run` prints: the segments, replay or spool, the files of each
step and the estimated temporary disk ([docs/chaining.md](chaining.md) section 4). It reads
the header of a file source and no record. With `qc` arguments it is what follows.

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
bammetrics CollectJumpingLibraryMetrics I=lane1.bam I=lane2.bam O=out.txt MQ=20 T=100
```

The Picard face takes `KEY=VALUE`, `--KEY VALUE`, `--KEY=VALUE` and `-K VALUE` with long
or short names, writes `O=` (standard output without it) and exits 0 or 1. It is
byte-identical to 0.12, apart from `--help` alone, which now prints the verb help, and an
option given twice. An option Picard takes once, given twice (`SP=a SP=b`, `O=a O=b`), exits
1 with Picard's `Option 'SPECIES' cannot be specified more than once.`; 0.12 kept the last
value. A list option (`METRIC_ACCUMULATION_LEVEL`, `BAIT_INTERVALS`, `PROGRAM`) adds each
value. The verbs keep the last value of an option given twice. `ASSUME_SORTED=false` exits
1 on a BAM whose header is not coordinate sorted, as Picard does.
`METRIC_ACCUMULATION_LEVEL` (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`; `null`
clears the default) works where Picard has it. Where the output differs from Picard on
purpose is listed in KNOWN_DIFFERENCES.md.

## Chaining examples

The worked examples of docs/05 section 2.5 that work today, written against the committed
fixtures so each runs from the repository root as shown (`tests/data/insert_size.bam` has
no `@RG SM`, so the sample is the BAM name). Output is what 0.30.0 prints; long lines are
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
insert_size  FAIL  preset clinical-panel v3 (6f0d59)  6 collectors, 1 BAM read (10 records)
FAIL          hs.mean-target-coverage   Mean target coverage is 0.377778, below the fail limit of 100
FAIL          hs.target-20x             Target bases at 20x is 0, below the fail limit of 0.95
FAIL          hs.target-100x            Target bases at 100x is 0, below the fail limit of 0.6
FAIL          hs.zero-coverage-targets  Targets with zero coverage is 0.2, above the fail limit of 0.03
FAIL          hs.gc-dropout             GC dropout is 20.588235, above the fail limit of 10
FAIL          gaps.reportable-range     Target bases at MIN_DEPTH is 0, below the fail limit of 0.95
FAIL          gaps.targets-complete     Targets without a gap is 0, below the fail limit of 0.9
+2 more, see out/NA12878.report.html
PASS          9 checks
wrote 11 files to out/NA12878.*, report out/NA12878.report.html
use --fail-on fail or bammetrics check to gate
```

The two lines past the 12-line budget are `WARN sex.inferred` and `NOT MEASURED
hs.evenness`. The 11 files add `NA12878.coverage_gaps`, `NA12878.coverage_gaps.bed` and
`NA12878.sex_check_metrics` to the run of 0.17.0.

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
2 samples, preset clinical-panel v3 (6f0d59), bammetrics 0.30.0
FAIL        NA1  7 fail, 1 warn, 9 pass, 1 not measured
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

**9. MarkDuplicates to CollectDuplicateMetrics.**

```
$ bammetrics chain tests/data/markdup.bam -o out/s + mark-duplicates + duplicate-metrics + write
bammetrics chain: step 1 (mark-duplicates): 36 records in, 36 out, 16 changed
bammetrics chain: 3 steps finished, 3 files written
$ ls out
s.duplicate_metrics
s.marked.bam
s.marked.duplication_metrics
```

One process, no pipe and no second read of the marked BAM: `duplicate-metrics` gets the
marked records in memory. The same with two processes and a pipe, as docs/05 first wrote
it:

```bash
bammetrics mark-duplicates tests/data/markdup.bam --metrics out/s.markdup_metrics -o - \
  | tee out/s.marked.bam | bammetrics duplicate-metrics - -o out/s.dup
```

The marker reads its input twice, so the input must be a file (standard input is copied
to `--tmp-dir` first). `duplicate-metrics` reads the flags it set: its pair and duplicate
counts equal the marker's, and its `READ_PAIR_OPTICAL_DUPLICATES` is 0, as in Picard.

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
`tools` and `presets` run the command's code (0.15.0). `chain(source, steps, ...)` runs
the `chain` verb, with `step(verb, **options)` for a step that has options. Every function
and type is in [docs/python.md](python.md).

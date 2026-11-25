# Exon coverage (CollectExonCoverageMetrics)

The depth of every exon of a gene annotation, wave W07 of
[docs/07](07_full_parity_program.md). Picard has no such tool. The depth is the one
`CollectHsMetrics` counts, so an exon's mean and minimum equal the `mean_coverage` and
`min_coverage` that `CollectHsMetrics PER_TARGET_COVERAGE` prints for the same interval.

| | |
|---|---|
| Picard class name | `CollectExonCoverageMetrics` |
| Verb | `exon-coverage` |
| Python | `collect_exon_coverage_metrics(bam, annotation=..., ...)` |
| Metrics classes | `bammetrics.qc.ExonCoverageMetrics` (`OUTPUT`), `bammetrics.qc.GeneCoverageMetrics` (`SUMMARY_OUTPUT`) |
| `CollectMultipleMetrics`-style names | `PREFIX.exon_coverage_metrics`, `PREFIX.gene_coverage_metrics` |
| Source | `crates/core/src/exon_coverage.rs`, `crates/core/src/gene_model.rs` |
| Oracle tier | A: `CollectHsMetrics PER_TARGET_COVERAGE` of Picard 3.5.0 |

```bash
bammetrics exon-coverage sample.bam --annotation gencode.v46.gtf.gz --gene-list panel.txt \
    --canonical --summary-output sample.gene_coverage_metrics -o sample.exon_coverage_metrics
bammetrics CollectExonCoverageMetrics I=sample.bam ANNOTATION=genes.bed GENE=BRCA1 MIN_DEPTH=20 MIN_DEPTH=50 O=out.txt
```

## 1. Annotation

`ANNOTATION` is read line by line, plain or gzip (one or more members, BGZF too), told by
its first two bytes. The format is told by the first line that is not blank, a `#` comment,
or a `track` or `browser` line. A 9-field line is GTF or GFF3 when its third field is text,
so a BED9 line with an integer name and score is BED. A byte order mark before the first
line is dropped.

| Format | Told by | Exons | Names |
|---|---|---|---|
| GTF 2.2 | 9 fields, `key "value";` attributes | `exon` lines (`CDS` with `FEATURE=CDS`) | `gene_id` and `transcript_id` required; `gene_name`, `exon_number` and each `tag` read |
| GFF3 | 9 fields, `key=value` attributes, or a `##gff-version 3` line | `exon` or `CDS` lines; each `Parent` is a transcript | the transcript and gene from the `Parent` chain, in any line order; `Name` of the gene as the gene name; `rank` as the exon number; `%XX` escapes decoded |
| BED12 | 12 fields | each block; with `FEATURE=CDS` the part of each block inside `thickStart` to `thickEnd` | the name is the gene and the transcript |
| BED | 3 to 11 fields | each line, a region | the name (or `chrom:start-end` without one) is the gene and the transcript; lines of one name are numbered in file order |

Coordinates are 1-based and closed in every row: a BED start is read as `start + 1`. A GTF
or GFF3 position is 1 to 2,147,483,647, as htsjdk keeps it. A file stops at a `##FASTA`
line. Every malformed line is an error that names the file and the line: a field count, an
integer, a strand, a quote, a block outside its line, a feature without a `gene_id` or
`transcript_id`, a GFF3 `Parent` no line declares, text that is not UTF-8.

Size caps: a line holds at most 1 MiB, an annotation at most 5,000,000 features of the kind
read (GENCODE 46 holds 1.65 million exon lines), a BED12 line at most 100,000 blocks, a
`GENE_LIST` file at most 16 MiB. The 5,000,000 cap also bounds the genes, the transcripts
and the GFF3 parent features: an `exon` line of a `FEATURE=CDS` run adds a gene and a
transcript and no feature. A `GENE_LIST` byte order mark is dropped. A BED12 block count and its size and offset lists must
agree before any block is placed.

### 1.1 Selection

| Option | Keeps |
|---|---|
| `GENE` (repeated) | the genes of that `gene_name` or `gene_id`; an Ensembl ID matches with or without its `.version` |
| `GENE_LIST` | the same, one name per line (the first word; `#` lines skipped), added to `GENE` |
| `TRANSCRIPT` (repeated) | the transcripts of that ID, with or without the `.version` |
| `CANONICAL` | one transcript per gene: tagged `Ensembl_canonical`, else `MANE_Select` (or `MANE Select`), else the most CDS bases, else the most exon bases, else the lowest ID |
| `FEATURE` | `exon` (default) or `CDS` |

The filters combine: a transcript is kept when every given filter keeps it. A requested gene
or transcript the annotation lacks is a usage error that names each one (the first with the
closest name). A selection that keeps no feature is a usage error. `CANONICAL` on a BED file
is a usage error; a BED file has no transcripts.

### 1.2 Merging and numbering

Overlapping features of one transcript merge into one exon, which keeps the lowest exon
number written. Abutting features stay apart. Features of different transcripts never merge:
each transcript reports its own exons, so a base two transcripts share is in a row of each.
An exon without a written number is numbered along its transcript from the 5' end: on the
minus strand the exon with the highest start is 1.

### 1.3 Contigs

Each annotation contig is looked up in the BAM header by its own name, then with `chr` added
or removed, then for a mitochondrial name the other three of `chrM`, `MT`, `M` and `chrMT`
(`genome_build::aliases`, the table `doctor`, `sex-check` and `mt-copy-number` read). A contig read under another name is
listed on standard error (`1 annotation contig read under the BAM header's name: 1 as
chr1`). An exon on a contig the header lacks under every name gets a row with empty depth
cells and is not counted; standard error names the contigs. The verb's and the chain's `-q`
drop both notes, and the Python function prints neither. When no selected exon is on the
header the run fails, as it does for an exon that ends past its contig's `@SQ` length: the
annotation is for another build.

## 2. Depth

The tool runs the `ALL_READS` unit of `hs_metrics::CaptureMetrics` with the placed exons as
both baits and targets and reads its depth arrays (`CaptureMetrics::into_target_depths`), as
`CollectCoverageGaps` does ([docs/06](06_clinical_checks.md) section 1.1). The filters are
the Hs ones: secondary, unmapped, QC-failed and duplicate records add nothing, an adapter read
at MAPQ 0 adds nothing, a read under `MINIMUM_MAPPING_QUALITY` (default 20) adds nothing, a
base under `MINIMUM_BASE_QUALITY` (default 20) or at quality 2 or below adds nothing, and with
`CLIP_OVERLAPPING_READS` (default true) the overlap of two mates counts once. Supplementary
records count. With `INCLUDE_INDELS` deleted bases count as covered.

Overlapping exons of different transcripts are one merged Hs target; each exon reads its own
bases from that target's array, so its depth does not depend on the other transcripts.

`COVERAGE_CAP` caps the Hs histograms and the base quality counts. It never caps the per-base
depth that `PER_TARGET_COVERAGE` and this tool read, so it is not an option.

## 3. Rows

`OUTPUT` holds one `ExonCoverageMetrics` row per selected exon, in the file order of each
exon's first feature: the gene and transcript names, `EXON_NUMBER`, `CHROM` (the BAM
header's name), `START`, `END`, `STRAND`, `LENGTH`, `MEAN_DEPTH` (the depth sum over
`LENGTH`), `MEDIAN_DEPTH` (the mean of the two middle depths for an even length),
`MINIMUM_DEPTH`, `MAXIMUM_DEPTH`, then `BASES_<N>X` for each `MIN_DEPTH` and `PCT_<N>X` for
each. `MIN_DEPTH` is repeated (at most 16 values, each at least 1); given values replace the
default 10, 20 and 30, sorted and without repeats.

`SUMMARY_OUTPUT` holds `GeneCoverageMetrics` rows: one `GENE` row per gene and contig, then a
`TRANSCRIPT` row per transcript of it. A row's bases are the union of its exons, so a base two
transcripts share counts once in the `GENE` row; introns are not counted. `EXON_COUNT` counts
distinct exon intervals.

`SUMMARY_OUTPUT` is Picard text in every face: `--format` applies to `OUTPUT` alone, as it
does to the further files of every verb. `convert` and Python `read_table` read the summary
file. The verb refuses a `SUMMARY_OUTPUT` that names `-o` or a file `--format` writes beside
it, and makes a missing directory above it before the first read, as it does for `-o`.

The glossary rows `BASES_<N>X` and `PCT_<N>X` stand for every `MIN_DEPTH`
(`glossary::column_matches`), so `--format json` and Python `read_table` type a `PCT_45X`
column as a float even when every cell prints `0`.

## 4. Correctness

No Picard tool computes per-exon depth from an annotation. Three checks hold the numbers:

- Planted truth. `crates/core/tests/exon_coverage.rs` writes a SAM file whose depth at every
  base is known (a duplicate, a MAPQ 10 read, a quality 2 read, a secondary and a QC-failed
  record that add nothing; a 10-base deletion; a pair whose 20-base overlap counts once) and
  a GTF with two genes, two overlapping transcripts, a one-base exon on the deletion, a
  single-exon gene on the minus strand under the alias `2`, and a gene on `chrUn`. Every cell
  of every exon, gene and transcript row is computed in the test from the planted intervals.
  The same exons as BED12 lines with unordered blocks give the same rows.
  `crates/core/src/gene_model/tests.rs` plants the reader cases (GTF, GFF3 with parents after
  children, BED12 offsets and thick parts, BED, gzip, canonical choice, every error).
- Tier A. `tests/concordance/test_exon_coverage.py` writes a GTF from each input's targets
  (genes of three exons, every other gene on the minus strand, a shifted second transcript, a
  one-base exon, a gene off the BAM), splits the exons into layers where no two touch (Picard
  merges overlapping and abutting targets), and runs `CollectHsMetrics PER_TARGET_COVERAGE`
  of Picard 3.5.0 and of BAMmetrics once per layer. Each exon's `MEAN_DEPTH` and
  `MINIMUM_DEPTH` must print as `mean_coverage` and `min_coverage`. Ten cases: the foundation
  corpus, the exome slice at the defaults, `MQ=0 Q=0`, `CLIP_OVERLAPPING_READS=false` and
  `INCLUDE_INDELS=true Q=30`, the `hs_bams.py` BAMs `deep`, `first_empty` and `unknown`, and
  corpus tiers t1 and t2 (`MQ=1`): 2,674 cells, all equal.
- Cross-check. `samtools bedcov -Q MAPQ -j` sums the depth of every base of the reads its
  flags keep (unmapped, secondary, QC-failed and duplicate dropped). Its mean is never below
  `MEAN_DEPTH`. With `MQ=0 Q=0 CLIP_OVERLAPPING_READS=false` the two differ only by the bases
  at quality 2 or below and the MAPQ 0 adapter reads: at most 0.55% per exon on the exome
  slice, 0.48% on the foundation corpus, 0 on `deep` (the test allows 1%). At the defaults the
  gap is the base quality filter and the mate overlap: up to 15.2% on the exome slice and
  8.2% on the foundation corpus. mosdepth, the oracle docs/07 names, is not installed
  (docs/07 section 6, Q9).

## 5. Scale

Measured on 2026-11-17 on the Mac Studio (M1 Ultra), `--threads 4`, the HG002 exome BAM
(10 GB, hs37d5, contigs `1` to `22`, `X`, `Y`, `MT`) against GENCODE v19 (chr names, so all
25 contigs are read under the alias):

| Selection | Rows | Wall time | Peak RSS |
|---|---:|---:|---:|
| every exon of every transcript | 1,196,293 | 93.2 s | 3.36 GiB |
| `--canonical --feature CDS` | 197,363 | 71.9 s | 1.07 GiB |

The depth arrays hold 8 bytes per merged exon base, and every row is held until the end of
the pass.

## 6. Verdict rule and plot

No shipped preset has an exon rule, so `qc` never runs this tool and its files are
unchanged. `qc` has no annotation input, and per-exon thresholds are set per assay. A lab
writes the rule docs/07 names, `exon.min-pct20`, in its own preset; `check` and `report`
read the two files like any metrics file. The verdict takes the worst row the scope
selects, so this rule fails on the least covered gene. The values are starting values, not
validated limits:

<!-- exon preset start -->
```toml
[preset]
name = "exon-panel"
version = 1
note = "starting values"

[[rule]]
id = "exon.min-pct20"
label = "Least covered gene at 20x"
class = "GeneCoverageMetrics"
metric = "PCT_20X"
select = "LEVEL=GENE"
warn = "< 0.98"
fail = "< 0.95"
```
<!-- exon preset end -->

`crates/cli/tests/exon_coverage.rs` runs this preset through `check`. No rule names a plot,
so the report draws no per-exon plot.

## 7. Departures from docs/07

- No shipped `exon.min-pct20` rule (section 3.4): `qc` output stays byte for byte, and the
  rule lives in a lab's preset (section 6 above).
- No per-exon plot family (section 3.4): no shipped rule names one.
- `COVERAGE_CAP` is not an option: it does not change the per-target depth.
- The cross-check is `samtools bedcov`; mosdepth (W07 row) is not installed.
- `doctor` and `MtCopyNumberMetrics` of the W07 row, and the D8 and D13 changes that go
  with `doctor`, are other parts of the wave.

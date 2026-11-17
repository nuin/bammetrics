# doctor: a pre-flight check of one BAM, CRAM or SAM file

`bammetrics doctor` checks one alignment file before a long run. It reads the header, the
last bytes of the file and the first 10,000 records, never the whole file. Each question
is one check line with a status (`pass`, `warn`, `fail` or `skip`), the value found and
the reason. On the 16 GB NA12878 CRAM with its GRCh38 reference it takes 0.35 s.

Source: `crates/core/src/doctor/` (checks), `crates/core/src/genome_build.rs` (build
table), `crates/cli/src/doctor.rs` (options, forms, exit code), `crates/py/src/doctor.rs`
(Python). Plan: docs/07, wave W07, decisions D8 and D13.

```bash
bammetrics doctor in.bam
bammetrics doctor in.cram -R ref.fa --expect-build hg38 --format json
bammetrics doctor in.bam --sample-sheet samples.csv --sample NA12878 --fail-on warn
```

## 1. Options

| Option | Meaning |
|---|---|
| `FILE` | The BAM, CRAM or SAM file (plain, bgzip or gzip). A regular file: doctor reads its end, which a pipe does not have. |
| `-R`, `--reference FASTA` | The FASTA with its `.fai`. A CRAM is decoded with it. Its contig names, lengths and MD5s are checked against `@SQ`, for a BAM too. |
| `--records N` | Records read from the start (default 10,000). 0 reads the header only. |
| `--check-md5` | Compute the MD5 of each `@SQ` contig with an `M5` from the FASTA. Default: read the `.dict` beside it. |
| `--expect-build B` | The build the header should name (section 3). |
| `--expect-sex male\|female` | A warning at most (section 2, `expected_sex`). |
| `--sample-sheet FILE --sample NAME` | Check the sheet's row for NAME against the file (section 4). Each needs the other. |
| `--format text\|json\|jsonl` | The form on standard output (section 5). |
| `--fail-on warn\|fail` | The status that exits 10 (default `fail`). |
| `--threads N` | CRAM container decoders and MD5 workers (default: the cores, at most 8). |

## 2. Checks

The checks come in this order in every form. A check that does not apply is `skip`, with
the reason.

| Check | Pass | Warn | Fail |
|---|---|---|---|
| `file` | BAM, CRAM 3.0 or 3.1, SAM, SAM (bgzip), SAM (gzip); the `@HD VN` | | none of those |
| `eof` | BAM or bgzip: the 28-byte BGZF end-of-file block; CRAM: the 38-byte EOF container; SAM: a final line end | SAM text without a final line end | no end-of-file block (the last block is whole), cut short (the last block is not), no CRAM EOF container |
| `header` | the header lines by type | | the header cannot be read |
| `references` | `@SQ` count and total length; no `@SQ` and no placed record (an unaligned file) | | an `@SQ` without `SN` or a number `LN`, a name twice, a BAM binary reference list that differs from the `@SQ` lines, placed records under no `@SQ` |
| `build` | GRCh37, hg19, GRCh38, hg38 or T2T-CHM13, or `unknown`; the contigs by kind | | primary chromosomes of two builds; `--expect-build` not met |
| `read_groups` | every `@RG` has ID, SM, LB and PL | no `@RG`; an `@RG` without LB or PL | duplicate IDs; an `@RG` without ID or SM |
| `sort_order` | the records read agree with `SO` | the records are in coordinate order but `SO` is absent, `unknown` or `unsorted` | `SO:coordinate` over records out of coordinate order; `SO:queryname` over records out of name order |
| `records` | records read (at most `--records`) | the file holds no record | a read error (a cut record) |
| `record_read_groups` | every `RG` tag is in the header | records without `RG` while the header has `@RG` | records whose `RG` the header lacks |
| `paired` | all paired or none | paired and single-end mixed | |
| `mapped` | at least 50% of primary records mapped; an unaligned file | under 50% | |
| `duplicates` | some flagged, or none and an `@PG` of a duplicate marker | none flagged and no such `@PG`: duplicates were probably not marked | |
| `read_length` | the mean over primary records with bases | no record has bases | |
| `qualities` | the range, Phred+33 | no record has qualities | Phred+64 suspected |
| `md_nm` | MD and NM on every mapped primary record (skip for a CRAM) | some lack them | |
| `reference` | with `-R`: every `@SQ` contig in the `.fai` at its length, the `M5`s equal; without `-R`: a CRAM that decodes without one | without `-R`: a CRAM that needs one | no `.fai`; a contig absent or of another length; an `M5` that differs |
| `index` | `.bai` or `.csi` (BAM), `.crai` (CRAM) beside the file, not older | no index; an index older than the file | |
| `sample_sheet` | the row agrees | | no row, two rows, or a column that disagrees |
| `expected_sex` | male with a Y contig, or female | male expected and no Y contig | |

Thresholds, all starting values:

- `mapped` warns under 50% of the primary records read.
- `qualities` suspects Phred+64 when no quality is under Q31 and some is over Q45. Phred+64
  text read as Phred+33 is 31 too high: Illumina 1.3 to 1.7 qualities Q0 to Q40 become Q31
  to Q71.
- `duplicates` takes an `@PG` whose `ID` or `PN` contains `markduplicates`, `markdup`,
  `samblaster`, `sambamba`, `bammarkduplicates` or `dedup` as a duplicate marker.
- `index` warns on an index older than the file. A download or a copy can write the two in
  that order, so the check does not fail.
- `expected_sex` never fails. 10,000 records of a coordinate-sorted file all lie on the
  first contig, which says nothing about X and Y. `bammetrics sex-check` infers sex from
  every read (docs/06 section 2).

A coordinate-sorted file has its unmapped reads at the end, so `mapped`, `paired` and
`duplicates` describe the start of the first contig. The records are counted, never
judged one by one: a check that reads records reports the first record that breaks a rule
and the counts.

### 2.1 What is read

- The kind, from the first bytes, as every bammetrics reader tells it (`cram::sniff`).
- The end: the last 64 KiB and 28 bytes of a BGZF file, the last 38 bytes of a CRAM, the
  last byte of plain SAM text. A BGZF block that starts in that window and ends at the last
  byte is whole.
- The header text, raw: a BAM header up to its binary reference list, a CRAM's SAM header
  container, the `@` lines of SAM text. No parser refuses it, so a header with two `@RG`
  lines of one ID is reported, not rejected. At most 64 MiB.
- The first records: a BAM's records from the same stream; a CRAM's through the CRAM
  decoder with `-R`; SAM text through the SAM reader at `SILENT`. A CRAM read against a
  FASTA whose contigs differ from `@SQ` is not decoded: the `reference` check names the
  contig.
- With `-R`: the `.fai` (at most 64 MiB) and the `.dict` beside the FASTA (`ref.dict` for
  `ref.fa`, `.fasta`, `.fna`, each with an optional `.gz`; then `ref.fa.dict`).

Nothing is written. A missing file, a directory or a pipe exits 1.

### 2.2 MD5

The `M5` of an `@SQ` line is the MD5 of the contig's bases in upper case, without line
ends. Computing it reads the whole FASTA: 3.3 GB for GRCh38 takes 2.5 s on 8 threads.
So by default doctor compares `M5` with the `M5` of the same contig in the `.dict` beside
the FASTA, which `CreateSequenceDictionary` and `samtools dict` write. Without a `.dict`
the check says the MD5s were not compared. `--check-md5` computes them, one contig per
worker, 1 MiB at a time.

## 3. Builds

`crates/core/src/genome_build.rs` holds the lengths of chromosomes 1 to 22, X and Y in
GRCh37, GRCh38 and T2T-CHM13 v2.0 (the UCSC `chrom.sizes` of hg19, hg38 and hs1). No two
builds share a length for any of them, so each primary chromosome present names its
build. No file or network is read.

- A name is compared without `chr`; `MT` is `M`.
- The file is of a build when every primary chromosome present has that build's length.
  One chromosome is enough; a header with chromosome 21 alone is named.
- Most names with `chr` gives `hg19` or `hg38`, without it `GRCh37` or `GRCh38`.
  T2T-CHM13 has one name.
- Alt (`_alt`), decoy (`_decoy`, `hs37d5`), HLA (`HLA-`), mitochondrial and other contigs
  are counted and never decide the build.
- Primary chromosomes of two builds fail the check: the reason names how many agree with
  the likeliest build and each that does not.
- No primary chromosome of a known length gives `unknown` (another species, or a test
  reference), which passes.

`--expect-build` and the sheet's `build` column take `GRCh37`, `b37`, `hs37d5`, `hg19`,
`GRCh38`, `hg38`, `T2T-CHM13`, `T2T`, `CHM13`, `chm13v2`, `hs1`, in any case. `hg19` and
`hg38` expect `chr` names, `GRCh37` and `GRCh38` names without it; a T2T name takes
either. A mismatch fails `build` and names what was found:
`expected hg38, found GRCh38 (contig names without chr)`.

## 4. Sample sheet

`--sample-sheet FILE --sample NAME` checks one sample per call (docs/07 D13). The sheet is
CSV or TSV (a tab in the first line makes it TSV), at most 16 MiB, read by
`crates/formats/src/sheet.rs`. Its header row names the columns in any case: `sample` is
required, `bam`, `sex`, `build` and `library` are compared when present and the cell is
not empty or `NA`. Other columns are ignored.

```
sample,bam,sex,build,library
NA12878,NA12878.bam,female,hg38,lib1
```

`sample_sheet` fails when:

- no row or two rows name the sample;
- no `@RG` has `SM` NAME;
- the `build` cell is no build name, or the header is of another build;
- the `bam` cell names another file (resolved against the sheet's directory, compared by
  device and inode);
- no `@RG` has `LB` equal to the `library` cell;
- the `sex` cell is not male, female (`m`, `f`, `xy`, `xx`) or `unknown`.

The `sex` cell feeds `expected_sex` when `--expect-sex` is not given. A sheet without a
`sample` column, or that cannot be read, exits 1.

## 5. Output and exit codes

Text: the file line (path, kind, status, counts), one line per check, the exit line.

```
tests/data/doctor/cut.bam  BAM  FAIL  2 fail, 1 warn, 13 pass, 3 skip
  PASS  file                BAM, @HD VN:1.6
  FAIL  eof                 cut short                        the last BGZF block is incomplete: ...
  FAIL  records             162 records read                 record 163: the file ends inside the record
  ...
exit 10: a check reaches --fail-on fail
```

`--format json` prints one object, schema `bammetrics.doctor/1`: `schema`, `path`, `kind`,
`status`, `exit_code`, `records` (the limit), `counts` (`fail`, `warn`, `pass`, `skip`) and
`checks`, each with `check`, `status`, `value` and `reason`. `--format jsonl` prints one
`check` event per line, then a `doctor` event with the status, exit code and counts.
Both go through `crates/formats/src/json.rs`.

| Code | Meaning |
|---|---|
| 0 | no check reaches `--fail-on` |
| 1 | the file, the `-R` FASTA or the sample sheet could not be read |
| 2 | usage error, nothing ran |
| 10 | a check fails, or warns with `--fail-on warn` |

The codes are those of `check` (docs/cli.md, "Exit codes"). Code 11 means an incomplete
verdict there; doctor has no verdict and never exits 11.

## 6. Python

```python
import bammetrics as bm

result = bm.doctor("in.cram", reference="ref.fa", expect_build="hg38")
result.status            # "pass", "warn" or "fail"
result.exit_code         # 0 or 10, under fail_on
result["eof"].status     # one DoctorCheck: name, status, value, reason
[c.name for c in result.checks if c.status == "fail"]
result.to_json()         # the text of --format json
```

The keywords are the command's long options in snake case. A failed check is a value. A
file that cannot be read raises `InputError`, a bad keyword `UsageError`. `Doctor` and
`DoctorCheck` pickle.

## 7. Tests

No other tool answers these questions in this form, so correctness rests on planted
defects. `tests/data/make_doctor_fixtures.py` writes `tests/data/doctor/`: `good.bam` and
`good.sam`, then one file per defect, each different from `good.bam` in one way (a cut
BAM, no end-of-file block, coordinate over unsorted records, an RG the header lacks, two
`@RG` with one ID, Phred+64 qualities, hg19, GRCh38, T2T and mixed headers, an unaligned
file, a CRAM and a FASTA of other lengths, a sample sheet in CSV and TSV).
`crates/cli/tests/doctor.rs` checks the status of each check on each of them, and on files
it makes: a CRAM without its EOF container, SAM text without a final line end, a BAM cut
inside its header, a header naming a contig twice, mostly unmapped records without bases
or qualities, an index older and one newer than its BAM, a FASTA with one base changed and
a `.dict` with a wrong `M5`. A test fails when a check has no pass case or no warn or fail
case. `tests/golden/doctor_*.txt` hold the text, JSON and JSON Lines forms;
`tests/python/test_doctor.py` the Python face. The build table was compared with the
`.fai` of GRCh38 and hs37d5 on the scale data.

## 8. Measurements

`/usr/bin/time -l`, release build, Apple M1 Ultra (20 cores), warm cache, 2026-11-17:

| Run | Real | Peak RSS |
|---|---|---|
| NA12878.final.cram (16 GB, CRAM 3.0) `-R GRCh38.fa`, median of 5 | 0.35 s | 666 MB |
| the same with `--threads 1` | 0.34 s | 557 MB |
| the same with `--check-md5` (3,366 contigs, 3.3 GB) | 2.52 s | 673 MB |
| the same without `-R` (records skipped) | 0.04 s | 85 MB |
| HG002.exome.bam (BAM) `-R hs37d5.fa` | 0.03 s | 7 MB |

The CRAM decoder loads the first contig whole (chr1, 249 MB), which sets the memory.

# In-process chaining

A chain runs several steps on one BAM in one process. Records move from step to step in
memory. No shell pipe and no intermediate file is used. This page describes what waves W03
to W05 of docs/07_full_parity_program.md (section 2) built: the stage layer and its
executor (W03, 0.21.0), then the `chain` verb, its plan file, its plan text and the Python
`chain` and `step` (W04, 0.22.0), then the sort barrier, the spool, the planner and SAM
text in and out (W05, 0.23.0). Sections 22 to 32 hold chain 2 (0.32.0): steps without
`+`, `-o` naming the BAM, the config file and the recipes (22 to 26), the manifest,
`--resume`, `--keep` and `qc` as the last step (27 to 30), every tool that reads a BAM as
a step (31), and the three together (32).

```bash
# 0.32.0: a bare word starts a step, -o names the BAM, the reference comes from the config file
bammetrics chain in.bam mark-duplicates insert-size -o out/NA12878.bam
bammetrics prep in.bam -o out/NA12878.bam qc --fail-on fail --resume    # a recipe, the verdict, resumable
```

The long form puts a `+` in front of every step and names each file. Every form of this
page written with `+` runs as it did at 0.31.0, byte for byte:

```bash
# MarkDuplicates, four collectors on the marked records and the marked BAM, one process
bammetrics chain in.bam -R ref.fa -o out/NA12878 \
    + mark-duplicates \
    + alignment-summary + insert-size + quality-yield + duplicate-metrics \
    + write --index
```

This writes `out/NA12878.marked.duplication_metrics`, the four metrics files
`out/NA12878.alignment_summary_metrics`, `.insert_size_metrics`, `.quality_yield_metrics`
and `.duplicate_metrics`, and `out/NA12878.marked.bam` with `out/NA12878.marked.bai`. The
four metrics files equal, byte for byte after the `StringHeader` lines, the files the four
tools write when they read the marked BAM. The BAM holds the records the `MarkDuplicates`
tool writes.

## 1. Steps

A chain is a source and a list of steps. Each step sees the records the step before it
emits.

| Step | What it is | Today |
|---|---|---|
| source | One BAM, CRAM or SAM text file, or a pipe holding one | the positional argument, `-` for standard input |
| transform | A tool that changes the records later steps see. A barrier needs every record before one goes on; a streaming transform edits each batch as it passes (section 16) | `mark-duplicates`, `sort-sam`, `revert-original-base-qualities-and-add-mate-cigar` (barriers); `add-comments-to-bam`, `replace-sam-header`, `add-o-a-tag`, `clean-sam`, `set-nm-md-and-uq-tags`, `set-nm-and-uq-tags` (streaming); `add-or-replace-read-groups` (streaming, a barrier with a `--sort-order` that differs from the stream's; section 17); `reorder-sam` (a barrier; section 19); `fix-mate-information` (a barrier), `filter-sam-reads` (streaming, a barrier when it sorts; section 20); `downsample-sam` (streaming; section 21) |
| tap | A one-pass tool on the records of its position. It changes nothing | 20 verbs: `insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics`, `gc-bias`, `quality-score-distribution`, `mean-quality-by-cycle`, `base-distribution-by-cycle`, `duplicate-metrics`, `rna-seq`, `sequencing-artifact`, `coverage-gaps`, `sex-check`, `mt-copy-number`, `exon-coverage`, `adapter-content`, `duplicate-by-tile`, `telomere`, `quality-yield-s-n-v-q`, `rrbs` |
| write | A tap that writes the records of its position as a BAM, as SAM text, or as a CRAM against `-R` | `write [PATH]` |
| split | A sink that writes the records of its position over BAMs in a directory (section 18) | `split-sam-by-library`, `split-sam-by-number-of-reads` |
| materialize | A tool with neither a stage nor a collector, run on a file of the records of its position: the source itself, or a copy in the temporary directory (section 31). One that writes a BAM hands it on as the next stream | 17 verbs: `wgs`, `raw-wgs`, `wgs-with-non-zero-coverage`, `targeted-pcr`, `library-complexity`, `multiple-metrics`, `oxo-g`, `bam-index-stats`, `index-depth`, `check-duplicate-marking`, `jumping-library`, `view-sam`, `calculate-read-group-checksum`, `umi-prevalence`; `merge-sam-files`, `gather-bam-files`, `position-based-downsample-sam` (each leaves a stream) |

A tap or `write` does not change the stream, so the step after it sees the same records.
`in.bam + tap + mark-duplicates + tap + write` gives the first tap the source records and
the second tap and the writer the marked records.

Every registered tool that reads a BAM is a step at any position, 53 of the 74 (section
31), but the four named below. The 17 tools that read no BAM (interval, VCF, FASTA and text tools,
`compare-metrics`, `fifo-buffer`) are refused with exit 2 and the reason. So is
`build-bam-index`: it names `write --index`, which writes the `.bai` of the BAM a chain
writes. So is `sam-format-converter` (wave W15): it writes the records it reads in the
format its name gives, which is what `write out.sam`, `write out.bam` and `write
out.cram` (with `-R`) do. The refusal names `write`.

`fifo-buffer` is never a step. `FifoBuffer` copies the bytes of a pipe between two
commands. The steps of a chain hand records on in memory, so a chain has no pipe for it
to buffer. The chain refuses it with exit 2 before the first read, in every spelling, on
the command line, in a plan file, in `plan` and in Python:

```
$ bammetrics chain in.bam + mark-duplicates + fifo-buffer --buffer-size 1024 + write -
bammetrics chain: step 2 (fifo-buffer): FifoBuffer buffers the pipe between two commands and a chain has no pipe: it cannot be a step
  leave the step out, or put the buffer on a pipe outside the chain: bammetrics chain BAM + write - | bammetrics fifo-buffer | ...
  the steps: bammetrics chain --help
```

On a pipe into or out of a chain the verb works as on any pipe:
`bammetrics chain in.bam + mark-duplicates + write - | bammetrics fifo-buffer | samtools
flagstat -`, or `producer | bammetrics fifo-buffer | bammetrics chain - + insert-size`.
`crates/cli/tests/chain.rs` checks both: the bytes behind the buffer equal the BAM the
chain writes to a file.

`compare-s-a-ms` and `check-terminator-block` are never steps, each with its own reason
(exit 2). `CompareSAMs` reads two files and a chain has one stream: write the stream with
`write out.bam`, then run `bammetrics compare-s-a-ms BAM out.bam`. `CheckTerminatorBlock`
reads the last bytes of a file and a stream has none: run it on the file the chain wrote.
A spooled copy of the stream gives neither tool its input, so no later wave makes either
a step.

## 2. The chain verb

```
bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]
bammetrics chain [OPTIONS] BAM STEP [STEP OPTIONS] [STEP ...] [OPTIONS]      (section 22)
bammetrics chain --plan FILE [OPTIONS] [BAM]
bammetrics chain BAM --recipe NAME [OPTIONS] [STEP.option=value ...] [qc [QC OPTIONS]]   (section 25)
```

This section describes the `+` form. The chain's options and the source come first. Every `+` starts a step. The first word of
a step is its verb: a tool verb (`insert-size`), its Picard class
(`CollectInsertSizeMetrics`), or `write`. The words after it are that step's options.

A step is parsed by its own verb's flags (`Verb::flags`). So every spelling a verb takes
works in a step: `--minimum-pct 0.1`, `--MINIMUM_PCT 0.1`, `MINIMUM_PCT=0.1`, `-M 0.1`,
`--no-x` for a boolean. An unknown option gets the closest name of that verb's options.
`bammetrics chain BAM + VERB --help` lists the options of one step. A value that is a lone
`+` is written `--option=+`.

### 2.1 Chain options

| Option | Meaning |
|---|---|
| `BAM`, `-i`, `I=` | The source: a BAM, or a CRAM with `-R`. `-` reads standard input |
| `-o`, `--output-prefix`, `O=` | The prefix of every step that names no file of its own. Default `./` plus the source name without its extension; with a `write PATH` step, `PATH` without its `.bam` or `.sam` (`write d/out.bam` gives `d/out`); else with a split step's `-o DIR`, `DIR/` plus the source name. Standard error says which. A directory (an existing one, or a path ending in `/`) gets the source name as the file part |
| `-R`, `--reference`, `R=` | Reference FASTA with a `.fai`. It decodes a CRAM source. It goes to every step that reads a reference and names none |
| `-@`, `--threads` | BGZF threads of each read and of `write` (default: cores, at most 8; at most 256 may be asked for) |
| `--tmp-dir`, `TMP_DIR=` | Where a transform spills and a piped source is copied |
| `--format LIST` | `picard` first. Each further format is written beside every metrics file |
| `--also LIST` | `json`, `jsonl`, `tsv`, `csv`, `parquet`, `xlsx`, `docx`, `sql` beside every metrics file; `markdown` as `PREFIX.report.md` |
| `--sample NAME` | Sample name of the json and jsonl forms (default: the first `@RG SM`, else the source name) |
| `--sort-memory BYTES` | What the barriers of the chain hold: a `sort-sam` step spills a sorted run at it, less 64 MiB and the read ends of a later `mark-duplicates`. A number, or one with K, M, G or T (1024 based), at most 18446744073709551615 bytes (exit 2 past it). Default 2G (section 15) |
| `--force-space` | Run when the estimated temporary and output bytes exceed the free space of their file system (section 15) |
| `--validation-stringency` | `STRICT` (default), `LENIENT` or `SILENT` for a SAM text source. A step refuses it |
| `--resume` | Skip each step whose record in `PREFIX.manifest.json` matches and whose transforms in front are skipped too (section 27). Standard input as the source is a usage error |
| `--plan FILE` | Read the source, the options and the steps from a plan file (section 3) |
| `--recipe NAME` | Run a named chain; `STEP.option=value` sets an option of one of its steps, and a last `qc` step may follow (section 25) |
| `--no-config` | Leave the user config file unread (section 24) |
| `--dry-run` | Print the plan and run nothing (section 4) |
| `-q`, `-v` | No notes and no summary; or each command line and each file written |

The source, `--output-prefix`, `--threads`, `--sample`, `--also`, `--plan`, `--dry-run`,
`-q` and `-v` belong to the chain. Given after the first `+` they are a usage error that
says so: `step 1 (insert-size): --also belongs to the chain: give it before the first +`.
One exception: `--dry-run` after the last step, when that step has no flag of that
spelling, is the chain's (`chain BAM + quality-yield + write x.bam --index --dry-run`).
`write` has none of the chain's options, so every one of them gets that error there.

An empty prefix (`-o ''`) is a usage error, as in `qc`. It would name hidden files in the
working directory (`.insert_size_metrics`).

### 2.2 Step options

- A tap takes every option of its verb. `-o FILE` names its file (`-o -` writes it to
  standard output). On a verb that writes several files (`gc-bias`,
  `sequencing-artifact`) `-o` is a prefix for that step alone.
- `mark-duplicates` takes every option of the tool except `OUTPUT`, `COMPRESSION_LEVEL`,
  `CREATE_INDEX` and `CREATE_MD5_FILE`. Those belong to `write`. `--metrics-file` names
  its metrics file. `--as NAME` names the stream it leaves (default `marked`).
- `sort-sam` takes `--sort-order` (required), `--max-records-in-ram` (a record bound beside
  `--sort-memory`; none by default in a chain), `--tmp-dir` and `-R`. `OUTPUT` and the
  options of the BAM file belong to `write`. Its stream is `sorted`. A default name an
  earlier transform took gets `-2`, `-3` (`sorted-2`); an `--as` name that is taken is a
  usage error.
- `write [PATH]` takes `--index` (`CREATE_INDEX`), `--md5` (`CREATE_MD5_FILE`) and
  `--compression-level N` (0 to 9, default 5). `-` as the path is standard output, so a
  pipe to another tool still works: `bammetrics chain in.bam + mark-duplicates + write -
  | samtools flagstat -`. The file is a BAM, SAM text for a name ending in `.sam`, or a
  CRAM against the chain's `-R` for a name ending in `.cram`, lower case as htsjdk tells
  them (`out.SAM` and `out.CRAM` are BAMs). `--output-format bam|sam|cram` says which
  whatever the name is, so `write - --output-format sam` prints SAM text. A CRAM needs
  `-R` and takes no `--index`: either is a usage error before the first read
  (KNOWN_DIFFERENCES.md, CRAM output). A name ending in `.fastq` or `.fq` (with `.gz`
  too) is refused with the wave that adds it, W22. `.sam.gz` is refused: SAM text is
  written uncompressed. Section 14 has the details.
- `--index` writes the `.bai` htsjdk names: `out.bam` gives `out.bai`, `out.dat` gives
  `out.dat.bai`. It writes one only for a coordinate-sorted stream. For any other sort
  order `write` writes no index and says so on standard error. The exit code stays 0.
  The `.bai` is the file Picard `BuildBamIndex` writes for the written BAM, byte for byte
  (W05, `crates/core/src/bai.rs`). `build-bam-index` is no step: `chain BAM +
  build-bam-index` exits 2 and names `write --index`.
- A file option with an empty value (`-o ''`, `--metrics-file ''`, `write ''`) is a usage
  error: `step 1 (gc-bias): -o names no file: its value is empty`.
- `--keep` on a transform writes the stream it leaves as `PREFIX.<stream>.bam` (section
  28). A materialized tool that writes a BAM takes it too. On any other step it is a
  usage error.
- `qc` is the last step: the verdict and report over the taps of its stream (section 29).

### 2.3 Files and names

Names follow `qc`. The extensions come from `multiple_metrics::file_extensions`, the one
table `qc`, `collect` and the verbs read. A stream has a name: `in` for the source, the
label of the transform in front otherwise.

| Step | File |
|---|---|
| A tap on the last stream | `PREFIX.<extension>`, the name `qc` writes |
| A tap in front of a transform | `PREFIX.<stream>.<extension>`: `PREFIX.in.insert_size_metrics` on the source |
| `mark-duplicates` | `PREFIX.<stream>.duplication_metrics`: `PREFIX.marked.duplication_metrics` |
| `write` without a path | `PREFIX.<stream>.bam`: `PREFIX.marked.bam`, or `PREFIX.in.bam` in a chain without a transform |
| A step with its own `-o`, or `write PATH` | That name |

```
$ bammetrics chain tests/data/markdup.bam -o out/s --also csv \
      + quality-yield \
      + mark-duplicates --as dedup --remove-duplicates \
      + insert-size \
      + quality-yield -o out/after.txt \
      + write out/kept.bam --md5
bammetrics chain: step 2 (mark-duplicates): 36 records in, 19 out, 16 changed
bammetrics chain: 5 steps finished, 10 files written
$ ls out
after.txt
after.txt.csv
kept.bam
kept.bam.md5
s.dedup.duplication_metrics
s.dedup.duplication_metrics.csv
s.in.quality_yield_metrics
s.in.quality_yield_metrics.csv
s.insert_size_metrics
s.insert_size_metrics.csv
```

Rules:

- Every output needs its own file. The outputs are each tap's files, the metrics file of
  a transform, the BAM of `write` with its `.bai` and `.md5`, and each further format
  beside a metrics file. Two outputs that are one file (links and `..` resolved, hard
  links by device and inode) are a usage error:
  `step 1 (quality-yield) and step 2 (insert-size) both write x.out: every output of a chain needs its own file`.
  A device or a pipe is exempt: several steps may name `/dev/null`.
- Two taps of one tool at one position write one name under the prefix, so they are that
  usage error. Give one of them `-o`.
- A tool whose `OUTPUT` is a prefix (`sequencing-artifact`) names its files when it
  finishes. If one of them is the file of another step, or a file the chain reads, that
  tool fails alone, writes none of its files, and the other file stays.
- `rrbs` is the second tap whose `OUTPUT` is a prefix. It names both files
  (`PREFIX.rrbs_summary_metrics`, `PREFIX.rrbs_detail_metrics`) before the first read,
  so they get the checks of every other output.
- `create-baf-regress-file`, `create-verify-i-d-intensity-contamination-file` and
  `arrays-variant-calling` read a text file or a VCF and no BAM. A chain refuses them as
  steps with exit 2.
- `scatter-intervals-by-ns`, `non-n-fasta-size` and `create-sequence-dictionary` read a
  FASTA and no BAM. A chain refuses them as steps with exit 2 and the same reason.
- An output that is a file the chain reads (the source, the reference and its `.fai`, an
  interval list of a step, the plan file of `--plan`) is a usage error. Nothing is read
  and nothing is written. The same file through a symbolic link, a hard link or `..`
  counts: `step 1 (insert-size): plan.toml is an input; writing it would destroy it`.
- Missing directories above an output are created, as every verb creates them.
- At most one step writes standard output (`write -`, or `-o -` on a single-file tap).

### 2.4 Formats

The Picard text format is the default and its bytes do not change. `--format picard,json`
writes each metrics file and `FILE.json` beside it, from one set of header lines, through
`crates/formats`. A step's own `--format` replaces the chain's list for that step. The list
starts with `picard`, on the chain and on a step: a chain names Picard files. `--also
json|jsonl|tsv|csv|parquet|xlsx|docx|sql` adds the same siblings. `--also markdown` writes `PREFIX.report.md`
with every metrics table the run wrote, in step order. `--also junit` and `--also multiqc`
print a verdict and a chain has none: they are a usage error that names `bammetrics check
DIR --format junit`. A step that writes standard output writes the Picard file alone.

The json, tsv and csv siblings are written from the values in memory, as `--format
picard,json` of a tool verb writes them, so a double keeps its full precision. `qc --also`
reads each Picard file back and so writes the 6 fraction digits of that file.

On a chain `--format` names the formats of the metrics files. `qc`, and `plan` with `qc`
arguments, give the same flag another meaning: `--format text|jsonl` there picks how
progress and the plan are printed. So `bammetrics plan BAM --format jsonl + STEP` is a
chain and is refused with exit 2: `--format jsonl starts with jsonl: a chain names Picard
files`. A chain has no event stream in this version. What it gives in its place:

- Standard error carries the lines of section 2.5. `-v` adds each step's command line and
  each file written.
- `--dry-run` and `plan` print the plan as text (section 4).
- `bm.chain` returns `Run.stages`: one entry per step with its status, error, files and
  record counts (section 5.3).
- `--also jsonl` writes each metrics file as JSON Lines beside the Picard file. That is a
  file format, with no events in it.

`qc --format jsonl` stays the event stream of a `qc` run.

### 2.5 What it prints, and exit codes

Standard output is empty unless a step writes there. Standard error gets one line per
transform with its record counts, each note, each failed step, and one summary line:

```
bammetrics chain: step 1 (mark-duplicates): 36 records in, 36 out, 16 changed
bammetrics chain: 6 steps finished, 7 files written
```

A message names a step by its number and verb: `step 3 (insert-size)`.

| Code | When | Found |
|---|---|---|
| 0 | Every step finished | |
| 2 | A usage error: an unknown step or option, a chain option on a step, a value that is no number, a missing required option, an empty prefix or file name, a name `write` cannot write yet, two outputs in one file, an output that is an input, a stream name another transform took, a bad `--sort-memory`, a bad plan file (one that is not UTF-8 text included) | Before the first read. No file is written and no directory is made |
| 1 | An input that is not there or cannot be read (the source, `-R`, the plan file) | Before the first read |
| 1 | A file system with less free space than the estimate needs there, without `--force-space` | Before the first read. No directory is made |
| 1 | An output that cannot be written (its directory cannot be made, the file is read only) | Before the first read, as the Picard face exits for its `OUTPUT` |
| 1 | A step failed | During the run. The other steps still write their files (section 9) |
| 10, 11 | The qc step's verdict reached its `--fail-on`, or is incomplete, and every step finished | After the run, as `check` gives them (section 29) |

```
$ bammetrics chain in.bam + quality-yield + insert-size --minmum-pct 0.1
bammetrics chain: step 2 (insert-size): unknown option --minmum-pct
  did you mean --minimum-pct (Picard MINIMUM_PCT)?
  every option: bammetrics help insert-size
$ bammetrics chain in.bam + inser-size
bammetrics chain: step 1: unknown step "inser-size"
  did you mean insert-size?
  a step is a tool verb (bammetrics tools lists them) or write
$ bammetrics chain in.bam + write out/x.cram
bammetrics chain: step 1 (write): out/x.cram: a CRAM is written against the reference: give the chain -R FASTA
```

What is checked before the first read: the step name, every option name, every value of
an int or float option, required options, files a step reads, `STOP_AFTER`, every option
of `mark-duplicates` (`stage::check_stage_options`), the names and the writability of
every output. Section 2.6 adds six checks.

### 2.6 More that is refused before the first read

The 0.22.0 release review found these. Tests: `crates/cli/tests/chain_safety.rs` (10
tests), `tests/python/test_chain_safety.py` and 13 golden records. The tests show "before
the first read" on a source cut after 150 bytes: its header block is whole, the block of
its records is not, and a run that reads it exits 1 with `truncated BAM`.

| What | The message says | Exit |
|---|---|---|
| A text value only the tool can judge on a tap (`--level NOPE`, `--expected-sex NOPE`) | The tool's own words under the step's label: `step 2 (insert-size): METRIC_ACCUMULATION_LEVEL "NOPE": expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP` | 2 |
| A number the tool refuses without a header (`sequencing-artifact --context-size 99`, `--maximum-insert-size 10`, `gc-bias --scan-window-size 0`, `quality-yield --flow-mode`) | The tool's own words: `step 1 (sequencing-artifact): CONTEXT_SIZE=99: too many contexts to count` | 2 |
| `--threads` over 256, or `threads` over 256 in a plan file | `100000000 is more than the maximum of 256 threads` | 2 |
| An output through a directory that does not exist with a `..` behind it (`missing/../in.bam`) | `missing does not exist, so the .. behind it names no directory` | 2 |
| A plan file over 1 MiB (1,048,576 bytes) | `not a plan file: it is larger than 1 MiB` | 2 |
| An output that is a symbolic link, pipe, socket or device another user owns | `owned by another user; not written in place`, the words and the exit code of the tool verb | 1 |

Values. `bammetrics_core::tap_options::check` holds one row per one-pass tool (20). Each
row calls the option parser that the tool's own collector calls, with no header and no
file opened. A unit test fails when a one-pass tool has no row. In 0.22.0 the tap failed
alone when it was built, after the transform in front of it had read the whole source,
and the chain exited 1. A value the tool judges against the header of its stream, or
against a file it reads, is still found when the tap is built, and that tap fails alone
with exit 1: a sort order it requires (`--no-assume-sorted` on a queryname stream), an
interval list with a contig the header lacks. The reason: the header of a piped source
can be read once, and the header behind a transform exists only after the transform ran.

Three rows run more than the parser, because the tool checks more without a header.
`sequencing-artifact` runs the stand-alone validation (`CONTEXT_SIZE` and the insert
sizes not negative, the maximum not under the minimum, each `CONTEXTS_TO_PRINT` of the
implied length) and the context count (at most 4^11, `CONTEXT_SIZE` 5). The verb ran
both in `run_tool`. The collector a chain builds ran the count alone, when the tap was
built (exit 1), and never the validation: `chain ... + sequencing-artifact
--maximum-insert-size 10` wrote its five files where the verb exits 2.
`CollectMultipleMetrics` still skips the validation, as Picard does. `quality-yield
--flow-mode` is refused by `quality_yield::checked_options`, which the collector calls
too. `gc-bias --scan-window-size` under 1 is refused by the option parser: 0 was a
division by zero and a negative window an index past the contig, a panic with exit 101
on the verb, the Picard face and the chain. Picard 3.5.0 throws on both
(`tests/concordance/test_gc_bias.py::test_a_window_under_one_base_fails`).

Threads. `pass::check_threads` is the one check and `pass::MAX_THREADS` is 256. Every
place that parses a thread count asks it: the tool verbs, `chain`, the plan file, `qc`,
`plan`, `collect`, `threads=` on the Picard face (exit 1 there, as every error of that
face) and the Python `chain`, `qc` and `collect`. In 0.22.0 and before,
`--threads 100000000` started threads for about two minutes and then panicked inside
`noodles-bgzf` with exit 101, on every verb. It now exits 2 in 0.04 s. The maximum runs:
`quality-yield` on the 10.3 GB HG002 exome BAM takes 33.1 s at 8 threads (7.8 MB peak)
and 30.4 s at 256 (48.8 MB peak), with the same metrics (Mac Studio, 20 cores).

Paths. A run makes the directories that are missing above its outputs. `missing/../in.bam`
names no file before the run and names `in.bam` once `missing` is made, so neither the
one-file-twice check nor the input check saw the file. In 0.22.0 `write missing/../in.bam`
failed at write time with exit 1 and left `missing/` behind, and a tap with that `-o`
replaced the source BAM with its metrics file. The tool verbs had the same hole
(`insert-size in.bam -o missing/../in.bam`) and refuse the path with the same words.
`qc`, `plan`, `collect` and the Python `qc` and `collect(metrics=[...])` refuse such a
prefix too (`qc in.bam -o missing/../q`). Their extensions are fixed, so no input could
be overwritten: the run wrote `q.*` and left `missing/` empty. A
`..` behind a directory that exists resolves and is checked as any spelling. Two
spellings of one new file under a directory that does not exist yet (`link/new/a.txt`
and `new/a.txt`, `link` leading to the working directory) are one file: the deepest
directory that exists is resolved (`mark_duplicates::same_file_key`).

Directories. `cleanup::MadeDirs` records each directory the run makes above an output. A
run that ends with exit 1, before the first read or after it, removes those that are
still empty, deepest first. A directory that was there before, or that holds a file,
stays. A signal removes them under the same rule. A run that finishes keeps every one.
The tool verbs, `collect` and `qc` hold a `MadeDirs` too: one that fails on its read
(`insert-size cut.bam -o made/deeper/x.txt`) leaves no `made/`. The Picard face makes no
directory.

Files read whole. `crates/cli/src/chain` reads two files whole, each through
`chain::read_capped`, which stops one byte past its cap. The plan file holds at most
1 MiB (`PLAN_FILE_CAP`). `--also markdown` reads back each metrics file the run wrote, at
most 64 MiB each (`REPORT_FILE_CAP`); over that the report fails with exit 1 after every
other file is written.

The help of a step that is no verb names the step by its number, as a run does:
`chain BAM + quality-yield + inser-size --help` prints `step 2: unknown step "inser-size"`.

## 3. The plan file

`bammetrics chain --plan FILE` reads the source, the options and the steps from a TOML
file. It is the same chain: `chain --plan FILE` and the command line it stands for write
the same files and print the same plan.

```toml
# tests/data/chain_markdup.toml
[chain]
source = "tests/data/markdup.bam"
reference = "tests/data/ref_100k.fa"

[[step]]
verb = "mark-duplicates"

[[step]]
verb = "alignment-summary"

[[step]]
verb = "insert-size"

[[step]]
verb = "quality-yield"

[[step]]
verb = "duplicate-metrics"

[[step]]
verb = "write"
index = true
```

```bash
bammetrics chain --plan tests/data/chain_markdup.toml -o out/sample
# the same chain:
bammetrics chain tests/data/markdup.bam -R tests/data/ref_100k.fa -o out/sample \
    + mark-duplicates \
    + alignment-summary + insert-size + quality-yield + duplicate-metrics \
    + write --index
```

The second example is a chain without a transform: four collectors on one read, each
metrics file also as JSON.

```toml
# tests/data/chain_taps.toml
[chain]
source = "tests/data/insert_size.bam"
format = ["picard", "json"]
sample = "S1"

[[step]]
verb = "quality-yield"
use-original-qualities = false

[[step]]
verb = "alignment-summary"
level = ["SAMPLE", "LIBRARY"]
stop-after = 8

[[step]]
verb = "insert-size"
minimum-pct = 0.1

[[step]]
verb = "duplicate-metrics"
```

```bash
bammetrics chain --plan tests/data/chain_taps.toml -o out/sample
# the same chain:
bammetrics chain tests/data/insert_size.bam -o out/sample --format picard,json --sample S1 \
    + quality-yield --no-use-original-qualities \
    + alignment-summary --level SAMPLE --level LIBRARY --stop-after 8 \
    + insert-size --minimum-pct 0.1 \
    + duplicate-metrics
```

The format:

- One `[chain]` table. Its 13 keys: `source`, `output-prefix`, `reference`, `threads`,
  `tmp-dir`, `format`, `also`, `sample`, `quiet`, `verbose`, `sort-memory`,
  `force-space`, `validation-stringency`. `format` and `also` take a
  string (`"picard,json"`) or a list of strings. An empty `output-prefix` is an error.
- One `[[step]]` table per step, in chain order. `verb` is required: a tool verb, a
  Picard class, or `write`. Every other key is the long name of an option of that verb
  (`minimum-pct`), or its Picard name (`MINIMUM_PCT`). A boolean option takes `true` or
  `false`, an int option an integer, a float option a number, a text or path option a
  string. One number also stands for the text of a text option (`program-record-id = 7`).
  A repeatable option takes a list; a list of a text option holds strings only. The keys
  of `write` are `output`, `index`, `md5`, `compression-level` and `output-format`.
- A path is read as the command line reads it: relative to the working directory of the
  command. It is never read relative to the directory of the plan file.
- The plan file is a file the chain reads. A step that names it as an output is a usage
  error (section 2.3).
- An option on the command line replaces the file's: `chain --plan FILE -o out/other
  other.bam` runs the file's steps on another source under another prefix. A `+` step
  beside `--plan` is a usage error.
- An unknown table, an unknown key, a key given twice and a value of the wrong type are
  usage errors with the file, the line and the closest name. Two spellings of one option
  in one step are that option given twice:
  `plan.toml:7: step 1 (insert-size): MINIMUM_PCT is given twice: line 6 spells it minimum-pct`.
  A chain option in a `[[step]]` says where it goes:
  `plan.toml:5: step 1 (insert-size): also belongs to the chain: write also under [chain]`.
- A file that cannot be read exits 1. A file that is read and is no plan exits 2; a file
  that is not UTF-8 text is one of those.

```
$ bammetrics chain --plan tests/data/chain_bad_key.toml
bammetrics chain: tests/data/chain_bad_key.toml:8: step 1 (insert-size): unknown key minmum-pct
  did you mean minimum-pct (Picard MINIMUM_PCT)?
  every option: bammetrics chain BAM + insert-size --help
```

The reader is `bammetrics_core::toml::parse_tables`, the subset of TOML the verdict
presets use (it moved from `crates/verdict/src/table.rs` in this wave). It gained one-line
arrays of scalars for repeatable options. No dependency was added. Every file it reads is
valid TOML. It does not read all of TOML:

| Read | Not read |
|---|---|
| `# comment`, on a line of its own or after a value | |
| `[table]` and `[[table]]` headers | Dotted names (`[a.b]`), inline tables (`{ a = 1 }`) |
| Bare keys: letters, digits, `_`, `-` | Quoted and dotted keys |
| Basic strings in double quotes on one line, with every TOML escape: `\"`, `\\`, `\b`, `\t`, `\n`, `\f`, `\r`, `\uXXXX`, `\UXXXXXXXX` | Literal strings in single quotes, multi-line strings (`"""`, `'''`) |
| Integers, decimals (`0.5`, `1e3`), `true`, `false` | `inf`, `nan`, a decimal too large for a double, dates and times, hexadecimal, octal and binary integers, `_` between digits |
| One-line arrays of those values (`["a", "b"]`), a trailing comma allowed | Arrays over several lines, arrays in arrays |
| A UTF-8 byte order mark in front of the first line | Any other encoding |

A value in a form of the right column is an error that names the form:

```
bammetrics chain: plan.toml:3: also: expected a quoted string, a number, true or false, got "[" (multi-line arrays are not read: write the array on one line)
bammetrics chain: plan.toml:2: source: expected a quoted string, a number, true or false, got "'in.bam'" (literal strings in single quotes are not read: use double quotes)
```

## 4. The plan

`chain ... --dry-run` prints what the chain would do and runs nothing. `bammetrics plan`
with a chain's arguments prints the same text: `plan BAM + STEP ...`, `plan --plan FILE`,
or `plan chain ...`. (`plan` with `qc` arguments, which hold no `+` and no `--plan`, is
still the dry run of `qc`.)

```
$ bammetrics chain tests/data/markdup.bam -R tests/data/ref_100k.fa -o out/sample --dry-run \
      + mark-duplicates \
      + alignment-summary + insert-size + quality-yield + duplicate-metrics \
      + write --index
chain plan for tests/data/markdup.bam
  source    tests/data/markdup.bam (BAM, 1,158 bytes, sort order coordinate; a file, it can be read again)
  prefix    out/sample (from -o)
  reference tests/data/ref_100k.fa
segment 1  one read of tests/data/markdup.bam
  step 1  mark-duplicates    transform, first read  takes every record, then finds the duplicate sets
segment 2  a second read of tests/data/markdup.bam through step 1 (replay: nothing is spooled)
  step 1  mark-duplicates    stream marked          out/sample.marked.duplication_metrics
  step 2  alignment-summary  tap on marked          out/sample.alignment_summary_metrics
  step 3  insert-size        tap on marked          out/sample.insert_size_metrics
  step 4  quality-yield      tap on marked          out/sample.quality_yield_metrics
  step 5  duplicate-metrics  tap on marked          out/sample.duplicate_metrics
  step 6  write              BAM of marked          out/sample.marked.bam, out/sample.marked.bai
barrier    step 1 mark-duplicates: replays the source; about 16 records by the file size, no spill under MAX_RECORDS_IN_RAM 4,000,000; memory about 100.0 MiB
temp disk  none expected
output     about 1.1 KiB written by write
memory     about 100.0 MiB at the peak, a fixed 100.0 MiB for the chain beside its barriers included (a small input takes less)
nothing ran: this is a plan
```

The plan lists:

- The segments, and what each reads. A transform that reads its input twice
  (`mark-duplicates`) replays a file source: the file is read a second time and nothing
  is spooled. A piped source is copied to the temporary directory once (the spool) and
  the copy is read twice.
- Each step with its role, the stream it reads and every file it writes, the further
  formats included. A tool whose `OUTPUT` is a prefix is listed as `PREFIX.*`.
- One `barrier` line per transform (section 15): whether it holds its records, replays the
  source, replays a copy of the piped source, or spools its input; the bytes it spills or
  the estimated record count with the bound under which nothing spills; its memory. The
  record count comes from the file size: 69 compressed bytes per record for a BAM (102
  when its header says `queryname`), 21 for a CRAM, 434 for SAM text (the HG002 exome BAM
  and the NA12878 genome CRAM of docs/benchmarks.md). A file that compresses better gets a low count and a low spill estimate: a 22 MB BAM of 3,000,000
  identical-pattern records (7.3 bytes each) is estimated at 318,938 records. `mark-duplicates` above
  `MAX_RECORDS_IN_RAM` records spills 59 bytes per record, the larger of the two peaks
  measured there (46 and 59): 8.2 GiB for the exome BAM (measured 6.5 GiB) and 41.3 GiB
  for the genome CRAM (measured 42.1 GiB).
- The temporary disk (every spill and spool, as if all were held at once), the bytes
  `write` makes, and the peak memory. A pipe has no size to read, so those read
  `unknown`.
- A `space` line for a file system with less free space than the chain needs there. The
  run refuses then, before the first read (section 15).

The plan reads the header of a file source, so a file that is no BAM exits 1 there. It
reads no record and nothing of a pipe.

## 5. The Rust API

### 5.1 The stage layer

`bammetrics_core::stage` (`crates/core/src/stage.rs`, `stage/exec.rs`,
`stage/sinks.rs`).

| Item | Role |
|---|---|
| `Source` | Opens the stream. `read(threads, bytes, consume)` hands `consume` the `Frame` and the `Records`. `replayable()` is true for a regular file. `spool(temp)` copies a pipe. `file()` names the file a lowered chain reads |
| `FileSource` | A BAM, CRAM or SAM text by its first bytes, read through `pass::with_counted_stream` |
| `Records` | `next_batch(&mut Vec<RawRecord>)`: at most 1024 records per batch (`BATCH_RECORDS`) |
| `Frame` | The stream at one position: `kind`, `order`, `header`, `cram`, `label`, `reads` (the files a writer must not replace) |
| `Transform` | `apply(batch)` in place, `finish(out)` for held records, `header(up)`, `counters()`, `products()` |
| `Barrier` | `push(batch)`, then `seal()` gives `Sealed { header, next }` |
| `Next` | `Records(...)` when the barrier holds the records (a sort). `Replay(transform)` when the source is read again through a transform (`MarkDuplicates`) |
| `Sink` | `accept_batch(&[RawRecord])`, `done()`, `finish()`. One thread per sink |
| `CollectorSink` | The tap adapter: any `pass::Collector` as a `Sink` |
| `Emit` | Where a producer puts each record. `BamWriter` and `Vec<RawRecord>` implement it |
| `Step` | `tap`, `tool_tap(tool, inv)`, `tool_stage(tool, inv)`, `stage`, `sink`, `write(WriteOptions, threads)` |
| `ChainPlan` | `source`, `steps`, `threads`, `tmp_dir`, `lower` |
| `run_chain(plan, progress)` | One `Result<Products>` per step, in step order. `Err` for a plan refused before the first read or a source that cannot be read in the first segment |
| `Products` | `output` (metrics the caller writes), `files` (a writer's renamed files), `counters` (a stage's records in, out, changed), `notes` (lines for standard error) |
| `halt(error)`, `Halt` | `Stopped` or `NotRun` for the result of a step another step halted, `None` for a step's own failure |
| `WriteOptions` | `path`, `format` (`WriteFormat::Bam`, `Sam` or `Cram`, from the name or `--output-format`), `compression_level`, `create_index`, `create_md5`. `check()` refuses a name of a later wave and an index beside SAM text. `files()` lists the BAM, SAM or CRAM, the `.bai` and the `.md5` |
| `later_kind(path)` | The kind and wave of a name `write` cannot make yet: `("FASTQ", "W22")` |
| `Outputs` | The files of a chain's steps in a hash map, keyed by device and inode or by the resolved path. `claim(label, path)` is a usage error when a step already writes the file. `owner(path)` |
| `check_distinct_outputs(outputs)`, `writer_of(outputs, path)` | The same check over a slice |
| `stage_of(tool)`, `label_of(tool)`, `check_stage_options(tool, inv)` | The stage function of a registry tool, the default name of its stream (`marked`), and its options checked before a run. One tool: `MarkDuplicates` |

`crates/core/src/bam_writer.rs` holds `BamWriter` (the header byte for byte, BGZF on
worker threads, htsjdk's order check, an `AtomicFile` under a temporary name) and
`finish` (the `.bai` and `.md5`, then the renames). `MarkDuplicates` and `write` share
it. The `.bai` comes from `bai::index_bam`: the finished BAM is read once more on the
writer's BGZF workers, so its offsets are the ones of the file as written.

### 5.2 The verb as a library

`bammetrics_cli::chain` (`crates/cli/src/chain.rs` and `chain/`). A caller that is not the
command line builds a `Spec`, calls `prepare`, then `execute`. The Python module spells the
argv of its keywords and calls `parse`, so one parser builds the `Spec` of a command line
and of a Python call (section 5.3).

| Item | Role |
|---|---|
| `Spec`, `StepSpec` | A chain as asked for: `source`, `reference`, `prefix`, `threads`, `tmp_dir`, `formats`, `also`, `sample`, `quiet`, `verbose`, `sort_memory`, `force_space`, `stringency`, `plan_file`, `stage_executor`, `free_space`, `steps`. A step is its verb and its arguments as text (`StepSpec::new("insert-size", &["--minimum-pct=0.1"])`) |
| `parse(args)` | The command line into `Asked::Chain { spec, dry_run }` or `Asked::Help(text)` (`chain/args.rs`) |
| `read_plan(path)`, `parse_plan(text, file)` | A plan file into a `Spec` (`chain/planfile.rs`) |
| `prepare(&Spec) -> Result<Chain, VerbError>` | Every usage check. Nothing is read or written (`chain/prepare.rs`) |
| `Chain` | `source`, `source_path`, `reference`, `prefix`, `threads`, `steps`, `report`, `replayable()`, `segments()` |
| `Planned` | One step: `number`, `verb`, `tool`, `role` (`Role::Transform`, `Tap`, `Write`), `label`, `reads`, `leaves`, `segment`, `files` (`PlannedFile { path, kind }`), `late_prefix`, `stdout`, `command` (the `StringHeader` command line) |
| `plan_text(&Chain)` | The text of section 4 (`chain/describe.rs`) |
| `execute(Chain) -> Result<Outcome, VerbError>` | Runs the chain and writes every file. It prints nothing. `Err` for what stops the chain before its first record (`chain/run.rs`) |
| `Outcome` | `chain`, `steps`, `files` (the chain's own, `PREFIX.report.md`), `code` (0 or 1) |
| `StepOutcome` | One per step: `number`, `status` (`StepStatus::Ok`, `Failed`, `Stopped`, `NotRun`), `error`, `files` (`Written { path, kind }`), `counters` (a transform's `records_in`, `records_out`, `changed`), `notes` |
| `summary(&Outcome)` | The lines the command prints on standard error |
| `global_flags()`, `step_flags(verb)`, `write_flags()` | The option tables, as `grammar::Flag`s: the source of `--help`, of the completions, and of a caller's keyword list |
| `sample_name(&Chain)` | `--sample`, else the first `@RG SM` of a file source, else the source name: the Python `Run.sample` |
| `Chain::last_stream()` | The stream after the last transform, `in` without one: a tap on it writes the names `qc` writes |
| `run(args)`, `run_plan(args)` | The command and its `plan` form: the exit code |

`crates/cli/tests/chain.rs::the_library_entry_returns_one_outcome_per_step` pins this
entry.

### 5.3 Python

`bammetrics.chain(source, steps, output=, reference=, threads=, tmp_dir=, format=, also=,
sample=, resume=, plan=, dry_run=)` and `bammetrics.step(verb, **options)` run this verb in the
process. docs/python.md is their reference.

```python
import bammetrics as bm

run = bm.chain("in.bam", [
    bm.step("mark_duplicates"),
    "alignment_summary", "insert_size", "quality_yield", "duplicate_metrics",
    bm.step("write", index=True),
], reference="ref.fa", output="out/NA12878")
run.stages[0].records_in            # a transform's record counts
run.metrics["insert-size"]          # out/NA12878.insert_size_metrics, parsed
```

`crates/cli/src/chain/python.rs` holds the tables both sides share:

| Item | Role |
|---|---|
| `chain_keywords()`, `chain_keyword(name)`, `source_keyword()` | The keywords of `chain`: every flag of `global_flags()` in snake case (`output` for `--output-prefix`), without the source, `--quiet`, `--verbose` and `--help` |
| `step_table(name)` | A `StepTable` for a verb named in snake or kebab case, by its Picard class or by its function name: `verb`, `role`, `keywords`. An unknown name is a usage error with the closest step |
| `StepTable::keyword(name)` | One keyword, or the usage error with the closest name. `input` (except on `merge_sam_files`, `gather_bam_files` and `jumping_library`, where it names the further inputs) and `output` of `mark_duplicates` are refused with the reason |
| `Keyword` | `name`, `long` (the flag it spells), `ty`, `many` (`One`, `Repeated`: one occurrence per item, `Joined`: one comma list) |
| `step_names()` | The 55 names `step` takes, snake case: the 53 tool verbs in `bammetrics tools` order, then `write` and `qc` |
| `metric_key(chain, step, path)` | The key of a metrics file in `Run.metrics`: the verb, `STREAM.` in front for a tap that does not read the last stream, `.FILE` after it for a further file |
| `stubs()` | The `bammetrics.pyi` entries: `_StepVerb`, `step`, three `chain` overloads |

A step's keywords are those of the tool's Python function
(`bammetrics_core::registry::python_keywords`, the list the generated registry stub
prints), each matched to its flag in `step_flags(verb)`, plus `format`, `panel` on
`hs_metrics` (the verb flag `--panel`, which the error for missing intervals names) and,
on a transform, `as_`. The test `no_keyword_of_a_step_verb_is_lost` fails when a keyword
of a function has no flag.

`step` and `chain` check each value when the call is made (`value_text` in
`crates/py/src/chain.rs`): the type of the option, a float that is NaN or infinite, and a
list of a text option, which holds `str` items only. Each is a `UsageError` naming the
keyword. One number alone stands for the text of a text option, as in a plan file.

`crates/py/src/chain.rs` holds `Step`, `Stage`, `chain` and `step`. A call builds the argv
(`--long=VALUE` per keyword, `+ VERB --long=VALUE ...` per step), then runs `parse`,
`prepare` and `execute` with the GIL released, reads every metrics file back, and returns
the `Run`. `dry_run=True` returns `plan_text`. `tests/python/test_chain.py` (77 tests)
checks chained equals staged through the Python functions, both plan files against their
inline steps, the bytes of a Python call against the command's, `Run.stages`,
`Run.metrics`, one error per misuse, pickling and the GIL release.

## 6. Segments and threads

The executor cuts the steps at barriers. Each part is a segment.

- One segment is one read of its input on the calling thread. That thread also runs the
  transforms and feeds the barrier.
- Every tap and every sink has its own thread. It takes batches of at most 1024 records
  through a queue of 8 batches. A batch is shared by the sinks of one position
  (`Arc<Vec<RawRecord>>`). A position in front of a transform copies its batch once,
  and only when a sink listens there.
- The reader keeps one reference to each batch in flight (`Pool`). It reads the next
  records into the buffers of the oldest batch no sink holds any more. A segment
  allocates no record in steady state, and no sink thread frees one.
- BGZF blocks are inflated on `threads` workers. `write` compresses on `threads`
  workers.
- Segments run in turn. Peak memory is the largest segment, not the sum.

`in.bam + tap A + mark-duplicates + tap B + write` is two segments. Segment 1 reads
`in.bam`, feeds tap A and the first read of `MarkDuplicates`. The barrier is sealed.
Segment 2 reads `in.bam` again through the second read of `MarkDuplicates` and feeds
tap B and the writer.

A chain of taps alone is not run by this executor. It is lowered to
`pass::run_isolated`, the read `qc` uses, so it costs what `qc` costs and writes the bytes
`qc` writes (section 10). `qc` itself is unchanged. `--no-lower` (not in `--help`) forces
the stage executor for such a chain, to measure it.

## 7. MarkDuplicates as a barrier

The tool has two reads. In a chain they become the two halves of a barrier.

- First read: `build::Builder::push` takes each record of the stream and makes the read
  ends. The tool calls it from its own read of the inputs. The barrier calls it per
  batch.
- Seal: `mark::run` finds the duplicate sets. `write::output_header` adds the `@PG`
  line. The barrier returns `Next::Replay` with the `write::Copier`.
- Second read: the executor opens the source again. Each batch goes through
  `Copier::copy`, which sets or clears the duplicate flag, edits the tags and drops
  removed records. Each kept record goes to an `Emit`: the BAM writer for the tool,
  the batch of the next step for a chain.

The file source is read twice. Nothing is spooled for a regular file. A pipe is copied
to `--tmp-dir` once, as the tool does. Read ends spill to `--tmp-dir` as in the tool
(`MAX_RECORDS_IN_RAM`).

A step after the barrier is built from the sealed header, so it sees the new `@PG` line.

Behind another transform the input of `mark-duplicates` is a stream no file holds. Its
first read then also writes every record it takes to a spool in `--tmp-dir`, and the
second read replays the spool (section 15). W04 refused that chain; W05 runs it.

## 8. Options in the middle of a chain

`STOP_AFTER` and `ASSUME_SORTED` belong to the tap that carries them.

- `--stop-after N` on a tap: the tap takes the first N records its position emits. After
  `mark-duplicates --remove-duplicates` those are the first N records that were kept,
  the same N records the tool reads from the written BAM.
- It never shortens the stream. The steps after the tap get every record.
- The read stops early only in a segment with no stage, when every tap of it is done.
  A lowered chain stops at the largest `STOP_AFTER` when every tap has one.
- `--no-assume-sorted` on a tap checks the header of its position. After a barrier
  that is the sealed header.

`crates/cli/tests/chain.rs::stop_after_in_the_middle_of_a_chain_belongs_to_its_tap` pins
this against the staged runs. `mark-duplicates` has no `STOP_AFTER`.

## 9. Failure

| What fails | Result |
|---|---|
| A tap (cannot be built, fails on a record, fails to finish, panics) | That tap has its error. Every other step finishes |
| A transform on a record, or a transform or `write` that cannot be built | Its segment stops. The failing step has its error. Every other step of the segment is `stopped: <step> failed: <cause>`, taps in front of it included, because they saw part of the stream. Every later step is `not run: <step> failed: <cause>` |
| A barrier at its seal | The barrier has its error. The taps of its segment keep their results: they saw the whole stream. Later steps are `not run` |
| `write` on a record | The segment stops as for a transform when the reader sees the failure between two batches. When the stream ended first, the taps keep their results. Later steps are `not run` either way |
| The source in the first segment | `run_chain` returns the error, as `pass::run_isolated` does. The chain exits 1 |
| The source in a later segment | The steps of that segment are `stopped`, later steps `not run` |

A piped source in front of a transform is read from its copy in the temporary directory.
An error of that read names the source as given (`/dev/stdin: not a BAM or CRAM file`),
never the copy (`FileSource::read`, `Error::naming_input`). The `MarkDuplicates` tool
names its piped `INPUT` the same way.

Steps of an earlier segment keep their results. A written BAM is renamed into place only
when its segment finished. A failed or stopped writer leaves no file and no temporary
file. A tap that failed on its own before its segment stopped keeps its own error.

The summary says which case it was. A tap that fails reads `failed alone`. A transform or
`write` that halted other steps reads `failed and stopped its segment`. The last line
counts the steps by how they ended. The exit code is 1.

```
$ bammetrics chain tests/data/duplicate_marking.bam -o out/s \
      + quality-yield --no-assume-sorted + mark-duplicates + insert-size + write
bammetrics chain: step 1 (quality-yield): failed alone: tests/data/duplicate_marking.bam should be coordinate sorted but the header says the sort order is queryname; pass ASSUME_SORTED=true to read it anyway
bammetrics chain: step 2 (mark-duplicates): 13 records in, 13 out, 5 changed
bammetrics chain: 3 of 4 steps finished, 1 failed; 3 files written
$ bammetrics chain tests/data/insert_size.bam -o out/s + quality-yield + mark-duplicates + insert-size + write
bammetrics chain: step 1 (quality-yield): stopped: step 2 (mark-duplicates) failed: Cannot invoke ... is null (read r0)
bammetrics chain: step 2 (mark-duplicates): failed and stopped its segment: Cannot invoke ... is null (read r0)
bammetrics chain: step 3 (insert-size): not run: step 2 (mark-duplicates) failed: Cannot invoke ... is null (read r0)
bammetrics chain: step 4 (write): not run: step 2 (mark-duplicates) failed: Cannot invoke ... is null (read r0)
bammetrics chain: 0 of 4 steps finished, 1 failed, 1 stopped, 2 not run; 0 files written
```

Tests: `crates/core/src/stage/tests.rs` (synthetic stages, each row of the table) and
`crates/cli/tests/chain.rs::a_failing_tap_fails_alone_and_a_failing_transform_stops_its_segment`
(a tap that asks for a coordinate-sorted stream on a queryname-sorted BAM,
`mark-duplicates` on a BAM without a read group).

## 10. What equals what

### 10.1 Chained equals staged

Staged means `MarkDuplicates I=in.bam O=marked.bam M=dup.txt`, then each tool on
`marked.bam`. Chained means one `chain` command.

What is equal:

- Each tap's metrics file after its `StringHeader` lines, byte for byte.
- The `mark-duplicates` metrics file after its `StringHeader` lines.
- Every record of the written BAM, byte for byte, and every header line with the `CL`
  field of `@PG` set aside.
- With `PROGRAM_GROUP_COMMAND_LINE` given, the whole BAM, its `.bai` and its `.md5`,
  written to a file or to standard output.

`crates/cli/tests/chain.rs` checks this for 8 taps (`alignment-summary`, `insert-size`,
`quality-yield`, `duplicate-metrics`, `quality-score-distribution`,
`mean-quality-by-cycle`, `base-distribution-by-cycle`, `gc-bias` with its summary file) on
`tests/data/markdup.bam` (2 option sets), the query-name sorted
`tests/data/duplicate_marking.bam`, corpus t1 (2 option sets) and the foundation corpus
(2 option sets, one spilling with `MAX_RECORDS_IN_RAM=500`). The corpus cases are
skipped when `tests/corpus/out` is absent.

### 10.2 A chain without a transform equals qc

`chain in.bam -o P + quality-yield + alignment-summary + insert-size + duplicate-metrics`
writes the four metrics files `qc in.bam -o P` writes: the same names and the same bytes,
apart from the `Started on` line, which carries the clock. The command line in the
`StringHeader` is the tool's: `INPUT`, the file options, `REFERENCE_SEQUENCE` when the
tool reads one or the source is a CRAM, then the step's options.
`crates/cli/tests/chain.rs::a_chain_without_a_transform_equals_qc` checks a BAM, a BAM
with `-R` and a CRAM. The golden suite checks the records of `chain_taps` and
`qc_generic` against each other.

### 10.3 The plan file equals the command line

`crates/cli/tests/chain.rs::the_plan_file_equals_the_command_line` runs both examples of
section 3 both ways: the same file names, the same bytes apart from the clock line, the
same BAM and index, the same `--dry-run` text, the same `plan` text. The golden suite
asserts that the records of `chain_plan_file` and `chain_mark_duplicates`,
`chain_taps_plan_file` and `chain_taps_formats`, `chain_dry_run_plan_file` and
`chain_dry_run` hold the same exit code, text and files.

### 10.4 What differs on purpose (also in KNOWN_DIFFERENCES.md)

- The `StringHeader` command line. A tap writes `INPUT=<the chain's source>` and, behind
  a transform, ` CHAIN=<the Picard classes of the transforms in front, joined by +>`. A
  tap on the source writes the line the tool alone writes.
- The `@PG` `CL` of the BAM is the step's command line, which has no `OUTPUT`.
- `mark-duplicates` in a chain writes no BAM. `write` does.

## 11. Measurements

Mac Studio, Apple M1 Ultra, 20 cores, 64 GB. Release build, 8 threads, `/usr/bin/time -l`.
Other builds and test runs used the machine during every run.

### 11.1 A chain without a transform against `qc` (the W03 gate, 5 percent)

Input: chr20, 21 and 22 of the HG002 exome, 573 MB, 8.4 M records
(`~/data/bammetrics-scale/waves/w03/chr20_21_22.bam`). Four collectors on one read
(quality yield, alignment summary, insert size, duplicate metrics), the pass of
`qc --preset wgs`. Seven runs each, interleaved.

W03, through the hidden `__chain` entry:

| Run | Median wall | Range | Against `qc` | Median max RSS |
|---|---:|---|---:|---:|
| `qc --preset wgs` (also writes the verdict and the page) | 6.49 s | 6.48 to 6.72 | | 44 MiB |
| `__chain`, lowered to `pass::run_isolated` | 6.44 s | 6.40 to 6.50 | -0.8 % | 42 MiB |
| `__chain LOWER=false`, the stage executor | 6.65 s | 6.54 to 6.95 | +2.5 % | 19 MiB |

W04, through the `chain` verb:

| Run | Median wall | Range | Against `qc` | Median max RSS |
|---|---:|---|---:|---:|
| `qc --preset wgs` (also writes the verdict and the page) | 6.48 s | 6.44 to 6.70 | | 44 MiB |
| `chain`, lowered to `pass::run_isolated` | 6.44 s | 6.38 to 6.98 | -0.6 % | 42 MiB |
| `chain --no-lower`, the stage executor | 6.53 s | 6.40 to 6.85 | +0.8 % | 20 MiB |

Both forms are inside the gate in both waves. The kill rule of docs/07 does not apply. A
chain of taps is lowered anyway: it is the same code as `qc`. The four files of the lowered
chain and of `qc` are equal apart from the `Started on` line and the prefix in `OUTPUT=`.

### 11.2 One chain against the staged sequence

W03, on the HG002 exome (10,339,495,604 bytes, 150,538,907 records, coordinate sorted),
one run each, default options, `PROGRAM_GROUP_COMMAND_LINE` given on both sides so the
BAMs can be compared whole:

| Run | Wall | Max RSS |
|---|---:|---:|
| `MarkDuplicates` to a file (staged, step 1) | 308.3 s | 800 MiB |
| The four collectors on the marked BAM, one shared read (staged, step 2) | 106.0 s | 95 MiB |
| Staged sequence, both steps | 414.3 s | 800 MiB |
| `__chain`: `MarkDuplicates` + 4 collectors + `write` | 316.0 s | 807 MiB |

The chain takes 23.7 percent less wall time than the staged sequence and 7 MiB more
memory at its peak. It takes 2.5 percent longer than `MarkDuplicates` alone and also runs
the four collectors. Both BAMs are 10,210,246,600 bytes and `cmp` finds no difference.
The five metrics files are equal after their `StringHeader` lines. The chain reads the
10 GB input twice and the marked BAM never; the staged sequence reads the input twice and
the marked BAM once.

W04, the same comparison through the `chain` verb on the 573 MB slice, three runs each,
interleaved:

| Run | Median wall | Median max RSS |
|---|---:|---:|
| `MarkDuplicates` to a file (staged, step 1) | 15.69 s | 643 MiB |
| The four collectors on the marked BAM, a chain of taps (staged, step 2) | 6.49 s | 42 MiB |
| Staged sequence, both steps | 22.18 s | 643 MiB |
| `chain`: `mark-duplicates` + 4 collectors + `write` | 16.49 s | 645 MiB |

The chain takes 25.7 percent less wall time than the staged sequence and 5.1 percent more
than `MarkDuplicates` alone. Both BAMs are 584,413,820 bytes and `cmp` finds no
difference. The four metrics files are equal after their `StringHeader` lines. `chain
--dry-run` on the slice took under 0.01 s and 4 MiB. It estimated 8,306,326 records (the
file holds 8.4 M) and 467.4 MiB of temporary disk; the peak of that run was not measured.

### 11.3 Batch size and the pool

The first build of the executor used batches of 4096 records and allocated every record.
On the exome it took 345.2 s for the chain of 11.2 (the tool 300.8 s, the staged sequence
408.4 s, same BAM bytes). Two changes followed, each measured on the slice with 4
interleaved runs, medians:

| Build | `MarkDuplicates` + `write` | + 4 collectors | Stage executor, no stage |
|---|---:|---:|---:|
| 4096 records, every record allocated | 18.9 s | 20.4 s | 6.73 s |
| 4096 records, pool | 18.0 s | 19.3 s | 6.70 s |
| 1024 records, pool (as built) | 15.4 s | 17.0 s | 6.6 s |
| 512 records, pool | 15.9 s | 17.0 s | 6.76 s |

The rows come from three sessions with different machine load, so compare within a pair:
no pool against pool (first two rows), 4096 against 512 (17.0 s against 15.2 s for
`write`, 18.7 s against 17.0 s with the collectors), 1024 against 512 (last two rows).
The `MarkDuplicates` tool alone took 16.2 to 16.9 s in those sessions. With 512 records
the system time of the chain rose from 12 s to 16 s for the same wall time.

### 11.4 The Python call against the command

W04, `bammetrics.chain` against `bammetrics chain` on the 573 MB slice: `mark-duplicates`,
four collectors (quality yield, alignment summary, insert size, duplicate metrics) and
`write`, 8 threads, `PROGRAM_GROUP_COMMAND_LINE` given on both sides. Three runs each,
interleaved, each in its own process under `/usr/bin/time -l`. The Python time includes
the interpreter start, the import and reading the five metrics files back.

| Run | Median wall | Range | Median max RSS |
|---|---:|---|---:|
| `bammetrics chain` | 16.04 s | 15.95 to 16.24 | 646 MiB |
| `bammetrics.chain(...)` | 16.13 s | 15.94 to 16.20 | 658 MiB |

The Python call takes 0.6 percent longer and 12 MiB more. Both BAMs are 584,413,868 bytes
and `cmp` finds no difference. The five metrics files are equal after their `StringHeader`
lines. `Run.stages[0]` reports 8,445,233 records in, 8,445,233 out and 422 changed.

## 12. Departures from docs/07 section 2

W03:

| docs/07 | As built | Reason |
|---|---|---|
| `Source::open` returns `Box<dyn Records>` | `Source::read` takes a closure | The CRAM decoder (`cram::decode`) runs its workers in a thread scope and hands the stream to a closure. An owned reader cannot leave that scope |
| `Tool` gains `stage`, `source`, `sink`, `tier` fields | `stage::stage_of(tool)`, a table by name | A new field breaks every `Tool { ... }` literal, including those on the two branches built at the same time. One tool has a stage |
| `Transform::requires()` checked at plan time | A stage is built from the `Frame` of its position and fails there with the tool's own message | The planner is W05. `MarkDuplicates` accepts two orders, which one `Order` value cannot say |
| `Worker` and `Fanout` of `pass.rs` reused for the taps | `SinkWorker` in `stage/exec.rs`, the same shape | `pass::Worker` carries `bam::Record` batches. A sink takes `RawRecord` batches and needs the done, failed and stopped signals |
| Batches of at most 4096 records | 1024 records, and a pool that reuses the buffers of finished batches | Measured on the exome slice (section 11): with 4096 and no pool the chain `MarkDuplicates` + `write` was slower than the tool. A batch of 4096 exome records is 1.4 MB, one of 1024 is 360 KB (opinion: the transform reads the smaller batch back from a nearer cache). 512 gave the same wall time and more system time |
| One reader thread and one transform thread | The transforms run on the reader thread | The tool does the same. A second thread adds one queue hop per batch. Measured on the exome: the chain takes 316 s, the tool alone 308 s (section 11) |
| 4 deflate workers per writer | `threads` workers, as the tool uses | The written BAM is then the tool's BAM byte for byte at any thread count, and compression is the slowest part of the chain |
| `Barrier::seal` on a `Send` trait | `Transform` and `Barrier` are not `Send` | They run on the calling thread and borrow the chain's `TempDir` |
| `run_chain(&ChainPlan, &dyn Progress) -> Result<Vec<Products>>` | `run_chain(ChainPlan, Option<OnProgress>) -> Result<Vec<Result<Products>>>` | A step can fail alone, so each step has its own result. `OnProgress` is the callback `qc` already has |

W04:

| docs/07 | As built | Reason |
|---|---|---|
| The verb is `chain`, with `run` as an alias or not (Q2) | `chain` only | The owner's decision on Q2 |
| `write` picks BAM, SAM, CRAM or FASTQ by extension | BAM only. A `.sam`, `.cram`, `.fastq` or `.fq` name is a usage error naming W05, W15 or W22 | The SAM text writer is W05 (section 14). The CRAM writer shipped in W15, 0.33.0: `write r.cram` against the chain's `-R` (section 23). FASTQ output is W22 |
| Two taps of one tool at one position is a usage error | It is one when both are named by the prefix, since they write one file. With `-o` on one of them the chain runs | A step's own `-o` is allowed, and two option sets of one tool on one stream is a use (`quality-yield` with and without original qualities) |
| The `# INPUT=` line names the chain source and the step list | `INPUT=<source>` plus ` CHAIN=<transforms in front>`, as W03 wrote it. A tap on the source writes the tool's own line | A chain without a transform must write the bytes `qc` writes (section 10.2). The steps behind a tap do not change its records, so the line names only the transforms in front |
| A transform's kept product is `PREFIX.<label>.bam` | The same, written by `write` without a path. The transform's own metrics file is `PREFIX.<label>.duplication_metrics` | docs/07 names no file for a transform's metrics. The label says which step wrote it |
| `plan` prints segments, replay or spool per barrier, temp disk, and what `--resume` skips | Segments, replay or spool, the files of each step, the temp disk estimate. No `--resume` line | Resume per step came with W14 (manifest 2, section 27) |
| The planner refuses when free space is below the estimate | The plan prints the estimate and refuses nothing | The planner and its `statvfs` call are W05 |
| `--also` on `chain` as on `qc` | `json`, `jsonl`, `tsv`, `csv`, `parquet`, `xlsx`, `docx`, `sql`, `markdown`. `junit` and `multiqc` are a usage error | They print a verdict. `qc` as the last step of a chain came with W14 (section 29) |
| `--format a,b` | `picard` must be first, on the chain and on a step | A chain names Picard files under a prefix, as the multi-file verbs do. A first format that is not `picard` would put JSON under a Picard name |
| `parse_tables` moves to core | It moved, and it reads one-line arrays of scalars as well | A repeatable option of a step (`level`, `comment`) needs a list. A preset still takes no list |
| `chain --plan plan.toml` reads the same steps | It reads the source and the chain's options too, and the command line may replace them | One file then describes one run, and one file can run on many sources (`--plan FILE other.bam -o other`) |
| `plan` gains the chain form | `plan` with a `+` step, with `--plan FILE`, or with the word `chain` first | `plan` is a clap verb with the arguments of `qc`. A `qc` command line holds no `+` and no `--plan`, so the two cannot be mistaken |
| Text-valued options are checked before the first read | They are, since the 0.22.0 review: `tap_options::check` runs each one-pass tool's option parser without a header (section 2.6). A value judged against the header or a file still fails its tap when it is built | A collector parses its options and then builds from the header of its position. The parser alone needs no header and opens no file, so no reference or interval list is read twice |

W04, Python (docs/07 section 2.2, the `crates/py` paragraph):

| docs/07 | As built | Reason |
|---|---|---|
| `crates/py/src/verbs.rs` holds `bm.chain` | `crates/py/src/chain.rs` | `verbs.rs` has 1,138 lines. The chain classes and functions are one file of 632 |
| `bm.chain(..., resume=)` | No `resume` keyword | Resume per step came with W14 (`resume=`, section 27). A keyword that does nothing would be a lie in the signature |
| `Run.stages` holds the transform counters | One `Stage` per step: `status`, `error`, `files`, `reads`, `leaves`, `segment`, and `records_in`, `records_out`, `changed` for a transform (`None` for a tap and `write`) | A failed run needs the status of every step, and `Products::counters` exists for a stage only. A chain of taps lowers to the read of `qc`, which counts no records |
| `bm.step` validates against the generated registry signature | It does, through `registry::python_keywords`. `input` is refused except on `merge-sam-files`, `gather-bam-files` and `jumping-library` (the files after the stream), `output` is the step's own `-o` (a prefix for a tool that writes several files), `mark_duplicates` takes no `output`, and `format` and `as_` are added | A step reads the stream, plus the further inputs of those three, a transform writes no BAM in a chain, and `--format` and `--as` are step options with no row in the registry |
| Steps are named `"sort_sam"`, `"insert_size"` | Snake case, and also kebab case, the Picard class and the function name | `Run.metrics` and `collect(metrics=[...])` already spell verbs with a dash, and `tools()` gives the class |
| The keywords of `bm.chain` are `output`, `reference`, `resume` | `output`, `reference`, `threads`, `tmp_dir`, `format`, `also`, `sample`, `plan`, `dry_run` | Every option of the verb that is not about printing. The table is read from `global_flags()`, so a new option becomes a keyword |
| No word on `dry_run` | `dry_run=True` returns the plan text as a str | The command prints that text. The stub has three overloads, so `chain(...)` types as `Run[None]` and `chain(..., dry_run=True)` as `str` |
| No word on a failed step | `InputError` after the other steps wrote their files, with `.run` | `qc` and `collect` do the same for a failed collector |
| The Python module builds a `Spec` (section 5.2 as W04 part A wrote it) | It builds the argv and calls `parse` | `qc`, `check` and `report` already work that way. `--plan` with its replacement rule, the `+` refusal beside a plan and every option check then have one implementation |

W04, `fifo-buffer` (docs/07 section 2.2, the `chain.rs` paragraph):

| docs/07 | As built | Reason |
|---|---|---|
| `fifo-buffer` copies stdin to stdout and prints that chaining removed the pipe it buffered | It does, as registry tool 33 with a Picard face (`FifoBuffer`) and the verb. The line is on standard error and `-q` drops it. The Picard face never prints it | The Picard face prints what Picard prints. Standard output carries the copied bytes |
| No word on `fifo-buffer` inside a chain | Refused as a step: a usage error with its own reason (section 1) | It moves bytes between two commands and writes no record. A chain has no pipe between its steps. The generic refusal named W14, which will not make it a step |
| No word on the help of a step that cannot be one | `chain BAM + VERB --help` for such a verb prints the refusal and exits 2 | It printed the options of `wgs` and of `fifo-buffer` under the words `a tap` |

Added beyond the list of docs/07: `Emit` takes `&mut RawRecord` so a batch can take the
record without a copy. `Transform::products` carries the `MarkDuplicates` metrics.
`Frame` has `cram` (htsjdk rebuilds `MD` and `NM` for CRAM records), `label` and
`reads`. `Sink::done` carries `STOP_AFTER`. `stage::Outputs` makes the one-file-twice
check linear in the number of steps.

## 13. What comes next

| Wave | Adds |
|---|---|
| W05, done | The sort barrier (`Next::Records`), the planner (`Frame` per position, `requires`, the free-space refusal), the spool for a replay behind another transform, `--sort-memory`, SAM output of `write`; `write --index` writes Picard's `.bai` of the written BAM (`bai.rs`); `BuildBamIndex` and `make memgate` |
| W09 on | Each BAM-writing tool registers its stage in `stage_of` and its stream name in `label_of` |
| W14, done | `qc` as the last step of a chain, manifest 2, resume per step (`resume=` in Python), kept products (sections 27 to 30). The spool for tools that read their input their own way landed as the materialize fallback (section 31). Left: a progress line for a resumed run |
| W15, done | `write` to a CRAM against the chain's `-R` (`cram_writer.rs`) |
| W22 | `write` to FASTQ |
| W46 | Variant streams (`Kind`), pileup and capture sinks |

## 14. W05 part A: SAM text

W05 part A adds SAM text on both ends of a chain and the record orders the sort barrier
will use. The sort barrier, the planner and `write --index` through the BAI writer are
the rest of W05.

- The source may be SAM text, plain or in gzip members, a file or a pipe. It is told by
  its first bytes (`cram::Kind`) and parsed by `sam_text::Reader` into the records a BAM
  holds, at `VALIDATION_STRINGENCY=STRICT`. A SAM text file can be read again, so
  `mark-duplicates` replays it as it replays a BAM.
- `write PATH` writes SAM text when `PATH` ends in `.sam` (`WriteFormat::of`), and
  `--output-format bam|sam` (`output_format=` in Python, `output-format = "sam"` in a
  plan file) says which whatever the name is. `write - --output-format sam` prints SAM
  text on standard output. The plan lists the file as kind `sam`.
- The file is the one Picard `SamFormatConverter` writes for the same records, byte for
  byte (`tests/concordance/test_view_sam.py`): `sam_text::header_text(header, false)`,
  which is `@HD VN:1.6` first, then one `sam_text::write_record` line per record, tags in
  htsjdk's order. `--md5` writes `PATH.md5`. `--index` with SAM text is a usage error.
  `--compression-level` changes nothing for SAM text.
- `crates/core/src/sam_order.rs` holds the orders a sort barrier needs (`Order`,
  `CoordinateKey`, `QuerynameKey`, `DuplicateKey`), each checked against Picard `SortSam`
  (`tests/concordance/test_sam_order.py`).

| docs/07 | As built | Reason |
|---|---|---|
| `write` picks BAM, SAM, CRAM or FASTQ by extension | BAM and SAM by extension, and `--output-format` | Standard output has no extension. CRAM shipped in W15, 0.33.0 (section 23, the `r.cram` row). FASTQ stays W22 |
| SAM text reader and writer in one infrastructure row with header-merge and the BAI reader and writer | The reader, the writer and the comparators only | Part A. Header merge exists in `mark_duplicates/inputs.rs` since 0.18.0; BAI and the sort barrier are parts B and C |
| A chain carries `VALIDATION_STRINGENCY` | The source is read at `STRICT` | A chain has one source and no option of the source's own; the tools and verbs carry the option |
| The SAM file of `write` has the header of the stream, as the BAM has | The SAM file has htsjdk's `SAMFileWriter` header (`VN:1.6`) | Picard is the oracle for SAM text. The BAM of `write` keeps the stream's `VN`, which htsjdk also rewrites; changing that would change the bytes 0.22.0 writes |

## 15. W05 part B: the sort barrier, replay or spool, the planner

```bash
# Sort by name, mark duplicates, sort back, write the BAM and its index: one process
bammetrics chain in.bam -o out/s --sort-memory 4G \
    + sort-sam --sort-order queryname \
    + mark-duplicates \
    + sort-sam --sort-order coordinate \
    + write out/s.bam --index
```

### 15.1 sort-sam as a step

- `sort-sam` is a barrier that holds its records (`Next::Records`). Its first read pushes
  each record into a `sam_sort::RecordSorter`. The seal sorts what is held, merges the
  spilled runs, and hands the records on in order under the sorted header (`SO` set,
  `@HD VN:1.6`, as the tool writes it). The next segment reads them from the sorter.
- The buffer holds at most `--sort-memory` bytes (default 2 GiB) less what is resident
  beside it (`planner::sort_budget`): 64 MiB for the run merge, a spool writer and the
  batches in flight, and the read ends of the largest `mark-duplicates` step after it
  (`mark_duplicates::held_bytes`, 640 MB at `MAX_RECORDS_IN_RAM` 4,000,000), at least a
  quarter of `--sort-memory`. The allocator keeps the pages of the sort's records resident
  after the sort seals, so that step's lists come on top of them. A record counts its
  bytes, the item and two index slots of the sort (`sam_sort::ITEM_BYTES`). The sort of
  a run permutes the records in place, so no second vector of the run exists. The step's
  own `--max-records-in-ram` adds a record bound. Runs spill to `--tmp-dir` as BGZF level
  3, are merged at most 64 at a time with ties to the earlier run, and are removed once
  read, at a failure and at a signal (`spill::TempDir`).
- A run in memory is sorted by the port of Java's `TimSort` (`crates/core/src/timsort.rs`),
  as htsjdk sorts it. For a queryname order that is not total (a record flagged both
  first and second of its pair) the order is Java's or Java's
  `Comparison method violates its general contract!`, never a panic.
- A CRAM source gets the `MD` and `NM` tags htsjdk rebuilds when the step has a
  reference (`-R` reaches it).
- Standard error gets `step N (sort-sam): R records in, R out, 0 changed` and, when it
  spilled, `spilled B bytes in sorted runs, L merge levels`.

### 15.2 Replay or spool

docs/07 section 2.3: a barrier that reads its input twice (`mark-duplicates`) replays the
source when the source is a file and nothing stands in front of the barrier in its
segment. Otherwise it spools.

| Where the replaying barrier is | Its second read |
|---|---|
| First stage, file source | The source again (replay: nothing is spooled) |
| First stage, piped source | A copy of the pipe made before the first read |
| Behind another transform | The spool of its own input: every record its first read took, written to `--tmp-dir` as BGZF level 3 (`stage/spool.rs`), read back and removed after the replay |

So any number of transforms may follow one another. A second `mark-duplicates` runs on
the records the first one marked (its stream is `marked-2`), and standard error says
`spooled B bytes of its input to replay them`.

### 15.3 The planner and the free-space refusal

`crates/cli/src/chain/planner.rs` estimates, before the first read, per barrier: how it
hands its stream on, the bytes it spills or spools, and its memory. The record count
comes from the source's size (section 4): 69 bytes per record for a BAM, 102 for a BAM
whose header says `queryname` (since 0.23.0), 21 for a CRAM, 434 for SAM text. The
constants were measured on 8,445,233 records of the HG002 exome (docs/benchmarks.md, Sort
barrier and planner):

| Constant | Value | Measured |
|---|---:|---|
| Bytes of one record in memory | 350 | 2,951,869,900 inflated bytes over 8,445,233 records |
| Resident bytes per byte a sort holds | 1.25 | 4,390 MB resident for 3,526 MB held |
| Resident memory beside the barriers | 100 MiB | 425 MB at `--sort-memory 256M` less 1.25 times 268 MB |
| A queryname run or spool, per record | 105 bytes | 871,834,880, 877,285,505 and 879,914,157 bytes |
| A coordinate run, per record | 76 bytes | 616,683,279 and 639,427,369 bytes |
| `write` BAM per record, coordinate or other | 70 or 102 bytes | 584,011,111 and 856,813,333 bytes |
| `mark-duplicates` memory per record held | 200 bytes | 646 MiB (8.4 M records), 800 MiB (150.5 M) |

A barrier's memory adds the memory of the barrier before it: a sort that held its
records still holds them while the next barrier reads them, and the allocator keeps what
the barrier before freed. The two-barrier chain above peaked at 1.65 GiB with
`--sort-memory 1G`, more than either sort alone. The 1.25 resident bytes per byte held
were measured before the run sort permuted in place; since 0.23.0 a sort holding 1.44
GB peaked at 1,636 MiB (1.12 with the base), so the estimate errs high.

Before the first read `execute` adds the temporary bytes (every spill and spool, as if
all were held at once) and the bytes of each `write` file, by file system (the device of
`--tmp-dir` and of each output's directory), and asks each for its free space (one
`statvfs` call). When a file system has less free space than the chain needs there, the
chain exits 1 and makes nothing:

```
$ bammetrics chain in.bam -o out/s + sort-sam -SO queryname + write
bammetrics chain: not enough free space in out: the chain needs about 54.9 GiB there (out/s.sorted.bam) and 40.2 GiB are free
  --force-space runs the chain anyway
```

`--force-space` runs it. The plan prints a `space` line for such a file system and does
not refuse. A pipe has no size, so its estimate is unknown and nothing is refused.

### 15.4 Departures from docs/07

| docs/07 | As built | Reason |
|---|---|---|
| Spilled runs as BGZF level 1 (section 2.1) | Level 3 | zlib-rs level 1 spilled a queryname sort of 8.4 M records as 1,142,443,761 bytes in 16.7 s, level 3 as 871,834,880 bytes in 17.5 s; a coordinate sort 820,851,148 against 616,683,279 bytes in 8.1 against 8.9 s |
| The planner in `crates/core/src/plan.rs` | `crates/cli/src/chain/planner.rs` | The planned steps (`Planned`) and their Picard arguments live in the cli crate; the core keeps the executor and the spool |
| Replay when every earlier transform is deterministic | Replay only when the barrier reads the source itself; spool otherwise | Re-running earlier transforms needs each to be replayable; a spool of the barrier's own input is one rule for every stage. Measured cost: 879,914,157 bytes and the time to write them for 8.4 M records |
| `Transform::requires()` checked at plan time | Not built | `MarkDuplicates` accepts two orders and a sort accepts any; the stage still checks its input when it is built |
| `--sort-memory` sets the buffer | It does, and the step's `MAX_RECORDS_IN_RAM` bounds it too when given | The tool's option; Picard's `SortSam` counts records only, so the tool keeps that bound alone |
| Free-space refusal in `plan` | `plan` prints the shortage; the run refuses | A plan runs nothing; it says what the run would do |
| No word on a free-space override in tests | A hidden `--free-space BYTES` (not in `--help`) replaces the file system's answer | The golden suite records the refusal on any machine |
| One `libc::statvfs` call, dependency decision in W05 | `libc` 0.2 in the cli crate, already in the lock file through `signal-hook` and `clap` | No new crate enters the tree |
| A sort spill under 1.2 times the input (section 6) | A queryname run spills 1.52 times a coordinate BAM, a coordinate run 1.08 times | Queryname order puts records of distant positions side by side; the queryname BAM at level 5 is already 1.49 times the coordinate one, and no BGZF level reaches 1.2 (docs/benchmarks.md) |
| The planner's estimate within 25 percent of the measured spill (section 6, Q5) | Within 2.2 percent for the sort spill on both scale inputs, after one fix | At 69 bytes per record the 2,870,148,655-byte queryname slice counted 41.6 M records for 29.3 M and the sort spill estimate was 51 percent over. A queryname BAM now counts at 102 bytes per record, the bytes `write` makes in that order (docs/benchmarks.md, Scale gate of wave W05) |
| The scale gate on the exome (section 6) | On chromosomes 1 to 5 of the exome, 3,183,979,932 bytes | The planner gives the whole exome 29.4 GiB of temporary files and 9.8 GiB of output, over the 20 GB a run may hold on this machine |
| The header merge lifted for `SortSam` | Lifted to `crates/core/src/header_merge.rs` and used by `MarkDuplicates`; `SortSam` takes one `INPUT` and merges nothing | Picard 3.5.0 `SortSam` has one `INPUT`; `MergeSamFiles` (W13) is the next user |

## 16. W09: record-editing transforms

`add-comments-to-bam` and `replace-sam-header` are streaming transforms. They are built on
`crates/core/src/record_edit.rs`, the shape every record-editing tool of waves W09 and
W10 uses. A streaming transform edits each batch of records on the reading thread and
hands it on in the same read. It starts no segment and holds no record.

```bash
bammetrics chain in.bam -o out/s \
    + add-comments-to-bam -C "library prep v2" \
    + replace-sam-header --header new_header.sam \
    + insert-size \
    + mark-duplicates \
    + write --index
```

The plan lists both edits in segment 1, the segment they read. Their streams are named
`commented` and `reheadered`, so the tap on `insert-size` writes
`out/s.reheadered.insert_size_metrics`. `mark-duplicates` reads the edited stream, which no
file holds, so it spools its input (section 15.2). The run prints one line per edit:
`step 1 (add-comments-to-bam): 36 records in, 36 out, 0 changed`.

| Tool | Default stream | What the step does |
|---|---|---|
| `add-comments-to-bam` | `commented` | Adds one `@CO` line per `--comment` to the header of its stream. The source may be SAM text or a CRAM; the tool on its own copies a BAM only, as Picard does |
| `replace-sam-header` | `reheadered` | Puts the header of `--header` in place of its stream's. Each record keeps its reference index, as Picard's block copy of a `.bam` keeps it. The new `SO` must be `unsorted` or the stream's |
| `add-o-a-tag` | `oa_tagged` | Appends each record's alignment to its `OA` tag, or with `-L` only on the records over an interval. Every record leaves with the `bin` htsjdk's writer computes, as the tool on its own writes it |
| `clean-sam` | `cleaned` | Soft-clips each alignment past the end of its reference and the `MC` of a mate past the end, sets MAPQ 0 on unmapped records, and gives each record the `bin` htsjdk computes. The chain reads a SAM text source at its `--validation-stringency`, so a source with lines past the end of their reference needs `LENIENT`; the tool alone reads `STRICT` as `LENIENT`. Chained equals staged in `crates/cli/tests/clean_sam.rs` |
| `set-nm-md-and-uq-tags` | `tagged` | Sets `MD`, `NM` and `UQ` on each mapped record of a coordinate-sorted stream against the chain's `-R` (section 16.1) |
| `set-nm-and-uq-tags` | `tagged` | Picard's deprecated name of `set-nm-md-and-uq-tags`: the same step |

The rules of a record-editing step:

- The options of the BAM a tool writes on its own (`-o`, `--compression-level`,
  `--create-index`, `--create-md5-file`) belong to the `write` step. Each exits 2 before
  the first read: `CREATE_INDEX belongs to the write step in a chain`.
- A header without `@HD` leaves the step with `@HD VN:1.6`, as htsjdk holds it.
- Records of a CRAM source get `MD` and `NM` when the chain has `-R`, as the `sort-sam`
  step gives them.
- The step's counters are records in, records out and records changed. A tool that drops
  or adds records says so in them.

Chained equals staged (`crates/cli/tests/record_edit.rs`): each tool first, in the middle
and last, beside `mark-duplicates` and `sort-sam`, the two after each other, on the
fixtures, corpus t1 and the foundation corpus. The records of the written BAM are compared
as bytes. One case differs from the staged tool on purpose: `replace-sam-header` at the
first position of a SAM text or CRAM source keeps each reference index, where Picard's
tool on such a file keeps each reference name. With one dictionary in both headers the two
agree (`a_sam_source_reheaders_as_the_staged_tool`).

### 16.1 set-nm-md-and-uq-tags

`set-nm-md-and-uq-tags` (and `set-nm-and-uq-tags`, Picard's deprecated name for the same
tool) is a streaming transform that reads the chain's `-R`. It sets `MD`, `NM` and `UQ`
on each mapped record of its stream, which must be coordinate sorted, so the usual place
is after `sort-sam -SO coordinate`:

```bash
bammetrics chain reads.queryname.bam -R ref.fa -o out/s \
    + sort-sam -SO coordinate \
    + set-nm-md-and-uq-tags \
    + insert-size \
    + write out/tagged.bam --index
```

| Tool | Default stream | What the step does |
|---|---|---|
| `set-nm-md-and-uq-tags`, `set-nm-and-uq-tags` | `tagged` | `MD`, `NM` and `UQ` of each mapped record against `-R`, as the tool computes them; `--is-bisulfite-sequence`, `--set-only-uq` |

- Without `-R` the step exits 2 before the first read: `missing --reference-sequence`.
- A stream whose header is not `SO:coordinate` fails the step as the tool fails: `Input
  must be coordinate-sorted for this program to run. Found: queryname`.
- The tool on its own writes every record through htsjdk's writer, which gives each
  record the `bin` it computes. The step gives each record that `bin` too
  (`RecordEdit::rewrites`), so the written BAM equals the staged tool's byte for byte.
- `changed` counts the records whose tags the step changed: a second step after the first
  changes none.

Chained equals staged (`crates/cli/tests/set_nm_md_uq.rs`): the step first, after
`sort-sam`, last, alone, twice and under the old name, beside `mark-duplicates` and
`sort-sam`, with each option, on the fixtures, SAM text, a CRAM, corpus t1 and the
foundation corpus.

## 17. W10: steps that sort

`add-or-replace-read-groups` and `revert-original-base-qualities-and-add-mate-cigar` are
record-editing steps whose tool on its own writes through Picard's `presorted=false`
writer, which sorts at close. In a chain that sort makes the step a barrier that holds its
records, as `sort-sam` holds them (section 15.1): the plan lists it with a `barrier` line,
its spill and its memory, and the next segment reads the records it hands on.
`add-or-replace-read-groups` sorts only when its `--sort-order` differs from the order the
stream's header names, as the tool decides from its input's header. The plan is made
before the header is read, so it lists the step as a barrier either way; a step whose
stream is already in that order streams the records with the order check of Picard's
presorted writer and starts no segment.

```bash
bammetrics chain in.bam -o out/s \
    + add-or-replace-read-groups --rglb lib1 --rgpl ILLUMINA --rgpu unit1 --rgsm NA12878 \
    + revert-original-base-qualities-and-add-mate-cigar \
    + mark-duplicates \
    + write out/s.bam --index
```

| Tool | Default stream | What the step does |
|---|---|---|
| `add-or-replace-read-groups` | `rg` | Replaces every `@RG` line by one read group and sets `RG` on each record, streaming. With `--sort-order queryname`, `coordinate` or `duplicate` that differs from the stream's order it holds its records and hands them on sorted; the stream's own order goes through as it comes, its order checked; `--sort-order unknown` exits 2 |
| `revert-original-base-qualities-and-add-mate-cigar` | `reverted` | Holds its records: `QUAL` from `OQ` and the clip of `clean-sam` as they come, a queryname sort, the mate fields, `MQ` and `MC` of each pair, then `--sort-order` (default the stream's). It never skips its stream, so `--max-records-to-examine` exits 2 |

- The spilled runs go to the chain's `--tmp-dir`. Each sort holds at most
  `--max-records-in-ram` records when the step gives it, and its share of
  `--sort-memory` (the reverting step splits it between its two sorts).
- Every record leaves with the `bin` htsjdk's writer computes, set after the mate pass, as
  the tool on its own writes it.

Chained equals staged (`crates/cli/tests/read_groups.rs`): each step first, in the middle
and last, streaming and sorting, beside `mark-duplicates` and `sort-sam`, on
`tests/data/markdup.bam`, `tests/data/revert.bam`, corpus t1 and the foundation corpus,
and `--sort-order coordinate` on t1 and the foundation corpus, which samtools sorted with
its own order of tied records.
`test_add_or_replace_read_groups.py` and
`test_revert_original_base_qualities_and_add_mate_cigar.py` also compare each step's
written BAM with Picard's file of the tool.

## 18. W10: split steps

A split writes the records of its position over several BAMs in a directory: one per
library (`split-sam-by-library`) or files of about one size (`split-sam-by-number-of-reads`).
It is a sink of the executor, as `write` is: it changes nothing, so the step after it sees
the same records, and its failure stops its segment.

```bash
bammetrics chain sample.bam -o out/s \
    + mark-duplicates \
    + split-sam-by-library -o out/libraries --create-index
```

- The step's own `-o DIR` names the directory; the run makes it when missing. Without `-o`
  `split-sam-by-library` writes into the current directory, as the tool does, and
  `split-sam-by-number-of-reads` exits 2. A chain with no `-o` of its own and no `write
  PATH` names its metrics files under the step's directory: `chain s.bam insert-size
  split-sam-by-library -o s1` writes `s1/s.insert_size_metrics`.
- The files are BAMs named `.bam`, whatever the source: `LIBRARY.bam` and `unknown.bam`, or
  `OUT_PREFIX_0001.bam` on. `--create-index`, `--create-md5-file` and `--compression-level`
  belong to the step.
- A count split needs `--total-reads`: the stream is read once, and Picard refuses a stream
  without it (exit 2 before the first read).
- At most 64 files are open at once, whatever the number of libraries. Every file is
  renamed into place after the last record, with its `.bai` and `.md5`. A split that fails
  leaves none of its files; a chain that does not finish removes the directories it made
  that are still empty.
- The names a split makes are known only once the stream is read. The plan claims
  `DIR/<LIBRARY>.bam` or `DIR/OUT_PREFIX_<N>.bam` for the step, so two splits with one
  naming rule in one directory are a usage error before the first read. Each file the
  split makes is claimed when it is made against every planned file of the chain and the
  files of the other splits, and against every file the chain reads: a `write` into the
  directory under a library's name fails the split.

Chained equals staged (`crates/cli/tests/split_sam.rs`, `test_split_sam_by_library.py`,
`test_split_sam_by_number_of_reads.py`): `mark-duplicates` then a library split on the
fixtures and the foundation corpus, `sort-sam -SO queryname` then a count split, a split
beside a tap and `write`; every file has the staged run's name, header lines and record
bytes.

## 19. W13: ReorderSam, and MergeSamFiles outside the chain

`reorder-sam` is a record-editing step whose tool on its own writes through Picard's
`presorted=false` writer in the input's order (section 17). In a chain it moves each
record to the dictionary of `--sequence-dictionary` and holds the records when the
stream's header names a sort (`coordinate`, `queryname`, `duplicate`), then hands them on
sorted under the new dictionary; an `unsorted` stream goes through as it comes. The plan
lists it as a barrier either way: a sort of one pass in the stream's order ("takes every
record and sorts them in the stream order"), with the memory and spill of a `sort-sam`
step under `--sort-memory`. The dictionary is read when the step is built, before
the first record.

```bash
bammetrics chain in.bam -o out/s \
    + reorder-sam --sequence-dictionary hg38.dict -S \
    + mark-duplicates \
    + write out/s.bam --index
```

| Tool | Default stream | What the step does |
|---|---|---|
| `reorder-sam` | `reordered` | New `@SQ` lines, each record's reference and mate reference moved by name; with `-S` a record on a contig the dictionary lacks becomes unplaced and its mate loses `MC`; then the sort of the stream's order |

`merge-sam-files` is a materialized step since wave W14 (section 31): the stream is its
first `INPUT`, and the files after the verb are the further ones.

Chained equals staged (`crates/cli/tests/reorder_sam.rs`): the step first, in the middle
and last, beside `mark-duplicates` and `sort-sam`, on `tests/data/markdup.bam`, corpus t1
and the foundation corpus, with a contig added and with contigs dropped.
`test_reorder_sam.py` also compares the written BAM with Picard's file of the tool.

## 20. W13: mate fixing and filtering steps

`fix-mate-information` holds its records as `revert-original-base-qualities-and-add-mate-cigar`
does (section 17): a queryname sort, the mate fields of each name, then `--sort-order`
(default the stream's). `filter-sam-reads` streams the records its `--filter` keeps, or
holds and sorts them when `--sort-order` names another order than the stream's.

```bash
bammetrics chain in.bam -o out/s \
    + fix-mate-information \
    + filter-sam-reads --filter includePairedIntervals --interval-list panel.interval_list \
    + mark-duplicates \
    + write out/s.bam --index
```

| Tool | Default stream | What the step does |
|---|---|---|
| `fix-mate-information` | `fixed` | Holds its records: a queryname sort, the mate reference, start, strand, `MQ`, `MC` and `TLEN` of each pair, then `--sort-order` (default the stream's). It sorts by name even when the stream says `SO:queryname`, where the tool on its own takes the records as they come; `--assume-sorted`, `--sort-order unknown` and a second input exit 2 |
| `filter-sam-reads` | `filtered` | Keeps the records `--filter` keeps, streaming; with a `--sort-order` that differs from the stream's order it holds them and hands them on sorted. The pair filters (`includeAligned`, `excludeAligned`) need a stream in queryname order, `includePairedIntervals` one in coordinate order. `--write-reads-files` exits 2 |

Chained equals staged (`crates/cli/tests/mate_filter.rs`): each step first, in the middle
and last, beside `mark-duplicates` and `sort-sam`, streaming and sorting, on
`tests/data/mate_filter.bam`, its queryname form, `tests/data/revert.bam`, corpus t1 and
the foundation corpus.

## 21. W13: downsample-sam

`downsample-sam` is a streaming transform (stream name `downsampled`): it keeps the
records of a random part of the read names at its position and hands on the stream under
the header with its `@PG` line, as `DownsampleSam` writes it on its own.

| Step | Stream | What it does |
|---|---|---|
| `downsample-sam` | `downsampled` | ConstantMemory decides each record as it passes. HighAccuracy and Chained hold one chunk of records (`ceil(1 / ACCURACY)` new names) and hand it on once decided, so the steps after it see the same records in the same order. The seed is resolved against the header of the stream at its position |

- `--metrics-file` on the step writes the QualityYieldMetrics of the records it keeps.
- `-o`, `--compression-level`, `--create-index` and `--create-md5-file` belong to `write`.
- `position-based-downsample-sam` reads its input twice: since wave W14 it is a
  materialized step that runs on a file of its stream (section 31).

Chained equals staged (`crates/cli/tests/downsample.rs`): each strategy first, in the
middle and last beside `mark-duplicates` and `sort-sam` on the fixture, and first on t1;
the records as bytes, the header lines and every metrics file equal.

## 22. Steps without +

The owner's goal of 2026-11-22: the chain command line is too long. Both forms below run
the same chain and write the same files, byte for byte (`chain_short_form` and
`chain_short_form_long` in `tests/golden`):

```bash
# The + form, as before
bammetrics chain fixed.bam -R hs37d5.fa -o c/r \
    + reorder-sam --sequence-dictionary hs37d5.dict + insert-size + write c/r.bam --index

# The short form: -R gives the reference, and the step's dictionary is its .dict
bammetrics chain fixed.bam reorder-sam insert-size -R hs37d5.fa -o c/r.bam

# The same, with the line reference = "hs37d5.fa" in the config file (section 24)
bammetrics chain fixed.bam reorder-sam insert-size -o c/r.bam
```

The rules, applied to each word in turn (`crates/cli/src/chain/bare.rs`):

- `+` starts a step, and the word after it is the step's verb, as before.
- A bare word that names a step verb, a Picard class, `write` or `qc` starts a step.
- An option takes its values with it. `--as insert-size` keeps `insert-size` as the value
  of `--as`, so a value that names a verb stays a value.
- The first bare word in front of the first step is the source, unless `-i`, `--input`
  or `I=` gave it.
- A bare word one or two letters away from a verb is an unknown step, named by its
  number, with the closest verb (`step 2: unknown step "quality-yeld"`, `did you mean
  quality-yield?`). A verb the word starts ranks first: `ins` is `insert-size`. A word
  with a `/` or a `.` is never taken for a misspelled verb.
- After `--`, every word up to the next `+` is positional.
- The chain's options go in front of the first step, or after the last step when a bare
  word started that step. There an option goes to the chain when the step has no flag of
  its spelling, or when the step's flag means what the chain's does: `-o`, `--format`,
  `--threads`, the reference and `TMP_DIR`. A step option of another meaning stays with
  the step: `downsample-sam -R 7` is its `RANDOM_SEED`. Two steps keep their `-o` in
  last place, since it names what they write: `write` without a path (`write -o
  out.bam`) and a split step (`split-sam-by-library -o DIR`). Any other step's own `-o`
  in last place needs the `+` form. In the middle of the chain a step keeps its own
  `-o`, `--format` and `-R`: `insert-size --format picard,json alignment-summary`
  writes JSON for `insert-size` alone.
- `qc` in last place keeps its own options (`--preset`, `--report`, `--fail-on`, its
  `--also` of verdict forms) and hands on the rest: `qc --preset wgs -o r.bam --resume`
  gives `-o` and `--resume` to the chain (section 32).
- A bare word after `merge-sam-files`, `gather-bam-files` or `jumping-library` that is
  no verb is a further input of that step (`merge-sam-files second.bam`), and so is its
  `--input`, in the middle and in last place.
- A word that names a verb is a step, so a `write` path that is also a verb name is
  written `./NAME`.

A command line in the new form (a bare word started a step, or it has no `+` and no
`--plan`) gets its errors in its own terms, each exit 2 before the first read:

| Command line | Error |
|---|---|
| `chain insert-size in.bam` | `insert-size names a step: the source comes first` |
| `chain in.bam bogus` | `"bogus" is no step verb, and the chain takes one source` |
| `chain in.bam a.bam insert-size` | `one source expected, got 2`, hint `the source comes first, then the steps` |
| `chain in.bam insert-size mark` | `step 1 (insert-size): "mark" is no step verb, and a step takes no positional value`, `did you mean mark-duplicates?` |
| `chain in.bam insert-size -q alignment-summary` | `step 1 (insert-size): -q belongs to the chain: give it before the first step or after the last` |
| `chain in.bam mark-duplicates -o x.bam --md5` | `step 1 (mark-duplicates): --md5 is an option of write`, with the `write OUT.bam --md5` form |
| `chain -o a in.bam insert-size -o b` | `-o is given twice` |

The closest verb of a stray word is the shortest verb it starts, else the closest by
spelling, else the shortest verb that holds it, else the verb that shares its longest
start (4 letters or more). A command line with a `+` in front of every step comes back
unchanged and keeps the parser's errors, so every old form runs as before. Both forms mix: `chain in.bam + mark-duplicates quality-yield write
m.bam`. `plan` tells a chain from a `qc` plan by a `+`, `--plan` or `--recipe`; for the
short form write `bammetrics plan chain BAM STEP ...`.

## 23. -o names the result

`-o` ending in `.bam` (any case), `.sam` or `.cram` (lower case) on a chain with a
transform and no `write` step adds `write OUT --index` as the last step (`write OUT` for
`.sam` and `.cram`: SAM text has no index, and a CRAM index is refused as `CREATE_INDEX`
refuses it), and the prefix of every metrics file is the name without its extension:
`-o r.bam` writes `r.bam`, `r.bai` and `r.insert_size_metrics`
(`crates/cli/src/chain/resolve.rs`). A materialized tool that writes a BAM
(`merge-sam-files`, `gather-bam-files`, `position-based-downsample-sam`) counts as a
transform here. A `qc` step stays last: the write goes in front of it, and the manifest,
the verdict and the report take the same prefix (`r.manifest.json`, `r.verdict.json`,
`r.report.html`; section 32).

| `-o` | The chain | What it does |
|---|---|---|
| `out/s` (no extension) | any | The prefix, as before |
| `r.bam` | a transform, no `write` | `write r.bam --index` last, prefix `r` |
| `r.sam` | a transform, no `write` | `write r.sam` last, prefix `r` |
| `r.bam` | a `write` step | The prefix `r.bam`, as before |
| `r.bam` | taps only | Exit 2: `-o PREFIX names the metrics files` |
| `r.cram` | a transform, no `write` | `write r.cram` last, prefix `r`: a CRAM against the chain's reference (`-R`, else `BAMMETRICS_REFERENCE`, else `reference` of the config file); exit 2 before the first read without one |
| `r.SAM`, `r.Sam`, `r.CRAM` | a transform, no `write` | Exit 2: the write step, as htsjdk's writer, writes SAM text or a CRAM only for a name ending in lower-case `.sam` or `.cram`, and a BAM under any other name |

Two notes go to standard error before the first read. A chain of transforms that name
no file (`chain in.bam reorder-sam -R ref.fa -o out/s`: no tap, no `write`, no split, no metrics
file) gets `no step writes the records: -o OUT.bam or a write step writes them` on its
last step. A `.bam` or `.sam` prefix beside a `write` step of the
chain's own gets a note on that step: the prefix names the metrics files
(`-o r.bam` gives `r.bam.insert_size_metrics`).

The added `write` is the `write r.bam --index` of the long form. A stream that is not
sorted by coordinate gets that step's note (`CREATE_INDEX: no index written, the sort
order is queryname ...`) and no `.bai`.

## 24. Defaults from the config file and the reference

The user config file is `$XDG_CONFIG_HOME/bammetrics/config.toml` (an absolute
`XDG_CONFIG_HOME`), else `$HOME/.config/bammetrics/config.toml`:

```toml
reference = "~/refs/hs37d5.fa"
threads = 8
tmp_dir = "/scratch"
sort_memory = "4G"
```

- The file is read through `bammetrics_core::toml` with a cap of 64 KiB. Every key is
  checked; an unknown one is an error with the line and the closest name.
- `~/` is the home directory. A relative path is relative to the directory of the file.
- `BAMMETRICS_REFERENCE` replaces the file's reference. `-R` replaces both.
- `--no-config` leaves the file unread; the variable still counts.
- The file must be a regular file. A FIFO, a device or a directory named `config.toml`
  is a usage error before a byte is read; the file is opened without blocking.
- A reference from the file or the variable that is not a file fails a chain only when a
  step reads it: a tap whose `REFERENCE_SEQUENCE` or `REFERENCE` is a value option, a
  transform that requires one (`set-nm-md-and-uq-tags`), a `reorder-sam` without its
  dictionary, or any step on a CRAM source. That is an input error (exit 1) naming where
  the path came from, with the hint `--no-config` for the file and `unset
  BAMMETRICS_REFERENCE` for the variable. Any other chain runs without it. `qc`, `plan`
  and `collect` leave such a path out with a note on standard error (`bammetrics qc:
  note: /x.fa not found (the reference of ...): not used`).
- `bammetrics config` prints each value and where it came from (`crates/cli/src/config.rs`).

| Verb | Takes from the config |
|---|---|
| `chain` (and every recipe) | `reference`, `threads`, `tmp_dir`, `sort_memory` |
| `qc`, `plan` | `reference` |
| `collect` | `reference` |
| a tool verb | `reference`, when the tool reads one, or the input is a CRAM |
| the Picard face (`KEY=VALUE` after a class) | nothing: it stays Picard's command line |

A step option that takes the reference itself gets the chain's `-R` when the step names
none, as before (section 2.1): `REFERENCE_SEQUENCE` of `alignment-summary`,
`hs-metrics`, `gc-bias`, `sequencing-artifact`, `set-nm-md-and-uq-tags`,
`set-nm-and-uq-tags` and every record-editing step (for a CRAM source), and `REFERENCE`
of `rrbs`. One step option names a dictionary:

| Step | Option | Default |
|---|---|---|
| `reorder-sam` | `--sequence-dictionary` | The reference's `.dict` (its FASTA extension replaced, as htsjdk finds it), else its `.fai` |

The step's notes say which file it took:

```text
bammetrics chain: step 1 (reorder-sam): SEQUENCE_DICTIONARY: the dictionary is the reference's index ref.fa.fai (no ref.dict beside the reference): names and lengths only
```

A dictionary read from a `.fai` holds `SN` and `LN` alone; the `.dict` adds `M5` and `UR`.
Picard 3.5.0 refuses a `.fai` as `SEQUENCE_DICTIONARY`; bammetrics reads one
(KNOWN_DIFFERENCES.md). A `reorder-sam` step with no dictionary and no reference is a
usage error that names the ways to give one (`--sequence-dictionary`, `--from-bam`,
`--dict`, `-R`, `reference` in the config file, `BAMMETRICS_REFERENCE`). The step's own `-R` comes first, then the chain's: in
the middle of a bare chain `reorder-sam -R ref.fa insert-size` takes the dictionary of
`ref.fa`. The default applies in a chain only; the `reorder-sam` verb still needs its
`--sequence-dictionary`, as before.

## 25. Recipes

A recipe is a named chain in the plan file format (section 3), with a `[recipe]` table
and no source. Three are compiled in from `crates/cli/recipes/`; the user's are read from
the `recipes` directory beside the config file (`~/.config/bammetrics/recipes/NAME.toml`,
64 KiB each), and one with the name of a shipped recipe replaces it.

| Recipe | Steps | `-o` |
|---|---|---|
| `prep` | `fix-mate-information`, `mark-duplicates`, `sort-sam --sort-order coordinate`, `quality-yield`, `alignment-summary`, `insert-size`, `duplicate-metrics` (the taps of qc's generic preset), then the BAM with its index | the BAM: `-o out.bam` |
| `qc-only` | the four taps of qc's generic preset | the prefix |
| `downsample` | `downsample-sam --strategy ConstantMemory`; needs `downsample-sam.probability`; then the BAM with its index | the BAM: `-o out.bam` |

```bash
bammetrics recipes                                     # one line each, with its steps
bammetrics prep in.bam -o out.bam                      # the shorthand verb
bammetrics chain in.bam --recipe prep -o out.bam       # the same chain
bammetrics downsample in.bam -o d.bam downsample-sam.probability=0.1
bammetrics chain in.bam --recipe prep -o out.bam --dry-run   # its plan
```

```toml
[recipe]
description = "one line for bammetrics recipes"
output = "bam"                                 # optional: -o must name the BAM
requires = ["downsample-sam.probability"]      # optional

[[step]]
verb = "downsample-sam"
strategy = "ConstantMemory"
```

- A word `STEP.option=value` sets that option on every step of the verb `STEP` in the
  recipe, in place of the recipe's value. `STEP` is a verb or its Picard class, `option`
  any spelling of the step's option. A recipe without such a step, or a step without such
  an option, is a usage error with the closest name.
- Any word `STEP.option=value` whose two parts are plain words is an override, so a
  misspelled step gets the closest step of the recipe (`recipe prep has no step
  mark-dups`, `did you mean mark-duplicates?`). A `STEP.option` word without `=` whose
  step names a verb exits 2 with `needs =VALUE`. An error in a step an override set names
  the word beside the recipe line (`recipe prep:21 with sort-sam.sort-order=a=b`).
- A recipe gives the steps: a step on the command line beside `--recipe` exits 2, and so
  does `--plan` beside it. A second `--recipe` exits 2. One step may follow a recipe:
  `qc`, the last step of any chain (`bammetrics prep in.bam -o out.bam qc --preset
  generic`, section 32). A recipe file may also end in a `qc` step, and `qc.option=value`
  sets its options; a second `qc` after such a recipe exits 2.
- A user recipe is a file the chain reads, and so is the config file: a step that would
  write either exits 2 (`is an input; writing it would destroy it`).
- `STEP.option=value` without a recipe exits 2.
- `bammetrics recipes` lists a user recipe that cannot be read as `NAME  (yours) not
  usable: REASON` and goes on with the rest. A user file that is not a regular file is
  not listed. A user recipe named like a verb (`sort-sam`) is listed with `bammetrics
  sort-sam runs the verb`: run it with `chain BAM --recipe sort-sam`.
- `bammetrics NAME --help` prints the recipe's description, its steps with their options
  and its required options.
- `output = "bam"` makes `-o` name a BAM; the implicit `write` of section 23 writes it.

## 26. Departures from docs/07 (the simpler grammar)

The decision rows of the simpler grammar, D160 to D174 of docs/07. Rows D199 to D205 hold
what changed when the grammar met the manifest and the materialized steps (section 32):

| Row | docs/07 or the task | As built | Reason |
|---|---|---|---|
| D160 | Chain options "after the last step keep their meaning" | In the last bare step an option goes to the chain when the step has no flag of its spelling or the step's flag has the chain's meaning (`-o`, `--format`, `--threads`, the reference, `TMP_DIR`); others stay with the step | The owner's example puts `-o` after `insert-size`, which has its own `-o`. `downsample-sam -R` is `RANDOM_SEED`, not the reference |
| D161 | "indexed when the stream is coordinate sorted" | The added step is `write OUT --index`: an unsorted stream gets that step's note and no `.bai` | The sort order is known when the stream reaches `write`. The new form then equals the long form byte for byte, notes included |
| D162 | The -o rule for every chain | A chain with its own `write` keeps `-o r.bam` as its prefix | Today's forms stay unchanged byte for byte |
| D163 | "the dictionary read from its .fai" | `ReorderSam` reads a `.fai` as `SEQUENCE_DICTIONARY` (names and lengths); Picard 3.5.0 exits on one | No file is written for the dictionary. The project has the precedent of `CollectWgsMetrics` on a `.fai` without a `.dict` |
| D164 | Config `threads`, `tmp_dir`, `sort_memory` for every verb | Those three apply to `chain` and the recipes; `qc`, `plan`, `collect` and the tool verbs take `reference` | `qc --threads` has a clap default, and its Python keywords come from the clap arguments; the tool verbs keep their documented defaults |
| D165 | `--no-config` on every verb | `qc`, `plan`, `collect` and the tool verbs take it out of their arguments before parsing; `chain --help` lists it | Their `--help` texts and golden records stay unchanged |
| D166 | Override "STEP.option=value or the step repeated" | `STEP.option=value` | One form; a step repeated after a recipe exits 2 |
| D167 | Config paths | `~/` is `$HOME`; a relative path is relative to the file's directory | The file names the same reference from any working directory |
| D168 | Recipe names | A user recipe replaces a shipped one of its name | A user changes `prep` without a new name |
| D169 | Config reference on the tool verbs | Only when the tool reads the reference, or the input is a CRAM | A BAM run of a tool that ignores `REFERENCE_SEQUENCE` prints no new note |
| D170 | Python `config=False` | On `chain` and `recipe` only; Python `qc` reads no config | The task names chain and recipe; a notebook's `qc` keeps its explicit keywords |
| D171 | D160: a step's own `-o` in last place needs the `+` form | `write` without a path and a split step keep their `-o` in last place | `write -o out.bam` sent `-o` to the prefix and wrote `out.bam.marked.bam`; `split-sam-by-library -o DIR` wrote into `.` |
| D172 | A config or variable reference that is not a file exits 1 | Only a chain with a step that reads it exits 1; `qc`, `plan` and `collect` leave it out with a note | A stale config failed every chain, `insert-size` alone included |
| D173 | `-o` ending in `.sam` in any case | `.SAM` exits 2 | The write step writes SAM text for lower-case `.sam` only, as htsjdk's `SAMFileWriterFactory`; `-o r.SAM` wrote BGZF bytes with no `.bai` |
| D174 | The dictionary default for every `-R` user | In a chain only | The `reorder-sam` verb keeps its command line and golden records; its error names `--from-bam` and `--dict`, and so does the chain's |

Tests: `crates/cli/src/chain/bare.rs`, `recipes.rs`, `resolve.rs` and
`crates/cli/src/config.rs` (unit), `crates/cli/tests/chain_grammar.rs` (each new form
against its long form, byte for byte; the reference order; the Picard face untouched;
recipes, user recipes and the 64 KiB cap; a stale reference, a FIFO config file, a
broken user recipe, the last writer's `-o`), the 41 records of `CHAIN_GRAMMAR_CASES` in
`crates/cli/tests/golden.rs`, and `tests/python/test_chain.py`.

## 27. W14: the manifest and `--resume`

A chain writes `PREFIX.manifest.json` (schema `bammetrics.manifest/2`,
[docs/schema/manifest-2.json](schema/manifest-2.json)) when it resumes, keeps a stream or
ends in `qc`. A chain that asks for none of these writes the files it wrote before.

```bash
bammetrics chain in.bam -o out/s --resume \
    + quality-yield + mark-duplicates --keep \
    + sort-sam -SO queryname + insert-size + write
# Edit step 3 and run again: steps 1 and 2 are skipped, the read starts at
# out/s.marked.bam, steps 3, 4 and 5 run.
```

The manifest holds the source's fingerprint (size, mtime, the FNV-1a 64 of its first and
last 64 KiB), the reference's, the validation stringency, when the run started, its wall
time and exit code. Then one record per step:

| Member | What it is |
|---|---|
| `id` | The stream of a transform or of a materialized tool that writes a BAM (`marked`, `merged`), else `VERB@STREAM` (`insert-size@marked`), `#2` for the second of one verb on one stream; `qc` for the qc step |
| `upstream` | The ids of the transforms in front of the step (a materialized tool that writes a BAM counts as one) |
| `argv` | The canonical argv: a tool's Picard arguments (with `KEEP=PATH` for a kept stream), `write`'s and `qc`'s options as `KEY=VALUE` |
| `options_hash` | FNV-1a 64 of the argv, the tool version, the formats, `--sample` and the stringency |
| `inputs` | Every file the step reads by name, with its fingerprint: the reference and its `.fai`, dictionaries, interval and read lists |
| `outputs` | Every file it wrote, with its kind and its fingerprint after the write. A `write` step's file is `bam`, `sam` or `cram` by its extension, its companions `index` and `md5` |
| `status`, `error`, `skipped` | `ok`, `failed`, `stopped` or `not run`; `skipped` is true when the run kept the record of the run before |
| `started`, `seconds` | The start and wall time of the run that ran the step |

Then `kept`: each kept stream with its transform's id and its BAM's fingerprint.

`--resume` skips a step when the record under its id has the same verb and hash, ended
`ok`, read files with the same fingerprints, and every file it wrote is there with its
fingerprint; and when every transform in front of it is skipped too. A tap in front of
the changed step is not upstream of it, so it stays skipped. The qc step and a step that
writes standard output always run. Every step runs when there is no manifest, when it is
the manifest of a `qc` run (schema 1), when it is damaged (a note says why, the run goes
on), or when the tool version, the source, the reference or the stringency changed.

The steps that run share reads. The read starts at the last transform in front of the
first step that runs whose stream is in a BAM the run before wrote and that is still
there with its fingerprint: the stream it kept (`--keep`), or the BAM a `write` step of
that stream wrote, when that `write` step is skipped in this run. So `-o out.bam` is a
place to resume from without `--keep` (decision D201). A `write` step that wrote SAM text
or a CRAM is skipped as a BAM's is, and the read never starts at its file. Without such a BAM the read starts
at the source, and every transform between the start and the last step that runs runs
again to feed the steps after it. A run where every step is skipped reads nothing. The
plan names the start: `skipped; its kept stream out/s.marked.bam starts the read`, or
`skipped; the BAM a write step wrote of its stream, out/r.bam, starts the read`.

`plan --resume` and `chain --resume --dry-run` print one line per step after the plan:

```
resume     from out/s.manifest.json
resume     step 1 quality-yield: skipped: its record and files match
resume     step 2 mark-duplicates: skipped; its kept stream out/s.marked.bam starts the read
resume     step 3 sort-sam: runs: its options changed
resume     step 4 insert-size: runs: step 3 in front of it runs
resume     step 5 write: runs: step 3 in front of it runs
```

A skipped step prints `skipped: the manifest proves it done, its files stay`, and the
summary counts it: `5 steps finished, 2 of them skipped by --resume, 3 files written`.

One file has one name. `--resume` compares paths resolved: absolute, links followed, `.`
and `..` taken out. The manifest records the run's working directory (`directory`), and
its other paths are as the run named them, relative to it. So `m.bam`, `./m.bam`,
`$PWD/m.bam`, a link to it and `../m.bam` from a subdirectory name one source, and an
`-o` spelled another way names one prefix. The options hash reads each path option
resolved. Another file under the source's name prints `--resume: the run before read
m.bam; this run reads other.bam: every step runs`.

Fingerprints. The source and the reference are fingerprinted before the first read, and
each output as soon as the run wrote it. A source replaced while a run reads it therefore
shows as changed to the next `--resume`, which reruns every step.

One run per prefix. A chain that writes a manifest holds an advisory lock (`flock`) on
`PREFIX.manifest.json.lock` from before its first read to its end, and removes the file
then (on SIGINT, SIGTERM and SIGHUP too). A second run on the same prefix exits 2 before
it reads: `another chain is writing out/s.manifest.json (it holds
out/s.manifest.json.lock): one run per output prefix`. A lock file a killed run left is
no lock: `flock` ends with the process.

File safety. The manifest is read under a 16 MiB cap through the JSON reader with a
nesting limit of 64, so a damaged file is a note, never a crash. Its directory is checked
before the first read by making and removing the temporary file it is written under; a
directory that takes no new file fails with exit 1 and names the directory: `the manifest
out/s.manifest.json: cannot write in out: Permission denied`. A step that runs replaces its
files atomically; a skipped step's files are not opened for writing.

After a resumed run that finished, each file the manifest of the run before lists that
this run neither wrote nor kept is removed, and named: `--resume: removed out/s.sorted.bai:
the run before wrote it and this run does not`. It goes only while it is a regular file
(no link is followed) with the fingerprint recorded, the chain does not read it, and it is
a `PREFIX.*` file in the prefix's directory or the `.bai` or `.md5` of a BAM the same step
wrote in this run. A planned file this run did not write is stale too: the `.bai` of a
kept stream that is no longer coordinate sorted goes. These stay, each named with why:

- the files of a `write` or split step this run does not have (`--resume: left in place:
  deliver/final.bam (a write or split step's file; this run has no such step)`), so
  dropping `write` to rerun faster keeps the BAM;
- a file outside the prefix's files, which only a hand-edited manifest can list;
- a file whose recorded kind is none a chain writes, or a `bam`, `kept`, `sam`, `index`
  or `md5` record whose name lacks that kind's extension (`--resume: left in place:
  out/s.insert_size_metrics (its kind is none a chain writes under that name)`);
- a file changed since the run before wrote it.

The manifest is opened once without following a link, and its owner, mode and text come
from that one descriptor. A manifest that is a symbolic link, is not a regular file, is
owned by another user or may be written by group or others is not read, so every step
runs and no file is removed: `out/s.manifest.json has mode 664: group or others may write
it, so it is not read: every step runs`.

A run that fails still writes its manifest, on purpose: its record of the failed step makes
the next `--resume` rerun that step, and the steps that finished are skipped. It removes no
file of the run before.

A run killed with SIGINT, SIGTERM or SIGHUP leaves no temporary file and the files of the
run before as they were (exit 130 on SIGINT). A run killed with SIGKILL cannot clean up:
its temporary files stay as `.NAME.PID.N.HASH.tmp` beside the prefix, up to the size of
the file being written. The next `--resume` names each one whose process is gone
(`--resume: left in place: out/.s.marked.bam.4571.0.e6d139e96494597c.tmp (a temporary file
of a run that was killed; remove it)`) and leaves it: another program may own a name of
that form.

Tests: `crates/cli/tests/chain_resume.rs` (the gate: an option of step 3 of five reruns
steps 3, 4 and 5 and every file equals a fresh run; the feed without a kept stream; a
changed source; a damaged manifest; a hand-edited file left alone; a failed step that
reruns), `crates/cli/tests/chain_resume_safety.rs` (a respelled source or prefix, a
source swapped during a run, the lock, what a resumed run removes and leaves, a grafted
record, a record whose kind does not fit its name, a linked or group-writable manifest, an
unwritable directory, a killed run's temporary file, the examples of
docs/cli.md and the README), the golden records `chain_resume_*` and
`plan_chain_resume_option_change`.

## 28. W14: kept streams

`--keep` on a transform writes the stream it leaves as `PREFIX.<stream>.bam`, with
`PREFIX.<stream>.bai` when the stream is coordinate sorted. A materialized tool that
writes a BAM takes it too (`merge-sam-files --keep` writes `PREFIX.merged.bam`). The file
is the BAM `write` writes at that position, byte for byte, and its index Picard's
`BuildBamIndex` index of it. `--keep` on the last transform beside `-o out.bam` writes the
same records twice (`PREFIX.<stream>.bam` and `out.bam`); `--resume` starts at either.
No note is printed for a stream that is not coordinate sorted: no index is the rule there.
A plan file writes `keep = true` in the step's table; Python passes `keep=True` to `step`.

The kept BAM is a planned file of its transform: it is claimed, checked for writability and
against the inputs before the first read, so `+ mark-duplicates --keep + write` (both
`PREFIX.marked.bam`) is a usage error: `step 1 (mark-duplicates) --keep and step 2 (write)
both write out/s.marked.bam: give write its own path, or drop --keep`. `--keep` before the
first `+` is a usage error that names the transform step, and `keep = true` on a tap in a
plan file gets the message of `--keep` on a tap. The plan lists it as `PREFIX.marked.bam (kept)` and
its index as `(when coordinate sorted)`, and the output line counts its bytes: `output
about 2.7 KiB written by write and the kept streams`.

A kept stream read by `--resume` in place of a CRAM source keeps the CRAM's frame
(`FileSource::decoded_from_cram`), so the stages after it rebuild `MD` and `NM` as they did.

## 29. W14: `qc` as the last step

```bash
bammetrics chain in.bam -o out/s + mark-duplicates + write out/s.bam + qc --preset generic
```

The qc step evaluates the preset over the metrics files of the taps on its stream and
writes `PREFIX.verdict.json` and `PREFIX.report.html` (`--report pdf` adds
`PREFIX.report.pdf`) as `bammetrics qc` writes them for the same taps on the same records:
the same `verdict.json`, and the same page apart from the date and each metrics file's
command line. A collector the preset needs that no tap of the stream runs is added as a
tap, built as the qc verb builds it (its file `PREFIX.<extension>`, its levels, its
reference) and numbered after the qc step: `step 3 quality-yield tap on marked, for qc`.
The qc step keeps the number it was given. A transform's own metrics file
(`PREFIX.marked.duplication_metrics`) is no tap's, so it does not enter the verdict.

| Option | Meaning |
|---|---|
| `--preset NAME\|FILE` | The rules (default `generic`) |
| `--report FORM[,FORM]` | `html` (default), `pdf`, or both |
| `--also FORMAT[,FORMAT]` | The forms of the verdict: `markdown`, `junit`, `multiqc`. A per-file format is the chain's `--also` and a usage error here |
| `--fail-on warn\|fail` | The chain exits 10 when the verdict reaches it, 11 when it is incomplete; 1 still wins when a step failed |
| `--hs BAITS TARGETS`, `--panel FILE`, `--ref-flat FILE`, `--strand`, `--ribosomal-intervals`, `--sex` | As for `qc`: the capture, RNA and sex check collectors |
| `--all-plots`, `--date`, `--pdf-font`, `--page-size` | As for `qc` |

The source, `-o`, `-R`, `--threads`, `--sample` and `--resume` belong to the chain. A step
after `qc` is a usage error. The step prints its verdict on standard error: `step 2 (qc):
verdict pass: 6 checks, 6 pass, 0 warn, 0 fail, 0 not measured`. A class no step can write
is named in the plan (`qc         WgsMetrics not measured: needs -R FASTA`); a class a tap
of the stream writes is measured, so an `hs-metrics` step in front gives `HsMetrics`.
The `sex-check` tap the qc step adds counts reads per base over each whole contig: it does
not take `TARGET_INTERVALS` from an `hs-metrics` step, where `qc --hs` gives it the
targets, so its `X_READS_PER_BASE` and `Y_READS_PER_BASE` differ from the verb's.

`CollectWgsMetrics` has no one-pass collector. Under a qc step it reads a file of the
stream on its own after the chain: the source when no transform is in front, else the
stream the last transform keeps (`--keep`), else the BAM a `write` step writes of that
stream (the BAM `-o out.bam` names). With no such file it is a materialized step of the
chain's read (section 31): it runs on a copy of the stream in the temporary directory,
the plan lists it as `materialize on marked, for qc` and counts its disk. A `wgs` step
the chain already runs on the stream is not added twice, and the verdict reads its file
(decision D202). A `--panel` BED is written as
`PREFIX.panel.interval_list` from the source's dictionary after the chain, as `qc` writes
it; a piped source with a BED panel is a usage error.

The chain writes `PREFIX.manifest.json` with schema 2; the qc verb keeps writing schema 1.
`--resume` always reruns the qc step, so the verdict and the report are rewritten over
the files of the steps it skips.

Tests: `crates/cli/tests/chain_resume.rs` (the qc step after `mark-duplicates` against
`qc` on the BAM the chain wrote), the golden records `chain_qc_step*` (the four metrics
files and `verdict.json` of `chain_qc_step` equal those of `qc_generic`),
`tests/python/test_chain.py::test_qc_is_the_last_step`.

## 30. W14: departures from docs/07

| Item | docs/07 | Built | Why |
|---|---|---|---|
| When manifest 2 is written | Every chain and qc run | A chain that resumes, keeps a stream or ends in `qc`; `qc` keeps schema 1 | The golden records of every earlier chain form stay byte for byte |
| `Key.upstream` | A hash of every earlier step's key | The ids of the transforms in front, and a step is skipped only when they are | A hash over every earlier step reruns the taps after a changed tap; the gate wants taps in front of step 3 kept |
| Timing | Per step | The start and wall time of the run that ran the step | The steps of a chain share reads and run at once |
| `qc.rs` | `Plan.bam` becomes a chain, `Read` gains `Tap(n)` | The qc verb is unchanged; the qc step reuses its collector planning (`qc::plan_chain_step`) | qc golden records stay byte for byte |
| Kept product kind | A kept product may be CRAM | BAM, with its index when coordinate sorted | CRAM writing shipped in W15, 0.33.0, as `write r.cram` (section 23); `--keep` writes a BAM (section 28) |
| Materialize fallback | A tool with only `run` reads a spooled copy of the stream | Built in section 31. The qc step's `CollectWgsMetrics` reads the source, a kept stream or the BAM a `write` step wrote, else it is a materialized step (D202) | No copy is written when a file already holds the stream |
| Paths and fingerprints | The manifest's fingerprints | Paths compared resolved against the recorded `directory`; the source fingerprinted before the read, each output when written | A respelled path reran every step; a source swapped during a run was recorded with its new fingerprint |
| Concurrent runs | Not covered | One run per prefix through `PREFIX.manifest.json.lock`; a second exits 2 | Two runs on one prefix wrote one kept BAM and a manifest that skipped over the other run's records |
| What a resumed run removes | Every file the old manifest lists and this run does not write | `PREFIX.*` files and the companions of a BAM the same step rewrote; a dropped `write` step's files and any file elsewhere are left and named | A dropped `write` lost its BAM; a hand-edited manifest removed a file of another prefix |

## 31. W14: every registered tool at any position, the materialize fallback

Every registered tool that reads a BAM runs at any position of a chain: 53 of the 71
tools at 0.32.0, 53 of the 74 at 0.33.0, whose three new tools are refused as steps with
their reasons (section 3; docs/07 D224 and D233). A tool with a stage is a transform, a one-pass tool is a tap, a split tool is a
split. A tool with neither a stage nor a collector is materialized (docs/07 section
2.3): it runs on a file that holds the records of its position.

```bash
# wgs and library-complexity behind fix-mate-information, the BAM merged with a second
# lane, then insert size on the merged stream and the merged BAM written
bammetrics chain lane1.bam -o out/s \
    + fix-mate-information \
    + wgs -R ref.fa + library-complexity \
    + merge-sam-files lane2.bam \
    + insert-size \
    + write out/s.bam --index
```

### 31.1 What materialize does

- On the source with no stage in front, and the source a regular file, the tool reads
  the source itself: nothing is copied.
- Anywhere else the records of the position are written to a BAM in the chain's
  temporary directory (`--tmp-dir`, or the system one): BGZF level 3, by the writer
  `write` uses, with htsjdk's order check. `bam-index-stats` and `index-depth` get its
  `.bai` when the stream is coordinate sorted. `calculate-read-group-checksum` reads the
  header alone: the file holds no record.
- The tool then runs through its own `run` on that file, so its output is the output of
  `chain ... + write x.bam` and then the tool on `x.bam`. The file is removed when the
  tool is done, or with the temporary directory.
- A tool that writes a BAM (`merge-sam-files`, `gather-bam-files`,
  `position-based-downsample-sam`) is a barrier. Its stream is its first `INPUT`, and
  files after the verb (or `--input FILE`, `input = [...]` in a plan file, `input=[...]`
  in Python) are the further inputs of `merge-sam-files`, `gather-bam-files` and
  `jumping-library`. The tool writes its BAM in the temporary directory at level 3, and
  that BAM is the stream of the next segment, named `merged`, `gathered` or
  `downsampled` (`--as NAME` renames it, as on a transform; `merge-sam-files` takes its
  `ASSUME_SORTED` as `--assume-sorted`, as `mark-duplicates` does). Its `-o`, `--compression-level`, `--create-index` and
  `--create-md5-file` belong to `write` (exit 2).
- A materialized tool that writes no BAM fails alone, as a tap does. A materialized
  tool that writes a BAM fails as a barrier: the later steps are not run.
- A materialized tool that writes a BAM is a transform to the rest of the chain: it takes
  `--keep` (`PREFIX.merged.bam`), `--resume` starts a read at its kept stream and reruns
  the steps behind it when it runs, `-o out.bam` writes its records, and a `qc` step
  reads the taps of its stream. A materialized tool that writes no BAM is a tap to the
  qc step: the verdict reads its metrics file (`wgs` in front of `qc --preset wgs`).
  Section 32 holds the cases.
- `check-duplicate-marking` writes its file and then fails (exit 1) when a record
  disagrees, as Picard exits 1.
- The command line the tool writes (`@PG CL`, the `StringHeader`) names the step's input
  as the chain planned it, never the file in the temporary directory.
- Each materialized step writes its own copy. Three such steps on one stream write three
  copies, so the disk and the compression grow with the count of steps; the plan lists
  each copy.
- A tool that spills (`check-duplicate-marking`, `umi-prevalence`, the sort of
  `merge-sam-files`) spills under `--tmp-dir` too, unless the step gives its own
  `--tmp-dir`. The planner counts the copy and the BAM a tool writes, not the tool's own
  spill.
- Before the first read, as for every step: the tool's own option parser judges every
  value (`check-duplicate-marking --mode NOPE`, `multiple-metrics --level NOPE`, exit 2);
  each further input must be a BAM, CRAM (with the chain's `-R`) or SAM file, and
  `merge-sam-files` takes at most 1,000 inputs, the stream among them; a
  `bam-index-stats` or `index-depth` step that reads the source itself needs the
  source's `.bai` or `.csi` (`index-depth` also refuses a CRAM there); the files every
  `multiple-metrics` program writes are claimed, so a tap that writes one of them
  (`insert-size`) is a usage error; the `.bai` or `.csi` a step reads beside its inputs is
  an input no step may write; and `--tmp-dir` must be a directory the user may write, or a missing one whose nearest existing parent is.
- `bam-index-stats` and `index-depth` on a stream that is not coordinate sorted fail
  alone when the stream starts: its copy gets no `.bai`. `sort-sam --sort-order
  coordinate` in front gives it one.
- An error of the copy names the stream (`markdup.bam + step 1 (gather-bam-files)`). A
  stream whose records break the order its header says fails the next materialized step
  with `the stream ... is not in the order its header says`.

### 31.2 Names

A materialized tap writes `PREFIX<extension>` on the last stream and
`PREFIX.<stream><extension>` in front of a transform, as a tap does, or its own `-o
FILE`. `-o -` prints `view-sam` and `bam-index-stats` on standard output.

| Verb | Extension |
|---|---|
| `wgs` | `.wgs_metrics` (the name `qc` writes) |
| `multiple-metrics` | the prefix itself: each program adds its extension |
| `raw-wgs` | `.raw_wgs_metrics` |
| `wgs-with-non-zero-coverage` | `.wgs_non_zero_coverage_metrics` |
| `targeted-pcr` | `.targeted_pcr_metrics` |
| `library-complexity` | `.library_complexity_metrics` |
| `oxo-g` | `.oxog_metrics` |
| `bam-index-stats` | `.bam_index_stats` |
| `index-depth` | `.index_depth_metrics` |
| `check-duplicate-marking` | `.duplicate_marking_check` |
| `jumping-library` | `.jumping_library_metrics` |
| `view-sam` | `.view.sam` |
| `calculate-read-group-checksum` | `.read_group_md5` |
| `umi-prevalence` | `.umi_prevalence_metrics` |

### 31.3 The plan

`--dry-run` and `plan` show each materialized step and its disk, and the free-space
refusal counts it in the temporary directory:

```
$ bammetrics chain tests/data/markdup.bam -o out/s --dry-run + raw-wgs -R ref.fa \
    + add-comments-to-bam --comment c + raw-wgs -R ref.fa + merge-sam-files tests/data/markdup.bam \
    + insert-size
segment 1  one read of tests/data/markdup.bam
  step 1  raw-wgs              materialize on in         out/s.in.raw_wgs_metrics
  step 2  add-comments-to-bam  stream commented          no file
  step 3  raw-wgs              materialize on commented  out/s.commented.raw_wgs_metrics
  step 4  merge-sam-files      materialize, first read   writes its input to a file, then runs on it
segment 2  a read of the BAM step 4 wrote (materialize)
  step 4  merge-sam-files      stream merged             no file
  step 5  insert-size          tap on merged             out/s.insert_size_metrics
materialize step 4 merge-sam-files: writes its input to a file, about 1.2 KiB; writes its BAM, about 2.3 KiB
materialize step 1 raw-wgs: reads the source file itself, nothing is copied
materialize step 3 raw-wgs: writes its stream to a file, about 1.2 KiB
temp disk  about 4.7 KiB in the system temporary directory
```

The copy takes the bytes per record of a spool (76 in coordinate order, 105 in another,
section 15.3) times the records the source size gives. The BAM of a tool that writes one
takes the same per record plus the size of every further input. The planner sums every
copy as if all were held at once. The spill of a tool that sorts on its own
(`check-duplicate-marking`, `umi-prevalence`, `merge-sam-files`) is not counted
(decision D198).

### 31.4 Every tool and its role

| Count | Role |
|---|---|
| 20 | tap |
| 14 | transform |
| 2 | split |
| 17 | materialize (3 of them leave a stream) |
| 18 | refused: 17 read no BAM, and `build-bam-index` |

| Verb | Picard class | Role in a chain |
|---|---|---|
| `insert-size` | `CollectInsertSizeMetrics` | tap |
| `alignment-summary` | `CollectAlignmentSummaryMetrics` | tap |
| `quality-yield` | `CollectQualityYieldMetrics` | tap |
| `hs-metrics` | `CollectHsMetrics` | tap |
| `gc-bias` | `CollectGcBiasMetrics` | tap |
| `quality-score-distribution` | `QualityScoreDistribution` | tap |
| `mean-quality-by-cycle` | `MeanQualityByCycle` | tap |
| `base-distribution-by-cycle` | `CollectBaseDistributionByCycle` | tap |
| `duplicate-metrics` | `CollectDuplicateMetrics` | tap |
| `library-complexity` | `EstimateLibraryComplexity` | materialize |
| `multiple-metrics` | `CollectMultipleMetrics` | materialize |
| `accumulate-quality-yield` | `AccumulateQualityYieldMetrics` | refused: reads no BAM |
| `wgs` | `CollectWgsMetrics` | materialize |
| `raw-wgs` | `CollectRawWgsMetrics` | materialize |
| `wgs-with-non-zero-coverage` | `CollectWgsMetricsWithNonZeroCoverage` | materialize |
| `targeted-pcr` | `CollectTargetedPcrMetrics` | materialize |
| `rna-seq` | `CollectRnaSeqMetrics` | tap |
| `variant-calling` | `CollectVariantCallingMetrics` | refused: reads no BAM |
| `accumulate-variant-calling` | `AccumulateVariantCallingMetrics` | refused: reads no BAM |
| `genotype-concordance` | `GenotypeConcordance` | refused: reads no BAM |
| `bed-to-interval-list` | `BedToIntervalList` | refused: reads no BAM |
| `compare-metrics` | `CompareMetrics` | refused: reads no BAM |
| `create-sequence-dictionary` | `CreateSequenceDictionary` | refused: reads no BAM |
| `mark-duplicates` | `MarkDuplicates` | transform, barrier, replays |
| `sequencing-artifact` | `CollectSequencingArtifactMetrics` | tap |
| `convert-sequencing-artifact-to-oxo-g` | `ConvertSequencingArtifactToOxoG` | refused: reads no BAM |
| `oxo-g` | `CollectOxoGMetrics` | materialize |
| `coverage-gaps` | `CollectCoverageGaps` | tap |
| `sex-check` | `CollectSexCheckMetrics` | tap |
| `bam-index-stats` | `BamIndexStats` | materialize |
| `check-duplicate-marking` | `CheckDuplicateMarking` | materialize |
| `jumping-library` | `CollectJumpingLibraryMetrics` | materialize |
| `build-bam-index` | `BuildBamIndex` | refused: `write --index` writes the `.bai` |
| `fifo-buffer` | `FifoBuffer` | refused: a chain has no pipe |
| `view-sam` | `ViewSam` | materialize |
| `sort-sam` | `SortSam` | transform, barrier, holds |
| `mt-copy-number` | `CollectMtCopyNumberMetrics` | tap |
| `exon-coverage` | `CollectExonCoverageMetrics` | tap |
| `adapter-content` | `CollectAdapterContentMetrics` | tap |
| `duplicate-by-tile` | `CollectDuplicateByTileMetrics` | tap |
| `telomere` | `CollectTelomereMetrics` | tap |
| `index-depth` | `CollectIndexDepthMetrics` | materialize |
| `add-comments-to-bam` | `AddCommentsToBam` | transform, streams |
| `replace-sam-header` | `ReplaceSamHeader` | transform, streams |
| `add-o-a-tag` | `AddOATag` | transform, streams |
| `clean-sam` | `CleanSam` | transform, streams |
| `set-nm-md-and-uq-tags` | `SetNmMdAndUqTags` | transform, streams |
| `set-nm-and-uq-tags` | `SetNmAndUqTags` | transform, streams |
| `add-or-replace-read-groups` | `AddOrReplaceReadGroups` | transform, streams (a barrier when it sorts) |
| `revert-original-base-qualities-and-add-mate-cigar` | `RevertOriginalBaseQualitiesAndAddMateCigar` | transform, barrier, holds |
| `split-sam-by-library` | `SplitSamByLibrary` | split |
| `split-sam-by-number-of-reads` | `SplitSamByNumberOfReads` | split |
| `gather-bam-files` | `GatherBamFiles` | materialize, leaves `gathered` |
| `calculate-read-group-checksum` | `CalculateReadGroupChecksum` | materialize |
| `quality-yield-s-n-v-q` | `CollectQualityYieldMetricsSNVQ` | tap |
| `umi-prevalence` | `CollectUmiPrevalenceMetrics` | materialize |
| `rrbs` | `CollectRrbsMetrics` | tap |
| `create-baf-regress-file` | `CreateBafRegressMetricsFile` | refused: reads no BAM |
| `create-verify-i-d-intensity-contamination-file` | `CreateVerifyIDIntensityContaminationMetricsFile` | refused: reads no BAM |
| `arrays-variant-calling` | `CollectArraysVariantCallingMetrics` | refused: reads no BAM |
| `interval-list-tools` | `IntervalListTools` | refused: reads no BAM |
| `interval-list-to-bed` | `IntervalListToBed` | refused: reads no BAM |
| `vcf-to-interval-list` | `VcfToIntervalList` | refused: reads no BAM |
| `scatter-intervals-by-ns` | `ScatterIntervalsByNs` | refused: reads no BAM |
| `non-n-fasta-size` | `NonNFastaSize` | refused: reads no BAM |
| `merge-sam-files` | `MergeSamFiles` | materialize, leaves `merged` |
| `reorder-sam` | `ReorderSam` | transform, barrier, holds |
| `fix-mate-information` | `FixMateInformation` | transform, barrier, holds |
| `filter-sam-reads` | `FilterSamReads` | transform, streams (a barrier when it sorts) |
| `downsample-sam` | `DownsampleSam` | transform, streams |
| `position-based-downsample-sam` | `PositionBasedDownsampleSam` | materialize, leaves `downsampled` |

### 31.5 What equals what

`crates/cli/tests/chain_positions.rs` walks `TOOLS`. Each of the 53 tools that are steps
runs first (`SOURCE + TOOL`), in the middle (`SOURCE + C1 + TOOL + C2`) and last
(`SOURCE + C1 + C2 + TOOL`), where C1 and C2 are `add-comments-to-bam` steps, on a
fixture with its options. Each run equals the staged run: the stream in front written by
`write`, the tool's own verb on that file, then C2 on the tool's BAM. Metrics files are
equal after their `StringHeader` lines, BAMs in their header lines (the `CL` of `@PG`
aside) and every record as bytes, text files whole. The test fails when a registered
tool has no row in its table of steps and none in its table of refusals, and each of the
18 refusals must exit 2 with its reason.

### 31.6 Departures from docs/07

The decision rows D190 to D198 of docs/07 hold each departure with its reason:

- D190: a table of materialized tools (`stage::materialize`), no `Tool` field; the tool
  reads the source itself when it can.
- D191: "every registered tool runs at any position" means the 53 tools that read a BAM
  and are no `build-bam-index`.
- D192: the W09 to W13 tools without a stage (`gather-bam-files`,
  `calculate-read-group-checksum`, `umi-prevalence`, `merge-sam-files`,
  `position-based-downsample-sam`) are materialized, not given a stage.
- D193: `view-sam` is a materialized step, no longer refused.
- D194: `check-duplicate-marking` fails after its file.
- D195: the names under the prefix of the tools `qc` does not name.
- D196: the plan counts a copy at the bytes per record of a spool, and a merge's BAM
  adds the size of its further inputs; the steps after a merge keep the source's record
  count.
- D197: a second `INPUT` of `fix-mate-information` and `--write-reads-files` of
  `filter-sam-reads` stay refused, with a new reason.
- D198: one copy per materialized step, and the tool's own spill is not in the plan.

## 32. 0.32.0: the short form, the manifest and the materialized steps together

Sections 22 to 26 (the simpler grammar), 27 to 30 (the manifest) and 31 (materialize)
were built on three branches. This section holds the rules where they meet, each with
its test in `crates/cli/tests/chain_combined.rs` and its golden record. The decision rows
are D199 to D205 of docs/07.

```bash
# Before (0.31.0): every step behind a +, the write step named, a prefix beside it
bammetrics chain fixed.bam -R hs37d5.fa -o c/r \
    + reorder-sam --sequence-dictionary hs37d5.dict + insert-size + write c/r.bam --index

# 0.32.0: -R gives the reference, and the step's dictionary is its .dict
bammetrics chain fixed.bam reorder-sam insert-size -R hs37d5.fa -o r.bam

# Run it again after an edit: only what changed runs
bammetrics chain fixed.bam reorder-sam insert-size --minimum-pct 0.2 -R hs37d5.fa -o r.bam --resume

# With the line reference = "hs37d5.fa" in ~/.config/bammetrics/config.toml, -R may go
bammetrics chain fixed.bam reorder-sam insert-size -o r.bam

# The recipe, the verdict and the report in one command
bammetrics prep sample.bam -o out.bam qc --preset generic --fail-on fail --resume
```

### 32.1 qc and --resume in the short form

- `qc` is a step word (`chain BAM mark-duplicates insert-size qc -o r.bam`). In last
  place it keeps `--preset`, `--report`, `--fail-on`, its `--also` of verdict forms and
  the capture, RNA and sex options. Every other option after it is the chain's: `-o`,
  `-R`, `--threads`, `--sample`, `--resume`, `--dry-run`.
- `--resume` after the last step is the chain's. `--keep` after a transform's verb is
  that step's.
- `-o r.bam` on a chain that ends in `qc` adds `write r.bam --index` in front of the qc
  step (decision D200). The prefix is `r`: `r.insert_size_metrics`, `r.manifest.json`,
  `r.verdict.json`, `r.report.html`.
- `qc --index` in last place exits 2 and names the form: `... write OUT.bam --index qc`.
  A step after `qc` exits 2 (`step 2 (qc) is the last step: no step may follow it`).
- `chain - STEP ... --resume` exits 2: a pipe has no fingerprint.
- The reference from the config file or `BAMMETRICS_REFERENCE` is the chain's reference
  in the manifest. It is fingerprinted before the first read, as `-R` is, so another
  reference in the config file reruns every step.

### 32.2 A resumed run starts at the BAM it wrote

`-o r.bam` writes the last stream. A resumed run whose `write` step is skipped (its
record matches and `r.bam` and `r.bai` hold their fingerprints) starts its read at
`r.bam`, as it would at a kept stream (decision D201):

```
$ bammetrics chain fixed.bam reorder-sam insert-size --minimum-pct 0.2 -R hs37d5.fa -o r.bam --resume --dry-run
resume     from r.manifest.json
resume     step 1 reorder-sam: skipped; the BAM a write step wrote of its stream, r.bam, starts the read
resume     step 2 insert-size: runs: its options changed
resume     step 3 write: skipped: its record and files match
```

- The BAM is the one `--keep` would write at that position: the same writer on the same
  records. A `write` step counts when it reads the stream of the transform in question,
  writes a BAM (no SAM text, no standard output) and is skipped in this run.
- `r.bam` changed by another hand has another fingerprint. The `write` step then runs
  again and the transform in front runs to feed it, from the source.
- A kept stream comes first when both are there. `--keep` on the last transform beside
  `-o r.bam` writes the records twice and is not needed for `--resume`.
- A resumed run opens no file of a skipped step for writing. `r.bam` is read, never
  written, in the run above.

### 32.3 Recipes with qc and --resume

- One step may follow a recipe on the command line: `qc` (decision D199). `bammetrics
  prep in.bam -o out.bam qc --preset generic` runs the seven steps of `prep`, writes
  `out.bam` and `out.bai`, then the verdict and the page over the four metrics files the
  recipe's taps wrote. The verdict equals the verdict of `bammetrics qc out.bam`.
- Any other step beside a recipe exits 2, as before. A recipe file may end in a `qc`
  step; a second `qc` after it exits 2 (`recipe NAME ends in a qc step`).
- `qc.option=value` sets an option of a recipe's own qc step (`qc.preset=wgs`).
  `STEP.keep=true` keeps the stream of a recipe's transform.
- `--resume` works on a recipe as on its chain. An override is an option of its step, so
  `bammetrics prep in.bam -o out.bam --resume insert-size.minimum-pct=0.2` after a first
  run skips seven of the eight steps and reads `out.bam` for `insert-size`. An override
  of a transform (`mark-duplicates.remove-duplicates=true`) reruns that step and every
  step behind it.
- A user recipe and the config file are inputs of the chain (decision D204). A step that
  would write either exits 2 before the first read.
- Python: `bm.recipe("prep", bam, output="out.bam", qc=True, resume=True)`, or
  `qc=bm.step("qc", preset="generic", fail_on="fail")`. `Run.exit_code` holds 10 or 11
  from the qc step; no exception is raised for them.

### 32.4 Materialized steps in the short form

- A bare word after `merge-sam-files`, `gather-bam-files` or `jumping-library` that names
  no verb is a further input of that step: `chain a.bam merge-sam-files b.bam insert-size
  -o m.bam`. `--input b.bam` there is the step's too, in the middle and in last place
  (decision D205). A file named like a verb is written `./NAME`.
- `-o m.bam` on a chain whose only change is a materialized tool that writes a BAM
  (`chain a.bam merge-sam-files b.bam -o m.bam`) adds the `write`: the merge changes the
  records (decision D203).
- `--keep` on `merge-sam-files`, `gather-bam-files` and `position-based-downsample-sam`
  writes `PREFIX.merged.bam`, `PREFIX.gathered.bam` or `PREFIX.downsampled.bam`, and
  `--resume` starts there. Such a step has its stream as its id in the manifest and is
  upstream of the steps behind it, as a transform is. A further input that changed
  reruns it (`runs: b.bam changed`) and every step behind it.
- A materialized tool that writes no BAM (`wgs`, `library-complexity`) has a record of
  its own: unchanged, it is skipped; with a changed option it runs alone, on the kept
  stream or the written BAM when one starts the read.

### 32.5 The qc step and tools with no collector

`qc --preset wgs` needs `CollectWgsMetrics`, which reads a file. The qc step picks the
first of these (decision D202):

| The chain | What `CollectWgsMetrics` reads | The plan |
|---|---|---|
| a `wgs` step on the qc step's stream | nothing more: the verdict reads that step's file | `step 2  wgs  materialize on marked` |
| no transform in front, the source a file | the source, after the chain | `own read, after the chain` |
| `--keep` on the last transform | `PREFIX.<stream>.bam`, after the chain | `own read, after the chain` |
| a `write` step of the stream to a BAM (`-o out.bam`) | that BAM, after the chain | `own read, after the chain` |
| none of these | a copy of the stream in the temporary directory, in the chain's read | `materialize on marked, for qc` |

All five give the metrics `wgs` gives for the BAM the chain wrote, and one verdict
(`the_qc_step_measures_wgs_behind_a_transform`). On the manifest branch alone the last
row was `WgsMetrics not measured`.

### 32.6 What equals what

`crates/cli/tests/chain_combined.rs`, each against the long form or a fresh run:

| Test | What it holds |
|---|---|
| `the_owners_example_resumes` | `chain BAM reorder-sam insert-size -o r.bam --resume` with the config reference: the second run reads nothing; a changed `insert-size` option reruns that step alone from `r.bam` and leaves `r.bam` untouched; a replaced `r.bam` starts no read; the files equal the long form |
| `a_recipe_resumes_after_an_option_change` | `prep --resume`, then with `insert-size.minimum-pct=0.2`: 7 of 8 steps skipped, only that file and the manifest rewritten, the files of a fresh run; an override of a transform reruns every step |
| `qc_is_the_last_step_of_the_prep_recipe` | `prep BAM -o p.bam qc` equals the nine steps by hand and the verdict of `qc p.bam`; another step and a second `qc` exit 2 |
| `a_user_recipe_ends_in_qc` | A recipe file with a `qc` step, `qc.date=...`, the write in front of qc; a step that writes the recipe or the config file exits 2 |
| `a_merge_in_the_short_form_keeps_and_resumes` | `merge-sam-files FILE --keep` in the short form, the resumed read at `m.merged.bam`, the long form |
| `the_qc_step_measures_wgs_behind_a_transform` | The ways of 32.5 against `wgs` on the written BAM |
| `a_materialized_step_resumes` | `library-complexity` skipped, `wgs` rerun alone |
| `the_short_form_refuses_what_the_long_form_refuses` | The usage errors of 32.1, nothing written |

Golden records: the 9 of `CHAIN_COMBINED_CASES` and the 4 of `COMBINED_RESUME_CASES`
(`chain_short_form_qc`, `recipe_prep_qc`, `chain_short_form_merge`, the two plans
`chain_plan_qc_wgs_*`, four errors, `chain_short_form_resume_skips`,
`chain_short_form_resume_option_change`, its plan, `recipe_prep_resume_option_change`).
`tests/python/test_chain.py` holds the same cases through `chain`, `step` and `recipe`.

## 33. W16: a VCF source, the variant stream

A chain whose source is a VCF carries the source's record lines from step to step
(`crates/core/src/stage/variant.rs`). The source is told by its text: a regular file,
plain, gzip or BGZF, that starts with `#`. htsjdk reads the `##` lines in any order, so
`##fileformat` need not be first; a text that names no version fails when its header is
read, with the codec's message. A line ends at `\n`, `\r` or `\r\n`, as htsjdk's line
reader ends it. Standard input stays a stream of alignments.

```bash
# Copy a VCF as htsjdk writes it, with its .tbi: Picard VcfFormatConverter's output
bammetrics chain calls.vcf.gz -o out/calls.vcf.gz

# Two taps and the copy, one read
bammetrics chain calls.vcf.gz variant-calling --dbsnp dbsnp.vcf.gz \
    genotype-concordance --truth-vcf truth.vcf.gz --truth-sample NA12878 \
    -o out/calls.vcf.gz
```

### 33.1 Steps of a variant stream

| Step | What it does |
|---|---|
| `write [PATH]` | The lines of its position as a VCF (`crates/core/src/vcf_writer.rs`). A name ending in `.gz`, `.gzip`, `.bgz` or `.bgzf`, in any case, is BGZF; a name ending in lower-case `.vcf` is plain text; any other name (`x.bam`, `x.VCF`, `x.vcf.gz.tbi`) exits 2 before the first read, as htsjdk cannot type it; `.bcf` exits 2. `--index` writes the `.tbi` of a BGZF VCF; `--md5` writes `PATH.md5`; `--compression-level` 0 to 9 (default 5). Without a path, `PREFIX.<stream>.vcf.gz` |
| `variant-calling` | `CollectVariantCallingMetrics` as a tap: `PREFIX.variant_calling_detail_metrics` and `PREFIX.variant_calling_summary_metrics`. The tool alone reads `INPUT` whole and Picard needs it indexed; the tap reads the stream in its order and needs no index |
| `genotype-concordance` | `GenotypeConcordance` as a tap: the stream is the call VCF, `--truth-vcf` is read in step with it. `--call-vcf` on the step exits 2. `INTERVALS` still needs both files indexed, as the tool checks |

`-o r.vcf.gz` (also `.vcf.bgz` and `.vcf`, any case) on a chain whose source is a VCF
and that has no `write` step adds `write r.vcf.gz --index` as the last step, `write
r.vcf` for a plain name, and the prefix of every metrics file is `r`. A chain with no
step and `-o r.vcf.gz` is the copy above. On a VCF source any name htsjdk writes as
BGZF (`.gz`, `.gzip`, `.bgz`, `.bgzf`, any case) is a VCF too: `-o r.gz` writes `r.gz`
and `r.gz.tbi` with the prefix `r` (decision D244).

Refused before the first read, exit 2:

| Command line | Error |
|---|---|
| `chain calls.vcf insert-size` | `step 1 (insert-size): CollectInsertSizeMetrics reads alignments, and this chain's stream holds VCF records (the source is a VCF): it cannot be a step of this chain` |
| `chain in.bam variant-calling --dbsnp d.vcf` | `step 1 (variant-calling): CollectVariantCallingMetrics reads VCF records, and this chain's stream holds alignments (the source is a BAM, CRAM or SAM file): it cannot be a step of this chain` |
| `chain calls.vcf vcf-to-interval-list` | `step 1 (vcf-to-interval-list): VcfToIntervalList reads a VCF file of its own: it cannot be a step; ...` |
| `chain in.bam mark-duplicates -o r.vcf.gz` | `-o r.vcf.gz names a VCF, and this chain's source is no VCF: its stream holds alignments` |
| `chain calls.vcf -o r.bam` | `-o r.bam names a BAM, SAM or CRAM, and this chain's source is a VCF: its stream holds VCF records` |
| `chain calls.vcf write r.vcf --index` | the index of a plain VCF is tribble's `.idx`, which this writer does not write (decision D238) |
| `chain calls.vcf write -` | a VCF is written to a file |
| `chain calls.vcf ... qc` | qc reads the metrics of the taps of a BAM |
| `chain notes.txt -o r.vcf` | `-o r.vcf names a VCF, and the source notes.txt is no VCF: its text does not start with #, as a VCF header does`; an empty file `it is empty`, a BCF `it is a BCF, and bammetrics reads VCF text only` |

A record the codec or the encoder refuses (an `INFO` key or `FILTER` name with no header
line, a negative `QUAL`, a sample with more values than keys when the samples are out
of sorted order, a `GL` with some values `.`, a field a bare `\r` cuts short) fails the
`write` step with exit 1 and leaves no file, as Picard fails
(`tests/concordance/test_vcf_writer.py::test_error_parity`). The `.tbi` needs the
records sorted by contig and start: `Sequence c1 added out of order` or `Features added
out of order` fail the step, and nothing is renamed into place. A record that starts or
ends past 536,903,680 (the 32,770 16 kb windows of a `.tbi`) fails the index, as htsjdk
fails there with an `ArrayIndexOutOfBoundsException`; the plain copy writes it.

The reader holds a header of at most 64 MiB (`MAX_HEADER_BYTES`) and a line of at most
16 MiB (`MAX_LINE_BYTES`); htsjdk reads any length. A longer header or line fails the
read with exit 1 and leaves no file (KNOWN_DIFFERENCES.md).

### 33.2 What the writer writes

The writer is htsjdk 5.0.0's `VCFWriter` for the record of `VCFCodec`:

- The header: `##fileformat=VCFv4.2` first, whatever version the input names, then every
  kept line once in sorted order (contig lines in their own order), then the column
  line. A standard `INFO` or `FORMAT` line whose count or type differs from htsjdk's is
  replaced by htsjdk's (`FT` with `Number=1` becomes `Number=.`, `Genotype-level
  filter`: the rewrite the inventory found in every case). An `INFO` or `FORMAT` line
  loses keys other than `ID`, `Number`, `Type`, `Description`, `Source` and `Version`;
  one without `Description` gets `Not provided in original VCF header`; a line with no
  `=` is dropped.
- Each record: `POS` as an `int`, `REF` and the alternates upper case (symbolic ones as
  given), `QUAL` as `%.2f` with `.00` cut (missing for `.`, `-1` and `-10`), the
  `FILTER` names sorted once, the `INFO` keys sorted with a flag's value dropped, and
  the genotypes as read. When the header names its samples out of sorted order htsjdk
  decodes them, and so does the writer: the keys `GT` first then sorted, `GQ` at most
  99 and absent when it rounds to -1, `FT` `PASS` for an unfiltered genotype, a missing
  value one `.` whatever the key's `Number`, trailing missing values dropped. `GL`
  becomes `PL` (`10 x` the distance to the largest likelihood, rounded) until a `PL`
  value of any sample of the record is read; `GL` itself is not written. A `QUAL` of
  `-0` prints `0`.
- The text is ISO-8859-1, as `VCFWriter` writes it: a character up to U+00FF as its
  one byte, any other as `?`. The reader decodes UTF-8, so `é` in the input is the byte
  `E9` in the output, and a byte that is no UTF-8 is `?`.

The memory: the source is read in batches of at most 1024 lines and 4 MiB of text
(`stage::variant::BATCH_BYTES`), at most 8 batches queued per tap or writer; a writer
holds one line and its BGZF workers' blocks; a record whose samples are out of sorted
order is decoded into one array per record, with no allocation per sample. The header is
held once, shared by every step, its lines as written in one buffer: a header of `n`
lines of `b` bytes costs about `2b + 100n` at its peak. Measured: 1,000,000 records
under 64 MiB; 20,000 and 100,000 samples out of sorted order 61 and 133 MiB; 3,150,000
`##contig` lines (65 MB, the largest header read) 463 MiB; 1,000,000 `##contig` lines
199 MiB (`benches/memgate.toml`, `chain vcf *`). docs/07 section 5 sets 500 MiB.

### 33.3 What equals what

| Test | What it holds |
|---|---|
| `tests/concordance/test_vcf_writer.py::test_round_trip` | every VCF of `variant_vcfs.py` and 16 corner files (GL, missing values, GQ and QUAL, an empty sample cell, ISO-8859-1 text, bare `\r` line ends, `##fileformat` not first among them), plain and `.vcf.gz`, through `chain X -o OUT` and Picard `VcfFormatConverter`: header lines (date, command and inventoried volatile lines aside) and every record line equal; each `.tbi` answers every whole-contig and 200-base `tabix` query as Picard's does |
| `test_headers_match_the_inventory` | each case of `tests/concordance/vcf_headers.json` whose output header adds and removes no line: the chain's header of its input header equals the recorded output |
| `test_error_parity`, `test_unsorted_records_fail_the_index` | 13 refused inputs fail both tools with no file left; unsorted records fail the index in both and copy equal without it |
| `test_index_limits_fail_in_both`, `test_compressed_names_are_bgzf` | a record past the windows of a `.tbi`, by its start or its `END`, fails the index in both; `.gz`, `.bgz` and `.GZ` names are BGZF with a `.tbi` in both |
| `test_variant_calling_tap`, `test_genotype_concordance_tap` | the taps write the bodies the tools write alone |
| `crates/core/src/stage/variant.rs` tests | taps see the lines of their position, a transform's counters, a tap fails alone, a failing writer stops the chain and leaves no file |
| `crates/cli/tests/vcf_chain.rs` | `--resume` skips a VCF `write` and reruns a changed tap; the manifest records kinds `vcf`, `index` and `md5` |

Golden records: the 14 of `VCF_WRITER_CASES` on `tests/data/vcf_writer`
(`tests/data/make_vcf_writer_fixtures.py`).

### 33.4 Departures from docs/07

The decision rows D236 to D245 of docs/07 section 7. In short: the variant stream runs
on an executor of the same shape beside `exec.rs` (no barrier yet); its item is the
record line, which the writer needs to pass genotypes through as htsjdk does; the
plain VCF's `.idx` is refused; `-o r.vcf.gz` copies a VCF with no transform; the `.tbi`
is built from the finished file read back.

## 34. W16: the VCF tools as steps of a variant stream

`make-sites-only-vcf` is a transform and `split-vcfs` a split of a variant stream;
`make-vcf-sample-name-map` is no step.

```bash
bammetrics chain calls.vcf.gz make-sites-only-vcf -o sites.vcf.gz
bammetrics chain calls.vcf.gz make-sites-only-vcf --sample S1 split-vcfs --strict false -o out/r
bammetrics chain calls.vcf.gz split-vcfs --strict false --snp-output s.vcf.gz --indel-output i.vcf.gz variant-calling --dbsnp dbsnp.vcf.gz -o out/r
```

- `make-sites-only-vcf` edits each line as it passes (`stage::variant::VARIANT_TRANSFORMS`,
  `vcf::SampleSubset`) and leaves the stream `sites` (or `--as NAME`). The steps after it
  read the header without sample columns, or with those of `--sample`. It starts no
  segment. `-o`, `--create-index` and `--compression-level` on the step exit 2: the VCF is
  the `write` step's, or the chain's `-o r.vcf.gz`.
- `split-vcfs` is a sink of the stream at its position (`VARIANT_SPLITS`): the SNP and
  indel VCFs, with a `.tbi` beside each `.vcf.gz`, renamed into place as one set when the
  stream ends. Without `--snp-output` and `--indel-output` they are `PREFIX.snp.vcf.gz`
  and `PREFIX.indel.vcf.gz` under the step's `-o` or the chain's prefix (`PREFIX.STREAM`
  when a transform comes after the step). Its failure stops the chain, as a `write` step's does, so
  `--strict` (the default) on a record of another type exits 1 and leaves no file. Its
  note counts the records `--strict false` dropped.
- `--keep` on a step of a variant stream exits 2: a kept stream is a BAM.
- After a transform each step but `write` reads the stream as a reader of the transform's
  own file would read it: each line encoded as the writer encodes it, the header and the
  lines through ISO-8859-1 and back as UTF-8 (a character above U+00FF is `?`, a Latin-1
  letter the replacement character), the header at `VCFv4.2`. A `write` step right after
  a transform writes what the transform's own file holds. So `chain X make-sites-only-vcf
  split-vcfs` equals `MakeSitesOnlyVcf` then `SplitVcfs` for any header text, a `VCFv4.3`
  input included (docs/07 D252).
- `make-sites-only-vcf` exits 2 under the chain's `--validation-stringency` other than
  `STRICT`: alone the tool writes keys the header lacks there, and the steps after it check
  every key (D254).
- On a BAM stream all three exit 2 naming the stream's kind; `make-vcf-sample-name-map`
  exits 2 on either kind: it reads the headers of several files.
- The manifest records the split's files as kinds `vcf` and `index`; `--resume` skips the
  step when its files and options are unchanged, and runs the transform in front of any
  step that runs.

### 34.1 What equals what

| Test | What it holds |
|---|---|
| `tests/concordance/test_split_vcfs.py` | every VCF of `variant_vcfs.py` and `types.vcf`, plain and gzip, into plain and `.vcf.gz` outputs: record and header lines equal to Picard's, `.tbi` queries equal, plain files byte for byte; `STRICT` fails both on the same type; every option; 6 error cases; `chain SRC make-sites-only-vcf split-vcfs` equals Picard `SplitVcfs` on `MakeSitesOnlyVcf`'s file |
| `tests/concordance/test_make_sites_only_vcf.py` | the same inputs with every `SAMPLE` shape (each, all, a missing name, `null`), haploid no-calls, lazy and eager genotypes, 6 error cases; `chain X make-sites-only-vcf -o r.vcf.gz` equals Picard |
| `tests/concordance/test_make_vcf_sample_name_map.py` | the map byte for byte for relative, absolute, `file://`, repeated and 40 inputs; the duplicate-sample warning; 10 error cases |
| `crates/cli/tests/vcf_tools_chain.rs` | chained equals staged for both steps, `--resume` over a split with the transform in front rerun, a strict split with no file left, the refusals |

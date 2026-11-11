# In-process chaining

A chain runs several steps on one BAM in one process. Records move from step to step in
memory. No shell pipe and no intermediate file is used. This page describes what waves W03
and W04 of docs/07_full_parity_program.md (section 2) built: the stage layer and its
executor (W03, 0.21.0), then the `chain` verb, its plan file, its plan text and the Python
`chain` and `step` (W04, 0.22.0).

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
| source | One BAM or CRAM file, or a pipe holding one | the positional argument, `-` for standard input |
| transform | A tool that changes the records later steps see. A barrier needs every record before one goes on | `mark-duplicates` (a barrier) |
| tap | A one-pass tool on the records of its position. It changes nothing | 13 verbs: `insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics`, `gc-bias`, `quality-score-distribution`, `mean-quality-by-cycle`, `base-distribution-by-cycle`, `duplicate-metrics`, `rna-seq`, `sequencing-artifact`, `coverage-gaps`, `sex-check` |
| write | A tap that writes the records of its position as a BAM | `write [PATH]` |

A tap or `write` does not change the stream, so the step after it sees the same records.
`in.bam + tap + mark-duplicates + tap + write` gives the first tap the source records and
the second tap and the writer the marked records.

A tool without a collector and without a stage (`wgs`, `targeted-pcr`,
`library-complexity`, `multiple-metrics`, every tool that reads no BAM) cannot be a step.
W14 spools the stream for those.

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

## 2. The chain verb

```
bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]
bammetrics chain --plan FILE [OPTIONS] [BAM]
```

The chain's options and the source come first. Every `+` starts a step. The first word of
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
| `-o`, `--output-prefix`, `O=` | The prefix of every step that names no file of its own. Default `./` plus the source name without its extension, and standard error says so. A directory (an existing one, or a path ending in `/`) gets the source name as the file part |
| `-R`, `--reference`, `R=` | Reference FASTA with a `.fai`. It decodes a CRAM source. It goes to every step that reads a reference and names none |
| `-@`, `--threads` | BGZF threads of each read and of `write` (default: cores, at most 8; at most 256 may be asked for) |
| `--tmp-dir`, `TMP_DIR=` | Where a transform spills and a piped source is copied |
| `--format LIST` | `picard` first. Each further format is written beside every metrics file |
| `--also LIST` | `json`, `jsonl`, `tsv`, `csv` beside every metrics file; `markdown` as `PREFIX.report.md` |
| `--sample NAME` | Sample name of the json and jsonl forms (default: the first `@RG SM`, else the source name) |
| `--plan FILE` | Read the source, the options and the steps from a plan file (section 3) |
| `--dry-run` | Print the plan and run nothing (section 4) |
| `-q`, `-v` | No notes and no summary; or each command line and each file written |

The source, `--output-prefix`, `--threads`, `--sample`, `--also`, `--plan`, `--dry-run`,
`-q` and `-v` belong to the chain. Given after the first `+` they are a usage error that
says so: `step 1 (insert-size): --also belongs to the chain: give it before the first +`.
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
- `write [PATH]` takes `--index` (`CREATE_INDEX`), `--md5` (`CREATE_MD5_FILE`) and
  `--compression-level N` (0 to 9, default 5). `-` as the path is standard output, so a
  pipe to another tool still works: `bammetrics chain in.bam + mark-duplicates + write -
  | samtools flagstat -`. The file is always a BAM. A name ending in `.sam`, `.cram`,
  `.fastq` or `.fq` (with `.gz` too) is refused with the wave that adds that kind: SAM is
  W05, CRAM is W15, FASTQ is W22.
- `--index` writes `<name>.bai` only for a coordinate-sorted stream. For any other sort
  order `write` writes no index and says so on standard error. The exit code stays 0.
- A file option with an empty value (`-o ''`, `--metrics-file ''`, `write ''`) is a usage
  error: `step 1 (gc-bias): -o names no file: its value is empty`.

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
json|jsonl|tsv|csv` adds the same siblings. `--also markdown` writes `PREFIX.report.md`
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
| 2 | A usage error: an unknown step or option, a chain option on a step, a value that is no number, a missing required option, an empty prefix or file name, a name `write` cannot write yet, two outputs in one file, an output that is an input, a second transform, a bad plan file (one that is not UTF-8 text included) | Before the first read. No file is written and no directory is made |
| 1 | An input that is not there or cannot be read (the source, `-R`, the plan file) | Before the first read |
| 1 | An output that cannot be written (its directory cannot be made, the file is read only) | Before the first read, as the Picard face exits for its `OUTPUT` |
| 1 | A step failed | During the run. The other steps still write their files (section 9) |

```
$ bammetrics chain in.bam + quality-yield + insert-size --minmum-pct 0.1
bammetrics chain: step 2 (insert-size): unknown option --minmum-pct
  did you mean --minimum-pct (Picard MINIMUM_PCT)?
  every option: bammetrics help insert-size
$ bammetrics chain in.bam + inser-size
bammetrics chain: step 1: unknown step "inser-size"
  did you mean insert-size?
  a step is a tool verb (bammetrics tools lists them) or write
$ bammetrics chain in.bam + write out/x.sam
bammetrics chain: step 1 (write): out/x.sam: bammetrics writes BAM only: SAM output is wave W05 of docs/07_full_parity_program.md
```

What is checked before the first read: the step name, every option name, every value of
an int or float option, required options, files a step reads, `STOP_AFTER`, every option
of `mark-duplicates` (`stage::check_stage_options`), the names and the writability of
every output. Section 2.6 adds five checks.

### 2.6 More that is refused before the first read

The 0.22.0 release review found these. Tests: `crates/cli/tests/chain_safety.rs` (7
tests), `tests/python/test_chain_safety.py` and 8 golden records.

| What | The message says | Exit |
|---|---|---|
| A text value only the tool can judge on a tap (`--level NOPE`, `--expected-sex NOPE`) | The tool's own words under the step's label: `step 2 (insert-size): METRIC_ACCUMULATION_LEVEL "NOPE": expected ALL_READS, SAMPLE, LIBRARY or READ_GROUP` | 2 |
| `--threads` over 256, or `threads` over 256 in a plan file | `100000000 is more than the maximum of 256 threads` | 2 |
| An output through a directory that does not exist with a `..` behind it (`missing/../in.bam`) | `missing does not exist, so the .. behind it names no directory` | 2 |
| A plan file over 1 MiB (1,048,576 bytes) | `not a plan file: it is larger than 1 MiB` | 2 |
| An output that is a symbolic link, pipe, socket or device another user owns | `owned by another user; not written in place`, the words and the exit code of the tool verb | 1 |

Values. `bammetrics_core::tap_options::check` holds one row per one-pass tool (13). Each
row calls the option parser that the tool's own collector calls, with no header and no
file opened. A unit test fails when a one-pass tool has no row. In 0.22.0 the tap failed
alone when it was built, after the transform in front of it had read the whole source,
and the chain exited 1. A value the tool judges against the header of its stream, or
against a file it reads, is still found when the tap is built, and that tap fails alone
with exit 1: a sort order it requires (`--no-assume-sorted` on a queryname stream), an
interval list with a contig the header lacks. The reason: the header of a piped source
can be read once, and the header behind a transform exists only after the transform ran.

Threads. `pass::check_threads` is the one check and `pass::MAX_THREADS` is 256. Every
place that parses a thread count asks it: the tool verbs, `chain`, the plan file, `qc`,
`plan`, `collect`, `threads=` on the Picard face (exit 1 there, as every error of that
face) and the Python `chain`, `qc` and `collect`. In 0.22.0 and before,
`--threads 100000000` started threads for about two minutes and then panicked inside
`noodles-bgzf` with exit 101, on every verb.

Paths. A run makes the directories that are missing above its outputs. `missing/../in.bam`
names no file before the run and names `in.bam` once `missing` is made, so neither the
one-file-twice check nor the input check saw the file. In 0.22.0 `write missing/../in.bam`
failed at write time with exit 1 and left `missing/` behind, and a tap with that `-o`
replaced the source BAM with its metrics file. The tool verbs had the same hole
(`insert-size in.bam -o missing/../in.bam`) and refuse the path with the same words. A
`..` behind a directory that exists resolves and is checked as any spelling. Two
spellings of one new file under a directory that does not exist yet (`link/new/a.txt`
and `new/a.txt`, `link` leading to the working directory) are one file: the deepest
directory that exists is resolved (`mark_duplicates::same_file_key`).

Directories. `cleanup::MadeDirs` records each directory the run makes above an output. A
run that ends with exit 1, before the first read or after it, removes those that are
still empty, deepest first. A directory that was there before, or that holds a file,
stays. A signal removes them under the same rule. A run that finishes keeps every one.

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

- One `[chain]` table. Its keys: `source`, `output-prefix`, `reference`, `threads`,
  `tmp-dir`, `format`, `also`, `sample`, `quiet`, `verbose`. `format` and `also` take a
  string (`"picard,json"`) or a list of strings. An empty `output-prefix` is an error.
- One `[[step]]` table per step, in chain order. `verb` is required: a tool verb, a
  Picard class, or `write`. Every other key is the long name of an option of that verb
  (`minimum-pct`), or its Picard name (`MINIMUM_PCT`). A boolean option takes `true` or
  `false`, an int option an integer, a float option a number, a text or path option a
  string. One number also stands for the text of a text option (`program-record-id = 7`).
  A repeatable option takes a list; a list of a text option holds strings only. The keys
  of `write` are `output`, `index`, `md5` and `compression-level`.
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
temp disk  none expected: about 16 records, and mark-duplicates keeps 4,000,000 in memory before it spills (MAX_RECORDS_IN_RAM)
nothing ran: this is a plan
```

The plan lists:

- The segments, and what each reads. A transform that reads its input twice
  (`mark-duplicates`) replays a file source: the file is read a second time and nothing
  is spooled. A piped source is copied to the temporary directory once (the spool) and
  the copy is read twice.
- Each step with its role, the stream it reads and every file it writes, the further
  formats included. A tool whose `OUTPUT` is a prefix is listed as `PREFIX.*`.
- The estimated temporary disk. The record count comes from the file size: 69 compressed
  bytes per record for a BAM, 21 for a CRAM (the HG002 exome BAM and the NA12878 genome
  CRAM of docs/benchmarks.md). Under `MAX_RECORDS_IN_RAM` records nothing spills. Above
  it the estimate is 59 bytes per record, the larger of the two peaks measured there
  (46 and 59). It gives 8.2 GiB for the exome BAM (measured 6.5 GiB) and 41.3 GiB for the
  genome CRAM (measured 42.1 GiB). A pipe has no size to read, so the plan says what
  spills and gives no number.

The plan reads the header of a file source, so a file that is no BAM exits 1 there. It
reads no record and nothing of a pipe. The planner of W05 adds the free-space refusal and
the sort order each step requires.

## 5. The Rust API

### 5.1 The stage layer

`bammetrics_core::stage` (`crates/core/src/stage.rs`, `stage/exec.rs`,
`stage/sinks.rs`).

| Item | Role |
|---|---|
| `Source` | Opens the stream. `read(threads, bytes, consume)` hands `consume` the `Frame` and the `Records`. `replayable()` is true for a regular file. `spool(temp)` copies a pipe. `file()` names the file a lowered chain reads |
| `FileSource` | A BAM or CRAM by its magic bytes, read through `pass::with_counted_stream` |
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
| `WriteOptions` | `path`, `compression_level`, `create_index`, `create_md5`. `check()` refuses a name of a later wave. `files()` lists the BAM, the `.bai` and the `.md5` |
| `later_kind(path)` | The kind and wave of a name `write` cannot make yet: `("SAM", "W05")` |
| `Outputs` | The files of a chain's steps in a hash map, keyed by device and inode or by the resolved path. `claim(label, path)` is a usage error when a step already writes the file. `owner(path)` |
| `check_distinct_outputs(outputs)`, `writer_of(outputs, path)` | The same check over a slice |
| `stage_of(tool)`, `label_of(tool)`, `check_stage_options(tool, inv)` | The stage function of a registry tool, the default name of its stream (`marked`), and its options checked before a run. One tool: `MarkDuplicates` |

`crates/core/src/bam_writer.rs` holds `BamWriter` (the header byte for byte, BGZF on
worker threads, htsjdk's order check, an `AtomicFile` under a temporary name) and
`finish` (the `.bai` and `.md5`, then the renames). `MarkDuplicates` and `write` share
it.

### 5.2 The verb as a library

`bammetrics_cli::chain` (`crates/cli/src/chain.rs` and `chain/`). A caller that is not the
command line builds a `Spec`, calls `prepare`, then `execute`. The Python module spells the
argv of its keywords and calls `parse`, so one parser builds the `Spec` of a command line
and of a Python call (section 5.3).

| Item | Role |
|---|---|
| `Spec`, `StepSpec` | A chain as asked for: `source`, `reference`, `prefix`, `threads`, `tmp_dir`, `formats`, `also`, `sample`, `quiet`, `verbose`, `steps`. A step is its verb and its arguments as text (`StepSpec::new("insert-size", &["--minimum-pct=0.1"])`) |
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
sample=, plan=, dry_run=)` and `bammetrics.step(verb, **options)` run this verb in the
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
| `StepTable::keyword(name)` | One keyword, or the usage error with the closest name. `input`, and `output` of `mark_duplicates`, are refused with the reason |
| `Keyword` | `name`, `long` (the flag it spells), `ty`, `many` (`One`, `Repeated`: one occurrence per item, `Joined`: one comma list) |
| `step_names()` | The 15 names `step` takes, snake case, in `bammetrics tools` order, `write` last |
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
the `Run`. `dry_run=True` returns `plan_text`. `tests/python/test_chain.py` (75 tests)
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

Limit: a transform that replays must sit directly behind a file source. A second
transform is refused by `prepare`, before the first read: `step 3 (mark-duplicates): it
reads its input twice and cannot follow step 1 (mark-duplicates): only a file source can
be read again in this version`. It is a usage error (exit 2). No step runs and no file is
written. `run_chain` refuses the same plan for a caller of the stage layer. W05 adds the
spool and the planner that lift this.

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
| `write` picks BAM, SAM, CRAM or FASTQ by extension | BAM only. A `.sam`, `.cram`, `.fastq` or `.fq` name is a usage error naming W05, W15 or W22 | The SAM text writer is W05, the CRAM writer W15, FASTQ output W22 |
| Two taps of one tool at one position is a usage error | It is one when both are named by the prefix, since they write one file. With `-o` on one of them the chain runs | A step's own `-o` is allowed, and two option sets of one tool on one stream is a use (`quality-yield` with and without original qualities) |
| The `# INPUT=` line names the chain source and the step list | `INPUT=<source>` plus ` CHAIN=<transforms in front>`, as W03 wrote it. A tap on the source writes the tool's own line | A chain without a transform must write the bytes `qc` writes (section 10.2). The steps behind a tap do not change its records, so the line names only the transforms in front |
| A transform's kept product is `PREFIX.<label>.bam` | The same, written by `write` without a path. The transform's own metrics file is `PREFIX.<label>.duplication_metrics` | docs/07 names no file for a transform's metrics. The label says which step wrote it |
| `plan` prints segments, replay or spool per barrier, temp disk, and what `--resume` skips | Segments, replay or spool, the files of each step, the temp disk estimate. No `--resume` line | Resume per step is W14 (manifest 2) |
| The planner refuses when free space is below the estimate | The plan prints the estimate and refuses nothing | The planner and its `statvfs` call are W05 |
| `--also` on `chain` as on `qc` | `json`, `jsonl`, `tsv`, `csv`, `markdown`. `junit` and `multiqc` are a usage error | They print a verdict. `qc` as the last step of a chain is W14 |
| `--format a,b` | `picard` must be first, on the chain and on a step | A chain names Picard files under a prefix, as the multi-file verbs do. A first format that is not `picard` would put JSON under a Picard name |
| `parse_tables` moves to core | It moved, and it reads one-line arrays of scalars as well | A repeatable option of a step (`level`, `comment`) needs a list. A preset still takes no list |
| `chain --plan plan.toml` reads the same steps | It reads the source and the chain's options too, and the command line may replace them | One file then describes one run, and one file can run on many sources (`--plan FILE other.bam -o other`) |
| `plan` gains the chain form | `plan` with a `+` step, with `--plan FILE`, or with the word `chain` first | `plan` is a clap verb with the arguments of `qc`. A `qc` command line holds no `+` and no `--plan`, so the two cannot be mistaken |
| Text-valued options are checked before the first read | They are, since the 0.22.0 review: `tap_options::check` runs each one-pass tool's option parser without a header (section 2.6). A value judged against the header or a file still fails its tap when it is built | A collector parses its options and then builds from the header of its position. The parser alone needs no header and opens no file, so no reference or interval list is read twice |

W04, Python (docs/07 section 2.2, the `crates/py` paragraph):

| docs/07 | As built | Reason |
|---|---|---|
| `crates/py/src/verbs.rs` holds `bm.chain` | `crates/py/src/chain.rs` | `verbs.rs` has 1,138 lines. The chain classes and functions are one file of 610 |
| `bm.chain(..., resume=)` | No `resume` keyword | Resume per step is W14 (manifest 2). A keyword that does nothing would be a lie in the signature |
| `Run.stages` holds the transform counters | One `Stage` per step: `status`, `error`, `files`, `reads`, `leaves`, `segment`, and `records_in`, `records_out`, `changed` for a transform (`None` for a tap and `write`) | A failed run needs the status of every step, and `Products::counters` exists for a stage only. A chain of taps lowers to the read of `qc`, which counts no records |
| `bm.step` validates against the generated registry signature | It does, through `registry::python_keywords`. `input` is refused, `output` is the step's own `-o` (a prefix for a tool that writes several files), `mark_duplicates` takes no `output`, and `format` and `as_` are added | A step has no input of its own, a transform writes no BAM in a chain, and `--format` and `--as` are step options with no row in the registry |
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
| W05 | The sort barrier (`Next::Records`), the planner (`Frame` per position, `requires`, the free-space refusal), the spool for a replay behind another transform, `--sort-memory`, SAM output of `write`, `write --index` through the BAI writer |
| W09 on | Each BAM-writing tool registers its stage in `stage_of` and its stream name in `label_of` |
| W14 | `qc` as the last step of a chain, manifest 2, resume per step (`resume=` in Python), kept products, the spool for tools that read their input their own way, a progress line |
| W15, W22 | `write` to CRAM, FASTQ |
| W46 | Variant streams (`Kind`), pileup and capture sinks |

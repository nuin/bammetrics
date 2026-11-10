# In-process chaining

A chain runs several steps on one BAM in one process. Records move from step to step in
memory. No shell pipe and no intermediate file is used. This page describes wave W03 of
docs/07_full_parity_program.md (section 2): the stage layer, the executor, `MarkDuplicates` as
the first stage, the tap adapter and the `write` step. The public `chain` verb is wave
W04 (section 11).

```bash
# MarkDuplicates, four collectors on the marked records and the marked BAM, one process
bammetrics __chain I=in.bam \
    STEP=MarkDuplicates M=out.duplication_metrics \
    STEP=CollectAlignmentSummaryMetrics O=out.alignment_summary_metrics \
    STEP=CollectInsertSizeMetrics O=out.insert_size_metrics \
    STEP=CollectQualityYieldMetrics O=out.quality_yield_metrics \
    STEP=CollectDuplicateMetrics M=out.duplicate_metrics \
    STEP=write O=out.marked.bam CREATE_INDEX=true
```

The four metrics files equal, byte for byte after the `StringHeader` lines, the files the
four tools write when they read `out.marked.bam`. The BAM holds the records the
`MarkDuplicates` tool writes.

## 1. Steps

A chain is a source and a list of steps. Each step sees the records the step before it
emits.

| Step | What it is | In W03 |
|---|---|---|
| source | One BAM or CRAM file, or a pipe holding one | `I=` |
| tap | A one-pass tool on the records of its position. It changes nothing | every registry tool with a collector (13 tools) |
| stage | A tool that changes the stream. A transform edits batches in place. A barrier needs every record before one goes on | `MarkDuplicates` (barrier) |
| sink | A writer of the records of its position | `write` (BAM) |

A tap or a sink does not change the stream, so the step after it sees the same records.
`MarkDuplicates I=in.bam + tap + tap + write` therefore gives every tap and the writer
the marked records.

A tool without a collector (`CollectWgsMetrics`, `CollectTargetedPcrMetrics`,
`EstimateLibraryComplexity`, `CollectMultipleMetrics`) cannot be a step yet. W14 spools
the stream for those.

## 2. The test entry

`bammetrics __chain` is hidden. `tools`, `help` and the completions do not list it. It
exists so the tests and the measurements can run a chain from the command line. W04
replaces it.

```
bammetrics __chain I=in.bam [R=ref.fa] [threads=N] [TMP_DIR=dir] \
    STEP=<Picard class or write> [KEY=VALUE ...] \
    STEP=...
```

- Options before the first `STEP=` belong to the chain: `INPUT` (`I`),
  `REFERENCE_SEQUENCE` (`R`, for a CRAM source), `threads` (default: the core
  count, at most 8), `TMP_DIR` (where a barrier spills). `LOWER=false` forces the stage executor
  for a chain of taps (section 9).
- Every `STEP=` starts a step. The options after it are that step's, in the tool's own
  Picard spelling. The chain's `INPUT` is given to every tool step. The chain's `R=` is
  given to each tool that reads the reference, unless the step names its own.
- A tap needs its file option (`O=`, or `M=` for `CollectDuplicateMetrics`).
- `MarkDuplicates` takes every option of the tool except `OUTPUT`, `COMPRESSION_LEVEL`,
  `CREATE_INDEX` and `CREATE_MD5_FILE`. Those belong to `write`. `METRICS_FILE` is
  optional. Its record counts (in, out, changed) go to standard error.
- `write` takes `OUTPUT` (`O`), `COMPRESSION_LEVEL` (0 to 9, default 5),
  `CREATE_INDEX` and `CREATE_MD5_FILE`. The file is always BAM. A name ending in
  `.sam` or `.cram` is refused. A path that is the source or the reference is refused.

Exit codes: 0, 1 when a step failed, 2 for a usage error found before the run. A failed
step prints `bammetrics __chain: <step>: <error>`. The other steps still write their
files.

## 3. The Rust API

`bammetrics_core::stage` (`crates/core/src/stage.rs`, `stage/exec.rs`,
`stage/sinks.rs`).

| Item | Role |
|---|---|
| `Source` | Opens the stream. `read(threads, bytes, consume)` hands `consume` the `Frame` and the `Records`. `replayable()` is true for a regular file. `spool(temp)` copies a pipe. `file()` names the file a lowered chain reads |
| `FileSource` | A BAM or CRAM by its magic bytes, read through `pass::with_counted_stream` |
| `Records` | `next_batch(&mut Vec<RawRecord>)`: at most 4096 records per batch |
| `Frame` | The stream at one position: `kind`, `order`, `header`, `cram`, `label`, `reads` (the files a writer must not replace) |
| `Transform` | `apply(batch)` in place, `finish(out)` for held records, `header(up)`, `counters()`, `products()` |
| `Barrier` | `push(batch)`, then `seal()` gives `Sealed { header, next }` |
| `Next` | `Records(...)` when the barrier holds the records (a sort). `Replay(transform)` when the source is read again through a transform (`MarkDuplicates`) |
| `Sink` | `accept_batch(&[RawRecord])`, `done()`, `finish()`. One thread per sink |
| `CollectorSink` | The tap adapter: any `pass::Collector` as a `Sink` |
| `Emit` | Where a producer puts each record. `BamWriter` and `Vec<RawRecord>` implement it |
| `Step` | `tap`, `tool_tap(tool, inv)`, `tool_stage(tool, inv)`, `stage`, `sink`, `write(WriteOptions, threads)` |
| `ChainPlan` | `source`, `steps`, `threads`, `tmp_dir`, `lower` |
| `run_chain(plan, progress)` | One `Result<Products>` per step, in step order |
| `Products` | `output` (metrics the caller writes), `files` (a writer's renamed files), `counters` (a stage's records in, out, changed) |
| `stage_of(tool)` | The stage function of a registry tool. One row: `MarkDuplicates` |

`crates/core/src/bam_writer.rs` holds `BamWriter` (the header byte for byte, BGZF on
worker threads, htsjdk's order check, an `AtomicFile` under a temporary name) and
`finish` (the `.bai` and `.md5`, then the renames). `MarkDuplicates` and `write` share
it.

## 4. Segments and threads

The executor cuts the steps at barriers. Each part is a segment.

- One segment is one read of its input on the calling thread. That thread also runs the
  transforms and feeds the barrier.
- Every tap and every sink has its own thread. It takes batches of at most 4096 records
  through a queue of 8 batches. A batch is shared by the sinks of one position
  (`Arc<Vec<RawRecord>>`). A position in front of a transform copies its batch once,
  and only when a sink listens there.
- BGZF blocks are inflated on `threads` workers. `write` compresses on `threads`
  workers.
- Segments run in turn. Peak memory is the largest segment, not the sum.

`in.bam + tap A + MarkDuplicates + tap B + write` is two segments. Segment 1 reads
`in.bam`, feeds tap A and the first read of `MarkDuplicates`. The barrier is sealed.
Segment 2 reads `in.bam` again through the second read of `MarkDuplicates` and feeds
tap B and the writer.

A chain of taps alone is not run by this executor. It is lowered to
`pass::run_isolated`, the read `qc` uses, so `qc` and every golden record keep their
bytes and their speed. `qc` itself is unchanged in W03.

## 5. MarkDuplicates as a barrier

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
to `TMP_DIR` once, as the tool does. Read ends spill to `TMP_DIR` as in the tool
(`MAX_RECORDS_IN_RAM`).

A step after the barrier is built from the sealed header, so it sees the new `@PG` line.

Limit of W03: a stage that replays must sit directly behind a file source. A replaying
stage behind another stage is refused before any record is read (`reads its input
twice and cannot follow ...`). W05 adds the spool and the planner that lift this.

## 6. Options in the middle of a chain

`STOP_AFTER` and `ASSUME_SORTED` belong to the tap that carries them.

- `STOP_AFTER=N` on a tap: the tap takes the first N records its position emits. After
  `MarkDuplicates REMOVE_DUPLICATES=true` those are the first N records that were kept,
  the same N records the tool reads from the written BAM.
- It never shortens the stream. The steps after the tap get every record.
- The read stops early only in a segment with no stage, when every tap of it is done.
  A lowered chain stops at the largest `STOP_AFTER` when every tap has one.
- `ASSUME_SORTED=false` on a tap checks the header of its position. After a barrier
  that is the sealed header.

`crates/cli/tests/chain.rs::stop_after_in_the_middle_of_a_chain_belongs_to_its_tap` pins
this against the staged runs. `MarkDuplicates` has no `STOP_AFTER`.

## 7. Failure

| What fails | Result |
|---|---|
| A tap (cannot be built, fails on a record, fails to finish, panics) | That tap has its error. Every other step finishes |
| A stage on a record, or a stage or `write` that cannot be built | Its segment stops. The failing step has its error. Every other step of the segment is `stopped: <step> failed: <cause>`, taps in front of it included, because they saw part of the stream. Every later step is `not run: <step> failed: <cause>` |
| A barrier at its seal | The barrier has its error. The taps of its segment keep their results: they saw the whole stream. Later steps are `not run` |
| `write` on a record | The segment stops as for a stage when the reader sees the failure between two batches. When the stream ended first, the taps keep their results. Later steps are `not run` either way |
| The source in the first segment | `run_chain` returns the error, as `pass::run_isolated` does |
| The source in a later segment | The steps of that segment are `stopped`, later steps `not run` |

Steps of an earlier segment keep their results. A written BAM is renamed into place only
when its segment finished. A failed or stopped writer leaves no file and no temporary
file. A tap that failed on its own before its segment stopped keeps its own error.

Tests: `crates/core/src/stage/tests.rs` (synthetic stages, each row of the table) and
`crates/cli/tests/chain.rs::a_failing_tap_fails_alone_and_a_failing_stage_stops_its_segment`
(`CollectGcBiasMetrics` without a reference, `MarkDuplicates` on a BAM without a read
group).

## 8. Chained equals staged

Staged means `MarkDuplicates I=in.bam O=marked.bam M=dup.txt`, then each tool on
`marked.bam`. Chained means one `__chain` command.

What is equal:

- Each tap's metrics file after its `StringHeader` lines, byte for byte.
- The `MarkDuplicates` metrics file after its `StringHeader` lines.
- Every record of the written BAM, byte for byte, and every header line with the `CL`
  field of `@PG` set aside.
- With `PROGRAM_GROUP_COMMAND_LINE` given, the whole BAM, its `.bai` and its `.md5`.

`crates/cli/tests/chain.rs` checks this for 8 taps (`CollectAlignmentSummaryMetrics`,
`CollectInsertSizeMetrics`, `CollectQualityYieldMetrics`, `CollectDuplicateMetrics`,
`QualityScoreDistribution`, `MeanQualityByCycle`, `CollectBaseDistributionByCycle`,
`CollectGcBiasMetrics` with its summary file) on `tests/data/markdup.bam` (2 option
sets), the query-name sorted `tests/data/duplicate_marking.bam`, corpus t1 (2 option sets) and the foundation corpus
(2 option sets, one spilling with `MAX_RECORDS_IN_RAM=500`). The corpus cases are
skipped when `tests/corpus/out` is absent.

What differs on purpose (also in KNOWN_DIFFERENCES.md):

- The `StringHeader` command line. A tap writes `INPUT=<the chain's source>` and, behind
  a stage, ` CHAIN=<stages in front, joined by +>`. A tap on the source writes the line
  the tool alone writes.
- The `@PG` `CL` of the BAM is the step's command line, which has no `OUTPUT`.
- `MarkDuplicates` in a chain writes no BAM. `write` does.

## 9. Measurements

Mac Studio, Apple M1 Ultra, 20 cores, 64 GB. Release build, 8 threads, `/usr/bin/time -l`.
Other builds ran on the machine during both runs.

Transform-free chain against `qc` (the W03 gate, 5 percent). Input: chr20, 21 and 22 of
the HG002 exome, 573 MB, 8.4 M records. Four collectors on one read (quality yield,
alignment summary, insert size, duplicate metrics), the pass of `qc --preset wgs`.
Seven runs each, interleaved.

| Run | Median wall | Range | Against `qc` | Max RSS |
|---|---:|---|---:|---:|
| `qc --preset wgs` (also writes the verdict and the page) | 6.52 s | 6.43 to 6.88 | | 43 MiB |
| `__chain`, lowered to `pass::run_isolated` | 6.43 s | 6.36 to 6.48 | -1.4 % | 42 MiB |
| `__chain LOWER=false`, the stage executor | 6.59 s | 6.55 to 6.91 | +1.1 % | 42 MiB (78 MiB in run 1) |

Both forms are inside the gate. The kill rule of docs/07 does not apply. A chain of taps
is lowered anyway: it is the same code as `qc`, and it was 2.5 percent faster than the
stage executor here. The four files of the three runs are equal after their
`StringHeader` lines.

Chain against the staged sequence on the HG002 exome (10,339,495,604 bytes, 150,538,907
records, coordinate sorted). One run each, default options, `PROGRAM_GROUP_COMMAND_LINE`
given on both sides so the BAMs can be compared whole.

| Run | Wall | Max RSS |
|---|---:|---:|
| `MarkDuplicates` to a file (staged, step 1) | 300.8 s | 764 MiB |
| The four collectors on the marked BAM, one shared read (staged, step 2) | 107.7 s | 72 MiB |
| Staged sequence, both steps | 408.4 s | 764 MiB |
| `__chain`: `MarkDuplicates` + 4 collectors + `write` | 345.2 s | 770 MiB |

The chain takes 15.5 percent less wall time than the staged sequence and 6 MiB more
memory at its peak. Both BAMs are 10,210,246,600 bytes and `cmp` finds no difference.
The five metrics files are equal after their `StringHeader` lines. The chain reads the
10 GB input twice and the marked BAM never; the staged sequence reads the input twice and
the marked BAM once. The chain's second segment took about 44 s longer than the tool's
second read. The tool reuses one record buffer; a chain allocates each record of a batch
(opinion: that is the cost; an arena batch in W05 would show it).

`MarkDuplicates` alone was measured at 299 s and 783 MiB before this wave
(docs/benchmarks.md), so the tool did not slow down.

## 10. Departures from docs/07 section 2

| docs/07 | As built | Reason |
|---|---|---|
| `Source::open` returns `Box<dyn Records>` | `Source::read` takes a closure | The CRAM decoder (`cram::decode`) runs its workers in a thread scope and hands the stream to a closure. An owned reader cannot leave that scope |
| `Tool` gains `stage`, `source`, `sink`, `tier` fields | `stage::stage_of(tool)`, a table by name | A new field breaks every `Tool { ... }` literal, including those on the two branches built at the same time. One tool has a stage |
| `Transform::requires()` checked at plan time | A stage is built from the `Frame` of its position and fails there with the tool's own message | The planner is W05. `MarkDuplicates` accepts two orders, which one `Order` value cannot say |
| `Worker` and `Fanout` of `pass.rs` reused for the taps | `SinkWorker` in `stage/exec.rs`, the same shape | `pass::Worker` carries `bam::Record` batches. A sink takes `RawRecord` batches and needs the done, failed and stopped signals |
| One reader thread and one transform thread | The transforms run on the reader thread | The tool does the same. A second thread adds one queue hop per batch. Measured on the exome: the chain takes 345 s, the tool alone 301 s, the staged sequence 408 s (section 9) |
| 4 deflate workers per writer | `threads` workers, as the tool uses | The written BAM is then the tool's BAM byte for byte at any thread count, and compression is the slowest part of the chain |
| `Barrier::seal` on a `Send` trait | `Transform` and `Barrier` are not `Send` | They run on the calling thread and borrow the chain's `TempDir` |
| `write` picks BAM, SAM, CRAM or FASTQ by extension, `--index` | BAM only, `CREATE_INDEX`, `CREATE_MD5_FILE` | SAM text is W05, the CRAM writer W15, FASTQ W22 |
| The `# INPUT=` line names the source and the step list | `INPUT=<source>` plus ` CHAIN=<stages>` behind a stage | The step grammar is W04. A tap on the source keeps today's line |
| `run_chain(&ChainPlan, &dyn Progress) -> Result<Vec<Products>>` | `run_chain(ChainPlan, Option<OnProgress>) -> Result<Vec<Result<Products>>>` | A step can fail alone, so each step has its own result. `OnProgress` is the callback `qc` already has |

Added beyond the list of docs/07: `Emit` takes `&mut RawRecord` so a batch can take the
record without a copy. `Transform::products` carries the `MarkDuplicates` metrics.
`Frame` has `cram` (htsjdk rebuilds `MD` and `NM` for CRAM records), `label` and
`reads`. `Sink::done` carries `STOP_AFTER`.

## 11. What comes next

| Wave | Adds |
|---|---|
| W04 | The `chain` verb with the `+` grammar and did-you-mean per step, `plan` for a chain, the TOML plan file, Python `bm.chain` and `bm.step`, `Run.stages` from `Products::counters`, the final form of the `StringHeader` line, product naming (`PREFIX.<label>.<ext>`). `__chain` and `crates/cli/src/chain.rs` are replaced. `Step::tool_tap`, `Step::tool_stage`, `Step::write` and `run_chain` stay |
| W05 | The sort barrier (`Next::Records`), the planner (`Frame` per position, `requires`, free-space check), the spool for a replay behind another stage, `--sort-memory`, `write --index` through the BAI writer |
| W09 on | Each BAM-writing tool registers its stage in `stage_of` |
| W14 | `qc` as the last step of a chain, manifest 2, resume per step, kept products, the spool for tools that read their input their own way |
| W15, W22 | `write` to CRAM, FASTQ |
| W46 | Variant streams (`Kind`), pileup and capture sinks |

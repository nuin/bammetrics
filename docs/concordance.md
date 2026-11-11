# Concordance corpus and harness

How the synthetic corpus is built and how every case is compared with Picard 3.5.0. The validation model and exit criteria are in [04_picard_parity_plan.md](04_picard_parity_plan.md) section 5.

## Regenerate the test BAM

```bash
uv run --no-project --with pysam python tests/make_bam.py   # writes tests/data/insert_size.bam
uv run --no-project --with pysam python tests/make_orientation_bam.py   # writes tests/data/orientation.bam
python3 tests/data/make_ref.py                               # writes tests/data/ref_100k.fa and .fai
```

`tests/data/targets.interval_list` (5 targets over the fixture reads) and the two
`quality_yield_*.txt` files (`CollectQualityYieldMetrics` on each BAM) complete the
committed fixtures.

## Golden suite

```bash
make golden          # cargo test --release -p bammetrics-cli --test golden
make golden-update   # rewrite tests/golden after an intended command-line change
```

`crates/cli/tests/golden.rs` runs every documented command form from the repository
root: each registry tool in Picard syntax with and without `O=`, the five 0.12
subcommands (verbs since 0.13.0) with and without `-O`, the verb grammar's spellings,
formats, prefixes and info verbs, `check`, `report`, `qc` and `plan`, and the documented
error paths, 215 cases on the committed fixtures and 22 more on corpus t1 when
`GOLDEN_CORPUS=1` is set and `tests/corpus/out` exists. A case that writes into `.`
(`collect` without `-o`) runs in its output directory and records the fixture paths as
`ROOT/...`. One
record per case under `tests/golden/` holds the argv, the exit code, stdout, stderr and
every file written (the output directory reads `OUT`, the `Started on` header line
`<time>`). Any byte change fails, with the first differing line named. The records were
first written by the 0.6.0 binary, so a later command-line change has to keep them or
show its change in the diff of `tests/golden`. `BAMMETRICS_BIN=path make golden` runs
another binary against the same records. `make concordance` runs the suite before the
pytest cases.

## Corpus and harness

```bash
make corpus             # tests/corpus/out: reference, intervals, t0-t3 BAMs (11 MB, ~30 s)
make reference-metrics  # Picard 3.5.0 output: tests/corpus/out/picard/<case>.<metric>.txt
make concordance        # cargo build --release, then pytest tests/concordance
python3 tests/data/make_corpus.py out_dir   # seeded multi-sample corpus: 2 samples x 2 libraries x 3 read groups
```

`tests/concordance/compare.py` has `compare_metrics(tool, bam, picard_args)`: it runs Picard and
bammetrics on one BAM and asserts every column of every row, and the histogram section,
equal as printed strings. `CLAUDE.md` lists the steps for adding a tool.

Needs samtools, bwa, `picard` (3.5.0, checked) and uv. Nothing in
`tests/corpus/out/` is committed; `generate.py` rebuilds it from one seed. Records
are identical across runs. BAM bytes differ only in `@PG` lines (paths, and the
random ID suffixes `samtools merge` adds).

The corpus (`tests/corpus/generate.py`):

- `reference.fa` (+ `.fai`, `.dict`, bwa index): random chr1 300 kb, chr2 200 kb,
  chrM 16 kb. A 3 kb exact copy and a 2 kb copy with 2% divergence of chr1 sit in
  chr2, so some reads get MAPQ 0 or intermediate MAPQ and secondary records.
- `targets.interval_list` / `targets.bed`: 120 targets of 100-400 bp on chr1 and chr2.
  `baits.*`: targets padded by 50 bp and merged (110 baits). Interval lists carry the
  `.dict` header.
- t0: tiny pysam-built BAMs, one Picard rule each: even count with a half-integer
  median, odd count, single pair (SD `?`), mode tie, all TANDEM, RF at exactly 5%
  (dropped) and at 5.3% (kept), a wide spread with a long tail, the record filters
  (duplicate, secondary, supplementary, cross-contig pair with non-zero TLEN), and
  no pairs. The committed fixtures in `tests/data/` are t0 too.
- t1: 20k clean FR pairs, 2x150, two libraries (libA, libB) in sample SYN1,
  70% of fragments on target, binned qualities.
- t2: 50k pairs with everything t1 lacks: 2% adapter read-through, 5% reads with a
  1-3 bp indel, qualities 2-41 and Q2 tails, 3% RF and 1% TANDEM pairs, 2%
  cross-contig chimeras, 2% split reads (supplementary records), 1% unmapped
  mates, 5% PCR duplicates marked by `picard MarkDuplicates`, 1% QC-fail flags,
  `bwa mem -a` secondary records.
- t3: 10k single-end reads.
- t0_hs: pysam-built pairs and reads on the corpus targets for `CollectHsMetrics`:
  overlapping mates in every clipping shape htsjdk handles (equal starts, `=`/`X`
  elements, an insertion, a trailing soft clip, a deletion at the clip, a mate on
  another contig), duplicates, an adapter read at MAPQ 0, MAPQ 0 and 10, Q2 and Q3
  bases, QC-fail, secondary, supplementary and mate-unmapped records, near-bait and
  off-bait reads, indels, and an 8-deep stack. `hs_targets.interval_list` holds the
  first targets unsorted, duplicated, overlapping, abutting, empty, and on a contig
  outside the dictionary.
- t0_asm, t0_asm_cycles, t0_asm_empty: pysam-built BAMs on the corpus reference
  for `CollectAlignmentSummaryMetrics`: one record per rule (mismatches, IUPAC and
  `=` read bases, clips, indels, `XN` noise, adapter reads at MAPQ 0 and unmapped,
  chimeras by insert size, contig, orientation and SA tag, MQ gating, QC-fail,
  secondary and supplementary records, a read off the end of chrM), the BAD_CYCLES
  rules, and an empty BAM (one all-zero UNPAIRED row).

Each library is aligned with `bwa mem -R`, run through `samtools fixmate -m`
(so every record has `MC`), sorted, merged, then duplicate-marked.

The harness (`tests/concordance/`):

- `cases.py` lists every case: BAM, metric, Picard arguments and the matching
  bammetrics arguments. `run_picard.py` and the tests both read it.
- `run_picard.py` runs Picard per case and records in `manifest.json` whether Picard
  wrote a file. `CollectInsertSizeMetrics` will not start without `Rscript`, so a
  no-op one in `target/fakeR` goes first on `PATH`; `H=/dev/null` is accepted.
- `test_concordance.py` runs `target/release/bammetrics` per case and compares
  with `picard_metrics.compare_files`: every `METRICS CLASS` field, the column
  order, and every `HISTOGRAM` cell. Rows pair by each metric's key columns
  (`Metric.key_columns`: `PAIR_ORIENTATION` plus `SAMPLE`, `LIBRARY`, `READ_GROUP`
  for insert size, `CATEGORY` for alignment summary; quality yield and HS metrics
  have one row).
- Comparison rules: integer literals compare exactly; other numbers compare with
  `math.isclose(rel_tol=1e-9, abs_tol=1e-6)` after our value is printed the way
  htsjdk `FormatUtil` prints doubles (at most 6 decimals, trailing zeros dropped).
  Picard's printed value is already rounded to 6 decimals, so rounding ours first
  keeps a 7th-decimal difference from failing. Empty fields equal only empty
  fields; `?` (FormatUtil's NaN and infinity) equals `?` or NaN only. Every
  mismatching field is reported with both values.
- Native runs repeat at `--threads 1` and `--threads 8`. `test_collect` runs
  `bammetrics collect` on t1, t2, t3 and t0_asm at both thread counts and checks each
  file it writes against the Picard reference of the matching single-tool case.
- `PER_BASE_COVERAGE` cases (`t1_pb`, `t2_pb`, `t0_hs_pb`, `t0_hs_pb_messy`, each at the
  defaults and with `MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0`, plus `t1_pb_q0_gz`)
  compare the whole per-base file byte for byte (after gunzip for the `.gz` case).
  `test_pipeline_coverage_command` runs the clinical pipeline's coverage command argument
  for argument on t1 and t2.
- The run ends with one line per case: fields matched / total, and whether the
  file body is byte-identical to Picard's.

## Oracle tiers

Each tool is compared with Picard 3.5.0 at one tier.

| Tier | Meaning | Harness call |
|---|---|---|
| A | Picard's file byte for byte after the `StringHeader` lines | `compare_metrics`, `compare_outputs`, `compare_prefixed` |
| B | Picard's output decoded record by record, `@PG` normalised | `compare_bam`, `compare_vcf` |
| C | Picard on synthetic data with a written tolerance | per tool |
| D | no real-data oracle: round trip and parser tests | per tool |
| E | not a Picard tool: planted truth, then a named tool within a stated tolerance | per tool |

Metrics files, interval lists, `.dict` files and standard output are tier A. BAM, CRAM and
VCF outputs are tier B (next two sections).

## Tools that write reads or variants

`compare_bam(tool, inputs, picard_args, bam_args=None, workdir=None, outcome="equal",
input_key="I", output_key="O", metrics_key=None, suffix=".bam", reference=None, label=None,
record_bytes=False)` in `tests/concordance/compare.py` runs Picard and bammetrics with
`O=<dir>/<side><suffix>` and reads both files back with `samtools view --no-PG`, one character
per byte.

- Header lines other than `@PG` must be equal, in order. The `DT` field of an `@RG` line is
  cut first: htsjdk rewrites it in the JVM's time zone.
- The `@PG` lines must carry the same `ID` and `PP`, in order. `VN` and `CL` name the program
  that wrote the file and are not compared.
- Every record must be equal as `samtools view` prints it: every field, every tag, file order.
  The text is split at `\n` only. `str.splitlines` also splits at `\r`, `\x0b`, `\x0c`, `\x1c`
  to `\x1e` and `\x85`, so a tag value holding one of them would be cut in two and the byte
  itself never compared.
- Each BAM must end with the 28-byte BGZF end-of-file block. `samtools view` exits 0 without
  it and prints every record, so a writer that never wrote it would pass otherwise. The same
  check runs on a bgzip VCF and under `record_bytes`.
- The failure names the first differing header line or record by its index, from 0, with both
  lines.
- `metrics_key="M"` also compares a metrics file as `compare_metrics` does.
- `outcome="fail"` requires both tools to fail, bammetrics to leave no output, and neither to
  write a metrics file. `outcome="same"` lets Picard decide between the two.
- `record_bytes=True` also compares each record as bytes after BGZF inflation. That catches
  what `samtools view` does not print: the `bin` field and the width of an integer tag.
- `suffix=".cram"` with `reference=` reads a CRAM back. It returns `Written` (`outcome`,
  Picard's `header` and `records`, the parsed `metrics`, `cells`).

`tests/concordance/test_mark_duplicates.py` runs its 160 cases through `compare_bam`.

`compare_vcf(tool, inputs, picard_args, bam_args=None, workdir=None, outcome="equal",
input_key="I", output_key="O", suffix=".vcf", drop=(), label=None)` does the same for a VCF,
plain or gzip. The header lines must be equal in order, without the date and command lines
(`##fileDate=`, any key ending in `CommandLine` or `Command`) and without the keys in `drop`.
Every record line must be equal, split at `\n` only. `drop` takes the lines one Picard writer adds that name the
run, for example `drop=["liftOverTime"]`. Picard 3.5.0 `SortVcf` and `MakeSitesOnlyVcf` add no
such line and write `##FORMAT` before `##contig`. No registered tool writes a VCF yet, so
`compare_vcf` is tested on stand-in tools only.

## BAM bytes

Decision (2026-11-10, wave W02 of the full parity program plan): a tool that writes a BAM is
compared at tier B,
decoded records. BAM bytes are not compared and bammetrics does not try to write Picard's
bytes. The same holds for CRAM and for a gzip VCF.

The spike that decided it is `make deflater-spike`
(`benches/deflater_spike.py` and the package `benches/deflater_spike`, which is outside the
workspace, so no shipped crate depends on it). It runs `MarkDuplicates` three times per input:
Picard at its defaults (`COMPRESSION_LEVEL=5`), Picard with `USE_JDK_DEFLATER=true`, and
bammetrics. It then inflates every BGZF block Picard wrote, deflates it again with each
candidate at each level, and counts the blocks that come out byte for byte.

Run of 2026-11-10: Mac Studio M1 Ultra, OpenJDK 27.0.1 aarch64, Picard 3.5.0, zlib-rs 0.6.8,
miniz_oxide 0.8.9, libdeflater 1.26.1, on `tests/data/markdup.bam` (36 records) and corpus
`t2.bam` (101,099 records).

Which deflater Picard used, from its log:

| Picard run | Log line | Library |
|---|---|---|
| defaults | `Intel Deflater not supported, using Java.util.zip.Deflater` | `/usr/lib/libz.1.dylib`, zlib 1.2.12 |
| `USE_JDK_DEFLATER=true` | `libdeflate is available; using libdeflate for DEFLATE compression` | libdeflate, loaded by htsjdk |

The Intel GKL library in `picard.jar` is x86_64 only (`incompatible architecture (have
'x86_64', need 'arm64')`), so it was not measured. On x86_64 it is a third stream.

Whole files:

| Input | File | Bytes | Data blocks | Largest block | Against Picard defaults |
|---|---|---:|---:|---:|---|
| markdup.bam | Picard defaults | 2,022 | 1 | 9,461 | equal |
| markdup.bam | Picard `USE_JDK_DEFLATER=true` | 2,000 | 1 | 9,443 | differ at byte 16 |
| markdup.bam | bammetrics | 1,355 | 1 | 7,331 | differ at byte 16 |
| t2.bam | Picard defaults | 8,375,031 | 585 | 65,498 | equal |
| t2.bam | Picard `USE_JDK_DEFLATER=true` | 8,414,929 | 585 | 65,498 | differ at byte 16 |
| t2.bam | bammetrics | 8,180,379 | 585 | 65,495 | differ at byte 16 |

Blocks each candidate reproduces, at its best level:

| Picard stream | Candidate | Levels tried | markdup.bam | t2.bam |
|---|---|---|---:|---:|
| defaults | zlib-rs 0.6.8 | 1 to 9 | 0 of 1 | 0 of 585 |
| defaults | miniz_oxide 0.8.9 | 1 to 10 | 0 of 1 | 0 of 585 |
| defaults | C libdeflate 1.26.1 | 1 to 12 | 0 of 1 | 0 of 585 |
| defaults | C zlib 1.2.12, memLevel 8 | 1 to 9 | 1 of 1 (level 5) | 585 of 585 (level 5) |
| `USE_JDK_DEFLATER=true` | zlib-rs 0.6.8 | 1 to 9 | 0 of 1 | 0 of 585 |
| `USE_JDK_DEFLATER=true` | miniz_oxide 0.8.9 | 1 to 10 | 0 of 1 | 0 of 585 |
| `USE_JDK_DEFLATER=true` | C zlib 1.2.12, memLevel 8 | 1 to 9 | 0 of 1 | 0 of 585 |
| `USE_JDK_DEFLATER=true` | C libdeflate 1.26.1 | 1 to 12 | 1 of 1 (level 5) | 585 of 585 (level 5) |

Why tier B:

1. Neither pure-Rust deflater writes a single Picard block. zlib-rs is the deflater under
   noodles-bgzf 0.52, so every BAM bammetrics writes goes through it.
2. Picard has no one byte stream. Picard at its defaults and with `USE_JDK_DEFLATER=true`
   differ by 39,898 bytes on t2, on one machine. Two runs of the same command write the same
   blocks. The default stream depends on the platform (Intel GKL on x86_64, the machine's libz
   elsewhere) and on the libz and libdeflate versions.
3. Both streams are reachable only through a C library pinned to Picard's version (libz
   1.2.12 here, libdeflate at level 5), which the shipped crates do not link today.
4. Bytes would still differ. htsjdk cuts blocks at 65,498 bytes and noodles-bgzf at 65,495.
   A tool that writes a `@PG` line writes its own `VN` and `CL`, so the first block differs
   in content (9,461 against 7,331 inflated bytes on markdup.bam).

What tier B still proves, and what it does not:

- Records after inflation, bammetrics against Picard defaults: 101,099 of 101,099 byte-equal on
  t2, 34 of 36 on markdup.bam. `compare_bam(record_bytes=True)` is that check.
- `test_mark_duplicates.py` with `record_bytes=True` fails 93 of its 160 tests, all on one
  field: `bin` of a record without a reference (`refID` -1). The inputs carry 4680 there
  (`reg2bin(-1, 0)`, as samtools and pysam write it; two foundation records carry 4683),
  Picard writes 0, bammetrics copies the input. A probe planted `bin` 1 on two placed records
  of markdup.bam: Picard wrote 4681, bammetrics 1. So htsjdk computes `bin` for every record
  it writes and bammetrics copies it, which KNOWN_DIFFERENCES.md already lists under
  MarkDuplicates. `samtools view` does not print `bin`. `MarkDuplicates` therefore keeps
  `record_bytes` off; a later writer that computes `bin` can turn it on.
- A BAI that Picard writes beside its BAM (`CREATE_INDEX=true`) holds virtual offsets into
  Picard's compressed blocks. It cannot be byte-equal to an index of a BAM bammetrics wrote.
  `BuildBamIndex` on one given BAM can be tier A. It is: the next section.

## BAM index

A `.bai` is tier A for a given BAM. `tests/concordance/test_build_bam_index.py` runs
`BuildBamIndex` through both tools on the same file and compares the two indexes with
`==` on their bytes: the six coordinate-sorted BAMs of `tests/data`, the foundation
corpus, every coordinate-sorted tier BAM (t0 to t3), the exome slice and 26 BAMs that
`index_bams.py` writes byte by byte. 107 tests, 68 indexes compared.

`index_bams.py` packs each record with `struct` and each BGZF block with `zlib`, because
an index depends on where the blocks end and pysam and samtools choose that themselves.
Its cases: blocks of 250 to 5000 bytes, records that span blocks, a header that shares
its block with records, bins met again many blocks later (95 bins with more than one
chunk over the module), a wrong bin field in every record, placed unmapped reads, reads
on the edges of the 16384-base windows, unplaced reads only, no record, no reference,
references without reads before, between and after the one with reads, a reference of
2^29 + 40000 bases with reads up to position 2^29, a bin of each of the six levels, a
file without its EOF block, a reference of length 0 (last, between two references with
reads, alone).

What running Picard 3.5.0 showed, each pinned by a test:

1. htsjdk computes the bin of every record and never reads the bin field
   (`stored_bin_zero`, `stored_bin_wrong`).
2. Without `OUTPUT` the index goes into the current directory, under the file name of
   `INPUT`. The help says beside the BAM (`test_default_output_name`).
3. At the default `VALIDATION_STRINGENCY=STRICT` Picard refuses the foundation corpus (a
   supplementary flag on an unmapped read). The corpora run with `SILENT` on both sides,
   as in every other module. The three `lenient` BAMs show that bammetrics, which
   validates nothing, writes the index Picard writes under `SILENT`.
4. An end past position 2^29 is `Read position too high for BAI bin indexing.`, so no
   bin above 37448 and no window above 32767 is ever written.
5. An empty BGZF block in the middle of a file ends htsjdk's read: Picard indexes 101 of
   300 records and exits 0. bammetrics indexes all 300 (KNOWN_DIFFERENCES.md).
6. Picard's writer and its reader disagree on one offset. `MarkDuplicates
   CREATE_INDEX=true` on corpus t2 ends the last chunk at `8374551:1514`, inside the
   last block; `BuildBamIndex` of that BAM ends it at `8375239:0`, the EOF block. A BAM
   that ends with unplaced reads shows no difference
   (`test_index_and_md5_beside_the_output`).
7. htsjdk reads a header with `@SQ LN:0` and Picard indexes the file. noodles holds a
   length as a non-zero number and its header reader refuses it, so `bai/header.rs`
   reads the header for `BuildBamIndex` (`empty_reference_between`,
   `empty_reference_last`, `empty_reference_only`). Picard at `STRICT` refuses a record
   placed on such a reference and indexes it under `SILENT`, as bammetrics does
   (`records_on_an_empty_reference`). samtools 1.x prints no header with `LN:0` and
   finds no region of the file; `samtools idxstats` reads the bammetrics index.

16 BAMs both tools refuse (`index_bams.refused`): each of the five headers that are not
`SO:coordinate`, references out of order, five alignments past 2^29, a position without
a reference, a start past the reference, a reference the header lacks, a truncated file,
a record whose `block_size` is 0xFFFFFF00 in a file of 227 bytes (Picard: `Invalid
record length: -256`; bammetrics: `truncated BAM`). Picard leaves an empty `OUTPUT`;
bammetrics leaves none and keeps the file that was there.
`test_a_block_size_past_the_file_costs_the_bytes_present` runs `build-bam-index`,
`mark-duplicates` and `chain + write --index` on that file and fails above 64 MiB of
peak RSS: each takes 5 MiB. A buffer of `block_size` bytes took 4.3 GB.

samtools reads each index as its own (`cross_check`): `samtools idxstats` and `samtools
view -c` on four regions per reference give the same numbers through the bammetrics
index and through `samtools index`. The two files differ: samtools writes 195,448 bytes
for the 573 MB exome slice, Picard and bammetrics 346,128.

A BAM bammetrics wrote holds its own compressed bytes (tier B), and its index is tier A
against Picard's `BuildBamIndex` of that file: `chain BAM + write OUT --index` on five
corpora at compression levels 1 and 5, and `MarkDuplicates CREATE_INDEX=true` on five
inputs (`test_mark_duplicates.py`). The `.md5` beside it is the MD5 of the BAM, 32 hex
digits and no line end, as Picard writes it.

`tests/data/bai` holds the six fixture indexes Picard wrote (`tests/data/make_bai.py`).
`crates/cli/tests/build_bam_index.rs` compares both faces and eight spellings with them,
so `cargo test` holds Picard's bytes without a JVM; `test_fixtures` fails when those files
go stale. The golden suite records the index of a committed BAM as hex.

On real data (docs/benchmarks.md): the index of the HG002 exome BAM, 7,186,968 bytes, is
equal to Picard's under `cmp`, and so is the index of its chr20 to chr22 slice.

## Picard queue

One Picard JVM runs at a time per user and machine. `run_picard.picard_env()` puts
`target/picard-queue/<key>/picard` first on `PATH`. `<key>` is 12 hex characters of the
SHA-256 of the script's text, which holds the interpreter path. A script at one path never
changes, so two harness processes started from different interpreters (a uv environment and
`.venv`) never replace each other's script. That script runs
`tests/concordance/picard_lock.py`, which:

1. takes an exclusive `flock` on `~/.cache/bammetrics/picard.lock` (file mode 0600,
   never through a symbolic link);
2. writes its pid, the time and the command into that file;
3. waits while `/bin/ps -eo pid=,args=` lists a `java` process running `picard.jar`, which
   covers a harness older than the lock and a Picard run started by hand;
4. replaces itself with the real `picard`. The JVM inherits the descriptor, so the kernel
   releases the lock when the JVM exits, however it exits.

Two pytest processes, two worktrees or two agents therefore never run two JVMs at once.
`JAVA_OPTS`, the arguments, standard streams and the exit code are the real command's. The
Homebrew `picard` wrapper and the CI wrapper both expand `JAVA_OPTS`. A cached result
(`target/fuzz/cache`, the `run_picard.py` manifest) starts no process and takes no lock.

```bash
python3 tests/concordance/picard_lock.py --status          # free, or who holds the lock
python3 tests/concordance/picard_lock.py SortSam --help    # one queued Picard call
```

| Variable | Effect |
|---|---|
| `PICARD_LOCK_FILE` | another lock file |
| `PICARD_LOCK_TIMEOUT` | seconds to wait, then exit 75; unset waits without limit; a value that is not a finite number from 0 up exits 2 with one line |
| `PICARD_LOCK_WAIT_FOREIGN` | `0` skips step 3 |

`benches/bench.py` and `benches/wgs_bench.py` hold the same lock around each timed Picard run
(`picard_lock.held`), so the wait is outside the measured time. The queue adds 123 ms per call
on the Mac Studio (10 ms direct, 133 ms queued, mean of 50 calls of a no-op `picard`; 78 ms of
it is `/bin/ps`). `make picard-fields` reads class files with `java` for about 1 s and is not
queued. `tests/harness/test_picard_lock.py` fails when a file under `tests`, `benches` or
`tools` starts `picard` without `env=picard_env()`. It reads a list literal that starts with
`"picard"` as a command, except a list of two or more output format names
(`format=["picard", "json"]` of the Python `chain`).

## Registration rule

`crates/cli/tests/registration.rs` fails while a tool in `TOOLS` lacks any of:

1. a verb in `VERBS`;
2. an argument list in `crates/core/tests/glossary.rs`, or a line in its `EXEMPT` table. The
   test reads the keys of that file's `args()` map and the first field of each `EXEMPT` line.
   A name in a comment or in an argument value does not count;
3. a concordance module that names it: `tests/concordance/test_<verb>.py` or
   `test_<picard_name_in_snake_case>.py`, a line in `CONCORDANCE` (tool, shared module), or a
   line in `NO_PICARD_ORACLE` with the reason;
4. its name in `docs/cli.md`;
5. a row of the README tool table that starts with its name.

One failure lists every gap of every tool with the edit that closes it. `BuildBamIndex`
(wave W05 of docs/07) is tool 34, with `tests/concordance/test_build_bam_index.py` as its
module and an argument list in the glossary test (it returns rows without a class). At
0.22.0 the 33 tools
pass: 13 have a module named after them, 18 are mapped in `CONCORDANCE` to a shared module
(`test_read_level.py` 6, `test_wgs_metrics.py` 3, `test_variant_tools.py` 3, six others 1
each), and `CollectCoverageGaps` and `CollectSexCheckMetrics` are in `NO_PICARD_ORACLE`
(docs/06, not Picard tools). `FifoBuffer` (wave W04 of docs/07) is tool 33. Its module is
`tests/concordance/test_fifo_buffer.py`. It is the one line of `EXEMPT`: it copies
standard input to standard output and writes no metrics column.

## Spec script

```bash
python3 tools/spec_agent.py SortSam               # tab separated
python3 tools/spec_agent.py SortSam --markdown    # a table
python3 tools/spec_agent.py --from-help saved.txt
```

`tools/spec_agent.py` prints the option table of `picard <Tool> --help` as a builder's
starting spec: `KEY`, `SHORT`, `TYPE` (Picard's Java type), `OPT` (the registry row kind it
suggests: bool, int, float, path, text), `DEFAULT`, `NOTES` and `DOC`. `NOTES` holds
`required`, `repeated`, `advanced`, `values A|B`, `excludes X`, `common` for a name the
registry's `COMMON` rows take, and `parser` for the four options of Picard's command line
parser. `values` comes from `Possible values: {A, B}`, or from the `A (what it does)` lines
Picard prints for an enum with a text per value (`SORT_ORDER` of `SortSam`, `MODE` of
`CheckDuplicateMarking`). A name Picard does not know exits 1 with one line, for example
`spec_agent: SortBam: not a Picard tool; Picard suggests SortSam`. Picard's colour codes are
dropped. One JVM runs, through the queue, with `JAVA_OPTS=-Xmx1g` unless `JAVA_OPTS` is set.
The table says what Picard accepts. What an option does to the output comes from the Picard
3.5.0 source.

## Dependency checks

```bash
cargo install cargo-deny --locked    # once; 0.20.2 was used
make deny                            # cargo deny check
```

`deny.toml` checks the 105 packages of `Cargo.lock`. Result at 0.21.0: advisories ok, bans
ok, licenses ok, sources ok. 0.22.0 adds no package: its `Cargo.lock` differs from the
0.21.0 one in the six workspace versions only. `cargo-deny` was not rerun for 0.22.0.

- Licences allowed: MIT, Apache-2.0, Unicode-3.0 (unicode-ident), Zlib (zlib-rs),
  bzip2-1.0.6 (libbz2-rs-sys under noodles-cram) and Apache-2.0 WITH LLVM-exception
  (target-lexicon, a build dependency of pyo3). The last two were the only findings of the
  first run. Both are permissive, so they are allowed with the reason in the file.
- Bans: two versions of one crate, a wildcard version, and seven network or TLS crates
  (docs/05 section 4.9: no network call).
- Sources: crates.io only, no git dependency.
- Advisories: the RustSec database, yanked crates denied, no exception.

cargo-deny is not installed on the CI runners, so `make deny` is not a CI step. A duplicate
version that a later dependency brings in needs a `skip` entry with its reason.

## Harness tests

```bash
make harness-test    # pytest tests/harness: 79 tests, about 12 s, needs samtools
```

`tests/harness` tests the harness with stand-in tools, so it starts no JVM and needs no release
build: the queue (`test_picard_lock.py`), `compare_bam` and `compare_vcf` on files with one
planted difference each (`test_compare_bam.py`: one flag, one tag, one base quality, one `@PG`
ID, one `bin` field, one separator byte in a tag, a missing end-of-file block), and the spec
script on the `--help` text of `MarkDuplicates` and `SortSam` (`test_spec_agent.py`). It is a
step of `make concordance`, before the pytest run.

## Python module tests

```bash
make python-test     # maturin develop --release, then pytest tests/python
```

`tests/python` runs against the maturin-built module. `test_cli_parity.py` compares every
field the module returns on corpus BAM t1 with the CLI's Picard-format output, and
`test_api.py` checks `to_dict()`, the error types and that `bammetrics.pyi` lists every field.
Plain `cargo build` of the cdylib is not a usable module on macOS; use maturin.

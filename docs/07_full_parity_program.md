# 07 Full parity program: every Picard tool, chaining, formats, reports, more analyses

Status: program plan for the maintainer. Date: 2026-11-09. Base: main at 0.17.0 (23 tools), release/0.18.0 in flight.

Sources: the 91-tool survey (`picard <Tool> --help` on 3.5.0, one JVM at a time), the chaining design, the format and PDF spikes, the beyond-Picard survey, three draft plans and three judge reviews. Plan 1 (foundations first) won. Its grafts from the other two and the gaps the judges listed are closed here.

## 0. Answer

- 48 waves (W00 to W47), 572 agents, versions 0.18.0 to 1.7.0. 1.0.0 is the parity gate at W40.
- All 91 unbuilt Picard tools land by 0.57.0. 15 analyses beyond Picard land in W07, W08 and W41 to W47.
- Chaining ships at 0.22.0 (`chain`), with sort at 0.23.0, resume at 0.32.0, VCF streams at 1.6.0.
- Formats core ships at 0.20.0. PDF at 0.24.0. Parquet, XLSX, DOCX at 0.39.0.
- Every foundation lands before the first tool that needs it, with one proving tool. No tool is ported file-to-file and refactored later.
- 12 tools have no real-data oracle. They ship labelled `best effort` in `tools` output and KNOWN_DIFFERENCES.md.

Decisions the owner must take, by the wave each one blocks:

| # | Decision | Blocks |
|---|---|---|
| Q1 | Agent cap: 40 total is read as per workflow (14 at most each). Confirm. | W02 |
| Q2 | Verb name: `chain`, with `run` as an alias or not. | W04 |
| Q3 | krilla for PDF: +3.8 MB behind feature `pdf`, rust-version 1.88 to 1.92. | W06 |
| Q4 | Push `feat/check-duplicate-marking` and `feat/clinical-checks`; locate `feat/jumping-library-metrics`. | W01 |
| Q5 | Scale gate on a CRAM slice (section 6), since no 30x BAM fits. | W05 |
| Q6 | `JAVASCRIPT_FILE` rejected on FilterSamReads and FilterVcf. | W13 |
| Q7 | docs/05 4.9 edit: a stateless cohort report over directories passed in is open. | W45 |
| Q8 | Site sets from 1000 Genomes and gnomAD, not somalier. | W24 |
| Q9 | Public downloads (Picard 3.5.0 test resources, GIAB HG002 chr20, 1000 Genomes BAMs) and oracle installs (section 6). | W35, W42 |
| Q10 | Run-folder and array tools ship best effort without real instrument files. | W35 |

## 1. Goal

The owner's words: "I want ALL modules, and more, also want daisy chaining of analysis without pipes, new output formats (including the original one) and reports, HTML, pdfs, etc".

What "all" means. Picard 3.5.0 has 121 programs (docs/04 section 1.2).

| Set | Count | Content |
|---|---:|---|
| On main, 0.17.0 | 23 | The 23 `TOOLS` rows in `crates/core/src/registry.rs` |
| In flight | 7 | MarkDuplicates, CollectSequencingArtifactMetrics, ConvertSequencingArtifactToOxoG, CollectOxoGMetrics, BamIndexStats, CheckDuplicateMarking, CollectJumpingLibraryMetrics |
| To build | 91 | The survey list: 65 docs/04 `out`, 26 `later` |

Branch state in this clone: `origin/feat/mark-duplicates`, `origin/release/0.18.0`, `origin/feat/bam-index-stats`, `origin/feat/oxog-metrics`, `origin/feat/sequencing-artifact-metrics` exist. `feat/check-duplicate-marking` and `feat/clinical-checks` are local worktrees only. `feat/jumping-library-metrics` is absent. W01 recounts 23 + 7 + 91 = 121 against docs/04 before the first tool wave (Q4). `feat/clinical-checks` (coverage gaps, sex check, docs/06) is not a Picard tool and is counted under "and more".

"And more": 15 analyses, ranked by demand in the beyond-Picard survey. doctor, mtDNA copy number, per-exon coverage (W07). Adapter content on BAM, duplicates by tile, telomere content, index-only depth (W08). FASTQ QC (W41). VCF QC, stratified benchmarking (W42). Sample identity, contamination (W43). Relatedness, cfDNA fragmentomics (W44). Haplotype-aware benchmarking (W47). MultiQC export is a format (W02). The cohort report is W45.

"Daisy chaining without pipes": every step of a chain runs in one process, records flow step to step in memory, no shell pipe, no intermediate file unless a step cannot stream (section 2).

"Output formats including the original one": Picard text stays the default for every metrics tool, including the 15 analyses (class names `bammetrics.<group>.<X>Metrics`). JSON, JSON Lines, TSV, CSV, Markdown, Parquet, XLSX, SQL, JUnit and MultiQC custom content come from one typed table (section 3). `convert` reads Picard text, JSON or Parquet and writes any of them, so every output exists in the original form. Tools that write BAM, CRAM, VCF, FASTQ or interval lists write the same file kind Picard writes. `convert` does not apply to them. `write` in a chain picks the kind by extension.

"Reports, HTML, pdfs, etc": the HTML page and index stay. PDF, Markdown, DOCX and a cohort report are added (section 3). Each new metrics tool gets a glossary row per column, a verdict rule where a threshold exists, a plot family only when a rule names it.

## 2. In-process chaining

A chain is a line of steps. It is cut into segments at barriers (sort, mark-duplicates). One segment is one reader, one fused transform thread, and a fan-out to taps (collectors) and writers. Records move as batches of at most 4096 `RawRecord` (`bam_raw.rs`, the record bytes without `block_size`). A transform-free chain lowers to `pass::run_isolated`, so `qc` keeps its speed and its golden bytes.

### 2.1 Types to add, by file

`crates/core/src/stage.rs` (new, W03):

```rust
pub trait Source: Send {
    fn open(&self, threads: usize) -> Result<(sam::Header, Box<dyn Records>)>;
    fn replayable(&self) -> bool;               // file true, stdin false
}
pub trait Records: Send { fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool>; }
pub trait Transform: Send {                     // 1 batch in, 1 batch out
    fn requires(&self) -> Order;                // checked at plan time
    fn header(&self, up: &sam::Header) -> Result<sam::Header>;
    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()>;
    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()>;
    fn counters(&self) -> Counters;             // in, out, changed
}
pub trait Barrier: Send {                       // sort, mark-duplicates phase one, merge
    fn push(&mut self, batch: &[RawRecord]) -> Result<()>;
    fn seal(self: Box<Self>) -> Result<Sealed>;
}
pub struct Sealed { pub header: sam::Header, pub next: Next }
pub enum Next { Records(Box<dyn Records>), Replay(Box<dyn Transform>) }
pub trait Sink: Send {                          // taps and writers
    fn accept_batch(&mut self, b: &[RawRecord]) -> Result<()>;
    fn finish(self: Box<Self>) -> Result<Products>;
}
pub struct Frame { pub kind: Kind, pub order: Order, pub header: sam::Header }
pub struct Segment { /* source or sealed barrier, transforms, sinks, barrier */ }
pub fn run_chain(plan: &ChainPlan, progress: &dyn Progress) -> Result<Vec<Products>>;
```

`CollectorSink` wraps the existing `Collector` (`pass.rs:55`): one reused `bam::Record` filled from the `RawRecord` bytes through `Reader<&[u8]>::read_record`, then `accept`. None of the 23 collectors change. `Worker` (`pass.rs:399`) and `Fanout` (`pass.rs:455`) are reused for the tap fan-out.

`crates/core/src/registry.rs` (W03): `Tool` gains optional fields beside `collector`: `stage: Option<fn(&Invocation, &sam::Header) -> Result<Stage>>` (`Stage::Transform` or `Stage::Barrier`), `source`, `sink`. A tool with only `run` (wgs, hs, library complexity, rna) is still chainable: the planner spools the stream to a temp BAM in `TMP_DIR`, runs `Tool::run`, deletes it, and `plan` prints that line with its disk cost (W14). A `tier: Tier` field (section 5) prints in `tools`.

`crates/core/src/mark_duplicates/` (W03): `write.rs` `Writer` becomes an `Emit` trait so a channel replaces the file. `build.rs` `run` becomes `Barrier::push`. Replay reopens the file source and reruns the streaming prefix. Copier logic is untouched. `test_mark_duplicates.py`, the fuzz run and the golden cases gate it.

`crates/core/src/spill.rs` (W05): `Spill` for `RawRecord` with a coordinate or queryname key, runs as BGZF level 1, merge keeps ties in run order. `--sort-memory` (default 2 GiB) sets the buffer.

`crates/core/src/plan.rs` (new, W05): the planner. Tracks `Frame` per position, rejects a tap whose `requires` the stream does not meet before any read, inserts nothing silently, chooses replay (file source, every earlier transform deterministic) or spool, estimates temp disk, refuses when free space is below the estimate (one `libc::statvfs` call, dependency decision in W05).

`crates/cli/src/chain.rs` (new, W04): `chain` verb. Global options first, then the source, then steps split at `+`. Each step is parsed by its own `Verb::flags` (`verbs.rs:44`), so did-you-mean works per step. `write [PATH]` is a tap; extension picks BAM, SAM, CRAM, FASTQ; `--index` writes a BAI. `-` reads stdin or writes stdout, so pipes to bwa and samtools keep working. `chain --plan plan.toml` reads the same steps from a file through `parse_tables` (moved from `crates/verdict/src/table.rs` to core in W04). `plan` prints segments, replay or spool per barrier, temp disk, and what `--resume` skips. `fifo-buffer` copies stdin to stdout and prints that chaining removed the pipe it buffered (FifoBuffer).

`crates/cli/src/qc.rs` (W14): `Plan.bam: PathBuf` becomes a chain. `Read` gains `Tap(n)` beside `Pass`, `Own`, `Header`. `collectors` yields taps for the final position. `chain ... + qc --preset P` expands to taps, verdict and report. `qc` with no transform is one segment and today's bytes.

`crates/cli/src/manifest.rs` (W14): schema `bammetrics.manifest/2`. `Key` gains `upstream`, a hash of every earlier step's key. `StepRecord::resumes` stays per step. Output kind `bam` with the existing `Fingerprint`. `/1` records stay valid.

`crates/py/src/verbs.rs` (W04): `bm.chain(source, [bm.step("sort_sam", sort_order="coordinate"), "insert_size", ...], output=, reference=, resume=)` returns the existing `Run`, GIL released as in `collect_run`. `Run.stages` holds the transform counters. `bm.step` validates against the generated registry signature.

### 2.2 Rules

- Names. Each transform has a default label (`sorted`, `marked`, `fixed`, `rg`, `cleaned`, `downsampled`, `filtered`, `merged`), `--as NAME` overrides. Kept product `PREFIX.<label>.bam`. Tap on the final stream `PREFIX.<ext>`, the same as `qc`. Tap earlier `PREFIX.<label>.<ext>`, `in` for the source. Two taps of one tool at one position is a usage error (exit 2).
- Headers. A tap after a barrier is built from the sealed header, so new `@PG` and `SO` lines are visible to it. Chained equals staged compares BAM records decoded and headers with `@PG` `CL` and `PP` normalised; every other header line must be equal. Metrics headers: the `# INPUT=` line names the chain source and the step list; a file-to-file run writes today's line.
- Failure. A failing tap fails alone, as `run_isolated` does. A failing transform or writer aborts its segment; taps in that segment are marked failed with the cause, because they saw a prefix. Finished segments keep their files. Every written product is atomic (temp name, then rename) and a failed one is removed. Exit code follows docs/05 2.4. `--resume` restarts at the last kept product whose key chain matches; a sort reruns from its source, because its runs live in `TMP_DIR` and are removed on drop.
- `STOP_AFTER` and `ASSUME_SORTED` apply per tap through `Accumulator::done`. Mid-chain `STOP_AFTER` is pinned by a staged comparison in W03.
- Threads. One `--threads` budget (default min(cores, 8)) split by the planner: BGZF inflate workers, one reader, one transform thread, one thread per tap, 4 deflate workers per writer, one sort worker. Queues bounded at 8 batches, so a segment holds 12 to 25 MB in flight. Memory is the maximum over segments, because they run in turn.
- Disk at 30x (derived): about 620 M records, about 220 GB uncompressed, a 70 GB BAM. A spool near BAM size does not fit on this machine. The planner prefers replay, `plan` prints the estimate, `plan` refuses when free space is short. A kept product may be CRAM (`write out.cram`) to halve the size.
- Variant streams (W46) use the same executor with `vcf::Variant` as the item. Cross-kind steps (FASTQ to BAM, BAM to FASTQ) are sources and sinks, never transforms.

### 2.3 Increments

| Wave | Content | Gate |
|---|---|---|
| W03 | stage.rs, executor, MarkDuplicates barrier, tap adapter, `write` | Chained equals staged on corpus t1; transform-free chain within 5 percent of `qc` on chr21.bam, else transform-free chains keep `run_isolated` (kill rule) |
| W04 | `chain`, `+`, TOML plan, Python, FifoBuffer | Both documented examples run; `plan` prints segments and disk |
| W05 | sort barrier, planner, replay or spool, memgate | Sort chain RSS under `--sort-memory` plus 1 GiB on the scale slice |
| W14 | manifest 2, resume, kept products, `qc` as terminal step, materialize fallback | Edit step 3 reruns 3 and later only; `qc` golden unchanged |
| W46 | VCF kind, streaming pileup and capture sinks | Wgs and hs mid-chain with no temp BAM, equal to staged |

Each BAM-writing tool from W09 on registers its `stage` on the day it lands.

## 3. Output and report formats

### 3.1 Model

`Table { class, columns: Vec<(name, Ty)>, rows }` in a new `crates/formats` (`bammetrics-formats`, std only by default). `Ty` is Int, Float or Str from `glossary::lookup(class, column)` (487 columns: 299 float, 144 int, 44 string). Reason: on the NA12878 fixture 63 of 165 non-empty cells infer a different type from their text than the glossary declares (MEDIAN_INSERT_SIZE prints 310, so `parse_value` makes it Int; Picard's type is double). Without glossary types, polars frames from two files have different dtypes. A class the glossary lacks falls back to inference, and a test fails when any registered tool emits a column with no glossary row.

Rows stay rows (ALL_READS, SAMPLE, LIBRARY, READ_GROUP). Histograms are a second table with the key column and one column per series. StringHeader lines travel as `headers[]`. Structured formats carry full f64; Picard text prints at most 6 fraction digits. Lossless set: Picard text, Parquet, Arrow (empty is null, `?` is NaN). Lossy: JSON, CSV, SQLite, XLSX (null covers both). A golden round trip, Picard text to Parquet to Picard text, runs on every fixture.

### 3.2 Format table

| Format | Wave | Cost | Note |
|---|---|---|---|
| Picard text | exists | 0 | Default. Byte-identical after the StringHeader lines. Every new tool writes it |
| JSON | W02 | 0 | Additive: `schema: "bammetrics.metrics/1"`, `headers[]`, column types; `docs/schema/metrics-1.json`; golden JSON records are a reviewed diff |
| JSON Lines | W02 | 0 | One record per row, a `histogram` discriminator for histogram rows |
| TSV | W02 | 0 | Fix: escape tab, newline, backslash in string cells |
| CSV | W02 | under 5 KB code | RFC 4180, hand-written, no crate |
| Markdown | W02 | under 5 KB | GFM tables per class, verdict on top, 65,536 character cap |
| JUnit XML | W02 | under 3 KB | `check`: one testcase per rule |
| MultiQC custom content | W02 | under 3 KB | `PREFIX.verdict_mqc.json`. Picard-class files already parse: MultiQC 1.35 output on 8 tools equals Picard's |
| `convert` verb | W02 | small | Picard text, JSON or Parquet in; any format out |
| `--format a,b`, `--also` on `qc` and `chain` | W02 | small | Sibling files under the prefix, listed in the manifest |
| Parquet | W21 | +0.53 MB, feature `parquet` on in release | Low-level `SerializedFileWriter`, snappy, one file per class plus a long file (sample, class, row_key, column, value) |
| XLSX | W21 | 0 crates | Hand-written OOXML in a stored zip, one sheet per class; openpyxl reads it |
| DOCX | W21 | 0 crates | Tables and text, no plots; pandoc and python-docx read it |
| SQL dump | W21 | 0 crates | `sqlite3 < file.sql`; native `.sqlite` behind feature `sqlite` (+1.55 MB) in W46 |
| Arrow IPC, `to_polars`, `to_pandas` | W46 | +0.46 MB, wheel only | `read_table(path)` with glossary dtypes ships in W02 at no Rust cost |

### 3.3 PDF route

krilla 0.8.2 with krilla-svg in a new `crates/pdf`, cargo feature `pdf`. Evidence from stripped release probes (lto fat, codegen-units 1, the method of `crates/cli/tests/budgets.rs`); base probe 633,760 bytes:

| Route | Growth | Crates | Time | PDF | Plots |
|---|---|---:|---|---|---|
| krilla, no SVG | +1.97 MB | 53 | 17 ms | 35 KB | by hand |
| krilla + krilla-svg (chosen) | +3.81 MB | 97 | 25 ms | 47 KB, 28 KB with 70 KB subset fonts | vector, text selectable |
| svg2pdf + pdf-writer | +3.06 MB | 93 | n/a | 55 KB | vector, no font subsetting or tagging |
| typst 0.15.1 as a library | +29.0 MB | 346 | 20 to 50 ms | 71 KB | good, MPL-2.0 crates in tree |
| Chrome headless | 0 | 0 | 2.0 s warm | 366 KB | exact, needs a browser, bytes differ run to run |
| wkhtmltopdf | 0 | 0 | 12 s | n/a | plots lost (no CSS variables) |
| printpdf html | n/a | 181 | 2 s | 41 pages | plots lost as text |

krilla is pure Rust, byte-identical across runs (needed by the golden suite), 8.5 MB peak RSS, validates PDF/A and PDF/UA-1, tags content, writes outlines. Costs: about 1,000 lines of layout code (wrap, tables, page breaks, tags), three Latin subsets of DejaVu (70 KB), `--pdf-font PATH` for non-Latin sample names, rust-version 1.92 (Q3). One change in `crates/report`: `plots::figure_standalone` replaces `var(--x)` with the print palette (`theme::token(Theme::Print, name)`), because usvg ignores CSS variables.

Kill rule (W06): layout code over 1,500 lines or growth over 4 MB keeps PDF as a documented Chrome recipe with feature `pdf` off. Budget rows added to docs/05 4.8: PDF growth under 4 MB, render under 100 ms per sample, under 150 KB per sample, each failing at 2x. The 300 KB report row keeps its probes.

PDF content follows docs/05 4.2: header, verdict band, failing checks each with its plot, check table, provenance, an appendix with every table (the print CSS hides `#tables` today), an outline entry per section. A later krilla-only backend for `svg.rs` `Chart` (about 12 primitives) saves about 1.8 MB (estimate).

### 3.4 Reports

HTML: unchanged. Markdown and plain text: from the same `Verdict`. Cohort (W45, Q7): built only from the directories passed in. HTML adds a strip plot per key metric, MAD outliers, a sample by rule heat table, inline fragments for non-passing samples only; a 100-sample file stays under 1.5 MB. PDF adds page 1 the sorted index, page 2 the distribution plots, a section per non-passing sample, bookmarks per sample.

Coverage for new tools: each of the 15 analyses gets glossary rows, a KNOWN_DIFFERENCES section naming its oracle tier, and a rule in the verdict language where the survey names one (`exon.min-pct20`, `contamination.freemix`, `identity.lod`, `vcf.titv`, `mtcn`). A plot family is added only when a rule names it (coverage for exon and gaps, insert-size for fragmentomics).

## 4. Shared infrastructure

Each piece is built once, in the wave that first needs it, with a proving tool from its own group. A foundation with no consumer in its wave is not merged.

| Piece | Wave | Lines | Tools unblocked | Proving tool |
|---|---|---:|---:|---|
| bam-writer, spill-sort, md-nm | W00 (in flight) | 0 | 27 | MarkDuplicates |
| Table model, format writers, `convert` | W02 | 1,200 | all metrics tools | every registered tool |
| compare_bam, deflater spike, Picard queue, memgate, spec agent | W02, W05 | harness | 40 BAM, CRAM, gz tools | MarkDuplicates output |
| stage.rs, executor, planner | W03, W05 | 1,500 | all | MarkDuplicates, SortSam |
| sam-comparators, sam-text-reader and writer, header-merge, bai-reader and writer | W05 | 950 | 10 | SortSam, ViewSam, BuildBamIndex |
| crates/pdf | W06 | 1,000 | reports | fixture pages |
| GTF and BED12 reader, into_target_depths | W07 | 400 | per-exon, gaps | ExonCoverageMetrics |
| interval-ops (union, intersect, subtract, pad, scatter) | W12 | 500 | 5 | IntervalListTools |
| java-random reuse, murmur3 | W13 | 150 | 2 | DownsampleSam |
| cram-writer | W15 | 400 | 3 | SamFormatConverter |
| vcf-writer, VCF header inventory, Variant stream kind | W16 | 700 | 18 | SplitVcfs |
| vcf-sort-merge, ped-reader | W18 | 500 | 5 | SortVcf |
| chain-reader, fasta-writer, gff-io | W19 | 650 | 7 | LiftoverVcf |
| bcf-codec (noodles `bcf` feature) | W20 | 1,000 | 1 plus BCF input for all VCF tools | VcfFormatConverter |
| Parquet writer, OOXML zip | W21 | 600 | formats | fixtures |
| fastq-io, read-clipping, bfq-writer | W22 | 700 | 5 | SamToFastq |
| adapter-clip (ClippingUtility) | W23 | 300 | 2 | MarkIlluminaAdapters |
| mate-cigar-engine, umi-graph | W24, W25 | 1,300 | 4 | MarkDuplicatesWithMateCigar |
| haplotype-map, fingerprint-model | W26 | 1,900 | 9 | CheckFingerprint |
| mersenne-twister, chi-squared, regularized_gamma | W28 | 400 | 3 | CalculateFingerprintMetrics |
| sam-validator, lenient-bam-reader, raw-bytes BAM writer | W29 | 3,000 | 1 | ValidateSamFile |
| bcl-reader (BCL, CBCL, filter, locs, clocs, RunInfo.xml), interop-reader | W35 | 1,200 | 6 | CheckIlluminaDirectory |
| barcode-matcher, basecalls-converter | W36, W37 | 1,150 | 4 | ExtractIlluminaBarcodes |
| bpm-egt-reader, gtc-reader, adpc-writer | W38, W39 | 1,000 | 5 | CompareGtcFiles |
| site sets (GRCh37, GRCh38), FASTQ to Record adapter | W43, W41 | data | 4 analyses | IdentityMetrics |

Java behavior stays in `java.rs` and `commons_math.rs`. Already on main and reused: `cram.rs` (read), `adapter.rs`, `java.rs`, `optical.rs`, `md5.rs`, `gzip.rs`, `intervals.rs`, `levels.rs`, `locus.rs`, `vcf.rs`.

Registration hot spots. Every tool edits `TOOLS` in `registry.rs`, `VERBS` in `verbs.rs`, `glossary.tsv`, `bammetrics.pyi` and `tests/golden`. Rule from W02 on: builders write the tool file and its tests only, plus a `registration.md` snippet in the PR. The wave integrator applies the five-file edits in branch order, regenerates the stub with `UPDATE_PYI=1`, and runs `make golden-update` once. Builders never touch those five files. Before a wave starts, the integrator runs `bammetrics tools` on main and lists `origin/feat/*` to confirm no other workflow has built a listed tool.

## 5. Wave table

Wave template, as 0.17.0 ran: per branch one builder (own worktree), two lenses (A oracle cases, B surface: options, stubs, glossary, formats, chain stage, docs), one fixer. Per wave one spec agent (Haiku: `picard X --help` plus Picard 3.5.0 source notes per tool, before builders start), one integrator (holds the Picard queue slot), one reviewer (Opus on W03, W05, W23, W26 to W33, else Sonnet). Agents per wave = 4 x branches + 2 (the spec agent replaces the plan-1 docs count). Three branches 14, two branches 10, one branch 6. Peak concurrency 6.

Oracle tiers. A: Picard 3.5.0 output byte-identical after the StringHeader lines. B: Picard output decoded record by record (BAM, CRAM, FASTQ.gz, BCF, VCF with a header rule per tool), `@PG` normalised. C: Picard on synthetic data with a documented tolerance. D: no real-data oracle, round trip and parser-derived tests only, `best effort` in `tools`. E: beyond Picard, planted truth then a named tool within a stated tolerance. The W02 deflater spike (Picard with `USE_JDK_DEFLATER=true COMPRESSION_LEVEL=5` against our writer, `cmp`) decides once whether BAM tools are A or B.

Base gate, every wave: fmt, clippy `-D warnings`, `cargo test --workspace`, `make golden` changed only by added or reviewed records, concordance at the tool's tier with 3 distinct non-default values per column, a glossary row per emitted column, the format round trip, chained equals staged for every tool with a stage, `make python-test`, `make wheel-check`, `make budgets`, `make memgate` for every writer, version bump, KNOWN_DIFFERENCES.md and CLAUDE.md updated, one PR per branch, none to main, no AI reference in any commit.

| Wave | Version | Title | Items | Needs | Agents | Oracle and data | Acceptance beyond the base gate |
|---|---|---|---|---|---:|---|---|
| W00 | 0.18.0 | Release MarkDuplicates | review and tag release/0.18.0 | branch | 2 | A metrics, B BAM; markdup corpus, fuzz | 0.18.0 tagged, `test_mark_duplicates.py` and fuzz pass |
| W01 | 0.19.0 | Merge train | bam-index-stats, oxog-metrics, sequencing-artifact-metrics, then the pushed check-duplicate-marking, jumping-library-metrics, clinical-checks (Q4); tool recount; docs/04 rows for the 7 | W00 | 6 | each branch's own suite rerun on the merged head | 30 tools in `tools`, one golden rewrite reviewed, clinical-panel preset carries the docs/06 rules |
| W02 | 0.20.0 | Table model, formats core, harness | A: `crates/formats`, Table, CSV, Markdown, TSV escape, JSON schema, JSONL, JUnit, verdict_mqc. B: `--format` list, `--also`, `convert`, `read_table`. C: compare_bam, deflater spike, Picard queue through `agent-job`, cargo-deny, spec-agent script, registration rule | W01 | 14 | round trip on `tests/data/picard`; MultiQC 1.35 parses the mqc file | `convert` round trip identical on every fixture; compare_bam fails on a planted difference; tier A or B decided for BAM; JUnit renders in a CI viewer |
| W03 | 0.21.0 | Stage executor, MarkDuplicates barrier | A: stage.rs, executor. B: barrier plus replay, `Emit`. C: tap adapter, `write` sink | W00 | 14 | chained equals staged on corpus t1 | within 5 percent of `qc` on chr21.bam or the kill rule applies; mid-chain STOP_AFTER pinned |
| W04 | 0.22.0 | `chain` verb, plan file, Python | A: `chain`, `+`, did-you-mean, `plan`. B: `--plan` TOML, `parse_tables` to core. C: `bm.chain`, `bm.step`, `Run.stages`, FifoBuffer | W03, W02 | 14 | self; golden cases for chain text and errors | both documented examples run; `plan` prints segments and disk; FifoBuffer exits 0 |
| W05 | 0.23.0 | Record IO, sort barrier, planner, memgate | A: SAM text reader and writer, comparators with a total-order test. B: header-merge, Spill for RawRecord, planner, free-space refusal. C: bai reader and writer, `write --index`, `make memgate`. Tools SortSam, ViewSam, BuildBamIndex | W04 | 14 | A for BAI and SAM text, B for BAM; shuffled corpus; htsjdk tie order on random keys | BAI byte-identical; sort chain on the scale slice (section 6) with RSS under `--sort-memory` plus 1 GiB, spill under 1.2x input |
| W06 | 0.24.0 | PDF report | A: `figure_standalone`. B: `crates/pdf`, fonts, `--pdf-font`. C: `report --pdf`, `qc --report pdf`, PDF/UA-1 test, budget rows | W02 | 14 | self; veraPDF or a structure test; two runs `cmp` | byte-stable, 3 pages for the fixture, growth under 4 MB, kill rule checked first |
| W07 | 0.25.0 | Beyond A | doctor (header, last BGZF block, first 10,000 records, build detection, `--sample-sheet`), MtCopyNumberMetrics, ExonCoverageMetrics with GTF and BED12 reader | W01, W02 | 14 | E; per-exon mean and min tier A against `CollectHsMetrics PER_TARGET_COVERAGE`; planted defects; samtools coverage | doctor under 1 s on the scale CRAM header; exit codes as `check` |
| W08 | 0.26.0 | Beyond B | AdapterContentMetrics, DuplicateByTileMetrics, TelomereMetrics, index-only depth | W01, W02 | 14 | E; planted reads; samtools idxstats | 4 analyses match planted truth; tolerances written |
| W09 | 0.27.0 | BAM edit transforms 1 | AddCommentsToBam, AddOATag, CleanSam, SetNmMdAndUqTags, SetNmAndUqTags, ReplaceSamHeader | W05 | 14 | A or B per W02; corpus with tags stripped; reads off contig ends | 6 tools equal; each runs mid-chain |
| W10 | 0.28.0 | BAM edit transforms 2, split sinks | AddOrReplaceReadGroups, RevertOriginalBaseQualitiesAndAddMateCigar, SplitSamByLibrary, SplitSamByNumberOfReads, GatherBamFiles, CalculateReadGroupChecksum | W09 | 14 | B; GatherBamFiles A by block copy; corpus has 2 libraries, 3 read groups | fan-out sinks in the executor; checksum equal on any header |
| W11 | 0.29.0 | Metrics tools, no new infrastructure | CollectQualityYieldMetricsSNVQ, CollectUmiPrevalenceMetrics, CollectRrbsMetrics, CreateBafRegressMetricsFile, CreateVerifyIDIntensityContaminationMetricsFile, CollectArraysVariantCallingMetrics | W02 | 14 | A; qa qc qg qt tags, RX and BQ tags, C to T reads added to the corpus; hand-written text inputs; synthetic array VCF | 6 tools equal; three have glossary rows from scratch |
| W12 | 0.30.0 | Interval foundation | interval-ops; IntervalListToBed, IntervalListTools, ScatterIntervalsByNs, VcfToIntervalList, NonNFastaSize | W05 | 14 | A; Picard-made lists; reference with N runs; scatter tree compared as a set | every `ACTION` of IntervalListTools equal; `qc --panel` takes their output |
| W13 | 0.31.0 | Sort, merge, filter, downsample | MergeSamFiles, FilterSamReads (`JAVASCRIPT_FILE` rejected, Q6), FixMateInformation, ReorderSam, DownsampleSam, PositionBasedDownsampleSam | W05, W12 | 14 | B; corpus split by RG and merged with id collisions; reordered reference; `RANDOM_SEED`, three strategies, 50 seeds | seeded records equal; merge RSS under the W05 numbers |
| W14 | 0.32.0 | Chain 2 | manifest 2, per-step resume, kept products, `qc` as terminal step, materialize fallback, stage on W09 to W13 tools, schema freeze (section 6) | W13 | 14 | self; `qc` golden unchanged | edit step 3 reruns 3 and later only; every registered tool runs at any position |
| W15 | 0.33.0 | CRAM writer, convert, compare, terminator | cram-writer; SamFormatConverter, CompareSAMs, CheckTerminatorBlock | W13 | 10 | B; corpus both ways; damaged BAM tails | six conversions decoded equal; exit 100 on a damaged tail |
| W16 | 0.34.0 | VCF writer, variant stream kind | Haiku inventory: all 14 Picard VCF writers run once, added header lines recorded in `tests/concordance/vcf_headers.json`; vcf-writer; Variant kind; SplitVcfs, MakeSitesOnlyVcf, MakeVcfSampleNameMap | W14 | 14 | B with the header rule, A for the TSV; `variant_vcfs.py` | header lines match the inventory; variant stream runs in `chain` |
| W17 | 0.35.0 | VCF transforms | RenameSampleInVcf, UpdateVcfSequenceDictionary, FixVcfHeader, FilterVcf (`JAVASCRIPT_FILE` rejected), GatherVcfs | W16 | 14 | B; `calls3.vcf` split by contig | 5 tools equal apart from inventoried lines |
| W18 | 0.36.0 | VCF sort and merge, PED | vcf-sort-merge, ped-reader; SortVcf, MergeVcfs, FindMendelianViolations, CombineGenotypingArrayVcfs, MergePedIntoVcf | W17 | 14 | B; synthetic trio and PED; hand-written PED, MAP, zCall | SortVcf on 10 M records under 1 GiB; the 3.5.0 `VCF_DIR` bug reproduced or rejected as docs/04 says |
| W19 | 0.37.0 | Chain files, FASTA, GFF | chain-reader, fasta-writer, gff-io; LiftOverIntervalList, LiftoverVcf, ExtractSequences, NormalizeFasta, SortGff | W18, W12 | 14 | A, B; synthetic chain with strand flips on two references; synthetic GFF3 | liftover equal including the reject file |
| W20 | 0.38.0 | BCF | bcf-codec, BCF input for every VCF tool; VcfFormatConverter | W17 | 6 | B after BGZF decompression | BCF read by every VCF tool; both directions equal |
| W21 | 0.39.0 | Formats 2 | Parquet and the long file, XLSX, DOCX, SQL dump | W02, W06 | 14 | round trip; DuckDB, polars, openpyxl, pandoc | Parquet types equal the glossary; text to Parquet to text identical |
| W22 | 0.40.0 | FASTQ io, uBAM conversion | fastq-io, read-clipping, bfq-writer, FASTQ comparator; FastqToSam, SamToFastq, SamToFastqWithTags, RevertSam, BamToBfq | W13 | 14 | A for FASTQ text and bfq, B for gz and BAM; FASTQ cut from the corpus; `OUTPUT_BY_READGROUP` | 5 tools equal; FASTQ reader feeds the four read-level collectors |
| W23 | 0.41.0 | Adapter clipping, MergeBamAlignment | adapter-clip; MarkIlluminaAdapters; MergeBamAlignment (34 options, two branches) | W22 | 10 | B; uBAM from RevertSam, aligned BAM from bwa on the corpus; adapters inserted by script | every option equal on 3 aligned inputs; XT tags equal |
| W24 | 0.42.0 | Mate-CIGAR engine | mate-cigar-engine; MarkDuplicatesWithMateCigar, SimpleMarkDuplicatesWithMateCigar; site-set licence spike (Q8) | W13 | 10 | B and A for metrics; markdup BAMs with MC | equal on the markdup corpus and `make fuzz SEED=1 N=200` |
| W25 | 0.43.0 | UMI graph, replicate metrics | umi-graph in Java hash order; UmiAwareMarkDuplicatesWithMateCigar, CollectIndependentReplicateMetrics | W24 | 10 | B, A; UMI families at edit distance 1; replicate sets of size 2 and 3 hand-built (tier D for the alpha tool) | UMI inference equal; replicate tool labelled if two Picard runs differ |
| W26 | 0.44.0 | Haplotype map, fingerprint model | haplotype-map, fingerprint-model, synthetic map generator; CheckFingerprint, ConvertHaplotypeDatabaseToVcf, LiftOverHaplotypeMap | W19, W20 | 14 | C; synthetic map on the corpus reference, synthetic people sampled to VCF and BAM | 100 seeded cases, LOD columns equal at 6 digits, log10 summation order ported |
| W27 | 0.45.0 | Crosscheck family | CrosscheckFingerprints, CrosscheckReadGroupFingerprints, ClusterCrosscheckMetrics | W26 | 10 | C; `NUM_THREADS=1` exact, higher at 1e-5 relative; hash-order rows by `java.rs` | row order equal under HashMap order for 3 library counts |
| W28 | 0.46.0 | Seeded fingerprint tools | mersenne-twister, chi-squared; ExtractFingerprint, IdentifyContaminant, CalculateFingerprintMetrics | W26 | 14 | C; mixes at a known rate; 50 seeds | seeded sampling equal; VCF equal apart from header lines |
| W29 | 0.47.0 | ValidateSamFile 1 | lenient-bam-reader, rule framework, raw-bytes BAM writer (one broken record per rule); record, flag, CIGAR and mate rules | W15 | 14 | A; one Picard run per broken file | each rule's message equal |
| W30 | 0.48.0 | ValidateSamFile 2 | header, RG, PG, sort and tag rules | W29 | 14 | A | each rule's message equal |
| W31 | 0.49.0 | ValidateSamFile 3 | CRAM and index rules, `MODE`, `IGNORE`, `MAX_OUTPUT`, `INDEX_VALIDATION_STRINGENCY`, SUMMARY and VERBOSE | W30 | 14 | A | about 70 rules covered; rules 3.5.0 cannot trigger are listed |
| W32 | 0.50.0 | CollectSamErrorMetrics 1 | pileup sink, error core, 5 stratifiers | W17 | 10 | A; corpus plus `variant_vcfs.py` sites | 5 stratifier files equal |
| W33 | 0.51.0 | CollectSamErrorMetrics 2 | remaining stratifiers, `FILE_EXTENSION` outputs, hash order | W32 | 10 | A | every stratifier equal |
| W34 | 0.52.0 | BaitDesigner, Flow quality yield | BaitDesigner (two Picard runs first to confirm determinism); CollectQualityYieldMetricsFlow | W12 | 10 | A; soft-masked reference; D for Flow (no public Ultima BAM) | BaitDesigner equal; Flow labelled best effort |
| W35 | 0.53.0 | BCL reader, run-folder checks | data spike (Q9): Picard 3.5.0 `testdata/picard/illumina`, nf-core bcl test data; bcl-reader, interop-reader; CheckIlluminaDirectory, CollectIlluminaLaneMetrics, CollectHiSeqXPfFailMetrics; synthetic run-folder writer by a second agent from public specs | W11 | 14 | A on fetched fixtures, D otherwise | tier A where fixtures exist; CBCL without a fixture best effort |
| W36 | 0.54.0 | Barcode matcher | CollectIlluminaBasecallingMetrics, ExtractIlluminaBarcodes | W35 | 10 | as W35 | equal on every available folder |
| W37 | 0.55.0 | Basecalls converter | IlluminaBasecallsToFastq, IlluminaBasecallsToSam | W36, W22 | 10 | as W35; gz after decompression | 96-barcode run under 2 GiB RSS |
| W38 | 0.56.0 | Array file readers | data spike (`testdata/picard/arrays`); bpm-egt-reader, gtc-reader; BpmToNormalizationManifestCsv, CompareGtcFiles, CreateExtendedIlluminaManifest | W19 | 14 | A on fixtures, D otherwise | round trip passes; `tools` shows the label |
| W39 | 0.57.0 | Array VCF and adpc | adpc-writer; GtcToVcf, VcfToAdpc | W38 | 10 | A for adpc bytes, B for VCF | adpc `cmp` equal |
| W40 | 1.0.0 | Parity gate | 121-tool registry check, fuzz N=500 on BAM writers, docs/04 and docs/05 edits (section 7), README tool table with tier, KNOWN_DIFFERENCES tier table, PyPI release | W39 | 6 | all | 121 registered, 0 `out`, no `both-err` rise |
| W41 | 1.1.0 | FASTQ QC | FastqMetrics, overrepresented sequences (cap 100,000 prefixes), adapter content by cycle, per-tile quality, preset `fastq` | W22 | 10 | A through FastqToSam for the four read-level tools; E against FastQC and fastp | four collectors byte-equal; throughput published before any claim |
| W42 | 1.2.0 | VCF QC, stratified benchmarking | VariantQcMetrics, StratifiedConcordanceMetrics (size bins, left-align, no haplotype matching) | W20 | 10 | E; bcftools stats exact on counts; hap.py or vcfeval within 0.1 percent on HG002 chr20 (Q9) | counts exact; F1 within 0.1 percent |
| W43 | 1.3.0 | Identity and contamination | site sets (Q8), IdentityMetrics and SiteCounts sketch, ContaminationMetrics | W28 | 14 | E; `samtools view -s` mixes at 0, 1, 2, 5, 10 percent; VerifyBamID2 | within 0.005 of planted at 30x; INCOMPLETE under 5,000 sites |
| W44 | 1.4.0 | Relatedness, cfDNA fragmentomics | RelatednessMetrics over explicit sketch paths, FragmentMetrics, end-motif file | W43 | 10 | E; somalier within 0.02; insert histogram tier A by reuse | explicit paths only, nothing discovered or stored |
| W45 | 1.5.0 | Cohort report | cohort HTML and PDF | W21, Q7 | 10 | self; 100 synthetic samples | HTML under 1.5 MB; PDF bookmarks per sample |
| W46 | 1.6.0 | Streaming sinks, Python frames, SQLite | VCF stream kind, `locus.rs` stream sink for wgs and hs, Arrow IPC in the wheel, `to_polars`, `to_pandas`, feature `sqlite` | W21, W18 | 14 | chained equals staged | wgs and hs mid-chain with no temp BAM; `to_polars` dtypes equal the glossary |
| W47 | 1.7.0 | Haplotype-aware benchmarking | variant-path matching for stratified concordance | W42 | 10 | E; vcfeval on HG002 chr20 | within 0.1 percent of vcfeval |

Totals: 48 waves, 572 agents, 91 Picard tools (W04 1, W05 3, W09 6, W10 6, W11 6, W12 5, W13 6, W15 3, W16 3, W17 5, W18 5, W19 5, W20 1, W22 5, W23 2, W24 2, W25 2, W26 3, W27 3, W28 3, W29 to W31 1, W32 to W33 1, W34 2, W35 3, W36 2, W37 2, W38 3, W39 2), 15 analyses. Inside W09 to W13 and W22 the order follows nf-core sarek and GATK usage: SortSam, MergeSamFiles, FixMateInformation, AddOrReplaceReadGroups, SetNmMdAndUqTags, FastqToSam, SamToFastq, RevertSam before the rest.

Budget cuts if Q1 lands below 572: first the 12 best-effort tools (W34 Flow, W35 to W39, 70 agents), then W47 and W45 (20), then the fingerprint family W26 to W28 (38). Nothing in W02 to W14 is cut.

## 6. Machine and agent schedule

Machine: Mac Studio M1 Ultra, 20 cores, 64 GB, 77 GB free of 1.8 TB at plan time (15 percent used, the rest is on other volumes). `.claude/worktrees` holds about 32 GB.

- Agents: at most 8 at once, 14 per workflow (Q1). Builders run `cargo test`, `cargo clippy` and `make golden` only. Concordance, fuzz, memgate and scale runs happen once per wave in the integrator slot.
- Picard: one JVM at a time. The slot is the existing `agent-job --queue bammetrics "<task>"` (pueue group `agents`, parallelism 1); W02 makes `run_picard.py` refuse to start while another `picard` process is listed by `/bin/ps`. Four other workflows run now, so a wave waits for the queue. `JAVA_OPTS=-Xmx1g` for `--help` and metrics tools, `-Xmx4g` for Hs, WGS and fingerprint runs, `-Xmx8g` for MarkDuplicates, MergeBamAlignment and ValidateSamFile on the exome. Call `picard`, never bare `java`.
- JVM hours (opinion): a tool wave runs about 90 Picard cases at 10 to 40 s each, under 1 hour with the `target/fuzz/cache` hit rate above 50 percent on a rerun. W29 to W31 run about 70 one-record files each, under 1 hour. W23 and W27 run under 2 hours.
- Disk: floor 25 GB free. Before each wave `git worktree remove` for merged worktrees and `cargo clean` in the rest; a wave of 3 worktrees uses about 12 GB of targets plus 6 GB of corpora. Picard cache pruned to 5 GB, `TMP_DIR` for spill tests at most 20 GB, deleted at wave end.
- Scale data: `~/data/bammetrics-scale/genome` holds `NA12878.final.cram` (16 GB), `GRCh38.fa` (3.3 GB) and their indexes; the exome BAM (10 GB) and two markdup BAMs (10 GB each) sit beside it. No 30x BAM exists and one (60 to 70 GB) does not fit. Scale gate (Q5): W05 writes a chr20 plus chr21 BAM slice from the CRAM with `samtools view` (about 4 GB), runs each writer on it and on the exome, records peak RSS, wall time and spill bytes with `/usr/bin/time -l` into docs/benchmarks.md (`make memgate`), and states the 30x number as a linear extrapolation. Targets: streaming writer RSS at most 1.0 GiB, sort or merge RSS at most `--sort-memory` plus 1 GiB, spill at most 1.2x input, VCF writer at most 500 MiB. Scale runs happen at W05, W14, W22 and W40 only, alone, after the wave merges. No other workflow runs on `~/data` meanwhile.
- Oracle installs (Q9): picard, samtools, bwa and bcftools exist. mosdepth, somalier, VerifyBamID2, FastQC, fastp, hap.py or rtg vcfeval, veraPDF, DuckDB and polars do not. W06 needs veraPDF (JVM) or falls back to a structure test plus qpdf. W07 needs mosdepth. W41 to W44 need the rest. The W06 and W40 spec agents list the installs; the owner approves them as one item.
- Release shape: one minor per wave. PyPI and bioconda publish at the end of each phase (0.24.0, 0.33.0, 0.39.0, 0.46.0, 0.57.0, 1.0.0), so the stub changes six times before 1.0.0. Stability rule from W14: the plan-file keys, `bammetrics.metrics/1`, `bammetrics.manifest/2` and the Python `chain` and `step` signatures change only additively until 1.0.0.
- Python surface: 91 functions and 15 analyses added. Budgets from W14: wheel under 25 MB with features `pdf` and `parquet`, import under 300 ms, `UPDATE_PYI` owned by the integrator.
- Rules that hold: no `/dev/tcp` probes, an unreachable host is reported. No patient data enters the repo. Long unattended waves go through `agent-job`.
- Wall time (opinion): one wave per working day when the Picard queue is free, so about 10 weeks; W29 to W33 and W35 to W39 are the likely slips.

## 7. Changes in docs/04 and docs/05

Each change is one reviewed edit in the wave named.

| # | Place | Today | Change | Wave |
|---|---|---|---|---|
| D1 | docs/04 1.1 counts, 1.2 table | 19 core, 37 later, 65 out | The 65 `out` rows and 26 unbuilt `later` rows are scheduled: decision `build (docs/07 Wnn)`, covering tool kept as a note | W01 (pointer), W40 (row by row) |
| D2 | docs/04 after 1.2, rule 1 | Harness-prep tools never built as Rust ports | Reversed. SortSam, BuildBamIndex, FixMateInformation, SetNmMdAndUqTags, AddOrReplaceReadGroups, CleanSam, DownsampleSam, FastqToSam are built in W05, W09, W13, W22 | W05 |
| D3 | docs/04 after 1.2, rule 2 | Out-of-scope tools get one README line | README table of 121 tools with wave, tier and covering tool | W40 |
| D4 | docs/04 R13, Flow and Arrays "verify only if the decision changes" | out | Built in W34 and W11 | W11 |
| D5 | docs/04 5.5 and 5.9 claim levels | A to C2 | Tier letters A to E added; tier D text for the 12 best-effort tools | W02 |
| D6 | docs/04 R3 | Picard source unread for three tools | Each wave's spec agent reads the 3.5.0 source for every tool first | W02 |
| D7 | docs/05 1 non-goals | "new metric tools beyond the bed-to-interval-list step" | Removed. "A workflow engine" stays: a chain is a line of steps with taps, no branching, loops, scatter or cluster | W04 |
| D8 | docs/05 2.1 dispatch | Picard class, grammar verbs, clap verbs | `chain`, `convert`, `doctor`, `fifo-buffer` added | W02, W04, W07 |
| D9 | docs/05 3 layers | Layer 3 is the fixed `qc` graph | Layers 1 (pipes) and 2 (prefix files) stay. Layer 3 becomes a preset of `chain`. Example 9 becomes one `chain` command. Stage table gains the `stage` column | W04, W14 |
| D10 | docs/05 4.8 budgets | 300 KB report growth | Row kept. New rows: PDF growth under 4 MB, render under 100 ms, under 150 KB per sample, fail at 2x; data assets: identity set 70 KB embedded, contamination set 400 KB beside the binary | W06, W43 |
| D11 | docs/05 4.9 | "cohort and lab-wide views" hosted | See the boundary below (Q7) | W45 |
| D12 | docs/05 8 decision 4 | hand-written TOML reader in verdict | Holds; `parse_tables` moves to core | W04 |
| D13 | docs/05 8 decision 6 | no sample sheet runner | Holds; `doctor --sample-sheet` checks one sample per call | W07 |
| D14 | docs/05 8 decision 10 | MarkDuplicates next after 0.17.0 | Closed by W00 | W00 |
| D15 | CLAUDE.md "no new dependency" for verdict and report | holds | Holds. krilla, parquet, rusqlite, arrow live in `crates/pdf` and `crates/formats` behind features; rust-version 1.92 for `pdf` (Q3); cargo-deny in CI | W02, W06 |
| D16 | docs/04 line on JS filters | n/a | `JAVASCRIPT_FILE` on FilterSamReads and FilterVcf is `Opt::rejected` (Q6); both listed as not fully concordant | W13, W17 |
| D17 | ROADMAP.md Phase 0 | bammetrics done in 0.2.0 | Pointer to docs/07 and the 121-tool goal | this PR |
| D18 | KNOWN_DIFFERENCES.md | per tool | A tier table and one section per tier B to E tool and per analysis | every wave |

The docs/05 4.9 boundary. Stays open: a verdict for one run; an index or cohort page over the directories passed on one command line, one invocation, no stored state; identity and relatedness over sketch paths given explicitly. Stays hosted: stored cohorts, runs over time, cross-run comparison, baselines, audit trail, users, comments, uploads. The report and verdict crates keep their rules: no read of other run directories, no history store, no user identity, no network call, no comment field. A cohort-trained classifier score (cfDNA) and index-depth cohort normalisation are hosted.

## 8. Risks and the wave that retires each

| Risk | Retired by | Proof or fallback |
|---|---|---|
| Picard source unread; in and out types, determinism and sizes come from `--help` and memory | every wave, spec agent | One `--help` plus source note per tool before the builder starts; two Picard runs confirm determinism for BaitDesigner (W34) and the replicate tool (W25) |
| BAM bytes differ from Picard (Intel deflater claim, unverified on this Mac) | W02 | Deflater spike with `USE_JDK_DEFLATER=true COMPRESSION_LEVEL=5`; tier A or B fixed before W09 |
| PDF route too large or too much layout code | W06 | Kill rule: over 1,500 lines or 4 MB keeps PDF a documented recipe, feature off |
| Chaining slower than `qc` | W03 | Kill rule: over 5 percent on chr21.bam keeps `run_isolated` for transform-free chains |
| Genome-scale memory and disk; no 30x BAM fits | W05, W14, W22, W40 | memgate on the CRAM slice and the exome, linear extrapolation, free-space refusal in `plan` |
| VCF header lines Picard adds are unknown for 14 tools | W16 | Haiku inventory of all 14 writers before any VCF writer is built |
| Merge conflicts in the five registration files, 3 branches per wave, 4 other workflows | W02 | Integrator owns the five files; builders ship a snippet; pre-wave check of main and `origin/feat/*` |
| Fingerprint float order and hash-order rows (5 tools) | W26, W27 | log10 summation order ported; `NUM_THREADS=1` exact, else 1e-5 relative; `java.rs` orders |
| 12 tools without real instrument files | W35, W38 | Data spike for Picard 3.5.0 test resources first; synthetic writer by a second agent; tier D label in `tools` |
| Site-set licence | W24 | 1000 Genomes and gnomAD (CC0) only; somalier list not copied |
| Oracle tools not installed; public downloads needed | W06, W40 | Install list approved as one item (Q9); each E-tier acceptance names its tool |
| Agent budget: 572 against a standing cap of 40 total | W02 | Q1; cut order in section 5 |
| Three of the seven in-flight branches are not on origin | W01 | Q4; W01 merges what is pushed and recounts |
| Line estimates 30 percent either way; ValidateSamFile and CollectSamErrorMetrics may split again | W29 to W33 | Three and two workflows reserved; Opus review |
| Disk at the floor from worktree growth | every wave | Worktree removal and `cargo clean` before each wave, hard stop at 25 GB free |
| JSON changes rewrite golden records | W02 | Additive schema; one reviewed `make golden-update` |

## 9. Open questions for the owner

1. Agent cap: is 40 per workflow the reading, so the program may spend 572 agents over 48 workflows? If the budget is smaller, section 5 names what is cut first.
2. Verb name: `chain` for the inline form and `chain --plan FILE` for the file. Add `run` as an alias?
3. Accept krilla (+3.81 MB behind feature `pdf`, 97 crates, rust-version 1.88 to 1.92) with the W06 kill rule?
4. Push `feat/check-duplicate-marking` and `feat/clinical-checks` to origin, and say where `feat/jumping-library-metrics` lives, before W01?
5. Scale gate on a chr20 plus chr21 BAM slice of `NA12878.final.cram` and the exome, with a linear extrapolation to 30x, since a 30x BAM does not fit in 77 GB? Or provide an external disk?
6. FilterSamReads and FilterVcf reject `JAVASCRIPT_FILE`. Confirm.
7. Approve the docs/05 4.9 boundary in section 7 (stateless cohort page open, stored cohorts hosted) before W45?
8. Site sets for identity and contamination from 1000 Genomes and gnomAD, with the somalier list excluded. Confirm.
9. Approve the downloads (Picard 3.5.0 test resources, GIAB HG002 chr20, 1000 Genomes BAMs) and the oracle installs (mosdepth, somalier, VerifyBamID2, FastQC, fastp, rtg vcfeval or hap.py, veraPDF, DuckDB, polars)?
10. Run-folder and array tools ship best effort when no real files exist. Confirm, or supply real files that never enter the repo.
11. PyPI and bioconda releases at phase ends only (six before 1.0.0), or at every minor?
12. Who reviews the clinical thresholds the analyses introduce (contamination 0.02 and 0.05, LOD plus or minus 5, Ti/Tv ranges), and in which wave?

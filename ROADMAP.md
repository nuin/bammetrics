# BAMmetrics platform — execution roadmap (money + reputation)

The committed plan. One open-core platform (BAMmetrics QC engine + vcfclick variant query), SequenceServer-style: a free core and a citable paper build **reputation**; a paid cloud makes **money**. The two goals are sequenced, not traded off — reputation is the distribution engine that sells the paid layer for you, which matters because building and publishing are your strengths and selling is not.

## The core principle

In this model you never have to "sell." You **build** (a free tool) and **publish** (a paper), both of which you're good at. Those create reputation; reputation creates adoption; adoption is the funnel a small paying subset converts through. Order is fixed: reputation precedes revenue. You will not be paid before you are known.

## The one product (no more ideas)

A unified genomics platform with two capabilities that match the workflow — **QC your data, then query it**:
- **BAMmetrics** — Picard-concordant alignment/QC + variant-call benchmarking (Rust + PyO3, on noodles). Already v0.1-verified.
- **vcfclick** — correctness-first cohort variant query with an NL/SQL interface. Already operational.

Library underneath (permissive, adoption), application + cloud on top (AGPL, monetization).

---

## Phase 0b: full Picard parity and more  *(program plan: docs/07)*

Owner direction (2026-11-09): every one of the 121 Picard 3.5.0 tools (23 on main at 0.17.0, 7 in flight, 91 to build), analyses chained in one process without pipes or temp files (`chain`), every output in Picard text and in JSON, CSV, Parquet, XLSX and more, reports in HTML, PDF and DOCX, and 15 analyses beyond Picard (identity, contamination, FASTQ QC, per-exon coverage, doctor, VCF QC, benchmarking, relatedness, cfDNA fragmentomics). `docs/07_full_parity_program.md` holds the 48 waves (0.18.0 to 1.7.0, parity gate at 1.0.0), the machine schedule, the docs/04 and docs/05 changes and the open questions. The 4.9 boundary of docs/05 holds: stored cohorts, history and cross-run comparison stay in the hosted product below.

## Phase 0 — Finish the core  *(BAMmetrics done in 0.2.0; vcfclick open)*

- [x] BAMmetrics: `CollectInsertSizeMetrics`, `CollectAlignmentSummaryMetrics`,
  `CollectQualityYieldMetrics` and `CollectHsMetrics`, concordant with Picard 3.5.0 on 90 cases.
- [x] `noodles-fasta` reference path (`.fai`, one contig at a time, with `Reference::subsequence`).
- [x] Picard-format output files and Picard `KEY=VALUE` command-line syntax, a drop-in for a pipeline.
- [x] One-pass `collect`, threaded BGZF decoding, benchmark: 14.4x on a 4-million-read BAM.
- [x] Python API with row objects, `to_dict()`, type stubs, `py.typed`; wheel from `maturin build --release`.
- [x] CI workflow (fmt, clippy, tests, concordance, wheel), MIT OR Apache-2.0 licence files, CHANGELOG.
- [ ] vcfclick: settle the engine question — **DuckDB** for the open/embeddable core (conda-packageable; chDB is not), chDB reserved for cloud scale.

Shipped since: every accumulation level, `CollectHsMetrics` `PER_TARGET_COVERAGE` and
`THEORETICAL_SENSITIVITY_OUTPUT`, `CollectWgsMetrics` (0.8.0), `CollectVariantCallingMetrics`
(0.11.0), the `qc`, `check` and `report` verbs (0.12.0 to 0.16.0), `MarkDuplicates` (0.18.0),
`CollectSequencingArtifactMetrics`, `CollectOxoGMetrics`, `ConvertSequencingArtifactToOxoG`,
`BamIndexStats`, `CheckDuplicateMarking`, `CollectJumpingLibraryMetrics` and the two clinical
checks `CollectCoverageGaps` and `CollectSexCheckMetrics` (0.19.0), the output formats
(`convert`, `--format` lists, `qc --also`) and the harness for tools that write reads
(0.20.0), the stage layer of in-process chaining (0.21.0), and the `chain` verb with
its plan file, its plan text and the Python `chain` and `step`, with `FifoBuffer` as
registry tool 33 (0.22.0, wave W04 of docs/07; docs/chaining.md, docs/python.md), SAM
text in and out with `ViewSam`, `SortSam` with the sort barrier, the spool and the planner
of a chain, `BuildBamIndex` with Picard's `.bai` behind every index written, and `make
memgate` with the scale gate of docs/07 section 6 (tools 34 to 36, 0.23.0, wave W05),
the PDF report behind the cargo feature `pdf`: `report --pdf`, `qc --report pdf`,
`--pdf-font` and the Python keywords, tagged PDF/UA-1 (0.24.0, wave W06; docs/pdf.md),
and `bammetrics doctor`, `CollectMtCopyNumberMetrics` and `CollectExonCoverageMetrics`
with the GTF, GFF3 and BED12 reader (tools 37 and 38, 0.25.0, wave W07; docs/doctor.md,
docs/06 section 6, docs/exon_coverage.md), and `CollectAdapterContentMetrics`,
`CollectDuplicateByTileMetrics`, `CollectTelomereMetrics` and `CollectIndexDepthMetrics`
(tools 39 to 42, 0.26.0, wave W08; docs/06 sections 7 to 10, docs/telomere.md), and the
record-editing shape with `AddCommentsToBam`, `ReplaceSamHeader`, `AddOATag`, `CleanSam`,
`SetNmMdAndUqTags` and `SetNmAndUqTags`, each a streaming chain transform (tools 43 to 48,
0.27.0, wave W09; docs/chaining.md section 16), and `AddOrReplaceReadGroups` and
`RevertOriginalBaseQualitiesAndAddMateCigar` with the writer's sort in the record-editing
shape, `SplitSamByLibrary` and `SplitSamByNumberOfReads` as fan-out sinks of a chain,
`GatherBamFiles` by BGZF block copy and `CalculateReadGroupChecksum` (tools 49 to 54,
0.28.0, wave W10; docs/chaining.md sections 17 and 18), and `CollectQualityYieldMetricsSNVQ`
and `CollectRrbsMetrics` as chain taps, `CollectUmiPrevalenceMetrics` over the duplicate-set
sort, and `CreateBafRegressMetricsFile`, `CreateVerifyIDIntensityContaminationMetricsFile`
and `CollectArraysVariantCallingMetrics` on text and arrays VCF input (tools 55 to 60,
0.29.0, wave W11), and the interval layer `interval_ops` with `IntervalListTools`,
`IntervalListToBed`, `VcfToIntervalList`, `ScatterIntervalsByNs` and `NonNFastaSize`
(tools 61 to 65, 0.30.0, wave W12), and `MergeSamFiles` with merged sequence
dictionaries, `ReorderSam`, `FixMateInformation`, `FilterSamReads`, `DownsampleSam` with
Picard's seeds and `PositionBasedDownsampleSam`, four of them chain transforms (tools 66
to 71, 0.31.0, wave W13; docs/chaining.md sections 19 to 21), and chain 2 with no new
tool (0.32.0, wave W14; docs/chaining.md sections 22 to 32): steps without `+`
(`bammetrics chain fixed.bam reorder-sam insert-size -o r.bam`), `-o` naming the BAM, the
user config file, the recipes `prep`, `qc-only` and `downsample`, manifest 2 with
`--resume` per step and `--keep`, `qc` as the last step, and every tool that reads a BAM
as a step at any position (53 of the 71) through the materialize fallback, and the CRAM
writer behind every `OUTPUT` and `write` step ending in `.cram`, `SamFormatConverter`,
`CompareSAMs` and `CheckTerminatorBlock` (tools 72 to 74, 0.33.0, wave W15; docs/07 D206
to D235), and the VCF writer with the variant stream of `chain`, the header inventory of
21 Picard VCF tools, `SplitVcfs`, `MakeSitesOnlyVcf` and `MakeVcfSampleNameMap` (tools 75 to
77, 0.34.0, wave W16; docs/chaining.md sections 33 and 34, docs/07 D236 to D256).
CHANGELOG.md has each.

Remaining for BAMmetrics:

- bioconda recipe and a PyPI release (first push of the branch runs the CI workflow).
- The Picard tools still to build, further formats and DOCX reports: the program plan
  `docs/07_full_parity_program.md` (waves W17 onward);
  docs/04 section 4 holds the P9 queue.

Exit criterion: a `pip install`-able library that produces correct, hand-verified numbers.

## Phase 1 — Reputation  *(free; this is the flywheel)*

The entire goal of this phase is to become the known name, at near-zero marketing cost, using only building and publishing.

1. **Ship the library permissively (MIT/Apache) to PyPI + bioconda.** It rides noodles, so it packages cleanly. Permissive license = it becomes a transitive dependency (the cyvcf2 path to real install numbers). bioconda is non-negotiable for reach; the library qualifies because it has no chDB dependency.
2. **Ship the quick-scripting CLI** on the same core — the seqkit-style top-of-funnel utility that gets installed everywhere.
3. **Write the paper — the reputation engine.** JOSS for the library ("Picard-concordant BAM metrics in Rust"); an Application Note later for the platform. SequenceServer's 2019 paper, cited weekly, is what made it ubiquitous. This is the single highest-leverage action and it's pure strength for you.
4. **Visibility through strengths, not selling:** the "Bioinformatics in a Nutshell" Substack, a FOSDEM-bioinformatics or AMP talk, an honest README and concordance report. The concordance report against pinned Picard *is* the credibility artifact.

Exit criterion: installs accruing, first citations, name recognition in the clinical/research-bioinformatics niche.

## Phase 2 — Money  *(paid; on reputation's back)*

1. **Hosted cloud:** QC dashboards per run/pipeline + variant query + collaboration + audit trails + the managed NL/MCP interface + no-server convenience + benchmarking-on-every-pipeline-change (recurring → sticky). AGPL the application core (the monetization lever that pushes commercial users to the cloud/commercial license). DuckDB for the embeddable tier, chDB for scale.
2. **Sell to the paying subset your reputation reached:** biotech, pharma, CROs, sequencing core facilities — research, non-clinical, so no employer COI.
3. **Incorporate** when cloud revenue justifies it (SequenceServer → Pragmatic Genomics Ltd). Sole proprietor is fine to ship the free core; the paid layer benefits from the liability shield and vendor credibility.
4. **Target:** ~$1,500/year × ~67 customers = $100k ARR — reachable because you need few customers and they find you through the phase-1 reputation, not through marketing.

## Phase 3 — High-value expansion  *(optional, later)*

The clinical validation-report tier (the regulated, highest-willingness-to-pay play), once incorporated and AHS/APL COI-cleared. Same engine, premium tier. Don't start here — arrive here.

---

## License logic (it follows the goal)

- **Library** → MIT/Apache: you want it to be everyone's free dependency (adoption = reputation).
- **Application + cloud** → AGPL: you want commercial users to pay (network-copyleft is the lever).

## Honest constraints

- **Timeline is compounding, not fast.** SequenceServer is a ~decade arc (2013 paper → 2021 commercial arm). Reputation matures, then revenue follows. $100k is reachable on this model but it is the slow-flywheel route.
- **The binding risk now is not the idea or the build — it's idea-hopping.** The idea is settled; the build you can do. The failure mode is starting over. Commit to Phase 1.
- **Willingness to pay comes from the commercial subset** (biotech/pharma/CRO/core facilities), never from academics — academics are the reputation, not the revenue.

## The single next action

Finish BAMmetrics to a shippable, MIT-licensed library, get it on PyPI + bioconda, and draft the JOSS paper. The paper is the reputation engine everything else hangs on, and it's the part you're uniquely equipped to execute. Everything in Phase 2 depends on Phase 1 existing first.

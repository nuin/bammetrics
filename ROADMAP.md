# bammetrics platform — execution roadmap (money + reputation)

The committed plan. One open-core platform (bammetrics QC engine + vcfclick variant query), SequenceServer-style: a free core and a citable paper build **reputation**; a paid cloud makes **money**. The two goals are sequenced, not traded off — reputation is the distribution engine that sells the paid layer for you, which matters because building and publishing are your strengths and selling is not.

## The core principle

In this model you never have to "sell." You **build** (a free tool) and **publish** (a paper), both of which you're good at. Those create reputation; reputation creates adoption; adoption is the funnel a small paying subset converts through. Order is fixed: reputation precedes revenue. You will not be paid before you are known.

## The one product (no more ideas)

A unified genomics platform with two capabilities that match the workflow — **QC your data, then query it**:
- **bammetrics** — Picard-concordant alignment/QC + variant-call benchmarking (Rust + PyO3, on noodles). Already v0.1-verified.
- **vcfclick** — correctness-first cohort variant query with an NL/SQL interface. Already operational.

Library underneath (permissive, adoption), application + cloud on top (AGPL, monetization).

---

## Phase 0b: full Picard parity and more  *(program plan: docs/07)*

Owner direction (2026-11-09): every one of the 121 Picard 3.5.0 tools (23 on main at 0.17.0, 7 in flight, 91 to build), analyses chained in one process without pipes or temp files (`chain`), every output in Picard text and in JSON, CSV, Parquet, XLSX and more, reports in HTML, PDF and DOCX, and 15 analyses beyond Picard (identity, contamination, FASTQ QC, per-exon coverage, doctor, VCF QC, benchmarking, relatedness, cfDNA fragmentomics). `docs/07_full_parity_program.md` holds the 48 waves (0.18.0 to 1.7.0, parity gate at 1.0.0), the machine schedule, the docs/04 and docs/05 changes and the open questions. The 4.9 boundary of docs/05 holds: stored cohorts, history and cross-run comparison stay in the hosted product below.

## Phase 0 — Finish the core  *(bammetrics done in 0.2.0; vcfclick open)*

- [x] bammetrics: `CollectInsertSizeMetrics`, `CollectAlignmentSummaryMetrics`,
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
its plan file and plan text (unreleased, wave W04 of docs/07; docs/chaining.md).
CHANGELOG.md has each.

Remaining for bammetrics:

- bioconda recipe and a PyPI release (first push of the branch runs the CI workflow).
- The Picard tools still to build, the Python surface of `chain`, the sort barrier and
  the planner, further formats and PDF reports: the program plan
  `docs/07_full_parity_program.md` (W04 part C, then W05 onward); docs/04 section 4
  holds the P9 queue.

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

Finish bammetrics to a shippable, MIT-licensed library, get it on PyPI + bioconda, and draft the JOSS paper. The paper is the reputation engine everything else hangs on, and it's the part you're uniquely equipped to execute. Everything in Phase 2 depends on Phase 1 existing first.

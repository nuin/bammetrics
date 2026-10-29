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

## Phase 0 — Finish the core  *(bammetrics done in 0.2.0; vcfclick open)*

- [x] bammetrics: `CollectInsertSizeMetrics`, `CollectAlignmentSummaryMetrics`,
  `CollectQualityYieldMetrics` and `CollectHsMetrics`, concordant with Picard 3.5.0 on 90 cases.
- [x] `noodles-fasta` reference path (`.fai`, one contig at a time, with `Reference::subsequence`).
- [x] Picard-format output files and Picard `KEY=VALUE` command-line syntax, a drop-in for a pipeline.
- [x] One-pass `collect`, threaded BGZF decoding, benchmark: 14.4x on a 4-million-read BAM.
- [x] Python API with row objects, `to_dict()`, type stubs, `py.typed`; wheel from `maturin build --release`.
- [x] CI workflow (fmt, clippy, tests, concordance, wheel), MIT OR Apache-2.0 licence files, CHANGELOG.
- [ ] vcfclick: settle the engine question — **DuckDB** for the open/embeddable core (conda-packageable; chDB is not), chDB reserved for cloud scale.

Remaining for bammetrics:

- `CollectVariantCallingMetrics` (the variant-concordance engine).
- `MarkDuplicates` (optical duplicates, duplicate scoring, library complexity; reuses `duplication`).
- `CollectWgsMetrics` (reuses `intervals::OverlapIndex` and the `pass` reader).
- `METRIC_ACCUMULATION_LEVEL` `LIBRARY` and `READ_GROUP`: needs the read-group to library map and
  per-level histogram labels. The corpus BAMs t1 and t2 already carry two libraries.
- `CollectHsMetrics` `PER_TARGET_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT` (`PER_BASE_COVERAGE` is done in 0.3.0).
- bioconda recipe and a PyPI release (first push of the branch runs the CI workflow).

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

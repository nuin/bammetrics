# Picard parity status at 0.34.0

This page is the state of the Picard port at 0.34.0 for a maintainer. It is a snapshot of origin/fix/0.34.0-review. Recount with the commands in section 6 before quoting a number.

## 1. Summary

- Registry rows: 77 (`TOOLS` in `crates/core/src/registry.rs`).
- Picard tools ported: 69 of 121. The other 8 rows are not Picard programs (section 2.2).
- Picard tools remaining: 52. Every one has a wave in docs/07 (W17 to W39).
- Waves left to 1.0.0: 23 (W17 to W20, W22 to W40). W21 shipped early as 0.31.0.
- Waves beyond Picard: 7 (W41 to W47, 1.1.0 to 1.7.0).
- Registry after W40: 129 rows (77 plus 52), 121 of them Picard.
- Release state: 0.34.0 (wave W16, 77 tools, decisions D236 to D256) is unmerged. PR #64 targets main. PR #65 targets release/0.34.0. The last tag is v0.33.0. No PyPI or bioconda release exists. `CHANGELOG.md` has an empty Unreleased heading above 0.34.0.
- Chain coverage: 53 of 77 tools are steps on a BAM stream. 24 are refused there with a reason (`crates/cli/tests/chain_positions.rs`, `assert_eq!((CASES.len(), REFUSED.len()), (53, 24))`).

## 2. Shipped

### 2.1 Registry tools

Chain role is the role on a BAM stream. T is transform, tap is tap, M is materialized (runs on a file of the stream), S is split sink, R is refused. Oracle is whether Picard 3.5.0 produces the expected output.

| Picard name | Verb | Version | Chain role | Picard oracle |
|---|---|---|---|---|
| CollectInsertSizeMetrics | insert-size | 0.1.0 | tap | yes |
| CollectAlignmentSummaryMetrics | alignment-summary | 0.2.0 | tap | yes |
| CollectQualityYieldMetrics | quality-yield | 0.2.0 | tap | yes |
| CollectHsMetrics | hs-metrics | 0.2.0 | tap | yes |
| QualityScoreDistribution | quality-score-distribution | 0.6.0 | tap | yes |
| MeanQualityByCycle | mean-quality-by-cycle | 0.6.0 | tap | yes |
| CollectBaseDistributionByCycle | base-distribution-by-cycle | 0.6.0 | tap | yes |
| CollectDuplicateMetrics | duplicate-metrics | 0.6.0 | tap | yes |
| AccumulateQualityYieldMetrics | accumulate-quality-yield | 0.6.0 | R (no BAM) | yes |
| CollectGcBiasMetrics | gc-bias | 0.6.0 | tap | yes |
| CollectMultipleMetrics | multiple-metrics | 0.6.0 | M | yes |
| EstimateLibraryComplexity | library-complexity | 0.6.0 | M | yes |
| CollectWgsMetrics | wgs | 0.8.0 | M | yes |
| CollectRawWgsMetrics | raw-wgs | 0.8.0 | M | yes |
| CollectWgsMetricsWithNonZeroCoverage | wgs-with-non-zero-coverage | 0.8.0 | M | yes |
| BedToIntervalList | bed-to-interval-list | 0.8.0 | R (no BAM) | yes |
| CollectTargetedPcrMetrics | targeted-pcr | 0.9.0 | M | yes |
| CollectRnaSeqMetrics | rna-seq | 0.10.0 | tap | yes |
| CollectVariantCallingMetrics | variant-calling | 0.11.0 | R on BAM, tap on VCF stream | yes |
| AccumulateVariantCallingMetrics | accumulate-variant-calling | 0.11.0 | R (no BAM) | yes |
| GenotypeConcordance | genotype-concordance | 0.11.0 | R on BAM, tap on VCF stream | yes |
| CompareMetrics | compare-metrics | 0.17.0 | R (no BAM) | yes |
| CreateSequenceDictionary | create-sequence-dictionary | 0.17.0 | R (no BAM) | yes |
| MarkDuplicates | mark-duplicates | 0.18.0 | T (replays barrier) | yes |
| CollectSequencingArtifactMetrics | sequencing-artifact | 0.19.0 | tap | yes |
| ConvertSequencingArtifactToOxoG | convert-sequencing-artifact-to-oxo-g | 0.19.0 | R (no BAM) | yes |
| CollectOxoGMetrics | oxo-g | 0.19.0 | M | yes |
| BamIndexStats | bam-index-stats | 0.19.0 | M (indexed BAM) | yes |
| CheckDuplicateMarking | check-duplicate-marking | 0.19.0 | M | yes |
| CollectJumpingLibraryMetrics | jumping-library | 0.19.0 | M | yes |
| FifoBuffer | fifo-buffer | 0.22.0 | R (a chain has no pipe) | yes |
| ViewSam | view-sam | 0.23.0 | M | yes |
| SortSam | sort-sam | 0.23.0 | T (sort barrier) | yes |
| BuildBamIndex | build-bam-index | 0.23.0 | R (`write --index` writes the .bai) | yes |
| AddCommentsToBam | add-comments-to-bam | 0.27.0 | T (streams) | yes |
| ReplaceSamHeader | replace-sam-header | 0.27.0 | T (streams) | yes |
| AddOATag | add-o-a-tag | 0.27.0 | T (streams) | yes |
| CleanSam | clean-sam | 0.27.0 | T (streams) | yes |
| SetNmMdAndUqTags | set-nm-md-and-uq-tags | 0.27.0 | T (streams) | yes |
| SetNmAndUqTags | set-nm-and-uq-tags | 0.27.0 | T (streams) | yes |
| AddOrReplaceReadGroups | add-or-replace-read-groups | 0.28.0 | T (barrier when it sorts) | yes |
| RevertOriginalBaseQualitiesAndAddMateCigar | revert-original-base-qualities-and-add-mate-cigar | 0.28.0 | T (holds) | yes |
| SplitSamByLibrary | split-sam-by-library | 0.28.0 | S | yes |
| SplitSamByNumberOfReads | split-sam-by-number-of-reads | 0.28.0 | S | yes |
| GatherBamFiles | gather-bam-files | 0.28.0 | M (leaves `gathered`) | yes |
| CalculateReadGroupChecksum | calculate-read-group-checksum | 0.28.0 | M (header only) | yes |
| CollectQualityYieldMetricsSNVQ | quality-yield-s-n-v-q | 0.29.0 | tap | yes |
| CollectUmiPrevalenceMetrics | umi-prevalence | 0.29.0 | M | yes |
| CollectRrbsMetrics | rrbs | 0.29.0 | tap | yes |
| CreateBafRegressMetricsFile | create-baf-regress-file | 0.29.0 | R (no BAM) | yes |
| CreateVerifyIDIntensityContaminationMetricsFile | create-verify-i-d-intensity-contamination-file | 0.29.0 | R (no BAM) | yes |
| CollectArraysVariantCallingMetrics | arrays-variant-calling | 0.29.0 | R (no BAM) | yes |
| IntervalListTools | interval-list-tools | 0.30.0 | R (no BAM) | yes |
| IntervalListToBed | interval-list-to-bed | 0.30.0 | R (no BAM) | yes |
| VcfToIntervalList | vcf-to-interval-list | 0.30.0 | R (no BAM, not a VCF step) | yes |
| ScatterIntervalsByNs | scatter-intervals-by-ns | 0.30.0 | R (no BAM) | yes |
| NonNFastaSize | non-n-fasta-size | 0.30.0 | R (no BAM) | yes |
| MergeSamFiles | merge-sam-files | 0.31.0 | M (leaves `merged`) | yes |
| ReorderSam | reorder-sam | 0.31.0 | T (holds) | yes |
| FixMateInformation | fix-mate-information | 0.31.0 | T (holds) | yes |
| FilterSamReads | filter-sam-reads | 0.31.0 | T (barrier when it sorts) | yes |
| DownsampleSam | downsample-sam | 0.31.0 | T (streams) | yes |
| PositionBasedDownsampleSam | position-based-downsample-sam | 0.31.0 | M (leaves `downsampled`) | yes |
| SamFormatConverter | sam-format-converter | 0.33.0 | R (`write` converts by name) | yes |
| CompareSAMs | compare-s-a-ms | 0.33.0 | R (two files, D233) | yes |
| CheckTerminatorBlock | check-terminator-block | 0.33.0 | R (last bytes of a file, D233) | yes |
| SplitVcfs | split-vcfs | 0.34.0 | R on BAM, split sink on VCF stream | yes |
| MakeSitesOnlyVcf | make-sites-only-vcf | 0.34.0 | R on BAM, transform on VCF stream | yes |
| MakeVcfSampleNameMap | make-vcf-sample-name-map | 0.34.0 | R (no step on any stream) | yes |
| CollectCoverageGaps | coverage-gaps | 0.19.0 | tap | no |
| CollectSexCheckMetrics | sex-check | 0.19.0 | tap | no |
| CollectMtCopyNumberMetrics | mt-copy-number | 0.25.0 | tap | no |
| CollectExonCoverageMetrics | exon-coverage | 0.25.0 | tap | no |
| CollectAdapterContentMetrics | adapter-content | 0.26.0 | tap | no |
| CollectDuplicateByTileMetrics | duplicate-by-tile | 0.26.0 | tap | no |
| CollectTelomereMetrics | telomere | 0.26.0 | tap | no |
| CollectIndexDepthMetrics | index-depth | 0.26.0 | M (indexed BAM) | no |

Role counts on a BAM stream: 14 transforms, 20 taps, 2 splits, 17 materialized, 24 refused. The sum is 77. The last 8 rows are not Picard programs.

### 2.2 Checks with no Picard program

8 registry rows have no Picard counterpart. Their truth comes from planted BAMs and outside tools.

- CollectCoverageGaps (docs/06 section 1): runs of target bases under `MIN_DEPTH`, plus a `GAPS_FILE` BED. `PCT_REPORTABLE` equals Hs `PCT_TARGET_BASES_NX`.
- CollectSexCheckMetrics (docs/06 section 2): reads per base on X, Y and the autosomes against `EXPECTED_SEX`. Thresholds are starting values.
- CollectMtCopyNumberMetrics (docs/06 section 6): 2 x chrM depth over autosomal depth. Truth: planted BAMs and `samtools depth`.
- CollectExonCoverageMetrics (docs/exon_coverage.md): per exon and gene depth from GTF, GFF3, BED12. `MEAN_DEPTH` and `MINIMUM_DEPTH` equal Hs `PER_TARGET_COVERAGE` on 2,674 cells.
- CollectAdapterContentMetrics (docs/06 section 7): FastQC's Adapter Content on a BAM. Oracle: FastQC 0.13.0 on `samtools fastq`.
- CollectDuplicateByTileMetrics (docs/06 section 8): duplicates per lane and tile with `TILE_OUTLIER`. The `ALL` row equals MarkDuplicates and CollectDuplicateMetrics totals.
- CollectTelomereMetrics (docs/telomere.md): TelSeq formula. TelSeq itself did not run (D29).
- CollectIndexDepthMetrics (docs/06 section 10): depth from the `.bai` or `.csi` alone. Oracle: `samtools idxstats`. A CRAM is refused (D23, D32).

`bammetrics doctor` (0.25.0, docs/doctor.md) is a clap verb outside `TOOLS`. It runs 19 pre-flight checks on one BAM, CRAM or SAM.

### 2.3 Platform features

- Picard face (`KEY=VALUE` and `--KEY VALUE`) and one kebab verb per tool, with `--help`, completions and `tools`.
- Python API through PyO3 and maturin: one function per tool, typed stubs, `read_metrics`, `read_table`, `qc`, `check`, `report`, `collect`, `chain`, `step`, `recipe`, `doctor`.
- `qc` and `plan` (0.12.0): the collectors a preset needs on one read, manifest, `--resume`, `--also`, exit codes 10 and 11.
- Verdict presets (0.7.0): generic, clinical-panel, rna-seq, wgs. Each is marked "starting values". `verdict.json` schema is `bammetrics.verdict/1`.
- `check` and `report`: HTML page and run index (axe-core clean, 0.16.0). PDF/UA-1 report behind cargo feature `pdf` (0.24.0, docs/pdf.md).
- Output formats (0.20.0): Picard text, JSON, TSV, CSV, JSON Lines, Markdown, JUnit, MultiQC, `convert`. Parquet, XLSX, DOCX and SQL since 0.31.0 (W21).
- `chain` (0.22.0, W04): `+` grammar, TOML plan file, planner, free-space refusal. Stage layer since 0.21.0 (W03).
- Chain grammar 2 (0.32.0, W14): steps without `+`, `-o` naming the BAM, config file, recipes (`prep`, `qc-only`, `downsample`), manifest 2, `--resume`, `--keep`, `qc` as last step, materialize fallback for 17 tools.
- CRAM reader (0.17.0, CRAM 3.0 and 3.1, needs `R=`). CRAM writer (0.33.0, W15, CRAM 3.0, by `.cram` name).
- SAM text in and out (0.23.0, W05), read as htsjdk reads it.
- VCF writer and VCF stream of a chain (0.34.0, W16, docs/chaining.md sections 33 and 34).
- Gates: golden suite (917 records at 0.33.0, more at 0.34.0), Picard concordance, fuzzer, `make memgate`, `make budgets`, VCF header inventory over 21 Picard VCF tools.

## 3. Gaps in shipped tools

### 3.1 Options not implemented

| Tool | Missing |
|---|---|
| GenotypeConcordance | `OUTPUT_VCF`. Unknown option, exit 1. |
| MarkDuplicates | Flow-based options (accepted at defaults only). `READ_NAME_REGEX` runs on the Rust regex crate, not Java. |
| CollectOxoGMetrics | `CONTEXT_SIZE` above 4 without `CONTEXTS`. Exit 1. |
| CollectSequencingArtifactMetrics | `CONTEXT_SIZE` above 5. |
| FilterSamReads | `FILTER=includeJavascript`, `JAVASCRIPT_FILE` (Q6, D16). FilterVcf will refuse it the same way in W17. |
| CollectWgsMetrics family | `THEORETICAL_SENSITIVITY_OUTPUT`, `ALLELE_FRACTION`. Exit 1 naming the option. |
| CollectMultipleMetrics | `PROGRAM=RnaSeqMetrics` (`NOT_IMPLEMENTED` in `multiple_metrics.rs`). |
| EstimateLibraryComplexity | `ElcHashBasedDuplicatesFinder` for groups of 100 or more. No spill past `MAX_RECORDS_IN_RAM`. |
| CheckDuplicateMarking | `MAX_RECORDS_IN_RAM` and `TMP_DIR` do nothing on an unsorted header. |
| All collectors | `FLOW_MODE=true` exits 1. `CHART` and `HISTOGRAM_FILE` draw no plot. |
| CRAM output | CRAM 3.1, `CREATE_INDEX`, `.crai`. Block copies refuse a `.cram` name. |
| VCF writer, SplitVcfs, MakeSitesOnlyVcf | Tribble `.idx` for a plain `.vcf` (D238), CSI index, BCF, `.md5` from the two tools. |
| Picard face | `--arguments_file` is not read. |

Tools never built: FindMendelianViolations and the 51 others of section 4.

### 3.2 Differences on purpose

`KNOWN_DIFFERENCES.md` is the list. It has 39 `##` sections. Counts by area:

- Output file, Options, Charts, Python, Inputs, Validation, Numbers: 7 sections.
- Picard behaviour reproduced (HashMap order, shared Hs histograms, WGS locus caps): 1 section, about 370 lines.
- Failures and Threads: 2 sections.
- Per tool or tool pair: 20 sections (MarkDuplicates, SortSam, SamFormatConverter, the W09 to W13 tools, CompareSAMs and CheckTerminatorBlock, ReorderSam).
- Formats and streams (Chains, SAM text, ViewSam, CRAM output, Interval lists, VCF writer and the variant stream): 6 sections.
- Not in Picard, Not implemented, Fuzzing: 3 sections. The fuzz run of seeds 1 to 3 gave 0 failures in 30,000 commands.

## 4. Remaining

Each wave runs under the base gate: fmt, clippy, tests, golden, concordance, round trip, chained equals staged, python-test, wheel-check, budgets, memgate, version bump. The exit gate column is the extra gate. Truth tiers: A is Picard text equal, B is Picard BAM or VCF decoded equal, C is Picard model output equal on synthetic input, D is best effort, E is an outside tool.

### 4.1 Waves W17 to W40 (52 Picard tools)

| Wave | Version | Tools | Depends on | Truth tier | Exit gate | Data needed |
|---|---|---|---|---|---|---|
| W17 | 0.35.0 | RenameSampleInVcf, UpdateVcfSequenceDictionary, FixVcfHeader, FilterVcf, GatherVcfs (5) | W16 | B | 5 tools equal apart from inventoried header lines | Exists: `calls3.vcf`, `variant_vcfs.py`, `vcf_headers.json`. Missing: multi-shard VCF set (split by script). Port the tribble `.idx` writer (D238). |
| W18 | 0.36.0 | SortVcf, MergeVcfs, FindMendelianViolations, CombineGenotypingArrayVcfs, MergePedIntoVcf (5) | W17 | B | SortVcf on 10 M records under 1 GiB RSS; the 3.5.0 `VCF_DIR` bug reproduced or rejected | Missing: PED, MAP, zCall files, 10 M record VCF, trio VCF with GT, GQ, DP, AD, PL. All synthetic. |
| W19 | 0.37.0 | LiftOverIntervalList, LiftoverVcf, ExtractSequences, NormalizeFasta, SortGff (5) | W18, W12 | A for text, B for VCF | Liftover equal, reject file included | Missing: chain files, second reference, GFF3. All synthetic. |
| W20 | 0.38.0 | VcfFormatConverter (1) | W17 | B after BGZF inflate | BCF read by every VCF tool, both directions equal | Missing: any BCF file (made by Picard and bcftools). |
| W22 | 0.40.0 | FastqToSam, SamToFastq, SamToFastqWithTags, RevertSam, BamToBfq (5) | W13 | A for FASTQ and bfq, B for gz and BAM | 5 tools equal; FASTQ reader feeds the 4 read-level collectors; scale run | Missing: FASTQ files (cut from the corpus). |
| W23 | 0.41.0 | MarkIlluminaAdapters, MergeBamAlignment (2) | W22 | B | Every option equal on 3 aligned inputs; XT tags equal | Missing: uBAM, bwa-aligned BAMs (made in the wave). |
| W24 | 0.42.0 | MarkDuplicatesWithMateCigar, SimpleMarkDuplicatesWithMateCigar (2) | W13 | B for BAM, A for metrics | Equal on the markdup corpus and `make fuzz SEED=1 N=200` | Exists: markdup corpus. Q8 licence spike for 1000 Genomes and gnomAD site sets. |
| W25 | 0.43.0 | UmiAwareMarkDuplicatesWithMateCigar, CollectIndependentReplicateMetrics (2) | W24 | B and A, D for the replicate tool | UMI inference equal; replicate tool labelled if two Picard runs differ | Missing: UMI families at edit distance 1, replicate sets of 2 and 3 (hand built). |
| W26 | 0.44.0 | CheckFingerprint, ConvertHaplotypeDatabaseToVcf, LiftOverHaplotypeMap (3) | W19, W20 | C | 100 seeded cases, LOD equal at 6 digits, log10 summation order ported | Missing: haplotype map and people genotypes. All synthetic. |
| W27 | 0.45.0 | CrosscheckFingerprints, CrosscheckReadGroupFingerprints, ClusterCrosscheckMetrics (3) | W26 | C | Row order equal under HashMap order for 3 library counts | Missing: multi-library fingerprint inputs (synthetic). |
| W28 | 0.46.0 | ExtractFingerprint, IdentifyContaminant, CalculateFingerprintMetrics (3) | W26 | C | Seeded sampling equal (static `Random(42)`); VCF equal apart from header lines | Missing: mixed samples at known contamination (planted). |
| W29 | 0.47.0 | ValidateSamFile, part 1 | W15 | A | Record, flag, CIGAR and mate rule messages equal | Missing: one broken record per rule (about 70 files over W29 to W31). |
| W30 | 0.48.0 | ValidateSamFile, part 2 | W29 | A | Each header, RG, PG, sort and tag message equal | Missing: broken header, RG, PG, sort files. |
| W31 | 0.49.0 | ValidateSamFile, part 3 | W30 | A | About 70 rules covered; rules 3.5.0 cannot trigger listed | Missing: broken CRAM and index files. |
| W32 | 0.50.0 | CollectSamErrorMetrics, part 1 | W17 | A | 5 stratifier files equal | Exists: corpus, `variant_vcfs.py`. Needs an indexed VCF. |
| W33 | 0.51.0 | CollectSamErrorMetrics, part 2 | W32 | A | Every stratifier equal | As W32. docs/04 lists about 30 stratifiers. |
| W34 | 0.52.0 | BaitDesigner, CollectQualityYieldMetricsFlow (2) | W12 | A for BaitDesigner, D for Flow | BaitDesigner equal; Flow labelled best effort | Missing: any public Ultima BAM; soft-masked reference (made). |
| W35 | 0.53.0 | CheckIlluminaDirectory, CollectIlluminaLaneMetrics, CollectHiSeqXPfFailMetrics (3) | W11 | A on fetched fixtures, D otherwise | Tier A where fixtures exist; CBCL without a fixture best effort | Missing: any run folder, BCL, CBCL, filter, locs, RunInfo.xml, InterOp file. Spike on Picard 3.5.0 testdata and nf-core bcl data. |
| W36 | 0.54.0 | CollectIlluminaBasecallingMetrics, ExtractIlluminaBarcodes (2) | W35 | as W35 | Equal on every available folder | As W35, with barcodes. |
| W37 | 0.55.0 | IlluminaBasecallsToFastq, IlluminaBasecallsToSam (2) | W36, W22 | as W35 | 96-barcode run under 2 GiB RSS | As W35, plus a synthetic 96-barcode run. |
| W38 | 0.56.0 | BpmToNormalizationManifestCsv, CompareGtcFiles, CreateExtendedIlluminaManifest (3) | W19 | A on fixtures, D otherwise | Round trip passes; `tools` shows the best-effort label | Missing: any BPM, EGT, GTC file. `tests/data/arrays` holds text and VCF only. |
| W39 | 0.57.0 | GtcToVcf, VcfToAdpc (2) | W38 | A for adpc, B for VCF | adpc `cmp` equal | As W38. |
| W40 | 1.0.0 | none (gate) | W39 | all | 121 Picard tools registered; `fuzz N=500` on BAM writers; docs/04 and docs/05 edits (D1, D3); README table with wave, tier and covering tool; PyPI release | None new. Scale run. |

Tools by wave sum to 52: 5+5+5+1+5+2+2+2+3+3+3+1+1+2+3+2+2+3+2.

Agent budget in docs/07: W17 to W40 without W21 sum to 266, W41 to W47 sum to 78. Total 344.

### 4.2 Waves W41 to W47 (beyond Picard)

| Wave | Version | Content | Depends on | Truth tier | Exit gate |
|---|---|---|---|---|---|
| W41 | 1.1.0 | FASTQ QC: FastqMetrics, overrepresented sequences (cap 100,000 prefixes), adapter content by cycle, per-tile quality, preset fastq | W22 | A through FastqToSam, E against FastQC and fastp (not installed) | 4 collectors byte equal; throughput published first |
| W42 | 1.2.0 | VCF QC: VariantQcMetrics, StratifiedConcordanceMetrics (size bins, left-align) | W20, Q9 | E: `bcftools stats` exact, hap.py or rtg vcfeval on HG002 chr20 (not installed) | Counts exact, F1 within 0.1 percent |
| W43 | 1.3.0 | Identity and contamination: site sets (identity 70 KB, contamination 400 KB, D10), IdentityMetrics, ContaminationMetrics | W28, Q8, item 12 | E: `samtools view -s` mixes at 0, 1, 2, 5, 10 percent; VerifyBamID2 (not installed) | Within 0.005 of planted at 30x; INCOMPLETE under 5,000 sites |
| W44 | 1.4.0 | Relatedness and cfDNA fragmentomics: RelatednessMetrics, FragmentMetrics, end-motif file | W43 | E: somalier within 0.02 (not installed) | Explicit paths only, nothing stored |
| W45 | 1.5.0 | Cohort report, HTML and PDF, from directories on one command line | W21 (shipped), Q7 | self checked, 100 synthetic samples | HTML under 1.5 MB; PDF bookmark per sample |
| W46 | 1.6.0 | Streaming sinks, Python frames, SQLite. Shipped early: VCF stream (W16), Parquet read (D152), `to_pandas`, `to_polars`. Left: `locus.rs` stream sink, Arrow IPC (+0.46 MB), native `.sqlite` (+1.55 MB), variant barrier from SortVcf | W21 (shipped), W18 | chained equals staged | wgs and hs mid-chain with no temp BAM; `to_polars` dtypes equal glossary |
| W47 | 1.7.0 | Haplotype-aware benchmarking: variant-path matching | W42 | E: vcfeval on HG002 chr20 | Within 0.1 percent of vcfeval |

Out of the open crates by the docs/07 section 7 boundary: stored cohorts, runs over time, baselines, audit trail, users, uploads, a cohort-trained cfDNA score. The somalier site list is not copied (Q8).

## 5. Blockers

### 5.1 Data that does not exist

The tree has no BCF, GTC, BPM, EGT, BCL, CBCL, chain, GFF, PED or FASTQ file under `tests/data`.

- W17: multi-shard VCF set.
- W18: PED, MAP, zCall files, 10 M record VCF, trio VCF.
- W19: chain files, second reference, GFF3.
- W20: BCF files.
- W22: FASTQ files.
- W23: uBAM and bwa-aligned BAMs.
- W25 to W28: UMI families, replicate sets, haplotype maps, mixed samples.
- W29 to W31: about 70 one-record broken files.
- W34: an Ultima BAM.
- W35 to W39: run folders, GTC, BPM, EGT files.
- W41 to W47: no oracle installed for fastp, mosdepth, somalier, VerifyBamID2, hap.py or vcfeval, veraPDF. FastQC 0.13.0 was used for W08.

Real-data check gaps: 6 scale memgate rows are skipped (`TIER=scale make memgate` passes 41, skips 6, inputs under `~/data/bammetrics-scale` were removed). No real VCF exists under `~/data`. The 2 W16 scale rows read synthetic data from `benches/vcf_scale.py` (1,000,000 records, 57.1 and 15.3 MiB peak RSS).

### 5.2 Open questions

| Id | Question | State |
|---|---|---|
| Q1 | Agent cap per workflow | Applied as 14 per workflow, 8 at once. The global guide says 40 total. Unconfirmed. |
| Q2 | Verb name `chain` | Answered: `chain` only. |
| Q3 | krilla for PDF | Answered: accepted behind feature `pdf`. Shipped 0.24.0. |
| Q4 | Push feat branches | Resolved in W01 (0.19.0). No written answer. |
| Q5 | Scale gate on a CRAM slice | Answered and done (W05). |
| Q6 | `JAVASCRIPT_FILE` refused | Answered: refused (D16). FilterSamReads done, FilterVcf in W17. |
| Q7 | docs/05 4.9 stateless cohort page | Open. Blocks W45. |
| Q8 | Site sets from 1000 Genomes and gnomAD | Open. Blocks the W24 licence spike and W43. |
| Q9 | Public downloads and oracle installs | Answered as one item. Installs partial: mosdepth absent, TelSeq did not run. W35 and W42 downloads not started. |
| Q10 | Run-folder and array tools ship best effort | Open for W35. |
| Item 11 | PyPI and bioconda at phase ends only | Open. No release made. Phase ends: 0.46.0, 0.57.0, 1.0.0 are named in docs/07. |
| Item 12 | Reviewer of clinical thresholds (0.02 and 0.05 contamination, LOD plus or minus 5) | Open. Blocks W26, W28, W43. |

### 5.3 Open problems

- 0.34.0 is unmerged: PR #65 into release/0.34.0, then PR #64 into main. A fix PR comes before each release merge. No CI minutes exist, so local checks are the gate.
- Clinical-panel preset review is open. Gap and sex rules and all thresholds of generic, wgs and rna-seq are starting values.
- D255 (no sink finishes until every sink took every batch) does not cover a sink that fails while closing its files.
- CRAM writer departures D206 to D223 stay. One input was written at scale (HG002 exome). Output is 16.6 percent larger and 3.7 times slower than samtools' CRAM. CRAM bytes are never compared with Picard's.
- BAM bytes are never compared: zlib-rs and miniz_oxide reproduce 0 of 585 Picard blocks. `.tbi` bytes differ from Picard's and are checked by `tabix` queries.
- The VCF header inventory skips GtcToVcf (no GTC, EGT or BPM file) and BCF outputs (no codec until W20).
- VCF reader limits: header over 64 MiB and line over 16 MiB fail (D245). The noodles reader differs from htsjdk on some files (D92). BCF is refused everywhere. Standard input is never a VCF source (D241).
- Java HashMap tree bins are not emulated in `java.rs` for MakeVcfSampleNameMap, CollectOxoGMetrics and EstimateLibraryComplexity. Inputs that form one exit 1.
- Stale lines in `KNOWN_DIFFERENCES.md`: Not implemented still lists CRAM output (shipped 0.33.0), and line 539 of Inputs still says SAM text input exits 1 in every tool (shipped 0.23.0).
- docs/07 drift: section 0 still places VCF streams at 1.6.0 and formats at 0.39.0 (both shipped). Slot 0.39.0 has no wave. docs/07 says 12 best-effort tools, but W34 to W39 hold 13 tools.

## 6. Sources and recounts

Drawn from:

- `docs/04_picard_parity_plan.md`, `docs/07_full_parity_program.md` (sections 0, 4 to 9, D1 to D256).
- `ROADMAP.md`, `README.md`, `CHANGELOG.md`, `KNOWN_DIFFERENCES.md`, project guide `CLAUDE.md`.
- `crates/core/src/registry.rs`, `crates/cli/src/verbs.rs`, `crates/cli/tests/chain_positions.rs`, `crates/core/src/stage/materialize.rs`.
- `docs/chaining.md`, `docs/concordance.md`, `docs/benchmarks.md`.

Recount commands:

```bash
git fetch -q origin
git grep -c '&crate::.*::TOOL,' origin/fix/0.34.0-review -- crates/core/src/registry.rs   # 77
cargo test -p bammetrics-cli --test registration                      # every tool has a verb, glossary row, concordance module, doc line
cargo test -p bammetrics-cli --test chain_positions                   # asserts (53, 24)
git grep -n 'REFUSED' crates/cli/tests/chain_positions.rs
git grep -nE '^\| [A-Z][A-Za-z]+ \|' origin/fix/0.34.0-review -- docs/04_picard_parity_plan.md   # tool rows of docs/04
```

To recount the remaining 52, list the docs/04 tool names, remove every `Tool.name` in `TOOLS` and the 8 non-Picard rows, and match the rest against the Items cells of docs/07 section 5.

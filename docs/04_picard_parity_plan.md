# Picard parity plan

Status: draft for the maintainer and the agents who implement each tool. Target: field-for-field output identical to Picard 3.5.0. Date: 2026-10-29.

Sources: the per-tool inventory (49 tools, each run on Picard 3.5.0), a triage of the other 74 tools, two validation designs, `CLAUDE.md`, `ROADMAP.md` and the v0.1 implementation guide. The two tool lists overlapped by 2 tools, so the union is 121. "Verified" means the column list or behaviour was seen in a Picard 3.5.0 run. Two tools are not verified by a run and are marked. Nothing here was changed in the repo.

## 0. How to read this

- A "decision" is one of core (build for 1.0), later (build after 1.0 or as a cheap by-product), out (not built, with a reason).
- Counts, from the 121 programs in Picard 3.5.0: **19 core, 37 later, 65 out**.
- The inventory covered 49 tools: 17 core, 25 later, 7 out. The triage covered 74: 2 core, 12 later, 60 out. Two out tools were counted in both inputs, so 49 + 74 = 123 and the union is 121.
- The tool list was diffed against `picard` with no arguments (ANSI stripped): 121 names, and the table in section 1.2 has 121 rows. No name is missing, extra or duplicated. The diff is done.
- Column lists are verified by Picard 3.5.0 runs (spot-checked for InsertSize, AlignmentSummary, QualityYield, BaseDistributionByCycle, CollectDuplicateMetrics, CollectWgsMetrics and GcBias). Accumulation-level behaviour (READ_GROUP value, `unknown` rows, extra all-zero rows) is not verified: every run so far used ALL_READS only. See Q-001 in section 5.7.
- Repo state: the plan is rebased on `origin/main` (version 0.3.0). The local branch `feat/picard-metrics` has 1 commit not on main (merge base `af8b521`). `origin/main` has 23 commits not on the branch. Main already holds the Cargo workspace and clap CLI (M1), `metrics_file.rs` (M2), `intervals.rs` (M15), `sensitivity.rs` (M16), `collect.rs` (a one-pass multi-collector, M9 in part), `alignment_summary.rs`, `quality_yield.rs`, `hs_metrics.rs`, `duplication.rs`, `tests/concordance/` (90 cases, Picard-syntax CLI), `KNOWN_DIFFERENCES.md` and a CHANGELOG. Hs, including HS_LIBRARY_SIZE, HET_SNP_SENSITIVITY and PER_BASE_COVERAGE, shipped without the Wgs pileup (M14). ROADMAP on main says "bammetrics done in 0.2.0". Work that main already holds is marked "verify against section 5" in section 4, not "build". Open work on main: accumulation levels (main accepts ALL_READS only), PER_TARGET_COVERAGE, THEORETICAL_SENSITIVITY_OUTPUT, STOP_AFTER, and every tool not on main.

## 1. Scope

### 1.1 Decision counts

| Decision | Inventory tools | Triage tools | Overlap counted twice | Total |
|---|---:|---:|---:|---:|
| core | 17 | 2 | 0 | 19 |
| later | 25 | 12 | 0 | 37 |
| out | 7 | 60 | 2 | 65 |
| all | 49 | 74 | 2 | 121 |

### 1.2 Table of all tools

The "Topic (ours)" column is our own grouping, not Picard's. Picard 3.5.0 has 10 groups: Diagnostics and Quality Control 47, Read Data Manipulation 31, Variant Manipulation 12, Genotyping Arrays Manipulation 9, Base Calling 7, Reference 6, Intervals Manipulation 4, Other 2, Variant Evaluation and Refinement 2, Variant Filtering 1. For example CollectGcBiasMetrics, MarkDuplicates and the fingerprint tools all sit in Diagnostics and Quality Control. "Covered elsewhere by" names a tool that gives the same or a similar answer. "Partial" means the numbers differ in definition, so it is not a substitute for concordance.

| Tool | Topic (ours) | Decision | Reason | Covered elsewhere by |
|---|---|---|---|---|
| CollectAlignmentSummaryMetrics | Metrics | core | Most-run alignment QC, next milestone in ROADMAP. L | samtools stats, flagstat (partial) |
| CollectInsertSizeMetrics | Metrics | core | In v0.1 for ALL_READS. Needs levels, widths, options. M | samtools stats IS (partial) |
| CollectQualityYieldMetrics | Metrics | core | Single-pass counters, no reference. S | fastp, seqkit stats (partial) |
| CollectQualityYieldMetricsSNVQ | Metrics | later | Needs Ultima qa/qc/qg/qt tags. About 150 lines on top of QualityYield. S | none |
| CollectQualityYieldMetricsFlow | Metrics | out | Experimental, Ultima only, needs FlowBasedRead port and flow fixtures. L | none |
| QualityScoreDistribution | Metrics | core | About 100 lines of integer counts. S | none |
| MeanQualityByCycle | Metrics | core | Under 200 lines, one cycle-offset trap. S | FastQC (partial) |
| CollectBaseDistributionByCycle | Metrics | core | About 150 lines of integer counts. S | FastQC (partial) |
| CollectWgsMetrics | Metrics | core | Most-run coverage QC, owns the pileup engine. L | mosdepth, samtools coverage (partial) |
| CollectRawWgsMetrics | Metrics | core | Preset of the WGS engine plus a class name. S | none |
| CollectWgsMetricsWithNonZeroCoverage | Metrics | later | About 60 lines after WGS. Chart needs R, skipped. S | none |
| CollectHsMetrics | Metrics | core | Main clinical capture metric, 70 columns. XL | mosdepth plus scripts (partial) |
| CollectTargetedPcrMetrics | Metrics | later | Same collector as Hs with other defaults. M | none |
| MarkDuplicates | Duplicates | later | Most-used tool, writes BAMs, needs ReadEnds. XL | samtools markdup (different algorithm) |
| MarkDuplicatesWithMateCigar | Duplicates | out | Gives different flags than MarkDuplicates on the same input. L | none |
| SimpleMarkDuplicatesWithMateCigar | Duplicates | out | Test aid in the Picard source. M | none |
| UmiAwareMarkDuplicatesWithMateCigar | Duplicates | later | UMI graph, fails on multi-library input in Picard. L | fgbio, UMI-tools (partial) |
| CollectDuplicateMetrics | Duplicates | core | Pure flag counting, builds the DuplicationMetrics block. S | samtools stats (partial) |
| EstimateLibraryComplexity | Duplicates | core | Works on unaligned input, no BAM output. M | none |
| CheckDuplicateMarking | Duplicates | later | Validator, exit code is the output. S | none |
| CollectUmiPrevalenceMetrics | Duplicates | later | Small output, but inverted barcode filter to copy. M | none |
| CollectIndependentReplicateMetrics | Duplicates | out | Alpha research tool, niche protocol. L | none |
| CollectGcBiasMetrics | Metrics | core | Standard WGS bias QC, reference windows plus read starts. M | Picard only (deepTools computeGCBias, partial) |
| CollectSequencingArtifactMetrics | Metrics | later | 5 output files, pro/con rules. FFPE and OxoG use. L | none |
| CollectOxoGMetrics | Metrics | later | Needs pileup, Java hash row order. M | none |
| ConvertSequencingArtifactToOxoG | Metrics | later | Text conversion, needs the artifact tool. S | none |
| CollectSamErrorMetrics | Metrics | later | About 30 stratifiers, research tool. XL | none |
| CollectRrbsMetrics | Metrics | later | RRBS niche, level quirks. M | Bismark reports (partial) |
| CollectRnaSeqMetrics | Metrics | core | Mainstream RNA-seq QC (nf-core/rnaseq). L | RSeQC, Qualimap (partial) |
| CollectJumpingLibraryMetrics | Metrics | later | Mate-pair libraries are rare. S | none |
| CollectMultipleMetrics | Metrics | core | Most users call this wrapper. Single pass is a speed argument. M | none |
| BamIndexStats | Metrics | later | Same numbers as samtools idxstats. S | samtools idxstats |
| ValidateSamFile | Diagnostics | out | Defined by htsjdk parser behaviour over about 70 rules. XL | samtools quickcheck (partial) |
| CompareMetrics | Diagnostics | core | The oracle comparator. Its semantics go into the harness. S | none |
| CompareSAMs | Diagnostics | later | Pipeline regression check. Quirks to copy. M | none |
| CalculateReadGroupChecksum | Diagnostics | later | 40 lines, verified algorithm. S | none |
| CheckTerminatorBlock | Diagnostics | later | 50 lines, pre-flight check. S | noodles-bgzf, samtools quickcheck |
| CollectVariantCallingMetrics | Variant | core | Variant-benchmarking pillar of the roadmap. M | bcftools stats (partial) |
| AccumulateVariantCallingMetrics | Variant | later | Lossy shard merge, copy the inversion. S | none |
| AccumulateQualityYieldMetrics | Metrics | later | Needs QualityYield first. S | none |
| GenotypeConcordance | Variant | core | The benchmarking tool people cite. L | hap.py, vcfeval (different method) |
| FindMendelianViolations | Variant | later | Deterministic, needs VCF layer. M | bcftools +mendelian2 (partial) |
| CollectArraysVariantCallingMetrics | Arrays | out | Broad GtcToVcf output only. M | none |
| CrosscheckFingerprints | Fingerprint | later | Identity check, needs pileup and likelihood model. XL | NGSCheckMate, somalier (different method) |
| CheckFingerprint | Fingerprint | later | Reuses Crosscheck core. L | none |
| ExtractFingerprint | Fingerprint | later | VCF output only. M | none |
| CalculateFingerprintMetrics | Fingerprint | later | Needs Mersenne Twister and chi-squared. L | none |
| ClusterCrosscheckMetrics | Fingerprint | later | Java hash order. S | none |
| IdentifyContaminant | Fingerprint | out | Niche, VCF output only. M | VerifyBamID |
| BedToIntervalList | Intervals | core | Hs, TargetedPcr and Wgs take interval_list. Harness needs it. | none |
| IntervalListToBed | Intervals | later | Reverse of the parser. | none |
| IntervalListTools | Intervals | later | Implement only the ops metrics use (merge, pad, uniquify). | bedtools |
| LiftOverIntervalList | Intervals | out | Liftover, no metrics. | UCSC liftOver, CrossMap |
| CreateSequenceDictionary | Reference | core | Metrics need contig lengths. Harness prep. | samtools dict |
| ExtractSequences | Reference | out | Region FASTA subsetting. | samtools faidx |
| NonNFastaSize | Reference | later | One number, WGS computes it internally. | seqkit stats |
| NormalizeFasta | Reference | out | FASTA rewrap. | seqkit seq -w |
| ScatterIntervalsByNs | Reference | out | Pipeline scatter helper. | bedtools |
| BaitDesigner | Reference | out | Legacy bait design. | none |
| FixMateInformation | Read data | later | Harness input prep (MC tag). | samtools fixmate -m |
| SetNmMdAndUqTags | Read data | later | Harness prep. Metrics use the reference, not tags. | samtools calmd |
| SetNmAndUqTags | Read data | out | Deprecated. | samtools calmd |
| BuildBamIndex | Read data | later | Harness prep. | samtools index |
| SortSam | Read data | later | Harness prep. | samtools sort |
| SamFormatConverter | Read data | out | Format conversion. | samtools view |
| AddCommentsToBam | Read data | out | Header edit. | samtools reheader |
| AddOATag | Read data | out | No metric reads OA. | none |
| AddOrReplaceReadGroups | Read data | later | Harness fixture builder for per-RG rows. | samtools addreplacerg |
| BamToBfq | Read data | out | Obsolete aligner format. | none |
| CleanSam | Read data | later | Harness may need it to make fixtures valid. | none |
| DownsampleSam | Read data | later | Harness shrinks real BAMs once. | samtools view -s |
| PositionBasedDownsampleSam | Read data | out | Tile-based downsampling for duplicate studies. | none |
| FastqToSam | Read data | later | Harness builds unaligned fixtures. | samtools import, fgbio |
| FilterSamReads | Read data | out | Subsetting. | samtools view |
| GatherBamFiles | Read data | out | Concatenation. | samtools cat |
| MergeSamFiles | Read data | out | Merge. | samtools merge |
| MergeBamAlignment | Read data | out | uBAM merge step. | fgbio ZipperBams |
| ReorderSam | Read data | out | Contig reorder. | none |
| ReplaceSamHeader | Read data | out | Header swap. | samtools reheader |
| RevertOriginalBaseQualitiesAndAddMateCigar | Read data | out | Legacy pipeline step. | samtools fixmate (MC only) |
| RevertSam | Read data | out | BAM to uBAM. | none |
| SamToFastq | Read data | out | Export. | samtools fastq |
| SamToFastqWithTags | Read data | out | Export with tags. | samtools fastq -T |
| SplitSamByLibrary | Read data | out | File split. | samtools split |
| SplitSamByNumberOfReads | Read data | out | File split. | samtools split |
| ViewSam | Diagnostics | out | Prints records. | samtools view |
| CollectHiSeqXPfFailMetrics | Diagnostics | out | Needs basecall directory. | none |
| CrosscheckReadGroupFingerprints | Diagnostics | out | Deprecated, CrosscheckFingerprints covers it. | CrosscheckFingerprints |
| ConvertHaplotypeDatabaseToVcf | Diagnostics | out | Haplotype map conversion. | none |
| LiftOverHaplotypeMap | Diagnostics | out | Haplotype map liftover. | none |
| CheckIlluminaDirectory | Base calling | out | Reads BCL folders. | none |
| CollectIlluminaBasecallingMetrics | Base calling | out | Reads BCL. | bcl-convert reports |
| CollectIlluminaLaneMetrics | Base calling | out | Reads InterOp. | Illumina InterOp |
| ExtractIlluminaBarcodes | Base calling | out | Demultiplexing. | bcl-convert |
| IlluminaBasecallsToFastq | Base calling | out | BCL conversion. | bcl-convert |
| IlluminaBasecallsToSam | Base calling | out | BCL conversion. | bcl-convert |
| MarkIlluminaAdapters | Base calling | out | Adapter tagging in uBAM pipelines. | fastp, cutadapt |
| BpmToNormalizationManifestCsv | Arrays | out | Array manifest. | none |
| CombineGenotypingArrayVcfs | Arrays | out | Array VCF merge. | bcftools merge |
| CompareGtcFiles | Arrays | out | GTC diff. | none |
| CreateBafRegressMetricsFile | Arrays | out | Wraps bafRegress. | none |
| CreateExtendedIlluminaManifest | Arrays | out | Manifest builder. | none |
| CreateVerifyIDIntensityContaminationMetricsFile | Arrays | out | Wraps VerifyIDIntensity. | none |
| GtcToVcf | Arrays | out | GTC to VCF. | bcftools +gtc2vcf |
| MergePedIntoVcf | Arrays | out | zCall merge. | bcftools |
| VcfToAdpc | Arrays | out | Array export. | none |
| FifoBuffer | Other | out | Stream buffer. | mbuffer, pv |
| SortGff | Other | out | GFF sort. | gff3sort |
| FilterVcf | Variant | out | Hard filtering. | bcftools filter |
| FixVcfHeader | Variant | out | Header edit. | bcftools reheader |
| GatherVcfs | Variant | out | Scatter-gather. | bcftools concat |
| LiftoverVcf | Variant | out | Liftover. | CrossMap |
| MakeSitesOnlyVcf | Variant | out | Drops genotypes. | bcftools view -G |
| MakeVcfSampleNameMap | Variant | out | TSV helper. | bcftools query -l |
| MergeVcfs | Variant | out | VCF merge. | bcftools merge |
| RenameSampleInVcf | Variant | out | Rename. | bcftools reheader -s |
| SortVcf | Variant | out | Sort. | bcftools sort |
| SplitVcfs | Variant | out | SNP/indel split. | bcftools view -v |
| UpdateVcfSequenceDictionary | Variant | out | Header contig update. | bcftools reheader |
| VcfFormatConverter | Variant | out | VCF/BCF conversion. | bcftools view |
| VcfToIntervalList | Variant | later | Source of Hs-style targets from VCF sites. | bcftools query plus BedToIntervalList |

Rules that follow from the table:
- "later" harness-prep tools (SortSam, BuildBamIndex, FixMateInformation, SetNmMdAndUqTags, AddOrReplaceReadGroups, CleanSam, DownsampleSam, FastqToSam) are never built as Rust ports. The harness calls samtools or Picard for them. Their "later" status means "prepare fixtures with the covering tool".
- Out-of-scope tools get one line in the README so users know the covering tool.

## 2. Feature inventory

Format per tool: inputs, options that change numbers, metrics class and columns, histograms, traps, complexity. All column lists were seen in a Picard 3.5.0 run unless marked UNVERIFIED. Level behaviour (accumulation levels beyond ALL_READS) is not verified (see Q-001). Numbers print with 6 decimals, trailing zeros dropped, integer-valued doubles without a point, NaN as `?` (see section 3, writer).

Common to BAM tools: `VALIDATION_STRINGENCY` has no effect on counts when SILENT. Chart options (`H=`, `CHART_OUTPUT`) need R and are not implemented, the harness uses a no-op `Rscript` shim so Picard finishes. Option presence is not uniform: `ASSUME_SORTED` and `STOP_AFTER` exist on QualityYield and RnaSeq, but CollectHsMetrics has neither and CollectWgsMetrics has `STOP_AFTER` only. Option presence is therefore taken only from the generated `tools/<Tool>.yaml`, never from prose. The accumulation-level label printed in the level cell differs per tool: most tools emit `ALL_READS`, but CollectGcBiasMetrics prints `All Reads` (with a space) and puts `ALL` in READS_USED. The comparator normalises level names per tool in the row key. Levels SAMPLE, LIBRARY, READ_GROUP rows follow ALL_READS when requested. The claim that the READ_GROUP value is the RG `PU` field, `unknown` when PU or RG is missing, is unconfirmed (Q-001).

### 2.1 Core tools

#### CollectInsertSizeMetrics (M, in v0.1 partly)
- Inputs: BAM/SAM/CRAM. No reference.
- Number-changing options: `DEVIATIONS=10.0`, `HISTOGRAM_WIDTH`, `MIN_HISTOGRAM_WIDTH`, `MINIMUM_PCT=0.05` (float, valid 0 to 0.5), `INCLUDE_DUPLICATES=false`, `METRIC_ACCUMULATION_LEVEL=[ALL_READS]`.
- Class `picard.analysis.InsertSizeMetrics` (23 columns: 19 value columns plus PAIR_ORIENTATION, SAMPLE, LIBRARY, READ_GROUP as key columns): MEDIAN_INSERT_SIZE, MODE_INSERT_SIZE, MEDIAN_ABSOLUTE_DEVIATION, MIN_INSERT_SIZE, MAX_INSERT_SIZE, MEAN_INSERT_SIZE, STANDARD_DEVIATION, READ_PAIRS, PAIR_ORIENTATION, WIDTH_OF_10_PERCENT, WIDTH_OF_20_PERCENT, WIDTH_OF_30_PERCENT, WIDTH_OF_40_PERCENT, WIDTH_OF_50_PERCENT, WIDTH_OF_60_PERCENT, WIDTH_OF_70_PERCENT, WIDTH_OF_80_PERCENT, WIDTH_OF_90_PERCENT, WIDTH_OF_95_PERCENT, WIDTH_OF_99_PERCENT, SAMPLE, LIBRARY, READ_GROUP. Rows: orientation (FR, RF, TANDEM) times accumulation unit. `MINIMUM_PCT=0` emits empty RF and TANDEM rows of zeros.
- Histogram: `insert_size` by `<prefix>.fr_count`, `.rf_count`, `.tandem_count`. Prefix is `All_Reads`, or the sample, library or PU name.
- Traps:
  - Only second-of-pair records count. The record's own TLEN sign and mate start decide orientation.
  - TLEN==0, secondary, supplementary and mate-unmapped records are skipped (verified: a second-of-pair with TLEN 0 and a secondary record both dropped, READ_PAIRS fell from 5 to 3).
  - QC-fail reads and MAPQ are not filtered. Duplicates are excluded by default.
  - Mean and SD come from the histogram trimmed to median + DEVIATIONS * MAD. SD uses n-1.
  - Median and MAD average the two middle values on even counts. Mode ties go to the lowest key.
  - Window loop for WIDTH_OF_X uses `low >= min-1 || high <= max+1`, truncates the median to int, and reports `(int)(high-low)+1`.
  - `total * pct` uses the double of the float 0.05f.
  - If every orientation is dropped, Picard writes no file and logs a warning. Picard exits 1 without Rscript when a chart is requested or the shim is missing.
  - Open gaps in the repo: levels, histogram output, mate end from `MC` not used for forward reads. Verify the current state on main before listing gaps.

#### CollectAlignmentSummaryMetrics (L)
- Inputs: BAM/SAM/CRAM. Reference (+ .dict/.fai) only for PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_HQ_MEDIAN_MISMATCHES. Without it these print 0.
- Options: `METRIC_ACCUMULATION_LEVEL`, `MAX_INSERT_SIZE=100000`, `EXPECTED_PAIR_ORIENTATIONS=[FR]`, `ADAPTER_SEQUENCE` (6 Illumina strings), `IS_BISULFITE_SEQUENCED=false`, `COLLECT_ALIGNMENT_INFORMATION=true`. Hard-coded: MAPQ threshold 20, base quality threshold 20.
- Class `picard.analysis.AlignmentSummaryMetrics` (36 columns, the 35 value columns plus CATEGORY first): CATEGORY, TOTAL_READS, PF_READS, PCT_PF_READS, PF_NOISE_READS, PF_READS_ALIGNED, PCT_PF_READS_ALIGNED, PF_ALIGNED_BASES, PF_HQ_ALIGNED_READS, PF_HQ_ALIGNED_BASES, PF_HQ_ALIGNED_Q20_BASES, PF_HQ_MEDIAN_MISMATCHES, PF_MISMATCH_RATE, PF_HQ_ERROR_RATE, PF_INDEL_RATE, MEAN_READ_LENGTH, SD_READ_LENGTH, MEDIAN_READ_LENGTH, MAD_READ_LENGTH, MIN_READ_LENGTH, MAX_READ_LENGTH, MEAN_ALIGNED_READ_LENGTH, READS_ALIGNED_IN_PAIRS, PCT_READS_ALIGNED_IN_PAIRS, PF_READS_IMPROPER_PAIRS, PCT_PF_READS_IMPROPER_PAIRS, BAD_CYCLES, STRAND_BALANCE, PCT_CHIMERAS, PCT_ADAPTER, PCT_SOFTCLIP, PCT_HARDCLIP, AVG_POS_3PRIME_SOFTCLIP_LENGTH, SAMPLE, LIBRARY, READ_GROUP. Rows: CATEGORY (FIRST_OF_PAIR, SECOND_OF_PAIR, PAIR, UNPAIRED) times unit.
- Histogram: READ_LENGTH by PAIRED_TOTAL_LENGTH_COUNT, PAIRED_ALIGNED_LENGTH_COUNT, UNPAIRED_TOTAL_LENGTH_COUNT, UNPAIRED_ALIGNED_LENGTH_COUNT. Always written.
- Traps:
  - PF_READS == 0 with TOTAL_READS > 0 throws. No row is written.
  - Rows: FIRST/SECOND/PAIR only if FIRST_OF_PAIR.TOTAL_READS > 0. UNPAIRED if it has reads or there are no first-of-pair reads. Order FIRST, SECOND, PAIR, UNPAIRED.
  - SD_READ_LENGTH prints `?` (NaN) for a single-read row (verified: UNPAIRED row with TOTAL_READS 1).
  - Unmapped reads land in the UNPAIRED_ALIGNED_LENGTH_COUNT histogram bin at key 0 (verified: key 0 gives UNPAIRED_ALIGNED_LENGTH_COUNT 1; key 10 gives UNPAIRED_TOTAL_LENGTH_COUNT 1).
  - Adapter test only for unmapped reads and MAPQ 0 reads. Match is any 16-mer of the adapters (and reverse complements) at read start, with at most 1 mismatch.
  - MEAN_ALIGNED_READ_LENGTH averages over all PF reads, unmapped reads count as 0.
  - PCT_SOFTCLIP and PCT_HARDCLIP divide by the PF base total from the read-length histogram. AVG_POS_3PRIME_SOFTCLIP_LENGTH divides by reads with a 3' soft clip, 0 on zero denominator.
  - 3' soft clip counts only clips after the first non-clip element in strand-oriented CIGAR.
  - PF_HQ_ALIGNED_READS needs MAPQ >= 20 and not QC-fail. MAPQ 0 reads still add PF_ALIGNED_BASES.
  - No reference means mismatch is false for every base.
  - Mismatch histograms include zero-mismatch reads. PF_HQ_MEDIAN_MISMATCHES uses the Histogram median.
  - Chimera: mate-mapped path uses MQ (absent MQ counts toward the denominator) and `ChimeraUtil.isChimeric` (TLEN above max, different contig, unexpected orientation, or SA tag). Orientation uses the record's own TLEN.
  - BAD_CYCLES: cycles where N or mismatch fraction of TOTAL_READS is at least 0.8. PAIR BAD_CYCLES = first + second.
  - PF_NOISE_READS tests XN equals integer 1.
  - Implementation guide section 2.4 omits BAD_CYCLES, PCT_SOFTCLIP, PCT_HARDCLIP, AVG_POS_3PRIME_SOFTCLIP_LENGTH, MEAN_ALIGNED_READ_LENGTH and the aligned-length histogram. Update the guide.

#### CollectQualityYieldMetrics (S)
- Options: `USE_ORIGINAL_QUALITIES=true` (OQ if present, else QUAL, per read), `INCLUDE_SECONDARY_ALIGNMENTS=false`, `INCLUDE_SUPPLEMENTAL_ALIGNMENTS=false`. `FLOW_MODE=true` throws.
- Class `CollectQualityYieldMetrics$QualityYieldMetrics` (11 columns): TOTAL_READS, PF_READS, READ_LENGTH, TOTAL_BASES, PF_BASES, Q20_BASES, PF_Q20_BASES, Q30_BASES, PF_Q30_BASES, Q20_EQUIVALENT_YIELD, PF_Q20_EQUIVALENT_YIELD. One row, no levels.
- Traps: yield = sum of qualities / 20 with long truncation. READ_LENGTH = TOTAL_BASES / TOTAL_READS as integer division. Unmapped, duplicate and QC-fail reads count in TOTAL. Hard-clipped bases are not counted.

#### QualityScoreDistribution (S)
- Options: `ALIGNED_READS_ONLY`, `PF_READS_ONLY`, `INCLUDE_NO_CALLS=false`.
- No metrics class. Histogram `java.lang.Byte`: QUALITY by COUNT_OF_Q, COUNT_OF_OQ. OQ column only if any read has OQ.
- Traps: N bases dropped from both counts. OQ counted only for reads that carry OQ. Row set is the union of qualities with a nonzero count. No METRICS CLASS block (verified).

#### MeanQualityByCycle (S)
- Options: `ALIGNED_READS_ONLY`, `PF_READS_ONLY`.
- Histogram `java.lang.Integer`: CYCLE by MEAN_QUALITY, MEAN_ORIGINAL_QUALITY (OQ column only if any read has OQ; missing cells print 0). No METRICS CLASS block (verified).
- Traps: second-of-pair cycles are offset by the highest first-read cycle with a positive quality sum. An all-Q0 cycle disappears and shifts the offset. Reverse-strand reads are flipped. Unpaired reads count as first read.

#### CollectBaseDistributionByCycle (S)
- Options: `ALIGNED_READS_ONLY`, `PF_READS_ONLY`.
- Class `picard.analysis.BaseDistributionByCycleMetrics` (7 columns): READ_END, CYCLE, PCT_A, PCT_C, PCT_G, PCT_T, PCT_N.
- Traps: any non-ACGT symbol counts as N. PCT is 0 to 100. READ_END 2 cycles are offset by the last nonzero first-read cycle. Zero-count cycles are omitted.

#### CollectDuplicateMetrics (S)
- Options: the metrics output option is `M=` (`METRICS_FILE`), not `O=`. Other options (`ASSUME_SORTED`, `STOP_AFTER`) come from `tools/CollectDuplicateMetrics.yaml`. Reads flags only.
- Class `picard.sam.DuplicationMetrics` (10 columns): LIBRARY, UNPAIRED_READS_EXAMINED, READ_PAIRS_EXAMINED, SECONDARY_OR_SUPPLEMENTARY_RDS, UNMAPPED_READS, UNPAIRED_READ_DUPLICATES, READ_PAIR_DUPLICATES, READ_PAIR_OPTICAL_DUPLICATES, PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE. One row per library.
- Traps: pair counts are accumulated per read then halved. OPTICAL is always 0, so ESTIMATED_LIBRARY_SIZE differs from MarkDuplicates. Library from RG LB, `Unknown` when absent. Library size empty when there are no duplicates. Needs the Lander-Waterman bisection solver.

#### EstimateLibraryComplexity (M)
- Inputs: one or more BAMs, mapped or unmapped.
- Options: `MIN_IDENTICAL_BASES=5`, `MAX_DIFF_RATE=0.03`, `MIN_MEAN_QUALITY=20`, `MAX_GROUP_RATIO=500`, `MAX_READ_LENGTH=0`, `MIN_GROUP_COUNT=2`, `OPTICAL_DUPLICATE_PIXEL_DISTANCE=100`, `READ_NAME_REGEX`, barcode tags.
- Class `DuplicationMetrics` (same columns as above). Histogram `duplication_group_count` (Integer), one value column per library.
- Traps: only histogram bins with at least MIN_GROUP_COUNT groups enter READ_PAIRS_EXAMINED. Mean quality is integer division. Reads are reverse-complemented back to original orientation. Greedy clustering depends on the sort order of the pair comparator. Other metric columns are 0, library size may be empty. Optical duplicate counting needs M18.

#### CollectGcBiasMetrics (M)
- Inputs: BAM and reference (needed in practice).
- Options: `SCAN_WINDOW_SIZE=100`, `METRIC_ACCUMULATION_LEVEL`, `ALSO_IGNORE_DUPLICATES=false`, `IS_BISULFITE_SEQUENCED=false`, `MINIMUM_GENOME_FRACTION=1e-5`.
- Two files. `GcBiasDetailMetrics` (11 columns): ACCUMULATION_LEVEL, READS_USED, GC, WINDOWS, READ_STARTS, MEAN_BASE_QUALITY, NORMALIZED_COVERAGE, ERROR_BAR_WIDTH, SAMPLE, LIBRARY, READ_GROUP. 101 rows per level and READS_USED value, zero rows kept. `GcBiasSummaryMetrics` (15 columns): ACCUMULATION_LEVEL, READS_USED, WINDOW_SIZE, TOTAL_CLUSTERS, ALIGNED_READS, AT_DROPOUT, GC_DROPOUT, GC_NC_0_19, GC_NC_20_39, GC_NC_40_59, GC_NC_60_79, GC_NC_80_100, SAMPLE, LIBRARY, READ_GROUP.
- Level label: ACCUMULATION_LEVEL prints `All Reads` (with a space, not `ALL_READS`) and READS_USED prints `ALL` (verified: summary row `All Reads  ALL  100  6  10 ...`, detail file has 101 `All Reads` rows). A comparator key built from `ALL_READS` would mismatch every row. Registered as a quirk (Q-020).
- Traps: window loop starts at index 1 and stops before refLen - windowSize, so edge windows are skipped. Windows with 5 or more N are excluded. GC% truncates by integer division and N counts in the denominator. Reverse reads use alignment end minus window size. READ_GROUP keyed by PU, two RGs without PU collapse into one `unknown` row (unconfirmed, Q-001). UNIQUE rows exist only with ALSO_IGNORE_DUPLICATES. TOTAL_CLUSTERS counts unpaired and first-of-pair reads, mapped or not (verified: 5 first-of-pair plus 1 unmapped unpaired gives 6, ALIGNED_READS 10). Open: secondary and supplementary handling and the MINIMUM_GENOME_FRACTION use were not confirmed. Read `CollectGcBiasMetrics.java` and `GcBiasMetricsCollector` before coding.

#### CollectWgsMetrics (L), CollectRawWgsMetrics (S)
- Inputs: coordinate-sorted BAM/CRAM, reference with .fai and .dict.
- Options (WGS defaults / Raw defaults): `MINIMUM_MAPPING_QUALITY` 20 / 0, `MINIMUM_BASE_QUALITY` 20 / 3, `COVERAGE_CAP` 250 / 100000, `LOCUS_ACCUMULATION_CAP` 100000 / 200000, `COUNT_UNPAIRED=false`, `INTERVALS`, `READ_LENGTH=150` (theoretical sensitivity only), `SAMPLE_SIZE=10000`, `ALLELE_FRACTION` list, `INCLUDE_BQ_HISTOGRAM=false`, `USE_FAST_ALGORITHM=false`, `THEORETICAL_SENSITIVITY_OUTPUT`. WgsMetrics has `STOP_AFTER` but no `ASSUME_SORTED`.
- Class `WgsMetrics` (Raw: `CollectRawWgsMetrics$RawWgsMetrics`), 32 columns: GENOME_TERRITORY, MEAN_COVERAGE, SD_COVERAGE, MEDIAN_COVERAGE, MAD_COVERAGE, PCT_EXC_ADAPTER, PCT_EXC_MAPQ, PCT_EXC_DUPE, PCT_EXC_UNPAIRED, PCT_EXC_BASEQ, PCT_EXC_OVERLAP, PCT_EXC_CAPPED, PCT_EXC_TOTAL, PCT_1X, PCT_5X, PCT_10X, PCT_15X, PCT_20X, PCT_25X, PCT_30X, PCT_40X, PCT_50X, PCT_60X, PCT_70X, PCT_80X, PCT_90X, PCT_100X, FOLD_80_BASE_PENALTY, FOLD_90_BASE_PENALTY, FOLD_95_BASE_PENALTY, HET_SNP_SENSITIVITY, HET_SNP_Q. One row.
- Histogram: `coverage` by `high_quality_coverage_count`, `unfiltered_baseq_count` (keys 0 to cap; third column populated only with INCLUDE_BQ_HISTOGRAM). Raw writes about 100k lines.
- Other output: theoretical sensitivity file (ALLELE_FRACTION, THEORETICAL_SENSITIVITY, THEORETICAL_SENSITIVITY_Q, LOG_ODDS_THRESHOLD, SAMPLE_SIZE).
- Traps:
  - PCT_EXC_* denominator is high-quality bases plus all excluded bases, not GENOME_TERRITORY.
  - MEAN_COVERAGE divides by GENOME_TERRITORY including zero-depth positions.
  - Exclusion order decides which reason a base is charged to. Overlap rule drops the later observation of a template, so read order at a locus matters.
  - FOLD_xx_BASE_PENALTY prints `?` for NaN.
  - Depth counts aligned M/=/X bases only.
  - LOCUS_ACCUMULATION_CAP drops reads silently and they are not counted as capped.
  - HET_SNP_SENSITIVITY and HET_SNP_Q use floating-point convolution. Last-digit drift is likely. Needs a registered tolerance or exact port of the arithmetic order. Raw took 0.76 min on a 3 kb genome because of the sensitivity step.
  - CRAM and USE_FAST_ALGORITHM=true take other code paths. SD_COVERAGE denominator (n vs n-1) must be confirmed on the corpus.
  - Unfiltered histograms skip quality-2 bases.
  - Raw class line must read `picard.analysis.CollectRawWgsMetrics$RawWgsMetrics`.

#### CollectHsMetrics (XL)
- Inputs: BAM, `BAIT_INTERVALS`, `TARGET_INTERVALS` (interval_list with SQ header). Reference optional, needed for `%gc` and dropout.
- Options: `MINIMUM_MAPPING_QUALITY=20`, `MINIMUM_BASE_QUALITY=20`, `CLIP_OVERLAPPING_READS=true`, `NEAR_DISTANCE=250`, `INCLUDE_INDELS=false`, `COVERAGE_CAP=200`, `SAMPLE_SIZE=10000`, `ALLELE_FRACTION`, `METRIC_ACCUMULATION_LEVEL`, `BAIT_SET_NAME`, `PER_TARGET_COVERAGE`, `PER_BASE_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT`. CollectHsMetrics has neither `STOP_AFTER` nor `ASSUME_SORTED` (verified from `-h`).
- Class `picard.analysis.directed.HsMetrics`: BAIT_SET, BAIT_TERRITORY, BAIT_DESIGN_EFFICIENCY, ON_BAIT_BASES, NEAR_BAIT_BASES, OFF_BAIT_BASES, PCT_SELECTED_BASES, PCT_OFF_BAIT, ON_BAIT_VS_SELECTED, MEAN_BAIT_COVERAGE, PCT_USABLE_BASES_ON_BAIT, PCT_USABLE_BASES_ON_TARGET, FOLD_ENRICHMENT, HS_LIBRARY_SIZE, HS_PENALTY_10X, HS_PENALTY_20X, HS_PENALTY_30X, HS_PENALTY_40X, HS_PENALTY_50X, HS_PENALTY_100X, TARGET_TERRITORY, GENOME_SIZE, TOTAL_READS, PF_READS, PF_BASES, PF_UNIQUE_READS, PF_UQ_READS_ALIGNED, PF_BASES_ALIGNED, PF_UQ_BASES_ALIGNED, ON_TARGET_BASES, PCT_PF_READS, PCT_PF_UQ_READS, PCT_PF_UQ_READS_ALIGNED, MEAN_TARGET_COVERAGE, MEDIAN_TARGET_COVERAGE, MAX_TARGET_COVERAGE, MIN_TARGET_COVERAGE, ZERO_CVG_TARGETS_PCT, PCT_EXC_DUPE, PCT_EXC_ADAPTER, PCT_EXC_MAPQ, PCT_EXC_BASEQ, PCT_EXC_OVERLAP, PCT_EXC_OFF_TARGET, FOLD_80_BASE_PENALTY, PCT_TARGET_BASES_1X, _2X, _10X, _20X, _30X, _40X, _50X, _100X, _250X, _500X, _1000X, _2500X, _5000X, _10000X, _25000X, _50000X, _100000X (16 columns, each prefixed `PCT_TARGET_BASES`), AT_DROPOUT, GC_DROPOUT, HET_SNP_SENSITIVITY, HET_SNP_Q, SAMPLE, LIBRARY, READ_GROUP.
- Histogram: `coverage_or_base_quality` by `high_quality_coverage_count`, `unfiltered_baseq_count`.
- Other outputs: PER_TARGET_COVERAGE tsv (chrom, start, end, length, name, %gc, mean_coverage, normalized_coverage, min_normalized_coverage, max_normalized_coverage, min_coverage, max_coverage, pct_0x, read_count). PER_BASE_COVERAGE tsv (chrom, pos, target, coverage). Theoretical sensitivity tsv.
- Traps:
  - Territory uses merged intervals, but ZERO_CVG_TARGETS_PCT uses the unmerged target list.
  - GENOME_SIZE is the sum of SQ lengths in the bait file header.
  - BAIT_SET is the bait file basename without extension.
  - ON/NEAR/OFF bait counts are taken before dupe, adapter, MAPQ filters and overlap clipping. ON_TARGET_BASES is after them.
  - PCT_EXC_* use PF_BASES_ALIGNED as denominator (WGS uses another). PCT_EXC_OFF_TARGET counts only high-quality aligned bases not on target.
  - Supplementary adds bases but not reads. PF_UQ_READS_ALIGNED excludes duplicates and secondary.
  - HS_LIBRARY_SIZE uses the Lander-Waterman solver. HS_PENALTY prints -1 when not computable. FOLD_80 prints `?` for NaN.
  - AT_DROPOUT and GC_DROPOUT use 101 integer GC bins and need the reference.
  - Rows order ALL_READS, SAMPLE, LIBRARY, READ_GROUP. INCLUDE_INDELS changes the on-target definition.

#### CollectRnaSeqMetrics (L)
- Inputs: BAM, `REF_FLAT` (UCSC refFlat), `STRAND_SPECIFICITY` (NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND, no default). Optional but number-changing: `RIBOSOMAL_INTERVALS`.
- Options: `END_BIAS_BASES=100`, `MINIMUM_LENGTH=500`, `RRNA_FRAGMENT_PERCENTAGE=0.8`, `IGNORE_SEQUENCE`, `METRIC_ACCUMULATION_LEVEL`. Without Rscript the real tool writes nothing when a chart is requested.
- Class `picard.analysis.RnaSeqMetrics`: PF_BASES, PF_ALIGNED_BASES, RIBOSOMAL_BASES, CODING_BASES, UTR_BASES, INTRONIC_BASES, INTERGENIC_BASES, IGNORED_READS, CORRECT_STRAND_READS, INCORRECT_STRAND_READS, NUM_R1_TRANSCRIPT_STRAND_READS, NUM_R2_TRANSCRIPT_STRAND_READS, NUM_UNEXPLAINED_READS, PCT_R1_TRANSCRIPT_STRAND_READS, PCT_R2_TRANSCRIPT_STRAND_READS, PCT_RIBOSOMAL_BASES, PCT_CODING_BASES, PCT_UTR_BASES, PCT_INTRONIC_BASES, PCT_INTERGENIC_BASES, PCT_MRNA_BASES, PCT_USABLE_BASES, PCT_CORRECT_STRAND_READS, MEDIAN_CV_COVERAGE, MEDIAN_5PRIME_BIAS, MEDIAN_3PRIME_BIAS, MEDIAN_5PRIME_TO_3PRIME_BIAS, SAMPLE, LIBRARY, READ_GROUP.
- Histogram: `normalized_position` (0 to 100, 101 rows) with one `<prefix>normalized_coverage` column per level.
- Traps:
  - RIBOSOMAL_BASES and PCT_RIBOSOMAL_BASES are blank, not 0, without RIBOSOMAL_INTERVALS.
  - PCT_USABLE_BASES divides by PF_BASES, other PCT_ base fractions by PF_ALIGNED_BASES.
  - rRNA test: largest single interval overlap over fragment length, at least 0.8, removes the whole record from every other counter.
  - Reads overlapping two genes classify bases but are excluded from strand counts.
  - Mate end comes from MC, else read length.
  - Transcript pick per gene: highest mean coverage, ties to first seen, Gene transcript order is a HashMap. Top-1000 cutoff index is n-1001.
  - With STRAND=NONE, CORRECT and INCORRECT stay 0.
  - Duplicate flag and MAPQ are not filtered. refFlat is 0-based half-open.
  - READ_GROUP without PU gives one `unknown` row plus an all-zero row labelled with the header read group (unconfirmed, Q-001).

#### CollectMultipleMetrics (M)
- Dispatcher. Options: `PROGRAM` (default five, verified: CollectAlignmentSummaryMetrics, CollectBaseDistributionByCycle, CollectInsertSizeMetrics, MeanQualityByCycle, QualityScoreDistribution; also CollectGcBiasMetrics, RnaSeqMetrics, CollectSequencingArtifactMetrics, CollectQualityYieldMetrics), `METRIC_ACCUMULATION_LEVEL`, `EXTRA_ARGUMENT=<PROGRAM>::<arg>`, `INTERVALS`, `DB_SNP`, `INCLUDE_UNPAIRED`, `IGNORE_SEQUENCE`, `FILE_EXTENSION`, `DO_NOT_CREATE_PLOTS`.
- Output files by fixed extension: `.alignment_summary_metrics`, `.insert_size_metrics`, `.base_distribution_by_cycle_metrics`, `.quality_by_cycle_metrics`, `.quality_distribution_metrics`, `.quality_yield_metrics`, `.gc_bias.detail_metrics`, `.gc_bias.summary_metrics`, `.rna_metrics`, and the five artifact files.
- Traps: file names must match exactly. Legacy `NAME::arg=value` form works in 3.5.0, `--NAME value` fails to parse. `METRIC_ACCUMULATION_LEVEL` applies to all programs at once. INPUT, REFERENCE_SEQUENCE, ASSUME_SORTED and STOP_AFTER are silently ignored in EXTRA_ARGUMENT. STOP_AFTER is int here. Insert size on tiny inputs can error through MINIMUM_PCT.

#### CompareMetrics (S)
- Inputs: two metrics files. Options: `KEY` (repeatable), `METRIC_ALLOWABLE_RELATIVE_CHANGE=NAME:value`, `METRICS_TO_IGNORE`, `METRICS_NOT_REQUIRED`, `IGNORE_HISTOGRAM_DIFFERENCES`, `OUTPUT`, `OUTPUT_TABLE`. Exit 0 equal, 1 different.
- Class `CompareMetrics$MetricComparisonDifferences`: KEY, METRIC, VALUE1, VALUE2.
- Traps: without KEY the row label is the zero-based index. Zero value1 gives max-double relative change. Tolerance is strict less-than. Both-NaN is equal. Histogram equality is all or nothing. Row-count mismatch returns before any class check.
- Role here: implement the same semantics in the harness (section 5) and optionally as a CLI. It is not a user feature first.

#### CollectVariantCallingMetrics (M)
- Inputs: VCF (indexed if TARGET_INTERVALS or THREAD_COUNT > 1), `DBSNP` VCF (bgzip and tabix worked), sequence dictionary. Options: `GVCF_INPUT=false`, `TARGET_INTERVALS`, `THREAD_COUNT=1`.
- Detail class `VariantCallingDetailMetrics` (file `.variant_calling_detail_metrics`): SAMPLE_ALIAS, HET_HOMVAR_RATIO, PCT_GQ0_VARIANTS, TOTAL_GQ0_VARIANTS, TOTAL_HET_DEPTH, TOTAL_SNPS, NUM_IN_DB_SNP, NOVEL_SNPS, FILTERED_SNPS, PCT_DBSNP, DBSNP_TITV, NOVEL_TITV, TOTAL_INDELS, NOVEL_INDELS, FILTERED_INDELS, PCT_DBSNP_INDELS, NUM_IN_DB_SNP_INDELS, DBSNP_INS_DEL_RATIO, NOVEL_INS_DEL_RATIO, TOTAL_MULTIALLELIC_SNPS, NUM_IN_DB_SNP_MULTIALLELIC, TOTAL_COMPLEX_INDELS, NUM_IN_DB_SNP_COMPLEX_INDELS, SNP_REFERENCE_BIAS, NUM_SINGLETONS.
- Summary class `VariantCallingSummaryMetrics`: TOTAL_SNPS, NUM_IN_DB_SNP, NOVEL_SNPS, FILTERED_SNPS, PCT_DBSNP, DBSNP_TITV, NOVEL_TITV, TOTAL_INDELS, NOVEL_INDELS, FILTERED_INDELS, PCT_DBSNP_INDELS, NUM_IN_DB_SNP_INDELS, DBSNP_INS_DEL_RATIO, NOVEL_INS_DEL_RATIO, TOTAL_MULTIALLELIC_SNPS, NUM_IN_DB_SNP_MULTIALLELIC, TOTAL_COMPLEX_INDELS, NUM_IN_DB_SNP_COMPLEX_INDELS, SNP_REFERENCE_BIAS, NUM_SINGLETONS. One row.
- Traps: PCT_DBSNP and PCT_DBSNP_INDELS are computed in float, other ratios in double. TITV and INS_DEL stay 0, not NaN, on zero denominator. HET_HOMVAR_RATIO and PCT_GQ0_VARIANTS print `?` on zero. dbSNP match is by contig and start only, SNP and indel in separate bitsets. MNP and symbolic sites fall in no category. GQ0 check does not exclude no-calls. Detail FILTERED_* is per sample, summary per site. The variant-type classifier (SNP, MNP, indel, complex, symbolic, mixed) is the main risk and GenotypeConcordance reuses it.

#### GenotypeConcordance (L)
- Inputs: `TRUTH_VCF`, `CALL_VCF`, sample names. Options: `MIN_GQ=0`, `MIN_DP=0`, `IGNORE_FILTER_STATUS`, `MISSING_SITES_HOM_REF` (needs intervals), `INTERVALS`, `INTERSECT_INTERVALS=true`, `OUTPUT_ALL_ROWS`, `OUTPUT_VCF`, `USE_VCF_INDEX`.
- Summary class: VARIANT_TYPE, TRUTH_SAMPLE, CALL_SAMPLE, HET_SENSITIVITY, HET_PPV, HET_SPECIFICITY, HOMVAR_SENSITIVITY, HOMVAR_PPV, HOMVAR_SPECIFICITY, VAR_SENSITIVITY, VAR_PPV, VAR_SPECIFICITY, GENOTYPE_CONCORDANCE, NON_REF_GENOTYPE_CONCORDANCE. Rows SNP, INDEL.
- Detail class: VARIANT_TYPE, TRUTH_SAMPLE, CALL_SAMPLE, TRUTH_STATE, CALL_STATE, COUNT, CONTINGENCY_VALUES (non-zero rows unless OUTPUT_ALL_ROWS).
- Contingency class: VARIANT_TYPE, TRUTH_SAMPLE, CALL_SAMPLE, TP_COUNT, TN_COUNT, FP_COUNT, FN_COUNT, EMPTY_COUNT.
- Optional output: annotated VCF with INFO CONC_ST, samples named `truth` and `call`.
- Traps: one scheme cell can add to two categories (for example FP,TN). Zero denominators print `?`. Filtered truth is always treated as missing. LOW_GQ beats LOW_DP, absent GQ/DP (-1) never triggers. One-allele genotypes become NO_CALL. MISSING_SITES_HOM_REF counts every uncovered interval base as TN. Indel allele normalization decides matches. About 60 scheme cells (GA4GH default).

### 2.2 Later tools

Columns for these tools are verified by run, except where marked. Entries are shorter because they are not 1.0 work.

| Tool | Inputs and options that change numbers | Class and columns | Histograms | Main traps | Cx |
|---|---|---|---|---|---|
| CollectQualityYieldMetricsSNVQ | Tags qa, qc, qg, qt on every non-excluded read, else NPE. `ALTERNATE_QUALITY_ATTRIBUTE`, `INCLUDE_BQ_HISTOGRAM` | `...SNVQ$QualityYieldMetrics`: TOTAL_READS, PF_READS, READ_LENGTH, TOTAL_BASES, PF_BASES, Q20_BASES, PF_Q20_BASES, Q30_BASES, PF_Q30_BASES, Q40_BASES, PF_Q40_BASES, PCT_Q20_BASES, PCT_Q30_BASES, PCT_Q40_BASES, PCT_PF_Q20_BASES, PCT_PF_Q30_BASES, PCT_PF_Q40_BASES, TOTAL_SNVQ, PF_SNVQ, Q20_SNVQ, PF_Q20_SNVQ, Q30_SNVQ, PF_Q30_SNVQ, Q40_SNVQ, PF_Q40_SNVQ, PCT_Q20_SNVQ, PCT_Q30_SNVQ, PCT_Q40_SNVQ, PCT_PF_Q20_SNVQ, PCT_PF_Q30_SNVQ, PCT_PF_Q40_SNVQ | With BQ option: KEY by BQ_COUNT, SNVQ_COUNT, READ_INDEX_MEAN_BQ, READ_INDEX_MEAN_SNVQ | N or lowercase base gives 4 SNVQ observations. Q40 cascades into Q20 and Q30. Missing tag crashes: decide error vs mimic. No Q20_EQUIVALENT_YIELD | S |
| CollectWgsMetricsWithNonZeroCoverage | As WGS plus `CHART_OUTPUT` (run aborts without it and without Rscript) | `...$WgsMetricsWithNonZeroCoverage`: CATEGORY + the 32 WgsMetrics columns. Rows WHOLE_GENOME, NON_ZERO_REGIONS | `coverage`, `count_WHOLE_GENOME`, `count_NON_ZERO_REGIONS` | PCT_EXC_* copied unscaled to the second row. Each row has its own NaN rule. HET_SNP differs per row | S |
| CollectTargetedPcrMetrics | `AMPLICON_INTERVALS`, `TARGET_INTERVALS`. Defaults MAPQ 1, BQ 0, `CLIP_OVERLAPPING_READS=false`, `CUSTOM_AMPLICON_SET_NAME` | `TargetedPcrMetrics`: CUSTOM_AMPLICON_SET, AMPLICON_TERRITORY, ON_AMPLICON_BASES, NEAR_AMPLICON_BASES, OFF_AMPLICON_BASES, PCT_AMPLIFIED_BASES, PCT_OFF_AMPLICON, ON_AMPLICON_VS_SELECTED, MEAN_AMPLICON_COVERAGE, FOLD_ENRICHMENT, PF_SELECTED_PAIRS, PF_SELECTED_UNIQUE_PAIRS, ON_TARGET_FROM_PAIR_BASES, then TARGET_TERRITORY through HET_SNP_Q, SAMPLE, LIBRARY, READ_GROUP exactly as Hs (no library size, no HS_PENALTY, no BAIT_* columns) | As Hs | Per-tool defaults on a shared collector. PF_SELECTED_PAIRS counts first-of-pair with mapped mate on an amplicon window | M |
| MarkDuplicates | Coordinate-sorted BAM. `DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES`, `OPTICAL_DUPLICATE_PIXEL_DISTANCE=100`, `READ_NAME_REGEX`, `REMOVE_DUPLICATES`, `TAGGING_POLICY`, `TAG_DUPLICATE_SET_MEMBERS`, barcode tags, `MAX_OPTICAL_DUPLICATE_SET_SIZE=300000`, `CLEAR_DT` | `DuplicationMetrics` (10 columns, section 2.1) | `set_size` (Double keys) by all_sets, optical_sets, non_optical_sets | Pair counts halved. PERCENT_DUPLICATION = (unpaired dups + 2 pair dups) / (unpaired + 2 pairs). Score ties go to first ReadEnds. 5' unclipped position needs clips and MC. Writes a BAM with 0x400, PG record, optional DT/DI/DS tags. Optical sets split by strand | XL |
| UmiAwareMarkDuplicatesWithMateCigar | RX on every read, MC, `MAX_EDIT_DISTANCE_TO_JOIN=1`, `UMI_TAG_NAME`, `DUPLEX_UMI`, `ALLOW_MISSING_UMIS`, `UMI_METRICS_FILE` | `DuplicationMetrics` plus `UmiMetrics`: LIBRARY, MEAN_UMI_LENGTH, OBSERVED_UNIQUE_UMIS, INFERRED_UNIQUE_UMIS, OBSERVED_BASE_ERRORS, DUPLICATE_SETS_IGNORING_UMI, DUPLICATE_SETS_WITH_UMI, OBSERVED_UMI_ENTROPY, INFERRED_UMI_ENTROPY, UMI_BASE_QUALITIES, PCT_UMI_WITH_N | `set_size` as above | Picard throws "UMIs of differing lengths" with two libraries. Throws if any read lacks RX. Graph order affects inferred UMI | L |
| CheckDuplicateMarking | `MODE` (ALL, PRIMARY_ONLY, PRIMARY_MAPPED_ONLY, PRIMARY_PROPER_PAIR_ONLY), `OUTPUT` | none | none | Exit code 1 and ERROR line are the output. Queryname sort | S |
| CollectUmiPrevalenceMetrics | Indexed BAM. `MINIMUM_MQ=30`, `BARCODE_TAG=RX`, `BARCODE_BQ=BQ`, `MINIMUM_BARCODE_BQ=30`, `FILTER_UNPAIRED_READS` | none | numUmis by duplicateSets (Integer) | Barcode quality filter is inverted in source: keeps low or missing BQ. Copy as is | M |
| CollectSequencingArtifactMetrics | Reference. `CONTEXT_SIZE=1`, `MINIMUM_QUALITY_SCORE=20`, `MINIMUM_MAPPING_QUALITY=30`, `MINIMUM_INSERT_SIZE=60`, `MAXIMUM_INSERT_SIZE=600`, `INCLUDE_DUPLICATES`, `INCLUDE_NON_PF_READS`, `INCLUDE_UNPAIRED`, `TANDEM_READS`, `USE_OQ=true`, `DB_SNP`, `INTERVALS`, `CONTEXTS_TO_PRINT` | 5 files. `PreAdapterDetailMetrics`: SAMPLE_ALIAS, LIBRARY, REF_BASE, ALT_BASE, CONTEXT, PRO_REF_BASES, PRO_ALT_BASES, CON_REF_BASES, CON_ALT_BASES, ERROR_RATE, QSCORE. `PreAdapterSummaryMetrics` and `BaitBiasSummaryMetrics`: SAMPLE_ALIAS, LIBRARY, REF_BASE, ALT_BASE, TOTAL_QSCORE, WORST_CXT, WORST_CXT_QSCORE, WORST_PRE_CXT, WORST_PRE_CXT_QSCORE, WORST_POST_CXT, WORST_POST_CXT_QSCORE, ARTIFACT_NAME. `BaitBiasDetailMetrics`: SAMPLE_ALIAS, LIBRARY, REF_BASE, ALT_BASE, CONTEXT, FWD_CXT_REF_BASES, FWD_CXT_ALT_BASES, REV_CXT_REF_BASES, REV_CXT_ALT_BASES, FWD_ERROR_RATE, REV_ERROR_RATE, ERROR_RATE, QSCORE. `ErrorSummaryMetrics`: REF_BASE, ALT_BASE, SUBSTITUTION, REF_COUNT, ALT_COUNT, SUBSTITUTION_RATE | none | Pro/con rules live in `ArtifactCounter.java` (not read). QSCORE caps at 100, error floor 1e-10, rounding rule unconfirmed. Library order from a HashMap (observed L1, L2). Row sort order unconfirmed. Missing sample or library becomes UnknownSample or UnknownLibrary | L |
| CollectOxoGMetrics | Reference. `CONTEXT_SIZE`, `CONTEXTS`, `DB_SNP`, `INTERVALS`, `INCLUDE_NON_PF_READS=true`, insert 60 to 600, MAPQ 30, BQ 20, `USE_OQ` | `CollectOxoGMetrics$CpcgMetrics`: SAMPLE_ALIAS, LIBRARY, CONTEXT, TOTAL_SITES, TOTAL_BASES, REF_NONOXO_BASES, REF_OXO_BASES, REF_TOTAL_BASES, ALT_NONOXO_BASES, ALT_OXO_BASES, OXIDATION_ERROR_RATE, OXIDATION_Q, C_REF_REF_BASES, G_REF_REF_BASES, C_REF_ALT_BASES, G_REF_ALT_BASES, C_REF_OXO_ERROR_RATE, C_REF_OXO_Q, G_REF_OXO_ERROR_RATE, G_REF_OXO_Q | none | Row order follows Java HashSet iteration and interleaves libraries by context. No zero-denominator guard on C_REF/G_REF rates. SAMPLE_ALIAS for multi-sample unverified | M |
| ConvertSequencingArtifactToOxoG | `INPUT_BASE` or `PRE_ADAPTER_IN` + `BAIT_BIAS_IN`, reference | `CpcgMetrics` (same 20 columns) | none | TOTAL_SITES always 0. Use counts, not rounded rates. Order is library then context. Field mapping unread: read the source | S |
| CollectSamErrorMetrics | Reference, indexed VCF, `ERROR_METRICS` (27 defaults), `MIN_BASE_Q`, `MIN_MAPPING_Q`, `PRIOR_Q=30`, `PROBABILITY`, `MAX_LOCI`, `LONG_HOMOPOLYMER`, `LOCATION_BIN_SIZE` | Three classes. `BaseErrorMetric`: ERROR_BASES, Q_SCORE, COVARIATE, TOTAL_BASES. `OverlappingErrorMetric`: NUM_BASES_WITH_OVERLAPPING_READS, NUM_DISAGREES_WITH_REFERENCE_ONLY, DISAGREES_WITH_REFERENCE_ONLY_Q, NUM_DISAGREES_WITH_REF_AND_MATE, DISAGREES_WITH_REF_AND_MATE_ONLY_Q, NUM_THREE_WAYS_DISAGREEMENT, THREE_WAYS_DISAGREEMENT_ONLY_Q, COVARIATE, TOTAL_BASES. `IndelErrorMetric`: NUM_INSERTIONS, NUM_INSERTED_BASES, INSERTIONS_Q, NUM_DELETIONS, NUM_DELETED_BASES, DELETIONS_Q, ERROR_BASES, Q_SCORE, COVARIATE, TOTAL_BASES. One file per entry (26 files seen) | none | About 30 stratifiers. Q formula in `BaseCalculator` unread. Sampling uses Random(42). VCF masking skips all non-filtered records. If ever built: BASE_QUALITY, CYCLE, READ_ORDINALITY, ONE_BASE_PADDED_CONTEXT, GC_CONTENT first | XL |
| CollectRrbsMetrics | Reference. `C_QUALITY_THRESHOLD=20`, `NEXT_BASE_QUALITY_THRESHOLD=10`, `MINIMUM_READ_LENGTH=5`, `MAX_MISMATCH_RATE=0.1`, `SEQUENCE_NAMES`, levels | `RrbsSummaryMetrics`: READS_ALIGNED, NON_CPG_BASES, NON_CPG_CONVERTED_BASES, PCT_NON_CPG_BASES_CONVERTED, CPG_BASES_SEEN, CPG_BASES_CONVERTED, PCT_CPG_BASES_CONVERTED, MEAN_CPG_COVERAGE, MEDIAN_CPG_COVERAGE, READS_WITH_NO_CPG, READS_IGNORED_SHORT, READS_IGNORED_MISMATCHES, SAMPLE, LIBRARY, READ_GROUP. `RrbsCpgDetailMetrics`: SEQUENCE_NAME, POSITION, TOTAL_SITES, CONVERTED_SITES, PCT_CONVERTED, SAMPLE, LIBRARY, READ_GROUP | none | Default mismatch rate dropped 691 of 1600 bisulfite reads in one test. Extra all-zero row with `?` at some level combinations. Non-CpG counted even with no CpG | M |
| CollectJumpingLibraryMetrics | Paired reads only (unpaired throws). `MINIMUM_MAPPING_QUALITY=0`, `TAIL_LIMIT=10000`, `CHIMERA_KB_MIN=100000` | `JumpingLibraryMetrics`: JUMP_PAIRS, JUMP_DUPLICATE_PAIRS, JUMP_DUPLICATE_PCT, JUMP_LIBRARY_SIZE, JUMP_MEAN_INSERT_SIZE, JUMP_STDEV_INSERT_SIZE, NONJUMP_PAIRS, NONJUMP_DUPLICATE_PAIRS, NONJUMP_DUPLICATE_PCT, NONJUMP_LIBRARY_SIZE, NONJUMP_MEAN_INSERT_SIZE, NONJUMP_STDEV_INSERT_SIZE, CHIMERIC_PAIRS, FRAGMENTS, PCT_JUMPS, PCT_NONJUMPS, PCT_CHIMERAS | none written | Check order supersized, tandem, cross-chrom. Strict `>` against max(mode, 100000). Hidden SAMPLE_FOR_MODE=50000 | S |
| BamIndexStats | BAM with .bai | none. Stdout lines `chr1 length=<TAB>6000<TAB>Aligned= 185<TAB>Unaligned= 0`, then `NoCoordinateCount= 1` | none | Exact spacing. BAI only, no CSI or CRAM. Header contig order | S |
| CompareSAMs | Two files. `LENIENT_*`, `LOW_MQ_THRESHOLD=3`, `COMPARE_MQ` | `SamComparisonMetric`: LEFT_FILE, RIGHT_FILE, MAPPINGS_MATCH, MAPPINGS_DIFFER, UNMAPPED_BOTH, UNMAPPED_LEFT, UNMAPPED_RIGHT, MISSING_LEFT, MISSING_RIGHT, DUPLICATE_MARKINGS_DIFFER, ARE_EQUAL | With COMPARE_MQ: BIN by VALUE, string keys | Match is contig plus start only. Strand never counts (htsjdk compares a record with itself). Paths are absolute. STRICT aborts on odd mate fields | M |
| CalculateReadGroupChecksum | Header only | One-line file, 32 hex characters, no trailing newline | none | ID excluded. Read groups sorted PU, LB, DT, SM, CN, PL, DS, ID. No separators. Output beside the input | S |
| CheckTerminatorBlock | BGZF file | Prints HAS_TERMINATOR_BLOCK, HAS_HEALTHY_LAST_BLOCK or DEFECTIVE. Exit 100 on DEFECTIVE | none | Three outcomes. Tail check without full decompression | S |
| AccumulateVariantCallingMetrics | Shard prefixes | Same two classes as CollectVariantCallingMetrics | none | Lossy inversion of printed ratios (numHomVar = round(TOTAL_SNPS / (ratio + 1))). Reads rounded text. Row order from a HashMap (needs M22) | S |
| AccumulateQualityYieldMetrics | Several QualityYield files | Same class as QualityYield | none | READ_LENGTH integer division. Fails on mixed USE_ORIGINAL_QUALITIES | S |
| FindMendelianViolations | VCF with GT, GQ, DP, AD, PL, PED. `MIN_GQ=30`, `MIN_DP=0`, `MIN_HET_FRACTION=0.3`, `SKIP_CHROMS`, `FEMALE_CHROMS`, `MALE_CHROMS`, `PSEUDO_AUTOSOMAL_REGIONS` | `MendelianViolationMetrics`: FAMILY_ID, MOTHER, FATHER, OFFSPRING, OFFSPRING_SEX, NUM_VARIANT_SITES, NUM_DIPLOID_DENOVO, NUM_HOMVAR_HOMVAR_HET, NUM_HOMREF_HOMVAR_HOM, NUM_HOM_HET_HOM, NUM_HAPLOID_DENOVO, NUM_HAPLOID_OTHER, NUM_OTHER, TOTAL_MENDELIAN_VIOLATIONS | none | NUM_VARIANT_SITES incremented before the haploid-het skip. De novo GQ test uses PL[0], missing PL throws. Default PAR list mixes b37 and hg38. `VCF_DIR` threw "duplicate sample names" in 3.5.0 | M |
| CrosscheckFingerprints | BAM/VCF, haplotype map. `CROSSCHECK_BY`, `CROSSCHECK_MODE`, `LOD_THRESHOLD`, `CALCULATE_TUMOR_AWARE_RESULTS`, `LOSS_OF_HET_RATE`, `MAX_EFFECT_OF_EACH_HAPLOTYPE_BLOCK` | `CrosscheckMetric`: LEFT_GROUP_VALUE, RIGHT_GROUP_VALUE, RESULT, DATA_TYPE, LOD_SCORE, LOD_SCORE_TUMOR_NORMAL, LOD_SCORE_NORMAL_TUMOR, LEFT_RUN_BARCODE, LEFT_LANE, LEFT_MOLECULAR_BARCODE_SEQUENCE, LEFT_LIBRARY, LEFT_SAMPLE, LEFT_FILE, RIGHT_RUN_BARCODE, RIGHT_LANE, RIGHT_MOLECULAR_BARCODE_SEQUENCE, RIGHT_LIBRARY, RIGHT_SAMPLE, RIGHT_FILE | none | Group values are Java `String.hashCode` strings. Identity is PU, not RG ID. FILE columns are absolute URIs. Row order from hash maps. LOCUS_MAX_READS sampling uses a shared Random(42) | XL |
| CheckFingerprint | BAM/VCF, genotypes VCF, haplotype map | `FingerprintingSummaryMetrics`: READ_GROUP, SAMPLE, LL_EXPECTED_SAMPLE, LL_RANDOM_SAMPLE, LOD_EXPECTED_SAMPLE, HAPLOTYPES_WITH_GENOTYPES, HAPLOTYPES_CONFIDENTLY_CHECKED, HAPLOTYPES_CONFIDENTLY_MATCHING, HET_AS_HOM, HOM_AS_HET, HOM_AS_OTHER_HOM. `FingerprintingDetailMetrics`: READ_GROUP, SAMPLE, SNP, SNP_ALLELES, CHROM, POSITION, EXPECTED_GENOTYPE, OBSERVED_GENOTYPE, LOD, OBS_A, OBS_B | none | READ_GROUP is PU, blank when absent. Float-precision log10, summation order matters. Exit codes 1 and 2 | L |
| ExtractFingerprint | BAM, map, reference | VCF with AD:DP:PL | none | Subsampling via static Random(42). PL from capped-evidence model | M |
| CalculateFingerprintMetrics | BAM/VCF, map. `CALCULATE_BY`, `GENOTYPE_LOD_THRESHOLD=3.0`, `NUMBER_OF_SAMPLING=100` | `FingerprintMetrics`: SAMPLE_ALIAS, SOURCE, INFO, HAPLOTYPES, HAPLOTYPES_WITH_EVIDENCE, DEFINITE_GENOTYPES, NUM_HOM_ALLELE1, NUM_HOM_ALLELE2, NUM_HOM_ANY, NUM_HET, EXPECTED_HOM_ALLELE1, EXPECTED_HOM_ALLELE2, EXPECTED_HET, CHI_SQUARED_PVALUE, LOG10_CHI_SQUARED_PVALUE, CROSS_ENTROPY_LOD, HET_CHI_SQUARED_PVALUE, LOG10_HET_CHI_SQUARED_PVALUE, HET_CROSS_ENTROPY_LOD, HOM_CHI_SQUARED_PVALUE, LOG10_HOM_CHI_SQUARED_PVALUE, HOM_CROSS_ENTROPY_LOD, LOD_SELF_CHECK, DISCRIMINATORY_POWER | none | Needs bit-exact commons-math MersenneTwister and permutation, and a chi-squared p-value matching to 6 decimals | L |
| ClusterCrosscheckMetrics | CrosscheckMetric file, `LOD_THRESHOLD` | `ClusteredCrosscheckMetric`: CLUSTER, CLUSTER_SIZE, then the 19 CrosscheckMetric columns | none | Cluster ids and order from Java hash iteration. Strict `>` threshold | S |

### 2.3 Later tools from triage (harness and interval helpers)

- BedToIntervalList (core): parse BED, write interval_list with `@SQ` header from a dictionary. Needed by Hs, TargetedPcr, Wgs `INTERVALS`. Both sides must read identical interval_lists.
- IntervalListToBed, IntervalListTools (merge, pad, uniquify, subtract only as the metrics need them), VcfToIntervalList, NonNFastaSize: small, built as library functions, not as CLI clones.
- CreateSequenceDictionary (core for the harness): `picard CreateSequenceDictionary` or `samtools dict`. The library reads `.dict` or `.fai` for contig lengths.
- Harness prep only: FixMateInformation (MC tag, via `samtools fixmate -m`), SetNmMdAndUqTags (`samtools calmd`), BuildBamIndex, SortSam, AddOrReplaceReadGroups, CleanSam, DownsampleSam, FastqToSam. Never ported.

### 2.4 Complexity summary

| Size | Tools |
|---|---|
| S | QualityYield, QualityScoreDistribution, MeanQualityByCycle, BaseDistributionByCycle, CollectDuplicateMetrics, RawWgs, CompareMetrics, Accumulate* (2), SNVQ, NonZeroCoverage, CheckDuplicateMarking, ConvertSequencingArtifactToOxoG, Jumping, BamIndexStats, ReadGroupChecksum, CheckTerminatorBlock, ClusterCrosscheck |
| M | InsertSize (rest), EstimateLibraryComplexity, GcBias, CollectMultipleMetrics, CollectVariantCallingMetrics, TargetedPcr, OxoG, Rrbs, UmiPrevalence, CompareSAMs, FindMendelianViolations, ExtractFingerprint |
| L | AlignmentSummary, WgsMetrics, RnaSeq, GenotypeConcordance, SequencingArtifact, UmiAware, CheckFingerprint, CalculateFingerprintMetrics |
| XL | HsMetrics, MarkDuplicates, CrosscheckFingerprints, CollectSamErrorMetrics |

## 3. Shared machinery and what it unblocks

Status on `origin/main` 0.3.0: M1 (workspace), M2 (metrics file), M15 (`intervals.rs`) and M16 (`sensitivity.rs`) exist. M9 exists in part (`collect.rs`). Alignment summary, quality yield, Hs and duplication code exist, so M6, M7, M8, M10, M12 and M13 are likely present in some form. Each one is marked "verify against section 5" until the harness confirms it. M14 (pileup) is not built and Hs shipped without it.

| # | Machinery | Core scope | Unblocks |
|---|---|---|---|
| M1 | Cargo workspace with one `core` crate, `py` depends on it by path | Removes the two byte-identical copies of `histogram.rs` and `insert_size.rs`. Done on main | Every tool |
| M2 | Metrics file writer and reader (header lines, `## METRICS CLASS`, `## HISTOGRAM`, Java double format, `?` for NaN, empty cell for null, `-1` conventions) | Reader parses its own and Picard text. Done on main | Every tool, CompareMetrics, Accumulate* |
| M3 | `Histogram<K>` extensions: median, MAD, mode, mean, SD, sum, `trimByWidth`, `trimByTailLimit`, Integer, Double, Byte and String keys | Exists for i64 keys | InsertSize, AlignmentSummary, RnaSeq, Jumping, Wgs, Hs, Rrbs |
| M4 | BAM/SAM/CRAM reader with header RG map and auxiliary tags (MC, MQ, SA, XN, OQ, RG, NM, BQ, RX) | CRAM needs reference | Every BAM tool |
| M5 | Accumulation-level framework (ALL_READS, SAMPLE, LIBRARY, READ_GROUP, PU key, `unknown`, `All_Reads.` histogram prefixes, per-tool level label) | One generic collector wrapper. Main accepts ALL_READS only | InsertSize, AlignmentSummary, GcBias, Hs, TargetedPcr, RnaSeq, Rrbs, Artifacts, CollectMultipleMetrics |
| M6 | Pair orientation (`SamPairUtil.getPairOrientation`) with MC mate end, and `ChimeraUtil` | Shared with insert size | InsertSize, AlignmentSummary, Jumping |
| M7 | Reference FASTA reader with `.fai` and `.dict`, cached uppercase contig, N detection | noodles-fasta | AlignmentSummary, GcBias, Wgs, Hs, Artifacts, OxoG, Rrbs, SamError |
| M8 | Alignment-block walker (CIGAR to ref and read offsets, strand flip, clips) | | AlignmentSummary, GcBias, Rrbs, Artifacts, Wgs, Hs |
| M9 | Single-pass multi-collector trait (`accept_record`, `finish`) | Needed so one read serves many collectors. `collect.rs` on main covers part | CollectMultipleMetrics, speed claim |
| M10 | Adapter utility (16-mer matcher, 1 mismatch) | | AlignmentSummary, Wgs, Hs, TargetedPcr |
| M11 | Quality decode with OQ fallback | | QualityYield, QualityScoreDistribution, MeanQualityByCycle, Artifacts |
| M12 | Per-cycle accumulator | | MeanQualityByCycle, BaseDistributionByCycle, AlignmentSummary BAD_CYCLES |
| M13 | Lander-Waterman library size solver (`estimateLibrarySize`) | Exact port including `Math.exp` order | CollectDuplicateMetrics, MarkDuplicates, EstimateLibraryComplexity, Hs, Jumping |
| M14 | Locus pileup engine (sorted multi-contig, filters with per-reason base counters, overlap rule, caps) | The largest new piece. Not needed by Hs | Wgs, Raw, NonZero, TargetedPcr (if it uses the pileup), OxoG, SamError, fingerprints |
| M15 | Interval list reader, merger, overlap detector, SQ-header genome size | BedToIntervalList. Exists on main | Wgs `INTERVALS`, Hs, TargetedPcr, RnaSeq (rRNA), Artifacts |
| M16 | Theoretical sensitivity model | Floating-point order must match. Exists on main | Wgs family, Hs, TargetedPcr |
| M17 | refFlat reader and gene, transcript, exon model | | RnaSeq |
| M18 | Read-name tile/x/y parser and optical duplicate finder | | MarkDuplicates, EstimateLibraryComplexity |
| M19 | ReadEnds sort with spill, duplicate scoring, BAM writer with tag edits | | MarkDuplicates, UmiAware |
| M20 | VCF/BCF reader with index, `VariantContext` type classifier, dbSNP position bitsets | | CollectVariantCallingMetrics, GenotypeConcordance, Mendelian, SamError, Artifacts (dbSNP mask) |
| M21 | GenotypeConcordance scheme table and allele normalization | | GenotypeConcordance |
| M22 | Java compatibility: `String.hashCode`, `HashMap` order, `java.util.Random`, MersenneTwister | Only for tools whose output order or sampling depends on them | AccumulateVariantCallingMetrics, OxoG, Artifacts, Crosscheck family, SamError |
| M23 | Harness: oracle runner, comparator, report (section 5) | Extends `tests/concordance/` on main | Every tool's exit criterion |

Dependency order by tool:

| Tool | Needs beyond M1 to M4 |
|---|---|
| InsertSize | M3, M5, M6 |
| QualityYield, QualityScoreDistribution | M11 |
| MeanQualityByCycle, BaseDistributionByCycle | M11, M12 |
| CollectDuplicateMetrics | M5-lite (library only), M13 |
| AlignmentSummary | M5, M6, M7, M8, M10, M12 |
| GcBias | M5, M7, M8 |
| EstimateLibraryComplexity | M13, M18, sort |
| CollectMultipleMetrics | M9 and every collector it dispatches |
| WgsMetrics, Raw | M7, M10, M14, M15, M16 |
| HsMetrics | M5, M13, M15, M16 (and M10). Shipped on main without M14 |
| RnaSeq | M5, M6, M15, M17 |
| CollectVariantCallingMetrics | M20 |
| GenotypeConcordance | M20, M21, M15 |
| AccumulateVariantCallingMetrics | M20, M22 |
| MarkDuplicates | M13, M18, M19 |

Scheduling check: a CI script reads this table and section 4 and fails if any M-item a tool needs is not built or scheduled in a phase at or before the tool's phase. This check would have caught M18 (EstimateLibraryComplexity) and M22 (AccumulateVariantCallingMetrics), which the earlier draft left unscheduled.

## 4. Phased build order

Start state: `origin/main` at 0.3.0 (see section 0). Items main already holds are marked "verify against section 5": the harness runs them, and they count as done only when they reach claim level C2. Open work on main: accumulation levels (ALL_READS only today), PER_TARGET_COVERAGE, THEORETICAL_SENSITIVITY_OUTPUT, STOP_AFTER, and every tool not on main. Extend `tests/concordance/` (90 cases, Picard-syntax CLI) and `KNOWN_DIFFERENCES.md` instead of adding a parallel quirks file and a second harness. Release numbers continue from 0.3.0.

Release rule: a phase ends in a tagged release. A tool is listed "concordant" in that release only if it reaches claim level C2 (section 5.9). A tool at C1 ships as "experimental".

| Phase | Release | Work | Tools that must pass the exit criterion | Milestone check |
|---|---|---|---|---|
| P0 Verify and levels | 0.3.x | Verify against section 5: M1, M2, M9 in part, M15, M16, InsertSize, AlignmentSummary, QualityYield, CollectDuplicateMetrics, Hs (all 70 columns, HS_LIBRARY_SIZE, HET_SNP_*, PER_BASE_COVERAGE). Build M5 accumulation levels, InsertSize widths, histogram output, MC mate end. Extend the harness: oracle lock, cell-level output, report v1, `fields.yaml` completeness test, level-matrix probe (2 LB, 2 PU, one RG without PU), the 36-column AlignmentSummary check with an unmapped read and a single-read category, `STOP_AFTER` | InsertSize, AlignmentSummary, QualityYield, CollectDuplicateMetrics, Hs at ALL_READS | Report shows the side-by-side table for the tools on main. Q-001 confirmed or corrected |
| P1 Alignment summary levels | 0.4.0 | Levels for AlignmentSummary, bisulfite flag, adapters, chimera, reference input. Guide section 2.4 updated | CollectAlignmentSummaryMetrics at all levels | All 36 columns and the 4-column read-length histogram match on T0 to T3 |
| P2 Read-level tools | 0.5.0 | M11, M12. Build QualityScoreDistribution, MeanQualityByCycle, BaseDistributionByCycle. AccumulateQualityYield as a free add | The three new tools plus QualityYield and CollectDuplicateMetrics | Histogram-only files (no METRICS CLASS block) pass byte-for-byte after header strip |
| P3 Wrapper and complexity | 0.6.0 | Finish M9 trait, CollectMultipleMetrics with the default five plus QualityYield. GcBias (two files, `All Reads` label). M18 optical duplicate and tile parser. EstimateLibraryComplexity. Fuzzer and shrinker in harness. **Submit JOSS here (opinion)**: four tools plus the wrapper make a citable library, and the single-pass speed claim is measurable | CollectMultipleMetrics, GcBias, EstimateLibraryComplexity | File names and contents of every Multiple output match Picard. 10,000 fuzz examples clean |
| P4 Coverage | 0.7.0 | M14 pileup. BedToIntervalList. WgsMetrics, RawWgsMetrics, NonZeroCoverage (later, 60 lines, free). M15 and M16 already on main, verify | CollectWgsMetrics, CollectRawWgsMetrics | HET_SNP_* either exact or inside a registered bound. T4 scale run passes |
| P5 Capture | 0.8.0 | Hs levels, PER_TARGET_COVERAGE, THEORETICAL_SENSITIVITY_OUTPUT, then TargetedPcr on the same collector | CollectHsMetrics (TargetedPcr at C1 then C2) | Exome T2d slice matches. Dropout columns match with reference |
| P6 RNA | 0.9.0 | M17. CollectRnaSeqMetrics with ribosomal intervals and the 101-row histogram | CollectRnaSeqMetrics | nf-core-style stranded RNA slice matches |
| P7 Variants | 0.10.0 | M20, M21, M22. CollectVariantCallingMetrics, GenotypeConcordance, then AccumulateVariantCallingMetrics (needs M22; move it to the later queue instead if M22 is not ready) | CollectVariantCallingMetrics, GenotypeConcordance | Scheme table cells all hit by the corpus. Zero `UNTESTED` columns |
| P8 Release 1.0 | 1.0.0 | CompareMetrics CLI. Stability gate (two clean releases). Full validation pack. Register frozen | All 19 core tools at C2, a subset at C3 | Section 5.9 criteria |
| P9 Later queue | 1.x | In order: MarkDuplicates (metrics and flags, needs M13, M18, M19), SequencingArtifact + ConvertToOxoG + OxoG (need M22), Jumping, Rrbs, BamIndexStats, ReadGroupChecksum, CheckTerminatorBlock, FindMendelianViolations, SNVQ, UmiPrevalence, UmiAware (M19), fingerprint family (M22), SamErrorMetrics | Each at C2 before its release note | Each later tool re-reads the Picard source for its unread parts first (ArtifactCounter, BaseCalculator, ConvertSequencingArtifactToOxoG) |

Note on core count: the 19 core tools include CompareMetrics, BedToIntervalList and CreateSequenceDictionary. These are harness and library pieces. They reach the exit criterion by being used by the harness without error, not by a Picard-vs-ours cell comparison. CompareMetrics is validated by the planted-fault test (section 5.7).

Order rationale:
- P0 first because main already ships tools that have no claim level yet. Verification costs less than new code and fixes the starting line.
- P1 and P2 follow because levels (M5) touch every tool and the S tools prove the histogram-only file shape cheaply.
- P3 holds M18 because EstimateLibraryComplexity needs the optical duplicate finder.
- P4 (pileup) is independent of Hs on main, so it follows the verified tools and is not a prerequisite for P5.
- P7 is last of the core because the VCF layer shares nothing with the BAM path.
- MarkDuplicates leads the later queue because it is the most used tool, but it needs M19 and BAM writing, which doubles the surface to validate.

## 5. Validation model

One model, merged from two designs. The unit is one **cell**: (tool, input, accumulation level, row key, field) or (…, histogram, bin). Each cell holds the Picard 3.5.0 printed string, the bammetrics printed string, a rule and a status. The central artifact is the side-by-side report rendered from one `compare.jsonl`, produced by the existing `tests/concordance/` runner (extended, not replaced).

### 5.1 Principles

- Differential testing. Same input, same arguments, two implementations, one comparator.
- Picard 3.5.0 is the specification. Where Picard's behaviour is odd, bammetrics copies it and the quirks register records it.
- Compare printed text. Picard prints doubles with at most 6 decimals. The printed string is the contract.
- Oracle is checked too. On T0, expected values are also derived by hand or by an independent Python reference. The report shows hand, Picard and bammetrics for those cases.
- Fail on unknown. A new column, row, bin or file with no policy entry fails the gate.
- Test the tester. Planted faults must fail the harness.
- No vacuous pass. A case where both tools fail is not agreement, and a case with zero compared cells is a failure (5.5, 5.9).

### 5.2 Corpus tiers

Raw public data is never committed. The repo keeps manifests (URL, size, sha256, slicing command, accession, licence). No clinical or AHS data enters any tier. Real slices are re-headered: `@PG` stripped, `@RG` sample names replaced.

| Tier | Content | Source | Size | Runs |
|---|---|---|---|---|
| T0 | About 60 hand-built BAMs, one Picard rule each: even/odd counts, half-integer median, single pair, mode tie, 4.9% vs 5.1% orientation cutoff, MC present/absent, dovetail, MAPQ 19/20, base quality 19/20, IUPAC and `=` bases, QC-fail, secondary, supplementary, TLEN 0, mate-unmapped, one unmapped read, one single-read category (SD_READ_LENGTH `?`), empty and header-only BAM, and the level-matrix probe (2 LB, 2 PU, one RG without PU). One case per trap in section 2 | SAM text, `samtools view -b`, `sort`, `index`, `faidx`, `picard CreateSequenceDictionary`. No pysam | Under 20 KB each | Every PR |
| T0h | Hand-computed answers for 3 BAMs (5-pair fixture: median 310, MAD 10, SD 8.164966, W50 21) | By hand, checked once against Picard | 3 files | Every PR |
| T1 | Seeded generator grid: read length 36/75/100/150/301, insert mean and SD, 1/2/4 RGs with and without PU, error and indel rates, MAPQ, flag and orientation mix. 240 seeds, 24 in the PR subset | Generator regenerated from seed and generator hash | 1 to 5 MB each | PR subset, nightly all |
| T2a | GIAB HG002 2x250 and HG001 2x100 or 2x150: chr21 1 Mb slice, chr6 MHC 1 Mb slice, chrM | NIST/GIAB, ENA | About 25 MB each | Nightly, release |
| T2b | One 1000 Genomes high-coverage CRAM slice (CRAM path) | 1000G | About 25 MB | Nightly, release |
| T2c | ENCODE paired-end stranded total RNA-seq, STAR, GRCh38: N ops, MAPQ 0/1/3/255, no MC | ENCODE portal (check licence per accession) | About 40 MB | Release; needed from P6 |
| T2d | Public exome slice with targets (GIAB HG001 exome or 1000G exome with pulldown BED), plus GENCODE CDS exons as a licence-clean target set | ENA, GIAB, 1000G | About 50 MB | Needed from P5 |
| T3 | Adversarial: Picard-marked duplicates, secondary and supplementary, QC-fail, clips, I/D/N/P, `=`/`X`, lowercase and IUPAC, zero-read, header-only, single-read, missing MC, missing RG, unsorted. SAM, BAM and CRAM variants | Mutations of T1 and T2 | Under 10 MB each | Every PR (small set), nightly (full) |
| T4 | Differential fuzz BAMs from a `BamSpec` (section 5.4) | Hypothesis, seeded | 200 per PR, 10,000 nightly, each under 20 KB | PR, nightly |
| T4r | Shrunk regressions, one SAM text file per bug plus Picard output | Fuzzer output, committed | Under 5 KB each | Every PR |
| T5 | Scale: one full 30x chromosome (about 20M reads) | 1000G or GIAB | About 1.5 GB | Release |
| T6 | Variant inputs: truth and call VCF pairs (GIAB truth sets), multi-sample trio VCF with PED, dbSNP slice (bgzip + tabix) | GIAB, 1000G | Under 50 MB | Needed from P7 |
| T7 | Format edge cases: SAM and CRAM, no index, unsorted, no `@RG`, one LB over several RG, one SM over several LB | Generator | 5 MB | Every PR |
| T8 | Picard and htsjdk test BAMs at tag 3.5.0 | `git sparse-checkout` | Under 20 MB | Release |

Storage: T0 to T3, T4r and T7 live in the repo or as release assets keyed by sha256. Regenerated files are compared by sha256 and a change is flagged, because aligner output changes with its version. Slices and Picard outputs go to a Zenodo deposit so the paper carries a DOI.

Independence count for a 95% claim: zero failures in n independent draws bounds the pass rate at 1 - 0.05^(1/n). That is 95.05% at n=59 and 99.0% at n=299. T0 cases are hand-targeted and T1 seeds come from one generator, so the cases are not all independent draws. The bound is computed on T1 and T2 draws only, and the statistical wording is dropped for T0. Target at least 299 (tool, input) cases per tool by C2 (T0 about 60, T1 240, T2 about 12). Cases where both Picard and bammetrics fail are counted separately and never toward the 299. The bound applies to the generator distribution.

### 5.3 Oracle and pinning

- Jar: `picard.jar` 3.5.0. One design recorded sha256 `b7d97861c3a54ba5a421f5a317f38382f955803862d30ef4aca2bcdc54943631` (Homebrew `picard-tools` 3.5.0 on the dev Mac). Verify before writing `oracle.lock`. Also record the bundled htsjdk version and `java -version`.
- JDK: Temurin 17 with the exact build in `oracle.lock`. The dev Mac runs OpenJDK 27.0.1. A nightly job runs JDK 17 and 21 and diffs the outputs byte for byte. A difference is an oracle entry, because `Math.exp`-type calls may differ in the last ulp.
- Fixed environment: `LC_ALL=C`, `-Duser.language=en -Duser.country=US -Dfile.encoding=UTF-8`, `VALIDATION_STRINGENCY=SILENT`, fixed command-line template per tool, no-op `Rscript` shim on PATH. One probe under `de_DE` decides whether the flags matter.
- Entry point is always `picard` (plain `java` is broken here).
- `oracle.lock` lists, per (tool, input, variant): input sha256, exact command, output body sha256 (from `## METRICS CLASS` on), wall time.
- Qualification: run each Picard command 3 times. Different bytes below the header flag the input `ORACLE_FAULT` and exclude it until explained.
- Stored outputs `oracle/picard-3.5.0/<input_id>/<Tool>.<variant>.txt` are committed. PR CI reads them and needs no JVM.
- Picard failures are recorded: exit code and first exception line. bammetrics must also exit non-zero and write no file. Messages are not compared. Examples: PF_READS == 0, unindexed VCF, missing SNVQ tags. These cases check error parity only.
- Header lines (`# Started on:`, command line) are stripped before parsing.
- Upgrade policy: a new Picard version gets `oracle/picard-X.Y.Z/`. CI emits a Picard-vs-Picard delta. Claims are per pinned version.

### 5.4 Fuzz and shrinking

- Hypothesis builds a `BamSpec` (not bytes): 1 to 3 contigs of 1 to 5 kb with N runs, lowercase and IUPAC; up to 40 templates with orientation class, read length from {1, 36, 100, 150, 1000}, CIGAR over `M I D S H N = X P` under validity rules, MAPQ from {0, 1, 19, 20, 21, 60, 255}, base quality from {0, 2, 19, 20, 21, 40, 93}, valid flag combinations, tags `MC` (correct, missing, wrong), `MQ`, `SA`, `XN`, `OQ`, `RG`, `NM`, `MD`; file-level sort order, header variants, empty and single-end files. Integers are edge-biased at thresholds and 5% cutoffs.
- Per example: write SAM text, convert with samtools, run Picard (cached by sha256 of BAM bytes + args + jar sha + JVM major), run bammetrics, compare. At most 8 JVMs at once.
- Shrinking: Hypothesis shrinks the spec, then `ddmin` on SAM lines (header and order kept, mates removed together). Output `tests/regressions/<sha8>.sam`, `.args`, `.picard.txt`. The case list globs the folder.
- Picard-free properties (5,000 per PR): thread invariance (1 vs 8), surface invariance (CLI, Python, `collect`), format invariance (BAM, SAM, CRAM), duplication (BAM plus a renamed copy doubles counts, ratios unchanged), level additivity (SAMPLE counts equal the sum over LIBRARY rows), irrelevant records (unmapped mate-less reads leave insert size unchanged), order (permuting equal-coordinate records changes nothing in order-insensitive tools).
- Crash fuzz (nightly, no oracle): `cargo-fuzz` on the readers. Property: no panic, no hang over 10 s.

### 5.5 Comparison rules

Rows match by key, never by position. Key = metric-specific columns (CATEGORY, PAIR_ORIENTATION, CONTEXT, etc.) plus SAMPLE, LIBRARY, READ_GROUP. Level names are normalised per tool before keying (GcBias prints `All Reads` and READS_USED `ALL`). The field policy lives in `tests/concordance/fields.yaml`: one entry per (tool, column) with type, rule and quirk ids. A completeness test fails if any Picard column in any reference file lacks an entry. Option and level matrices come from `picard <Tool> -h`, extracted by script into `tests/concordance/tools/<Tool>.yaml`.

| Value type | Rule | Status on failure |
|---|---|---|
| Integer, count, min/max | Decimal string equality | DIFF |
| Double | bammetrics formats with the Java `FormatUtil` port (6 decimals, trailing zeros dropped). Printed strings must be equal. The report also stores abs diff and relative diff. ULP distance is computed only from unrounded doubles and is `n/a` for printed-string DIFFs | DIFF |
| NaN, infinity | Both `?` | DIFF |
| Null versus zero | `""` equals only `""`. `0` is not `""`. `-1` sentinels compared as text | DIFF |
| String, enum | Exact, case-sensitive | DIFF |
| Row set | Same key set. Row order is also checked (`ROW_ORDER`) except for tools whose Picard order comes from Java hash iteration (OxoG, Artifacts, Accumulate*, Crosscheck family). For those, order is compared as a quirk entry | ROW_ONLY_*, ROW_ORDER |
| Column set | Same names, same order, same `METRICS CLASS` name | COLUMN_* |
| Histogram | Header key class (`Integer`, `Double`, `Byte`, `String`), column names, key list in order, every cell as text. A bin absent on one side is a failure. Picard fills 0 itself where the union of keys across levels applies, so the zero cells are compared as printed. The report also shows the bin sum against READ_PAIRS or the equivalent total | BIN_ONLY_*, DIFF |
| Extra files | Same rules per file. PER_BASE_COVERAGE byte-for-byte after gunzip | FILE_* |
| BAM output | `samtools view` text per record, `@PG` excluded | RECORD_DIFF |
| Errors | Exit zero or non-zero, file present or absent. Cases where both sides fail are tallied as `BOTH-ERR`, never as agreement toward a gate count. Per tool and per N/A option, at least k cases (k set in `fields.yaml`, default 3) must have Picard exit 0 with a non-empty body | ERROR_MISMATCH, VACUOUS |
| Cells compared | A case fails when the compared cell count is 0 or differs from the Picard cell count in the reference files | VACUOUS |
| Charts, header lines | Not compared | none |

Decisions on points where the two designs differed:
- No tolerance at the gate. Design 1's MATCH-ROUNDED is dropped because both sides print through the same formatter, so rounding differences are real differences. A tolerance exists only as a quirks-register entry with a bound (expected users: HET_SNP_SENSITIVITY, HET_SNP_Q, HS_LIBRARY_SIZE if the JDK `exp` differs).
- Absent bin is not zero (Design 2). It catches histogram trimming bugs.
- Two reports are added to the comparator: a rounding-boundary list (unrounded Python values within 1e-9 of a 6th-decimal tie) and a degenerate-field list (a column that is all `0` or `""` across the corpus is `UNTESTED`, which is how `PCT_ADAPTER = 0` on a synthetic tier would surface), with the exempt list described in 5.9.
- Second comparator: `picard CompareMetrics` runs on every case. Disagreement on pass or fail is a harness bug.
- Each `compare.jsonl` row records the bammetrics binary version and the input sha256.

### 5.6 Triage

Every non-MATCH cell gets one disposition in `tests/concordance/triage/<TOOL>-<n>.yaml` (cell filter, class, evidence, owner, status). Decision order:

1. ORACLE_FAULT. Picard fails, differs across 3 runs or depends on environment. Fix the oracle job. Exclude the input. Log O-nnn.
2. CORPUS_FAULT. Input invalid for the tool. Fix the generator.
3. RULE. The comparator rule is wrong. Fix the rule, rerun T0 probes.
4. QUIRK. The Picard 3.5.0 or htsjdk source (tag cited) shows the odd output is intended. Reproduce it, add Q-nnn, add a case that exercises it. Status MATCH-Q, a pass.
5. DEVIATION. bammetrics differs on purpose. Allowed only as an explicit error or behind `--picard-compat=false`. Record D-nnn with rationale, quantified impact, approver, review date, bound. Headline fields cannot carry a DEVIATION. The claim covers compat mode.
6. BUG. Everything else. Add the failing case to T0 first (or commit the shrunk regression), then fix. Open BUG blocks release.

| Status | Meaning | Gate |
|---|---|---|
| MATCH | Printed strings equal | pass |
| MATCH-Q | Equal by reproducing a registered quirk | pass |
| DEVIATION | Registered, inside bound | pass inside bound |
| DIFF | Unequal, open | fail |
| MISSING-OURS, EXTRA-OURS | Column or row only on one side | fail |
| ORACLE-ERR | Picard aborted or drifted | input excluded, listed |
| BOTH-ERR | Picard and bammetrics both failed | counted separately, not toward any case total |
| VACUOUS | Zero cells compared, or count differs from Picard | fail |
| N/A | Option or level rejected with exit 1 | tool not concordant for that option |
| SKIP | Not run on this platform | listed |

### 5.7 Quirks register

The register extends `KNOWN_DIFFERENCES.md` on main (structured entries rendered into it and into the report), not a parallel `quirks.yaml`. Entry: id, tool, fields, rule, Picard source reference (tag 3.5.0), case ids, bound (deviation only). A quirk with no case that exercises it fails. A deviation without an entry fails.

Seed entries from this plan:
- Q-001 READ_GROUP level value is the PU field, `unknown` when absent, with an extra all-zero row for some tools (all level tools). Unconfirmed: all 7 oracle runs used ALL_READS only, and main's `KNOWN_DIFFERENCES.md` says bammetrics produces no LIBRARY or READ_GROUP rows. Do not mark confirmed until the level-matrix probe (2 LB, 2 PU, one RG without PU) is in the T0 oracle run.
- Q-002 InsertSize mean and SD come from the trimmed histogram.
- Q-003 Mode ties go to the lowest key.
- Q-004 InsertSize min and max are histogram keys.
- Q-005 MC-missing orientation fallback.
- Q-006 AlignmentSummary adapter check applies only to unmapped and MAPQ 0 reads.
- Q-007 AlignmentSummary MEAN_ALIGNED_READ_LENGTH counts unmapped reads as 0.
- Q-008 Chimera denominator includes pairs with absent MQ.
- Q-009 QualityYield READ_LENGTH is integer division.
- Q-010 MeanQualityByCycle second-read offset uses the last positive-sum cycle.
- Q-011 WGS PCT_EXC_* denominator.
- Q-012 Hs ZERO_CVG_TARGETS_PCT uses unmerged targets.
- Q-013 UMI prevalence barcode filter inversion.
- Q-014 `?` for NaN and `-1` for non-computable HS_PENALTY.
- Q-015 HET_SNP floating-point tolerance bound (candidate DEVIATION).
- Q-016 Java hash-order row output in OxoG, Artifacts, Accumulate*, Crosscheck.
- Q-017 AccumulateVariantCallingMetrics lossy inversion.
- Q-018 CompareSAMs strand bug.
- Q-019 InsertSize skips records with TLEN 0, secondary, supplementary or unmapped mate.
- Q-020 GcBias prints the level as `All Reads` and READS_USED as `ALL`.
- Q-021 AlignmentSummary SD_READ_LENGTH prints `?` for a single-read row, and unmapped reads fall in the UNPAIRED_ALIGNED_LENGTH_COUNT bin at key 0.
- Entries Q-002 to Q-005 come from the implementation guide and need confirmation against the 3.5.0 source on T0. Q-006 to Q-021 come from the inventory runs and the correction runs.

### 5.8 CI, performance and tooling

Harness contract (extends the existing `tests/concordance/` runner and generator, does not replace them):
- The runner returns cell records and asserts nothing inside.
- Rules load from `fields.yaml`. A pytest option `--concordance-out DIR` writes `compare.jsonl`.
- `tools/make_oracle.py` regenerates oracle outputs. `tools/report.py` renders the report. A `make repro CASE=<id>` prints both commands and both outputs.
- Planted faults in `tests/concordance/faults/` (applied by patch): population SD, mean on the untrimmed histogram, mode ties to the highest key, MAD without the 2x key, count first-of-pair, ignore MC, count second-of-pair records with TLEN 0 or secondary or supplementary flags. Each later tool adds its own. The harness must fail on every fault.
- The workspace (M1) is on main, so the earlier CI step that diffs two crate copies of `histogram.rs` and `insert_size.rs` applies only to the old branch and is dropped.

| Job | Trigger | Content | Gate |
|---|---|---|---|
| unit | Every PR | fmt, clippy, `cargo test`, jar sha256 check | blocks |
| concordance-pr | Every PR | T0, T0h, T4r, T7, 24 T1 seeds, T3 small set, 200 derandomised fuzz examples, 5,000 Picard-free properties, threads 1 and 8, ubuntu and macos-14 against committed oracle | blocks on any DIFF, VACUOUS or missing waiver. Target under 10 min |
| oracle-refresh | Weekly, on `oracle.lock` change | Temurin 17 and 21, Picard by sha256, regenerate, diff | fails on drift |
| concordance-full | Nightly (seed logged) | All T1 seeds, T2, T3, 10,000 fuzz, `cargo-fuzz` 30 min. Failures shrunk and filed as issues | blocks release |
| mutation | Weekly | `cargo-mutants` on metric modules, `cargo llvm-cov` | gate 5 |
| scale and perf | Release | T5 and perf on a dedicated runner (Mac Studio through `agent-job --queue`, public data only) | blocks release |
| negative-control | Release | Planted faults | blocks release |

Every job uploads `compare.jsonl`. A PR comment shows the matrix. The README badge reads `summary.json`.

Performance side by side (`benches/bench.py` extended):
- 1 warm-up and 5 timed runs on a warm cache. Report min, median, max wall time.
- Peak RSS from `/usr/bin/time -l` (macOS, bytes) or `-v` (Linux, KB), normalised to MB.
- Thread peak sampled every 50 ms (`/bin/ps -M` or `/proc/<pid>/status`).
- bammetrics at 1 and 8 threads. Picard with default JVM flags and with `-Xmx4g`.
- Inputs: T1, a 2M-pair BAM, a 20M-pair BAM (nightly). Each row carries host, CPU, cores, OS, commit sha, Picard sha, bammetrics version and "outputs identical".
- Gate: dedicated runner only; fail on more than 25% median wall regression against the last release. Shared runners report but never gate.
- `perf.json` is the only source for a speed claim. Known measurements: Picard 0.74 s and 138 MB peak RSS on a 5-pair BAM (JVM floor); bammetrics 0.00 s and 3.6 MB.

### 5.9 Exit criterion and claim levels

A tool is **concordant with Picard 3.5.0** at level C2 when all hold on one tagged commit:

1. Every cell on T0, T0h, T1, T2, T3, T7 and T8 is MATCH or MATCH-Q, at every accumulation level the tool supports. Zero DIFF, zero VACUOUS, zero unregistered DEVIATION, zero open BUG, zero DEVIATION on headline fields.
2. Every option that can change a number (from `tools/<Tool>.yaml`) has at least one non-default case, or exits 1 and is listed N/A. For the tool and for each N/A option, at least k cases have Picard exit 0 with a non-empty body, so a tool that errors on every input cannot pass.
3. Every concordance trap in section 2 has a T0 case that fails under its planted fault.
4. At least 299 (tool, input) cases pass, T2 included. BOTH-ERR cases do not count. The 95% bound is stated on T1 and T2 draws only (5.2).
5. Every output column has at least 3 distinct non-default values in the corpus. Columns that fail are `UNTESTED` and block the claim, except columns on the explicit exempt list in `fields.yaml`: key columns (SAMPLE, LIBRARY and READ_GROUP at ALL_READS are blank by design) and structural constants or enums (PAIR_ORIENTATION, BAIT_SET, CATEGORY, METRICS CLASS name). Each exempt entry carries a reason. An exemption without a reason fails the completeness test.
6. Thread, surface and format invariance pass on T1 to T3.
7. Differential fuzz: 10,000 examples over 3 nightly seeds with zero failures.
8. Every case compared at least one cell, and the compared count equals Picard's cell count (the VACUOUS rule).

Level C1 is gates 1 to 3 on T0 and T1 only (label "experimental"). Level C3 adds: T5 scale run passes on a recorded date with recorded hashes; precision (10 repeat runs byte-identical; identical across Linux x86_64 and macOS arm64, 1 vs N threads, input formats, record order); `cargo-mutants` kill rate at least 90% and line coverage at least 95% on the tool's modules; negative control passes; two consecutive releases pass with an unchanged `oracle.lock`; the validation report is signed with name and date. C3 is the level offered to a lab as the start of its own verification.

Claim wording: "N of M printed cells identical to Picard 3.5.0 on corpus sha X". Never "bit-for-bit" or "100% accurate".

Reportable range, first declaration: short reads of length 36 to 301, BAM, SAM and CRAM with reference, coordinate, queryname and unsorted order, insert size up to 100 kb. Out-of-range inputs (long reads, colour space) produce an explicit error, never a number. Long-read slices (GIAB ONT or HiFi) are tested as error or match cases.

Lab terms map as follows: reference method is the pinned jar; accuracy is cell concordance; repeatability is the 10-run test; intermediate precision is the cross-platform test; analytical specificity is the filter cases (secondary, supplementary, duplicate, QC-fail, unmapped); interference is the missing MC, RG or PU and odd-CIGAR cases. Cite ISO 15189:2022 clause 7.3, CLIA 42 CFR 493.1253 and the AMP/CAP NGS bioinformatics guidance from the current texts when the pack is written. They are not quoted here.

### 5.10 Side-by-side report

Canonical output `compare.jsonl`, one record per cell:

```json
{"tool":"CollectInsertSizeMetrics","case":"t0_fixture","level":"ALL_READS",
 "row":{"PAIR_ORIENTATION":"FR","SAMPLE":"","LIBRARY":"","READ_GROUP":""},
 "section":"metric","column":"STANDARD_DEVIATION","picard":"8.164966",
 "bammetrics":"8.164966","status":"MATCH","quirk":null,"triage":null,"ulp":null,
 "bammetrics_version":"0.3.0","input_sha256":"<sha>"}
```

`make report` renders `report.html` (single file, data inlined), `report.md`, `summary.json`, `quirks.md`, `perf.md`. Run metadata at the top: git sha, Picard jar sha256, JDK, OS, corpus manifest sha256, seed.

Layout:
- Matrix: tool by tier, columns cells, MATCH, MATCH-Q, DEVIATION, DIFF, MISSING, N/A, BOTH-ERR, VACUOUS.
- Field scoreboard: tool by field, cells, matches, max abs diff, worst input.
- Per case: one table per row (Column, Picard, bammetrics, Status, Note), then histograms. Matching bins collapse by default with a toggle. For a numeric DIFF the first differing character is highlighted.
- Triage index with links to open items.
- Markdown version: failing and quirk rows in full, MATCH rows collapsed to a count, capped at 50 rows with a link to the HTML, so it fits a PR comment (65,536 characters).

Mock, using the real 5-pair fixture `tests/data/insert_size.bam` (insert sizes 300, 310, 310, 320, 500). On a rebuilt fixture run, all 19 value columns and the 3-bin histogram (300:1, 310:2, 320:1, bin 500 trimmed) matched Picard 3.5.0, so the passing rows below are real. The failure rows are invented and marked.

```
CollectInsertSizeMetrics | t0_fixture | ALL_READS | PAIR_ORIENTATION=FR SAMPLE="" LIBRARY="" READ_GROUP=""
Picard 3.5.0 (sha b7d97861) vs bammetrics (commit abc1234)

Column                      Picard      bammetrics   Status
MEDIAN_INSERT_SIZE          310         310          MATCH
MODE_INSERT_SIZE            310         310          MATCH
MEDIAN_ABSOLUTE_DEVIATION   10          10           MATCH
MIN_INSERT_SIZE             300         300          MATCH
MAX_INSERT_SIZE             500         500          MATCH
MEAN_INSERT_SIZE            310         310          MATCH-Q  Q-002 trimmed histogram
STANDARD_DEVIATION          8.164966    8.164966     MATCH
READ_PAIRS                  5           5            MATCH
WIDTH_OF_10..40_PERCENT     1           1            MATCH x4
WIDTH_OF_50..80_PERCENT     21          21           MATCH x4
WIDTH_OF_90_PERCENT         381         381          MATCH
WIDTH_OF_95_PERCENT         381         381          MATCH
WIDTH_OF_99_PERCENT         381         381          MATCH
Row total: 19 of 19 value columns pass (4 key columns: PAIR_ORIENTATION, SAMPLE, LIBRARY, READ_GROUP)

HISTOGRAM java.lang.Integer  (insert_size | All_Reads.fr_count)
Bin    Picard   bammetrics   Status
300    1        1            MATCH
310    2        2            MATCH
320    1        1            MATCH
(500 absent on both sides: histogram trimmed)    3 of 3 bins MATCH

Invented failure states, to show the statuses:
STANDARD_DEVIATION   8.164966   8.164967   DIFF         ulp n/a (printed strings), first diff at char 8
HS_LIBRARY_SIZE      3207741    3207742    DEVIATION    D-004 bound +-1, JDK 21 vs 17
WIDTH_OF_99_PERCENT  381        (missing)  MISSING-OURS
READ_GROUP row rg2   -          present    EXTRA-OURS
```

Baseline for the CLI on the old branch (from the 400-pair probe, 2 libraries, 2 PUs): 95 cells, 8 MATCH, 1 DIFF (SD printed with 4 decimals), 86 MISSING-OURS. Concordance 8.4%. That number is stale for `origin/main`. Re-measure on main at P0 and use that as the starting line.

### 5.11 Feeding the JOSS concordance report

- Paper Table 1: tool by tier from `summary.json` (cases, cells compared, matched, quirk, deviation, N/A).
- Quirks table straight from the register with Picard source locations.
- Fuzz summary: examples run, bugs the fuzzer found, shrunk size distribution.
- Performance table from `perf.json`.
- Appendix: Picard-vs-Picard delta for each oracle bump.
- Archive per release on Zenodo: `report.html`, `compare.jsonl`, corpus (T0 to T4r, T7, hashes for the rest), `oracle.lock`, and the `make paper-repro` command. The paper cites the DOI. Every number in the paper or README traces to a `compare.jsonl` line or a `perf.json` row.
- README block between `<!-- concordance:start -->` and `<!-- concordance:end -->`, rewritten at release.
- A reviewer can run T0 and T1 without a JVM, against the committed oracle outputs.
- Validation pack under `docs/validation/`, generated where possible: plan, per-tool specification, test protocol, traceability matrix (Picard field to rule to test ids to status), per-release signed report, quirk and deviation registers, change control (a Picard, noodles or rustc bump or a metric-module edit triggers T0 to T2 in full), SBOM from `Cargo.lock` and `pyproject.toml`, reproduction instructions.

## 6. Risks and open questions

Risks:

| # | Risk | Effect | Mitigation |
|---|---|---|---|
| R1 | `HET_SNP_SENSITIVITY`, `HS_LIBRARY_SIZE` and any `Math.exp` or `log10` use differ in the last digit across JDKs or in Rust | Printed 6th decimal flips | Port arithmetic order exactly. JDK 17 vs 21 nightly. Register a bound only after observing a drift |
| R2 | Java hash and set order drives row order and cluster ids (OxoG, Artifacts, Accumulate*, Crosscheck, Cluster) | Row order not reproducible | Compare by key. Emulate order only if users diff files byte for byte (M22) |
| R3 | Picard source unread for ArtifactCounter, BaseCalculator, ConvertSequencingArtifactToOxoG, GcBias secondary handling | Wrong port | Read the 3.5.0 source before each later tool starts |
| R4 | Picard needs Rscript for some tools even when no chart is asked | No oracle output | Shim on PATH. Never compare charts. State that PDF output is not implemented |
| R5 | Resolved: `origin/main` (0.3.0) is ahead of the local branch by 23 commits | Duplicate work if the branch is used as the base | Rebase the plan and the work on `origin/main` (done in section 4) |
| R6 | noodles rejects or repairs records that htsjdk accepts | Different counts on malformed input | T3 cases. Explicit error for out-of-range input |
| R7 | Synthetic corpus gives zero values (PCT_ADAPTER = 0) | False sense of coverage | Degenerate-field report and `UNTESTED` gate, with a reasoned exempt list |
| R8 | License of ENCODE and vendor BEDs | Cannot redistribute | Ship manifests only. Use GENCODE CDS targets |
| R9 | Memory and time of T5 and 100k-line Raw histograms | CI cost | Run on the Mac Studio queue, not shared runners |
| R10 | Picard behaviours that look like bugs (inverted UMI filter, strand compare in CompareSAMs, NPE on missing tags) | Copying them is wrong for users | Compat mode default on. Document each in the register. Missing tags: explicit error is a DEVIATION (D-nnn) |
| R11 | Core set is 19 tools but three are harness pieces | Count looks larger than the user-facing work (16 tools) | State both numbers in release notes |
| R12 | Done: the 121-tool list was diffed against `picard` with no arguments, no missing, extra or duplicate names | none | none |
| R13 | Unverified columns: CollectQualityYieldMetricsFlow and CollectArraysVariantCallingMetrics (both out of scope) | None for the plan | Verify only if the decision changes |
| R14 | Gate can pass with no real comparison (both tools fail, zero cells compared, non-independent case counts, unreachable column rule) | False concordance claim | BOTH-ERR and VACUOUS statuses, minimum Picard-success cases per tool and option, bound on T1 and T2 only, reasoned exempt-column list (5.5, 5.9) |
| R15 | A machinery item is needed by a tool but built in no phase (M18 and M22 were missed in the first draft) | Tool blocked at its phase | Scheduling check in section 3 |

Open questions for the maintainer: see the list returned with this document.

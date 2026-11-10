"""Type stubs for the bammetrics Python module (Picard-concordant BAM QC metrics)."""
from os import PathLike
from typing import Any, Generic, Literal, Sequence, TypeVar, TypedDict, overload

_Path = str | PathLike[str]
_Paths = _Path | Sequence[_Path]

_MetricsClass = TypedDict("_MetricsClass", {"class": str, "columns": list[str], "rows": list[dict[str, Any]]})
_MetricsHistogram = TypedDict("_MetricsHistogram", {"key_class": str, "columns": list[str], "rows": list[dict[str, Any]]})
_MetricsDoc = TypedDict("_MetricsDoc", {"headers": list[str], "classes": list[_MetricsClass], "histogram": _MetricsHistogram | None})

class _VerdictPreset(TypedDict):
    name: str
    version: int
    hash: str
    source: str
    note: str

class _VerdictCounts(TypedDict):
    pass_: int
    warn: int
    fail: int
    not_measured: int

class _VerdictFile(TypedDict):
    path: str
    classes: list[str]

class _VerdictCheck(TypedDict):
    rule: str
    label: str
    metric: str
    level: str
    select: str | None
    unit: str | None
    value: float | None
    warn: str | None
    fail: str | None
    status: str
    reason: str | None
    plot: str | None

class _VerdictDict(TypedDict):
    """verdict.json (schema bammetrics.verdict/1) as dicts; see docs/verdict.md.

    status is "pass", "warn", "fail" or "incomplete"; a check's status is
    "pass", "warn", "fail" or "not_measured". counts has the key "pass"
    (spelled pass_ above because pass is a Python keyword)."""
    schema: str
    tool_version: str
    sample: str
    preset: _VerdictPreset
    status: str
    counts: _VerdictCounts
    files: list[_VerdictFile]
    checks: list[_VerdictCheck]

class _BamIndexRow(TypedDict):
    """One row of bam_index_stats(): a reference of the BAM header with its length and
    the aligned and unaligned record counts of the index, then a last row whose
    reference is "*" and whose unaligned count is Picard's NoCoordinateCount."""
    reference: str
    length: int
    aligned: int
    unaligned: int

class BammetricsError(Exception):
    """Base of UsageError, InputError and VerdictError."""

class UsageError(BammetricsError, ValueError):
    """The arguments cannot be acted on and nothing ran (the command's exit 2)."""

class InputError(BammetricsError, OSError):
    """An input could not be read, or a collector failed (the command's exit 1).
    After a failed collector, .run holds the run with the other files; it is
    None otherwise."""
    run: Run[Any] | None

class VerdictError(BammetricsError):
    """A verdict reached fail_on (the command's exit 10)."""
    verdicts: list[Verdict]
    run: Run[Verdict] | None
    code: int
    def __init__(self, message: str, verdicts: list[Verdict], run: Run[Verdict] | None = None) -> None: ...

class Tool(str):
    """One row of `bammetrics tools`; the string is the Picard class name."""
    picard_class: str
    verb: str
    summary: str
    def __new__(cls, picard_class: str, verb: str, summary: str) -> Tool: ...

class PresetInfo(str):
    """One row of `bammetrics presets`; the string is the preset name."""
    name: str
    version: int
    hash: str
    rules: int
    note: str
    def __new__(cls, name: str, version: int, hash: str, rules: int, note: str) -> PresetInfo: ...

class Check:
    """One check of a verdict. status is "pass", "warn", "fail" or "not_measured";
    text is the cell's text for a rule with an "is TEXT" bound (value is then None);
    limit is the fail bound of a fail, else the warn bound (or the fail bound
    when the rule has none); why is the reason in plain words, None for a pass."""
    rule: str
    status: str
    value: float | None
    text: str | None
    limit: str | None
    why: str | None
    label: str
    metric: str
    level: str
    unit: str | None
    warn: str | None
    fail: str | None

class Verdict:
    """A sample's verdict, as PREFIX.verdict.json holds it."""
    sample: str
    status: Literal["pass", "warn", "fail", "incomplete"]
    checks: list[Check]
    counts: dict[str, int]
    preset: dict[str, Any]
    def to_json(self) -> str: ...
    def to_dict(self) -> _VerdictDict: ...
    def exit_code(self, fail_on: Literal["warn", "fail"] = "fail") -> int:
        """The exit code bammetrics check gives this verdict: 0, 10 or 11."""

_V = TypeVar("_V", Verdict, None)

class Run(Generic[_V]):
    """A qc run (Run[Verdict]), or a collect(metrics=[...]) run (Run[None]: verdict
    and manifest None).

    files lists every file written or kept, the manifest last. metrics maps the
    verb of each tool (insert-size, hs-metrics.BAIT_SET for one of several hs
    runs, gc-bias.summary for a further file) to the parsed file, as
    read_metrics returns it. In a notebook the run shows the report body."""
    sample: str
    prefix: str
    files: list[str]
    manifest: str | None
    exit_code: int
    errors: dict[str, str]
    metrics: dict[str, _MetricsDoc]
    verdict: _V
    def _repr_html_(self) -> str | None: ...

class InsertSizeMetrics:
    """One CollectInsertSizeMetrics row; fields are Picard's columns in lower case."""
    pair_orientation: str
    read_pairs: int
    min_insert_size: int
    max_insert_size: int
    median_insert_size: float
    mode_insert_size: int
    median_absolute_deviation: float
    mean_insert_size: float
    standard_deviation: float
    width_of_10_percent: int
    width_of_20_percent: int
    width_of_30_percent: int
    width_of_40_percent: int
    width_of_50_percent: int
    width_of_60_percent: int
    width_of_70_percent: int
    width_of_80_percent: int
    width_of_90_percent: int
    width_of_95_percent: int
    width_of_99_percent: int
    sample: str | None
    library: str | None
    read_group: str | None
    def to_dict(self) -> dict[str, Any]: ...
    def __repr__(self) -> str: ...

class AlignmentSummaryMetrics:
    """One CollectAlignmentSummaryMetrics row (a CATEGORY); fields are Picard's columns in lower case."""
    category: str
    total_reads: int
    pf_reads: int
    pct_pf_reads: float
    pf_noise_reads: int
    pf_reads_aligned: int
    pct_pf_reads_aligned: float
    pf_aligned_bases: int
    pf_hq_aligned_reads: int
    pf_hq_aligned_bases: int
    pf_hq_aligned_q20_bases: int
    pf_hq_median_mismatches: float
    pf_mismatch_rate: float
    pf_hq_error_rate: float
    pf_indel_rate: float
    mean_read_length: float
    sd_read_length: float
    median_read_length: float
    mad_read_length: float
    min_read_length: float
    max_read_length: float
    mean_aligned_read_length: float
    reads_aligned_in_pairs: int
    pct_reads_aligned_in_pairs: float
    pf_reads_improper_pairs: int
    pct_pf_reads_improper_pairs: float
    bad_cycles: int
    strand_balance: float
    pct_chimeras: float
    pct_adapter: float
    pct_softclip: float
    pct_hardclip: float
    avg_pos_3prime_softclip_length: float
    def to_dict(self) -> dict[str, Any]: ...
    def __repr__(self) -> str: ...

class QualityYieldMetrics:
    """The CollectQualityYieldMetrics row; fields are Picard's columns in lower case."""
    total_reads: int
    pf_reads: int
    read_length: int
    total_bases: int
    pf_bases: int
    q20_bases: int
    pf_q20_bases: int
    q30_bases: int
    pf_q30_bases: int
    q20_equivalent_yield: int
    pf_q20_equivalent_yield: int
    def to_dict(self) -> dict[str, Any]: ...
    def __repr__(self) -> str: ...

class HsMetrics:
    """The CollectHsMetrics row; fields are Picard's columns in lower case. hs_library_size and the read-group columns are None when Picard leaves them empty."""
    bait_set: str
    bait_territory: int
    bait_design_efficiency: float
    on_bait_bases: int
    near_bait_bases: int
    off_bait_bases: int
    pct_selected_bases: float
    pct_off_bait: float
    on_bait_vs_selected: float
    mean_bait_coverage: float
    pct_usable_bases_on_bait: float
    pct_usable_bases_on_target: float
    fold_enrichment: float
    hs_library_size: int | None
    hs_penalty_10x: float
    hs_penalty_20x: float
    hs_penalty_30x: float
    hs_penalty_40x: float
    hs_penalty_50x: float
    hs_penalty_100x: float
    target_territory: int
    genome_size: int
    total_reads: int
    pf_reads: int
    pf_bases: int
    pf_unique_reads: int
    pf_uq_reads_aligned: int
    pf_bases_aligned: int
    pf_uq_bases_aligned: int
    on_target_bases: int
    pct_pf_reads: float
    pct_pf_uq_reads: float
    pct_pf_uq_reads_aligned: float
    mean_target_coverage: float
    median_target_coverage: float
    max_target_coverage: int
    min_target_coverage: int
    zero_cvg_targets_pct: float
    pct_exc_dupe: float
    pct_exc_adapter: float
    pct_exc_mapq: float
    pct_exc_baseq: float
    pct_exc_overlap: float
    pct_exc_off_target: float
    fold_80_base_penalty: float
    pct_target_bases_1x: float
    pct_target_bases_2x: float
    pct_target_bases_10x: float
    pct_target_bases_20x: float
    pct_target_bases_30x: float
    pct_target_bases_40x: float
    pct_target_bases_50x: float
    pct_target_bases_100x: float
    pct_target_bases_250x: float
    pct_target_bases_500x: float
    pct_target_bases_1000x: float
    pct_target_bases_2500x: float
    pct_target_bases_5000x: float
    pct_target_bases_10000x: float
    pct_target_bases_25000x: float
    pct_target_bases_50000x: float
    pct_target_bases_100000x: float
    at_dropout: float
    gc_dropout: float
    het_snp_sensitivity: float
    het_snp_q: float
    sample: str | None
    library: str | None
    read_group: str | None
    high_quality_coverage: dict[int, float]
    unfiltered_base_quality: dict[int, float]
    def to_dict(self) -> dict[str, Any]: ...
    def __repr__(self) -> str: ...

def insert_size(
    path: _Path,
    minimum_pct: float = 0.05,
    deviations: float = 10.0,
    histogram_width: int | None = None,
    min_histogram_width: int | None = None,
    include_duplicates: bool = False,
    metric_accumulation_level: Sequence[str] | None = None,
    reference: _Path | None = None,
) -> list[InsertSizeMetrics]:
    """CollectInsertSizeMetrics: one row per pair orientation and accumulation unit.

    metric_accumulation_level adds ALL_READS, SAMPLE, LIBRARY or READ_GROUP to the
    default ALL_READS; "null" clears the list, as in Picard."""

def tools() -> list[Tool]:
    """The tool verbs as `bammetrics tools` lists them: each a str equal to its
    Picard class, with verb and summary."""

def read_metrics(path: _Path) -> _MetricsDoc:
    """Read a Picard metrics file written by Picard or bammetrics.

    headers are the StringHeader texts; classes has one entry per
    ## METRICS CLASS block with its rows as dicts keyed by column; histogram
    is the ## HISTOGRAM section (key_class, columns, rows) or None. Cells are
    typed from their text: int, float, str or None for an empty field; "?"
    is nan."""

def presets() -> list[PresetInfo]:
    """The shipped verdict presets as `bammetrics presets` lists them (generic,
    clinical-panel, wgs, rna-seq): each a str equal to its name, with version,
    hash, rules (the rule count) and note."""

def evaluate(preset: str | PathLike[str], paths: Sequence[_Path], sample: str | None = None) -> _VerdictDict:
    """Evaluate a preset over one sample's metrics files and return the verdict.

    preset is a shipped name from presets() or the path of a preset file
    (docs/verdict.md). paths are Picard metrics files written by bammetrics
    or Picard. A bad preset raises ValueError; an unreadable file OSError."""

# BEGIN registry tools: generated from crates/core/src/registry.rs (UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs)
# Each returns its rows as dicts keyed by Picard's column names; output= also writes the files.

def collect_insert_size_metrics(
    input: _Path,
    *,
    deviations: float | None = None,
    histogram_width: int | None = None,
    minimum_histogram_width: int | None = None,
    minimum_pct: float | None = None,
    include_duplicates: bool | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectInsertSizeMetrics: one row per pair orientation and accumulation unit, keyed
    by Picard's column names; the metrics file adds the insert size histogram."""

def collect_alignment_summary_metrics(
    input: _Path,
    *,
    reference_sequence: _Path | None = None,
    max_insert_size: int | None = None,
    expected_pair_orientations: str | Sequence[str] | None = None,
    adapter_sequence: str | Sequence[str] | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    is_bisulfite_sequenced: bool | None = None,
    collect_alignment_information: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
) -> list[dict[str, Any]]:
    """CollectAlignmentSummaryMetrics: one row per CATEGORY of each accumulation unit,
    keyed by Picard's column names; the metrics file adds the read length histograms."""

def collect_quality_yield_metrics(
    input: _Path,
    *,
    use_original_qualities: bool | None = None,
    include_secondary_alignments: bool | None = None,
    include_supplemental_alignments: bool | None = None,
    flow_mode: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectQualityYieldMetrics: one row of read, base and Q20/Q30 totals."""

def quality_score_distribution(
    input: _Path,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    include_no_calls: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """QualityScoreDistribution: the histogram, one row per quality: QUALITY, COUNT_OF_Q
    and COUNT_OF_OQ (only when a counted base has OQ)."""

def mean_quality_by_cycle(
    input: _Path,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """MeanQualityByCycle: the histogram, one row per cycle: CYCLE, MEAN_QUALITY and
    MEAN_ORIGINAL_QUALITY (only when a record has OQ; 0 at a cycle without an OQ mean)."""

def collect_base_distribution_by_cycle(
    input: _Path,
    *,
    aligned_reads_only: bool | None = None,
    pf_reads_only: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectBaseDistributionByCycle: one row per read end and cycle (READ_END, CYCLE,
    PCT_A, PCT_C, PCT_G, PCT_T, PCT_N on a 0 to 100 scale)."""

def collect_duplicate_metrics(
    input: _Path,
    *,
    output: _Path | None = None,
    stop_after: int | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectDuplicateMetrics: one row per library, from the duplicate flags already in
    the file; the metrics file goes to METRICS_FILE."""

def accumulate_quality_yield_metrics(
    input: _Paths,
    *,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """AccumulateQualityYieldMetrics: QualityYieldMetrics files (a list of paths) merged
    into one row."""

def collect_gc_bias_metrics(
    input: _Path,
    *,
    summary_output: _Path,
    reference_sequence: _Path,
    scan_window_size: int | None = None,
    minimum_genome_fraction: float | None = None,
    is_bisulfite_sequenced: bool | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    also_ignore_duplicates: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
) -> list[dict[str, Any]]:
    """CollectGcBiasMetrics: the detail rows (Picard's OUTPUT file), one row per GC bin,
    accumulation unit and READS_USED value; SUMMARY_OUTPUT gets the GcBiasSummaryMetrics
    file."""

def collect_multiple_metrics(
    input: _Path,
    *,
    output: _Path,
    reference_sequence: _Path | None = None,
    stop_after: int | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    file_extension: str | None = None,
    program: str | Sequence[str] | None = None,
    extra_argument: str | Sequence[str] | None = None,
    intervals: _Path | None = None,
    db_snp: _Path | None = None,
    include_unpaired: bool | None = None,
    assume_sorted: bool | None = None,
) -> list[dict[str, Any]]:
    """CollectMultipleMetrics: every PROGRAM on one read of the BAM, each written to OUTPUT
    plus its Picard extension (plus FILE_EXTENSION); Python returns one row per file:
    file and rows."""

def estimate_library_complexity(
    input: _Paths,
    *,
    min_identical_bases: int | None = None,
    max_diff_rate: float | None = None,
    min_mean_quality: int | None = None,
    max_group_ratio: int | None = None,
    barcode_tag: str | None = None,
    read_one_barcode_tag: str | None = None,
    read_two_barcode_tag: str | None = None,
    max_read_length: int | None = None,
    min_group_count: int | None = None,
    read_name_regex: str | None = None,
    optical_duplicate_pixel_distance: int | None = None,
    max_optical_duplicate_set_size: int | None = None,
    output: _Path | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """EstimateLibraryComplexity: one row per library (DuplicationMetrics columns) from
    read sequences alone; the input may be unaligned, or several files. The metrics file
    adds the duplication_group_count histogram."""

def collect_wgs_metrics(
    input: _Path,
    *,
    reference_sequence: _Path,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    coverage_cap: int | None = None,
    locus_accumulation_cap: int | None = None,
    stop_after: int | None = None,
    include_bq_histogram: bool | None = None,
    count_unpaired: bool | None = None,
    sample_size: int | None = None,
    intervals: _Path | None = None,
    use_fast_algorithm: bool | None = None,
    read_length: int | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectWgsMetrics: one row of the 32 WgsMetrics columns; the metrics file adds the
    coverage histogram."""

def collect_raw_wgs_metrics(
    input: _Path,
    *,
    reference_sequence: _Path,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    coverage_cap: int | None = None,
    locus_accumulation_cap: int | None = None,
    stop_after: int | None = None,
    include_bq_histogram: bool | None = None,
    count_unpaired: bool | None = None,
    sample_size: int | None = None,
    intervals: _Path | None = None,
    use_fast_algorithm: bool | None = None,
    read_length: int | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectRawWgsMetrics: CollectWgsMetrics with MINIMUM_MAPPING_QUALITY 0,
    MINIMUM_BASE_QUALITY 3, COVERAGE_CAP 100000 and LOCUS_ACCUMULATION_CAP 200000 as
    defaults, and the RawWgsMetrics class name."""

def collect_wgs_metrics_with_non_zero_coverage(
    input: _Path,
    *,
    reference_sequence: _Path,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    coverage_cap: int | None = None,
    locus_accumulation_cap: int | None = None,
    stop_after: int | None = None,
    include_bq_histogram: bool | None = None,
    count_unpaired: bool | None = None,
    sample_size: int | None = None,
    intervals: _Path | None = None,
    use_fast_algorithm: bool | None = None,
    read_length: int | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectWgsMetricsWithNonZeroCoverage: two rows, CATEGORY WHOLE_GENOME and
    NON_ZERO_REGIONS (the zero-depth loci left out), with the WgsMetrics columns."""

def bed_to_interval_list(
    input: _Path,
    *,
    sequence_dictionary: _Path,
    sort: bool | None = None,
    unique: bool | None = None,
    drop_missing_contigs: bool | None = None,
    keep_length_zero_intervals: bool | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """BedToIntervalList: one row per interval written (contig, start, end, strand, name);
    OUTPUT gets the interval list."""

def collect_hs_metrics(
    input: _Path,
    *,
    target_intervals: _Path | Sequence[_Path],
    reference_sequence: _Path | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    per_target_coverage: _Path | None = None,
    per_base_coverage: _Path | None = None,
    near_distance: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    clip_overlapping_reads: bool | None = None,
    include_indels: bool | None = None,
    coverage_cap: int | None = None,
    sample_size: int | None = None,
    theoretical_sensitivity_output: _Path | None = None,
    allele_fraction: float | Sequence[float] | None = None,
    bait_intervals: _Path | Sequence[_Path],
    bait_set_name: str | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectHsMetrics: one row of the 69 HsMetrics columns per accumulation unit; the
    metrics file adds the coverage and base quality histograms."""

def collect_targeted_pcr_metrics(
    input: _Path,
    *,
    target_intervals: _Path | Sequence[_Path],
    reference_sequence: _Path | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    per_target_coverage: _Path | None = None,
    per_base_coverage: _Path | None = None,
    near_distance: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    clip_overlapping_reads: bool | None = None,
    include_indels: bool | None = None,
    coverage_cap: int | None = None,
    sample_size: int | None = None,
    theoretical_sensitivity_output: _Path | None = None,
    allele_fraction: float | Sequence[float] | None = None,
    amplicon_intervals: _Path,
    custom_amplicon_set_name: str | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectTargetedPcrMetrics: one row of the 62 TargetedPcrMetrics columns per
    accumulation unit, from the CollectHsMetrics collector with MINIMUM_MAPPING_QUALITY
    1, MINIMUM_BASE_QUALITY 0 and CLIP_OVERLAPPING_READS false as defaults."""

def collect_rna_seq_metrics(
    input: _Path,
    *,
    ref_flat: _Path,
    ribosomal_intervals: _Path | None = None,
    strand_specificity: str,
    minimum_length: int | None = None,
    ignore_sequence: str | Sequence[str] | None = None,
    rrna_fragment_percentage: float | None = None,
    metric_accumulation_level: str | Sequence[str] | None = None,
    end_bias_bases: int | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectRnaSeqMetrics: one row of the 30 RnaSeqMetrics columns per accumulation unit;
    the metrics file adds the normalized_position histogram, one normalized_coverage
    column per unit."""

def collect_variant_calling_metrics(
    input: _Path,
    *,
    output: _Path,
    dbsnp: _Path,
    target_intervals: _Path | None = None,
    sequence_dictionary: _Path | None = None,
    gvcf_input: bool | None = None,
    thread_count: int | None = None,
) -> list[dict[str, Any]]:
    """CollectVariantCallingMetrics: the input is a VCF (plain or bgzip). OUTPUT is a
    prefix: PREFIX.variant_calling_summary_metrics and
    PREFIX.variant_calling_detail_metrics; Python returns one row per file: file and
    rows."""

def accumulate_variant_calling_metrics(
    input: _Paths,
    *,
    output: _Path,
) -> list[dict[str, Any]]:
    """AccumulateVariantCallingMetrics: the prefixes of CollectVariantCallingMetrics shards
    merged as Picard merges them; OUTPUT is the prefix of the two files. Python returns
    one row per file: file and rows."""

def genotype_concordance(
    truth_vcf: _Path,
    call_vcf: _Path,
    *,
    output: _Path,
    truth_sample: str | None = None,
    call_sample: str | None = None,
    intervals: _Path | Sequence[_Path] | None = None,
    intersect_intervals: bool | None = None,
    min_gq: int | None = None,
    min_dp: int | None = None,
    output_all_rows: bool | None = None,
    use_vcf_index: bool | None = None,
    missing_sites_hom_ref: bool | None = None,
    ignore_filter_status: bool | None = None,
) -> list[dict[str, Any]]:
    """GenotypeConcordance: the truth and call VCFs compared by position. OUTPUT is a
    prefix: PREFIX.genotype_concordance_summary_metrics, _detail_metrics and
    _contingency_metrics; Python returns one row per file: file and rows."""

def compare_metrics(
    input: _Paths,
    *,
    output: _Path | None = None,
    output_table: _Path | None = None,
    metrics_to_ignore: str | Sequence[str] | None = None,
    metrics_not_required: str | Sequence[str] | None = None,
    metric_allowable_relative_change: str | Sequence[str] | None = None,
    ignore_histogram_differences: bool | None = None,
    key: str | Sequence[str] | None = None,
    verbosity: str | None = None,
) -> list[dict[str, Any]]:
    """CompareMetrics: one row per difference Picard reports between two metrics files:
    KEY, METRIC, VALUE1 and VALUE2 for a value, DIFFERENCE (the line Picard logs) for
    every one; no row when the files are equal. The command line exits 1 when they are
    not."""

def create_sequence_dictionary(
    reference: _Path,
    *,
    genome_assembly: str | None = None,
    uri: str | None = None,
    species: str | None = None,
    truncate_names_at_whitespace: bool | None = None,
    num_sequences: int | None = None,
    alt_names: _Path | None = None,
    create_md5_file: bool | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CreateSequenceDictionary: one row per FASTA sequence (name, length, md5, assembly,
    uri, species, alternate_names); OUTPUT gets the .dict."""

def mark_duplicates(
    input: _Paths,
    *,
    output: _Path,
    metrics_file: _Path,
    remove_duplicates: bool | None = None,
    remove_sequencing_duplicates: bool | None = None,
    assume_sorted: bool | None = None,
    assume_sort_order: str | None = None,
    duplicate_scoring_strategy: str | None = None,
    program_record_id: str | None = None,
    program_group_version: str | None = None,
    program_group_command_line: str | None = None,
    program_group_name: str | None = None,
    add_pg_tag_to_reads: bool | None = None,
    comment: str | Sequence[str] | None = None,
    barcode_tag: str | None = None,
    read_one_barcode_tag: str | None = None,
    read_two_barcode_tag: str | None = None,
    tag_duplicate_set_members: bool | None = None,
    tagging_policy: str | None = None,
    clear_dt: bool | None = None,
    duplex_umi: bool | None = None,
    molecular_identifier_tag: str | None = None,
    read_name_regex: str | None = None,
    optical_duplicate_pixel_distance: int | None = None,
    max_optical_duplicate_set_size: int | None = None,
    compression_level: int | None = None,
    create_index: bool | None = None,
    create_md5_file: bool | None = None,
    tmp_dir: _Path | Sequence[_Path] | None = None,
    max_records_in_ram: int | None = None,
    flow_mode: bool | None = None,
    flow_dup_strategy: str | None = None,
    flow_effective_quality_threshold: int | None = None,
    flow_q_is_known_end: bool | None = None,
    flow_skip_first_n_flows: int | None = None,
    flow_unpaired_end_uncertainty: int | None = None,
    flow_unpaired_start_uncertainty: int | None = None,
    flow_use_end_in_unpaired_reads: bool | None = None,
    flow_use_unpaired_clipped_end: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """MarkDuplicates: one row per library (DuplicationMetrics) plus the set_size
    histograms; OUTPUT gets every record with the duplicate flag set or cleared,
    METRICS_FILE the metrics."""

def collect_sequencing_artifact_metrics(
    input: _Path,
    *,
    output: _Path,
    reference_sequence: _Path,
    intervals: _Path | None = None,
    db_snp: _Path | None = None,
    minimum_quality_score: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_unpaired: bool | None = None,
    include_duplicates: bool | None = None,
    include_non_pf_reads: bool | None = None,
    tandem_reads: bool | None = None,
    use_oq: bool | None = None,
    context_size: int | None = None,
    contexts_to_print: str | Sequence[str] | None = None,
    file_extension: str | None = None,
    stop_after: int | None = None,
    assume_sorted: bool | None = None,
) -> list[dict[str, Any]]:
    """CollectSequencingArtifactMetrics: OUTPUT is a prefix:
    PREFIX.pre_adapter_detail_metrics, .pre_adapter_summary_metrics,
    .bait_bias_detail_metrics, .bait_bias_summary_metrics and .error_summary_metrics
    (plus FILE_EXTENSION); Python returns one row per file: file and rows."""

def convert_sequencing_artifact_to_oxo_g(
    *,
    output: _Path | None = None,
    input_base: _Path | None = None,
    pre_adapter_in: _Path | None = None,
    bait_bias_in: _Path | None = None,
    oxog_out: _Path | None = None,
) -> list[dict[str, Any]]:
    """ConvertSequencingArtifactToOxoG: one CpcgMetrics row per library and C-centred
    context, read from the CollectSequencingArtifactMetrics detail files of INPUT_BASE;
    the file is OUTPUT_BASE.oxog_metrics. Python returns one row per file: file and
    rows."""

def collect_oxo_g_metrics(
    input: _Path,
    *,
    reference_sequence: _Path,
    intervals: _Path | None = None,
    db_snp: _Path | None = None,
    minimum_quality_score: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_insert_size: int | None = None,
    maximum_insert_size: int | None = None,
    include_non_pf_reads: bool | None = None,
    use_oq: bool | None = None,
    context_size: int | None = None,
    contexts: str | Sequence[str] | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectOxoGMetrics: one CpcgMetrics row per sequence context and library: reference
    and alternate base counts at C and G sites split by read and strand, the oxidation
    error rate and its Phred score."""

def collect_coverage_gaps(
    input: _Path,
    *,
    target_intervals: _Path | Sequence[_Path],
    gaps_file: _Path | None = None,
    min_depth: int | None = None,
    minimum_mapping_quality: int | None = None,
    minimum_base_quality: int | None = None,
    clip_overlapping_reads: bool | None = None,
    include_indels: bool | None = None,
    count_unpaired: bool | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectCoverageGaps: one row of the 10 CoverageGapMetrics columns: the target bases
    under MIN_DEPTH as maximal runs, with the depth CollectHsMetrics writes per base;
    the metrics file adds a histogram of gap lengths, and GAPS_FILE a BED of the gaps.
    Not a Picard tool (docs/06)."""

def collect_sex_check_metrics(
    input: _Path,
    *,
    expected_sex: str | None = None,
    minimum_mapping_quality: int | None = None,
    target_intervals: _Path | Sequence[_Path] | None = None,
    x_contig: str | None = None,
    y_contig: str | None = None,
    autosomes: str | Sequence[str] | None = None,
    stop_after: int | None = None,
    output: _Path | None = None,
    assume_sorted: bool | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CollectSexCheckMetrics: one row of the 12 SexCheckMetrics columns: reads per base on
    the X contig, the Y contig and the autosomes, their ratios, the sex they imply
    (FEMALE, MALE or UNKNOWN) and whether it matches EXPECTED_SEX. Not a Picard tool
    (docs/06)."""

def bam_index_stats(
    input: _Path,
    *,
    reference_sequence: _Path | None = None,
) -> list[_BamIndexRow]:
    """BamIndexStats: one row per reference of the BAM header (reference, length, aligned,
    unaligned) read from the .bai or .csi beside INPUT, then a * row whose unaligned
    count is the NoCoordinateCount; the command prints the text Picard prints."""

def check_duplicate_marking(
    input: _Path,
    *,
    output: _Path | None = None,
    mode: str | None = None,
    verbosity: str | None = None,
    reference_sequence: _Path | None = None,
) -> list[dict[str, Any]]:
    """CheckDuplicateMarking: one row per record whose duplicate flag differs from the
    first record of its queryname (read_name, contig, start); no row when every
    queryname agrees. OUTPUT gets the queryname of each such record, one per line. The
    command line exits 1 when any is found."""

# END registry tools

def alignment_summary(
    path: _Path,
    reference: _Path | None = None,
    max_insert_size: int = 100000,
    expected_pair_orientations: Sequence[str] | None = None,
    adapter_sequences: Sequence[str] | None = None,
    collect_alignment_information: bool = True,
) -> list[AlignmentSummaryMetrics]:
    """CollectAlignmentSummaryMetrics at ALL_READS: one row per CATEGORY."""

def quality_yield(
    path: _Path,
    use_original_qualities: bool = True,
    include_secondary_alignments: bool = False,
    include_supplemental_alignments: bool = False,
    reference: _Path | None = None,
) -> QualityYieldMetrics:
    """CollectQualityYieldMetrics. reference decodes a CRAM."""

def hs_metrics(
    path: _Path,
    bait_intervals: _Paths,
    target_intervals: _Paths,
    reference: _Path | None = None,
    bait_set_name: str | None = None,
    minimum_mapping_quality: int = 20,
    minimum_base_quality: int = 20,
    clip_overlapping_reads: bool = True,
    include_indels: bool = False,
    coverage_cap: int = 200,
    near_distance: int = 250,
    sample_size: int = 10000,
    per_base_coverage: _Path | None = None,
) -> HsMetrics:
    """CollectHsMetrics. Interval files are interval lists or BED.

    per_base_coverage writes Picard's PER_BASE_COVERAGE file there.
    """

_Strand = Literal["NONE", "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"]
_OldMetric = Literal["quality_yield", "alignment_summary", "insert_size", "hs_metrics"]

@overload
def collect(
    path: _Path,
    metrics: Sequence[_OldMetric] | None = None,
    reference: _Path | None = None,
    hs: Sequence[tuple[_Paths, _Paths]] | None = None,
    output_prefix: _Path | None = None,
    threads: int | None = None,
    hs_per_base: bool = False,
) -> dict[str, Any]:
    """Several metrics in one pass over the BAM.

    Keys are the metric names; values are what the single-metric functions
    return, and a list of HsMetrics for "hs_metrics".
    """

@overload
def collect(
    path: _Path,
    metrics: Sequence[_CollectorVerb],
    reference: _Path | None = None,
    hs: tuple[_Path, _Path] | Sequence[tuple[_Path, _Path]] | None = None,
    *,
    output: _Path | None = None,
    threads: int | None = None,
    panel: _Path | None = None,
    ref_flat: _Path | None = None,
    strand: _Strand | None = None,
    ribosomal_intervals: _Path | None = None,
) -> Run[None]:
    """The collectors qc runs for verb names (insert-size, gc-bias, hs-metrics, ...),
    on one read, files under output (default ./<bam stem>). Taken when a metric
    name holds a dash or a keyword after * is given."""

@overload
def collect(
    path: _Path,
    metrics: Sequence[str] | None = None,
    reference: _Path | None = None,
    hs: tuple[_Path, _Path] | Sequence[tuple[_Path, _Path]] | None = None,
    *,
    output: _Path,
    threads: int | None = None,
    panel: _Path | None = None,
    ref_flat: _Path | None = None,
    strand: _Strand | None = None,
    ribosomal_intervals: _Path | None = None,
) -> Run[None]: ...

@overload
def collect(
    path: _Path,
    metrics: Sequence[str] | None = None,
    reference: _Path | None = None,
    hs: tuple[_Path, _Path] | Sequence[tuple[_Paths, _Paths]] | None = None,
    output_prefix: _Path | None = None,
    threads: int | None = None,
    hs_per_base: bool = False,
    *,
    output: _Path | None = None,
    panel: _Path | None = None,
    ref_flat: _Path | None = None,
    strand: _Strand | None = None,
    ribosomal_intervals: _Path | None = None,
) -> dict[str, Any] | Run[None]:
    """Names held in a variable of type list[str]: the form follows the names at run time."""

# BEGIN verbs: generated from crates/cli/src/python.rs (UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs)
# Keywords are the command's long flags in snake case; None leaves an option at the command's default.

# The metric names with a dash: collect(metrics=[...]) given one runs the collectors as qc does.
_CollectorVerb = Literal[
    "quality-yield",
    "alignment-summary",
    "insert-size",
    "duplicate-metrics",
    "hs-metrics",
    "gc-bias",
    "rna-seq",
    "sex-check",
    "coverage-gaps",
]

def qc(
    bam: _Path,
    *,
    output: _Path | None = None,
    reference: _Path | None = None,
    preset: str | None = "generic",
    hs: tuple[_Path, _Path] | Sequence[tuple[_Path, _Path]] | None = None,
    panel: _Path | None = None,
    ref_flat: _Path | None = None,
    strand: Literal["NONE", "FIRST_READ_TRANSCRIPTION_STRAND", "SECOND_READ_TRANSCRIPTION_STRAND"] | None = None,
    ribosomal_intervals: _Path | None = None,
    sex: Literal["male", "female"] | None = None,
    sample: str | None = None,
    fail_on: Literal["warn", "fail"] | None = None,
    resume: bool = False,
    force: bool = False,
    threads: int | None = None,
    all_plots: bool = False,
    date: str | None = None,
    report: bool = True,
) -> Run[Verdict]:
    """bammetrics qc: the collectors the preset needs on one read of the BAM, the verdict
    and the page, the files the command writes. report=False leaves out
    PREFIX.report.html. A failed verdict is a value; VerdictError only when fail_on is
    given and reached. An incomplete verdict returns the Run (its exit_code is 11). A
    failed collector raises InputError (its .run holds the rest)."""

def check(
    inputs: _Path | Sequence[_Path],
    *,
    preset: str | None = "generic",
    fail_on: Literal["warn", "fail"] | None = "fail",
) -> Verdict | list[Verdict]:
    """bammetrics check: verdicts from directories of metrics files, verdict.json files or
    qc prefixes. One path with one sample gives a Verdict, else a list. VerdictError
    when a verdict reaches fail_on (fail_on=None never raises). An incomplete verdict is
    returned (its exit_code(fail_on) is 11)."""

def report(
    dirs: _Path | Sequence[_Path],
    *,
    preset: str | None = "generic",
    output: _Path | None = None,
    all_plots: bool = False,
    date: str | None = None,
) -> list[str]:
    """bammetrics report: SAMPLE.verdict.json and SAMPLE.report.html per sample and
    index.html, under output (default: the one directory given). Returns the paths
    written."""

# END verbs

"""Type stubs for the bammetrics Python module (Picard-concordant BAM QC metrics)."""
from os import PathLike
from typing import Any, Sequence, TypedDict

_Path = str | PathLike[str]
_Paths = _Path | Sequence[_Path]

MetricsClass = TypedDict("MetricsClass", {"class": str, "columns": list[str], "rows": list[dict[str, Any]]})
MetricsHistogram = TypedDict("MetricsHistogram", {"key_class": str, "columns": list[str], "rows": list[dict[str, Any]]})
MetricsDoc = TypedDict("MetricsDoc", {"headers": list[str], "classes": list[MetricsClass], "histogram": MetricsHistogram | None})

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
) -> list[InsertSizeMetrics]:
    """CollectInsertSizeMetrics: one row per pair orientation and accumulation unit.

    metric_accumulation_level adds ALL_READS, SAMPLE, LIBRARY or READ_GROUP to the
    default ALL_READS; "null" clears the list, as in Picard."""

def tools() -> list[str]:
    """Picard names of the tools in the registry."""

def read_metrics(path: _Path) -> MetricsDoc:
    """Read a Picard metrics file written by Picard or bammetrics.

    headers are the StringHeader texts; classes has one entry per
    ## METRICS CLASS block with its rows as dicts keyed by column; histogram
    is the ## HISTOGRAM section (key_class, columns, rows) or None. Cells are
    typed from their text: int, float, str or None for an empty field; "?"
    is nan."""

def collect_insert_size_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectInsertSizeMetrics through the registry: one dict per row, keyed by
    Picard's column names. Options are Picard's, lower case: minimum_pct=0.1,
    metric_accumulation_level=["SAMPLE"]. output=PATH also writes the metrics file."""

def collect_alignment_summary_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectAlignmentSummaryMetrics through the registry: one dict per row (a
    CATEGORY of one accumulation unit), keyed by Picard's column names. Options
    are Picard's, lower case: reference_sequence=PATH, is_bisulfite_sequenced=True,
    metric_accumulation_level=["READ_GROUP"], adapter_sequence=[...] (adds to
    Picard's default list). output=PATH also writes the metrics file."""

def collect_quality_yield_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectQualityYieldMetrics through the registry: one dict. Options:
    use_original_qualities, include_secondary_alignments,
    include_supplemental_alignments, stop_after, assume_sorted."""

def quality_score_distribution(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """QualityScoreDistribution. Picard's file has no metrics class, so the
    result is the histogram: one dict per quality, keys QUALITY, COUNT_OF_Q and
    COUNT_OF_OQ (only when a counted base has OQ). Options: aligned_reads_only,
    pf_reads_only, include_no_calls, stop_after, assume_sorted."""

def mean_quality_by_cycle(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """MeanQualityByCycle. The result is the histogram: one dict per cycle,
    keys CYCLE, MEAN_QUALITY and MEAN_ORIGINAL_QUALITY (only when a record has
    OQ; 0 at a cycle without an OQ mean). Options: aligned_reads_only,
    pf_reads_only, stop_after, assume_sorted."""

def collect_base_distribution_by_cycle(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectBaseDistributionByCycle: one dict per read end and cycle (READ_END,
    CYCLE, PCT_A, PCT_C, PCT_G, PCT_T, PCT_N on a 0 to 100 scale). Options:
    aligned_reads_only, pf_reads_only, stop_after, assume_sorted."""

def collect_duplicate_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectDuplicateMetrics from the duplicate flags already in the file: one
    dict per library. output=PATH writes the file Picard writes to METRICS_FILE."""

def accumulate_quality_yield_metrics(input: _Paths, **options: Any) -> list[dict[str, Any]]:
    """AccumulateQualityYieldMetrics: QualityYieldMetrics files (a list of paths)
    merged into one dict. output=PATH writes the merged file."""

def collect_gc_bias_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectGcBiasMetrics: the detail rows (Picard's OUTPUT file), one dict per
    GC bin, accumulation unit and READS_USED value. reference_sequence=PATH is
    required. summary_output=PATH writes the GcBiasSummaryMetrics file and
    output=PATH the detail file. Options: scan_window_size, is_bisulfite_sequenced,
    also_ignore_duplicates, metric_accumulation_level, minimum_genome_fraction
    (accepted, unused by Picard 3.5.0), stop_after, assume_sorted."""

def collect_multiple_metrics(input: _Path, **options: Any) -> list[dict[str, Any]]:
    """CollectMultipleMetrics: every PROGRAM on one read of the BAM, each written
    to output + its Picard extension (+ file_extension). output=PREFIX is
    required. Returns one dict per file written: file and rows. Options:
    program=[...] (adds to Picard's default five; "null" clears), reference_sequence,
    metric_accumulation_level, file_extension, extra_argument=["PROGRAM::KEY=VALUE"],
    stop_after, assume_sorted."""

def estimate_library_complexity(input: _Paths, **options: Any) -> list[dict[str, Any]]:
    """EstimateLibraryComplexity: one dict per library (DuplicationMetrics columns)
    from read sequences alone; the input may be unaligned, or a list of files.
    Options: min_identical_bases, max_diff_rate, min_mean_quality, max_group_ratio,
    max_read_length, min_group_count, read_name_regex, optical_duplicate_pixel_distance,
    max_optical_duplicate_set_size, barcode_tag, read_one_barcode_tag,
    read_two_barcode_tag. output=PATH writes the file with the
    duplication_group_count histogram."""

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
) -> QualityYieldMetrics:
    """CollectQualityYieldMetrics."""

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

def collect(
    path: _Path,
    metrics: Sequence[str] | None = None,
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


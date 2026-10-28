"""The concordance cases: which BAM, which metric, which arguments.

One case is one Picard run and one bammetrics run on the same BAM with the
same settings. Both run_picard.py and test_concordance.py read CASES, so a
case added here gets a reference file and a test.
"""
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out"
FIXTURES = ROOT / "tests" / "data"
PICARD_DIR = CORPUS / "picard"
OURS_DIR = ROOT / "target" / "concordance"
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"
PICARD_VERSION = "3.5.0"


@dataclass(frozen=True)
class Metric:
    picard_tool: str
    cli: str
    picard_args: tuple = ()
    # Columns that identify a row (the rest are compared).
    key_columns: tuple = ()


METRICS = {
    # H= is required; with the no-op Rscript shim nothing is drawn.
    "insert_size": Metric("CollectInsertSizeMetrics", "insert-size", ("H=/dev/null",),
                          ("PAIR_ORIENTATION", "SAMPLE", "LIBRARY", "READ_GROUP")),
    "alignment_summary": Metric("CollectAlignmentSummaryMetrics", "alignment-summary", (),
                                ("CATEGORY",)),
    # One row, no key column.
    "quality_yield": Metric("CollectQualityYieldMetrics", "quality-yield"),
}


@dataclass(frozen=True)
class Case:
    name: str
    bam: Path
    metric: str
    picard_args: tuple = ()
    our_args: tuple = field(default=())

    @property
    def id(self):
        return f"{self.name}.{self.metric}"

    @property
    def picard_txt(self):
        return PICARD_DIR / f"{self.id}.txt"

    @property
    def ours_txt(self):
        return OURS_DIR / f"{self.id}.txt"


def _insert_size_cases():
    t0 = ["t0_even_half_median", "t0_odd", "t0_single_pair", "t0_mode_tie", "t0_all_tandem",
          "t0_cutoff_at", "t0_cutoff_over", "t0_spread", "t0_filters", "t0_no_pairs"]
    cases = [
        Case("t0_fixture", FIXTURES / "insert_size.bam", "insert_size"),
        Case("t0_orientation", FIXTURES / "orientation.bam", "insert_size"),
        Case("t0_orientation_m0", FIXTURES / "orientation.bam", "insert_size",
             ("M=0",), ("-M", "0")),
    ]
    cases += [Case(name, CORPUS / f"{name}.bam", "insert_size") for name in t0]
    cases += [
        Case("t0_single_pair_m0", CORPUS / "t0_single_pair.bam", "insert_size",
             ("M=0",), ("-M", "0")),
        Case("t0_spread_w300", CORPUS / "t0_spread.bam", "insert_size",
             ("W=300",), ("-W", "300")),
        Case("t0_spread_mw3000", CORPUS / "t0_spread.bam", "insert_size",
             ("MW=3000",), ("--min-histogram-width", "3000")),
        Case("t0_spread_dev2", CORPUS / "t0_spread.bam", "insert_size",
             ("DEVIATIONS=2.5",), ("--deviations", "2.5")),
        Case("t0_filters_dups", CORPUS / "t0_filters.bam", "insert_size",
             ("INCLUDE_DUPLICATES=true",), ("--include-duplicates",)),
        Case("t1", CORPUS / "t1.bam", "insert_size"),
        Case("t2", CORPUS / "t2.bam", "insert_size"),
        Case("t2_m0", CORPUS / "t2.bam", "insert_size", ("M=0",), ("-M", "0")),
        Case("t2_dups", CORPUS / "t2.bam", "insert_size",
             ("INCLUDE_DUPLICATES=true",), ("--include-duplicates",)),
        Case("t3", CORPUS / "t3.bam", "insert_size"),
    ]
    return cases


REFERENCE = CORPUS / "reference.fa"


def _alignment_summary_cases():
    """Every BAM with and without R= where the BAM header matches the corpus
    reference. The insert-size t0 BAMs have their own @SQ lengths, so they
    run without R= only."""
    ref = ((f"R={REFERENCE}",), ("-R", str(REFERENCE)))
    # t0_asm has a read running off the end of chrM, which Picard's default
    # STRICT validation rejects before any metric is computed.
    silent = "VALIDATION_STRINGENCY=SILENT"
    cases = [Case(name, CORPUS / f"{name}.bam", "alignment_summary")
             for name in ("t0_filters", "t0_all_tandem", "t0_no_pairs", "t0_spread")]
    cases.append(Case("t0_fixture", FIXTURES / "insert_size.bam", "alignment_summary"))
    for name in ("t0_asm", "t0_asm_cycles", "t0_asm_empty", "t1", "t2", "t3"):
        bam = CORPUS / f"{name}.bam"
        extra = (silent,) if name == "t0_asm" else ()
        cases.append(Case(name, bam, "alignment_summary", extra))
        cases.append(Case(f"{name}_ref", bam, "alignment_summary", extra + ref[0], ref[1]))
    t2, t0 = CORPUS / "t2.bam", CORPUS / "t0_asm.bam"
    adapter = "AGATCGGAAGAGCACACGTCTGAACTCC"
    cases += [
        Case("t2_max_insert_400", t2, "alignment_summary",
             ref[0] + ("MAX_INSERT_SIZE=400",), ref[1] + ("--max-insert-size", "400")),
        # Picard appends to the default [FR]; null empties the list first.
        Case("t2_orient_add_rf", t2, "alignment_summary",
             ref[0] + ("EXPECTED_PAIR_ORIENTATIONS=RF",),
             ref[1] + ("--expected-pair-orientation", "RF")),
        Case("t2_orient_tandem_only", t2, "alignment_summary",
             ref[0] + ("EXPECTED_PAIR_ORIENTATIONS=null", "EXPECTED_PAIR_ORIENTATIONS=TANDEM"),
             ref[1] + ("--expected-pair-orientation", "null",
                       "--expected-pair-orientation", "TANDEM")),
        Case("t2_no_alignment_info", t2, "alignment_summary",
             ref[0] + ("COLLECT_ALIGNMENT_INFORMATION=false",),
             ref[1] + ("--collect-alignment-information", "false")),
        Case("t0_asm_no_adapters", t0, "alignment_summary",
             (silent, "ADAPTER_SEQUENCE=null"), ("--adapter-sequence", "null")),
        Case("t0_asm_own_adapter", t0, "alignment_summary",
             (silent, "ADAPTER_SEQUENCE=null", f"ADAPTER_SEQUENCE={adapter}"),
             ("--adapter-sequence", "null", "--adapter-sequence", adapter)),
    ]
    return cases


def _quality_yield_cases():
    """Every BAM at the defaults, then each option on the BAMs where it matters.
    t0_qy carries OQ tags of every length, a record without QUAL, and secondary
    and supplementary records; t1 to t3 carry OQ on a quarter of their reads."""
    silent = "VALIDATION_STRINGENCY=SILENT"
    names = ["t0_qy", "t0_asm", "t0_asm_empty", "t0_filters", "t0_no_pairs", "t1", "t2", "t3"]
    cases = [Case(name, CORPUS / f"{name}.bam", "quality_yield",
                  (silent,) if name == "t0_asm" else ()) for name in names]
    cases.append(Case("t0_fixture", FIXTURES / "insert_size.bam", "quality_yield"))
    options = {
        "no_oq": (("USE_ORIGINAL_QUALITIES=false",), ("--use-original-qualities", "false")),
        "oq_flag": (("OQ=false",), ("--use-original-qualities", "false")),
        "secondary": (("INCLUDE_SECONDARY_ALIGNMENTS=true",),
                      ("--include-secondary-alignments",)),
        "supplemental": (("INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true",),
                         ("--include-supplemental-alignments",)),
        "both": (("INCLUDE_SECONDARY_ALIGNMENTS=true", "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true"),
                 ("--include-secondary-alignments", "--include-supplemental-alignments")),
    }
    for name in ("t0_qy", "t2"):
        for label, (picard_args, our_args) in options.items():
            cases.append(Case(f"{name}_{label}", CORPUS / f"{name}.bam", "quality_yield",
                              picard_args, our_args))
    return cases


CASES = _insert_size_cases() + _alignment_summary_cases() + _quality_yield_cases()

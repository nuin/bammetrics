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


CASES = _insert_size_cases()

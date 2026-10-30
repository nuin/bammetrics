"""Hand-built edge-case BAMs (edge_bams.py) against Picard 3.5.0, each with several option sets."""
import pytest

from compare import compare_metrics
from conftest import RESULTS
from edge_bams import CASES, build

INSERT_SIZE = "CollectInsertSizeMetrics"
OPTION_SETS = {
    "defaults": [],
    "minimum_pct_0": ["M=0"],
    "minimum_pct_50": ["M=0.5"],
    "width_0": ["W=0"],
    "width_1": ["W=1"],
    "width_100": ["W=100"],
    "deviations_0": ["DEVIATIONS=0"],
    "deviations_1": ["DEVIATIONS=1"],
    "min_width_500": ["MW=500"],
    "include_duplicates": ["INCLUDE_DUPLICATES=true"],
    "read_group_level": ["LEVEL=ALL_READS", "LEVEL=READ_GROUP"],
}


@pytest.fixture(scope="module")
def edge_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("edge")


@pytest.mark.parametrize("options", OPTION_SETS)
@pytest.mark.parametrize("case", CASES)
def test_edge_case(edge_dir, case, options):
    bam = edge_dir / f"{case}.bam"
    if not bam.exists():
        build(edge_dir, case)
    # Picard writes no file when every row falls under MINIMUM_PCT or the input
    # has no pairs, so only cases marked as writing one must compare cells.
    cells = compare_metrics(INSERT_SIZE, bam, ["H=/dev/null", *OPTION_SETS[options]],
                            workdir=edge_dir)
    if CASES[case][2] and options in ("defaults", "minimum_pct_0"):
        assert cells > 0
    RESULTS.append((f"edge {case} {options}", cells, cells, "every cell the same printed string"))

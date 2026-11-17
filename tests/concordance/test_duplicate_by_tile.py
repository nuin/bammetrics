"""CollectDuplicateByTileMetrics (docs/06 section 7) on BAMs Picard 3.5.0
MarkDuplicates marked with TAGGING_POLICY=All. No Picard tool writes per-tile
duplicate rates, so the check is the reconciliation of docs/06 section 7.4:
the ALL row equals the sum of Picard's library rows in UNPAIRED_READS_EXAMINED,
READ_PAIRS_EXAMINED, UNMAPPED_READS, UNPAIRED_READ_DUPLICATES,
READ_PAIR_DUPLICATES and READ_PAIR_OPTICAL_DUPLICATES, and the TILE rows add
up to their LANE rows. Picard runs through the queue (run_picard.picard_env).
"""
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, ROOT
from markdup_bams import build_all as build_markdup
from run_picard import picard_env

PAIRS = {
    "UNPAIRED_READS_EXAMINED": "UNPAIRED_READS",
    "READ_PAIRS_EXAMINED": "PAIRS",
    "UNMAPPED_READS": "UNMAPPED",
    "UNPAIRED_READ_DUPLICATES": "UNPAIRED_DUPLICATES",
    "READ_PAIR_DUPLICATES": "PAIR_DUPLICATES",
    "READ_PAIR_OPTICAL_DUPLICATES": "OPTICAL_DUPLICATES",
}
COUNTS = ("READS", "UNMAPPED", "UNPAIRED_READS", "PAIRS", "UNPAIRED_DUPLICATES",
          "PAIR_DUPLICATES", "OPTICAL_DUPLICATES")


def metrics_rows(path):
    """The rows of the first `## METRICS CLASS` block, as dicts of text."""
    lines = Path(path).read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("## METRICS CLASS"))
    columns = lines[start + 1].split("\t")
    rows = []
    for line in lines[start + 2:]:
        if not line:
            break
        rows.append(dict(zip(columns, line.split("\t"))))
    return rows


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    out = tmp_path_factory.mktemp("md")
    built = build_markdup(out)
    return {"fixture": ROOT / "tests" / "data" / "markdup.bam",
            "markdup": built["markdup"], "markdup_duplex": built["markdup_duplex"],
            "markdup_second": built["markdup_second"]}


@pytest.mark.parametrize("name", ["fixture", "markdup", "markdup_duplex", "markdup_second"])
def test_all_row_reconciles_with_picard_mark_duplicates(inputs, tmp_path, name):
    marked = tmp_path / "picard.bam"
    metrics = tmp_path / "picard.txt"
    subprocess.run(["picard", "MarkDuplicates", f"I={inputs[name]}", f"O={marked}",
                    f"M={metrics}", "TAGGING_POLICY=All", "VALIDATION_STRINGENCY=SILENT"],
                   check=True, capture_output=True, env=picard_env())
    ours = tmp_path / "tiles.txt"
    subprocess.run([str(BAMMETRICS), "duplicate-by-tile", str(marked), "-o", str(ours), "-q"],
                   check=True)
    rows = metrics_rows(ours)
    everything = rows[-1]
    assert everything["UNIT"] == "ALL"
    picard = metrics_rows(metrics)
    for theirs, column in PAIRS.items():
        assert int(everything[column]) == sum(int(r[theirs]) for r in picard), (name, column)
    lanes = [r for r in rows if r["UNIT"] == "LANE"]
    for lane in lanes:
        tiles = [r for r in rows if r["UNIT"] == "TILE" and r["FLOWCELL"] == lane["FLOWCELL"]
                 and r["LANE"] == lane["LANE"]]
        for column in COUNTS:
            assert int(lane[column]) == sum(int(t[column]) for t in tiles), (name, column)
    (unparsed,) = [r for r in rows if r["UNIT"] == "UNPARSED"]
    for column in COUNTS:
        assert int(everything[column]) == int(unparsed[column]) + sum(
            int(lane[column]) for lane in lanes), (name, column)

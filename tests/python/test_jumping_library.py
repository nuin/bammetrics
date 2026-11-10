"""collect_jumping_library_metrics: the one JumpingLibraryMetrics row, several inputs,
options and the rejected ASSUME_SORTED. The counts are Picard 3.5.0's on the seeded
foundation corpus."""
import math
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out" / "foundation"
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")


@pytest.fixture(scope="module")
def corpus():
    if not (CORPUS / "corpus.bam").exists():
        subprocess.run(["python3", str(ROOT / "tests" / "data" / "make_corpus.py"), str(CORPUS)],
                       check=True, capture_output=True)
    return CORPUS


COLUMNS = ["JUMP_PAIRS", "JUMP_DUPLICATE_PAIRS", "JUMP_DUPLICATE_PCT", "JUMP_LIBRARY_SIZE",
           "JUMP_MEAN_INSERT_SIZE", "JUMP_STDEV_INSERT_SIZE", "NONJUMP_PAIRS",
           "NONJUMP_DUPLICATE_PAIRS", "NONJUMP_DUPLICATE_PCT", "NONJUMP_LIBRARY_SIZE",
           "NONJUMP_MEAN_INSERT_SIZE", "NONJUMP_STDEV_INSERT_SIZE", "CHIMERIC_PAIRS",
           "FRAGMENTS", "PCT_JUMPS", "PCT_NONJUMPS", "PCT_CHIMERAS"]


def test_one_row_with_picards_columns(corpus, tmp_path):
    out = tmp_path / "jump.txt"
    rows = bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"), output=str(out))
    assert len(rows) == 1 and list(rows[0]) == COLUMNS
    row = rows[0]
    assert (row["JUMP_PAIRS"], row["JUMP_DUPLICATE_PAIRS"], row["JUMP_LIBRARY_SIZE"]) == (533, 13, 10748)
    assert (row["NONJUMP_PAIRS"], row["CHIMERIC_PAIRS"], row["FRAGMENTS"]) == (5169, 584, 157)
    assert row["PCT_JUMPS"] == pytest.approx(533 / (533 + 5169 + 584))
    assert "## METRICS CLASS\tpicard.analysis.JumpingLibraryMetrics" in out.read_text()


def test_several_inputs_add_up(corpus):
    one = bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"))[0]
    fixture = bammetrics.collect_jumping_library_metrics(FIXTURE)[0]
    both = bammetrics.collect_jumping_library_metrics([str(corpus / "corpus.bam"), FIXTURE])[0]
    assert both["NONJUMP_PAIRS"] == one["NONJUMP_PAIRS"] + fixture["NONJUMP_PAIRS"]
    assert fixture["JUMP_PAIRS"] == 0 and math.isnan(fixture["JUMP_MEAN_INSERT_SIZE"])


def test_options_change_the_counts(corpus):
    base = bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"))[0]
    strict = bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"),
                                                        minimum_mapping_quality=20, tail_limit=2)[0]
    assert strict["NONJUMP_PAIRS"] < base["NONJUMP_PAIRS"]
    short = bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"), chimera_kb_min=500)[0]
    assert short["CHIMERIC_PAIRS"] > base["CHIMERIC_PAIRS"]
    assert short["NONJUMP_PAIRS"] + short["JUMP_PAIRS"] + short["CHIMERIC_PAIRS"] == \
        base["NONJUMP_PAIRS"] + base["JUMP_PAIRS"] + base["CHIMERIC_PAIRS"]
    with pytest.raises(ValueError, match="ASSUME_SORTED"):
        bammetrics.collect_jumping_library_metrics(str(corpus / "corpus.bam"), assume_sorted=True)

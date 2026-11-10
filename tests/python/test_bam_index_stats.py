"""bam_index_stats(): one row per reference of the header, then the `*` row, typed."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "tests" / "data" / "insert_size.bam"
CRAM = ROOT / "tests" / "data" / "insert_size.cram"
REF = ROOT / "tests" / "data" / "ref_100k.fa"
ROWS = [
    {"reference": "chr1", "length": 100000, "aligned": 10, "unaligned": 0},
    {"reference": "*", "length": 0, "aligned": 0, "unaligned": 0},
]


def test_rows_from_the_index():
    rows = bammetrics.bam_index_stats(str(BAM))
    assert rows == ROWS
    assert all(isinstance(row["aligned"], int) for row in rows)


def test_cram_is_counted_with_the_reference():
    assert bammetrics.bam_index_stats(str(CRAM), reference_sequence=str(REF)) == ROWS
    with pytest.raises(OSError, match="REFERENCE_SEQUENCE"):
        bammetrics.bam_index_stats(str(CRAM))


def test_no_index_and_output_raise(tmp_path):
    copy = tmp_path / "noindex.bam"
    copy.write_bytes(BAM.read_bytes())
    with pytest.raises(OSError, match="No index for bam file"):
        bammetrics.bam_index_stats(str(copy))
    with pytest.raises(ValueError, match="Unrecognized option"):
        bammetrics.bam_index_stats(str(BAM), output=str(tmp_path / "stats.txt"))

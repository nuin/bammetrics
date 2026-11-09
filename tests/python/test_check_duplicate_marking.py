"""check_duplicate_marking(): Picard's CheckDuplicateMarking as rows, one per record
whose duplicate flag differs from the first record of its queryname. An empty list
means every queryname agrees; the Picard exit code is not raised."""
from pathlib import Path

import pytest

import bammetrics

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
MARKED = str(DATA / "duplicate_marking.bam")
SORTED = str(DATA / "insert_size.bam")
CRAM = str(DATA / "insert_size.cram")
REF = str(DATA / "ref_100k.fa")


def triples(rows):
    return [(r["read_name"], r["contig"], r["start"]) for r in rows]


def test_bad_records_are_rows_in_file_order():
    assert triples(bammetrics.check_duplicate_marking(MARKED)) == [
        ("r2", "chr1", 1060), ("r4", "chr1", 1060), ("r4", "chr1", 1070),
        ("r5", None, 1030), ("r6", "chr1", 1080), ("r7", None, 0),
    ]


def test_mode_filters_the_records_checked():
    assert len(bammetrics.check_duplicate_marking(MARKED, mode="PRIMARY_ONLY")) == 4
    assert len(bammetrics.check_duplicate_marking(MARKED, mode="PRIMARY_MAPPED_ONLY")) == 2


def test_output_lists_one_queryname_per_bad_record(tmp_path):
    out = tmp_path / "bad.txt"
    rows = bammetrics.check_duplicate_marking(MARKED, output=str(out), mode="PRIMARY_MAPPED_ONLY")
    assert out.read_text() == "r2\nr4\n"
    assert len(rows) == 2


def test_a_coordinate_sorted_file_is_sorted_and_agrees(tmp_path):
    out = tmp_path / "bad.txt"
    assert bammetrics.check_duplicate_marking(SORTED, output=str(out)) == []
    assert out.read_text() == ""
    assert bammetrics.check_duplicate_marking(CRAM, reference_sequence=REF) == []


def test_errors_raise():
    with pytest.raises(ValueError, match="'NOPE' is not a valid value for Mode."):
        bammetrics.check_duplicate_marking(MARKED, mode="NOPE")
    with pytest.raises(Exception, match="Inappropriate call if not paired read"):
        bammetrics.check_duplicate_marking(MARKED, mode="PRIMARY_PROPER_PAIR_ONLY")
    with pytest.raises(ValueError, match="ASSUME_SORTED"):
        bammetrics.check_duplicate_marking(MARKED, assume_sorted=True)

"""compare_s_a_ms() and check_terminator_block(): Picard's CompareSAMs and
CheckTerminatorBlock as rows, on the fixtures of
tests/data/make_compare_sams_fixtures.py and the doctor fixtures. Picard's own
answers for these files are in tests/concordance. Run: pytest tests/python after
maturin develop."""
import os
from pathlib import Path

import pytest

import bammetrics

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
LEFT = str(DATA / "compare_sams_left.bam")
RIGHT = str(DATA / "compare_sams_right.bam")
OTHER = str(DATA / "compare_sams_other.bam")
COUNTS = ["MAPPINGS_MATCH", "MAPPINGS_DIFFER", "UNMAPPED_BOTH", "UNMAPPED_LEFT", "UNMAPPED_RIGHT",
          "MISSING_LEFT", "MISSING_RIGHT", "DUPLICATE_MARKINGS_DIFFER"]


def counts(row):
    return [row[name] for name in COUNTS]


def test_the_row_of_two_files(tmp_path):
    out = tmp_path / "compared.txt"
    (row,) = bammetrics.compare_s_a_ms([LEFT, RIGHT], output=str(out))
    assert counts(row) == [10, 2, 1, 0, 0, 0, 0, 4]
    assert (row["LEFT_FILE"], row["RIGHT_FILE"], row["ARE_EQUAL"]) == (LEFT, RIGHT, False)
    assert "## METRICS CLASS\tpicard.sam.SamComparisonMetric\n" in out.read_text()
    (other,) = bammetrics.compare_s_a_ms([OTHER, LEFT])
    assert counts(other) == [9, 1, 1, 1, 0, 1, 1, 2]


def test_lenient_options_make_the_two_files_equal(tmp_path):
    out = tmp_path / "compared.txt"
    (row,) = bammetrics.compare_s_a_ms(
        [LEFT, RIGHT], lenient_header=True, lenient_dup=True, lenient_low_mq_alignment=True,
        lenient_unknown_mq_alignment=True, compare_mq=True, max_records_in_ram=1,
        tmp_dir=str(tmp_path), output=str(out))
    assert counts(row) == [12, 0, 1, 0, 0, 0, 0, 0] and row["ARE_EQUAL"] is True
    assert "## HISTOGRAM\tjava.lang.String\nBIN\tVALUE\n0,0\t1\n2,3\t1\n255,255\t1\n60,60\t10\n" in out.read_text()
    assert [p.name for p in tmp_path.iterdir()] == ["compared.txt"]
    (strict,) = bammetrics.compare_s_a_ms([LEFT, RIGHT], lenient_low_mq_alignment=True,
                                          low_mq_threshold=2)
    assert strict["MAPPINGS_DIFFER"] == 2


def test_a_cram_against_its_bam():
    bam, cram = str(DATA / "insert_size.bam"), str(DATA / "insert_size.cram")
    (row,) = bammetrics.compare_s_a_ms([bam, cram], reference_sequence=str(DATA / "ref_100k.fa"))
    assert row["MAPPINGS_MATCH"] == 10 and row["ARE_EQUAL"] is True


def test_compare_errors_raise(tmp_path):
    with pytest.raises(ValueError, match="exactly two input files"):
        bammetrics.compare_s_a_ms([LEFT])
    with pytest.raises(ValueError, match="exactly two input files"):
        bammetrics.compare_s_a_ms([LEFT, RIGHT, OTHER])
    with pytest.raises(OSError):
        bammetrics.compare_s_a_ms([LEFT, str(tmp_path / "none.bam")])
    copy = tmp_path / "copy.bam"
    copy.write_bytes(Path(LEFT).read_bytes())
    with pytest.raises(Exception, match="is an input"):
        bammetrics.compare_s_a_ms([str(copy), RIGHT], output=str(copy))
    assert copy.read_bytes() == Path(LEFT).read_bytes()


@pytest.mark.parametrize("name,expected", [
    ("insert_size.bam", "HAS_TERMINATOR_BLOCK"),
    ("doctor/no_eof.bam", "HAS_HEALTHY_LAST_BLOCK"),
    ("doctor/cut.bam", "DEFECTIVE"),
    ("arrays/arrays.vcf.gz", "HAS_TERMINATOR_BLOCK"),
    ("insert_size.cram", "DEFECTIVE"),
])
def test_the_end_of_a_file(name, expected):
    path = str(DATA / name)
    assert bammetrics.check_terminator_block(path) == [{"file": path, "termination": expected}]


def test_terminator_errors_raise(tmp_path):
    with pytest.raises(Exception, match="Cannot read non-existent file"):
        bammetrics.check_terminator_block(str(tmp_path / "none.bam"))
    with pytest.raises(Exception, match="because it is a directory"):
        bammetrics.check_terminator_block(str(tmp_path))
    pipe = tmp_path / "pipe.bam"
    os.mkfifo(pipe)
    with pytest.raises(ValueError, match="not a regular file"):
        bammetrics.check_terminator_block(str(pipe))

"""Python binding of CollectRrbsMetrics on the bisulfite fixture: one entry
per file, the detail rows made as the file is written and filled for Python,
equal to the rows `read_metrics` reads back from both files. Run: pytest
tests/python after maturin develop."""
from pathlib import Path

import bammetrics

HERE = Path(__file__).parents[1]
BAM = str(HERE / "data" / "rrbs.bam")
REFERENCE = str(HERE / "data" / "ref_100k.fa")


def test_both_files_return_the_rows_they_hold(tmp_path):
    prefix = tmp_path / "rrbs"
    files = bammetrics.collect_rrbs_metrics(BAM, metrics_file_prefix=str(prefix),
                                            reference=REFERENCE,
                                            metric_accumulation_level=["ALL_READS", "SAMPLE"])
    names = [Path(f["file"]).name for f in files]
    assert names == ["rrbs.rrbs_summary_metrics", "rrbs.rrbs_detail_metrics"]
    for entry in files:
        written = bammetrics.read_metrics(entry["file"])["classes"][0]["rows"]
        assert len(entry["rows"]) == len(written) > 0
        for returned, read in zip(entry["rows"], written):
            assert list(returned) == list(read)
            for column, value in returned.items():
                if isinstance(value, float):
                    assert abs(value - read[column]) < 1e-6, column
                else:
                    assert value == read[column], column

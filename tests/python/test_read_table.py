"""read_table: columns as lists with the glossary's dtypes, from Picard text and
from the JSON form; qc(also=...) writes the files the command writes."""
import json
import math
import pickle
import subprocess
import sys

import pytest

import bammetrics

PICARD = "tests/data/picard/NA12878/NA12878"


def test_columns_carry_the_glossary_types():
    table = bammetrics.read_table(PICARD + ".insert_size_metrics")
    assert table.name == "picard.analysis.InsertSizeMetrics"
    assert table.tool == "CollectMultipleMetrics"
    assert len(table) == 1
    # MEDIAN_INSERT_SIZE prints 310: a float by the glossary, an int by its text.
    assert table.dtypes["MEDIAN_INSERT_SIZE"] == "float"
    assert table.columns["MEDIAN_INSERT_SIZE"] == [310.0]
    assert isinstance(table.columns["MEDIAN_INSERT_SIZE"][0], float)
    assert table.dtypes["READ_PAIRS"] == "int"
    assert table.columns["READ_PAIRS"] == [5]
    assert table.dtypes["PAIR_ORIENTATION"] == "string"
    assert table.columns["SAMPLE"] == [None]
    by_text = bammetrics.read_metrics(PICARD + ".insert_size_metrics")["classes"][0]["rows"][0]
    assert isinstance(by_text["MEDIAN_INSERT_SIZE"], int)
    assert table.rows()[0]["STANDARD_DEVIATION"] == pytest.approx(8.164966)
    assert list(table.columns) == list(table.dtypes)
    assert table.headers[1].startswith("Started on: ")
    assert repr(table) == "Table('picard.analysis.InsertSizeMetrics', 1 rows, 23 columns)"


def test_the_histogram_is_a_second_table():
    table = bammetrics.read_table(PICARD + ".insert_size_metrics")
    histogram = table.histogram
    assert histogram.name == "java.lang.Integer"
    assert histogram.dtypes == {"insert_size": "int", "All_Reads.fr_count": "float"}
    assert histogram.columns["insert_size"][:2] == [300, 310]
    assert histogram.histogram is None
    same = bammetrics.read_table(PICARD + ".insert_size_metrics", table="histogram")
    assert same.columns == histogram.columns
    # A file with a histogram only gives the histogram.
    quality = bammetrics.read_table(PICARD + ".quality_distribution_metrics")
    assert quality.name == "java.lang.Byte"
    assert quality.dtypes["QUALITY"] == "int"
    assert quality.histogram is None
    by_class = bammetrics.read_table(PICARD + ".insert_size_metrics", table="InsertSizeMetrics")
    assert by_class.columns == table.columns


def test_nan_and_empty_cells(tmp_path):
    path = tmp_path / "x.gc_bias.summary_metrics"
    path.write_text(
        "## htsjdk.samtools.metrics.StringHeader\n# CollectGcBiasMetrics INPUT=x.bam\n\n"
        "## METRICS CLASS\tpicard.analysis.GcBiasSummaryMetrics\n"
        "ACCUMULATION_LEVEL\tWINDOW_SIZE\tAT_DROPOUT\tGC_DROPOUT\tSAMPLE\n"
        "All Reads\t100\t?\t\t\n\n\n"
    )
    table = bammetrics.read_table(path)
    assert table.dtypes["AT_DROPOUT"] == "float"
    assert math.isnan(table.columns["AT_DROPOUT"][0])
    assert table.columns["GC_DROPOUT"] == [None]
    assert table.columns["SAMPLE"] == [None]
    assert table.columns["WINDOW_SIZE"] == [100]


def test_the_json_form_reads_the_same(tmp_path):
    out = tmp_path / "insert.json"
    # Needs target/release/bammetrics (`cargo build --release`), as test_mirror.py does.
    done = subprocess.run(["target/release/bammetrics", "convert",
                           PICARD + ".insert_size_metrics", "-o", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert json.loads(out.read_text())["schema"] == "bammetrics.metrics/1"
    from_json = bammetrics.read_table(out)
    from_text = bammetrics.read_table(PICARD + ".insert_size_metrics")
    assert from_json.columns == from_text.columns
    assert from_json.dtypes == from_text.dtypes
    assert from_json.histogram.columns == from_text.histogram.columns
    assert from_json.headers == from_text.headers


def test_errors():
    with pytest.raises(bammetrics.InputError, match="nope.txt"):
        bammetrics.read_table("tests/data/nope.txt")
    with pytest.raises(bammetrics.InputError, match="not a metrics file"):
        bammetrics.read_table("tests/data/ref_100k.fa")
    with pytest.raises(bammetrics.UsageError, match="has no table 'HsMetrics'|has no table \"HsMetrics\""):
        bammetrics.read_table(PICARD + ".insert_size_metrics", table="HsMetrics")


def test_a_table_pickles():
    table = bammetrics.read_table(PICARD + ".insert_size_metrics")
    # The bytes unpickled are the ones this test just made.
    copy = pickle.loads(pickle.dumps(table))
    assert copy.columns == table.columns
    assert copy.dtypes == table.dtypes
    assert copy.histogram.columns == table.histogram.columns
    assert copy.tool == table.tool


def test_to_pandas_takes_the_dtypes():
    pandas = pytest.importorskip("pandas")
    frame = bammetrics.read_table(PICARD + ".insert_size_metrics").to_pandas()
    assert str(frame["MEDIAN_INSERT_SIZE"].dtype) == "float64"
    assert str(frame["READ_PAIRS"].dtype) == "Int64"
    assert str(frame["PAIR_ORIENTATION"].dtype).startswith("string")
    assert frame["READ_PAIRS"][0] == 5
    assert pandas.isna(frame["SAMPLE"][0])
    assert frame.shape == (1, 23)


def test_to_polars_takes_the_dtypes():
    polars = pytest.importorskip("polars")
    frame = bammetrics.read_table(PICARD + ".insert_size_metrics").to_polars()
    assert frame.schema["MEDIAN_INSERT_SIZE"] == polars.Float64
    assert frame.schema["READ_PAIRS"] == polars.Int64
    assert frame.schema["PAIR_ORIENTATION"] == polars.String
    assert frame["SAMPLE"][0] is None
    assert frame.shape == (1, 23)


def csv_with_nan(tmp_path):
    """The CSV form of a table with NaN, -inf and an empty cell."""
    source = tmp_path / "x.gc_bias.summary_metrics"
    source.write_text(
        "\n## METRICS CLASS\tpicard.analysis.GcBiasSummaryMetrics\n"
        "ACCUMULATION_LEVEL\tWINDOW_SIZE\tAT_DROPOUT\tGC_DROPOUT\n"
        "All Reads\t100\t?\t1.5\nSample\t100\t-?\t\n\n\n"
    )
    out = tmp_path / "x.csv"
    done = subprocess.run(["target/release/bammetrics", "convert", str(source),
                           "--format", "csv", "-o", str(out)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert out.read_bytes() == (b"ACCUMULATION_LEVEL,WINDOW_SIZE,AT_DROPOUT,GC_DROPOUT\r\n"
                                b"All Reads,100,NaN,1.5\r\nSample,100,-inf,\r\n")
    return out


def test_pandas_reads_the_csv_form(tmp_path):
    pandas = pytest.importorskip("pandas")
    frame = pandas.read_csv(csv_with_nan(tmp_path))
    assert str(frame["AT_DROPOUT"].dtype) == "float64"
    assert math.isnan(frame["AT_DROPOUT"][0])
    assert frame["AT_DROPOUT"][1] == -math.inf
    assert frame["WINDOW_SIZE"].tolist() == [100, 100]


def test_polars_reads_the_csv_form(tmp_path):
    polars = pytest.importorskip("polars")
    frame = polars.read_csv(csv_with_nan(tmp_path))
    assert frame.schema["AT_DROPOUT"] == polars.Float64
    assert math.isnan(frame["AT_DROPOUT"][0])
    assert frame["AT_DROPOUT"][1] == -math.inf
    assert frame["GC_DROPOUT"][1] is None


def test_frames_are_not_imported_with_the_module():
    code = ("import sys, bammetrics; bammetrics.read_table(%r); "
            "print('pandas' in sys.modules, 'polars' in sys.modules)"
            % (PICARD + ".insert_size_metrics"))
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert done.stdout.split() == ["False", "False"]


def test_qc_also_writes_the_listed_files(tmp_path):
    prefix = tmp_path / "s"
    run = bammetrics.qc("tests/data/insert_size.bam", output=prefix, also=["json", "markdown"])
    written = {path[len(str(prefix)):] for path in run.files}
    assert ".insert_size_metrics.json" in written
    assert ".report.md" in written
    table = bammetrics.read_table(str(prefix) + ".insert_size_metrics.json")
    assert table.sample == run.sample
    assert table.dtypes["MEDIAN_INSERT_SIZE"] == "float"
    # A str is one format or a list with commas, as on the command line.
    again = bammetrics.qc("tests/data/insert_size.bam", output=prefix, also="junit,multiqc",
                          force=True)
    written = {path[len(str(prefix)):] for path in again.files}
    assert {".junit.xml", ".verdict_mqc.json"} <= written
    assert ".report.md" not in written
    with pytest.raises(bammetrics.UsageError, match="also"):
        bammetrics.qc("tests/data/insert_size.bam", output=prefix, also="parquet", force=True)

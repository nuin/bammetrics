"""CollectAdapterContentMetrics (docs/06 section 7) from Python: the planted
BAM of tests/data/make_adapter_bams.py gives its planted reads, percentages,
mean starts and by-cycle curve exactly; FastQC on `samtools fastq` of the
same BAM gives the same curve; the committed fixture gives its rows by
library. Run: pytest tests/python after maturin develop."""
import importlib.util
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
FIXTURE = str(ROOT / "data" / "adapter_content.bam")
SAMTOOLS = shutil.which("samtools")
FASTQC = shutil.which("fastqc")
LAST_CYCLE = 100 - 12 + 1


def generator():
    spec = importlib.util.spec_from_file_location("make_adapter_bams", ROOT / "data" / "make_adapter_bams.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GEN = generator()
TRUTH = GEN.expected(GEN.PLANTED)


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    out = tmp_path_factory.mktemp("adapter")
    subprocess.run(
        [sys.executable, str(ROOT / "data" / "make_adapter_bams.py"), str(out)],
        check=True, capture_output=True,
    )
    bam = out / "adapter_planted.bam"
    metrics = out / "adapter.txt"
    rows = bammetrics.collect_adapter_content_metrics(str(bam), output=str(metrics))
    return bam, rows, bammetrics.read_metrics(str(metrics))["histogram"]["rows"]


def test_planted_rows_are_exact(planted):
    _, rows, _ = planted
    assert len(rows) == 3 * len(GEN.KEYS)
    for row in rows:
        reads, starts = TRUTH[row["READ_END"]]
        found = starts[row["ADAPTER"]]
        assert row["KMER"] == GEN.KEYS[row["ADAPTER"]]
        assert (row["READS"], row["READS_WITH_ADAPTER"]) == (reads, len(found)), row
        assert row["PCT_READS_WITH_ADAPTER"] == 100 * len(found) / reads
        assert row["PCT_AT_LAST_CYCLE"] == row["PCT_READS_WITH_ADAPTER"]
        if found:
            assert row["MEAN_ADAPTER_START"] == sum(s + 1 for s in found) / len(found)
        else:
            assert math.isnan(row["MEAN_ADAPTER_START"])


def test_read_through_follows_the_planted_inserts(planted):
    _, rows, _ = planted
    through = {r["READ_END"]: r["PAIR_READ_THROUGH"] for r in rows}
    # 20 of the 200 proper pairs have an insert of 30 under 100-base reads.
    assert through["FIRST"] == through["SECOND"] == 0.1
    assert math.isnan(through["UNPAIRED"])


def test_planted_curve_is_exact(planted):
    _, _, table = planted
    assert [r["READ_CYCLE"] for r in table] == list(range(1, LAST_CYCLE + 1))
    for end, (reads, starts) in TRUTH.items():
        for adapter, found in starts.items():
            column = f"All_Reads.{end}.{adapter}"
            for row in table:
                want = 100 * sum(1 for s in found if s < row["READ_CYCLE"]) / reads
                assert row[column] == want, (column, row["READ_CYCLE"])


def fastqc_curves(fastq, out):
    """FastQC's Adapter Content table of `fastq`, ungrouped, by adapter."""
    subprocess.run([FASTQC, "--nogroup", "--extract", "-q", "--outdir", str(out), str(fastq)],
                   check=True, capture_output=True)
    text = (out / f"{fastq.stem}_fastqc" / "fastqc_data.txt").read_text()
    block = text.split(">>Adapter Content")[1].split(">>END_MODULE")[0].splitlines()[1:]
    names = block[0].lstrip("#").split("\t")[1:]
    curves = {name: {} for name in names}
    for line in block[1:]:
        cells = line.split("\t")
        for name, value in zip(names, cells[1:]):
            curves[name][int(cells[0])] = float(value)
    return curves


@pytest.mark.skipif(not (SAMTOOLS and FASTQC), reason="needs samtools and fastqc")
def test_fastqc_gives_the_same_curve(planted, tmp_path):
    bam, _, table = planted
    files = {e: tmp_path / f"{e.lower()}.fq" for e in ("FIRST", "SECOND", "UNPAIRED")}
    subprocess.run([SAMTOOLS, "fastq", "-F", "0x900", "-1", str(files["FIRST"]),
                    "-2", str(files["SECOND"]), "-0", str(files["UNPAIRED"]),
                    "-s", str(tmp_path / "singletons.fq"), str(bam)],
                   check=True, capture_output=True)
    for end, fastq in files.items():
        for adapter, curve in fastqc_curves(fastq, tmp_path).items():
            assert len(curve) == LAST_CYCLE
            for row in table:
                ours = row[f"All_Reads.{end}.{adapter}"]
                assert abs(ours - curve[row["READ_CYCLE"]]) < 5e-7, (end, adapter, row["READ_CYCLE"])


def test_fixture_rows_by_library():
    rows = bammetrics.collect_adapter_content_metrics(
        FIXTURE, metric_accumulation_level="LIBRARY", adapter_sequence=["U=AGATCGGAAGAG"])
    shape = [(r["LIBRARY"], r["READ_END"], r["READS"], r["READS_WITH_ADAPTER"]) for r in rows]
    assert shape == [
        (None, "FIRST", 12, 4), (None, "SECOND", 12, 4), (None, "UNPAIRED", 8, 4),
        ("L1", "FIRST", 8, 4), ("L1", "SECOND", 8, 4), ("L1", "UNPAIRED", 0, 0),
        ("L2", "FIRST", 4, 0), ("L2", "SECOND", 4, 0), ("L2", "UNPAIRED", 8, 4),
    ]


def test_bad_adapter_is_a_usage_error():
    with pytest.raises(bammetrics.UsageError, match="must be A, C, G or T"):
        bammetrics.collect_adapter_content_metrics(FIXTURE, adapter_sequence="x=ACGTNACGTACGT")

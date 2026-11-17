"""CollectMtCopyNumberMetrics (docs/06 section 6) from Python: the planted BAMs
of tests/data/make_mt_bams.py give their copy number exactly, the committed
fixture gives 100 in all and 120 and 80 per library, and the depths equal
samtools depth on the fixture and on corpus t1 (within the 0.5 percent of
docs/06 section 6.5). Run: pytest tests/python after maturin develop."""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
FIXTURE = str(ROOT / "data" / "mt_copy_number.bam")
NO_MT = str(ROOT / "data" / "insert_size.bam")
T1 = ROOT / "corpus" / "out" / "t1.bam"
SAMTOOLS = shutil.which("samtools")
# samtools depth with the tool's filters: base and mapping quality 20, the
# four flags sex_check.rs drops plus supplementary records.
DEPTH = ["depth", "-a", "-q", "20", "-Q", "20", "-G", "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY"]


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    out = tmp_path_factory.mktemp("mt")
    subprocess.run(
        [sys.executable, str(ROOT / "data" / "make_mt_bams.py"), str(out)],
        check=True, capture_output=True,
    )
    return out


@pytest.mark.parametrize(
    ("name", "mt_contig", "mt_bases", "copy_number"),
    [("mt_cn100", "chrM", 50_000, 100.0), ("mt_cn1000", "MT", 500_000, 1000.0),
     ("mt_cn0", "chrM", 0, 0.0)],
)
def test_planted_copy_number_is_exact(planted, name, mt_contig, mt_bases, copy_number):
    row = bammetrics.collect_mt_copy_number_metrics(str(planted / f"{name}.bam"))[0]
    assert row["MT_CONTIG"] == mt_contig
    assert (row["MT_LENGTH"], row["MT_BASES"]) == (500, mt_bases)
    assert (row["AUTOSOMAL_LENGTH"], row["AUTOSOMAL_BASES"]) == (4000, 8000)
    assert row["AUTOSOMAL_MEAN_DEPTH"] == 2.0
    assert row["MT_COPY_NUMBER"] == copy_number
    assert row["READS_USED"] == 80 + mt_bases // 100


def test_no_mitochondrial_contig_is_an_error(planted):
    for bam in (planted / "mt_none.bam", NO_MT):
        with pytest.raises(bammetrics.UsageError, match="no mitochondrial contig"):
            bammetrics.collect_mt_copy_number_metrics(str(bam))


def test_fixture_rows_by_library():
    rows = bammetrics.collect_mt_copy_number_metrics(
        FIXTURE, metric_accumulation_level="LIBRARY")
    assert [(r["LIBRARY"], r["MT_COPY_NUMBER"]) for r in rows] == [
        (None, 100.0), ("L1", 120.0), ("L2", 80.0)]


def depth_sums(bam, contigs):
    """Bases and positions per contig from one samtools depth -a run (the
    fixture has no index, so no region is asked for)."""
    text = subprocess.run([SAMTOOLS, *DEPTH, str(bam)],
                          check=True, capture_output=True, text=True).stdout
    sums = {contig: (0, 0) for contig in contigs}
    for line in text.splitlines():
        contig, _, depth = line.split("\t")
        if contig in sums:
            bases, positions = sums[contig]
            sums[contig] = (bases + int(depth), positions + 1)
    return sums


def check_against_samtools(bam, mt, autosomes, tolerance):
    row = bammetrics.collect_mt_copy_number_metrics(str(bam), autosomes=autosomes)[0]
    sums = depth_sums(bam, [mt, *autosomes])
    mt_depth = sums[mt][0] / sums[mt][1]
    auto_depth = sum(s for s, _ in (sums[a] for a in autosomes)) / sum(
        n for _, n in (sums[a] for a in autosomes))
    assert row["MT_MEAN_DEPTH"] == pytest.approx(mt_depth, rel=tolerance)
    assert row["AUTOSOMAL_MEAN_DEPTH"] == pytest.approx(auto_depth, rel=tolerance)
    assert row["MT_COPY_NUMBER"] == pytest.approx(2 * mt_depth / auto_depth, rel=tolerance)
    return row, sums


@pytest.mark.skipif(SAMTOOLS is None, reason="needs samtools")
def test_fixture_depths_equal_samtools_depth():
    row, sums = check_against_samtools(FIXTURE, "chrM", ["chr1", "chr2"], 0)
    assert (row["MT_BASES"], row["MT_LENGTH"]) == sums["chrM"]


@pytest.mark.skipif(SAMTOOLS is None or not T1.exists(), reason="needs samtools and make corpus")
def test_corpus_t1_within_half_a_percent_of_samtools_depth():
    check_against_samtools(T1, "chrM", ["chr1", "chr2"], 0.005)

"""CollectIndexDepthMetrics (docs/06 section 10) against samtools: no Picard program
exists, so MAPPED and UNMAPPED_PLACED must equal samtools idxstats exactly on the
foundation corpus, the exome slice, every corpus tier BAM and the HG002 exome when
present; ESTIMATED_MEAN_DEPTH must stay within the tolerance of docs/06 section 10.4 of
samtools coverage meandepth on corpus t1 (and on the HG002 exome with
BAMMETRICS_SCALE=1, a whole read of 10 GB).

    pytest tests/concordance/test_index_depth.py
"""
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from cases import CORPUS, ROOT

BIN = ROOT / "target" / "release" / "bammetrics"
SAMTOOLS = shutil.which("samtools")
HG002 = Path(os.environ.get("BAMMETRICS_HG002",
                            Path.home() / "data/bammetrics-scale/exome/HG002.exome.bam"))
# docs/06 section 10.4: the estimate against samtools coverage meandepth per contig.
T1_TOLERANCE = 0.005
# On a real exome the estimate counts duplicates (4.9% of the HG002 mapped records),
# secondary records and soft-clipped bases that samtools coverage drops: it runs 5.5% to
# 7.0% above on the 25 primary contigs of HG002.
HG002_TOLERANCE = 0.10

pytestmark = pytest.mark.skipif(SAMTOOLS is None, reason="samtools is not on PATH")


def tier_bams():
    bams = [CORPUS / "foundation" / "corpus.bam", CORPUS / "exome" / "exome.bam"]
    bams += sorted(CORPUS.glob("t*.bam"))
    return [b for b in bams if Path(f"{b}.bai").exists()]


def index_depth(bam):
    """The rows of `bammetrics index-depth --format tsv` as dicts of text."""
    text = subprocess.run([str(BIN), "index-depth", str(bam), "--format", "tsv"],
                          check=True, capture_output=True, text=True).stdout
    lines = text.splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def idxstats(bam):
    text = subprocess.run([SAMTOOLS, "idxstats", str(bam)],
                          check=True, capture_output=True, text=True).stdout
    return [line.split("\t") for line in text.splitlines()]


def check_counts(bam):
    rows = index_depth(bam)
    contigs = [(r["CONTIG"], r["LENGTH"], r["MAPPED"], r["UNMAPPED_PLACED"])
               for r in rows if r["CONTIG"]]
    expected = idxstats(bam)
    assert contigs == [tuple(row) for row in expected if row[0] != "*"], bam.name
    star = [row[3] for row in expected if row[0] == "*"]
    assert rows[-1]["NO_COORDINATE"] == star[0], bam.name
    return len(contigs)


def test_counts_equal_samtools_idxstats(corpus, exome):
    bams = tier_bams()
    assert len(bams) >= 20, bams
    assert sum(check_counts(bam) for bam in bams) > 0


@pytest.mark.skipif(not HG002.exists(), reason="no HG002 exome under ~/data")
def test_hg002_counts_equal_samtools_idxstats():
    assert check_counts(HG002) == 86


def coverage(bam):
    """samtools coverage at its defaults: meandepth per contig."""
    text = subprocess.run([SAMTOOLS, "coverage", str(bam)],
                          check=True, capture_output=True, text=True).stdout
    return {row[0]: float(row[6]) for row in (line.split("\t") for line in text.splitlines()[1:])}


def check_depth(bam, tolerance, contigs=None):
    rows = {r["CONTIG"]: float(r["ESTIMATED_MEAN_DEPTH"])
            for r in index_depth(bam) if r["CONTIG"]}
    depths = coverage(bam)
    checked = 0
    for contig, mean in depths.items():
        if (contigs and contig not in contigs) or mean == 0:
            continue
        assert abs(rows[contig] - mean) / mean <= tolerance, (contig, rows[contig], mean)
        checked += 1
    return checked


def test_t1_depth_against_samtools_coverage(corpus):
    assert check_depth(CORPUS / "t1.bam", T1_TOLERANCE) == 3


@pytest.mark.skipif(not (HG002.exists() and os.environ.get("BAMMETRICS_SCALE")),
                    reason="BAMMETRICS_SCALE=1 and the HG002 exome run the whole read")
def test_hg002_depth_against_samtools_coverage():
    primary = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}
    assert check_depth(HG002, HG002_TOLERANCE, primary) == 25

"""CollectIndexDepthMetrics (docs/06 section 7) from Python: the planted BAMs of
tests/data/make_index_depth_bams.py give their counts, read length, sex call and
copy number exactly through a .bai and through a .csi, and MAPPED and
UNMAPPED_PLACED equal samtools idxstats on every tests/data BAM with an index.
Run: pytest tests/python after maturin develop."""
import math
import shutil
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
DATA = ROOT / "data"
sys.path.insert(0, str(DATA))
import make_index_depth_bams as plant  # noqa: E402


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    out = tmp_path_factory.mktemp("index_depth")
    subprocess.run([sys.executable, str(DATA / "make_index_depth_bams.py"), str(out)],
                   check=True, capture_output=True)
    return out


def split(rows):
    """Contig rows by name, group rows by group, the ALL row."""
    contigs = {r["CONTIG"]: r for r in rows if r["CONTIG"] is not None}
    groups = {r["GROUP"]: r for r in rows if r["CONTIG"] is None}
    return contigs, groups, groups.pop("ALL")


def per_base(mapped, length):
    return mapped / length if length else math.nan


def same(a, b):
    return (math.isnan(a) and math.isnan(b)) or a == b


@pytest.mark.parametrize("index", ["bai", "csi"])
@pytest.mark.parametrize("name", list(plant.PLANTED))
def test_planted_counts_are_exact(planted, name, index):
    contigs, no_coordinate, _, read_length, sex, copy_number = plant.PLANTED[name]
    bam = planted / f"{name}.bam" if index == "bai" else planted / "csi" / f"{name}.bam"
    rows = bammetrics.collect_index_depth_metrics(str(bam))
    by_contig, groups, all_row = split(rows)
    assert list(by_contig) == [c[0] for c in contigs]
    auto = [c for c in contigs if bammetrics_group(c[0]) == "AUTOSOME"]
    auto_per_base = per_base(sum(c[2] for c in auto), sum(c[1] for c in auto))
    for contig, length, mapped, unmapped in contigs:
        row = by_contig[contig]
        assert (row["LENGTH"], row["MAPPED"], row["UNMAPPED_PLACED"]) == (length, mapped, unmapped)
        assert row["READ_LENGTH"] == read_length
        depth = (math.nan if length == 0 else 0.0 if mapped == 0
                 else math.nan if read_length is None else mapped * read_length / length)
        assert same(row["ESTIMATED_MEAN_DEPTH"], depth), contig
        assert same(row["RATIO_TO_AUTOSOMES"], per_base(mapped, length) / auto_per_base
                    if auto_per_base else math.nan), contig
    assert all_row["NO_COORDINATE"] == no_coordinate
    assert all_row["MAPPED"] == sum(c[2] for c in contigs)
    assert all_row["UNMAPPED_PLACED"] == sum(c[3] for c in contigs)
    assert all_row["INFERRED_SEX"] == sex
    if copy_number is None:
        assert math.isnan(all_row["MT_COPY_NUMBER"])
    else:
        assert all_row["MT_COPY_NUMBER"] == copy_number
    assert set(groups) == {"AUTOSOME", "X", "Y", "MT", "ALT", "DECOY", "HLA", "OTHER"}


def bammetrics_group(name):
    """The group of a planted contig name, as docs/06 section 7.2 lists them."""
    bare = name.removeprefix("chr")
    if bare in {str(i) for i in range(1, 23)}:
        return "AUTOSOME"
    return {"X": "X", "Y": "Y", "M": "MT", "MT": "MT"}.get(bare, "OTHER")


def test_read_length_sample(planted):
    bam = str(planted / "mixed_lengths.bam")
    assert bammetrics.collect_index_depth_metrics(bam)[0]["READ_LENGTH"] == 137.5
    rows = bammetrics.collect_index_depth_metrics(bam, read_length_sample=10)
    assert rows[0]["READ_LENGTH"] == 100.0
    assert rows[0]["ESTIMATED_MEAN_DEPTH"] == 40 * 100.0 / 8192
    rows = bammetrics.collect_index_depth_metrics(bam, read_length=250)
    assert rows[0]["READ_LENGTH"] == 250.0


def test_male_fixture_groups():
    rows = bammetrics.collect_index_depth_metrics(str(DATA / "index_depth.bam"))
    _, groups, all_row = split(rows)
    assert {g: (r["LENGTH"], r["MAPPED"]) for g, r in groups.items()} == {
        "AUTOSOME": (12288, 384), "X": (8192, 128), "Y": (4096, 64), "MT": (1024, 512),
        "ALT": (2048, 16), "DECOY": (1024, 8), "HLA": (4096, 8), "OTHER": (2048, 4)}
    assert groups["X"]["RATIO_TO_AUTOSOMES"] == 0.5
    assert groups["Y"]["RATIO_TO_AUTOSOMES"] == 0.5
    assert (all_row["INFERRED_SEX"], all_row["MT_COPY_NUMBER"]) == ("MALE", 32.0)
    assert all_row["REASON"] is None


def test_refusals(tmp_path):
    with pytest.raises(bammetrics.BammetricsError, match=r"\.crai\) holds no read counts"):
        bammetrics.collect_index_depth_metrics(str(DATA / "insert_size.cram"))
    with pytest.raises(bammetrics.BammetricsError, match="no .bai or .csi beside the BAM"):
        bammetrics.collect_index_depth_metrics(str(DATA / "markdup.bam"))
    with pytest.raises(bammetrics.UsageError):
        bammetrics.collect_index_depth_metrics(str(DATA / "index_depth.bam"),
                                               read_length_sample=0)


def indexed_fixtures(tmp_path):
    """Every tests/data BAM with an index: beside it, or its Picard .bai in
    tests/data/bai copied next to a copy of the BAM."""
    bams = [DATA / "insert_size.bam", DATA / "index_depth.bam"]
    for bai in sorted((DATA / "bai").glob("*.bai")):
        bam = tmp_path / f"{bai.stem}.bam"
        shutil.copyfile(DATA / f"{bai.stem}.bam", bam)
        shutil.copyfile(bai, tmp_path / f"{bai.stem}.bam.bai")
        bams.append(bam)
    return bams


def idxstats(bam):
    """samtools idxstats: {name: (length, mapped, unmapped)} and the * count."""
    rows = [line.split("\t") for line in pysam.idxstats(str(bam)).splitlines()]
    counts = {name: (int(length), int(m), int(u)) for name, length, m, u in rows if name != "*"}
    star = [int(u) for name, _, _, u in rows if name == "*"]
    return counts, star[0] if star else 0


def test_counts_equal_samtools_idxstats(tmp_path):
    bams = indexed_fixtures(tmp_path)
    assert len(bams) == 8
    for bam in bams:
        counts, star = idxstats(bam)
        by_contig, _, all_row = split(bammetrics.collect_index_depth_metrics(str(bam)))
        assert {c: (r["LENGTH"], r["MAPPED"], r["UNMAPPED_PLACED"])
                for c, r in by_contig.items()} == counts, bam.name
        assert all_row["NO_COORDINATE"] == star, bam.name

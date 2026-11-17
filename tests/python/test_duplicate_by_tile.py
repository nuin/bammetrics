"""CollectDuplicateByTileMetrics (docs/06 section 7) from Python.

The planted BAMs of tests/data/make_tile_bams.py give every count of every
row exactly. The ALL row equals the CollectDuplicateMetrics and MarkDuplicates
totals (READ_PAIR_DUPLICATES, UNPAIRED_READ_DUPLICATES,
READ_PAIR_OPTICAL_DUPLICATES and the examined counts) on tests/data/markdup.bam,
the markdup_bams.py BAMs and the complexity_bams.py BAMs after bammetrics
MarkDuplicates TAGGING_POLICY=All. Run: pytest tests/python after maturin develop.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests" / "data"
FIXTURE = str(DATA / "duplicate_by_tile.bam")
MARKDUP = str(DATA / "markdup.bam")
SAMTOOLS = shutil.which("samtools")
COUNTS = ("READS", "UNMAPPED", "UNPAIRED_READS", "PAIRS", "UNPAIRED_DUPLICATES",
          "PAIR_DUPLICATES", "DUPLICATES", "OPTICAL_DUPLICATES")
EXTENT = ("MIN_X", "MAX_X", "MIN_Y", "MAX_Y")
# The OUTLIER_LIMIT tolerance of docs/06 section 7.4: a sum of doubles.
LIMIT_TOLERANCE = 1e-12


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    out = tmp_path_factory.mktemp("tiles")
    subprocess.run([sys.executable, str(DATA / "make_tile_bams.py"), str(out)],
                   check=True, capture_output=True)
    return out


def by_unit(rows):
    tiles = {f"{r['LANE']}:{r['TILE']}": r for r in rows if r["UNIT"] == "TILE"}
    lanes = {r["LANE"]: r for r in rows if r["UNIT"] == "LANE"}
    (unparsed,) = [r for r in rows if r["UNIT"] == "UNPARSED"]
    (everything,) = [r for r in rows if r["UNIT"] == "ALL"]
    return tiles, lanes, unparsed, everything


def median(values):
    values = sorted(values)
    n = len(values)
    return values[n // 2] if n % 2 else (values[n // 2 - 1] + values[n // 2]) / 2


def test_planted_counts_are_exact(planted):
    truth = json.loads((planted / "truth.json").read_text())
    rows = bammetrics.collect_duplicate_by_tile_metrics(str(planted / "tile_truth.bam"))
    assert [r["UNIT"] for r in rows] == ["TILE"] * 6 + ["LANE"] + ["TILE"] * 6 + ["LANE"] + [
        "UNPARSED", "ALL"]
    tiles, lanes, unparsed, everything = by_unit(rows)
    assert set(tiles) == set(truth["tiles"])
    for key, want in truth["tiles"].items():
        row = tiles[key]
        assert {k: row[k] for k in COUNTS + EXTENT} == want, key
        assert row["FLOWCELL"] == "HTILEFCXX"
        examined = want["UNPAIRED_READS"] + 2 * want["PAIRS"]
        assert row["PCT_DUPLICATION"] == want["DUPLICATES"] / examined
        assert row["PCT_OPTICAL"] == want["OPTICAL_DUPLICATES"] / want["PAIRS"]
        assert row["TILE_OUTLIER"] == ("Y" if key == truth["hot"] else "N"), key
    for lane, want in truth["lanes"].items():
        assert {k: lanes[lane][k] for k in COUNTS} == want, lane
        pcts = [tiles[k]["PCT_DUPLICATION"] for k in tiles if k.startswith(lane + ":")]
        center = median(pcts)
        mad = median([abs(p - center) for p in pcts])
        assert abs(lanes[lane]["OUTLIER_LIMIT"] - (center + 3 * 1.4826 * mad)) < LIMIT_TOLERANCE
        assert lanes[lane]["TILE_OUTLIER"] is None
        assert lanes[lane]["MIN_X"] == min(tiles[k]["MIN_X"] for k in tiles if k[0] == lane)
    assert {k: unparsed[k] for k in COUNTS} == truth["unparsed"]
    assert all(unparsed[k] is None for k in EXTENT + ("LANE", "TILE"))
    assert {k: everything[k] for k in COUNTS} == truth["all"]
    assert everything["NOTES"] == (
        f"{truth['unparsed']['READS']} reads have names without a lane and tile under "
        "READ_NAME_REGEX; they are on the UNPARSED row")


def test_a_lower_outlier_sd_flags_more_tiles(planted):
    rows = bammetrics.collect_duplicate_by_tile_metrics(
        str(planted / "tile_truth.bam"), outlier_sd=0.5)
    flagged = sorted(f"{r['LANE']}:{r['TILE']}" for r in rows if r["TILE_OUTLIER"] == "Y")
    assert "1:1104" in flagged and len(flagged) > 1


def test_without_dt_tags_optical_is_empty_and_noted(planted):
    truth = json.loads((planted / "truth.json").read_text())
    rows = bammetrics.collect_duplicate_by_tile_metrics(str(planted / "tile_untagged.bam"))
    tiles, _, _, everything = by_unit(rows)
    assert all(r["OPTICAL_DUPLICATES"] is None and r["PCT_OPTICAL"] is None for r in rows)
    want = {k: v for k, v in truth["all"].items() if k != "OPTICAL_DUPLICATES"}
    assert {k: everything[k] for k in want} == want
    assert "TAGGING_POLICY=All" in everything["NOTES"]
    assert "never marked" not in everything["NOTES"]
    assert tiles[truth["hot"]]["TILE_OUTLIER"] == "Y"


def test_an_unmarked_file_is_noted(planted):
    rows = bammetrics.collect_duplicate_by_tile_metrics(str(planted / "tile_unmarked.bam"))
    everything = rows[-1]
    assert everything["DUPLICATES"] == 0 and everything["PAIR_DUPLICATES"] == 0
    assert everything["NOTES"].startswith(
        "no record carries the duplicate flag and no @PG line names a duplicate marker")


def test_committed_fixture_is_the_planted_truth(planted):
    fixture = bammetrics.collect_duplicate_by_tile_metrics(FIXTURE)
    fresh = bammetrics.collect_duplicate_by_tile_metrics(str(planted / "tile_truth.bam"))
    assert fixture == fresh


def totals(rows):
    """Library rows of a DuplicationMetrics file summed."""
    keys = ("UNPAIRED_READS_EXAMINED", "READ_PAIRS_EXAMINED", "UNMAPPED_READS",
            "UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES", "READ_PAIR_OPTICAL_DUPLICATES")
    return {k: sum(r[k] for r in rows) for k in keys}


def reconcile(bam, workdir):
    """MarkDuplicates TAGGING_POLICY=All on `bam`; the ALL row of the marked BAM
    against the MarkDuplicates metrics and CollectDuplicateMetrics on it."""
    marked = workdir / (Path(bam).stem + ".marked.bam")
    metrics = workdir / (Path(bam).stem + ".dup.txt")
    md = totals(bammetrics.mark_duplicates(
        str(bam), output=str(marked), metrics_file=str(metrics), tagging_policy="All"))
    cdm = totals(bammetrics.collect_duplicate_metrics(str(marked)))
    everything = bammetrics.collect_duplicate_by_tile_metrics(str(marked))[-1]
    assert everything["UNIT"] == "ALL"
    ours = {
        "UNPAIRED_READS_EXAMINED": everything["UNPAIRED_READS"],
        "READ_PAIRS_EXAMINED": everything["PAIRS"],
        "UNMAPPED_READS": everything["UNMAPPED"],
        "UNPAIRED_READ_DUPLICATES": everything["UNPAIRED_DUPLICATES"],
        "READ_PAIR_DUPLICATES": everything["PAIR_DUPLICATES"],
        "READ_PAIR_OPTICAL_DUPLICATES": everything["OPTICAL_DUPLICATES"],
    }
    if ours["READ_PAIR_OPTICAL_DUPLICATES"] is None:
        # A BAM with no duplicate gets no DT tag, and MarkDuplicates writes
        # no @PG line into a BAM without records: the column stays empty.
        assert md["READ_PAIR_DUPLICATES"] == md["UNPAIRED_READ_DUPLICATES"] == 0
        ours["READ_PAIR_OPTICAL_DUPLICATES"] = 0
    # CollectDuplicateMetrics reads no optical duplicate: its column is 0.
    assert ours == dict(cdm, READ_PAIR_OPTICAL_DUPLICATES=ours["READ_PAIR_OPTICAL_DUPLICATES"])
    for key in ("UNPAIRED_READ_DUPLICATES", "READ_PAIR_DUPLICATES",
                "READ_PAIR_OPTICAL_DUPLICATES"):
        assert ours[key] == md[key], (bam, key)
    return ours


def test_fixture_reconciles_with_mark_duplicates(tmp_path):
    ours = reconcile(MARKDUP, tmp_path)
    assert ours["READ_PAIR_OPTICAL_DUPLICATES"] == 6
    assert ours["READ_PAIR_DUPLICATES"] == 7


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    if SAMTOOLS is None:
        pytest.skip("samtools is not on PATH")
    sys.path.insert(0, str(ROOT / "tests" / "concordance"))
    import complexity_bams
    import markdup_bams
    out = tmp_path_factory.mktemp("generated")
    (out / "markdup").mkdir()
    (out / "complexity").mkdir()
    md = markdup_bams.build_all(out / "markdup")
    cx = complexity_bams.build_all(out / "complexity")
    return md, cx


MARKDUP_CASES = ("markdup", "markdup_queryname", "markdup_second", "markdup_duplex",
                 "markdup_duplex_queryname", "markdup_text")


@pytest.mark.parametrize("name", MARKDUP_CASES)
def test_markdup_bams_reconcile(generated, tmp_path, name):
    ours = reconcile(generated[0][name], tmp_path)
    assert ours["READ_PAIR_DUPLICATES"] > 0


# complexity_bams.py BAMs MarkDuplicates refuses, as Picard does: a record
# without RG, and samtools' name order, which htsjdk's writer rejects.
REFUSED = {"illumina": "getReadGroup", "unaligned_queryname": "out of order"}


def test_complexity_bams_reconcile(generated, tmp_path):
    compared = []
    for name, bam in sorted(generated[1].items()):
        if name in REFUSED:
            with pytest.raises(bammetrics.BammetricsError, match=REFUSED[name]):
                reconcile(bam, tmp_path)
            # The flags already in the file still reconcile.
            cdm = totals(bammetrics.collect_duplicate_metrics(str(bam)))
            everything = bammetrics.collect_duplicate_by_tile_metrics(str(bam))[-1]
            assert everything["PAIR_DUPLICATES"] == cdm["READ_PAIR_DUPLICATES"]
            assert everything["UNPAIRED_DUPLICATES"] == cdm["UNPAIRED_READ_DUPLICATES"]
            continue
        reconcile(bam, tmp_path)
        compared.append(name)
    assert len(compared) == len(generated[1]) - len(REFUSED), compared

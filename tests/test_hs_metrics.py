"""Python binding of hs_metrics on the committed fixture (five FR pairs of 100
bases on chr1, 100000 bp, Q40; the first pair at 1001-1100 and 1201-1300),
and on the corpus t1 BAM, against Picard 3.5.0's own numbers for it
(`make corpus`). Run: pytest tests/ after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

HERE = Path(__file__).parent
FIXTURE = str(HERE / "data" / "insert_size.bam")
CORPUS = HERE / "corpus" / "out"


def write_panel(directory):
    path = directory / "panel.interval_list"
    path.write_text("@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n")
    return str(path)


def test_fixture_fields_and_histograms(tmp_path):
    panel = write_panel(tmp_path)
    m = bammetrics.hs_metrics(FIXTURE, panel, [panel], sample_size=100)
    assert m.bait_set == "panel"
    assert (m.total_reads, m.on_bait_bases, m.near_bait_bases, m.off_bait_bases) == (10, 100, 100, 800)
    assert m.mean_target_coverage == 1.0
    assert m.hs_library_size is None
    assert m.sample is None
    assert m.high_quality_coverage[1] == 100
    assert m.unfiltered_base_quality[40] == 100
    fields = m.as_dict()
    assert list(fields)[:3] == ["bait_set", "bait_territory", "bait_design_efficiency"]
    assert len(fields) == 69
    with pytest.raises(AttributeError):
        m.no_such_field


def test_bed_and_options(tmp_path):
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t1000\t1100\tt1\n")
    m = bammetrics.hs_metrics(FIXTURE, str(bed), str(bed), bait_set_name="kit",
                              near_distance=0, sample_size=100)
    assert m.bait_set == "kit"
    assert m.near_bait_bases == 0


def test_missing_intervals_and_file(tmp_path):
    with pytest.raises(ValueError):
        bammetrics.hs_metrics(FIXTURE, [], [])
    with pytest.raises(OSError):
        bammetrics.hs_metrics("/nonexistent.bam", write_panel(tmp_path), write_panel(tmp_path))


@pytest.mark.skipif(not (CORPUS / "t1.bam").exists(), reason="corpus not generated")
def test_t1_matches_picard():
    m = bammetrics.hs_metrics(str(CORPUS / "t1.bam"), str(CORPUS / "baits.interval_list"),
                              str(CORPUS / "targets.interval_list"),
                              reference=str(CORPUS / "reference.fa"))
    assert (m.bait_territory, m.target_territory, m.on_bait_bases) == (40783, 29316, 3303305)
    assert m.median_target_coverage == 72
    assert round(m.fold_80_base_penalty, 6) == 1.229678
    assert round(m.het_snp_sensitivity, 6) == 0.994345
    assert m.het_snp_q == 22

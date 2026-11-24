"""SplitVcfs, MakeSitesOnlyVcf and MakeVcfSampleNameMap from Python (wave W16): the
registry functions, their rows, and the two steps of a variant chain."""
from pathlib import Path

import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
CALLS = DATA / "vcf_writer" / "calls.vcf"
TYPES = DATA / "vcf_tools" / "types.vcf"


def records(path):
    import gzip
    data = Path(path).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return [line for line in data.decode().splitlines() if not line.startswith("#")]


def test_split_vcfs_returns_its_counts(tmp_path):
    rows = bm.split_vcfs(TYPES, snp_output=tmp_path / "s.vcf.gz",
                         indel_output=tmp_path / "i.vcf.gz", strict=False)
    # A `*` alternate has the reference's length: htsjdk types it SNP.
    assert rows == [{"snp_records": 5, "indel_records": 3, "other_records": 4}]
    assert len(records(tmp_path / "s.vcf.gz")) == rows[0]["snp_records"]
    assert (tmp_path / "i.vcf.gz.tbi").exists()
    with pytest.raises(bm.BammetricsError, match="Found a record with type"):
        bm.split_vcfs(TYPES, snp_output=tmp_path / "a.vcf.gz", indel_output=tmp_path / "b.vcf.gz")
    assert not (tmp_path / "a.vcf.gz").exists()


def test_make_sites_only_vcf_keeps_the_samples_asked_for(tmp_path):
    rows = bm.make_sites_only_vcf(CALLS, output=tmp_path / "s.vcf", sample=["S2"],
                                  create_index=False)
    assert rows == [{"records": len(records(CALLS))}]
    header = [line for line in (tmp_path / "s.vcf").read_text().splitlines()
              if line.startswith("#CHROM")][0]
    assert header.endswith("\tFORMAT\tS2")


def test_make_vcf_sample_name_map_rows_are_its_lines(tmp_path):
    truth, other = DATA / "vcf_writer" / "truth.vcf", DATA / "vcf_tools" / "u.vcf"
    rows = bm.make_vcf_sample_name_map([str(truth), str(other)], output=tmp_path / "m.tsv")
    text = (tmp_path / "m.tsv").read_text()
    assert text == "".join(f"{r['input']}\t{r['sample']}\n" for r in rows)
    assert sorted(r["sample"] for r in rows) == ["T", "U"]
    with pytest.raises(bm.BammetricsError, match="single sample"):
        bm.make_vcf_sample_name_map([str(DATA / "vcf_writer" / "call.vcf")],
                                    output=tmp_path / "x.tsv")


def test_the_steps_of_a_variant_chain(tmp_path):
    run = bm.chain(CALLS, [bm.step("make_sites_only_vcf", sample=["S1"]),
                           bm.step("split_vcfs", strict=False)],
                   output=tmp_path / "r")
    assert [(s.verb, s.role, s.status) for s in run.stages] == [
        ("make-sites-only-vcf", "transform", "ok"), ("split-vcfs", "split", "ok")]
    assert [Path(f).name for f in run.stages[1].files] == [
        "r.snp.vcf.gz", "r.snp.vcf.gz.tbi", "r.indel.vcf.gz", "r.indel.vcf.gz.tbi"]
    with pytest.raises(bm.UsageError, match="sample map"):
        bm.chain(CALLS, ["make_vcf_sample_name_map"], output=tmp_path / "x")

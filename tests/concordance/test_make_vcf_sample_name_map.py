"""MakeVcfSampleNameMap against Picard 3.5.0 (tier A: the TSV byte for byte).

Picard 3.5.0 writes one line per distinct INPUT: the input as given (relative, absolute,
with `./` or `..`, a file:// URI), a tab and the sample, in the iteration order of a
`HashMap<String, String>(INPUT.size())` keyed by the input text. A sample in two inputs
warns and keeps both lines. Inputs: single-sample VCFs cut from truth.vcf (plain, gzip
with and without .tbi) under names chosen to collide in the hash table and to resize it.
Error parity: no sample, two samples, a missing file, a file that is no VCF, an empty
file, an OUTPUT whose directory is missing: both fail and bammetrics writes no file.
"""
import gzip
import os
import shutil
import subprocess

import pytest

from conftest import RESULTS
from vcf_tool_cases import bammetrics, corpus, ours, picard


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("sample_map")
    corpus(root)
    text = (root / "truth.vcf").read_text()
    names = [f"s{i}.vcf" for i in range(40)] + ["a.vcf", "b b.vcf", "Aa.vcf", "BB.vcf",
                                                 "ünï.vcf", "x.vcf.gz", "y.vcf.gz"]
    sub = root / "m"
    sub.mkdir()
    for name in names:
        sample = "S" + name.split(".")[0].replace(" ", "_")
        body = text.replace("\tT\n", f"\t{sample}\n", 1)
        if name.endswith(".gz"):
            (sub / name).write_bytes(gzip.compress(body.encode()))
        else:
            (sub / name).write_text(body)
    with open(sub / "z.vcf.gz", "wb") as out:
        subprocess.run(["bgzip", "-c", str(root / "truth.vcf")], stdout=out, check=True)
    subprocess.run(["tabix", "-p", "vcf", str(sub / "z.vcf.gz")], check=True)
    shutil.copy(root / "truth.vcf", sub / "dup.vcf")
    return root


def run_both(cwd, inputs, tmp):
    sides = {}
    old = os.getcwd()
    os.chdir(cwd)
    try:
        for side, run in (("picard", picard), ("ours", ours)):
            out = tmp / f"{side}.tsv"
            sides[side] = (run("MakeVcfSampleNameMap", *[f"I={i}" for i in inputs], f"O={out}"), out)
    finally:
        os.chdir(old)
    return sides


def assert_same(sides, label):
    (p, pout), (o, oout) = sides["picard"], sides["ours"]
    if p.returncode != 0:
        assert o.returncode != 0, f"{label}: Picard failed, bammetrics did not:\n{p.stderr[-600:]}"
        assert not oout.exists(), f"{label}: bammetrics left {oout.name}"
        return 0
    assert o.returncode == 0, f"{label}: {o.stderr}"
    assert pout.read_bytes() == oout.read_bytes(), (label, pout.read_bytes(), oout.read_bytes())
    return len(pout.read_bytes().splitlines())


SHAPES = {
    "six": ["m/s0.vcf", "m/s1.vcf", "m/a.vcf", "m/Aa.vcf", "m/BB.vcf", "m/x.vcf.gz"],
    "one": ["m/s5.vcf"],
    "resized": [f"m/s{i}.vcf" for i in range(40)],
    "unicode": ["m/ünï.vcf", "m/y.vcf.gz"],
    "path_forms": ["./m/s1.vcf", "m/../m/s2.vcf", "m//s3.vcf", "m/z.vcf.gz"],
    "repeated": ["m/s1.vcf", "m/s1.vcf", "m/s2.vcf", "m/s1.vcf"],
    "duplicate_sample": ["m/dup.vcf", "truth.vcf", "m/z.vcf.gz"],
}


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_relative_inputs(data, tmp_path, shape):
    lines = assert_same(run_both(data, SHAPES[shape], tmp_path), shape)
    RESULTS.append((f"MakeVcfSampleNameMap {shape}", lines, lines, "file equal"))


def test_absolute_and_uri_inputs(data, tmp_path):
    inputs = [str(data / "m" / "s7.vcf"), f"file://{data / 'm' / 's8.vcf'}", "m/s9.vcf"]
    lines = assert_same(run_both(data, inputs, tmp_path), "absolute")
    RESULTS.append(("MakeVcfSampleNameMap absolute and file://", lines, lines, "file equal"))


def test_duplicate_samples_warn(data, tmp_path):
    sides = run_both(data, SHAPES["duplicate_sample"], tmp_path)
    assert_same(sides, "warning")
    for side in ("picard", "ours"):
        assert "Duplicate sample: T. Sample was found in both m/dup.vcf and in truth.vcf" in \
            sides[side][0].stderr, side


@pytest.mark.parametrize("case", ["no_sample", "two_samples", "missing", "not_a_vcf", "empty",
                                  "missing_directory", "second_bad",
                                  "space_in_name", "bad_escape", "other_scheme"])
def test_error_parity(data, tmp_path, case):
    inputs = {
        "no_sample": ["sites_only.vcf.gz"],
        "two_samples": ["call.vcf"],
        "missing": ["m/none.vcf"],
        "not_a_vcf": ["ref.dict"],
        "empty": ["empty.vcf"],
        "missing_directory": ["m/s1.vcf"],
        "second_bad": ["m/s1.vcf", "call.vcf"],
        "space_in_name": ["m/b b.vcf"],
        "bad_escape": ["m/s1%zz.vcf"],
        "other_scheme": ["gs://bucket/s1.vcf"],
    }[case]
    (data / "empty.vcf").write_text("")
    target = tmp_path / ("nodir" if case == "missing_directory" else "")
    sides = run_both(data, inputs, target)
    assert sides["picard"][0].returncode != 0, f"{case}: Picard succeeded"
    assert sides["ours"][0].returncode != 0, f"{case}: bammetrics succeeded"
    assert not sides["ours"][1].exists(), case
    if case in ("no_sample", "two_samples"):
        n = "0" if case == "no_sample" else "2"
        message = f"was expected to contain a single sample but actually contained {n} samples."
        assert message in sides["picard"][0].stderr and message in sides["ours"][0].stderr
    RESULTS.append((f"MakeVcfSampleNameMap error {case}", 1, 1, "both fail"))


def test_not_a_chain_step(data, tmp_path):
    done = bammetrics("chain", str(data / "truth.vcf.gz"), "make-vcf-sample-name-map", "-o",
                      str(tmp_path / "x"))
    assert done.returncode == 2 and "sample map" in done.stderr, done.stderr


def test_bad_verbosity_fails_both(data, tmp_path):
    """Picard parses VERBOSITY as a LogLevel before the tool runs."""
    for value in ("bogus", "debug"):
        sides = {}
        for side, run in (("picard", picard), ("ours", ours)):
            out = tmp_path / f"{side}.{value}.tsv"
            sides[side] = (run("MakeVcfSampleNameMap", f"I={data / 'truth.vcf'}", f"O={out}",
                               f"VERBOSITY={value}"), out)
        assert sides["picard"][0].returncode != 0 and sides["ours"][0].returncode != 0, value
        assert not sides["ours"][1].exists(), value
    RESULTS.append(("MakeVcfSampleNameMap VERBOSITY=bogus", 2, 2, "both fail"))


def test_verb_reads_a_file_uri(data, tmp_path):
    """The verb passes a file: URI to the tool, as the Picard face does."""
    uri = f"file://{data / 'm' / 's8.vcf'}"
    verb = bammetrics("make-vcf-sample-name-map", uri, "-o", str(tmp_path / "v.tsv"))
    assert verb.returncode == 0, verb.stderr
    p = picard("MakeVcfSampleNameMap", f"I={uri}", f"O={tmp_path / 'p.tsv'}")
    assert p.returncode == 0, p.stderr[-600:]
    assert (tmp_path / "v.tsv").read_bytes() == (tmp_path / "p.tsv").read_bytes()
    RESULTS.append(("MakeVcfSampleNameMap verb file://", 1, 1, "file equal"))

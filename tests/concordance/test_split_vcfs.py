"""SplitVcfs against Picard 3.5.0.

Truth: compare_vcf's rule with the header rule of vcf_headers.json, on every VCF of
variant_vcfs.py and types.vcf (vcf_tool_cases.py), plain and gzip input, plain and
.vcf.gz output, every option. A .vcf.gz output's .tbi answers every tabix query as
Picard's does. A plain output gets no .idx (docs/07 decision W16-3): Picard writes one
under the default CREATE_INDEX, and CREATE_INDEX=true given for a plain name exits 2 here.
Error parity: both fail, and bammetrics leaves no output (Picard leaves the files it
began). The chain's split step writes the files the tool writes.
"""
import shutil
from pathlib import Path

import pytest

from conftest import RESULTS
from test_vcf_writer import contigs_of, tabix_regions
from compare import vcf_lines
from vcf_tool_cases import (BAD_OUTPUT_SUFFIXES, BGZF_OUTPUT_SUFFIXES, INPUTS,
                            REVIEW_CONTIGS, bammetrics, compare_files, corpus, files_in, ours,
                            picard, review_inputs)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return corpus(tmp_path_factory.mktemp("split_vcfs"))


def run_split(data, tmp, name, suffix, *args):
    sides = {}
    for side, run in (("picard", picard), ("ours", ours)):
        directory = tmp / side
        directory.mkdir()
        done = run("SplitVcfs", f"I={data / name}", f"SNP_OUTPUT={directory / ('snp' + suffix)}",
                   f"INDEL_OUTPUT={directory / ('indel' + suffix)}", *args)
        sides[side] = (done, directory)
    return sides


def assert_equal_split(sides, suffix, label):
    (p, pd), (o, od) = sides["picard"], sides["ours"]
    assert p.returncode == 0, f"{label}: Picard failed:\n{p.stderr[-1500:]}"
    assert o.returncode == 0, f"{label}: bammetrics failed:\n{o.stderr}"
    pairs = [(pd / f"{kind}{suffix}", od / f"{kind}{suffix}") for kind in ("snp", "indel")]
    problems = compare_files(pairs)
    assert not problems, f"{label}:\n  " + "\n  ".join(problems[:20])
    lines = sum(len(vcf_lines(b)) for _, b in pairs)
    if suffix == ".vcf.gz" and shutil.which("tabix"):
        for a, b in pairs:
            if Path(f"{a}.tbi").exists():
                assert Path(f"{b}.tbi").exists(), f"{label}: no .tbi"
                contigs = contigs_of(vcf_lines(b))
                assert tabix_regions(a, contigs) == tabix_regions(b, contigs), label
    return lines


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("suffix", [".vcf.gz", ".vcf"])
def test_not_strict(data, tmp_path, name, suffix):
    args = ["STRICT=false"] + (["CREATE_INDEX=false"] if suffix == ".vcf" else [])
    sides = run_split(data, tmp_path, name, suffix, *args)
    if sides["picard"][0].returncode != 0:
        assert sides["ours"][0].returncode != 0, f"{name}: Picard failed, bammetrics did not"
        assert files_in(sides["ours"][1]) == [], f"{name}: bammetrics left files"
        RESULTS.append((f"SplitVcfs {name}{suffix} both fail", 1, 1, "both fail"))
        return
    lines = assert_equal_split(sides, suffix, f"{name}{suffix}")
    RESULTS.append((f"SplitVcfs {name} -> {suffix}", lines, lines, "every line equal"))


@pytest.mark.parametrize("name", INPUTS)
def test_strict(data, tmp_path, name):
    """STRICT (the default): a record neither SNP nor indel fails both."""
    sides = run_split(data, tmp_path, name, ".vcf.gz")
    (p, _), (o, od) = sides["picard"], sides["ours"]
    if p.returncode != 0:
        assert o.returncode != 0, f"{name}: Picard failed, bammetrics did not"
        assert files_in(od) == [], f"{name}: bammetrics left {files_in(od)}"
        found = [line for line in p.stderr.splitlines() if "Found a record with type" in line]
        if found:
            kind = found[0].rsplit(" ", 1)[1]
            assert f"Found a record with type {kind}" in o.stderr, (found, o.stderr)
        RESULTS.append((f"SplitVcfs STRICT {name} both fail", 1, 1, "both fail"))
        return
    lines = assert_equal_split(sides, ".vcf.gz", f"strict {name}")
    RESULTS.append((f"SplitVcfs STRICT {name}", lines, lines, "every line equal"))


def test_plain_outputs_take_no_idx(data, tmp_path):
    """The default CREATE_INDEX: Picard writes .idx beside plain outputs; bammetrics
    writes the same VCFs, no .idx, and says so on standard error."""
    sides = run_split(data, tmp_path, "truth.vcf.gz", ".vcf", "STRICT=false")
    lines = assert_equal_split(sides, ".vcf", "default index, plain")
    assert (sides["picard"][1] / "snp.vcf.idx").exists()
    assert files_in(sides["ours"][1]) == ["indel.vcf", "snp.vcf"]
    assert "W16-3" in sides["ours"][0].stderr
    refused = ours("SplitVcfs", f"I={data / 'truth.vcf'}", f"SNP_OUTPUT={tmp_path / 'x.vcf'}",
                   f"INDEL_OUTPUT={tmp_path / 'y.vcf'}", "STRICT=false", "CREATE_INDEX=true")
    assert refused.returncode == 1 and not (tmp_path / "x.vcf").exists(), refused.stderr
    RESULTS.append(("SplitVcfs plain outputs", lines, lines, "VCFs equal, no .idx"))


def test_mixed_output_kinds(data, tmp_path):
    sides = {}
    for side, run in (("picard", picard), ("ours", ours)):
        d = tmp_path / side
        d.mkdir()
        sides[side] = (run("SplitVcfs", f"I={data / 'calls3.vcf.gz'}", f"SNP_OUTPUT={d / 'snp.vcf'}",
                           f"INDEL_OUTPUT={d / 'indel.vcf.gz'}", "STRICT=false",
                           "CREATE_INDEX=false"), d)
    assert sides["picard"][0].returncode == 0 and sides["ours"][0].returncode == 0
    problems = compare_files([(sides["picard"][1] / n, sides["ours"][1] / n)
                              for n in ("snp.vcf", "indel.vcf.gz")])
    assert not problems, problems
    RESULTS.append(("SplitVcfs plain and gzip outputs", 2, 2, "every line equal"))


@pytest.mark.parametrize("args", [
    ["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"], ["CREATE_MD5_FILE=true"],
    ["VALIDATION_STRINGENCY=LENIENT"], ["MAX_RECORDS_IN_RAM=10"], ["TMP_DIR=."],
    ["VERBOSITY=DEBUG"], ["VERBOSITY=ERROR"],
], ids=lambda a: a[0])
def test_options(data, tmp_path, args):
    sides = run_split(data, tmp_path, "calls3.vcf.gz", ".vcf.gz", "STRICT=false", *args)
    lines = assert_equal_split(sides, ".vcf.gz", " ".join(args))
    if "CREATE_MD5_FILE=true" in args:
        # Picard's VCF writers do not read it: neither side writes an .md5.
        assert not list(sides["picard"][1].glob("*.md5"))
        assert not list(sides["ours"][1].glob("*.md5"))
    if "VERBOSITY=DEBUG" in args:
        line = "Found 5 records that didn't match SNP or INDEL"
        assert line in sides["picard"][0].stderr and line in sides["ours"][0].stderr
    RESULTS.append((f"SplitVcfs {' '.join(args)}", lines, lines, "every line equal"))


def test_plain_outputs_hold_picards_bytes(data, tmp_path):
    """SplitVcfs adds no header line that names the run: a plain output is Picard's
    file byte for byte."""
    sides = run_split(data, tmp_path, "calls3.vcf", ".vcf", "STRICT=false",
                      "CREATE_INDEX=false")
    assert_equal_split(sides, ".vcf", "bytes")
    for kind in ("snp", "indel"):
        p = (sides["picard"][1] / f"{kind}.vcf").read_bytes()
        o = (sides["ours"][1] / f"{kind}.vcf").read_bytes()
        assert p == o, kind
    RESULTS.append(("SplitVcfs plain bytes", 2, 2, "files equal"))


def test_sequence_dictionary_indexes_a_vcf_without_contigs(data, tmp_path):
    sides = run_split(data, tmp_path, "nocontigs.vcf", ".vcf.gz", "STRICT=false",
                      f"SEQUENCE_DICTIONARY={data / 'ref.dict'}")
    lines = assert_equal_split(sides, ".vcf.gz", "SEQUENCE_DICTIONARY")
    RESULTS.append(("SplitVcfs SEQUENCE_DICTIONARY", lines, lines, "every line equal"))


@pytest.mark.parametrize("case", ["no_dictionary", "missing_input", "undeclared_filter",
                                  "bad_compression", "missing_dictionary"])
def test_error_parity(data, tmp_path, case):
    name, args = "calls3.vcf.gz", ["STRICT=false"]
    if case == "no_dictionary":
        name = "nocontigs.vcf"
    elif case == "missing_input":
        name = "missing.vcf"
    elif case == "undeclared_filter":
        bad = tmp_path / "bad.vcf"
        bad.write_text((data / "calls3.vcf").read_text().replace("\tPASS\t", "\tnope\t", 1))
        name = bad
    elif case == "bad_compression":
        args += ["COMPRESSION_LEVEL=x"]
    elif case == "missing_dictionary":
        name = "nocontigs.vcf"
        args += [f"SEQUENCE_DICTIONARY={tmp_path / 'none.dict'}"]
    sides = run_split(data, tmp_path, name, ".vcf.gz", *args)
    assert sides["picard"][0].returncode != 0, f"{case}: Picard succeeded"
    assert sides["ours"][0].returncode != 0, f"{case}: bammetrics succeeded"
    assert files_in(sides["ours"][1]) == [], f"{case}: {files_in(sides['ours'][1])}"
    RESULTS.append((f"SplitVcfs error {case}", 1, 1, "both fail"))


def test_bammetrics_refuses_one_file_for_both(data, tmp_path):
    done = ours("SplitVcfs", f"I={data / 'calls3.vcf.gz'}", f"SNP_OUTPUT={tmp_path / 'x.vcf.gz'}",
                f"INDEL_OUTPUT={tmp_path / 'x.vcf.gz'}", "STRICT=false")
    assert done.returncode == 1 and files_in(tmp_path) == [], done.stderr
    done = ours("SplitVcfs", f"I={data / 'calls3.vcf'}", f"SNP_OUTPUT={data / 'calls3.vcf'}",
                f"INDEL_OUTPUT={tmp_path / 'y.vcf.gz'}", "STRICT=false")
    assert done.returncode == 1, done.stderr
    # Picard reads any text but true as false (KNOWN_DIFFERENCES.md); bammetrics exits 1.
    done = ours("SplitVcfs", f"I={data / 'calls3.vcf'}", f"SNP_OUTPUT={tmp_path / 's.vcf.gz'}",
                f"INDEL_OUTPUT={tmp_path / 'i.vcf.gz'}", "STRICT=maybe")
    assert done.returncode == 1 and "expected true or false" in done.stderr, done.stderr


def test_chain_split_equals_the_tool(data, tmp_path):
    """chain SRC make-sites-only-vcf split-vcfs writes what the two tools write staged,
    and split-vcfs alone what Picard writes."""
    src = data / "calls3.vcf.gz"
    chained = bammetrics("chain", str(src), "-q", "make-sites-only-vcf", "split-vcfs",
                         "--strict", "false", "-o", str(tmp_path / "c"))
    assert chained.returncode == 0, chained.stderr
    staged = ours("MakeSitesOnlyVcf", f"I={src}", f"O={tmp_path / 's.vcf.gz'}")
    assert staged.returncode == 0, staged.stderr
    p = picard("SplitVcfs", f"I={tmp_path / 's.vcf.gz'}", f"SNP_OUTPUT={tmp_path / 'p.snp.vcf.gz'}",
               f"INDEL_OUTPUT={tmp_path / 'p.indel.vcf.gz'}", "STRICT=false")
    assert p.returncode == 0, p.stderr[-1000:]
    problems = compare_files([(tmp_path / f"p.{k}.vcf.gz", tmp_path / f"c.{k}.vcf.gz")
                              for k in ("snp", "indel")])
    assert not problems, problems
    alone = bammetrics("chain", str(src), "-q", "split-vcfs", "--strict", "false", "-o",
                       str(tmp_path / "a"))
    assert alone.returncode == 0, alone.stderr
    p = picard("SplitVcfs", f"I={src}", f"SNP_OUTPUT={tmp_path / 'q.snp.vcf.gz'}",
               f"INDEL_OUTPUT={tmp_path / 'q.indel.vcf.gz'}", "STRICT=false")
    problems = compare_files([(tmp_path / f"q.{k}.vcf.gz", tmp_path / f"a.{k}.vcf.gz")
                              for k in ("snp", "indel")])
    assert not problems, problems
    RESULTS.append(("SplitVcfs chain", 4, 4, "chained equals Picard staged"))


@pytest.fixture(scope="module")
def review(tmp_path_factory):
    return review_inputs(tmp_path_factory.mktemp("split_vcfs_review"))


# The W16 tools review: (input, suffix, arguments, Picard's outcome). Both succeed with
# equal files, or both fail and bammetrics leaves no file.
REVIEW_CASES = [
    ("v43.vcf", ".vcf", ["CREATE_INDEX=false"], 1),
    ("nonascii3.vcf", ".vcf", ["CREATE_INDEX=false"], 1),
    ("nonascii1.vcf", ".vcf", ["CREATE_INDEX=false"], 0),
    ("undeclared_info.vcf", ".vcf", ["CREATE_INDEX=false", "VALIDATION_STRINGENCY=LENIENT"], 1),
    ("extra_columns.vcf", ".vcf", ["CREATE_INDEX=false"], 0),
    ("unsorted_extra.vcf", ".vcf", ["CREATE_INDEX=false"], 0),
    ("v42.vcf", ".vcf", ["CREATE_INDEX=false", "SEQUENCE_DICTIONARY=nosq.dict"], 0),
    ("v42.vcf", ".vcf.gz", ["SEQUENCE_DICTIONARY=nosq.dict"], 1),
    ("v42.vcf", ".vcf", ["CREATE_INDEX=false", "SEQUENCE_DICTIONARY=dup.dict"], 1),
    ("v42.vcf", ".vcf", ["CREATE_INDEX=false", "SEQUENCE_DICTIONARY=big.dict"], 1),
    ("v42.vcf", ".vcf", ["CREATE_INDEX=false", "SEQUENCE_DICTIONARY=badname.dict"], 0),
    ("v42.vcf", ".vcf.gz", ["SEQUENCE_DICTIONARY=ok.dict"], 0),
    ("contig0.vcf", ".vcf", ["CREATE_INDEX=false", "SEQUENCE_DICTIONARY=ok.dict"], 0),
] + [(f"contig{i}.vcf", ".vcf", ["CREATE_INDEX=false"], 1 if i < 4 else 0)
     for i in range(len(REVIEW_CONTIGS))]


@pytest.mark.parametrize("case", range(len(REVIEW_CASES)))
def test_review_cases(review, tmp_path, case):
    name, suffix, args, picard_exit = REVIEW_CASES[case]
    args = ["STRICT=false"] + [a.replace("DICTIONARY=", f"DICTIONARY={review}/") for a in args]
    sides = run_split(review, tmp_path, name, suffix, *args)
    p, o = sides["picard"][0], sides["ours"][0]
    assert (p.returncode != 0) == (picard_exit != 0), p.stderr[-800:]
    if picard_exit:
        assert o.returncode != 0, f"{name} {args}: Picard failed, bammetrics did not"
        assert files_in(sides["ours"][1]) == [], f"{name}: bammetrics left files"
        RESULTS.append((f"SplitVcfs review {name} {' '.join(args)}", 1, 1, "both fail"))
        return
    lines = assert_equal_split(sides, suffix, f"{name} {args}")
    RESULTS.append((f"SplitVcfs review {name} {' '.join(args)}", lines, lines,
                    "every line equal"))


@pytest.mark.parametrize("suffix", BAD_OUTPUT_SUFFIXES + BGZF_OUTPUT_SUFFIXES)
def test_output_names(review, tmp_path, suffix):
    sides = run_split(review, tmp_path, "v42.vcf", suffix, "CREATE_INDEX=false")
    p, o = sides["picard"][0], sides["ours"][0]
    good = suffix in BGZF_OUTPUT_SUFFIXES
    assert (p.returncode == 0) == good, p.stderr[-800:]
    if not good:
        assert o.returncode != 0 and files_in(sides["ours"][1]) == [], suffix
        RESULTS.append((f"SplitVcfs SNP_OUTPUT=x{suffix}", 1, 1, "both fail"))
        return
    problems = compare_files([(sides["picard"][1] / f"{k}{suffix}", sides["ours"][1] / f"{k}{suffix}")
                              for k in ("snp", "indel")])
    assert not problems, problems
    RESULTS.append((f"SplitVcfs SNP_OUTPUT=x{suffix}", 2, 2, "every line equal"))


@pytest.mark.parametrize("version", [1, 3])
def test_chain_sites_split_equals_picard_staged_for_any_header_text(review, tmp_path, version):
    """chain X make-sites-only-vcf split-vcfs equals Picard MakeSitesOnlyVcf then
    SplitVcfs, for Latin-1 and wider header text and a VCFv4.3 input."""
    src = review / f"nonascii{version}.vcf"
    chained = bammetrics("chain", str(src), "-q", "make-sites-only-vcf", "split-vcfs",
                         "--strict", "false", "--snp-output", str(tmp_path / "c.snp.vcf"),
                         "--indel-output", str(tmp_path / "c.indel.vcf"))
    assert chained.returncode == 0, chained.stderr
    p = picard("MakeSitesOnlyVcf", f"I={src}", f"O={tmp_path / 'p.vcf'}", "CREATE_INDEX=false")
    assert p.returncode == 0, p.stderr[-800:]
    p = picard("SplitVcfs", f"I={tmp_path / 'p.vcf'}", f"SNP_OUTPUT={tmp_path / 'p.snp.vcf'}",
               f"INDEL_OUTPUT={tmp_path / 'p.indel.vcf'}", "STRICT=false", "CREATE_INDEX=false")
    assert p.returncode == 0, p.stderr[-800:]
    for kind in ("snp", "indel"):
        assert (tmp_path / f"p.{kind}.vcf").read_bytes().split(b"\n#CHROM")[1] == \
            (tmp_path / f"c.{kind}.vcf").read_bytes().split(b"\n#CHROM")[1], kind
    assert not compare_files([(tmp_path / f"p.{k}.vcf", tmp_path / f"c.{k}.vcf")
                              for k in ("snp", "indel")])
    RESULTS.append((f"SplitVcfs chain sites non-ASCII v4.{version}", 2, 2,
                    "chained equals Picard staged"))

"""MakeSitesOnlyVcf against Picard 3.5.0.

Truth: compare_vcf's rule with the header rule of vcf_headers.json on every VCF of
variant_vcfs.py and types.vcf (vcf_tool_cases.py), plain and gzip input, plain and
.vcf.gz output, every SAMPLE shape (none, each sample, several, a name the input lacks,
null clearing the list), the writer options. The .tbi answers every tabix query as
Picard's does; a plain output gets no .idx (docs/07 decision W16-3). Error parity: both
fail and bammetrics leaves no file. The chain's transform writes what the tool writes.
"""
import shutil

import pytest

from compare import vcf_lines
from conftest import RESULTS
from test_vcf_writer import contigs_of, tabix_regions
from vcf_tool_cases import (BAD_OUTPUT_SUFFIXES, BGZF_OUTPUT_SUFFIXES, INPUTS, REVIEW_CONTIGS,
                            bammetrics, compare_files, corpus, files_in, ours, picard,
                            review_inputs)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return corpus(tmp_path_factory.mktemp("sites_only"))


def samples_of(path):
    header = [line for line in vcf_lines(path) if line.startswith("#CHROM")][0]
    return header.split("\t")[9:]


def run(data, tmp, name, suffix, *args):
    sides = {}
    for side, tool in (("picard", picard), ("ours", ours)):
        out = tmp / f"{side}{suffix}"
        sides[side] = (tool("MakeSitesOnlyVcf", f"I={data / name}", f"O={out}", *args), out)
    return sides


def assert_equal(sides, label):
    (p, pout), (o, oout) = sides["picard"], sides["ours"]
    if p.returncode != 0:
        assert o.returncode != 0, f"{label}: Picard failed, bammetrics did not"
        assert not oout.exists(), f"{label}: bammetrics left {oout.name}"
        return 0
    assert o.returncode == 0, f"{label}: bammetrics failed:\n{o.stderr}"
    problems = compare_files([(pout, oout)])
    assert not problems, f"{label}:\n  " + "\n  ".join(problems[:20])
    if str(oout).endswith(".gz") and shutil.which("tabix") and pout.with_name(pout.name + ".tbi").exists():
        contigs = contigs_of(vcf_lines(oout))
        assert tabix_regions(pout, contigs) == tabix_regions(oout, contigs), label
    return len(vcf_lines(oout))


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("suffix", [".vcf.gz", ".vcf"])
def test_sites_only(data, tmp_path, name, suffix):
    args = ["CREATE_INDEX=false"] if suffix == ".vcf" else []
    lines = assert_equal(run(data, tmp_path, name, suffix, *args), f"{name}{suffix}")
    RESULTS.append((f"MakeSitesOnlyVcf {name} -> {suffix}", lines, lines,
                    "every line equal" if lines else "both fail"))


@pytest.mark.parametrize("name", INPUTS)
def test_each_sample_and_all_and_a_missing_one(data, tmp_path, name):
    samples = samples_of(data / name)
    shapes = [[s] for s in samples] + [list(reversed(samples)), samples[:1] + ["NOPE"],
                                       ["NOPE"], ["ZZ", "null"] + samples[-1:]]
    total = 0
    for i, shape in enumerate(shapes):
        case = tmp_path / str(i)
        case.mkdir()
        args = [f"SAMPLE={s}" for s in shape] + ["CREATE_INDEX=false"]
        total += assert_equal(run(data, case, name, ".vcf", *args), f"{name} {shape}")
    RESULTS.append((f"MakeSitesOnlyVcf SAMPLE {name}", total, total, "every line equal"))


def test_haploid_missing_sample_is_one_no_call(data, tmp_path):
    lines = assert_equal(run(data, tmp_path, "haploid_call.vcf", ".vcf", "SAMPLE=NOPE",
                             "SAMPLE=C", "CREATE_INDEX=false"), "haploid")
    RESULTS.append(("MakeSitesOnlyVcf haploid no-call", lines, lines, "every line equal"))


@pytest.mark.parametrize("args", [
    ["COMPRESSION_LEVEL=1"], ["CREATE_MD5_FILE=true"], ["VALIDATION_STRINGENCY=SILENT"],
    ["MAX_RECORDS_IN_RAM=5"], ["S=S2", "S=S1"],
], ids=lambda a: a[0])
def test_options(data, tmp_path, args):
    sides = run(data, tmp_path, "calls3.vcf.gz", ".vcf.gz", *args)
    lines = assert_equal(sides, " ".join(args))
    if "CREATE_MD5_FILE=true" in args:
        assert not list(tmp_path.glob("*.md5")), "Picard's VCF writers write no .md5"
    RESULTS.append((f"MakeSitesOnlyVcf {' '.join(args)}", lines, lines, "every line equal"))


def test_plain_output_is_picards_bytes_and_has_no_idx(data, tmp_path):
    sides = run(data, tmp_path, "calls3.vcf", ".vcf")
    assert_equal(sides, "default index")
    assert sides["picard"][1].read_bytes() == sides["ours"][1].read_bytes()
    assert (tmp_path / "picard.vcf.idx").exists() and not (tmp_path / "ours.vcf.idx").exists()
    RESULTS.append(("MakeSitesOnlyVcf plain bytes", 1, 1, "file equal, no .idx"))


@pytest.mark.parametrize("case", ["no_dictionary", "missing_input", "undeclared_filter",
                                  "bad_genotype_sorted", "bad_genotype_unsorted",
                                  "bad_compression"])
def test_error_parity(data, tmp_path, case):
    name, args = "calls3.vcf.gz", []
    if case == "no_dictionary":
        name = "nocontigs.vcf"
    elif case == "missing_input":
        name = "missing.vcf"
    elif case == "undeclared_filter":
        bad = tmp_path / "bad.vcf"
        bad.write_text((data / "calls3.vcf").read_text().replace("\tPASS\t", "\tnope\t", 1))
        name = bad
    elif case.startswith("bad_genotype"):
        # An allele index past the alleles: htsjdk decodes the genotypes when it cuts
        # them to SAMPLE, or as it reads a file whose samples are out of order.
        text = (data / "types.vcf").read_text().replace("\t0/1:30:10\t", "\t0/7:30:10\t")
        if case == "bad_genotype_sorted":
            text = text.replace("\tS2\tS1", "\tS1\tS2")
            args = ["SAMPLE=S1"]
        bad = tmp_path / "bad.vcf"
        bad.write_text(text)
        name = bad
    elif case == "bad_compression":
        args = ["COMPRESSION_LEVEL=x"]
    sides = run(data, tmp_path, name, ".vcf.gz", *args)
    assert sides["picard"][0].returncode != 0, f"{case}: Picard succeeded"
    assert sides["ours"][0].returncode != 0, f"{case}: bammetrics succeeded"
    assert not sides["ours"][1].exists(), case
    RESULTS.append((f"MakeSitesOnlyVcf error {case}", 1, 1, "both fail"))


def test_unsorted_samples_untouched_genotypes_still_fail_when_decoded(data, tmp_path):
    """Sorted sample names and no SAMPLE: the genotypes stay unread, so a bad one
    passes in both."""
    text = (data / "types.vcf").read_text().replace("\t0/1:30:10\t", "\t0/7:30:10\t")
    bad = tmp_path / "bad.vcf"
    bad.write_text(text.replace("\tS2\tS1", "\tS1\tS2"))
    lines = assert_equal(run(data, tmp_path, bad, ".vcf", "CREATE_INDEX=false"), "lazy")
    assert lines
    RESULTS.append(("MakeSitesOnlyVcf lazy genotypes", lines, lines, "every line equal"))


@pytest.mark.parametrize("name", ["calls3.vcf.gz", "truth.vcf", "gvcf.vcf.gz", "types.vcf.gz"])
def test_chain_equals_the_tool(data, tmp_path, name):
    """chain X make-sites-only-vcf -o r.vcf.gz writes Picard's file."""
    chained = bammetrics("chain", str(data / name), "-q", "make-sites-only-vcf", "-o",
                         str(tmp_path / "r.vcf.gz"))
    p = picard("MakeSitesOnlyVcf", f"I={data / name}", f"O={tmp_path / 'p.vcf.gz'}")
    assert chained.returncode == 0 and p.returncode == 0, (chained.stderr, p.stderr[-800:])
    problems = compare_files([(tmp_path / "p.vcf.gz", tmp_path / "r.vcf.gz")])
    assert not problems, problems
    contigs = contigs_of(vcf_lines(tmp_path / "r.vcf.gz"))
    if shutil.which("tabix"):
        assert tabix_regions(tmp_path / "p.vcf.gz", contigs) == tabix_regions(tmp_path / "r.vcf.gz", contigs)
    with_samples = bammetrics("chain", str(data / name), "-q", "make-sites-only-vcf",
                              *[a for s in samples_of(data / name)[:1] for a in ("--sample", s)],
                              "-o", str(tmp_path / "s.vcf.gz"))
    p = picard("MakeSitesOnlyVcf", f"I={data / name}", f"O={tmp_path / 'q.vcf.gz'}",
               *[f"SAMPLE={s}" for s in samples_of(data / name)[:1]])
    assert with_samples.returncode == 0 and p.returncode == 0
    assert not compare_files([(tmp_path / "q.vcf.gz", tmp_path / "s.vcf.gz")])
    assert files_in(tmp_path).count("r.vcf.gz.tbi") == 1
    RESULTS.append((f"MakeSitesOnlyVcf chain {name}", 2, 2, "chained equals Picard"))


@pytest.fixture(scope="module")
def review(tmp_path_factory):
    return review_inputs(tmp_path_factory.mktemp("sites_only_review"))


# The W16 tools review: each case runs both faces and holds them to each other (both
# succeed with equal files, or both fail and bammetrics leaves no file). The third item
# is Picard's outcome, so a case cannot pass by both tools failing where Picard works.
REVIEW_CASES = [
    ("undeclared_info.vcf", ["VALIDATION_STRINGENCY=LENIENT"], 0),
    ("undeclared_info.vcf", ["VALIDATION_STRINGENCY=SILENT", "SAMPLE=S1"], 0),
    ("undeclared_info.vcf", [], 1),
    ("undeclared_filter.vcf", ["VALIDATION_STRINGENCY=SILENT"], 0),
    ("undeclared_flag.vcf", ["VALIDATION_STRINGENCY=LENIENT"], 0),
    ("undeclared_format.vcf", ["VALIDATION_STRINGENCY=LENIENT", "SAMPLE=S1"], 0),
    ("undeclared_format.vcf", ["SAMPLE=S1"], 1),
    ("extra_columns.vcf", ["SAMPLE=S1"], 0),
    ("extra_columns.vcf", [], 0),
    ("unsorted_extra.vcf", [], 0),
    ("unsorted_extra.vcf", ["SAMPLE=S1"], 0),
    ("v43.vcf", [], 0),
    ("nonascii1.vcf", [], 0),
    ("nonascii3.vcf", ["SAMPLE=S2"], 0),
    ("v42.vcf", ["VERBOSITY=bogus"], 1),
    ("v42.vcf", ["VERBOSITY=debug"], 1),
    ("v42.vcf", ["SAMPLE=a\nb"], 1),
] + [(f"contig{i}.vcf", [], 1 if i < 4 else 0) for i in range(len(REVIEW_CONTIGS))]


@pytest.mark.parametrize("case", range(len(REVIEW_CASES)))
def test_review_cases(review, tmp_path, case):
    name, args, picard_exit = REVIEW_CASES[case]
    sides = run(review, tmp_path, name, ".vcf", "CREATE_INDEX=false", *args)
    assert (sides["picard"][0].returncode != 0) == (picard_exit != 0), sides["picard"][0].stderr[-800:]
    lines = assert_equal(sides, f"{name} {args}")
    RESULTS.append((f"MakeSitesOnlyVcf review {name} {' '.join(args)}", lines, lines,
                    "every line equal" if picard_exit == 0 else "both fail"))


@pytest.mark.parametrize("suffix", BAD_OUTPUT_SUFFIXES + BGZF_OUTPUT_SUFFIXES)
def test_output_names(review, tmp_path, suffix):
    """htsjdk types an output by its name: .vcf (lower case) is text, a gzip suffix in any
    case BGZF, any other name fails before a file is made."""
    sides = run(review, tmp_path, "v42.vcf", suffix, "CREATE_INDEX=false")
    assert (sides["picard"][0].returncode == 0) == (suffix in BGZF_OUTPUT_SUFFIXES)
    lines = assert_equal(sides, suffix)
    RESULTS.append((f"MakeSitesOnlyVcf O=x{suffix}", lines, lines,
                    "every line equal" if lines else "both fail"))

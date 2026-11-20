"""SplitSamByNumberOfReads against Picard 3.5.0.

Picard counts the records of a first read unless TOTAL_READS_IN_INPUT is given, takes
SPLIT_TO_N_FILES or `ceil(total / SPLIT_TO_N_READS)` files of `ceil(total / files)` records,
goes on past that number while the record name stays the same, and writes
`OUT_PREFIX_%04d.<input extension>` through `makeSAMOrBAMWriter` (a BAM under a `.cram`
name for a CRAM input). Every case runs both tools into a directory of their own and
compares every file as test_split_sam_by_library.py does (tier B with `record_bytes=True`).

Inputs: `split_bams.py` (`runs`: names of 1 to 6 records with secondary, supplementary and
unmapped ones, as BAM and SAM text; `queryname`, `unsorted`, `empty`), the foundation
corpus queryname sorted and as written (coordinate: Picard's warning), t1 to t3, the exome
slice, the fixture BAMs, SAM text, gzip SAM, CRAM 3.0 of the foundation corpus and
`insert_size.cram`. Options: SPLIT_TO_N_FILES (2, 3, 7, more files than records),
SPLIT_TO_N_READS (2, 10, 1000), TOTAL_READS_IN_INPUT (equal, smaller, larger),
OUT_PREFIX (a space, `null`, non-ASCII), COMPRESSION_LEVEL, CREATE_INDEX, CREATE_MD5_FILE,
MAX_RECORDS_IN_RAM, REFERENCE_SEQUENCE, VALIDATION_STRINGENCY. Error parity: both counts,
neither, SPLIT_TO_N_FILES=1, a negative total, a pipe without the total, ASSUME_SORTED, a
missing input, a missing output directory, an OUT_PREFIX through a missing directory. The
chain with TOTAL_READS_IN_INPUT against Picard. The distinct-value gate wants three file
counts and three records-per-file sizes.
"""
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from test_add_comments_to_bam import run_tool
from test_split_sam_by_library import (SEEN, SILENT, compare_split, diff_file,  # noqa: F401
                                       records_of, split)

TOOL = "SplitSamByNumberOfReads"
REF = FIXTURES / "ref_100k.fa"

# (files written, records) of every case that compared equal.
COUNTS = []


def compare(source, args, label, outcome="equal"):
    files = compare_split(TOOL, source, args, label, outcome)
    if outcome == "equal":
        COUNTS.append(SEEN[-1])
    return files


@pytest.fixture(scope="module")
def queryname(corpus, tmp_path_factory):
    out = tmp_path_factory.mktemp("qname") / "corpus.qname.bam"
    subprocess.run(["samtools", "sort", "-n", "--no-PG", "-o", str(out), str(corpus / "corpus.bam")],
                   check=True)
    return out


@pytest.mark.parametrize("args", [["N_FILES=2"], ["N_FILES=3"], ["N_FILES=7"], ["N_READS=1000"],
                                  ["N_READS=10"], ["N_FILES=3", "TOTAL_READS=12918"],
                                  ["N_FILES=3", "TOTAL_READS=5000"],
                                  ["N_READS=1000", "TOTAL_READS=40000"]])
def test_the_queryname_corpus(queryname, args):
    compare(queryname, [*args, *SILENT], f"foundation queryname {' '.join(args)}")


@pytest.mark.parametrize("args", [["N_FILES=4"], ["N_READS=50"], ["N_READS=2"], ["N_FILES=2000"]])
def test_runs_of_names(split, args):
    compare(split["runs"], args, f"runs {' '.join(args)}")


@pytest.mark.parametrize("name", ["runs_sam", "queryname", "unsorted", "empty"])
def test_split_inputs(split, name):
    compare(split[name], ["N_FILES=3", *SILENT], name)


@pytest.mark.parametrize("name", ["t1", "t2", "t3", "foundation", "exome", "markdup",
                                  "duplicate_marking", "oxog"])
def test_every_input(corpus, exome, name):
    source = {"t1": CORPUS / "t1.bam", "t2": CORPUS / "t2.bam", "t3": CORPUS / "t3.bam",
              "foundation": corpus / "corpus.bam", "exome": exome / "exome.bam",
              "markdup": FIXTURES / "markdup.bam", "oxog": FIXTURES / "oxog.bam",
              "duplicate_marking": FIXTURES / "duplicate_marking.bam"}[name]
    compare(source, ["N_FILES=3", *SILENT], name)


def test_sam_text_and_cram(corpus, split, tmp_path):
    compare(FIXTURES / "markdup.sam", ["N_READS=10"], "markdup.sam")
    compare(split["names_sam_gz"], ["N_FILES=2", *SILENT], "gzip SAM")
    compare(FIXTURES / "insert_size.cram", ["N_FILES=2", f"R={REF}"], "insert_size.cram")
    cram = tmp_path / "corpus.cram"
    reference = corpus / "ref.fa"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(reference), "-o", str(cram),
                    str(corpus / "corpus.bam")], check=True)
    files = compare(cram, ["N_FILES=3", f"R={reference}", "CREATE_INDEX=true", *SILENT],
                    "foundation CRAM 3.0")
    assert "shard_0001.cram" in files and "shard_0001.cram.bai" in files, files


@pytest.mark.parametrize("args", [["OUT_PREFIX=part a"], ["OUT_PREFIX=null"],
                                  ["OUT_PREFIX=café"], ["COMPRESSION_LEVEL=0"],
                                  ["COMPRESSION_LEVEL=9"], ["CREATE_INDEX=true"],
                                  ["CREATE_MD5_FILE=true"], ["MAX_RECORDS_IN_RAM=1"],
                                  [f"R={REF}"], ["VALIDATION_STRINGENCY=LENIENT"]])
def test_options(args):
    compare(FIXTURES / "markdup.bam", ["N_FILES=3", *args], " ".join(args))


@pytest.mark.parametrize("case", ["both", "none", "one_file", "negative_total", "assume_sorted",
                                  "missing_input", "prefix_directory"])
def test_errors_as_picard(tmp_path, case):
    source, args = FIXTURES / "markdup.bam", ["N_FILES=3"]
    args = {"both": ["N_FILES=3", "N_READS=5"], "none": [], "one_file": ["N_FILES=1"],
            "negative_total": ["N_FILES=3", "TOTAL_READS=-1"],
            "assume_sorted": ["N_FILES=3", "ASSUME_SORTED=true"],
            "prefix_directory": ["N_FILES=3", "OUT_PREFIX=nodir/x"]}.get(case, args)
    if case == "missing_input":
        source = tmp_path / "nope.bam"
    compare_split(TOOL, source, args, f"error {case}", "fail")


def test_a_pipe_needs_the_total(tmp_path):
    """Without TOTAL_READS_IN_INPUT both refuse a stream; with it both split it."""
    source = FIXTURES / "markdup.bam"
    for side in ("picard", "bammetrics"):
        out = tmp_path / side
        out.mkdir()
        with open(source, "rb") as stdin:
            if side == "picard":
                done = run_tool(side, TOOL, ["I=/dev/stdin", f"O={out}", "N_FILES=2"])
            else:
                done = subprocess.run([str(BAMMETRICS), TOOL, "I=/dev/stdin", f"O={out}",
                                       "N_FILES=2"], stdin=stdin, capture_output=True, text=True)
        assert done.returncode == 1, f"{side}: {done.stderr[-500:]}"
        assert not list(out.iterdir()), side
    out = tmp_path / "piped"
    out.mkdir()
    with open(source, "rb") as stdin:
        done = subprocess.run([str(BAMMETRICS), TOOL, "I=/dev/stdin", f"O={out}", "N_FILES=2",
                               "TOTAL_READS=36"], stdin=stdin, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    staged = tmp_path / "staged"
    staged.mkdir()
    run_tool("bammetrics", TOOL, [f"I={source}", f"O={staged}", "N_FILES=2"])
    names = sorted(p.name for p in out.iterdir())
    assert names == sorted(p.name for p in staged.iterdir()) == ["shard_0001.bam", "shard_0002.bam"]
    for name in names:
        assert not diff_file(staged / name, out / name), name
    RESULTS.append((f"{TOOL} pipe", 36, 36, "both refuse without the total; the total splits it"))


def test_the_step_in_a_chain_equals_picard(queryname, tmp_path):
    total = records_of(queryname)
    picard, chained = tmp_path / "picard", tmp_path / "chained"
    picard.mkdir()
    done = run_tool("picard", TOOL, [f"I={queryname}", f"O={picard}", "N_FILES=5", *SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    chain = subprocess.run([str(BAMMETRICS), "chain", str(queryname), "+",
                            "split-sam-by-number-of-reads", "-o", str(chained), "--n-files", "5",
                            "--total-reads", str(total)], capture_output=True, text=True)
    assert chain.returncode == 0, chain.stderr
    names = sorted(p.name for p in picard.iterdir())
    assert names == sorted(p.name for p in chained.iterdir())
    problems = [p for name in names for p in diff_file(picard / name, chained / name)]
    assert not problems, "\n".join(problems[:20])
    RESULTS.append((f"{TOOL} chain on the queryname corpus", total, total,
                    f"{len(names)} files equal Picard's"))


def test_values_are_distinct():
    """The cases write 1, 2, 3, 7 and more files from inputs of at least three sizes."""
    assert len({files for files, _ in COUNTS}) >= 3, COUNTS
    assert len({records for _, records in COUNTS}) >= 3, COUNTS

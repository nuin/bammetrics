"""PositionBasedDownsampleSam against Picard 3.5.0.

Picard reads the input twice: the first read finds each tile's x and y extent from the
read names (`ReadNameParser`), the second keeps the records whose name falls in a circle
of area FRACTION on its tile (outside one of area 1 - FRACTION above 0.5), clears their
duplicate flag and writes them through a presorted `SAMFileWriter` with an `@PG` line.
Tier B through `compare.compare_bam` with `record_bytes=True`.

Inputs: the markdup_bams and complexity_bams BAMs (Illumina names of 5 and 7 fields on
several tiles, names Picard cannot parse: 6 fields, a last field that is no number, plain
names), the MarkDuplicates fixture, the foundation corpus (no location in its names: one
tile -1) and t1; FRACTION 0.1, 0.5 and 0.9 on each. Options: REMOVE_DUPLICATE_INFORMATION
false, READ_NAME_REGEX (a regex with three groups, one whose y group is no int, null),
STOP_AFTER (0, a few records, null), FRACTION 0 and 1, COMPRESSION_LEVEL, CREATE_INDEX,
CREATE_MD5_FILE, SAM text in and out, a CRAM. ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS on
Picard's own output, which fails without it. Error parity: FRACTION missing, below 0,
above 1, not a number, a READ_NAME_REGEX that does not compile, ASSUME_SORTED, a missing
input, standard input (bammetrics refuses it, Picard reads /dev/stdin three times and
fails). The distinct-value gate wants three kept counts.
"""
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from markdup_bams import build_all as build_markdup
from test_add_comments_to_bam import compare, index_and_md5, run_tool

TOOL = "PositionBasedDownsampleSam"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    bams = {"fixture": FIXTURES / "markdup.bam", "t1": CORPUS / "t1.bam"}
    for prefix, build in (("md", build_markdup), ("elc", build_complexity)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = bam
    return bams


NAMES = ["fixture", "md_markdup", "md_markdup_second", "md_markdup_queryname", "elc_illumina",
         "elc_plain_names", "elc_unaligned", "t1"]


def run_case(source, args, label, outcome="equal", **kwargs):
    written = compare(TOOL, source, [*args, SILENT], label, outcome=outcome, **kwargs)
    if written.outcome == "equal":
        SEEN.append(len(written.records))
    return written


@pytest.mark.parametrize("fraction", ["0.1", "0.5", "0.9"])
@pytest.mark.parametrize("name", NAMES)
def test_fractions(built, name, fraction):
    if name not in built:
        pytest.skip(f"{name} not built")
    run_case(built[name], [f"FRACTION={fraction}"], f"{name} FRACTION={fraction}")


@pytest.mark.parametrize("fraction", ["0.1", "0.9"])
def test_foundation(corpus, fraction):
    run_case(corpus / "corpus.bam", [f"FRACTION={fraction}"], f"foundation F={fraction}")


OPTIONS = [
    ["F=0.3", "REMOVE_DUPLICATE_INFORMATION=false"],
    ["F=0.3", r"READ_NAME_REGEX=[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+)"],
    ["F=0.3", r"READ_NAME_REGEX=.*:([0-9]+):([0-9]+):([^:]+)"],
    ["F=0.3", "READ_NAME_REGEX=null"],
    ["F=0.3", "STOP_AFTER=0"],
    ["F=0.3", "STOP_AFTER=15"],
    ["F=0.3", "STOP_AFTER=null"],
    ["F=0"],
    ["F=1"],
    ["F=0.5", "COMPRESSION_LEVEL=1", "MAX_RECORDS_IN_RAM=5"],
]


@pytest.mark.parametrize("args", OPTIONS, ids=lambda a: " ".join(a))
def test_options(built, args):
    run_case(built["md_markdup"], args, f"md_markdup {' '.join(args)}")


def test_index_and_md5(tmp_path):
    index_and_md5(TOOL, [f"I={FIXTURES / 'markdup.bam'}", "F=0.5"], tmp_path, "out.bai")


def test_sam_text_and_cram(built, tmp_path):
    source = built["md_markdup"]
    sam = tmp_path / "in.sam"
    subprocess.run(["samtools", "view", "-h", "--no-PG", "-o", str(sam), str(source)], check=True)
    run_case(sam, ["F=0.4"], "SAM text in")
    run_case(source, ["F=0.4"], "SAM text out", suffix=".sam")
    reference = source.parent / "ref.fa"
    cram = tmp_path / "in.cram"
    subprocess.run(["samtools", "view", "-C", "--no-PG", "-T", str(reference), "-o", str(cram),
                    str(source)], check=True)
    run_case(cram, ["F=0.4", f"R={reference}"], "CRAM", reference=reference)


def test_a_second_run(built, tmp_path):
    first = tmp_path / "first.bam"
    done = run_tool("picard", TOOL, [f"I={built['md_markdup']}", f"O={first}", "F=0.8", SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    run_case(first, ["F=0.5"], "downsampled before", outcome="fail")
    run_case(first, ["F=0.5", "ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS=true"],
             "downsampled before, allowed")


ERRORS = [
    [],
    ["F=-0.1"],
    ["F=1.1"],
    ["F=half"],
    ["F=0.5", "READ_NAME_REGEX=(["],
    ["F=0.5", "ASSUME_SORTED=true"],
    ["F=0.5", "STOP_AFTER=x"],
]


@pytest.mark.parametrize("args", ERRORS, ids=lambda a: " ".join(a) or "no FRACTION")
def test_error_parity(args):
    run_case(FIXTURES / "markdup.bam", args, " ".join(args) or "no FRACTION", outcome="fail")


def test_missing_input(tmp_path):
    run_case(tmp_path / "missing.bam", ["F=0.5"], "missing input", outcome="fail")


def test_standard_input_is_refused(tmp_path):
    out = tmp_path / "out.bam"
    with open(FIXTURES / "markdup.bam", "rb") as source:
        done = subprocess.run([str(BAMMETRICS), TOOL, "I=/dev/stdin",
                               f"O={out}", "F=0.5"], stdin=source, capture_output=True, text=True)
    assert done.returncode != 0 and "standard input" in done.stderr, done.stderr
    assert not out.exists()


def test_values_are_distinct():
    if not SEEN:
        pytest.skip("no case ran")
    assert len(set(SEEN)) >= 3, sorted(set(SEEN))
    RESULTS.append((f"{TOOL} distinct kept counts", len(set(SEEN)), len(set(SEEN)), ""))

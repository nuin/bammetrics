"""CollectQualityYieldMetricsSNVQ against Picard 3.5.0: the one
QualityYieldMetrics row and the INCLUDE_BQ_HISTOGRAM section on the
snvq_umi_bams inputs (qa, qc, qg and qt tags on every base kind, PF and
QC-failed, secondary and supplementary, records without bases), every
option, CRAM and SAM text input (lower-case bases, which htsjdk reads as
upper case), and error parity: a missing SNVQ tag (Picard's
NullPointerException), a tag too short where a base reads it, a character
outside ! to ~, a tag of another type, bases without QUAL, a missing or bad
alternate quality tag, and the corpus BAMs, which carry no SNVQ tag.

The last test checks that every column took at least 3 distinct
non-default values over all Picard files compared.
"""
import subprocess
from collections import defaultdict

import pysam
import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_metrics
from conftest import RESULTS
from run_picard import picard_env
from snvq_umi_bams import build_all

TOOL = "CollectQualityYieldMetricsSNVQ"
SEEN = []
COLUMNS = 31


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build_all(tmp_path_factory.mktemp("snvq"))


def compare(bam, args, label, histogram=False):
    cells = compare_metrics(TOOL, bam, args, workdir=bam.parent, expect_file=True, seen=SEEN)
    if not histogram:
        assert cells == COLUMNS, f"{TOOL} {label}: {cells} cells compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))
    return cells


INPUTS = ["snvq_main", "snvq_short_unread", "snvq_secondary_bare", "snvq_empty"]


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, name):
    compare(built[name], [], name)


OPTIONS = {
    "histogram": ["INCLUDE_BQ_HISTOGRAM=true"],
    "secondary": ["INCLUDE_SECONDARY_ALIGNMENTS=true"],
    "supplemental": ["INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true"],
    "both": ["INCLUDE_SECONDARY_ALIGNMENTS=true", "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true"],
    "alternate": ["ALTERNATE_QUALITY_ATTRIBUTE=XQ"],
    "alternate_short": ["AQA=XQ", "INCLUDE_BQ_HISTOGRAM=true"],
    "alternate_null": ["ALTERNATE_QUALITY_ATTRIBUTE=null"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_50": ["STOP_AFTER=50"],
    "stop_after_400": ["STOP_AFTER=400", "INCLUDE_BQ_HISTOGRAM=true"],
    "assume_sorted_false": ["ASSUME_SORTED=false"],
    "everything": ["INCLUDE_BQ_HISTOGRAM=true", "INCLUDE_SECONDARY_ALIGNMENTS=true",
                   "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true", "AQA=XQ"],
    "stringency_silent": ["VALIDATION_STRINGENCY=SILENT"],
}


@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, option):
    args = OPTIONS[option]
    histogram = any("HISTOGRAM=true" in a for a in args)
    compare(built["snvq_main"], args, f"snvq_main {option}", histogram)


@pytest.mark.parametrize("args", [[], ["INCLUDE_BQ_HISTOGRAM=true", "AQA=XQ"]],
                         ids=["defaults", "histogram_alternate"])
def test_golden_fixture(args):
    """tests/data/snvq.bam, the fixture of the golden records."""
    compare(FIXTURES / "snvq.bam", args, f"fixture {' '.join(args)}", histogram=bool(args))


def test_histogram_of_an_empty_file_is_left_out(built):
    compare(built["snvq_empty"], ["INCLUDE_BQ_HISTOGRAM=true"], "snvq_empty histogram")


def test_bare_secondary_passes_by_default_and_fails_when_counted(built, tmp_path):
    compare(built["snvq_secondary_bare"], ["INCLUDE_BQ_HISTOGRAM=true"], "bare", True)
    both_fail(built["snvq_secondary_bare"], ["INCLUDE_SECONDARY_ALIGNMENTS=true"], tmp_path,
              "no qa tag")


def test_short_unread_tag_with_histogram(built):
    compare(built["snvq_short_unread"], ["INCLUDE_BQ_HISTOGRAM=true"], "short unread", True)


# --- CRAM and SAM text ----------------------------------------------------------------

def test_cram_input(built, tmp_path):
    cram = tmp_path / "snvq_main.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(built["ref"]), "-o", str(cram),
                    str(built["snvq_main"])], check=True)
    compare(cram, [f"R={built['ref']}"], "snvq_main as CRAM")
    compare(cram, [f"R={built['ref']}", "INCLUDE_BQ_HISTOGRAM=true"], "CRAM histogram", True)


def lower_case_sam(bam, out):
    """`bam` as SAM text with every other base of each read in lower case."""
    text = pysam.view("-h", str(bam))
    lines = []
    for line in text.splitlines():
        if not line.startswith("@"):
            fields = line.split("\t")
            if fields[9] != "*":
                fields[9] = "".join(b.lower() if i % 2 else b for i, b in enumerate(fields[9]))
            line = "\t".join(fields)
        lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


def test_sam_text_input_with_lower_case_bases(built, tmp_path):
    sam = lower_case_sam(built["snvq_main"], tmp_path / "lower.sam")
    compare(sam, [], "lower-case SAM text")
    compare(sam, ["INCLUDE_BQ_HISTOGRAM=true", "VALIDATION_STRINGENCY=LENIENT"],
            "lower-case SAM text histogram", True)


# --- error parity -----------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"O={picard_out}", *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


FAILURES = {
    "missing_tag": ("snvq_missing_tag", [], "no qt tag"),
    "short_tag_read": ("snvq_short_read_tag", [], "of the qc tag"),
    "space_in_tag": ("snvq_space", [], "read index"),
    "int_tag": ("snvq_int_tag", [], "Value for tag qa is not a String"),
    "no_qual": ("snvq_no_qual", [], "quality string length does not match"),
    "missing_alternate": ("snvq_missing_alt", ["AQA=XQ"], "no XQ tag"),
    "bad_alternate": ("snvq_bad_alt", ["ALTERNATE_QUALITY_ATTRIBUTE=XQ"], "Invalid fastq"),
    "alternate_name": ("snvq_main", ["AQA=XQZ"], "length() == 2"),
    "alternate_absent": ("snvq_main", ["AQA=ZZ"], "no ZZ tag"),
    "level": ("snvq_main", ["METRIC_ACCUMULATION_LEVEL=SAMPLE"], "METRIC_ACCUMULATION_LEVEL"),
    "flow_mode": ("snvq_main", ["FLOW_MODE=true"], "FLOW_MODE"),
    "missing_input": ("snvq_main", ["I=/nonexistent/x.bam"], None),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    name, args, message = FAILURES[case]
    both_fail(built[name], args, tmp_path, message)


@pytest.mark.parametrize("name", ["t1", "t2", "foundation"])
def test_corpus_without_snvq_tags_fails_in_both(corpus, tmp_path, name):
    bam = corpus / "corpus.bam" if name == "foundation" else CORPUS / f"{name}.bam"
    both_fail(bam, [], tmp_path, "no qa tag")


def test_alternate_on_corpus_fails_in_both(corpus, tmp_path):
    both_fail(CORPUS / "t1.bam", ["AQA=OQ"], tmp_path)


# --- every column exercised ------------------------------------------------------------

DEFAULTS = {"0", "", "?"}


def test_every_column_takes_three_non_default_values():
    if len(SEEN) < 20:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    columns = []
    for parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns:
            if column not in columns:
                columns.append(column)
        for key, row in parsed.histogram.items():
            if parsed.histogram_columns:
                values[parsed.histogram_columns[0]].add(str(key))
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.histogram_columns:
            if column not in columns:
                columns.append(column)
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

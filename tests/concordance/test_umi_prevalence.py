"""CollectUmiPrevalenceMetrics against Picard 3.5.0: the numUmis by
duplicateSets histogram on the snvq_umi_bams UMI inputs (duplicate sets of 1
to 50 pairs, UMIs with errors, BQ above and below the threshold, missing and
empty BQ, low MAPQ, fragments and mate-unmapped reads, two libraries), the
same records unsorted and name sorted (Picard sorts them itself), the
MarkDuplicates fixture and corpus, the tier and foundation corpora (no RX
tag: an empty histogram), every option, CRAM and SAM text input, and error
parity: a pair without MC, an RX or BQ tag of another type, a BQ character
past `~`, and the option values Picard fails on.

The barcode quality filter is Picard's inverted one (KNOWN_DIFFERENCES.md):
a record with BQ is kept only when a quality is below MINIMUM_BARCODE_BQ.
"""
import subprocess
from collections import defaultdict

import pysam
import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_metrics
from conftest import RESULTS
from markdup_bams import build_all as build_markdup
from run_picard import picard_env
from snvq_umi_bams import build_all

TOOL = "CollectUmiPrevalenceMetrics"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build_all(tmp_path_factory.mktemp("umi"))


@pytest.fixture(scope="module")
def markdup(tmp_path_factory):
    return build_markdup(tmp_path_factory.mktemp("umi_markdup"))


def compare(bam, args, label, expect_sets=True):
    cells = compare_metrics(TOOL, bam, args, workdir=bam.parent, expect_file=True, seen=SEEN)
    if expect_sets:
        assert cells > 0, f"{TOOL} {label}: no set counted"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))
    return cells


INPUTS = ["umi_main", "umi_unsorted", "umi_queryname", "umi_one_no_mc"]


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, name):
    compare(built[name], [], name)


@pytest.mark.parametrize("args", [[], ["FUR=false", "BQ=100"], ["MQ=10"]],
                         ids=["defaults", "unpaired_kept", "mq_10"])
def test_golden_fixture(args):
    """tests/data/umi_prevalence.bam, the fixture of the golden records."""
    compare(FIXTURES / "umi_prevalence.bam", args, f"fixture {' '.join(args)}")


def test_empty_file(built):
    compare(built["umi_empty"], [], "umi_empty", expect_sets=False)


OPTIONS = {
    "mq_0": ["MINIMUM_MQ=0"],
    "mq_29": ["MQ=29"],
    "mq_31": ["MINIMUM_MQ=31"],
    "mq_61": ["MINIMUM_MQ=61"],
    "mq_negative": ["MINIMUM_MQ=-1"],
    "bq_0": ["MINIMUM_BARCODE_BQ=0"],
    "bq_2": ["BQ=2"],
    "bq_31": ["MINIMUM_BARCODE_BQ=31"],
    "bq_41": ["MINIMUM_BARCODE_BQ=41"],
    "bq_100": ["MINIMUM_BARCODE_BQ=100"],
    "bq_negative": ["MINIMUM_BARCODE_BQ=-1"],
    "unpaired_kept": ["FILTER_UNPAIRED_READS=false"],
    "unpaired_kept_short": ["FUR=false", "MQ=0"],
    "barcode_tag_rg": ["BARCODE_TAG=RG"],
    "barcode_tag_mc": ["BARCODE_TAG=MC", "FILTER_UNPAIRED_READS=false"],
    "barcode_tag_absent": ["BARCODE_TAG=ZZ"],
    "bq_tag_absent": ["BARCODE_BQ=ZZ"],
    "bq_tag_rx": ["BARCODE_BQ=RX", "MINIMUM_BARCODE_BQ=33"],
    "progress_1": ["PROGRESS_STEP_INTERVAL=1"],
    "progress_negative": ["PROGRESS_STEP_INTERVAL=-5"],
    "records_in_ram_1": ["MAX_RECORDS_IN_RAM=1"],
    "records_in_ram_7": ["MAX_RECORDS_IN_RAM=7"],
    "records_in_ram_100": ["MAX_RECORDS_IN_RAM=100", "FUR=false"],
    "stringency_silent": ["VALIDATION_STRINGENCY=SILENT"],
    "everything": ["MQ=10", "BQ=35", "FUR=false", "MAX_RECORDS_IN_RAM=13"],
}


@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, option):
    compare(built["umi_main"], OPTIONS[option], f"umi_main {option}",
            expect_sets=option not in ("barcode_tag_absent", "mq_61"))


@pytest.mark.parametrize("option", ["unpaired_kept", "records_in_ram_7", "bq_41"])
def test_unsorted_options(built, option):
    compare(built["umi_unsorted"], OPTIONS[option], f"umi_unsorted {option}")


def test_tmp_dir(built, tmp_path):
    compare(built["umi_main"], [f"TMP_DIR={tmp_path}", "MAX_RECORDS_IN_RAM=5"], "umi_main tmp")
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("args", [[], ["FUR=false", "MQ=0", "BQ=100"]],
                         ids=["defaults", "everything_kept"])
def test_markdup_corpus(markdup, args):
    for name in ("markdup", "markdup_duplex", "markdup_second"):
        compare(markdup[name], args, f"{name} {' '.join(args)}", expect_sets=False)


@pytest.mark.parametrize("args", [[], ["FUR=false", "MQ=0", "BQ=100"]],
                         ids=["defaults", "everything_kept"])
def test_markdup_fixture(args):
    compare(FIXTURES / "markdup.bam", args, f"markdup fixture {' '.join(args)}",
            expect_sets=False)


@pytest.mark.parametrize("name", ["t1", "t2", "t3", "foundation"])
def test_corpus_without_umis(corpus, name):
    bam = corpus / "corpus.bam" if name == "foundation" else CORPUS / f"{name}.bam"
    # Picard's STRICT reader refuses the foundation corpus; bammetrics reads a BAM unvalidated.
    silent = ["VALIDATION_STRINGENCY=SILENT"] if name == "foundation" else []
    compare(bam, silent, name, expect_sets=False)
    compare(bam, ["BARCODE_TAG=RG", "FUR=false", "MQ=0", *silent], f"{name} RG as the UMI",
            expect_sets=name == "foundation")


# --- CRAM and SAM text ------------------------------------------------------------------

def test_cram_input(built, tmp_path):
    cram = tmp_path / "umi_main.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(built["ref"]), "-o", str(cram),
                    str(built["umi_main"])], check=True)
    compare(cram, [f"R={built['ref']}"], "umi_main as CRAM")
    compare(cram, [f"R={built['ref']}", "FUR=false"], "umi_main as CRAM, unpaired kept")


def test_sam_text_input(built, tmp_path):
    sam = tmp_path / "umi_main.sam"
    sam.write_text(pysam.view("-h", str(built["umi_unsorted"])))
    compare(sam, [], "umi_unsorted as SAM text")
    compare(sam, ["FUR=false", "VALIDATION_STRINGENCY=LENIENT"], "SAM text, unpaired kept")


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
    "no_mc": ("umi_no_mc", [], "Mate CIGAR (Tag MC) not found"),
    "int_rx": ("umi_int_rx", [], "Value for tag RX is not a String"),
    "int_bq": ("umi_int_bq", [], "Value for tag BQ is not a String"),
    "bad_bq": ("umi_bad_bq", [], "Invalid fastq character"),
    "mq_null": ("umi_main", ["MINIMUM_MQ=null"], "MINIMUM_MQ"),
    "bq_null": ("umi_main", ["MINIMUM_BARCODE_BQ=null"], "MINIMUM_BARCODE_BQ"),
    "barcode_tag_long": ("umi_main", ["BARCODE_TAG=RXX"], "length() == 2"),
    "barcode_tag_null": ("umi_main", ["BARCODE_TAG=null"], "null"),
    "barcode_bq_long": ("umi_main", ["BARCODE_BQ=B"], "length() == 2"),
    "progress_zero": ("umi_main", ["PROGRESS_STEP_INTERVAL=0"], "PROGRESS_STEP_INTERVAL"),
    "records_in_ram_zero": ("umi_main", ["MAX_RECORDS_IN_RAM=0"], "MAX_RECORDS_IN_RAM"),
    "assume_sorted": ("umi_main", ["ASSUME_SORTED=true"], "ASSUME_SORTED"),
    "stop_after": ("umi_main", ["STOP_AFTER=5"], "STOP_AFTER"),
    "mq_beyond_int": ("umi_main", ["MINIMUM_MQ=2147483648"], "MINIMUM_MQ"),
    "missing_input": ("umi_main", ["I=/nonexistent/x.bam"], None),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    name, args, message = FAILURES[case]
    both_fail(built[name], args, tmp_path, message)


@pytest.mark.parametrize("args", [["MINIMUM_BARCODE_BQ=null"], ["BARCODE_TAG=RXX"],
                                  ["PROGRESS_STEP_INTERVAL=0"]],
                         ids=["bq_null", "tag_long", "progress_zero"])
def test_values_picard_never_reads_pass_on_an_empty_file(built, args):
    compare(built["umi_empty"], args, f"umi_empty {' '.join(args)}", expect_sets=False)


def test_cram_without_reference_fails_in_both(built, tmp_path):
    cram = tmp_path / "umi_main.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(built["ref"]), "-o", str(cram),
                    str(built["umi_main"])], check=True)
    both_fail(cram, [], tmp_path)


# --- every column exercised -------------------------------------------------------------

def test_every_column_takes_three_non_default_values():
    if len(SEEN) < 20:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    for parsed in SEEN:
        for key, row in parsed.histogram.items():
            values["numUmis"].add(str(key))
            for column, value in row.items():
                if value not in {"0", ""}:
                    values[column].add(value)
    assert set(values) == {"numUmis", "duplicateSets"}
    short = {c: sorted(v) for c, v in values.items() if len(v) < 3}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    RESULTS.append((f"{TOOL} distinct values", 2, 2,
                    f"numUmis {len(values['numUmis'])}, duplicateSets "
                    f"{len(values['duplicateSets'])}"))

"""EstimateLibraryComplexity against Picard 3.5.0: the DuplicationMetrics
rows per library and the duplication_group_count histogram, on aligned and
unaligned input, one or several inputs, every option that changes a number,
and error parity. Every case runs Picard with VALIDATION_STRINGENCY=SILENT.

Inputs: the complexity_bams BAMs (Illumina read names, duplicate sets with
optical copies, barcodes, two read-group sets), the read_level_bams and
asm_bams BAMs, the foundation corpus and the bwa tiers (plain read names:
every pair goes to library Unknown).

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared, except the columns in
FEWER_THAN_THREE with the reason.
"""
import subprocess
from collections import defaultdict

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_metrics
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

TOOL = "EstimateLibraryComplexity"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    bams = {}
    for prefix, build in (("asm", build_asm), ("rl", build_read_level), ("elc", build_complexity)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = bam
    return bams


def bam_path(name, built, corpus):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "fixture":
        return FIXTURES / "insert_size.bam"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam"


# Inputs where Picard writes a file with neither row nor bin (no pair passes;
# t3 is single-end), so the byte comparison of the file body is the check.
NO_CELLS = {"rl_empty", "elc_elc_empty", "asm_empty", "asm_unpaired_only", "t0_asm_empty",
            "asm_single_unmapped", "asm_all_qcfail", "rl_all_qcfail", "t3"}


def compare(bam, args, label, extra_inputs=()):
    cells = compare_metrics(TOOL, bam, [SILENT, *args, *[f"I={i}" for i in extra_inputs]],
                            workdir=bam.parent, expect_file=True, seen=SEEN)
    assert cells > 0 or label.split()[0] in NO_CELLS, f"{TOOL} {label}: no cell compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))
    return cells


# --- every input at the defaults -------------------------------------------------

INPUTS = ["elc_illumina", "elc_unaligned", "elc_unaligned_queryname", "elc_second_file",
          "elc_plain_names", "elc_elc_empty", "rl_reads", "rl_single_library", "rl_all_duplicates",
          "rl_star_fields", "rl_queryname", "rl_empty", "rl_all_qcfail", "asm_mix", "asm_unpaired_only",
          "asm_single_unmapped", "asm_all_qcfail", "asm_empty", "foundation", "fixture", "t0_asm",
          "t0_asm_empty", "t1", "t2", "t3"]


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, corpus, name):
    compare(bam_path(name, built, corpus), [], name)


def test_two_inputs(built):
    compare(built["elc_illumina"], [], "elc_illumina with elc_second_file", [built["elc_second_file"]])


def test_same_input_twice(built):
    compare(built["elc_illumina"], [], "elc_illumina twice", [built["elc_illumina"]])


def test_three_inputs(built, corpus):
    compare(built["elc_unaligned"], [], "elc_unaligned with two more",
            [built["elc_second_file"], corpus / "corpus.bam"])


# --- every option that changes a number --------------------------------------------

OPTIONS = {
    "identical_3": ["MIN_IDENTICAL_BASES=3"],
    "identical_8": ["MIN_IDENTICAL_BASES=8"],
    "identical_12": ["MIN_IDENTICAL_BASES=12"],
    "diff_rate_0": ["MAX_DIFF_RATE=0"],
    "diff_rate_0.1": ["MAX_DIFF_RATE=0.1"],
    "diff_rate_0.5": ["MAX_DIFF_RATE=0.5"],
    "mean_quality_0": ["MIN_MEAN_QUALITY=0"],
    "mean_quality_30": ["MIN_MEAN_QUALITY=30"],
    "mean_quality_25": ["MIN_MEAN_QUALITY=25"],
    "group_ratio_1": ["MAX_GROUP_RATIO=1"],
    "group_ratio_2": ["MAX_GROUP_RATIO=2"],
    "read_length_20": ["MAX_READ_LENGTH=20"],
    "read_length_50": ["MAX_READ_LENGTH=50"],
    "read_length_6": ["MAX_READ_LENGTH=6"],
    "group_count_1": ["MIN_GROUP_COUNT=1"],
    "group_count_3": ["MIN_GROUP_COUNT=3"],
    "group_count_50": ["MIN_GROUP_COUNT=50"],
    "pixels_10": ["OPTICAL_DUPLICATE_PIXEL_DISTANCE=10"],
    "pixels_2500": ["OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"],
    "pixels_0": ["OPTICAL_DUPLICATE_PIXEL_DISTANCE=0"],
    "regex_null": ["READ_NAME_REGEX=null"],
    "regex_custom": ["READ_NAME_REGEX=[^:]+:[^:]+:[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+)"],
    "regex_five": ["READ_NAME_REGEX=[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+)"],
    "regex_default_text": ["READ_NAME_REGEX=<optimized capture of last three ':' separated fields "
                           "as numeric values>"],
    "set_size_2": ["MAX_OPTICAL_DUPLICATE_SET_SIZE=2"],
    "set_size_minus_1": ["MAX_OPTICAL_DUPLICATE_SET_SIZE=-1"],
    "barcode": ["BARCODE_TAG=BC"],
    "read_one_barcode": ["READ_ONE_BARCODE_TAG=BX"],
    "read_two_barcode": ["READ_TWO_BARCODE_TAG=BX"],
    "all_barcodes": ["BARCODE_TAG=BC", "READ_ONE_BARCODE_TAG=BX", "READ_TWO_BARCODE_TAG=BX"],
    "barcode_absent_tag": ["BARCODE_TAG=ZZ"],
    "everything": ["MIN_IDENTICAL_BASES=4", "MAX_DIFF_RATE=0.08", "MIN_MEAN_QUALITY=15",
                   "MAX_GROUP_RATIO=3", "MAX_READ_LENGTH=40", "MIN_GROUP_COUNT=1",
                   "OPTICAL_DUPLICATE_PIXEL_DISTANCE=300"],
    "reference_ignored": [f"R={CORPUS / 'reference.fa'}"],
}
OPTION_INPUTS = ["elc_illumina", "elc_unaligned", "elc_plain_names", "rl_reads", "foundation"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, corpus, option, name):
    compare(bam_path(name, built, corpus), OPTIONS[option], f"{name} {option}")


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"O={picard_out}", SILENT, *args, "QUIET=true",
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
    "assume_sorted": (["ASSUME_SORTED=true"], "ASSUME_SORTED"),
    "stop_after": (["STOP_AFTER=5"], "STOP_AFTER"),
    "read_length_below_seed": (["MAX_READ_LENGTH=3"], "MAX_READ_LENGTH"),
    "seed_zero": (["MIN_IDENTICAL_BASES=0"], "MIN_IDENTICAL_BASES"),
    "bad_regex": (["READ_NAME_REGEX=("], "READ_NAME_REGEX"),
    "missing_input": (["I=/nonexistent/x.bam"], None),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    args, message = FAILURES[case]
    both_fail(built["elc_illumina"], args, tmp_path, message)


def test_barcode_tag_that_is_not_a_string_fails(built, tmp_path):
    """XN is an integer tag on some asm_mix reads: Picard casts it to String."""
    both_fail(built["asm_mix"], ["BARCODE_TAG=XN", "MIN_MEAN_QUALITY=0"], tmp_path)


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {
    "UNPAIRED_READS_EXAMINED": "always 0: EstimateLibraryComplexity counts pairs only",
    "SECONDARY_OR_SUPPLEMENTARY_RDS": "always 0: secondary and supplementary records are skipped",
    "UNMAPPED_READS": "always 0: alignment is not read",
    "UNPAIRED_READ_DUPLICATES": "always 0: pairs only",
}


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
        for key, row in parsed.histogram.items():
            values[parsed.histogram_columns[0]].add(key)
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns + parsed.histogram_columns:
            if column not in columns:
                columns.append(column)
    short = {c: sorted(values[c]) for c in columns
             if len(values[c]) < 3 and c not in FEWER_THAN_THREE}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

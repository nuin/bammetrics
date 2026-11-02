"""Read-level tools against Picard 3.5.0: QualityScoreDistribution,
MeanQualityByCycle, CollectBaseDistributionByCycle, CollectQualityYieldMetrics,
CollectDuplicateMetrics and AccumulateQualityYieldMetrics.

Each case runs Picard and `bammetrics <Tool>` with the same KEY=VALUE
arguments; compare_metrics checks every cell and the whole file after the
StringHeader lines byte for byte. Inputs: the bwa tiers t0 to t3
(tests/corpus/out), the foundation corpus, the insert-size fixture and the
hand-built BAMs of read_level_bams.py. Every case runs Picard with
VALIDATION_STRINGENCY=SILENT, because several BAMs name read groups missing
from the header or carry flag combinations STRICT validation rejects.

Charts: QualityScoreDistribution, MeanQualityByCycle and
CollectBaseDistributionByCycle draw a PDF with R when CHART is given.
run_picard.picard_env puts a no-op Rscript first on PATH, so Picard's R call
succeeds and writes no PDF; bammetrics accepts CHART and writes no PDF either.
The CHART cases check the text file is unchanged.

The last tests check that every column took at least 3 distinct non-default
values over all Picard files compared, except the columns in FEWER_THAN_THREE
with the reason.

Needs picard, samtools and pysam; run `make concordance`.
"""
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_metrics, run_picard_tool
from conftest import RESULTS
from read_level_bams import build_all
from run_picard import picard_env

SILENT = "VALIDATION_STRINGENCY=SILENT"
QSD, MQBC, CBDBC = "QualityScoreDistribution", "MeanQualityByCycle", "CollectBaseDistributionByCycle"
QY, DUP, ACC = "CollectQualityYieldMetrics", "CollectDuplicateMetrics", "AccumulateQualityYieldMetrics"
BY_CYCLE = [QSD, MQBC, CBDBC]
OUTPUT_KEY = {DUP: "M"}
SEEN = defaultdict(list)

TIERS = ["t0_qy", "t0_asm", "t0_asm_cycles", "t0_asm_empty", "t0_filters", "t0_no_pairs",
         "t0_all_tandem", "t1", "t2", "t3"]
BUILT = ["reads", "reads_no_oq", "q0_cycle", "second_only", "all_qcfail", "star_fields",
         "single_library", "no_duplicates", "all_duplicates", "queryname", "empty",
         "no_read_groups"]


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build_all(tmp_path_factory.mktemp("read_level"))


def bam_path(name, built, corpus):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "fixture":
        return FIXTURES / "insert_size.bam"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam"


# Cases where Picard's file has no metrics row and no histogram bin, so the
# byte comparison of the file body is the whole check.
NO_CELLS = {
    "empty", "t0_asm_empty", "all_qcfail pf_only", "all_qcfail aligned_pf",
    "all_qcfail no_calls_pf_short", "all_qcfail everything",
    (MQBC, "reads stop_after_1"),  # the first record has quality 0 at every base
}


def compare(tool, bam, args, label):
    """AccumulateQualityYieldMetrics takes no VALIDATION_STRINGENCY."""
    extra = [] if tool == ACC else [SILENT]
    cells = compare_metrics(tool, bam, [*extra, *args], workdir=bam.parent, expect_file=True,
                            seen=SEEN[tool], output_key=OUTPUT_KEY.get(tool, "O"))
    assert cells > 0 or label in NO_CELLS or (tool, label) in NO_CELLS, (
        f"{tool} {label}: no cell compared")
    RESULTS.append((f"{tool} {label}", cells, cells, "every cell and the file body identical"))
    return cells


# --- every input at the defaults -------------------------------------------------

INPUTS = TIERS + ["foundation", "fixture"] + BUILT
# t0_qy carries OQ tags shorter than the read: Picard throws an index error.
FAILS_IN_PICARD = {(QSD, "t0_qy"): "OQ is shorter"}


@pytest.mark.parametrize("tool", [*BY_CYCLE, QY, DUP])
@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, corpus, tmp_path, tool, name):
    if (tool, name) in FAILS_IN_PICARD:
        both_fail(tool, bam_path(name, built, corpus), [], tmp_path, FAILS_IN_PICARD[tool, name])
        return
    compare(tool, bam_path(name, built, corpus), [], name)


# --- every option that changes a number --------------------------------------------

FILTERS = {
    "aligned_only": ["ALIGNED_READS_ONLY=true"],
    "pf_only": ["PF_READS_ONLY=true"],
    "aligned_pf": ["ALIGNED_READS_ONLY=true", "PF_READS_ONLY=true"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_137": ["STOP_AFTER=137"],
    "assume_sorted_false": ["ASSUME_SORTED=false"],
}
QSD_ONLY = {
    "no_calls": ["INCLUDE_NO_CALLS=true"],
    "no_calls_pf_short": ["INCLUDE_NO_CALLS=true", "PF=true"],
    "everything": ["INCLUDE_NO_CALLS=true", "ALIGNED_READS_ONLY=true", "PF_READS_ONLY=true"],
}
OPTION_INPUTS = ["reads", "reads_no_oq", "t2", "foundation", "all_qcfail", "star_fields"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", FILTERS)
@pytest.mark.parametrize("tool", BY_CYCLE)
def test_filter_options(built, corpus, tool, option, name):
    compare(tool, bam_path(name, built, corpus), FILTERS[option], f"{name} {option}")


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", QSD_ONLY)
def test_include_no_calls(built, corpus, option, name):
    compare(QSD, bam_path(name, built, corpus), QSD_ONLY[option], f"{name} {option}")


@pytest.mark.parametrize("tool", BY_CYCLE)
def test_chart_is_accepted_and_changes_nothing(built, tmp_path, tool):
    compare(tool, built["reads"], [f"CHART={tmp_path / 'chart.pdf'}"], "reads CHART")


QY_OPTIONS = {
    "no_oq": ["USE_ORIGINAL_QUALITIES=false"],
    "oq_short_name": ["OQ=false"],
    "secondary": ["INCLUDE_SECONDARY_ALIGNMENTS=true"],
    "supplemental": ["INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true"],
    "both": ["INCLUDE_SECONDARY_ALIGNMENTS=true", "INCLUDE_SUPPLEMENTAL_ALIGNMENTS=true"],
    "flow_mode_false": ["FLOW_MODE=false"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_137": ["STOP_AFTER=137"],
    "assume_sorted_false": ["ASSUME_SORTED=false"],
}


@pytest.mark.parametrize("name", ["reads", "reads_no_oq", "star_fields", "t0_qy", "t2", "foundation"])
@pytest.mark.parametrize("option", QY_OPTIONS)
def test_quality_yield_options(built, corpus, option, name):
    compare(QY, bam_path(name, built, corpus), QY_OPTIONS[option], f"{name} {option}")


DUP_OPTIONS = {
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_137": ["STOP_AFTER=137"],
    "stop_after_400": ["STOP_AFTER=400"],
    "assume_sorted_false": ["ASSUME_SORTED=false"],
    "metrics_file_long_name": [],
}


@pytest.mark.parametrize("name", ["reads", "single_library", "foundation", "t2"])
@pytest.mark.parametrize("option", DUP_OPTIONS)
def test_duplicate_metrics_options(built, corpus, option, name):
    compare(DUP, bam_path(name, built, corpus), DUP_OPTIONS[option], f"{name} {option}")


# --- AccumulateQualityYieldMetrics -----------------------------------------------------

SHARDS = {
    "two_tiers": [("t1", []), ("t2", [])],
    "same_file_twice": [("reads", []), ("reads", [])],
    "options_mixed": [("reads", ["OQ=false"]), ("star_fields", []), ("t0_qy", ["INCLUDE_SECONDARY_ALIGNMENTS=true"])],
    "one_file": [("t3", [])],
    "with_empty": [("empty", []), ("t0_asm_empty", []), ("reads_no_oq", [])],
    "only_empty": [("empty", [])],
}


@pytest.mark.parametrize("shards", SHARDS)
def test_accumulate_quality_yield(built, corpus, tmp_path, shards):
    """Picard QualityYieldMetrics files of several inputs, merged by both tools."""
    files = []
    for i, (name, args) in enumerate(SHARDS[shards]):
        out = tmp_path / f"shard{i}.txt"
        run_picard_tool(QY, bam_path(name, built, corpus), out, [SILENT, *args])
        files.append(out)
    compare(ACC, files[0], [f"I={f}" for f in files[1:]], f"shards {shards}")


# --- error parity ---------------------------------------------------------------------

def both_fail(tool, bam, args, tmp_path, message=None):
    key = OUTPUT_KEY.get(tool, "O")
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        run_picard_tool(tool, bam, picard_out, [SILENT, *args], key)
    ours = subprocess.run([str(BAMMETRICS), tool, f"I={bam}", f"{key}={ours_out}", *args],
                          capture_output=True, text=True)
    assert ours.returncode != 0, ours.stdout
    if message:
        assert message in ours.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{tool} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


@pytest.mark.parametrize("tool", [*BY_CYCLE, QY, DUP])
def test_assume_sorted_false_on_query_name_order_fails(built, tmp_path, tool):
    both_fail(tool, built["queryname"], ["ASSUME_SORTED=false"], tmp_path, "coordinate sorted")


def test_flow_mode_fails(built, tmp_path):
    both_fail(QY, built["reads"], ["FLOW_MODE=true"], tmp_path, "FLOW_MODE")


def test_oq_shorter_than_qual_fails_quality_score_distribution(built, tmp_path):
    both_fail(QSD, built["oq_short"], [], tmp_path, "OQ is shorter")


@pytest.mark.parametrize("tool", [MQBC, CBDBC, QY, DUP])
def test_oq_shorter_than_qual_elsewhere(built, tool):
    compare(tool, built["oq_short"], [], "oq_short")


def test_duplicate_metrics_takes_no_output_option(built, tmp_path):
    ours = subprocess.run([str(BAMMETRICS), DUP, f"I={built['reads']}", f"O={tmp_path / 'x.txt'}"],
                          capture_output=True, text=True)
    assert ours.returncode != 0 and "option O is not implemented" in ours.stderr


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {
    (CBDBC, "READ_END"): "1 or 2 by definition",
    (DUP, "READ_PAIR_OPTICAL_DUPLICATES"): "always 0: CollectDuplicateMetrics detects no optical duplicates",
}


@pytest.mark.parametrize("tool", [*BY_CYCLE, QY, DUP, ACC])
def test_every_column_takes_three_non_default_values(tool):
    seen = SEEN[tool]
    if len(seen) < 5:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    columns = []
    for parsed in seen:
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
             if len(values[c]) < 3 and (tool, c) not in FEWER_THAN_THREE}
    assert not short, f"{tool}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{tool} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))


def test_picard_is_on_path():
    out = subprocess.run(["picard", QSD, "--version"], capture_output=True, text=True,
                         env=picard_env())
    assert "3.5.0" in out.stdout + out.stderr

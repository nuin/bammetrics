"""CollectJumpingLibraryMetrics against Picard 3.5.0: the one JumpingLibraryMetrics
row on the jumping_bams corpus, the foundation corpus, corpus t1 and the two
fixture BAMs, every option that changes a number, two and three INPUT files,
CRAM input, and error parity.

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared.
"""
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_metrics
from conftest import RESULTS
from jumping_bams import build_all
from run_picard import picard_env

TOOL = "CollectJumpingLibraryMetrics"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build_all(tmp_path_factory.mktemp("jump"))


def bam_path(name, built, corpus):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name in ("insert_size", "orientation"):
        return FIXTURES / f"{name}.bam"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam"


def compare(bam, args, label, extra_inputs=()):
    cells = compare_metrics(TOOL, bam, [*args, *[f"I={i}" for i in extra_inputs]],
                            workdir=bam.parent, expect_file=True, seen=SEEN)
    assert cells == 17, f"{TOOL} {label}: {cells} cells compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))


# --- every input at the defaults -------------------------------------------------

INPUTS = ["main", "second", "dupes", "no_jumps", "big", "empty", "unmapped_only", "foundation",
          "t1", "t2", "insert_size", "orientation"]


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, corpus, name):
    compare(bam_path(name, built, corpus), [], name)


# --- every option that changes a number --------------------------------------------

OPTIONS = {
    "mq_1": ["MINIMUM_MAPPING_QUALITY=1"],
    "mq_20": ["MQ=20"],
    "mq_60": ["MINIMUM_MAPPING_QUALITY=60"],
    "mq_61": ["MINIMUM_MAPPING_QUALITY=61"],
    "tail_1": ["TAIL_LIMIT=1"],
    "tail_2": ["T=2"],
    "tail_3": ["TAIL_LIMIT=3"],
    "tail_10": ["TAIL_LIMIT=10"],
    "tail_0": ["TAIL_LIMIT=0"],
    "tail_negative": ["TAIL_LIMIT=-1"],
    "chimera_0": ["CHIMERA_KB_MIN=0"],
    "chimera_negative": ["CHIMERA_KB_MIN=-5"],
    "chimera_1000": ["CHIMERA_KB_MIN=1000"],
    "chimera_2999": ["CHIMERA_KB_MIN=2999"],
    "chimera_3000": ["CHIMERA_KB_MIN=3000"],
    "chimera_5000": ["CHIMERA_KB_MIN=5000"],
    "chimera_99999": ["CHIMERA_KB_MIN=99999"],
    "chimera_100001": ["CHIMERA_KB_MIN=100001"],
    "everything": ["MQ=20", "T=2", "CHIMERA_KB_MIN=5000"],
    "reference_ignored": [f"R={CORPUS / 'reference.fa'}"],
    "stringency_silent": ["VALIDATION_STRINGENCY=SILENT"],
}
OPTION_INPUTS = ["main", "second", "dupes", "no_jumps", "foundation"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, corpus, option, name):
    compare(bam_path(name, built, corpus), OPTIONS[option], f"{name} {option}")


@pytest.mark.parametrize("option", ["chimera_1000", "chimera_2000", "chimera_2500", "mq_20"])
def test_big_mode_sample(built, option):
    """The mode comes from the first 50000 RF pairs: 2000 here, 2600 over every pair."""
    args = {"chimera_2000": ["CHIMERA_KB_MIN=2000"], "chimera_2500": ["CHIMERA_KB_MIN=2500"]}
    compare(built["big"], args.get(option, OPTIONS[option]), f"big {option}")


# --- several inputs ---------------------------------------------------------------

def test_two_inputs(built):
    compare(built["main"], [], "main with second", [built["second"]])


def test_two_inputs_with_options(built):
    compare(built["main"], ["MQ=20", "T=2", "CHIMERA_KB_MIN=1000"], "main with second, options",
            [built["second"]])


def test_big_first_of_two_samples_25000(built):
    """25000 pairs per file: the first region of `big` (mode 1800) sets the chimera size."""
    compare(built["big"], ["CHIMERA_KB_MIN=1000"], "big with main", [built["main"]])


def test_big_second_of_two(built):
    compare(built["main"], ["CHIMERA_KB_MIN=1000"], "main with big", [built["big"]])


def test_same_input_twice(built):
    compare(built["dupes"], [], "dupes twice", [built["dupes"]])


def test_three_inputs(built, corpus):
    compare(built["no_jumps"], ["CHIMERA_KB_MIN=1000"], "no_jumps, foundation, second",
            [corpus / "corpus.bam", built["second"]])


def test_empty_and_unmapped_inputs_beside_a_full_one(built):
    compare(built["empty"], [], "empty, unmapped_only, main",
            [built["unmapped_only"], built["main"]])


# --- CRAM -----------------------------------------------------------------------

def test_cram_input(built, tmp_path):
    cram = tmp_path / "main.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(built["ref"]), "-o", str(cram),
                    str(built["main"])], check=True)
    compare(cram, [f"R={built['ref']}", "CHIMERA_KB_MIN=5000"], "main as CRAM")
    compare(cram, [f"R={built['ref']}"], "main as CRAM with the BAM", [built["second"]])


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None, extra_inputs=()):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    inputs = [f"I={bam}", *[f"I={i}" for i in extra_inputs]]
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, *inputs, f"O={picard_out}", *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, *inputs, f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


FAILURES = {
    "assume_sorted": ("main", ["ASSUME_SORTED=true"], "ASSUME_SORTED"),
    "stop_after": ("main", ["STOP_AFTER=5"], "STOP_AFTER"),
    "level": ("main", ["METRIC_ACCUMULATION_LEVEL=SAMPLE"], "METRIC_ACCUMULATION_LEVEL"),
    "queryname_sorted": ("queryname", [], "sort order"),
    "unpaired_record": ("unpaired", [], "paired"),
    "mq_tag_not_an_integer": ("bad_mq", [], "MQ tag"),
    "missing_input": ("main", ["I=/nonexistent/x.bam"], None),
    "bad_tail_limit": ("main", ["TAIL_LIMIT=ten"], "TAIL_LIMIT"),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    name, args, message = FAILURES[case]
    both_fail(built[name], args, tmp_path, message)


def test_a_failing_second_input_fails_the_run(built, tmp_path):
    both_fail(built["main"], [], tmp_path, "sort order", extra_inputs=[built["queryname"]])


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", "", "?", "-0"}


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
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

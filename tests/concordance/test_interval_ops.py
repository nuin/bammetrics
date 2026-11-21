"""The interval layer crates/core/src/interval_ops against Picard 3.5.0 IntervalListTools.

The example `interval_ops` (crates/core/examples/interval_ops.rs) runs the steps of
`IntervalListTools.doWork` over the layer: each input read and padded, the inputs reduced
by the action (`concatenate`, or `intersection` for INTERSECT), the action, `sorted`,
`invert`, `unique_intervals`, `break_at_band_multiples`, the `@PG` and `@CO` lines, then
one list or a scatter directory. Each case runs Picard and the example on the same
arguments and compares every line of every file written, the `@PG` line by its `ID` and
`PN` (its `CL` is each program's own command line); a scatter is compared as the set of
its files, names and lines. Standard output (`OUTPUT_VALUE`) is compared too. Where
Picard fails, the example must fail and leave no output.

Inputs come from interval_lists.py (seeded), from Picard itself (a UNION it wrote) and
from the foundation corpus (its bait and target lists, written by Picard's
BedToIntervalList).
"""
import subprocess

import pytest

import interval_lists
from cases import ROOT
from conftest import RESULTS
from run_picard import picard_env

EXAMPLE = ROOT / "target" / "release" / "examples" / "interval_ops"
ACTIONS = ("CONCAT", "UNION", "INTERSECT", "SUBTRACT", "SYMDIFF", "OVERLAPS")
SECOND = {"SUBTRACT", "SYMDIFF", "OVERLAPS"}
MODES = ("INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
         "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW", "INTERVAL_COUNT",
         "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER")
OPTION_SETS = {
    "defaults": [],
    "pad0": ["PADDING=0"],
    "pad7": ["PADDING=7"],
    "pad_negative": ["PADDING=-5"],
    "unique": ["UNIQUE=true"],
    "unsorted": ["SORT=false"],
    "unsorted_unique": ["SORT=false", "UNIQUE=true"],
    "unique_no_abutting": ["UNIQUE=true", "DONT_MERGE_ABUTTING=true"],
    "invert": ["INVERT=true"],
    "invert_no_abutting": ["INVERT=true", "DONT_MERGE_ABUTTING=true"],
    "bands": ["BREAK_BANDS_AT_MULTIPLES_OF=100"],
    "bands_unique": ["BREAK_BANDS_AT_MULTIPLES_OF=37", "UNIQUE=true"],
    "pad7_unique_comments": ["PADDING=7", "UNIQUE=true", "COMMENT=a note", "COMMENT=@COkept"],
}
SCATTER_COUNTS = (1, 2, 3, 7, 50)
# SCATTER_CONTENT values per mode: bases for the first three, intervals for the others.
SCATTER_CONTENTS = [(mode, content) for mode in MODES[:3] for content in (100, 1000)] + \
    [(mode, content) for mode in MODES[3:] for content in (7, 30)]
# Cases that matched, by kind, for test_cases_are_covered.
SEEN = {"action": set(), "mode": set(), "error": 0}


@pytest.fixture(scope="module")
def example():
    subprocess.run(["cargo", "build", "--release", "-p", "bammetrics-core", "--example", "interval_ops"],
                   cwd=ROOT, check=True, capture_output=True)
    return EXAMPLE


@pytest.fixture(scope="module")
def lists(tmp_path_factory):
    paths = interval_lists.write(tmp_path_factory.mktemp("lists"), seed=1)
    made = tmp_path_factory.mktemp("picard_made") / "union.interval_list"
    done = subprocess.run(["picard", "IntervalListTools", f"I={paths['a']}", f"I={paths['c']}", "ACTION=UNION",
                           f"O={made}", "QUIET=true", "VERBOSITY=ERROR"],
                          env=picard_env(), capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-500:]
    paths["picard_union"] = made
    return paths


def comparable(path):
    """The lines of `path`, each `@PG` line without its `CL` field."""
    lines = []
    for line in path.read_text().split("\n"):
        if line.startswith("@PG"):
            line = "\t".join(f for f in line.split("\t") if not f.startswith("CL:"))
        lines.append(line)
    return lines


def tree(path):
    """Every file under `path` by its relative name, `None` when nothing is there."""
    if not path.exists():
        return None
    if path.is_file():
        return {".": comparable(path)}
    return {str(f.relative_to(path)): comparable(f) for f in sorted(path.rglob("*")) if f.is_file()}


def run(example, args, tmp_path):
    picard_out, ours_out = tmp_path / "picard" / "out", tmp_path / "ours" / "out"
    picard_out.parent.mkdir()
    ours_out.parent.mkdir()
    picard = subprocess.run(["picard", "IntervalListTools", *args, f"O={picard_out}", "QUIET=true",
                             "VERBOSITY=ERROR"], env=picard_env(), capture_output=True, text=True)
    ours = subprocess.run([str(example), *args, f"O={ours_out}"], capture_output=True, text=True)
    return picard, ours, picard_out, ours_out


def compare(example, args, tmp_path, label):
    """Both succeed with the same files and standard output, or both fail and the
    example leaves nothing. Returns True when both succeeded."""
    picard, ours, picard_out, ours_out = run(example, args, tmp_path)
    if picard.returncode != 0:
        assert ours.returncode != 0, f"{label}: Picard failed, the example did not:\n{picard.stderr[-800:]}"
        assert not ours_out.exists(), f"{label}: the example failed and left {ours_out}"
        SEEN["error"] += 1
        RESULTS.append((f"interval_ops {label}", 0, 0, "both exit non-zero, no output"))
        return False
    assert ours.returncode == 0, f"{label}: the example failed:\n{ours.stderr}"
    expected, got = tree(picard_out), tree(ours_out)
    assert expected == got, f"{label}:\nPicard: {expected}\nexample: {got}"
    assert picard.stdout == ours.stdout, f"{label}: standard output {picard.stdout!r} != {ours.stdout!r}"
    lines = sum(len(v) for v in (expected or {}).values())
    RESULTS.append((f"interval_ops {label}", lines, lines, f"{len(expected or {})} files identical"))
    return True


def action_args(action, lists, first=("a",), second=("b",)):
    args = [f"ACTION={action}"] + [f"I={lists[name]}" for name in first]
    if action in SECOND:
        args += [f"SI={lists[name]}" for name in second]
    return args


@pytest.mark.parametrize("option_set", list(OPTION_SETS))
@pytest.mark.parametrize("action", ACTIONS)
def test_action(example, lists, tmp_path, action, option_set):
    """Every ACTION under every option set: INPUT a, b and c for the actions that take
    one input option, INPUT a and SECOND_INPUT b and c for the others."""
    first = ("a",) if action in SECOND else ("a", "b", "c")
    args = action_args(action, lists, first, ("b", "c")) + OPTION_SETS[option_set]
    if compare(example, args, tmp_path, f"{action} {option_set}"):
        SEEN["action"].add(action)


INPUT_CASES = {
    "sorted_header_unique": (["I={sorted_header}", "SORT=false", "UNIQUE=true"]),
    "sorted_header_sorted": (["I={sorted_header}"]),
    "extra_header": (["I={extra_header}", "I={a}", "COMMENT=more"]),
    "extra_header_subtract": (["ACTION=SUBTRACT", "I={extra_header}", "SI={a}"]),
    "empty_concat": (["I={empty}"]),
    "empty_invert": (["I={empty}", "INVERT=true"]),
    "empty_subtract": (["ACTION=SUBTRACT", "I={a}", "SI={empty}"]),
    "empty_from": (["ACTION=SUBTRACT", "I={empty}", "SI={a}"]),
    "empty_intersect": (["ACTION=INTERSECT", "I={a}", "I={empty}"]),
    "empty_symdiff": (["ACTION=SYMDIFF", "I={empty}", "SI={b}"]),
    "empty_overlaps": (["ACTION=OVERLAPS", "I={a}", "SI={empty}"]),
    "bare": (["I={bare}"]),
    "bare_invert": (["I={bare}", "INVERT=true"]),
    "bare_with_a": (["I={bare}", "I={a}"]),
    "other_dictionary": (["I={a}", "I={other_dictionary}"]),
    "other_dictionary_subtract": (["ACTION=SUBTRACT", "I={a}", "SI={other_dictionary}"]),
    "picard_made_subtract": (["ACTION=SUBTRACT", "I={picard_union}", "SI={b}"]),
    "picard_made_overlaps": (["ACTION=OVERLAPS", "I={b}", "SI={picard_union}", "SORT=false"]),
    "picard_made_symdiff": (["ACTION=SYMDIFF", "I={picard_union}", "SI={a}", "UNIQUE=true"]),
    "repeat_input": (["ACTION=INTERSECT", "I={a}", "I={a}"]),
    "non_utf8": (["I={non_utf8}", "UNIQUE=true"]),
    "count_bases": (["I={a}", "UNIQUE=true", "OUTPUT_VALUE=BASES"]),
    "count_intervals": (["I={a}", "I={b}", "OUTPUT_VALUE=INTERVALS"]),
    "subtract_no_second": (["ACTION=SUBTRACT", "I={a}"]),
    "concat_with_second": (["I={a}", "SI={b}"]),
    "scatter_zero": (["I={a}", "SCATTER_COUNT=0"]),
    "bands_negative": (["I={a}", "BREAK_BANDS_AT_MULTIPLES_OF=-1"]),
    "count_output_none": (["I={a}", "COUNT_OUTPUT=/dev/null"]),
    "empty_scatter": (["I={empty}", "SCATTER_COUNT=3"]),
    "invert_scatter": (["I={a}", "INVERT=true", "SCATTER_COUNT=4"]),
    "bands_scatter": (["I={a}", "BREAK_BANDS_AT_MULTIPLES_OF=250", "SCATTER_COUNT=5",
                       "SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION"]),
    "scatter_content_zero": (["I={a}", "SCATTER_CONTENT=0"]),
}


@pytest.mark.parametrize("case", list(INPUT_CASES))
def test_inputs(example, lists, tmp_path, case):
    """Lists Picard made, header shapes, empty lists, mismatched dictionaries and the
    command-line checks."""
    args = [a.format(**lists) for a in INPUT_CASES[case]]
    compare(example, args, tmp_path, case)


@pytest.mark.parametrize("name", list(interval_lists.BAD))
def test_bad_list(example, lists, tmp_path, name):
    """Each file htsjdk's reader refuses fails both, and the example writes nothing."""
    assert not compare(example, [f"I={lists[name]}"], tmp_path, f"bad {name}")


@pytest.mark.parametrize("count", SCATTER_COUNTS)
@pytest.mark.parametrize("mode", MODES)
def test_scatter(example, lists, tmp_path, mode, count):
    """Every SUBDIVISION_MODE at SCATTER_COUNT 1, 2, 3, 7 and 50 over INPUT a and b: the
    directory of temp_NNNN_of_N/scattered.interval_list files, including the empty lists
    and the fewer lists a mode makes."""
    args = [f"I={lists['a']}", f"I={lists['b']}", f"SUBDIVISION_MODE={mode}", f"SCATTER_COUNT={count}",
            "OUTPUT_VALUE=INTERVALS"]
    if compare(example, args, tmp_path, f"scatter {mode} {count}"):
        SEEN["mode"].add(mode)


@pytest.mark.parametrize("count", (50, 200))
@pytest.mark.parametrize("mode", MODES)
def test_scatter_fewer(example, lists, tmp_path, mode, count):
    """More lists asked than list b can fill: the lists each mode makes before the
    intervals run out, and the last one with the rest."""
    args = [f"I={lists['b']}", f"SUBDIVISION_MODE={mode}", f"SCATTER_COUNT={count}"]
    compare(example, args, tmp_path, f"scatter fewer {mode} {count}")


@pytest.mark.parametrize(("mode", "content"), SCATTER_CONTENTS)
def test_scatter_content(example, lists, tmp_path, mode, content):
    """SCATTER_CONTENT in place of SCATTER_COUNT, on the unique union of a and b."""
    args = [f"I={lists['a']}", f"I={lists['b']}", "ACTION=UNION", f"SUBDIVISION_MODE={mode}",
            f"SCATTER_CONTENT={content}"]
    compare(example, args, tmp_path, f"scatter content {mode} {content}")


@pytest.mark.parametrize("action", ACTIONS)
def test_corpus_lists(example, corpus, tmp_path, action):
    """The foundation corpus's bait and target lists, as Picard's BedToIntervalList
    wrote them, under every action."""
    paths = {"a": corpus / "target.interval_list", "b": corpus / "bait.interval_list"}
    compare(example, action_args(action, paths), tmp_path, f"corpus {action}")


def test_cases_are_covered():
    """Every action and every scatter mode produced a file equal to Picard's, and the
    error cases ran."""
    assert SEEN["action"] == set(ACTIONS), SEEN
    assert SEEN["mode"] == set(MODES), SEEN
    assert SEEN["error"] >= len(interval_lists.BAD), SEEN

"""CreateVerifyIDIntensityContaminationMetricsFile against Picard 3.5.0 on hand-written
VerifyIDIntensity outputs.

Each valid input is converted by both tools and the file under the prefix compared cell by
cell and byte for byte after the header lines (compare_prefixed). Each invalid input makes
both tools exit non-zero with no file left: a line outside Picard's three patterns, or a
line the pattern lets through whose number Java cannot read. The last test is the
distinct-value gate.
"""
import random
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS
from compare import compare_prefixed
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CreateVerifyIDIntensityContaminationMetricsFile"
HEAD = "ID  %Mix     LLK          LLK0"
RULE = "-" * 38
SEEN = []


def rows(seed, n):
    rng = random.Random(seed)
    return [f"{i}  {rng.uniform(0, 0.3):.5f}  {rng.uniform(-90000, -1000):.4f}  "
            f"{rng.uniform(-90000, -1000):.3f}" for i in range(n)]


VALID = {
    "one_row": "\n".join([HEAD, RULE, "0   0.0034  -1234.5  -1240.25"]) + "\n",
    "several_rows": "\n".join([HEAD, RULE, *rows(1, 12)]) + "\n",
    "tabs_and_trailing_spaces": "\n".join(["ID\t%Mix\tLLK\tLLK0  ", "-", "7\t.5\t-.25\t5  ",
                                           "8 \t 10   12 -0"]) + "\n",
    "crlf": "\r\n".join([HEAD, RULE, *rows(2, 3)]) + "\r\n",
    "cr": "\r".join([HEAD, RULE, *rows(3, 3)]),
    "header_and_rule_only": HEAD + "\n" + RULE + "\n",
    "header_and_rule_no_line_end": HEAD + "\n" + RULE,
    "leading_zeros_and_signs": "\n".join([HEAD, RULE, "007 000.5 -0 -00", "2147483647 1 1 1",
                                          "3 0 0 -0.0"]) + "\n",
    "many_rows": "\n".join([HEAD, RULE, *rows(4, 120)]) + "\n",
}

INVALID = {
    "empty": "",
    "header_only": HEAD + "\n",
    "wrong_header_case": "ID %MIX LLK LLK0\n" + RULE + "\n",
    "header_extra_word": HEAD + " X\n" + RULE + "\n",
    "header_leading_space": " " + HEAD + "\n" + RULE + "\n",
    "rule_with_spaces": HEAD + "\n--- ---\n",
    "rule_empty": HEAD + "\n\n",
    "rule_equals": HEAD + "\n=====\n",
    "blank_line": HEAD + "\n" + RULE + "\n0 0.1 -1 -2\n\n",
    "leading_space": HEAD + "\n" + RULE + "\n 0 0.1 -1 -2\n",
    "five_fields": HEAD + "\n" + RULE + "\n0 0.1 -1 -2 9\n",
    "three_fields": HEAD + "\n" + RULE + "\n0 0.1 -1\n",
    "mix_ends_in_point": HEAD + "\n" + RULE + "\n0 1. -1 -2\n",
    "mix_negative": HEAD + "\n" + RULE + "\n0 -0.1 -1 -2\n",
    "scientific": HEAD + "\n" + RULE + "\n0 1e-5 -1 -2\n",
    "llk_dash_inside": HEAD + "\n" + RULE + "\n0 0.1 1-2 -2\n",
    "llk0_double_minus": HEAD + "\n" + RULE + "\n0 0.1 -1 --2.5\n",
    "id_negative": HEAD + "\n" + RULE + "\n-1 0.1 -1 -2\n",
    "id_overflow": HEAD + "\n" + RULE + "\n2147483648 0.1 -1 -2\n",
    "nan": HEAD + "\n" + RULE + "\n0 NaN -1 -2\n",
    "unicode_digits": HEAD + "\n" + RULE + "\n١ 0.1 -1 -2\n",
}


def write(tmp_path, name, text):
    path = tmp_path / f"{name}.txt"
    path.write_bytes(text.encode("utf-8"))
    return path


@pytest.mark.parametrize("case", sorted(VALID))
def test_valid_input(case, tmp_path):
    path = write(tmp_path, case, VALID[case])
    cells = compare_prefixed(TOOL, [f"I={path}"], workdir=tmp_path, seen=SEEN)
    RESULTS.append((f"{TOOL} {case}", cells, cells, "every file identical"))


def both_fail(args, tmp_path, message=""):
    picard = subprocess.run(["picard", TOOL, *args, f"O={tmp_path / 'picard'}", "QUIET=true"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    ours = subprocess.run([str(BAMMETRICS), TOOL, *args, f"O={tmp_path / 'bammetrics'}"],
                          capture_output=True, text=True)
    assert ours.returncode != 0, ours.stdout
    assert message in ours.stderr, ours.stderr
    assert not list(tmp_path.glob("picard.*")) and not list(tmp_path.glob("bammetrics.*"))
    RESULTS.append((f"{TOOL} fails: {' '.join(args)[-50:]}", 0, 0, "both exit non-zero, no file"))


@pytest.mark.parametrize("case", sorted(INVALID))
def test_invalid_input(case, tmp_path):
    path = write(tmp_path, case, INVALID[case])
    both_fail([f"I={path}"], tmp_path)


@pytest.mark.parametrize("args", [
    ["ASSUME_SORTED=true"], ["GA4GH_CLIENT_SECRETS=x"], ["NO_SUCH_OPTION=1"],
])
def test_bad_option(args, tmp_path):
    path = write(tmp_path, "ok", VALID["one_row"])
    both_fail([f"I={path}", *args], tmp_path)


def test_missing_input(tmp_path):
    both_fail([f"I={tmp_path / 'absent.txt'}"], tmp_path, "absent.txt")


DEFAULTS = {"0", ""}
EXPECTED_FILES = len(VALID)


def test_every_column_takes_three_non_default_values(request):
    if len(SEEN) != EXPECTED_FILES:
        if request.config.getoption("keyword") or any("::" in a for a in request.config.args):
            pytest.skip("run the whole module")
        pytest.fail(f"{len(SEEN)} of {EXPECTED_FILES} Picard files compared")
    values = defaultdict(set)
    columns = []
    for _, parsed in SEEN:
        columns = parsed.columns or columns
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3}
    assert not short, f"columns with fewer than 3 distinct non-default values: {short}"
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns), "columns checked"))

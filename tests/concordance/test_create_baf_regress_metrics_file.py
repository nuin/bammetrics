"""CreateBafRegressMetricsFile against Picard 3.5.0 on hand-written bafRegress outputs.

Each valid input is converted by both tools and the file under the prefix compared cell by
cell and byte for byte after the header lines (compare_prefixed). Each invalid input makes
both tools exit non-zero with no file left. The last test is the distinct-value gate.
"""
import os
import random
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS
from compare import compare_prefixed
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CreateBafRegressMetricsFile"
HEADER = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom"
SEEN = []


def rows(seed, n):
    """`n` bafRegress lines with numbers in the forms bafRegress prints."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        estimate = rng.choice([f"{rng.uniform(-0.01, 0.2):.6f}", f"{rng.uniform(0, 1e-3):.4e}", "0"])
        stderr = f"{rng.uniform(1e-5, 0.01):.6g}"
        tval = f"{rng.uniform(-5, 80):.5f}"
        pval = rng.choice([f"{rng.random():.6f}", f"{rng.uniform(1, 9):.3f}e-{rng.randint(5, 300)}", "0"])
        callrate = f"{rng.uniform(0.9, 1):.5f}"
        out.append("\t".join([f"S{seed}_{i}", estimate, stderr, tval, pval, callrate,
                              str(rng.randint(0, 700000))]))
    return out


VALID = {
    "three_rows": "\n".join([HEADER, "NA12878\t0.0012\t0.0003\t4.1\t2.1e-5\t0.995\t123456",
                             "NA12891\t-0.0004\t0.00021\t-1.9\t0.97\t0.9981\t130001",
                             "NA12892\t0.15\t0.002\t75.3\t0\t0.93\t99999"]) + "\n",
    "crlf": "\r\n".join([HEADER, *rows(1, 4)]) + "\r\n",
    "cr": "\r".join([HEADER, *rows(2, 4)]) + "\r",
    "no_final_line_end": "\n".join([HEADER, *rows(3, 3)]),
    "nan_and_infinity": "\n".join([HEADER, "a\tNaN\tInfinity\t-Infinity\t1\t-0\t0",
                                   "b\t+0.5\t-NaN\t+Infinity\t1e-400\t0.0\t-1",
                                   "c\t1\t2\t3\t-1\t4\t5"]) + "\n",
    "number_forms": "\n".join([HEADER, "a\t1E5\t2.5e+3\t0x1p-3\t1d\t3f\t+7",
                               "b\t.5\t5.\t1e-300\t4.9e-324\t1.7976931348623157E308\t2147483647",
                               "c\t0X1.8P1\t1e0\t-.25\t100\t0.1F\t-2147483648"]) + "\n",
    "spaces": "\n".join([HEADER, "x  1   2\t\t3 4 5 6  ", " 1 2 3 4 5 6",
                         "y\t1\x0b2\x0c3 4 5 6"]) + "\n",
    "header_only": HEADER + "\n",
    "header_no_line_end": HEADER,
    "utf8_names": "\n".join([HEADER, "échantillon\t0.1\t0.2\t0.3\t0.4\t0.5\t6",
                             "样本 一\t0.1\t0.2\t0.3\t0.4\t0.5\t6"]) + "\n",
    "many_rows": "\n".join([HEADER, *rows(7, 150)]) + "\n",
}

INVALID = {
    "empty": "",
    "blank_first_line": "\n" + HEADER + "\n",
    "header_extra_column": HEADER + "\textra\n",
    "header_missing_column": HEADER.rsplit("\t", 1)[0] + "\n",
    "header_spaces": HEADER.replace("\t", " ") + "\n",
    "header_trailing_space": HEADER + " \n",
    "header_upper_case": HEADER.replace("Nhom", "NHOM") + "\n",
    "byte_order_mark": "﻿" + HEADER + "\n",
    "row_missing_column": HEADER + "\nS\t1\t2\t3\t4\t5\n",
    "row_extra_column": HEADER + "\nS\t1\t2\t3\t4\t5\t6\t7\n",
    "empty_line": HEADER + "\nS\t1\t2\t3\t4\t5\t6\n\nT\t1\t2\t3\t4\t5\t6\n",
    "blank_line": HEADER + "\n   \n",
    "bad_number": HEADER + "\nS\t1\tabc\t3\t4\t5\t6\n",
    "nhom_decimal": HEADER + "\nS\t1\t2\t3\t4\t5\t6.5\n",
    "nhom_overflow": HEADER + "\nS\t1\t2\t3\t4\t5\t3000000000\n",
    "nhom_nan": HEADER + "\nS\t1\t2\t3\t4\t5\tNaN\n",
    "comma_decimal": HEADER + "\nS\t0,5\t2\t3\t4\t5\t6\n",
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


def test_bytes_that_are_not_utf8(tmp_path):
    """A byte Java's decoder replaces becomes U+FFFD in both; an encoded surrogate (ED A0
    80) is one U+FFFD in Java."""
    path = tmp_path / "latin1.txt"
    path.write_bytes(HEADER.encode() + b"\nS\xe9\xff\t1\t2\t3\t4\t5\t6\n"
                     + b"T\xed\xa0\x80\xc0\x80\t1\t2\t3\t4\t5\t6\n")
    cells = compare_prefixed(TOOL, [f"I={path}"], workdir=tmp_path, seen=SEEN)
    RESULTS.append((f"{TOOL} latin1", cells, cells, "every file identical"))


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
    ["ASSUME_SORTED=true"], ["GA4GH_CLIENT_SECRETS=x"], ["DBSNP=x"], ["NO_SUCH_OPTION=1"],
])
def test_bad_option(args, tmp_path):
    path = write(tmp_path, "ok", VALID["three_rows"])
    both_fail([f"I={path}", *args], tmp_path)


def test_missing_input_and_output(tmp_path):
    both_fail([f"I={tmp_path / 'absent.txt'}"], tmp_path, "absent.txt")
    path = write(tmp_path, "ok", VALID["three_rows"])
    for args in ([f"I={path}"], [f"O={tmp_path / 'x'}"]):
        picard = subprocess.run(["picard", TOOL, *args], env=picard_env(), capture_output=True,
                                text=True)
        ours = subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True)
        assert picard.returncode != 0 and ours.returncode != 0, args
        assert "is required" in ours.stderr, ours.stderr
    picard = subprocess.run(["picard", TOOL, f"I={path}", f"O={tmp_path / 'no' / 'p'}"],
                            env=picard_env(), capture_output=True, text=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={path}", f"O={tmp_path / 'no' / 'p'}"],
                          capture_output=True, text=True)
    assert picard.returncode != 0 and ours.returncode != 0
    assert not (tmp_path / "no").exists()


def test_prefix_with_a_trailing_separator(tmp_path):
    """java.io.File drops the trailing /: O=dir/ writes dir.bafregress_metrics beside dir."""
    path = write(tmp_path, "ok", VALID["three_rows"])
    for side in ("picard", "ours"):
        (tmp_path / side).mkdir()
    picard = subprocess.run(["picard", TOOL, f"I={path}", f"O={tmp_path / 'picard'}/"],
                            env=picard_env(), capture_output=True, text=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={path}", f"O={tmp_path / 'ours'}/"],
                          capture_output=True, text=True)
    assert picard.returncode == 0 and ours.returncode == 0, ours.stderr
    for side in ("picard", "ours"):
        assert (tmp_path / f"{side}.bafregress_metrics").is_file()
        assert not list((tmp_path / side).iterdir())


def test_output_hard_linked_to_the_input(tmp_path):
    """bammetrics alone: an OUTPUT that is a hard link to INPUT is refused and left whole."""
    target = tmp_path / "x.bafregress_metrics"
    target.write_text(VALID["three_rows"])
    link = tmp_path / "y.txt"
    os.link(target, link)
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={link}", f"O={tmp_path / 'x'}"],
                          capture_output=True, text=True)
    assert ours.returncode != 0 and "is an input" in ours.stderr, ours.stderr
    assert target.read_text() == VALID["three_rows"]
    assert os.stat(target).st_ino == os.stat(link).st_ino


def test_ignored_common_options(tmp_path):
    path = write(tmp_path, "ok", VALID["three_rows"])
    args = [f"I={path}", "R=ref.fa", "VALIDATION_STRINGENCY=LENIENT", "COMPRESSION_LEVEL=1",
            "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "MAX_RECORDS_IN_RAM=10",
            f"TMP_DIR={tmp_path}", "USE_JDK_DEFLATER=true"]
    cells = compare_prefixed(TOOL, args, workdir=tmp_path, seen=SEEN)
    RESULTS.append((f"{TOOL} common options", cells, cells, "every file identical"))


DEFAULTS = {"0", ""}
EXPECTED_FILES = len(VALID) + 2


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

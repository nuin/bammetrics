"""bammetrics vs Picard 3.5.0, field by field, on every case in cases.py.

Run with `make concordance` (builds the corpus, the Picard references and
target/release/bammetrics first). bammetrics runs through its CLI, so the
Python module is not needed.
"""
import json
import subprocess

import pytest

from cases import BAMMETRICS, CASES, METRICS, OURS_DIR, PICARD_DIR
from conftest import RESULTS
from picard_metrics import compare_files, parse_picard_metrics

MANIFEST = PICARD_DIR / "manifest.json"


def picard_manifest():
    if not MANIFEST.exists():
        pytest.skip("no Picard references: run `make reference-metrics`")
    return json.loads(MANIFEST.read_text())


def run_ours(case, picard_syntax=False):
    """Run bammetrics on a case: native flags, or the exact arguments Picard got."""
    if not BAMMETRICS.exists():
        pytest.fail(f"{BAMMETRICS} missing: run `cargo build --release`")
    OURS_DIR.mkdir(parents=True, exist_ok=True)
    case.ours_txt.unlink(missing_ok=True)
    metric = METRICS[case.metric]
    if picard_syntax:
        cmd = [str(BAMMETRICS), metric.picard_tool, f"I={case.bam}", f"O={case.ours_txt}",
               *metric.picard_args, *case.picard_args, "VERBOSITY=ERROR", "QUIET=true",
               "VALIDATION_STRINGENCY=SILENT"]
    else:
        cmd = [str(BAMMETRICS), metric.cli, str(case.bam), "-O", str(case.ours_txt),
               *case.our_args]
    subprocess.run(cmd, check=True, capture_output=True)


def check(case, picard_syntax):
    entry = picard_manifest().get(case.id)
    if entry is None:
        pytest.fail(f"{case.id} has no Picard reference: run `make reference-metrics`")
    run_ours(case, picard_syntax)
    if not entry["written"]:
        # Picard wrote no metrics file; bammetrics must not write one either.
        assert not case.ours_txt.exists(), "Picard wrote no file, bammetrics did"
        return entry, None
    assert case.ours_txt.exists(), "Picard wrote a file, bammetrics did not"
    picard = parse_picard_metrics(case.picard_txt)
    ours = parse_picard_metrics(case.ours_txt)
    result = compare_files(picard, ours, METRICS[case.metric].key_columns)
    assert result.ok, "\n".join(result.mismatches)
    # Whole file from `## METRICS CLASS` on: header row, row count, histogram, blank lines.
    assert body(case.picard_txt) == body(case.ours_txt), "values match but the text differs"
    assert headers(case.ours_txt)[1].startswith("# Started on: "), "no Started on header"
    return entry, result


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_concordance(case):
    entry, result = check(case, picard_syntax=False)
    if result is None:
        RESULTS.append((case.id, 0, 0, "no file from either (no pairs kept)"))
    else:
        RESULTS.append((case.id, result.matched, result.total,
                        "byte-identical from METRICS CLASS on"))


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_picard_syntax(case):
    """`bammetrics <PicardTool> I= O= ...` with Picard's own arguments gives the same file."""
    check(case, picard_syntax=True)


def headers(path):
    """The StringHeader text lines (`# ...`) before `## METRICS CLASS`."""
    text = path.read_text()
    return [line for line in text[:text.find("## METRICS CLASS")].splitlines()
            if line.startswith("# ")]


def body(path):
    """The file from `## METRICS CLASS` on; the StringHeader lines differ by design."""
    text = path.read_text()
    return text[text.find("## METRICS CLASS"):]

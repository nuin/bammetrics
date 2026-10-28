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


def run_ours(case):
    if not BAMMETRICS.exists():
        pytest.fail(f"{BAMMETRICS} missing: run `cargo build --release`")
    OURS_DIR.mkdir(parents=True, exist_ok=True)
    case.ours_txt.unlink(missing_ok=True)
    metric = METRICS[case.metric]
    subprocess.run([str(BAMMETRICS), metric.cli, str(case.bam), "-O", str(case.ours_txt),
                    *case.our_args], check=True, capture_output=True)


@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_concordance(case):
    entry = picard_manifest().get(case.id)
    if entry is None:
        pytest.fail(f"{case.id} has no Picard reference: run `make reference-metrics`")
    run_ours(case)
    if not entry["written"]:
        # Picard wrote no metrics file; bammetrics must not write one either.
        assert not case.ours_txt.exists(), "Picard wrote no file, bammetrics did"
        RESULTS.append((case.id, 0, 0, "no file from either (no pairs kept)"))
        return
    assert case.ours_txt.exists(), "Picard wrote a file, bammetrics did not"
    picard = parse_picard_metrics(case.picard_txt)
    ours = parse_picard_metrics(case.ours_txt)
    result = compare_files(picard, ours, METRICS[case.metric].key_columns)
    if not result.ok:
        note = f"{len(result.mismatches)} mismatches"
    elif body(case.picard_txt) == body(case.ours_txt):
        note = "byte-identical from METRICS CLASS on"
    else:
        note = "values match, text differs"
    RESULTS.append((case.id, result.matched, result.total, note))
    assert result.ok, "\n".join(result.mismatches)


def body(path):
    """The file from `## METRICS CLASS` on; the StringHeader lines differ by design."""
    text = path.read_text()
    return text[text.find("## METRICS CLASS"):]

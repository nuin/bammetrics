"""bammetrics vs Picard 3.5.0, field by field, on every case in cases.py.

Run with `make concordance` (builds the corpus, the Picard references and
target/release/bammetrics first). bammetrics runs through its CLI, so the
Python module is not needed.
"""
import json
import subprocess
from pathlib import Path

import pytest

from cases import (BAMMETRICS, CASES, CASES_BY_ID, COLLECT_CASES, METRICS, OURS_DIR,
                   PICARD_DIR)
from conftest import RESULTS
from picard_metrics import compare_files, parse_picard_metrics

MANIFEST = PICARD_DIR / "manifest.json"


def picard_manifest():
    if not MANIFEST.exists():
        pytest.skip("no Picard references: run `make reference-metrics`")
    return json.loads(MANIFEST.read_text())


# BGZF worker counts the native runs use; the metrics must not depend on them.
THREADS = (1, 8)


def run_ours(case, picard_syntax=False, threads=None):
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
        if threads is not None:
            cmd += ["--threads", str(threads)]
    subprocess.run(cmd, check=True, capture_output=True)


def check(case, picard_syntax, threads=None):
    entry = reference(case)
    run_ours(case, picard_syntax, threads)
    return entry, compare_output(case, entry, case.ours_txt)


def reference(case):
    """The manifest entry of a case's Picard reference."""
    entry = picard_manifest().get(case.id)
    if entry is None:
        pytest.fail(f"{case.id} has no Picard reference: run `make reference-metrics`")
    return entry


def compare_output(case, entry, ours_txt):
    """Compare one file bammetrics wrote with the Picard reference of `case`;
    None when neither wrote a file."""
    if not entry["written"]:
        # Picard wrote no metrics file; bammetrics must not write one either.
        assert not ours_txt.exists(), "Picard wrote no file, bammetrics did"
        return None
    assert ours_txt.exists(), "Picard wrote a file, bammetrics did not"
    picard = parse_picard_metrics(case.picard_txt)
    ours = parse_picard_metrics(ours_txt)
    result = compare_files(picard, ours, METRICS[case.metric].key_columns)
    assert result.ok, "\n".join(result.mismatches)
    # Whole file from `## METRICS CLASS` on: header row, row count, histogram, blank lines.
    assert body(case.picard_txt) == body(ours_txt), "values match but the text differs"
    assert headers(ours_txt)[1].startswith("# Started on: "), "no Started on header"
    return result


@pytest.mark.parametrize("threads", THREADS)
@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_concordance(case, threads):
    _, result = check(case, picard_syntax=False, threads=threads)
    label = f"{case.id} --threads {threads}"
    if result is None:
        RESULTS.append((label, 0, 0, "no file from either (no pairs kept)"))
    else:
        RESULTS.append((label, result.matched, result.total,
                        "byte-identical from METRICS CLASS on"))


@pytest.mark.parametrize("threads", THREADS)
@pytest.mark.parametrize("collect", COLLECT_CASES, ids=[c.name for c in COLLECT_CASES])
def test_collect(collect, threads):
    """One `collect` pass writes, per metric, the file the single-tool case wrote."""
    if not BAMMETRICS.exists():
        pytest.fail(f"{BAMMETRICS} missing: run `cargo build --release`")
    collect.prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs = [Path(f"{collect.prefix}{suffix}") for suffix, _ in collect.expected]
    for path in outputs:
        path.unlink(missing_ok=True)
    subprocess.run([str(BAMMETRICS), "collect", str(collect.bam), "-O", str(collect.prefix),
                    "--threads", str(threads), *collect.our_args],
                   check=True, capture_output=True)
    for (suffix, case_id), path in zip(collect.expected, outputs):
        case = CASES_BY_ID[case_id]
        result = compare_output(case, reference(case), path)
        label = f"collect {collect.name}{suffix} --threads {threads}"
        if result is None:
            RESULTS.append((label, 0, 0, f"no file from either ({case_id})"))
        else:
            RESULTS.append((label, result.matched, result.total, f"= {case_id}"))


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

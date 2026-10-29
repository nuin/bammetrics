"""bammetrics vs Picard 3.5.0, field by field, on every case in cases.py.

Run with `make concordance` (builds the corpus, the Picard references and
target/release/bammetrics first). bammetrics runs through its CLI, so the
Python module is not needed.
"""
import gzip
import json
import subprocess
from pathlib import Path

import pytest

from cases import (BAMMETRICS, CASES, CASES_BY_ID, COLLECT_CASES, CORPUS, METRICS, OURS_DIR,
                   PICARD_DIR, REFERENCE)
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
        if case.per_base:
            cmd.append(f"PER_BASE_COVERAGE={case.ours_per_base}")
    else:
        cmd = [str(BAMMETRICS), metric.cli, str(case.bam), "-O", str(case.ours_txt),
               *case.our_args]
        if case.per_base:
            cmd += ["--per-base-coverage", str(case.ours_per_base)]
        if threads is not None:
            cmd += ["--threads", str(threads)]
    if case.per_base:
        case.ours_per_base.unlink(missing_ok=True)
    subprocess.run(cmd, check=True, capture_output=True)


def check(case, picard_syntax, threads=None):
    entry = reference(case)
    run_ours(case, picard_syntax, threads)
    result = compare_output(case, entry, case.ours_txt)
    if case.per_base:
        compare_per_base(case, case.ours_per_base)
    return entry, result


def per_base_text(path):
    """A PER_BASE_COVERAGE file as text, decompressed when its name says gzip."""
    data = path.read_bytes()
    return gzip.decompress(data) if path.name.endswith(".gz") else data


def compare_per_base(case, ours):
    """The whole per-base file, byte for byte (after gunzip for .gz names)."""
    assert ours.exists(), "Picard wrote a per-base file, bammetrics did not"
    picard, mine = per_base_text(case.picard_per_base), per_base_text(ours)
    assert picard.count(b"\n") == mine.count(b"\n"), "per-base line counts differ"
    assert picard == mine, "per-base files differ"
    if ours.name.endswith(".gz"):
        assert ours.read_bytes()[:2] == b"\x1f\x8b", "not gzip"
    PER_BASE_LINES[case.id] = picard.count(b"\n")


# Lines in each per-base file checked, for the summary.
PER_BASE_LINES = {}


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
    if case.per_base:
        lines = PER_BASE_LINES[case.id]
        RESULTS.append((f"{label} per-base", lines, lines, "lines, whole file byte-identical"))


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
        if suffix.endswith(".per_base_coverage"):
            compare_per_base(case, path)
            lines = PER_BASE_LINES[case.id]
            RESULTS.append((f"collect {collect.name}{suffix} --threads {threads}", lines, lines,
                            f"per-base lines byte-identical = {case_id}"))
            continue
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


PIPELINE_CASES = [c for c in CASES if c.name in ("t1_pb_q0", "t2_pb_q0")]


@pytest.mark.parametrize("case", PIPELINE_CASES, ids=[c.id for c in PIPELINE_CASES])
def test_pipeline_coverage_command(case):
    """The clinical pipeline's coverage command, argument for argument, with
    bammetrics in place of `java -jar picard.jar`."""
    entry = reference(case)
    panel = CORPUS / "targets.interval_list"
    OURS_DIR.mkdir(parents=True, exist_ok=True)
    summary = OURS_DIR / f"{case.name}.summary.partial.out"
    per_base = OURS_DIR / f"{case.name}.nucl.partial.out"
    for path in (summary, per_base):
        path.unlink(missing_ok=True)
    subprocess.run([str(BAMMETRICS), "CollectHsMetrics", f"BI={panel}", f"I={case.bam}",
                    f"PER_BASE_COVERAGE={per_base}", "MINIMUM_MAPPING_QUALITY=0",
                    "MINIMUM_BASE_QUALITY=0", f"TARGET_INTERVALS={panel}", f"OUTPUT={summary}",
                    f"R={REFERENCE}", "QUIET=true"], check=True, capture_output=True)
    result = compare_output(case, entry, summary)
    compare_per_base(case, per_base)
    RESULTS.append((f"pipeline command {case.name}", result.matched, result.total,
                    f"summary and {PER_BASE_LINES[case.id]} per-base lines byte-identical"))


def headers(path):
    """The StringHeader text lines (`# ...`) before `## METRICS CLASS`."""
    text = path.read_text()
    return [line for line in text[:text.find("## METRICS CLASS")].splitlines()
            if line.startswith("# ")]


def body(path):
    """The file from `## METRICS CLASS` on; the StringHeader lines differ by design."""
    text = path.read_text()
    return text[text.find("## METRICS CLASS"):]

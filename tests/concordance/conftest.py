"""Collects one concordance line per case and prints the table at the end."""
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from cases import CORPUS, ROOT  # noqa: E402

RESULTS = []

MAKE_CORPUS = ROOT / "tests" / "data" / "make_corpus.py"
MAKE_EXOME = ROOT / "tests" / "data" / "make_exome.py"
FOUNDATION = CORPUS / "foundation"
EXOME = CORPUS / "exome"


def generate(out, seed=None, script=MAKE_CORPUS):
    """Run a tests/data generator script into `out`."""
    cmd = [sys.executable, str(script), str(out)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    subprocess.run(cmd, check=True, capture_output=True)


def generated(out, bam, script):
    """`out` with `bam` written by `script`, regenerated when the script is newer."""
    if not (out / bam).exists() or (out / bam).stat().st_mtime < script.stat().st_mtime:
        generate(out, script=script)
    return out


@pytest.fixture(scope="session")
def corpus():
    """Directory of the generated multi-sample corpus (ref.fa, corpus.bam, interval lists)."""
    return generated(FOUNDATION, "corpus.bam", MAKE_CORPUS)


@pytest.fixture(scope="session")
def exome():
    """Directory of the exome-like slice (ref.fa, exome.bam, targets, baits, gc_targets)."""
    return generated(EXOME, "exome.bam", MAKE_EXOME)


def pytest_terminal_summary(terminalreporter):
    if not RESULTS:
        return
    terminalreporter.section("concordance (fields matched / total)")
    for case_id, matched, total, note in RESULTS:
        terminalreporter.write_line(f"{case_id:40s} {matched:6d} / {total:<6d} {note}")

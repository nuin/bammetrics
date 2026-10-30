"""Collects one concordance line per case and prints the table at the end."""
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))

from cases import CORPUS, ROOT  # noqa: E402

RESULTS = []

MAKE_CORPUS = ROOT / "tests" / "data" / "make_corpus.py"
FOUNDATION = CORPUS / "foundation"


def generate(out, seed=None):
    """Run tests/data/make_corpus.py into `out`."""
    cmd = [sys.executable, str(MAKE_CORPUS), str(out)]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    subprocess.run(cmd, check=True, capture_output=True)


@pytest.fixture(scope="session")
def corpus():
    """Directory of the generated multi-sample corpus (ref.fa, corpus.bam, interval lists)."""
    if not (FOUNDATION / "corpus.bam").exists() or (
            (FOUNDATION / "corpus.bam").stat().st_mtime < MAKE_CORPUS.stat().st_mtime):
        generate(FOUNDATION)
    return FOUNDATION


def pytest_terminal_summary(terminalreporter):
    if not RESULTS:
        return
    terminalreporter.section("concordance (fields matched / total)")
    for case_id, matched, total, note in RESULTS:
        terminalreporter.write_line(f"{case_id:40s} {matched:6d} / {total:<6d} {note}")

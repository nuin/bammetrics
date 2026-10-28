"""Collects one concordance line per case and prints the table at the end."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

RESULTS = []


def pytest_terminal_summary(terminalreporter):
    if not RESULTS:
        return
    terminalreporter.section("concordance (fields matched / total)")
    for case_id, matched, total, note in RESULTS:
        terminalreporter.write_line(f"{case_id:40s} {matched:6d} / {total:<6d} {note}")

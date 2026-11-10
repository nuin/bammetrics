"""Tests of the harness itself: no Picard and no release binary, only samtools.

`fake_tools` gives each test a `picard` and a `bammetrics` that copy files, so
`compare_bam`, `compare_vcf` and the Picard queue run end to end in seconds.
"""
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "concordance"))
sys.path.insert(0, str(ROOT / "tools"))

# A stand-in for either tool: `TOOL KEY=VALUE...` copies `I=` (or the file
# FAKE_PLANTED names, for the side FAKE_PLANT names) to `O=` and FAKE_METRICS
# to `M=`. `FAIL=<side>` or `FAIL=both` exits 1 before any copy.
FAKE_TOOL = '''#!{python}
import os, shutil, sys
side = "{side}"
args = dict(a.split("=", 1) for a in sys.argv[2:] if "=" in a)
if args.get("FAIL") in (side, "both"):
    sys.exit(1)
planted = os.environ.get("FAKE_PLANTED") if os.environ.get("FAKE_PLANT") == side else None
shutil.copy(planted or args["I"], args["O"])
if "M" in args:
    shutil.copy(os.environ["FAKE_METRICS"], args["M"])
'''


def write_tool(path, side):
    path.write_text(FAKE_TOOL.format(python=sys.executable, side=side))
    path.chmod(0o755)
    return path


@pytest.fixture
def queue_env(tmp_path, monkeypatch):
    """A lock file of this test's own and no wait for JVMs of other agents."""
    monkeypatch.setenv("PICARD_LOCK_FILE", str(tmp_path / "picard.lock"))
    monkeypatch.setenv("PICARD_LOCK_WAIT_FOREIGN", "0")
    monkeypatch.delenv("PICARD_LOCK_TIMEOUT", raising=False)
    monkeypatch.delenv("PICARD_LOCK_HELD", raising=False)
    return tmp_path / "picard.lock"


@pytest.fixture
def fake_tools(tmp_path, monkeypatch, queue_env):
    """`picard` (first on PATH, behind the queue shim) and `compare.BAMMETRICS`
    replaced by the copying stand-in."""
    import compare

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    write_tool(bin_dir / "picard", "picard")
    monkeypatch.setattr(compare, "BAMMETRICS", write_tool(bin_dir / "bammetrics", "bammetrics"))
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return bin_dir

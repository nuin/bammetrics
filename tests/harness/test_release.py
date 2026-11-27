"""The release run shims (scripts/linux-bammetrics) and the BAMMETRICS_BIN
override of the concordance harness. A stand-in `docker` records its argv."""
import importlib
import os
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import ROOT

SHIM = ROOT / "scripts" / "linux-bammetrics"
DOCKER = """#!/bin/sh
for a in "$@"; do echo "$a"; done > "$DOCKER_LOG"
exit "${DOCKER_EXIT:-0}"
"""


def release_env():
    values = {}
    for line in (ROOT / "scripts" / "release.env").read_text().splitlines():
        if "=" in line and not line.startswith("#"):
            key, _, value = line.partition("=")
            values[key] = value.strip('"')
    return values


@pytest.fixture
def shim_env(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    docker = bin_dir / "docker"
    docker.write_text(DOCKER)
    docker.chmod(0o755)
    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("DOCKER_LOG", str(tmp_path / "argv"))
    monkeypatch.setenv("BM_LINUX_ARCH", "arm64")
    monkeypatch.setenv("BM_LINUX_BIN", "/Users/x/bammetrics")
    return work


def run_shim(cwd, *args):
    return subprocess.run([str(SHIM), *args], cwd=cwd, capture_output=True, text=True)


def test_linux_shim_mounts_cwd_and_temp(shim_env, tmp_path):
    result = run_shim(shim_env, "tools", "--format", "json")
    assert result.returncode == 0, result.stderr
    argv = (tmp_path / "argv").read_text().split("\n")[:-1]
    assert argv[:3] == ["run", "--rm", "-i"]
    assert argv[argv.index("--platform") + 1] == "linux/arm64"
    assert argv[argv.index("--user") + 1] == f"{os.getuid()}:{os.getgid()}"
    assert argv[argv.index("-w") + 1] == str(shim_env.resolve())
    mounts = [argv[i + 1] for i, a in enumerate(argv) if a == "-v"]
    for mount in ("/Users:/Users", "/private:/private", "/private/var/folders:/var/folders"):
        assert mount in mounts
    envs = [argv[i + 1] for i, a in enumerate(argv) if a == "-e"]
    for name in ("TZ", "XDG_CONFIG_HOME", "BAMMETRICS_REFERENCE"):
        assert name in envs
    assert release_env()["RUNNER_IMAGE"] in argv
    assert argv[-4:] == ["/Users/x/bammetrics", "tools", "--format", "json"]


def test_linux_shim_passes_exit_code(shim_env, monkeypatch):
    monkeypatch.setenv("DOCKER_EXIT", "11")
    assert run_shim(shim_env, "x").returncode == 11


@pytest.mark.parametrize("name", ["BM_LINUX_ARCH", "BM_LINUX_BIN"])
def test_linux_shim_refuses_unset(shim_env, monkeypatch, name):
    monkeypatch.delenv(name)
    result = run_shim(shim_env, "x")
    assert result.returncode == 2
    assert name in result.stderr


def test_cases_reads_bammetrics_bin(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "tests" / "concordance"))
    monkeypatch.setenv("BAMMETRICS_BIN", "/x/y")
    cases = importlib.reload(importlib.import_module("cases"))
    assert cases.BAMMETRICS == Path("/x/y")
    monkeypatch.delenv("BAMMETRICS_BIN")
    cases = importlib.reload(cases)
    assert cases.BAMMETRICS.parts[-3:] == ("target", "release", "bammetrics")

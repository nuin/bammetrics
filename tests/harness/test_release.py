"""The release run shims (scripts/linux-bammetrics) and the BAMMETRICS_BIN
override of the concordance harness. A stand-in `docker` records its argv."""
import importlib
import os
import shutil
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


# scripts/release.sh, planned with --dry-run in a throwaway repository.

RELEASE = ROOT / "scripts" / "release.sh"
STUB = """#!/bin/sh
case "$(basename "$0") $*" in
  "rustup target list"*) printf 'aarch64-apple-darwin\\nx86_64-apple-darwin\\n' ;;
esac
exit 0
"""
LINUX = ("x86_64-unknown-linux-musl", "aarch64-unknown-linux-musl")


def git(repo, *args):
    subprocess.run(
        ["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t", *args],
        check=True, capture_output=True,
    )


@pytest.fixture
def release_repo(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    (repo / "scripts").mkdir(parents=True)
    for name in ("release.sh", "release.env"):
        (repo / "scripts" / name).write_text((ROOT / "scripts" / name).read_text())
    (repo / "scripts" / "release.sh").chmod(0o755)
    (repo / "Cargo.toml").write_text('[workspace.package]\nversion = "1.2.3"\n')
    (repo / "README.md").write_text("readme\n")
    (repo / "CHANGELOG.md").write_text("# Changelog\n\n## Unreleased\n\n## 1.2.3\n\nNotes.\n")
    git(repo, "init", "-q", "-b", "main")
    git(repo, "add", ".")
    git(repo, "commit", "-q", "-m", "one")
    (repo / "README.md").write_text("readme two\n")
    git(repo, "commit", "-q", "-am", "two")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name in ("docker", "lipo", "rustup", "gh"):
        (bin_dir / name).write_text(STUB)
        (bin_dir / name).chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return repo


def release(repo, *args, env=None):
    return subprocess.run(
        [str(repo / "scripts" / "release.sh"), *args],
        cwd=repo, capture_output=True, text=True, env={**os.environ, **(env or {})},
    )


def test_dirty_tree_refused(release_repo):
    (release_repo / "README.md").write_text("changed\n")
    result = release(release_repo, "--snapshot", "--dry-run")
    assert result.returncode == 2
    assert "README.md" in result.stderr
    assert not (release_repo / "dist").exists()


def test_tag_must_point_at_head(release_repo):
    result = release(release_repo, "--tag", "--dry-run")
    assert result.returncode == 2
    assert "v1.2.3" in result.stderr
    git(release_repo, "tag", "v1.2.3", "HEAD~1")
    result = release(release_repo, "--tag", "--dry-run")
    assert result.returncode == 2
    assert "v1.2.3" in result.stderr


def test_publish_needs_tag_mode(release_repo):
    result = release(release_repo, "--publish", "--snapshot", "--dry-run")
    assert result.returncode == 2
    assert "--publish" in result.stderr


def test_unknown_option_refused(release_repo):
    result = release(release_repo, "--nope", "--dry-run")
    assert result.returncode == 2
    assert "--nope" in result.stderr


def test_missing_tool_named(release_repo, monkeypatch):
    bin_dir = release_repo.parent / "bin"
    (bin_dir / "lipo").unlink()
    for tool in ("bash", "env", "git", "grep", "sort", "tr", "cut", "dirname", "cat", "id",
                 "shasum", "python3", "awk", "head", "ls", "sh", "uname", "mkdir", "rm"):
        found = shutil.which(tool)
        if found:
            (bin_dir / tool).symlink_to(found)
    monkeypatch.setenv("PATH", str(bin_dir))
    result = release(release_repo, "--snapshot", "--dry-run")
    assert result.returncode == 2
    assert "lipo" in result.stderr


def index_of(lines, *needles):
    for i, line in enumerate(lines):
        if all(n in line for n in needles):
            return i
    raise AssertionError(f"no line with {needles}:\n" + "\n".join(lines))


def test_dry_run_plan(release_repo):
    result = release(release_repo, "--snapshot", "--dry-run")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    steps = []
    for target, platform in zip(LINUX, ("linux/amd64", "linux/arm64")):
        steps.append(index_of(
            lines, "--platform " + platform,
            "cargo auditable build --release --locked -p bammetrics-cli --features pdf --target " + target,
        ))
    for target in ("aarch64-apple-darwin", "x86_64-apple-darwin"):
        steps.append(index_of(
            lines, "rustup run 1.100.0 cargo build", "--features pdf", "--target " + target))
    steps.append(index_of(lines, "lipo -create"))
    for target in (*LINUX, "universal-apple-darwin"):
        steps.append(index_of(lines, "tar ", f"bammetrics-1.2.3-{target}.tar.gz"))
    steps.append(index_of(lines, "SHA256SUMS"))
    assert steps == sorted(steps) and len(set(steps)) == len(steps)
    assert not any("gh release" in line for line in lines)
    assert (release_repo / "dist").exists() is False


def test_dry_run_publish_line(release_repo):
    git(release_repo, "tag", "v1.2.3")
    result = release(release_repo, "--tag", "--publish", "--dry-run")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert "dist/v1.2.3" in result.stdout
    assert index_of(lines, "gh release create v1.2.3") > index_of(lines, "SHA256SUMS")


def test_blocked_target_skipped(release_repo):
    result = release(
        release_repo, "--snapshot", "--dry-run",
        env={"BLOCKED_TARGETS": "aarch64-unknown-linux-musl"},
    )
    assert result.returncode == 0, result.stderr
    assert "aarch64-unknown-linux-musl" not in result.stdout
    assert result.stderr.count("aarch64-unknown-linux-musl") == 1
    assert "x86_64-unknown-linux-musl" in result.stdout

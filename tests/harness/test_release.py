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
    for name in ("release.sh", "release.env", "sbom_pin.py"):
        (repo / "scripts" / name).write_text((ROOT / "scripts" / name).read_text())
    (repo / "scripts" / "release.sh").chmod(0o755)
    (repo / "Cargo.toml").write_text('[workspace.package]\nversion = "1.2.3"\n')
    (repo / "README.md").write_text("readme\n")
    (repo / "LICENSE-MIT").write_text("mit\n")
    (repo / "LICENSE-APACHE").write_text("apache\n")
    (repo / "CHANGELOG.md").write_text("# Changelog\n\n## Unreleased\n\n## 1.2.3\n\nNotes.\n")
    git(repo, "init", "-q", "-b", "main")
    git(repo, "add", ".")
    git(repo, "commit", "-q", "-m", "one")
    (repo / "README.md").write_text("readme two\n")
    git(repo, "commit", "-q", "-am", "two")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name in ("docker", "lipo", "rustup", "gh", "file", "arch"):
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
    steps.append(index_of(lines, "cargo cyclonedx"))
    steps.append(index_of(lines, "shasum -a 256"))
    steps.append(index_of(lines, "brew_formula.py"))
    assert steps == sorted(steps) and len(set(steps)) == len(steps)
    assert not any("gh release" in line for line in lines)
    assert (release_repo / "dist").exists() is False


def test_dry_run_publish_line(release_repo):
    git(release_repo, "tag", "v1.2.3")
    result = release(release_repo, "--tag", "--publish", "--dry-run")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    assert "dist/v1.2.3" in result.stdout
    assert index_of(lines, "gh release create v1.2.3") > index_of(lines, "shasum -a 256")


def test_blocked_target_skipped(release_repo):
    result = release(
        release_repo, "--snapshot", "--dry-run",
        env={"BLOCKED_TARGETS": "aarch64-unknown-linux-musl"},
    )
    assert result.returncode == 0, result.stderr
    plan = [l for l in result.stdout.splitlines() if "brew_formula.py" not in l]
    assert not any("aarch64-unknown-linux-musl" in l for l in plan)
    assert result.stderr.count("aarch64-unknown-linux-musl") == 1
    assert "x86_64-unknown-linux-musl" in result.stdout


def test_sbom_line_in_dry_run(release_repo):
    result = release(release_repo, "--snapshot", "--dry-run")
    assert result.returncode == 0, result.stderr
    lines = result.stdout.splitlines()
    sbom = index_of(lines, "cargo cyclonedx")
    assert "--spec-version 1.5" in lines[sbom]
    assert "crates/cli/Cargo.toml" in lines[sbom]
    assert "bammetrics-1.2.3.cdx.json" in "\n".join(lines)
    assert sbom < index_of(lines, "shasum -a 256") < index_of(lines, "brew_formula.py")


def test_sums_list_only_planned_files(release_repo, tmp_path):
    (tmp_path / "bin" / "docker").write_text(
        '#!/bin/sh\ncase "$*" in *cyclonedx*) echo \'{"metadata": {}}\' > "$SBOM_STUB" ;; esac\nexit 0\n')
    sha = subprocess.run(["git", "-C", str(release_repo), "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    out = release_repo / "dist" / f"1.2.3-{sha}"
    out.mkdir(parents=True)
    touch = out / "bammetrics-1.2.3.cdx.json"
    (out / "bammetrics-1.2.3-aarch64-unknown-linux-musl.tar.gz").write_text("stale")
    (out / "bammetrics-0.0.1.cdx.json").write_text("stale")
    result = release(
        release_repo, "--snapshot",
        env={"BLOCKED_TARGETS": " ".join((*LINUX, "universal-apple-darwin")), "SBOM_STUB": str(touch)},
    )
    assert result.returncode == 0, result.stderr
    assert not list(out.glob("*.tar.gz"))
    assert not (out / "bammetrics-0.0.1.cdx.json").exists()
    sums = (out / "SHA256SUMS").read_text()
    assert "aarch64" not in sums and "0.0.1" not in sums
    assert sums.split()[1:] == ["bammetrics-1.2.3.cdx.json"]


# scripts/brew_formula.py

BREW = ROOT / "scripts" / "brew_formula.py"
TARGETS = ("universal-apple-darwin", "x86_64-unknown-linux-musl", "aarch64-unknown-linux-musl")


def sums_file(tmp_path, targets=TARGETS, version="0.34.2"):
    lines = [f"{str(i) * 64}  bammetrics-{version}-{t}.tar.gz" for i, t in enumerate(targets, 1)]
    path = tmp_path / "SHA256SUMS"
    path.write_text("\n".join(lines) + "\n")
    return path


def brew(tmp_path, sums, *extra):
    out = tmp_path / "bammetrics.rb"
    result = subprocess.run(
        [sys.executable, str(BREW), "0.34.2", str(sums), str(out), *extra],
        capture_output=True, text=True)
    return result, out


def test_formula_from_sums(tmp_path):
    result, out = brew(tmp_path, sums_file(tmp_path))
    assert result.returncode == 0, result.stderr
    text = out.read_text()
    base = "https://github.com/nuin/bammetrics/releases/download/v0.34.2/bammetrics-0.34.2-"
    for want in (
        "class Bammetrics < Formula", 'version "0.34.2"',
        'homepage "https://github.com/nuin/bammetrics"',
        "on_macos do", "on_linux do", "on_intel do", "on_arm do",
        f'url "{base}universal-apple-darwin.tar.gz"', "1" * 64,
        f'url "{base}x86_64-unknown-linux-musl.tar.gz"', "2" * 64,
        f'url "{base}aarch64-unknown-linux-musl.tar.gz"', "3" * 64,
        'bin.install "bammetrics"',
        "test do", 'shell_output("#{bin}/bammetrics --version")', "0.34.2",
    ):
        assert want in text, want
    assert 'license any_of: ["MIT", "Apache-2.0"]' in text


def test_formula_missing_target(tmp_path):
    result, out = brew(tmp_path, sums_file(tmp_path, TARGETS[:2]))
    assert result.returncode == 2
    assert "aarch64-unknown-linux-musl" in result.stderr
    assert not out.exists()


def test_formula_skips_blocked_target(tmp_path):
    result, out = brew(tmp_path, sums_file(tmp_path, TARGETS[:2]),
                       "--skip", "aarch64-unknown-linux-musl")
    assert result.returncode == 0, result.stderr
    text = out.read_text()
    assert "aarch64" not in text and "on_intel do" in text and "x86_64-unknown-linux-musl" in text


EXPECTED_SERIAL = "d4ca3c70-334d-5979-a8e3-59e98d8aaf5c"


def test_sbom_pin_sets_timestamp_and_serial(tmp_path):
    sys.path.insert(0, str(ROOT / "scripts"))
    sbom_pin = importlib.import_module("sbom_pin")
    sha = "e8e34200000000000000000000000000000000aa"
    source = (
        '{\n  "bomFormat": "CycloneDX",\n  "serialNumber": "urn:uuid:random",\n'
        '  "metadata": {\n    "timestamp": "2026-11-27T01:57:27.097447295Z"\n  },\n'
        '  "components": []\n}'
    )
    text = sbom_pin.pin(source, sha, "1795000000")
    assert text == source.replace("urn:uuid:random", "urn:uuid:" + EXPECTED_SERIAL).replace(
        "2026-11-27T01:57:27.097447295Z", "2026-11-18T11:06:40Z"
    )
    assert sbom_pin.pin(source, sha, "1795000000") == text
    path = tmp_path / "x.cdx.json"
    path.write_text(source)
    assert sbom_pin.main(["sbom_pin.py", str(path), sha, "1795000000"]) == 0
    assert path.read_text() == text


# Review fixes: no host path in artifacts, licences, publish checks, pins.

def planned(repo, *args, env=None):
    result = release(repo, *args, "--dry-run", env=env)
    assert result.returncode == 0, result.stderr
    return result.stdout.splitlines()


def test_sbom_container_mounts_repo_at_build(release_repo):
    lines = planned(release_repo, "--snapshot")
    sbom = lines[index_of(lines, "cargo cyclonedx") - 0]
    assert "/build" in sbom
    assert f"{release_repo}:{release_repo}" not in sbom
    assert "-w /build" in sbom


def test_leak_check_in_plan_before_sums(release_repo):
    lines = planned(release_repo, "--snapshot")
    assert index_of(lines, "leak check") < index_of(lines, "shasum -a 256")


def stub_sbom(release_repo, tmp_path, text):
    sha = subprocess.run(["git", "-C", str(release_repo), "rev-parse", "--short", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    out = release_repo / "dist" / f"1.2.3-{sha}"
    out.mkdir(parents=True)
    (tmp_path / "bin" / "docker").write_text(
        f'#!/bin/sh\ncase "$*" in *cyclonedx*) echo \'{text}\' > "$SBOM_STUB" ;; esac\nexit 0\n')
    env = {"BLOCKED_TARGETS": " ".join((*LINUX, "universal-apple-darwin")),
           "SBOM_STUB": str(out / "bammetrics-1.2.3.cdx.json")}
    return out, env


def test_leak_check_fires_on_repo_path(release_repo, tmp_path):
    out, env = stub_sbom(release_repo, tmp_path,
                         '{"metadata": {}, "ref": "path+file://%s/crates/x"}' % release_repo)
    result = release(release_repo, "--snapshot", env=env)
    assert result.returncode == 1
    assert "cdx.json" in result.stderr and str(release_repo) in result.stderr
    assert not (out / "SHA256SUMS").exists()


def test_leak_check_fires_on_home(release_repo, tmp_path):
    out, env = stub_sbom(release_repo, tmp_path, '{"metadata": {}, "p": "/home/zed/x"}')
    result = release(release_repo, "--snapshot", env={**env, "HOME": "/home/zed"})
    assert result.returncode == 1
    assert "/home/zed" in result.stderr


def test_leak_check_passes_clean(release_repo, tmp_path):
    out, env = stub_sbom(release_repo, tmp_path, '{"metadata": {}}')
    result = release(release_repo, "--snapshot", env=env)
    assert result.returncode == 0, result.stderr
    assert (out / "SHA256SUMS").exists()
    info = (out / "BUILDINFO").read_text()
    assert "cargo_auditable " + release_env()["CARGO_AUDITABLE_VERSION"] in info
    assert "cargo_cyclonedx " + release_env()["CARGO_CYCLONEDX_VERSION"] in info


def test_tarballs_get_both_licences(release_repo):
    lines = planned(release_repo, "--snapshot")
    assert any("rm -f" in l and "SHA256SUMS" in l and "BUILDINFO" in l for l in lines)
    text = (ROOT / "scripts" / "release.sh").read_text()
    assert "LICENSE-MIT LICENSE-APACHE" in text


def test_tool_versions_pinned(release_repo):
    lines = planned(release_repo, "--snapshot")
    env = release_env()
    auditable = lines[index_of(lines, "cargo auditable build")]
    cyclonedx = lines[index_of(lines, "cargo cyclonedx")]
    assert "cargo install --locked cargo-auditable@" + env["CARGO_AUDITABLE_VERSION"] in auditable
    assert "cargo install --locked cargo-cyclonedx@" + env["CARGO_CYCLONEDX_VERSION"] in cyclonedx
    assert "exit 1" in auditable and "exit 1" in cyclonedx


def test_dry_run_publish_plan_verifies_tag_and_drafts(release_repo):
    git(release_repo, "tag", "v1.2.3")
    lines = planned(release_repo, "--tag", "--publish")
    ls = index_of(lines, "git ls-remote --tags origin")
    assert ls < index_of(lines, "--platform")
    create = lines[index_of(lines, "gh release create v1.2.3")]
    assert "--draft" in create and "--verify-tag" in create
    assert index_of(lines, "gh release edit v1.2.3", "--draft=false") > index_of(lines, "gh release create")


def test_publish_refuses_when_origin_lacks_tag(release_repo, tmp_path):
    git(release_repo, "tag", "v1.2.3")
    bare = tmp_path / "origin.git"
    subprocess.run(["git", "init", "-q", "--bare", str(bare)], check=True)
    git(release_repo, "remote", "add", "origin", str(bare))
    result = release(release_repo, "--tag", "--publish")
    assert result.returncode == 2
    assert "v1.2.3" in result.stderr and "origin" in result.stderr


def test_publish_refuses_when_origin_tag_is_another_commit(release_repo, tmp_path):
    git(release_repo, "tag", "v1.2.3")
    bare = tmp_path / "origin.git"
    subprocess.run(["git", "init", "-q", "--bare", str(bare)], check=True)
    git(release_repo, "remote", "add", "origin", str(bare))
    git(release_repo, "push", "-q", "origin", "HEAD~1:refs/tags/v1.2.3")
    result = release(release_repo, "--tag", "--publish")
    assert result.returncode == 2
    assert "not" in result.stderr


def test_empty_release_notes_refused(release_repo):
    (release_repo / "CHANGELOG.md").write_text("# Changelog\n\n## 1.2.3\n\n## 1.2.2\n\nOld.\n")
    git(release_repo, "commit", "-q", "-am", "empty notes")
    git(release_repo, "tag", "v1.2.3")
    result = release(release_repo, "--tag", "--publish", "--dry-run")
    assert result.returncode == 2
    assert "notes" in result.stderr


def test_untracked_file_under_crates_refused(release_repo):
    (release_repo / "crates" / "core").mkdir(parents=True)
    (release_repo / "crates" / "core" / "stray.rs").write_text("x\n")
    result = release(release_repo, "--snapshot", "--dry-run")
    assert result.returncode == 2
    assert "stray.rs" in result.stderr


def test_stale_sums_removed_at_start(release_repo, tmp_path):
    out, env = stub_sbom(release_repo, tmp_path, '{"metadata": {}, "p": "%s"}' % release_repo)
    (out / "SHA256SUMS").write_text("old\n")
    (out / "BUILDINFO").write_text("old\n")
    result = release(release_repo, "--snapshot", env=env)
    assert result.returncode == 1
    assert not (out / "SHA256SUMS").exists()
    assert (out / "BUILDINFO").read_text() != "old\n"

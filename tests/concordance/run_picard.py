"""Write Picard reference metrics for every case in cases.py.

Usage: python tests/concordance/run_picard.py [--force]

Output: tests/corpus/out/picard/<case>.<metric>.txt (and <case>.<metric>.per_base.txt
for PER_BASE_COVERAGE cases), plus manifest.json that
records, per case, the Picard command and whether Picard wrote a file
(CollectInsertSizeMetrics writes none when no orientation is kept). A case
is rerun only when its BAM is newer than the manifest entry or --force is
given. Picard must be the pinned PICARD_VERSION.

CollectInsertSizeMetrics refuses to run without Rscript, so a no-op Rscript
in target/fakeR goes first on PATH; it is never asked to draw anything
useful because H=/dev/null.
"""
import json
import os
import subprocess
import sys

from cases import CASES, METRICS, PICARD_DIR, PICARD_VERSION, ROOT

FAKE_R = ROOT / "target" / "fakeR"


def picard_env():
    FAKE_R.mkdir(parents=True, exist_ok=True)
    rscript = FAKE_R / "Rscript"
    rscript.write_text("#!/bin/sh\nexit 0\n")
    rscript.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{FAKE_R}{os.pathsep}{env['PATH']}"
    return env


def picard_version(env):
    out = subprocess.run(["picard", "CollectInsertSizeMetrics", "--version"],
                         capture_output=True, text=True, env=env)
    text = out.stdout + out.stderr
    for line in text.splitlines():
        if line.startswith("Version:"):
            return line.split(":", 1)[1].strip()
    raise SystemExit(f"cannot read the Picard version from: {text!r}")


def main():
    force = "--force" in sys.argv[1:]
    env = picard_env()
    version = picard_version(env)
    if version != PICARD_VERSION:
        raise SystemExit(f"Picard {version} found, the references are pinned to {PICARD_VERSION}")
    PICARD_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = PICARD_DIR / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    for case in CASES:
        metric = METRICS[case.metric]
        cmd = ["picard", metric.picard_tool, f"I={case.bam}", f"O={case.picard_txt}",
               *metric.picard_args, *case.picard_args, "VERBOSITY=ERROR", "QUIET=true"]
        if case.per_base:
            cmd.append(f"PER_BASE_COVERAGE={case.picard_per_base}")
        entry = manifest.get(case.id)
        if (not force and entry and entry["command"] == cmd[1:]
                and entry["bam_mtime"] == case.bam.stat().st_mtime_ns
                and (case.picard_txt.exists() or not entry["written"])
                and (not case.per_base or case.picard_per_base.exists())):
            continue
        case.picard_txt.unlink(missing_ok=True)
        if case.per_base:
            case.picard_per_base.unlink(missing_ok=True)
        subprocess.run(cmd, check=True, env=env, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        manifest[case.id] = {
            "command": cmd[1:],
            "picard_version": version,
            "bam_mtime": case.bam.stat().st_mtime_ns,
            "written": case.picard_txt.exists(),
        }
        print(f"{case.id}: {'written' if case.picard_txt.exists() else 'no file (as Picard does)'}")
        manifest_path.write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    print(f"{len(CASES)} cases, Picard {version}, references in {PICARD_DIR}")


if __name__ == "__main__":
    main()

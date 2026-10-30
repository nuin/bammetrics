"""Time Picard and bammetrics on the benchmark BAM (benches/make_bam.py).

Usage:
    uv run --no-project python benches/bench.py [--reps N] [--bench-dir DIR] [--coverage-only]

Runs, each under /usr/bin/time -l (wall seconds, peak RSS):
  picard      CollectQualityYieldMetrics, CollectAlignmentSummaryMetrics R=,
              CollectInsertSizeMetrics, CollectHsMetrics R= twice (baits and
              targets, then targets as both), one JVM per tool as the
              clinical pipeline runs them
  bammetrics  the same five runs, one process each
  collect     bammetrics collect: all five metrics in one pass, at the
              default thread count and at --threads 1
  coverage    the clinical pipeline's coverage command (CollectHsMetrics with
              the targets as baits, R=, MINIMUM_MAPPING_QUALITY=0,
              MINIMUM_BASE_QUALITY=0 and PER_BASE_COVERAGE) in Picard syntax,
              Picard and bammetrics; --coverage-only runs just this

Each bammetrics file is compared with Picard's from `## METRICS CLASS` on;
the table says whether the text is identical. Wall time is the median of
--reps runs (default 3); peak RSS is the largest. Results go to stdout as a
Markdown table and to DIR/results.json.
"""
import argparse
import json
import os
import platform
import re
import statistics
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"
FAKE_R = ROOT / "target" / "fakeR"


def picard_env():
    """CollectInsertSizeMetrics needs an Rscript on PATH; a no-op one is enough."""
    FAKE_R.mkdir(parents=True, exist_ok=True)
    rscript = FAKE_R / "Rscript"
    rscript.write_text("#!/bin/sh\nexit 0\n")
    rscript.chmod(0o755)
    env = dict(os.environ)
    env["PATH"] = f"{FAKE_R}{os.pathsep}{env['PATH']}"
    return env


def timed(cmd, env):
    """Wall seconds and peak RSS in MB of one run, from /usr/bin/time -l."""
    out = subprocess.run(["/usr/bin/time", "-l", *cmd], env=env, capture_output=True, text=True)
    if out.returncode:
        raise SystemExit(f"failed: {' '.join(cmd)}\n{out.stderr[-2000:]}")
    wall = float(re.search(r"([\d.]+) real", out.stderr).group(1))
    rss = int(re.search(r"(\d+)\s+maximum resident set size", out.stderr).group(1))
    return wall, rss / 1e6


def measure(cmd, env, reps):
    runs = [timed(cmd, env) for _ in range(reps)]
    return statistics.median(w for w, _ in runs), max(r for _, r in runs)


def body(path):
    text = Path(path).read_text()
    return text[text.find("## METRICS CLASS"):]


def versions():
    def first_line(cmd):
        out = subprocess.run(cmd, capture_output=True, text=True)
        return (out.stdout + out.stderr).strip().splitlines()[0]

    model = first_line(["sysctl", "-n", "hw.model"])
    cores = first_line(["sysctl", "-n", "hw.ncpu"])
    perf = first_line(["sysctl", "-n", "hw.perflevel0.physicalcpu"])
    eff = first_line(["sysctl", "-n", "hw.perflevel1.physicalcpu"])
    memory = int(first_line(["sysctl", "-n", "hw.memsize"])) // 2**30
    return {
        "machine": f"{model}, {cores} cores ({perf} performance, {eff} efficiency), {memory} GB",
        "os": f"macOS {platform.mac_ver()[0]}",
        "bammetrics": first_line([str(BAMMETRICS), "--version"]),
        "rustc": first_line(["rustc", "--version"]),
        "picard": first_line(["picard", "CollectHsMetrics", "--version"]),
        "java": first_line(["java", "-version"]) if _has_java() else _picard_java(),
    }


def _has_java():
    return subprocess.run(["java", "-version"], capture_output=True).returncode == 0


def _picard_java():
    """The JVM the Homebrew picard wrapper runs (its JAVA_HOME default)."""
    wrapper = Path(subprocess.run(["which", "picard"], capture_output=True,
                                  text=True).stdout.strip()).read_text()
    home = re.search(r'JAVA_HOME:-([^}"]+)', wrapper).group(1)
    out = subprocess.run([f"{home}/bin/java", "-version"], capture_output=True, text=True)
    return out.stderr.strip().splitlines()[0]


def coverage_run(d, out, env, reps):
    """The pipeline's coverage step, argument for argument, with `bammetrics`
    in place of `picard`. Identical means the summary from `## METRICS CLASS`
    on and the whole per-base file."""
    targets, ref, bam = d / "targets.interval_list", d / "reference.fa", d / "bench.bam"
    walls = {}
    for tool, exe in (("picard", "picard"), ("bammetrics", str(BAMMETRICS))):
        cmd = [exe, "CollectHsMetrics", f"BI={targets}", f"I={bam}",
               f"PER_BASE_COVERAGE={out / f'{tool}.nucl.partial.out'}",
               "MINIMUM_MAPPING_QUALITY=0", "MINIMUM_BASE_QUALITY=0",
               f"TARGET_INTERVALS={targets}", f"OUTPUT={out / f'{tool}.summary.partial.out'}",
               f"R={ref}", "QUIET=true", "VALIDATION_STRINGENCY=SILENT", "VERBOSITY=ERROR"]
        walls[tool] = measure(cmd, env, reps)
        print(f"coverage command {tool}: {walls[tool][0]:.2f} s", flush=True)
    per_base = [(out / f"{t}.nucl.partial.out").read_bytes() for t in ("picard", "bammetrics")]
    summary = [body(out / f"{t}.summary.partial.out") for t in ("picard", "bammetrics")]
    return {"picard_s": walls["picard"][0], "picard_mb": walls["picard"][1],
            "bammetrics_s": walls["bammetrics"][0], "bammetrics_mb": walls["bammetrics"][1],
            "per_base_lines": per_base[0].count(b"\n"),
            "identical": per_base[0] == per_base[1] and summary[0] == summary[1]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--bench-dir", type=Path, default=ROOT / "target" / "bench")
    parser.add_argument("--coverage-only", action="store_true")
    args = parser.parse_args()
    d = args.bench_dir
    bam, ref = d / "bench.bam", d / "reference.fa"
    baits, targets = d / "baits.interval_list", d / "targets.interval_list"
    out = d / "out"
    out.mkdir(exist_ok=True)
    env = picard_env()
    silent = ["VALIDATION_STRINGENCY=SILENT", "VERBOSITY=ERROR", "QUIET=true"]

    # (label, Picard tool and arguments, bammetrics subcommand and arguments)
    tools = [
        ("CollectQualityYieldMetrics", ["CollectQualityYieldMetrics"],
         ["quality-yield"], ".quality_yield_metrics"),
        ("CollectAlignmentSummaryMetrics R=", ["CollectAlignmentSummaryMetrics", f"R={ref}"],
         ["alignment-summary", "-R", str(ref)], ".alignment_summary_metrics"),
        ("CollectInsertSizeMetrics", ["CollectInsertSizeMetrics", "H=/dev/null"],
         ["insert-size"], ".insert_size_metrics"),
        ("CollectHsMetrics R= (baits, targets)",
         ["CollectHsMetrics", f"R={ref}", f"BAIT_INTERVALS={baits}", f"TARGET_INTERVALS={targets}"],
         ["hs-metrics", "-R", str(ref), "--bait-intervals", str(baits),
          "--target-intervals", str(targets)], ".baits.hs_metrics"),
        ("CollectHsMetrics R= (targets as baits)",
         ["CollectHsMetrics", f"R={ref}", f"BAIT_INTERVALS={targets}",
          f"TARGET_INTERVALS={targets}"],
         ["hs-metrics", "-R", str(ref), "--bait-intervals", str(targets),
          "--target-intervals", str(targets)], ".targets.hs_metrics"),
    ]
    coverage = coverage_run(d, out, env, args.reps)
    if args.coverage_only:
        print(json.dumps(coverage, indent=2))
        return
    rows = []
    for i, (label, picard, ours, suffix) in enumerate(tools):
        picard_txt, ours_txt = out / f"picard.{i}.txt", out / f"ours.{i}.txt"
        p_wall, p_rss = measure(["picard", *picard, f"I={bam}", f"O={picard_txt}", *silent],
                                env, args.reps)
        b_wall, b_rss = measure([str(BAMMETRICS), *ours, str(bam), "-O", str(ours_txt)],
                                env, args.reps)
        rows.append({"tool": label, "suffix": suffix, "picard_s": p_wall, "picard_mb": p_rss,
                     "bammetrics_s": b_wall, "bammetrics_mb": b_rss,
                     "identical": body(picard_txt) == body(ours_txt),
                     "picard_txt": str(picard_txt)})
        print(f"{label}: picard {p_wall:.2f} s, bammetrics {b_wall:.2f} s", flush=True)

    collect_args = ["collect", str(bam), "-R", str(ref),
                    "--hs", str(baits), str(targets), "--hs", str(targets), str(targets)]
    collect = {}
    for label, threads in (("default", []), ("threads_1", ["--threads", "1"])):
        prefix = out / f"collect_{label}"
        collect[label] = measure([str(BAMMETRICS), *collect_args, "-O", str(prefix), *threads],
                                 env, args.reps)
        for row in rows:
            row[f"collect_{label}_identical"] = (
                body(row["picard_txt"]) == body(f"{prefix}{row['suffix']}"))
        print(f"collect {label}: {collect[label][0]:.2f} s", flush=True)

    result = {"versions": versions(), "reps": args.reps, "rows": rows,
              "collect": {k: {"wall_s": w, "peak_mb": m} for k, (w, m) in collect.items()},
              "coverage": coverage,
              "bench": json.loads((d / "bench.json").read_text()),
              "bam_mb": bam.stat().st_size / 1e6}
    (d / "results.json").write_text(json.dumps(result, indent=2) + "\n")

    print("\n| run | Picard wall s | Picard peak MB | bammetrics wall s | bammetrics peak MB "
          "| speed-up | text identical |")
    print("|---|---:|---:|---:|---:|---:|---|")
    for r in rows:
        print(f"| {r['tool']} | {r['picard_s']:.2f} | {r['picard_mb']:.0f} "
              f"| {r['bammetrics_s']:.2f} | {r['bammetrics_mb']:.0f} "
              f"| {r['picard_s'] / r['bammetrics_s']:.1f}x | {r['identical']} |")
    p_total = sum(r["picard_s"] for r in rows)
    b_total = sum(r["bammetrics_s"] for r in rows)
    print(f"| all five, separate runs | {p_total:.2f} | | {b_total:.2f} | "
          f"| {p_total / b_total:.1f}x | |")
    for label, (wall, rss) in collect.items():
        same = all(r[f"collect_{label}_identical"] for r in rows)
        print(f"| collect, {label} | | | {wall:.2f} | {rss:.0f} | {p_total / wall:.1f}x | {same} |")
    c = coverage
    print(f"| coverage command (PER_BASE_COVERAGE) | {c['picard_s']:.2f} | {c['picard_mb']:.0f} "
          f"| {c['bammetrics_s']:.2f} | {c['bammetrics_mb']:.0f} "
          f"| {c['picard_s'] / c['bammetrics_s']:.1f}x | {c['identical']} |")
    print()
    for key, value in result["versions"].items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    main()

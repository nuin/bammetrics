"""T4 scale run for the coverage tools: a 2 Mb reference with 200,000 pairs, Picard
versus bammetrics wall time and peak RSS.

Usage:
    uv run --no-project --with pysam --with numpy python benches/wgs_bench.py [--reps N] [--pairs N]

Builds target/bench_wgs/ (not committed) with the corpus simulator on two contigs
of 1.2 Mb and 0.8 Mb: messy 2x150 pairs (adapters, indels, Q2 tails, RF and
TANDEM pairs, chimeras, split reads, bwa mem -a secondary records, duplicates
marked by Picard MarkDuplicates, QC-fail flags, OQ tags), then runs under
/usr/bin/time -l:

  CollectWgsMetrics                   default pileup
  CollectWgsMetrics USE_FAST_ALGORITHM=true
  CollectRawWgsMetrics                100,000-bin depth histogram and sensitivity
  CollectWgsMetricsWithNonZeroCoverage

Each bammetrics file is compared with Picard's after the StringHeader lines. Wall
time is the median of --reps runs (default 3); peak RSS the largest. Results go
to stdout as a Markdown table and to target/bench_wgs/results.json.
"""
import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "corpus"))
sys.path.insert(0, str(ROOT / "benches"))
import bench  # noqa: E402
import generate  # noqa: E402

SEED = 20261102
CONTIGS = [("chr1", 1_200_000), ("chr2", 800_000)]
TARGETS_PER_CONTIG = {"chr1": 240, "chr2": 160}
DUP_RATE = 0.08


def make_targets(rng):
    targets = []
    for contig, length in CONTIGS:
        pos = 2_000
        for _ in range(TARGETS_PER_CONTIG[contig]):
            size = int(rng.integers(100, 401))
            targets.append((contig, pos, pos + size))
            pos += size + int(rng.integers(800, 4_000))
        assert pos < length - 2_000, contig
    return targets


def build(out, n_pairs):
    out.mkdir(parents=True, exist_ok=True)
    generate.CONTIGS = CONTIGS
    rng = np.random.default_rng(SEED)
    reference = generate.make_reference(rng, out)
    targets = make_targets(rng)
    generate.write_intervals(out, "targets", targets)
    sim = generate.Simulator(rng, reference, targets, messy=True)
    generate.make_paired_tier(sim, out, "bench", n_pairs, dup_rate=DUP_RATE, all_alignments=True)
    records = subprocess.run(["samtools", "view", "-c", str(out / "bench.bam")], capture_output=True,
                             text=True, check=True).stdout.strip()
    summary = {"seed": SEED, "pairs": n_pairs, "contigs": CONTIGS, "records": int(records)}
    (out / "bench.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def body(path):
    lines = Path(path).read_bytes().split(b"\n")
    kept, skip = [], False
    for line in lines:
        if skip:
            skip = False
            continue
        if line == b"## htsjdk.samtools.metrics.StringHeader":
            skip = True
            continue
        kept.append(line)
    return b"\n".join(kept)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=3)
    parser.add_argument("--pairs", type=int, default=200_000)
    parser.add_argument("--bench-dir", type=Path, default=ROOT / "target" / "bench_wgs")
    args = parser.parse_args()
    d = args.bench_dir
    summary = (json.loads((d / "bench.json").read_text()) if (d / "bench.json").exists()
               else build(d, args.pairs))
    bam, ref = d / "bench.bam", d / "reference.fa"
    out = d / "out"
    out.mkdir(exist_ok=True)
    env = bench.picard_env()
    silent = ["VALIDATION_STRINGENCY=SILENT", "VERBOSITY=ERROR", "QUIET=true"]
    runs = [
        ("CollectWgsMetrics", ["CollectWgsMetrics"]),
        ("CollectWgsMetrics USE_FAST_ALGORITHM=true", ["CollectWgsMetrics", "USE_FAST_ALGORITHM=true"]),
        ("CollectRawWgsMetrics", ["CollectRawWgsMetrics"]),
        ("CollectWgsMetricsWithNonZeroCoverage", ["CollectWgsMetricsWithNonZeroCoverage"]),
    ]
    rows = []
    for label, (tool, *options) in runs:
        files = {}
        walls = {}
        for side, exe, extra in (("picard", "picard", silent), ("bammetrics", str(bench.BAMMETRICS), [])):
            files[side] = out / f"{side}.{tool}{'.fast' if options else ''}.txt"
            cmd = [exe, tool, f"I={bam}", f"O={files[side]}", f"R={ref}", *options, *extra]
            walls[side] = bench.measure(cmd, env, args.reps)
            print(f"{label} {side}: {walls[side][0]:.2f} s, {walls[side][1]:.0f} MB", flush=True)
        identical = body(files["picard"]) == body(files["bammetrics"])
        rows.append({"run": label, "picard_s": walls["picard"][0], "picard_mb": walls["picard"][1],
                     "bammetrics_s": walls["bammetrics"][0], "bammetrics_mb": walls["bammetrics"][1],
                     "identical": identical})
    results = {"setup": bench.versions(), "bam": summary, "reps": args.reps, "rows": rows}
    (d / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print()
    print("| run | Picard wall s | Picard peak MB | bammetrics wall s | bammetrics peak MB | speed-up | text identical |")
    print("|---|---:|---:|---:|---:|---:|---|")
    for r in rows:
        print(f"| {r['run']} | {r['picard_s']:.2f} | {r['picard_mb']:.0f} | {r['bammetrics_s']:.2f} | "
              f"{r['bammetrics_mb']:.0f} | {r['picard_s'] / r['bammetrics_s']:.1f}x | "
              f"{'yes' if r['identical'] else 'NO'} |")


if __name__ == "__main__":
    main()

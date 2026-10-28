"""Generate the benchmark BAM under target/bench/ (not committed).

Usage:
    uv run --no-project --with pysam --with numpy python benches/make_bam.py [N_PAIRS] [OUT_DIR]

Reuses the corpus simulator (tests/corpus/generate.py) on a larger genome:
two 6 and 4 Mb contigs plus chrM, 2000 targets of 100-400 bp, baits padded
by 50 bp, 70% of fragments on target. Reads are the messy t2 kind (adapters,
indels, Q2 tails, RF/TANDEM, chimeras, split reads, secondary alignments
from bwa mem -a), aligned with bwa, duplicates marked by Picard
MarkDuplicates, QC-fail flags and OQ tags added. 2 million pairs by default.
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "corpus"))
import generate  # noqa: E402

SEED = 20261029
CONTIGS = [("chr1", 6_000_000), ("chr2", 4_000_000), ("chrM", 16_000)]
TARGETS_PER_CONTIG = {"chr1": 1200, "chr2": 800}
DUP_RATE = 0.08


def make_targets(rng):
    """Sorted, non-overlapping 100-400 bp targets, gaps 800-4000 bp, and every
    tenth gap under 2 * BAIT_PAD so that neighbouring baits merge."""
    targets = []
    for contig, length in CONTIGS:
        count = TARGETS_PER_CONTIG.get(contig, 0)
        pos = 2_000
        for i in range(count):
            size = int(rng.integers(100, 401))
            targets.append((contig, pos, pos + size))
            gap = int(rng.integers(20, 90)) if i % 10 == 9 else int(rng.integers(800, 4_000))
            pos += size + gap
        assert pos < length - 2_000, contig
    return targets


def main():
    n_pairs = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "target" / "bench"
    out.mkdir(parents=True, exist_ok=True)
    # The simulator and make_reference read the module-level contig list.
    generate.CONTIGS = CONTIGS
    rng = np.random.default_rng(SEED)
    reference = generate.make_reference(rng, out)
    targets = make_targets(rng)
    generate.write_intervals(out, "targets", targets)
    pad = generate.BAIT_PAD
    baits = generate.merge([(c, max(s - pad, 0), e + pad) for c, s, e in targets])
    generate.write_intervals(out, "baits", baits)
    sim = generate.Simulator(rng, reference, targets, messy=True)
    generate.make_paired_tier(sim, out, "bench", n_pairs, dup_rate=DUP_RATE, all_alignments=True)
    summary = {"seed": SEED, "pairs": n_pairs, "targets": len(targets), "baits": len(baits),
               "contigs": CONTIGS}
    (out / "bench.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()

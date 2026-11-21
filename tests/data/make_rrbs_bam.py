"""Write tests/data/rrbs.bam: 80 seeded bisulfite reads on chr1 of
tests/data/ref_100k.fa for CollectRrbsMetrics, coordinate sorted, with
three read groups (samples S1 and S2, libraries L1, L2 and L3, platform
units differing from the IDs) and records without RG.

Reads start at the first CG after a seeded position, 3 to 75 bases.
A forward read has C to T at every C except a CpG C kept with probability
0.6; a reverse read has G to A the same way. Every tenth read is not
converted, every ninth has 8 sequencing errors (over MAX_MISMATCH_RATE),
every eleventh is shorter than 5 bases, some bases have quality 2, and one
read is a duplicate, one secondary, one QC-failed and two unmapped.

    python3 tests/data/make_rrbs_bam.py

Needs samtools.
"""
import random
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
REF = DATA / "ref_100k.fa"
BAM = DATA / "rrbs.bam"
GROUPS = [("a", "S1", "L1", "PU.a"), ("b", "S1", "L2", "PU.b"), ("c", "S2", "L3", "PU.c")]


def reference():
    lines = REF.read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">")).upper()


def convert(rng, ref, start, length, reverse, unconverted):
    out = []
    for pos in range(start, start + length):
        base = ref[pos]
        target = "G" if reverse else "C"
        if unconverted or base != target:
            out.append(base)
            continue
        cpg = ref[pos - 1] == "C" if reverse else ref[pos + 1] == "G"
        keep = cpg and rng.random() < 0.6
        out.append(base if keep else ("A" if reverse else "T"))
    return out


def main():
    rng = random.Random(61)
    ref = reference()
    lines = []
    for i in range(80):
        start = ref.index("CG", rng.randint(1000, 1300))
        length = rng.choice((3, 4)) if i % 11 == 0 else rng.choice((30, 50, 75))
        reverse = i % 2 == 1
        seq = convert(rng, ref, start, length, reverse, unconverted=i % 10 == 0)
        if i % 9 == 0:
            for at in rng.sample(range(length), min(8, length)):
                seq[at] = rng.choice([b for b in "ACGT" if b != seq[at]])
        qual = "".join("#" if rng.random() < 0.08 else "I" for _ in range(length))
        flag = 16 if reverse else 0
        flag |= {13: 0x400, 27: 0x100, 41: 0x200}.get(i, 0)
        group = GROUPS[i % 4][0] if i % 4 < 3 else None
        tags = [f"RG:Z:{group}"] if group else []
        lines.append("\t".join([f"r{i:03d}", str(flag), "chr1", str(start + 1), "60",
                                f"{length}M", "*", "0", "0", "".join(seq), qual, *tags]))
    for i in range(2):
        lines.append("\t".join([f"u{i}", "4", "*", "0", "0", "*", "*", "0", "0", "ACGTACGTAC",
                                "IIIIIIIIII"]))
    header = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000"]
    header += [f"@RG\tID:{g}\tSM:{s}\tLB:{lb}\tPU:{pu}" for g, s, lb, pu in GROUPS]
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp) / "rrbs.sam"
        sam.write_text("\n".join(header + lines) + "\n")
        subprocess.run(["samtools", "sort", "--no-PG", "-o", str(BAM), str(sam)], check=True)


if __name__ == "__main__":
    main()

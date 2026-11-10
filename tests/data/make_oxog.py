"""Write tests/data/oxog.bam: 60 seeded read pairs on chr1 of
tests/data/ref_100k.fa for CollectOxoGMetrics, with two read groups (one
sample, libraries L1 and L2) and oxidation damage: about 4% of the reference
G bases of read 1 are read as T and of the reference C bases of read 2 as A,
both counted as sequenced, plus a 0.5% background of other substitutions.
Some bases have quality 2; every fifth read carries an OQ tag.

    python3 tests/data/make_oxog.py

Needs samtools.
"""
import random
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
REF = DATA / "ref_100k.fa"
BAM = DATA / "oxog.bam"
READ_LENGTH = 100
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reference():
    lines = REF.read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">")).upper()


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def damage(rng, as_read, first):
    """Oxidation on the read as sequenced: G to T in read 1, C to A in read 2."""
    out = []
    for base in as_read:
        if first and base == "G" and rng.random() < 0.04:
            base = "T"
        elif not first and base == "C" and rng.random() < 0.04:
            base = "A"
        elif rng.random() < 0.005:
            base = rng.choice([b for b in "ACGT" if b != base])
        out.append(base)
    return "".join(out)


def qualities(rng):
    return "".join("#" if rng.random() < 0.05 else "I" for _ in range(READ_LENGTH))


def records(rng, ref):
    lines = []
    for i in range(60):
        start = 1000 + rng.randrange(0, 2000)
        insert = rng.randrange(150, 400)
        mate = start + insert - READ_LENGTH
        first_reverse = rng.random() < 0.5
        group = "rgA" if i % 2 == 0 else "rgB"
        starts = (start, mate) if not first_reverse else (mate, start)
        for first in (True, False):
            pos = starts[0] if first else starts[1]
            reverse = first_reverse if first else not first_reverse
            segment = ref[pos - 1:pos - 1 + READ_LENGTH]
            as_read = revcomp(segment) if reverse else segment
            read = damage(rng, as_read, first)
            seq = revcomp(read) if reverse else read
            flag = 0x1 | 0x2 | (0x40 if first else 0x80)
            flag |= 0x10 if reverse else 0
            flag |= 0x20 if (not reverse) else 0
            mate_pos = starts[1] if first else starts[0]
            tlen = insert if pos == min(starts) else -insert
            tags = [f"RG:Z:{group}", f"MC:Z:{READ_LENGTH}M"]
            if (2 * i + first) % 5 == 0:
                tags.append("OQ:Z:" + qualities(rng))
            lines.append("\t".join([f"p{i}", str(flag), "chr1", str(pos), "60", f"{READ_LENGTH}M",
                                    "=", str(mate_pos), str(tlen), seq, qualities(rng), *tags]))
    return lines


def main():
    rng = random.Random(8)
    header = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000",
              "@RG\tID:rgA\tSM:S1\tLB:L1", "@RG\tID:rgB\tSM:S1\tLB:L2"]
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp, "oxog.sam")
        sam.write_text("\n".join(header + records(rng, reference())) + "\n")
        subprocess.run(["samtools", "sort", "--no-PG", "-o", str(BAM), str(sam)], check=True)
    print(f"{BAM.name}: {BAM.stat().st_size} bytes")


if __name__ == "__main__":
    main()

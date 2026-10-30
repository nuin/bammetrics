"""Small hand-built BAMs for insert-size edge cases: single pair, even and odd
counts, mode ties, one orientation only, duplicates, missing MC, filtered
records, an empty file. Built from SAM text with samtools."""
import random
import subprocess

L = 50
HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:100000\n@SQ\tSN:c2\tLN:100000\n"
RG_A = "@RG\tID:a\tSM:s\tLB:l\tPU:p\n"
RG_B = "@RG\tID:b\tSM:s2\tLB:l2\tPU:q\n"


def pair(name, pos, mpos, tlen, f1, f2, ref="c1", mref="=", mc=True, rg="a", mapq=60):
    tags = f"RG:Z:{rg}" + (f"\tMC:Z:{L}M" if mc else "")
    seq, qual = "A" * L, "I" * L
    return [f"{name}\t{f1}\t{ref}\t{pos}\t{mapq}\t{L}M\t{mref}\t{mpos}\t{tlen}\t{seq}\t{qual}\t{tags}",
            f"{name}\t{f2}\t{ref}\t{mpos}\t{mapq}\t{L}M\t{mref}\t{pos}\t{-tlen}\t{seq}\t{qual}\t{tags}"]


def fr(name, pos, insert, extra=0, mc=True):
    return pair(name, pos, pos + insert - L, insert, 97 | extra, 145 | extra, mc=mc)


def fr_set(prefix, inserts, extra=0):
    return sum([fr(f"{prefix}{i}", 100 + i * 1000, n, extra) for i, n in enumerate(inserts)], [])


def _filters():
    recs = fr("ok1", 100, 200) + fr("ok2", 1100, 220)
    recs += fr("sec", 2100, 500, 256) + fr("sup", 3100, 500, 2048) + fr("qcf", 4100, 500, 512)
    recs += pair("mq0", 5100, 5350, 333, 97, 145, mapq=0)
    recs += pair("ic", 6100, 6200, 0, 97, 145, mref="c2")
    recs += pair("z", 7100, 7150, 0, 97, 145)
    # Mate unmapped: the mate is placed at the mapped read's position with no CIGAR.
    recs += [f"mu\t105\tc1\t8100\t60\t{L}M\t=\t8100\t0\t{'A' * L}\t{'I' * L}\tRG:Z:a\tMC:Z:*",
             f"mu\t133\tc1\t8100\t0\t*\t=\t8100\t0\t{'A' * L}\t{'I' * L}\tRG:Z:a"]
    return recs


def _tlen_mismatch():
    recs = fr("m1", 100, 300)
    recs[0] = recs[0].replace("\t300\t", "\t290\t", 1)
    return recs + fr("m2", 1000, 300)


def _wide():
    rng = random.Random(1)
    sizes = [150] * 30 + [160] * 20 + [170] * 10 + [800, 5000, 20000, 40000]
    return sum([fr(f"w{i}", 100 + i * 100, rng.choice(sizes)) for i in range(60)], [])


def _orientation(prefix, count, flags, gap):
    return sum([pair(f"{prefix}{i}", 100 + i * 1000, 100 + i * 1000 + gap + i, 250 + i, *flags)
                for i in range(count)], [])


# name -> (read group header lines, records, whether Picard writes a metrics file with defaults)
CASES = {
    "single": (RG_A, fr("p1", 100, 300), True),
    "all_duplicates": (RG_A, fr_set("d", [200, 201, 202, 203, 204], 1024), False),
    "mode_tie": (RG_A, fr_set("t", [200, 200, 300, 300, 400, 250]), True),
    "even_count": (RG_A, fr_set("e", [100, 110, 120, 131]), True),
    "odd_count": (RG_A, fr_set("o", [100, 110, 121, 135, 500]), True),
    "only_tandem": (RG_A, _orientation("x", 4, (65, 129), 200), True),
    "only_tandem_reverse": (RG_A, _orientation("y", 4, (81, 145), 200), True),
    "only_rf": (RG_A, _orientation("r", 3, (113, 177), 200), True),
    "no_mc": (RG_A, sum([fr(f"n{i}", 100 + i * 1000, 150 + i * 7, mc=False) for i in range(8)], []), True),
    "dovetail_no_mc": (RG_A, pair("dv", 200, 150, 100, 97, 145, mc=False)
                       + pair("dv2", 5000, 4980, 70, 97, 145, mc=False), True),
    "dovetail_mc": (RG_A, pair("dv", 200, 150, 100, 97, 145) + pair("dv2", 5000, 4980, 70, 97, 145), True),
    "filtered_records": (RG_A, _filters(), True),
    "tlen_mismatch": (RG_A, _tlen_mismatch(), True),
    "wide_spread": (RG_A, _wide(), True),
    "two_read_groups": (RG_A + RG_B, sum([pair(f"g{i}", 100 + i * 500, 250 + i * 500 + i, 200 + i, 97, 145,
                                               rg="a" if i % 2 else "b") for i in range(12)], []), True),
    "tiny_inserts": (RG_A, pair("tn", 100, 100, 50, 97, 145) + pair("tn2", 300, 300, 50, 97, 145), True),
    "reverse_first": (RG_A, sum([pair(f"v{i}", 100 + i * 500, 250 + i * 500, 200, 81, 161)
                                 for i in range(5)], []), True),
    "unpaired_only": (RG_A, [f"u1\t0\tc1\t100\t60\t{L}M\t*\t0\t0\t{'A' * L}\t{'I' * L}\tRG:Z:a"], False),
    "empty": (RG_A, [], False),
}


def build(out_dir, name):
    """Write and index `name`.bam under out_dir; return its path."""
    rg, records, _ = CASES[name]
    records = sorted(records, key=lambda r: (r.split("\t")[2], int(r.split("\t")[3])))
    sam = out_dir / f"{name}.sam"
    sam.write_text(HEADER + rg + "".join(r + "\n" for r in records))
    bam = out_dir / f"{name}.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam

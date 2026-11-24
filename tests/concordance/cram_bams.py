"""Inputs of test_cram_writer.py: SAM text written line by line, BAMs through samtools.

    python3 tests/concordance/cram_bams.py OUT

`reference(directory)` writes `cram.fa` (+ `.fai`, `.dict`): `e1` (600 bases, 101 to 140 lower
case, every IUPAC code from 201, an N run from 301), `e2` (150 bases, for reads past a contig
end), `e3` (3,000 bases) and `e4` (40,000 bases, for a CIGAR of 70,000 operators).

`edges(directory)`: one coordinate-sorted BAM whose records each store something a CRAM
holds in a way of its own. Bases: a match over lower-case reference, ACGT mismatches, read
N, read IUPAC codes and `=`, reads over the IUPAC reference bases (equal and not) and over
the N run, a read 20 bases past the end of e2. CIGAR: one- and five-base insertions, a
deletion, a skip, soft and hard clips, padding, an insertion first, 70,000 operators (a `CG`
tag in the BAM). No qualities, with and without a mismatch; no bases (a secondary record
with its CIGAR). Tags: every type and every `B` subtype, an empty string, `RG` first, last,
absent and naming a read group the header lacks, `MD` and `NM` right, wrong and absent.
Mates: a proper pair, a pair whose TLEN is not the span, a pair at one start, a mapped read
with an unmapped mate placed on it, a mate on another contig. Flags: secondary, supplementary, duplicate, QC fail, reverse. Names: 254
characters, punctuation, one name three times. MAPQ 0 and 255. Unplaced records last, one
without bases.

`no_hd(directory)` is the same records under a header without `@HD`; `empty(directory)` a
header and no record; `slices(directory)` 14,320 records that make every kind of slice (300
on e1 and 20 on e2, too few to end a slice; 12,000 on e3; 500 on e4; 1,500 unplaced);
`shuffled(directory)` those records in random order under `SO:unsorted`, so every slice
holds several references. `overhang(directory)`: 28 records at the end of e2, each with bases
past the last base of the contig: forward and reverse, with and without qualities, the
bases past the end N, ACGT, IUPAC codes, `=` or a mix, an insertion and a soft clip past the
end, a deletion over the end, a read whose first base is the last base of the contig, and a
read that ends on it. `overhang_clean(directory)` is the 27 of them CleanSam takes.
`overhang_shuffled(directory)` is those records among 1,200 others in
random order, so they fall in multi-reference slices, where the writer reads the bases
under each record through the `.fai`.

`limits(directory)` holds what a CRAM does not give back, whoever writes it: `=` and `X`
operators and neighbouring `M` operators (read back as one `M`), a CIGAR and a mapping
quality on an unmapped record, and an unpaired record with mate fields set (htslib 1.25
reads its mate reference as `*`, htsjdk 5.0.0 its mate reference and start as `*` and 0).
samtools' own CRAM of that file is the truth for it. A mapped record with CIGAR `*` is not
here: `samtools view` reads it from SAM text as unmapped.
"""
import random
import subprocess
import sys
from pathlib import Path

IUPAC = "RYKMSWBDHVN"
SEED = 20261123
PUNCTUATION = "!#$%&'()+,-./:;<=>?[]^_{|}~"


def contigs():
    rng = random.Random(SEED)
    e1 = [rng.choice("ACGT") for _ in range(600)]
    for at in range(100, 140):
        e1[at] = e1[at].lower()
    for offset, code in enumerate(IUPAC):
        e1[200 + offset] = code
    for at in range(300, 320):
        e1[at] = "N"
    return [("e1", "".join(e1))] + [
        (name, "".join(rng.choice("ACGT") for _ in range(length)))
        for name, length in (("e2", 150), ("e3", 3000), ("e4", 40000))]


def reference(directory):
    path = Path(directory) / "cram.fa"
    if not path.exists():
        with open(path, "w") as out:
            for name, bases in contigs():
                out.write(f">{name}\n")
                for at in range(0, len(bases), 60):
                    out.write(bases[at:at + 60] + "\n")
        subprocess.run(["samtools", "faidx", str(path)], check=True)
        subprocess.run(["samtools", "dict", "-o", str(path.with_suffix(".dict")), str(path)],
                       check=True, capture_output=True)
    return path


def header(order="coordinate", hd=True):
    lines = [f"@HD\tVN:1.6\tSO:{order}"] if hd else []
    lines += [f"@SQ\tSN:{name}\tLN:{len(bases)}" for name, bases in contigs()]
    lines += ["@RG\tID:g1\tSM:s1\tLB:l1\tPL:ILLUMINA", "@RG\tID:g2\tSM:s1\tLB:l2\tPL:ILLUMINA",
              "@PG\tID:aligner\tPN:aligner\tVN:1", "@CO\ta comment"]
    return "".join(line + "\n" for line in lines)


def to_bam(sam, bam):
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


class Records:
    """SAM lines kept with their place, written in coordinate order, unplaced last."""

    def __init__(self):
        self.contigs = dict(contigs())
        self.order = {name: at for at, name in enumerate(self.contigs)}
        self.lines = []
        self.rng = random.Random(SEED + 1)

    def bases(self, contig, pos, length):
        return self.contigs[contig][pos - 1:pos - 1 + length].upper()

    def quals(self, length):
        return "".join(chr(33 + self.rng.randrange(2, 41)) for _ in range(length))

    def add(self, name, flag, contig, pos, cigar, seq, qual=None, tags=(), mapq=60,
            mate=("*", 0, 0)):
        if qual is None:
            qual = "*" if seq == "*" else self.quals(len(seq))
        fields = [name, str(flag), contig, str(pos), str(mapq), cigar, mate[0], str(mate[1]),
                  str(mate[2]), seq, qual, *tags]
        key = (self.order.get(contig, len(self.order)), pos, len(self.lines))
        self.lines.append((key, "\t".join(fields)))

    def text(self, order="coordinate", hd=True, shuffle=None):
        lines = [line for _, line in sorted(self.lines)]
        if shuffle is not None:
            random.Random(shuffle).shuffle(lines)
        return header(order, hd) + "".join(line + "\n" for line in lines)


def replace(bases, changes):
    bases = list(bases)
    for at, base in changes.items():
        bases[at] = base
    return "".join(bases)


def other(base):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}.get(base, "A")


def edge_records():
    r = Records()
    b = r.bases
    r.add("match_lower", 0, "e1", 95, "60M", b("e1", 95, 60), tags=["RG:Z:g1", "NM:i:0", "MD:Z:60"])
    seq = b("e1", 20, 50)
    r.add("mm_acgt", 0, "e1", 20, "50M",
          replace(seq, {at: other(seq[at]) for at in (0, 7, 8, 30, 49)}), tags=["RG:Z:g2"])
    r.add("read_n", 16, "e1", 30, "40M", replace(b("e1", 30, 40), {2: "N", 20: "N", 39: "N"}))
    r.add("read_iupac", 0, "e1", 40, "40M",
          replace(b("e1", 40, 40), {1: "R", 5: "Y", 6: "M", 30: "K", 39: "B"}))
    r.add("read_eq", 0, "e1", 50, "30M", replace(b("e1", 50, 30), {0: "=", 10: "=", 29: "="}))
    r.add("over_iupac_same", 0, "e1", 195, "30M", b("e1", 195, 30))
    r.add("over_iupac_diff", 0, "e1", 195, "30M",
          replace(b("e1", 195, 30), {5: "A", 6: "C", 7: "R", 8: "N", 15: "G"}))
    r.add("over_n_run", 0, "e1", 295, "40M",
          replace(b("e1", 295, 40), {at: "ACGT"[at % 4] for at in range(5, 20)}))
    r.add("n_on_n", 0, "e1", 300, "30M", b("e1", 300, 30))
    r.add("ins1", 0, "e1", 400, "10M1I20M", b("e1", 400, 10) + "T" + b("e1", 410, 20))
    r.add("ins5", 0, "e1", 402, "10M5I20M", b("e1", 402, 10) + "ACGTN" + b("e1", 412, 20))
    r.add("del", 0, "e1", 404, "15M3D15M", b("e1", 404, 15) + b("e1", 422, 15))
    r.add("skip", 0, "e1", 406, "10M100N10M", b("e1", 406, 10) + b("e1", 516, 10))
    r.add("softs", 0, "e1", 408, "5S20M6S", "ACGTA" + b("e1", 408, 20) + "TTGGCC")
    r.add("hards", 0, "e1", 410, "3H20M4H", b("e1", 410, 20))
    r.add("pad", 0, "e1", 412, "10M1I2P1I10M", b("e1", 412, 10) + "GC" + b("e1", 422, 10))
    r.add("soft_hard", 16, "e1", 414, "2H3S10M4S5H", "NNN" + b("e1", 414, 10) + "ACGT")
    r.add("ins_first", 0, "e1", 416, "3I20M", "TTT" + b("e1", 416, 20))
    r.add("del_then_ins", 0, "e1", 418, "10M2D2I10M", b("e1", 418, 10) + "AA" + b("e1", 430, 10))
    r.add("noqual", 0, "e1", 430, "30M", b("e1", 430, 30), qual="*")
    seq = b("e1", 432, 30)
    r.add("noqual_mm", 0, "e1", 432, "30M", replace(seq, {3: other(seq[3]), 9: "R", 12: "N"}),
          qual="*")
    r.add("noseq", 256, "e1", 434, "3S10M2I5M1H", "*", mapq=3)
    r.add("tags_all", 0, "e1", 440, "25M", b("e1", 440, 25), tags=[
        "XA:A:x", "Xc:i:-5", "XC:i:200", "Xs:i:-3000", "XS:i:60000", "Xi:i:-100000",
        "XI:i:3000000000", "Xf:f:3.5", "XZ:Z:hello world", "XE:Z:", "XH:H:1AE301",
        "Bc:B:c,-1,2", "BC:B:C,1,200", "Bs:B:s,-300,5", "BS:B:S,60000", "Bi:B:i,-70000",
        "BI:B:I,3000000000", "Bf:B:f,1.5,-2"])
    r.add("rg_first", 0, "e1", 442, "25M", b("e1", 442, 25), tags=["RG:Z:g1", "XX:i:1", "YY:Z:y"])
    r.add("rg_middle", 0, "e1", 444, "25M", b("e1", 444, 25), tags=["XX:i:1", "RG:Z:g2", "YY:Z:y"])
    r.add("rg_unknown", 0, "e1", 446, "25M", b("e1", 446, 25), tags=["XX:i:1", "RG:Z:nope"])
    r.add("mdnm_right", 0, "e1", 448, "25M", b("e1", 448, 25), tags=["NM:i:0", "MD:Z:25", "XX:i:2"])
    r.add("mdnm_wrong", 0, "e1", 450, "25M", b("e1", 450, 25), tags=["MD:Z:3A21", "NM:i:9"])
    r.add("md_only", 0, "e1", 452, "25M", b("e1", 452, 25), tags=["MD:Z:25"])
    r.add("pair_ok", 99, "e1", 460, "30M", b("e1", 460, 30), mate=("=", 520, 90),
          tags=["MC:Z:30M", "RG:Z:g1"])
    r.add("pair_ok", 147, "e1", 520, "30M", b("e1", 520, 30), mate=("=", 460, -90),
          tags=["MC:Z:30M", "RG:Z:g1"])
    r.add("pair_tlen_off", 99, "e1", 462, "30M", b("e1", 462, 30), mate=("=", 522, 5000))
    r.add("pair_tlen_off", 147, "e1", 522, "30M", b("e1", 522, 30), mate=("=", 462, -1))
    r.add("pair_one_start", 163, "e1", 464, "30M", b("e1", 464, 30), mate=("=", 464, 30))
    r.add("pair_one_start", 83, "e1", 464, "30M", b("e1", 464, 30), mate=("=", 464, -30))
    r.add("mate_unmapped", 73, "e1", 466, "30M", b("e1", 466, 30), mate=("=", 466, 0))
    r.add("mate_unmapped", 133, "e1", 466, "*", "ACGTRYACGTNN", mapq=0, mate=("=", 466, 0))
    r.add("mate_far", 65, "e1", 468, "30M", b("e1", 468, 30), mate=("e3", 100, 0))
    r.add("mate_far", 129, "e3", 100, "30M", b("e3", 100, 30), mate=("e1", 468, 0))
    r.add("secondary", 256, "e1", 472, "30M", b("e1", 472, 30))
    r.add("supplementary", 2048, "e1", 474, "10M20H", b("e1", 474, 10),
          tags=["SA:Z:e3,200,+,10S20M,60,0;"])
    r.add("duplicate", 1024, "e1", 476, "30M", b("e1", 476, 30))
    r.add("qc_fail", 512 + 16, "e1", 478, "30M", b("e1", 478, 30))
    r.add("n" * 254, 0, "e1", 480, "30M", b("e1", 480, 30))
    r.add(PUNCTUATION, 0, "e1", 482, "30M", b("e1", 482, 30))
    for at in (484, 486, 488):
        r.add("thrice", 0, "e1", at, "30M", b("e1", at, 30))
    r.add("mapq_0", 0, "e1", 490, "30M", b("e1", 490, 30), mapq=0)
    r.add("mapq_255", 0, "e1", 492, "30M", b("e1", 492, 30), mapq=255)
    r.add("off_end", 0, "e2", 141, "30M", b("e2", 141, 10) + "ACGTNACGTAACGTACGTAC")
    r.add("off_end_noqual", 16, "e2", 145, "20M", b("e2", 145, 6) + "NNNNACGTACGTAC", qual="*")
    r.add("last_base", 0, "e2", 121, "30M", b("e2", 121, 30))
    seq = "".join(b("e4", 1 + at, 1) + "ACGT"[at % 4] for at in range(35000))
    r.add("long_cigar", 0, "e4", 1, "1M1I" * 35000, seq, qual="I" * len(seq))
    r.add("unplaced_pair", 77, "*", 0, "*", "ACGTNRYACGT", mapq=0)
    r.add("unplaced_pair", 141, "*", 0, "*", "TTTTGGGGCC", mapq=0, tags=["RG:Z:g1", "XX:i:-1"])
    r.add("unplaced_no_seq", 4, "*", 0, "*", "*", mapq=0)
    r.add("unplaced_no_qual", 4, "*", 0, "*", "ACGT=N", qual="*", mapq=0)
    return r


def write(directory, name, text):
    sam = Path(directory) / f"{name}.sam"
    sam.write_text(text)
    return to_bam(sam, Path(directory) / f"{name}.bam")


def edges(directory):
    return write(directory, "edges", edge_records().text())


def no_hd(directory):
    return write(directory, "no_hd", edge_records().text(hd=False))


def empty(directory):
    return write(directory, "empty", header())


def slice_records():
    r = Records()
    rng = random.Random(SEED + 2)

    def read(name, contig, pos, flag=0, mate=("*", 0, 0), group="g1"):
        length = rng.choice((36, 50, 75))
        length = min(length, len(r.contigs[contig]) - pos + 1)
        seq = r.bases(contig, pos, length)
        if rng.random() < 0.3:
            at = rng.randrange(length)
            seq = replace(seq, {at: other(seq[at])})
        cigar = f"{length}M"
        if rng.random() < 0.1 and length > 30:
            seq = seq[:10] + "G" + seq[10:-1]
            cigar = f"10M1I{length - 11}M"
        r.add(name, flag, contig, pos, cigar, seq, tags=[f"RG:Z:{group}"], mate=mate)

    for at in range(300):
        read(f"a{at:05d}", "e1", 330 + at % 200)
    for at in range(20):
        read(f"b{at:05d}", "e2", 1 + at * 3, group="g2")
    for at in range(6000):
        left, right = 1 + rng.randrange(2500), 0
        right = left + rng.randrange(50, 300)
        right = min(right, 2900)
        span = right + 75 - left
        read(f"c{at:05d}", "e3", left, 99, ("=", right, span))
        read(f"c{at:05d}", "e3", right, 147, ("=", left, -span), group="g2")
    for at in range(500):
        read(f"d{at:05d}", "e4", 1 + rng.randrange(39000))
    for at in range(750):
        for flag in (77, 141):
            r.add(f"u{at:05d}", flag, "*", 0, "*", "".join(rng.choice("ACGTN") for _ in range(40)),
                  mapq=0, tags=["RG:Z:g1"])
    return r


def slices(directory):
    return write(directory, "slices", slice_records().text())


def shuffled(directory):
    return write(directory, "shuffled", slice_records().text(order="unsorted", shuffle=SEED + 3))


def limits(directory):
    r = Records()
    b = r.bases
    r.add("eq_x", 0, "e3", 10, "10=1X19=", replace(b("e3", 10, 30), {10: other(b("e3", 20, 1))}))
    r.add("m_then_m", 0, "e3", 20, "10M5M15M", b("e3", 20, 30))
    r.add("unmapped_cigar", 4, "e3", 40, "30M", b("e3", 40, 30), mapq=37)
    r.add("unpaired_mate_fields", 0, "e3", 50, "30M", b("e3", 50, 30), mate=("e2", 77, 123))
    return write(directory, "limits", r.text())


# What lies past the last base of e2, by name.
TAILS = (("n", "NNNNNNNN"), ("acgt", "ACGTACGT"), ("iupac", "RYKMSWBD"), ("eq", "=A=C=G=T"),
         ("mix", "NA=RNCGT"))


def overhang_records(r=None, deletion=True):
    r = r or Records()
    b = r.bases
    end = len(r.contigs["e2"])
    for name, tail in TAILS:
        for flag in (0, 16):
            for quals in (True, False):
                label = f"over_{name}_{'r' if flag else 'f'}_{'q' if quals else 'noq'}"
                on = 6 if quals else 9
                r.add(label, flag, "e2", end - on + 1, f"{on + len(tail)}M", b("e2", end - on + 1, on) + tail,
                      qual=None if quals else "*")
    r.add("over_insertion", 0, "e2", end - 4, "7M3I4M", b("e2", end - 4, 5) + "AC" + "GGG" + "TTAA")
    r.add("over_soft_clip", 0, "e2", end - 4, "8M5S", b("e2", end - 4, 5) + "ACG" + "TTTTT")
    if deletion:
        r.add("over_deletion", 0, "e2", end - 4, "3M6D4M", b("e2", end - 4, 3) + "ACGT")
    r.add("starts_on_last_base", 16, "e2", end, "10M", b("e2", end, 1) + "NACGTRACG")
    r.add("ends_on_last_base", 0, "e2", end - 19, "20M", b("e2", end - 19, 20))
    r.add("over_mismatch_before", 0, "e2", end - 9, "16M",
          replace(b("e2", end - 9, 10), {4: other(b("e2", end - 5, 1))}) + "NNACGT")
    r.add("over_two_bases", 0, "e2", end - 27, "30M", b("e2", end - 27, 28) + "NN")
    r.add("over_all_past_but_one", 16, "e2", end, "30M", b("e2", end, 1) + "ACGTN" * 5 + "ACGT",
          qual="*")
    return r


def overhang(directory):
    return write(directory, "overhang", overhang_records().text())


def overhang_clean(directory):
    """`overhang` without the deletion over the contig end, which bammetrics' CleanSam
    refuses (the clip gives a CIGAR that covers another number of bases than the read has;
    Picard writes that record and no reader opens the file)."""
    return write(directory, "overhang_clean", overhang_records(deletion=False).text())


def overhang_shuffled(directory):
    r = Records()
    rng = random.Random(SEED + 4)
    for at in range(600):
        pos = 1 + rng.randrange(2900)
        r.add(f"f{at:04d}", 0, "e3", pos, "50M", r.bases("e3", pos, 50))
        pos = 1 + rng.randrange(500)
        r.add(f"g{at:04d}", 16, "e1", pos, "50M", r.bases("e1", pos, 50))
    overhang_records(r)
    return write(directory, "overhang_shuffled", r.text(order="unsorted", shuffle=SEED + 5))


def write_all(directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    reference(directory)
    return {make.__name__: make(directory)
            for make in (edges, no_hd, empty, slices, shuffled, overhang, overhang_clean,
                         overhang_shuffled, limits)}


if __name__ == "__main__":
    for name, path in write_all(sys.argv[1]).items():
        print(name, path)

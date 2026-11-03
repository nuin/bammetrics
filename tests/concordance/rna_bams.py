"""Hand-built RNA-seq BAMs for CollectRnaSeqMetrics, on a reference of their own
(c1 60 kb, c2 30 kb, cR 8 kb of rRNA, cM 5 kb to ignore).

`build_all(out)` writes ref.fa (.fai, .dict), `rna.refFlat`, `rrna.interval_list`
and one BAM per entry of CASES:

  stranded    a dUTP-style library: the second read lies on the transcript
              strand (SECOND_READ_TRANSCRIPTION_STRAND), a tenth of the pairs
              flipped, so FIRST_READ_TRANSCRIPTION_STRAND counts them the other
              way round.
  unstranded  the same layout with the strand of each pair drawn at random.
  empty       the header and no record.

The gene model (TRANSCRIPTS): GA on the + strand of c1 with a 4-exon and a
3-exon transcript (one exon skipped); GB on the - strand with a coding and a
non-coding transcript; GC on + overlapping GB's last exon, so a read there
touches two genes and counts in no strand column; GD, 401 bases, under the
default MINIMUM_LENGTH; GE on c2 with 5'-biased coverage against GA's 3' bias;
GF single-exon on -; GG with no read; GH whose two lines disagree on strand
and GI with an exon count that is not its number of starts (both dropped by
Picard); GX on a contig the BAM lacks (its line is skipped).

Reads are 100 bases (a few 75), drawn in transcript coordinates and spliced
back to the genome, so reads crossing an exon boundary carry N operators.
Beyond the exonic pairs: intronic and intergenic pairs, pairs inside the rRNA
intervals and pairs half inside one (ribosomal at RRNA_FRAGMENT_PERCENTAGE 0.5,
not at 0.8), pairs in the c1 rRNA interval, reads on cM, pairs reaching past
the gene ends, RF pairs, mate-unmapped and cross-contig pairs, pairs without
MC, pairs with TLEN 0, unpaired reads on both strands, secondary,
supplementary, QC-failed and duplicate records, soft clips, an insertion and
a deletion, and unmapped reads at the end. Read groups: four over two samples
and three libraries, and records without RG.
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, FIRST, MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER, QCFAIL, REVERSE,
                      SECOND, SECONDARY, SUPPLEMENTARY, UNMAPPED, Sam, run)

CONTIGS = [("c1", 60_000), ("c2", 30_000), ("cR", 8_000), ("cM", 5_000)]
GROUPS = [("r1", {"SM": "S1", "LB": "L1", "PU": "P1"}),
          ("r2", {"SM": "S1", "LB": "L1", "PU": "P2"}),
          ("r3", {"SM": "S1", "LB": "L2", "PU": "P3"}),
          ("r4", {"SM": "S2", "LB": "L3", "PU": "P4"})]
RG_CHOICES = ["r1", "r1", "r2", "r3", "r4", None]

# (gene, transcript, contig, strand, exons as 1-based closed (start, end), coding (start, end) or None)
TRANSCRIPTS = [
    ("GA", "GA-1", "c1", "+", [(1000, 1600), (2500, 3100), (4800, 5400), (8200, 9000)], (1300, 8600)),
    ("GA", "GA-2", "c1", "+", [(1000, 1600), (4800, 5400), (8200, 9000)], (1300, 8600)),
    ("GB", "GB-1", "c1", "-", [(12000, 12800), (15000, 15900), (19000, 20000)], (12300, 19500)),
    ("GB", "GB-2", "c1", "-", [(15000, 15900), (19000, 20000)], None),
    ("GC", "GC-1", "c1", "+", [(19500, 20500), (23000, 24000)], (19800, 23800)),
    ("GD", "GD-1", "c1", "-", [(30000, 30400)], (30050, 30350)),
    ("GE", "GE-1", "c2", "+", [(2000, 2900), (6000, 7000), (11000, 12000)], (2200, 11700)),
    ("GE", "GE-2", "c2", "+", [(2000, 2900), (11000, 12000)], (2200, 11700)),
    ("GF", "GF-1", "c2", "-", [(15000, 16000)], (15200, 15800)),
    ("GG", "GG-1", "c2", "+", [(20000, 20600), (20800, 21000)], (20100, 20900)),
]
# (contig, 1-based start, end, name)
RRNA = [("cR", 1, 3000, "rRNA_5S"), ("cR", 3500, 6000, "rRNA_18S"), ("c1", 40000, 41000, "rRNA_c1")]
READ = 100


def write_reference(out, seed=23):
    rng = random.Random(seed)
    ref = {}
    with open(out / "ref.fa", "w") as fasta:
        for name, length in CONTIGS:
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            ref[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, length, 70):
                fasta.write(seq[i:i + 70] + "\n")
    run(["samtools", "faidx", str(out / "ref.fa")])
    run(["picard", "CreateSequenceDictionary", f"R={out / 'ref.fa'}", f"O={out / 'ref.dict'}",
         "QUIET=true", "VERBOSITY=ERROR"])
    return ref


def refflat_line(gene, transcript, contig, strand, exons, coding):
    tx_start, tx_end = exons[0][0] - 1, exons[-1][1]
    cds_start, cds_end = (coding[0] - 1, coding[1]) if coding else (tx_end, tx_end)
    return "\t".join([gene, transcript, contig, strand, str(tx_start), str(tx_end), str(cds_start),
                      str(cds_end), str(len(exons)), "".join(f"{s - 1}," for s, _ in exons),
                      "".join(f"{e}," for _, e in exons)])


def write_annotation(out):
    lines = [refflat_line(*t) for t in TRANSCRIPTS]
    lines.append(refflat_line("GH", "GH-1", "c2", "+", [(25000, 25500)], (25100, 25400)))
    lines.append(refflat_line("GH", "GH-2", "c2", "-", [(25000, 25500)], (25100, 25400)))
    lines.append(refflat_line("GI", "GI-1", "c2", "+", [(26000, 26500)], (26100, 26400)))
    lines.append("GI\tGI-2\tc2\t+\t25999\t26500\t26099\t26400\t2\t25999,\t26500,")
    lines.append(refflat_line("GX", "GX-1", "chrUn", "+", [(100, 900)], (200, 800)))
    (out / "rna.refFlat").write_text("# refFlat\n" + "\n".join(lines) + "\n\n")
    header = "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
    (out / "rrna.interval_list").write_text(
        "@HD\tVN:1.6\tSO:coordinate\n" + header
        + "".join(f"{c}\t{s}\t{e}\t+\t{n}\n" for c, s, e, n in RRNA))


class RnaSam(Sam):
    """`Sam` with this module's contigs in the header."""

    def write(self, out, name, groups):
        header = "@HD\tVN:1.6\tSO:coordinate\n"
        header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
        for rg_id, fields in groups:
            header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
        sam = out / f"{name}.sam"
        sam.write_text(header + "".join(line + "\n" for line in self.lines))
        bam = out / f"{name}.bam"
        run(["samtools", "sort", "-o", str(bam), str(sam)])
        run(["samtools", "index", str(bam)])
        return bam


def transcript_length(exons):
    return sum(e - s + 1 for s, e in exons)


def spliced(sam, contig, exons, tpos, length):
    """Genomic start, CIGAR and sequence of a read over transcript bases
    `tpos` to `tpos + length - 1` (1-based in the concatenated exons)."""
    pieces, offset, start = [], 0, None
    remaining, t = length, tpos
    for s, e in exons:
        size = e - s + 1
        if t > offset + size:
            offset += size
            continue
        gpos = s + (t - offset - 1)
        take = min(remaining, e - gpos + 1)
        if start is None:
            start = gpos
        pieces.append((gpos, take))
        remaining -= take
        t += take
        offset += size
        if remaining == 0:
            break
    cigar, seq = "", ""
    for i, (gpos, take) in enumerate(pieces):
        if i > 0:
            gap = gpos - (pieces[i - 1][0] + pieces[i - 1][1])
            cigar += f"{gap}N"
        cigar += f"{take}M"
        seq += sam.ref_seq(contig, gpos, take)
    return start, cigar, seq


def exonic_pairs(sam, rng, name, contig, exons, strand, count, stranded, bias, prefix, flip=0.1):
    """`count` FR pairs inside the transcript. With `stranded` the second read
    lies on the transcript strand (a `flip` fraction the other way); otherwise
    the pair's strand is a coin flip. `bias` draws positions towards the 3'
    end (1), the 5' end (-1) or evenly (0)."""
    length = transcript_length(exons)
    for i in range(count):
        insert = rng.randint(180, min(320, length))
        u = rng.random()
        if bias > 0:
            u = u ** 0.5
        elif bias < 0:
            u = u ** 2
        if strand == "-":
            u = 1 - u
        tpos = 1 + int(u * (length - insert))
        tpos = max(1, min(tpos, length - insert + 1))
        left = spliced(sam, contig, exons, tpos, READ)
        right = spliced(sam, contig, exons, tpos + insert - READ, READ)
        if stranded:
            second_on_transcript = rng.random() >= flip
            first_left = (strand == "-") == second_on_transcript
        else:
            first_left = rng.random() < 0.5
        rg = rng.choice(RG_CHOICES)
        mapq = (rng.choice((60, 60, 255, 3, 1)), rng.choice((60, 60, 255, 0)))
        flags = PAIRED | PROPER
        if i % 17 == 0:
            flags |= DUPLICATE
        r1, r2 = (left, right) if first_left else (right, left)
        sam.pair(f"{prefix}{name}{i}", contig, r1[0], r2[0], (r1[2], r2[2]), flags=flags,
                 reverse=(not first_left, first_left), mapq=mapq, rg=rg, cigars=(r1[1], r2[1]))
        if i % 23 == 0:
            drop_mate_cigar(sam)


def drop_mate_cigar(sam):
    """Strip MC from the last pair written."""
    for j in (-2, -1):
        fields = [f for f in sam.lines[j].split("\t") if not f.startswith("MC:Z:")]
        sam.lines[j] = "\t".join(fields)


def plain_pair(sam, name, contig, pos, insert, rg=None, reverse=(False, True), flags=PAIRED | PROPER,
               mapq=(60, 60), seqs=None, cigars=(None, None), length=READ):
    seqs = seqs or (sam.ref_seq(contig, pos, length), sam.ref_seq(contig, pos + insert - length, length))
    sam.pair(name, contig, pos, pos + insert - length, seqs, flags=flags, reverse=reverse, mapq=mapq,
             rg=rg, cigars=cigars)


def gene_pairs(sam, rng, stranded):
    bias = {"GA": 1, "GE": -1}
    counts = {"GA-1": 60, "GA-2": 30, "GB-1": 50, "GB-2": 20, "GC-1": 40, "GD-1": 30, "GE-1": 50, "GE-2": 30,
              "GF-1": 40}
    for gene, transcript, contig, strand, exons, _ in TRANSCRIPTS:
        if transcript not in counts:
            continue
        exonic_pairs(sam, rng, transcript, contig, exons, strand, counts[transcript], stranded,
                     bias.get(gene, 0), "e")


def other_reads(sam, rng, stranded):
    # Intronic pairs of GA (1601..2499) and GB (12801..14999), intergenic pairs.
    for i in range(12):
        plain_pair(sam, f"intronA{i}", "c1", rng.randint(1650, 2150), rng.randint(200, 300),
                   rg=rng.choice(RG_CHOICES))
        plain_pair(sam, f"intronB{i}", "c1", rng.randint(12850, 14600), rng.randint(200, 300),
                   rg=rng.choice(RG_CHOICES))
    for i in range(30):
        plain_pair(sam, f"inter{i}", "c1", rng.randint(45000, 55000), rng.randint(200, 400),
                   rg=rng.choice(RG_CHOICES))
    # rRNA: pairs inside the cR intervals, half inside one, and in the c1 interval.
    for i in range(40):
        contig, start, end = ("cR", 1, 3000) if i % 2 else ("cR", 3500, 6000)
        plain_pair(sam, f"rrna{i}", contig, rng.randint(start, end - 400), rng.randint(200, 350),
                   rg=rng.choice(RG_CHOICES), reverse=(True, False) if i % 3 == 0 else (False, True))
    for i in range(10):
        # Fragment 2801..3200 overlaps rRNA_5S by 200 of 400 bases: ribosomal at 0.5, not at 0.8.
        plain_pair(sam, f"rrna_half{i}", "cR", 2801, 400, rg=rng.choice(RG_CHOICES))
    for i in range(10):
        plain_pair(sam, f"rrna_c1_{i}", "c1", rng.randint(40000, 40700), rng.randint(200, 300),
                   rg=rng.choice(RG_CHOICES))
    # Reads on the contig IGNORE_SEQUENCE names, paired and unpaired.
    for i in range(12):
        plain_pair(sam, f"mito{i}", "cM", rng.randint(100, 4500), rng.randint(200, 300),
                   rg=rng.choice(RG_CHOICES))
    for i in range(4):
        pos = rng.randint(100, 4800)
        sam.single(f"mito_u{i}", REVERSE if i % 2 else 0, "cM", pos, sam.ref_seq("cM", pos, READ), rg="r2")
    # Pairs over GB's last exon and GC's first: two genes, no strand count.
    for i in range(15):
        plain_pair(sam, f"two{i}", "c1", rng.randint(19500, 19800), rng.randint(200, 300),
                   rg=rng.choice(RG_CHOICES))
    # Unpaired reads in GA exon 4, GE exon 1 and GF, on both strands.
    for i in range(20):
        contig, pos = [("c1", rng.randint(8200, 8900)), ("c2", rng.randint(2000, 2800)),
                       ("c2", rng.randint(15000, 15900))][i % 3]
        sam.single(f"single{i}", REVERSE if i % 2 else 0, contig, pos, sam.ref_seq(contig, pos, READ),
                   rg=rng.choice(RG_CHOICES))
    # Pairs reaching past a gene end (left of GA, right of GF), and RF pairs inside GA.
    for i in range(10):
        plain_pair(sam, f"past{i}", "c1", 950 + i, 300, rg=rng.choice(RG_CHOICES))
        plain_pair(sam, f"pastF{i}", "c2", 15800 + i, 300, rg=rng.choice(RG_CHOICES))
        plain_pair(sam, f"rf{i}", "c1", rng.randint(8200, 8600), 300, rg=rng.choice(RG_CHOICES),
                   reverse=(True, False), flags=PAIRED)
    # Mate unmapped, cross-contig, TLEN 0.
    for i in range(6):
        pos = rng.randint(4800, 5300)
        seq = sam.ref_seq("c1", pos, READ)
        rg = rng.choice(RG_CHOICES)
        tag = f"RG:Z:{rg}" if rg else None
        sam.record(f"half{i}", PAIRED | FIRST | MATE_UNMAPPED, "c1", pos, 60, f"{READ}M", "=", pos, 0, seq,
                   None, (tag,))
        sam.record(f"half{i}", PAIRED | SECOND | UNMAPPED, "c1", pos, 0, None, "=", pos, 0, seq, None, (tag,))
    for i in range(5):
        pos = rng.randint(1000, 1500)
        sam.pair(f"cross{i}", "c1", pos, 6000 + i * 10, (sam.ref_seq("c1", pos, READ), sam.ref_seq("c2", 6000 + i * 10, READ)),
                 flags=PAIRED, mate_contig="c2", rg=rng.choice(RG_CHOICES))
    for i in range(4):
        pos = rng.randint(15000, 15800)
        seq1, seq2 = sam.ref_seq("c2", pos, READ), sam.ref_seq("c2", pos + 150, READ)
        sam.record(f"tlen0_{i}", PAIRED | FIRST | MATE_REVERSE, "c2", pos, 60, f"{READ}M", "=", pos + 150, 0,
                   seq1, None, (f"MC:Z:{READ}M", "RG:Z:r1"))
        sam.record(f"tlen0_{i}", PAIRED | SECOND | REVERSE, "c2", pos + 150, 60, f"{READ}M", "=", pos, 0,
                   seq2, None, (f"MC:Z:{READ}M", "RG:Z:r1"))
    # Secondary, supplementary, QC-failed records; soft clips, an insertion, a deletion.
    for i in range(6):
        pos = rng.randint(8200, 8900)
        seq = sam.ref_seq("c1", pos, READ)
        sam.single(f"sec{i}", SECONDARY | (REVERSE if i % 2 else 0), "c1", pos, seq, mapq=0, rg="r1")
        sam.single(f"supp{i}", SUPPLEMENTARY | (REVERSE if i % 2 else 0), "c1", pos, seq[:60],
                   cigar="60M40H", rg="r3")
        plain_pair(sam, f"qcfail{i}", "c1", rng.randint(1000, 1400), 200, rg="r4", flags=PAIRED | PROPER | QCFAIL)
    for i in range(8):
        pos = rng.randint(6000, 6800)
        soft = 5 + i
        clipped = "A" * soft + sam.ref_seq("c2", pos, READ - soft)
        plain_pair(sam, f"clip{i}", "c2", pos, 250, rg=rng.choice(RG_CHOICES),
                   seqs=(clipped, sam.ref_seq("c2", pos + 150, READ)), cigars=(f"{soft}S{READ - soft}M", None))
    pos = 11100
    ins = sam.ref_seq("c2", pos, 40) + "GG" + sam.ref_seq("c2", pos + 40, 58)
    plain_pair(sam, "indel_ins", "c2", pos, 260, rg="r2", seqs=(ins, sam.ref_seq("c2", pos + 160, READ)),
               cigars=("40M2I58M", None))
    dele = sam.ref_seq("c2", pos + 20, 50) + sam.ref_seq("c2", pos + 73, 50)
    plain_pair(sam, "indel_del", "c2", pos + 20, 260, rg="r2",
               seqs=(dele, sam.ref_seq("c2", pos + 180, READ)), cigars=("50M3D50M", None))
    # 75-base pairs in GE-1 exon 2 (only in GE-1) and unmapped reads.
    for i in range(10):
        pos = rng.randint(6000, 6700)
        plain_pair(sam, f"short{i}", "c2", pos, 220, rg=rng.choice(RG_CHOICES), length=75)
    for i in range(10):
        rg = rng.choice(RG_CHOICES)
        sam.single(f"unmapped{i}", UNMAPPED, None, 0, "".join(rng.choice("ACGT") for _ in range(READ)), rg=rg)
    if stranded:
        for i in range(3):
            plain_pair(sam, f"extra{i}", "c1", rng.randint(2500, 2900), 220, rg="r1")


def stranded(sam):
    rng = sam.rng
    gene_pairs(sam, rng, True)
    other_reads(sam, rng, True)


def unstranded(sam):
    rng = sam.rng
    gene_pairs(sam, rng, False)
    other_reads(sam, rng, False)


def empty(sam):
    pass


CASES = {
    "stranded": (stranded, GROUPS),
    "unstranded": (unstranded, GROUPS),
    "empty": (empty, GROUPS),
}


def build_all(out: Path, seed=31):
    """Reference, refFlat, rRNA intervals and every CASES BAM in `out`;
    returns {name: bam path}."""
    ref = write_reference(out)
    write_annotation(out)
    bams = {}
    for i, (name, (builder, groups)) in enumerate(CASES.items()):
        sam = RnaSam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups)
    return bams


if __name__ == "__main__":
    import sys
    for name, path in build_all(Path(sys.argv[1])).items():
        print(name, path)

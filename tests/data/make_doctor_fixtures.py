"""Write tests/data/doctor/: the planted defects of `bammetrics doctor`
(docs/doctor.md). Each file differs from good.bam in one way, known by
construction, so each check has a file it passes and a file it does not.

    python3 tests/data/make_doctor_fixtures.py

BAM bytes are written here with struct and zlib (BGZF blocks of at most
4,096 input bytes, so a cut file still holds whole blocks before the cut).
The CRAM goes through samtools (view -C, --no-PG) on tests/data/ref_100k.fa.

good.bam      hg38 header (chr1 to chrY, chrM), SO:coordinate, one @RG with
              SM, LB and PL, @PG MarkDuplicates; 300 paired records on chr1 in
              coordinate order, every 10th pair flagged duplicate, MD and NM on
              every record, qualities Q2 to Q40, read length 100.
good.sam      the same header and records as SAM text.
no_eof.bam    good.bam without its 28-byte end-of-file block.
cut.bam       good.bam cut in the middle of a block.
unsorted.bam  SO:coordinate over the records of good.bam shuffled.
rg_missing.bam    every 7th record carries RG:Z:rg2, which no @RG declares.
rg_duplicate.bam  two @RG lines with ID rg1.
phred64.bam   every quality 31 higher (Q33 to Q71): Phred+64 read as Phred+33.
hg19.bam      the header of hg19 (chr names, chrM of 16,571, chrUn contigs).
grch38.bam    GRCh38 without chr, with alt, decoy and HLA contigs.
t2t.bam       T2T-CHM13 v2.0 (chr names).
mixed.bam     hg38 with chr2 of its GRCh37 length.
unaligned.bam no @SQ, 200 unmapped unpaired records, no @RG.
ref.cram      40 pairs on chr1 of ref_100k.fa as CRAM 3.0 (samtools adds M5).
other.fa      chr1 of 1,000 bases: a reference of other lengths (with .fai).
sheet.csv, sheet.tsv  the sample sheet: NA12878 agrees with good.bam, NA12891
              disagrees (other sample, build and file), DUP has two rows.
"""
import random
import re
import shutil
import struct
import subprocess
import tempfile
import zlib
from pathlib import Path

DATA = Path(__file__).resolve().parent
OUT = DATA / "doctor"
SEED = 7
BLOCK = 4096
EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

# Primary chromosome lengths: GRCh37, GRCh38, T2T-CHM13 v2.0 (UCSC chrom.sizes).
PRIMARY = [
    ("1", 249250621, 248956422, 248387328), ("2", 243199373, 242193529, 242696752),
    ("3", 198022430, 198295559, 201105948), ("4", 191154276, 190214555, 193574945),
    ("5", 180915260, 181538259, 182045439), ("6", 171115067, 170805979, 172126628),
    ("7", 159138663, 159345973, 160567428), ("8", 146364022, 145138636, 146259331),
    ("9", 141213431, 138394717, 150617247), ("10", 135534747, 133797422, 134758134),
    ("11", 135006516, 135086622, 135127769), ("12", 133851895, 133275309, 133324548),
    ("13", 115169878, 114364328, 113566686), ("14", 107349540, 107043718, 101161492),
    ("15", 102531392, 101991189, 99753195), ("16", 90354753, 90338345, 96330374),
    ("17", 81195210, 83257441, 84276897), ("18", 78077248, 80373285, 80542538),
    ("19", 59128983, 58617616, 61707364), ("20", 63025520, 64444167, 66210255),
    ("21", 48129895, 46709983, 45090682), ("22", 51304566, 50818468, 51324926),
    ("X", 155270560, 156040895, 154259566), ("Y", 59373566, 57227415, 62460029),
]
GRCH37, GRCH38, T2T = 1, 2, 3


def dictionary(column, chr_names, mito):
    prefix = "chr" if chr_names else ""
    contigs = [(prefix + row[0], row[column]) for row in PRIMARY]
    contigs.append((mito[0], mito[1]))
    return contigs


HG38 = dictionary(GRCH38, True, ("chrM", 16569))
HG19 = dictionary(GRCH37, True, ("chrM", 16571)) + [
    ("chr1_gl000191_random", 106433), ("chrUn_gl000220", 161802)]
GRCH38_NOCHR = dictionary(GRCH38, False, ("MT", 16569)) + [
    ("KI270706.1", 175055), ("chr6_GL000250v2_alt", 4672374),
    ("chrUn_JTFH01000001v1_decoy", 25139), ("HLA-A*01:01:01:01", 3503)]
T2T_DICT = dictionary(T2T, True, ("chrM", 16569))
MIXED = [(name, PRIMARY[1][GRCH37] if name == "chr2" else length) for name, length in HG38]

RG = "@RG\tID:rg1\tSM:NA12878\tLB:lib1\tPL:ILLUMINA\tPU:unit1\n"
PG = "@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tVN:3.5.0\n"
CIGAR_OPS = "MIDNSHP=X"
SEQ_CODES = "=ACMGRSVTWYHKDBN"


def header(contigs, so="coordinate", rg=RG, pg=PG):
    text = f"@HD\tVN:1.6\tSO:{so}\n"
    text += "".join(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in contigs)
    return text + rg + pg


def reg2bin(beg, end):
    end -= 1
    for shift, offset in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return offset + (beg >> shift)
    return 0


class Record:
    """One alignment, written as BAM bytes or as a SAM line."""

    def __init__(self, name, flag, ref, pos, mapq, cigar, mate_ref, mate_pos, tlen, seq, qual,
                 tags):
        self.name, self.flag, self.ref, self.pos, self.mapq = name, flag, ref, pos, mapq
        self.cigar, self.mate_ref, self.mate_pos, self.tlen = cigar, mate_ref, mate_pos, tlen
        self.seq, self.qual, self.tags = seq, qual, tags

    def bam(self):
        ops = [(int(n), CIGAR_OPS.index(op)) for n, op in self.cigar]
        span = sum(n for n, op in ops if op in (0, 2, 3, 7, 8)) or 1
        bin_ = reg2bin(self.pos, self.pos + span) if self.ref >= 0 else 4680
        name = self.name.encode() + b"\0"
        packed = bytearray()
        for i in range(0, len(self.seq), 2):
            hi = SEQ_CODES.index(self.seq[i]) << 4
            lo = SEQ_CODES.index(self.seq[i + 1]) if i + 1 < len(self.seq) else 0
            packed.append(hi | lo)
        tags = b""
        for tag, kind, value in self.tags:
            if kind == "Z":
                tags += tag.encode() + b"Z" + value.encode() + b"\0"
            else:
                tags += tag.encode() + b"i" + struct.pack("<i", value)
        body = struct.pack("<iiBBHHHiiii", self.ref, self.pos, len(name), self.mapq, bin_,
                           len(ops), self.flag, len(self.seq), self.mate_ref, self.mate_pos,
                           self.tlen)
        body += name + b"".join(struct.pack("<I", n << 4 | op) for n, op in ops)
        body += bytes(packed) + bytes(self.qual) + tags
        return struct.pack("<i", len(body)) + body

    def sam(self, names):
        rname = names[self.ref] if self.ref >= 0 else "*"
        if self.mate_ref < 0:
            rnext = "*"
        elif self.mate_ref == self.ref:
            rnext = "="
        else:
            rnext = names[self.mate_ref]
        cigar = "".join(f"{n}{op}" for n, op in self.cigar) or "*"
        qual = "".join(chr(q + 33) for q in self.qual)
        tags = "".join(f"\t{tag}:{kind}:{value}" for tag, kind, value in self.tags)
        return (f"{self.name}\t{self.flag}\t{rname}\t{self.pos + 1}\t{self.mapq}\t{cigar}\t"
                f"{rnext}\t{self.mate_pos + 1}\t{self.tlen}\t{self.seq}\t{qual}{tags}\n")


def bgzf(payload):
    out = bytearray()
    for start in range(0, len(payload), BLOCK):
        chunk = payload[start:start + BLOCK]
        deflate = zlib.compressobj(6, zlib.DEFLATED, -15)
        data = deflate.compress(chunk) + deflate.flush()
        out += struct.pack("<4BIBBHBBHH", 0x1F, 0x8B, 8, 4, 0, 0, 0xFF, 6, ord("B"), ord("C"), 2,
                           len(data) + 25)
        out += data + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    return bytes(out)


def bam_bytes(text, contigs, records, eof=True):
    raw = text.encode()
    payload = b"BAM\1" + struct.pack("<i", len(raw)) + raw + struct.pack("<i", len(contigs))
    for name, length in contigs:
        payload += struct.pack("<i", len(name) + 1) + name.encode() + b"\0"
        payload += struct.pack("<i", length)
    payload += b"".join(r.bam() for r in records)
    return bgzf(payload) + (EOF if eof else b"")


def write_bam(name, text, contigs, records, eof=True):
    (OUT / name).write_bytes(bam_bytes(text, contigs, records, eof))


def pairs(rng, count, start, quality_offset=0, rg="rg1", rg_every=0, read=100):
    """`count` FR pairs on reference 0 from `start`, in coordinate order."""
    records = []
    bases = "ACGT"
    for i in range(count):
        pos = start + i * 37
        mate = pos + 150 + rng.randrange(50)
        dup = 0x400 if i % 10 == 9 else 0
        name = f"pair{i:04d}"
        for first, (here, there) in enumerate(((pos, mate), (mate, pos))):
            seq = "".join(rng.choice(bases) for _ in range(read))
            qual = [rng.randrange(2, 41) + quality_offset for _ in range(read)]
            flag = (0x1 | 0x2 | dup | (0x20 | 0x40 if first == 0 else 0x10 | 0x80))
            tlen = (mate + read - pos) * (1 if first == 0 else -1)
            group = "rg2" if rg_every and (2 * i + first) % rg_every == 0 else rg
            tags = [("MD", "Z", str(read)), ("NM", "i", 0), ("RG", "Z", group)]
            records.append(Record(name, flag, 0, here, 60, [(read, "M")], 0, there, tlen, seq,
                                  qual, tags))
    records.sort(key=lambda r: (r.pos, r.name))
    return records


def ref_bases():
    lines = (DATA / "ref_100k.fa").read_text().splitlines()
    return "".join(line for line in lines if not line.startswith(">"))


def write_fasta(path, name, bases, width=60):
    body = "".join(bases[i:i + width] + "\n" for i in range(0, len(bases), width))
    path.write_text(f">{name}\n{body}")
    offset = len(name) + 2
    Path(str(path) + ".fai").write_text(f"{name}\t{len(bases)}\t{offset}\t{width}\t{width + 1}\n")


def beyond_dictionary():
    """A header with one contig and a second record on contig number 7."""
    small = [("chr1", 100_000)]
    records = [Record(f"r{i}", 0, ref, 999, 60, [(10, "M")], -1, -1, 0, "ACGTACGTAC", [30] * 10,
                      [("RG", "Z", "rg1")]) for i, ref in enumerate((0, 7))]
    return header(small), small, records


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    rng = random.Random(SEED)
    good = pairs(rng, 150, 10_000)
    hg38 = header(HG38)
    write_bam("good.bam", hg38, HG38, good)
    names = [name for name, _ in HG38]
    (OUT / "good.sam").write_text(hg38 + "".join(r.sam(names) for r in good))
    write_bam("no_eof.bam", hg38, HG38, good, eof=False)
    whole = bam_bytes(hg38, HG38, good)
    (OUT / "cut.bam").write_bytes(whole[: len(whole) * 3 // 5])
    shuffled = list(good)
    random.Random(SEED).shuffle(shuffled)
    write_bam("unsorted.bam", hg38, HG38, shuffled)
    write_bam("rg_missing.bam", hg38, HG38, pairs(random.Random(SEED), 150, 10_000, rg_every=7))
    write_bam("rg_duplicate.bam", header(HG38, rg=RG + RG.replace("lib1", "lib2")), HG38, good)
    phred64 = pairs(random.Random(SEED), 150, 10_000, quality_offset=31)
    write_bam("phred64.bam", hg38, HG38, phred64)
    few = good[:20]
    for name, contigs in (("hg19.bam", HG19), ("grch38.bam", GRCH38_NOCHR), ("t2t.bam", T2T_DICT),
                          ("mixed.bam", MIXED)):
        write_bam(name, header(contigs), contigs, few)
    unaligned = [Record(f"u{i:04d}", 0x4, -1, -1, 0, [], -1, -1, 0, "ACGT" * 25, [30] * 100, [])
                 for i in range(200)]
    write_bam("unaligned.bam", "@HD\tVN:1.6\tSO:unsorted\n", [], unaligned)
    write_bam("refid_big.bam", *beyond_dictionary())

    small = [("chr1", 100_000)]
    ref_pairs = pairs(random.Random(SEED), 40, 1_000)
    with tempfile.TemporaryDirectory() as tmp:
        bam = Path(tmp, "ref.bam")
        bam.write_bytes(bam_bytes(header(small), small, ref_pairs))
        cram = OUT / "ref.cram"
        subprocess.run(["samtools", "view", "--no-PG", "-C", "--output-fmt-option", "version=3.0",
                        "-T", str(DATA / "ref_100k.fa"), "-o", str(cram), str(bam)], check=True)
        # samtools writes the absolute path of -T as UR; the fixture names the file alone.
        text = subprocess.run(["samtools", "view", "-H", "--no-PG", str(cram)], check=True,
                              capture_output=True, text=True).stdout
        text = "".join(f"UR:ref_100k.fa" if field.startswith("UR:") else field
                       for field in re.split(r"(\t|\n)", text))
        Path(tmp, "header.sam").write_text(text)
        subprocess.run(["samtools", "reheader", "--no-PG", "-i", str(Path(tmp, "header.sam")),
                        str(cram)], check=True)
    write_fasta(OUT / "other.fa", "chr1", ref_bases()[:1000])

    rows = [
        ("sample", "bam", "sex", "build", "library"),
        ("NA12878", "good.bam", "female", "hg38", "lib1"),
        ("NA12891", "hg19.bam", "male", "GRCh38", "lib9"),
        ("DUP", "good.bam", "", "", ""),
        ("DUP", "good.bam", "", "", ""),
    ]
    (OUT / "sheet.csv").write_text("".join(",".join(row) + "\n" for row in rows))
    (OUT / "sheet.tsv").write_text("".join("\t".join(row) + "\n" for row in rows))


if __name__ == "__main__":
    main()

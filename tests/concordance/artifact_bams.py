"""Inputs for CollectSequencingArtifactMetrics built from the foundation corpus.

    build_all(out, foundation) -> {name: (bam, reference)}

* decorated: the corpus BAM on a copy of its reference with chr1 2001 to 4000
  in lower case, an R every 500 bases of chr2 and an N every 700 bases of chr1.
  Every 13th record fails vendor QC, every 5th carries an OQ tag 15 above or
  below its qualities, every 23rd mapped record has MAPQ 255, every 29th is
  made unpaired, and every 17th and 19th change read group: to a header line
  without LB (UnknownLibrary, sample S3) and to an ID the header lacks (also
  UnknownLibrary).
* no_read_groups: the corpus BAM without its @RG lines; the records keep their
  RG tags, so a counted record has a library missing from the header.
* unplaced_middle: the corpus BAM with an unmapped pair without a reference
  between the chr1 and chr2 records. Picard stops reading there.
* sparse: one pair whose forward read mismatches every 7th base, so contexts
  with an error rate of 1 print QSCORE 0 (never -0).
* oq_below_bang: the corpus BAM with an OQ of spaces on every 7th record,
  which Picard rejects (Invalid fastq character).

    without_contigs(vcf, out) -> the VCF without ##contig lines, plus a POS 0 record

    write_dbsnp(reference, out) -> (plain VCF, bgzip VCF with .tbi)

writes a sites-only VCF on the reference's dictionary: a SNP every 97 bases,
a 4-base deletion every 1009, an insertion every 1511, an MNP every 2003 and a
filtered SNP every 3001, over the first 30 kb of each contig.
"""
import shutil
import subprocess
from pathlib import Path

import pysam

PAIRED_ONLY = 0x1 | 0x2 | 0x8 | 0x20 | 0x40 | 0x80


def decorate_reference(source: Path, out: Path) -> Path:
    """The reference with a lower-case stretch, R and N bases; .fai and .dict beside it."""
    fasta = pysam.FastaFile(str(source))
    path = out / "decorated.fa"
    with open(path, "w") as handle:
        for name in fasta.references:
            seq = list(fasta.fetch(name))
            if name == "chr1":
                for i in range(2000, min(4000, len(seq))):
                    seq[i] = seq[i].lower()
                for i in range(350, len(seq), 700):
                    seq[i] = "N"
            if name == "chr2":
                for i in range(250, len(seq), 500):
                    seq[i] = "R"
            text = "".join(seq)
            handle.write(f">{name}\n")
            for i in range(0, len(text), 60):
                handle.write(text[i:i + 60] + "\n")
    pysam.faidx(str(path))
    subprocess.run(["samtools", "dict", "-o", str(out / "decorated.dict"), str(path)], check=True,
                   capture_output=True)
    return path


def decorated_bam(source: Path, out: Path) -> Path:
    path = out / "decorated.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        header = bam.header.to_dict()
        header["RG"].append({"ID": "S3.noLB", "SM": "S3", "PL": "ILLUMINA", "PU": "FC301.1"})
        with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
            for i, r in enumerate(bam):
                if i % 13 == 0:
                    r.flag |= 0x200
                if i % 5 == 0 and r.query_qualities is not None:
                    shift = 15 if i % 10 == 0 else -15
                    oq = [max(0, min(60, q + shift)) for q in r.query_qualities]
                    r.set_tag("OQ", "".join(chr(q + 33) for q in oq))
                if i % 23 == 0 and not r.is_unmapped:
                    r.mapping_quality = 255
                if i % 29 == 0:
                    r.flag &= ~PAIRED_ONLY
                    r.next_reference_id = -1
                    r.next_reference_start = -1
                    r.template_length = 0
                    if r.has_tag("MC"):
                        r.set_tag("MC", None)
                if i % 17 == 0:
                    r.set_tag("RG", "S3.noLB")
                elif i % 19 == 0:
                    r.set_tag("RG", "ghost")
                sink.write(r)
    return path


def no_read_groups_bam(source: Path, out: Path) -> Path:
    path = out / "no_read_groups.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        header = bam.header.to_dict()
        header.pop("RG", None)
        with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
            for r in bam:
                sink.write(r)
    return path


def unplaced_middle_bam(source: Path, out: Path) -> Path:
    path = out / "unplaced_middle.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        header = bam.header.to_dict()
        records = list(bam)
    first_chr2 = next(i for i, r in enumerate(records) if r.reference_id == 1)
    with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
        for i, r in enumerate(records):
            if i == first_chr2:
                for flag in (77, 141):
                    u = pysam.AlignedSegment(sink.header)
                    u.query_name = "unplaced"
                    u.flag = flag
                    u.reference_id = -1
                    u.reference_start = -1
                    u.next_reference_id = -1
                    u.next_reference_start = -1
                    u.query_sequence = "ACGT" * 25
                    u.query_qualities = [30] * 100
                    u.set_tag("RG", records[0].get_tag("RG"))
                    sink.write(u)
            sink.write(r)
    return path


def sparse_bam(source: Path, reference: Path, out: Path) -> Path:
    """One FR pair on chr1 whose forward read mismatches every 7th base: with
    CONTEXT_SIZE=2 most mismatch contexts hold that one base, so their error
    rate is 1 and QSCORE prints 0."""
    path = out / "sparse.bam"
    seq = pysam.FastaFile(str(reference)).fetch("chr1", 10000, 10300).upper()
    swap = {"A": "C", "C": "G", "G": "T", "T": "A"}
    with pysam.AlignmentFile(str(source)) as bam:
        header = bam.header.to_dict()
    group = header["RG"][0]["ID"]
    with pysam.AlignmentFile(str(path), "wb", header=header) as sink:
        for flag, start, tlen in ((99, 0, 300), (147, 200, -300)):
            bases = list(seq[start:start + 100])
            if flag == 99:
                for i in range(3, 97, 7):
                    bases[i] = swap.get(bases[i], "A")
            r = pysam.AlignedSegment(sink.header)
            r.query_name = "sparse"
            r.flag = flag
            r.reference_id = 0
            r.reference_start = 10000 + start
            r.mapping_quality = 60
            r.cigarstring = "100M"
            r.next_reference_id = 0
            r.next_reference_start = 10200 if flag == 99 else 10000
            r.template_length = tlen
            r.query_sequence = "".join(bases)
            r.query_qualities = [30] * 100
            r.set_tag("RG", group)
            r.set_tag("MC", "100M")
            sink.write(r)
    return path


def oq_below_bang_bam(source: Path, out: Path) -> Path:
    """The corpus BAM with an OQ of spaces (below '!') on every 7th record."""
    path = out / "oq_below_bang.bam"
    with pysam.AlignmentFile(str(source)) as bam:
        with pysam.AlignmentFile(str(path), "wb", template=bam) as sink:
            for i, r in enumerate(bam):
                if i % 7 == 0 and r.query_length:
                    r.set_tag("OQ", " " * r.query_length)
                sink.write(r)
    return path


def build_all(out: Path, foundation: Path):
    """{name: (bam, reference)} of the BAMs."""
    out.mkdir(parents=True, exist_ok=True)
    corpus, reference = foundation / "corpus.bam", foundation / "ref.fa"
    return {
        "decorated": (decorated_bam(corpus, out), decorate_reference(reference, out)),
        "no_read_groups": (no_read_groups_bam(corpus, out), reference),
        "unplaced_middle": (unplaced_middle_bam(corpus, out), reference),
        "sparse": (sparse_bam(corpus, reference, out), reference),
        "oq_below_bang": (oq_below_bang_bam(corpus, out), reference),
    }


def without_contigs(vcf: Path, out: Path) -> Path:
    """`vcf` without its ##contig lines and with a POS 0 record first: Picard's
    artifact tool reads both."""
    lines = vcf.read_text().splitlines()
    header = [line for line in lines if line.startswith("#") and not line.startswith("##contig")]
    records = [line for line in lines if not line.startswith("#")]
    contig = records[0].split("\t", 1)[0]
    path = out / f"{vcf.stem}.no_contigs.vcf"
    path.write_text("\n".join([*header, f"{contig}\t0\trs0\tA\tC\t.\tPASS\t.", *records]) + "\n")
    return path


def write_dbsnp(reference: Path, out: Path, span=30000):
    """(plain VCF, bgzip VCF indexed by tabix) of dbSNP-like sites on `reference`."""
    fasta = pysam.FastaFile(str(reference))
    lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=LowQual,Description="Low quality">']
    lines += [f"##contig=<ID={n},length={fasta.get_reference_length(n)}>" for n in fasta.references]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    for name in fasta.references:
        seq = fasta.fetch(name).upper()
        end = min(span, len(seq) - 10)
        sites = {}
        for pos in range(97, end, 97):
            ref = seq[pos - 1]
            sites[pos] = (ref, "T" if ref != "T" else "C", "PASS")
        for pos in range(1009, end, 1009):
            sites[pos] = (seq[pos - 1:pos + 3], seq[pos - 1], "PASS")
        for pos in range(1511, end, 1511):
            sites[pos] = (seq[pos - 1], seq[pos - 1] + "AC", "PASS")
        for pos in range(2003, end, 2003):
            ref = seq[pos - 1:pos + 1]
            sites[pos] = (ref, "".join("T" if b != "T" else "C" for b in ref), "PASS")
        for pos in range(3001, end, 3001):
            ref = seq[pos - 1]
            sites[pos] = (ref, "G" if ref != "G" else "A", "LowQual")
        for pos in sorted(sites):
            ref, alt, filt = sites[pos]
            if "N" in ref or "N" in alt:
                continue
            lines.append(f"{name}\t{pos}\trs{pos}\t{ref}\t{alt}\t.\t{filt}\t.")
    stem = f"{reference.parent.name}_{reference.stem}"
    plain = out / f"{stem}.dbsnp.vcf"
    plain.write_text("\n".join(lines) + "\n")
    gz = out / f"{stem}.dbsnp.vcf.gz"
    shutil.copy(plain, out / "tmp.vcf")
    pysam.tabix_index(str(out / "tmp.vcf"), preset="vcf", force=True)
    (out / "tmp.vcf.gz").rename(gz)
    (out / "tmp.vcf.gz.tbi").rename(gz.with_name(gz.name + ".tbi"))
    return plain, gz

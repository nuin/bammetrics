"""Write BAMs whose index counts are known by construction, for
CollectIndexDepthMetrics (docs/06 section 7).

The BAM bytes come from the struct writer of tests/concordance/index_bams.py, which
keeps a reference of length 0 in the header. samtools (through pysam) writes each
.bai and .csi. Every planted mapped read is 100M; lengths and counts are powers of
two, so every reads-per-base ratio is exact in binary. Each contig of PLANTED gives
(name, length, mapped records, placed unmapped records). On chr1 of `male` the
mapped records include 4 duplicates, 2 QC failures, 3 secondary records without a
sequence and 3 supplementary records of 40 bases: the index counts them all, the
read length sample (mapped primary records with a sequence) skips the last two kinds.

    python3 tests/data/make_index_depth_bams.py OUT        # PLANTED, .bai beside each,
                                                           # OUT/csi/ the same with .csi
    python3 tests/data/make_index_depth_bams.py --fixture  # tests/data/index_depth.bam,
                                                           # .bam.bai and .bam.csi

Needs pysam.
"""
import argparse
import shutil
import sys
from pathlib import Path

import pysam

DATA = Path(__file__).resolve().parent
sys.path.insert(0, str(DATA.parent / "concordance"))
import index_bams as ib  # noqa: E402

READ = 100
DUPLICATE, QC_FAIL, SECONDARY, SUPPLEMENTARY = 0x400, 0x200, 0x100, 0x800

MALE = [
    ("chr1", 8192, 256, 7),
    ("chr2", 4096, 128, 0),
    ("chrX", 8192, 128, 0),
    ("chrY", 4096, 64, 0),
    ("chrM", 1024, 512, 0),
    ("chr1_KI270706v1_alt", 2048, 16, 0),
    ("chrUn_JTFH01000001v1_decoy", 1024, 8, 0),
    ("HLA-A*01:01:01:01", 4096, 8, 0),
    ("chrEBV", 1024, 4, 0),
    ("chrUn_GL000220v1", 1024, 0, 0),
    ("zero", 0, 0, 0),
]

# name: (contigs, records without a coordinate, read lengths in file order or None for
# all READ, expected READ_LENGTH at the default sample, expected INFERRED_SEX,
# expected MT_COPY_NUMBER or None for `?`).
PLANTED = {
    "male": (MALE, 9, None, 100.0, "MALE", 32.0),
    # No mitochondrial contig; X at the autosomal depth, no Y read.
    "female_no_mt": ([("1", 8192, 256, 3), ("2", 8192, 256, 0), ("X", 4096, 128, 0),
                      ("Y", 2048, 0, 0)], 0, None, 100.0, "FEMALE", None),
    # References and no record: depth 0, no read length.
    "empty": ([("chr1", 4096, 0, 0), ("chrX", 4096, 0, 0)], 0, None, None, "UNKNOWN", None),
    # No reference: every record lacks a coordinate.
    "unaligned": ([], 5, None, None, "UNKNOWN", None),
    # Ten 100-base reads, then thirty 150-base reads: 100 over a sample of 10, 137.5
    # over all.
    "mixed_lengths": ([("chr1", 8192, 40, 0), ("chrM", 1024, 0, 0)], 0,
                      [100] * 10 + [150] * 30, 137.5, "UNKNOWN", 0.0),
}


def contig_records(ref, name, length, mapped, unmapped, lengths, special):
    """`mapped` records tiled over the contig, then `unmapped` placed unmapped records.
    With `special`, the last 12 mapped records carry the flags the index counts and the
    sample partly skips."""
    flags = [0] * mapped
    if special:
        flags[-12:] = [DUPLICATE] * 4 + [QC_FAIL] * 2 + [SECONDARY] * 3 + [SUPPLEMENTARY] * 3
    out = []
    span = max(length - 200, 1)
    for i in range(mapped):
        pos = 1 + (i * 37) % span
        read = lengths[i] if lengths else READ
        flag = flags[i]
        if flag == SECONDARY:
            out.append(ib.record(f"{name}_s{i}", ref, pos, cigar=f"{READ}M", flag=flag,
                                 seq_len=0))
        elif flag == SUPPLEMENTARY:
            out.append(ib.record(f"{name}_p{i}", ref, pos, cigar="40M60H", flag=flag))
        else:
            out.append(ib.record(f"{name}_r{i}", ref, pos, cigar=f"{read}M", flag=flag))
    out += [ib.unmapped(f"{name}_u{i}", ref, 1 + i * 11) for i in range(unmapped)]
    return out


def write(path, contigs, no_coordinate, lengths=None, special=False):
    references = [(name, length) for name, length, _, _ in contigs]
    records = []
    for ref, (name, length, mapped, unmapped) in enumerate(contigs):
        records += contig_records(ref, name, length, mapped, unmapped, lengths,
                                  special and ref == 0)
    records.sort(key=ib.position)
    records += ib.unplaced(no_coordinate)
    path.write_bytes(ib.bam(references, records))


def index(bam, csi):
    """samtools index: `x.bam.bai`, or `x.bam.csi` with `csi`."""
    pysam.index("-c" if csi else "-b", str(bam))


def write_planted(out):
    out.mkdir(parents=True, exist_ok=True)
    (out / "csi").mkdir(exist_ok=True)
    for name, (contigs, no_coordinate, lengths, *_) in PLANTED.items():
        bam = out / f"{name}.bam"
        write(bam, contigs, no_coordinate, lengths, special=name == "male")
        index(bam, csi=False)
        copy = out / "csi" / f"{name}.bam"
        shutil.copyfile(bam, copy)
        index(copy, csi=True)


def write_fixture():
    bam = DATA / "index_depth.bam"
    contigs, no_coordinate, *_ = PLANTED["male"]
    write(bam, contigs, no_coordinate, special=True)
    index(bam, csi=False)
    index(bam, csi=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--fixture", action="store_true")
    args = parser.parse_args()
    if args.fixture:
        write_fixture()
    elif args.out:
        write_planted(args.out)
    else:
        parser.error("give OUT or --fixture")


if __name__ == "__main__":
    main()

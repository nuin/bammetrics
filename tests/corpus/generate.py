"""Generate the synthetic concordance corpus under tests/corpus/out/.

Usage:
    uv run --no-project --with pysam --with numpy python tests/corpus/generate.py [OUT_DIR]

Needs samtools, bwa and picard (3.5.0) on PATH. Everything is derived from
SEED, so two runs write the same reads, intervals and reference.

Tiers:
  t0_*  tiny hand-built BAMs with known histogram shapes (written by pysam)
  t1    20k clean FR pairs, two libraries, 70% of fragments on target
  t2    50k messy pairs: adapters, indels, Q2 tails, RF and TANDEM pairs,
        cross-contig chimeras, split reads (supplementary), unmapped mates,
        repeats (low MAPQ, secondary with bwa mem -a), PCR duplicates marked
        by picard MarkDuplicates, QC-fail flags
  t3    10k single-end reads
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pysam

SEED = 20261028
READ_LEN = 150
CONTIGS = [("chr1", 300_000), ("chr2", 200_000), ("chrM", 16_000)]
N_TARGETS = 120
BAIT_PAD = 50
ON_TARGET = 0.7
SAMPLE = "SYN1"
LIBRARIES = [("rgA", "libA"), ("rgB", "libB")]
# Illumina TruSeq adapter starts: read 1 runs into ADAPTER_R1 when the
# fragment is shorter than the read, read 2 into ADAPTER_R2.
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
# Two repeat copies: an exact one (MAPQ 0, secondary records) and a 2%
# diverged one (intermediate MAPQ). (source contig, start, length, dest contig, start).
REPEATS = [("chr1", 50_000, 3_000, "chr2", 100_000, 0.0),
           ("chr1", 150_000, 2_000, "chr2", 20_000, 0.02)]

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = np.array(list("ACGT"))


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def run(cmd, **kwargs):
    subprocess.run(cmd, check=True, **kwargs)


def random_seq(rng, n):
    return "".join(BASES[rng.integers(0, 4, n)])


# --- reference and intervals -------------------------------------------------

def make_reference(rng, out):
    seqs = {name: list(random_seq(rng, length)) for name, length in CONTIGS}
    for src, s_start, length, dst, d_start, divergence in REPEATS:
        copy = seqs[src][s_start:s_start + length]
        for i in np.nonzero(rng.random(length) < divergence)[0]:
            copy[i] = BASES[(list("ACGT").index(copy[i]) + rng.integers(1, 4)) % 4]
        seqs[dst][d_start:d_start + length] = copy
    fasta = out / "reference.fa"
    with fasta.open("w") as f:
        for name, _ in CONTIGS:
            seq = "".join(seqs[name])
            f.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")
    run(["samtools", "faidx", str(fasta)])
    dict_path = out / "reference.dict"
    run(["samtools", "dict", "-u", "reference.fa", "-o", str(dict_path), str(fasta)])
    run(["bwa", "index", str(fasta)], stderr=subprocess.DEVNULL)
    return {name: "".join(s) for name, s in seqs.items()}


def make_targets(rng):
    """120 sorted, non-overlapping targets of 100-400 bp on chr1 and chr2.

    Gaps are mostly wide, but every tenth gap is under 2 * BAIT_PAD so the
    padded baits of neighbouring targets merge.
    """
    targets = []
    for contig, length, count in (("chr1", 300_000, 80), ("chr2", 200_000, 40)):
        pos = 2_000
        for i in range(count):
            size = int(rng.integers(100, 401))
            targets.append((contig, pos, pos + size))  # 0-based half open
            gap = int(rng.integers(20, 90)) if i % 10 == 9 else int(rng.integers(800, 2_800))
            pos += size + gap
        assert pos < length - 2_000, contig
    return targets


def merge(intervals):
    merged = []
    for contig, start, end in intervals:
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1] = (contig, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((contig, start, end))
    return merged


def write_intervals(out, name, intervals):
    header = (out / "reference.dict").read_text()
    with (out / f"{name}.interval_list").open("w") as il, (out / f"{name}.bed").open("w") as bed:
        il.write(header)
        for i, (contig, start, end) in enumerate(intervals, 1):
            label = f"{name[:-1]}_{i:03d}"
            il.write(f"{contig}\t{start + 1}\t{end}\t+\t{label}\n")
            bed.write(f"{contig}\t{start}\t{end}\t{label}\n")


# --- read simulation -----------------------------------------------------------

class Simulator:
    def __init__(self, rng, reference, targets, messy):
        self.rng = rng
        self.ref = reference
        self.targets = targets
        self.messy = messy
        self.names = [name for name, _ in CONTIGS]
        lengths = np.array([length for _, length in CONTIGS], dtype=float)
        self.weights = lengths / lengths.sum()

    def locus(self, frag_len):
        """Fragment start: ON_TARGET of fragments overlap a target."""
        rng = self.rng
        if rng.random() < ON_TARGET:
            contig, start, end = self.targets[rng.integers(len(self.targets))]
            centre = int(rng.integers(start, end))
            pos = centre - int(rng.integers(0, frag_len))
        else:
            contig = self.names[rng.choice(len(self.names), p=self.weights)]
            pos = int(rng.integers(0, len(self.ref[contig]) - frag_len))
        pos = min(max(pos, 0), len(self.ref[contig]) - frag_len)
        return contig, pos

    def fragment(self, contig, pos, frag_len, reverse):
        seq = self.ref[contig][pos:pos + frag_len]
        return revcomp(seq) if reverse else seq

    def qualities(self):
        rng = self.rng
        if self.messy and rng.random() < 0.1:
            q = rng.integers(2, 42, READ_LEN)
        else:
            # Binned NovaSeq-style qualities, dropping toward the 3' end.
            q = np.full(READ_LEN, 37)
            drop = int(rng.integers(90, 160))
            if drop < READ_LEN:
                q[drop:] = rng.choice([12, 23, 37], READ_LEN - drop, p=[0.3, 0.5, 0.2])
            q[rng.random(READ_LEN) < 0.02] = 23
        if self.messy and rng.random() < 0.08:
            q[READ_LEN - int(rng.integers(10, 41)):] = 2  # Illumina Q2 "B tail"
        return q

    def read(self, seq):
        """Pad or cut to READ_LEN, add an indel sometimes, then errors by quality."""
        rng = self.rng
        if self.messy and rng.random() < 0.05:
            at = int(rng.integers(20, 130))
            size = int(rng.integers(1, 4))
            if rng.random() < 0.5:
                seq = seq[:at] + random_seq(rng, size) + seq[at:]
            else:
                seq = seq[:at] + seq[at + size:]
        if len(seq) < READ_LEN:
            seq += random_seq(rng, READ_LEN - len(seq))
        seq = np.array(list(seq[:READ_LEN]))
        q = self.qualities()
        err = rng.random(READ_LEN) < 10 ** (-q / 10)
        seq[err] = BASES[rng.integers(0, 4, int(err.sum()))]
        return "".join(seq), "".join(chr(33 + int(x)) for x in q)

    def pair(self, contig, pos, frag_len, reverse, kind):
        frag = self.fragment(contig, pos, frag_len, reverse)
        if kind == "adapter":
            r1 = frag + ADAPTER_R1
            r2 = revcomp(frag) + ADAPTER_R2
        elif kind == "rf":
            r1, r2 = revcomp(frag[:READ_LEN]), frag[-READ_LEN:]
        elif kind == "tandem":
            r1, r2 = frag[:READ_LEN], frag[-READ_LEN:]
        elif kind == "cross_contig":
            other = self.names[(self.names.index(contig) + 1) % len(self.names)]
            o_pos = int(self.rng.integers(0, len(self.ref[other]) - frag_len))
            r1, r2 = frag[:READ_LEN], revcomp(self.ref[other][o_pos:o_pos + frag_len][-READ_LEN:])
        elif kind == "unmapped_mate":
            r1, r2 = frag[:READ_LEN], random_seq(self.rng, READ_LEN)
        elif kind == "split":
            other_pos = int(self.rng.integers(0, len(self.ref[contig]) - 200))
            r1 = frag[:90] + self.ref[contig][other_pos:other_pos + 60]
            r2 = revcomp(frag[-READ_LEN:])
        else:
            r1, r2 = frag[:READ_LEN], revcomp(frag[-READ_LEN:])
        return self.read(r1), self.read(r2)

    def kind(self):
        if not self.messy:
            return "normal"
        x = self.rng.random()
        for kind, p in (("adapter", 0.02), ("rf", 0.03), ("tandem", 0.01),
                        ("cross_contig", 0.02), ("unmapped_mate", 0.01), ("split", 0.02)):
            if x < p:
                return kind
            x -= p
        return "normal"

    def frag_len(self, kind):
        if kind == "adapter":
            return int(self.rng.integers(60, 140))
        return int(np.clip(self.rng.normal(300, 60), 160, 700))

    def paired(self, tier, n_pairs, dup_rate):
        """Yield (library index, name, (seq1, qual1), (seq2, qual2))."""
        rng = self.rng
        history = []
        for i in range(n_pairs):
            lib = int(rng.integers(len(LIBRARIES)))
            if history and rng.random() < dup_rate:
                lib, contig, pos, frag_len, reverse, kind = history[rng.integers(len(history))]
            else:
                kind = self.kind()
                frag_len = self.frag_len(kind)
                contig, pos = self.locus(frag_len)
                reverse = bool(rng.random() < 0.5)
                history.append((lib, contig, pos, frag_len, reverse, kind))
            name = f"{tier}:{i}"
            if self.messy and rng.random() < 0.01:
                name += ":qcfail"
            r1, r2 = self.pair(contig, pos, frag_len, reverse, kind)
            yield lib, name, r1, r2

    def single(self, tier, n_reads):
        for i in range(n_reads):
            frag_len = self.frag_len("normal")
            contig, pos = self.locus(frag_len)
            reverse = bool(self.rng.random() < 0.5)
            yield int(self.rng.integers(len(LIBRARIES))), f"{tier}:{i}", self.read(
                self.fragment(contig, pos, frag_len, reverse)[:READ_LEN])


# --- alignment -------------------------------------------------------------------

def write_fastq(path, records):
    with path.open("w") as f:
        for name, (seq, qual) in records:
            f.write(f"@{name}\n{seq}\n+\n{qual}\n")


def align(out, tier, per_library_fastqs, all_alignments):
    """bwa mem per library, fixmate -m (MC and ms tags), sort, merge,
    MarkDuplicates, then QC-fail flags from the read names."""
    work = out / f"{tier}.work"
    reference = out / "reference.fa"
    sorted_bams = []
    for (rg, lib), fastqs in zip(LIBRARIES, per_library_fastqs):
        read_group = f"@RG\\tID:{rg}\\tSM:{SAMPLE}\\tLB:{lib}\\tPL:ILLUMINA"
        bwa = ["bwa", "mem", "-t", "4", "-K", "10000000", "-R", read_group]
        if all_alignments:
            bwa.append("-a")
        bwa += [str(reference)] + [str(p) for p in fastqs]
        fixmate = work / f"{lib}.fixmate.bam"
        sorted_bam = work / f"{lib}.sorted.bam"
        with subprocess.Popen(bwa, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL) as p:
            run(["samtools", "fixmate", "-m", "-", str(fixmate)], stdin=p.stdout)
        if p.returncode:
            raise subprocess.CalledProcessError(p.returncode, bwa)
        run(["samtools", "sort", "-o", str(sorted_bam), str(fixmate)])
        sorted_bams.append(sorted_bam)
    merged = work / "merged.bam"
    run(["samtools", "merge", "-f", "-o", str(merged)] + [str(b) for b in sorted_bams])
    marked = work / "marked.bam"
    run(["picard", "MarkDuplicates", f"I={merged}", f"O={marked}",
         f"M={out / f'{tier}.duplicate_metrics.txt'}", "VALIDATION_STRINGENCY=SILENT",
         "QUIET=true", "VERBOSITY=ERROR"], stderr=subprocess.DEVNULL)
    final = out / f"{tier}.bam"
    with pysam.AlignmentFile(str(marked)) as src, \
            pysam.AlignmentFile(str(final), "wb", template=src) as dst:
        for r in src:
            if r.query_name.endswith(":qcfail"):
                r.flag |= 0x200
            dst.write(r)
    run(["samtools", "index", str(final)])
    shutil.rmtree(work)


def make_paired_tier(sim, out, tier, n_pairs, dup_rate, all_alignments):
    work = out / f"{tier}.work"
    work.mkdir(exist_ok=True)
    reads = [([], []) for _ in LIBRARIES]
    for lib, name, r1, r2 in sim.paired(tier, n_pairs, dup_rate):
        reads[lib][0].append((name, r1))
        reads[lib][1].append((name, r2))
    fastqs = []
    for (_, lib), (r1s, r2s) in zip(LIBRARIES, reads):
        paths = (work / f"{lib}_R1.fq", work / f"{lib}_R2.fq")
        write_fastq(paths[0], r1s)
        write_fastq(paths[1], r2s)
        fastqs.append(paths)
    align(out, tier, fastqs, all_alignments)


def make_single_tier(sim, out, tier, n_reads):
    work = out / f"{tier}.work"
    work.mkdir(exist_ok=True)
    reads = [[] for _ in LIBRARIES]
    for lib, name, r in sim.single(tier, n_reads):
        reads[lib].append((name, r))
    fastqs = []
    for (_, lib), records in zip(LIBRARIES, reads):
        path = work / f"{lib}.fq"
        write_fastq(path, records)
        fastqs.append((path,))
    align(out, tier, fastqs, all_alignments=False)


# --- t0: hand-built histogram shapes ------------------------------------------------

PAIRED, PROPER, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x10, 0x20
FIRST, SECOND, SECONDARY, DUPLICATE, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x400, 0x800
T0_LEN = 100


def t0_pair(name, insert, orientation="FR", pos=1000, extra=0, mc=True, contigs=(0, 0)):
    """Two records of one pair with |TLEN| = insert. Read 2 is the counted read.

    FR: read 1 forward at pos, read 2 reverse ending at pos + insert - 1.
    RF: read 2 reverse at pos, read 1 forward starting at pos + insert - 100.
    TANDEM: both forward.
    """
    far = pos + insert - T0_LEN
    if orientation == "FR":
        r1, r2 = (pos, False), (far, True)
    elif orientation == "RF":
        r1, r2 = (far, False), (pos, True)
    else:
        r1, r2 = (pos, False), (far, False)
    records = []
    for mine, mate, which, contig, mate_contig in ((r1, r2, FIRST, contigs[0], contigs[1]),
                                                   (r2, r1, SECOND, contigs[1], contigs[0])):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.flag = PAIRED | PROPER | which | extra
        a.flag |= REVERSE if mine[1] else 0
        a.flag |= MATE_REVERSE if mate[1] else 0
        a.reference_id = contig
        a.reference_start = mine[0]
        a.mapping_quality = 60
        a.cigarstring = f"{T0_LEN}M"
        a.query_sequence = "A" * T0_LEN
        a.query_qualities = pysam.qualitystring_to_array("I" * T0_LEN)
        a.next_reference_id = mate_contig
        a.next_reference_start = mate[0]
        a.template_length = insert if mine[0] < mate[0] or (
            mine[0] == mate[0] and which == FIRST) else -insert
        if mc:
            a.set_tag("MC", f"{T0_LEN}M")
        records.append(a)
    return records


def fr_pairs(prefix, inserts, orientation="FR", start=1000):
    records = []
    for i, insert in enumerate(inserts):
        records += t0_pair(f"{prefix}{i}", insert, orientation, pos=start + i * 2_000)
    return records


def t0_cases():
    """Name -> records. Each case pins one Picard histogram or filter rule."""
    cases = {
        # Even count, median 101.5: half-integer MAD keys and trim width int(101.5 + 10 * MAD).
        "t0_even_half_median": fr_pairs("e", [100, 101, 102, 105, 140, 141]),
        # Odd count: median is the middle item.
        "t0_odd": fr_pairs("o", [100, 150, 151, 152, 300]),
        # One pair: trimmed SD is 0/0, printed '?'.
        "t0_single_pair": fr_pairs("s", [250]),
        # Mode tie: lowest key wins.
        "t0_mode_tie": fr_pairs("m", [200, 200, 220, 220, 230]),
        "t0_all_tandem": fr_pairs("t", [300, 310, 310, 330, 400, 1_000], "TANDEM"),
        # RF at exactly 5% (1 of 20): dropped by the float MINIMUM_PCT.
        "t0_cutoff_at": fr_pairs("a", [300] * 19) + fr_pairs("ar", [500], "RF", start=90_000),
        # RF at 1 of 19 (5.3%): kept.
        "t0_cutoff_over": fr_pairs("b", [300] * 18) + fr_pairs("br", [500], "RF", start=90_000),
        # Wide spread: every WIDTH_OF_X_PERCENT differs, long tail trimmed from mean/SD.
        "t0_spread": fr_pairs("w", [150 + 7 * i for i in range(40)] + [2_000, 5_000]),
    }
    # Filters: these records must be skipped, leaving the 5 FR pairs.
    filtered = fr_pairs("f", [300, 310, 310, 320, 330])
    filtered += t0_pair("dup", 900, pos=20_000, extra=DUPLICATE)
    filtered += t0_pair("sec", 900, pos=22_000, extra=SECONDARY)
    filtered += t0_pair("sup", 900, pos=24_000, extra=SUPPLEMENTARY)
    # Cross-contig pair with a non-zero TLEN (bwa writes 0; other aligners may not).
    filtered += t0_pair("xc", 777, pos=26_000, contigs=(0, 1))
    cases["t0_filters"] = filtered
    # Single-end reads only: no pairs, Picard writes no metrics file.
    se = []
    for i in range(4):
        a = pysam.AlignedSegment()
        a.query_name = f"se{i}"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1_000 + 300 * i
        a.mapping_quality = 60
        a.cigarstring = f"{T0_LEN}M"
        a.query_sequence = "A" * T0_LEN
        a.query_qualities = pysam.qualitystring_to_array("I" * T0_LEN)
        se.append(a)
    cases["t0_no_pairs"] = se
    return cases


def make_t0(out):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 200_000}, {"SN": "chr2", "LN": 200_000}]}
    for name, records in t0_cases().items():
        records.sort(key=lambda r: (r.reference_id, r.reference_start, r.query_name, r.flag))
        path = out / f"{name}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=header) as f:
            for r in records:
                f.write(r)
        pysam.index(str(path))


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent / "out")
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    reference = make_reference(rng, out)
    targets = make_targets(rng)
    write_intervals(out, "targets", targets)
    baits = merge([(c, max(s - BAIT_PAD, 0), e + BAIT_PAD) for c, s, e in targets])
    write_intervals(out, "baits", baits)
    make_t0(out)
    make_paired_tier(Simulator(rng, reference, targets, messy=False), out, "t1",
                     20_000, dup_rate=0.0, all_alignments=False)
    make_paired_tier(Simulator(rng, reference, targets, messy=True), out, "t2",
                     50_000, dup_rate=0.05, all_alignments=True)
    make_single_tier(Simulator(rng, reference, targets, messy=False), out, "t3", 10_000)
    summary = {"seed": SEED, "targets": len(targets), "baits": len(baits)}
    (out / "corpus.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()

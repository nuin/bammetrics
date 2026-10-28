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
  t0_qy hand-built reads for CollectQualityYieldMetrics: OQ tags of every
        length, no QUAL, quality thresholds, QC-fail, secondary, supplementary
  t0_hs hand-built reads for CollectHsMetrics: overlap clipping, filters and
        bait distances; hs_targets.interval_list holds unsorted, duplicate,
        overlapping, abutting and empty targets

t1 to t3 also carry an OQ tag on a quarter of the reads, with the QUAL field
lowered the way base recalibration does, so the OQ and QUAL totals differ.
"""
import json
import shutil
import zlib
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


# Fraction of records (1 in OQ_EVERY, picked by a hash of the read name so that
# both mates and all alignments of a read agree) that get an OQ tag.
OQ_EVERY = 4
RECALIBRATION_SHIFT = 8


def add_original_qualities(record):
    """Move QUAL into OQ and lower QUAL by RECALIBRATION_SHIFT (floor Q2)."""
    quals = record.query_qualities
    if quals is None or zlib.crc32(record.query_name.encode()) % OQ_EVERY:
        return
    record.set_tag("OQ", pysam.qualities_to_qualitystring(quals))
    record.query_qualities = pysam.qualitystring_to_array(
        "".join(chr(33 + max(2, q - RECALIBRATION_SHIFT)) for q in quals))


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
            add_original_qualities(r)
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


# --- t0 for CollectAlignmentSummaryMetrics: one record per Picard rule ------------

QCFAIL, UNMAPPED, MATE_UNMAPPED = 0x200, 0x4, 0x8
# Starts with a 16-base window of the INDEXED 3' adapter (AdapterUtility).
ADAPTER_READ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"


def substitute(seq, positions, bases="ACGT"):
    """Replace seq[i] by the next base in `bases` that differs from it."""
    seq = list(seq)
    for i in positions:
        seq[i] = next(b for b in bases if b != seq[i])
    return "".join(seq)


class AsmBuilder:
    """pysam records on the corpus reference for the alignment-summary t0 BAMs."""

    def __init__(self, reference):
        self.ref = reference
        self.contig_ids = {name: i for i, (name, _) in enumerate(CONTIGS)}
        self.records = []

    def ref_seq(self, contig, pos, length):
        return self.ref[contig][pos:pos + length]

    def make(self, name, flag, contig, pos, cigar, seq, mapq=60, quals=None, tags=(),
             mate=None, tlen=0):
        a = pysam.AlignedSegment()
        a.query_name = name
        a.flag = flag
        if contig is None:
            a.reference_id, a.reference_start = -1, -1
        else:
            a.reference_id, a.reference_start = self.contig_ids[contig], pos
        a.mapping_quality = 0 if flag & UNMAPPED else mapq
        if cigar and not flag & UNMAPPED:
            a.cigarstring = cigar
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(quals or "I" * len(seq))
        if mate is not None:
            mate_contig, mate_pos = mate
            a.next_reference_id = -1 if mate_contig is None else self.contig_ids[mate_contig]
            a.next_reference_start = mate_pos
            a.template_length = tlen
        for tag, value in tags:
            a.set_tag(tag, value)
        self.records.append(a)
        return a

    def pair(self, name, contig, pos1, pos2, len1=100, len2=100, flags=PAIRED | PROPER,
             reverse=(False, True), mapq=(60, 60), mq_tags=True, mate_contig=None,
             seqs=(None, None), cigars=(None, None), tags=((), ())):
        """Two mapped mates; TLEN from outer ends when on one contig, else 0."""
        mate_contig = mate_contig or contig
        cigars = tuple(c or f"{n}M" for c, n in zip(cigars, (len1, len2)))
        seqs = tuple(s or self.ref_seq(c, p, n) for s, c, p, n in
                     zip(seqs, (contig, mate_contig), (pos1, pos2), (len1, len2)))
        if contig == mate_contig:
            span = max(pos1 + len1, pos2 + len2) - min(pos1, pos2)
            tlens = (span, -span) if pos1 <= pos2 else (-span, span)
        else:
            tlens = (0, 0)
        ends = ((contig, pos1), (mate_contig, pos2))
        for i, which in enumerate((FIRST, SECOND)):
            j = 1 - i
            flag = flags | which
            flag |= REVERSE if reverse[i] else 0
            flag |= MATE_REVERSE if reverse[j] else 0
            extra = [("MC", cigars[j])]
            if mq_tags:
                extra.append(("MQ", mapq[j]))
            self.make(name, flag, ends[i][0], ends[i][1], cigars[i], seqs[i], mapq=mapq[i],
                      tags=extra + list(tags[i]), mate=ends[j], tlen=tlens[i])

    def write(self, path):
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": n, "LN": length} for n, length in CONTIGS]}
        placed = [r for r in self.records if r.reference_id >= 0]
        unplaced = [r for r in self.records if r.reference_id < 0]
        placed.sort(key=lambda r: (r.reference_id, r.reference_start, r.query_name, r.flag))
        with pysam.AlignmentFile(str(path), "wb", header=header) as f:
            for r in placed + unplaced:
                f.write(r)
        pysam.index(str(path))


def t0_asm(reference):
    """Every AlignmentSummaryMetrics rule on a few reads: four categories,
    mismatches, IUPAC and '=' read bases, clips, indels, noise, adapters,
    chimeras, MQ gating, QC-fail, secondary and supplementary records."""
    b = AsmBuilder(reference)
    r = b.ref_seq
    # Unpaired reads.
    b.make("u_clean", 0, "chr1", 1_000, "100M", r("chr1", 1_000, 100))
    quals = "I" * 50 + "+" * 50  # Q40 then Q10
    b.make("u_mismatch", 0, "chr1", 2_000, "100M",
           substitute(r("chr1", 2_000, 100), [3, 40, 70]), quals=quals)
    b.make("u_low_mapq", REVERSE, "chr1", 3_000, "100M",
           substitute(r("chr1", 3_000, 100), [10, 90]), mapq=10)
    iupac = list(r("chr1", 4_000, 100))
    iupac[5], iupac[6], iupac[7], iupac[8] = "R", "N", "=", "n"
    b.make("u_iupac", 0, "chr1", 4_000, "100M", "".join(iupac))
    b.make("u_softclip", 0, "chr1", 5_000, "5S90M5S",
           "TTTTT" + r("chr1", 5_000, 90) + "GGGGG")
    b.make("u_softclip_rev", REVERSE, "chr1", 6_000, "10S85M5S",
           "C" * 10 + r("chr1", 6_000, 85) + "A" * 5)
    b.make("u_hardclip", 0, "chr1", 7_000, "10H90M", r("chr1", 7_000, 90))
    seq = r("chr1", 8_000, 40) + "TT" + r("chr1", 8_040, 30) + r("chr1", 8_073, 28)
    b.make("u_indel", 0, "chr1", 8_000, "40M2I30M3D28M", seq)
    b.make("u_skip", 0, "chr1", 8_500, "50M100N50M",
           r("chr1", 8_500, 50) + r("chr1", 8_650, 50))
    b.make("u_eq_x", 0, "chr1", 8_800, "20=1X79=",
           substitute(r("chr1", 8_800, 100), [20]))
    b.make("u_sa_hq", 0, "chr1", 9_000, "60M40S", r("chr1", 9_000, 60) + "A" * 40,
           tags=[("SA", "chr2,5000,+,60S40M,60,0;")])
    b.make("u_sa_low", 0, "chr1", 9_500, "100M", r("chr1", 9_500, 100), mapq=5,
           tags=[("SA", "chr2,6000,+,100M,5,0;")])
    b.make("u_noise", 0, "chr1", 10_000, "100M", r("chr1", 10_000, 100), tags=[("XN", 1)])
    b.make("u_noise_two", 0, "chr1", 10_100, "100M", r("chr1", 10_100, 100), tags=[("XN", 2)])
    b.make("u_noise_text", 0, "chr1", 10_200, "100M", r("chr1", 10_200, 100),
           tags=[("XN", "1")])
    b.make("u_qcfail", QCFAIL, "chr1", 11_000, "100M",
           substitute(r("chr1", 11_000, 100), [0, 1]) [:50] + "N" * 50)
    b.make("u_adapter_mq0", 0, "chr1", 12_000, f"{len(ADAPTER_READ)}M", ADAPTER_READ, mapq=0)
    b.make("u_adapter_mq0_rev", REVERSE, "chr1", 12_500, f"{len(ADAPTER_READ)}M",
           revcomp(ADAPTER_READ), mapq=0)
    b.make("u_adapter_mq60", 0, "chr1", 13_000, f"{len(ADAPTER_READ)}M", ADAPTER_READ)
    b.make("u_short", 0, "chr1", 14_000, "36M", r("chr1", 14_000, 36))
    # Runs 40 bases past the end of chrM: only the first 60 are compared.
    b.make("u_past_end", 0, "chrM", 15_940, "100M", r("chrM", 15_940, 60) + "C" * 40)
    b.make("u_secondary", SECONDARY, "chr1", 15_000, "100M", "N" * 100)
    b.make("u_supplementary", SUPPLEMENTARY, "chr2", 5_000, "60H40M", r("chr2", 5_000, 40))
    b.make("u_unmapped_adapter", UNMAPPED, None, 0, None, ADAPTER_READ)
    b.make("u_unmapped_one_error", UNMAPPED, None, 0, None, substitute(ADAPTER_READ, [4]))
    b.make("u_unmapped_two_errors", UNMAPPED, None, 0, None, substitute(ADAPTER_READ, [4, 9]))
    b.make("u_unmapped_n", UNMAPPED | REVERSE, None, 0, None, "NNA" + "ACGT" * 20)
    # Pairs.
    b.pair("p_fr", "chr1", 20_000, 20_200)
    b.pair("p_lengths", "chr1", 21_000, 21_250, len1=120, len2=75)
    b.pair("p_low_mapq", "chr1", 22_000, 22_200, mapq=(60, 10))
    b.pair("p_no_mq_tag", "chr1", 23_000, 23_200, mapq=(60, 10), mq_tags=False)
    b.pair("p_rf", "chr1", 24_000, 24_200, reverse=(True, False), flags=PAIRED)
    b.pair("p_tandem", "chr1", 25_000, 25_200, reverse=(False, False), flags=PAIRED)
    b.pair("p_far", "chr1", 26_000, 186_000, flags=PAIRED)
    b.pair("p_cross", "chr1", 27_000, 9_000, mate_contig="chr2", flags=PAIRED)
    b.pair("p_sa", "chr1", 28_000, 28_200, tags=([("SA", "chr2,7000,+,100M,60,0;")], ()))
    b.pair("p_improper", "chr1", 29_000, 29_300, flags=PAIRED)
    b.pair("p_dup", "chr1", 30_000, 30_200, flags=PAIRED | PROPER | DUPLICATE)
    b.pair("p_qcfail", "chr1", 31_000, 31_200, flags=PAIRED | PROPER | QCFAIL,
           seqs=("N" * 100, None))
    b.pair("p_blocks", "chr1", 32_000, 32_300, cigars=(None, "30M2D70M"),
           seqs=(None, substitute(r("chr1", 32_300, 30) + r("chr1", 32_332, 70), [5, 35])))
    b.pair("p_clipped", "chr1", 33_000, 33_200, cigars=("3H97M", "90M10S"),
           seqs=(r("chr1", 33_000, 97), r("chr1", 33_200, 90) + "G" * 10), len1=97, len2=90)
    # Half-mapped pair: read 1 mapped, read 2 unmapped and placed at read 1.
    b.make("p_half", PAIRED | FIRST | MATE_UNMAPPED, "chr2", 10_000, "100M",
           r("chr2", 10_000, 100), mate=("chr2", 10_000))
    b.make("p_half", PAIRED | SECOND | UNMAPPED, "chr2", 10_000, None, ADAPTER_READ,
           mate=("chr2", 10_000))
    return b


def t0_asm_cycles(reference):
    """BAD_CYCLES: a cycle is bad when no-calls plus mismatches there reach
    80% of TOTAL_READS. Five pairs and three unpaired reads.

    FIRST_OF_PAIR: N at cycle 5 in all five, N at cycle 51 in four -> 2.
    SECOND_OF_PAIR: reverse 30M2D70M reads with a mismatch at read index 0
    (cycle 100) in all five, and a mismatch at read index 5 in three and 35 in
    two. Picard takes the cycle from the offset in the alignment block, so
    index 35 (offset 5 in the second block) also lands on cycle 95: 2 bad.
    PAIR prints 2 + 2. UNPAIRED: an N at the first cycle of a QC-failed, an
    unmapped and a mapped read -> 1."""
    b = AsmBuilder(reference)
    r = b.ref_seq
    for i in range(5):
        pos = 40_000 + 1_000 * i
        first = list(r("chr1", pos, 100))
        first[4] = "N"
        if i < 4:
            first[50] = "N"
        second = r("chr1", pos + 300, 30) + r("chr1", pos + 332, 70)
        second = substitute(second, [0, 5 if i < 3 else 35])
        b.pair(f"c{i}", "chr1", pos, pos + 300, cigars=(None, "30M2D70M"),
               seqs=("".join(first), second))
    b.make("cu_qcfail", QCFAIL, "chr2", 1_000, "100M", "N" + r("chr2", 1_001, 99))
    b.make("cu_mapped", 0, "chr2", 2_000, "100M", "N" + r("chr2", 2_001, 99))
    b.make("cu_unmapped", UNMAPPED, None, 0, None, "N" + r("chr2", 3_001, 99))
    return b


def t0_qy(reference):
    """CollectQualityYieldMetrics rules on a few reads: Q20 and Q30 edges,
    OQ longer, shorter, empty and different from QUAL, a record without QUAL,
    and the records the default settings skip or keep."""
    b = AsmBuilder(reference)
    r = b.ref_seq

    def qstring(values):
        return "".join(chr(33 + v) for v in values)

    edges = [0, 1, 19, 20, 21, 29, 30, 31, 41, 60, 93]
    b.make("q_plain", 0, "chr1", 1_000, "100M", r("chr1", 1_000, 100))
    b.make("q_edges", 0, "chr1", 2_000, "44M", r("chr1", 2_000, 44), quals=qstring(edges * 4))
    b.make("q_oq_differs", REVERSE, "chr1", 3_000, "100M", r("chr1", 3_000, 100),
           tags=[("OQ", qstring(edges * 9 + [30]))])
    b.make("q_oq_short", 0, "chr1", 4_000, "100M", r("chr1", 4_000, 100),
           tags=[("OQ", qstring([35] * 60))])
    b.make("q_oq_long", 0, "chr1", 5_000, "100M", r("chr1", 5_000, 100),
           tags=[("OQ", qstring([25] * 120))])
    b.make("q_oq_empty", 0, "chr1", 6_000, "100M", r("chr1", 6_000, 100),
           quals=qstring([22] * 100), tags=[("OQ", "")])
    b.make("q_no_qual", 0, "chr1", 7_000, "100M", r("chr1", 7_000, 100)).query_qualities = None
    b.make("q_qcfail", QCFAIL, "chr1", 8_000, "100M", r("chr1", 8_000, 100),
           quals=qstring([40] * 100), tags=[("OQ", qstring([10] * 100))])
    b.make("q_secondary", SECONDARY, "chr1", 9_000, "100M", r("chr1", 9_000, 100),
           quals=qstring([33] * 100))
    b.make("q_supplementary", SUPPLEMENTARY, "chr2", 5_000, "60H40M", r("chr2", 5_000, 40),
           quals=qstring([28] * 40))
    b.make("q_unmapped", UNMAPPED, None, 0, None, "ACGT" * 15, quals=qstring([20] * 60))
    b.make("q_short", 0, "chr1", 10_000, "36M", r("chr1", 10_000, 36), quals=qstring([30] * 36))
    b.pair("q_pair", "chr1", 20_000, 20_200, len1=120, len2=75)
    b.pair("q_pair_qcfail", "chr1", 21_000, 21_200, flags=PAIRED | PROPER | QCFAIL)
    return b


def t0_hs(reference, targets):
    """CollectHsMetrics rules on a few reads placed on the corpus targets:
    overlapping mates clipped as htsjdk does (equal starts, = and X elements,
    insertions, trailing soft clips and deletions at the clip, a mate on
    another contig), the duplicate, adapter, MAPQ and base-quality filters,
    QC-fail, secondary, supplementary and mate-unmapped records, near-bait
    and off-bait reads, and a pile-up deeper than a small COVERAGE_CAP."""
    b = AsmBuilder(reference)
    r = b.ref_seq
    t = [start for contig, start, _ in targets if contig == "chr1"]
    # Overlapping mates: the left end loses the overlap.
    b.pair("p_overlap", "chr1", t[0], t[0] + 50)
    # Same start: the first end keeps everything, the second is dropped.
    b.pair("p_same_start", "chr1", t[1], t[1])
    # = and X elements straddling the mate start are clipped whole.
    b.pair("p_seq_match", "chr1", t[2], t[2] + 30, cigars=("40=1X59=", None),
           seqs=(substitute(r("chr1", t[2], 100), [40]), None))
    # An insertion before the clip point and a trailing soft clip.
    b.pair("p_insertion", "chr1", t[3], t[3] + 47, len1=95, cigars=("45M5I50M", None),
           seqs=(r("chr1", t[3], 45) + "ACGTA" + r("chr1", t[3] + 45, 50), None))
    b.pair("p_soft_clip", "chr1", t[4], t[4] + 60, len1=90, cigars=("90M10S", None),
           seqs=(r("chr1", t[4], 90) + "A" * 10, None))
    # The mate starts inside a deletion.
    b.pair("p_deletion", "chr1", t[5], t[5] + 55, len1=110, cigars=("50M10D50M", None),
           seqs=(r("chr1", t[5], 50) + r("chr1", t[5] + 60, 50), None))
    # A mate on chr2 at a larger position still clips (contigs are not compared).
    b.pair("p_other_contig", "chr1", t[6], t[6] + 40, mate_contig="chr2",
           flags=PAIRED)
    # Duplicates and unique first ends near the baits give HS_LIBRARY_SIZE.
    for i in range(4):
        b.pair(f"p_dup{i}", "chr1", t[7] + 10 * i, t[7] + 150,
               flags=PAIRED | PROPER | (DUPLICATE if i % 2 else 0))
    # Adapter at MAPQ 0, MAPQ 0 without adapter, MAPQ 10.
    b.make("u_adapter", 0, "chr1", t[8], "100M", ADAPTER_READ + r("chr1", t[8] + 64, 36), mapq=0)
    b.make("u_mapq0", 0, "chr1", t[8] + 20, "100M", r("chr1", t[8] + 20, 100), mapq=0)
    b.make("u_mapq10", REVERSE, "chr1", t[8] + 40, "100M", r("chr1", t[8] + 40, 100), mapq=10)
    # Base qualities 2, 3, 10 and 40.
    quals = "#" * 10 + "$" * 10 + "+" * 30 + "I" * 50
    b.make("u_low_quality", 0, "chr1", t[9], "100M", r("chr1", t[9], 100), quals=quals)
    # Partly off target, then near-bait and off-bait reads.
    b.make("u_edge", 0, "chr1", t[10] - 60, "100M", r("chr1", t[10] - 60, 100))
    b.make("u_near", 0, "chr1", t[11] - 250, "100M", r("chr1", t[11] - 250, 100))
    b.make("u_off", 0, "chr1", t[11] - 1_000, "100M", r("chr1", t[11] - 1_000, 100))
    # Flags: QC-fail, secondary, supplementary (hard clipped), mate unmapped.
    b.make("u_qcfail", QCFAIL, "chr1", t[12], "100M", r("chr1", t[12], 100))
    b.make("u_secondary", SECONDARY, "chr1", t[12] + 10, "100M", r("chr1", t[12] + 10, 100))
    b.make("u_supplementary", SUPPLEMENTARY, "chr1", t[12] + 20, "50H50M",
           r("chr1", t[12] + 20, 50))
    b.make("p_mate_unmapped", PAIRED | FIRST | MATE_UNMAPPED, "chr1", t[13], "100M",
           r("chr1", t[13], 100), mate=("chr1", t[13]))
    b.make("p_mate_unmapped", PAIRED | SECOND | UNMAPPED, "chr1", t[13], None,
           r("chr1", t[13], 100), mate=("chr1", t[13]))
    # Indels on target for INCLUDE_INDELS.
    b.make("u_indels", 0, "chr1", t[14], "30M3I30M4D37M",
           r("chr1", t[14], 30) + "TTT" + r("chr1", t[14] + 30, 30) + r("chr1", t[14] + 64, 37))
    # Eight reads on one spot: depth 8 against COVERAGE_CAP=5.
    for i in range(8):
        b.make(f"u_stack{i}", 0, "chr1", t[15], "100M", r("chr1", t[15], 100),
               quals="5" * 50 + "I" * 50)
    return b


def write_hs_intervals(out, targets):
    """Targets for the t0_hs merge rules, in a deliberately unsorted list:
    two copies of one target, an overlapping and an abutting pair, an empty
    interval and a contig missing from the dictionary (skipped by htsjdk)."""
    (c0, s0, e0), (_, s1, e1), (_, s2, e2) = targets[0], targets[1], targets[2]
    rows = [
        (c0, s2 + 1, e2), (c0, s0 + 1, e0), (c0, s0 + 1, e0),
        (c0, s1 + 1, s1 + 100), (c0, s1 + 50, e1),
        (c0, e1 + 1, e1 + 30), (c0, e2 + 500, e2 + 499), ("chrUn", 1, 10),
    ]
    header = (out / "reference.dict").read_text()
    with (out / "hs_targets.interval_list").open("w") as f:
        f.write(header)
        for i, (contig, start, end) in enumerate(rows, 1):
            f.write(f"{contig}\t{start}\t{end}\t+\ths_{i}\n")


def make_t0_asm(out, reference):
    t0_asm(reference).write(out / "t0_asm.bam")
    t0_asm_cycles(reference).write(out / "t0_asm_cycles.bam")
    AsmBuilder(reference).write(out / "t0_asm_empty.bam")
    t0_qy(reference).write(out / "t0_qy.bam")


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
    make_t0_asm(out, reference)
    t0_hs(reference, targets).write(out / "t0_hs.bam")
    write_hs_intervals(out, targets)
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

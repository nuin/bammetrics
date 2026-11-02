"""Differential fuzzer: seeded random BAMs through Picard 3.5.0 and bammetrics
(parity plan section 5.4), with a shrinker for the failures.

    python tests/concordance/fuzz.py --seed 1 --n 200 [--jobs 8] [--out target/fuzz]

One example is a BamSpec drawn from `random.Random(seed * 1_000_003 + index)`:
1 to 3 contigs of 1 to 5 kb with N runs, lower-case and IUPAC bases; 0 to 3
read groups with or without SM, LB and PU; up to 40 templates, paired or
single, with read lengths from {1, 36, 100, 150, 1000}, CIGARs over
M I D S H N = X P, MAPQ from {0, 1, 19, 20, 21, 60, 255}, base qualities from
{0, 2, 19, 20, 21, 40, 93}, flag combinations, tags MC (right, wrong or
missing), MQ, SA, XN, OQ, RG, NM, MD, BC and BX, Illumina-style or plain
read names with clustered tile coordinates, duplicate sets with copies of
the sequence, unmapped mates and templates, and empty files.

Each example runs three commands on both tools, with options drawn from a
small menu: CollectMultipleMetrics with every program (seven files),
CollectDuplicateMetrics and EstimateLibraryComplexity. Picard output is
cached by the sha256 of the BAM bytes and the arguments under
`<out>/cache`. A mismatch is any file present on one side only, any file
differing after the header lines, or one tool failing where the other
succeeds; a case where both fail is counted apart (BOTH-ERR), never as
agreement.

A failing example is shrunk with ddmin over its templates (mates stay
together, the header stays) to the fewest templates that still fail, and
written to tests/concordance/regressions/<sha8>.sam with <sha8>.fa (the
reference), <sha8>.args (the failing command) and <sha8>.picard.txt
(Picard's files). test_regressions.py replays every case in that folder.
"""
import argparse
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from cases import BAMMETRICS, ROOT  # noqa: E402
from compare import body  # noqa: E402
from run_picard import picard_env  # noqa: E402

REGRESSIONS = Path(__file__).parent / "regressions"
PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
FIRST, SECOND, SECONDARY, QCFAIL, DUPLICATE, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x200, 0x400, 0x800
READ_LENGTHS = (1, 36, 100, 150, 1000)
MAPQS = (0, 1, 19, 20, 21, 60, 255)
QUALITIES = (0, 2, 19, 20, 21, 40, 93)
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
SILENT = "VALIDATION_STRINGENCY=SILENT"


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


# --- the spec ------------------------------------------------------------------------

class Template:
    def __init__(self, name, records):
        self.name = name
        self.records = records  # SAM lines without the name


class BamSpec:
    def __init__(self, contigs, groups, templates, sort_order="coordinate"):
        self.contigs = contigs        # [(name, sequence)]
        self.groups = groups          # [(id, {field: value})]
        self.templates = templates
        self.sort_order = sort_order

    def header(self):
        lines = [f"@HD\tVN:1.6\tSO:{self.sort_order}"]
        lines += [f"@SQ\tSN:{name}\tLN:{len(seq)}" for name, seq in self.contigs]
        for rg_id, fields in self.groups:
            lines.append("@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()))
        return "\n".join(lines) + "\n"

    def sam_text(self):
        text = self.header()
        for template in self.templates:
            for record in template.records:
                text += template.name + "\t" + record + "\n"
        return text

    def subset(self, templates):
        return BamSpec(self.contigs, self.groups, templates, self.sort_order)


def random_sequence(rng, length):
    bases = [rng.choice("ACGT") for _ in range(length)]
    for _ in range(rng.choice((0, 0, 1, 2))):
        start = rng.randrange(length)
        for i in range(start, min(length, start + rng.randint(1, 40))):
            bases[i] = "N"
    for _ in range(rng.choice((0, 0, 3))):
        i = rng.randrange(length)
        bases[i] = rng.choice("RYKMSWBDHVrykmswbdhv")
    if rng.random() < 0.3:
        start = rng.randrange(length)
        for i in range(start, min(length, start + rng.randint(1, 200))):
            bases[i] = bases[i].lower()
    return "".join(bases)


def random_cigar(rng, read_length, max_span):
    """A CIGAR consuming `read_length` query bases, spanning at most
    `max_span` reference bases. Returns (cigar, reference span)."""
    if read_length == 0:
        return "*", 0
    ops = []
    remaining = read_length
    left_hard = rng.choice((0, 0, 0, 3))
    right_hard = rng.choice((0, 0, 0, 5))
    left_soft = rng.choice((0, 0, 0, 1, 10)) if remaining > 1 else 0
    left_soft = min(left_soft, remaining - 1)
    remaining -= left_soft
    right_soft = rng.choice((0, 0, 0, 2, 7)) if remaining > 1 else 0
    right_soft = min(right_soft, remaining - 1)
    remaining -= right_soft
    middle = []
    span = 0
    first = True
    while remaining > 0:
        kind = rng.choice("MMMMM=XIDNP") if not first else rng.choice("MMM=X")
        first = False
        if kind in "M=X":
            n = rng.randint(1, remaining)
            middle.append((n, kind))
            remaining -= n
            span += n
        elif kind == "I":
            n = rng.randint(1, remaining)
            middle.append((n, "I"))
            remaining -= n
        elif kind == "D":
            middle.append((rng.randint(1, 5), "D"))
            span += middle[-1][0]
        elif kind == "N":
            middle.append((rng.randint(1, 300), "N"))
            span += middle[-1][0]
        else:
            middle.append((rng.randint(1, 3), "P"))
    while middle and middle[-1][1] in "IDNP":
        last = middle.pop()
        if last[1] == "I":
            remaining += last[0]
        elif last[1] in "DN":
            span -= last[0]
    if remaining:
        middle.append((remaining, "M"))
        span += remaining
    if span > max_span:
        return None, None
    if left_hard:
        ops.append((left_hard, "H"))
    if left_soft:
        ops.append((left_soft, "S"))
    ops += middle
    if right_soft:
        ops.append((right_soft, "S"))
    if right_hard:
        ops.append((right_hard, "H"))
    merged = []
    for n, kind in ops:
        if merged and merged[-1][1] == kind:
            merged[-1] = (merged[-1][0] + n, kind)
        else:
            merged.append((n, kind))
    return "".join(f"{n}{k}" for n, k in merged), span


def read_name(rng, index, cluster):
    tile, x, y = cluster
    form = rng.choice(("7", "7", "5", "plain", "6"))
    x += rng.choice((0, 0, 1, 50, 99, 100, 101, 2500, 10_000))
    y += rng.choice((0, 0, 1, 100, 3000))
    if form == "7":
        return f"M{index % 2}:1:FC0:{1 + index % 2}:{tile}:{x}:{y}"
    if form == "5":
        return f"M0:1:{tile}:{x}:{y}"
    if form == "6":
        return f"M0:1:FC0:{tile}:{x}:{y}"
    return f"r{index}"


def random_quals(rng, length):
    kind = rng.random()
    if kind < 0.2:
        q = rng.choice(QUALITIES)
        return "".join(chr(33 + q) for _ in range(length))
    return "".join(chr(33 + rng.choice(QUALITIES)) for _ in range(length))


def random_tags(rng, mapped, mate_cigar, mate_mapq, rg, seq):
    tags = []
    if rg is not None:
        tags.append(f"RG:Z:{rg}")
    if mate_cigar is not None:
        choice = rng.choice(("right", "right", "missing", "wrong"))
        if choice == "right":
            tags.append(f"MC:Z:{mate_cigar}")
        elif choice == "wrong":
            tags.append(f"MC:Z:{rng.randint(1, 50)}M")
    if mate_mapq is not None and rng.random() < 0.7:
        tags.append(f"MQ:i:{mate_mapq}")
    if mapped and rng.random() < 0.15:
        tags.append("SA:Z:c1,100,+,10M,60,0;")
    if rng.random() < 0.1:
        tags.append(f"XN:i:{rng.choice((1, 2))}")
    # An OQ shorter than the read makes Picard's QualityScoreDistribution
    # throw (a both-fail case the read-level tests cover), so OQ is as long
    # as the read or one base longer.
    if rng.random() < 0.3:
        tags.append("OQ:Z:" + random_quals(rng, len(seq) + rng.choice((0, 0, 0, 1))))
    if mapped and rng.random() < 0.3:
        tags.append(f"NM:i:{rng.randint(0, 5)}")
    if rng.random() < 0.2:
        tags.append(f"BC:Z:{rng.choice(('ACGT', 'TTTT'))}")
    if rng.random() < 0.2:
        tags.append(f"BX:Z:{rng.choice(('AAAA', 'CCCC'))}")
    return tags


def random_spec(rng):
    contigs = [(f"c{i + 1}", random_sequence(rng, rng.choice((1000, 1500, 2000, 5000))))
               for i in range(rng.choice((1, 1, 2, 3)))]
    groups = []
    for i in range(rng.choice((0, 1, 1, 2, 3))):
        fields = {}
        if rng.random() < 0.85:
            fields["SM"] = rng.choice(("S1", "S2"))
        if rng.random() < 0.8:
            fields["LB"] = rng.choice(("L1", "L2", "L3"))
        if rng.random() < 0.8:
            fields["PU"] = rng.choice(("P1", "P2", "P3"))
        groups.append((f"rg{i}", fields))
    group_ids = [g for g, _ in groups] + [None] + (["ghost"] if rng.random() < 0.3 else [])
    templates = []
    count = rng.choice((0, 1, 2, 5, 10, 20, 40))
    clusters = [(rng.choice((1101, 1102)), rng.randint(1000, 10_000), rng.randint(1000, 10_000))
                for _ in range(max(1, count // 3))]
    sequences = {}
    index = 0
    for _ in range(count):
        index += 1
        copies = rng.choice((1, 1, 1, 2, 3))
        contig, seq = rng.choice(contigs)
        length1 = rng.choice([n for n in READ_LENGTHS if n <= len(seq) - 10])
        length2 = rng.choice([n for n in READ_LENGTHS if n <= len(seq) - 10])
        pos1 = rng.randint(1, max(1, len(seq) - length1 - 400))
        pos2 = min(len(seq) - length2, pos1 + rng.choice((0, 1, 50, 200, 300, 350)))
        paired = rng.random() < 0.75
        key = (contig, pos1, pos2)
        if key not in sequences:
            s1 = seq[pos1 - 1:pos1 - 1 + length1].upper().replace("N", rng.choice("ACGT"))
            s2 = seq[pos2 - 1:pos2 - 1 + length2].upper().replace("N", rng.choice("ACGT"))
            sequences[key] = (s1, s2)
        s1, s2 = sequences[key]
        for copy in range(copies):
            index += 1
            name = read_name(rng, index, rng.choice(clusters))
            rg = rng.choice(group_ids)
            mutated1 = list(s1)
            for _ in range(rng.choice((0, 0, 1, 3))):
                k = rng.randrange(len(mutated1))
                mutated1[k] = rng.choice("ACGTN")
            mutated1 = "".join(mutated1)
            dup = DUPLICATE if copy > 0 and rng.random() < 0.5 else 0
            qcfail = QCFAIL if rng.random() < 0.1 else 0
            if not paired:
                flag = rng.choice((0, REVERSE)) | dup | qcfail
                if rng.random() < 0.1:
                    flag |= UNMAPPED
                flag |= rng.choice((0, 0, 0, 0, SECONDARY, SUPPLEMENTARY))
                if rng.random() < 0.05:
                    flag |= SECOND
                templates.append(Template(name, [single_record(rng, flag, contig, pos1, mutated1, rg,
                                                               len(seq))]))
                continue
            reverse1 = rng.random() < 0.5
            reverse2 = not reverse1 if rng.random() < 0.8 else reverse1
            unmapped = rng.choice((None, None, None, None, 1, 2, "both"))
            records = pair_records(rng, contig, pos1, pos2, mutated1, s2, reverse1, reverse2, rg,
                                   unmapped, dup | qcfail, len(seq))
            templates.append(Template(name, records))
    return BamSpec(contigs, groups, templates)


def single_record(rng, flag, contig, pos, seq, rg, contig_length):
    seq_text = revcomp(seq) if flag & REVERSE else seq
    if flag & UNMAPPED:
        tags = random_tags(rng, False, None, None, rg, seq_text)
        return "\t".join([str(flag), "*", "0", "0", "*", "*", "0", "0", seq_text,
                          random_quals(rng, len(seq_text)), *tags])
    cigar, _ = random_cigar(rng, len(seq), contig_length - pos + 1)
    if cigar is None:
        cigar = f"{len(seq)}M" if pos + len(seq) - 1 <= contig_length else "*"
        if cigar == "*":
            flag |= UNMAPPED
            return single_record(rng, flag, contig, pos, seq, rg, contig_length)
    tags = random_tags(rng, True, None, None, rg, seq_text)
    if flag & SUPPLEMENTARY and rng.random() < 0.5:
        half = max(1, len(seq_text) // 2)
        seq_text = seq_text[:half]
        cigar = f"{len(seq) - half}H{half}M"
        if pos + half - 1 > contig_length:
            cigar = f"{len(seq) - half}H{contig_length - pos + 1}M{half - (contig_length - pos + 1)}S"
    return "\t".join([str(flag), contig, str(pos), str(rng.choice(MAPQS)), cigar, "*", "0", "0",
                      seq_text, random_quals(rng, len(seq_text)), *tags])


def pair_records(rng, contig, pos1, pos2, s1, s2, reverse1, reverse2, rg, unmapped, extra_flag,
                 contig_length):
    ends = []
    for i, (pos, seq, reverse) in enumerate(((pos1, s1, reverse1), (pos2, s2, reverse2))):
        mapped = unmapped not in ("both", i + 1)
        cigar, span = (random_cigar(rng, len(seq), contig_length - pos + 1) if mapped else ("*", 0))
        if mapped and cigar is None:
            mapped = False
            cigar, span = "*", 0
        ends.append({"pos": pos, "seq": seq, "reverse": reverse, "mapped": mapped, "cigar": cigar,
                     "span": span, "mapq": rng.choice(MAPQS)})
    records = []
    both_mapped = ends[0]["mapped"] and ends[1]["mapped"]
    if both_mapped:
        left = min(ends[0]["pos"], ends[1]["pos"])
        right = max(ends[0]["pos"] + ends[0]["span"], ends[1]["pos"] + ends[1]["span"])
        tlen = right - left
    for i, end in enumerate(ends):
        mate = ends[1 - i]
        flag = PAIRED | (FIRST if i == 0 else SECOND) | extra_flag
        if both_mapped and rng.random() < 0.8:
            flag |= PROPER
        if not end["mapped"]:
            flag |= UNMAPPED
        if not mate["mapped"]:
            flag |= MATE_UNMAPPED
        if end["reverse"]:
            flag |= REVERSE
        if mate["reverse"]:
            flag |= MATE_REVERSE
        seq_text = revcomp(end["seq"]) if end["reverse"] else end["seq"]
        if end["mapped"]:
            rname, pos, mapq, cigar = contig, end["pos"], end["mapq"], end["cigar"]
        elif mate["mapped"]:
            rname, pos, mapq, cigar = contig, mate["pos"], 0, "*"
        else:
            rname, pos, mapq, cigar = "*", 0, 0, "*"
        if mate["mapped"]:
            rnext, pnext = "=" if rname == contig else contig, mate["pos"]
        elif end["mapped"]:
            rnext, pnext = "=", end["pos"]
        else:
            rnext, pnext = "*", 0
        if both_mapped:
            this_tlen = tlen if end["pos"] <= mate["pos"] else -tlen
            if end["pos"] == mate["pos"]:
                this_tlen = tlen if i == 0 else -tlen
            if rng.random() < 0.05:
                this_tlen = 0
        else:
            this_tlen = 0
        tags = random_tags(rng, end["mapped"], mate["cigar"] if mate["mapped"] else None,
                           mate["mapq"] if mate["mapped"] else None, rg, seq_text)
        records.append("\t".join([str(flag), rname, str(pos), str(mapq), cigar, rnext, str(pnext),
                                  str(this_tlen), seq_text, random_quals(rng, len(seq_text)), *tags]))
    return records


# --- running --------------------------------------------------------------------------

def write_inputs(spec, directory):
    """ref.fa (+ .fai, .dict) and example.bam (+ .bai) in `directory`."""
    directory.mkdir(parents=True, exist_ok=True)
    fasta = directory / "ref.fa"
    with open(fasta, "w") as out:
        for name, seq in spec.contigs:
            out.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")
    subprocess.run(["samtools", "faidx", str(fasta)], check=True, capture_output=True)
    with open(directory / "ref.dict", "w") as out:
        out.write("@HD\tVN:1.6\n")
        for name, seq in spec.contigs:
            out.write(f"@SQ\tSN:{name}\tLN:{len(seq)}\n")
    sam = directory / "example.sam"
    sam.write_text(spec.sam_text())
    bam = directory / "example.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True, capture_output=True)
    subprocess.run(["samtools", "index", str(bam)], check=True, capture_output=True)
    return bam, fasta


def random_commands(rng, bam, fasta):
    """Three (tool, arguments) sets for one example; `O` or `M` is added at run time."""
    multiple = [f"R={fasta}", "PROGRAM=CollectGcBiasMetrics", "PROGRAM=CollectQualityYieldMetrics"]
    for level in rng.sample(("SAMPLE", "LIBRARY", "READ_GROUP"), rng.choice((0, 0, 1, 2, 3))):
        multiple.append(f"LEVEL={level}")
    extras = [
        "EXTRA_ARGUMENT=CollectInsertSizeMetrics::MINIMUM_PCT=0",
        "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=300",
        "EXTRA_ARGUMENT=CollectInsertSizeMetrics::INCLUDE_DUPLICATES=true",
        "EXTRA_ARGUMENT=CollectGcBiasMetrics::SCAN_WINDOW_SIZE=50",
        "EXTRA_ARGUMENT=CollectGcBiasMetrics::ALSO_IGNORE_DUPLICATES=true",
        "EXTRA_ARGUMENT=CollectGcBiasMetrics::IS_BISULFITE_SEQUENCED=true",
        "EXTRA_ARGUMENT=CollectAlignmentSummaryMetrics::MAX_INSERT_SIZE=300",
        "EXTRA_ARGUMENT=CollectAlignmentSummaryMetrics::IS_BISULFITE_SEQUENCED=true",
        "EXTRA_ARGUMENT=QualityScoreDistribution::INCLUDE_NO_CALLS=true",
        "EXTRA_ARGUMENT=MeanQualityByCycle::ALIGNED_READS_ONLY=true",
        "EXTRA_ARGUMENT=CollectBaseDistributionByCycle::PF_READS_ONLY=true",
        "EXTRA_ARGUMENT=CollectQualityYieldMetrics::OQ=false",
    ]
    multiple += rng.sample(extras, rng.choice((0, 1, 2, 3)))
    if rng.random() < 0.1:
        multiple.append(f"STOP_AFTER={rng.choice((1, 3, 17))}")
    duplicate = []
    complexity = []
    if rng.random() < 0.5:
        complexity.append(f"MIN_IDENTICAL_BASES={rng.choice((3, 5, 8))}")
    if rng.random() < 0.5:
        complexity.append(f"MAX_DIFF_RATE={rng.choice((0, 0.03, 0.1))}")
    if rng.random() < 0.5:
        complexity.append(f"MIN_MEAN_QUALITY={rng.choice((0, 20, 30))}")
    if rng.random() < 0.5:
        complexity.append(f"OPTICAL_DUPLICATE_PIXEL_DISTANCE={rng.choice((10, 100, 2500))}")
    if rng.random() < 0.3:
        complexity.append(f"MIN_GROUP_COUNT={rng.choice((1, 2))}")
    if rng.random() < 0.2:
        complexity.append("READ_NAME_REGEX=null")
    if rng.random() < 0.2:
        complexity.append(f"MAX_READ_LENGTH={rng.choice((20, 50))}")
    if rng.random() < 0.2:
        complexity.append(f"{rng.choice(('BARCODE_TAG=BC', 'READ_ONE_BARCODE_TAG=BX'))}")
    return [("CollectMultipleMetrics", "O", multiple), ("CollectDuplicateMetrics", "M", duplicate),
            ("EstimateLibraryComplexity", "O", complexity)]


def cache_key(bam, tool, args):
    digest = hashlib.sha256(bam.read_bytes())
    digest.update(json.dumps([tool, args]).encode())
    digest.update(b"picard-3.5.0")
    return digest.hexdigest()


def run_picard(tool, bam, key, args, out_dir, cache):
    """Picard's files for one command, from the cache when seen before.
    Returns {suffix: bytes} or None when Picard failed."""
    entry = cache / cache_key(bam, tool, args)
    if not (entry / "status.json").exists():
        work = out_dir / "picard_work"
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        result = subprocess.run(["picard", tool, f"I={bam}", f"{key}={work / 'out'}", SILENT, *args,
                                 "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(),
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        entry.mkdir(parents=True, exist_ok=True)
        for path in work.glob("out*"):
            shutil.copy(path, entry / path.name)
        (entry / "status.json").write_text(json.dumps({"exit": result.returncode}))
    status = json.loads((entry / "status.json").read_text())
    if status["exit"] != 0:
        return None
    return {p.name[3:]: p.read_bytes() for p in entry.glob("out*")}


def run_bammetrics(tool, bam, key, args, out_dir):
    work = out_dir / "bammetrics_work"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    result = subprocess.run([str(BAMMETRICS), tool, f"I={bam}", f"{key}={work / 'out'}", *args],
                            capture_output=True, text=True)
    if result.returncode != 0:
        return None, result.stderr
    return {p.name[3:]: p.read_bytes() for p in work.glob("out*")}, ""


def compare_command(tool, key, args, bam, out_dir, cache):
    """None when the two tools agree, 'both-err' when both failed, else a
    description of the mismatch."""
    picard = run_picard(tool, bam, key, args, out_dir, cache)
    ours, stderr = run_bammetrics(tool, bam, key, args, out_dir)
    if picard is None and ours is None:
        return "both-err"
    if picard is None:
        return f"{tool}: Picard failed, bammetrics wrote {sorted(ours)}"
    if ours is None:
        return f"{tool}: bammetrics failed ({stderr.strip()[:200]}), Picard wrote {sorted(picard)}"
    if set(picard) != set(ours):
        return f"{tool}: files Picard {sorted(picard)} bammetrics {sorted(ours)}"
    for suffix in sorted(picard):
        if strip_headers(picard[suffix]) != strip_headers(ours[suffix]):
            return f"{tool}: file {suffix!r} differs"
    return None


def strip_headers(data):
    lines = data.split(b"\n")
    kept, skip = [], False
    for line in lines:
        if skip:
            skip = False
            continue
        if line == b"## htsjdk.samtools.metrics.StringHeader":
            skip = True
            continue
        kept.append(line)
    return b"\n".join(kept)


def run_example(spec, commands_seed, out_dir, cache):
    """Every command of one example. Returns (mismatches, both_err_count)."""
    bam, fasta = write_inputs(spec, out_dir)
    rng = random.Random(commands_seed)
    mismatches, both_err = [], 0
    for tool, key, args in random_commands(rng, bam, fasta):
        args = [a.replace("R=" + str(fasta), f"R={fasta}") for a in args]
        outcome = compare_command(tool, key, args, bam, out_dir, cache)
        if outcome == "both-err":
            both_err += 1
        elif outcome is not None:
            mismatches.append((tool, key, args, outcome))
    return mismatches, both_err


# --- shrinking ------------------------------------------------------------------------

def shrink(spec, failing, out_dir, cache):
    """ddmin over the templates: the smallest subset for which `failing`
    still reports a mismatch (any mismatch of the same tool and arguments)."""
    tool, key, args, _ = failing

    def fails(templates):
        reduced = spec.subset(templates)
        bam, fasta = write_inputs(reduced, out_dir / "shrink")
        args_here = [a if not a.startswith("R=") else f"R={fasta}" for a in args]
        outcome = compare_command(tool, key, args_here, bam, out_dir / "shrink", cache)
        return outcome not in (None, "both-err")

    templates = list(spec.templates)
    n = 2
    while len(templates) >= 2:
        chunk = max(1, len(templates) // n)
        subsets = [templates[i:i + chunk] for i in range(0, len(templates), chunk)]
        reduced = False
        for subset in subsets:
            complement = [t for t in templates if t not in subset]
            if complement and fails(complement):
                templates, n, reduced = complement, max(n - 1, 2), True
                break
        if not reduced:
            if n >= len(templates):
                break
            n = min(n * 2, len(templates))
    return spec.subset(templates)


def save_regression(spec, failing, picard_files):
    REGRESSIONS.mkdir(exist_ok=True)
    text = spec.sam_text()
    sha8 = hashlib.sha256(text.encode()).hexdigest()[:8]
    (REGRESSIONS / f"{sha8}.sam").write_text(text)
    with open(REGRESSIONS / f"{sha8}.fa", "w") as out:
        for name, seq in spec.contigs:
            out.write(f">{name}\n{seq}\n")
    tool, key, args, outcome = failing
    (REGRESSIONS / f"{sha8}.args").write_text(json.dumps(
        {"tool": tool, "output_key": key, "args": args, "outcome": outcome,
         "reference": [[name, len(seq)] for name, seq in spec.contigs]}, indent=2) + "\n")
    with open(REGRESSIONS / f"{sha8}.picard.txt", "w") as out:
        for suffix, data in sorted((picard_files or {}).items()):
            out.write(f"### {suffix}\n")
            out.write(data.decode(errors="replace"))
    return sha8


# --- main -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--n", type=int, default=200)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--out", type=Path, default=ROOT / "target" / "fuzz")
    parser.add_argument("--no-shrink", action="store_true")
    options = parser.parse_args()
    cache = options.out / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    started = time.time()
    failures = []
    both_err_total = 0
    templates_total = 0

    def one(index):
        rng = random.Random(options.seed * 1_000_003 + index)
        spec = random_spec(rng)
        out_dir = options.out / f"seed{options.seed}" / f"{index:05d}"
        mismatches, both_err = run_example(spec, options.seed * 7_919 + index, out_dir, cache)
        if not mismatches:
            shutil.rmtree(out_dir, ignore_errors=True)
        return index, spec, mismatches, both_err

    with ThreadPoolExecutor(max_workers=options.jobs) as pool:
        for index, spec, mismatches, both_err in pool.map(one, range(options.n)):
            both_err_total += both_err
            templates_total += len(spec.templates)
            if mismatches:
                failures.append((index, spec, mismatches))
                print(f"example {index}: {mismatches[0][3]}", flush=True)
    for index, spec, mismatches in failures:
        failing = mismatches[0]
        out_dir = options.out / f"seed{options.seed}" / f"{index:05d}"
        reduced = spec if options.no_shrink else shrink(spec, failing, out_dir, cache)
        bam, fasta = write_inputs(reduced, out_dir / "final")
        tool, key, args, _ = failing
        args = [a if not a.startswith("R=") else f"R={fasta}" for a in args]
        picard_files = run_picard(tool, bam, key, args, out_dir / "final", cache)
        sha8 = save_regression(reduced, failing, picard_files)
        print(f"example {index}: shrunk {len(spec.templates)} -> {len(reduced.templates)} templates, "
              f"saved tests/concordance/regressions/{sha8}.sam")
    elapsed = time.time() - started
    print(f"seed {options.seed}: {options.n} examples, {templates_total} templates, "
          f"{len(failures)} failures, {both_err_total} both-err commands, {elapsed:.0f} s")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

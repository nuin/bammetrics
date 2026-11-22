"""Inputs for FixMateInformation and FilterSamReads (wave W13), cut from the foundation
corpus and seeded, so every run writes the same bytes.

`write(corpus_dir, out)` writes:

  damaged.bam             the corpus with one record kept per (name, first or second,
                          primary): Picard fails on a name with two primary first reads.
                          Then, seeded, mate fields damaged on a third of the records: `MC`
                          removed or wrong, `MQ` removed, the mate start moved, `TLEN` 0 or
                          negated, the mate-reverse flag flipped, the mate-unmapped flag
                          set on a mapped mate or cleared on an unmapped one; one end
                          dropped from 3% of the pairs (orphans). Coordinate order, the
                          corpus header (`SO:coordinate`).
  damaged.queryname.bam   the same records in htsjdk's queryname order (bammetrics SortSam,
                          which equals Picard's).
  damaged.natural.bam     the same records in samtools' natural name order under
                          `SO:queryname`: Picard takes it as sorted.
  primary.queryname.bam   the primary records of whole pairs and the unpaired records, in
                          queryname order: the pair filters of FilterSamReads read it.
  split_N.bam             the damaged records of read groups `N`, `N+3`, ... (N 0 to 2),
                          each file's read groups renamed `g0`, `g1`, ... so every ID
                          collides across files with other fields; coordinate order. Also
                          `split_N.queryname.bam`.
  tagged.bam              the corpus records with tags of every type: `XZ:Z`, `XA:A`,
                          `XI:i` (several widths), `XF:f`, `XH:H`, `XB:B`, values cycling so
                          a value repeats on many records.
  names_*.txt             read lists: corpus names one per line, with words after them,
                          CRLF and CR line ends, blank and space-only lines, a leading
                          space, a `#` line, names absent from the BAM, repeated names, a
                          gzip copy, and names beyond ASCII (`names_latin1.txt`, for
                          `latin1.sam`, UTF-8 SAM text whose names htsjdk holds as the
                          characters E9 and FC; a name above U+00FF matches none).
"""
import array
import gzip
import random
import subprocess
import sys
from pathlib import Path

import pysam

SEED = 20260131
TAG_VALUES = ["alpha", "beta", "gamma", "delta"]


def run(*cmd):
    subprocess.run([str(c) for c in cmd], check=True)


def bammetrics():
    root = Path(__file__).resolve().parents[2]
    return root / "target" / "release" / "bammetrics"


def one_primary_per_end(records):
    """Keep the first primary first and second record of each name, every other record."""
    seen, kept = set(), []
    for record in records:
        if not record.is_secondary and not record.is_supplementary and record.is_paired:
            end = (record.query_name, record.is_read1, record.is_read2)
            if end in seen:
                continue
            seen.add(end)
        kept.append(record)
    return kept


def damage(record, rng):
    """One of the mate fields of `record` made wrong."""
    kind = rng.randrange(8)
    if kind == 0 and record.has_tag("MC"):
        record.set_tag("MC", None)
    elif kind == 1:
        record.set_tag("MC", "7M3S")
    elif kind == 2 and record.has_tag("MQ"):
        record.set_tag("MQ", None)
    elif kind == 3:
        record.next_reference_start = max(0, record.next_reference_start + rng.randint(-50, 50))
    elif kind == 4:
        record.template_length = 0
    elif kind == 5:
        record.template_length = -record.template_length
    elif kind == 6 and record.is_paired:
        record.mate_is_reverse = not record.mate_is_reverse
    elif kind == 7 and record.is_paired:
        record.mate_is_unmapped = not record.mate_is_unmapped


def write_records(path, header, records):
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for record in records:
            out.write(record)


def queryname(source, out):
    """htsjdk's queryname order, through bammetrics SortSam (equal to Picard's)."""
    run(bammetrics(), "SortSam", f"I={source}", f"O={out}", "SORT_ORDER=queryname")


def write_damaged(corpus_bam, out, rng):
    with pysam.AlignmentFile(str(corpus_bam)) as bam:
        header = bam.header.to_dict()
        records = one_primary_per_end(list(bam))
    names = sorted({r.query_name for r in records if r.is_paired})
    orphans = {n: rng.choice((True, False)) for n in rng.sample(names, len(names) * 3 // 100)}
    kept = []
    for record in records:
        if record.query_name in orphans and record.is_paired and \
                record.is_read1 == orphans[record.query_name] and not record.is_secondary \
                and not record.is_supplementary:
            continue
        if rng.random() < 0.33:
            damage(record, rng)
        kept.append(record)
    write_records(out / "damaged.bam", header, kept)
    queryname(out / "damaged.bam", out / "damaged.queryname.bam")
    run("samtools", "sort", "-n", "--no-PG", "-o", out / "damaged.natural.bam",
        out / "damaged.bam")
    pairs = {}
    for record in kept:
        if not record.is_secondary and not record.is_supplementary and record.is_paired:
            pairs.setdefault(record.query_name, []).append(record)
    whole = {n for n, ends in pairs.items() if len(ends) == 2}
    primary = [r for r in kept if not r.is_secondary and not r.is_supplementary
               and (not r.is_paired or r.query_name in whole)]
    write_records(out / "primary.bam", header, primary)
    queryname(out / "primary.bam", out / "primary.queryname.bam")
    return header, kept


def write_split(header, records, out):
    groups = [rg["ID"] for rg in header.get("RG", [])]
    for part in range(3):
        mine = groups[part::3]
        renamed = {old: f"g{i}" for i, old in enumerate(mine)}
        part_header = dict(header)
        part_header["RG"] = [dict(rg, ID=renamed[rg["ID"]]) for rg in header["RG"]
                             if rg["ID"] in renamed]
        kept = []
        for record in records:
            group = record.get_tag("RG") if record.has_tag("RG") else None
            if group in renamed:
                copy = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(part_header))
                copy.query_name = record.query_name
                copy.flag = record.flag
                copy.reference_id = record.reference_id
                copy.reference_start = record.reference_start
                copy.mapping_quality = record.mapping_quality
                copy.cigartuples = record.cigartuples
                copy.next_reference_id = record.next_reference_id
                copy.next_reference_start = record.next_reference_start
                copy.template_length = record.template_length
                copy.query_sequence = record.query_sequence
                copy.query_qualities = record.query_qualities
                copy.set_tags([(t, v) for t, v in record.get_tags() if t != "RG"]
                              + [("RG", renamed[group])])
                kept.append(copy)
        path = out / f"split_{part}.bam"
        write_records(path, part_header, kept)
        queryname(path, out / f"split_{part}.queryname.bam")


def write_tagged(corpus_bam, out, rng):
    with pysam.AlignmentFile(str(corpus_bam)) as bam:
        header = bam.header.to_dict()
        records = list(bam)
    for i, record in enumerate(records):
        record.set_tag("XZ", TAG_VALUES[i % 4], "Z")
        record.set_tag("XA", "ab"[i % 2], "A")
        record.set_tag("XI", [7, 300, 70000, -2][i % 4], "i")
        record.set_tag("XF", 1.5 * (i % 3), "f")
        record.set_tag("XH", "1AE3" if i % 2 else "FF00", "H")
        record.set_tag("XB", array.array("h", [1, 2, i % 5]))
        if i % 7 == 0:
            record.set_tag("XZ", None)
    write_records(out / "tagged.bam", header, records)


def write_lists(corpus_bam, out, rng):
    with pysam.AlignmentFile(str(corpus_bam)) as bam:
        names = sorted({r.query_name for r in bam})
    picked = rng.sample(names, len(names) // 4)
    (out / "names_plain.txt").write_text("".join(f"{n}\n" for n in picked))
    lines = []
    for i, name in enumerate(picked):
        lines.append([f"{name}", f"{name}\tafter a tab", f"{name} words after",
                      f" {name}", "", "   ", f"#{name}", f"{name}\r", f"absent_{i}"][i % 9])
    lines += picked[:20]
    (out / "names_kinds.txt").write_text("\n".join(lines) + "\n")
    (out / "names_crlf.txt").write_bytes("\r\n".join(picked[:200]).encode() + b"\r\n")
    (out / "names_cr.txt").write_bytes("\r".join(picked[:200]).encode())
    with gzip.open(out / "names_plain.txt.gz", "wb") as handle:
        handle.write("".join(f"{n}\n" for n in picked).encode())
    (out / "names_absent.txt").write_text("nobody\nnot_here\n")
    (out / "names_latin1.txt").write_text("café\nmüller x\nĀbc\n", "utf-8")
    sam = ("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000\n"
           "café\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
           "müller\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
           "plain\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")
    (out / "latin1.sam").write_bytes(sam.encode("utf-8"))


def write(corpus_dir, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    corpus_bam = Path(corpus_dir) / "corpus.bam"
    header, kept = write_damaged(corpus_bam, out, rng)
    write_split(header, kept, out)
    write_tagged(corpus_bam, out, rng)
    write_lists(corpus_bam, out, rng)
    return out


if __name__ == "__main__":
    write(sys.argv[1], sys.argv[2])

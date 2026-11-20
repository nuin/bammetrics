"""Hand-built BAMs for CollectQualityYieldMetricsSNVQ and
CollectUmiPrevalenceMetrics (docs/07 wave W11), written with pysam from a
seed.

`build_all(out)` writes a reference of its own (u1 60 kb, u2 20 kb; ref.fa
with .fai) and one BAM per entry of the returned dict, coordinate sorted
and indexed unless the name says otherwise.

SNVQ inputs carry the four SNV quality tags qa, qc, qg and qt (Phred+33,
one character per base, any value from ! to ~):

* `snvq_main`: 600 records, paired and unpaired, forward and reverse,
  mapped and unmapped; bases A, C, G, T, N, IUPAC codes and `=`; read
  lengths 1 to 151 and a few with no bases at all (SEQ `*`, empty tags);
  qualities 0 to 45 and some 93; PF and QC-failed records; secondary and
  supplementary records; tags longer than the read; an `XQ` alternate
  quality tag on every record.
* `snvq_short_unread`: a qa tag shorter than an all-A read, never read.
* `snvq_secondary_bare`: a secondary and a supplementary record without
  any SNVQ tag, which pass unless their INCLUDE option is set.
* `snvq_empty`: a header and no record.
* The failing inputs, one defect each: `snvq_missing_tag` (no qt),
  `snvq_short_read_tag` (qc shorter than a read whose late base needs it),
  `snvq_space` (a space in qg where it is read), `snvq_int_tag`
  (`qa:i:30`), `snvq_no_qual` (bases without QUAL), `snvq_missing_alt`
  (one record without XQ, failing under ALTERNATE_QUALITY_ATTRIBUTE=XQ),
  `snvq_bad_alt` (a space in XQ).

UMI inputs carry `RX` (UMI) and `BQ` (its base qualities, a space between
the two halves of a duplex UMI) on paired reads with `MC`:

* `umi_main`: two read groups in two libraries; duplicate sets of 1 to 50
  pairs at one position and orientation, with UMIs drawn from a few true
  ones plus copies with one or two errors; BQ above and below 30, missing
  and empty; MAPQ 0 to 60; soft clips that move the unclipped 5' end;
  fragments and mate-unmapped reads beside the pairs, both strands, so
  FILTER_UNPAIRED_READS=false folds them into the pairs' sets; unmapped,
  secondary, supplementary, QC-failed and duplicate-flagged records;
  records without RX; pairs across the two contigs.
* `umi_unsorted` and `umi_queryname`: the same records, shuffled under
  `SO:unsorted` and name sorted. Picard sorts them itself.
* `umi_one_no_mc`: one kept record, without MC; nothing compares it, so
  both tools write one set of one UMI.
* `umi_empty`: a header and no record.
* The failing inputs: `umi_no_mc` (pairs without MC), `umi_int_rx`
  (`RX:i:5` on a kept record), `umi_int_bq` (`BQ:i:30`), `umi_bad_bq`
  (a BQ character past `~`, written as raw BAM bytes by pysam).

Repeated tags, as htsjdk's `BinaryTagCodec.readTags` lists them
(`bam_raw::HtsjdkTags`): `snvq_repeated` (qa, qc and XQ given twice, the
last value read), `snvq_repeated_lost` (a second qa replaces the append
point, so qc, qg and qt are lost and both tools fail), `umi_repeated` (RX
given twice, the last UMI counted; RG given twice before RX, so RX is lost
and the read left out).
"""
import random

import pysam

CONTIGS = [("u1", 60_000), ("u2", 20_000)]
GROUPS = [{"ID": "ra", "SM": "S", "LB": "La", "PU": "P1"},
          {"ID": "rb", "SM": "S", "LB": "Lb", "PU": "P2"}]
PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
FIRST, SECOND, SECONDARY, QCFAIL, DUPLICATE, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x200, 0x400, 0x800
SNVQ_TAGS = ("qa", "qc", "qg", "qt")
BASES = "ACGT"


def write_reference(out, seed=71):
    rng = random.Random(seed)
    ref = {}
    with open(out / "ref.fa", "w") as fasta:
        for name, length in CONTIGS:
            seq = "".join(rng.choice(BASES) for _ in range(length))
            ref[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, length, 70):
                fasta.write(seq[i:i + 70] + "\n")
    pysam.faidx(str(out / "ref.fa"))
    return ref


def header(order="coordinate", groups=GROUPS):
    return {"HD": {"VN": "1.6", "SO": order},
            "SQ": [{"SN": n, "LN": length} for n, length in CONTIGS],
            "RG": [dict(g) for g in groups]}


def phred_text(rng, length, low=0, high=93):
    return "".join(chr(33 + rng.randint(low, high)) for _ in range(length))


def segment(name, flag, contig, start, bases, quals, cigar=None, mapq=60, mate=None, tags=()):
    """A record; `mate` is (contig, start) or None. `start` is 0-based."""
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    contig_id = [n for n, _ in CONTIGS].index(contig) if contig else -1
    a.reference_id = contig_id
    a.reference_start = start if contig else -1
    unmapped = flag & UNMAPPED
    a.mapping_quality = 0 if unmapped else mapq
    a.cigartuples = None if unmapped or not bases else (cigar or [(0, len(bases))])
    if bases:
        a.query_sequence = bases
        if quals is not None:
            a.query_qualities = pysam.qualitystring_to_array(quals)
    if mate:
        a.next_reference_id = [n for n, _ in CONTIGS].index(mate[0]) if mate[0] else -1
        a.next_reference_start = mate[1]
    else:
        a.next_reference_id = -1
        a.next_reference_start = -1
    for tag, value, *kind in tags:
        if value is not None:
            a.set_tag(tag, value, *kind)
    return a


def write(out, name, records, order="coordinate", groups=GROUPS):
    """Write `records` as `name.bam`, sorted and indexed for a coordinate order."""
    path = out / f"{name}.bam"
    raw = out / f"{name}.unsorted.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header(order, groups)) as bam:
        for r in records:
            bam.write(r)
    if order == "coordinate":
        pysam.sort("-o", str(path), str(raw))
        pysam.index(str(path))
    elif order == "queryname":
        pysam.sort("-n", "-o", str(path), str(raw))
    else:
        raw.replace(path)
    raw.unlink(missing_ok=True)
    return path


# --- SNVQ ----------------------------------------------------------------------

def snvq_bases(rng, length):
    pool = "ACGTACGTACGTACGTNNRY="
    return "".join(rng.choice(pool) for _ in range(length))


def snvq_quals(rng, length):
    return "".join(chr(33 + rng.choice((0, 2, 10, 19, 20, 25, 29, 30, 35, 39, 40, 41, 45, 93)))
                   for _ in range(length))


def snvq_tags(rng, length, extra=0):
    tags = [(tag, phred_text(rng, length + extra, *rng.choice(((0, 93), (15, 50), (35, 60)))))
            for tag in SNVQ_TAGS]
    return tags + [("XQ", snvq_quals(rng, length))]


def snvq_main(rng, ref):
    records = []
    for i in range(600):
        length = rng.choice((1, 2, 17, 50, 75, 100, 151))
        bases = snvq_bases(rng, length)
        flag = rng.choice((0, 0, REVERSE, PAIRED | FIRST, PAIRED | SECOND | REVERSE,
                           PAIRED | FIRST | MATE_UNMAPPED, UNMAPPED, UNMAPPED | PAIRED | SECOND))
        if rng.random() < 0.15:
            flag |= QCFAIL
        if not flag & UNMAPPED:
            if rng.random() < 0.05:
                flag |= SECONDARY
            elif rng.random() < 0.05:
                flag |= SUPPLEMENTARY
        if rng.random() < 0.05:
            flag |= DUPLICATE
        contig, size = rng.choice(CONTIGS)
        start = rng.randint(0, size - 200)
        mate = (contig, start + 100) if flag & PAIRED else None
        if flag & UNMAPPED and not flag & PAIRED and rng.random() < 0.5:
            contig = None
        extra = rng.choice((0, 0, 0, 3))
        records.append(segment(f"s{i}", flag, contig, start, bases, snvq_quals(rng, length),
                               mate=mate,
                               tags=[*snvq_tags(rng, length, extra), ("RG", rng.choice(GROUPS)["ID"])]))
    for i in range(3):
        # No bases: SEQ and QUAL `*`, empty tags.
        records.append(segment(f"empty{i}", UNMAPPED, None, 0, "", None,
                               tags=[(t, "") for t in (*SNVQ_TAGS, "XQ")]))
    return records


def snvq_variants(rng):
    """The single-defect inputs: name -> records."""
    def read(name, bases, tags=None, flag=0, quals=True, drop=()):
        tags = dict(tags or [(t, "I" * len(bases)) for t in (*SNVQ_TAGS, "XQ")])
        for tag in drop:
            tags.pop(tag)
        return segment(name, flag, "u1", 100, bases, "I" * len(bases) if quals else None,
                       tags=list(tags.items()))
    good = [read(f"g{i}", snvq_bases(rng, 30)) for i in range(5)]
    short_qa = {"qa": "II", "qc": "I" * 10, "qg": "I" * 10, "qt": "I" * 10, "XQ": "I" * 10}
    short_qc = {"qa": "I" * 10, "qc": "IIII", "qg": "I" * 10, "qt": "I" * 10, "XQ": "I" * 10}
    space = {"qa": "I" * 10, "qc": "I" * 10, "qg": "III II" + "I" * 4, "qt": "I" * 10,
             "XQ": "I" * 10}
    bad_alt = {**{t: "I" * 10 for t in SNVQ_TAGS}, "XQ": "II II" + "I" * 5}
    int_qa = segment("int", 0, "u1", 100, "ACGTACGTAC", "I" * 10,
                     tags=[("qa", 30, "i"), *[(t, "I" * 10) for t in SNVQ_TAGS[1:]]])
    return {
        "snvq_short_unread": [*good, read("allA", "A" * 10, short_qa)],
        "snvq_secondary_bare": [*good, read("sec", "ACGT", flag=SECONDARY, drop=(*SNVQ_TAGS, "XQ")),
                                read("sup", "ACGT", flag=SUPPLEMENTARY, drop=(*SNVQ_TAGS, "XQ"))],
        "snvq_empty": [],
        "snvq_missing_tag": [*good, read("noqt", "ACGT", drop=("qt",))],
        "snvq_short_read_tag": [*good, read("shortqc", "AAAAAAAAGA", short_qc)],
        "snvq_space": [*good, read("space", "ACGTACGTAC", space)],
        "snvq_int_tag": [*good, int_qa],
        "snvq_no_qual": [*good, read("noqual", "ACGT", quals=False)],
        "snvq_missing_alt": [*good, read("noxq", "ACGT", drop=("XQ",))],
        "snvq_bad_alt": [*good, read("badxq", "ACGTACGTAC", bad_alt)],
    }


# --- UMI -----------------------------------------------------------------------

def mutate_umi(rng, umi, errors):
    out = list(umi)
    for _ in range(errors):
        i = rng.randrange(len(out))
        if out[i] != "-":
            out[i] = rng.choice([b for b in BASES if b != out[i]])
    return "".join(out)


def bq_for(rng, umi):
    """BQ: high, low somewhere, empty or absent, a space where RX has `-`."""
    kind = rng.choice(("high", "high", "low", "low", "absent", "empty"))
    if kind == "absent":
        return None
    if kind == "empty":
        return ""
    quals = [rng.randint(30, 41) for _ in umi.replace("-", "")]
    if kind == "low":
        quals[rng.randrange(len(quals))] = rng.randint(2, 29)
    text = "".join(chr(33 + q) for q in quals)
    if "-" in umi:
        half = umi.index("-")
        text = text[:half] + " " + text[half:]
    return text


def pair_records(rng, name, contig, pos1, pos2, flags1, rg, umi, mapqs=(60, 60), clip=0,
                 length=40, mate_contig=None, mc=True, extra_flags=0):
    """Read 1 at pos1 and read 2 at pos2 (0-based); flags1 says the strand of read 1."""
    mate_contig = mate_contig or contig
    reverse1 = bool(flags1 & REVERSE)
    cigar1 = [(4, clip), (0, length - clip)] if clip else [(0, length)]
    cigar2 = [(0, length)]
    text = {1: "".join(f"{n}{'SM'[op == 0]}" for op, n in cigar1), 2: f"{length}M"}
    bq = bq_for(rng, umi)
    out = []
    for which, contig_i, pos, mate_pos, cigar, mapq, reverse, mate_reverse in (
            (FIRST, contig, pos1, (mate_contig, pos2), cigar1, mapqs[0], reverse1, not reverse1),
            (SECOND, mate_contig, pos2, (contig, pos1), cigar2, mapqs[1], not reverse1, reverse1)):
        flag = PAIRED | PROPER | which | extra_flags
        flag |= REVERSE if reverse else 0
        flag |= MATE_REVERSE if mate_reverse else 0
        bases = "".join(rng.choice(BASES) for _ in range(length))
        mate_cigar = text[2] if which == FIRST else text[1]
        out.append(segment(name, flag, contig_i, pos, bases, "I" * length, cigar=cigar,
                           mapq=mapq, mate=mate_pos,
                           tags=[("RG", rg), ("MC", mate_cigar if mc else None), ("RX", umi),
                                 ("BQ", bq)]))
    return out


def umi_main(rng):
    records = []
    truths = ["".join(rng.choice(BASES) for _ in range(8)) for _ in range(12)]
    duplex = [f"{t[:4]}-{t[4:]}" for t in truths[:4]]
    pos = 1_000
    for s in range(70):
        contig = "u1" if s % 5 else "u2"
        size = rng.choice((1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 13, 20, 30, 50))
        rg = rng.choice(GROUPS)["ID"]
        pool = rng.sample(truths + duplex, rng.choice((1, 1, 2, 3, 5)))
        pos = pos + rng.randint(150, 700) if contig == "u1" else rng.randint(500, 18_000)
        insert = rng.randint(150, 400)
        flags1 = rng.choice((0, REVERSE))
        for c in range(size):
            umi = mutate_umi(rng, rng.choice(pool), rng.choice((0, 0, 0, 0, 1, 2)))
            mapqs = (rng.choice((60, 60, 60, 30, 29, 0)), rng.choice((60, 60, 40, 10)))
            clip = rng.choice((0, 0, 0, 3))
            extra = QCFAIL if rng.random() < 0.05 else (DUPLICATE if c and rng.random() < 0.5 else 0)
            mate_contig = "u2" if contig == "u1" and s % 7 == 0 else contig
            pos2 = pos + insert if mate_contig == contig else rng.randint(100, 15_000)
            records += pair_records(rng, f"set{s}_{c}", contig, pos, pos2, flags1, rg, umi,
                                    mapqs, clip, mate_contig=mate_contig, extra_flags=extra)
        # Fragments and mate-unmapped reads at the set's 5' end, both strands.
        for f in range(rng.choice((0, 0, 1, 2))):
            umi = rng.choice(pool)
            strand = rng.choice((0, REVERSE))
            start = pos
            flag = strand if rng.random() < 0.5 else PAIRED | FIRST | MATE_UNMAPPED | strand
            mate = (contig, start) if flag & PAIRED else None
            records.append(segment(f"frag{s}_{f}", flag, contig, start,
                                   "".join(rng.choice(BASES) for _ in range(40)), "I" * 40,
                                   mate=mate, tags=[("RG", rg), ("RX", umi), ("BQ", bq_for(rng, umi))]))
            if flag & PAIRED:
                records.append(segment(f"frag{s}_{f}", PAIRED | SECOND | UNMAPPED, contig, start,
                                       "".join(rng.choice(BASES) for _ in range(40)), "I" * 40,
                                       mate=(contig, start), tags=[("RG", rg), ("RX", umi)]))
    # Records every filter drops, and records without RX.
    for i in range(30):
        kind = rng.choice((UNMAPPED | PAIRED | FIRST | MATE_UNMAPPED, SECONDARY | PAIRED | FIRST,
                           SUPPLEMENTARY | PAIRED | FIRST, PAIRED | FIRST))
        start = rng.randint(1_000, 50_000)
        tags = [("RG", "ra"), ("MC", "40M"), ("RX", None if kind == PAIRED | FIRST else "ACGTACGT")]
        records.append(segment(f"drop{i}", kind | MATE_REVERSE, "u1", start,
                               "".join(rng.choice(BASES) for _ in range(40)), "I" * 40,
                               mate=("u1", start + 200), tags=tags))
    return records


def umi_variants(rng):
    def pair(name, umi="ACGTACGT", mc=True, rx_kind=None, bq=None, bq_kind=None):
        recs = pair_records(rng, name, "u1", 500, 700, 0, "ra", umi, mc=mc)
        for r in recs:
            r.set_tag("BQ", None)
            if rx_kind:
                r.set_tag("RX", 5, rx_kind)
            if bq is not None:
                r.set_tag("BQ", bq, bq_kind)
        return recs
    good = [r for i in range(4) for r in pair(f"g{i}")]
    lone = segment("lone", PAIRED | FIRST | MATE_REVERSE, "u1", 300, "A" * 40, "I" * 40,
                   mate=("u1", 500), tags=[("RG", "ra"), ("RX", "ACGT")])
    return {
        "umi_one_no_mc": [lone],
        "umi_empty": [],
        "umi_no_mc": [*good, *pair("nomc", mc=False)],
        "umi_int_rx": [*good, *pair("intrx", rx_kind="i")],
        "umi_int_bq": [*good, *pair("intbq", bq=30, bq_kind="i")],
        "umi_bad_bq": [*good, *pair("badbq", bq="++\x7f+", bq_kind="Z")],
    }


def with_tags(record, tags):
    """`record` with its tags replaced by `tags` in this order, repeats kept."""
    record.set_tags([(tag, value, "Z") for tag, value in tags])
    return record


def repeated_tag_variants(rng):
    """Records with a repeated tag: name -> records."""
    def snvq(name, tags):
        return with_tags(segment(name, 0, "u1", 100, "ACGT", "IIII"), tags)
    good = [snvq(f"g{i}", [(t, "IIII") for t in (*SNVQ_TAGS, "XQ")]) for i in range(3)]
    low, high = "!!!!", "IIII"
    repeated = [
        snvq("appended", [("qa", low), ("qc", high), ("qg", high), ("qt", high), ("XQ", low),
                          ("qa", high), ("XQ", high)]),
        snvq("inserted", [("qt", high), ("qg", low), ("qc", low), ("qa", low), ("qa", high),
                          ("qc", high), ("XQ", high)]),
    ]
    lost = snvq("lost", [("qa", low), ("qa", high), ("qc", high), ("qg", high), ("qt", high),
                         ("XQ", high)])
    umi = []
    for i, (first, last) in enumerate((("AAAA", "CCCC"), ("GGGG", "CCCC"), ("TTTT", "CCCC"))):
        for record in pair_records(rng, f"rep{i}", "u1", 800, 1_000, 0, "ra", first):
            umi.append(with_tags(record, [("RG", "ra"), ("MC", "40M"), ("RX", first),
                                          ("RX", last)]))
    for record in pair_records(rng, "lostrx", "u1", 800, 1_000, 0, "ra", "GGGG"):
        umi.append(with_tags(record, [("RG", "ra"), ("MC", "40M"), ("RG", "ra"),
                                      ("RX", "GGGG")]))
    return {
        "snvq_repeated": [*good, *repeated],
        "snvq_repeated_lost": [*good, lost],
        "umi_repeated": umi,
    }


def build_all(out, seed=97):
    """Write every input under `out`; returns name -> path, plus `ref`."""
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    ref = write_reference(out)
    built = {"ref": out / "ref.fa"}
    built["snvq_main"] = write(out, "snvq_main", snvq_main(rng, ref))
    for name, records in snvq_variants(rng).items():
        built[name] = write(out, name, records)
    main = umi_main(rng)
    built["umi_main"] = write(out, "umi_main", main)
    shuffled = list(main)
    rng.shuffle(shuffled)
    built["umi_unsorted"] = write(out, "umi_unsorted", shuffled, order="unsorted")
    built["umi_queryname"] = write(out, "umi_queryname", main, order="queryname")
    for name, records in umi_variants(rng).items():
        built[name] = write(out, name, records)
    # A generator of its own, so the inputs above keep their records.
    for name, records in repeated_tag_variants(random.Random(seed + 1)).items():
        built[name] = write(out, name, records)
    return built


if __name__ == "__main__":
    import sys
    from pathlib import Path
    for name, path in build_all(Path(sys.argv[1])).items():
        print(name, path)

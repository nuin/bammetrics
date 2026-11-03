"""The VCF corpus for CollectVariantCallingMetrics, AccumulateVariantCallingMetrics
and GenotypeConcordance (plan tier T6, synthetic).

build_all(out) writes into `out` and returns {name: path}:

  truth.vcf / truth.vcf.gz (+ .tbi)   one sample T: one record per (truth state, call
                                       state) cell of the GenotypeConcordance scheme that
                                       the states can reach, SNPs on c1 and indels on c2,
                                       plus the odd sites (MNP, symbolic, SNP against an
                                       indel, two records at one position, a c3-only site,
                                       deletions and insertions whose reference allele is
                                       longer on one side)
  call.vcf / call.vcf.gz (+ .tbi)     samples C (the call side of each cell) and C2 (hom-ref
                                       everywhere, so the sample choice changes the numbers)
  calls3.vcf / calls3.vcf.gz (+.tbi)  samples S1, S2, S3 for CollectVariantCallingMetrics:
                                       transitions and transversions in and out of dbSNP,
                                       multiallelic SNPs, insertions, deletions, complex
                                       indels, MNPs, a symbolic and a mixed record, filtered
                                       sites, hom-ref-only and no-call-only sites, GQ 0,
                                       AD present, absent and unparsable, singletons
  dbsnp.vcf / dbsnp.vcf.gz (+ .tbi)   sites only: SNPs, an MNP, a multiallelic SNP, an
                                       insertion, a five-base deletion, a filtered site
  gvcf.vcf.gz (+ .tbi)                one sample G1 with <NON_REF> alleles and END blocks
  sites_only.vcf.gz (+ .tbi)          calls3 without genotype columns
  nocontigs.vcf                       truth.vcf without ##contig lines
  calls3_nocontigs.vcf.gz (+ .tbi)    calls3 without ##contig lines
  mismatch_truth.vcf, mismatch_call.vcf  one site whose reference alleles share no prefix
  haploid_call.vcf                    one site whose call genotype has one allele
  targets.interval_list               c1:1-2500 and c2:1000-2000
  wide.interval_list                  c1:2000-5000, c2:1-1500 and c3:1-1000
  left.interval_list, right.interval_list  c1 and c2 whole, the Accumulate shards
  other.interval_list                 the intervals of another dictionary
  ref.dict                            the dictionary of the three contigs
"""
import random
from itertools import product
from pathlib import Path

import pysam

CONTIGS = [("c1", 5000), ("c2", 3000), ("c3", 1000)]

TRUTH_STATES = ["MISSING", "HOM_REF", "HET_REF_VAR1", "HET_VAR1_VAR2", "HOM_VAR1", "NO_CALL",
                "LOW_GQ", "LOW_DP", "VC_FILTERED", "GT_FILTERED", "IS_MIXED"]
CALL_STATES = ["MISSING", "HOM_REF", "HET_REF_VAR1", "HET_VAR1_VAR2", "HOM_VAR1", "HET_REF_VAR2",
               "HET_REF_VAR3", "HET_VAR1_VAR3", "HET_VAR3_VAR4", "HOM_VAR2", "HOM_VAR3", "NO_CALL",
               "LOW_GQ", "LOW_DP", "VC_FILTERED", "GT_FILTERED", "IS_MIXED"]

# Cells the GA4GH scheme marks NA: the call's extra alleles need the truth alleles.
NA_CELLS = (
    {(t, "HET_REF_VAR2") for t in TRUTH_STATES if t not in ("HET_REF_VAR1", "HOM_VAR1")}
    | {(t, "HET_REF_VAR3") for t in TRUTH_STATES if t != "HET_VAR1_VAR2"}
    | {(t, "HET_VAR1_VAR3") for t in TRUTH_STATES if t != "HET_VAR1_VAR2"}
    | {(t, "HET_VAR3_VAR4") for t in TRUTH_STATES[5:]}
    | {(t, "HOM_VAR2") for t in TRUTH_STATES if t not in ("HET_REF_VAR1", "HET_VAR1_VAR2", "HOM_VAR1")}
    | {(t, "HOM_VAR3") for t in TRUTH_STATES if t != "HET_VAR1_VAR2"}
    | {("MISSING", "MISSING")}
)

# Cells no pair of records can produce, with the reason (KNOWN_DIFFERENCES.md lists them).
UNREACHABLE = {
    ("MISSING", "HOM_REF"): "both sub-contexts are NO_VARIATION, so the pair is not classified",
    ("MISSING", "IS_MIXED"): "a NO_VARIATION truth with a MIXED call is not classified",
    ("HOM_REF", "MISSING"): "both sub-contexts are NO_VARIATION",
    ("HOM_REF", "HOM_REF"): "both sub-contexts are NO_VARIATION",
    ("HOM_REF", "IS_MIXED"): "a NO_VARIATION truth with a MIXED call is not classified",
    ("IS_MIXED", "MISSING"): "a MIXED truth with a NO_VARIATION call is not classified",
    ("IS_MIXED", "HOM_REF"): "a MIXED truth with a NO_VARIATION call is not classified",
    ("IS_MIXED", "IS_MIXED"): "two MIXED records are not classified",
    ("MISSING", "HET_VAR3_VAR4"): "without truth alleles the call's alleles take positions 1 and 2",
    ("HOM_REF", "HET_VAR3_VAR4"): "without truth alleles the call's alleles take positions 1 and 2",
    ("HET_VAR1_VAR2", "HOM_VAR2"): "a call allele at position 2 makes Picard re-add the truth alleles swapped, so the call reads HOM_VAR1",
}

# The cell the hom-ref scheme fills from the interval bases no record covered.
MISSING_HOM_ONLY = {("MISSING", "MISSING")}

REACHABLE = {(t, c) for t in TRUTH_STATES for c in CALL_STATES} - NA_CELLS - set(UNREACHABLE)

HEADER_LINES = [
    "##fileformat=VCFv4.2",
    '##FILTER=<ID=LowQual,Description="Low quality">',
    '##FILTER=<ID=q10,Description="Quality below 10">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">',
]


def header(samples, contigs=True):
    lines = list(HEADER_LINES)
    if contigs:
        lines[1:1] = [f"##contig=<ID={name},length={length}>" for name, length in CONTIGS]
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples]
    if not samples:
        columns = columns[:8]
    return lines + ["\t".join(columns)]


def record(contig, pos, ref, alts, genotypes, filt="PASS", info=".", keys=("GT", "GQ", "DP", "AD", "FT"),
           ident="."):
    """One VCF line. `genotypes` is a list of dicts with GT, GQ, DP, AD, FT (missing keys
    print `.`); with no genotypes the line has no FORMAT column."""
    fields = [contig, str(pos), ident, ref, ",".join(alts) if alts else ".", "50", filt, info]
    if genotypes:
        fields.append(":".join(keys))
        for g in genotypes:
            values = [str(g[k]) if k in g else "." for k in keys]
            # Trailing missing values are dropped, as callers write them.
            while len(values) > 1 and values[-1] == ".":
                values.pop()
            fields.append(":".join(values))
    return "\t".join(fields)


def gt(call, gq=40, dp=20, ad=None, ft=None):
    g = {"GT": call, "GQ": gq, "DP": dp}
    if ad is not None:
        g["AD"] = ad
    if ft is not None:
        g["FT"] = ft
    return g


# --- GenotypeConcordance: one record per reachable scheme cell -------------------------


def truth_recipe(state, alleles):
    """(alts, genotype dict, FILTER) of the truth sample for a state. `alleles` gives the
    variant alleles v1 and v2 and the mixed-type allele."""
    v1, v2, mixed = alleles
    return {
        "HOM_REF": ([v1], gt("0/0"), "PASS"),
        "HET_REF_VAR1": ([v1], gt("0/1"), "."),
        "HET_VAR1_VAR2": ([v1, v2], gt("1/2"), "PASS"),
        "HOM_VAR1": ([v1], gt("1/1"), "PASS"),
        "NO_CALL": ([v1], gt("./1"), "PASS"),
        "LOW_GQ": ([v1], gt("0/1", gq=5, dp=3), "PASS"),
        "LOW_DP": ([v1], gt("0/1", gq=50, dp=3), "PASS"),
        "VC_FILTERED": ([v1], gt("0/1"), "LowQual"),
        "GT_FILTERED": ([v1], gt("0/1", ft="lowq"), "PASS"),
        "IS_MIXED": ([v1, mixed], gt("1/2"), "PASS"),
    }[state]


def call_recipe(state, truth_state, alleles):
    """(alts, genotype dict, FILTER, FORMAT keys) of the call sample for a state against the
    truth state. v3 and v4 are alleles the truth never carries."""
    v1, v2, v3, v4, mixed = alleles
    plain = ("GT", "GQ", "DP", "AD", "FT")
    if state == "HET_VAR3_VAR4":
        alts = [v3, v4] if truth_state == "HET_VAR1_VAR2" else [v2, v3]
        return alts, gt("1/2"), "PASS", plain
    if state == "NO_CALL":
        # A half call keeps the record's variant type; a full no-call drops it.
        call = "./." if truth_state in ("HET_REF_VAR1", "HET_VAR1_VAR2", "HOM_VAR1") else "./1"
        return [v1], gt(call), "PASS", plain
    return {
        "HOM_REF": ([v1], gt("0/0"), "PASS", plain),
        "HET_REF_VAR1": ([v1], gt("0/1"), "PASS", plain),
        "HET_VAR1_VAR2": ([v1, v2], gt("1/2"), "PASS", plain),
        # GQ and DP absent: neither threshold can fire.
        "HOM_VAR1": ([v1], {"GT": "1/1"}, "PASS", ("GT",)),
        "HET_REF_VAR2": ([v2], gt("0/1"), "PASS", plain),
        "HET_REF_VAR3": ([v3], gt("0/1"), "PASS", plain),
        "HET_VAR1_VAR3": ([v1, v3], gt("1/2"), "PASS", plain),
        "HOM_VAR2": ([v2], gt("1/1"), "PASS", plain),
        "HOM_VAR3": ([v3], gt("1/1"), "PASS", plain),
        "LOW_GQ": ([v1], gt("0/1", gq=5, dp=3), "PASS", plain),
        "LOW_DP": ([v1], gt("0/1", gq=50, dp=3), "PASS", plain),
        "VC_FILTERED": ([v1], gt("0/1"), "LowQual", plain),
        "GT_FILTERED": ([v1], gt("0/1", ft="lowq"), "PASS", plain),
        "IS_MIXED": ([v1, mixed], gt("1/2"), "PASS", plain),
    }[state]


SNP_ALLELES = ("A", ("G", "T", "AT"), ("G", "T", "C", "N", "AT"))
INDEL_ALLELES = ("A", ("AT", "ATT", "G"), ("AT", "ATT", "ATTT", "ATTTT", "G"))


def cell_records(contig, first_pos, ref, truth_alleles, call_alleles):
    """Truth and call lines for every reachable cell, one position each, ten bases apart."""
    truth, call = [], []
    pos = first_pos
    c2 = gt("0/0", gq=99, dp=30)
    for truth_state, call_state in product(TRUTH_STATES, CALL_STATES):
        if (truth_state, call_state) not in REACHABLE:
            continue
        if truth_state != "MISSING":
            alts, g, filt = truth_recipe(truth_state, truth_alleles)
            truth.append(record(contig, pos, ref, alts, [g], filt))
        if call_state != "MISSING":
            alts, g, filt, keys = call_recipe(call_state, truth_state, call_alleles)
            c2_here = {"GT": "0/0"} if keys == ("GT",) else c2
            call.append(record(contig, pos, ref, alts, [g, c2_here], filt, keys=keys))
        pos += 10
    return truth, call, pos


def odd_records():
    """Sites outside the cell product: unclassified pairs, normalization, a second contig."""
    c2 = gt("0/0", gq=99, dp=30)
    truth = [
        # A SNP against an indel at the same position: not classified.
        record("c1", 20, "A", ["G"], [gt("0/1")]),
        # An MNP: not classified.
        record("c1", 30, "AC", ["GT"], [gt("0/1")]),
        # A symbolic deletion with END: not classified.
        record("c1", 40, "A", ["<DEL>"], [gt("0/1")], info="END=45"),
        # Two truth records at one position pair with the call's one and nothing.
        record("c1", 50, "A", ["G"], [gt("0/1")]),
        record("c1", 50, "A", ["AT"], [gt("0/1")]),
        # A deletion whose reference is longer on the call side.
        record("c2", 20, "AT", ["A"], [gt("0/1")]),
        # An insertion whose reference is longer on the call side: the splice
        # gives ACT against ATC.
        record("c2", 30, "A", ["AT"], [gt("1/1")]),
        # A deletion whose reference is longer on the truth side.
        record("c2", 40, "ATCG", ["AG"], [gt("0/1")]),
        # Truth only, on c3 (outside targets.interval_list).
        record("c3", 100, "C", ["T"], [gt("0/1")]),
        # A hom-ref truth against a filtered call: an EMPTY cell.
        record("c1", 60, "A", ["G"], [gt("0/0")]),
    ]
    call = [
        record("c1", 20, "A", ["AT"], [gt("0/1"), c2]),
        record("c1", 30, "AC", ["GT"], [gt("0/1"), c2]),
        record("c1", 40, "A", ["G"], [gt("0/1"), c2]),
        record("c1", 50, "A", ["G"], [gt("0/1"), c2]),
        record("c2", 20, "ATC", ["AC"], [gt("0/1"), c2]),
        record("c2", 30, "AC", ["ATC"], [gt("1/1"), c2]),
        record("c2", 40, "AT", ["A"], [gt("0/1"), c2]),
        record("c1", 60, "A", ["G"], [gt("0/1"), c2], "q10"),
        # Call only, before the first cell site, with a GQ of 0.
        record("c1", 70, "A", ["C"], [gt("0/1", gq=0, dp=1), c2]),
    ]
    return truth, call


def concordance_files(out):
    truth, call = odd_records()
    t, c, _ = cell_records("c1", 100, *SNP_ALLELES)
    truth += t
    call += c
    t, c, _ = cell_records("c2", 100, *INDEL_ALLELES)
    truth += t
    call += c
    truth.sort(key=sort_key)
    call.sort(key=sort_key)
    write_vcf(out / "truth.vcf", header(["T"]), truth, index=True)
    write_vcf(out / "call.vcf", header(["C", "C2"]), call, index=True)
    write_vcf(out / "nocontigs.vcf", header(["T"], contigs=False), truth)
    write_vcf(out / "mismatch_truth.vcf", header(["T"]), [record("c1", 10, "AT", ["A"], [gt("0/1")])])
    write_vcf(out / "mismatch_call.vcf", header(["C"]), [record("c1", 10, "AG", ["A"], [gt("0/1")])])
    write_vcf(out / "haploid_call.vcf", header(["C"]),
              [record("c1", 10, "A", ["G"], [{"GT": "1", "GQ": 30}], keys=("GT", "GQ"))])


def sort_key(line):
    fields = line.split("\t")
    return ([c for c, _ in CONTIGS].index(fields[0]), int(fields[1]))


# --- CollectVariantCallingMetrics ----------------------------------------------------


def calling_records(rng):
    """calls3.vcf: three samples over c1 and c2 (none on c3)."""
    lines = []

    def hom_ref():
        return gt("0/0", gq=rng.choice([30, 60, 99]), dp=rng.randint(10, 40))

    def het(ad=True, gq=None):
        r, a = rng.randint(2, 30), rng.randint(2, 30)
        g = gt("0/1", gq=rng.choice([0, 20, 45, 99]) if gq is None else gq, dp=r + a)
        if ad:
            g["AD"] = f"{r},{a}"
        return g

    def hom_var(gq=None):
        a = rng.randint(5, 40)
        g = gt("1/1", gq=rng.choice([0, 30, 99]) if gq is None else gq, dp=a, ad=f"0,{a}")
        return g

    def no_call():
        return {"GT": "./.", "GQ": 0}

    # Biallelic SNPs: transitions and transversions, varied genotype mixes.
    snp_alts = {"A": ["G", "C", "T"], "C": ["T", "A", "G"], "G": ["A", "C", "T"], "T": ["C", "A", "G"]}
    pos = 100
    for i in range(60):
        ref = rng.choice("ACGT")
        alt = snp_alts[ref][0] if i % 3 else snp_alts[ref][1 + i % 2]
        mix = i % 7
        genotypes = {
            0: [het(), hom_ref(), hom_ref()],  # singleton
            1: [het(), het(), hom_var()],
            2: [hom_var(), hom_ref(), hom_ref()],
            3: [het(ad=False), het(), no_call()],
            4: [no_call(), no_call(), het()],  # singleton with no-calls
            5: [hom_ref(), het(gq=0), hom_ref()],  # singleton, GQ 0
            6: [hom_var(), hom_var(), het()],
        }[mix]
        filt = "LowQual" if i % 11 == 10 else "PASS"
        lines.append(record("c1", pos, ref, [alt], genotypes, filt))
        pos += 7
    # Multiallelic SNPs.
    for i in range(8):
        genotypes = [gt("1/2", gq=50, dp=20, ad="0,9,11"), het(), hom_ref()] if i % 2 else \
            [hom_ref(), gt("0/2", gq=40, dp=18, ad="9,0,9"), hom_ref()]
        lines.append(record("c1", 600 + 9 * i, "A", ["G", "T"], genotypes, "LowQual" if i == 7 else "."))
    # Insertions and deletions, some filtered.
    pos = 1000
    for i in range(30):
        if i % 2:
            ref, alts = "A", ["A" + "T" * (1 + i % 3)]
        else:
            ref, alts = "A" + "C" * (1 + i % 4), ["A"]
        genotypes = {0: [het(), hom_ref(), hom_ref()], 1: [hom_var(), het(), hom_ref()],
                     2: [no_call(), hom_var(), hom_var()], 3: [het(), het(), het()],
                     4: [hom_ref(), hom_ref(), het(gq=0)]}[i % 5]
        filt = "q10" if i % 9 == 8 else "PASS"
        lines.append(record("c1", pos, ref, alts, genotypes, filt))
        pos += 11
    # Complex indels: unequal lengths with different first bases, and multiallelic indels.
    for i, (ref, alts) in enumerate([("AT", ["GCC"]), ("A", ["AT", "ATT"]), ("ACG", ["T"]), ("ACT", ["A", "AC"]),
                                     ("AT", ["GCC"]), ("A", ["AT", "ATT"])]):
        genotypes = [hom_ref(), gt("1/2", gq=33, dp=14), hom_ref()] if len(alts) > 1 else [het(), hom_ref(), hom_var()]
        lines.append(record("c1", 1500 + 13 * i, ref, alts, genotypes, "LowQual" if i == 5 else "PASS"))
    # MNPs, a symbolic record, mixed records: no bucket, but singletons count.
    lines.append(record("c1", 1700, "AC", ["GT"], [het(), hom_ref(), hom_ref()]))
    lines.append(record("c1", 1710, "A", ["<DEL>"], [het(ad=False), hom_ref(), hom_ref()], info="END=1720"))
    lines.append(record("c1", 1730, "A", ["G", "AT"], [gt("1/2", gq=20, dp=10), hom_ref(), het()]))
    lines.append(record("c1", 1740, "A", ["G", "AT"], [hom_ref(), hom_ref(), gt("0/1", gq=0, dp=10)]))
    # Every genotype hom-ref, every genotype no-call, and an ALT of `.`.
    lines.append(record("c1", 1800, "A", ["G"], [hom_ref(), hom_ref(), hom_ref()]))
    lines.append(record("c1", 1810, "A", ["G"], [no_call(), no_call(), no_call()]))
    lines.append(record("c1", 1820, "A", [], [hom_ref(), hom_ref(), hom_ref()]))
    # An AD that does not parse, an AD missing on one het, a GQ with a fraction.
    lines.append(record("c1", 1830, "A", ["G"], [gt("0/1", gq=20, dp=9, ad="5,."), het(), hom_ref()]))
    lines.append(record("c1", 1840, "C", ["T"], [{"GT": "0/1", "GQ": "20.5", "DP": 9}, hom_ref(), hom_ref()]))
    # c2: a second contig with SNPs and indels, in and out of dbSNP.
    pos = 100
    for i in range(40):
        if i % 4 == 3:
            ref, alts = ("AC", ["A"]) if i % 8 == 3 else ("C", ["CGA"])
        else:
            ref = rng.choice("ACGT")
            alts = [snp_alts[ref][i % 3]]
        genotypes = [het(), hom_var(), hom_ref()] if i % 2 else [hom_ref(), het(), het()]
        lines.append(record("c2", pos, ref, alts, genotypes, "LowQual" if i % 13 == 12 else "PASS"))
        pos += 17
    # c2 also carries multiallelic SNPs, complex indels, filtered indels and an
    # insertion-heavy block, so the c1 and c2 shards differ in every column.
    for i in range(3):
        genotypes = [gt("1/2", gq=50, dp=20, ad="0,9,11"), hom_ref(), het()] if i else [hom_ref(), hom_ref(), gt("0/2", gq=40, dp=18, ad="9,0,9")]
        lines.append(record("c2", 900 + 9 * i, "G", ["A", "C"], genotypes))
    for i, (ref, alts) in enumerate([("GT", ["ACC"]), ("G", ["GA", "GAA"])]):
        lines.append(record("c2", 950 + 13 * i, ref, alts, [het(), hom_var(), hom_ref()]))
    for i in range(6):
        ref, alts = ("T", ["TG"]) if i < 5 else ("TGA", ["T"])
        lines.append(record("c2", 1200 + 11 * i, ref, alts, [het(), hom_ref(), het()], "q10" if i % 3 == 2 else "PASS"))
    lines.sort(key=sort_key)
    return lines


def dbsnp_records(calls):
    """dbSNP sites: half of the calls3 SNPs and indels by position, some with a different
    type at the same position, plus a five-base deletion, an MNP and a filtered site."""
    lines = []
    for i, line in enumerate(calls):
        contig, pos, _, ref, alt = line.split("\t")[:5]
        if i % 2:
            continue
        alts = alt.split(",")
        if alt == "." or alts[0].startswith("<"):
            continue
        if i % 6 == 0 and len(ref) == 1 and len(alts[0]) == 1:
            # A dbSNP indel where the call is a SNP: the SNP bit stays clear.
            lines.append(record(contig, int(pos), ref, [ref + "A"], [], ident=f"rs{i}"))
        elif i % 6 == 2 and (len(ref) > 1 or len(alts[0]) > 1):
            # A dbSNP SNP where the call is an indel.
            lines.append(record(contig, int(pos), ref[0], ["G" if ref[0] != "G" else "A"], [], ident=f"rs{i}"))
        else:
            lines.append(record(contig, int(pos), ref, alts[:1], [], ident=f"rs{i}"))
    lines.append(record("c1", 998, "ACCCCC", ["A"], [], ident="rs_del"))  # spans 998 to 1003
    lines.append(record("c1", 1697, "GTA", ["CAT"], [], ident="rs_mnp"))
    lines.append(record("c2", 117, "A", ["G"], [], "LowQual", ident="rs_filtered"))
    lines.append(record("c1", 1831, "A", ["T", "C"], [], ident="rs_multi"))
    lines.sort(key=sort_key)
    return lines


def gvcf_records():
    def g(call, gq, dp, ad):
        return {"GT": call, "GQ": gq, "DP": dp, "AD": ad}

    lines = [
        record("c1", 1, "A", ["<NON_REF>"], [g("0/0", 50, 20, "20,0")], info="END=99"),
        record("c1", 100, "A", ["G", "<NON_REF>"], [g("0/1", 40, 20, "9,11,0")]),
        record("c1", 101, "C", ["<NON_REF>"], [g("0/0", 50, 20, "20,0")], info="END=199"),
        record("c1", 200, "C", ["T", "<NON_REF>"], [g("1/1", 60, 22, "0,22,0")]),
        record("c1", 210, "G", ["A", "T", "<NON_REF>"], [g("1/2", 30, 20, "0,9,11,0")]),
        record("c1", 220, "T", ["TA", "<NON_REF>"], [g("0/1", 35, 18, "8,10,0")]),
        record("c1", 230, "TAC", ["T", "<NON_REF>"], [g("1/1", 45, 15, "0,15,0")]),
        record("c1", 240, "A", ["G", "<NON_REF>"], [g("0/0", 10, 5, "5,0,0")]),
        record("c1", 250, "A", ["C", "<NON_REF>"], [g("0/1", 0, 3, "1,2,0")], "LowQual"),
        record("c1", 260, "A", ["AT", "<NON_REF>"], [g("0/1", 0, 4, "2,2,0")], "LowQual"),
        record("c1", 270, "A", ["<NON_REF>"], [g("./.", 0, 0, "0,0")], info="END=300"),
        record("c1", 1500, "AT", ["GCC", "<NON_REF>"], [g("0/1", 25, 12, "6,6,0")]),
        record("c2", 100, "A", ["G", "<NON_REF>"], [g("0/1", 40, 20, "9,11,0")]),
        record("c2", 117, "T", ["C", "<NON_REF>"], [g("1/1", 40, 20, "0,20,0")]),
        record("c2", 151, "C", ["CGA", "<NON_REF>"], [g("0/1", 40, 20, "9,11,0")]),
    ]
    return lines


# --- files -------------------------------------------------------------------------------


def write_vcf(path, header_lines, records, index=False):
    text = "\n".join(header_lines + records) + "\n"
    path.write_text(text)
    if index:
        gz = str(path) + ".gz"
        pysam.tabix_compress(str(path), gz, force=True)
        pysam.tabix_index(gz, preset="vcf", force=True)
        return Path(gz)
    return path


def interval_list(path, intervals, contigs=CONTIGS):
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in contigs]
    lines += [f"{c}\t{s}\t{e}\t+\t{n}" for c, s, e, n in intervals]
    path.write_text("\n".join(lines) + "\n")
    return path


def build_all(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(7)
    concordance_files(out)
    calls = calling_records(rng)
    write_vcf(out / "calls3.vcf", header(["S1", "S2", "S3"]), calls, index=True)
    write_vcf(out / "sites_only.vcf", header([]), ["\t".join(line.split("\t")[:8]) for line in calls],
              index=True)
    write_vcf(out / "dbsnp.vcf", header([]), dbsnp_records(calls), index=True)
    write_vcf(out / "gvcf.vcf", header(["G1"]), gvcf_records(), index=True)
    write_vcf(out / "calls3_nocontigs.vcf", header(["S1", "S2", "S3"], contigs=False), calls, index=True)
    interval_list(out / "targets.interval_list", [("c1", 1, 2500, "left_half"), ("c2", 1000, 2000, "middle")])
    interval_list(out / "wide.interval_list", [("c1", 2000, 5000, "right"), ("c2", 1, 1500, "front"),
                                               ("c3", 1, 1000, "all")])
    interval_list(out / "left.interval_list", [("c1", 1, 5000, "c1")])
    interval_list(out / "right.interval_list", [("c2", 1, 3000, "c2")])
    interval_list(out / "other.interval_list", [("chr1", 1, 100, "x")], contigs=[("chr1", 1000)])
    (out / "ref.dict").write_text("@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{n}\tLN:{l}\n" for n, l in CONTIGS))
    names = ["truth.vcf", "truth.vcf.gz", "call.vcf", "call.vcf.gz", "calls3.vcf", "calls3.vcf.gz",
             "sites_only.vcf.gz", "dbsnp.vcf", "dbsnp.vcf.gz", "gvcf.vcf.gz", "nocontigs.vcf",
             "calls3_nocontigs.vcf.gz", "mismatch_truth.vcf", "mismatch_call.vcf", "haploid_call.vcf",
             "targets.interval_list", "wide.interval_list", "left.interval_list", "right.interval_list",
             "other.interval_list", "ref.dict"]
    return {name: out / name for name in names}


if __name__ == "__main__":
    import sys

    for name, path in build_all(sys.argv[1]).items():
        print(name, path)

"""Synthetic arrays VCFs as GtcToVcf writes them, for CollectArraysVariantCallingMetrics.

build(out, name, seed, ...) writes `out/name.vcf`, `out/name.vcf.gz` and its `.tbi`; the
plain file gets the Tribble `.idx` Picard needs from index_plain(). dbsnp() writes the
sites-only dbSNP VCF of a record list. Every name below is a header key, FORMAT or INFO
field, filter or control that Picard 3.5.0 reads (picard.arrays.illumina.InfiniumVcfFields,
ArraysControlInfo).

The records of one sample (or several) on the autosomes 1 and 2, X, Y and MT: every
genotype form (hom-ref, het, hom-var, phased, multiallelic het, mixed, no-call, haploid),
a GTA that agrees, says ./. (a zCall rescue), is . or is missing, SNPs, multiallelic SNPs,
insertions, deletions, MNPs and sites with no alternate, and the filters PASS, ., DUPE,
ZEROED_OUT_ASSAY, TRIALLELIC, FAIL_REF and their pairs.
"""
import random
import subprocess
from pathlib import Path

import pysam

from run_picard import picard_env

CONTIGS = [("1", 40000), ("2", 30000), ("X", 30000), ("Y", 12000), ("MT", 16569)]

CONTROLS = ["DNP(High)", "DNP(Bgnd)", "Biotin(High)", "Biotin(Bgnd)", "Extension(A)",
            "Extension(T)", "Extension(C)", "Extension(G)", "TargetRemoval", "Hyb(High)",
            "Hyb(Medium)", "Hyb(Low)", "String(PM)", "String(MM)", "NSB(Bgnd)Red",
            "NSB(Bgnd)Purple", "NSB(Bgnd)Blue", "NSB(Bgnd)Green", "NP(A)", "NP(T)", "NP(C)",
            "NP(G)", "Restore"]
CATEGORIES = {"DNP": "Staining", "Biotin": "Staining", "Extension": "Extension",
              "TargetRemoval": "TargetRemoval", "Hyb": "Hybridization", "String": "Stringency",
              "NSB": "Non-SpecificBinding", "NP": "Non-Polymorphic", "Restore": "Restoration"}

REQUIRED = ["sampleAlias", "arrayType", "autocallGender", "autocallVersion", "autocallDate",
            "imagingDate", "extendedIlluminaManifestVersion", "clusterFile", "p95Green", "p95Red",
            "scannerName"]
OPTIONAL = ["pipelineVersion", "analysisVersionNumber", "expectedGender", "fingerprintGender",
            "gtcCallRate", "zcallVersion", "zcallThresholds"]

FORMAT = ["GT", "IGC", "X", "Y", "NORMX", "NORMY", "R", "THETA", "BAF", "LRR", "GTA", "GTZ"]
FORMAT_LINES = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=IGC,Number=1,Type=Float,Description="Illumina GenCall Confidence Score">',
    '##FORMAT=<ID=X,Number=1,Type=Integer,Description="Raw X intensity">',
    '##FORMAT=<ID=Y,Number=1,Type=Integer,Description="Raw Y intensity">',
    '##FORMAT=<ID=NORMX,Number=1,Type=Float,Description="Normalized X intensity">',
    '##FORMAT=<ID=NORMY,Number=1,Type=Float,Description="Normalized Y intensity">',
    '##FORMAT=<ID=R,Number=1,Type=Float,Description="Normalized R value">',
    '##FORMAT=<ID=THETA,Number=1,Type=Float,Description="Normalized Theta value">',
    '##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B Allele Frequency">',
    '##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R Ratio">',
    '##FORMAT=<ID=GTA,Number=1,Type=String,Description="Illumina Autocall Genotype">',
    '##FORMAT=<ID=GTZ,Number=1,Type=String,Description="Illumina zCall Genotype">',
]
INFO_LINES = [
    '##INFO=<ID=ALLELE_A,Number=1,Type=String,Description="A allele">',
    '##INFO=<ID=ALLELE_B,Number=1,Type=String,Description="B allele">',
    '##INFO=<ID=BEADSET_ID,Number=1,Type=Integer,Description="Bead set ID">',
    '##INFO=<ID=ILLUMINA_STRAND,Number=1,Type=String,Description="Probe strand">',
    '##INFO=<ID=SOURCE,Number=1,Type=String,Description="Probe source">',
]
FILTER_LINES = [
    '##FILTER=<ID=DUPE,Description="Duplicate assays at position">',
    '##FILTER=<ID=TRIALLELIC,Description="Tri-allelic assay">',
    '##FILTER=<ID=FAIL_REF,Description="Assay failed to map to reference">',
    '##FILTER=<ID=ZEROED_OUT_ASSAY,Description="Assay Zeroed out (marked as uncallable)">',
]
FILTERS = ["PASS"] * 14 + [".", ".", "DUPE", "DUPE", "ZEROED_OUT_ASSAY", "ZEROED_OUT_ASSAY",
                           "TRIALLELIC", "FAIL_REF", "DUPE;TRIALLELIC",
                           "ZEROED_OUT_ASSAY;TRIALLELIC", "TRIALLELIC;FAIL_REF"]


def default_values(seed):
    """The header values of one chip well, a little different per seed."""
    rng = random.Random(seed)
    return [
        ("arrayType", rng.choice(["GSA-24v3-0_A1", "PsychChip_v1-1_15073391_A1", "MEG_AllofUs_20002558X351448_A2"])),
        ("extendedIlluminaManifestVersion", rng.choice(["1.3", "1.5", "2.0"])),
        ("chipWellBarcode", f"2041262900{seed:02d}_R0{1 + seed % 8}C01"),
        ("analysisVersionNumber", str(1 + seed % 5)),
        ("sampleAlias", f"NA{12878 + seed}"),
        ("expectedGender", rng.choice(["Female", "Male", "Unknown", "NotReported"])),
        ("fingerprintGender", rng.choice(["F", "M", "U", "N"])),
        ("autocallGender", rng.choice(["Female", "Male", "Unknown", "F", "M", "U"])),
        ("autocallDate", f"0{1 + seed % 9}/1{seed % 10}/20{10 + seed % 15} {seed % 24}:{seed % 60:02d}"),
        ("imagingDate", f"{1 + seed % 12:02d}/{1 + seed % 28:02d}/20{10 + seed % 15} "
                        f"{1 + seed % 12:02d}:{seed % 60:02d}:{(7 * seed) % 60:02d} {'AM' if seed % 2 else 'PM'}"),
        ("clusterFile", f"cluster_{seed}.egt"),
        ("manifestFile", f"manifest_{seed}.bpm"),
        ("extendedManifestFile", f"manifest_{seed}.extended.csv"),
        ("autocallVersion", rng.choice(["3.0.0", "2.0.1", "1.6.2.2"])),
        ("zcallVersion", rng.choice(["1.0.0.0", "1.1", "1.0.0.3"])),
        ("zcallThresholds", f"thresholds_{seed}.7.txt"),
        ("p95Red", str(9000 + 37 * seed)),
        ("p95Green", str(7000 + 29 * seed)),
        ("scannerName", rng.choice(["N370", "N580", "N1043"])),
        ("pipelineVersion", f"1.{seed % 7}.{seed % 3}"),
        ("gtcCallRate", f"0.9{seed % 10}{(3 * seed) % 10}"),
    ]


def control_lines(seed):
    rng = random.Random(seed + 100)
    lines = []
    for name in CONTROLS:
        category = next(v for k, v in CATEGORIES.items() if name.startswith(k))
        lines.append(f"##{name}={name}|{category}|{rng.randint(100, 60000)}|{rng.randint(100, 60000)}")
    return lines


def header_lines(values, samples, contigs=CONTIGS, controls=None):
    """`values` is a list of (key, value) other lines in order; `controls` replaces the
    control lines when given."""
    lines = ["##fileformat=VCFv4.2", *FILTER_LINES, *FORMAT_LINES, *INFO_LINES]
    lines += [f"##{key}={value}" for key, value in values]
    lines += controls if controls is not None else control_lines(0)
    lines += [f"##contig=<ID={name},length={length},assembly=GRCh37>" for name, length in contigs]
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                            *samples]))
    return lines


def site_alleles(rng):
    """(ref, alts) of one site."""
    kind = rng.choices(["snp", "multi", "ins", "del", "mnp", "none"], [60, 6, 8, 8, 3, 3])[0]
    bases = "ACGT"
    ref = rng.choice(bases)
    others = [b for b in bases if b != ref]
    if kind == "snp":
        return ref, [rng.choice(others)]
    if kind == "multi":
        return ref, rng.sample(others, 2)
    if kind == "ins":
        return ref, [ref + rng.choice(bases) * rng.randint(1, 3)]
    if kind == "del":
        return ref + "".join(rng.choice(bases) for _ in range(rng.randint(1, 3))), [ref]
    if kind == "mnp":
        return ref + rng.choice(bases), [rng.choice(others) + rng.choice(bases)]
    return ref, []


def genotype(rng, n_alleles, haploid):
    if haploid:
        return rng.choice(["0", "1", "."] if n_alleles > 1 else ["0", "."])
    calls = ["0/0", "./.", "0/0"]
    if n_alleles > 1:
        calls += ["0/1", "0/1", "1/1", "0|1", "1|0", "./1", "1/1"]
    if n_alleles > 2:
        calls += ["1/2", "0/2", "2/2"]
    return rng.choice(calls)


def sample_value(rng, gt):
    """One sample's FORMAT values: the GTA agrees with the GT, says ./. (a zCall rescue),
    is ., or is missing with GTZ."""
    numbers = [f"{rng.random():.4f}", str(rng.randint(0, 30000)), str(rng.randint(0, 30000)),
               f"{rng.random() * 3:.3f}", f"{rng.random() * 3:.3f}", f"{rng.random() * 4:.3f}",
               f"{rng.random():.4f}", f"{rng.random():.4f}", f"{rng.uniform(-2, 2):.4f}"]
    gta = rng.choices([gt, "./.", "."], [6, 2, 1])[0]
    values = [gt, *numbers, gta, gt]
    if rng.random() < 0.08:
        values = values[:-2]
    return ":".join(values)


def records(seed, samples=1, per_contig=None):
    """Record lines of `samples` samples on every contig, sorted by position: 40 to 68
    sites per contig by the seed."""
    rng = random.Random(seed)
    per_contig = 40 + 7 * (seed % 5) if per_contig is None else per_contig
    lines = []
    for contig, length in CONTIGS:
        positions = sorted(rng.sample(range(10, length - 10), per_contig))
        haploid = contig in ("Y", "MT")
        for i, pos in enumerate(positions):
            ref, alts = site_alleles(rng)
            n_alleles = 1 + len(alts)
            filt = rng.choice(FILTERS)
            info = f"ALLELE_A={ref};ALLELE_B={alts[0] if alts else ref};BEADSET_ID={rng.randint(1, 999)}"
            values = [sample_value(rng, genotype(rng, n_alleles, haploid and rng.random() < 0.8))
                      for _ in range(samples)]
            lines.append("\t".join([contig, str(pos), f"rs{seed}{contig}{i}", ref,
                                    ",".join(alts) or ".", ".", filt, info, ":".join(FORMAT),
                                    *values]))
    return lines


def dbsnp_lines(seed, record_lines, contigs=CONTIGS):
    """A sites-only dbSNP VCF: about half the SNP sites, some indel sites, an indel record at a
    SNP site, a SNP record at an indel site, and sites of its own."""
    rng = random.Random(seed + 1)
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={name},length={length},assembly=GRCh37>" for name, length in contigs]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    sites = []
    for line in record_lines:
        contig, pos, _, ref, alt = line.split("\t")[:5]
        alts = [] if alt == "." else alt.split(",")
        roll = rng.random()
        snp = alts and len(ref) == 1 and all(len(a) == 1 for a in alts)
        if snp and roll < 0.5:
            sites.append((contig, int(pos), ref, alts[0]))
        elif snp and roll < 0.6:
            sites.append((contig, int(pos), ref, ref + "T"))
        elif not snp and roll < 0.4:
            sites.append((contig, int(pos), ref, alts[0] if alts else ref + "T"))
        elif not snp and roll < 0.55:
            sites.append((contig, int(pos), ref[0], "G" if ref[0] != "G" else "C"))
    for contig, length in contigs:
        for _ in range(5):
            sites.append((contig, rng.randint(10, length - 10), "A", "C"))
    for i, (contig, pos, ref, alt) in enumerate(sorted(sites, key=lambda s: (
            [c for c, _ in contigs].index(s[0]), s[1]))):
        lines.append(f"{contig}\t{pos}\trs{900000 + i}\t{ref}\t{alt}\t.\t.\t.")
    return lines


def write(path, lines):
    """`path` (plain) and `path.gz` with its tabix index."""
    path = Path(path)
    path.write_text("\n".join(lines) + "\n")
    gz = str(path) + ".gz"
    pysam.tabix_compress(str(path), gz, force=True)
    pysam.tabix_index(gz, preset="vcf", force=True)
    return Path(gz)


def index_plain(vcf, out):
    """A plain VCF with the Tribble .idx htsjdk requires, written by Picard SortVcf."""
    subprocess.run(["picard", "SortVcf", f"I={vcf}", f"O={out}", "CREATE_INDEX=true",
                    "VERBOSITY=ERROR", "QUIET=true"], check=True, env=picard_env(),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return Path(out)


def build(out, name, seed, values=None, samples=("SAMPLE",), controls=None, contigs=CONTIGS,
          record_lines=None):
    """`out/name.vcf.gz` (with .tbi) and `out/name.vcf`; returns the .gz path."""
    values = default_values(seed) if values is None else values
    if record_lines is None:
        record_lines = records(seed, samples=len(samples))
    controls = control_lines(seed) if controls is None else controls
    lines = header_lines(values, list(samples), contigs=contigs, controls=controls) + record_lines
    return write(Path(out) / f"{name}.vcf", lines)

"""CollectOxoGMetrics against Picard 3.5.0: every CpcgMetrics cell, every
option, CONTEXT_SIZE 0 to 4, CONTEXTS lists, INTERVALS, DB_SNP, and error
parity.

Inputs: the bwa tiers t1 to t3 with the corpus reference, the foundation
corpus, the exome slice, the oxog_bams BAMs (six read groups with per-library
oxidation, an IUPAC and soft-masked reference) and the hand-built BAMs of
asm_bams, read_level_bams and complexity_bams whose records all carry a
header read group. Every case runs Picard with VALIDATION_STRINGENCY=SILENT.
DB_SNP files are written here from each reference: SNPs, deletions,
insertions, MNPs and a filtered record, with and without ##contig lines.

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared. It fails when a comparison above did
not hand over its Picard file, and skips only under -k, --deselect or a node
id.
"""
import random
import subprocess
from collections import defaultdict

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, REFERENCE
from compare import compare_metrics
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from oxog_bams import build_all as build_oxog
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

TOOL = "CollectOxoGMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Every hand-built BAM, each with its reference."""
    bams = {}
    for prefix, build in (("ox", build_oxog), ("asm", build_asm), ("rl", build_read_level),
                          ("elc", build_complexity)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = (bam, bam.parent / "ref.fa")
    return bams


def inputs(name, built, corpus, exome):
    """(bam, reference) of a named input."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name == "exome":
        return exome / "exome.bam", exome / "ref.fa"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam", REFERENCE


def read_fasta(path):
    """{contig: upper-case sequence} in file order."""
    contigs, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                contigs[name] = "".join(parts).upper()
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name:
        contigs[name] = "".join(parts).upper()
    return contigs


def write_dbsnp(reference, out, seed=3, contig_lines=True, extra_contig=None):
    """A sites-only VCF over `reference`: about one record per 60 bases, half
    SNPs, the rest deletions, insertions, MNPs and filtered SNPs."""
    rng = random.Random(seed)
    lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=LowQual,Description="Low quality">']
    contigs = read_fasta(reference)
    if contig_lines:
        lines += [f"##contig=<ID={name},length={len(seq)}>" for name, seq in contigs.items()]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    def base(seq, pos):
        b = seq[pos - 1]
        return b if b in "ACGT" else "N"

    records = list(contigs.items())
    if extra_contig:
        records.append((extra_contig, "ACGT" * 100))
    for name, seq in records:
        pos = rng.randint(1, 30)
        while pos + 6 < len(seq):
            ref = base(seq, pos)
            kind = rng.random()
            filt = "PASS"
            if kind < 0.5:
                alleles = (ref, "A" if ref != "A" else "C")
            elif kind < 0.7:
                alleles = ("".join(base(seq, p) for p in range(pos, pos + 4)), ref)
            elif kind < 0.85:
                alleles = (ref, ref + "TT")
            elif kind < 0.95:
                two = "".join(base(seq, p) for p in range(pos, pos + 2))
                alleles = (two, "".join("A" if b != "A" else "C" for b in two))
            else:
                alleles, filt = (ref, "G" if ref != "G" else "T"), "LowQual"
            lines.append(f"{name}\t{pos}\trs{len(lines)}\t{alleles[0]}\t{alleles[1]}\t50\t{filt}\t.")
            pos += rng.randint(20, 100)
    out.write_text("\n".join(lines) + "\n")
    return out


def compare(bam, ref, args, label):
    cells = compare_metrics(TOOL, bam, [SILENT, f"R={ref}", *args], workdir=bam.parent,
                            expect_file=True, seen=SEEN)
    assert cells > 0, f"{label}: no cell compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))
    return cells


# --- every input at the defaults -------------------------------------------------

# The other hand-built BAMs hold records without an RG tag, or with one the
# header lacks, at counted bases: both tools fail on them (below).
TIERS = ["t1", "t2", "t3"]
BUILT = ["ox_damage", "ox_single_library", "ox_rg_free_low_quality", "asm_unpaired_only",
         "asm_no_pf_group", "asm_all_qcfail", "asm_clips_and_flags", "asm_empty", "rl_q0_cycle",
         "rl_second_only", "rl_star_fields", "rl_oq_short", "rl_empty", "elc_gc_skew",
         "elc_unaligned"]
INPUTS = TIERS + ["foundation", "exome"] + BUILT


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, corpus, exome, name):
    bam, ref = inputs(name, built, corpus, exome)
    compare(bam, ref, [], name)


# Unpaired reads and the insert size filter: with both bounds 0 the filter is gone.
PERMISSIVE = ["MIN_INS=0", "MAX_INS=0", "MQ=0", "Q=0"]
PERMISSIVE_INPUTS = ["asm_unpaired_only", "asm_clips_and_flags", "rl_q0_cycle", "rl_second_only",
                     "elc_gc_skew", "ox_single_library"]


@pytest.mark.parametrize("name", PERMISSIVE_INPUTS)
def test_permissive(built, corpus, exome, name):
    bam, ref = inputs(name, built, corpus, exome)
    compare(bam, ref, PERMISSIVE, f"{name} permissive")


# --- options ---------------------------------------------------------------------------

OPTIONS = {
    "context_size_0": ["CONTEXT_SIZE=0"],
    "context_size_2": ["CONTEXT_SIZE=2"],
    "contexts": ["CONTEXTS=ACA", "CONTEXTS=TCG", "CONTEXTS=GCT", "CONTEXTS=ACA"],
    "contexts_null": ["CONTEXTS=ACA", "CONTEXTS=null", "CONTEXTS=GCG", "CONTEXTS=CCC"],
    "contexts_size_2": ["CONTEXT_SIZE=2", "CONTEXTS=AACAA", "CONTEXTS=GTCGA", "CONTEXTS=TTCTT"],
    "contexts_size_0": ["CONTEXT_SIZE=0", "CONTEXTS=C"],
    "q_0": ["Q=0"],
    "q_30": ["MINIMUM_QUALITY_SCORE=30"],
    "mq_0": ["MQ=0"],
    "mq_60": ["MINIMUM_MAPPING_QUALITY=60"],
    "insert_0_0": ["MIN_INS=0", "MAX_INS=0"],
    "insert_0_200": ["MIN_INS=0", "MAX_INS=200"],
    "insert_100_300": ["MINIMUM_INSERT_SIZE=100", "MAXIMUM_INSERT_SIZE=300"],
    "insert_250_250": ["MIN_INS=250", "MAX_INS=250"],
    "non_pf_false": ["NON_PF=false"],
    "use_oq_false": ["USE_OQ=false"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_0": ["STOP_AFTER=0"],
    "stop_after_negative": ["STOP_AFTER=-5"],
    "stop_after_700": ["STOP_AFTER=700"],
    "all_loose": ["MQ=0", "Q=0", "MIN_INS=0", "MAX_INS=0", "USE_OQ=false", "CONTEXT_SIZE=2"],
}
OPTION_INPUTS = ["foundation", "ox_damage", "t1", "exome"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, corpus, exome, option, name):
    bam, ref = inputs(name, built, corpus, exome)
    compare(bam, ref, OPTIONS[option], f"{name} {option}")


# Contexts the reference holds besides ACGT: N runs, IUPAC codes, and a
# lower-case context that never matches the upper-cased reference.
ODD_CONTEXTS = ["CONTEXTS=NCA", "CONTEXTS=ACN", "CONTEXTS=NCN", "CONTEXTS=RCG", "CONTEXTS=YCS",
                "CONTEXTS=aCg", "CONTEXTS=ACA"]


@pytest.mark.parametrize("extra", [[], ["MIN_INS=0", "MAX_INS=0", "Q=0"]])
def test_iupac_and_lower_case_contexts(built, corpus, exome, extra):
    bam, ref = inputs("ox_damage", built, corpus, exome)
    compare(bam, ref, [*ODD_CONTEXTS, *extra], f"ox_damage odd contexts {extra}")


@pytest.mark.parametrize("name, size", [("foundation", 3), ("ox_damage", 3), ("t1", 4)])
def test_large_context_sizes(built, corpus, exome, name, size):
    bam, ref = inputs(name, built, corpus, exome)
    compare(bam, ref, [f"CONTEXT_SIZE={size}"], f"{name} context_size_{size}")


# --- intervals and dbSNP --------------------------------------------------------------------

@pytest.fixture(scope="module")
def damage_intervals(built, tmp_path_factory):
    """An interval list over the oxog_bams reference: overlapping intervals
    (uniqued), the first and last bases of contigs, and a single base."""
    bam, ref = built["ox_damage"]
    directory = tmp_path_factory.mktemp("oxog_intervals")
    header = (ref.parent / "ref.dict").read_text()
    body = ["c1\t100\t5000\t+\ta", "c1\t4000\t8000\t+\tb", "c2\t1\t300\t+\tc", "c2\t10000\t10000\t-\td",
            "c3\t5800\t6000\t+\te"]
    path = directory / "damage.interval_list"
    path.write_text(header + "\n".join(body) + "\n")
    return path


def interval_cases(corpus, exome, damage_intervals):
    return {
        "foundation": corpus / "target.interval_list",
        "foundation_baits": corpus / "bait.interval_list",
        "t1": CORPUS / "baits.interval_list",
        "t3": CORPUS / "baits.interval_list",
        "exome": exome / "targets.interval_list",
        "ox_damage": damage_intervals,
    }


INTERVAL_NAMES = ["foundation", "foundation_baits", "t1", "t3", "exome", "ox_damage"]


@pytest.mark.parametrize("name", INTERVAL_NAMES)
@pytest.mark.parametrize("extra", [[], ["STOP_AFTER=300", "CONTEXT_SIZE=2"]])
def test_intervals(built, corpus, exome, damage_intervals, name, extra):
    intervals = interval_cases(corpus, exome, damage_intervals)[name]
    bam, ref = inputs(name.replace("_baits", ""), built, corpus, exome)
    compare(bam, ref, [f"INTERVALS={intervals}", *extra], f"{name} intervals {extra}")


DBSNP_NAMES = ["foundation", "t1", "exome", "ox_damage"]


@pytest.mark.parametrize("name", DBSNP_NAMES)
def test_db_snp(built, corpus, exome, tmp_path, name):
    bam, ref = inputs(name, built, corpus, exome)
    vcf = write_dbsnp(ref, tmp_path / "dbsnp.vcf")
    compare(bam, ref, [f"DB_SNP={vcf}"], f"{name} db_snp")


def test_db_snp_without_contig_lines_and_intervals(built, corpus, exome, damage_intervals, tmp_path):
    """No ##contig line: Picard compares no dictionary."""
    bam, ref = inputs("ox_damage", built, corpus, exome)
    vcf = write_dbsnp(ref, tmp_path / "dbsnp.vcf", seed=9, contig_lines=False)
    compare(bam, ref, [f"DB_SNP={vcf}", f"INTERVALS={damage_intervals}", "MIN_INS=0", "MAX_INS=0"],
            "ox_damage db_snp no contigs intervals")


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"O={picard_out}", SILENT, *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


def test_without_reference_fails(built, tmp_path):
    bam, _ = built["ox_damage"]
    both_fail(bam, [], tmp_path, "REFERENCE_SEQUENCE")


@pytest.mark.parametrize("name", ["rl_no_read_groups", "asm_no_read_groups"])
def test_without_read_groups_fails(built, tmp_path, name):
    bam, ref = built[name]
    both_fail(bam, [f"R={ref}"], tmp_path, "requires a read group entry")


@pytest.mark.parametrize("name", ["asm_mix", "asm_some_without_rg", "rl_reads", "elc_illumina"])
def test_counted_base_without_read_group_fails(built, tmp_path, name):
    """Picard's getReadGroup() is null for a record without RG or with an ID
    the header lacks: a NullPointerException at the first such base that
    passes MINIMUM_QUALITY_SCORE."""
    bam, ref = built[name]
    both_fail(bam, [f"R={ref}"], tmp_path, "no read group in the header")


def test_short_oq_fails_once_reached(built, tmp_path):
    """rl_oq_short has an OQ shorter than its read: the defaults never reach a
    base past its end (test_defaults), the permissive filters do, and both
    tools fail; USE_OQ=false reads QUAL and completes."""
    bam, ref = built["rl_oq_short"]
    both_fail(bam, [f"R={ref}", *PERMISSIVE], tmp_path, "OQ shorter than the read")
    compare(bam, ref, [*PERMISSIVE, "USE_OQ=false"], "rl_oq_short permissive no OQ")


def test_low_quality_bases_never_look_up_the_read_group(built, tmp_path):
    """ox_rg_free_low_quality passes at the default (test_defaults) and fails at Q=0."""
    bam, ref = built["ox_rg_free_low_quality"]
    both_fail(bam, [f"R={ref}", "Q=0"], tmp_path, "no read group in the header")
    both_fail(bam, [f"R={ref}", "Q=2"], tmp_path, "no read group in the header")


@pytest.mark.parametrize("args, message", [
    (["CONTEXTS=AC"], "Context AC is not 3 long as implied by CONTEXT_SIZE=1"),
    (["CONTEXTS=AAA"], "Middle base of context sequence AAA must be C"),
    (["CONTEXTS=AcA"], "Middle base of context sequence AcA must be C"),
    (["CONTEXT_SIZE=2", "CONTEXTS=ACA"], "is not 5 long"),
    (["CONTEXT_SIZE=-1"], "CONTEXT_SIZE=-1"),
    (["CONTEXT_SIZE=-1", "CONTEXTS=C"], "is not -1 long"),
    (["MIN_INS=-1"], "MINIMUM_INSERT_SIZE cannot be negative"),
    (["MAX_INS=-5", "MIN_INS=0"], "MAXIMUM_INSERT_SIZE cannot be negative"),
    (["MIN_INS=700"], "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE"),
    (["ASSUME_SORTED=true"], "Unrecognized option: ASSUME_SORTED"),
    (["METRIC_ACCUMULATION_LEVEL=SAMPLE"], "METRIC_ACCUMULATION_LEVEL"),
    (["Q=abc"], "MINIMUM_QUALITY_SCORE=abc"),
])
def test_validation_fails_in_both(built, tmp_path, args, message):
    bam, ref = built["ox_damage"]
    both_fail(bam, [f"R={ref}", *args], tmp_path, message)


@pytest.mark.parametrize("name", ["rl_queryname", "elc_unaligned_queryname"])
def test_queryname_order_fails(built, tmp_path, name):
    bam, ref = built[name]
    both_fail(bam, [f"R={ref}"], tmp_path)


def test_reference_contigs_must_match_the_header(built, corpus, tmp_path):
    bam, _ = built["ox_damage"]
    both_fail(bam, [f"R={corpus / 'ref.fa'}"], tmp_path, "in the BAM header")


def test_db_snp_dictionary_must_match(built, corpus, tmp_path):
    bam, ref = built["ox_damage"]
    vcf = write_dbsnp(corpus / "ref.fa", tmp_path / "other.vcf")
    both_fail(bam, [f"R={ref}", f"DB_SNP={vcf}"], tmp_path, "does not match")


def test_db_snp_contig_missing_from_the_bam_fails(built, tmp_path):
    """Without ##contig lines nothing is compared, and a record on a contig
    the BAM lacks is a NullPointerException in Picard."""
    bam, ref = built["ox_damage"]
    vcf = write_dbsnp(ref, tmp_path / "extra.vcf", contig_lines=False, extra_contig="zz")
    both_fail(bam, [f"R={ref}", f"DB_SNP={vcf}"], tmp_path, "not in the sequence dictionary")


def test_context_size_5_is_refused(built, tmp_path):
    """Picard writes 1,048,576 contexts per library in an order its HashMap
    tree bins decide; bammetrics exits 1 naming the option
    (KNOWN_DIFFERENCES.md). A short CONTEXTS list of that size runs."""
    bam, ref = built["ox_damage"]
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={tmp_path / 'b.txt'}", f"R={ref}",
                             "CONTEXT_SIZE=5"], capture_output=True, text=True)
    assert result.returncode == 1 and "CONTEXT_SIZE=5" in result.stderr
    assert not (tmp_path / "b.txt").exists()
    compare(bam, ref, ["CONTEXT_SIZE=5", "CONTEXTS=AAAAACAAAAA", "CONTEXTS=GGTACCGTACG"],
            "ox_damage context_size_5 contexts")


def test_bed_intervals_are_read_by_bammetrics_only(corpus, tmp_path):
    """Picard 3.5.0 rejects a BED as INTERVALS; bammetrics converts it as
    BedToIntervalList would (KNOWN_DIFFERENCES.md)."""
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t1000\t3000\ta\nchr2\t100\t900\tc\n")
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={corpus / 'corpus.bam'}", f"O={tmp_path / 'p.txt'}",
                        f"R={corpus / 'ref.fa'}", SILENT, f"INTERVALS={bed}", "QUIET=true", "VERBOSITY=ERROR"],
                       check=True, env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    converted = tmp_path / "t.interval_list"
    subprocess.run([str(BAMMETRICS), "BedToIntervalList", f"I={bed}", f"SD={corpus / 'ref.dict'}",
                    f"O={converted}"], check=True, capture_output=True)
    cells = compare_metrics(TOOL, corpus / "corpus.bam", [SILENT, f"R={corpus / 'ref.fa'}", f"INTERVALS={converted}"],
                            bam_args=[f"R={corpus / 'ref.fa'}", f"INTERVALS={bed}"], workdir=tmp_path,
                            expect_file=True)
    assert cells > 0


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Picard files the comparisons above hand to SEEN when every one of them runs.
EXPECTED_FILES = (len(INPUTS) + len(PERMISSIVE_INPUTS) + len(OPTIONS) * len(OPTION_INPUTS) + 2 + 3
                  + len(INTERVAL_NAMES) * 2 + len(DBSNP_NAMES) + 1 + 1 + 1)


def narrowed(config):
    """True when -k, --deselect or a node id selected part of the module."""
    return bool(config.getoption("keyword") or config.getoption("deselect")
                or any("::" in arg for arg in config.args))


def test_every_column_takes_three_non_default_values(request):
    if len(SEEN) != EXPECTED_FILES:
        message = f"{len(SEEN)} of {EXPECTED_FILES} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(f"{message}: a comparison above failed or EXPECTED_FILES is stale")
    values = defaultdict(set)
    columns = []
    for parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns:
            if column not in columns:
                columns.append(column)
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3}
    assert not short, f"columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

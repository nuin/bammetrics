"""CollectWgsMetrics, CollectRawWgsMetrics and CollectWgsMetricsWithNonZeroCoverage
against Picard 3.5.0: the metrics row, the coverage histogram and the base-quality
histogram, every option that changes a number, both pileup algorithms, and error
parity.

Inputs: the bwa tiers t0 to t3 with the corpus reference, the foundation corpus with
its own reference, and the asm_bams, read_level_bams and complexity_bams BAMs, each
set with its reference. Every case runs Picard with VALIDATION_STRINGENCY=SILENT.
INTERVALS cases take an interval list that bammetrics BedToIntervalList wrote from a
BED, so both tools read the same file.

CollectWgsMetricsWithNonZeroCoverage with USE_FAST_ALGORITHM=true is an error in
both tools once a record reaches a locus (Picard's cast fails); only the two empty
BAMs compare files with it.

The last tests check that every column took at least 3 distinct non-default values
over all Picard files compared, except the columns in FEWER_THAN_THREE with the
reason. They fail when a comparison above did not hand over its Picard file, and
skip only under -k, --deselect or a node id.
"""
import subprocess
from collections import defaultdict

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, REFERENCE
from compare import compare_metrics
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

WGS = "CollectWgsMetrics"
RAW = "CollectRawWgsMetrics"
NON_ZERO = "CollectWgsMetricsWithNonZeroCoverage"
TOOLS = [WGS, RAW, NON_ZERO]
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = {tool: [] for tool in TOOLS}


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Every hand-built BAM, each set with its reference."""
    bams = {}
    for prefix, build in (("asm", build_asm), ("rl", build_read_level), ("elc", build_complexity)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = (bam, bam.parent / "ref.fa")
    return bams


@pytest.fixture(scope="module")
def interval_lists(tmp_path_factory, corpus):
    """Interval lists written by bammetrics BedToIntervalList: the tier baits from
    baits.bed, and a hand-written BED of the foundation corpus with an interval
    on each contig, an abutting pair and a zero-length one."""
    directory = tmp_path_factory.mktemp("intervals")
    lists = {}
    foundation_bed = directory / "foundation.bed"
    foundation_bed.write_text("chr1\t1000\t3000\ta\nchr1\t3000\t3500\tb\nchr1\t20000\t20000\tz\n"
                              "chr2\t100\t900\tc\nchrM\t0\t12000\td\n")
    for name, bed, dictionary in (("tiers", CORPUS / "baits.bed", CORPUS / "reference.dict"),
                                  ("foundation", foundation_bed, corpus / "ref.dict")):
        out = directory / f"{name}.interval_list"
        subprocess.run([str(BAMMETRICS), "BedToIntervalList", f"I={bed}", f"SD={dictionary}", f"O={out}"],
                       check=True, capture_output=True)
        lists[name] = out
    return lists


def inputs(name, built, corpus):
    """(bam, reference) of a named input."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam", REFERENCE


def compare(tool, bam, ref, args, label):
    cells = compare_metrics(tool, bam, [SILENT, f"R={ref}", *args], workdir=bam.parent,
                            expect_file=True, seen=SEEN[tool])
    assert cells > 0, f"{tool} {label}: no cell compared"
    RESULTS.append((f"{tool} {label}", cells, cells, "every cell and the file body identical"))
    return cells


# --- every input at the defaults -------------------------------------------------

# t0_asm has an unpaired read running past the end of chrM: the default pileup
# drops its loci past the contig, the fast one fails in both tools once
# COUNT_UNPAIRED=true lets the read through (below). t0_filters has
# two @SQ lines against a three-contig reference: both tools fail (below). An
# input without a base above Q2 (t0_asm_empty, asm_empty, rl_empty, elc_unaligned,
# asm_all_qcfail, asm_single_unmapped, rl_all_duplicates) costs Picard about 50 s:
# every draw of the sensitivity sampler is rejected 600 times.
TIERS = ["t0_asm", "t0_asm_cycles", "t0_asm_empty", "t0_hs", "t0_qy", "t1", "t2", "t3"]
BUILT = ["asm_mix", "asm_unpaired_only", "asm_no_pf_group", "asm_all_qcfail", "asm_single_unmapped",
         "asm_clips_and_flags", "asm_no_read_groups", "asm_some_without_rg", "asm_empty", "rl_reads",
         "rl_star_fields", "rl_all_duplicates", "rl_empty", "rl_no_read_groups", "rl_second_only",
         "elc_illumina", "elc_unaligned", "elc_gc_skew"]
INPUTS = TIERS + ["foundation"] + BUILT


@pytest.mark.parametrize("name", INPUTS)
def test_wgs_defaults(built, corpus, name):
    bam, ref = inputs(name, built, corpus)
    compare(WGS, bam, ref, [], name)


# CollectRawWgsMetrics spends most of its time in the sensitivity simulation
# (a 100,001-bin depth histogram), so it runs on fewer inputs.
RAW_INPUTS = ["t1", "t3", "foundation", "asm_mix", "asm_clips_and_flags", "rl_reads", "rl_star_fields",
              "elc_illumina"]


@pytest.mark.parametrize("name", RAW_INPUTS)
def test_raw_defaults(built, corpus, name):
    bam, ref = inputs(name, built, corpus)
    compare(RAW, bam, ref, [], name)


NON_ZERO_INPUTS = ["t0_asm", "t1", "t3", "foundation", "asm_mix", "asm_unpaired_only",
                   "asm_clips_and_flags", "rl_reads", "rl_star_fields", "rl_empty", "elc_illumina",
                   "elc_gc_skew"]


@pytest.mark.parametrize("name", NON_ZERO_INPUTS)
def test_non_zero_defaults(built, corpus, name):
    bam, ref = inputs(name, built, corpus)
    compare(NON_ZERO, bam, ref, [], name)


# --- options ---------------------------------------------------------------------------

OPTIONS = {
    "mq_0": ["MQ=0"],
    "mq_60": ["MINIMUM_MAPPING_QUALITY=60"],
    "q_0": ["Q=0"],
    "q_3": ["Q=3"],
    "q_30": ["MINIMUM_BASE_QUALITY=30"],
    "cap_1": ["CAP=1"],
    "cap_5": ["COVERAGE_CAP=5"],
    "cap_50": ["CAP=50"],
    "cap_1000": ["CAP=1000"],
    "locus_cap_2": ["CAP=2", "LOCUS_ACCUMULATION_CAP=2"],
    "locus_cap_5": ["CAP=5", "LOCUS_ACCUMULATION_CAP=5"],
    "locus_cap_below_cap": ["LOCUS_ACCUMULATION_CAP=1"],  # raised to COVERAGE_CAP
    "count_unpaired": ["COUNT_UNPAIRED=true"],
    "unpaired_mq0_q0": ["COUNT_UNPAIRED=true", "MQ=0", "Q=0"],
    "bq_histogram": ["INCLUDE_BQ_HISTOGRAM=true"],
    "sample_100": ["SAMPLE_SIZE=100"],
    "sample_1": ["SAMPLE_SIZE=1"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_1765": ["STOP_AFTER=1765"],
    "stop_after_100000": ["STOP_AFTER=100000"],
    "stop_after_0": ["STOP_AFTER=0"],  # reads everything
    "fast": ["USE_FAST_ALGORITHM=true"],
    "fast_bq_unpaired": ["USE_FAST_ALGORITHM=true", "INCLUDE_BQ_HISTOGRAM=true", "COUNT_UNPAIRED=true"],
    "fast_q0_cap_3": ["USE_FAST_ALGORITHM=true", "Q=0", "CAP=3"],
    "fast_stop_after_1765": ["USE_FAST_ALGORITHM=true", "STOP_AFTER=1765"],
    "fast_stop_after_30000": ["USE_FAST_ALGORITHM=true", "STOP_AFTER=30000"],
    # The default READ_LENGTH=150 window (300,000 loci) slides once on t1's chr1;
    # 1000 never slides, and anything under the read length is an error (below).
    "fast_read_length_1000": ["USE_FAST_ALGORITHM=true", "READ_LENGTH=1000"],
    "read_length_1_default_algorithm": ["READ_LENGTH=1"],  # unused by the default pileup
}
OPTION_INPUTS = ["foundation", "asm_mix", "rl_reads", "t1", "elc_illumina"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_wgs_options(built, corpus, option, name):
    bam, ref = inputs(name, built, corpus)
    compare(WGS, bam, ref, OPTIONS[option], f"{name} {option}")


RAW_OPTIONS = {
    "mq_20_q_20": ["MQ=20", "Q=20"],
    "cap_10_sample_100": ["CAP=10", "SAMPLE_SIZE=100"],
    "cap_300_bq": ["CAP=300", "INCLUDE_BQ_HISTOGRAM=true", "COUNT_UNPAIRED=true"],
    "fast_sample_100": ["USE_FAST_ALGORITHM=true", "SAMPLE_SIZE=100"],
    "stop_after_5000": ["STOP_AFTER=5000", "SAMPLE_SIZE=100"],
}


@pytest.mark.parametrize("name", ["foundation", "asm_mix", "t1"])
@pytest.mark.parametrize("option", RAW_OPTIONS)
def test_raw_options(built, corpus, option, name):
    bam, ref = inputs(name, built, corpus)
    compare(RAW, bam, ref, RAW_OPTIONS[option], f"{name} {option}")


# USE_FAST_ALGORITHM=true is an error in both tools whenever a record reaches a
# locus (test_non_zero_fast_algorithm_fails_in_both); the two empty BAMs run it.
NON_ZERO_OPTIONS = {
    "bq_chart": ["INCLUDE_BQ_HISTOGRAM=true", "CHART=/dev/null"],
    "mq_0_q_0_cap_5": ["MQ=0", "Q=0", "CAP=5"],
    "unpaired_sample_100": ["COUNT_UNPAIRED=true", "SAMPLE_SIZE=100"],
    "stop_after_1765": ["STOP_AFTER=1765", "CHART_OUTPUT=/dev/null"],
}
NON_ZERO_OPTION_INPUTS = ["foundation", "asm_mix", "rl_reads", "t1"]
NON_ZERO_FAST_EMPTY_INPUTS = ["rl_empty", "asm_empty"]


@pytest.mark.parametrize("name", NON_ZERO_OPTION_INPUTS)
@pytest.mark.parametrize("option", NON_ZERO_OPTIONS)
def test_non_zero_options(built, corpus, option, name):
    bam, ref = inputs(name, built, corpus)
    compare(NON_ZERO, bam, ref, NON_ZERO_OPTIONS[option], f"{name} {option}")


@pytest.mark.parametrize("name", NON_ZERO_FAST_EMPTY_INPUTS)
def test_non_zero_fast_algorithm_on_an_empty_bam(built, corpus, name):
    """No record reaches a locus, so Picard's collector never casts one and
    both tools write a file with every depth 0."""
    bam, ref = inputs(name, built, corpus)
    compare(NON_ZERO, bam, ref, ["USE_FAST_ALGORITHM=true"], f"{name} fast")


# --- intervals ---------------------------------------------------------------------------

INTERVAL_CASES = {
    "intervals": [],
    "intervals_fast": ["USE_FAST_ALGORITHM=true"],
    "intervals_stop_after_100": ["STOP_AFTER=100"],
    "intervals_fast_stop_after_100": ["USE_FAST_ALGORITHM=true", "STOP_AFTER=100"],
    "intervals_bq_unpaired": ["INCLUDE_BQ_HISTOGRAM=true", "COUNT_UNPAIRED=true", "MQ=0"],
}
FAST_INTERVAL_CASES = ["intervals_fast", "intervals_fast_stop_after_100"]
# (tool, case) pairs that compare files: the NonZero fast cases fail in both tools.
FOUNDATION_INTERVAL_CASES = [(tool, case) for tool in TOOLS for case in INTERVAL_CASES
                             if not (tool == NON_ZERO and case in FAST_INTERVAL_CASES)]
TIER_INTERVAL_INPUTS = ["t1", "t3", "t0_hs"]
TIER_INTERVAL_CASES = ["intervals", "intervals_fast", "intervals_bq_unpaired"]


@pytest.mark.parametrize("tool, case", FOUNDATION_INTERVAL_CASES)
def test_foundation_intervals(corpus, interval_lists, tool, case):
    args = [f"INTERVALS={interval_lists['foundation']}", *INTERVAL_CASES[case]]
    if tool == RAW:
        args.append("SAMPLE_SIZE=100")
    compare(tool, corpus / "corpus.bam", corpus / "ref.fa", args, f"foundation {case}")


@pytest.mark.parametrize("name", TIER_INTERVAL_INPUTS)
@pytest.mark.parametrize("case", TIER_INTERVAL_CASES)
def test_tier_intervals(interval_lists, name, case):
    args = [f"INTERVALS={interval_lists['tiers']}", *INTERVAL_CASES[case]]
    compare(WGS, CORPUS / f"{name}.bam", REFERENCE, args, f"{name} {case}")


def test_picard_interval_lists_are_read_too(corpus):
    compare(WGS, corpus / "corpus.bam", corpus / "ref.fa",
            [f"INTERVALS={corpus / 'target.interval_list'}", "USE_FAST_ALGORITHM=true"],
            "foundation picard target list fast")
    compare(WGS, corpus / "corpus.bam", corpus / "ref.fa", [f"INTERVALS={corpus / 'bait.interval_list'}"],
            "foundation picard bait list")


# --- error parity ---------------------------------------------------------------------

def both_fail(tool, bam, args, tmp_path, message=None):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", tool, f"I={bam}", f"O={picard_out}", SILENT, *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), tool, f"I={bam}", f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{tool} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


def test_without_reference_fails(built, tmp_path):
    bam, _ = built["rl_reads"]
    both_fail(WGS, bam, [], tmp_path, "REFERENCE_SEQUENCE")


def test_assume_sorted_is_not_an_option(built, tmp_path):
    """Picard names the option as given: `Unrecognized option: AS`."""
    bam, ref = built["rl_reads"]
    both_fail(WGS, bam, [f"R={ref}", "ASSUME_SORTED=true"], tmp_path, "Unrecognized option: ASSUME_SORTED")
    both_fail(RAW, bam, [f"R={ref}", "AS=false"], tmp_path, "Unrecognized option: AS")


def test_queryname_order_fails(built, tmp_path):
    bam, ref = built["rl_queryname"]
    both_fail(WGS, bam, [f"R={ref}"], tmp_path, "coordinate sorted")
    both_fail(WGS, bam, [f"R={ref}", "USE_FAST_ALGORITHM=true"], tmp_path, "coordinate sorted")


def test_coverage_cap_and_sample_size_must_be_positive(built, tmp_path):
    bam, ref = built["rl_reads"]
    both_fail(WGS, bam, [f"R={ref}", "CAP=0"], tmp_path, "Coverage cap must be positive")
    both_fail(NON_ZERO, bam, [f"R={ref}", "SAMPLE_SIZE=0"], tmp_path, "Sample size is required")


@pytest.mark.parametrize("read_length", [0, 1, 2])
def test_fast_algorithm_window_too_small_fails(built, corpus, tmp_path, read_length):
    """READ_LENGTH sizes the fast collector's counter window (READ_LENGTH x
    2000 loci) and the window slides only when fewer than READ_LENGTH loci are
    left, so a 150-base read near the end of a 4,000-locus window is an index
    error in both tools."""
    for name in ("foundation", "rl_reads", "t1"):
        bam, ref = inputs(name, built, corpus)
        both_fail(WGS, bam, [f"R={ref}", "USE_FAST_ALGORITHM=true", f"READ_LENGTH={read_length}"], tmp_path,
                  "out of counter bounds")


def test_fast_algorithm_on_a_read_past_the_contig_end(tmp_path):
    """t0_asm's u_past_end is unpaired and runs 40 bases past chrM. At the
    default COUNT_UNPAIRED=false the pair filter drops it and both tools agree;
    with COUNT_UNPAIRED=true the fast collector reads the reference at every
    base of the block and Picard throws an index error."""
    compare(WGS, CORPUS / "t0_asm.bam", REFERENCE, ["USE_FAST_ALGORITHM=true"], "t0_asm fast")
    both_fail(WGS, CORPUS / "t0_asm.bam", [f"R={REFERENCE}", "USE_FAST_ALGORITHM=true", "COUNT_UNPAIRED=true"],
              tmp_path, "past the end of the contig")


def test_reference_contigs_must_match_the_header(built, corpus, tmp_path):
    bam, _ = built["rl_reads"]
    both_fail(WGS, bam, [f"R={corpus / 'ref.fa'}"], tmp_path,
              "contig chr1 (60000 bp) in the reference, c1 (30000 bp) in the BAM header")
    both_fail(WGS, CORPUS / "t0_filters.bam", [f"R={REFERENCE}"], tmp_path, "in the BAM header")


def test_interval_list_dictionary_must_match(built, corpus, interval_lists, tmp_path):
    bam, ref = built["rl_reads"]
    both_fail(WGS, bam, [f"R={ref}", f"INTERVALS={interval_lists['foundation']}"], tmp_path,
              "contig chr1 (60000 bp) in the interval list, c1 (30000 bp) in the BAM header")


@pytest.mark.parametrize("name, case", [(name, "fast") for name in NON_ZERO_OPTION_INPUTS]
                         + [("foundation", case) for case in FAST_INTERVAL_CASES])
def test_non_zero_fast_algorithm_fails_in_both(built, corpus, interval_lists, tmp_path, name, case):
    """Picard 3.5.0 keeps the default collector for this tool whatever
    USE_FAST_ALGORITHM says, but still reads with the edge iterator, so the
    collector's cast of the first event fails (KNOWN_DIFFERENCES.md)."""
    bam, ref = inputs(name, built, corpus)
    args = [f"R={ref}"]
    if case in FAST_INTERVAL_CASES:
        args += [f"INTERVALS={interval_lists['foundation']}", *INTERVAL_CASES[case]]
    else:
        args.append("USE_FAST_ALGORITHM=true")
    both_fail(NON_ZERO, bam, args, tmp_path, "USE_FAST_ALGORITHM=true")


def test_theoretical_sensitivity_output_is_not_implemented(corpus, tmp_path):
    """Picard writes the second file; bammetrics exits 1 naming the option
    (KNOWN_DIFFERENCES.md). ALLELE_FRACTION only feeds that file and is
    rejected the same way."""
    picard_out = tmp_path / "picard.txt"
    subprocess.run(["picard", WGS, f"I={corpus / 'corpus.bam'}", f"O={picard_out}", f"R={corpus / 'ref.fa'}",
                    SILENT, f"THEORETICAL_SENSITIVITY_OUTPUT={tmp_path / 'ts.txt'}", "ALLELE_FRACTION=0.1",
                    "QUIET=true", "VERBOSITY=ERROR"], check=True, env=picard_env(),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    assert picard_out.exists() and (tmp_path / "ts.txt").exists()
    for option in (f"THEORETICAL_SENSITIVITY_OUTPUT={tmp_path / 'b.ts.txt'}", "ALLELE_FRACTION=0.1"):
        result = subprocess.run([str(BAMMETRICS), WGS, f"I={corpus / 'corpus.bam'}", f"O={tmp_path / 'b.txt'}",
                                 f"R={corpus / 'ref.fa'}", option], capture_output=True, text=True)
        assert result.returncode != 0 and option.split("=")[0] in result.stderr
        assert not (tmp_path / "b.txt").exists()


def test_bed_intervals_are_read_by_bammetrics_only(corpus, tmp_path):
    """Picard 3.5.0 rejects a BED as INTERVALS; bammetrics converts it as
    BedToIntervalList would (KNOWN_DIFFERENCES.md)."""
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t1000\t3000\ta\nchr2\t100\t900\tc\n")
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", WGS, f"I={corpus / 'corpus.bam'}", f"O={tmp_path / 'p.txt'}",
                        f"R={corpus / 'ref.fa'}", SILENT, f"INTERVALS={bed}", "QUIET=true", "VERBOSITY=ERROR"],
                       check=True, env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    converted = tmp_path / "t.interval_list"
    subprocess.run([str(BAMMETRICS), "BedToIntervalList", f"I={bed}", f"SD={corpus / 'ref.dict'}",
                    f"O={converted}"], check=True, capture_output=True)
    cells = compare_metrics(WGS, corpus / "corpus.bam", [SILENT, f"R={corpus / 'ref.fa'}", f"INTERVALS={converted}"],
                            bam_args=[SILENT, f"R={corpus / 'ref.fa'}", f"INTERVALS={bed}"], workdir=tmp_path,
                            expect_file=True)
    assert cells > 0


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {
    (NON_ZERO, "CATEGORY"): "WHOLE_GENOME or NON_ZERO_REGIONS by definition",
}
# Picard files the comparisons above hand to SEEN when every one of them runs.
EXPECTED_FILES = {
    WGS: (len(INPUTS) + len(OPTIONS) * len(OPTION_INPUTS) + len(INTERVAL_CASES)
          + len(TIER_INTERVAL_INPUTS) * len(TIER_INTERVAL_CASES) + 2 + 1),
    RAW: len(RAW_INPUTS) + len(RAW_OPTIONS) * 3 + len(INTERVAL_CASES),
    NON_ZERO: (len(NON_ZERO_INPUTS) + len(NON_ZERO_OPTIONS) * len(NON_ZERO_OPTION_INPUTS)
               + len(NON_ZERO_FAST_EMPTY_INPUTS) + len(INTERVAL_CASES) - len(FAST_INTERVAL_CASES)),
}


def narrowed(config):
    """True when -k, --deselect or a node id selected part of the module."""
    return bool(config.getoption("keyword") or config.getoption("deselect")
                or any("::" in arg for arg in config.args))


@pytest.mark.parametrize("tool", TOOLS)
def test_every_column_takes_three_non_default_values(request, tool):
    files = SEEN[tool]
    if len(files) != EXPECTED_FILES[tool]:
        message = f"{tool}: {len(files)} of {EXPECTED_FILES[tool]} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(f"{message}: a comparison above failed or EXPECTED_FILES is stale")
    values = defaultdict(set)
    columns = []
    for parsed in files:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns:
            if column not in columns:
                columns.append(column)
    short = {c: sorted(values[c]) for c in columns
             if len(values[c]) < 3 and (tool, c) not in FEWER_THAN_THREE}
    assert not short, f"{tool}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{tool} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

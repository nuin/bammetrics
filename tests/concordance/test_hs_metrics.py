"""CollectHsMetrics and CollectTargetedPcrMetrics against Picard 3.5.0: the metrics
row at every accumulation level, the two histograms, every option that changes a
number, PER_TARGET_COVERAGE, PER_BASE_COVERAGE, THEORETICAL_SENSITIVITY_OUTPUT
with ALLELE_FRACTION, and error parity.

Inputs: the bwa tiers t0_hs, t1, t2 and t3 with the corpus baits and targets, the
foundation corpus, the exome slice (tests/data/make_exome.py: 150 exons, tiled baits,
a GC-skewed target subset), the hs_bams.py BAMs (level shapes, adapter reads, the
deep targets) and a few asm_bams and read_level_bams BAMs on the hs_bams panel.
Every case runs Picard with VALIDATION_STRINGENCY=SILENT. Interval lists for the
hand-built BAMs are written by bammetrics BedToIntervalList from the BEDs, so both
tools read the same file.

Levels: Picard 3.5.0 shares the four histograms between units and throws for every
unit that finishes after a unit with a counted base (KNOWN_DIFFERENCES.md), so the
level cases that compare files are the shapes that complete (hs_bams.py); the
others are error parity.

The last tests check that every column took at least 3 distinct non-default values
over all Picard files compared, except the columns in FEWER_THAN_THREE with the
reason.
"""
import gzip
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, REFERENCE
from compare import body, compare_metrics, diff_metrics
from conftest import RESULTS
from hs_bams import build_all as build_hs
from picard_metrics import parse_picard_metrics
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

HS = "CollectHsMetrics"
PCR = "CollectTargetedPcrMetrics"
TOOLS = [HS, PCR]
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = {tool: [] for tool in TOOLS}
# Picard's sensitivity sampler costs 600 rejected draws per sample on a unit
# without a counted base; a small SAMPLE_SIZE keeps those cases quick.
SMALL = "SAMPLE_SIZE=100"


def bed_to_interval_list(bed, dictionary, out):
    subprocess.run([str(BAMMETRICS), "BedToIntervalList", f"I={bed}", f"SD={dictionary}", f"O={out}"],
                   check=True, capture_output=True)
    return out


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Every hand-built BAM as (bam, reference, baits, targets); the hs_bams panel
    and baits serve the asm and read-level BAMs, which share its reference."""
    bams = {}
    hs_dir = tmp_path_factory.mktemp("hs")
    hs = build_hs(hs_dir)
    panel = bed_to_interval_list(hs_dir / "hs_panel.bed", hs["first_empty"], hs_dir / "panel.interval_list")
    baits = bed_to_interval_list(hs_dir / "hs_baits.bed", hs["first_empty"], hs_dir / "baits.interval_list")
    deep = bed_to_interval_list(hs_dir / "deep_panel.bed", hs["deep"], hs_dir / "deep.interval_list")
    for name, bam in hs.items():
        lists = (deep, deep) if name == "deep" else (baits, panel)
        bams[f"hs_{name}"] = (bam, hs_dir / "ref.fa", *lists)
    for prefix, build in (("asm", build_asm), ("rl", build_read_level)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = (bam, directory / "ref.fa", baits, panel)
    return bams


def inputs(name, built, corpus, exome, targets=None):
    """(bam, reference, baits, targets) of a named input; `targets` picks another
    target list of the same corpus (`hs_targets` for the tiers, `gc_targets` for
    the exome)."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa", corpus / "bait.interval_list", corpus / "target.interval_list"
    if name == "exome":
        return (exome / "exome.bam", exome / "ref.fa", exome / "baits.interval_list",
                exome / f"{targets or 'targets'}.interval_list")
    if name in built:
        return built[name]
    return (CORPUS / f"{name}.bam", REFERENCE, CORPUS / "baits.interval_list",
            CORPUS / f"{targets or 'targets'}.interval_list")


def probe_args(tool, baits, targets):
    key = "BI" if tool == HS else "AI"
    return [f"{key}={baits}", f"TI={targets}"]


def compare(tool, bam, args, label, bam_args=None):
    cells = compare_metrics(tool, bam, [SILENT, *args], bam_args=None if bam_args is None else [SILENT, *bam_args],
                            workdir=bam.parent, expect_file=True, seen=SEEN[tool])
    assert cells > 0, f"{tool} {label}: no cell compared"
    RESULTS.append((f"{tool} {label}", cells, cells, "every cell and the file body identical"))
    return cells


# --- every input at the defaults, with and without the reference ----------------------

TIERS = ["t0_hs", "t1", "t2", "t3"]
BUILT = ["hs_first_empty", "hs_unknown", "hs_no_coverage", "hs_deep", "asm_mix", "asm_clips_and_flags",
         "asm_no_read_groups", "asm_some_without_rg", "rl_reads"]
INPUTS = TIERS + ["foundation", "exome"] + BUILT


@pytest.mark.parametrize("name", INPUTS)
def test_hs_defaults(built, corpus, exome, name):
    bam, _, baits, targets = inputs(name, built, corpus, exome)
    small = [SMALL] if name.startswith("hs_") else []
    compare(HS, bam, [*probe_args(HS, baits, targets), *small], name)


REFERENCE_INPUTS = ["t0_hs", "t1", "t2", "t3", "foundation", "exome", "hs_first_empty", "hs_deep", "asm_mix",
                    "rl_reads"]


@pytest.mark.parametrize("name", REFERENCE_INPUTS)
def test_hs_with_reference(built, corpus, exome, name):
    """R= adds AT_DROPOUT and GC_DROPOUT."""
    bam, ref, baits, targets = inputs(name, built, corpus, exome)
    small = [SMALL] if name.startswith("hs_") else []
    compare(HS, bam, [*probe_args(HS, baits, targets), f"R={ref}", *small], f"{name} R=")


OTHER_TARGETS = [("t0_hs", "hs_targets"), ("t2", "hs_targets"), ("exome", "gc_targets")]


@pytest.mark.parametrize("name, targets", OTHER_TARGETS)
def test_hs_other_target_lists(built, corpus, exome, name, targets):
    """hs_targets: unsorted, duplicate, overlapping, abutting and empty targets;
    gc_targets: the exons under 35% or over 65% GC."""
    bam, ref, baits, target_list = inputs(name, built, corpus, exome, targets)
    compare(HS, bam, [*probe_args(HS, baits, target_list), f"R={ref}"], f"{name} {targets}")


# --- options ------------------------------------------------------------------------------

OPTIONS = {
    "mq_0": ["MQ=0"],
    "mq_60": ["MINIMUM_MAPPING_QUALITY=60"],
    "q_0": ["Q=0"],
    "q_30": ["MINIMUM_BASE_QUALITY=30"],
    "mq_0_q_0": ["MQ=0", "Q=0"],
    "no_clip": ["CLIP_OVERLAPPING_READS=false"],
    "clip": ["CLIP_OVERLAPPING_READS=true"],
    "indels": ["INCLUDE_INDELS=true"],
    "cap_20": ["COVERAGE_CAP=20"],
    "cap_5": ["covMax=5"],
    "near_10": ["NEAR_DISTANCE=10"],
    "near_0": ["NEAR_DISTANCE=0"],
    "near_-5": ["NEAR_DISTANCE=-5"],
    "sample_100": [SMALL],
    "named": ["BAIT_SET_NAME=kit_v1"],
    "n_short": ["N=panel.v3"],
}
OPTION_INPUTS = ["t2", "exome", "t0_hs"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_hs_options(built, corpus, exome, option, name):
    bam, ref, baits, targets = inputs(name, built, corpus, exome)
    compare(HS, bam, [*probe_args(HS, baits, targets), f"R={ref}", *OPTIONS[option]], f"{name} {option}")


DEEP_OPTIONS = {"default": [], "mq_0": ["MQ=0"], "q_0": ["Q=0"], "mq_0_q_0": ["MQ=0", "Q=0"], "cap_1000": ["COVERAGE_CAP=1000"]}


@pytest.mark.parametrize("option", DEEP_OPTIONS)
def test_hs_deep_targets(built, corpus, exome, option):
    """Three 40-base targets whose depth falls from 160,000x along the target
    (hs_bams.py): PCT_TARGET_BASES_500X to _100000X count a different number of
    bases at each filter."""
    bam, ref, baits, targets = inputs("hs_deep", built, corpus, exome)
    compare(HS, bam, [*probe_args(HS, baits, targets), f"R={ref}", SMALL, *DEEP_OPTIONS[option]], f"deep {option}")


def test_hs_two_bait_files(corpus, tmp_path):
    """BAIT_INTERVALS repeats; BAIT_SET joins the file stems."""
    bam, ref, baits, targets = inputs("t0_hs", {}, corpus, None)
    compare(HS, bam, [f"BI={baits}", f"BAIT_INTERVALS={targets}", f"TI={targets}", f"R={ref}"], "t0_hs two bait files")


# --- the extra files ------------------------------------------------------------------------

def compare_with_file(tool, bam, args, option, label, suffix=".txt", gzipped=False, tmp_path=None):
    """compare_metrics with a second file named by `option`, compared byte for byte
    (after decompression when `gzipped`); its lines count as cells."""
    picard_file, ours_file = tmp_path / f"picard{suffix}", tmp_path / f"ours{suffix}"
    cells = compare(tool, bam, [*args, f"{option}={picard_file}"], label, bam_args=[*args, f"{option}={ours_file}"])
    assert picard_file.exists() and ours_file.exists(), f"{label}: {option} not written by both tools"
    read = gzip.open if gzipped else open
    with read(picard_file, "rb") as p, read(ours_file, "rb") as o:
        picard_bytes, ours_bytes = p.read(), o.read()
    assert picard_bytes == ours_bytes, f"{label}: {option} differs"
    lines = picard_bytes.count(b"\n")
    RESULTS.append((f"{tool} {label} {option}", lines, lines, "lines identical"))
    return cells + lines


PER_TARGET_CASES = [("t1", "targets", []), ("t2", "hs_targets", []), ("t0_hs", "hs_targets", []),
                    ("foundation", None, []), ("exome", "targets", []), ("exome", "gc_targets", []),
                    ("exome", "targets", ["MQ=0", "Q=0"]), ("hs_deep", None, [SMALL]),
                    ("hs_first_empty", None, [SMALL])]


@pytest.mark.parametrize("name, targets, extra", PER_TARGET_CASES)
def test_hs_per_target_coverage(built, corpus, exome, tmp_path, name, targets, extra):
    bam, ref, baits, target_list = inputs(name, built, corpus, exome, targets)
    compare_with_file(HS, bam, [*probe_args(HS, baits, target_list), f"R={ref}", *extra], "PER_TARGET_COVERAGE",
                      f"{name} {targets or ''} {' '.join(extra)}".strip(), tmp_path=tmp_path)


PER_BASE_CASES = [("t1", "targets", ["MQ=0", "Q=0"], ".txt"), ("t2", "hs_targets", [], ".txt"),
                  ("exome", "gc_targets", [], ".txt"), ("exome", "targets", [], ".txt.gz"),
                  ("hs_unknown", None, [SMALL], ".txt")]


@pytest.mark.parametrize("name, targets, extra, suffix", PER_BASE_CASES)
def test_hs_per_base_coverage(built, corpus, exome, tmp_path, name, targets, extra, suffix):
    bam, ref, baits, target_list = inputs(name, built, corpus, exome, targets)
    compare_with_file(HS, bam, [*probe_args(HS, baits, target_list), f"R={ref}", *extra], "PER_BASE_COVERAGE",
                      f"{name} {targets or ''} {suffix}".strip(), suffix=suffix, gzipped=suffix.endswith(".gz"),
                      tmp_path=tmp_path)


def compare_sensitivity_file(tool, bam, args, label, tmp_path):
    """compare_metrics plus the THEORETICAL_SENSITIVITY_OUTPUT metrics file, parsed and
    compared cell by cell and byte for byte after the header lines."""
    picard_file, ours_file = tmp_path / "picard.ts.txt", tmp_path / "ours.ts.txt"
    option = "THEORETICAL_SENSITIVITY_OUTPUT"
    cells = compare(tool, bam, [*args, f"{option}={picard_file}"], label, bam_args=[*args, f"{option}={ours_file}"])
    assert picard_file.exists() and ours_file.exists(), f"{label}: {option} not written by both tools"
    picard, ours = parse_picard_metrics(picard_file), parse_picard_metrics(ours_file)
    assert picard.metrics_class == "picard.analysis.TheoreticalSensitivityMetrics"
    assert (picard.metrics_class, picard.columns, picard.rows) == (ours.metrics_class, ours.columns, ours.rows)
    assert body(picard_file) == body(ours_file), f"{label}: {option} body differs"
    extra = len(picard.rows) * len(picard.columns)
    RESULTS.append((f"{tool} {label} {option}", extra, extra, "every cell identical"))
    return cells + extra


SENSITIVITY_CASES = [
    ("t0_hs", [SMALL]),
    ("t1", []),
    ("t2", ["ALLELE_FRACTION=null", "ALLELE_FRACTION=0.3", "ALLELE_FRACTION=0.05", SMALL]),
    ("foundation", ["COVERAGE_CAP=20", SMALL]),
    ("exome", ["ALLELE_FRACTION=0.15", "SAMPLE_SIZE=1000"]),
    ("hs_first_empty", ["LEVEL=null", "LEVEL=SAMPLE", SMALL]),
]


@pytest.mark.parametrize("name, extra", SENSITIVITY_CASES)
def test_hs_theoretical_sensitivity_output(built, corpus, exome, tmp_path, name, extra):
    bam, _, baits, targets = inputs(name, built, corpus, exome)
    compare_sensitivity_file(HS, bam, [*probe_args(HS, baits, targets), *extra], f"{name} {' '.join(extra)}", tmp_path)


# --- accumulation levels --------------------------------------------------------------------

LEVEL_CASES = [
    ("hs_first_empty", ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("hs_first_empty", ["LEVEL=null", "LEVEL=LIBRARY"]),
    ("hs_first_empty", ["LEVEL=null", "METRIC_ACCUMULATION_LEVEL=READ_GROUP"]),
    ("hs_unknown", ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("hs_unknown", ["LEVEL=null", "LEVEL=READ_GROUP"]),
    ("hs_no_coverage", ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"]),
    ("hs_no_coverage", ["LEVEL=null", "LEVEL=READ_GROUP", "LEVEL=SAMPLE"]),
    ("exome", ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("t1", ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("t3", ["LEVEL=null", "LEVEL=SAMPLE"]),
]


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("name, levels", LEVEL_CASES)
def test_levels_that_complete(built, corpus, exome, tool, name, levels):
    """One unit with coverage, after any empty units; or every unit empty."""
    bam, ref, baits, targets = inputs(name, built, corpus, exome)
    small = [SMALL] if name.startswith("hs_") else []
    compare(tool, bam, [*probe_args(tool, baits, targets), f"R={ref}", *levels, *small], f"{name} {' '.join(levels)}")


# --- CollectTargetedPcrMetrics -------------------------------------------------------------

PCR_INPUTS = ["t0_hs", "t1", "t2", "t3", "foundation", "exome", "hs_first_empty", "hs_deep", "asm_mix", "rl_reads"]


@pytest.mark.parametrize("name", PCR_INPUTS)
def test_pcr_defaults(built, corpus, exome, name):
    bam, _, baits, targets = inputs(name, built, corpus, exome)
    small = [SMALL] if name.startswith("hs_") else []
    compare(PCR, bam, [*probe_args(PCR, baits, targets), *small], name)


PCR_OPTIONS = {
    "reference": [],
    "hs_defaults": ["MQ=20", "Q=20", "CLIP_OVERLAPPING_READS=true"],
    "clip": ["CLIP_OVERLAPPING_READS=true"],
    "mq_60": ["MINIMUM_MAPPING_QUALITY=60"],
    "q_30": ["Q=30"],
    "indels": ["INCLUDE_INDELS=true"],
    "cap_20": ["COVERAGE_CAP=20"],
    "near_0": ["NEAR_DISTANCE=0"],
    "near_-5": ["NEAR_DISTANCE=-5"],
    "named": ["CUSTOM_AMPLICON_SET_NAME=amp_v2"],
    "n_short": ["N=amp.v3"],
    "sample_100": [SMALL],
}
PCR_OPTION_INPUTS = ["t2", "exome"]


@pytest.mark.parametrize("name", PCR_OPTION_INPUTS)
@pytest.mark.parametrize("option", PCR_OPTIONS)
def test_pcr_options(built, corpus, exome, option, name):
    bam, ref, baits, targets = inputs(name, built, corpus, exome)
    compare(PCR, bam, [*probe_args(PCR, baits, targets), f"R={ref}", *PCR_OPTIONS[option]], f"{name} {option}")


PCR_DEEP_OPTIONS = {"mq_20": ["MQ=20"], "q_20": ["Q=20"], "mq_20_q_20": ["MQ=20", "Q=20"]}


@pytest.mark.parametrize("option", PCR_DEEP_OPTIONS)
def test_pcr_deep_targets(built, corpus, exome, option):
    """The deep ladder under the filters CollectTargetedPcrMetrics leaves off by
    default (MAPQ 1, base quality 0): a fifth of the reads and a seventh of the
    bases drop out, so every PCT_TARGET_BASES_*X column moves."""
    bam, ref, baits, targets = inputs("hs_deep", built, corpus, exome)
    compare(PCR, bam, [*probe_args(PCR, baits, targets), f"R={ref}", SMALL, *PCR_DEEP_OPTIONS[option]],
            f"deep {option}")


@pytest.mark.parametrize("name, targets", [("exome", "gc_targets"), ("t0_hs", "hs_targets"), ("hs_deep", None)])
def test_pcr_per_target_coverage(built, corpus, exome, tmp_path, name, targets):
    bam, ref, baits, target_list = inputs(name, built, corpus, exome, targets)
    small = [SMALL] if name.startswith("hs_") else []
    compare_with_file(PCR, bam, [*probe_args(PCR, baits, target_list), f"R={ref}", *small], "PER_TARGET_COVERAGE",
                      f"{name} {targets or ''}".strip(), tmp_path=tmp_path)


def test_pcr_per_base_and_sensitivity(built, corpus, exome, tmp_path):
    bam, ref, baits, targets = inputs("t2", built, corpus, exome)
    compare_with_file(PCR, bam, [*probe_args(PCR, baits, targets), f"R={ref}"], "PER_BASE_COVERAGE", "t2",
                      tmp_path=tmp_path)
    compare_sensitivity_file(PCR, bam, [*probe_args(PCR, baits, targets), "ALLELE_FRACTION=0.1", SMALL], "t2 af 0.1",
                             tmp_path)


# --- error parity --------------------------------------------------------------------------

def both_fail(tool, bam, args, tmp_path, message, picard_message, output_written=False):
    """Both tools exit non-zero with `message` (bammetrics) and `picard_message`
    (Picard) in their output. Neither writes OUTPUT, or with `output_written`
    both do and the files agree cell for cell (Picard writes the metrics file
    before it computes the sensitivities)."""
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    picard = subprocess.run(["picard", tool, f"I={bam}", f"O={picard_out}", SILENT, *args, "QUIET=true",
                             "VERBOSITY=ERROR"], env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr
    result = subprocess.run([str(BAMMETRICS), tool, f"I={bam}", f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    assert message in result.stderr, result.stderr
    label = f"{tool} {bam.name} {' '.join(args)}"
    if not output_written:
        assert not picard_out.exists() and not ours_out.exists()
        RESULTS.append((label, 0, 0, "both exit non-zero, no file"))
        return
    assert picard_out.exists() and ours_out.exists(), f"{label}: OUTPUT not written by both tools"
    picard_metrics, ours_metrics = parse_picard_metrics(picard_out), parse_picard_metrics(ours_out)
    problems = diff_metrics(picard_metrics, ours_metrics)
    assert not problems, f"{label}:\n  " + "\n  ".join(problems[:40])
    assert body(picard_out) == body(ours_out), f"{label}: OUTPUT body differs"
    cells = (len(picard_metrics.rows) * len(picard_metrics.columns)
             + len(picard_metrics.histogram) * len(picard_metrics.histogram_columns))
    assert cells > 0, f"{label}: no cell compared"
    SEEN[tool].append(picard_metrics)
    RESULTS.append((label, cells, cells, "both exit non-zero, OUTPUT identical"))


def test_stop_after_and_assume_sorted_are_not_options(built, tmp_path):
    bam, _, baits, targets = built["hs_first_empty"]
    both_fail(HS, bam, [*probe_args(HS, baits, targets), "STOP_AFTER=10"], tmp_path, "Unrecognized option: STOP_AFTER",
              "Unrecognized option: STOP_AFTER")
    both_fail(PCR, bam, [*probe_args(PCR, baits, targets), "AS=false"], tmp_path, "Unrecognized option: AS",
              "Unrecognized option: AS")


LEVELS_THAT_THROW = [
    ("foundation", HS, ["LEVEL=SAMPLE"]),
    ("foundation", PCR, ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("t1", HS, ["LEVEL=null", "LEVEL=LIBRARY"]),
    ("exome", HS, ["LEVEL=null", "LEVEL=READ_GROUP"]),
    ("hs_first_empty", HS, ["LEVEL=SAMPLE"]),
    ("hs_first_empty", PCR, ["LEVEL=null", "LEVEL=SAMPLE", "LEVEL=LIBRARY"]),
]
THROW_MESSAGE = "numbers of bases in the base quality histogram and the coverage histogram are not equal"


@pytest.mark.parametrize("name, tool, levels", LEVELS_THAT_THROW)
def test_a_unit_after_a_covered_unit_fails_in_both(built, corpus, exome, tmp_path, name, tool, levels):
    """Picard's shared-histogram check throws in the second unit with a counted base."""
    bam, _, baits, targets = inputs(name, built, corpus, exome)
    both_fail(tool, bam, [*probe_args(tool, baits, targets), *levels, SMALL], tmp_path, THROW_MESSAGE, THROW_MESSAGE)


def test_per_target_coverage_needs_all_reads_and_a_reference(built, tmp_path):
    bam, ref, baits, targets = built["hs_first_empty"]
    both_fail(HS, bam, [*probe_args(HS, baits, targets), f"PER_TARGET_COVERAGE={tmp_path / 'pt.txt'}"], tmp_path,
              "Must supply REFERENCE_SEQUENCE", "Must supply REFERENCE_SEQUENCE")
    both_fail(PCR, bam, [*probe_args(PCR, baits, targets), f"R={ref}", "LEVEL=SAMPLE",
                         f"PER_TARGET_COVERAGE={tmp_path / 'pt.txt'}"], tmp_path, "ALL_READS", "ALL_READS")
    assert not (tmp_path / "pt.txt").exists()


def test_mapped_read_without_qualities_fails(built, tmp_path):
    """A `*` QUAL on a mapped read: Picard indexes an empty quality array."""
    bam, _, baits, targets = built["rl_star_fields"]
    both_fail(HS, bam, [*probe_args(HS, baits, targets), SMALL], tmp_path, "aligned base without a quality",
              "ArrayIndexOutOfBoundsException")


def test_negative_coverage_cap_fails(built, tmp_path):
    """COVERAGE_CAP under 0 keys the unfiltered depth histogram below 0, so the
    base-count check fails at the first unit with a counted base."""
    bam, _, baits, targets = built["hs_first_empty"]
    both_fail(HS, bam, [*probe_args(HS, baits, targets), "COVERAGE_CAP=-1", SMALL], tmp_path, THROW_MESSAGE, THROW_MESSAGE)


def test_amplicon_intervals_once_and_targets_required(built, tmp_path):
    """bammetrics prints Picard's argument parser messages, option name included."""
    bam, _, baits, targets = built["hs_first_empty"]
    both_fail(PCR, bam, [f"AI={baits}", f"AI={targets}", f"TI={targets}"], tmp_path,
              "Option 'AMPLICON_INTERVALS' cannot be specified more than once.",
              "Option 'AI' cannot be specified more than once.")
    targets_required = "Option 'TARGET_INTERVALS' must be specified at least 1 times."
    both_fail(HS, bam, [f"BI={baits}"], tmp_path, targets_required, targets_required)
    both_fail(PCR, bam, [f"AI={baits}"], tmp_path, targets_required, targets_required)
    baits_required = "Option 'BAIT_INTERVALS' must be specified at least 1 times."
    both_fail(HS, bam, [f"TI={targets}"], tmp_path, baits_required, baits_required)
    amplicons_required = "Option 'AMPLICON_INTERVALS' is required."
    both_fail(PCR, bam, [f"TI={targets}"], tmp_path, amplicons_required, amplicons_required)


def test_interval_list_of_another_dictionary_fails(built, corpus, tmp_path):
    """The foundation targets name chr1 (60,000 bp); the hs_bams header has c1 (30,000 bp)."""
    bam, _, baits, _ = built["hs_first_empty"]
    both_fail(HS, bam, [f"BI={baits}", f"TI={corpus / 'target.interval_list'}"], tmp_path,
              "contig chr1 (60000 bp) in the interval list, c1 (30000 bp) in the BAM header",
              "does not match sequence dictionary for (target intervals)")


def test_allele_fraction_outside_0_1_fails_with_the_sensitivity_file(built, tmp_path):
    """Picard writes OUTPUT, then throws computing the sensitivities: the metrics
    file exists and agrees, the sensitivity file does not."""
    bam, _, baits, targets = built["hs_first_empty"]
    sensitivity = tmp_path / "ts.txt"
    both_fail(HS, bam, [*probe_args(HS, baits, targets), "ALLELE_FRACTION=1.5", SMALL,
                        f"THEORETICAL_SENSITIVITY_OUTPUT={sensitivity}"], tmp_path,
              "Allele fractions must be between 0 and 1", "Allele fractions must be between 0 and 1",
              output_written=True)
    assert not sensitivity.exists()


def test_bed_intervals_are_read_by_bammetrics_only(built, tmp_path):
    """Picard 3.5.0 rejects a BED; bammetrics converts it as BedToIntervalList would.
    The BED is copied to the interval list's stem, so BAIT_SET is `panel` on both sides."""
    bam, ref, baits, targets = built["hs_first_empty"]
    bed = tmp_path / "panel.bed"
    shutil.copy(bam.parent / "hs_panel.bed", bed)
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", HS, f"I={bam}", f"O={tmp_path / 'p.txt'}", SILENT, f"BI={bed}", f"TI={bed}",
                        "QUIET=true", "VERBOSITY=ERROR"], check=True, env=picard_env(), stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    cells = compare_metrics(HS, bam, [SILENT, f"BI={targets}", f"TI={targets}", f"R={ref}", SMALL],
                            bam_args=[SILENT, f"BI={bed}", f"TI={bed}", f"R={ref}", SMALL], workdir=tmp_path,
                            expect_file=True, seen=SEEN[HS])
    assert cells > 0
    RESULTS.append((f"{HS} hs_first_empty BED", cells, cells, "every cell and the file body identical"))


# --- every column exercised ---------------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {}
EXPECTED_FILES = {
    # The last 1 of each: the BED case (Hs) and the written OUTPUT of the
    # ALLELE_FRACTION=1.5 case (Hs); none for CollectTargetedPcrMetrics.
    HS: (len(INPUTS) + len(REFERENCE_INPUTS) + len(OTHER_TARGETS) + len(OPTIONS) * len(OPTION_INPUTS)
         + len(DEEP_OPTIONS) + 1 + len(PER_TARGET_CASES) + len(PER_BASE_CASES) + len(SENSITIVITY_CASES)
         + len(LEVEL_CASES) + 1 + 1),
    PCR: (len(PCR_INPUTS) + len(PCR_OPTIONS) * len(PCR_OPTION_INPUTS) + len(PCR_DEEP_OPTIONS) + 3 + 2
          + len(LEVEL_CASES)),
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

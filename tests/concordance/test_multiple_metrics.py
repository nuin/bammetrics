"""CollectMultipleMetrics against Picard 3.5.0: every file of every PROGRAM,
named as Picard names it, with the default list and every option that
changes a file (PROGRAM, LEVEL, EXTRA_ARGUMENT, FILE_EXTENSION, STOP_AFTER,
REFERENCE_SEQUENCE), plus error parity. Every case runs Picard with
VALIDATION_STRINGENCY=SILENT.

compare_outputs runs both tools with O=<dir>/<side> and compares every file
either wrote: the set of file names must be the same and each file is
compared cell by cell and byte for byte after the header lines.
"""
import subprocess

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, REFERENCE
from compare import compare_outputs
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

TOOL = "CollectMultipleMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
ALL_PROGRAMS = ["PROGRAM=CollectGcBiasMetrics", "PROGRAM=CollectQualityYieldMetrics"]
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    asm_dir = tmp_path_factory.mktemp("asm")
    read_dir = tmp_path_factory.mktemp("read_level")
    elc_dir = tmp_path_factory.mktemp("complexity")
    bams = {}
    for prefix, build, directory in (("asm", build_asm, asm_dir), ("rl", build_read_level, read_dir),
                                     ("elc", build_complexity, elc_dir)):
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = (bam, bam.parent / "ref.fa")
    return bams


def inputs(name, built, corpus):
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam", REFERENCE


def compare(bam, args, label, expected_files=None):
    before = len(SEEN)
    cells = compare_outputs(TOOL, bam, [SILENT, *args], (("O", ""),), workdir=bam.parent, seen=SEEN)
    assert cells > 0, f"{TOOL} {label}: no cell compared"
    written = {suffix for suffix, _ in SEEN[before:]}
    if expected_files is not None:
        # The insert size file is absent when the input has no counted pair.
        assert set(expected_files) - {".insert_size_metrics"} <= written <= set(expected_files), (
            f"{label}: files {sorted(written)}")
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every file identical"))
    return written


DEFAULT_FILES = [".alignment_summary_metrics", ".base_distribution_by_cycle_metrics",
                 ".insert_size_metrics", ".quality_by_cycle_metrics", ".quality_distribution_metrics"]
ALL_FILES = DEFAULT_FILES + [".gc_bias.detail_metrics", ".gc_bias.summary_metrics",
                             ".quality_yield_metrics"]

# --- the default programs, with and without a reference --------------------------

# rl_star_fields has a mapped read without qualities: CollectAlignmentSummaryMetrics
# throws in Picard and fails in bammetrics (below).
INPUTS = ["t0_asm", "t0_filters", "t0_no_pairs", "t1", "t2", "t3", "foundation", "asm_mix",
          "asm_unpaired_only", "asm_clips_and_flags", "rl_reads", "elc_illumina", "elc_gc_skew"]


@pytest.mark.parametrize("name", INPUTS)
def test_defaults_without_reference(built, corpus, name):
    bam, _ = inputs(name, built, corpus)
    compare(bam, [], f"{name} defaults", DEFAULT_FILES)


# CollectGcBiasMetrics fails on t0_asm (a read past its contig), t2 (over 100
# secondary records without bases) and the insert-size t0 BAMs (their own @SQ
# lengths), so both tools fail there (below).
REFERENCE_INPUTS = [name for name in INPUTS if name not in ("t0_asm", "t2", "t0_filters", "t0_no_pairs")]


@pytest.mark.parametrize("name", REFERENCE_INPUTS)
def test_every_program_with_reference(built, corpus, name):
    bam, ref = inputs(name, built, corpus)
    compare(bam, [f"R={ref}", *ALL_PROGRAMS], f"{name} all programs", ALL_FILES)


# --- options -------------------------------------------------------------------------

def option_cases(ref):
    return {
        "levels": [f"R={ref}", *ALL_PROGRAMS, "LEVEL=SAMPLE", "LEVEL=LIBRARY"],
        "read_group_level": [f"R={ref}", *ALL_PROGRAMS, "METRIC_ACCUMULATION_LEVEL=READ_GROUP"],
        "null_level": [f"R={ref}", *ALL_PROGRAMS, "LEVEL=null", "LEVEL=LIBRARY"],
        "extension": [f"R={ref}", *ALL_PROGRAMS, "FILE_EXTENSION=.txt"],
        "ext_short": ["EXT=.tsv"],
        "subset": ["PROGRAM=null", "PROGRAM=QualityScoreDistribution", "PROGRAM=MeanQualityByCycle"],
        "subset_repeated": ["PROGRAM=null", "PROGRAM=CollectInsertSizeMetrics",
                            "PROGRAM=CollectInsertSizeMetrics", "PROGRAM=CollectQualityYieldMetrics"],
        "stop_after_50": [f"R={ref}", *ALL_PROGRAMS, "STOP_AFTER=50"],
        "stop_after_1": ["PROGRAM=CollectQualityYieldMetrics", "STOP_AFTER=1"],
        "extra_insert_width": ["EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200",
                               "EXTRA_ARGUMENT=CollectInsertSizeMetrics::MINIMUM_PCT=0.4"],
        "extra_gc_window": [f"R={ref}", "PROGRAM=CollectGcBiasMetrics",
                            "EXTRA_ARGUMENT=CollectGcBiasMetrics::SCAN_WINDOW_SIZE=50",
                            "EXTRA_ARGUMENT=CollectGcBiasMetrics::ALSO_IGNORE_DUPLICATES=true"],
        "extra_asm_bisulfite": [f"R={ref}", "EXTRA_ARGUMENT=CollectAlignmentSummaryMetrics::BS=true",
                                "EXTRA_ARGUMENT=CollectAlignmentSummaryMetrics::MAX_INSERT_SIZE=400"],
        "extra_qsd_no_calls": ["EXTRA_ARGUMENT=QualityScoreDistribution::INCLUDE_NO_CALLS=true",
                               "EXTRA_ARGUMENT=MeanQualityByCycle::PF_READS_ONLY=true",
                               "EXTRA_ARGUMENT=CollectBaseDistributionByCycle::ALIGNED_READS_ONLY=true"],
        "extra_quality_yield_oq": ["PROGRAM=CollectQualityYieldMetrics",
                                   "EXTRA_ARGUMENT=CollectQualityYieldMetrics::OQ=false"],
        # INPUT, R, ASSUME_SORTED and STOP_AFTER in an extra argument are dropped.
        "extra_run_level_ignored": ["PROGRAM=null", "PROGRAM=CollectQualityYieldMetrics",
                                    "EXTRA_ARGUMENT=CollectQualityYieldMetrics::STOP_AFTER=3",
                                    "EXTRA_ARGUMENT=CollectQualityYieldMetrics::ASSUME_SORTED=false",
                                    "EXTRA_ARGUMENT=CollectQualityYieldMetrics::INPUT=nothing.bam",
                                    "EXTRA_ARGUMENT=CollectQualityYieldMetrics::R=nothing.fa"],
        "extra_level_one_program": ["PROGRAM=null", "PROGRAM=CollectInsertSizeMetrics",
                                    "EXTRA_ARGUMENT=CollectInsertSizeMetrics::LEVEL=SAMPLE"],
        "no_plots": ["DO_NOT_CREATE_PLOTS=true", f"R={ref}", *ALL_PROGRAMS],
        "rna_options_dropped": ["INTERVALS=/dev/null", "DB_SNP=/dev/null", "REF_FLAT=/dev/null",
                                "IGNORE_SEQUENCE=c1", "INCLUDE_UNPAIRED=true"],
    }


OPTION_INPUTS = ["foundation", "asm_mix", "rl_reads", "t1"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", list(option_cases("x")))
def test_options(built, corpus, option, name):
    bam, ref = inputs(name, built, corpus)
    compare(bam, option_cases(ref)[option], f"{name} {option}")


def test_extra_output_override(built, tmp_path):
    """An OUTPUT in an extra argument moves that program's file."""
    bam, _ = built["rl_reads"]
    for side in ("picard", "bammetrics"):
        args = [TOOL, f"I={bam}", f"O={tmp_path / side}", "PROGRAM=null",
                "PROGRAM=CollectQualityYieldMetrics",
                f"EXTRA_ARGUMENT=CollectQualityYieldMetrics::OUTPUT={tmp_path / (side + '.moved')}"]
        if side == "picard":
            subprocess.run(["picard", *args, SILENT, "QUIET=true", "VERBOSITY=ERROR"], check=True,
                           env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.run([str(BAMMETRICS), *args], check=True, capture_output=True)
    assert sorted(p.name for p in tmp_path.glob("*")) == ["bammetrics.moved", "picard.moved"]
    body = lambda p: p.read_text().split("## METRICS CLASS", 1)[1]
    assert body(tmp_path / "picard.moved") == body(tmp_path / "bammetrics.moved")
    RESULTS.append((f"{TOOL} rl_reads extra OUTPUT", 1, 1, "moved file identical"))


def test_insert_size_file_absent_without_pairs(corpus):
    """No counted pair: neither tool writes the insert size file, the other
    files are compared."""
    written = compare(CORPUS / "t0_no_pairs.bam", [], "t0_no_pairs insert size absent")
    assert ".insert_size_metrics" not in written and len(written) == 4


def test_empty_bam(built):
    bam, ref = built["rl_empty"]
    compare(bam, [f"R={ref}", *ALL_PROGRAMS], "rl_empty all programs",
            [f for f in ALL_FILES if f != ".insert_size_metrics"])


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None):
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"O={tmp_path / 'p'}", SILENT, *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={tmp_path / 'b'}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not list(tmp_path.glob("b*")), "bammetrics wrote a file"
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


FAILURES = {
    "gc_bias_without_reference": (["PROGRAM=CollectGcBiasMetrics"], "REF"),
    "no_program": (["PROGRAM=null"], "No programs"),
    "unknown_program": (["PROGRAM=Bogus"], "Bogus"),
    "rna_seq_not_implemented": (["PROGRAM=RnaSeqMetrics"], "RnaSeqMetrics"),
    "artifacts_not_implemented": (["PROGRAM=CollectSequencingArtifactMetrics"],
                                  "CollectSequencingArtifactMetrics"),
    "extra_for_unrequested": (["PROGRAM=null", "PROGRAM=MeanQualityByCycle",
                               "EXTRA_ARGUMENT=CollectInsertSizeMetrics::HISTOGRAM_WIDTH=200"],
                              "wasn't requested"),
    "extra_unknown_program": (["EXTRA_ARGUMENT=Bogus::X=1"], "Bogus"),
    "extra_malformed": (["EXTRA_ARGUMENT=nothing"], "EXTRA_ARGUMENT"),
    "extra_unknown_option": (["EXTRA_ARGUMENT=CollectInsertSizeMetrics::BOGUS=1"], "BOGUS"),
    "extra_dashed": (["EXTRA_ARGUMENT=CollectInsertSizeMetrics::--HISTOGRAM_WIDTH 200"],
                     "Invalid argument"),
    "bad_level": (["LEVEL=BOGUS"], "BOGUS"),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    bam, _ = built["rl_reads"]
    args, message = FAILURES[case]
    both_fail(bam, args, tmp_path, message)


def test_assume_sorted_false_on_query_name_order_fails(built, tmp_path):
    bam, _ = built["rl_queryname"]
    both_fail(bam, ["ASSUME_SORTED=false"], tmp_path, "coordinate sorted")


@pytest.mark.parametrize("name", ["t0_asm", "t2", "t0_filters"])
def test_a_failing_program_fails_the_run(tmp_path, name):
    """CollectGcBiasMetrics throws on these inputs: no file of any program."""
    both_fail(CORPUS / f"{name}.bam", [f"R={REFERENCE}", "PROGRAM=CollectGcBiasMetrics"], tmp_path)


def test_a_read_without_qualities_fails_the_run(built, tmp_path):
    """CollectAlignmentSummaryMetrics throws on a mapped read with `*` for QUAL."""
    bam, _ = built["rl_star_fields"]
    both_fail(bam, [], tmp_path, "without base qualities")


def test_files_compared_cover_every_program():
    suffixes = {suffix for suffix, _ in SEEN}
    if len(SEEN) < 20:
        pytest.skip("run the whole module")
    assert suffixes >= set(ALL_FILES) | {".quality_yield_metrics.txt", ".insert_size_metrics.tsv"}

"""CollectGcBiasMetrics against Picard 3.5.0: the detail file (OUTPUT) and
the summary file (SUMMARY_OUTPUT), every accumulation level, every option
that changes a number, and error parity.

Inputs: the bwa tiers t0 to t3 with the corpus reference, the foundation
corpus with its own reference, and the asm_bams, read_level_bams and
complexity_bams BAMs, which share one reference. Every case runs Picard with
VALIDATION_STRINGENCY=SILENT.

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared, except the columns in
FEWER_THAN_THREE with the reason.
"""
import subprocess
from collections import defaultdict

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS, FIXTURES, REFERENCE
from compare import compare_outputs
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

TOOL = "CollectGcBiasMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
OUTPUTS = (("O", ".detail"), ("S", ".summary"))
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Every hand-built BAM, each set with its reference."""
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
    """(bam, reference) of a named input."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam", REFERENCE


# Inputs without an aligned read: Picard writes both files with neither class
# block nor row, so the byte comparison of the file bodies is the whole check.
NO_CELLS = {"t0_asm_empty", "asm_single_unmapped", "asm_empty", "rl_empty", "elc_elc_empty",
            "elc_unaligned"}


def compare(bam, ref, args, label):
    cells = compare_outputs(TOOL, bam, [SILENT, f"R={ref}", *args], OUTPUTS, workdir=bam.parent,
                            seen=SEEN)
    assert cells > 0 or label.split()[0] in NO_CELLS or "stop_after_1" in label, (
        f"{TOOL} {label}: no cell compared")
    RESULTS.append((f"{TOOL} {label}", cells, cells, "both files identical"))
    return cells


# --- every input at the defaults -------------------------------------------------

# t0_asm has a read running 40 bases past the end of chrM, t2 more than 100
# secondary records without bases, and the insert-size t0 BAMs their own @SQ
# lengths: both tools fail on them with the reference (below).
TIERS = ["t0_asm_cycles", "t0_asm_empty", "t1", "t3"]
BUILT = ["asm_mix", "asm_unpaired_only", "asm_no_pf_group", "asm_all_qcfail", "asm_single_unmapped",
         "asm_clips_and_flags", "asm_no_read_groups", "asm_some_without_rg", "rl_reads",
         "rl_star_fields", "rl_all_duplicates", "rl_empty", "rl_no_read_groups", "elc_illumina",
         "elc_unaligned", "elc_gc_skew"]


@pytest.mark.parametrize("name", TIERS + ["foundation"] + BUILT)
def test_defaults(built, corpus, name):
    bam, ref = inputs(name, built, corpus)
    compare(bam, ref, [], name)


# --- levels and options ------------------------------------------------------------

LEVELS = {
    "sample": ["LEVEL=SAMPLE"],
    "library": ["LEVEL=LIBRARY"],
    "read_group": ["LEVEL=READ_GROUP"],
    "all_four": ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"],
    "null_then_library": ["LEVEL=null", "LEVEL=LIBRARY"],
}
LEVEL_INPUTS = ["foundation", "asm_mix", "asm_some_without_rg", "rl_reads", "elc_illumina", "t1",
                "elc_gc_skew"]


@pytest.mark.parametrize("name", LEVEL_INPUTS)
@pytest.mark.parametrize("level", LEVELS)
def test_levels(built, corpus, level, name):
    bam, ref = inputs(name, built, corpus)
    compare(bam, ref, LEVELS[level], f"{name} {level}")


OPTIONS = {
    "window_50": ["WINDOW_SIZE=50"],
    "window_150": ["SCAN_WINDOW_SIZE=150"],
    "window_7000": ["WINDOW_SIZE=7000"],  # longer than contig c3, shorter than the others
    "window_70000": ["WINDOW_SIZE=70000"],  # longer than every contig: no window at all
    "unique": ["ALSO_IGNORE_DUPLICATES=true"],
    "bisulfite": ["BS=true"],
    "bisulfite_unique_levels": ["IS_BISULFITE_SEQUENCED=true", "ALSO_IGNORE_DUPLICATES=true",
                                "LEVEL=LIBRARY", "LEVEL=READ_GROUP"],
    "mgf": ["MGF=0.5"],  # accepted by both, read by neither
    "mgf_long": ["MINIMUM_GENOME_FRACTION=0"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_250": ["STOP_AFTER=250"],
    "chart": ["CHART=/dev/null"],
    "assume_sorted_true": ["ASSUME_SORTED=true"],
}
OPTION_INPUTS = ["foundation", "asm_mix", "rl_reads", "t1", "elc_illumina", "elc_gc_skew"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, corpus, option, name):
    bam, ref = inputs(name, built, corpus)
    compare(bam, ref, OPTIONS[option], f"{name} {option}")


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None):
    picard = [f"O={tmp_path / 'p.detail'}", f"S={tmp_path / 'p.summary'}"]
    ours = [f"O={tmp_path / 'b.detail'}", f"S={tmp_path / 'b.summary'}"]
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", *picard, SILENT, *args, "QUIET=true",
                        "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", *ours, *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr
    assert not list(tmp_path.glob("*.detail")) and not list(tmp_path.glob("*.summary"))
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


def test_without_reference_fails(built, tmp_path):
    bam, _ = built["rl_reads"]
    both_fail(bam, [], tmp_path, "REFERENCE_SEQUENCE")


def test_assume_sorted_false_on_query_name_order_fails(built, tmp_path):
    bam, ref = built["rl_queryname"]
    both_fail(bam, [f"R={ref}", "ASSUME_SORTED=false"], tmp_path, "coordinate sorted")


def test_without_summary_output_fails(built, tmp_path):
    bam, ref = built["rl_reads"]
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"O={tmp_path / 'p.detail'}", f"R={ref}", SILENT,
                        "QUIET=true", "VERBOSITY=ERROR"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={tmp_path / 'b.detail'}",
                             f"R={ref}"], capture_output=True, text=True)
    assert result.returncode != 0 and "SUMMARY_OUTPUT" in result.stderr


def test_read_past_the_contig_end_fails(tmp_path):
    """Picard throws an index error counting mismatches of a read that runs
    past the end of its contig."""
    both_fail(CORPUS / "t0_asm.bam", [f"R={REFERENCE}"], tmp_path, "past the end of the contig")


def test_reference_contigs_must_match_the_header(built, corpus, tmp_path):
    bam, _ = built["rl_reads"]
    both_fail(bam, [f"R={corpus / 'ref.fa'}"], tmp_path)
    both_fail(CORPUS / "t0_filters.bam", [f"R={REFERENCE}"], tmp_path)


def test_more_than_100_records_without_bases_fail(tmp_path):
    """t2 has over 100 secondary records with `*` for SEQ: each unit drops
    its first 100 and then Picard counts mismatches on an empty read."""
    both_fail(CORPUS / "t2.bam", [f"R={REFERENCE}"], tmp_path, "CIGAR longer than the sequence")


@pytest.mark.parametrize("window", [0, -5])
def test_a_window_under_one_base_fails(built, tmp_path, window):
    """Picard divides by zero at 0 and indexes past the contig below it.
    bammetrics panicked on both before it refused them."""
    bam, ref = built["rl_reads"]
    both_fail(bam, [f"R={ref}", f"SCAN_WINDOW_SIZE={window}"], tmp_path,
              f"SCAN_WINDOW_SIZE={window}: a window holds at least 1 base")


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {
    ("summary", "READS_USED"): "ALL or UNIQUE by definition",
    ("detail", "READS_USED"): "ALL or UNIQUE by definition",
}


@pytest.mark.parametrize("which", ["detail", "summary"])
def test_every_column_takes_three_non_default_values(which):
    files = [parsed for suffix, parsed in SEEN if suffix == f".{which}"]
    if len(files) < 10:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
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
             if len(values[c]) < 3 and (which, c) not in FEWER_THAN_THREE}
    assert not short, f"{TOOL} {which}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} {which} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))


def test_fixture_runs_without_read_groups(tmp_path):
    """The committed fixture has no @RG and no reference of its own: Picard
    needs one, so this only checks bammetrics rejects a missing reference."""
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={FIXTURES / 'insert_size.bam'}",
                             f"O={tmp_path / 'd'}", f"S={tmp_path / 's'}"], capture_output=True, text=True)
    assert result.returncode != 0 and "REFERENCE_SEQUENCE" in result.stderr

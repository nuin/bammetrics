"""CollectAlignmentSummaryMetrics against Picard 3.5.0 at every accumulation level.

Each case runs `picard CollectAlignmentSummaryMetrics` and `bammetrics
CollectAlignmentSummaryMetrics` with the same KEY=VALUE arguments and
compares every metrics cell and histogram cell as printed (compare_metrics).
Inputs: the bwa tiers t0 to t3 (tests/corpus/out), the 12-read-group
foundation corpus, and the hand-built BAMs of asm_bams.py.

The last test checks that the cases exercise each column: every column must
take at least 3 distinct non-default values over all Picard files compared,
except the columns listed in FEWER_THAN_THREE with the reason.

Needs picard, samtools and pysam; run `make concordance`.
"""
import subprocess
from collections import defaultdict

import pytest

from asm_bams import INDEXED_3P, PAIRED_3P, build_all
from cases import BAMMETRICS, CORPUS, REFERENCE
from compare import compare_metrics, run_picard_tool
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CollectAlignmentSummaryMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
ALL_LEVELS = ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"]
SEEN = []


def compare(bam, args, label):
    cells = compare_metrics(TOOL, bam, [SILENT, *args], workdir=bam.parent, expect_file=True,
                            seen=SEEN)
    assert cells > 0
    RESULTS.append((f"alignment_summary {label}", cells, cells, "every cell the same printed string"))
    return cells


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("asm")
    return out, build_all(out)


# --- tiers t0 to t3 and the foundation corpus, all four levels -----------------

# t0 BAMs built for CollectInsertSizeMetrics carry their own @SQ lengths, so
# they run without a reference.
TIERS = {
    "t0_asm": True, "t0_asm_cycles": True, "t0_asm_empty": True, "t0_filters": False,
    "t0_no_pairs": False, "t0_all_tandem": False, "t1": True, "t2": True, "t3": True,
}
TIER_CASES = [(name, ref) for name, has_ref in TIERS.items() for ref in (False, True)
              if has_ref or not ref]


@pytest.mark.parametrize("name,ref", TIER_CASES,
                         ids=[f"{n}{'_ref' if r else ''}" for n, r in TIER_CASES])
def test_tiers_at_all_levels(name, ref):
    args = ALL_LEVELS + ([f"R={REFERENCE}"] if ref else [])
    compare(CORPUS / f"{name}.bam", args, f"{name}{' R' if ref else ''} all levels")


@pytest.mark.parametrize("ref", [False, True], ids=["no_ref", "ref"])
def test_foundation_at_all_levels(corpus, ref):
    args = ALL_LEVELS + ([f"R={corpus / 'ref.fa'}"] if ref else [])
    compare(corpus / "corpus.bam", args, f"foundation{' R' if ref else ''} all levels")


LEVEL_SETS = {
    "sample": ["LEVEL=SAMPLE"],
    "library": ["METRIC_ACCUMULATION_LEVEL=LIBRARY"],
    "read_group": ["LEVEL=READ_GROUP"],
    "only_read_group": ["LEVEL=null", "LEVEL=READ_GROUP"],
    "only_sample_library": ["LEVEL=null", "LEVEL=LIBRARY", "LEVEL=SAMPLE"],
}


@pytest.mark.parametrize("levels", LEVEL_SETS)
def test_level_sets(built, corpus, levels):
    """Single levels, and lists without ALL_READS (no histogram is written)."""
    out, bams = built
    compare(corpus / "corpus.bam", LEVEL_SETS[levels], f"foundation {levels}")
    compare(bams["mix"], LEVEL_SETS[levels] + [f"R={out / 'ref.fa'}"], f"mix {levels}")


# --- every option that changes a number, at all four levels --------------------

OPTIONS = {
    "bisulfite": ["IS_BISULFITE_SEQUENCED=true"],
    "bisulfite_short_name": ["BS=true"],
    "adapters_cleared": ["ADAPTER_SEQUENCE=null"],
    "adapter_added": [f"ADAPTER_SEQUENCE={INDEXED_3P[:20]}TTTTTTTTTT"],
    "adapter_paired_only": ["ADAPTER_SEQUENCE=null", f"ADAPTER_SEQUENCE={PAIRED_3P}"],
    "max_insert_300": ["MAX_INSERT_SIZE=300"],
    "max_insert_0": ["MAX_INSERT_SIZE=0"],
    "orient_add_rf": ["EXPECTED_PAIR_ORIENTATIONS=RF"],
    "orient_tandem_only": ["EXPECTED_PAIR_ORIENTATIONS=null", "EXPECTED_PAIR_ORIENTATIONS=TANDEM"],
    "orient_all": ["EXPECTED_PAIR_ORIENTATIONS=RF", "EXPECTED_PAIR_ORIENTATIONS=TANDEM"],
    "stop_after_1": ["STOP_AFTER=1"],
    "stop_after_700": ["STOP_AFTER=700"],
    "stop_after_0": ["STOP_AFTER=0"],
    "no_alignment_info": ["COLLECT_ALIGNMENT_INFORMATION=false"],
}


@pytest.mark.parametrize("option", OPTIONS)
@pytest.mark.parametrize("ref", [False, True], ids=["no_ref", "ref"])
def test_options_on_mix(built, option, ref):
    out, bams = built
    args = ALL_LEVELS + OPTIONS[option] + ([f"R={out / 'ref.fa'}"] if ref else [])
    compare(bams["mix"], args, f"mix {option}{' R' if ref else ''}")


@pytest.mark.parametrize("option", ["bisulfite", "max_insert_300", "orient_all", "stop_after_700",
                                    "no_alignment_info", "adapters_cleared"])
def test_options_on_t2_and_foundation(corpus, option):
    compare(CORPUS / "t2.bam", ALL_LEVELS + OPTIONS[option] + [f"R={REFERENCE}"], f"t2 {option}")
    compare(corpus / "corpus.bam", ALL_LEVELS + OPTIONS[option] + [f"R={corpus / 'ref.fa'}"],
            f"foundation {option}")


# --- hand-built edge cases -----------------------------------------------------

EDGE = ["unpaired_only", "single_unmapped", "clips_and_flags", "no_read_groups",
        "some_without_rg", "empty"]


@pytest.mark.parametrize("variant", ["no_ref", "ref", "ref_bisulfite"])
@pytest.mark.parametrize("name", EDGE)
def test_edge_cases(built, name, variant):
    out, bams = built
    extra = {"no_ref": [], "ref": [f"R={out / 'ref.fa'}"],
             "ref_bisulfite": [f"R={out / 'ref.fa'}", "IS_BISULFITE_SEQUENCED=true"]}[variant]
    compare(bams[name], ALL_LEVELS + extra, f"{name} {variant}")


def test_unit_without_pf_reads_at_all_reads(built):
    """Group b fails QC throughout; at ALL_READS its reads are a minority."""
    _, bams = built
    compare(bams["no_pf_group"], [], "no_pf_group ALL_READS")


@pytest.mark.parametrize("name,levels", [
    ("no_pf_group", ["LEVEL=READ_GROUP"]),
    ("no_pf_group", ["LEVEL=null", "LEVEL=SAMPLE"]),
    ("all_qcfail", []),
])
def test_no_pf_reads_fails_like_picard(built, tmp_path, name, levels):
    """A unit whose reads all fail QC: Picard throws and writes no file."""
    _, bams = built
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    with pytest.raises(subprocess.CalledProcessError):
        run_picard_tool(TOOL, bams[name], picard_out, [SILENT, *levels])
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={bams[name]}", f"O={ours_out}", *levels],
                          capture_output=True, text=True)
    assert ours.returncode != 0 and "no PF_READS" in ours.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"alignment_summary {name} {' '.join(levels)}", 0, 0,
                    "both exit non-zero, no file"))


# --- every column exercised ----------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values in any input, and why.
FEWER_THAN_THREE = {}


def test_every_column_takes_three_non_default_values():
    if len(SEEN) < 20:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    # Metric cells by the deepest label a row carries: a row with READ_GROUP
    # set is a READ_GROUP row, and so on. An empty header bucket missing that
    # field counts one level up.
    cells_by_level = defaultdict(int)
    for parsed in SEEN:
        for row in parsed.rows:
            level = next((name for name in ("READ_GROUP", "LIBRARY", "SAMPLE") if row[name]),
                         "ALL_READS")
            cells_by_level[level] += len(row)
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
    columns = SEEN[0].columns
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3 and c not in FEWER_THAN_THREE}
    assert not short, f"columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append(("alignment_summary distinct values", len(columns), len(columns),
                    f"columns with at least 3 non-default values; fewest {fewest} "
                    f"{len(values[fewest])}"))
    for level, cells in sorted(cells_by_level.items()):
        RESULTS.append((f"alignment_summary metric cells at {level}", cells, cells, "compared"))


def test_picard_is_on_path():
    """The module needs the picard wrapper; fail loudly instead of skipping."""
    out = subprocess.run(["picard", TOOL, "--version"], capture_output=True, text=True,
                         env=picard_env())
    assert "3.5.0" in out.stdout + out.stderr

"""CollectVariantCallingMetrics, AccumulateVariantCallingMetrics and GenotypeConcordance
against Picard 3.5.0 on the variant_vcfs.py corpus: every file of every tool, every option
that changes a number at a non-default value, THREAD_COUNT invariance, error parity, the
scheme table coverage and the distinct-value gate per tool.

Every Picard file is compared cell by cell and byte for byte after the header lines
(compare_prefixed). The last tests check that every (truth state, call state) cell the
states can reach was counted by some site, and that every column of each tool took at
least 3 distinct non-default values, except the columns in FEWER_THAN_THREE with the
reason.
"""
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import pytest

from cases import BAMMETRICS
from compare import compare_prefixed, compare_written
from conftest import RESULTS
from run_picard import picard_env
from variant_vcfs import MISSING_HOM_ONLY, REACHABLE, build_all, header, record

CVCM = "CollectVariantCallingMetrics"
ACCUMULATE = "AccumulateVariantCallingMetrics"
GC = "GenotypeConcordance"
SEEN = {CVCM: [], ACCUMULATE: [], GC: []}


@pytest.fixture(scope="module")
def vcfs(tmp_path_factory):
    """{name: path} of the variant corpus, plus `dir`."""
    out = tmp_path_factory.mktemp("variants")
    files = build_all(out)
    files["dir"] = out
    return files


def compare(tool, args, label, bam_args=None, workdir=None):
    cells = compare_prefixed(tool, args, bam_args=bam_args, workdir=workdir, seen=SEEN[tool])
    assert cells > 0, f"{tool} {label}: no cell compared"
    RESULTS.append((f"{tool} {label}", cells, cells, "every file identical"))
    return cells


# --- CollectVariantCallingMetrics ----------------------------------------------------

CVCM_CASES = {
    "default": ("calls3.vcf.gz", "dbsnp.vcf.gz", []),
    "dbsnp_plain": ("calls3.vcf.gz", "dbsnp.vcf", []),
    "targets": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["TI=targets.interval_list"]),
    "wide": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["TARGET_INTERVALS=wide.interval_list"]),
    "left": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["TI=left.interval_list"]),
    "right": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["TI=right.interval_list"]),
    "sd_dict": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["SD=ref.dict"]),
    "sd_vcf": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["SEQUENCE_DICTIONARY=calls3.vcf.gz"]),
    "sd_interval_list": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["SD=left.interval_list", "TI=targets.interval_list"]),
    "threads2": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["THREAD_COUNT=2"]),
    "threads3_targets": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["THREAD_COUNT=3", "TI=targets.interval_list"]),
    "threads8": ("calls3.vcf.gz", "dbsnp.vcf.gz", ["THREAD_COUNT=8", "GVCF_INPUT=false"]),
    "gvcf": ("gvcf.vcf.gz", "dbsnp.vcf.gz", ["GVCF_INPUT=true"]),
    "gvcf_as_vcf": ("gvcf.vcf.gz", "dbsnp.vcf.gz", ["GVCF_INPUT=false"]),
    "gvcf_targets": ("gvcf.vcf.gz", "dbsnp.vcf.gz", ["GVCF_INPUT=true", "TI=targets.interval_list"]),
    "sites_only": ("sites_only.vcf.gz", "dbsnp.vcf.gz", []),
    "dbsnp_is_truth": ("calls3.vcf.gz", "truth.vcf.gz", []),
    "self_dbsnp": ("calls3.vcf.gz", "calls3.vcf.gz", ["TI=wide.interval_list"]),
}


def cvcm_args(vcfs, case):
    vcf, dbsnp, extra = CVCM_CASES[case]
    args = [f"I={vcfs[vcf]}", f"DBSNP={vcfs[dbsnp]}"]
    for arg in extra:
        key, value = arg.split("=", 1)
        args.append(f"{key}={vcfs[value]}" if value in vcfs else arg)
    return args


@pytest.mark.parametrize("case", sorted(CVCM_CASES))
def test_variant_calling_metrics(vcfs, case):
    compare(CVCM, cvcm_args(vcfs, case), case, workdir=vcfs["dir"])


def test_variant_calling_metrics_reads_an_unindexed_input(vcfs, tmp_path):
    """Picard queries INPUT through its index and fails without one; bammetrics reads the
    plain file and writes the files of the indexed run (KNOWN_DIFFERENCES.md)."""
    plain = [f"I={vcfs['calls3.vcf']}", f"DBSNP={vcfs['dbsnp.vcf.gz']}"]
    compare(CVCM, cvcm_args(vcfs, "default"), "plain input", bam_args=plain, workdir=vcfs["dir"])
    picard = subprocess.run(["picard", CVCM, *plain, f"O={tmp_path / 'p'}", "QUIET=true", "VERBOSITY=ERROR"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0
    assert "An index is required" in picard.stdout + picard.stderr
    assert not list(tmp_path.glob("p.*"))


# --- AccumulateVariantCallingMetrics -------------------------------------------------


@pytest.fixture(scope="module")
def shards(vcfs):
    """Picard's and bammetrics' CollectVariantCallingMetrics files for the left, right
    and whole-genome runs: {name: (picard prefix, bammetrics prefix)}."""
    out = vcfs["dir"] / "shards"
    out.mkdir()
    prefixes = {}
    for name in ("left", "right", "default"):
        args = cvcm_args(vcfs, name)
        picard, ours = out / f"picard_{name}", out / f"bammetrics_{name}"
        subprocess.run(["picard", CVCM, *args, f"O={picard}", "QUIET=true", "VERBOSITY=ERROR"], check=True,
                       env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(BAMMETRICS), CVCM, *args, f"O={ours}"], check=True, capture_output=True)
        prefixes[name] = (picard, ours)
    return prefixes


ACCUMULATE_CASES = {
    "two_shards": ["left", "right"],
    "three": ["left", "right", "default"],
    "single": ["default"],
    "twice": ["left", "left"],
    "right_only": ["right", "right", "right"],
    "left_thrice": ["left", "left", "left"],
    "skewed": ["left", "right", "right"],
}


@pytest.mark.parametrize("case", sorted(ACCUMULATE_CASES))
def test_accumulate_variant_calling_metrics(vcfs, shards, case):
    """Each tool merges its own shard files; the printed ratios are inverted back into
    counts, so `three` and `twice` double-count as Picard does."""
    names = ACCUMULATE_CASES[case]
    picard_args = [f"I={shards[n][0]}" for n in names]
    bam_args = [f"I={shards[n][1]}" for n in names]
    compare(ACCUMULATE, picard_args, case, bam_args=bam_args, workdir=vcfs["dir"])


def test_accumulate_cross_merges_picard_shards(vcfs, shards):
    """bammetrics on Picard's shard files gives Picard's merged files."""
    picard_args = [f"I={shards[n][0]}" for n in ("left", "right")]
    compare(ACCUMULATE, picard_args, "picard shards both", bam_args=picard_args, workdir=vcfs["dir"])


# --- GenotypeConcordance ---------------------------------------------------------------

GC_BASE = ["TV=truth.vcf", "CV=call.vcf", "CS=C"]
GC_INDEXED = ["TRUTH_VCF=truth.vcf.gz", "CALL_VCF=call.vcf.gz", "CALL_SAMPLE=C"]
GC_CASES = {
    "default": GC_BASE,
    "indexed": [*GC_INDEXED, "USE_VCF_INDEX=true"],
    "index_present_unused": [*GC_INDEXED],
    "intervals": [*GC_INDEXED, "INTERVALS=targets.interval_list"],
    "intervals_wide": [*GC_INDEXED, "INTERVALS=wide.interval_list"],
    "intervals_left": [*GC_INDEXED, "INTERVALS=left.interval_list"],
    "intersect": [*GC_INDEXED, "INTERVALS=targets.interval_list", "INTERVALS=wide.interval_list"],
    "intersect_explicit": [*GC_INDEXED, "INTERVALS=targets.interval_list", "INTERVALS=wide.interval_list",
                           "INTERSECT_INTERVALS=true"],
    "union": [*GC_INDEXED, "INTERVALS=targets.interval_list", "INTERVALS=wide.interval_list",
              "INTERSECT_INTERVALS=false"],
    "min_gq": [*GC_BASE, "MIN_GQ=20"],
    "min_gq_high": [*GC_BASE, "MIN_GQ=45"],
    "min_dp": [*GC_BASE, "MIN_DP=10"],
    "min_both": [*GC_BASE, "MIN_GQ=20", "MIN_DP=10"],
    "min_dp_high": [*GC_BASE, "MIN_DP=25"],
    "ignore_filter": [*GC_BASE, "IGNORE_FILTER_STATUS=true"],
    "ignore_filter_min": [*GC_BASE, "IGNORE_FILTER_STATUS=true", "MIN_GQ=20"],
    "all_rows": [*GC_BASE, "OUTPUT_ALL_ROWS=true"],
    "all_rows_min": [*GC_BASE, "OUTPUT_ALL_ROWS=true", "MIN_DP=10"],
    "missing_hom": [*GC_INDEXED, "INTERVALS=targets.interval_list", "MISSING_SITES_HOM_REF=true"],
    "missing_hom_union": [*GC_INDEXED, "INTERVALS=targets.interval_list", "INTERVALS=wide.interval_list",
                          "INTERSECT_INTERVALS=false", "MISSING_HOM=true"],
    "missing_hom_all_rows": [*GC_INDEXED, "INTERVALS=left.interval_list", "MISSING_HOM=true",
                             "OUTPUT_ALL_ROWS=true", "MIN_GQ=20"],
    "call_c2": ["TV=truth.vcf", "CV=call.vcf", "CS=C2"],
    "truth_explicit": ["TV=truth.vcf", "TS=T", "CV=call.vcf", "CS=C"],
    "swapped": ["TV=call.vcf", "TS=C", "CV=truth.vcf"],
    "swapped_c2": ["TV=call.vcf", "TRUTH_SAMPLE=C2", "CV=truth.vcf", "CS=T"],
    "swapped_min": ["TV=call.vcf", "TS=C", "CV=truth.vcf", "MIN_GQ=20", "MIN_DP=10"],
}


def gc_args(vcfs, case):
    args = []
    for arg in GC_CASES[case]:
        key, value = arg.split("=", 1)
        args.append(f"{key}={vcfs[value]}" if value in vcfs else arg)
    return args


@pytest.mark.parametrize("case", sorted(GC_CASES))
def test_genotype_concordance(vcfs, case):
    compare(GC, gc_args(vcfs, case), case, workdir=vcfs["dir"])


# --- error parity -----------------------------------------------------------------------


def both_fail(tool, args, tmp_path, message, picard_message):
    """Both tools exit non-zero with `message` (bammetrics) and `picard_message`
    (Picard) in their output, and neither writes a file under OUTPUT."""
    picard = subprocess.run(["picard", tool, *args, f"O={tmp_path / 'picard'}", "QUIET=true",
                             "VERBOSITY=ERROR"], env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr
    result = subprocess.run([str(BAMMETRICS), tool, *args, f"O={tmp_path / 'bammetrics'}"],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    assert message in result.stderr, result.stderr
    assert not list(tmp_path.glob("picard*")) and not list(tmp_path.glob("bammetrics*"))
    RESULTS.append((f"{tool} {' '.join(a.split('/')[-1] for a in args)[:60]}", 0, 0,
                    "both exit non-zero, no file"))


def test_dbsnp_is_required(vcfs, tmp_path):
    both_fail(CVCM, [f"I={vcfs['calls3.vcf.gz']}"], tmp_path, "Option 'DBSNP' is required.",
              "Option 'DBSNP' is required.")


def test_output_is_required(vcfs, tmp_path):
    for tool, args in ((CVCM, [f"I={vcfs['calls3.vcf.gz']}", f"DBSNP={vcfs['dbsnp.vcf.gz']}"]),
                       (GC, [f"TV={vcfs['truth.vcf']}", f"CV={vcfs['call.vcf']}", "CS=C"]),
                       (ACCUMULATE, [f"I={vcfs['dir'] / 'x'}"])):
        picard = subprocess.run(["picard", tool, *args], env=picard_env(), capture_output=True, text=True)
        assert "Option 'OUTPUT' is required." in picard.stdout + picard.stderr
        ours = subprocess.run([str(BAMMETRICS), tool, *args], capture_output=True, text=True)
        assert ours.returncode != 0 and "Option 'OUTPUT' is required." in ours.stderr


def test_dbsnp_dictionary_must_match(vcfs, tmp_path):
    other = tmp_path / "other_dbsnp.vcf"
    lines = [line for line in header([]) if not line.startswith("##contig")]
    lines.insert(1, "##contig=<ID=chr1,length=1000>")
    other.write_text("\n".join(lines + [record("chr1", 10, "A", ["G"], [])]) + "\n")
    both_fail(CVCM, [f"I={vcfs['calls3.vcf.gz']}", f"DBSNP={other}"], tmp_path,
              "does not match sequence dictionary for (INPUT)", "does not match sequence dictionary for (INPUT)")


def test_input_without_contig_lines(vcfs, tmp_path):
    """Picard compares the dbSNP dictionary with the input's null one."""
    both_fail(CVCM, [f"I={vcfs['calls3_nocontigs.vcf.gz']}", f"DBSNP={vcfs['dbsnp.vcf.gz']}"], tmp_path,
              "no ##contig line", "NullPointerException")


def test_gvcf_input_needs_one_sample(vcfs, tmp_path):
    both_fail(CVCM, [f"I={vcfs['calls3.vcf.gz']}", f"DBSNP={vcfs['dbsnp.vcf.gz']}", "GVCF_INPUT=true"], tmp_path,
              "Expected to have exactly 1 sample in a GVCF, found 3",
              "Expected to have exactly 1 sample in a GVCF, found 3")


def test_thread_count_must_be_positive(vcfs, tmp_path):
    both_fail(CVCM, [f"I={vcfs['calls3.vcf.gz']}", f"DBSNP={vcfs['dbsnp.vcf.gz']}", "THREAD_COUNT=0"], tmp_path,
              "Multithreading value must exceed 0.", "Multithreading value must exceed 0.")


def test_target_intervals_need_an_indexed_dbsnp(vcfs, tmp_path):
    """Picard queries the dbSNP file by interval; bammetrics reads it whole, so this is
    a Picard-only failure (KNOWN_DIFFERENCES.md): bammetrics writes the indexed run's files."""
    args = [f"I={vcfs['calls3.vcf.gz']}", f"DBSNP={vcfs['dbsnp.vcf']}", f"TI={vcfs['targets.interval_list']}"]
    picard = subprocess.run(["picard", CVCM, *args, f"O={tmp_path / 'p'}", "QUIET=true", "VERBOSITY=ERROR"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0 and "An index is required" in picard.stdout + picard.stderr
    compare(CVCM, cvcm_args(vcfs, "targets"), "plain dbSNP with targets", bam_args=args, workdir=vcfs["dir"])


def test_accumulate_needs_every_shard(vcfs, shards, tmp_path):
    both_fail(ACCUMULATE, [f"I={shards['left'][0]}", f"I={tmp_path / 'none'}"], tmp_path, "none",
              "Cannot read non-existent file")


def test_accumulate_summary_needs_one_row(vcfs, shards, tmp_path):
    prefix = tmp_path / "bad"
    shutil.copy(str(shards["left"][0]) + ".variant_calling_detail_metrics", f"{prefix}.variant_calling_detail_metrics")
    summary = Path(str(shards["left"][0]) + ".variant_calling_summary_metrics").read_text()
    row = summary.splitlines()[-3]
    (Path(f"{prefix}.variant_calling_summary_metrics")).write_text(summary.replace(row + "\n", row + "\n" + row + "\n"))
    both_fail(ACCUMULATE, [f"I={prefix}"], tmp_path, "Expected 1 row in the summary metrics file but saw 2",
              "Expected 1 row in the summary metrics file but saw 2")


def test_missing_hom_needs_intervals(vcfs, tmp_path):
    both_fail(GC, [*gc_args(vcfs, "default"), "MISSING_SITES_HOM_REF=true"], tmp_path,
              "You cannot use the MISSING_HOM option", "You cannot use the MISSING_HOM option")


def test_use_vcf_index_needs_the_index(vcfs, tmp_path):
    both_fail(GC, [*gc_args(vcfs, "default"), "USE_VCF_INDEX=true"], tmp_path,
              "The index file was not found for the TRUTH VCF", "The index file was not found for the TRUTH VCF")


def test_intervals_need_the_index(vcfs, tmp_path):
    both_fail(GC, [*gc_args(vcfs, "default"), f"INTERVALS={vcfs['targets.interval_list']}"], tmp_path,
              "The index file was not found for the CALL VCF", "The index file was not found for the CALL VCF")


def test_truth_sample_is_required_with_two_samples(vcfs, tmp_path):
    both_fail(GC, [f"TV={vcfs['call.vcf']}", f"CV={vcfs['truth.vcf']}"], tmp_path,
              "TRUTH_SAMPLE is required when the TRUTH_VCF has more than one sample",
              "TRUTH_SAMPLE is required when the TRUTH_VCF has more than one sample")


def test_call_sample_is_required_with_two_samples(vcfs, tmp_path):
    both_fail(GC, [f"TV={vcfs['truth.vcf']}", f"CV={vcfs['call.vcf']}"], tmp_path,
              "CALL_SAMPLE is required when the CALL_VCF has more than one sample",
              "CALL_SAMPLE is required when the CALL_VCF has more than one sample")


def test_unknown_sample(vcfs, tmp_path):
    both_fail(GC, [f"TV={vcfs['truth.vcf']}", f"CV={vcfs['call.vcf']}", "CS=X"], tmp_path,
              "does not contain genotypes for sample X", "does not contain genotypes for sample X")


def test_intervals_of_another_dictionary(vcfs, tmp_path):
    both_fail(GC, [*gc_args(vcfs, "indexed"), f"INTERVALS={vcfs['other.interval_list']}"], tmp_path,
              "Sequence dictionary for (provided intervals) does not match",
              "Sequence dictionary for (provided intervals) does not match")


def test_truth_without_contig_lines(vcfs, tmp_path):
    """Picard builds its position comparator from the truth header's dictionary."""
    both_fail(GC, [f"TV={vcfs['nocontigs.vcf']}", f"CV={vcfs['call.vcf']}", "CS=C"], tmp_path,
              "NullPointerException", "NullPointerException")


def test_truth_and_call_without_contig_lines(vcfs, tmp_path):
    call = tmp_path / "call_nocontigs.vcf"
    call.write_text("".join(line for line in vcfs["call.vcf"].open() if not line.startswith("##contig")))
    both_fail(GC, [f"TV={vcfs['nocontigs.vcf']}", f"CV={call}", "CS=C"], tmp_path, "no ##contig line",
              "NullPointerException")


def test_reference_alleles_must_share_a_prefix(vcfs, tmp_path):
    both_fail(GC, [f"TV={vcfs['mismatch_truth.vcf']}", f"CV={vcfs['mismatch_call.vcf']}"], tmp_path,
              "Ref alleles mismatch", "Ref alleles mismatch")


def test_haploid_genotype(vcfs, tmp_path):
    both_fail(GC, [f"TV={vcfs['mismatch_truth.vcf']}", f"CV={vcfs['haploid_call.vcf']}"], tmp_path,
              "does not have exactly 2 alleles", "does not have exactly 2 alleles")


def test_input_option_is_rejected(vcfs, tmp_path):
    both_fail(GC, [*gc_args(vcfs, "default"), f"I={vcfs['truth.vcf']}"], tmp_path, "Unrecognized option: I",
              "Unrecognized option: I")


# --- the scheme table ---------------------------------------------------------------------


def narrowed(config):
    """True when -k, --deselect or a node id selected part of the module."""
    return bool(config.getoption("keyword") or config.getoption("deselect")
                or any("::" in arg for arg in config.args))


def detail_cells():
    """{(truth state, call state): count} summed over every GenotypeConcordance detail file."""
    cells = defaultdict(int)
    for suffix, parsed in SEEN[GC]:
        if suffix.endswith("detail_metrics"):
            for row in parsed.rows:
                cells[(row["TRUTH_STATE"], row["CALL_STATE"])] += int(row["COUNT"])
    return cells


def test_every_reachable_scheme_cell_is_hit(request):
    if len(SEEN[GC]) != 3 * len(GC_CASES):
        message = f"{GC}: {len(SEEN[GC])} of {3 * len(GC_CASES)} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(message)
    hit = {cell for cell, count in detail_cells().items() if count > 0}
    expected = REACHABLE | MISSING_HOM_ONLY
    assert expected - hit == set(), f"cells no site reached: {sorted(expected - hit)}"
    assert hit - expected == set(), f"cells reached that variant_vcfs.py lists as unreachable: {sorted(hit - expected)}"
    RESULTS.append((f"{GC} scheme cells", len(hit), len(expected), "every reachable cell counted"))


# --- the distinct-value gate ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, per tool, and why.
FEWER_THAN_THREE = {
    CVCM: {},
    ACCUMULATE: {},
    GC: {
        "VARIANT_TYPE": "the row key: SNP or INDEL",
        "HET_SPECIFICITY": "Picard sets it to NaN (prints ?): no TN for a het truth",
        "HOMVAR_SPECIFICITY": "Picard sets it to NaN (prints ?): no TN for a hom-var truth",
    },
}
EXPECTED_FILES = {
    CVCM: 2 * (len(CVCM_CASES) + 2),
    ACCUMULATE: 2 * (len(ACCUMULATE_CASES) + 1),
    GC: 3 * len(GC_CASES),
}


@pytest.mark.parametrize("tool", [CVCM, ACCUMULATE, GC])
def test_every_column_takes_three_non_default_values(request, tool):
    if len(SEEN[tool]) != EXPECTED_FILES[tool]:
        message = f"{tool}: {len(SEEN[tool])} of {EXPECTED_FILES[tool]} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(f"{message}: a comparison above failed or EXPECTED_FILES is stale")
    values = defaultdict(set)
    columns = defaultdict(list)
    for suffix, parsed in SEEN[tool]:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[suffix, column].add(value)
        for column in parsed.columns:
            if column not in columns[suffix]:
                columns[suffix].append(column)
    short = {f"{suffix} {c}": sorted(values[suffix, c]) for suffix in columns for c in columns[suffix]
             if len(values[suffix, c]) < 3 and c not in FEWER_THAN_THREE[tool]}
    assert not short, f"{tool}: columns with fewer than 3 distinct non-default values: {short}"
    checked = sum(len(c) for c in columns.values())
    fewest = min(((suffix, c) for suffix in columns for c in columns[suffix]), key=lambda k: len(values[k]))
    RESULTS.append((f"{tool} distinct values", checked, checked,
                    f"columns checked; fewest {fewest[1]} {len(values[fewest])}"))

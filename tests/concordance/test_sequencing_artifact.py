"""CollectSequencingArtifactMetrics and ConvertSequencingArtifactToOxoG against
Picard 3.5.0: all five files and the OxoG file cell by cell, every option that
changes a number, CONTEXT_SIZE 0 to 2, INTERVALS and DB_SNP, and error parity.

Inputs: the foundation corpus (two samples, four libraries), corpus tiers t1 to
t3, the exome slice, and three BAMs artifact_bams.py builds from the foundation
corpus (decorated: QC-failed, OQ-tagged, MAPQ 255 and unpaired records, a read
group without LB, a lower-case and IUPAC reference; no_read_groups;
unplaced_middle). DB_SNP is a sites-only VCF artifact_bams.write_dbsnp writes on
each reference, plain and bgzip with a tabix index.

Each case runs both collectors with O=<dir>/<side>, then each side's converter
on its own five files, and compares every file. Where Picard's converter fails
(CONTEXTS_TO_PRINT leaves out a context it needs), bammetrics' must fail too.

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared, except the columns in FEWER_THAN_THREE
with the reason.
"""
import subprocess
from collections import defaultdict
from pathlib import Path

import pytest

from artifact_bams import build_all, write_dbsnp
from cases import BAMMETRICS, CORPUS
from compare import compare_written
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CollectSequencingArtifactMetrics"
CONVERT = "ConvertSequencingArtifactToOxoG"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []
DETAIL_FILES = (".pre_adapter_detail_metrics", ".bait_bias_detail_metrics")


@pytest.fixture(scope="module")
def built(tmp_path_factory, corpus):
    out = tmp_path_factory.mktemp("artifact")
    return build_all(out, corpus)


@pytest.fixture(scope="module")
def dbsnp(tmp_path_factory, corpus, exome):
    """{reference: (plain VCF, indexed bgzip VCF)} for every reference used."""
    out = tmp_path_factory.mktemp("dbsnp")
    references = [corpus / "ref.fa", CORPUS / "reference.fa", exome / "ref.fa"]
    return {ref: write_dbsnp(ref, out) for ref in references}


def inputs(name, built, corpus, exome):
    """(bam, reference, interval list) of a named input."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa", corpus / "target.interval_list"
    if name == "exome":
        return exome / "exome.bam", exome / "ref.fa", exome / "targets.interval_list"
    if name in built:
        bam, ref = built[name]
        return bam, ref, corpus / "target.interval_list"
    return CORPUS / f"{name}.bam", CORPUS / "reference.fa", CORPUS / "baits.interval_list"


def run(cmd, picard):
    """Run one command; Picard through its wrapper with the shim Rscript."""
    if picard:
        return subprocess.run(["picard", *cmd, "VERBOSITY=ERROR", "QUIET=true"], env=picard_env(),
                              capture_output=True, text=True)
    return subprocess.run([str(BAMMETRICS), *cmd], capture_output=True, text=True)


def convert_args(prefix, ext):
    """The converter's inputs: INPUT_BASE, or the two detail files when EXT renamed them."""
    if not ext:
        return [f"I={prefix}", f"O={prefix}"]
    return [f"PRE_ADAPTER_IN={prefix}{DETAIL_FILES[0]}{ext}",
            f"BAIT_BIAS_IN={prefix}{DETAIL_FILES[1]}{ext}", f"OXOG_OUT={prefix}.oxog_metrics"]


def compare(bam, ref, args, label, tmp, bam_args=None):
    """Both collectors and both converters; every file compared. Returns cells compared."""
    bam_args = args if bam_args is None else bam_args
    ext = next((a.split("=", 1)[1] for a in args if a.startswith(("EXT=", "FILE_EXTENSION="))), "")
    converted = []
    for side, side_args, picard in (("picard", args, True), ("bammetrics", bam_args, False)):
        prefix = tmp / side
        done = run([TOOL, f"I={bam}", f"O={prefix}", f"R={ref}", SILENT, *side_args], picard)
        assert done.returncode == 0, f"{side} {TOOL} {label}: {done.stderr[-2000:]}"
        done = run([CONVERT, *convert_args(prefix, ext)], picard)
        converted.append(done)
    picard_conv, ours_conv = converted
    assert (picard_conv.returncode == 0) == (ours_conv.returncode == 0), (
        f"{CONVERT} after {label}: Picard exit {picard_conv.returncode}, bammetrics "
        f"{ours_conv.returncode}: {ours_conv.stderr[-500:]}")
    cells = compare_written(TOOL, tmp, label, SEEN)
    assert cells > 0, f"{TOOL} {label}: no cell compared"
    note = "five files" + (" and the OxoG file" if picard_conv.returncode == 0 else
                           "; both converters fail")
    RESULTS.append((f"{TOOL} {label}", cells, cells, note + " identical"))
    return cells


INPUTS = ["foundation", "t1", "t2", "t3", "exome", "decorated", "unplaced_middle"]


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("context", ["0", "1", "2"])
def test_artifact_context_sizes(built, corpus, exome, name, context, tmp_path):
    bam, ref, _ = inputs(name, built, corpus, exome)
    compare(bam, ref, [f"CONTEXT_SIZE={context}"], f"{name} CONTEXT_SIZE={context}", tmp_path)


INTERVAL_INPUTS = ["foundation", "t1", "t3", "exome", "decorated"]


@pytest.mark.parametrize("name", INTERVAL_INPUTS)
def test_artifact_intervals(built, corpus, exome, name, tmp_path):
    bam, ref, intervals = inputs(name, built, corpus, exome)
    compare(bam, ref, [f"INTERVALS={intervals}"], f"{name} INTERVALS", tmp_path)


@pytest.mark.parametrize("name", INTERVAL_INPUTS)
def test_artifact_db_snp(built, corpus, exome, dbsnp, name, tmp_path):
    bam, ref, _ = inputs(name, built, corpus, exome)
    plain, _ = dbsnp[ref if name != "decorated" else corpus / "ref.fa"]
    compare(bam, ref, [f"DB_SNP={plain}"], f"{name} DB_SNP", tmp_path)


@pytest.mark.parametrize("name", INTERVAL_INPUTS)
def test_artifact_db_snp_and_intervals(built, corpus, exome, dbsnp, name, tmp_path):
    """With INTERVALS, the dbSNP VCF is read by interval through its index."""
    bam, ref, intervals = inputs(name, built, corpus, exome)
    _, indexed = dbsnp[ref if name != "decorated" else corpus / "ref.fa"]
    compare(bam, ref, [f"DB_SNP={indexed}", f"INTERVALS={intervals}", "CONTEXT_SIZE=2"],
            f"{name} DB_SNP INTERVALS", tmp_path)


# Options that change a number, each at a non-default value.
OPTIONS = {
    "q10": ["MINIMUM_QUALITY_SCORE=10"],
    "q30": ["Q=30"],
    "mq0": ["MINIMUM_MAPPING_QUALITY=0"],
    "mq50": ["MQ=50"],
    "min_ins_0": ["MINIMUM_INSERT_SIZE=0"],
    "min_ins_250": ["MIN_INS=250"],
    "max_ins_0": ["MAXIMUM_INSERT_SIZE=0"],
    "max_ins_300": ["MAX_INS=300"],
    "unpaired": ["INCLUDE_UNPAIRED=true"],
    "dupes": ["INCLUDE_DUPLICATES=true"],
    "non_pf": ["INCLUDE_NON_PF_READS=true"],
    "tandem": ["TANDEM_READS=true"],
    "no_oq": ["USE_OQ=false"],
    "contexts": ["CONTEXTS_TO_PRINT=ACA", "CONTEXTS_TO_PRINT=TGT", "CONTEXTS_TO_PRINT=GGC"],
    "contexts_oxog": ["CONTEXT_SIZE=0", "CONTEXTS_TO_PRINT=C", "CONTEXTS_TO_PRINT=G"],
    "contexts_null": ["CONTEXTS_TO_PRINT=AAA", "CONTEXTS_TO_PRINT=null", "CONTEXTS_TO_PRINT=CCC"],
    "ext": ["FILE_EXTENSION=.txt"],
    "stop_1000": ["STOP_AFTER=1000"],
    "stop_5000": ["STOP_AFTER=5000", "CONTEXT_SIZE=0"],
    "everything": ["Q=5", "MQ=0", "UNPAIRED=true", "DUPES=true", "NON_PF=true", "TANDEM=true",
                   "USE_OQ=false"],
}
OPTION_INPUTS = ["foundation", "decorated"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", sorted(OPTIONS))
def test_artifact_options(built, corpus, exome, name, option, tmp_path):
    bam, ref, _ = inputs(name, built, corpus, exome)
    compare(bam, ref, OPTIONS[option], f"{name} {option}", tmp_path)


TIER_OPTIONS = {
    "t2": ["TANDEM_READS=true", "MAX_INS=0"],
    "exome": ["Q=30", "MQ=20", "DUPES=true", "CONTEXT_SIZE=2"],
}


@pytest.mark.parametrize("name", sorted(TIER_OPTIONS))
def test_artifact_tier_options(built, corpus, exome, name, tmp_path):
    bam, ref, _ = inputs(name, built, corpus, exome)
    compare(bam, ref, TIER_OPTIONS[name], f"{name} options", tmp_path)


def test_converter_reads_picard_files(corpus, tmp_path):
    """bammetrics' converter on Picard's five files writes Picard's OxoG file."""
    picard = run([TOOL, f"I={corpus / 'corpus.bam'}", f"O={tmp_path / 'source'}",
                  f"R={corpus / 'ref.fa'}", SILENT], True)
    assert picard.returncode == 0, picard.stderr
    for side, is_picard in (("picard", True), ("bammetrics", False)):
        done = run([CONVERT, f"I={tmp_path / 'source'}", f"O={tmp_path / side}"], is_picard)
        assert done.returncode == 0, done.stderr
    cells = compare_written(CONVERT, tmp_path, "on Picard's files")
    RESULTS.append((f"{CONVERT} on Picard's files", cells, cells, "identical"))


def test_multiple_metrics_runs_the_artifact_program(corpus, dbsnp, tmp_path):
    """The program of CollectMultipleMetrics with the run's INTERVALS, DB_SNP,
    INCLUDE_UNPAIRED, FILE_EXTENSION and an EXTRA_ARGUMENT."""
    ref = corpus / "ref.fa"
    args = [f"I={corpus / 'corpus.bam'}", f"R={ref}", "PROGRAM=null",
            f"PROGRAM={TOOL}", "PROGRAM=CollectQualityYieldMetrics", "EXT=.txt",
            "INCLUDE_UNPAIRED=true", f"INTERVALS={corpus / 'target.interval_list'}",
            f"DB_SNP={dbsnp[ref][1]}", f"EXTRA_ARGUMENT={TOOL}::CONTEXT_SIZE=2",
            f"EXTRA_ARGUMENT={TOOL}::R={ref}", SILENT]
    for side, picard in (("picard", True), ("bammetrics", False)):
        done = run(["CollectMultipleMetrics", *args, f"O={tmp_path / side}"], picard)
        assert done.returncode == 0, done.stderr[-2000:]
    cells = compare_written("CollectMultipleMetrics", tmp_path, "artifact program")
    RESULTS.append(("CollectMultipleMetrics artifact program", cells, cells, "six files identical"))


# --- error parity -----------------------------------------------------------------------


def both_fail(cmd, tmp_path, message, picard_message):
    """Both tools exit non-zero with their message, and neither writes a file."""
    picard = run([*cmd, f"O={tmp_path / 'picard'}"], True)
    assert picard.returncode != 0, f"Picard completed: {cmd}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr[-2000:]
    ours = run([*cmd, f"O={tmp_path / 'ours'}"], False)
    assert ours.returncode != 0, ours.stdout
    assert message in ours.stderr, ours.stderr
    assert not list(tmp_path.glob("picard*")) and not list(tmp_path.glob("ours*"))
    shown = " ".join(a for a in cmd[1:] if "=" in a and "/" not in a and a != SILENT)
    RESULTS.append((f"{cmd[0]} {shown}"[:60], 0, 0, "both exit non-zero, no file"))


def collector(bam, ref, *args):
    reference = [f"R={ref}"] if ref else []
    return [TOOL, f"I={bam}", *reference, SILENT, *args]


def test_reference_is_required(corpus, tmp_path):
    both_fail(collector(corpus / "corpus.bam", None), tmp_path,
              "Option 'REFERENCE_SEQUENCE' is required.", "Option 'REFERENCE_SEQUENCE' is required.")


@pytest.mark.parametrize("args,message", [
    (["CONTEXT_SIZE=-1"], "CONTEXT_SIZE cannot be negative"),
    (["MIN_INS=-1"], "MINIMUM_INSERT_SIZE cannot be negative"),
    (["MAX_INS=-5"], "MAXIMUM_INSERT_SIZE cannot be negative"),
    (["MAX_INS=50", "MIN_INS=100"],
     "MAXIMUM_INSERT_SIZE cannot be less than MINIMUM_INSERT_SIZE unless set to 0"),
    (["CONTEXTS_TO_PRINT=AAAA", "CONTEXTS_TO_PRINT=AC"],
     "Context AC is not the length implied by CONTEXT_SIZE: 3"),
    (["CONTEXT_SIZE=0", "CONTEXTS_TO_PRINT=AAA"],
     "Context AAA is not the length implied by CONTEXT_SIZE: 1"),
])
def test_command_line_validation(corpus, tmp_path, args, message):
    both_fail(collector(corpus / "corpus.bam", corpus / "ref.fa", *args), tmp_path, message, message)


def test_no_read_groups(built, tmp_path):
    bam, ref = built["no_read_groups"]
    message = "Record contains library that is missing from header: UnknownLibrary"
    both_fail(collector(bam, ref), tmp_path, message, message)


def test_bed_intervals(corpus, tmp_path):
    bed = tmp_path / "targets.bed"
    bed.write_text("chr1\t1000\t2000\n")
    message = "Interval list file must contain header"
    both_fail(collector(corpus / "corpus.bam", corpus / "ref.fa", f"INTERVALS={bed}"), tmp_path,
              message, message)


def test_db_snp_with_intervals_needs_an_index(corpus, dbsnp, tmp_path):
    plain, _ = dbsnp[corpus / "ref.fa"]
    message = "An index is required, but none found"
    both_fail(collector(corpus / "corpus.bam", corpus / "ref.fa", f"DB_SNP={plain}",
                        f"INTERVALS={corpus / 'target.interval_list'}"), tmp_path, message, message)


def test_db_snp_of_another_dictionary(corpus, exome, dbsnp, tmp_path):
    plain, _ = dbsnp[exome / "ref.fa"]
    message = "does not match sequence dictionary for (INPUT)"
    both_fail(collector(corpus / "corpus.bam", corpus / "ref.fa", f"DB_SNP={plain}"), tmp_path,
              message, message)


def test_reference_of_another_dictionary(corpus, exome, tmp_path):
    both_fail(collector(corpus / "corpus.bam", exome / "ref.fa"), tmp_path, "in the BAM header",
              "Exception")


def test_queryname_sorted(corpus, tmp_path):
    bam = tmp_path / "byname.bam"
    subprocess.run(["samtools", "sort", "-n", "-o", str(bam), str(corpus / "corpus.bam")],
                   check=True, capture_output=True)
    both_fail(collector(bam, corpus / "ref.fa", "ASSUME_SORTED=false"), tmp_path,
              "should be coordinate sorted", "should be coordinate sorted")
    both_fail(collector(bam, corpus / "ref.fa"), tmp_path, "the input must be coordinate sorted",
              "Requesting earlier reference sequence")


@pytest.fixture(scope="module")
def picard_files(corpus, tmp_path_factory):
    """Picard's five files for the foundation corpus."""
    out = tmp_path_factory.mktemp("picard_files")
    done = run([TOOL, f"I={corpus / 'corpus.bam'}", f"O={out / 'run'}", f"R={corpus / 'ref.fa'}",
                SILENT], True)
    assert done.returncode == 0, done.stderr
    return out / "run"


def converter_both_fail(args, tmp_path, message, picard_message):
    picard = run([CONVERT, *args, f"OXOG_OUT={tmp_path / 'picard.oxog_metrics'}"], True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr[-2000:]
    ours = run([CONVERT, *args, f"OXOG_OUT={tmp_path / 'ours.oxog_metrics'}"], False)
    assert ours.returncode != 0, ours.stdout
    assert message in ours.stderr, ours.stderr
    assert not (tmp_path / "picard.oxog_metrics").exists()
    assert not (tmp_path / "ours.oxog_metrics").exists()
    shown = " ".join(a.split("=")[0] + "=" + Path(a.split("=", 1)[1]).name for a in args)
    RESULTS.append((f"{CONVERT} {shown}"[:60], 0, 0, "both exit non-zero, no file"))


def test_converter_needs_an_input(tmp_path):
    message = "Must specify either INPUT_BASE or PRE_ADAPTER_IN"
    converter_both_fail([], tmp_path, message, message)


def test_converter_swapped_files(picard_files, tmp_path):
    converter_both_fail([f"PRE_ADAPTER_IN={picard_files}.bait_bias_detail_metrics",
                         f"BAIT_BIAS_IN={picard_files}.pre_adapter_detail_metrics"], tmp_path,
                        "cannot be cast", "ClassCastException")


def test_converter_missing_file(picard_files, tmp_path):
    converter_both_fail([f"I={picard_files}.nope"], tmp_path, "nope", "nope")


def test_converter_empty_pre_adapter_file(picard_files, tmp_path):
    empty = tmp_path / "empty.pre_adapter_detail_metrics"
    lines = Path(f"{picard_files}.pre_adapter_detail_metrics").read_text().splitlines()
    header = [line for line in lines if line.startswith("#") or line.startswith("SAMPLE_ALIAS")
              or not line]
    empty.write_text("\n".join(header[:7]) + "\n")
    converter_both_fail([f"PRE_ADAPTER_IN={empty}",
                         f"BAIT_BIAS_IN={picard_files}.bait_bias_detail_metrics"], tmp_path,
                        "has no rows", "IndexOutOfBoundsException")


# --- the distinct-value gate ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {
    ("picard.analysis.artifacts.ErrorSummaryMetrics", "REF_BASE"):
        "a G or T reference is folded into C or A",
    ("picard.analysis.CollectOxoGMetrics$CpcgMetrics", "TOTAL_SITES"): "always 0",
}
CLASSES = 6
# Six Picard files per case; contexts_null leaves the converter without GGG, so
# those cases have five.
EXPECTED_FILES = (6 * (len(INPUTS) * 3 + 3 * len(INTERVAL_INPUTS) + len(OPTIONS) * len(OPTION_INPUTS)
                       + len(TIER_OPTIONS))
                  - len(OPTION_INPUTS))


def narrowed(config):
    """True when -k, --deselect or a node id selected part of the module."""
    return bool(config.getoption("keyword") or config.getoption("deselect")
                or any("::" in arg for arg in config.args))


def test_every_column_takes_three_non_default_values(request):
    if len(SEEN) != EXPECTED_FILES:
        message = f"{TOOL}: {len(SEEN)} of {EXPECTED_FILES} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(f"{message}: a comparison above failed or EXPECTED_FILES is stale")
    values = defaultdict(set)
    columns = set()
    for _, parsed in SEEN:
        for column in parsed.columns:
            columns.add((parsed.metrics_class, column))
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[(parsed.metrics_class, column)].add(value)
    short = {c: sorted(values[c]) for c in columns
             if len(values[c]) < 3 and c not in FEWER_THAN_THREE}
    assert not short, f"columns with fewer than 3 distinct non-default values: {short}"
    assert len({c for c, _ in columns}) == CLASSES
    fewest = min((c for c in columns if c not in FEWER_THAN_THREE), key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest[1]} {len(values[fewest])}"))

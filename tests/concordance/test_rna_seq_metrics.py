"""CollectRnaSeqMetrics against Picard 3.5.0: every column at every accumulation
level, the normalized coverage histogram, every option that changes a number, and
error parity.

Inputs: the rna_bams.py corpus (a stranded and an unstranded BAM over a refFlat
with multi-exon transcripts on both strands, overlapping genes, a transcript
under MINIMUM_LENGTH, dropped annotation lines and rRNA intervals; an empty BAM)
and the foundation corpus with its refflat.txt and target.interval_list as the
ribosomal intervals. Every case runs Picard with VALIDATION_STRINGENCY=SILENT.

The last test checks that every column took at least 3 distinct non-default
values over all Picard files compared, except the columns in FEWER_THAN_THREE
with the reason.
"""
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS
from compare import compare_metrics
from conftest import RESULTS
from rna_bams import RRNA, build_all
from run_picard import picard_env

TOOL = "CollectRnaSeqMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []
NONE, FIRST, SECOND = ("STRAND=" + s for s in ("NONE", "FIRST_READ_TRANSCRIPTION_STRAND",
                                              "SECOND_READ_TRANSCRIPTION_STRAND"))


@pytest.fixture(scope="module")
def rna(tmp_path_factory):
    """{name: bam} of the rna_bams corpus plus its annotation files."""
    out = tmp_path_factory.mktemp("rna")
    bams = build_all(out)
    bams["refflat"] = out / "rna.refFlat"
    bams["rrna"] = out / "rrna.interval_list"
    bams["dir"] = out
    return bams


def inputs(name, rna, corpus):
    """(bam, refFlat, ribosomal interval list) of a named input."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "refflat.txt", corpus / "target.interval_list"
    return rna[name], rna["refflat"], rna["rrna"]


def compare(bam, refflat, args, label, bam_args=None):
    args = [SILENT, f"REF_FLAT={refflat}", *args]
    bam_args = None if bam_args is None else [SILENT, f"REF_FLAT={refflat}", *bam_args]
    cells = compare_metrics(TOOL, bam, args, bam_args=bam_args, workdir=bam.parent, expect_file=True,
                            seen=SEEN)
    assert cells > 0, f"{TOOL} {label}: no cell compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every cell and the file body identical"))
    return cells


INPUTS = ["stranded", "unstranded", "empty", "foundation"]
STRANDS = [NONE, FIRST, SECOND]


@pytest.mark.parametrize("name", INPUTS)
@pytest.mark.parametrize("strand", STRANDS)
def test_rna_every_strand_specificity(rna, corpus, name, strand):
    bam, refflat, _ = inputs(name, rna, corpus)
    compare(bam, refflat, [strand], f"{name} {strand}")


@pytest.mark.parametrize("name", INPUTS)
def test_rna_ribosomal_intervals(rna, corpus, name):
    """RIBOSOMAL_INTERVALS fills RIBOSOMAL_BASES and PCT_RIBOSOMAL_BASES."""
    bam, refflat, rrna = inputs(name, rna, corpus)
    compare(bam, refflat, [SECOND, f"RIBOSOMAL_INTERVALS={rrna}"], f"{name} rRNA")


# Options that change a number, each at a non-default value; `rrna` adds the
# ribosomal intervals.
OPTIONS = {
    "min_200": (["MINIMUM_LENGTH=200"], False),
    "min_2500": (["MINIMUM_LENGTH=2500"], False),
    "min_50": (["MINIMUM_LENGTH=50", "END_BIAS_BASES=20"], False),
    "bias_50": (["END_BIAS_BASES=50"], False),
    "bias_300": (["END_BIAS_BASES=300"], False),
    "pct_0.5": (["RRNA_FRAGMENT_PERCENTAGE=0.5"], True),
    "pct_0.3": (["RRNA_FRAGMENT_PERCENTAGE=0.3"], True),
    "pct_0": (["RRNA_FRAGMENT_PERCENTAGE=0"], True),
    "pct_1": (["RRNA_FRAGMENT_PERCENTAGE=1"], True),
    "ignore_cM": (["IGNORE_SEQUENCE=cM"], False),
    "ignore_cM_cR": (["IGNORE_SEQUENCE=cM", "IGNORE_SEQUENCE=cR"], True),
    "ignore_levels": (["IGNORE_SEQUENCE=cM", "LEVEL=SAMPLE", "LEVEL=LIBRARY"], False),
    "stop_300": (["STOP_AFTER=300"], False),
    "stop_1500": (["STOP_AFTER=1500"], True),
    "chart": (["CHART=chart.pdf"], False),
    "level_sample": (["LEVEL=SAMPLE"], False),
    "level_library": (["LEVEL=LIBRARY"], True),
    "level_read_group": (["LEVEL=READ_GROUP"], False),
    "level_all": (["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"], True),
    "level_null_rg": (["LEVEL=null", "LEVEL=READ_GROUP", "MINIMUM_LENGTH=200"], False),
    "level_null_sl": (["LEVEL=null", "LEVEL=SAMPLE", "LEVEL=LIBRARY", "END_BIAS_BASES=200"], True),
}
OPTION_INPUTS = ["stranded", "unstranded"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", sorted(OPTIONS))
def test_rna_options(rna, corpus, name, option, tmp_path):
    args, with_rrna = OPTIONS[option]
    bam, refflat, rrna = inputs(name, rna, corpus)
    extra = [f"RIBOSOMAL_INTERVALS={rrna}"] if with_rrna else []
    args = [a.replace("chart.pdf", str(tmp_path / "chart.pdf")) for a in args]
    strand = FIRST if name == "stranded" else SECOND
    compare(bam, refflat, [strand, *extra, *args], f"{name} {option}")


FOUNDATION_OPTIONS = {
    "levels_rrna": ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP", "MINIMUM_LENGTH=100", "rrna"],
    "ignore": ["IGNORE_SEQUENCE=chrM", "MINIMUM_LENGTH=100", "END_BIAS_BASES=30"],
    "stop": ["STOP_AFTER=2000", "rrna", "RRNA_FRAGMENT_PERCENTAGE=0.2"],
}


@pytest.mark.parametrize("option", sorted(FOUNDATION_OPTIONS))
def test_rna_foundation_options(corpus, option):
    bam, refflat, rrna = inputs("foundation", None, corpus)
    args = [f"RIBOSOMAL_INTERVALS={rrna}" if a == "rrna" else a for a in FOUNDATION_OPTIONS[option]]
    compare(bam, refflat, [FIRST, *args], f"foundation {option}")


def test_rna_bed_ribosomal_intervals_are_read_by_bammetrics_only(rna, corpus, tmp_path):
    """bammetrics reads a BED as RIBOSOMAL_INTERVALS (converted as BedToIntervalList
    would); Picard 3.5.0 rejects it, so Picard gets the interval list."""
    bam, refflat, rrna = inputs("stranded", rna, corpus)
    bed = tmp_path / "rrna.bed"
    bed.write_text("".join(f"{c}\t{s - 1}\t{e}\t{n}\n" for c, s, e, n in RRNA))
    compare(bam, refflat, [SECOND, f"RIBOSOMAL_INTERVALS={rrna}"], "stranded rRNA BED",
            bam_args=[SECOND, f"RIBOSOMAL_INTERVALS={bed}"])
    picard = subprocess.run(["picard", TOOL, f"I={bam}", f"O={tmp_path / 'p.txt'}", f"REF_FLAT={refflat}",
                             SECOND, f"RIBOSOMAL_INTERVALS={bed}", "QUIET=true", "VERBOSITY=ERROR"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0
    assert "Interval list file must contain header" in picard.stderr + picard.stdout


# --- error parity -----------------------------------------------------------------------


def both_fail(bam, args, tmp_path, message, picard_message):
    """Both tools exit non-zero with `message` (bammetrics) and `picard_message`
    (Picard) in their output, and neither writes OUTPUT."""
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    picard = subprocess.run(["picard", TOOL, f"I={bam}", f"O={picard_out}", SILENT, *args, "QUIET=true",
                             "VERBOSITY=ERROR"], env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={ours_out}", *args],
                            capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    assert message in result.stderr, result.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(a for a in args if '=' in a)[:60]}", 0, 0,
                    "both exit non-zero, no file"))


def test_strand_specificity_is_required(rna, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}"], tmp_path,
              "Option 'STRAND_SPECIFICITY' is required.", "Option 'STRAND_SPECIFICITY' is required.")


def test_ref_flat_is_required(rna, tmp_path):
    both_fail(rna["stranded"], [NONE], tmp_path, "Option 'REF_FLAT' is required.",
              "Option 'REF_FLAT' is required.")


def test_bad_strand_value(rna, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}", "STRAND=BOTH"], tmp_path,
              "'BOTH' is not a valid value for StrandSpecificity.",
              "'BOTH' is not a valid value for StrandSpecificity.")


def test_unknown_ignore_sequence(rna, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}", NONE, "IGNORE_SEQUENCE=chrM"], tmp_path,
              "Unrecognized sequence chrM passed as argument to IGNORE_SEQUENCE",
              "Unrecognized sequence chrM passed as argument to IGNORE_SEQUENCE")


def test_zero_percentage_needs_ribosomal_intervals(rna, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}", NONE, "RRNA_FRAGMENT_PERCENTAGE=0"], tmp_path,
              "Must use a RIBOSOMAL_INTERVALS file if RRNA_FRAGMENT_PERCENTAGE = 0.0",
              "Must use a RIBOSOMAL_INTERVALS file if RRNA_FRAGMENT_PERCENTAGE = 0.0")


def test_negative_percentage_without_intervals_fails_in_both(rna, tmp_path):
    """Every fragment reaches a negative percentage, and Picard adds to a null
    RIBOSOMAL_BASES (NullPointerException)."""
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}", NONE, "RRNA_FRAGMENT_PERCENTAGE=-1"], tmp_path,
              "without RIBOSOMAL_INTERVALS", "NullPointerException")


def test_ribosomal_intervals_of_another_dictionary(rna, corpus, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={rna['refflat']}", NONE,
                                f"RIBOSOMAL_INTERVALS={corpus / 'target.interval_list'}"], tmp_path,
              "Sequence dictionaries differ", "Sequence dictionaries differ")


SHORT_LINE = "GZ\tGZ-1\tc1\t+\t100\t900\t200\t800\t1\t99,\n"


def test_ref_flat_with_a_short_last_line(rna, tmp_path):
    """Picard pads a later line with fewer fields than the first with nulls and
    fails splitting the null exon ends; bammetrics names the line."""
    refflat = tmp_path / "bad.refFlat"
    refflat.write_text(rna["refflat"].read_text() + SHORT_LINE)
    both_fail(rna["stranded"], [f"REF_FLAT={refflat}", NONE], tmp_path, "Wrong number of fields",
              "NullPointerException")


def test_ref_flat_with_a_short_first_line(rna, tmp_path):
    """The first line sets Picard's field count, so every later line overflows it."""
    refflat = tmp_path / "bad.refFlat"
    refflat.write_text(SHORT_LINE + rna["refflat"].read_text())
    both_fail(rna["stranded"], [f"REF_FLAT={refflat}", NONE], tmp_path, "Wrong number of fields",
              "Unexpected number of elements found when parsing file")


def test_missing_ref_flat_file(rna, tmp_path):
    both_fail(rna["stranded"], [f"REF_FLAT={tmp_path / 'none.refFlat'}", NONE], tmp_path,
              "none.refFlat", "none.refFlat")


def test_assume_sorted_false_on_a_queryname_sorted_bam(rna, tmp_path):
    bam = tmp_path / "byname.bam"
    subprocess.run(["samtools", "sort", "-n", "-o", str(bam), str(rna["stranded"])], check=True,
                   capture_output=True)
    both_fail(bam, [f"REF_FLAT={rna['refflat']}", NONE, "ASSUME_SORTED=false"], tmp_path,
              "should be coordinate sorted", "should be coordinate sorted")


# --- the distinct-value gate ----------------------------------------------------------

DEFAULTS = {"0", ""}
# Columns that cannot take 3 distinct non-default values, and why.
FEWER_THAN_THREE = {}
EXPECTED_FILES = (len(INPUTS) * len(STRANDS) + len(INPUTS) + len(OPTIONS) * len(OPTION_INPUTS)
                  + len(FOUNDATION_OPTIONS) + 1)


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
    columns = []
    histogram_rows = 0
    for parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns:
            if column not in columns:
                columns.append(column)
        histogram_rows += len(parsed.histogram)
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3 and c not in FEWER_THAN_THREE}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    assert histogram_rows > 0, "no normalized coverage histogram compared"
    assert len(columns) == 30 and SEEN[0].metrics_class == "picard.analysis.RnaSeqMetrics"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

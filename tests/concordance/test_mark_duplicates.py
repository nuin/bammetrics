"""MarkDuplicates against Picard 3.5.0: the DuplicationMetrics file cell by
cell and byte for byte from `## METRICS CLASS` on, and the BAM record by
record. Every case runs Picard with VALIDATION_STRINGENCY=SILENT.

The BAMs are compared as `samtools view --no-PG` text: the header lines
other than @PG must be equal, and so must every record, in file order, with
every field and every tag (the PG tag included). The @PG lines must carry
the same IDs and PP links; their VN and CL name the program that wrote
them. BAM bytes are never compared.

Inputs: the markdup_bams BAMs (optical sets, barcodes, UMIs, clips, every
record kind, a query-name copy, a second input with colliding read group
and program IDs), the complexity_bams, read_level_bams and asm_bams BAMs,
the foundation corpus, the exome slice, the bwa tiers t0 to t3, a CRAM.
Where Picard fails (a pair without RG, a file htsjdk's writer refuses),
bammetrics must fail too and neither may write a metrics file.

The last test checks that every metrics and histogram column took at
least 3 distinct non-default values over all Picard files compared.
"""
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pytest

from asm_bams import build_all as build_asm
from cases import BAMMETRICS, CORPUS
from compare import body, diff_metrics, first_difference
from complexity_bams import build_all as build_complexity
from conftest import RESULTS
from markdup_bams import build_all as build_markdup
from picard_metrics import parse_picard_metrics
from read_level_bams import build_all as build_read_level
from run_picard import picard_env

TOOL = "MarkDuplicates"
SILENT = "VALIDATION_STRINGENCY=SILENT"
SEEN = []


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    bams = {}
    for prefix, build in (("md", build_markdup), ("asm", build_asm),
                          ("rl", build_read_level), ("elc", build_complexity)):
        directory = tmp_path_factory.mktemp(prefix)
        for name, bam in build(directory).items():
            bams[f"{prefix}_{name}"] = bam
    reference = bams["md_markdup"].parent / "ref.fa"
    cram = bams["md_markdup"].with_suffix(".cram")
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(reference), "-o", str(cram),
                    str(bams["md_markdup"])], check=True)
    bams["md_markdup_cram"] = cram
    return bams


def bam_path(name, built, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    if name in built:
        return built[name]
    return CORPUS / f"{name}.bam"


def view(path, header=False):
    args = ["samtools", "view", "--no-PG", *(["-H"] if header else []), str(path)]
    return subprocess.run(args, check=True, capture_output=True).stdout.decode().splitlines()


def programs(lines):
    """The ID and PP of each @PG line."""
    out = []
    for line in lines:
        if line.startswith("@PG"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            out.append((fields.get("ID"), fields.get("PP")))
    return out


def run_both(inputs, args, tmp):
    picard_bam, ours_bam = tmp / "picard.bam", tmp / "bammetrics.bam"
    picard_txt, ours_txt = tmp / "picard.txt", tmp / "bammetrics.txt"
    ins = [f"I={i}" for i in inputs]
    picard = subprocess.run(["picard", TOOL, *ins, f"O={picard_bam}", f"M={picard_txt}", SILENT,
                             *args, "VERBOSITY=ERROR", "QUIET=true"], env=picard_env(),
                            capture_output=True, text=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, *ins, f"O={ours_bam}", f"M={ours_txt}", *args],
                          capture_output=True, text=True)
    return picard, ours, (picard_bam, ours_bam, picard_txt, ours_txt)


def compare(inputs, args, label, outcome="equal"):
    """Run both tools. `outcome` "equal": both succeed with equal output;
    "fail": both fail and neither writes a metrics file; "same": either of
    the two, decided by Picard. Returns the cells compared."""
    inputs = [inputs] if isinstance(inputs, Path) else list(inputs)
    with tempfile.TemporaryDirectory(dir=inputs[0].parent) as tmp:
        picard, ours, (picard_bam, ours_bam, picard_txt, ours_txt) = run_both(
            inputs, args, Path(tmp))
        picard_failed = picard.returncode != 0
        if outcome == "fail" or (outcome == "same" and picard_failed):
            assert picard_failed, f"{label}: Picard succeeded"
            assert ours.returncode != 0, f"{label}: Picard failed, bammetrics succeeded"
            assert not picard_txt.exists() and not ours_txt.exists(), f"{label}: metrics written"
            RESULTS.append((f"{TOOL} {label}", 0, 0, "both fail, no metrics file"))
            return 0
        assert not picard_failed, f"{label}: Picard failed:\n{picard.stderr[-2000:]}"
        assert ours.returncode == 0, f"{label}: bammetrics failed: {ours.stderr}"
        picard_metrics, ours_metrics = parse_picard_metrics(picard_txt), parse_picard_metrics(ours_txt)
        problems = diff_metrics(picard_metrics, ours_metrics)
        picard_body, ours_body = body(picard_txt), body(ours_txt)
        if not problems and picard_body != ours_body:
            problems.append("metrics file after the header lines differs: "
                            + first_difference(picard_body, ours_body))
        picard_header, ours_header = view(picard_bam, True), view(ours_bam, True)
        if [h for h in picard_header if not h.startswith("@PG")] != \
                [h for h in ours_header if not h.startswith("@PG")]:
            problems.append("BAM header lines other than @PG differ")
        if programs(picard_header) != programs(ours_header):
            problems.append(f"@PG IDs: Picard {programs(picard_header)}, "
                            f"bammetrics {programs(ours_header)}")
        picard_records, ours_records = view(picard_bam), view(ours_bam)
        if picard_records != ours_records:
            for i, (p, o) in enumerate(zip(picard_records, ours_records)):
                if p != o:
                    problems.append(f"record {i}:\n    Picard     {p}\n    bammetrics {o}")
                    break
            else:
                problems.append(f"record counts: Picard {len(picard_records)}, "
                                f"bammetrics {len(ours_records)}")
        assert not problems, f"{TOOL} {label}:\n  " + "\n  ".join(problems[:20])
    SEEN.append(picard_metrics)
    duplicates = sum(1 for r in picard_records if int(r.split("\t")[1]) & 0x400)
    cells = (len(picard_metrics.rows) * len(picard_metrics.columns)
             + len(picard_metrics.histogram) * len(picard_metrics.histogram_columns))
    RESULTS.append((f"{TOOL} {label}", cells, cells,
                    f"{len(picard_records)} records equal, {duplicates} flagged"))
    return cells


# --- every input at the defaults ---------------------------------------------------

# Picard's outcome decides: equal output, or both fail (a paired record
# without RG is a NullPointerException in Picard; samtools' natural name
# order is out of order for htsjdk's writer).
INPUTS = [
    "md_markdup", "md_markdup_queryname", "md_markdup_second", "md_markdup_second_queryname",
    "md_markdup_duplex", "md_markdup_duplex_queryname",
    "elc_illumina", "elc_unaligned", "elc_unaligned_queryname", "elc_second_file",
    "elc_plain_names", "elc_elc_empty", "rl_reads", "rl_single_library", "rl_all_duplicates",
    "rl_no_duplicates", "rl_star_fields", "rl_queryname", "rl_empty", "rl_all_qcfail",
    "rl_q0_cycle", "rl_second_only", "asm_mix", "asm_unpaired_only", "asm_single_unmapped",
    "asm_clips_and_flags", "asm_no_pf_group", "asm_all_qcfail", "asm_empty",
    "foundation", "exome", "t0_asm", "t0_all_tandem", "t0_filters", "t0_single_pair", "t0_spread",
    "t1", "t2", "t3",
]
# Inputs that must compare equal: a case Picard fails on would hide a regression.
MUST_SUCCEED = {"md_markdup", "md_markdup_queryname", "md_markdup_second",
                "md_markdup_second_queryname", "md_markdup_duplex", "md_markdup_duplex_queryname",
                "elc_unaligned", "elc_plain_names", "elc_second_file", "foundation", "exome",
                "t1", "t2", "t3"}


@pytest.mark.parametrize("name", INPUTS)
def test_defaults(built, corpus, exome, name):
    path = bam_path(name, built, corpus, exome)
    compare(path, [], name, "equal" if name in MUST_SUCCEED else "same")


def test_cram(built):
    reference = built["md_markdup"].parent / "ref.fa"
    compare(built["md_markdup_cram"], [f"R={reference}"], "md_markdup.cram")
    compare(built["md_markdup_cram"], [f"R={reference}", "CLEAR_DT=false",
                                       "PROGRAM_RECORD_ID=null"], "md_markdup.cram untouched")


# --- several inputs ------------------------------------------------------------------

def test_two_inputs(built):
    compare([built["md_markdup"], built["md_markdup_second"]], [], "markdup with second")
    compare([built["md_markdup_second"], built["md_markdup"]], ["TAG_DUPLICATE_SET_MEMBERS=true"],
            "second with markdup")


def test_two_queryname_inputs(built):
    compare([built["md_markdup_queryname"], built["md_markdup_second_queryname"]],
            ["TAGGING_POLICY=All"], "queryname with queryname")


def test_same_input_twice_puts_a_pair_twice(built):
    """The cross-contig pairs of one file given twice meet in Picard's
    CoordinateSortedPairInfoMap: both tools fail."""
    compare([built["md_markdup"], built["md_markdup"]], [], "markdup twice", "fail")


def test_inputs_in_two_orders_fail(built):
    compare([built["md_markdup"], built["md_markdup_queryname"]], [], "coordinate with queryname",
            "fail")


# --- every option --------------------------------------------------------------------

OPTIONS = {
    "remove_duplicates": ["REMOVE_DUPLICATES=true"],
    "remove_sequencing": ["REMOVE_SEQUENCING_DUPLICATES=true"],
    "tag_all": ["TAGGING_POLICY=All"],
    "tag_optical": ["TAGGING_POLICY=OpticalOnly"],
    "tag_all_remove": ["TAGGING_POLICY=All", "REMOVE_DUPLICATES=true"],
    "tag_optical_remove_sequencing": ["TAGGING_POLICY=OpticalOnly",
                                      "REMOVE_SEQUENCING_DUPLICATES=true"],
    "tag_all_remove_both": ["TAGGING_POLICY=All", "REMOVE_DUPLICATES=true",
                            "REMOVE_SEQUENCING_DUPLICATES=true"],
    "set_members": ["TAG_DUPLICATE_SET_MEMBERS=true"],
    "keep_dt": ["CLEAR_DT=false"],
    "keep_dt_no_pg_tag": ["CLEAR_DT=false", "ADD_PG_TAG_TO_READS=false"],
    "no_program": ["PROGRAM_RECORD_ID=null"],
    "program_collides": ["PROGRAM_RECORD_ID=pgA"],
    "program_fields": ["PROGRAM_GROUP_NAME=MD", "PROGRAM_GROUP_VERSION=9",
                       "PROGRAM_GROUP_COMMAND_LINE=cmd"],
    "no_pg_tag": ["ADD_PG_TAG_TO_READS=false"],
    "reference_length": ["DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH"],
    "random": ["DUPLICATE_SCORING_STRATEGY=RANDOM"],
    "regex_null": ["READ_NAME_REGEX=null"],
    "regex_seven": ["READ_NAME_REGEX=[^:]+:[^:]+:[^:]+:[^:]+:([0-9]+):([0-9]+):([0-9]+)"],
    "regex_partial": ["READ_NAME_REGEX=(?:.*:)?([0-9]+):([0-9]+):([0-9a-z]+)"],
    "pixels_2500": ["OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500"],
    "pixels_0": ["OPTICAL_DUPLICATE_PIXEL_DISTANCE=0"],
    "set_size_3": ["MAX_OPTICAL_DUPLICATE_SET_SIZE=3"],
    "set_size_off": ["MAX_OPTICAL_DUPLICATE_SET_SIZE=-1"],
    "umi": ["BARCODE_TAG=RX"],
    "umi_identifier": ["BARCODE_TAG=RX", "MOLECULAR_IDENTIFIER_TAG=MI"],
    "read_one_barcode": ["READ_ONE_BARCODE_TAG=BX"],
    "read_two_barcode": ["READ_TWO_BARCODE_TAG=BX"],
    "all_barcodes": ["BARCODE_TAG=BC", "READ_ONE_BARCODE_TAG=BX", "READ_TWO_BARCODE_TAG=BX"],
    "comments": ["COMMENT=hello", "COMMENT=world"],
    "assume_queryname": ["ASSUME_SORT_ORDER=queryname"],
    "small_sort_buffers": ["MAX_RECORDS_IN_RAM=7"],
    "compression_and_ignored": ["COMPRESSION_LEVEL=1", "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=10",
                                "SORTING_COLLECTION_SIZE_RATIO=0.1",
                                "MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=5"],
}
COORDINATE_ONLY = {"assume_coordinate": ["ASSUME_SORT_ORDER=coordinate"],
                   "assume_sorted": ["ASSUME_SORTED=true"]}
OPTION_INPUTS = ["md_markdup", "md_markdup_queryname"]


@pytest.mark.parametrize("name", OPTION_INPUTS)
@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, name, option):
    compare(built[name], OPTIONS[option], f"{name} {option}")


@pytest.mark.parametrize("option", COORDINATE_ONLY)
def test_assumed_coordinate_order(built, option):
    compare(built["md_markdup"], COORDINATE_ONLY[option], f"md_markdup {option}")


DUPLEX = {
    "duplex": ["BARCODE_TAG=RX", "DUPLEX_UMI=true"],
    "duplex_identifier": ["BARCODE_TAG=RX", "DUPLEX_UMI=true", "MOLECULAR_IDENTIFIER_TAG=MI"],
    "duplex_everything": ["BARCODE_TAG=RX", "DUPLEX_UMI=true", "MOLECULAR_IDENTIFIER_TAG=MI",
                          "TAG_DUPLICATE_SET_MEMBERS=true", "TAGGING_POLICY=All",
                          "READ_ONE_BARCODE_TAG=BX"],
}


@pytest.mark.parametrize("name", ["md_markdup_duplex", "md_markdup_duplex_queryname"])
@pytest.mark.parametrize("option", DUPLEX)
def test_duplex_umis(built, name, option):
    compare(built[name], DUPLEX[option], f"{name} {option}")


FOUNDATION_OPTIONS = {
    "tag_all_set_members": ["TAGGING_POLICY=All", "TAG_DUPLICATE_SET_MEMBERS=true"],
    "remove_sequencing": ["REMOVE_SEQUENCING_DUPLICATES=true", "TAGGING_POLICY=OpticalOnly"],
    "reference_length": ["DUPLICATE_SCORING_STRATEGY=TOTAL_MAPPED_REFERENCE_LENGTH"],
    "random_regex_null": ["DUPLICATE_SCORING_STRATEGY=RANDOM", "READ_NAME_REGEX=null"],
}


@pytest.mark.parametrize("option", FOUNDATION_OPTIONS)
def test_foundation_options(corpus, option):
    compare(corpus / "corpus.bam", FOUNDATION_OPTIONS[option], f"foundation {option}")


def test_exome_options(exome):
    compare(exome / "exome.bam", ["TAGGING_POLICY=All", "TAG_DUPLICATE_SET_MEMBERS=true",
                                  "MAX_RECORDS_IN_RAM=1000"], "exome tagged, small buffers")


# --- error parity --------------------------------------------------------------------

FAILURES = {
    "unsorted": ("md_markdup_unsorted", []),
    "assume_unsorted": ("md_markdup", ["ASSUME_SORT_ORDER=unsorted"]),
    "natural_name_order": ("md_markdup_natural", []),
    "queryname_assumed_coordinate": ("md_markdup_queryname", ["ASSUME_SORTED=true"]),
    "no_read_group": ("md_markdup_no_rg", []),
    "assume_both": ("md_markdup", ["ASSUME_SORTED=true", "ASSUME_SORT_ORDER=coordinate"]),
    "duplex_unpaired_read": ("md_markdup", ["BARCODE_TAG=RX", "DUPLEX_UMI=true"]),
    "umi_illegal_characters": ("md_markdup", ["BARCODE_TAG=PG"]),
    "barcode_not_a_string": ("md_markdup", ["BARCODE_TAG=XS"]),
    "barcode_tag_length": ("md_markdup", ["BARCODE_TAG=RXX"]),
    "bad_regex": ("md_markdup", ["READ_NAME_REGEX=("]),
    "bad_strategy": ("md_markdup", ["DUPLICATE_SCORING_STRATEGY=BEST"]),
    "bad_policy": ("md_markdup", ["TAGGING_POLICY=Some"]),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, case):
    name, args = FAILURES[case]
    compare(built[name], args, f"{name} {case}", "fail")


def test_missing_metrics_file_fails(built, tmp_path):
    bam = built["md_markdup"]
    picard = subprocess.run(["picard", TOOL, f"I={bam}", f"O={tmp_path / 'p.bam'}", SILENT,
                             "VERBOSITY=ERROR", "QUIET=true"], env=picard_env(),
                            capture_output=True, text=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"O={tmp_path / 'o.bam'}"],
                          capture_output=True, text=True)
    assert picard.returncode != 0 and ours.returncode != 0
    assert "METRICS_FILE" in ours.stderr
    assert not (tmp_path / "p.bam").exists() and not (tmp_path / "o.bam").exists()
    RESULTS.append((f"{TOOL} no METRICS_FILE", 0, 0, "both fail, no BAM"))


@pytest.mark.parametrize("option", ["FLOW_MODE=true", "FLOW_SKIP_FIRST_N_FLOWS=2"])
def test_flow_options_are_refused(built, tmp_path, option):
    """KNOWN_DIFFERENCES.md: the flow-based options are not implemented, so
    bammetrics refuses them where Picard runs."""
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={built['md_markdup']}",
                           f"O={tmp_path / 'o.bam'}", f"M={tmp_path / 'o.txt'}", option],
                          capture_output=True, text=True)
    assert ours.returncode != 0 and option.split("=")[0] in ours.stderr
    assert not (tmp_path / "o.bam").exists()
    RESULTS.append((f"{TOOL} {option}", 0, 0, "refused as not implemented"))


def test_missing_input_fails(tmp_path):
    compare(tmp_path / "nope.bam", [], "missing input", "fail")


# --- every column exercised ----------------------------------------------------------

DEFAULTS = {"0", "", "0.0"}


def test_every_column_takes_three_non_default_values():
    if len(SEEN) < 40:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    columns = []
    for parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for key, row in parsed.histogram.items():
            values[parsed.histogram_columns[0]].add(key)
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[column].add(value)
        for column in parsed.columns + parsed.histogram_columns:
            if column not in columns:
                columns.append(column)
    short = {c: sorted(values[c]) for c in columns if len(values[c]) < 3}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest} {len(values[fewest])}"))

"""RevertOriginalBaseQualitiesAndAddMateCigar against Picard 3.5.0.

Picard first reads up to MAX_RECORDS_TO_EXAMINE records and may decide the file needs
nothing (no `OQ`, `MC` on the first paired record with a mapped mate): it then writes no
file and exits 0. Otherwise each record gets CleanSam's clip and `QUAL` from `OQ`, the
records are sorted by name, each pair gets its mate information
(`SamPairUtil.SetMateInfoIterator`), and a writer that sorts writes them in SORT_ORDER.
Every written case goes through `compare.compare_bam` with `record_bytes=True` (tier B:
header lines, `@PG` ID and PP, every record as samtools prints it, then every record as
bytes after inflation: Picard reads with `EAGERLY_DECODE`, so every record is encoded
again, tags in htsjdk's order and widths).

Inputs: the `revert_bams.py` records (coordinate, queryname and unsorted BAM, SAM text,
CRAM 3.0 and 3.1), the fixture BAMs, `markdup.sam`, `insert_size.cram`, t1 to t3, the
foundation corpus and the exome slice (both hold two primary first reads under one name,
so both tools fail). Options: each SORT_ORDER, RESTORE_ORIGINAL_QUALITIES, MAX_RECORDS_IN_RAM
(runs of 20: one merge level, see `test_add_or_replace_read_groups.py`), TMP_DIR,
COMPRESSION_LEVEL, CREATE_INDEX and CREATE_MD5_FILE (on by default), SAM text output.
The early exit: files Picard skips (both exit 0, neither writes), MAX_RECORDS_TO_EXAMINE
0, 3 and 10000. Error parity: the `revert_bams.FAILS` files, a bad sort order,
`MAX_RECORDS_IN_RAM=0`, ASSUME_SORTED, a missing input, a missing output directory, a
reference or mate reference index outside the dictionary (`revert_bams.write_far`), the
default index and MD5 beside standard output. The step in a chain against Picard's file. The distinct-value gate wants 50 qualities
restored from `OQ`, 200 `MC` values set, 4 output orders.
"""
import subprocess

import pytest

import revert_bams
from cases import FIXTURES
from conftest import RESULTS
from test_add_comments_to_bam import (SOURCES, chain_against_picard, compare, index_and_md5,
                                      run_tool, source_path)
from test_add_or_replace_read_groups import standard_output

TOOL = "RevertOriginalBaseQualitiesAndAddMateCigar"
REF = FIXTURES / "ref_100k.fa"
RUN = ["MAX_RECORDS_TO_EXAMINE=0", "VALIDATION_STRINGENCY=SILENT"]

# What the equal cases changed: records whose QUAL came from OQ, MC values set, orders.
SEEN = {"restored": 0, "mc": 0, "orders": set()}


@pytest.fixture(scope="module")
def records(tmp_path_factory):
    return revert_bams.write(tmp_path_factory.mktemp("revert"))


def note(source, written, order):
    """Count the records whose QUAL came from OQ and the MC values set."""
    before = subprocess.run(["samtools", "view", str(source)], check=True, capture_output=True,
                            text=True).stdout.splitlines()
    SEEN["restored"] += sum("\tOQ:Z:" in line and not line.endswith("OQ:Z:") for line in before)
    SEEN["mc"] += sum("\tMC:Z:" in line for line in written.records)
    SEEN["orders"].add(order)


@pytest.mark.parametrize("kind", ["bam", "queryname", "unsorted", "sam"])
def test_the_generated_records(records, kind):
    written = compare(TOOL, records[kind], RUN, f"revert_bams {kind}")
    note(records["bam"], written, "input")
    compare(TOOL, records[kind], RUN, f"revert_bams {kind} into SAM text", suffix=".sam")


@pytest.mark.parametrize("order", ["unsorted", "queryname", "coordinate", "duplicate"])
def test_each_sort_order(records, order):
    """The generated records hold paired records flagged neither first nor second, which
    get no mate information and no `MC`: the duplicate order fails in both."""
    args = [*RUN, f"SORT_ORDER={order}"]
    written = compare(TOOL, records["bam"], args, f"revert_bams SORT_ORDER={order}",
                      "fail" if order == "duplicate" else "equal")
    if written.outcome == "equal":
        note(records["bam"], written, order)
    compare(TOOL, records["unsorted"], [*args, "MAX_RECORDS_IN_RAM=20"],
            f"revert_bams unsorted SORT_ORDER={order} in runs of 20",
            "fail" if order == "duplicate" else "equal")


def test_duplicate_order_where_every_mate_has_mc():
    written = compare(TOOL, SOURCES["oxog"], [*RUN, "SORT_ORDER=duplicate"],
                      "oxog SORT_ORDER=duplicate")
    note(SOURCES["oxog"], written, "duplicate")


def test_original_qualities_kept(records):
    compare(TOOL, records["bam"], [*RUN, "RESTORE_ORIGINAL_QUALITIES=false"],
            "revert_bams RESTORE_ORIGINAL_QUALITIES=false")


@pytest.mark.parametrize("version", ["3.0", "3.1"])
def test_cram(records, tmp_path, version):
    cram = tmp_path / f"revert.{version}.cram"
    options = ["--output-fmt-option", f"version={version}"]
    if version == "3.1":
        options += ["--output-fmt-option", "use_tok=0"]
    reference = records["reference"]
    subprocess.run(["samtools", "view", "--no-PG", "-C", *options, "-T", str(reference),
                    "-o", str(cram), str(records["bam"])], check=True)
    compare(TOOL, cram, [*RUN, f"R={reference}"], f"revert_bams CRAM {version}",
            reference=reference)


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    """markdup, duplicate_marking, the foundation corpus and the exome slice hold two
    primary first reads under one name: Picard decides, bammetrics follows."""
    source = source_path(name, corpus, exome)
    written = compare(TOOL, source, RUN, name, "same")
    if written.outcome == "equal":
        note(source, written, "input")


def test_sam_text_and_a_cram_fixture():
    compare(TOOL, FIXTURES / "markdup.sam", RUN, "markdup.sam", "same")
    compare(TOOL, FIXTURES / "insert_size.cram", [*RUN, f"R={REF}"], "insert_size.cram",
            reference=REF)


def skips(source, args, tmp_path, label):
    """Both exit 0 and neither writes OUTPUT, its index or its MD5."""
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        done = run_tool(side, TOOL, [f"I={source}", f"O={directory / 'out.bam'}", *args])
        assert done.returncode == 0, f"{label} {side}: {done.stderr[-800:]}"
        assert list(directory.iterdir()) == [], f"{label}: {side} wrote a file"
    RESULTS.append((f"{TOOL} skip {label}", 0, 0, "both skip the file and write nothing"))


def test_files_picard_skips(records, tmp_path):
    cases = [
        (SOURCES["oxog"], [], "oxog: MC on the first pair"),
        (SOURCES["t1"], ["VALIDATION_STRINGENCY=SILENT"], "t1"),
        (records["bam"], ["VALIDATION_STRINGENCY=SILENT"], "revert_bams"),
        (records["skip"], [], "unpaired records, no OQ, end of file"),
        (FIXTURES / "no_reads.bam", ["MAX_RECORDS_TO_EXAMINE=0"], "no record"),
        (records["skip"], ["MAX_RECORDS_TO_EXAMINE=5"], "five unpaired records of five"),
        (SOURCES["oxog"], ["RESTORE_ORIGINAL_QUALITIES=false", "SORT_ORDER=queryname"],
         "oxog without restoring"),
    ]
    for i, (source, args, label) in enumerate(cases):
        directory = tmp_path / f"case{i}"
        directory.mkdir()
        skips(source, args, directory, label)


def test_files_picard_reverts_after_a_look(records):
    """Fewer records examined than the file holds and no evidence: the file is
    reverted."""
    compare(TOOL, records["skip"], ["MAX_RECORDS_TO_EXAMINE=3"], "unpaired, 3 examined of 5")
    compare(TOOL, records["skip"], ["MAX_RECORDS_TO_EXAMINE=-1"], "unpaired, -1 examined")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=1", "SORT_ORDER=coordinate"],
                                  ["CREATE_INDEX=false", "CREATE_MD5_FILE=false"]])
def test_options(tmp_path, args):
    compare(TOOL, SOURCES["oxog"], [*RUN, f"TMP_DIR={tmp_path}", *args],
            f"oxog {' '.join(args)}")


def test_index_and_md5(records, tmp_path):
    index_and_md5(TOOL, [f"I={records['bam']}", *RUN], tmp_path, "out.bai")


def test_index_and_md5_by_default(records, tmp_path):
    """CREATE_INDEX and CREATE_MD5_FILE default to true: out.bai and out.bam.md5 beside a
    coordinate-sorted BAM, out.bam.md5 alone beside a queryname one."""
    for order, names in (("coordinate", ["out.bai", "out.bam", "out.bam.md5"]),
                         ("queryname", ["out.bam", "out.bam.md5"])):
        for side in ("picard", "bammetrics"):
            directory = tmp_path / f"{side}_{order}"
            directory.mkdir()
            done = run_tool(side, TOOL, [f"I={records['bam']}", f"O={directory / 'out.bam'}",
                                         *RUN, f"SORT_ORDER={order}"])
            assert done.returncode == 0, done.stderr[-800:]
            assert sorted(p.name for p in directory.iterdir()) == names, (side, order)
    RESULTS.append((f"{TOOL} default companions", 2, 2, "the same files beside each BAM"))


ERRORS = {
    "bad_order": ["SORT_ORDER=bogus"],
    "no_room": ["MAX_RECORDS_IN_RAM=0"],
    "assume_sorted": ["ASSUME_SORTED=true"],
}


@pytest.mark.parametrize("case", [*revert_bams.FAILS, "unknown_order"])
def test_failing_records(records, case):
    compare(TOOL, records["fails"][case], ["MAX_RECORDS_TO_EXAMINE=0"], f"error {case}", "fail")


def test_failing_records_while_looking(records):
    """An MC that is no string or no CIGAR fails while Picard looks at the first records."""
    for case in ("mc_integer", "mc_not_a_cigar"):
        compare(TOOL, records["fails"][case], [], f"error {case} while looking", "fail")


def test_an_unknown_input_order_of_one_record(records):
    compare(TOOL, records["unknown_one"], ["MAX_RECORDS_TO_EXAMINE=0"],
            "SO:unknown, one record")


@pytest.mark.parametrize("case", [*ERRORS, "missing_input", "missing_directory"])
def test_errors_as_picard(tmp_path, case):
    source = SOURCES["oxog"]
    args = [*RUN, *ERRORS.get(case, [])]
    if case == "missing_input":
        source = tmp_path / "nope.bam"
    if case == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"O={tmp_path / 'no' / 'out.bam'}",
                                         *args])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    compare(TOOL, source, args, f"error {case}", "fail")


@pytest.mark.parametrize("name", ["refid_far", "mate_refid_far", "refid_negative"])
def test_a_reference_index_outside_the_dictionary(tmp_path_factory, name):
    far = revert_bams.write_far(tmp_path_factory.mktemp("far"))
    compare(TOOL, far[name], RUN, name, "fail")


def test_default_companions_beside_standard_output(records, tmp_path):
    """CREATE_INDEX and CREATE_MD5_FILE default to true: a pipe gets the BAM alone, a
    redirect to a file fails in both."""
    standard_output(TOOL, [f"I={records['bam']}", *RUN], tmp_path, "default companions")


def test_the_step_in_a_chain_equals_picard(records, tmp_path):
    step = ["revert-original-base-qualities-and-add-mate-cigar"]
    cases = [
        (records["bam"], [], []),
        (records["unsorted"], ["--sort-order", "queryname"], ["SORT_ORDER=queryname"]),
        (SOURCES["oxog"], ["--restore-original-qualities", "false"],
         ["RESTORE_ORIGINAL_QUALITIES=false"]),
    ]
    for i, (source, extra, picard) in enumerate(cases):
        directory = tmp_path / f"case{i}"
        directory.mkdir()
        chain_against_picard(TOOL, [*step, *extra], source, [*RUN, *picard], directory)


def test_values_are_distinct():
    assert SEEN["restored"] >= 50, SEEN
    assert SEEN["mc"] >= 200, SEEN
    assert len(SEEN["orders"]) >= 4, SEEN

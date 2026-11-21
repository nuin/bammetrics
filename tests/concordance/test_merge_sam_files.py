"""MergeSamFiles against Picard 3.5.0.

Every case goes through `compare.compare_bam` with `record_bytes=True` (tier B: header
lines other than `@PG` with the `@RG` `DT` cut, `@PG` ID and PP, every record as samtools
prints it, then every record as bytes after inflation: the `bin` htsjdk computes, and the
tags of a record whose `RG` or `PG` was set in htsjdk's order and widths).

Inputs (`merge_bams.py`): the golden fixtures (`tests/data/make_merge_reorder_fixtures.py`); the
foundation corpus, the exome slice and the markdup fixture
split by read group and merged back; the same pieces with every read group and a `@PG`
renamed so that both collide; the pieces sorted by query name, under `SO:unsorted`,
without `@HD`, one by query name and the rest by coordinate; a first piece whose
dictionary lacks the last contig, one with a contig before and after its own, two whose
contigs cross, a first piece whose dictionary is reversed (with and without
MERGE_SEQUENCE_DICTIONARIES); a BAM with no record; one input;
64 inputs; pieces with `@CO` lines and COMMENT; SAM text and CRAM inputs. Every
SORT_ORDER against each input order; `SORT_ORDER=unsorted` takes the inputs in turn in
bammetrics and in `System.identityHashCode` order in Picard, so those cases compare the
header and the records as a multiset (docs/07 decision row); a sort of inputs that are
not in the order hands the sort the records in the same hash order, so records that tie
in queryname order (one name in two read groups of the markdup fixture) are compared up
to the order of the tie. `SORT_ORDER=duplicate` fails in both on a pair without `MC`. ASSUME_SORTED, INTERVALS on
indexed pieces, COMPRESSION_LEVEL, MAX_RECORDS_IN_RAM with spilled runs, TMP_DIR,
USE_THREADING, CREATE_INDEX and CREATE_MD5_FILE, SAM output. Error parity: dictionaries
that differ without MERGE_SEQUENCE_DICTIONARIES, crossed dictionaries with it, INTERVALS
without an index, CREATE_INDEX with a queryname order, a bad SORT_ORDER, `unknown` on two
inputs, MAX_RECORDS_IN_RAM=0 with a sort, ASSUME_SORTED on inputs out of order, a list
file as INPUT, a missing input. The foundation corpus runs at VALIDATION_STRINGENCY=SILENT
(Picard's STRICT refuses some of its records). The distinct-value gate wants 4 output
orders, renamed IDs and a merged dictionary among the cases that compared equal.
"""
import subprocess

import pytest

import merge_bams
from cases import BAMMETRICS, FIXTURES
from compare import compare_bam, sam_view, bam_header_lines
from conftest import RESULTS
from test_add_comments_to_bam import index_and_md5, run_tool

TOOL = "MergeSamFiles"
SILENT = ["VALIDATION_STRINGENCY=SILENT"]

# (output order, renamed IDs, merged dictionary) of every case that compared equal.
SEEN = []


@pytest.fixture(scope="module")
def pieces(tmp_path_factory, corpus, exome):
    directory = tmp_path_factory.mktemp("merge")
    return {
        "foundation": merge_bams.pieces(corpus / "corpus.bam", directory, "foundation",
                                        BAMMETRICS),
        "exome": merge_bams.pieces(exome / "exome.bam", directory, "exome", BAMMETRICS),
        "markdup": merge_bams.pieces(FIXTURES / "markdup.bam", directory, "markdup",
                                     BAMMETRICS),
    }


def note(args, written, merged=False):
    order = next((a.split("=", 1)[1] for a in args if a.split("=")[0] in ("SO", "SORT_ORDER")),
                 "coordinate")
    renamed = any(".1\t" in line or line.endswith(".1") for line in written.header
                  if line.startswith(("@RG", "@PG")))
    SEEN.append((order, renamed, merged))


def merge(inputs, args, label, outcome="equal", suffix=".bam", reference=None, merged=False):
    written = compare_bam(TOOL, inputs, [*SILENT, *args], outcome=outcome, label=label,
                          suffix=suffix, reference=reference, record_bytes=suffix == ".bam")
    if written.outcome == "fail":
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both fail"))
        return written
    note(args, written, merged)
    count = len(written.records)
    RESULTS.append((f"{TOOL} {label}", count, count, "records equal"))
    return written


def unordered(inputs, args, tmp_path, label):
    """`SORT_ORDER=unsorted`: the same header lines and the same records in any order."""
    lines = {}
    for side in ("picard", "bammetrics"):
        out = tmp_path / f"{side}.bam"
        done = run_tool(side, TOOL, [*(f"I={i}" for i in inputs), f"O={out}", *SILENT, *args])
        assert done.returncode == 0, f"{label} {side}: {done.stderr[-800:]}"
        lines[side] = (bam_header_lines(sam_view(out, True)), sorted(sam_view(out)))
    assert lines["picard"][0] == lines["bammetrics"][0], label
    assert lines["picard"][1] == lines["bammetrics"][1], label
    count = len(lines["picard"][1])
    SEEN.append(("unsorted", False, False))
    RESULTS.append((f"{TOOL} {label}", count, count, "records equal as a multiset"))


@pytest.mark.parametrize("source", ["foundation", "exome", "markdup"])
@pytest.mark.parametrize("kind", ["by_rg", "colliding", "with_comments", "empty", "many"])
def test_pieces_merge_back(pieces, source, kind):
    merge(pieces[source][kind], [], f"{source} {kind}")


@pytest.mark.parametrize("source", ["foundation", "markdup"])
def test_one_input(pieces, source):
    merge(pieces[source]["by_rg"][:1], [], f"{source} one input")


ORDERS = ["coordinate", "queryname", "duplicate"]
KINDS = ["by_rg", "queryname", "mixed", "unsorted_header", "no_hd", "colliding"]


def up_to_ties(inputs, args, tmp_path, label):
    """A sort of records that tie in queryname order (one name in two read groups): Picard
    hands them to its sort in `System.identityHashCode` order, bammetrics in input order
    (docs/07 decision row). The same header lines, the same records, and the same
    sequence of (name, flag) keys."""
    lines = {}
    for side in ("picard", "bammetrics"):
        out = tmp_path / f"{side}.bam"
        done = run_tool(side, TOOL, [*(f"I={i}" for i in inputs), f"O={out}", *SILENT, *args])
        assert done.returncode == 0, f"{label} {side}: {done.stderr[-800:]}"
        lines[side] = (bam_header_lines(sam_view(out, True)), sam_view(out))
    (p_header, p_records), (o_header, o_records) = lines["picard"], lines["bammetrics"]
    assert p_header == o_header, label
    assert sorted(p_records) == sorted(o_records), label
    keys = [[tuple(r.split("\t")[:2]) for r in records] for records in (p_records, o_records)]
    assert keys[0] == keys[1], label
    SEEN.append(("queryname", False, False))
    count = len(p_records)
    RESULTS.append((f"{TOOL} {label}", count, count, "records equal, ties in any order"))


def outcome_of(order):
    """Picard's duplicate comparator fails on a pair without `MC`; bammetrics too."""
    return "same" if order == "duplicate" else "equal"


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("kind", KINDS)
def test_sort_orders(pieces, tmp_path, order, kind):
    label = f"markdup {kind} SO={order}"
    if order == "queryname" and kind != "queryname":
        up_to_ties(pieces["markdup"][kind], [f"SO={order}"], tmp_path, label)
        return
    merge(pieces["markdup"][kind], [f"SO={order}"], label, outcome=outcome_of(order))


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("kind", ["queryname", "mixed"])
def test_sort_orders_on_the_corpus(pieces, order, kind):
    merge(pieces["foundation"][kind], [f"SO={order}"], f"foundation {kind} SO={order}",
          outcome=outcome_of(order))


@pytest.mark.parametrize("kind", ["by_rg", "queryname", "mixed"])
def test_unsorted_takes_every_record(pieces, tmp_path, kind):
    unordered(pieces["foundation"][kind], ["SO=unsorted"], tmp_path, f"foundation {kind} unsorted")


def test_unknown_order(pieces):
    merge(pieces["markdup"]["by_rg"], ["SO=unknown"], "markdup SO=unknown", outcome="fail")
    merge(pieces["markdup"]["by_rg"][:1], ["SO=unknown"], "markdup one input SO=unknown",
          outcome="same")


def test_assume_sorted(pieces):
    merge(pieces["markdup"]["unsorted_header"], ["AS=true"], "markdup unsorted header AS=true")
    merge(pieces["markdup"]["mixed"], ["AS=true"], "markdup mixed AS=true", outcome="same")
    merge(pieces["markdup"]["queryname"], ["AS=true", "SO=queryname"],
          "markdup queryname AS=true SO=queryname")


@pytest.mark.parametrize("source", ["foundation", "exome"])
@pytest.mark.parametrize("kind", ["subset_dict", "extra_contig"])
def test_dictionaries_merge(pieces, source, kind):
    merge(pieces[source][kind], [], f"{source} {kind} without MSD", outcome="fail")
    merge(pieces[source][kind], ["MERGE_SEQUENCE_DICTIONARIES=true"], f"{source} {kind} MSD",
          merged=True)
    merge(pieces[source][kind], ["MERGE_SEQUENCE_DICTIONARIES=true", "SO=queryname"],
          f"{source} {kind} MSD SO=queryname", merged=True)


@pytest.mark.parametrize("kind", ["crossed", "reordered_dict"])
def test_crossed_dictionaries_fail(pieces, kind):
    merge(pieces["foundation"][kind], ["MERGE_SEQUENCE_DICTIONARIES=true"],
          f"foundation {kind} MSD", outcome="fail")


def test_comments(pieces):
    merge(pieces["markdup"]["with_comments"], ["COMMENT=one", "CO=two words", "CO=@CO\tthree"],
          "markdup COMMENT")


@pytest.fixture(scope="module")
def indexed(pieces):
    for piece in pieces["foundation"]["by_rg"]:
        subprocess.run(["samtools", "index", str(piece)], check=True)
    return pieces["foundation"]["by_rg"]


def test_intervals_on_indexed_pieces(indexed, corpus):
    for name in ("target.interval_list", "bait.interval_list"):
        merge(indexed, [f"INTERVALS={corpus / name}"], f"foundation INTERVALS={name}")


def test_intervals_without_an_index(pieces, corpus):
    merge(pieces["markdup"]["by_rg"], [f"INTERVALS={corpus / 'target.interval_list'}"],
          "markdup INTERVALS without index", outcome="fail")


def test_sam_text_and_cram_inputs(pieces, corpus, tmp_path):
    ref = corpus / "ref.fa"
    first, second = pieces["foundation"]["by_rg"][:2]
    sam = merge_bams.write(merge_bams.text(first), tmp_path / "first.sam")
    cram = tmp_path / "second.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(ref), "-o", str(cram),
                    str(second)], check=True)
    merge([sam, cram, *pieces["foundation"]["by_rg"][2:]], [f"R={ref}"],
          "foundation SAM text and CRAM")
    merge(pieces["markdup"]["by_rg"], [], "markdup to SAM text", suffix=".sam")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["USE_THREADING=true"], ["MAX_RECORDS_IN_RAM=50"],
                                  ["MAX_RECORDS_IN_RAM=7", "SO=queryname"]],
                         ids=lambda a: " ".join(a))
def test_options(pieces, tmp_path, args):
    merge(pieces["foundation"]["mixed"], [*args, f"TMP_DIR={tmp_path}"],
          f"foundation mixed {' '.join(args)}")


def test_index_and_md5(pieces, tmp_path):
    index_and_md5(TOOL, [f"I={p}" for p in pieces["markdup"]["by_rg"]], tmp_path, "out.bai")


ERRORS = {
    "bad_order": ["SO=sideways"],
    "index_not_coordinate": ["SO=queryname", "CREATE_INDEX=true"],
    "no_room": ["SO=queryname", "MAX_RECORDS_IN_RAM=0"],
}


@pytest.mark.parametrize("case", sorted(ERRORS))
def test_errors_as_picard(pieces, case):
    merge(pieces["markdup"]["by_rg"], ERRORS[case], f"error {case}", outcome="fail")


def test_input_errors_as_picard(pieces, tmp_path):
    listing = tmp_path / "inputs.list"
    listing.write_text("".join(f"{p}\n" for p in pieces["markdup"]["by_rg"]))
    merge([listing], [], "error a list file", outcome="fail")
    merge([pieces["markdup"]["by_rg"][0], tmp_path / "missing.bam"], [], "error missing input",
          outcome="fail")


def test_the_golden_fixtures(tmp_path):
    """The inputs of the golden cases: the markdup and artifact BAMs (`rg1` collides),
    three fixture BAMs by query name, and the SAM text of another dictionary."""
    md, art = FIXTURES / "markdup.bam", FIXTURES / "artifact.bam"
    merge([md, art], ["COMMENT=lane 1"], "fixtures markdup artifact COMMENT")
    merge([md, FIXTURES / "oxog.bam", FIXTURES / "insert_size.bam"], ["SO=queryname"],
          "fixtures three inputs SO=queryname")
    other = FIXTURES / "merge_other_dict.sam"
    merge([md, other], ["MERGE_SEQUENCE_DICTIONARIES=true"], "fixtures other dictionary MSD",
          merged=True)
    merge([md, other], ["MERGE_SEQUENCE_DICTIONARIES=true"], "fixtures other dictionary SAM out",
          suffix=".sam", merged=True)
    merge([md, other], [], "fixtures other dictionary", outcome="fail")


def test_values_are_distinct():
    orders = {order for order, _, _ in SEEN}
    assert len(orders) >= 4, orders
    assert any(renamed for _, renamed, _ in SEEN)
    assert any(merged for _, _, merged in SEEN)

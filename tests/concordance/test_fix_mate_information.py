"""FixMateInformation against Picard 3.5.0.

Picard reads every INPUT (several merged by htsjdk's `SamFileHeaderMerger` and
`MergingSamRecordIterator`), sorts the records by name unless every input says
`SO:queryname` or ASSUME_SORTED is true, sets the mate fields of each name's records
(`SamPairUtil.SetMateInfoIterator`) and writes them through a writer presorted for
`queryname` and sorting for every other order. Every case goes through
`compare.compare_bam` (tier B: header lines, `@PG` ID and PP, every record as samtools
prints it); a BAM case also compares every record as bytes after inflation
(`record_bytes=True`: the `bin` field, tag order and integer widths).

Inputs: `mate_filter_bams.py` (the foundation corpus with one primary record per end
and mate fields damaged on a third of its records, orphans; in coordinate, htsjdk
queryname and samtools natural order; split by read group into three files whose read
group and program IDs collide), the `revert_bams.py` records (coordinate, queryname,
unsorted, SAM text, CRAM 3.0 and 3.1), the fixture BAMs, t1 to t3, the foundation corpus
and the exome slice as they are (two primary first reads under one name: Picard
decides). Options: each SORT_ORDER, ASSUME_SORTED, ADD_MATE_CIGAR, IGNORE_MISSING_MATES,
MAX_RECORDS_IN_RAM (runs of 20 and of 1), TMP_DIR, COMPRESSION_LEVEL, CREATE_INDEX,
CREATE_MD5_FILE, SAM text output, two and three INPUTs, a BAM without MD and NM merged
with a CRAM (only the CRAM's records get them). Without OUTPUT: the input is
rewritten in place and equals Picard's in-place file; bammetrics keeps the mode, a link,
and writes `x.bai` and `x.bam.md5` where Picard leaves them under its temporary name
(docs/07 decision rows). Error parity: two primary first reads, a missing mate under
IGNORE_MISSING_MATES=false, `SO:unknown` on two records, CREATE_INDEX with another order
than coordinate, two INPUTs without OUTPUT, a bad sort order, `MAX_RECORDS_IN_RAM=0`, a
missing input, a missing output directory, a reference index outside the dictionary. The
step in a chain against Picard's file. The distinct-value gate wants 200 MC values set,
50 records whose mate start changed, 4 output orders and 3 merged inputs.
"""
import shutil
import subprocess

import pytest

import mate_filter_bams
import revert_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_bam
from conftest import RESULTS
from test_add_comments_to_bam import (SOURCES, chain_against_picard, compare, index_and_md5,
                                      run_tool, source_path)

TOOL = "FixMateInformation"
REF = FIXTURES / "ref_100k.fa"
SILENT = ["VALIDATION_STRINGENCY=SILENT"]

SEEN = {"mc": 0, "moved": 0, "orders": set(), "inputs": 0}


@pytest.fixture(scope="module")
def bams(corpus, tmp_path_factory):
    return mate_filter_bams.write(corpus, tmp_path_factory.mktemp("mate_filter"))


@pytest.fixture(scope="module")
def records(tmp_path_factory):
    return revert_bams.write(tmp_path_factory.mktemp("revert"))


def view(path):
    return subprocess.run(["samtools", "view", str(path)], check=True, capture_output=True,
                          text=True).stdout.splitlines()


def note(source, written, order):
    """Count the MC values written and the records whose mate start moved."""
    SEEN["mc"] += sum("\tMC:Z:" in line for line in written.records)
    before = {}
    for line in view(source):
        fields = line.split("\t")
        before.setdefault((fields[0], fields[1], fields[3]), fields[7])
    for line in written.records:
        fields = line.split("\t")
        old = before.get((fields[0], fields[1], fields[3]))
        SEEN["moved"] += old is not None and old != fields[7]
    SEEN["orders"].add(order)


def merged(inputs, args, label, outcome="equal", suffix=".bam"):
    written = compare_bam(TOOL, inputs, args, args, outcome=outcome, label=label,
                          suffix=suffix, record_bytes=suffix == ".bam")
    RESULTS.append((f"{TOOL} {label}", len(written.records), len(written.records),
                    "both fail" if written.outcome == "fail" else "records equal"))
    return written


@pytest.mark.parametrize("kind", ["damaged", "damaged.queryname", "damaged.natural"])
def test_the_damaged_corpus(bams, kind):
    source = bams / f"{kind}.bam"
    outcome = "same" if kind == "damaged.natural" else "equal"
    written = compare(TOOL, source, SILENT, kind, outcome)
    if written.outcome == "equal":
        note(source, written, "input")
    compare(TOOL, source, SILENT, f"{kind} into SAM text", outcome, suffix=".sam")


@pytest.mark.parametrize("order", ["unsorted", "queryname", "coordinate", "duplicate"])
def test_each_sort_order(bams, order):
    """The duplicate order reads MC of every paired record with a mapped mate; Picard
    decides on the damaged records."""
    args = [*SILENT, f"SORT_ORDER={order}"]
    for kind in ("damaged", "damaged.queryname"):
        written = compare(TOOL, bams / f"{kind}.bam", args, f"{kind} SORT_ORDER={order}",
                          "same" if order == "duplicate" else "equal")
        if written.outcome == "equal":
            note(bams / f"{kind}.bam", written, order)


@pytest.mark.parametrize("args", [["ADD_MATE_CIGAR=false"], ["IGNORE_MISSING_MATES=false"],
                                  ["ASSUME_SORTED=true"],
                                  ["ASSUME_SORTED=true", "SORT_ORDER=queryname"],
                                  ["MAX_RECORDS_IN_RAM=20"], ["MAX_RECORDS_IN_RAM=1"],
                                  ["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"]])
def test_options(bams, tmp_path, args):
    """IGNORE_MISSING_MATES=false fails on the orphans; ASSUME_SORTED on a coordinate
    input sets the mates of adjacent records alone, and with SORT_ORDER=queryname
    Picard's presorted writer fails on the coordinate order."""
    run = [*SILENT, f"TMP_DIR={tmp_path}", *args]
    outcome = "same" if args[0] in ("IGNORE_MISSING_MATES=false", "ASSUME_SORTED=true") \
        else "equal"
    for kind in ("damaged", "damaged.queryname"):
        compare(TOOL, bams / f"{kind}.bam", run, f"{kind} {' '.join(args)}", outcome)


@pytest.mark.parametrize("order", ["", "queryname"])
def test_split_by_read_group_and_merged(bams, order):
    """Three inputs whose read group IDs `g0`, `g1`, ... and `@PG` IDs collide: htsjdk
    renames them; coordinate inputs are taken one after the other, queryname inputs
    merged by name."""
    for suffix in ("", ".queryname"):
        inputs = [bams / f"split_{i}{suffix}.bam" for i in range(3)]
        args = [*SILENT, *([f"SORT_ORDER={order}"] if order else [])]
        merged(inputs, args, f"split{suffix} {order or 'input order'}")
        merged(inputs[:2], args, f"two of split{suffix} {order or 'input order'}")
    SEEN["inputs"] = max(SEEN["inputs"], 3)


@pytest.mark.parametrize("kind", ["bam", "queryname", "unsorted", "sam"])
def test_the_revert_records(records, kind):
    for args in ([], ["SORT_ORDER=coordinate"], ["ADD_MATE_CIGAR=false"]):
        written = compare(TOOL, records[kind], [*SILENT, *args], f"revert_bams {kind} {args}",
                          "same")
        if written.outcome == "equal":
            note(records["bam"], written, "input")


@pytest.mark.parametrize("version", ["3.0", "3.1"])
def test_cram(bams, corpus, tmp_path, version):
    cram = tmp_path / f"damaged.{version}.cram"
    options = ["--output-fmt-option", f"version={version}"]
    if version == "3.1":
        options += ["--output-fmt-option", "use_tok=0"]
    reference = corpus / "ref.fa"
    subprocess.run(["samtools", "view", "--no-PG", "-C", *options, "-T", str(reference),
                    "-o", str(cram), str(bams / "damaged.bam")], check=True)
    compare(TOOL, cram, [*SILENT, f"R={reference}"], f"damaged CRAM {version}",
            reference=reference)


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    """markdup, the foundation corpus and the exome slice hold two primary first reads
    under one name: Picard decides, bammetrics follows."""
    source = source_path(name, corpus, exome)
    written = compare(TOOL, source, SILENT, name, "same")
    if written.outcome == "equal":
        note(source, written, "input")


def test_a_bam_and_a_cram_merged():
    """htsjdk rebuilds MD and NM as it decodes a CRAM: the 10 records of the CRAM get them,
    the 23 of a BAM without them stay without them. The two share no read group ID: htsjdk
    fails with a NullPointerException on a CRAM input whose read group ID collides."""
    bam, cram = FIXTURES / "mate_filter.bam", FIXTURES / "insert_size.cram"
    for inputs in ([bam, cram], [cram, bam]):
        for order in ("", "queryname"):
            args = [*SILENT, f"R={REF}", *([f"SORT_ORDER={order}"] if order else [])]
            written = merged(inputs, args, f"{inputs[0].name} and {inputs[1].name} "
                             f"{order or 'input order'}")
            assert sum("\tNM:i:" in line for line in written.records) == 10, written.records
            assert len(written.records) == 33


def test_sam_text_and_a_cram_fixture():
    compare(TOOL, FIXTURES / "markdup.sam", SILENT, "markdup.sam", "same")
    compare(TOOL, FIXTURES / "insert_size.cram", [*SILENT, f"R={REF}"], "insert_size.cram",
            reference=REF)


def test_index_and_md5(bams, tmp_path):
    index_and_md5(TOOL, [f"I={bams / 'damaged.bam'}", *SILENT], tmp_path, "out.bai")


def in_place(side, source, directory, args):
    target = directory / "fix.bam"
    shutil.copy(source, target)
    done = run_tool(side, TOOL, [f"I={target}", *SILENT, *args])
    assert done.returncode == 0, f"{side}: {done.stderr[-800:]}"
    return target


def test_in_place_equals_picard(bams, tmp_path):
    """Without OUTPUT both rewrite the input; the records are equal. Picard leaves the
    index and MD5 under its temporary name and warns; bammetrics names them after the
    input (a decision row)."""
    files = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        target = in_place(side, bams / "damaged.bam", directory, ["CREATE_INDEX=true"])
        files[side] = sorted(p.name for p in directory.iterdir())
        assert view(target), side
    assert view(tmp_path / "picard" / "fix.bam") == view(tmp_path / "bammetrics" / "fix.bam")
    assert files["bammetrics"] == ["fix.bai", "fix.bam"], files
    assert "fix.bam" in files["picard"], files
    RESULTS.append((f"{TOOL} in place", 1, 1, "the rewritten input equals Picard's"))


def test_in_place_keeps_mode_and_link(bams, tmp_path):
    target = tmp_path / "real.bam"
    shutil.copy(bams / "damaged.bam", target)
    target.chmod(0o640)
    link = tmp_path / "link.bam"
    link.symlink_to(target.name)
    done = run_tool("bammetrics", TOOL, [f"I={link}", *SILENT])
    assert done.returncode == 0, done.stderr
    assert link.is_symlink() and oct(target.stat().st_mode & 0o777) == oct(0o640)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["link.bam", "real.bam"]


def test_in_place_refuses_sam_text(tmp_path):
    """Picard writes BAM bytes over a .sam input; bammetrics refuses and leaves it."""
    target = tmp_path / "x.sam"
    shutil.copy(FIXTURES / "markdup.sam", target)
    before = target.read_bytes()
    done = run_tool("bammetrics", TOOL, [f"I={target}"])
    assert done.returncode == 1 and "give OUTPUT" in done.stderr, done.stderr
    verb = subprocess.run([str(BAMMETRICS), "fix-mate-information", str(target)],
                          capture_output=True, text=True)
    assert verb.returncode == 2, verb.stderr
    assert target.read_bytes() == before


ERRORS = {
    "bad_order": ["SORT_ORDER=bogus"],
    "no_room": ["MAX_RECORDS_IN_RAM=0"],
    "index_not_coordinate": ["CREATE_INDEX=true", "SORT_ORDER=queryname"],
    "missing_mates": ["IGNORE_MISSING_MATES=false"],
}


@pytest.mark.parametrize("case", [*ERRORS, "missing_input", "missing_directory",
                                  "two_inputs_in_place"])
def test_errors_as_picard(bams, tmp_path, case):
    source = bams / "damaged.bam"
    args = [*SILENT, *ERRORS.get(case, [])]
    if case == "missing_input":
        source = tmp_path / "nope.bam"
    if case in ("missing_directory", "two_inputs_in_place"):
        for side in ("picard", "bammetrics"):
            output = [] if case == "two_inputs_in_place" else \
                [f"O={tmp_path / 'no' / 'out.bam'}"]
            inputs = [f"I={source}"] * (2 if case == "two_inputs_in_place" else 1)
            done = run_tool(side, TOOL, [*inputs, *output, *args])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error {case}", 0, 0, "both fail"))
        return
    compare(TOOL, source, args, f"error {case}", "fail")


@pytest.mark.parametrize("case", ["two_first", "unknown_order"])
def test_failing_records(records, case):
    compare(TOOL, records["fails"][case], SILENT, f"error {case}", "fail")


@pytest.mark.parametrize("name", ["refid_far", "mate_refid_far", "refid_negative"])
def test_a_reference_index_outside_the_dictionary(tmp_path_factory, name):
    far = revert_bams.write_far(tmp_path_factory.mktemp("far"))
    compare(TOOL, far[name], SILENT, name, "fail")


def test_the_step_in_a_chain_equals_picard(bams, tmp_path):
    step = ["fix-mate-information"]
    cases = [
        (bams / "damaged.bam", [], []),
        (bams / "damaged.queryname.bam", ["--sort-order", "coordinate"],
         ["SORT_ORDER=coordinate"]),
        (bams / "damaged.bam", ["--add-mate-cigar", "false"], ["ADD_MATE_CIGAR=false"]),
        (CORPUS / "t1.bam", ["--sort-order", "queryname"], ["SORT_ORDER=queryname"]),
    ]
    for i, (source, extra, picard) in enumerate(cases):
        directory = tmp_path / f"case{i}"
        directory.mkdir()
        chain_against_picard(TOOL, [*step, *extra], source, [*SILENT, *picard], directory)


def test_values_are_distinct():
    assert SEEN["mc"] >= 200, SEEN
    assert SEEN["moved"] >= 50, SEEN
    assert len(SEEN["orders"]) >= 4, SEEN
    assert SEEN["inputs"] >= 3, SEEN

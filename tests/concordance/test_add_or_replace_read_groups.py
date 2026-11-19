"""AddOrReplaceReadGroups against Picard 3.5.0.

Picard replaces every `@RG` line by the read group its options give, sets `RG` on every
record and writes through a writer made with `presorted` true when the output order is the
input's (the order checked) and false otherwise (the records sorted at close). Every case
goes through `compare.compare_bam` with `record_bytes=True` (tier B: header lines, `@PG`
ID and PP, every record as samtools prints it, then every record as bytes after inflation:
every record is encoded again, tags in htsjdk's order and widths).

`compare_bam` cuts the `DT` field of `@RG` lines (htsjdk writes it in the JVM's time zone).
Here Picard runs under `TZ=UTC`, the zone bammetrics writes, and `same_read_group`
compares the whole `@RG` line, `DT` included, byte for byte.

Inputs: the fixture BAMs, `markdup.sam`, `insert_size.cram`, the foundation corpus (two
samples, two libraries, three read groups), t1 to t3, the exome slice and the
`revert_bams.py` records (coordinate, queryname and unsorted). Options: every read group
field, `RGDT` in 14 forms, `RGID=null`, each `SORT_ORDER` against each input order,
`MAX_RECORDS_IN_RAM`, `TMP_DIR`, COMPRESSION_LEVEL, CREATE_INDEX (`out.bai`),
CREATE_MD5_FILE, SAM text output. Error parity: a required option missing or `null`, a
value with a tab, a bad date, a bad integer, a bad sort order, `MAX_RECORDS_IN_RAM=0` with
a sort, `SORT_ORDER=unknown` on two records, records out of the header's order,
ASSUME_SORTED, a missing input, a missing output directory, a reference or mate reference
index outside the dictionary (`revert_bams.write_far`, at STRICT and SILENT), an index or
MD5 beside standard output. The step in a chain, streaming and sorting, against Picard's
file, `SORT_ORDER=coordinate` on the foundation corpus and t1 among them (samtools sorted
both, with its own order of tied records, which a presorted writer keeps). The distinct-value gate wants 3 read group IDs, 4
output orders and 10 `DT` values.
"""
import subprocess

import pytest

import revert_bams
from cases import BAMMETRICS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env
from test_add_comments_to_bam import (SOURCES, chain_against_picard, compare, index_and_md5,
                                      run_tool, source_path)

TOOL = "AddOrReplaceReadGroups"
REF = FIXTURES / "ref_100k.fa"
RG = ["RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1"]
SILENT = ["VALIDATION_STRINGENCY=SILENT"]

# (read group ID, output order, DT) of every case that compared equal: the gate's input.
SEEN = []


@pytest.fixture(autouse=True)
def utc(monkeypatch):
    """Picard writes `RGDT` in the JVM's time zone; bammetrics writes UTC."""
    monkeypatch.setenv("TZ", "UTC")


@pytest.fixture(scope="module")
def records(tmp_path_factory):
    return revert_bams.write(tmp_path_factory.mktemp("revert"))


def note(args, order="input"):
    rgid = next((a.split("=", 1)[1] for a in args if a.startswith("RGID=")), "1")
    dt = next((a.split("=", 1)[1] for a in args if a.startswith("RGDT=")), None)
    SEEN.append((rgid, order, dt))


def header_text(path):
    return subprocess.run(["samtools", "view", "--no-PG", "-H", str(path)], check=True,
                          capture_output=True, text=True).stdout


def same_read_group(source, args, tmp_path, label):
    """The `@RG` lines of both outputs, `DT` included, byte for byte."""
    lines = {}
    for side in ("picard", "bammetrics"):
        out = tmp_path / f"{side}.rg.bam"
        done = run_tool(side, TOOL, [f"I={source}", f"O={out}", *args])
        assert done.returncode == 0, f"{label} {side}: {done.stderr[-800:]}"
        lines[side] = [line for line in header_text(out).splitlines() if line.startswith("@RG")]
    assert lines["picard"] == lines["bammetrics"], f"{label}: {lines}"
    RESULTS.append((f"{TOOL} @RG {label}", 1, 1, lines["picard"][0]))
    return lines["picard"][0]


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    source = source_path(name, corpus, exome)
    compare(TOOL, source, [*RG, *SILENT], name)
    note(RG)


def test_sam_text_and_a_cram():
    compare(TOOL, FIXTURES / "markdup.sam", RG, "markdup.sam")
    compare(TOOL, FIXTURES / "markdup.bam", RG, "markdup.bam into SAM text", suffix=".sam")
    compare(TOOL, FIXTURES / "insert_size.cram", [*RG, f"R={REF}"], "insert_size.cram",
            reference=REF)


def test_the_generated_records(records):
    for kind in ("bam", "queryname", "unsorted"):
        compare(TOOL, records[kind], [*RG, *SILENT], f"revert_bams {kind}")
    compare(TOOL, records["sam"], [*RG, *SILENT], "revert_bams SAM text")


FIELDS = [
    ["RGID=7"],
    ["RGID=two words", "RGCN=Broad Institute", "RGDS=a description, with commas"],
    ["RGID=rg.1", "RGPI=350", "RGPG=bwa", "RGPM=HiSeq2500", "RGKS=TCAG", "RGFO=TACG"],
    ["RGPI=+5"],
    ["RGPI=-7"],
    ["RGCN=null", "RGDS=null"],
    ["RGCN=c", "RGDS=d", "RGDT=2020-01-02", "RGPI=1", "RGPG=p", "RGPM=m", "RGKS=k", "RGFO=f"],
]


@pytest.mark.parametrize("fields", FIELDS, ids=[" ".join(f) for f in FIELDS])
def test_every_field(tmp_path, fields):
    args = [*RG, *fields]
    compare(TOOL, SOURCES["markdup"], args, f"markdup {' '.join(fields)}")
    same_read_group(SOURCES["markdup"], args, tmp_path, " ".join(fields))
    note(args)


DATES = ["2020-01-02", "2020-01-02T03:04:05", "2020-01-02T03:04:05Z",
         "2020-01-02T03:04:05.123+05:30", "2020-01", "2020", "2020-01-02T03:04",
         "2020-01-02T03:04:05-0700", "2020-13-45", "1969-12-31T23:59:59Z",
         "2020-01-02T03:04:05.5Z", "2020-01-02T25:61:61", "2020-01-02T03:04+05",
         "1999-12-31T23:59:59.999-12:00", "1969-12-31T23:59:59.500Z",
         "1969-12-31T23:58:58.999Z"]


@pytest.mark.parametrize("date", DATES)
def test_run_dates(tmp_path, date):
    args = [*RG, f"RGDT={date}"]
    same_read_group(FIXTURES / "insert_size.bam", args, tmp_path, f"RGDT={date}")
    note(args)


def test_null_id_removes_each_rg_tag():
    """`RGID=null` is Java's null: `ID:null` in the header, no `RG` on any record."""
    written = compare(TOOL, SOURCES["oxog"], [*RG, "RGID=null"], "oxog RGID=null")
    assert all("\tRG:Z:" not in line for line in written.records)
    assert any(line.startswith("@RG\tID:null") for line in written.header)
    note(["RGID=null"])


ORDERS = ["unsorted", "queryname", "coordinate", "duplicate"]


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("name", ["markdup", "duplicate_marking", "oxog"])
def test_sort_orders(name, order):
    """Each order against a coordinate input and a queryname one. The duplicate order
    needs `MC` on every paired record with a mapped mate: Picard decides (markdup has a
    mate without it, and both fail)."""
    args = [*RG, f"SORT_ORDER={order}"]
    outcome = "same" if order == "duplicate" else "equal"
    written = compare(TOOL, SOURCES[name], args, f"{name} SORT_ORDER={order}", outcome)
    compare(TOOL, SOURCES[name], [*args, "MAX_RECORDS_IN_RAM=5"],
            f"{name} SORT_ORDER={order} in runs of 5", outcome)
    if written.outcome == "equal":
        note(args, order)


@pytest.mark.parametrize("order", ORDERS)
def test_sort_orders_on_generated_records(records, order):
    """The generated records hold pairs without `MC`, so the duplicate order fails in
    both. Runs of 20 keep the merge to one level: bammetrics merges 64 runs into one as
    they spill, which orders records whose queryname order is not total (paired records
    flagged neither first nor second) differently from htsjdk's one merge of every run."""
    for kind in ("bam", "unsorted"):
        args = [*RG, *SILENT, f"SORT_ORDER={order}", "MAX_RECORDS_IN_RAM=20"]
        compare(TOOL, records[kind], args, f"revert_bams {kind} SORT_ORDER={order}",
                "same" if order == "duplicate" else "equal")


def test_an_unknown_order(records, tmp_path):
    """`SORT_ORDER=unknown` over a coordinate input: htsjdk sorts with no comparator and
    fails on the second record; a file of one record is written."""
    compare(TOOL, SOURCES["markdup"], [*RG, "SORT_ORDER=unknown"], "markdup SORT_ORDER=unknown",
            "fail")
    one = tmp_path / "one.sam"
    text = header_text(SOURCES["markdup"])
    first = subprocess.run(["samtools", "view", str(SOURCES["markdup"])], check=True,
                           capture_output=True, text=True).stdout.splitlines()[0]
    one.write_text(text + first + "\n")
    compare(TOOL, one, [*RG, "SORT_ORDER=unknown"], "one record SORT_ORDER=unknown")
    compare(TOOL, records["unknown_one"], RG, "SO:unknown input, one record")
    note(["SORT_ORDER=unknown"], "unknown")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=0"], [f"R={REF}"],
                                  ["SORT_ORDER=queryname", "MAX_RECORDS_IN_RAM=1"],
                                  ["SORT_ORDER=coordinate"]])
def test_options(tmp_path, args):
    compare(TOOL, SOURCES["markdup"], [*RG, f"TMP_DIR={tmp_path}", *args],
            f"markdup {' '.join(args)}")


def test_index_and_md5(tmp_path):
    index_and_md5(TOOL, [f"I={SOURCES['markdup']}", *RG], tmp_path, "out.bai")


ERRORS = {
    "no_rglb": [a for a in RG if not a.startswith("RGLB")],
    "rglb_null": [*RG[1:], "RGLB=null"],
    "tab": [*RG, "RGCN=a\tb"],
    "non_ascii": [*RG, "RGSM=é"],
    "bad_date": [*RG, "RGDT=garbage"],
    "date_with_space": [*RG, "RGDT=2020-01-02 03:04"],
    "bad_integer": [*RG, "RGPI=abc"],
    "bad_order": [*RG, "SORT_ORDER=bogus"],
    "no_room": [*RG, "SORT_ORDER=queryname", "MAX_RECORDS_IN_RAM=0"],
    "assume_sorted": [*RG, "ASSUME_SORTED=true"],
}


@pytest.mark.parametrize("case", [*ERRORS, "out_of_order", "missing_input",
                                  "missing_directory"])
def test_errors_as_picard(tmp_path, case):
    source = SOURCES["markdup"]
    args = ERRORS.get(case, RG)
    if case == "out_of_order":
        # A coordinate header over records in another order: the presorted writer fails.
        text = subprocess.run(["samtools", "view", "-h", str(source)], check=True,
                              capture_output=True, text=True).stdout.splitlines(True)
        head = [line for line in text if line.startswith("@")]
        body = [line for line in text if not line.startswith("@")]
        source = tmp_path / "shuffled.sam"
        source.write_text("".join(head + body[::-1]))
    elif case == "missing_input":
        source = tmp_path / "nope.bam"
    if case == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"O={tmp_path / 'no' / 'out.bam'}", *RG])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    compare(TOOL, source, args, f"error {case}", "fail")


@pytest.fixture(scope="module")
def far(tmp_path_factory):
    return revert_bams.write_far(tmp_path_factory.mktemp("far"))


@pytest.mark.parametrize("stringency", ["STRICT", "SILENT"])
@pytest.mark.parametrize("name", ["refid_far", "mate_refid_far", "refid_negative"])
def test_a_reference_index_outside_the_dictionary(far, name, stringency):
    """htsjdk's writer fails on such a record whatever the stringency; bammetrics fails
    and leaves no output."""
    compare(TOOL, far[name], [*RG, f"VALIDATION_STRINGENCY={stringency}"],
            f"{name} {stringency}", "fail")


def run_bytes(tool, side, args, cwd, stdout=subprocess.PIPE):
    """One tool through Picard or bammetrics, standard output kept as bytes."""
    if side == "picard":
        env = picard_env()
        env.setdefault("JAVA_OPTS", "-Xmx8g")
        command = ["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"]
    else:
        env, command = None, [str(BAMMETRICS), tool, *args]
    return subprocess.run(command, env=env, stdout=stdout, stderr=subprocess.PIPE, cwd=cwd)


def standard_output(tool, args, tmp_path, label):
    """`O=/dev/stdout` with `args` (an index or MD5 asked). On a pipe htsjdk writes
    neither, warns, and both tools write the same records and exit 0. Redirected to a
    regular file, Picard fails opening `/dev/stdout.bai` or `.md5`; bammetrics refuses
    first and writes no byte."""
    texts = {}
    for side in ("picard", "bammetrics"):
        done = run_bytes(tool, side, ["O=/dev/stdout", *args], tmp_path)
        assert done.returncode == 0, f"{label} {side}: {done.stderr[-800:]}"
        bam = tmp_path / f"{side}.piped.bam"
        bam.write_bytes(done.stdout)
        texts[side] = subprocess.run(["samtools", "view", "--no-PG", str(bam)], check=True,
                                     capture_output=True, text=True).stdout
        bam.unlink()
    assert texts["picard"] == texts["bammetrics"], label
    assert texts["picard"], label
    for side in ("picard", "bammetrics"):
        redirected = tmp_path / f"{side}.redirected"
        with redirected.open("wb") as out:
            done = run_bytes(tool, side, ["O=/dev/stdout", *args], tmp_path, stdout=out)
        assert done.returncode != 0, f"{label} {side} redirected"
        if side == "bammetrics":
            assert redirected.stat().st_size == 0, f"{label}: bytes before the refusal"
        redirected.unlink()
    assert list(tmp_path.iterdir()) == [], label
    RESULTS.append((f"{tool} O=/dev/stdout {label}", 2, 2,
                    "a pipe gets the same records, a redirect fails in both"))


@pytest.mark.parametrize("companion", ["CREATE_INDEX=true", "CREATE_MD5_FILE=true"])
def test_companions_beside_standard_output(tmp_path, companion):
    standard_output(TOOL, [f"I={SOURCES['markdup']}", *RG, companion], tmp_path, companion)


def test_the_step_in_a_chain_equals_picard(corpus, tmp_path):
    flags = ["--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1"]
    cases = [
        (SOURCES["markdup"], [], []),
        (corpus / "corpus.bam", ["--rgid", "g7", "--rgdt", "2020-01-02"],
         ["RGID=g7", "RGDT=2020-01-02"]),
        (SOURCES["markdup"], ["--sort-order", "queryname"], ["SORT_ORDER=queryname"]),
        (SOURCES["duplicate_marking"], ["--sort-order", "coordinate"],
         ["SORT_ORDER=coordinate"]),
        (corpus / "corpus.bam", ["--sort-order", "coordinate"], ["SORT_ORDER=coordinate"]),
        (SOURCES["t1"], ["--sort-order", "coordinate"], ["SORT_ORDER=coordinate"]),
    ]
    for i, (source, extra, picard) in enumerate(cases):
        directory = tmp_path / f"case{i}"
        directory.mkdir()
        chain_against_picard(TOOL, ["add-or-replace-read-groups", *flags, *extra], source,
                             [*RG, *picard, *SILENT], directory)


def test_values_are_distinct():
    assert len({rgid for rgid, _, _ in SEEN}) >= 3, SEEN
    assert len({order for _, order, _ in SEEN}) >= 4, SEEN
    assert len({dt for _, _, dt in SEEN if dt}) >= 10, SEEN

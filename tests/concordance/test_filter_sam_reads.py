"""FilterSamReads against Picard 3.5.0.

Picard filters INPUT through htsjdk's `FilteringSamIterator` with the filter `FILTER`
names (`AlignedFilter` judging the pairs of a queryname input together, `ReadNameFilter`,
`TagFilter`, `IntervalKeepPairFilter`) and writes through a writer presorted for the input's
order and sorting for any other SORT_ORDER. Every case goes through `compare.compare_bam`
(tier B: header lines, `@PG` ID and PP, every record as samtools prints it; a BAM case also
every record as bytes after inflation).

Inputs: the foundation corpus, t1, the exome slice and the `mate_filter_bams.py` files
(`primary.queryname.bam`, primary records of whole pairs and unpaired records in htsjdk's
queryname order, for the pair filters; `tagged.bam`, tags of every type; read lists with
every name kind: words after the name, a leading space, blank lines, `#`, CRLF and CR line
ends, absent and repeated names, gzip, names beyond ASCII on `latin1.sam`), interval lists
written by IntervalListTools (padded, merged) and the corpus and exome target lists, SAM
text and CRAM input. Options: every FILTER but includeJavascript, TAG with values of
every type, several TAG_VALUEs, each SORT_ORDER, WRITE_READS_FILES (the `.reads` files
compared byte for byte), MAX_RECORDS_IN_RAM, COMPRESSION_LEVEL, CREATE_INDEX,
CREATE_MD5_FILE, SAM text output. Error parity: each FILTER without the option it needs,
an option for another FILTER, a TAG of three characters, the pair filters on a coordinate
input, includePairedIntervals on a queryname input, a mate without MC, an empty read list,
a bad FILTER, a missing input, a missing output directory. includeJavascript and
JAVASCRIPT_FILE are refused by bammetrics alone (docs/07 Q6, D16). The step in a chain
against Picard's file. The distinct-value gate wants every filter kept and dropped
records, 3 tag types and 3 read list kinds.
"""
import subprocess

import pytest

import mate_filter_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from test_add_comments_to_bam import chain_against_picard, compare, index_and_md5, run_tool

TOOL = "FilterSamReads"
SILENT = ["VALIDATION_STRINGENCY=SILENT"]

SEEN = {"filters": set(), "tags": set(), "lists": set()}


@pytest.fixture(scope="module")
def bams(corpus, tmp_path_factory):
    return mate_filter_bams.write(corpus, tmp_path_factory.mktemp("mate_filter"))


@pytest.fixture(scope="module")
def lists(corpus, exome, tmp_path_factory):
    """Interval lists: the corpus and exome targets, padded and merged by IntervalListTools
    (bammetrics, whose lists equal Picard's)."""
    out = tmp_path_factory.mktemp("lists")
    padded = out / "padded.interval_list"
    done = run_tool("bammetrics", "IntervalListTools",
                    [f"I={corpus / 'target.interval_list'}", f"O={padded}", "PADDING=150",
                     "UNIQUE=true"])
    assert done.returncode == 0, done.stderr
    return {"target": corpus / "target.interval_list", "bait": corpus / "bait.interval_list",
            "padded": padded, "exome": exome / "targets.interval_list"}


def records_of(path):
    return int(subprocess.run(["samtools", "view", "-c", str(path)], check=True,
                              capture_output=True, text=True).stdout)


def kept(source, args, label, outcome="equal", suffix=".bam", reference=None):
    written = compare(TOOL, source, [*SILENT, *args], label, outcome, reference, suffix)
    if written.outcome == "equal":
        total = records_of(source)
        if 0 < len(written.records) < total:
            SEEN["filters"].add(next(a for a in args if a.startswith("FILTER=")))
    return written


@pytest.mark.parametrize("filter_", ["includeReadList", "excludeReadList"])
@pytest.mark.parametrize("names", ["plain", "kinds", "crlf", "cr", "plain.txt.gz", "absent"])
def test_read_lists(corpus, bams, filter_, names):
    path = bams / (f"names_{names}" if names.endswith(".gz") else f"names_{names}.txt")
    written = kept(corpus / "corpus.bam", [f"FILTER={filter_}", f"READ_LIST_FILE={path}"],
                   f"{filter_} {path.name}")
    if written.outcome == "equal":
        SEEN["lists"].add(names)


def test_read_list_beyond_ascii(bams):
    for filter_ in ("includeReadList", "excludeReadList"):
        kept(bams / "latin1.sam", [f"FILTER={filter_}",
                                   f"READ_LIST_FILE={bams / 'names_latin1.txt'}"],
             f"{filter_} names beyond ASCII", suffix=".sam")


@pytest.mark.parametrize("tag,values", [("XZ", ["alpha"]), ("XZ", ["alpha", "gamma"]),
                                        ("XA", ["a"]), ("XI", ["7"]), ("XF", ["1.5"]),
                                        ("XH", ["1AE3"]), ("XB", ["1"]), ("RG", []),
                                        ("NO", ["x"])])
@pytest.mark.parametrize("filter_", ["includeTagValues", "excludeTagValues"])
def test_tags_of_every_type(bams, tag, values, filter_):
    """Only a string tag matches a value: htsjdk compares the tag's object with strings."""
    args = [f"FILTER={filter_}", f"TAG={tag}", *(f"TAG_VALUE={v}" for v in values)]
    written = kept(bams / "tagged.bam", args, " ".join(args))
    if written.outcome == "equal":
        SEEN["tags"].add(tag)


def test_tag_on_the_corpus(corpus):
    for filter_ in ("includeTagValues", "excludeTagValues"):
        with subprocess.Popen(["samtools", "view", "-H", str(corpus / "corpus.bam")],
                              stdout=subprocess.PIPE, text=True) as head:
            groups = [line.split("ID:")[1].split("\t")[0] for line in head.stdout
                      if line.startswith("@RG")]
        kept(corpus / "corpus.bam", [f"FILTER={filter_}", "TAG=RG",
                                     *(f"TAG_VALUE={g}" for g in groups[:2])],
             f"{filter_} RG {groups[:2]}")


@pytest.mark.parametrize("filter_", ["includeAligned", "excludeAligned"])
def test_aligned_pairs_together(bams, filter_):
    source = bams / "primary.queryname.bam"
    kept(source, [f"FILTER={filter_}"], f"{filter_} primary queryname")
    kept(source, [f"FILTER={filter_}", "SORT_ORDER=coordinate"],
         f"{filter_} primary queryname to coordinate")
    kept(source, [f"FILTER={filter_}"], f"{filter_} primary queryname into SAM text",
         suffix=".sam")


@pytest.mark.parametrize("filter_", ["includeAligned", "excludeAligned"])
def test_aligned_where_picard_decides(corpus, bams, filter_):
    """A coordinate input with pairs and a queryname input with supplementary records
    fail Picard's mate assertion; the damaged queryname file may or may not."""
    for source in (corpus / "corpus.bam", bams / "damaged.queryname.bam",
                   bams / "damaged.natural.bam", FIXTURES / "markdup.bam"):
        kept(source, [f"FILTER={filter_}"], f"{filter_} {source.name}", "same")


@pytest.mark.parametrize("name", ["target", "bait", "padded"])
def test_paired_intervals(corpus, lists, name):
    kept(corpus / "corpus.bam",
         ["FILTER=includePairedIntervals", f"INTERVAL_LIST={lists[name]}"],
         f"includePairedIntervals {name}")


def test_paired_intervals_on_the_exome_and_t1(exome, lists, tmp_path):
    kept(exome / "exome.bam",
         ["FILTER=includePairedIntervals", f"INTERVAL_LIST={lists['exome']}"],
         "includePairedIntervals exome", "same")
    names = subprocess.run(["samtools", "view", str(CORPUS / "t1.bam")], check=True,
                           capture_output=True, text=True).stdout.splitlines()
    listed = tmp_path / "t1_names.txt"
    listed.write_text("".join(line.split("\t")[0] + "\n" for line in names[::7]))
    for filter_ in ("includeReadList", "excludeReadList"):
        kept(CORPUS / "t1.bam", [f"FILTER={filter_}", f"READ_LIST_FILE={listed}"],
             f"t1 {filter_}")


def test_paired_intervals_after_fixing_mates(bams, lists, tmp_path):
    """FixMateInformation (Picard's) sets MC on every mapped mate of the whole pairs of
    `primary.bam`, so no mate lacks it."""
    fixed = tmp_path / "fixed.bam"
    done = run_tool("picard", "FixMateInformation",
                    [f"I={bams / 'primary.bam'}", f"O={fixed}", *SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    for name in ("padded", "target"):
        kept(fixed, ["FILTER=includePairedIntervals", f"INTERVAL_LIST={lists[name]}"],
             f"includePairedIntervals fixed primary records, {name}")


@pytest.mark.parametrize("order", ["unsorted", "queryname", "coordinate", "duplicate"])
def test_each_sort_order(corpus, bams, order):
    kept(bams / "tagged.bam", ["FILTER=excludeTagValues", "TAG=XZ", "TAG_VALUE=beta",
                               f"SORT_ORDER={order}"], f"excludeTagValues SORT_ORDER={order}",
         "same" if order == "duplicate" else "equal")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=20", "SORT_ORDER=queryname"]])
def test_options(bams, tmp_path, args):
    kept(bams / "tagged.bam", ["FILTER=includeTagValues", "TAG=XZ", "TAG_VALUE=alpha",
                               f"TMP_DIR={tmp_path}", *args], " ".join(args))


def test_sam_text_and_cram(corpus, bams, tmp_path):
    sam = tmp_path / "tagged.sam"
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(sam),
                    str(bams / "tagged.bam")], check=True)
    kept(sam, ["FILTER=includeTagValues", "TAG=XZ", "TAG_VALUE=delta"], "SAM text input")
    cram = tmp_path / "corpus.cram"
    reference = corpus / "ref.fa"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(reference), "-o", str(cram),
                    str(corpus / "corpus.bam")], check=True)
    kept(cram, ["FILTER=includeReadList", f"READ_LIST_FILE={bams / 'names_plain.txt'}",
                f"R={reference}"], "CRAM input", reference=reference)


def test_reads_files(corpus, bams, tmp_path):
    """WRITE_READS_FILES: `<basename>.reads` of INPUT and of OUTPUT beside OUTPUT."""
    texts = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        done = run_tool(side, TOOL, [f"I={corpus / 'corpus.bam'}", f"O={directory / 'kept.bam'}",
                                     *SILENT, "FILTER=includeReadList",
                                     f"READ_LIST_FILE={bams / 'names_plain.txt'}",
                                     "WRITE_READS_FILES=true", "SORT_ORDER=queryname"])
        assert done.returncode == 0, f"{side}: {done.stderr[-800:]}"
        texts[side] = {p.name: p.read_bytes() for p in directory.iterdir()
                       if p.suffix == ".reads"}
    assert sorted(texts["picard"]) == sorted(texts["bammetrics"]) == ["corpus.reads",
                                                                      "kept.reads"]
    for name in texts["picard"]:
        assert texts["picard"][name] == texts["bammetrics"][name], name
    RESULTS.append((f"{TOOL} WRITE_READS_FILES", 2, 2, ".reads files equal byte for byte"))


def test_index_and_md5(bams, tmp_path):
    index_and_md5(TOOL, [f"I={bams / 'tagged.bam'}", *SILENT, "FILTER=includeTagValues",
                         "TAG=XZ", "TAG_VALUE=alpha"], tmp_path, "out.bai")


ERRORS = {
    "no_filter": [],
    "bad_filter": ["FILTER=includeEverything"],
    "read_list_missing": ["FILTER=includeReadList"],
    "interval_list_missing": ["FILTER=includePairedIntervals"],
    "tag_missing": ["FILTER=includeTagValues", "TAG_VALUE=x"],
    "tag_for_another": ["FILTER=includeAligned", "TAG=RG"],
    "list_for_another": ["FILTER=excludeTagValues", "TAG=RG", "READ_LIST_FILE=@LIST"],
    "tag_three_chars": ["FILTER=includeTagValues", "TAG=RGX", "TAG_VALUE=A"],
    "pairs_on_coordinate": ["FILTER=includeAligned"],
    "intervals_on_queryname": ["FILTER=includePairedIntervals", "INTERVAL_LIST=@TARGET",
                               "SORT_ORDER=queryname"],
    "empty_read_list": ["FILTER=includeReadList", "READ_LIST_FILE=@EMPTY"],
    "bad_order": ["FILTER=includeAligned", "SORT_ORDER=bogus"],
}


@pytest.mark.parametrize("case", [*ERRORS, "missing_input", "missing_directory"])
def test_errors_as_picard(corpus, bams, tmp_path, case):
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    swap = {"@LIST": str(bams / "names_plain.txt"), "@EMPTY": str(empty),
            "@TARGET": str(corpus / "target.interval_list")}
    args = [*SILENT, *(a.split("=")[0] + "=" + swap.get(a.split("=", 1)[1], a.split("=", 1)[1])
                       for a in ERRORS.get(case, ["FILTER=includeAligned"]))]
    source = corpus / "corpus.bam"
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


def test_javascript_is_refused(corpus, tmp_path):
    for args in (["FILTER=includeJavascript"], ["FILTER=includeAligned",
                                                f"JAVASCRIPT_FILE={tmp_path / 'f.js'}"]):
        done = run_tool("bammetrics", TOOL, [f"I={corpus / 'corpus.bam'}",
                                             f"O={tmp_path / 'out.bam'}", *args])
        assert done.returncode == 1 and "JAVASCRIPT" in done.stderr.upper(), done.stderr
        verb = subprocess.run([str(BAMMETRICS), "filter-sam-reads", str(corpus / "corpus.bam"),
                               "-o", str(tmp_path / "out.bam"), "--filter", "includeJavascript"],
                              capture_output=True, text=True)
        assert verb.returncode == 2, verb.stderr
        assert not (tmp_path / "out.bam").exists()
    RESULTS.append((f"{TOOL} includeJavascript", 0, 0, "refused (docs/07 Q6)"))


def test_the_step_in_a_chain_equals_picard(corpus, bams, tmp_path):
    cases = [
        (bams / "tagged.bam", ["--filter", "includeTagValues", "--tag", "XZ", "--tag-value",
                               "gamma"], ["FILTER=includeTagValues", "TAG=XZ",
                                          "TAG_VALUE=gamma"]),
        (bams / "primary.queryname.bam", ["--filter", "excludeAligned"],
         ["FILTER=excludeAligned"]),
        (corpus / "corpus.bam", ["--filter", "includeReadList", "--read-list-file",
                                 str(bams / "names_kinds.txt"), "--sort-order", "queryname"],
         ["FILTER=includeReadList", f"READ_LIST_FILE={bams / 'names_kinds.txt'}",
          "SORT_ORDER=queryname"]),
    ]
    for i, (source, extra, picard) in enumerate(cases):
        directory = tmp_path / f"case{i}"
        directory.mkdir()
        chain_against_picard(TOOL, ["filter-sam-reads", *extra], source, [*SILENT, *picard],
                             directory)


def test_values_are_distinct():
    assert len(SEEN["filters"]) >= 6, SEEN
    assert len(SEEN["tags"]) >= 3, SEEN
    assert len(SEEN["lists"]) >= 3, SEEN

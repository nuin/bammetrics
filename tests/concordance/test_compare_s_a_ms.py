"""CompareSAMs against Picard 3.5.0: the exit code, standard output byte for byte (each
header difference as Picard's three lines, the two refusals, `SAM files match.` or
`SAM files differ.`) and the SamComparisonMetric file after its StringHeader lines byte
for byte, the COMPARE_MQ histogram included.

Pairs made from the foundation corpus by `compare_sams_bams.py`: the file against itself
(its 32 QC-failed copies share a name and a place with their original, so Picard counts
32 MISSING_LEFT and says the file differs from itself), strands flipped (never
compared), unmapped flags set on either side and on both, records moved, records
missing on each side, duplicate flags flipped, mapping qualities 0, 2, 3, 4 and 255 on
records that move (every LENIENT_* option and LOW_MQ_THRESHOLD), two duplicate markings
of the corpus by `bammetrics MarkDuplicates` with two scoring strategies (LENIENT_DUP, in
memory and spilled), every header field Picard compares changed, added and dropped in
both directions with and without LENIENT_HEADER, and the fields it does not compare (GO,
PN, PP, @CO). The same records queryname sorted (in Picard's key order and in samtools'
order), unsorted and shuffled, under SO:duplicate, without @HD, under two dictionaries
whose contigs are renamed or reordered, as SAM text, and as CRAM 3.0 against the BAM.
Hand-written duplicate sets for LENIENT_DUP. Corpus t1 and the exome slice against
themselves. Error parity: the right file ends first under queryname order (Picard
throws), SO:unknown, a pair without MC under LENIENT_DUP, bad option values, one input
and three, a missing file, an OUTPUT whose directory is missing. The verb gives the
Picard face's file and exit code. Header values and read names beyond ASCII (accented,
CJK, emoji, characters that share a low byte, characters whose low byte is a tab, a line
end or NUL) in SAM text, BAM and CRAM against each other. A boolean flag before the
files in Picard 3 syntax. An unmapped record with a mapped mate and no MC under
LENIENT_DUP and under SortSam's duplicate order: counted alone, failing on a tie.

Picard runs with VALIDATION_STRINGENCY=SILENT (its default refuses corpus records),
QUIET=true and JAVA_OPTS=-Xmx8g. Needs picard and samtools; run `make concordance`.
"""
import random
import subprocess

import pytest

import compare_sams_bams as make
import revert_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CompareSAMs"
SILENT = "VALIDATION_STRINGENCY=SILENT"
COUNTS = ["MAPPINGS_MATCH", "MAPPINGS_DIFFER", "UNMAPPED_BOTH", "UNMAPPED_LEFT", "UNMAPPED_RIGHT",
          "MISSING_LEFT", "MISSING_RIGHT", "DUPLICATE_MARKINGS_DIFFER"]
SEEN = []
# Every `X differs.` line the header cases printed.
PRINTED = set()
VERDICTS = ("SAM files match.", "SAM files differ.")

MAPQ = ([(0, 0, True)] * 6 + [(3, 3, True)] * 5 + [(3, 2, True)] * 4 + [(4, 3, True)] * 3
        + [(4, 4, True)] * 3 + [(255, 255, True)] * 7 + [(255, 60, True)] * 2 + [(0, 0, False)] * 4
        + [(255, 255, False)] * 3 + [(20, 30, False)] * 5 + [(2, 2, True)] * 3)
SPARE_SQ = "@SQ\tSN:spare\tLN:100\tM5:0123456789abcdef0123456789abcdef\tUR:file:/x.fa\tAS:asm1\tSP:human"
SPARE_RG = ("@RG\tID:spare\tSM:sx\tLB:lx\tDS:text\tPU:ux\tPI:300\tCN:centre"
            "\tDT:2020-01-02T03:04:05+0000\tPL:ILLUMINA")


def replace_last(text):
    """An edit of a header's last line of one kind."""
    return lambda lines: [*lines[:-1], text(lines[-1])]


def swap(old, new):
    return replace_last(lambda line: line.replace(old, new))


def drop_tag(tag):
    return replace_last(lambda line: "\t".join(c for c in line.split("\t") if not c.startswith(tag + ":")))


# name -> the edit of the base header, by the kind of line it changes.
HEADERS = {
    "vn": {"vn": "1.5"},
    "creator": {"hd_extra": "\tCR:maker"},
    "group_order": {"hd_extra": "\tGO:none"},
    "sq_more": {"sq": lambda lines: [*lines, "@SQ\tSN:more\tLN:50"]},
    "sq_name": {"sq": swap("SN:spare", "SN:other")},
    "sq_length": {"sq": swap("LN:100", "LN:101")},
    "sq_md5": {"sq": swap("M5:0123", "M5:9123")},
    "sq_no_md5": {"sq": drop_tag("M5")},
    "sq_uri": {"sq": swap("UR:file:/x.fa", "UR:file:/y.fa")},
    "sq_assembly": {"sq": swap("AS:asm1", "AS:asm2")},
    "sq_no_assembly": {"sq": drop_tag("AS")},
    "sq_species": {"sq": swap("SP:human", "SP:mouse")},
    "sq_all": {"sq": replace_last(lambda _: "@SQ\tSN:spare\tLN:100\tM5:ffff\tUR:u\tAS:a\tSP:s")},
    "rg_more": {"rg": lambda lines: [*lines, "@RG\tID:more\tSM:sx"]},
    "rg_id": {"rg": swap("ID:spare", "ID:other")},
    "rg_sample": {"rg": swap("SM:sx", "SM:sy")},
    "rg_library": {"rg": swap("LB:lx", "LB:ly")},
    "rg_no_library": {"rg": drop_tag("LB")},
    "rg_description": {"rg": swap("DS:text", "DS:other text")},
    "rg_unit": {"rg": swap("PU:ux", "PU:uy")},
    "rg_insert": {"rg": swap("PI:300", "PI:301")},
    "rg_centre": {"rg": swap("CN:centre", "CN:center")},
    "rg_date": {"rg": swap("DT:2020-01-02", "DT:2021-01-02")},
    "rg_no_date": {"rg": drop_tag("DT")},
    "rg_platform": {"rg": swap("PL:ILLUMINA", "PL:ELEMENT")},
    "rg_all": {"rg": replace_last(lambda _: "@RG\tID:spare\tSM:a\tLB:b\tDS:c\tPU:d\tPI:1\tCN:e\tPL:f")},
    "pg_fewer": {"pg": lambda lines: lines[:-1]},
    "pg_none": {"pg": lambda lines: []},
    "pg_id": {"pg": swap("ID:last", "ID:final")},
    "pg_version": {"pg": swap("VN:1.0", "VN:2.0")},
    "pg_no_version": {"pg": drop_tag("VN")},
    "pg_command": {"pg": swap("CL:last --x", "CL:last --y")},
    "pg_name": {"pg": swap("PN:last", "PN:lastly")},
    "pg_previous": {"pg": drop_tag("PP")},
}
# Fields Picard does not compare: these headers equal the base.
NOT_COMPARED = {"group_order", "pg_name", "pg_previous"}
# Dictionary differences: a coordinate-sorted pair is then not compared.
DICTIONARY = {"sq_more", "sq_name", "sq_length", "sq_md5", "sq_no_md5", "sq_uri", "sq_assembly",
              "sq_no_assembly", "sq_species", "sq_all"}
# What LENIENT_HEADER lets differ.
LENIENT_PASSES = {"sq_md5", "sq_no_md5", "sq_uri", "sq_assembly", "sq_no_assembly", "sq_species",
                  "sq_all", "pg_fewer", "pg_none", "pg_id", "pg_version", "pg_no_version",
                  "pg_command"}


def base_header(lines):
    """The corpus header with a spare @SQ, a spare @RG and a last @PG to edit."""
    last = "@PG\tID:last\tPN:last\tPP:samtools.3\tVN:1.0\tCL:last --x"
    return make.header(lines, sq=lambda sq: [*sq, SPARE_SQ], rg=lambda rg: [*rg, SPARE_RG],
                       pg=lambda pg: [*pg, last])


def corpus_files(directory, corpus):
    """The corpus and its edited copies, name -> path."""
    src = make.Source(corpus / "corpus.bam")
    head = src.header
    files = {"corpus": make.write(directory, "corpus", head, src.records)}
    files["strand"] = make.write(directory, "strand", head, src.strand(60))
    files["unmapped_a"] = make.write(directory, "unmapped_a", head, src.unmapped(25, "a"))
    files["unmapped_b"] = make.write(directory, "unmapped_b", head, src.unmapped(400, "b"))
    files["moved"] = make.write(directory, "moved", head, src.moved(40), sort=True)
    files["without_a"] = make.write(directory, "without_a", head, src.without(13, "a"))
    files["without_b"] = make.write(directory, "without_b", head, src.without(21, "b"))
    files["dups"] = make.write(directory, "dups", head, src.duplicates(30))
    left, right = src.mapq(MAPQ)
    files["mapq_left"] = make.write(directory, "mapq_left", head, left)
    files["mapq_right"] = make.write(directory, "mapq_right", head, right, sort=True)
    several = make.Source(files["moved"], seed=5)
    files["several"] = make.write(directory, "several", head, several.duplicates(17), sort=True)
    files["sam"] = make.write(directory, "corpus_text", head, src.records, sam=True)
    cram = directory / "corpus.cram"
    make.samtools("view", "--no-PG", "-C", "-T", corpus / "ref.fa", "--output-fmt-option",
                  "version=3.0", "-o", cram, files["corpus"])
    files["cram"] = cram
    for strategy in ("SUM_OF_BASE_QUALITIES", "RANDOM"):
        marked = directory / f"marked_{strategy.lower()}.bam"
        subprocess.run([str(BAMMETRICS), "MarkDuplicates", f"I={files['corpus']}", f"O={marked}",
                        f"M={directory / 'marked.txt'}", f"DUPLICATE_SCORING_STRATEGY={strategy}"],
                       check=True, capture_output=True)
        files[f"marked_{strategy.lower()}"] = marked
    return src, files


def order_files(directory, src, files):
    """The corpus records in the other sort orders."""
    by_name = make.queryname(src.records)
    qhead = make.header(src.header, so="queryname")
    files["q"] = make.write(directory, "q", qhead, by_name)
    middle = [line for i, line in enumerate(by_name[:-40]) if i % 97 != 5] + by_name[-40:]
    files["q_fewer"] = make.write(directory, "q_fewer", qhead, middle)
    files["q_short"] = make.write(directory, "q_short", qhead, by_name[:-9])
    files["q_edits"] = make.write(directory, "q_edits", qhead, make.queryname(src.duplicates(30)))
    natural = directory / "q_natural.bam"
    make.samtools("sort", "--no-PG", "-n", "-o", natural, files["corpus"])
    files["q_natural"] = natural
    shuffled = list(src.records)
    random.Random(3).shuffle(shuffled)
    uhead = make.header(src.header, so="unsorted")
    files["unsorted"] = make.write(directory, "unsorted", uhead, src.records)
    files["shuffled"] = make.write(directory, "shuffled", uhead, shuffled)
    files["shuffled_fewer"] = make.write(directory, "shuffled_fewer", uhead, shuffled[50:])
    nohd = make.header(src.header, drop_hd=True)
    files["no_hd"] = make.write(directory, "no_hd", nohd, shuffled)
    for order in ("duplicate", "unknown", "bogus"):
        ohead = make.header(src.header, so=order)
        files[order] = make.write(directory, f"so_{order}", ohead, src.records)
    renamed = make.header(uhead, sq=lambda sq: [line.replace("SN:chr2", "SN:chrTwo") for line in sq])
    files["renamed"] = make.write(directory, "renamed", renamed,
                                  [line.replace("\tchr2\t", "\tchrTwo\t") for line in shuffled])
    reordered = make.header(uhead, sq=lambda sq: sq[::-1])
    files["reordered"] = make.write(directory, "reordered", reordered, shuffled)
    return files


def header_files(directory, src, files):
    """One BAM per header edit, coordinate sorted and unsorted."""
    base = base_header(src.header)
    # Without the QC-failed copies, so the records of a file equal themselves.
    records = [line for line in src.records if not make.flag_of(line) & 0x200][:400]
    for order in ("coordinate", "unsorted"):
        ordered = make.header(base, so=order)
        files[f"base_{order}"] = make.write(directory, f"base_{order}", ordered, records)
        for name, edit in HEADERS.items():
            files[f"{name}_{order}"] = make.write(directory, f"{name}_{order}",
                                                  make.header(ordered, **edit), records)
    comment = [*make.header(base, so="coordinate"), "@CO\ta comment"]
    files["comment_coordinate"] = make.write(directory, "comment_coordinate", comment, records)
    return files


@pytest.fixture(scope="module")
def files(tmp_path_factory, corpus):
    directory = tmp_path_factory.mktemp("compare_sams")
    src, made = corpus_files(directory, corpus)
    order_files(directory, src, made)
    header_files(directory, src, made)
    for name, (left, right) in make.hand_sets().items():
        made[f"hand_{name}_left"] = make.write(directory, f"hand_{name}_left", make.HAND_HEADER, left)
        made[f"hand_{name}_right"] = make.write(directory, f"hand_{name}_right", make.HAND_HEADER,
                                                right)
    made["reference"] = corpus / "ref.fa"
    return made


def run(side, arguments, quiet=("QUIET=true",)):
    if side == "picard":
        env = picard_env()
        # Header values beyond ASCII are printed: the JVM writes them as UTF-8.
        env["JAVA_OPTS"] = "-Xmx8g -Dstdout.encoding=UTF-8"
        return subprocess.run(["picard", TOOL, *map(str, arguments), *quiet], env=env,
                              capture_output=True, text=True)
    return subprocess.run([str(BAMMETRICS), TOOL, *map(str, arguments)], capture_output=True,
                          text=True)


def body(path):
    """The metrics file after its StringHeader lines."""
    return path.read_text().split("\n\n", 1)[1]


def row_of(path):
    lines = body(path).splitlines()
    return dict(zip(lines[1].split("\t"), lines[2].split("\t")))


def compare(left, right, args, tmp_path, label):
    """Both tools on the pair: exit code, standard output and OUTPUT equal. Returns
    the exit code, the standard output lines and the row."""
    results, outputs = {}, {}
    for side in ("picard", "ours"):
        outputs[side] = tmp_path / f"{label}.{side}.txt".replace(" ", "_")
        results[side] = run(side, [left, right, f"O={outputs[side]}", SILENT, *args])
    picard, ours = results["picard"], results["ours"]
    done = picard.stdout.splitlines()[-1:] in ([VERDICTS[0]], [VERDICTS[1]])
    assert done, f"Picard {label}: {picard.stdout[-300:]} {picard.stderr[-1500:]}"
    assert ours.returncode == picard.returncode, (
        f"{label}: Picard exit {picard.returncode}, bammetrics exit {ours.returncode}\n"
        f"Picard: {picard.stdout}\nbammetrics: {ours.stdout}{ours.stderr[-800:]}")
    assert ours.stdout == picard.stdout, (
        f"{label}: standard output differs\nPicard: {picard.stdout}\nbammetrics: {ours.stdout}")
    expected, got = body(outputs["picard"]), body(outputs["ours"])
    assert got == expected, f"{label}: OUTPUT differs\nPicard:\n{expected}\nbammetrics:\n{got}"
    row = row_of(outputs["picard"])
    SEEN.append(tuple(int(row[column]) for column in COUNTS))
    cells = len(row) + expected.count("\n")
    RESULTS.append((f"{TOOL} {label}", cells, cells,
                    f"exit {picard.returncode}, stdout and OUTPUT identical"))
    return picard.returncode, picard.stdout.splitlines(), row


def both_fail(arguments, tmp_path, label, message=None):
    """Both tools exit non-zero without a verdict; bammetrics writes no OUTPUT."""
    ours_out = tmp_path / f"{label}.ours.txt".replace(" ", "_")
    picard_out = tmp_path / f"{label}.picard.txt".replace(" ", "_")
    picard = run("picard", [*arguments, f"O={picard_out}"])
    ours = run("ours", [*arguments, f"O={ours_out}"])
    assert picard.returncode != 0 and ours.returncode != 0, (picard.stderr[-500:], ours.stderr)
    for result in (picard, ours):
        assert not any(verdict in result.stdout for verdict in VERDICTS), result.stdout
    assert ours.stdout == picard.stdout, (picard.stdout, ours.stdout)
    if message:
        assert message in picard.stderr, f"Picard {label}: {picard.stderr[-800:]}"
        assert message in ours.stderr, f"bammetrics {label}: {ours.stderr}"
    assert not ours_out.exists() and not picard_out.exists()
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))


RECORD_PAIRS = [
    ("corpus", "corpus"), ("corpus", "strand"), ("strand", "corpus"), ("corpus", "unmapped_a"),
    ("unmapped_a", "corpus"), ("unmapped_a", "unmapped_b"), ("corpus", "moved"),
    ("moved", "corpus"), ("without_a", "without_b"), ("without_b", "without_a"),
    ("corpus", "dups"), ("corpus", "several"), ("several", "without_a"), ("corpus", "sam"),
    ("sam", "moved"), ("unmapped_b", "corpus"), ("unmapped_b", "unmapped_a"),
    ("unmapped_b", "without_a"),
]


@pytest.mark.parametrize("left,right", RECORD_PAIRS)
def test_records(files, tmp_path, left, right):
    compare(files[left], files[right], [], tmp_path, f"{left} {right}")


def test_the_strand_is_never_compared(files, tmp_path):
    """Picard's alignmentsMatch compares the left strand with itself (docs/04 Q-018):
    60 records on the other strand change no count."""
    _, _, same = compare(files["corpus"], files["corpus"], [], tmp_path, "strand base")
    _, _, flipped = compare(files["corpus"], files["strand"], [], tmp_path, "strand flipped")
    assert {c: same[c] for c in COUNTS} == {c: flipped[c] for c in COUNTS}
    assert flipped["MAPPINGS_DIFFER"] == "0"


def test_a_file_with_one_key_twice_differs_from_itself(files, tmp_path):
    code, lines, row = compare(files["corpus"], files["corpus"], [], tmp_path, "corpus itself")
    assert (code, lines, row["MISSING_LEFT"], row["ARE_EQUAL"]) == (1, [VERDICTS[1]], "32", "N")


MAPQ_ARGS = [
    [], ["LENIENT_LOW_MQ_ALIGNMENT=true"], ["LENIENT_UNKNOWN_MQ_ALIGNMENT=true"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LENIENT_UNKNOWN_MQ_ALIGNMENT=true"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=0"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=2"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=4"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=255"],
    ["LENIENT_LOW_MQ_ALIGNMENT=true", "LOW_MQ_THRESHOLD=-1"],
    ["LOW_MQ_THRESHOLD=60"], ["COMPARE_MQ=true"],
    ["COMPARE_MQ=true", "LENIENT_LOW_MQ_ALIGNMENT=true", "LENIENT_UNKNOWN_MQ_ALIGNMENT=true"],
]


@pytest.mark.parametrize("args", MAPQ_ARGS, ids=lambda a: " ".join(a) or "strict")
def test_mapping_quality_options(files, tmp_path, args):
    compare(files["mapq_left"], files["mapq_right"], args, tmp_path, "mapq " + " ".join(args))


def test_low_mq_threshold_is_inclusive(files, tmp_path):
    """Picard's help says below; the code tests at most."""
    differ = {}
    for threshold in (2, 3):
        args = ["LENIENT_LOW_MQ_ALIGNMENT=true", f"LOW_MQ_THRESHOLD={threshold}"]
        _, _, row = compare(files["mapq_left"], files["mapq_right"], args, tmp_path,
                            f"inclusive {threshold}")
        differ[threshold] = int(row["MAPPINGS_DIFFER"])
    assert differ[2] - differ[3] == 9  # the five (3, 3) and the four (3, 2) records


def test_compare_mq_histogram(files, tmp_path):
    _, _, _ = compare(files["mapq_left"], files["mapq_right"], ["COMPARE_MQ=true"], tmp_path,
                      "mq histogram")
    text = body(tmp_path / "mq_histogram.picard.txt")
    assert "## HISTOGRAM\tjava.lang.String\nBIN\tVALUE\n" in text
    keys = [line.split("\t")[0] for line in text.split("BIN\tVALUE\n")[1].splitlines() if line]
    assert keys == sorted(keys) and {"255,255", "3,2", "20,30", "255,60"} <= set(keys)
    # Text order, not number order.
    assert keys.index("255,255") < keys.index("3,2")


DUP_ARGS = [[], ["LENIENT_DUP=true"], ["LENIENT_DUP=true", "MAX_RECORDS_IN_RAM=50"],
            ["LENIENT_DUP=true", "LENIENT_HEADER=true"], ["LENIENT_HEADER=true"]]


@pytest.mark.parametrize("args", DUP_ARGS, ids=lambda a: " ".join(a) or "strict")
@pytest.mark.parametrize("pair", [("marked_sum_of_base_qualities", "marked_random"),
                                  ("marked_random", "marked_sum_of_base_qualities"),
                                  ("corpus", "dups"), ("dups", "marked_random")])
def test_duplicate_markings(files, tmp_path, pair, args):
    left, right = pair
    extra = [f"TMP_DIR={tmp_path}"] if "MAX_RECORDS_IN_RAM=50" in args else []
    compare(files[left], files[right], [*args, *extra], tmp_path,
            f"{left} {right} " + " ".join(args))
    assert not [p for p in tmp_path.iterdir() if p.name.startswith("bammetrics-")]


def test_lenient_dup_forgives_another_representative(files, tmp_path):
    pair = (files["marked_sum_of_base_qualities"], files["marked_random"])
    _, _, strict = compare(*pair, ["LENIENT_HEADER=true"], tmp_path, "marked strict")
    _, _, lenient = compare(*pair, ["LENIENT_HEADER=true", "LENIENT_DUP=true"], tmp_path,
                            "marked lenient")
    assert int(strict["DUPLICATE_MARKINGS_DIFFER"]) > 100
    assert int(lenient["DUPLICATE_MARKINGS_DIFFER"]) < int(strict["DUPLICATE_MARKINGS_DIFFER"])


@pytest.mark.parametrize("lenient", ["false", "true"])
@pytest.mark.parametrize("name", list(make.hand_sets()))
def test_hand_duplicate_sets(files, tmp_path, name, lenient):
    pair = (files[f"hand_{name}_left"], files[f"hand_{name}_right"])
    args = [f"LENIENT_DUP={lenient}"]
    compare(*pair, args, tmp_path, f"hand {name} {lenient}")
    compare(*pair[::-1], args, tmp_path, f"hand {name} reversed {lenient}")
    if lenient == "true":
        compare(*pair, [*args, "MAX_RECORDS_IN_RAM=1", f"TMP_DIR={tmp_path}"], tmp_path,
                f"hand {name} spilled")


def test_hand_swaps_count(files, tmp_path):
    counts = {}
    for name in ("swap", "pairs", "unmarked", "single"):
        pair = (files[f"hand_{name}_left"], files[f"hand_{name}_right"])
        for lenient in ("false", "true"):
            _, _, row = compare(*pair, [f"LENIENT_DUP={lenient}"], tmp_path,
                                f"count {name} {lenient}")
            counts[name, lenient] = int(row["DUPLICATE_MARKINGS_DIFFER"])
    assert (counts["swap", "false"], counts["swap", "true"]) == (2, 0)
    assert (counts["pairs", "false"], counts["pairs", "true"]) == (4, 0)
    assert counts["unmarked", "true"] == 2 and counts["single", "true"] == 1


@pytest.mark.parametrize("lenient", [False, True])
@pytest.mark.parametrize("order", ["coordinate", "unsorted"])
@pytest.mark.parametrize("name", list(HEADERS))
def test_headers(files, tmp_path, name, order, lenient):
    args = ["LENIENT_HEADER=true"] if lenient else []
    base, edited = files[f"base_{order}"], files[f"{name}_{order}"]
    label = f"{name} {order} {'lenient' if lenient else 'strict'}"
    code, lines, row = compare(base, edited, args, tmp_path, label)
    _, reversed_lines, _ = compare(edited, base, args, tmp_path, label + " reversed")
    PRINTED.update(line for line in [*lines, *reversed_lines] if line.endswith(" differs."))
    equal = name in NOT_COMPARED or (lenient and name in LENIENT_PASSES)
    assert (code == 0) == equal, lines
    refused = "Cannot compare coordinate-sorted SAM files because sequence dictionaries differ."
    strict_dictionary = name in DICTIONARY and not (lenient and name in LENIENT_PASSES)
    assert (refused in lines) == (order == "coordinate" and strict_dictionary)
    if order == "unsorted":
        assert row["MAPPINGS_MATCH"] != "0"


def test_comment_lines_are_not_compared(files, tmp_path):
    code, lines, _ = compare(files["base_coordinate"], files["comment_coordinate"], [], tmp_path,
                             "comment")
    assert (code, lines) == (0, [VERDICTS[0]])


def test_every_kind_of_header_line_is_printed(files, tmp_path):
    """Each label of SamComparison's reportDifference calls appeared in `test_headers`
    or appears here."""
    assert len(PRINTED) >= 20, "test_headers ran first"
    _, lines, _ = compare(files["corpus"], files["q"], [], tmp_path, "labels sort order")
    printed = PRINTED | {line for line in lines if line.endswith(" differs.")}
    wanted = ["File format version", "File creator", "Sort order", "Sort Order",
              "Length of sequence dictionaries", "Name of sequence record 4",
              "Length of sequence spare", "Species of sequence spare", "Assembly of sequence spare",
              "MD5 of sequence spare", "URI of sequence spare", "Number of read groups",
              "Read Group ID", "Sample for read group spare", "Library for read group spare",
              "DS for read group spare", "PU for read group spare", "PI for read group spare",
              "CN for read group spare", "DT for read group spare", "PL for read group spare",
              "Number of program records", "Program Name", "VN Program Record attribute",
              "CL Program Record attribute"]
    missing = [label for label in wanted if f"{label} differs." not in printed]
    assert not missing, (missing, sorted(printed))


ORDER_PAIRS = [
    ("q", "q"), ("q_fewer", "q"), ("q", "q_edits"), ("q_edits", "q_fewer"), ("q_short", "q"),
    ("q_natural", "q_natural"), ("q_natural", "q"), ("unsorted", "shuffled"),
    ("shuffled", "unsorted"), ("shuffled_fewer", "shuffled"), ("shuffled", "shuffled_fewer"),
    ("no_hd", "no_hd"), ("no_hd", "unsorted"), ("duplicate", "duplicate"),
    ("duplicate", "unsorted"), ("corpus", "q"), ("q", "unsorted"), ("corpus", "unknown"),
    ("shuffled", "renamed"), ("renamed", "shuffled"), ("shuffled", "reordered"),
    ("reordered", "renamed"),
]


@pytest.mark.parametrize("args", [[], ["LENIENT_HEADER=true", "COMPARE_MQ=true"]],
                         ids=["strict", "lenient"])
@pytest.mark.parametrize("left,right", ORDER_PAIRS)
def test_sort_orders(files, tmp_path, left, right, args):
    compare(files[left], files[right], args, tmp_path, f"{left} {right} {' '.join(args)}")


def test_reference_names_decide_a_match(files, tmp_path):
    """A contig renamed in one file differs at every record on it; a dictionary in
    another order, with the same names, matches."""
    _, _, renamed = compare(files["shuffled"], files["renamed"], [], tmp_path, "renamed counted")
    _, _, reordered = compare(files["shuffled"], files["reordered"], [], tmp_path,
                              "reordered counted")
    assert int(renamed["MAPPINGS_DIFFER"]) > 1000 and reordered["MAPPINGS_DIFFER"] == "0"


def test_cram_against_its_bam(files, tmp_path):
    """samtools writes M5 and UR into a CRAM header, so the strict run refuses the
    coordinate walk; LENIENT_HEADER compares the records, which equal the BAM's."""
    reference = [f"R={files['reference']}"]
    code, lines, _ = compare(files["corpus"], files["cram"], reference, tmp_path, "cram strict")
    assert code == 1 and "MD5 of sequence chr1 differs." in lines
    lenient = [*reference, "LENIENT_HEADER=true", "COMPARE_MQ=true"]
    _, _, cram = compare(files["corpus"], files["cram"], lenient, tmp_path, "cram lenient")
    _, _, bam = compare(files["corpus"], files["corpus"], ["LENIENT_HEADER=true"], tmp_path,
                        "cram bam")
    assert {c: cram[c] for c in COUNTS} == {c: bam[c] for c in COUNTS}
    compare(files["cram"], files["cram"], reference, tmp_path, "cram cram")
    compare(files["cram"], files["moved"], lenient, tmp_path, "cram moved")


@pytest.mark.parametrize("name,equal", [("t1.bam", True), ("exome/exome.bam", False)])
def test_a_file_against_itself(tmp_path, name, equal):
    """Corpus t1 equals itself. The exome slice, as the foundation corpus, gives a key
    twice at one place, so Picard counts MISSING_LEFT and says it differs."""
    bam = CORPUS / name
    if not bam.exists():
        pytest.skip(f"{name} not built (make corpus)")
    code, lines, row = compare(bam, bam, ["COMPARE_MQ=true"], tmp_path, f"itself {bam.stem}")
    expected = (0, [VERDICTS[0]], "Y") if equal else (1, [VERDICTS[1]], "N")
    assert (code, lines, row["ARE_EQUAL"]) == expected
    assert (row["MISSING_LEFT"] == "0") == equal


def test_the_committed_fixtures(tmp_path):
    """tests/data/compare_sams_{left,right,other}.bam, the inputs of the golden
    records: the counts tests/data/make_compare_sams_fixtures.py states, and the two
    committed files Picard reads as a CRAM and as SAM text."""
    left, right, other = (FIXTURES / f"compare_sams_{name}.bam" for name in ("left", "right", "other"))
    code, lines, row = compare(left, right, [], tmp_path, "fixtures strict")
    assert code == 1 and lines[0] == "VN Program Record attribute differs."
    assert (row["MAPPINGS_DIFFER"], row["DUPLICATE_MARKINGS_DIFFER"]) == ("2", "4")
    lenient = ["LENIENT_HEADER=true", "LENIENT_DUP=true", "LENIENT_LOW_MQ_ALIGNMENT=true",
               "LENIENT_UNKNOWN_MQ_ALIGNMENT=true", "COMPARE_MQ=true"]
    code, lines, _ = compare(left, right, lenient, tmp_path, "fixtures lenient")
    assert (code, lines) == (0, [VERDICTS[0]])
    _, _, row = compare(left, other, [], tmp_path, "fixtures other")
    assert [row[c] for c in COUNTS] == ["9", "1", "1", "0", "1", "1", "1", "2"]
    _, _, row = compare(other, left, ["LOW_MQ_THRESHOLD=2", "LENIENT_LOW_MQ_ALIGNMENT=true"],
                        tmp_path, "fixtures other reversed")
    assert row["UNMAPPED_LEFT"] == "1"
    reference = f"R={FIXTURES / 'ref_100k.fa'}"
    code, _, _ = compare(FIXTURES / "insert_size.bam", FIXTURES / "insert_size.cram",
                         [reference, "COMPARE_MQ=true"], tmp_path, "fixture cram")
    assert code == 0
    compare(FIXTURES / "markdup.sam", FIXTURES / "markdup.bam", ["LENIENT_DUP=true"], tmp_path,
            "fixture sam")
    compare(FIXTURES / "duplicate_marking.bam", FIXTURES / "duplicate_marking.bam", [], tmp_path,
            "fixture queryname")
    compare(FIXTURES / "markdup.bam", FIXTURES / "duplicate_marking.bam", [], tmp_path,
            "fixture sort orders")


def test_new_style_arguments(files, tmp_path):
    """Picard 3 syntax: --KEY VALUE, the files anywhere among the options."""
    outputs = {}
    for side in ("picard", "ours"):
        outputs[side] = tmp_path / f"new.{side}.txt"
        arguments = ["--LENIENT_DUP", "true", files["corpus"], "--OUTPUT", outputs[side],
                     files["dups"], "--VALIDATION_STRINGENCY", "SILENT", "--COMPARE_MQ"]
        result = run(side, arguments, quiet=("--QUIET", "true"))
        assert result.returncode == 1 and result.stdout == VERDICTS[1] + "\n", result.stderr[-500:]
    assert body(outputs["ours"]) == body(outputs["picard"])
    RESULTS.append((f"{TOOL} new style", 1, 1, "OUTPUT identical"))


def test_the_right_file_ending_first_fails_in_both(files, tmp_path):
    """Picard counts the left records that remain and then reads the current record of
    the ended left file, which throws; no result is written."""
    both_fail([files["q"], files["q_short"], SILENT], tmp_path, "right ends first")
    both_fail([files["q"], files["q_short"], SILENT, "LENIENT_HEADER=true", "COMPARE_MQ=true"],
              tmp_path, "right ends first lenient")


def test_unknown_sort_order_fails_in_both(files, tmp_path):
    message = "Unrecognized sort order (unknown) found."
    both_fail([files["unknown"], files["unknown"], SILENT], tmp_path, "so unknown", message)
    both_fail([files["bogus"], files["bogus"], SILENT], tmp_path, "so bogus", message)


def test_a_reference_index_outside_the_dictionary_fails_in_both(tmp_path):
    """htsjdk resolves the reference and mate reference names of each record it
    decodes: index 77, a mate index 99 and index -5 under two @SQ lines fail the read
    at any stringency."""
    revert_bams.write_far(tmp_path)
    for name, index in (("refid_far", 77), ("mate_refid_far", 99), ("refid_negative", -5)):
        bam = tmp_path / f"{name}.bam"
        message = f"Reference name for '{index}' not found in sequence dictionary."
        both_fail([bam, bam, SILENT], tmp_path, name, message)
        both_fail([bam, bam], tmp_path, f"{name} strict", message)


def test_lenient_dup_needs_mc_once_two_pairs_differ(tmp_path):
    """htsjdk reads MC when it compares two records: one differing pair alone is
    counted, two fail."""
    def records(first, second):
        tail = f"{'ACGT' * 5}\t{'I' * 20}\tRG:Z:g1"
        return [f"p\t{99 | first}\tc1\t100\t60\t20M\t=\t300\t220\t{tail}",
                f"p\t{147 | second}\tc1\t300\t60\t20M\t=\t100\t-220\t{tail}"]
    plain = make.write(tmp_path, "plain", make.HAND_HEADER, records(0, 0))
    one = make.write(tmp_path, "one", make.HAND_HEADER, records(1024, 0))
    two = make.write(tmp_path, "two", make.HAND_HEADER, records(1024, 1024))
    _, _, row = compare(plain, one, ["LENIENT_DUP=true"], tmp_path, "one pair without MC")
    assert row["DUPLICATE_MARKINGS_DIFFER"] == "1"
    both_fail([plain, two, SILENT, "LENIENT_DUP=true"], tmp_path, "two pairs without MC",
              "Mate CIGAR (Tag MC) not found")


BAD_ARGUMENTS = {
    "threshold text": ["LOW_MQ_THRESHOLD=x"],
    "threshold null": ["LOW_MQ_THRESHOLD=null"],
    "boolean text": ["LENIENT_DUP=maybe"],
    "unknown option": ["LENIENT_255_MQ_ALIGNMENT=true"],
    "assume sorted": ["ASSUME_SORTED=true"],
    "stringency": ["VALIDATION_STRINGENCY=LOUD"],
    "option twice": ["COMPARE_MQ=true", "COMPARE_MQ=true"],
}


@pytest.mark.parametrize("name", list(BAD_ARGUMENTS))
def test_bad_arguments(files, tmp_path, name):
    both_fail([files["corpus"], files["corpus"], *BAD_ARGUMENTS[name]], tmp_path, name)


def test_inputs_are_two_files(files, tmp_path):
    corpus = files["corpus"]
    both_fail([corpus], tmp_path, "one input")
    both_fail([corpus, corpus, corpus], tmp_path, "three inputs")
    both_fail([], tmp_path, "no input")
    both_fail([corpus, tmp_path / "none.bam"], tmp_path, "missing right")
    both_fail([tmp_path / "none.bam", corpus], tmp_path, "missing left")
    both_fail([f"I={corpus}", f"I={corpus}"], tmp_path, "inputs by name")
    both_fail([f"SAM_FILES={corpus}", f"SAM_FILES={corpus}"], tmp_path, "inputs by field name")


def test_a_file_that_is_no_sam_fails_here_and_differs_in_picard(files, tmp_path):
    """htsjdk opens any text as SAM text with an empty header and reads its lines
    only when the records are compared: against a coordinate-sorted file Picard prints
    the header differences and `SAM files differ.`. bammetrics tells the kind of a file
    from its first bytes and fails (KNOWN_DIFFERENCES.md). Both exit 1."""
    text = tmp_path / "not_a_bam.bam"
    text.write_text("not a BAM\n")
    picard = run("picard", [files["corpus"], text])
    ours = run("ours", [files["corpus"], text])
    assert picard.returncode == ours.returncode == 1
    assert picard.stdout.endswith(VERDICTS[1] + "\n") and ours.stdout == ""
    assert "not_a_bam.bam" in ours.stderr
    unsorted = run("picard", [files["unsorted"], text, SILENT])
    assert unsorted.returncode == 1 and not unsorted.stdout.endswith(VERDICTS[1] + "\n")
    RESULTS.append((f"{TOOL} not a SAM file", 1, 1, "both exit 1"))


def test_an_output_in_a_missing_directory_fails_in_both(files, tmp_path):
    missing = tmp_path / "none" / "out.txt"
    for side in ("picard", "ours"):
        result = run(side, [files["corpus"], files["corpus"], f"O={missing}", SILENT])
        assert result.returncode != 0, side
    assert not missing.parent.exists()
    RESULTS.append((f"{TOOL} missing output directory", 0, 0, "both exit non-zero"))


def test_no_output_writes_no_file_and_prints_the_verdict(files, tmp_path):
    for side in ("picard", "ours"):
        result = run(side, [files["hand_swap_left"], files["hand_swap_left"], SILENT])
        assert (result.returncode, result.stdout) == (0, VERDICTS[0] + "\n"), result.stderr[-500:]
    RESULTS.append((f"{TOOL} no OUTPUT", 1, 1, "stdout identical"))


def test_the_verb_gives_the_picard_face(files, tmp_path):
    """`compare-s-a-ms LEFT RIGHT -o FILE`: the same file body and exit code; the lines
    of standard output go to standard error."""
    face = tmp_path / "face.txt"
    verb = tmp_path / "verb.txt"
    pair = [str(files["corpus"]), str(files["mapq_right"])]
    picard_face = run("ours", [*pair, f"O={face}", "COMPARE_MQ=true", "LENIENT_DUP=true"])
    by_verb = subprocess.run([str(BAMMETRICS), "compare-s-a-ms", *pair, "-o", str(verb),
                              "--compare-mq", "--lenient-dup"], capture_output=True, text=True)
    assert by_verb.returncode == picard_face.returncode == 1
    assert body(verb) == body(face)
    assert by_verb.stdout == "" and "compare-s-a-ms: SAM files differ." in by_verb.stderr
    piped = subprocess.run([str(BAMMETRICS), "compare-s-a-ms", *pair, "-q", "--compare-mq",
                            "--lenient-dup"], capture_output=True, text=True)
    assert piped.stdout.split("\n\n", 1)[1] == body(face)
    assert piped.stderr == ""


def test_bammetrics_refuses_an_output_that_is_an_input(files, tmp_path):
    copy = tmp_path / "copy.bam"
    copy.write_bytes(files["hand_swap_left"].read_bytes())
    before = copy.read_bytes()
    result = run("ours", [copy, files["hand_swap_right"], f"O={copy}"])
    assert result.returncode != 0 and copy.read_bytes() == before, result.stderr


# Header values and read names beyond ASCII. htsjdk holds strings: a BAM's bytes one
# character each, SAM text and a CRAM decoded as UTF-8. samtools copies the UTF-8 bytes
# of SAM text into a BAM, so the BAM of `é` holds the two characters `Ã©`.
WIDE = {
    "latin": ("caf\u00e9", "caf\u00e8"),
    "cjk": ("\u6f22\u5b57", "\u6f22\u5b50"),
    "emoji": ("x\U0001F600", "x\U0001F601"),
    "low_byte": ("A", "\u0141"),
    "low_bytes": ("\u0141", "\u0241"),
    "tab_byte": ("a\u0109b", "ab"),
    "line_byte": ("a\u010ab", "ab"),
    "nul_byte": ("a\u0100b", "ab"),
}
WIDE_KINDS = [("sam", "sam"), ("sam", "bam"), ("bam", "sam"), ("sam", "cram"), ("cram", "sam"),
              ("cram", "bam"), ("bam", "bam")]


def wide_file(directory, name, header, records, kind):
    """SAM text written as UTF-8, or the BAM or CRAM samtools makes of it."""
    text = directory / f"{name}.sam"
    text.write_text("\n".join([*header, *records]) + "\n", encoding="utf-8")
    if kind == "sam":
        return text
    out = directory / f"{name}.{kind}"
    if kind == "bam":
        make.samtools("view", "--no-PG", "-b", "-o", out, text)
    else:
        make.samtools("view", "--no-PG", "-C", "-T", wide_reference(directory), "-o", out, text)
    return out


def wide_reference(directory):
    reference = directory / "wide_ref.fa"
    if not reference.exists():
        reference.write_text(">c1\n" + "ACGT" * 1250 + "\n")
        make.samtools("faidx", reference)
    return reference


def wide_header(value, order="unsorted"):
    return [f"@HD\tVN:1.6\tSO:{order}", "@SQ\tSN:c1\tLN:5000",
            f"@RG\tID:g1\tSM:{value}\tLB:lib\tDS:{value}\tPU:{value}\tCN:{value}",
            f"@PG\tID:p1\tVN:{value}\tCL:run {value}"]


def wide_record(name, flag=0, pos=100):
    return f"{name}\t{flag}\tc1\t{pos}\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:g1"


def same_outcome(left, right, args, tmp_path, label):
    """Both tools on the pair, whatever Picard does: the exit code and standard output
    equal, and OUTPUT equal or written by neither. Returns Picard's exit code and row."""
    results, outputs = {}, {}
    for side in ("picard", "ours"):
        outputs[side] = tmp_path / f"{label}.{side}.txt".replace(" ", "_").replace("/", "_")
        results[side] = run(side, [left, right, f"O={outputs[side]}", SILENT, *args])
    picard, ours = results["picard"], results["ours"]
    assert ours.returncode == picard.returncode, (
        f"{label}: Picard exit {picard.returncode}, bammetrics exit {ours.returncode}\n"
        f"Picard: {picard.stdout}{picard.stderr[-600:]}\nbammetrics: {ours.stdout}{ours.stderr[-600:]}")
    assert ours.stdout == picard.stdout, (
        f"{label}: standard output differs\nPicard: {picard.stdout}\nbammetrics: {ours.stdout}")
    assert outputs["ours"].exists() == outputs["picard"].exists(), f"{label}: {ours.stderr[-400:]}"
    row = None
    if outputs["picard"].exists():
        assert body(outputs["ours"]) == body(outputs["picard"]), label
        row = row_of(outputs["picard"])
    RESULTS.append((f"{TOOL} {label}", 1, 1, f"exit {picard.returncode}, stdout and OUTPUT identical"))
    return picard.returncode, row


@pytest.mark.parametrize("kinds", WIDE_KINDS, ids="/".join)
@pytest.mark.parametrize("name", list(WIDE))
def test_header_values_beyond_ascii(tmp_path, name, kinds):
    """SM, DS, PU, CN, VN and CL that differ in one character, and the same value on
    both sides, in each pair of file kinds. Two SAM text values that share their low
    bytes (`A` and U+0141) differ; a character whose low byte is a tab or a line end
    stays in its value."""
    ours, theirs = WIDE[name]
    reference = [f"R={wide_reference(tmp_path)}"]
    records = [wide_record("r1")]
    left = wide_file(tmp_path, "left", wide_header(ours), records, kinds[0])
    right = wide_file(tmp_path, "right", wide_header(theirs), records, kinds[1])
    same = wide_file(tmp_path, "same", wide_header(ours), records, kinds[1])
    label = f"header {name} {'/'.join(kinds)}"
    code, lines, _ = compare(left, right, reference, tmp_path, label)
    assert code == 1 and "DS for read group g1 differs." in lines
    if kinds == ("sam", "sam"):
        assert f"File 1: {ours}" in lines and f"File 2: {theirs}" in lines
    code, lines, _ = compare(left, same, reference, tmp_path, f"{label} one value")
    if kinds[0] == kinds[1]:
        assert code == 0, lines


@pytest.mark.parametrize("order", ["coordinate", "queryname", "unsorted"])
@pytest.mark.parametrize("name", ["latin", "cjk", "emoji", "low_byte", "low_bytes"])
def test_read_names_beyond_ascii(tmp_path, name, order):
    """Read names that differ in one character, and one name on both sides, in SAM text
    and BAM. A CRAM with such a name is left out: bammetrics refuses it for every tool
    (KNOWN_DIFFERENCES.md). Under queryname order Picard throws when the right file
    ends first, and bammetrics fails there too."""
    ours, theirs = WIDE[name]
    header = wide_header("s", order)
    for kinds in (("sam", "sam"), ("bam", "bam"), ("sam", "bam"), ("bam", "sam")):
        left = wide_file(tmp_path, "left", header, [wide_record(ours)], kinds[0])
        right = wide_file(tmp_path, "right", header, [wide_record(theirs)], kinds[1])
        same = wide_file(tmp_path, "same", header, [wide_record(ours)], kinds[1])
        label = f"names {name} {order} {'/'.join(kinds)}"
        code, row = same_outcome(left, right, [], tmp_path, label)
        assert code == 1
        if kinds == ("sam", "sam") and row is not None:
            assert (row["MISSING_LEFT"], row["MISSING_RIGHT"], row["MAPPINGS_MATCH"]) == ("1", "1", "0")
        code, row = same_outcome(left, same, [], tmp_path, f"{label} one name")
        if kinds[0] == kinds[1]:
            assert code == 0 and row["MAPPINGS_MATCH"] == "1"


def test_a_flag_before_the_files(files, tmp_path):
    """Picard 3 syntax: a boolean option takes the next word only when it is true or
    false, so `--LENIENT_DUP LEFT RIGHT` is a flag and two files."""
    left, right = files["hand_swap_left"], files["hand_swap_right"]
    forms = [["--LENIENT_DUP", left, right], [left, "--LENIENT_DUP", right],
             ["--LENIENT_DUP", "TRUE", left, right], ["--LENIENT_DUP", "False", left, right],
             ["--COMPARE_MQ", "--LENIENT_DUP", left, right], [left, right, "--LENIENT_DUP"]]
    for at, form in enumerate(forms):
        outputs = {}
        for side in ("picard", "ours"):
            outputs[side] = tmp_path / f"flag{at}.{side}.txt"
            result = run(side, [*form, "--OUTPUT", outputs[side], "--VALIDATION_STRINGENCY", "SILENT"],
                         quiet=("--QUIET", "true"))
            assert result.stdout.splitlines()[-1:] in ([VERDICTS[0]], [VERDICTS[1]]), (
                side, form[:2], result.stderr[-500:])
            outputs[side + " code"] = result.returncode
        assert outputs["ours code"] == outputs["picard code"], form[:2]
        assert body(outputs["ours"]) == body(outputs["picard"]), form[:2]
    RESULTS.append((f"{TOOL} a flag before the files", len(forms), len(forms), "OUTPUT identical"))


def unmapped_beside_a_mapped_mate(name, at, mark):
    """A pair whose second read is unmapped at its mate's place and has no MC."""
    tail = f"{'ACGT' * 5}\t{'I' * 20}\tRG:Z:g1"
    return [f"{name}\t73\tc1\t{at}\t60\t20M\t=\t{at}\t0\t{tail}\tMC:Z:20M",
            f"{name}\t{133 | mark}\tc1\t{at}\t0\t*\t=\t{at}\t0\t{tail}"]


def test_lenient_dup_scores_an_unmapped_record_only_on_a_tie(tmp_path):
    """htsjdk reads the mate CIGAR of an unmapped record only to score it, and scores
    a record only against one that ties with it in the file order: two such records at
    two places are counted, two at one place fail in both tools."""
    def pairs(second, mark):
        return [*unmapped_beside_a_mapped_mate("u", 100, mark),
                *unmapped_beside_a_mapped_mate("v", second, mark)]
    left = make.write(tmp_path, "apart_left", make.HAND_HEADER, pairs(200, 0))
    right = make.write(tmp_path, "apart_right", make.HAND_HEADER, pairs(200, 1024))
    _, _, row = compare(left, right, ["LENIENT_DUP=true"], tmp_path, "unmapped without MC apart")
    assert row["DUPLICATE_MARKINGS_DIFFER"] == "2"
    compare(left, right, [], tmp_path, "unmapped without MC apart strict")
    left = make.write(tmp_path, "tie_left", make.HAND_HEADER, pairs(100, 0))
    right = make.write(tmp_path, "tie_right", make.HAND_HEADER, pairs(100, 1024))
    both_fail([left, right, SILENT, "LENIENT_DUP=true"], tmp_path, "unmapped without MC at one place")


@pytest.mark.parametrize("label,places,fails", [
    ("alone", [100], False), ("apart", [200, 100], False), ("at one place", [100, 100], True)])
def test_the_duplicate_order_scores_an_unmapped_record_only_on_a_tie(tmp_path, label, places, fails):
    """The comparator LENIENT_DUP sorts with is SortSam's SORT_ORDER=duplicate: the same
    records through both tools' SortSam come out in one order, or fail in both."""
    records = [unmapped_beside_a_mapped_mate(f"r{at}", place, 0)[1] for at, place in enumerate(places)]
    source = make.write(tmp_path, "unsorted", ["@HD\tVN:1.6\tSO:unsorted", *make.HAND_HEADER[1:]], records)
    sorted_records = {}
    for side, command in (("picard", ["picard", "SortSam"]), ("ours", [str(BAMMETRICS), "SortSam"])):
        out = tmp_path / f"sorted.{side}.sam"
        env = picard_env()
        env["JAVA_OPTS"] = "-Xmx8g"
        result = subprocess.run([*command, f"I={source}", f"O={out}", "SORT_ORDER=duplicate", SILENT,
                                 "QUIET=true"], env=env, capture_output=True, text=True)
        assert (result.returncode != 0) == fails, (side, label, result.stderr[-600:])
        if not fails:
            sorted_records[side] = [line for line in out.read_text().splitlines() if not line.startswith("@")]
    if not fails:
        assert sorted_records["ours"] == sorted_records["picard"]
        assert len(sorted_records["ours"]) == len(places)
    RESULTS.append((f"SortSam duplicate order, unmapped without MC {label}", 1, 1,
                    "both fail" if fails else "records identical"))


def test_values_are_distinct():
    """The comparisons cover the counts: three distinct values of each column."""
    assert len(SEEN) >= 100, "the comparisons ran first"
    for index, column in enumerate(COUNTS):
        values = {seen[index] for seen in SEEN}
        assert len(values) >= 3, f"{column}: {sorted(values)}"

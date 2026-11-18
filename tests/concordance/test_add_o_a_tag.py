"""AddOATag against Picard 3.5.0.

Picard reads the input with a `SamReader` and writes every record through a presorted
`SAMFileWriter`, with the record's alignment appended to its `OA` tag
(`RNAME,POS,STRAND,CIGAR,MAPQ,NM;`, or `*,0,STRAND,*,255,;` for an unmapped record), or
only on the records over `INTERVAL_LIST`. Every case goes through `compare.compare_bam`
with `record_bytes=True` (tier B: the header lines, `@PG` ID and PP, every record as
samtools prints it, then every record as bytes after inflation: the `bin` field, the
place of `OA` among the tags and the width of each integer tag).

Inputs: the fixture BAMs, the foundation corpus, t1 to t3, the exome slice, each also with
its `NM` tags stripped and with `OA` already on its records (one entry on every third
record, two on every fifth); a hand-built BAM of reads off both contig ends, unplaced and
placed unmapped reads, secondary and supplementary records, CIGARs without a reference
base and `NM` of every scalar type; CRAM 3.0 of the fixture and the foundation corpus;
SAM text in, and SAM text out. INTERVAL_LIST: an interval list of each input's contigs
(intervals at both ends, one past the last read, an empty one, one on a contig the input
lacks), the corpus target list; a BED is refused by both. Options: COMPRESSION_LEVEL,
CREATE_INDEX (`out.bai`), CREATE_MD5_FILE, MAX_RECORDS_IN_RAM, REFERENCE_SEQUENCE,
VALIDATION_STRINGENCY. Error parity: a BED, a missing INTERVAL_LIST, records out of the
header's order, OVERWRITE_TAG and ASSUME_SORTED (no such options), a bad SAM line at
STRICT, a missing input and a missing output directory. bammetrics alone refuses an `NM`
or `OA` array, which Picard prints as a Java object reference, and a CRAM output. The step
in a chain against Picard's file. The distinct-value gate wants three entry counts, three
`NM` texts and three tagged counts.
"""
import subprocess
from pathlib import Path

import pysam
import pytest

from cases import CORPUS, FIXTURES
from conftest import RESULTS
from test_add_comments_to_bam import chain_against_picard, compare, index_and_md5, run_tool

TOOL = "AddOATag"
REF = FIXTURES / "ref_100k.fa"
SILENT = "VALIDATION_STRINGENCY=SILENT"

SOURCES = {
    "markdup": FIXTURES / "markdup.bam",
    "insert_size": FIXTURES / "insert_size.bam",
    "duplicate_marking": FIXTURES / "duplicate_marking.bam",
    "oxog": FIXTURES / "oxog.bam",
    "artifact": FIXTURES / "artifact.bam",
    "mt_copy_number": FIXTURES / "mt_copy_number.bam",
    "t1": CORPUS / "t1.bam",
    "t2": CORPUS / "t2.bam",
    "t3": CORPUS / "t3.bam",
}

# (OA entries per tagged record, NM text of the last entry, tagged records) of every
# case that compared equal: the gate's input.
SEEN = []


def source_path(name, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    return SOURCES[name]


def note(written):
    tagged = 0
    for line in written.records:
        oa = next((f[5:] for f in line.split("\t")[11:] if f.startswith("OA:Z:")), None)
        if oa is None:
            continue
        tagged += 1
        entries = oa.split(";")[:-1]
        SEEN.append((len(entries), entries[-1].split(",")[-1], None))
    SEEN.append((None, None, tagged))


def case(source, args, label, outcome="equal", suffix=".bam"):
    written = compare(TOOL, source, [*args, SILENT], label, outcome, suffix=suffix)
    if written.outcome == "equal":
        note(written)
    return written


def contigs(source):
    """(name, length) of each `@SQ` line of `source`."""
    header = subprocess.run(["samtools", "view", "--no-PG", "-H", str(source)], check=True,
                            capture_output=True, text=True).stdout
    out = []
    for line in header.splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            out.append((fields["SN"], int(fields["LN"])))
    return out


def interval_list(source, path):
    """An interval list of `source`'s contigs: both ends of each, a window a third of
    the way in, an empty interval, and one on a contig the input does not have."""
    lines = ["@HD\tVN:1.6"]
    seq = contigs(source)
    lines += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in seq]
    lines.append("@SQ\tSN:zz_absent\tLN:1000")
    for name, length in seq:
        third = max(1, length // 3)
        lines += [f"{name}\t1\t{min(length, 150)}\t+\tstart",
                  f"{name}\t{third}\t{min(length, third + length // 20)}\t-\tthird",
                  f"{name}\t{max(1, length - 100)}\t{length}\t+\tend",
                  f"{name}\t{min(length, third + 10)}\t{min(length, third + 9)}\t+\tempty"]
    lines.append("zz_absent\t1\t1000\t+\tabsent")
    path.write_text("\n".join(lines) + "\n")
    return path


def edited(source, path, strip_nm=False, add_oa=False):
    """`source` with its NM tags removed or an earlier OA on some records."""
    with pysam.AlignmentFile(str(source), "rb", check_sq=False) as bam, \
            pysam.AlignmentFile(str(path), "wb", template=bam) as out:
        for n, read in enumerate(bam.fetch(until_eof=True)):
            if strip_nm and read.has_tag("NM"):
                read.set_tag("NM", None)
            if add_oa and n % 3 == 0:
                read.set_tag("OA", "chrX,100,+,10M,60,1;", "Z")
            if add_oa and n % 5 == 0:
                read.set_tag("OA", "chrY,5,-,3S7M,0,;chrX,9,+,10M,3,0;", "Z")
            out.write(read)
    return path


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, tmp_path, name):
    source = source_path(name, corpus, exome)
    case(source, [], f"{name}")
    case(source, [f"L={interval_list(source, tmp_path / 'iv.interval_list')}"],
         f"{name} INTERVAL_LIST")
    case(edited(source, tmp_path / "no_nm.bam", strip_nm=True), [], f"{name} without NM")
    case(edited(source, tmp_path / "oa.bam", add_oa=True), [], f"{name} with OA")


def test_the_corpus_targets(corpus):
    case(corpus / "corpus.bam", [f"L={corpus / 'target.interval_list'}"],
         "foundation target.interval_list")


@pytest.fixture(scope="module")
def odd_shapes(tmp_path_factory):
    """Reads off both contig ends, unplaced and placed unmapped reads, secondary and
    supplementary records, CIGARs without a reference base, NM of each scalar type and
    earlier OA values of other types, as SAM text and as a BAM (both SO:coordinate)."""
    directory = tmp_path_factory.mktemp("shapes")
    seq, qual = "ACGTACGTAC", "IIIIIIIIII"

    def rec(name, flag, rname, pos, cigar, *tags, mapq=60):
        bases, quals = (seq, qual) if cigar != "*" or flag & 4 else ("*", "*")
        return "\t".join([name, str(flag), rname, str(pos), str(mapq), cigar, "*", "0", "0",
                          bases, quals, *tags])

    lines = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:c1\tLN:1000", "@SQ\tSN:c2\tLN:500",
             rec("before_start", 0, "c1", 1, "5S5M", "NM:i:0"),
             rec("nm", 0, "c1", 10, "10M", "NM:i:2"),
             rec("rev", 16, "c1", 20, "4M2I4M", mapq=30),
             rec("oa_one", 0, "c1", 30, "10M", "OA:Z:c2,5,+,10M,60,1;", "NM:i:0"),
             rec("oa_two", 256, "c1", 40, "5S5M", "OA:Z:c2,5,+,10M,60,1;c1,7,-,3S7M,10,;",
                 mapq=0),
             rec("nm_float", 2048, "c1", 50, "10M", "NM:f:1.5", mapq=1),
             rec("nm_text", 0, "c1", 60, "10M", "NM:Z:abc", mapq=1),
             rec("nm_char", 0, "c1", 70, "10M", "NM:A:x", mapq=1),
             rec("oa_int", 0, "c1", 80, "10M", "OA:i:5", mapq=1),
             rec("tag_order", 0, "c1", 90, "10M", "NM:i:300", "XA:i:1", "AA:Z:a", "zz:i:7"),
             rec("placed_unmapped", 20, "c1", 100, "*", "NM:i:1"),
             rec("insertion_only", 0, "c1", 120, "10I"),
             rec("deletion", 0, "c1", 990, "5M20D5M", "NM:i:20"),
             rec("off_end", 16, "c1", 995, "10M", "NM:i:5"),
             rec("c2_start", 0, "c2", 1, "10M", "NM:i:1000000"),
             rec("c2_off_end", 0, "c2", 498, "2M8S"),
             rec("unmapped", 4, "*", 0, "*", "NM:i:1"),
             rec("unmapped_rev", 20, "*", 0, "*", "OA:Z:c1,1,+,10M,3,;")]
    sam = directory / "shapes.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = directory / "shapes.bam"
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(bam), str(sam)], check=True)
    ivl = directory / "shapes.interval_list"
    ivl.write_text("@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n"
                   "c1\t15\t15\t+\tt1\nc1\t96\t98\t+\tt2\nc1\t1000\t1000\t+\tt3\n"
                   "c2\t1\t1\t+\tt4\nc2\t499\t500\t+\tt5\n")
    return {"sam": sam, "bam": bam, "intervals": ivl}


@pytest.mark.parametrize("kind", ["sam", "bam"])
@pytest.mark.parametrize("suffix", [".bam", ".sam"])
def test_odd_shapes(odd_shapes, kind, suffix):
    source = odd_shapes[kind]
    case(source, [], f"shapes {kind} to {suffix}", suffix=suffix)
    case(source, [f"L={odd_shapes['intervals']}"], f"shapes {kind} INTERVAL_LIST to {suffix}",
         suffix=suffix)


@pytest.fixture(scope="module")
def crams(corpus, tmp_path_factory):
    directory = tmp_path_factory.mktemp("cram")
    corpus_cram = directory / "corpus.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(corpus / "ref.fa"),
                    "-o", str(corpus_cram), str(corpus / "corpus.bam")], check=True)
    return {"insert_size": (FIXTURES / "insert_size.cram", REF),
            "foundation": (corpus_cram, corpus / "ref.fa")}


@pytest.mark.parametrize("name", ["insert_size", "foundation"])
def test_a_cram(crams, name):
    source, reference = crams[name]
    case(source, [f"R={reference}"], f"{name} CRAM")


def test_sam_text_in_and_out(corpus, tmp_path):
    corpus_sam = tmp_path / "corpus.sam"
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(corpus_sam),
                    str(corpus / "corpus.bam")], check=True)
    for source in (FIXTURES / "markdup.sam", corpus_sam):
        case(source, [], f"{source.name} to BAM")
        case(source, [], f"{source.name} to SAM", suffix=".sam")
    case(FIXTURES / "markdup.bam", [f"L={FIXTURES / 'targets.interval_list'}"],
         "markdup INTERVAL_LIST to SAM", suffix=".sam")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=1"], [f"R={REF}"],
                                  ["VALIDATION_STRINGENCY=LENIENT"]])
def test_options(args):
    if args[0].startswith("VALIDATION_STRINGENCY"):
        note(compare(TOOL, SOURCES["markdup"], args, f"markdup {args[0]}"))
    else:
        case(SOURCES["markdup"], args, f"markdup {' '.join(args)}")


def test_strict_is_the_default():
    """markdup.bam passes Picard's STRICT checks: the default run equals too."""
    note(compare(TOOL, SOURCES["markdup"], [], "markdup at STRICT"))


def test_index_and_md5(tmp_path):
    index_and_md5(TOOL, [f"I={SOURCES['markdup']}"], tmp_path, "out.bai")


@pytest.fixture(scope="module")
def bad_inputs(tmp_path_factory):
    directory = tmp_path_factory.mktemp("bad")
    unsorted = directory / "unsorted.sam"
    unsorted.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:1000\n"
                        "r1\t0\tc1\t50\t0\t2M\t*\t0\t0\tAC\tII\n"
                        "r2\t0\tc1\t5\t0\t2M\t*\t0\t0\tAC\tII\n")
    bad_line = directory / "bad_line.sam"
    bad_line.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000\n"
                        "r1\t0\tc1\t5\t60\t*\t*\t0\t0\tAC\tII\n")
    bed = directory / "targets.bed"
    bed.write_text("chr1\t1000\t1100\n")
    arrays = directory / "arrays.sam"
    arrays.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c1\tLN:1000\n"
                      "r1\t0\tc1\t5\t0\t2M\t*\t0\t0\tAC\tII\tNM:B:c,1,2\n"
                      "r2\t0\tc1\t5\t0\t2M\t*\t0\t0\tAC\tII\tOA:H:1A\n")
    return {"unsorted": unsorted, "bad_line": bad_line, "bed": bed, "arrays": arrays,
            "directory": directory}


@pytest.mark.parametrize("name", ["bed", "missing_intervals", "unsorted", "overwrite_tag",
                                  "assume_sorted", "strict", "missing_input"])
def test_errors_as_picard(bad_inputs, tmp_path, name):
    source, args = SOURCES["markdup"], [SILENT]
    if name == "bed":
        args.append(f"L={bad_inputs['bed']}")
    elif name == "missing_intervals":
        args.append(f"L={tmp_path / 'nope.interval_list'}")
    elif name == "unsorted":
        source = bad_inputs["unsorted"]
    elif name == "overwrite_tag":
        args.append("OVERWRITE_TAG=true")
    elif name == "assume_sorted":
        args.append("ASSUME_SORTED=true")
    elif name == "strict":
        source, args = bad_inputs["bad_line"], []
    elif name == "missing_input":
        source = tmp_path / "nope.bam"
    compare(TOOL, source, args, f"error {name}", "fail")


def test_a_missing_output_directory(tmp_path):
    for side in ("picard", "bammetrics"):
        done = run_tool(side, TOOL, [f"I={SOURCES['markdup']}", f"O={tmp_path / 'no' / 'out.bam'}"])
        assert done.returncode != 0, side
    assert not (tmp_path / "no").exists()
    RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))


def test_an_array_value_is_refused(bad_inputs, tmp_path):
    """Picard prints an NM or OA array as a Java object reference (`[B@1a2b3c`), which
    changes from run to run; bammetrics fails with the record's name and writes nothing."""
    picard_out, ours_out = tmp_path / "picard.sam", tmp_path / "ours.sam"
    args = [f"I={bad_inputs['arrays']}", SILENT]
    done = run_tool("picard", TOOL, [*args, f"O={picard_out}"])
    assert done.returncode == 0, done.stderr[-800:]
    assert "[B@" in picard_out.read_text()
    ours = run_tool("bammetrics", TOOL, [*args, f"O={ours_out}"])
    assert ours.returncode == 1 and "r1: its NM tag is an array" in ours.stderr, ours.stderr
    assert not ours_out.exists()
    RESULTS.append((f"{TOOL} array NM", 0, 0, "Picard prints an object reference, bammetrics fails"))


def test_a_cram_output_is_refused(tmp_path):
    out = tmp_path / "out.cram"
    done = run_tool("bammetrics", TOOL, [f"I={SOURCES['markdup']}", f"O={out}", f"R={REF}"])
    assert done.returncode == 1 and "CRAM output is wave" in done.stderr, done.stderr
    assert not out.exists()


def test_the_step_in_a_chain_equals_picard(corpus, tmp_path):
    for source in (SOURCES["markdup"], corpus / "corpus.bam"):
        directory = tmp_path / source.stem
        directory.mkdir()
        chain_against_picard(TOOL, ["add-o-a-tag"], source, [SILENT], directory)
    directory = tmp_path / "intervals"
    directory.mkdir()
    targets = FIXTURES / "targets.interval_list"
    chain_against_picard(TOOL, ["add-o-a-tag", "-L", str(targets)], SOURCES["markdup"],
                         [f"L={targets}"], directory)


def test_values_are_distinct():
    """Records carry one, two and three OA entries; NM prints empty, as an integer and in
    other forms; the cases tag at least three different numbers of records."""
    entries = {count for count, _, _ in SEEN if count is not None}
    nm = {text for _, text, _ in SEEN if text is not None}
    tagged = {count for _, _, count in SEEN if count is not None}
    assert len(entries) >= 3, entries
    assert "" in nm and len(nm) >= 3, nm
    assert len(tagged) >= 3, tagged

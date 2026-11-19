"""GatherBamFiles against Picard 3.5.0.

Every input name ends in `.bam`: Picard block copies (htsjdk
`BamFileIoUtils.gatherWithBlockCopying`), the first file whole, each later file from the
block of its first record, which htsjdk deflates again, then one end-of-file block. Tier A:
at COMPRESSION_LEVEL=0 (stored blocks) and with one input the two files are equal byte for
byte (`==` on the bytes). At level 1 to 9 every block is equal but the one deflated block
per later input with records, whose bytes come from Picard's zlib and bammetrics' zlib-rs;
those blocks inflate to equal bytes, and compare_bam with record_bytes=True passes.

Inputs: t1, the foundation corpus and the exome slice split by contig with samtools (the
pieces gathered back in order, and some in another order), the fixture BAMs, a BAM with no
record first and later, a BAM without its end-of-file block first and later, lists of
paths (nested, blank lines, spaces). A SAM or CRAM input sends both tools down Picard's
`gatherNormally` (tier B, records as bytes). CREATE_INDEX names the index `out.bai` beside
`out.bam` and beside `out.dat`, and follows htsjdk's `BAMIndexer` on the header order
(coordinate, unsorted, no `@HD`). Error parity: a file cut inside a block, text named `.bam`,
bgzip SAM named `.bam` after a BAM, a missing input, a list naming a missing file, a CRAM
without R, a missing output directory. bammetrics alone refuses an input with records
under another sequence dictionary, a block copy into a `.sam` name and an OUTPUT that is
an input (docs/07 decision rows).
"""
import struct
import subprocess
import zlib
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_bam
from conftest import RESULTS
from test_add_comments_to_bam import index_and_md5, index_by_sort_order, records_of, run_tool

TOOL = "GatherBamFiles"
T1 = CORPUS / "t1.bam"
REF = CORPUS / "reference.fa"


def blocks(path):
    """The BGZF blocks of `path`: raw bytes and inflated bytes."""
    data = Path(path).read_bytes()
    at, out = 0, []
    while at < len(data):
        size = struct.unpack_from("<H", data, at + 16)[0] + 1
        raw = data[at:at + size]
        out.append((raw, zlib.decompress(raw[18:-8], -15)))
        at += size
    return out


def split(bam, directory, name):
    """`bam` cut by contig into BAMs written by samtools, in header order."""
    contigs = [line.split("\t")[0] for line in subprocess.run(
        ["samtools", "idxstats", str(bam)], check=True, capture_output=True,
        text=True).stdout.splitlines()]
    parts = []
    for contig in [c for c in contigs if c != "*"]:
        part = directory / f"{name}_{contig}.bam"
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(part), str(bam), contig],
                       check=True)
        parts.append(part)
    unplaced = directory / f"{name}_unplaced.bam"
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(unplaced), str(bam), "*"],
                   check=True)
    return parts + [unplaced]


@pytest.fixture(scope="module")
def pieces(tmp_path_factory, corpus, exome):
    directory = tmp_path_factory.mktemp("gather")
    return {
        "t1": split(T1, directory, "t1"),
        "foundation": split(corpus / "corpus.bam", directory, "foundation"),
        "exome": split(exome / "exome.bam", directory, "exome"),
        "dir": directory,
    }


def both(args, directory, out_name="out.bam"):
    """Both tools on `args`, each into a directory of its own: (Picard, bammetrics) as
    (exit code, output path, stderr)."""
    result = {}
    for side in ("picard", "bammetrics"):
        side_dir = directory / side
        side_dir.mkdir(exist_ok=True)
        out = side_dir / out_name
        done = run_tool(side, TOOL, [*args, f"O={out}"])
        result[side] = (done.returncode, out, done.stderr)
    return result["picard"], result["bammetrics"]


def tier_a(inputs, directory, label, extra=(), gathered=None):
    """Picard's bytes at COMPRESSION_LEVEL=0, and at the default level every block but the
    deflated ones (at most one per later file of the `gathered` the lists unroll to), which
    inflate to Picard's."""
    gathered = len(inputs) if gathered is None else gathered
    args = [f"I={p}" for p in inputs] + list(extra)
    for level in ("0", "5"):
        level_dir = directory / f"level{level}"
        level_dir.mkdir()
        (prc, pout, perr), (orc, oout, oerr) = both([*args, f"COMPRESSION_LEVEL={level}"],
                                                    level_dir)
        assert prc == 0, perr[-800:]
        assert orc == 0, oerr[-800:]
        picard, ours = pout.read_bytes(), oout.read_bytes()
        if level == "0" or gathered == 1:
            assert picard == ours, f"{label}: level {level} not byte for byte"
            RESULTS.append((f"{TOOL} {label} level {level}", len(ours), len(ours), "bytes equal"))
            continue
        pb, ob = blocks(pout), blocks(oout)
        assert len(pb) == len(ob), label
        assert [p[1] for p in pb] == [o[1] for o in ob], f"{label}: inflated blocks differ"
        deflated = sum(p[0] != o[0] for p, o in zip(pb, ob))
        assert deflated <= gathered - 1, f"{label}: {deflated} blocks differ"
        RESULTS.append((f"{TOOL} {label} level 5", len(pb) - deflated, len(pb),
                        f"{deflated} blocks deflated again, inflate equal"))
    return records_of(oout)


@pytest.mark.parametrize("name", ["t1", "foundation", "exome"])
def test_pieces_gather_back(pieces, tmp_path, name):
    count = tier_a(pieces[name], tmp_path, f"{name} by contig")
    source = {"t1": T1}.get(name)
    if source is not None:
        assert count == records_of(source)


def test_pieces_in_another_order(pieces, tmp_path):
    parts = pieces["t1"]
    tier_a([parts[2], parts[0], parts[1]], tmp_path, "t1 reordered")


def test_records_compare_at_the_default_level(pieces):
    written = compare_bam(TOOL, pieces["foundation"], [], label="foundation pieces",
                          record_bytes=True)
    RESULTS.append((f"{TOOL} foundation records", len(written.records), len(written.records),
                    "records equal as bytes"))


@pytest.mark.parametrize("inputs", [
    ["markdup.bam"], ["markdup.bam", "markdup.bam"], ["insert_size.bam", "markdup.bam"],
    ["oxog.bam", "artifact.bam", "oxog.bam"], ["no_reads.bam", "insert_size.bam"],
    ["insert_size.bam", "no_reads.bam"], ["duplicate_marking.bam", "duplicate_marking.bam"],
])
def test_fixtures(tmp_path, inputs):
    tier_a([FIXTURES / name for name in inputs], tmp_path, " ".join(inputs))


def test_a_missing_end_of_file_block(pieces, tmp_path):
    parts = pieces["t1"]
    cut = tmp_path / "noeof.bam"
    cut.write_bytes(parts[1].read_bytes()[:-28])
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    tier_a([parts[0], cut], tmp_path / "a", "no EOF later")
    tier_a([cut, parts[0]], tmp_path / "b", "no EOF first")


def test_lists(pieces, tmp_path):
    parts = pieces["t1"]
    inner = tmp_path / "inner.txt"
    inner.write_text(f"{parts[2]}\r\n")
    outer = tmp_path / "parts.list"
    outer.write_text(f"   {parts[1]}  \n\n{inner}\n\t\n{parts[3]}\n")
    tier_a([parts[0], outer], tmp_path, "list of lists", gathered=4)


def test_index_and_md5(tmp_path):
    md = FIXTURES / "markdup.bam"
    index_and_md5(TOOL, [f"I={md}", f"I={md}"], tmp_path, "out.bai")


def test_the_index_of_a_dat_output(tmp_path):
    md = FIXTURES / "markdup.bam"
    (prc, pout, _), (orc, oout, _) = both([f"I={md}", "CREATE_INDEX=true"], tmp_path, "out.dat")
    assert prc == orc == 0
    names = [sorted(p.name for p in out.parent.iterdir()) for out in (pout, oout)]
    assert names[0] == names[1] == ["out.bai", "out.dat"]
    RESULTS.append((f"{TOOL} CREATE_INDEX out.dat", 1, 1, "out.bai in both"))


def test_the_index_follows_the_header_order_as_htsjdk(tmp_path):
    index_by_sort_order(TOOL, lambda bam: [f"I={bam}"], tmp_path, "out.bai")


def test_the_index_follows_the_first_header(tmp_path):
    """A queryname-sorted input after markdup.bam: the output takes the `SO:coordinate`
    header of the first, so both tools index the records (in no coordinate order)."""
    queryname = tmp_path / "queryname.bam"
    subprocess.run(["samtools", "sort", "-n", "--no-PG", "-o", str(queryname),
                    str(FIXTURES / "markdup.bam")], check=True)
    (prc, pout, perr), (orc, oout, oerr) = both(
        [f"I={FIXTURES / 'markdup.bam'}", f"I={queryname}", "CREATE_INDEX=true"], tmp_path)
    assert prc == 0, perr[-800:]
    assert orc == 0, oerr[-800:]
    indexes = [out.with_name("out.bai").read_bytes() for out in (pout, oout)]
    assert indexes[0] == indexes[1]
    RESULTS.append((f"{TOOL} CREATE_INDEX queryname after coordinate", 1, 1, "out.bai equal"))


def test_sam_and_cram_inputs_gather_record_by_record(pieces, tmp_path):
    parts = pieces["t1"]
    sam = tmp_path / "part.sam"
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(sam), str(parts[1])],
                   check=True)
    cram = tmp_path / "part.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(REF), "-o", str(cram),
                    str(parts[2])], check=True)
    silent = ["VALIDATION_STRINGENCY=SILENT"]
    for inputs, args, label in [
        ([parts[0], sam], silent, "BAM then SAM"),
        ([parts[0], sam, cram], [*silent, f"R={REF}"], "BAM, SAM, CRAM"),
        ([sam], silent, "SAM alone"),
    ]:
        written = compare_bam(TOOL, inputs, args, label=label, record_bytes=True)
        RESULTS.append((f"{TOOL} {label}", len(written.records), len(written.records),
                        "records equal as bytes"))
    written = compare_bam(TOOL, [parts[0], sam], silent, label="SAM text out", suffix=".sam")
    RESULTS.append((f"{TOOL} SAM text out", len(written.records), len(written.records),
                    "records equal"))


def error_inputs(pieces, directory):
    parts = pieces["t1"]
    cut = directory / "cut.bam"
    cut.write_bytes(parts[1].read_bytes()[:5000])
    text = directory / "text.bam"
    text.write_text("hello\n" * 10)
    bgzip_sam = directory / "bgzip_sam.bam"
    sam = subprocess.run(["samtools", "view", "--no-PG", "-h", str(parts[1])], check=True,
                         capture_output=True).stdout
    subprocess.run(["bgzip", "-c"], input=sam, stdout=bgzip_sam.open("wb"), check=True)
    cram = directory / "part.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(REF), "-o", str(cram),
                    str(parts[1])], check=True)
    missing_list = directory / "missing.list"
    missing_list.write_text(f"{directory / 'nope.bam'}\n")
    return {
        "cut": [parts[0], cut],
        "text": [parts[0], text],
        "bgzip SAM": [parts[0], bgzip_sam],
        "missing": [parts[0], directory / "nope.bam"],
        "missing in a list": [parts[0], missing_list],
        "CRAM without R": [parts[0], cram],
    }


@pytest.mark.parametrize("case", ["cut", "text", "bgzip SAM", "missing", "missing in a list",
                                  "CRAM without R"])
def test_errors_as_picard(pieces, tmp_path, case):
    inputs = error_inputs(pieces, tmp_path)[case]
    compare_bam(TOOL, inputs, ["VALIDATION_STRINGENCY=SILENT"], outcome="fail", label=case)
    RESULTS.append((f"{TOOL} error {case}", 0, 0, "both fail, no output"))


def test_a_missing_output_directory(tmp_path):
    md = FIXTURES / "markdup.bam"
    for side in ("picard", "bammetrics"):
        out = tmp_path / side / "missing" / "out.bam"
        done = run_tool(side, TOOL, [f"I={md}", f"O={out}"])
        assert done.returncode != 0, side
        assert not out.exists()


def test_bammetrics_refuses_what_picard_writes(pieces, tmp_path):
    """Picard writes these; bammetrics exits 1 on the Picard face and writes nothing
    (docs/07 decision rows)."""
    parts = pieces["t1"]
    md = FIXTURES / "markdup.bam"
    cases = {
        "another dictionary": ([f"I={md}", f"I={parts[0]}"], "out.bam", 1),
        "a first input with no record and another dictionary":
            ([f"I={FIXTURES / 'no_reads.bam'}", f"I={parts[0]}"], "out.bam", 1),
        "a .sam name": ([f"I={md}", f"I={md}"], "out.sam", 1),
    }
    for label, (args, name, code) in cases.items():
        directory = tmp_path / label.replace(" ", "_")
        directory.mkdir()
        (prc, _, perr), (orc, oout, _) = both(args, directory, name)
        assert prc == 0, f"{label}: {perr[-400:]}"
        assert orc == code and not oout.exists(), label
        RESULTS.append((f"{TOOL} refuses {label}", 1, 1, "Picard writes; bammetrics exits"))
    copy = tmp_path / "copy.bam"
    copy.write_bytes(md.read_bytes())
    done = subprocess.run([str(BAMMETRICS), TOOL, f"I={copy}", f"I={md}", f"O={copy}"],
                          capture_output=True, text=True)
    assert done.returncode == 1 and copy.read_bytes() == md.read_bytes()


def test_the_verb_and_a_list_on_standard_input(pieces, tmp_path):
    parts = pieces["t1"]
    out = tmp_path / "verb.bam"
    done = subprocess.run([str(BAMMETRICS), "gather-bam-files", str(parts[0]), "-i",
                           str(parts[1]), "--compression-level", "0", "-o", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    piped = tmp_path / "piped.bam"
    done = subprocess.run([str(BAMMETRICS), "gather-bam-files", "-", "-o", str(piped),
                           "--compression-level", "0"], input=f"{parts[0]}\n{parts[1]}\n",
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    assert out.read_bytes() == piped.read_bytes()

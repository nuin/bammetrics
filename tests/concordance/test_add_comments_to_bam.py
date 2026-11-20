"""AddCommentsToBam against Picard 3.5.0, and the helpers of the record-editing tools
(test_replace_sam_header.py imports them).

Picard block copies the BAM (htsjdk `BamFileIoUtils.reheaderBamFile`): a new header, then
the record bytes as they were. Every case goes through `compare.compare_bam` with
`record_bytes=True` (tier B: the header lines, `@PG` ID and PP, every record as samtools
prints it, then every record as bytes after inflation).

Inputs: the fixture BAMs, a BAM without an `@HD` line, a queryname-sorted BAM, the
foundation corpus, t1 to t3, the exome slice. Comments: plain, with a tab, with a carriage
return, with `@CO` and a tab in front (kept once), `@CO` alone, beyond ASCII (Latin-1
bytes in both). Options: COMPRESSION_LEVEL, CREATE_INDEX (`OUTPUT.bai`, also under
`SO:unsorted` and without `@HD`, none under another order), CREATE_MD5_FILE,
MAX_RECORDS_IN_RAM and REFERENCE_SEQUENCE (both ignored on a BAM). Error parity: a `.sam`
input, SAM text in another name, bgzip SAM, a CRAM, no COMMENT, a newline in a comment,
ASSUME_SORTED, a missing input and a missing output directory. The step in a chain against
Picard's file. The distinct-value gate wants three comment counts and three input sizes.
"""
import hashlib
import subprocess
from pathlib import Path

import pysam
import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_bam, diff_bam, diff_bam_bytes
from conftest import RESULTS
from run_picard import picard_env
from test_mark_duplicates import offsets_that_differ

TOOL = "AddCommentsToBam"
REF = FIXTURES / "ref_100k.fa"
EOF_BLOCK = 28

SOURCES = {
    "markdup": FIXTURES / "markdup.bam",
    "insert_size": FIXTURES / "insert_size.bam",
    "duplicate_marking": FIXTURES / "duplicate_marking.bam",
    "oxog": FIXTURES / "oxog.bam",
    "artifact": FIXTURES / "artifact.bam",
    "t1": CORPUS / "t1.bam",
    "t2": CORPUS / "t2.bam",
    "t3": CORPUS / "t3.bam",
}

# (comments added, records) of every case that compared equal: the gate's input.
SEEN = []


def source_path(name, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    return SOURCES[name]


def run_tool(side, tool, args, cwd=None):
    """One tool through Picard (side "picard") or bammetrics."""
    if side == "picard":
        env = picard_env()
        env.setdefault("JAVA_OPTS", "-Xmx8g")
        command = ["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"]
        return subprocess.run(command, env=env, capture_output=True, text=True, cwd=cwd)
    return subprocess.run([str(BAMMETRICS), tool, *args], capture_output=True, text=True, cwd=cwd)


def compare(tool, source, args, label, outcome="equal", reference=None, suffix=".bam"):
    """Both tools through compare_bam, records compared as bytes too."""
    written = compare_bam(tool, [source], args, args, outcome=outcome, label=label,
                          reference=reference, suffix=suffix, record_bytes=suffix == ".bam")
    detail = "both fail" if written.outcome == "fail" else f"{len(written.records)} records equal"
    RESULTS.append((f"{tool} {label}", len(written.records), len(written.records), detail))
    return written


def records_of(path):
    return int(subprocess.run(["samtools", "view", "-c", str(path)], check=True,
                              capture_output=True, text=True).stdout)


def index_and_md5(tool, args, tmp_path, index_name):
    """`tool` with CREATE_INDEX and CREATE_MD5_FILE through both tools, each into a
    directory of its own: the same file names in both, `index_name` among them. The
    bammetrics .bai is Picard's BuildBamIndex of the bammetrics BAM, byte for byte;
    Picard's own .bai equals its BuildBamIndex of its BAM except where the last record
    of the file ends (test_mark_duplicates.py). Each .md5 is the hex MD5 of its BAM."""
    names = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        done = run_tool(side, tool, [*args, f"O={directory / 'out.bam'}", "CREATE_INDEX=true",
                                     "CREATE_MD5_FILE=true"])
        assert done.returncode == 0, f"{side}: {done.stderr[-800:]}"
        names[side] = sorted(p.name for p in directory.iterdir())
        rebuilt = tmp_path / f"{side}.rebuilt.bai"
        index = run_tool("picard", "BuildBamIndex", [f"I={directory / 'out.bam'}", f"O={rebuilt}",
                                                     "VALIDATION_STRINGENCY=SILENT"])
        assert index.returncode == 0, index.stderr[-800:]
        written, bam = (directory / index_name).read_bytes(), (directory / "out.bam").read_bytes()
        if side == "bammetrics":
            assert written == rebuilt.read_bytes(), "the bammetrics .bai is not Picard's index"
        else:
            after_last_record = (len(bam) - EOF_BLOCK) << 16
            for picard_word, read_word in offsets_that_differ(written, rebuilt.read_bytes()):
                assert read_word == after_last_record and picard_word >> 16 < read_word >> 16
        md5 = hashlib.md5(bam).hexdigest()
        assert (directory / "out.bam.md5").read_text() == md5, side
    assert names["picard"] == names["bammetrics"] == sorted(["out.bam", index_name, "out.bam.md5"])
    RESULTS.append((f"{tool} CREATE_INDEX {' '.join(args)}", 1, 1,
                    f"{index_name} is Picard's index of the BAM; .md5 equal"))


def chain_against_picard(tool, step, source, picard_args, tmp_path, reference=None):
    """`chain SOURCE + STEP + write` against Picard's file of the tool on SOURCE: the
    same header lines, records and record bytes."""
    picard_out, ours_out = tmp_path / "picard.bam", tmp_path / "chained.bam"
    done = run_tool("picard", tool, [f"I={source}", f"O={picard_out}", *picard_args])
    assert done.returncode == 0, done.stderr[-800:]
    chain = subprocess.run([str(BAMMETRICS), "chain", str(source), "-o", str(tmp_path / "c"),
                            *(["-R", str(reference)] if reference else []), "+", *step, "+",
                            "write", str(ours_out)], capture_output=True, text=True)
    assert chain.returncode == 0, chain.stderr
    problems = diff_bam(picard_out, ours_out, reference)
    problems += [] if problems else diff_bam_bytes(picard_out, ours_out)
    assert not problems, "\n".join(problems[:20])
    count = records_of(picard_out)
    RESULTS.append((f"{tool} chain {' '.join(step)} on {Path(source).name}", count, count,
                    "chained equals Picard"))


COMMENTS = ["C=one", "C=two words", "COMMENT=@CO\tthree", "C=tab\there", "C=a\rb", "C=@CO",
            "C=é ü ñ"]


def sort_order_inputs(directory):
    """BAMs whose header order and record order disagree, for the index of a block copy:
    the markdup records under `SO:unsorted`, under no `@HD`, and sorted by name under
    `SO:queryname`; three records whose reference goes back (c2 then c1) under
    `SO:unsorted` and `SO:unknown`."""
    md = SOURCES["markdup"]
    header = subprocess.run(["samtools", "view", "-H", str(md)], check=True,
                            capture_output=True, text=True).stdout
    body = subprocess.run(["samtools", "view", str(md)], check=True, capture_output=True,
                          text=True).stdout
    no_hd = "".join(line for line in header.splitlines(True) if not line.startswith("@HD"))
    texts = {
        "unsorted_header": header.replace("SO:coordinate", "SO:unsorted") + body,
        "no_hd": no_hd + body,
    }
    back = "".join(f"r{i}\t0\t{ref}\t{pos}\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
                   for i, (ref, pos) in enumerate([("c2", 100), ("c1", 200), ("c1", 50)]))
    for order in ("unsorted", "unknown"):
        texts[f"back_{order}"] = (f"@HD\tVN:1.6\tSO:{order}\n@SQ\tSN:c1\tLN:5000\n"
                                  f"@SQ\tSN:c2\tLN:5000\n{back}")
    paths = {}
    for name, text in texts.items():
        paths[name] = directory / f"{name}.bam"
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(paths[name]), "-"],
                       input=text, text=True, check=True)
    paths["queryname"] = directory / "queryname.bam"
    subprocess.run(["samtools", "sort", "-n", "--no-PG", "-o", str(paths["queryname"]), str(md)],
                   check=True)
    return paths


def block_starts(bam):
    """The offset of each BGZF block of `bam`, in file order."""
    data, starts, at = bam.read_bytes(), [], 0
    while at < len(data):
        starts.append(at)
        at += int.from_bytes(data[at + 16:at + 18], "little") + 1
    return starts


def same_index_by_block(index, other, bam, other_bam):
    """Two .bai files equal once each virtual offset is read as (block number, offset
    in the block) of its own BAM. A block copy deflates the header block again, and
    zlib-rs and the JDK's zlib give it other lengths (docs/07 D63), so the compressed
    offsets of two indexes of the same records differ by that length."""
    index, other = index.read_bytes(), other.read_bytes()
    if len(index) != len(other):
        return False
    starts, other_starts = block_starts(bam), block_starts(other_bam)
    for word, other_word in offsets_that_differ(index, other):
        if word & 0xFFFF != other_word & 0xFFFF:
            return False
        if word >> 16 not in starts or other_word >> 16 not in other_starts:
            return False
        if starts.index(word >> 16) != other_starts.index(other_word >> 16):
            return False
    return True


def index_by_sort_order(tool, args_of, tmp_path, index_name):
    """CREATE_INDEX on a block copy follows htsjdk's `BAMIndexer`: the records are indexed
    under `SO:coordinate`, `SO:unsorted` or no `@HD` (the same .bai in both tools, block
    by block), no index under another order (both exit 0). A record the indexer refuses
    (its reference goes back) leaves Picard's run at exit 0 with an empty .bai; bammetrics
    exits 1 and writes nothing (KNOWN_DIFFERENCES.md)."""
    inputs = sort_order_inputs(tmp_path)
    for name, bam in inputs.items():
        result = {}
        for side in ("picard", "bammetrics"):
            directory = tmp_path / name / side
            directory.mkdir(parents=True)
            done = run_tool(side, tool, [*args_of(bam), f"O={directory / 'out.bam'}",
                                         "CREATE_INDEX=true"])
            result[side] = (done.returncode, sorted(p.name for p in directory.iterdir()))
        index = [tmp_path / name / side / index_name for side in ("picard", "bammetrics")]
        if name in ("unsorted_header", "no_hd"):
            assert result["picard"] == result["bammetrics"] == (0, sorted(["out.bam",
                                                                           index_name])), name
            bams = [path.parent / "out.bam" for path in index]
            assert same_index_by_block(*index, *bams), name
            detail = f"{index_name} equal by block"
        elif name in ("queryname", "back_unknown"):
            assert result["picard"] == result["bammetrics"] == (0, ["out.bam"]), name
            detail = "no index in either"
        else:
            assert result["picard"] == (0, sorted(["out.bam", index_name])), name
            assert index[0].stat().st_size == 0, name
            assert result["bammetrics"] == (1, []), name
            detail = "indexer refuses: Picard an empty .bai, bammetrics exit 1 and no file"
        RESULTS.append((f"{tool} CREATE_INDEX {name}", 1, 1, detail))


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    source = source_path(name, corpus, exome)
    for count in (1, 3, len(COMMENTS)):
        written = compare(TOOL, source, COMMENTS[:count], f"{name} {count} comments")
        SEEN.append((count, len(written.records)))


def header_only(tool, args, tmp_path):
    """`tool` on a BAM with a header and no record. Picard's block copy leaves out the
    BGZF end-of-file block (the input's last block is the one it skips as the
    terminator); bammetrics ends the file with it (KNOWN_DIFFERENCES.md). The header
    lines are equal and neither file holds a record."""
    outs = {}
    for side in ("picard", "bammetrics"):
        out = tmp_path / f"{side}.bam"
        done = run_tool(side, tool, [f"I={FIXTURES / 'no_reads.bam'}", f"O={out}", *args])
        assert done.returncode == 0, f"{side}: {done.stderr[-800:]}"
        outs[side] = out
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    assert not outs["picard"].read_bytes().endswith(eof)
    assert outs["bammetrics"].read_bytes().endswith(eof)
    headers = [subprocess.run(["samtools", "view", "--no-PG", "-H", str(out)], check=True,
                              capture_output=True, text=True).stdout for out in outs.values()]
    assert headers[0] == headers[1]
    assert [records_of(out) for out in outs.values()] == [0, 0]
    RESULTS.append((f"{tool} no_reads", 0, 0, "same header; Picard leaves out the EOF block"))


def test_a_bam_without_records(tmp_path):
    header_only(TOOL, ["C=x"], tmp_path)


def test_a_bam_without_hd_gets_htsjdk_hd(tmp_path):
    source = tmp_path / "no_hd.bam"
    text = subprocess.run(["samtools", "view", "--no-PG", "-h", str(SOURCES["markdup"])],
                          check=True, capture_output=True, text=True).stdout
    sam = tmp_path / "no_hd.sam"
    sam.write_text("".join(line for line in text.splitlines(True) if not line.startswith("@HD")))
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(source), str(sam)], check=True)
    compare(TOOL, source, ["C=x"], "a BAM without @HD")


def test_no_comment_copies_the_bam():
    """COMMENT is listed as required, and Picard 3.5.0 runs without one: the BAM is
    copied with its header as it is."""
    written = compare(TOOL, SOURCES["markdup"], [], "markdup without COMMENT")
    SEEN.append((0, len(written.records)))


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=1"], [f"R={REF}"]])
def test_options(args):
    compare(TOOL, SOURCES["markdup"], ["C=x", *args], f"markdup {' '.join(args)}")


def test_index_and_md5(tmp_path):
    index_and_md5(TOOL, [f"I={SOURCES['markdup']}", "C=x"], tmp_path, "out.bam.bai")


def test_no_index_beside_a_queryname_bam(tmp_path):
    """A queryname-sorted BAM gets no index from either tool (Picard's indexing thread
    fails and the run succeeds); the .md5 is written."""
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        done = run_tool(side, TOOL, [f"I={SOURCES['duplicate_marking']}", "C=x",
                                     f"O={directory / 'out.bam'}", "CREATE_INDEX=true",
                                     "CREATE_MD5_FILE=true"])
        assert done.returncode == 0, done.stderr[-800:]
        assert sorted(p.name for p in directory.iterdir()) == ["out.bam", "out.bam.md5"], side
    RESULTS.append((f"{TOOL} CREATE_INDEX queryname", 0, 0, "neither writes an index"))


def test_the_index_follows_the_header_order_as_htsjdk(tmp_path):
    index_by_sort_order(TOOL, lambda bam: [f"I={bam}", "C=x"], tmp_path, "out.bam.bai")


@pytest.fixture(scope="module")
def other_kinds(tmp_path_factory):
    """The markdup records as SAM text named .txt and as bgzip SAM."""
    directory = tmp_path_factory.mktemp("kinds")
    text = subprocess.run(["samtools", "view", "-h", str(SOURCES["markdup"])], check=True,
                          capture_output=True, text=True).stdout
    plain = directory / "markdup.txt"
    plain.write_text(text)
    zipped = directory / "markdup.sam.gz"
    pysam.tabix_compress(str(plain), str(zipped), force=True)
    return {"sam_txt": plain, "bgzip_sam": zipped}


@pytest.mark.parametrize("case", ["sam", "sam_txt", "bgzip_sam", "cram", "newline",
                                  "assume_sorted", "missing_input", "missing_directory"])
def test_errors_as_picard(other_kinds, tmp_path, case):
    source, args = SOURCES["markdup"], ["C=x"]
    if case == "sam":
        source = FIXTURES / "markdup.sam"
    elif case in other_kinds:
        source = other_kinds[case]
    elif case == "cram":
        source, args = FIXTURES / "insert_size.cram", ["C=x", f"R={REF}"]
    elif case == "newline":
        args = ["C=a\nb"]
    elif case == "assume_sorted":
        args = ["C=x", "ASSUME_SORTED=true"]
    elif case == "missing_input":
        source = tmp_path / "nope.bam"
    if case == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"O={tmp_path / 'no' / 'out.bam'}", "C=x"])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    compare(TOOL, source, args, f"error {case}", "fail")


def test_the_step_in_a_chain_equals_picard(corpus, tmp_path):
    for source in (SOURCES["markdup"], corpus / "corpus.bam"):
        directory = tmp_path / source.stem
        directory.mkdir()
        chain_against_picard(TOOL, ["add-comments-to-bam", "-C", "x", "-C", "@CO\ty"], source,
                             ["C=x", "C=@CO\ty"], directory)


def test_values_are_distinct():
    """The cases add 1, 3 and 8 comments to inputs of at least three sizes."""
    assert len({count for count, _ in SEEN}) >= 3, SEEN
    assert len({records for _, records in SEEN}) >= 3, SEEN

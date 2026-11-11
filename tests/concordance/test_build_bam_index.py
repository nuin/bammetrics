"""BuildBamIndex against Picard 3.5.0: the .bai, byte for byte, for the same BAM bytes.

Inputs: the six coordinate-sorted BAMs of tests/data, the foundation corpus, every
coordinate-sorted corpus tier BAM (t0 to t3), the exome slice, and the hand-written BAMs
of `index_bams.py` (BGZF blocks of 250 to 5000 bytes, records that span blocks, bins
met again in later blocks, a wrong bin field, unplaced reads only, no record, no
reference, references without reads, a reference longer than 512 Mbp with reads up to
position 2^29, a file without its EOF block, a reference of length 0).

Each index is also read by samtools: `samtools idxstats` and region counts through the
bammetrics index equal the ones through `samtools index`.

Error parity (`index_bams.refused`): a header that is not SO:coordinate (queryname,
unsorted, unknown, no SO, no @HD), references out of order, an alignment past 2^29, a
record outside its reference or the header, a truncated file, a record whose block_size
claims 4 GiB (bammetrics reads the bytes the file holds); a SAM, a CRAM, a missing
file, an OUTPUT whose directory is missing. Picard leaves an empty OUTPUT behind;
bammetrics leaves none and keeps an earlier file.

Where bammetrics differs on purpose (KNOWN_DIFFERENCES.md): records are not validated,
so the three `index_bams.lenient` BAMs are indexed as Picard indexes them under
VALIDATION_STRINGENCY=SILENT; an empty BGZF block in the middle of the file does not end
the read.

`write --index` of a chain and `MarkDuplicates CREATE_INDEX=true` write a BAM with the
bammetrics deflater: the .bai beside it equals the one Picard builds from that BAM.

`test_values_are_distinct` requires three distinct values of each count an index holds
and a bin of every level over the cases compared.
"""
import hashlib
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

import index_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "BuildBamIndex"
VERB = "build-bam-index"
# tests/data BAMs whose header says SO:coordinate.
SORTED_FIXTURES = ["artifact", "insert_size", "markdup", "no_reads", "orientation", "oxog"]
# The first bin of each level of the binning scheme, deepest last.
LEVEL_STARTS = (0, 1, 9, 73, 585, 4681)
METADATA_BIN = 37450
# The corpora hold records Picard's default STRICT validation refuses (a supplementary
# flag on an unmapped read), so both tools get SILENT, as in every other module.
SILENT = "VALIDATION_STRINGENCY=SILENT"
# Every Picard index compared, parsed: label -> index.
SEEN = {}


def env():
    e = picard_env()
    e.setdefault("JAVA_OPTS", "-Xmx8g")
    return e


def run_picard(args, quiet=True):
    extra = ["VERBOSITY=ERROR", "QUIET=true"] if quiet else []
    return subprocess.run(["picard", TOOL, *args, *extra], env=env(), capture_output=True)


def run_ours(args, face=TOOL, stdin=None):
    return subprocess.run([str(BAMMETRICS), face, *args], capture_output=True, stdin=stdin)


def parse(data):
    """A .bai as {"references": [{"bins": {bin: [(start, end)]}, "metadata": (first,
    last, aligned, unaligned) or None, "intervals": [offset]}], "no_coordinate": n}."""
    assert data[:4] == b"BAI\x01", data[:4]
    (n_ref,) = struct.unpack_from("<i", data, 4)
    at = 8
    references = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", data, at)
        at += 4
        bins, metadata = {}, None
        for _ in range(n_bin):
            number, n_chunk = struct.unpack_from("<Ii", data, at)
            at += 8
            chunks = [struct.unpack_from("<QQ", data, at + 16 * i) for i in range(n_chunk)]
            at += 16 * n_chunk
            if number == METADATA_BIN:
                metadata = (*chunks[0], *chunks[1])
            else:
                bins[number] = chunks
        (n_intv,) = struct.unpack_from("<i", data, at)
        at += 4
        intervals = list(struct.unpack_from(f"<{n_intv}Q", data, at))
        at += 8 * n_intv
        references.append({"bins": bins, "metadata": metadata, "intervals": intervals})
    (no_coordinate,) = struct.unpack_from("<Q", data, at)
    assert at + 8 == len(data)
    return {"references": references, "no_coordinate": no_coordinate}


def first_difference(a, b):
    """Where two indexes differ, for the assertion message."""
    if len(a) != len(b):
        return f"lengths {len(a)} and {len(b)}"
    at = next(i for i in range(len(a)) if a[i] != b[i])
    return f"byte {at}: {a[at:at + 16].hex()} and {b[at:at + 16].hex()}"


def compare(bam, label, tmp_path, extra=()):
    """Both tools on `bam` with `I=` and `O=`: exit 0 and the same bytes."""
    picard_bai, ours_bai = tmp_path / "picard.bai", tmp_path / "ours.bai"
    picard = run_picard([f"I={bam}", f"O={picard_bai}", *extra])
    ours = run_ours([f"I={bam}", f"O={ours_bai}", *extra])
    assert picard.returncode == 0, picard.stderr[-800:]
    assert ours.returncode == 0, ours.stderr
    assert ours.stdout == b"" and picard.stdout == b""
    expected, got = picard_bai.read_bytes(), ours_bai.read_bytes()
    assert got == expected, f"{TOOL} {label}: {first_difference(got, expected)}"
    SEEN[label] = parse(expected)
    RESULTS.append((f"{TOOL} {label}", len(expected), len(expected), ".bai identical"))
    return ours_bai


def beside(bam, index, directory):
    """`bam` linked into `directory` as x.bam with `index` beside it as x.bam.bai."""
    directory.mkdir()
    link = directory / "x.bam"
    link.symlink_to(Path(bam).resolve())
    shutil.copy(index, directory / "x.bam.bai")
    return link


def samtools(*args):
    return subprocess.run(["samtools", *args], check=True, capture_output=True).stdout


def regions(bam):
    """A few regions of each reference of `bam`: whole, a start, a middle window."""
    out = []
    for line in samtools("view", "-H", str(bam)).decode().splitlines():
        if not line.startswith("@SQ"):
            continue
        fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
        name, length = fields["SN"], int(fields["LN"])
        middle = length // 2
        out += [name, f"{name}:1-{min(length, 20000)}", f"{name}:{middle}-{min(length, middle + 40000)}",
                f"{name}:{max(1, length - 3000)}-{length}"]
    return out


def cross_check(bam, index, tmp_path, label):
    """samtools reads the bammetrics index of `bam` as it reads its own: the same
    idxstats rows and the same record count for every region tried."""
    ours = beside(bam, index, tmp_path / "with_ours")
    theirs_dir = tmp_path / "with_samtools"
    theirs_dir.mkdir()
    theirs = theirs_dir / "x.bam"
    theirs.symlink_to(Path(bam).resolve())
    samtools("index", str(theirs))
    assert samtools("idxstats", str(ours)) == samtools("idxstats", str(theirs)), label
    tried = regions(bam)
    for region in tried:
        assert samtools("view", "-c", str(ours), region) == samtools("view", "-c", str(theirs), region), (
            f"{label}: {region}")
    RESULTS.append((f"{TOOL} {label} samtools", len(tried), len(tried), "idxstats and region counts equal"))


def coordinate_sorted(bam):
    head = samtools("view", "-H", str(bam)).split(b"\n", 1)[0]
    return b"SO:coordinate" in head


def tier_bams():
    return sorted(p for p in CORPUS.glob("t*.bam") if coordinate_sorted(p))


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """The hand-written BAMs of index_bams.py, one file per case."""
    directory = tmp_path_factory.mktemp("index_bams")
    index_bams.write_all(directory)
    return directory


@pytest.mark.parametrize("name", SORTED_FIXTURES)
def test_fixtures(tmp_path, name):
    """tests/data/bai holds Picard's files for `cargo test`: they are still Picard's."""
    index = compare(FIXTURES / f"{name}.bam", name, tmp_path)
    assert index.read_bytes() == (FIXTURES / "bai" / f"{name}.bai").read_bytes(), (
        "tests/data/bai is stale: python3 tests/data/make_bai.py")
    cross_check(FIXTURES / f"{name}.bam", index, tmp_path, name)


def test_foundation(corpus, tmp_path):
    index = compare(corpus / "corpus.bam", "foundation", tmp_path, extra=[SILENT])
    cross_check(corpus / "corpus.bam", index, tmp_path, "foundation")


def test_exome(exome, tmp_path):
    index = compare(exome / "exome.bam", "exome", tmp_path, extra=[SILENT])
    cross_check(exome / "exome.bam", index, tmp_path, "exome")


@pytest.mark.parametrize("bam", tier_bams(), ids=lambda p: p.stem)
def test_tiers(bam, tmp_path):
    index = compare(bam, bam.stem, tmp_path, extra=[SILENT])
    cross_check(bam, index, tmp_path, bam.stem)


def test_tiers_t1_to_t3_are_present():
    names = {p.stem for p in tier_bams()}
    assert {"t1", "t2", "t3"} <= names, f"make corpus first: {sorted(names)}"


@pytest.mark.parametrize("name", sorted(index_bams.indexed()))
def test_hand_written(built, tmp_path, name):
    compare(built / f"{name}.bam", name, tmp_path)


@pytest.mark.parametrize("name", ["tiny_blocks", "interleaved_bins", "window_edges", "placed_unmapped",
                                  "with_unplaced_tail", "unplaced_only", "no_records",
                                  "middle_reference_only"])
def test_hand_written_indexes_work_in_samtools(built, tmp_path, name):
    """Region counts need the bin field samtools filters by to be right, so the cases
    with a wrong stored bin are left out."""
    bam = built / f"{name}.bam"
    index = tmp_path / "ours.bai"
    assert run_ours([f"I={bam}", f"O={index}"]).returncode == 0
    cross_check(bam, index, tmp_path, name)


def test_samtools_reads_the_index_beside_an_empty_reference(built, tmp_path):
    """samtools 1.x prints no header with `LN:0` and finds no region of such a BAM, so
    `cross_check` cannot run. `samtools idxstats` reads the bammetrics index: the row of
    the empty reference, the counts of the two others and the unplaced reads."""
    bam = built / "empty_reference_between.bam"
    index = tmp_path / "ours.bai"
    assert run_ours([f"I={bam}", f"O={index}"]).returncode == 0
    ours = beside(bam, index, tmp_path / "with_ours")
    assert samtools("idxstats", str(ours)).decode().splitlines() == [
        "chr1\t70000\t200\t0", "zero\t0\t0\t0", "chr3\t50000\t200\t0", "*\t0\t0\t5"]
    RESULTS.append((f"{TOOL} empty_reference_between samtools", 4, 4, "idxstats rows as written"))


def test_long_reference_reaches_the_last_bin_and_window(built, tmp_path):
    """A reference of 2^29 + 40000 bases: reads up to position 2^29 are indexed, into
    the last window bin and the last of 32768 windows."""
    compare(built / "long_reference.bam", "long_reference again", tmp_path)
    big = SEEN["long_reference again"]["references"][0]
    assert max(big["bins"]) == 37448
    assert len(big["intervals"]) == 32768
    assert {0, 1, 9, 73, 585} <= set(big["bins"])


@pytest.mark.parametrize("name", sorted(index_bams.refused()))
def test_refused_in_both(built, tmp_path, name):
    """Both tools fail. Picard leaves an empty OUTPUT, or none when it fails on the
    header; bammetrics writes nothing and keeps the file that was there."""
    _, picard_says, ours_says = index_bams.refused()[name]
    bam = built / f"{name}.bam"
    picard_bai, ours_bai = tmp_path / "picard.bai", tmp_path / "ours.bai"
    ours_bai.write_bytes(b"an earlier index")
    picard = run_picard([f"I={bam}", f"O={picard_bai}"])
    ours = run_ours([f"I={bam}", f"O={ours_bai}"])
    assert picard.returncode != 0 and picard_says.encode() in picard.stderr, picard.stderr[-800:]
    assert ours.returncode == 1 and ours_says.encode() in ours.stderr, ours.stderr
    assert not picard_bai.exists() or picard_bai.stat().st_size == 0
    assert ours_bai.read_bytes() == b"an earlier index"
    assert sorted(p.name for p in tmp_path.iterdir()) == sorted(
        ["ours.bai"] + (["picard.bai"] if picard_bai.exists() else []))
    RESULTS.append((f"{TOOL} {name}", 0, 0, "both exit non-zero"))


def test_a_block_size_past_the_file_costs_the_bytes_present(built, tmp_path):
    """bammetrics alone: a record whose block_size claims 0xFFFFFF00 bytes in a BAM of
    227 is read as the bytes the file holds. A buffer of block_size bytes took 4.3 GB."""
    bam = built / "block_size_past_the_file.bam"
    for face, args in [(VERB, [str(bam), "-o", str(tmp_path / "x.bai")]),
                       ("mark-duplicates", [str(bam), "-o", str(tmp_path / "x.bam"), "--metrics-file",
                                            str(tmp_path / "x.txt")]),
                       ("chain", [str(bam), "+", "write", str(tmp_path / "y.bam"), "--index"])]:
        process = subprocess.Popen([str(BAMMETRICS), face, *args], stderr=subprocess.DEVNULL)
        _, status, usage = os.wait4(process.pid, 0)
        process.returncode = os.waitstatus_to_exitcode(status)
        # ru_maxrss is bytes on macOS and kilobytes on Linux.
        peak = usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)
        assert process.returncode == 1, face
        assert peak < 64 << 20, f"{face}: peak RSS {peak} bytes"
    assert not any(tmp_path.iterdir())
    RESULTS.append((f"{TOOL} block_size past the file", 3, 3, "three commands under 64 MiB"))


@pytest.mark.parametrize("name", ["position_without_reference", "start_past_the_reference"])
def test_refused_by_the_indexer_under_silent(built, tmp_path, name):
    """With validation off Picard's indexer fails where the bammetrics one does."""
    bam = built / f"{name}.bam"
    picard = run_picard([f"I={bam}", f"O={tmp_path / 'picard.bai'}", "VALIDATION_STRINGENCY=SILENT"])
    ours = run_ours([f"I={bam}", f"O={tmp_path / 'ours.bai'}", "VALIDATION_STRINGENCY=SILENT"])
    assert picard.returncode != 0 and b"Exception creating BAM index for record" in picard.stderr
    assert ours.returncode == 1 and b"Exception creating BAM index for record" in ours.stderr
    assert not (tmp_path / "ours.bai").exists()
    RESULTS.append((f"{TOOL} {name} SILENT", 0, 0, "both exit non-zero in the indexer"))


@pytest.mark.parametrize("name", sorted(index_bams.lenient()))
def test_unvalidated_records_index_as_picard_silent(built, tmp_path, name):
    """Picard STRICT refuses the record. bammetrics validates nothing and writes the
    index Picard writes under VALIDATION_STRINGENCY=SILENT."""
    bam = built / f"{name}.bam"
    strict = run_picard([f"I={bam}", f"O={tmp_path / 'strict.bai'}"])
    assert strict.returncode != 0 and b"SAM validation error" in strict.stderr
    ours_default = tmp_path / "default.bai"
    assert run_ours([f"I={bam}", f"O={ours_default}"]).returncode == 0
    index = compare(bam, f"{name} SILENT", tmp_path, extra=["VALIDATION_STRINGENCY=SILENT"])
    assert index.read_bytes() == ours_default.read_bytes()


def test_empty_block_in_the_middle_does_not_end_the_read(built, tmp_path):
    """htsjdk reads an empty BGZF block as the end of the file: Picard indexes the
    records before it and exits 0. bammetrics indexes every record, as samtools does."""
    _, records, read_by_htsjdk = index_bams.EMPTY_BLOCK_IN_THE_MIDDLE
    bam = built / "empty_block_in_the_middle.bam"
    picard_bai, ours_bai = tmp_path / "picard.bai", tmp_path / "ours.bai"
    assert run_picard([f"I={bam}", f"O={picard_bai}"]).returncode == 0
    assert run_ours([f"I={bam}", f"O={ours_bai}"]).returncode == 0
    aligned = lambda path: parse(path.read_bytes())["references"][0]["metadata"][2]  # noqa: E731
    assert aligned(picard_bai) == read_by_htsjdk
    assert aligned(ours_bai) == records == int(samtools("view", "-c", str(bam)))
    cross_check(bam, ours_bai, tmp_path, "empty block in the middle")


def test_default_output_name(corpus, tmp_path):
    """Without OUTPUT Picard 3.5.0 writes into the current directory, under the file
    name of INPUT with .bai in place of .bam, or appended to another name. Its help
    says beside the BAM. The Picard face of bammetrics does what Picard does."""
    names = [("x.bam", "x.bai"), ("a.b.bam", "a.b.bai"), ("plain", "plain.bai"), ("x.BAM", "x.BAM.bai"),
             ("sp ace.bam", "sp ace.bai")]
    # Numbered directories: x.bam and x.BAM are one name on a case-insensitive disk.
    for number, (name, index) in enumerate(names):
        written = {}
        for side, command in [("picard", ["picard", TOOL]), ("ours", [str(BAMMETRICS), TOOL])]:
            data, cwd = tmp_path / f"{side}_{number}_data", tmp_path / f"{side}_{number}_cwd"
            data.mkdir()
            cwd.mkdir()
            shutil.copy(corpus / "corpus.bam", data / name)
            result = subprocess.run([*command, f"I={data / name}", SILENT, "QUIET=true", "VERBOSITY=ERROR"],
                                    env=env(), cwd=cwd, capture_output=True)
            assert result.returncode == 0, result.stderr[-800:]
            assert [p.name for p in data.iterdir()] == [name], side
            assert [p.name for p in cwd.iterdir()] == [index], side
            written[side] = (cwd / index).read_bytes()
        assert written["ours"] == written["picard"], name
    RESULTS.append((f"{TOOL} default name", 5, 5, "same file name, directory and bytes"))


def test_verb_forms(corpus, tmp_path):
    """`bammetrics build-bam-index BAM` writes beside the BAM, `-o FILE` to the file,
    `-o -` to standard output, and `-` reads the BAM from standard input."""
    expected = tmp_path / "picard.bai"
    assert run_picard([f"I={corpus / 'corpus.bam'}", f"O={expected}", SILENT]).returncode == 0
    expected = expected.read_bytes()
    bam = tmp_path / "x.bam"
    shutil.copy(corpus / "corpus.bam", bam)
    assert run_ours([str(bam)], face=VERB).returncode == 0
    assert (tmp_path / "x.bai").read_bytes() == expected
    out = tmp_path / "sub" / "named.bai"
    assert run_ours([str(bam), "-o", str(out)], face=VERB).returncode == 0
    assert out.read_bytes() == expected
    printed = run_ours([str(bam), "-o", "-"], face=VERB)
    assert printed.returncode == 0 and printed.stdout == expected
    with open(bam, "rb") as stdin:
        piped = run_ours(["-", "-o", "-"], face=VERB, stdin=stdin)
    assert piped.returncode == 0 and piped.stdout == expected, piped.stderr
    with open(bam, "rb") as stdin:
        unnamed = run_ours(["-"], face=VERB, stdin=stdin)
    assert unnamed.returncode == 2 and b"missing -o FILE" in unnamed.stderr
    RESULTS.append((f"{VERB} verb forms", 4, 4, "default name, -o, standard output and input identical"))


def test_options_change_nothing(corpus, tmp_path):
    """Every option Picard lists beside INPUT and OUTPUT: both tools accept it and the
    index is the same. CREATE_MD5_FILE writes no digest of an index in either."""
    options = [SILENT, "COMPRESSION_LEVEL=1", "CREATE_INDEX=true",
               "CREATE_MD5_FILE=true", "MAX_RECORDS_IN_RAM=10", f"TMP_DIR={tmp_path}",
               "USE_JDK_DEFLATER=true", "USE_JDK_INFLATER=true", f"R={corpus / 'ref.fa'}"]
    plain = compare(corpus / "corpus.bam", "foundation plain", tmp_path, extra=[SILENT]).read_bytes()
    directory = tmp_path / "options"
    directory.mkdir()
    index = compare(corpus / "corpus.bam", "foundation every option", directory, extra=options)
    assert index.read_bytes() == plain
    assert sorted(p.name for p in directory.iterdir()) == ["ours.bai", "picard.bai"]


@pytest.mark.parametrize("option", ["ASSUME_SORTED=true", "GA4GH_CLIENT_SECRETS=x", "STOP_AFTER=10",
                                    "SORT_ORDER=coordinate"])
def test_options_picard_does_not_have(corpus, tmp_path, option):
    args = [f"I={corpus / 'corpus.bam'}", f"O={tmp_path / 'x.bai'}", option]
    picard, ours = run_picard(args), run_ours(args)
    assert picard.returncode != 0 and ours.returncode != 0
    assert not (tmp_path / "x.bai").exists()
    RESULTS.append((f"{TOOL} {option}", 0, 0, "both exit non-zero"))


def test_inputs_that_are_no_bam_fail_in_both(corpus, tmp_path):
    """A SAM, a CRAM, a missing file and no INPUT at all."""
    bam, ref = corpus / "corpus.bam", corpus / "ref.fa"
    sam, cram = tmp_path / "x.sam", tmp_path / "x.cram"
    sam.write_bytes(samtools("view", "-h", str(bam)))
    subprocess.run(["samtools", "view", "-C", "-T", str(ref), "-o", str(cram), str(bam)], check=True)
    for label, args in [("SAM", [f"I={sam}"]), ("CRAM", [f"I={cram}", f"R={ref}"]),
                        ("missing", [f"I={tmp_path / 'none.bam'}"]), ("no INPUT", [])]:
        out = tmp_path / f"{label.replace(' ', '_')}.bai"
        picard, ours = run_picard([*args, f"O={out}"]), run_ours([*args, f"O={out}"])
        assert picard.returncode != 0, label
        assert ours.returncode != 0, label
        assert not out.exists() or out.stat().st_size == 0, label
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))
    picard_sam = run_picard([f"I={sam}", f"O={tmp_path / 'p.bai'}"])
    assert b"Input file must be bam file, not sam file" in picard_sam.stderr
    ours_cram = run_ours([f"I={cram}", f"O={tmp_path / 'o.bai'}"])
    assert b"Input file must be bam file, not sam file" in ours_cram.stderr


def test_unwritable_output_fails_in_both(corpus, tmp_path):
    bam = corpus / "corpus.bam"
    directory = tmp_path / "a_directory"
    directory.mkdir()
    for label, out in [("missing directory", tmp_path / "nowhere" / "x.bai"), ("a directory", directory)]:
        args = [f"I={bam}", f"O={out}", SILENT]
        picard, ours = run_picard(args), run_ours(args)
        assert picard.returncode != 0 and ours.returncode != 0, label
        RESULTS.append((f"{TOOL} OUTPUT {label}", 0, 0, "both exit non-zero"))
    assert not (tmp_path / "nowhere").exists() and not any(directory.iterdir())


def test_output_is_never_an_input(corpus, tmp_path):
    """bammetrics alone: OUTPUT naming INPUT, by name or by a hard link, is refused and
    the BAM is left whole."""
    bam = tmp_path / "x.bam"
    shutil.copy(corpus / "corpus.bam", bam)
    before = hashlib.md5(bam.read_bytes()).hexdigest()
    link = tmp_path / "link.bai"
    link.hardlink_to(bam)
    for out in (bam, link):
        ours = run_ours([f"I={bam}", f"O={out}"])
        assert ours.returncode != 0 and b"would destroy it" in ours.stderr
        assert hashlib.md5(bam.read_bytes()).hexdigest() == before


def chain_sources(corpus, exome):
    return {"foundation": corpus / "corpus.bam", "exome": exome / "exome.bam", "t1": CORPUS / "t1.bam",
            "t2": CORPUS / "t2.bam", "t3": CORPUS / "t3.bam"}


@pytest.mark.parametrize("source", ["foundation", "exome", "t1", "t2", "t3"])
@pytest.mark.parametrize("level", [1, 5])
def test_chain_write_index_is_the_index_of_the_written_bam(corpus, exome, tmp_path, source, level):
    """`chain BAM + write OUT --index` compresses with the bammetrics deflater, so the
    BAM bytes are its own: the .bai beside them is the one Picard builds from them."""
    bam = chain_sources(corpus, exome)[source]
    out = tmp_path / "written.bam"
    chain = run_ours([str(bam), "+", "write", str(out), "--index", "--md5",
                      "--compression-level", str(level)], face="chain")
    assert chain.returncode == 0, chain.stderr
    index = tmp_path / "written.bai"
    picard_bai = tmp_path / "picard.bai"
    picard = run_picard([f"I={out}", f"O={picard_bai}", SILENT])
    assert picard.returncode == 0, picard.stderr[-800:]
    expected = picard_bai.read_bytes()
    assert index.read_bytes() == expected, first_difference(index.read_bytes(), expected)
    assert (tmp_path / "written.bam.md5").read_text() == hashlib.md5(out.read_bytes()).hexdigest()
    label = f"chain write --index {source} level {level}"
    SEEN[label] = parse(expected)
    RESULTS.append((label, len(expected), len(expected), ".bai identical to Picard's of the written BAM"))


def test_values_are_distinct():
    """Three distinct values of every count an index holds, a bin of every level, bins
    with several chunks, and a reference with reads beside one without."""
    assert len(SEEN) >= 40, f"run the whole module: {len(SEEN)} indexes compared"
    references = [r for index in SEEN.values() for r in index["references"]]
    with_reads = [r for r in references if r["metadata"]]
    assert any(not r["bins"] for r in references) and with_reads
    counts = {
        "references": {len(index["references"]) for index in SEEN.values()},
        "no_coordinate": {index["no_coordinate"] for index in SEEN.values()},
        "aligned": {r["metadata"][2] for r in with_reads},
        "unaligned": {r["metadata"][3] for r in with_reads},
        "bins": {len(r["bins"]) for r in with_reads},
        "intervals": {len(r["intervals"]) for r in with_reads},
        "chunks in a bin": {len(chunks) for r in with_reads for chunks in r["bins"].values()},
    }
    for name, values in counts.items():
        assert len(values - {0}) >= 3, f"{name}: only {sorted(values)}"
    levels = set()
    for reference in with_reads:
        for number in reference["bins"]:
            levels.add(max(i for i, start in enumerate(LEVEL_STARTS) if number >= start))
    assert levels == set(range(6)), f"bin levels seen: {sorted(levels)}"
    multi = sum(1 for r in with_reads for chunks in r["bins"].values() if len(chunks) > 1)
    assert multi >= 50, f"only {multi} bins with more than one chunk"
    RESULTS.append((f"{TOOL} distinct values", len(SEEN), len(SEEN), f"{multi} bins of several chunks"))

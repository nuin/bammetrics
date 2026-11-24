"""SplitSamByLibrary against Picard 3.5.0, and the helpers of the split tools
(test_split_sam_by_number_of_reads.py imports them).

Picard opens one writer per library of the header before the first record, a writer for
`unknown` with the first record that has no library, and writes through presorted
`SAMFileWriter`s. Every case runs both tools into a directory of its own and requires the
same list of file names, then compares each file as `compare.compare_bam` compares one
(tier B with `record_bytes=True`: the header lines, the `@PG` IDs and PP, every record as
samtools prints it, then every record of a BAM as bytes after inflation). An index is
compared through `samtools idxstats` (Picard's writer-side `.bai` ends one offset off its
BuildBamIndex, test_mark_duplicates.py); an `.md5` must hold the MD5 of its own file.

Inputs: `split_bams.py` (300 libraries, odd library names, read groups without `LB`,
records without `RG` or with one the header lacks, headers without libraries, `VN:1.4`, SAM text with libraries above U+00FF,
queryname and unsorted headers), the fixture BAMs, the foundation corpus (four
libraries), t1 to t3, the exome slice, SAM text (plain and gzip). Picard reads the
corpora at SILENT: its record validation refuses them at STRICT (docs/07 D51). Options:
COMPRESSION_LEVEL, CREATE_INDEX, CREATE_MD5_FILE, MAX_RECORDS_IN_RAM, REFERENCE_SEQUENCE,
VALIDATION_STRINGENCY, the default OUTPUT. Error parity: no library (exit 2 in both),
records out of order, an `RG:i` tag, ASSUME_SORTED, a missing input and a missing output
directory. Departures, asserted as such: a CRAM input (Picard writes CRAM, bammetrics exits
2 saying so) and two libraries with one file name (Picard writes one over the other and
exits 0, bammetrics exits 1 and leaves no file). The chain against Picard's MarkDuplicates
then SplitSamByLibrary. The distinct-value gate wants three library counts and three input
sizes.
"""
import hashlib
import os
import subprocess
import tempfile
from pathlib import Path

import pytest

import split_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import diff_bam_bytes, diff_sam_lines, missing_eofs, sam_view
from conftest import RESULTS
from test_add_comments_to_bam import run_tool

TOOL = "SplitSamByLibrary"
SILENT = ["VALIDATION_STRINGENCY=SILENT"]
REF = FIXTURES / "ref_100k.fa"

# (files written, records) of every case that compared equal: the gate's input.
SEEN = []

SOURCES = {
    "markdup": FIXTURES / "markdup.bam",
    "oxog": FIXTURES / "oxog.bam",
    "artifact": FIXTURES / "artifact.bam",
    "duplicate_marking": FIXTURES / "duplicate_marking.bam",
    "t1": CORPUS / "t1.bam",
    "t2": CORPUS / "t2.bam",
    "t3": CORPUS / "t3.bam",
}


@pytest.fixture(scope="session")
def split(tmp_path_factory):
    return split_bams.write(tmp_path_factory.mktemp("split"))


def source_path(name, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    return SOURCES[name]


def records_of(path):
    return int(subprocess.run(["samtools", "view", "-c", str(path)], check=True,
                              capture_output=True, text=True).stdout)


def idxstats(path):
    return subprocess.run(["samtools", "idxstats", str(path)], check=True, capture_output=True,
                          text=True).stdout


def indexed(index):
    """The BAM an index belongs to: `x.bai` beside `x.bam`, `x.cram.bai` beside `x.cram`."""
    stem = index.with_name(index.name[:-4])
    return stem.with_name(stem.name + ".bam") if stem.with_name(stem.name + ".bam").exists() else stem


def diff_file(picard, ours):
    """Difference lines between one file of each side."""
    name = picard.name
    if name.endswith(".md5"):
        own = hashlib.md5(ours.with_name(name[:-4]).read_bytes()).hexdigest()
        return [] if ours.read_text() == own else [f"{name}: not the MD5 of its file"]
    if name.endswith(".bai"):
        same = idxstats(indexed(picard)) == idxstats(indexed(ours))
        return [] if same else [f"{name}: idxstats differ"]
    problems = missing_eofs(picard, ours)
    problems += diff_sam_lines(sam_view(picard, True), sam_view(picard, False),
                               sam_view(ours, True), sam_view(ours, False))
    if name.endswith(".sam"):
        if not problems and picard.read_bytes() != ours.read_bytes():
            problems.append("SAM text differs byte for byte")
    elif not problems:
        problems += diff_bam_bytes(picard, ours)
    return [f"{name}: {p}" for p in problems]


def compare_split(tool, source, args, label, outcome="equal", workdir=None, cwd_output=False):
    """Run both tools on `source` with `O=<dir>/<side>` and compare every file written.
    `outcome` "equal": both succeed with the same files; "fail": both fail and bammetrics
    leaves its directory empty. Returns the file names."""
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        dirs = {side: Path(tmp) / side for side in ("picard", "bammetrics")}
        done = {}
        for side, out in dirs.items():
            out.mkdir()
            output = [] if cwd_output else [f"O={out}"]
            done[side] = run_tool(side, tool, [f"I={source}", *output, *args],
                                  cwd=out if cwd_output else None)
        what = f"{tool} {label}"
        if outcome == "fail":
            assert done["picard"].returncode != 0, f"{what}: Picard succeeded"
            assert done["bammetrics"].returncode != 0, f"{what}: bammetrics succeeded"
            assert not os.listdir(dirs["bammetrics"]), f"{what}: bammetrics left files"
            RESULTS.append((what, 0, 0, "both fail"))
            return done
        for side in dirs:
            assert done[side].returncode == 0, f"{what} {side}: {done[side].stderr[-1500:]}"
        names = {side: sorted(os.listdir(out)) for side, out in dirs.items()}
        assert names["picard"] == names["bammetrics"], f"{what}: {names}"
        problems = []
        records = 0
        for name in names["picard"]:
            picard, ours = dirs["picard"] / name, dirs["bammetrics"] / name
            problems += diff_file(picard, ours)
            if not name.endswith((".bai", ".md5")):
                records += records_of(picard)
        assert not problems, f"{what}:\n  " + "\n  ".join(problems[:20])
        SEEN.append((len(names["picard"]), records))
        RESULTS.append((what, records, records, f"{len(names['picard'])} files equal"))
        return names["picard"]


@pytest.mark.parametrize("name", ["many", "names", "names_sam", "names_sam_gz", "old_header",
                                  "queryname", "unsorted", "empty"])
def test_split_inputs(split, name):
    """An `RG` the header lacks fails Picard's record validation at STRICT, so these run
    at SILENT; bammetrics validates no BAM record."""
    files = compare_split(TOOL, split[name], SILENT, name)
    if name == "many":
        assert len(files) == 300, len(files)


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    source = source_path(name, corpus, exome)
    compare_split(TOOL, source, SILENT, name)


def test_the_odd_names_are_picards(split):
    files = compare_split(TOOL, split["names"], SILENT, "odd names")
    assert "Lib_A.bam" in files and "lib_x.bam" in files and "sp.bam" in files
    # The UTF-8 bytes of the BAM header are two characters each, as htsjdk reads them.
    assert "cafÃ©.bam" in files, files
    assert "unknown.bam" in files and "no_reads.bam" in files
    sam = compare_split(TOOL, split["names_sam"], SILENT, "odd names in SAM text")
    assert "café.sam" in sam, sam
    wide = compare_split(TOOL, split["wide_sam"], SILENT, "libraries above U+00FF in SAM text")
    assert wide == sorted(["中文.sam", "\u3e2d文.sam", "中_文.sam"]), wide


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["CREATE_INDEX=true"], ["CREATE_MD5_FILE=true"],
                                  ["CREATE_INDEX=true", "CREATE_MD5_FILE=true"],
                                  ["MAX_RECORDS_IN_RAM=1"], [f"R={REF}"],
                                  ["VALIDATION_STRINGENCY=LENIENT"]])
def test_options(args):
    compare_split(TOOL, SOURCES["markdup"], args, " ".join(args))


def test_index_beside_queryname_files(split):
    files = compare_split(TOOL, split["queryname"], ["CREATE_INDEX=true"], "queryname index")
    assert not any(name.endswith(".bai") for name in files)


def test_the_default_output_is_the_current_directory():
    compare_split(TOOL, SOURCES["markdup"], [], "no OUTPUT", cwd_output=True)


def test_the_sam_text_fixture():
    compare_split(TOOL, FIXTURES / "markdup.sam", [], "markdup.sam")


def test_a_cram_input_is_wave_w15(corpus, tmp_path):
    """Picard writes a CRAM per library; bammetrics refuses before the first read."""
    cram = tmp_path / "corpus.cram"
    reference = corpus / "ref.fa"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(reference), "-o", str(cram),
                    str(corpus / "corpus.bam")], check=True)
    dirs = {side: tmp_path / side for side in ("picard", "bammetrics")}
    done = {}
    for side, out in dirs.items():
        out.mkdir()
        done[side] = run_tool(side, TOOL, [f"I={cram}", f"O={out}", f"R={reference}", *SILENT])
    assert done["picard"].returncode == 0, done["picard"].stderr[-800:]
    assert sorted(p.name for p in dirs["picard"].iterdir()) == [
        "LibA.cram", "LibB.cram", "LibC.cram", "LibD.cram"]
    assert done["bammetrics"].returncode == 1
    assert "Picard writes a CRAM per library" in done["bammetrics"].stderr
    assert not os.listdir(dirs["bammetrics"])
    verb = subprocess.run([str(BAMMETRICS), "split-sam-by-library", str(cram), "-R",
                           str(reference), "-o", str(tmp_path / "verb")], capture_output=True,
                          text=True)
    assert verb.returncode == 2 and not (tmp_path / "verb").exists()
    RESULTS.append((f"{TOOL} CRAM input", 0, 0, "Picard writes CRAM, bammetrics refuses (D58)"))


@pytest.mark.parametrize("case", ["safe_clash", "unknown_clash", "case_clash"])
def test_one_file_for_two_writers_fails(split, tmp_path, case):
    """Picard opens two writers on one file and exits 0; bammetrics leaves nothing."""
    (tmp_path / "Probe").mkdir()
    if case == "case_clash" and not (tmp_path / "probe").exists():
        pytest.skip("the file system tells AB from ab")
    dirs = {side: tmp_path / side for side in ("picard", "bammetrics")}
    for side, out in dirs.items():
        out.mkdir()
        done = run_tool(side, TOOL, [f"I={split[case]}", f"O={out}"])
        if side == "picard":
            assert done.returncode == 0, done.stderr[-800:]
        else:
            assert done.returncode == 1, done.stderr
    assert not os.listdir(dirs["bammetrics"])
    RESULTS.append((f"{TOOL} {case}", 0, 0, "Picard writes one file over the other (D59)"))


@pytest.mark.parametrize("case", ["no_rg", "no_lb"])
def test_no_library_is_exit_2_in_both(split, case):
    done = compare_split(TOOL, split[case], [], f"error {case}", "fail")
    assert done["picard"].returncode == done["bammetrics"].returncode == 2


@pytest.mark.parametrize("case", ["out_of_order", "rg_int", "assume_sorted", "missing_input",
                                  "missing_directory", "missing_rg_strict_sam"])
def test_errors_as_picard(split, tmp_path, case):
    source, args = SOURCES["markdup"], []
    if case in ("out_of_order", "rg_int"):
        source = split[case]
    elif case == "assume_sorted":
        args = ["ASSUME_SORTED=true"]
    elif case == "missing_input":
        source = tmp_path / "nope.bam"
    elif case == "missing_rg_strict_sam":
        source = split["names_sam"]
    if case == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"O={tmp_path / 'no' / 'dir'}"])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    compare_split(TOOL, source, args, f"error {case}", "fail")


def test_the_step_in_a_chain_equals_picard(corpus, tmp_path):
    """`chain SOURCE + mark-duplicates + split-sam-by-library -o DIR` against Picard's
    MarkDuplicates, then its SplitSamByLibrary on that file."""
    for source in (SOURCES["markdup"], corpus / "corpus.bam"):
        work = tmp_path / source.stem
        picard_dir, ours_dir = work / "picard", work / "chained"
        picard_dir.mkdir(parents=True)
        marked = work / "marked.bam"
        done = run_tool("picard", "MarkDuplicates", [f"I={source}", f"O={marked}",
                                                     f"M={work / 'm.txt'}", *SILENT])
        assert done.returncode == 0, done.stderr[-800:]
        done = run_tool("picard", TOOL, [f"I={marked}", f"O={picard_dir}", *SILENT])
        assert done.returncode == 0, done.stderr[-800:]
        chain = subprocess.run([str(BAMMETRICS), "chain", str(source), "-o", str(work / "c"),
                                "+", "mark-duplicates", "+", "split-sam-by-library", "-o",
                                str(ours_dir)], capture_output=True, text=True)
        assert chain.returncode == 0, chain.stderr
        names = sorted(os.listdir(picard_dir))
        assert names == sorted(os.listdir(ours_dir))
        problems = []
        for name in names:
            problems += diff_file(picard_dir / name, ours_dir / name)
        assert not problems, "\n".join(problems[:20])
        count = records_of(marked)
        RESULTS.append((f"{TOOL} chain mark-duplicates on {source.name}", count, count,
                        f"{len(names)} files equal Picard's"))


def test_values_are_distinct():
    """The cases write 1, 2, 4 and 301 files from inputs of at least three sizes."""
    assert len({files for files, _ in SEEN}) >= 3, SEEN
    assert len({records for _, records in SEEN}) >= 3, SEEN

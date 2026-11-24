"""SamFormatConverter against Picard 3.5.0: the six conversions of SAM text, BAM and CRAM,
each written by both tools from one source file.

Tier B (docs/concordance.md). Inputs: the foundation corpus, t1, the exome slice, every BAM
of tests/data that test_view_sam.py reads, and the BAMs of cram_bams.py. Each input is
taken in three forms, all written by samtools: the BAM, its SAM text (`samtools view -h`)
and a CRAM 3.0 against its reference. A fixture whose `@SQ` lines are not those of
`tests/data/ref_100k.fa` has no reference here, and only its SAM and BAM forms are used.
Both tools convert each form to the other two, reading at `VALIDATION_STRINGENCY=SILENT`
(Picard's default refuses the foundation corpus):

* to SAM text: the two files byte for byte, the `DT` field of `@RG` lines aside;
* to a BAM: the header and every record as samtools prints them, then every record as
  bytes after BGZF inflation (`compare.diff_bam_bytes`);
* to a CRAM: CRAM bytes are never compared. bammetrics' CRAM, decoded by samtools with
  `decode_md=0`, holds the header and every record, tags in order, of Picard's BAM of the
  same source (both tools hold that BAM equal). Picard's own CRAM, decoded the same way,
  holds the same fields, and the tags as a multiset without `MD` and `NM`: htsjdk's CRAM
  writer leaves both out and orders the tags its own way (docs/07 W15-4). That second
  truth is taken only where Picard's CRAM writer keeps the records (`PICARD_CRAM_LOSES`).

`reader_gaps` names the records the two readers give otherwise from one source, for any
writer, and `CRAM_LOSES` what no CRAM holds; the test fails when a set changes either way
(docs/07 W15-23). A CIGAR of 70,000 operators converts as in Picard in all six.

Round trips: each form converted by bammetrics to each other form and back holds the
records of the form it started from (samtools, the tags as a multiset: htsjdk orders
them by tag; a CRAM in the loop adds `MD` and `NM` as htsjdk decodes them and reads an
unsigned array signed, as `same_records` allows; an `H` tag comes back as `B:c`).

Also: Picard's own CRAM 3.1 as a source; the name that picks the format (`.BAM`,
`.CRAM`, `.sam.gz`, `.txt`, no extension: a BAM in both); every option (`CREATE_INDEX`,
`CREATE_MD5_FILE`, `COMPRESSION_LEVEL`, `MAX_RECORDS_IN_RAM`, `VALIDATION_STRINGENCY`
on bad SAM lines); headers without `@HD`, with `VN:1.5`, without `SO`. Error parity: a
CRAM out or in without a reference, `CREATE_INDEX` on a header that is not coordinate
sorted, records out of the header's order, `ASSUME_SORTED`, `OUTPUT` equal to `INPUT`
(Picard truncates it; bammetrics refuses, docs/07 W15-20).
"""
import gzip
import hashlib
import random
import re
import subprocess
from collections import Counter
from pathlib import Path

import pytest

import cram_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import bam_records, diff_bam, diff_bam_bytes, first_difference
from conftest import RESULTS
from run_picard import picard_env
from test_cram_writer import NO_REFERENCE, Decoded, header_lines, same_records, samtools_view
from test_view_sam import BAD, HEADER, VALID

TOOL = "SamFormatConverter"
VERB = "sam-format-converter"
SILENT = "VALIDATION_STRINGENCY=SILENT"
KINDS = ("sam", "bam", "cram")
CONVERSIONS = [(a, b) for a in KINDS for b in KINDS if a != b]
FIXTURE_REFERENCE = FIXTURES / "ref_100k.fa"
# ViewSam's fixtures: every BAM of tests/data but the one with `@SQ LN:0`.
FIXTURE_BAMS = sorted(p.stem for p in FIXTURES.glob("*.bam") if p.name != "index_depth.bam")
EDGES = ("edges", "no_hd", "empty", "slices", "shuffled", "overhang", "overhang_shuffled", "limits")
CORPORA = ("foundation", "t1", "exome")
INPUTS = (*CORPORA, *(f"data/{name}" for name in FIXTURE_BAMS), *(f"edge/{name}" for name in EDGES))

# Records the two readers give otherwise from one source, whoever writes the output
# (KNOWN_DIFFERENCES.md). The test fails when a set changes either way.
#   long_cigar (edges, no_hd): a CIGAR of 70,000 operators (`1M1I` repeated) that ends
#     the reference span of its slice. From a CRAM htsjdk rebuilds NM against that span
#     and leaves out the insertion after its last base (34,999); bammetrics counts it
#     (35,000). The CIGAR itself converts as in Picard in all six (docs/07 W15-23).
LONG_CIGAR = {"edge/edges", "edge/no_hd"}
# The record of each input with an `H` tag.
HEX_TAG = {"edge/edges": {"tags_all"}, "edge/no_hd": {"tags_all"}}


def reader_gaps(name, source_kind, target_kind):
    if name in LONG_CIGAR and source_kind == "cram":
        return {"long_cigar"}
    return set()


# Records a CRAM gives back otherwise under samtools whoever writes it (cram_writer.rs):
# the mapping quality of an unmapped record, `=` and `X` operators, neighbouring `M`
# operators, the CIGAR of an unmapped record, the mate fields of an unpaired record. Such
# an input is also held to samtools' own CRAM of Picard's BAM.
CRAM_LOSES = {
    "data/clean_sam": {"unmapped", "unplaced"},
    "edge/limits": {"eq_x", "m_then_m", "unmapped_cigar", "unpaired_mate_fields"},
}
# Inputs where Picard's CRAM is not a truth (test_cram_writer.py, docs/07 W15-3 to W15-5):
# a read group the header lacks (dropped), mates whose fields disagree with each other
# (htsjdk attaches them and rebuilds the fields: mate_filter, revert), reads past a contig
# end, what no CRAM holds (clean_sam, limits).
PICARD_CRAM_LOSES = {"edge/edges", "edge/no_hd", "edge/overhang", "edge/overhang_shuffled",
                     "edge/limits", "data/clean_sam", "data/mate_filter",
                     "data/mate_filter.queryname", "data/revert"}
SEEN = Counter()


def run(cmd, env=None):
    return subprocess.run([str(c) for c in cmd], capture_output=True, env=env)


def picard(*args):
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return run(["picard", TOOL, *args, "VERBOSITY=ERROR", "QUIET=true"], env)


def ours(*args, face=TOOL):
    return run([BAMMETRICS, face, *args])


def side_run(side, *args):
    return picard(*args) if side == "picard" else ours(*args)


@pytest.fixture(scope="module")
def work(tmp_path_factory):
    return tmp_path_factory.mktemp("sam_format_converter")


@pytest.fixture(scope="module")
def edge(work):
    directory = work / "cram_bams"
    return cram_bams.write_all(directory), cram_bams.reference(directory)


def fasta_contigs(fai):
    return [tuple(line.split("\t")[:2]) for line in Path(f"{fai}.fai").read_text().splitlines()]


def header_contigs(bam):
    out = run(["samtools", "view", "-H", bam]).stdout.decode("latin-1")
    return [(f["SN"], f["LN"]) for f in (
        dict(field.split(":", 1) for field in line.split("\t")[1:] if ":" in field)
        for line in out.split("\n") if line.startswith("@SQ"))]


def locate(name, corpus, exome, edge):
    """The BAM of input `name` and its reference, None where it has none."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name == "exome":
        return exome / "exome.bam", exome / "ref.fa"
    if name == "t1":
        return CORPUS / "t1.bam", CORPUS / "reference.fa"
    kind, stem = name.split("/")
    if kind == "edge":
        return edge[0][stem], edge[1]
    bam = FIXTURES / f"{stem}.bam"
    contigs = header_contigs(bam)
    fits = set(contigs) <= set(fasta_contigs(FIXTURE_REFERENCE))
    return bam, FIXTURE_REFERENCE if fits else None


class Input:
    """One input in its three forms, and every conversion of them, each run once."""

    def __init__(self, name, bam, reference, directory):
        self.name, self.bam, self.reference = name, bam, reference
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=True)
        self.forms = {"bam": bam}
        sam = directory / "source.sam"
        sam.write_bytes(run(["samtools", "view", "-h", "--no-PG", bam]).stdout)
        self.forms["sam"] = sam
        if reference is not None:
            cram = directory / "source.cram"
            made = run(["samtools", "view", "--no-PG", "-T", reference, "--output-fmt",
                        "cram,version=3.0", "-o", cram, bam])
            assert made.returncode == 0, made.stderr.decode()[-400:]
            self.forms["cram"] = cram
        self.done = {}

    def kinds(self):
        return [k for k in KINDS if k in self.forms]

    def args(self, source, out, *extra):
        reference = [f"R={self.reference}"] if self.reference is not None else []
        return [f"I={source}", f"O={out}", SILENT, *reference, *extra]

    def run(self, side, source_kind, target_kind):
        """The process and the output of one conversion, run once."""
        key = (side, source_kind, target_kind)
        if key not in self.done:
            out = self.directory / f"{side}.{source_kind}_to.{target_kind}"
            self.done[key] = side_run(side, *self.args(self.forms[source_kind], out)), out
        return self.done[key]

    def convert(self, side, source_kind, target_kind):
        done, out = self.run(side, source_kind, target_kind)
        assert done.returncode == 0, (
            f"{self.name} {source_kind} -> {target_kind} ({side}): {done.stderr.decode('latin-1')[:800]}")
        return out


INPUT_CACHE = {}


def the_input(name, corpus, exome, edge, work):
    if name not in INPUT_CACHE:
        bam, reference = locate(name, corpus, exome, edge)
        INPUT_CACHE[name] = Input(name, bam, reference, work / name.replace("/", "_"))
    return INPUT_CACHE[name]


DATE = re.compile(rb"\tDT:[^\t\n]*")


def without_dates(text):
    """SAM text with the DT field cut from @RG lines (htsjdk writes it in the JVM's zone)."""
    return b"".join(DATE.sub(b"", line) if line.startswith(b"@RG") else line
                    for line in text.splitlines(keepends=True))


def split_gaps(what, picard_items, ours_items, name_of, allow):
    """The two lists without the items `allow` names, after checking that each named item
    is in both, once, and differs."""
    if not allow:
        return picard_items, ours_items
    named = {}
    for side, items in (("picard", picard_items), ("bammetrics", ours_items)):
        named[side] = [item for item in items if name_of(item) in allow]
        assert sorted(map(name_of, named[side])) == sorted(allow), f"{what}: {side} lacks {allow}"
    for x, y in zip(named["picard"], named["bammetrics"]):
        assert x != y, f"{what}: {name_of(x)} no longer differs: drop it from the gaps"
    return ([item for item in picard_items if name_of(item) not in allow],
            [item for item in ours_items if name_of(item) not in allow])


def same_sam_text(what, picard_out, ours_out, allow=()):
    a = without_dates(picard_out.read_bytes()).split(b"\n")
    b = without_dates(ours_out.read_bytes()).split(b"\n")
    a, b = split_gaps(what, a, b, lambda line: line.split(b"\t", 1)[0].decode("latin-1"), allow)
    a, b = b"\n".join(a), b"\n".join(b)
    assert a == b, f"{what}: {first_difference(a, b)}"
    return a.count(b"\n")


def record_name(raw):
    return raw[32:32 + raw[8] - 1].decode("latin-1")


def same_bam(what, picard_out, ours_out, reference=None, allow=()):
    if allow:
        theirs, mine = samtools_view(picard_out, reference), samtools_view(ours_out, reference)
        assert header_without_dates(mine) == header_without_dates(theirs), f"{what}: header"
        same_records(what, theirs, mine, allow=allow)
        p, o = split_gaps(what, bam_records(picard_out), bam_records(ours_out), record_name, allow)
        problems = [f"record {at}: bytes differ" for at, (x, y) in enumerate(zip(p, o)) if x != y][:5]
    else:
        problems = diff_bam(picard_out, ours_out)
        if not problems:
            problems = diff_bam_bytes(picard_out, ours_out)
    assert not problems, f"{what}:\n  " + "\n  ".join(problems[:20])
    return int(run(["samtools", "view", "-c", ours_out]).stdout)


def without_md_nm(decoded):
    for record in decoded.records:
        record[11:] = sorted(tag for tag in record[11:] if tag[:2] not in ("MD", "NM"))
    return decoded


def header_without_dates(decoded):
    return [DATE.sub(b"", line.encode("latin-1")).decode("latin-1") for line in header_lines(decoded)]


def same_cram(what, item, source_kind, ours_out):
    """bammetrics' CRAM against Picard's BAM of the same source, then against Picard's CRAM."""
    picard_bam = item.convert("picard", source_kind, "bam")
    truth = samtools_view(picard_bam, item.reference)
    mine = samtools_view(ours_out, item.reference)
    assert ours_out.read_bytes()[:6] == b"CRAM\x03\x00", f"{what}: not a CRAM 3.0"
    assert header_without_dates(mine) == header_without_dates(truth), f"{what}: header"
    loses = CRAM_LOSES.get(item.name, set())
    same_records(f"{what}: bammetrics' CRAM against Picard's BAM", truth, mine, allow=loses)
    if loses:
        theirs = item.directory / f"samtools.{source_kind}.cram"
        made = run(["samtools", "view", "--no-PG", "-C", "-T", item.reference, "-o", theirs, picard_bam])
        assert made.returncode == 0, made.stderr.decode()[-300:]
        same_records(f"{what}: bammetrics' CRAM against samtools' CRAM of Picard's BAM",
                     samtools_view(theirs, item.reference), mine, multiset=True)
    if item.name not in PICARD_CRAM_LOSES:
        theirs = samtools_view(item.convert("picard", source_kind, "cram"), item.reference)
        assert header_without_dates(theirs) == header_without_dates(mine), f"{what}: header of Picard's CRAM"
        same_records(f"{what}: Picard's CRAM against bammetrics'", without_md_nm(theirs),
                     without_md_nm(samtools_view(ours_out, item.reference)))
    return len(truth.records)


@pytest.mark.parametrize("source_kind,target_kind", CONVERSIONS, ids=lambda c: c)
@pytest.mark.parametrize("name", INPUTS)
def test_conversion(name, source_kind, target_kind, corpus, exome, edge, work):
    item = the_input(name, corpus, exome, edge, work)
    if source_kind not in item.forms or (target_kind == "cram" and item.reference is None):
        pytest.skip(f"{name}: no reference, so no CRAM")
    what = f"{TOOL} {name} {source_kind} -> {target_kind}"
    ours_out = item.convert("bammetrics", source_kind, target_kind)
    allow = reader_gaps(name, source_kind, target_kind)
    if target_kind == "sam":
        count = same_sam_text(what, item.convert("picard", source_kind, "sam"), ours_out, allow)
    elif target_kind == "bam":
        count = same_bam(what, item.convert("picard", source_kind, "bam"), ours_out, item.reference, allow)
    else:
        count = same_cram(what, item, source_kind, ours_out)
    SEEN[source_kind, target_kind] += 1
    note = "equal to Picard's" + (f"; reader gaps {sorted(allow)}" if allow else "")
    RESULTS.append((what, count, count, note))


@pytest.mark.parametrize("source_kind,target_kind", CONVERSIONS, ids=lambda c: c)
@pytest.mark.parametrize("name", INPUTS)
def test_round_trip(name, source_kind, target_kind, corpus, exome, edge, work):
    """bammetrics from one form to another and back holds the records of the first."""
    item = the_input(name, corpus, exome, edge, work)
    if source_kind not in item.forms or target_kind not in item.forms:
        pytest.skip(f"{name}: no reference, so no CRAM")
    there = item.convert("bammetrics", source_kind, target_kind)
    back = item.directory / f"back.{source_kind}_{target_kind}.{source_kind}"
    done = ours(*item.args(there, back))
    assert done.returncode == 0, done.stderr.decode("latin-1")[-600:]
    first, last = samtools_view(item.forms[source_kind], item.reference), samtools_view(back, item.reference)
    assert [line for line in header_lines(last) if not line.startswith("@HD")] == [
        line for line in header_lines(first) if not line.startswith("@HD")], f"{name}: header"
    # A CRAM in the loop: MD and NM rebuilt, unsigned arrays read signed, and what no CRAM
    # holds lost on the way, unless the first form was a CRAM already.
    rules = ("rebuilt", "signed") if "cram" in (source_kind, target_kind) else ()
    allow = set(CRAM_LOSES.get(name, set())) if target_kind == "cram" else set()
    # From a CRAM samtools reads the mate start of an unpaired record that bammetrics,
    # as htsjdk, reads as 0.
    if source_kind == "cram" and "unpaired_mate_fields" in CRAM_LOSES.get(name, set()):
        allow.add("unpaired_mate_fields")
    # An `H` tag comes back as `B:c` over the same bytes: htsjdk holds it so, and Picard's
    # conversion writes it so (test_conversion holds the first step to Picard's file).
    allow |= HEX_TAG.get(name, set())
    # Tags come back in htsjdk's order (by tag value), as Picard writes them.
    used = same_records(f"{name} {source_kind} -> {target_kind} -> {source_kind}", first, last,
                        rules=rules, multiset=True, allow=allow)
    RESULTS.append((f"{TOOL} round trip {name} {source_kind}->{target_kind}", len(first.records),
                    len(first.records), f"MD/NM rebuilt {used['rebuilt']}, signed {used['signed']}"))


@pytest.mark.parametrize("target_kind", ("sam", "bam"))
def test_picards_cram_3_1_as_the_source(corpus, tmp_path, target_kind):
    bam, reference = corpus / "corpus.bam", corpus / "ref.fa"
    source = tmp_path / "picard.cram"
    assert picard(f"I={bam}", f"O={source}", f"R={reference}", SILENT).returncode == 0
    assert source.read_bytes()[:6] == b"CRAM\x03\x01"
    outs = {}
    for side in ("picard", "bammetrics"):
        outs[side] = tmp_path / f"{side}.{target_kind}"
        done = side_run(side, f"I={source}", f"O={outs[side]}", f"R={reference}", SILENT)
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-600:]}"
    what = f"{TOOL} Picard's CRAM 3.1 -> {target_kind}"
    compare = same_sam_text if target_kind == "sam" else same_bam
    count = compare(what, outs["picard"], outs["bammetrics"])
    RESULTS.append((what, count, count, "equal to Picard's"))


@pytest.mark.parametrize("name", ["out.txt", "out.BAM", "out.CRAM", "out.Sam", "out.sam.gz", "out"])
def test_the_name_decides_the_format(tmp_path, name):
    """`.sam` and `.cram`, lower case, name SAM text and a CRAM; any other name is a BAM in
    both, which htsjdk logs as `Unknown file extension, assuming BAM format`."""
    outs = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        outs[side] = directory / name
        done = side_run(side, f"I={FIXTURES / 'markdup.bam'}", f"O={outs[side]}", f"R={FIXTURE_REFERENCE}")
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-600:]}"
        assert outs[side].read_bytes()[:4] == b"\x1f\x8b\x08\x04", f"{side}: {name} is not a BAM"
    count = same_bam(f"{TOOL} O={name}", outs["picard"], outs["bammetrics"])
    RESULTS.append((f"{TOOL} O={name}", count, count, "a BAM in both"))


def both_with(tmp_path, source, out_name, *extra):
    """Both tools on `source`, each into a directory of its own under `tmp_path`:
    {side: (process, directory)}."""
    sides = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir(parents=True, exist_ok=True)
        done = side_run(side, f"I={source}", f"O={directory / out_name}", SILENT, *extra)
        sides[side] = (done, directory)
    return sides


def names(directory):
    return sorted(p.name for p in directory.iterdir())


@pytest.mark.parametrize("name", ["foundation", "markdup"])
def test_create_index_and_md5_beside_a_bam(corpus, tmp_path, name):
    """Both write out.bai and out.bam.md5. bammetrics' .bai is Picard BuildBamIndex of its
    BAM byte for byte; each .md5 is the hex MD5 of its BAM."""
    source = corpus / "corpus.bam" if name == "foundation" else FIXTURES / "markdup.bam"
    sides = both_with(tmp_path, source, "out.bam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true")
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-600:]}"
        assert names(directory) == ["out.bai", "out.bam", "out.bam.md5"], side
        bam = (directory / "out.bam").read_bytes()
        assert (directory / "out.bam.md5").read_bytes() == hashlib.md5(bam).hexdigest().encode(), side
    ours_dir = sides["bammetrics"][1]
    rebuilt = tmp_path / "rebuilt.bai"
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    index = run(["picard", "BuildBamIndex", f"I={ours_dir / 'out.bam'}", f"O={rebuilt}", SILENT,
                 "QUIET=true", "VERBOSITY=ERROR"], env)
    assert index.returncode == 0, index.stderr[-600:]
    assert (ours_dir / "out.bai").read_bytes() == rebuilt.read_bytes()
    count = same_bam(f"{TOOL} CREATE_INDEX {name}", sides["picard"][1] / "out.bam", ours_dir / "out.bam")
    SEEN["CREATE_INDEX"] += 1
    SEEN["CREATE_MD5_FILE"] += 1
    RESULTS.append((f"{TOOL} CREATE_INDEX {name}", count, count, ".bai is Picard's index of the BAM"))


def test_sam_text_gets_an_md5_and_no_index(tmp_path):
    sides = both_with(tmp_path, FIXTURES / "markdup.bam", "out.sam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true")
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-600:]}"
        assert names(directory) == ["out.sam", "out.sam.md5"], side
        text = (directory / "out.sam").read_bytes()
        assert (directory / "out.sam.md5").read_bytes() == hashlib.md5(text).hexdigest().encode(), side
    count = same_sam_text(f"{TOOL} SAM with CREATE_INDEX", sides["picard"][1] / "out.sam",
                          sides["bammetrics"][1] / "out.sam")
    RESULTS.append((f"{TOOL} SAM CREATE_INDEX", count, count, "no index in either, .md5 in both"))


def test_a_cram_gets_an_md5_and_its_index_is_refused(tmp_path):
    """A departure (docs/07 W15-2): Picard writes `<CRAM>.bai`; bammetrics refuses
    CREATE_INDEX with a CRAM and writes nothing."""
    base = (f"R={FIXTURE_REFERENCE}",)
    sides = both_with(tmp_path, FIXTURES / "markdup.bam", "out.cram", *base, "CREATE_MD5_FILE=true")
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-600:]}"
        assert names(directory) == ["out.cram", "out.cram.md5"], side
        cram = (directory / "out.cram").read_bytes()
        assert (directory / "out.cram.md5").read_bytes() == hashlib.md5(cram).hexdigest().encode(), side
    indexed = both_with(tmp_path / "indexed", FIXTURES / "markdup.bam", "out.cram", *base, "CREATE_INDEX=true")
    done, directory = indexed["picard"]
    assert done.returncode == 0 and names(directory) == ["out.cram", "out.cram.bai"], names(directory)
    done, directory = indexed["bammetrics"]
    assert done.returncode == 1 and b"no index is written beside a CRAM" in done.stderr, done.stderr
    assert names(directory) == []
    RESULTS.append((f"{TOOL} CRAM CREATE_INDEX", 1, 1, "MD5 in both; Picard's .cram.bai refused"))


@pytest.mark.parametrize("target", ["out.bam", "out.sam"])
def test_create_index_on_a_queryname_header_fails_in_both(tmp_path, target):
    """Picard throws after it has made OUTPUT and leaves it; bammetrics fails before."""
    sides = both_with(tmp_path, FIXTURES / "duplicate_marking.bam", target, "CREATE_INDEX=true")
    sentence = b"Can't CREATE_INDEX unless sort order is coordinate"
    for side, (done, _) in sides.items():
        assert done.returncode == 1 and sentence in done.stderr, f"{side}: {done.stderr[-400:]}"
    assert names(sides["bammetrics"][1]) == []
    assert names(sides["picard"][1]) == [target]
    RESULTS.append((f"{TOOL} CREATE_INDEX queryname {target}", 0, 0, "both exit 1"))


def test_options_that_change_nothing_or_the_bytes(corpus, tmp_path):
    """MAX_RECORDS_IN_RAM changes nothing; COMPRESSION_LEVEL changes the bytes of a BAM and
    not its records, in both tools."""
    bam = corpus / "corpus.bam"
    plain = tmp_path / "plain.bam"
    assert ours(f"I={bam}", f"O={plain}", SILENT).returncode == 0
    few = tmp_path / "few.bam"
    assert ours(f"I={bam}", f"O={few}", SILENT, "MAX_RECORDS_IN_RAM=1").returncode == 0
    assert few.read_bytes() == plain.read_bytes()
    sizes = {}
    for level in (0, 9):
        sides = both_with(tmp_path / f"level{level}", bam, "out.bam", f"COMPRESSION_LEVEL={level}")
        for side, (done, _) in sides.items():
            assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
        same_bam(f"{TOOL} COMPRESSION_LEVEL={level}", sides["picard"][1] / "out.bam",
                 sides["bammetrics"][1] / "out.bam")
        sizes[level] = (sides["bammetrics"][1] / "out.bam").stat().st_size
    assert sizes[9] < sizes[0], sizes
    SEEN["MAX_RECORDS_IN_RAM"] += 1
    SEEN["COMPRESSION_LEVEL"] += 1
    RESULTS.append((f"{TOOL} options", 3, 3, "MAX_RECORDS_IN_RAM same bytes, levels 0 and 9"))


@pytest.mark.parametrize("target", ["out.sam", "out.bam"])
@pytest.mark.parametrize("header", ["no_hd", "vn_1_5", "no_so", "unknown_so"])
def test_headers(tmp_path, header, target):
    """The header as htsjdk's writer writes it: `@HD VN:1.6`, `SO:unsorted` when absent."""
    head = {"no_hd": "", "vn_1_5": "@HD\tVN:1.5\tSO:unsorted\n", "no_so": "@HD\tVN:1.5\n",
            "unknown_so": "@HD\tVN:1.6\tSO:sorted\n"}[header]
    source = tmp_path / "in.sam"
    source.write_text(head + HEADER.split("\n", 1)[1] + "".join(line + "\n" for line in VALID))
    sides = both_with(tmp_path, source, target)
    for side, (done, _) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
    compare = same_sam_text if target.endswith(".sam") else same_bam
    count = compare(f"{TOOL} header {header} -> {target}", sides["picard"][1] / target,
                    sides["bammetrics"][1] / target)
    RESULTS.append((f"{TOOL} header {header} {target}", count, count, "equal to Picard's"))


BAD_LINES = BAD[:12] + BAD[-6:]


@pytest.mark.parametrize("stringency", ("STRICT", "LENIENT", "SILENT"))
@pytest.mark.parametrize("name,line", BAD_LINES, ids=[name for name, _ in BAD_LINES])
def test_bad_lines(tmp_path, name, line, stringency):
    """A bad SAM line at each stringency: both fail, or both write the same SAM text."""
    source = tmp_path / "bad.sam"
    source.write_text(HEADER + VALID[0] + "\n" + line + "\n")
    outs, codes = {}, {}
    for side in ("picard", "bammetrics"):
        outs[side] = tmp_path / f"{side}.sam"
        done = side_run(side, f"I={source}", f"O={outs[side]}", f"VALIDATION_STRINGENCY={stringency}")
        codes[side] = done.returncode
    assert (codes["picard"] == 0) == (codes["bammetrics"] == 0), f"{name} {stringency}: exit codes {codes}"
    if codes["picard"] == 0:
        same_sam_text(f"{TOOL} {name} {stringency}", outs["picard"], outs["bammetrics"])
    else:
        assert not outs["bammetrics"].exists(), f"{name} {stringency}: bammetrics left a file"
    SEEN["VALIDATION_STRINGENCY"] += 1
    RESULTS.append((f"{TOOL} bad line {name} {stringency}", 1, 1,
                    "both fail" if codes["picard"] else "same SAM text"))


@pytest.mark.parametrize("target", ["out.sam", "out.bam", "out.cram"])
@pytest.mark.parametrize("order", ["coordinate", "queryname"])
def test_records_out_of_the_headers_order_fail_in_both(tmp_path, order, target):
    source = tmp_path / "in.sam"
    records = [VALID[6], VALID[0]] if order == "coordinate" else [VALID[2], VALID[0]]
    source.write_text(f"@HD\tVN:1.6\tSO:{order}\n" + HEADER.split("\n", 1)[1]
                      + "".join(line + "\n" for line in records))
    reference = tmp_path / "ref.fa"
    write_handwritten_reference(reference)
    sides = both_with(tmp_path, source, target, f"R={reference}")
    for side, (done, _) in sides.items():
        assert done.returncode == 1, f"{side} exit {done.returncode}: {done.stderr[-400:]}"
    assert names(sides["bammetrics"][1]) == []
    RESULTS.append((f"{TOOL} {order} out of order {target}", 0, 0, "both exit 1"))


def write_handwritten_reference(path):
    """chr1 of ref_100k.fa and a chr2 of 5,000 bases: the dictionary of test_view_sam's HEADER."""
    bases = "".join(FIXTURE_REFERENCE.read_text().split("\n", 1)[1].split())
    with open(path, "w") as out:
        out.write(">chr1\n" + bases + "\n>chr2\n" + ("ACGTTGCA" * 625) + "\n")
    subprocess.run(["samtools", "faidx", str(path)], check=True)


def test_a_cram_without_a_reference(corpus, tmp_path):
    """Out: both exit 1 with htsjdk's sentence; Picard leaves an empty file, bammetrics none.
    In: both fail and bammetrics leaves no file."""
    bam = corpus / "corpus.bam"
    sides = both_with(tmp_path / "out", bam, "out.cram")
    for side, (done, _) in sides.items():
        assert done.returncode == 1 and NO_REFERENCE.encode() in done.stderr, f"{side}: {done.stderr[-400:]}"
    assert names(sides["bammetrics"][1]) == []
    assert (sides["picard"][1] / "out.cram").stat().st_size == 0
    cram = tmp_path / "in.cram"
    assert run(["samtools", "view", "--no-PG", "-C", "-T", corpus / "ref.fa", "-o", cram, bam]).returncode == 0
    sides = both_with(tmp_path / "in", cram, "out.sam")
    for side, (done, _) in sides.items():
        assert done.returncode != 0, f"{side} read a CRAM without its reference"
    assert names(sides["bammetrics"][1]) == []
    RESULTS.append((f"{TOOL} CRAM without R", 0, 0, "both fail, out and in"))


def test_assume_sorted_is_no_option(tmp_path):
    sides = both_with(tmp_path, FIXTURES / "markdup.bam", "out.sam", "ASSUME_SORTED=true")
    for side, (done, _) in sides.items():
        assert done.returncode != 0, f"{side} took ASSUME_SORTED"
    assert names(sides["bammetrics"][1]) == []
    RESULTS.append((f"{TOOL} ASSUME_SORTED", 0, 0, "both refuse"))


def test_output_equal_to_input(tmp_path):
    """Picard opens the output over its input: a BAM read in place fails and is left empty.
    bammetrics refuses before it reads (docs/07 W15-20) and leaves the input whole."""
    original = (FIXTURES / "markdup.bam").read_bytes()
    for side in ("picard", "bammetrics"):
        copy = tmp_path / f"{side}.bam"
        copy.write_bytes(original)
        done = side_run(side, f"I={copy}", f"O={copy}", SILENT)
        assert done.returncode == 1, f"{side}: exit {done.returncode}"
        if side == "picard":
            assert copy.stat().st_size == 0, "Picard no longer empties its input"
        else:
            assert copy.read_bytes() == original
            assert b"writing it would destroy it" in done.stderr
    RESULTS.append((f"{TOOL} OUTPUT=INPUT", 0, 0, "both exit 1; Picard empties the input"))


# The 28-byte BGZF end-of-file block.
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
TAILS = {
    "zero_1": b"\0", "zero_17": b"\0" * 17, "zero_18": b"\0" * 18, "zero_1000": b"\0" * 1000,
    "ff_10": b"\xff" * 10, "header_4": EOF_BLOCK[:4], "header_18": EOF_BLOCK[:18],
    "random_100": random.Random(15).randbytes(100), "gzip_member": gzip.compress(b"text\n"),
    "eof_twice": EOF_BLOCK,
}


@pytest.mark.parametrize("target", ["out.sam", "out.cram"])
@pytest.mark.parametrize("tail", sorted(TAILS))
def test_bytes_after_the_end_of_file_block(tmp_path, tail, target):
    """htsjdk reads an empty BGZF block as the end of the stream: the bytes after the
    end-of-file block are never read. Both convert every record (docs/07 W15-24)."""
    source = tmp_path / "in.bam"
    source.write_bytes((FIXTURES / "markdup.bam").read_bytes() + TAILS[tail])
    sides = both_with(tmp_path, source, target, f"R={FIXTURE_REFERENCE}")
    for side, (done, _) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
    picard_out, ours_out = (sides[side][1] / target for side in ("picard", "bammetrics"))
    if target == "out.sam":
        count = same_sam_text(f"{TOOL} tail {tail}", picard_out, ours_out)
    else:
        theirs = samtools_view(picard_out, FIXTURE_REFERENCE)
        count = len(theirs.records)
        same_records(f"{TOOL} tail {tail} CRAM", without_md_nm(theirs),
                     without_md_nm(samtools_view(ours_out, FIXTURE_REFERENCE)))
    assert count > 0
    RESULTS.append((f"{TOOL} bytes after the EOF block: {tail} -> {target}", count, count,
                    "equal to Picard's"))


@pytest.mark.parametrize("size", [400, 1000])
def test_a_bam_cut_inside_its_first_block(tmp_path, size):
    """KNOWN_DIFFERENCES.md: markdup.bam cut to `size` bytes. Picard exits 0 and writes the
    header alone; bammetrics exits 1 (`truncated BAM`) and writes nothing."""
    source = tmp_path / "cut.bam"
    source.write_bytes((FIXTURES / "markdup.bam").read_bytes()[:size])
    sides = both_with(tmp_path, source, "out.sam")
    done, directory = sides["picard"]
    assert done.returncode == 0, "Picard now fails on a cut BAM"
    assert all(line.startswith(b"@") for line in (directory / "out.sam").read_bytes().splitlines())
    done, directory = sides["bammetrics"]
    assert done.returncode == 1 and b"truncated BAM" in done.stderr, done.stderr[-300:]
    assert names(directory) == []
    RESULTS.append((f"{TOOL} BAM cut at {size} bytes", 0, 0, "Picard 0 and a header, bammetrics 1"))


HEAD_UTF8 = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
             "@RG\tID:g1\tSM:s\tDS:café über\n@CO\tcafé\n@CO\tcheck ✓\n")
BODY_UTF8 = "r1\t0\tchr1\t10\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:g1\n"


@pytest.mark.parametrize("source_kind,target_kind", CONVERSIONS, ids=lambda c: c)
def test_header_text_outside_ascii(tmp_path, source_kind, target_kind):
    """A header in UTF-8 (samtools' BAM and CRAM of the text hold its bytes). htsjdk holds
    a header one character per byte of a BAM, per character of SAM text or of a CRAM's
    UTF-8, writes SAM text and a BAM as the low byte of each character and a CRAM as UTF-8.
    Every header is Picard's, but one line: SAM text to a CRAM, where bammetrics holds `✓`
    (U+2713) as its low byte 0x13 and Picard stores its UTF-8 (docs/07 W15-13)."""
    sam = tmp_path / "in.sam"
    sam.write_bytes((HEAD_UTF8 + BODY_UTF8).encode())
    sources = {"sam": sam, "bam": tmp_path / "in.bam", "cram": tmp_path / "in.cram"}
    for kind in ("bam", "cram"):
        made = run(["samtools", "view", "--no-PG", "-T", FIXTURE_REFERENCE, "--output-fmt", kind,
                    "-o", sources[kind], sam])
        assert made.returncode == 0, made.stderr.decode()[-300:]
    target = f"out.{target_kind}"
    sides = both_with(tmp_path, sources[source_kind], target, f"R={FIXTURE_REFERENCE}")
    headers = {}
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
        headers[side] = run(["samtools", "view", "-H", "--no-PG", directory / target]).stdout
    picard_lines, ours_lines = (headers[side].split(b"\n") for side in ("picard", "bammetrics"))
    if (source_kind, target_kind) == ("sam", "cram"):
        differ = [(a, b) for a, b in zip(picard_lines, ours_lines) if a != b]
        assert differ == [(b"@CO\tcheck \xe2\x9c\x93", b"@CO\tcheck \x13")], differ
    else:
        assert ours_lines == picard_lines
    assert b"@RG\tID:g1\tSM:s\tDS:caf" in headers["bammetrics"]
    RESULTS.append((f"{TOOL} header outside ASCII {source_kind} -> {target_kind}", 1, 1,
                    "Picard's header" + (", but U+2713" if source_kind == "sam" and target_kind == "cram" else "")))


def test_a_repeated_read_group_id(tmp_path):
    """KNOWN_DIFFERENCES.md: a header with two `@RG` lines of one ID. Picard keeps both;
    bammetrics' header (noodles' `sam::Header`, read groups keyed by ID) holds one, so it
    exits 1 naming the line and writes nothing (docs/07 W15-25)."""
    source = tmp_path / "in.sam"
    source.write_text("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:rgA\tSM:s\n"
                      "@RG\tID:rgA\tSM:other\nr1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:rgA\n")
    for target in ("out.sam", "out.bam"):
        sides = both_with(tmp_path / target, source, target)
        done, directory = sides["picard"]
        assert done.returncode == 0, "Picard now refuses a repeated read group ID"
        text = run(["samtools", "view", "-H", "--no-PG", directory / target]).stdout
        assert text.count(b"@RG\tID:rgA") == 2
        done, directory = sides["bammetrics"]
        assert done.returncode == 1 and b"duplicate read group ID: rgA" in done.stderr, done.stderr
        assert names(directory) == []
    RESULTS.append((f"{TOOL} repeated @RG ID", 0, 0, "Picard keeps both lines, bammetrics exits 1"))


def test_a_record_without_bases_in_a_cram(tmp_path):
    """KNOWN_DIFFERENCES.md, docs/07 W15-26: a mapped record with SEQ `*`, `100M`. Picard
    reads its own CRAM of it as `0M`; Picard and samtools read bammetrics' CRAM as `100M`,
    and samtools reads Picard's as `100M` too."""
    source = tmp_path / "in.sam"
    source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
                      "r1\t0\tchr1\t49\t30\t100M\t*\t0\t0\t*\t*\tXA:Z:x\n"
                      "r2\t0\tchr1\t60\t30\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n")
    sides = both_with(tmp_path, source, "out.cram", f"R={FIXTURE_REFERENCE}")
    cigars = {}
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
        cram = directory / "out.cram"
        cigars[side] = ([r[5] for r in picard_view(cram).records],
                        [r[5] for r in samtools_view(cram, FIXTURE_REFERENCE).records])
    assert cigars == {"picard": (["0M", "8M"], ["100M", "8M"]),
                      "bammetrics": (["100M", "8M"], ["100M", "8M"])}, cigars
    RESULTS.append((f"{TOOL} SEQ * record to CRAM", 2, 2, "Picard reads its own as 0M"))


def picard_view(path):
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return Decoded("Picard", run(["picard", "ViewSam", f"I={path}", f"R={FIXTURE_REFERENCE}", SILENT,
                                  "QUIET=true", "VERBOSITY=ERROR"], env), header=True)


def test_a_cram_read_against_a_shorter_fasta(tmp_path):
    """KNOWN_DIFFERENCES.md: markdup.bam (reads at 1001 to 9550) as a CRAM, read with a FASTA
    whose chr1 is the first 19,620 bases. Picard converts it; bammetrics exits 1 naming the
    two lengths (docs/07 D52)."""
    bases = "".join(FIXTURE_REFERENCE.read_text().split("\n", 1)[1].split())[:19_620]
    short = tmp_path / "short.fa"
    short.write_text(">chr1\n" + "\n".join(bases[i:i + 60] for i in range(0, len(bases), 60)) + "\n")
    subprocess.run(["samtools", "faidx", str(short)], check=True)
    cram = tmp_path / "in.cram"
    made = run(["samtools", "view", "--no-PG", "-C", "-T", FIXTURE_REFERENCE, "-o", cram,
                FIXTURES / "markdup.bam"])
    assert made.returncode == 0, made.stderr.decode()[-300:]
    sides = both_with(tmp_path, cram, "out.sam", f"R={short}")
    assert sides["picard"][0].returncode == 0, "Picard now refuses the shorter FASTA"
    done, directory = sides["bammetrics"]
    assert done.returncode == 1 and b"(19620 bp) in the reference" in done.stderr, done.stderr
    assert names(directory) == []
    RESULTS.append((f"{TOOL} CRAM with a shorter FASTA", 0, 0, "Picard 0, bammetrics 1"))


@pytest.mark.parametrize("kind", ["embed_ref", "no_ref"])
def test_a_cram_with_its_own_reference_or_none(tmp_path, kind):
    """KNOWN_DIFFERENCES.md: samtools' CRAM of markdup.bam with `embed_ref=1`, read without
    `R=`, and with `no_ref=1`, read with it. The records are Picard's; `MD` and `NM` are
    not. htsjdk rebuilds both from an embedded reference and none without a reference;
    bammetrics rebuilds both from `R=` alone."""
    cram = tmp_path / "in.cram"
    made = run(["samtools", "view", "--no-PG", "-T", FIXTURE_REFERENCE, "--output-fmt",
                f"cram,{kind}=1", "-o", cram, FIXTURES / "markdup.bam"])
    assert made.returncode == 0, made.stderr.decode()[-300:]
    extra = [f"R={FIXTURE_REFERENCE}"] if kind == "no_ref" else []
    sides = both_with(tmp_path, cram, "out.sam", *extra)
    lines = {}
    for side, (done, directory) in sides.items():
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
        lines[side] = [line.split("\t") for line in (directory / "out.sam").read_text().splitlines()
                       if not line.startswith("@")]
    with_md = {side: sum(any(f.startswith("MD:Z:") for f in fields) for fields in records)
               for side, records in lines.items()}
    expected = {"picard": 33, "bammetrics": 0} if kind == "embed_ref" else {"picard": 0, "bammetrics": 33}
    assert with_md == expected, with_md
    strip = lambda records: [[f for f in fields if f[:3] not in ("MD:", "NM:")] for fields in records]
    assert strip(lines["picard"]) == strip(lines["bammetrics"])
    RESULTS.append((f"{TOOL} CRAM {kind}", len(lines["picard"]), len(lines["picard"]),
                    "records equal; MD and NM differ"))


def test_verb_forms(tmp_path):
    """The verb writes the file the Picard face writes; the chain refuses it as a step."""
    face, verb = tmp_path / "face.sam", tmp_path / "verb.sam"
    assert ours(f"I={FIXTURES / 'markdup.bam'}", f"O={face}").returncode == 0
    done = ours(str(FIXTURES / "markdup.bam"), "-o", str(verb), face=VERB)
    assert done.returncode == 0, done.stderr
    assert face.read_bytes() == verb.read_bytes()
    cram = tmp_path / "verb.cram"
    done = ours(str(FIXTURES / "markdup.sam"), "-o", str(cram), "-R", str(FIXTURE_REFERENCE),
                "--create-md5-file", face=VERB)
    assert done.returncode == 0, done.stderr
    assert sorted(p.name for p in tmp_path.iterdir() if p.name.startswith("verb.cram")) == [
        "verb.cram", "verb.cram.md5"]
    chained = run([BAMMETRICS, "chain", FIXTURES / "markdup.bam", "+", VERB])
    assert chained.returncode == 2 and b"what the write step does" in chained.stderr, chained.stderr
    RESULTS.append((f"{VERB} verb forms", 1, 1, "verb = Picard face; refused as a chain step"))


def test_cases_are_covered():
    """Every conversion compared on at least 20 inputs, every option exercised."""
    for conversion in CONVERSIONS:
        assert SEEN[conversion] >= 20, (conversion, SEEN[conversion])
    for option in ("CREATE_INDEX", "CREATE_MD5_FILE", "COMPRESSION_LEVEL", "MAX_RECORDS_IN_RAM",
                   "VALIDATION_STRINGENCY"):
        assert SEEN[option] >= 1, option

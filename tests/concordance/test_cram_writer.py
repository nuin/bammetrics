"""The CRAM bammetrics writes (crates/core/src/cram_writer.rs), read by samtools, by Picard
3.5.0 (htsjdk 5.0.0) and by bammetrics.

Tier B (docs/concordance.md): CRAM bytes are never compared. Each case writes one stream
twice with bammetrics, as a BAM and as a CRAM, and holds the CRAM to two truths.

1. The BAM of the same stream, per reader (`against_bam`).
   * samtools with `decode_md=0`: every line equal. No rule.
   * Picard `ViewSam`: every field equal, and every tag but two kinds, each counted:
     `MD` or `NM` the BAM line lacks (htsjdk rebuilds them on a CRAM), and a `B:C`, `B:S`
     or `B:I` array read as `B:c`, `B:s`, `B:i` (htsjdk's CRAM reader drops the unsigned
     flag; samtools reads the unsigned subtype from the same file).
   * bammetrics `ViewSam` against Picard `ViewSam`, both on the CRAM: equal, but for the
     records of `READER_GAPS`, which differ for any writer's CRAM.
2. samtools' own CRAM of that BAM (`against_samtools_cram`), under samtools with its
   default MD and NM rebuild and under Picard: every field equal, the tags as a multiset.
   A reader does with bammetrics' file what it does with htslib's.

A comparison fails when a reader exits non-zero or the two sides hold another number of
records. A rule is applied on the fields of a record and counted; a record that differs
in any other way is listed by name with its fields.

Inputs: the foundation corpus, t1, the exome slice and the BAMs of cram_bams.py (edges, a
header without @HD, no record, every kind of slice, a shuffled stream, reads over a contig
end sorted and shuffled, and what a CRAM does not give back). Writers: SortSam on every
input in two orders; MergeSamFiles, MarkDuplicates, AddOrReplaceReadGroups, CleanSam,
DownsampleSam and the chain `write` step on three. The nine other tools on
`record_edit::Out`, each on a fixture, read by samtools. Also Picard's own CRAM of each
tool against bammetrics', error parity without a reference, the options, and the
departures.

Picard's own CRAM is a truth only where its writer keeps the record. It is compared on
inputs without a repeated read name and whose records name read groups of the header. A
self-merge, where it sets the mate-unmapped flag on unmapped records, and a read group the
header lacks, which it drops, are each held to the two tools' BAMs instead.
"""
import hashlib
import subprocess
from collections import Counter

import pytest

import cram_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

SILENT = "VALIDATION_STRINGENCY=SILENT"
FIELDS = ("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL")
UNSIGNED = {"B:C": "B:c", "B:S": "B:s", "B:I": "B:i"}
NO_REFERENCE = "A reference must be supplied that includes the reference sequence for"

# Records bammetrics' ViewSam prints otherwise than Picard's, from a CRAM of any writer.
# Not the writer's: a reader gap. The test fails when the set changes either way.
#   long_cigar: a CIGAR of 70,000 operators that ends the reference span of its slice.
#     htsjdk rebuilds NM against that span and leaves out the insertion after its last
#     base (34,999); bammetrics counts it (35,000). The CIGAR prints as htsjdk prints it.
READER_GAPS = {"edges": {"long_cigar"}, "no_hd": {"long_cigar"}}


def run(cmd, env=None):
    return subprocess.run([str(c) for c in cmd], capture_output=True, env=env)


def picard(tool, *args):
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return run(["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"], env)


def ours(*args):
    return run([BAMMETRICS, *args])


def lines_of(done):
    return [line for line in done.stdout.decode("latin-1").split("\n") if line]


class Decoded:
    """What one reader printed for one file: its exit code and its record lines."""

    def __init__(self, who, done, header=False):
        self.who = who
        self.code = done.returncode
        self.error = done.stderr.decode("latin-1")[-600:]
        lines = lines_of(done)
        self.header = [line for line in lines if line.startswith("@")] if header else []
        self.records = [line.split("\t") for line in lines if not line.startswith("@")]


def samtools_view(path, reference, rebuild=False):
    cmd = ["samtools", "view", "--no-PG", "-h", "-T", reference]
    if not rebuild:
        cmd += ["--input-fmt-option", "decode_md=0"]
    return Decoded("samtools", run([*cmd, path]), header=True)


def picard_view(path, reference):
    return Decoded("Picard", picard("ViewSam", f"I={path}", f"R={reference}", SILENT), header=True)


def ours_view(path, reference):
    return Decoded("bammetrics", ours("ViewSam", f"I={path}", f"R={reference}"), header=True)


def rebuilt(expected_tags, tag):
    """An MD or NM the expected line lacks: a CRAM reader made it."""
    return tag[:2] in ("MD", "NM") and not any(t[:2] == tag[:2] for t in expected_tags)


def signed(tag):
    """`tag` as htsjdk's CRAM reader gives an unsigned array: the subtype signed over the
    same bytes. Returns None for any other tag."""
    name, kind = tag[:2], tag[3:6]
    if kind not in UNSIGNED:
        return None
    width = {"B:C": 8, "B:S": 16, "B:I": 32}[kind]
    values = [int(v) for v in tag[7:].split(",") if v]
    wrapped = [v - (1 << width) if v >= 1 << (width - 1) else v for v in values]
    return ",".join([f"{name}:{UNSIGNED[kind]}", *(str(v) for v in wrapped)])


def same_records(what, expected, got, rules=(), multiset=False, allow=()):
    """Assert `got` holds the records of `expected`. Both readers must have exited 0 and
    hold the same number of records. `rules` names the tag rules a record may use
    ("rebuilt", "signed"); `multiset` compares the tags without their order; `allow` names
    records that may differ and must. Returns a Counter of the rules used."""
    for side in (expected, got):
        assert side.code == 0, f"{what}: {side.who} exited {side.code}: {side.error}"
    assert len(expected.records) == len(got.records), (
        f"{what}: {expected.who} {len(expected.records)} records, {got.who} {len(got.records)}")
    used, problems, differing = Counter(), [], set()
    for x, y in zip(expected.records, got.records):
        if x == y:
            continue
        bad = [f"{FIELDS[i]} {x[i][:50]!r} -> {y[i][:50]!r}" for i in range(11) if x[i] != y[i]]
        want, have, record_rules = list(x[11:]), list(y[11:]), Counter()
        if "signed" in rules:
            for at, tag in enumerate(want):
                if signed(tag) is not None and tag not in have and signed(tag) in have:
                    want[at] = signed(tag)
                    record_rules["signed"] += 1
        if "rebuilt" in rules:
            kept = [tag for tag in have if not rebuilt(x[11:], tag)]
            record_rules["rebuilt"] += len(have) - len(kept)
            have = kept
        if (sorted(want) != sorted(have)) if multiset else (want != have):
            bad.append(f"tags {[t[:40] for t in want if t not in have]} -> "
                       f"{[t[:40] for t in have if t not in want]}"
                       + ("" if sorted(want) != sorted(have) else " (order)"))
        if bad:
            differing.add(x[0])
            if x[0] not in allow:
                problems.append(f"{x[0][:40]}: " + "; ".join(bad))
        else:
            used += record_rules
    assert not problems, (f"{what}: {len(problems)} of {len(expected.records)} records differ:\n  "
                          + "\n  ".join(problems[:12]))
    assert differing == set(allow), (
        f"{what}: expected exactly {sorted(allow)} to differ, got {sorted(differing)}")
    return used


def header_lines(decoded):
    return [line for line in decoded.header if not line.startswith("@PG")]


class Written:
    """One stream written by bammetrics as a BAM and as a CRAM, and samtools' CRAM of the
    BAM."""

    def __init__(self, tmp, label, reference, write):
        self.label, self.reference = label, reference
        self.bam, self.cram, self.theirs = tmp / "out.bam", tmp / "out.cram", tmp / "samtools.cram"
        for out in (self.bam, self.cram):
            done = write(out)
            assert done.returncode == 0, f"{label} -> {out.name}: {done.stderr.decode()[-600:]}"
        assert self.cram.read_bytes()[:6] == b"CRAM\x03\x00", f"{label}: not a CRAM 3.0"
        made = run(["samtools", "view", "--no-PG", "-C", "-T", reference, "-o", self.theirs, self.bam])
        assert made.returncode == 0, made.stderr.decode()[-300:]

    def against_bam(self, gaps=()):
        """Truth 1. Returns the records, and the rules Picard's comparison used."""
        from_bam = samtools_view(self.bam, self.reference)
        from_cram = samtools_view(self.cram, self.reference)
        assert header_lines(from_cram) == header_lines(from_bam), f"{self.label}: header (samtools)"
        same_records(f"{self.label}: samtools, CRAM against BAM", from_bam, from_cram)
        picard_bam = picard_view(self.bam, self.reference)
        picard_cram = picard_view(self.cram, self.reference)
        assert header_lines(picard_cram) == header_lines(picard_bam), f"{self.label}: header (Picard)"
        used = same_records(f"{self.label}: Picard, CRAM against BAM", picard_bam, picard_cram,
                            rules=("rebuilt", "signed"))
        same_records(f"{self.label}: on the CRAM, bammetrics against Picard", picard_cram,
                     ours_view(self.cram, self.reference), allow=gaps)
        return len(from_bam.records), used

    def against_samtools_cram(self, default_decode_fails=False):
        """Truth 2: each reader on bammetrics' CRAM and on samtools' CRAM of the BAM."""
        mine = samtools_view(self.cram, self.reference, rebuild=True)
        theirs = samtools_view(self.theirs, self.reference, rebuild=True)
        if default_decode_fails:
            # htslib cannot rebuild MD over a deletion past a contig end, in its own file.
            assert mine.code != 0 and theirs.code != 0, f"{self.label}: samtools read it"
            assert "extends beyond slice reference extents" in theirs.error
        else:
            same_records(f"{self.label}: samtools, ours against samtools' CRAM", theirs, mine,
                         multiset=True)
        same_records(f"{self.label}: Picard, ours against samtools' CRAM",
                     picard_view(self.theirs, self.reference), picard_view(self.cram, self.reference),
                     multiset=True)


@pytest.fixture(scope="session")
def edge(tmp_path_factory):
    """The BAMs of cram_bams.py and their reference."""
    directory = tmp_path_factory.mktemp("cram_bams")
    return cram_bams.write_all(directory), cram_bams.reference(directory)


def source(name, corpus, exome, edge):
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa"
    if name == "exome":
        return exome / "exome.bam", exome / "ref.fa"
    if name == "t1":
        return CORPUS / "t1.bam", CORPUS / "reference.fa"
    return edge[0][name], edge[1]


SORT_INPUTS = ("foundation", "exome", "t1", "edges", "no_hd", "empty", "slices", "shuffled",
               "overhang", "overhang_shuffled")


@pytest.mark.parametrize("order", ("coordinate", "queryname"))
@pytest.mark.parametrize("name", SORT_INPUTS)
def test_sort_sam_cram_reads_as_the_bam(name, order, corpus, exome, edge, tmp_path):
    bam, reference = source(name, corpus, exome, edge)
    written = Written(tmp_path, f"SortSam {name} {order}", reference, lambda out: ours(
        "SortSam", f"I={bam}", f"O={out}", f"SO={order}", f"R={reference}"))
    records, used = written.against_bam(gaps=READER_GAPS.get(name, ()))
    # The one shape default samtools cannot read back, from its own CRAM either.
    written.against_samtools_cram(default_decode_fails=name.startswith("overhang"))
    RESULTS.append((f"cram-writer SortSam {name} {order}", records, records,
                    f"3 readers; Picard rebuilt {used['rebuilt']} MD/NM, signed {used['signed']} arrays"))


def tool_args(tool, bam, out, reference):
    common = [f"I={bam}", f"O={out}", f"R={reference}", SILENT]
    return {
        "MergeSamFiles": ["MergeSamFiles", f"I={bam}", *common],
        "MarkDuplicates": ["MarkDuplicates", *common, f"M={out}.metrics", "ASSUME_SORT_ORDER=coordinate"],
        "AddOrReplaceReadGroups": ["AddOrReplaceReadGroups", *common, "RGLB=l", "RGPL=ILLUMINA",
                                   "RGPU=u", "RGSM=s"],
        "CleanSam": ["CleanSam", *common],
        "DownsampleSam": ["DownsampleSam", f"I={bam}", f"O={out}", f"REFERENCE_SEQUENCE={reference}",
                          SILENT, "P=0.5", "RANDOM_SEED=7"],
    }[tool]


TOOLS = ("MergeSamFiles", "MarkDuplicates", "AddOrReplaceReadGroups", "CleanSam", "DownsampleSam")


@pytest.mark.parametrize("tool", TOOLS)
@pytest.mark.parametrize("name", ("foundation", "slices", "overhang"))
def test_each_tool_writes_a_cram_that_reads_as_its_bam(tool, name, corpus, exome, edge, tmp_path):
    if (tool, name) == ("CleanSam", "overhang"):
        # CleanSam refuses the deletion over the contig end before it writes anything.
        name = "overhang_clean"
    bam, reference = source(name, corpus, exome, edge)
    written = Written(tmp_path, f"{tool} {name}", reference,
                      lambda out: ours(*tool_args(tool, bam, out, reference)))
    records, used = written.against_bam()
    RESULTS.append((f"cram-writer {tool} {name}", records, records,
                    f"3 readers; Picard rebuilt {used['rebuilt']} MD/NM"))


@pytest.mark.parametrize("name", ("foundation", "shuffled", "overhang_shuffled", "empty"))
def test_the_chain_write_step_keeps_the_stream(name, corpus, exome, edge, tmp_path):
    bam, reference = source(name, corpus, exome, edge)
    written = Written(tmp_path, f"chain write {name}", reference, lambda out: ours(
        "chain", bam, "-R", reference, "+", "write", out))
    records, _ = written.against_bam()
    # No sort in between: the CRAM holds the input's records in the input's order.
    same_records(f"chain write {name}: the input", samtools_view(bam, reference),
                 samtools_view(written.cram, reference))
    RESULTS.append((f"cram-writer chain write {name}", records, records, "3 readers; input order"))


FIXTURE_REFERENCE = FIXTURES / "ref_100k.fa"
# The other tools that write through `record_edit::Out`, each with arguments that work on a
# fixture of tests/data (the golden suite runs the same ones to a BAM).
OTHER_TOOLS = {
    "ReplaceSamHeader": [f"I={FIXTURES / 'markdup.sam'}", f"HEADER={FIXTURES / 'insert_size.bam'}"],
    "SetNmMdAndUqTags": [f"I={FIXTURES / 'markdup.bam'}"],
    "SetNmAndUqTags": [f"I={FIXTURES / 'insert_size.bam'}"],
    "RevertOriginalBaseQualitiesAndAddMateCigar": [f"I={FIXTURES / 'revert.bam'}"],
    "GatherBamFiles": [f"I={FIXTURES / 'markdup.sam'}"],
    "FixMateInformation": [f"I={FIXTURES / 'mate_filter.bam'}"],
    "FilterSamReads": [f"I={FIXTURES / 'mate_filter.bam'}", "FILTER=excludeTagValues", "TAG=XZ",
                       "TAG_VALUE=a", "TAG_VALUE=b"],
    "PositionBasedDownsampleSam": [f"I={FIXTURES / 'markdup.bam'}", "F=0.5"],
    "ReorderSam": [f"I={FIXTURES / 'markdup.bam'}", f"SD={FIXTURES / 'chr1.dict'}"],
}


def other_tool_cram_reads_as_its_bam(tool, tmp):
    """`tool` on its fixture to a BAM and to a CRAM; samtools reads both to the same
    header and records. Returns the number of records. No JVM."""
    outs = {}
    for kind in ("bam", "cram"):
        outs[kind] = tmp / f"out.{kind}"
        done = ours(tool, *OTHER_TOOLS[tool], f"O={outs[kind]}", f"R={FIXTURE_REFERENCE}", SILENT)
        assert done.returncode == 0, f"{tool} -> {kind}: {done.stderr.decode('latin-1')[-400:]}"
    assert outs["cram"].read_bytes()[:6] == b"CRAM\x03\x00", f"{tool}: not a CRAM 3.0"
    from_bam = samtools_view(outs["bam"], FIXTURE_REFERENCE)
    from_cram = samtools_view(outs["cram"], FIXTURE_REFERENCE)
    assert header_lines(from_cram) == header_lines(from_bam), f"{tool}: header"
    assert from_bam.records, f"{tool}: no record written"
    same_records(f"{tool}: samtools, CRAM against BAM", from_bam, from_cram)
    return len(from_bam.records)


@pytest.mark.parametrize("tool", sorted(OTHER_TOOLS))
def test_the_other_record_editing_tools_write_a_cram(tool, tmp_path):
    """Read by samtools alone: test_each_tool_writes_a_cram_that_reads_as_its_bam holds the
    writer to all three readers through five tools."""
    records = other_tool_cram_reads_as_its_bam(tool, tmp_path)
    RESULTS.append((f"cram-writer {tool} fixture", records, records, "samtools"))


def test_the_index_a_tool_makes_by_default_is_left_out_of_a_cram(tmp_path):
    """RevertOriginalBaseQualitiesAndAddMateCigar makes an index and an MD5 by default. With
    a CRAM it writes the MD5 and no index, where Picard writes `<CRAM>.bai`; an index asked
    for is refused, as for every tool."""
    tool = "RevertOriginalBaseQualitiesAndAddMateCigar"
    other_tool_cram_reads_as_its_bam(tool, tmp_path)
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "out.bai", "out.bam", "out.bam.md5", "out.cram", "out.cram.md5"]
    asked = tmp_path / "asked.cram"
    done = ours(tool, *OTHER_TOOLS[tool], f"O={asked}", f"R={FIXTURE_REFERENCE}", "CREATE_INDEX=true")
    assert done.returncode == 1 and b"no index is written beside a CRAM" in done.stderr, done.stderr
    assert not asked.exists()
    RESULTS.append(("cram-writer default index", 1, 1, "MD5 and no index beside a CRAM"))


def test_what_a_cram_does_not_give_back(edge, tmp_path):
    """Each loss named in cram_writer.rs, and the same loss from samtools' own CRAM."""
    bams, reference = edge
    bam = bams["limits"]
    cram, theirs = tmp_path / "ours.cram", tmp_path / "samtools.cram"
    assert ours("chain", bam, "-R", reference, "+", "write", cram).returncode == 0
    assert run(["samtools", "view", "--no-PG", "-C", "-T", reference, "-o", theirs, bam]).returncode == 0
    original = {r[0]: r for r in samtools_view(bam, reference).records}
    for who, view in (("samtools", samtools_view), ("Picard", picard_view)):
        mine = {r[0]: r for r in view(cram, reference).records}
        other = {r[0]: r for r in view(theirs, reference).records}
        assert set(mine) == set(original) == set(other), who
        # `=` and `X`, and neighbouring M operators, are one M.
        for name in ("eq_x", "m_then_m"):
            assert mine[name][5] == other[name][5] == "30M", (who, name, mine[name][5])
            assert mine[name][9] == original[name][9], (who, name)
        # An unmapped record keeps neither its CIGAR nor its mapping quality.
        assert mine["unmapped_cigar"][4:6] == other["unmapped_cigar"][4:6] == ["0", "*"], who
        # The mate fields of an unpaired record: the reader's rule, for either file.
        assert mine["unpaired_mate_fields"][6:9] == other["unpaired_mate_fields"][6:9], who
        assert mine["unpaired_mate_fields"][6] == "*", who
    assert original["unpaired_mate_fields"][6:9] == ["e2", "77", "123"]
    # bammetrics' reader gives `*` and 0, as htsjdk's, and keeps the template length.
    back = {r[0]: r for r in ours_view(cram, reference).records}
    assert back["unpaired_mate_fields"][6:9] == ["*", "0", "123"]
    RESULTS.append(("cram-writer limits", 4, 4, "each loss equal to samtools' CRAM"))


def halves(bam, tmp):
    """`bam` split in two, the first of each pair and every other record, so a merge of
    the two holds each record once. Both keep the whole header: `samtools view -r` would
    cut the @RG lines of both files to one, and leave records naming a read group their
    header lacks."""
    one, rest = tmp / "one.bam", tmp / "rest.bam"
    made = run(["samtools", "view", "--no-PG", "-b", "-f", "64", "-U", rest, "-o", one, bam])
    assert made.returncode == 0, made.stderr.decode()[-300:]
    head = lines_of(run(["samtools", "view", "-H", "--no-PG", bam]))
    counts = []
    for half in (one, rest):
        assert lines_of(run(["samtools", "view", "-H", "--no-PG", half])) == head, half.name
        counts.append(int(run(["samtools", "view", "-c", half]).stdout))
    assert min(counts) > 0 and sum(counts) == int(run(["samtools", "view", "-c", bam]).stdout), counts
    return one, rest


@pytest.mark.parametrize("tool", ("SortSam", *TOOLS))
def test_picards_own_cram_holds_the_same_records(tool, corpus, tmp_path):
    """Both tools write a CRAM; samtools reads both. Picard's writer leaves MD and NM out
    and orders tags its own way, so the tags are compared as a multiset without them.
    MergeSamFiles merges two halves of the corpus: Picard's CRAM writer is not faithful
    where a read name repeats (test_picards_cram_of_repeated_names_is_not_its_bam)."""
    bam, reference = corpus / "corpus.bam", corpus / "ref.fa"
    outs = {}
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.cram"
        if tool == "SortSam":
            args = ["SortSam", f"I={bam}", f"O={out}", "SO=queryname", f"R={reference}", SILENT]
        elif tool == "MergeSamFiles":
            one, rest = halves(bam, tmp_path)
            args = ["MergeSamFiles", f"I={one}", f"I={rest}", f"O={out}", f"R={reference}", SILENT]
        else:
            args = tool_args(tool, bam, out, reference)
        done = start(*args)
        assert done.returncode == 0, f"{tool} {side}: {done.stderr.decode('latin-1')[-600:]}"
        outs[side] = samtools_view(out, reference)

    def without_md_nm(decoded):
        for record in decoded.records:
            record[11:] = sorted(tag for tag in record[11:] if tag[:2] not in ("MD", "NM"))
        return decoded

    assert header_lines(outs["bammetrics"]) == header_lines(outs["picard"]), f"{tool}: header"
    same_records(f"{tool}: Picard's CRAM against bammetrics'", without_md_nm(outs["picard"]),
                 without_md_nm(outs["bammetrics"]))
    count = len(outs["picard"].records)
    RESULTS.append((f"cram-writer {tool} against Picard's CRAM", count, count, "tags without MD, NM"))


def test_picards_cram_of_repeated_names_is_not_its_bam(corpus, tmp_path):
    """A file merged with itself holds every read name twice over. Both tools write the
    same BAM, and bammetrics' CRAM reads as that BAM. Picard's own CRAM does not: samtools
    reads the mate-unmapped flag (0x8) on unmapped records whose BAM record lacks it. So
    Picard's CRAM is no truth for such a file; the BAM is."""
    bam, reference = corpus / "corpus.bam", corpus / "ref.fa"
    decoded = {}
    for side, start in (("picard", picard), ("bammetrics", ours)):
        for kind in ("bam", "cram"):
            out = tmp_path / f"{side}.{kind}"
            done = start("MergeSamFiles", f"I={bam}", f"I={bam}", f"O={out}", f"R={reference}", SILENT)
            assert done.returncode == 0, f"{side} {kind}: {done.stderr.decode('latin-1')[-400:]}"
            decoded[side, kind] = samtools_view(out, reference)
    truth = decoded["picard", "bam"]
    same_records("self-merge: the two BAMs", truth, decoded["bammetrics", "bam"])
    same_records("self-merge: bammetrics' CRAM against Picard's BAM", truth, decoded["bammetrics", "cram"])
    theirs = decoded["picard", "cram"]
    assert theirs.code == 0 and len(theirs.records) == len(truth.records)
    changed = [(x, y) for x, y in zip(truth.records, theirs.records) if x[:11] != y[:11]]
    assert changed, "Picard's CRAM of a self-merge now reads as its BAM: use it as a truth again"
    for x, y in changed:
        assert x[:1] + x[2:11] == y[:1] + y[2:11], f"{x[0]}: more than FLAG differs"
        assert int(y[1]) == int(x[1]) | 0x8 and int(x[1]) & 0x4, f"{x[0]}: FLAG {x[1]} -> {y[1]}"
    count = len(truth.records)
    RESULTS.append(("cram-writer self-merge", count, count,
                    f"ours = both BAMs; Picard's CRAM sets 0x8 on {len(changed)} records"))


def test_picards_cram_drops_a_read_group_the_header_lacks(edge, tmp_path):
    """A record whose RG names no @RG line. Picard's BAM keeps the tag and Picard's CRAM
    does not; bammetrics' CRAM keeps it, as its BAM and Picard's BAM do."""
    bams, reference = edge
    tags = {}
    for side, start in (("picard", picard), ("bammetrics", ours)):
        for kind in ("bam", "cram"):
            out = tmp_path / f"{side}.{kind}"
            done = start("SortSam", f"I={bams['edges']}", f"O={out}", "SO=coordinate",
                         f"R={reference}", SILENT)
            assert done.returncode == 0, f"{side} {kind}: {done.stderr.decode('latin-1')[-400:]}"
            found = [r for r in samtools_view(out, reference).records if r[0] == "rg_unknown"]
            assert len(found) == 1, (side, kind, len(found))
            tags[side, kind] = [tag for tag in found[0][11:] if tag.startswith("RG:")]
    assert tags["picard", "bam"] == tags["bammetrics", "bam"] == ["RG:Z:nope"], tags
    assert tags["bammetrics", "cram"] == ["RG:Z:nope"], tags
    assert tags["picard", "cram"] == [], "Picard's CRAM now keeps a read group its header lacks"
    RESULTS.append(("cram-writer RG not in the header", 1, 1, "kept, as in both BAMs; Picard's CRAM drops it"))


def test_a_placed_record_needs_the_reference_in_both(corpus, edge, tmp_path):
    """htsjdk fails at the first record with a reference sequence; unplaced records alone
    need no reference. bammetrics gives htsjdk's sentence and leaves no file."""
    bam = corpus / "corpus.bam"
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.cram"
        done = start("SortSam", f"I={bam}", f"O={out}", "SO=queryname", SILENT)
        assert done.returncode == 1, f"{side}: exit {done.returncode}"
        assert NO_REFERENCE in done.stderr.decode("latin-1"), side
    assert not (tmp_path / "bammetrics.cram").exists()
    assert (tmp_path / "picard.cram").stat().st_size == 0  # Picard leaves an empty file.
    unplaced = tmp_path / "unplaced.bam"
    made = run(["samtools", "view", "--no-PG", "-b", "-f", "12", "-o", unplaced, bam])
    assert made.returncode == 0
    count = int(run(["samtools", "view", "-c", unplaced]).stdout)
    assert count > 0
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.unplaced.cram"
        done = start("SortSam", f"I={unplaced}", f"O={out}", "SO=queryname", SILENT)
        assert done.returncode == 0, f"{side}: {done.stderr.decode('latin-1')[-400:]}"
    reference = corpus / "ref.fa"
    # Picard's CRAM writer orders the tags its own way: compared as a multiset.
    same_records("unplaced without a reference", samtools_view(tmp_path / "picard.unplaced.cram", reference),
                 samtools_view(tmp_path / "bammetrics.unplaced.cram", reference), multiset=True)
    RESULTS.append(("cram-writer no reference", count, count, "both fail placed, both write unplaced"))


def test_options_and_departures(corpus, tmp_path):
    bam, reference = corpus / "corpus.bam", corpus / "ref.fa"
    base = [f"I={bam}", "SO=coordinate", f"R={reference}", SILENT]
    # CREATE_MD5_FILE: the hex MD5 of the file, as Picard names it.
    out = tmp_path / "md5.cram"
    assert ours("SortSam", f"O={out}", *base, "CREATE_MD5_FILE=true").returncode == 0
    assert (tmp_path / "md5.cram.md5").read_text() == hashlib.md5(out.read_bytes()).hexdigest()
    # COMPRESSION_LEVEL changes the bytes and not the records; Picard ignores it for a CRAM.
    sizes = {}
    for level in (0, 9):
        out = tmp_path / f"level{level}.cram"
        assert ours("SortSam", f"O={out}", *base, f"COMPRESSION_LEVEL={level}").returncode == 0
        sizes[level] = out.stat().st_size
        same_records(f"COMPRESSION_LEVEL={level}", samtools_view(tmp_path / "md5.cram", reference),
                     samtools_view(out, reference))
    assert sizes[9] < sizes[0], sizes
    # A departure: Picard writes <CRAM>.bai, a BAI over slices; bammetrics refuses.
    theirs = tmp_path / "picard.cram"
    assert picard("SortSam", f"O={theirs}", *base, "CREATE_INDEX=true").returncode == 0
    assert (tmp_path / "picard.cram.bai").exists()
    refused = ours("SortSam", f"O={tmp_path / 'index.cram'}", *base, "CREATE_INDEX=true")
    assert refused.returncode == 1 and b"no index is written beside a CRAM" in refused.stderr
    assert not (tmp_path / "index.cram").exists()
    # The name decides, lower case, as in htsjdk: out.CRAM is a BAM in both.
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.CRAM"
        assert start("SortSam", f"O={out}", *base).returncode == 0, side
        assert out.read_bytes()[:4] == b"\x1f\x8b\x08\x04", side
    RESULTS.append(("cram-writer options", 5, 5, "md5, level, index refused, .CRAM is a BAM"))


def read_fasta(path):
    """[(name, bases)] of a FASTA."""
    contigs = []
    for line in path.read_text().split("\n"):
        if line.startswith(">"):
            contigs.append([line[1:].split()[0], []])
        elif line:
            contigs[-1][1].append(line)
    return [(name, "".join(lines)) for name, lines in contigs]


def write_fasta(path, contigs):
    with open(path, "w") as out:
        for name, bases in contigs:
            out.write(f">{name}\n")
            out.writelines(bases[at:at + 70] + "\n" for at in range(0, len(bases), 70))
    assert run(["samtools", "faidx", path]).returncode == 0
    return path


ABSENT = ("A reference must be supplied (reference sequence SAMSequenceRecord(name=chrM,length=12000,"
          "dict_index=2,assembly=null,alternate_names=[]) not found).")


def test_the_fasta_holds_the_headers_contigs_by_name(corpus, tmp_path):
    """htsjdk finds each contig of the header in the FASTA by its name. Both tools write
    against a FASTA in another order, one with more contigs, one whose contig is longer
    than its @SQ line and one that lacks a contig no record lies on; both fail at the first
    record on a contig the FASTA lacks, with the same sentence."""
    bam, reference = corpus / "corpus.bam", corpus / "ref.fa"
    contigs = read_fasta(reference)
    assert [name for name, _ in contigs] == ["chr1", "chr2", "chrM"]
    no_m = tmp_path / "no_m.bam"
    assert run(["samtools", "view", "--no-PG", "-b", "-e", 'rname != "chrM"', "-o", no_m, bam]).returncode == 0
    cases = {
        "order": (bam, contigs[::-1]),
        "more": (bam, [("chrFirst", "ACGT" * 500), *contigs, ("chrLast", "TTGCA" * 300)]),
        "longer": (bam, [(name, bases + "ACGTA" if name == "chr2" else bases) for name, bases in contigs]),
        "absent": (no_m, contigs[:2]),
    }
    total = 0
    for case, (source_bam, fasta_contigs) in cases.items():
        fasta = write_fasta(tmp_path / f"{case}.fa", fasta_contigs)
        outs = {}
        for side, start in (("picard", picard), ("bammetrics", ours)):
            outs[side] = tmp_path / f"{case}.{side}.cram"
            done = start("SortSam", f"I={source_bam}", f"O={outs[side]}", "SO=coordinate", f"R={fasta}", SILENT)
            assert done.returncode == 0, f"{case}, {side}: {done.stderr.decode('latin-1')[-400:]}"
        as_bam = tmp_path / f"{case}.bam"
        assert ours("SortSam", f"I={source_bam}", f"O={as_bam}", "SO=coordinate").returncode == 0
        from_bam = samtools_view(as_bam, fasta)
        assert from_bam.records
        same_records(f"{case}: samtools, CRAM against BAM", from_bam, samtools_view(outs["bammetrics"], fasta))
        same_records(f"{case}: Picard, CRAM against BAM", picard_view(as_bam, fasta),
                     picard_view(outs["bammetrics"], fasta), rules=("rebuilt", "signed"))
        # Picard's own CRAM against that FASTA holds the same records: its writer leaves MD
        # and NM out and orders the tags its own way.
        theirs, mine = samtools_view(outs["picard"], fasta), samtools_view(outs["bammetrics"], fasta)
        for decoded in (theirs, mine):
            for record in decoded.records:
                record[11:] = sorted(tag for tag in record[11:] if tag[:2] not in ("MD", "NM"))
        same_records(f"{case}: Picard's CRAM against bammetrics'", theirs, mine)
        total += len(from_bam.records)
    # A record on chrM, which the FASTA lacks.
    fasta = tmp_path / "absent.fa"
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"fails.{side}.cram"
        done = start("SortSam", f"I={bam}", f"O={out}", "SO=coordinate", f"R={fasta}", SILENT)
        assert done.returncode == 1, f"{side}: exit {done.returncode}"
        assert ABSENT in done.stderr.decode("latin-1"), f"{side}: {done.stderr.decode('latin-1')[-400:]}"
    assert not (tmp_path / "fails.bammetrics.cram").exists()
    RESULTS.append(("cram-writer FASTA by name", total, total, "order, more, longer, absent; both fail on a used absent contig"))


def test_a_tag_text_with_bytes_outside_printable_ascii(tmp_path):
    """A BAM whose `Z` value holds UTF-8 bytes or a control byte, and whose read name holds
    UTF-8 bytes (samtools makes it: both tools read SAM text one byte per character). Both
    tools write the CRAM. samtools reads bammetrics' CRAM to the bytes of its BAM, and
    bammetrics reads its CRAM as its BAM. Picard's own CRAM holds another `XZ` text."""
    head = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
    body = ("r1\t0\tchr1\t10\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tXZ:Z:café 中\tXY:Z:a\x01b\n"
            "café\t0\tchr1\t20\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tXZ:Z:plain\n")
    sam, bam = tmp_path / "text.sam", tmp_path / "text.bam"
    sam.write_bytes((head + body).encode())
    assert run(["samtools", "view", "--no-PG", "-b", "-o", bam, sam]).returncode == 0
    outs = {}
    for side, start in (("picard", picard), ("bammetrics", ours)):
        for kind in ("bam", "cram"):
            outs[side, kind] = tmp_path / f"{side}.{kind}"
            done = start("SortSam", f"I={bam}", f"O={outs[side, kind]}", "SO=coordinate",
                         f"R={FIXTURE_REFERENCE}", SILENT)
            assert done.returncode == 0, f"{side} {kind}: {done.stderr.decode('latin-1')[-400:]}"
    from_bam = samtools_view(outs["bammetrics", "bam"], FIXTURE_REFERENCE)
    from_cram = samtools_view(outs["bammetrics", "cram"], FIXTURE_REFERENCE)
    same_records("tag text: samtools, CRAM against BAM", from_bam, from_cram)
    # samtools prints the bytes; `lines_of` reads them one character per byte.
    expected = [field.encode().decode("latin-1") for field in ("XZ:Z:café 中", "XY:Z:a\x01b")]
    assert from_cram.records[0][11:] == expected, from_cram.records[0][11:]
    assert from_cram.records[1][0] == "café".encode().decode("latin-1")
    same_records("tag text: bammetrics, CRAM against BAM",
                 Decoded("bammetrics", ours("ViewSam", f"I={outs['bammetrics', 'bam']}"), header=True),
                 ours_view(outs["bammetrics", "cram"], FIXTURE_REFERENCE), rules=("rebuilt",))
    theirs = [tag for tag in samtools_view(outs["picard", "cram"], FIXTURE_REFERENCE).records[0][11:]
              if tag.startswith("XZ:")]
    assert theirs and theirs != expected[:1], "Picard's CRAM now keeps the bytes of the tag"
    RESULTS.append(("cram-writer tag text bytes", 2, 2, f"kept as in the BAM; Picard's CRAM holds {theirs[0]!r}"))


def test_a_mate_reference_the_header_lacks_fails_in_both(tmp_path):
    """A record whose mate reference index is not in the header: a CRAM of it is one no
    reader opens. Both tools fail, and bammetrics leaves no file."""
    import struct

    import index_bams

    forged = bytearray(index_bams.record("forged", 0, 100, flag=1))
    struct.pack_into("<i", forged, 24, 7)
    bam = tmp_path / "forged.bam"
    bam.write_bytes(index_bams.bam([("chr1", 100000)], [bytes(forged)], order="unsorted"))
    for side, start in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.cram"
        done = start("SortSam", f"I={bam}", f"O={out}", "SO=queryname", f"R={FIXTURE_REFERENCE}", SILENT)
        assert done.returncode == 1, f"{side}: exit {done.returncode}: {done.stderr.decode('latin-1')[-300:]}"
    assert b"record forged: mate reference index 7 is not in the header" in done.stderr
    assert not (tmp_path / "bammetrics.cram").exists()
    RESULTS.append(("cram-writer mate reference outside the header", 1, 1, "both fail"))


def test_two_records_picard_refuses_and_bammetrics_writes(tmp_path):
    """KNOWN_DIFFERENCES.md, CRAM output: a read that starts past the end of its contig and
    an unmapped record with a CIGAR. Picard's CRAM writer fails on each. bammetrics writes
    both, as samtools does: samtools reads its CRAM as it reads samtools' own."""
    head = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n"
    cases = {
        "past": "r1\t0\tchr1\t100001\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
        "cigar": "r1\t4\tchr1\t10\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n",
    }
    for case, line in cases.items():
        sam, bam = tmp_path / f"{case}.sam", tmp_path / f"{case}.bam"
        sam.write_text(head + line)
        assert run(["samtools", "view", "--no-PG", "-b", "-o", bam, sam]).returncode == 0
        theirs = picard("SortSam", f"I={bam}", f"O={tmp_path / f'{case}.picard.cram'}", "SO=queryname",
                        f"R={FIXTURE_REFERENCE}", SILENT)
        assert theirs.returncode == 1, f"{case}: Picard now writes it"
        out, own = tmp_path / f"{case}.cram", tmp_path / f"{case}.samtools.cram"
        done = ours("SortSam", f"I={bam}", f"O={out}", "SO=queryname", f"R={FIXTURE_REFERENCE}", SILENT)
        assert done.returncode == 0, done.stderr.decode("latin-1")[-300:]
        assert run(["samtools", "view", "--no-PG", "-C", "-T", FIXTURE_REFERENCE, "-o", own, bam]).returncode == 0
        same_records(f"{case}: samtools, ours against samtools' CRAM",
                     samtools_view(own, FIXTURE_REFERENCE), samtools_view(out, FIXTURE_REFERENCE))
    RESULTS.append(("cram-writer records Picard refuses", 2, 2, "written; read as samtools' own CRAM"))


def sam_to_bam(tmp_path, name, text):
    sam, bam = tmp_path / f"{name}.sam", tmp_path / f"{name}.bam"
    sam.write_bytes(text.encode())
    done = run(["samtools", "view", "--no-PG", "-b", "-o", bam, sam])
    assert done.returncode == 0, done.stderr.decode()[-300:]
    return bam


FIXTURE_HEAD = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"


def test_a_container_of_fewer_than_four_qualities(tmp_path):
    """rANS 4x8 of order 1 takes no block under 4 bytes, so a quality block of 1 to 3
    bytes is coded with order 0 (vendor/README.md, the eleventh change). Both tools write
    one read of 1, 2 or 3 bases, three reads of 1 base, and a tail container of one 1-base
    read after 10,000 reads of 3 bases; samtools and Picard read bammetrics' CRAM as the
    BAM."""
    line = "r{0}\t0\tchr1\t{1}\t60\t{2}M\t*\t0\t0\t{3}\t{4}\n"
    cases = {
        "1x1": [line.format(0, 100, 1, "A", "I")],
        "1x2": [line.format(0, 100, 2, "AC", "I#")],
        "1x3": [line.format(0, 100, 3, "ACG", "I#5")],
        "3x1": [line.format(i, 100 + i, 1, "A", "I") for i in range(3)],
        "tail": [line.format(i, 100 + i // 10, 3, "ACG", "I#5") for i in range(10_000)]
                + [line.format("t", 2000, 1, "A", "I")],
    }
    for case, lines in cases.items():
        bam = sam_to_bam(tmp_path, case, FIXTURE_HEAD + "".join(lines))
        for side, start in (("picard", picard), ("bammetrics", ours)):
            out = tmp_path / f"{case}.{side}.cram"
            done = start("SortSam", f"I={bam}", f"O={out}", "SO=coordinate", f"R={FIXTURE_REFERENCE}", SILENT)
            assert done.returncode == 0, f"{case}, {side}: {done.stderr.decode('latin-1')[-300:]}"
        mine = tmp_path / f"{case}.bammetrics.cram"
        same_records(f"{case}: samtools, CRAM against BAM", samtools_view(bam, FIXTURE_REFERENCE),
                     samtools_view(mine, FIXTURE_REFERENCE))
        same_records(f"{case}: Picard, CRAM against BAM", picard_view(bam, FIXTURE_REFERENCE),
                     picard_view(mine, FIXTURE_REFERENCE), rules=("rebuilt",))
    RESULTS.append(("cram-writer quality blocks under 4 bytes", len(cases), len(cases),
                    "written by both; read as the BAM"))


def test_an_operator_of_length_0_is_left_out(tmp_path):
    """docs/07 D223: bammetrics' CRAM leaves out every operator of length 0, so each
    CIGAR below reads back `10M` under samtools and Picard. Picard's CRAM leaves out `0I`
    and `0H` alone: samtools reads its `5M0D5M` as `5M5M`, its `0S10M` as `0S10M`, and
    fails on the slice of its `10M0S` (`Feature position 11 after end of read`)."""
    cigars = ("5M0I5M", "5M0D5M", "0S10M", "5M0N5M", "0H10M", "5M0P5M", "10M0S")
    picards = {"5M0I5M": "10M", "5M0D5M": "5M5M", "0S10M": "0S10M", "10M0S": None}
    read = {}
    for i, cigar in enumerate(cigars):
        line = f"z{i}\t0\tchr1\t{100 + i}\t60\t{cigar}\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n"
        bam = sam_to_bam(tmp_path, f"zero{i}", FIXTURE_HEAD + line)
        for side, start in (("picard", picard), ("bammetrics", ours)):
            if side == "picard" and cigar not in picards:
                continue
            out = tmp_path / f"zero{i}.{side}.cram"
            done = start("SortSam", f"I={bam}", f"O={out}", "SO=coordinate", f"R={FIXTURE_REFERENCE}", SILENT)
            assert done.returncode == 0, f"{cigar}, {side}: {done.stderr.decode('latin-1')[-300:]}"
            decoded = samtools_view(out, FIXTURE_REFERENCE)
            read[side, cigar] = [r[5] for r in decoded.records][0] if decoded.code == 0 else None
            if side == "bammetrics":
                assert [r[5] for r in picard_view(out, FIXTURE_REFERENCE).records] == ["10M"], cigar
    assert all(read["bammetrics", cigar] == "10M" for cigar in cigars), read
    assert {cigar: read["picard", cigar] for cigar in picards} == picards, read
    RESULTS.append(("cram-writer operators of length 0", len(cigars), len(cigars),
                    "left out; Picard's CRAM keeps 0D, 0N, 0P and 0S"))


def test_a_record_without_bases_over_16_mi_read_bases_fails(tmp_path):
    """docs/07 D223: SEQ `*` with a CIGAR over 16 Mi read bases fails bammetrics' CRAM
    writer, naming the record, and leaves no file. One of 1,001 read bases is written and
    samtools reads it back as the BAM."""
    fine = "ok\t0\tchr1\t100\t60\t1000S1M\t*\t0\t0\t*\t*\n"
    forged = "forged\t0\tchr1\t100\t60\t268435455S1M\t*\t0\t0\t*\t*\n"
    bam = sam_to_bam(tmp_path, "fine", FIXTURE_HEAD + fine)
    out = tmp_path / "fine.cram"
    done = ours("SortSam", f"I={bam}", f"O={out}", "SO=coordinate", f"R={FIXTURE_REFERENCE}")
    assert done.returncode == 0, done.stderr.decode()[-300:]
    same_records("SEQ *: samtools, CRAM against BAM", samtools_view(bam, FIXTURE_REFERENCE),
                 samtools_view(out, FIXTURE_REFERENCE))
    bam = sam_to_bam(tmp_path, "forged", FIXTURE_HEAD + fine + forged)
    out = tmp_path / "forged.cram"
    done = ours("SortSam", f"I={bam}", f"O={out}", "SO=coordinate", f"R={FIXTURE_REFERENCE}")
    assert done.returncode == 1
    assert b"record forged: SEQ * with a CIGAR over 268435456 read bases" in done.stderr, done.stderr
    assert not out.exists()
    RESULTS.append(("cram-writer SEQ * over 16 Mi bases", 1, 1, "refused, no file; Picard writes it"))


def test_header_text_outside_ascii(tmp_path):
    """docs/07 D218: the CRAM header holds the header as htsjdk holds it (one character
    per byte of the BAM, or per character of SAM text) encoded as UTF-8, and a CRAM
    header is read back as UTF-8. Picard's CRAM of the BAM and bammetrics' read alike under
    samtools, and Picard and bammetrics read each CRAM as the BAM. Converted from SAM text,
    both CRAMs hold the UTF-8 of the text, and both tools write the SAM text back. A
    character above U+00FF in SAM text is test_sam_format_converter.py's
    test_header_text_outside_ascii."""
    head = FIXTURE_HEAD + "@RG\tID:g1\tSM:s\tDS:café über\n@CO\tcafé\n"
    body = "r1\t0\tchr1\t10\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\tRG:Z:g1\n"
    bam = sam_to_bam(tmp_path, "header", head + body)
    sam = tmp_path / "header.sam"
    sam.write_text(head + body, encoding="utf-8")
    for source in (bam, sam):
        outs = {}
        for who, run_tool in (("ours", ours), ("picard", picard)):
            out = tmp_path / f"{who}.{source.suffix[1:]}.cram"
            done = run_tool("SortSam", f"I={source}", f"O={out}", "SO=coordinate",
                            f"R={FIXTURE_REFERENCE}")
            assert done.returncode == 0, done.stderr.decode()[-300:]
            outs[who] = out
        samtools_side = {who: header_lines(samtools_view(out, FIXTURE_REFERENCE))
                         for who, out in outs.items()}
        assert samtools_side["ours"] == samtools_side["picard"], source.name
        source_lines = header_lines(samtools_view(source, FIXTURE_REFERENCE))
        for who_reads, view in (("bammetrics", ours_view), ("Picard", picard_view)):
            for who, out in outs.items():
                read = header_lines(view(out, FIXTURE_REFERENCE))
                expected = header_lines(view(source, FIXTURE_REFERENCE))
                assert read == expected, f"{who_reads} reading {who}'s CRAM of {source.name}"
        if source == sam:
            assert "@CO\t" + "café".encode().decode("latin-1") in samtools_side["ours"]
        else:
            assert samtools_side["ours"] != source_lines
    RESULTS.append(("cram-writer header text outside ASCII", 1, 1,
                    "both CRAMs alike under samtools, each read as its source by both tools"))

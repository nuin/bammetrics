"""ViewSam against Picard 3.5.0: the SAM text on standard output, byte for byte.

Inputs: every BAM of tests/data, the foundation corpus, the corpus tiers under 64 MiB,
the exome slice, the foundation corpus and t1 as CRAM 3.0 and as SAM text (plain,
bgzip, gzip, standard input), and hand-written SAM files that hold every tag type,
every `*` form, CR and CRLF line ends and headers without `@HD`.

Options: the nine `ALIGNMENT_STATUS` by `PF_STATUS` pairs, `HEADER_ONLY`,
`RECORDS_ONLY`, `INTERVAL_LIST` on a BAM with its index (Picard queries the index),
on a BAM without one and on SAM text (Picard scans), with four interval lists.

Picard reads a BAM or a CRAM with `VALIDATION_STRINGENCY=SILENT` here, because
bammetrics validates no BAM record (KNOWN_DIFFERENCES.md). SAM text is read at
`STRICT`, `LENIENT` and `SILENT` by both tools: a bad line must fail in both or print
the same text in both.

Also here, since they share the SAM text writer and reader: `bammetrics chain BAM +
write x.sam` against Picard `SamFormatConverter`, whole files byte for byte, and each
metrics tool on SAM text against itself on the BAM.

Error parity: `HEADER_ONLY` with `RECORDS_ONLY`, `O=`, `ASSUME_SORTED=`, an enum
value in the wrong case, a missing input, a missing interval list, a CRAM without
`R=`, a text file that is no SAM.

`test_cases_are_covered` requires every tag type and every filter value in the
compared output.
"""
import gzip
import shutil
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "ViewSam"
VERB = "view-sam"
SILENT = "VALIDATION_STRINGENCY=SILENT"
STRINGENCIES = ("STRICT", "LENIENT", "SILENT")
TIER_BYTES = 64 << 20
# Tag types (`A`, `i`, `B:c`, ...) and options seen in compared Picard output.
SEEN_TAGS = set()
SEEN_OPTIONS = set()

HEADER = (
    "@HD\tVN:1.5\tSO:unsorted\n"
    "@SQ\tSN:chr1\tLN:100000\tM5:abc\tUR:file:/ref.fa\n"
    "@SQ\tSN:chr2\tLN:5000\n"
    "@RG\tID:rg1\tSM:s1\tLB:lib1\tPL:ILLUMINA\tPU:fc.1\n"
    "@RG\tPL:ILLUMINA\tID:rg2\tSM:s2\n"
    "@PG\tPN:bwa\tID:bwa\tVN:0.7.17\tCL:bwa mem -t 4 ref.fa r1.fq r2.fq\n"
    "@CO\ta comment\twith a tab\n"
)
# Records that are valid for htsjdk at STRICT: every tag type, the `*` forms, a flag
# in each spelling, a mate on another reference, clips, pads and skips.
VALID = [
    "r01\t99\tchr1\t100\t60\t5S20M\t=\t300\t225\tACGTACGTACGTACGTACGTACGTA\tIIIIIIIIIIIIIIIIIIIIIIIII\t"
    "RG:Z:rg1\tNM:i:0\tXA:A:q\tXF:f:1.5\tXZ:Z:a:b c\tXE:Z:\tAS:i:-7",
    "r01\t147\tchr1\t300\t60\t25M\t=\t100\t-225\tacgtnACGTN.=ACGTACGTACGTA\t!\"#$%&'()*+,-./0123456789\t"
    "RG:Z:rg1\tXB:B:c,-128,127\tXC:B:C,0,255\tXS:B:s,-32768,32767\tXT:B:S,0,65535\t"
    "XI:B:i,-2147483648,2147483647\tXJ:B:I,0,4294967295\tXG:B:f,0.1,-3,1e10,3.4028235E38\tXN:B:c",
    "r02\t0x41\tchr1\t500\t0\t3H10M2I5M3D5M4N5M1P5M2S\tchr2\t7\t0\tACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\t"
    "XH:H:1AE301\tXU:i:4294967295\tXL:i:-2147483648\tXK:i:65535\tXM:i:65536\tXO:i:255\tXP:i:256\tXQ:i:-129",
    "r02\t0201\tchr2\t7\t255\t10=2X\tchr1\t500\t0\tACGTACGTACGT\tIIIIIIIIIIII\tXF:f:NaN\tXG:f:Infinity\tXD:f:1e-10\tXW:f:100000000",
    "r03\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*",
    "r04\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\tRG:Z:rg2",
    "r05\tuUp1\t*\t0\t0\t*\t*\t0\t0\tAC\tII",
    "r05\t141\t*\t0\t0\t*\t*\t0\t0\tGT\t~~",
    "r06\t69\tchr2\t40\t0\t*\t=\t40\t0\tACGT\tIIII\tXZ:Z:placed unmapped",
    "r06\t137\tchr2\t40\t9\t4M\t=\t40\t0\tACGT\tIIII",
    "r07\t16\tchr1\t99990\t3\t11M\t*\t0\t0\tACGTACGTACG\tIIIIIIIIIII\tHI:i:2\tXX:Z:last base of chr1",
    "r08\t256\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXZ:Z:repeated\tXZ:Z:kept",
    "r09\t2048\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII",
    "r10\t512\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII",
    "r11\t1024\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII",
]
# One bad line each: (name, line). htsjdk fails some at every stringency.
BAD = [
    ("fields", "b01\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT"),
    ("flag_text", "b02\t12a\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("flag_letter", "b03\tq\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("pos_text", "b04\t0\tchr1\tseven\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("mapq_text", "b05\t0\tchr1\t7\tx\t4M\t*\t0\t0\tACGT\tIIII"),
    ("mpos_text", "b06\t0\tchr1\t7\t0\t4M\t*\tx\t0\tACGT\tIIII"),
    ("isize_text", "b07\t0\tchr1\t7\t0\t4M\t*\t0\tx\tACGT\tIIII"),
    ("no_rname_mapped", "b08\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"),
    ("pos_zero", "b09\t0\tchr1\t0\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("cigar_star_mapped", "b10\t0\tchr1\t7\t0\t*\t*\t0\t0\tACGT\tIIII"),
    ("unmapped_pos", "b11\t4\t*\t7\t0\t*\t*\t0\t0\tACGT\tIIII"),
    ("unmapped_mapq", "b12\t4\t*\t0\t9\t*\t*\t0\t0\tACGT\tIIII"),
    ("mate_unpaired", "b14\t0\tchr1\t7\t0\t4M\tchr1\t9\t0\tACGT\tIIII"),
    ("mate_missing", "b15\t65\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("mpos_zero", "b16\t65\tchr1\t7\t0\t4M\t=\t0\t0\tACGT\tIIII"),
    ("mpos_no_mate", "b17\t73\tchr1\t7\t0\t4M\t*\t9\t0\tACGT\tIIII"),
    ("isize_no_mate", "b18\t73\tchr1\t7\t0\t4M\t*\t0\t5\tACGT\tIIII"),
    ("qual_space", "b22\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tII I"),
    ("tag_fields", "b23\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX"),
    ("tag_type", "b24\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:q:1"),
    ("tag_int", "b25\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:i:1.5"),
    ("tag_int_range", "b26\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:i:4294967296"),
    ("tag_char", "b27\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:A:ab"),
    ("tag_float", "b28\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:f:x"),
    ("tag_hex", "b29\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:H:ABC"),
    ("tag_array_type", "b30\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:B:q,1"),
    ("tag_array_range", "b31\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXX:B:c,200"),
    ("tag_name", "b32\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tXXX:i:1"),
    ("rg_unknown", "b33\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:nope"),
    ("start_past_end", "b34\t4\tchr2\t9000\t0\t*\t*\t0\t0\tACGT\tIIII"),
    ("cigar_length", "b35\t0\tchr1\t7\t0\t5M\t*\t0\t0\tACGT\tIIII"),
    ("cigar_off_reference", "b38\t0\tchr2\t4999\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("cigar_adjacent_indels", "b39\t0\tchr1\t7\t0\t2M1I1I2M\t*\t0\t0\tACGTAC\tIIIIII"),
    ("cigar_clip_inside", "b40\t0\tchr1\t7\t0\t2M2S2M\t*\t0\t0\tACGTAC\tIIIIII"),
    ("proper_unpaired", "b41\t2\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("pair_no_end", "b42\t1\tchr1\t7\t0\t4M\t=\t9\t0\tACGT\tIIII"),
    ("unmapped_secondary", "b43\t260\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII"),
    ("empty_field", "b45\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII\t\tXX:i:1"),
]
# Bad lines htsjdk keeps as text under LENIENT or SILENT in a form a BAM record
# cannot hold: (name, line, what bammetrics prints for it, None when it fails).
# bammetrics reads them differently (KNOWN_DIFFERENCES.md); each is checked in
# test_lines_a_bam_cannot_hold.
UNHOLDABLE = [
    ("rname_unknown", "u01\t0\tchr9\t7\t0\t4M\t*\t0\t0\tACGT\tIIII", "u01\t0\t*\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("mate_unknown", "u02\t65\tchr1\t7\t0\t4M\tchr9\t9\t0\tACGT\tIIII",
     "u02\t65\tchr1\t7\t0\t4M\t*\t9\t0\tACGT\tIIII"),
    ("base", "u03\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACXT\tIIII", "u03\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACNT\tIIII"),
    ("mapq_256", "u04\t0\tchr1\t7\t256\t4M\t*\t0\t0\tACGT\tIIII", "u04\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIIII"),
    ("qual_length", "u05\t0\tchr1\t7\t0\t4M\t*\t0\t0\tACGT\tIII", None),
    ("qual_no_seq", "u06\t4\t*\t0\t0\t*\t*\t0\t0\t*\tI", None),
    ("cigar_text", "u07\t0\tchr1\t7\t0\tM4\t*\t0\t0\tACGT\tIIII", None),
    ("cigar_operator", "u08\t0\tchr1\t7\t0\t4Q\t*\t0\t0\tACGT\tIIII", None),
]


def env():
    e = picard_env()
    e.setdefault("JAVA_OPTS", "-Xmx8g")
    return e


def run_picard(args, tool=TOOL, stdin=None):
    return subprocess.run(["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"], env=env(),
                          capture_output=True, stdin=stdin)


def run_ours(args, face=TOOL, stdin=None):
    return subprocess.run([str(BAMMETRICS), face, *args], capture_output=True, stdin=stdin)


def note_tags(stdout):
    for line in stdout.split(b"\n"):
        if line.startswith(b"@"):
            continue
        for field in line.split(b"\t")[11:]:
            parts = field.split(b":", 2)
            if len(parts) == 3:
                kind = parts[1].decode()
                SEEN_TAGS.add(kind + ":" + parts[2][:1].decode() if kind == "B" else kind)


def first_difference(picard, ours):
    p, o = picard.split(b"\n"), ours.split(b"\n")
    for i, (a, b) in enumerate(zip(p, o)):
        if a != b:
            return f"line {i + 1}:\n  Picard:     {a[:300]!r}\n  bammetrics: {b[:300]!r}"
    return f"{len(p)} lines from Picard, {len(o)} from bammetrics"


def compare(label, args, picard_extra=(), ours_extra=(), stdin=None):
    """Both tools exit 0 and print identical standard output."""
    picard = run_picard([*args, *picard_extra], stdin=stdin() if stdin else None)
    ours = run_ours([*args, *ours_extra], stdin=stdin() if stdin else None)
    assert picard.returncode == 0, picard.stderr[-1500:].decode(errors="replace")
    assert ours.returncode == 0, ours.stderr.decode(errors="replace")
    assert picard.stdout == ours.stdout, f"{TOOL} {label}: {first_difference(picard.stdout, ours.stdout)}"
    note_tags(picard.stdout)
    SEEN_OPTIONS.update(a for a in args if not a.startswith(("I=", "R=", "INTERVAL_LIST=")))
    lines = picard.stdout.count(b"\n")
    RESULTS.append((f"{TOOL} {label}", lines, lines, "stdout identical"))
    return picard.stdout


def both_fail(label, args, picard_extra=(), stdin=None):
    picard = run_picard([*args, *picard_extra], stdin=stdin() if stdin else None)
    ours = run_ours(args, stdin=stdin() if stdin else None)
    assert picard.returncode != 0, f"{label}: Picard exits 0"
    assert ours.returncode != 0, f"{label}: bammetrics exits 0, Picard says {picard.stderr[-600:]!r}"
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))
    return picard, ours


def sam_text(bam, reference=None):
    cmd = ["samtools", "view", "-h", "--no-PG", str(bam)]
    if reference:
        cmd += ["-T", str(reference)]
    return subprocess.run(cmd, check=True, capture_output=True).stdout


def tier_bams():
    return sorted(p for p in CORPUS.glob("t*.bam") if p.stat().st_size < TIER_BYTES)


@pytest.fixture(scope="module")
def handwritten(tmp_path_factory):
    """The hand-written SAM: HEADER and the VALID records."""
    path = tmp_path_factory.mktemp("handwritten") / "valid.sam"
    path.write_text(HEADER + "".join(line + "\n" for line in VALID))
    return path


@pytest.fixture(scope="module")
def foundation_sam(corpus, tmp_path_factory):
    path = tmp_path_factory.mktemp("foundation") / "corpus.sam"
    path.write_bytes(sam_text(corpus / "corpus.bam"))
    return path


@pytest.mark.parametrize("bam", sorted(FIXTURES.glob("*.bam")), ids=lambda p: p.stem)
def test_fixtures(bam):
    compare(bam.name, [f"I={bam}"], picard_extra=[SILENT])


def test_foundation(corpus):
    compare("foundation", [f"I={corpus / 'corpus.bam'}"], picard_extra=[SILENT])


def test_exome(exome):
    compare("exome", [f"I={exome / 'exome.bam'}"], picard_extra=[SILENT])


@pytest.mark.parametrize("bam", tier_bams(), ids=lambda p: p.stem)
def test_tiers(bam):
    compare(bam.stem, [f"I={bam}"], picard_extra=[SILENT])


@pytest.mark.parametrize("alignment", ["Aligned", "Unaligned", "All"])
@pytest.mark.parametrize("pf", ["PF", "NonPF", "All"])
def test_status_filters(corpus, alignment, pf):
    compare(f"foundation {alignment} {pf}",
            [f"I={corpus / 'corpus.bam'}", f"ALIGNMENT_STATUS={alignment}", f"PF_STATUS={pf}"],
            picard_extra=[SILENT])


@pytest.mark.parametrize("option", ["HEADER_ONLY=true", "RECORDS_ONLY=true", "HEADER_ONLY=false",
                                    "RECORDS_ONLY=false"])
def test_header_and_records_only(corpus, handwritten, option):
    compare(f"foundation {option}", [f"I={corpus / 'corpus.bam'}", option], picard_extra=[SILENT])
    compare(f"handwritten {option}", [f"I={handwritten}", option])
    compare(f"markdup {option}", [f"I={FIXTURES / 'markdup.bam'}", option, "PF_STATUS=PF"],
            picard_extra=[SILENT])


def interval_lists(corpus, directory):
    """Four interval lists on the foundation reference: the corpus targets, one small
    interval, intervals that overlap, abut and come out of order, and none."""
    header = subprocess.run(["samtools", "view", "-H", str(corpus / "corpus.bam")], check=True,
                            capture_output=True, text=True).stdout
    sequences = [line for line in header.splitlines() if line.startswith("@SQ")]
    names = [dict(f.split(":", 1) for f in line.split("\t")[1:])["SN"] for line in sequences]
    lengths = [int(dict(f.split(":", 1) for f in line.split("\t")[1:])["LN"]) for line in sequences]
    head = "@HD\tVN:1.6\n" + "\n".join(sequences) + "\n"
    first, last = names[0], names[-1]
    lists = {
        "targets": (corpus / "target.interval_list").read_text(),
        "one": head + f"{first}\t1000\t1200\t+\tsmall\n",
        "mixed": head + (f"{last}\t500\t900\t-\tlate\n{first}\t3000\t3500\t+\ta\n{first}\t3400\t3600\t+\tb\n"
                         f"{first}\t3601\t3700\t+\tabuts\n{first}\t1\t50\t+\tstart\n"
                         f"{last}\t{lengths[-1] - 200}\t{lengths[-1]}\t+\tend\n"),
        "none": head,
    }
    paths = {}
    for name, text in lists.items():
        paths[name] = directory / f"{name}.interval_list"
        paths[name].write_text(text)
    return paths


@pytest.mark.parametrize("which", ["targets", "one", "mixed", "none"])
def test_interval_list(corpus, foundation_sam, tmp_path, which):
    """With the .bai beside the BAM Picard queries the index; without it, and on SAM
    text, it scans and stops after the last interval. The two paths differ on
    unmapped records placed with their mates."""
    intervals = interval_lists(corpus, tmp_path)[which]
    indexed = corpus / "corpus.bam"
    assert indexed.with_name("corpus.bam.bai").exists()
    bare = tmp_path / "bare.bam"
    shutil.copy(indexed, bare)
    by_index = compare(f"interval list {which}, index", [f"I={indexed}", f"INTERVAL_LIST={intervals}"],
                       picard_extra=[SILENT])
    by_scan = compare(f"interval list {which}, scan", [f"I={bare}", f"INTERVAL_LIST={intervals}"],
                      picard_extra=[SILENT])
    compare(f"interval list {which}, SAM", [f"I={foundation_sam}", f"INTERVAL_LIST={intervals}", SILENT])
    compare(f"interval list {which}, records", [f"I={indexed}", f"INTERVAL_LIST={intervals}",
                                                "RECORDS_ONLY=true", "ALIGNMENT_STATUS=Unaligned"],
            picard_extra=[SILENT])
    if which in ("targets", "mixed"):
        assert by_index != by_scan, "the index query and the scan are expected to differ here"


@pytest.mark.parametrize("source", ["foundation", "t1"])
def test_cram(corpus, tmp_path, source):
    bam, ref = {"foundation": (corpus / "corpus.bam", corpus / "ref.fa"),
                "t1": (CORPUS / "t1.bam", CORPUS / "reference.fa")}[source]
    cram = tmp_path / f"{source}.cram"
    subprocess.run(["samtools", "view", "-C", "--output-fmt-option", "version=3.0", "-T", str(ref),
                    "-o", str(cram), str(bam)], check=True)
    compare(f"{source} CRAM", [f"I={cram}", f"R={ref}"], picard_extra=[SILENT])
    both_fail(f"{source} CRAM without R", [f"I={cram}"], picard_extra=[SILENT])


@pytest.mark.parametrize("stringency", STRINGENCIES)
def test_handwritten_sam(handwritten, stringency):
    text = compare(f"handwritten {stringency}", [f"I={handwritten}", f"VALIDATION_STRINGENCY={stringency}"])
    assert text.count(b"\n") == HEADER.count("\n") + len(VALID)


def test_sam_of_the_corpora(corpus, foundation_sam, tmp_path):
    """SAM text of the foundation corpus and t1: the foundation corpus holds records
    htsjdk refuses at STRICT, so both fail there and agree at LENIENT and SILENT."""
    both_fail("foundation SAM STRICT", [f"I={foundation_sam}"])
    for stringency in ("LENIENT", "SILENT"):
        from_sam = compare(f"foundation SAM {stringency}",
                           [f"I={foundation_sam}", f"VALIDATION_STRINGENCY={stringency}"])
    from_bam = run_picard([f"I={corpus / 'corpus.bam'}", SILENT]).stdout
    assert from_sam == from_bam
    t1 = tmp_path / "t1.sam"
    t1.write_bytes(sam_text(CORPUS / "t1.bam"))
    compare("t1 SAM STRICT", [f"I={t1}"])


def test_compressed_sam_and_standard_input(handwritten, foundation_sam, tmp_path):
    expected = run_picard([f"I={handwritten}"]).stdout
    bgzipped = tmp_path / "valid.sam.gz"
    bgzipped.write_bytes(subprocess.run(["bgzip", "-c", str(handwritten)], check=True, capture_output=True).stdout)
    gzipped = tmp_path / "plain.sam.gz"
    gzipped.write_bytes(gzip.compress(handwritten.read_bytes()))
    for label, path in (("bgzip", bgzipped), ("gzip", gzipped)):
        assert compare(f"handwritten {label}", [f"I={path}"]) == expected
    assert compare("handwritten stdin", ["I=/dev/stdin"], stdin=lambda: open(handwritten, "rb")) == expected
    big = tmp_path / "corpus.sam.gz"
    big.write_bytes(subprocess.run(["bgzip", "-c", str(foundation_sam)], check=True, capture_output=True).stdout)
    compare("foundation bgzip SAM", [f"I={big}", SILENT])
    piped = run_ours(["-"], face=VERB, stdin=open(bgzipped, "rb"))
    assert piped.returncode == 0 and piped.stdout == expected, piped.stderr


@pytest.mark.parametrize("ending", ["\r\n", "\r"], ids=["crlf", "cr"])
def test_line_ends(tmp_path, ending):
    path = tmp_path / "ends.sam"
    path.write_bytes((HEADER + "".join(line + "\n" for line in VALID)).replace("\n", ending).encode())
    compare(f"line ends {ending!r}", [f"I={path}"])


@pytest.mark.parametrize("name,header", [
    ("no_hd", "@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:5000\n@RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:s\n"),
    ("hd_fields", "@HD\tSO:coordinate\tVN:1.0\tGO:none\tSS:coordinate:queryname\n@SQ\tSN:chr1\tLN:100000\n"
                  "@SQ\tSN:chr2\tLN:5000\n@RG\tID:rg1\tSM:s\n@RG\tID:rg2\tSM:s\n"),
    ("no_header", ""),
])
def test_headers(tmp_path, name, header):
    """Without `@HD` htsjdk writes `@HD VN:1.6`; `VN` comes first whatever the file says.
    Without any header the records name no reference htsjdk knows, so only the
    unmapped ones stay."""
    records = VALID if header else [line for line in VALID
                                    if line.split("\t")[2] == "*" and "RG:Z:" not in line]
    path = tmp_path / f"{name}.sam"
    path.write_text(header + "".join(line + "\n" for line in records))
    out = compare(f"header {name}", [f"I={path}"])
    assert out.startswith(b"@HD\tVN:")
    compare(f"header {name} only", [f"I={path}", "HEADER_ONLY=true"])


@pytest.mark.parametrize("name,line", [
    ("so", "@HD\tVN:1.6\tSO:sorted\n"),
    ("unknown_type", "@XX\tID:a\n"),
    ("rg_no_sm", "@RG\tID:rg3\n"),
    ("sq_no_ln", "@SQ\tSN:chr3\n"),
    ("key_value", "@PG\tID:p\tPN\n"),
    ("clash", "@PG\tID:p\tID:q\n"),
    ("long_key", "@PG\tID:p\tPNX:b\n"),
])
@pytest.mark.parametrize("stringency", STRINGENCIES)
def test_bad_header_lines(tmp_path, name, line, stringency):
    """A header line htsjdk refuses at STRICT and drops or repairs otherwise."""
    head = HEADER if name != "so" else HEADER.split("\n", 1)[1]
    path = tmp_path / f"{name}.sam"
    path.write_text(line + head + "".join(l + "\n" for l in VALID))
    args = [f"I={path}", f"VALIDATION_STRINGENCY={stringency}"]
    picard = run_picard(args)
    if picard.returncode != 0:
        both_fail(f"header {name} {stringency}", args)
    else:
        compare(f"header {name} {stringency}", args)


@pytest.mark.parametrize("name,line", BAD, ids=[name for name, _ in BAD])
def test_bad_lines(tmp_path, name, line):
    """One bad record line among good ones, at each stringency: both tools fail, or
    both print the same text."""
    path = tmp_path / f"{name}.sam"
    path.write_text(HEADER + VALID[0] + "\n" + line + "\n" + VALID[4] + "\n")
    outcomes = []
    for stringency in STRINGENCIES:
        args = [f"I={path}", f"VALIDATION_STRINGENCY={stringency}"]
        picard = run_picard(args)
        if picard.returncode != 0:
            both_fail(f"bad {name} {stringency}", args)
            outcomes.append("fails")
        else:
            compare(f"bad {name} {stringency}", args)
            outcomes.append("reads")
    assert outcomes[0] == "fails", f"{name}: htsjdk reads this line at STRICT"


@pytest.mark.parametrize("name,line,held", UNHOLDABLE, ids=[name for name, _, _ in UNHOLDABLE])
def test_lines_a_bam_cannot_hold(tmp_path, name, line, held):
    """Picard keeps these as text where it reads them (LENIENT, SILENT or both) and
    prints the line back. bammetrics holds BAM records: an unknown reference is read as
    none, an unknown base as N, a MAPQ above 255 by its low byte, and a QUAL of another
    length than SEQ or a CIGAR that does not parse is an error. Both fail at STRICT."""
    path = tmp_path / f"{name}.sam"
    path.write_text(HEADER + line + "\n")
    both_fail(f"unholdable {name} STRICT", [f"I={path}"])
    read_by_picard = 0
    for stringency in ("LENIENT", "SILENT"):
        args = [f"I={path}", f"VALIDATION_STRINGENCY={stringency}", "RECORDS_ONLY=true"]
        picard, ours = run_picard(args), run_ours(args)
        if picard.returncode != 0:
            assert ours.returncode != 0, f"{name} {stringency}"
            continue
        read_by_picard += 1
        assert picard.stdout == (line + "\n").encode()
        if held is None:
            assert ours.returncode != 0 and ours.stdout == b"", f"{name} {stringency}"
        else:
            assert ours.returncode == 0, ours.stderr
            assert ours.stdout == (held + "\n").encode(), f"{name} {stringency}"
    assert read_by_picard, f"{name}: Picard reads this line at no stringency"


def test_an_empty_file_is_an_empty_sam_for_picard_only(tmp_path):
    """htsjdk reads an empty file as SAM text with no header and no record; bammetrics
    names it as a file of no known kind."""
    path = tmp_path / "empty.sam"
    path.write_bytes(b"")
    picard, ours = run_picard([f"I={path}"]), run_ours([f"I={path}"])
    assert picard.returncode == 0 and picard.stdout == b"@HD\tVN:1.6\n"
    assert ours.returncode != 0 and b"not a BAM or CRAM file" in ours.stderr


def test_read_group_dates_stay_as_written(tmp_path):
    """htsjdk parses `@RG DT` and writes it back in the JVM's time zone; bammetrics
    keeps the text. Every other byte is the same."""
    path = tmp_path / "dated.sam"
    path.write_text(HEADER.replace("PU:fc.1", "PU:fc.1\tDT:2020-01-02T03:04:05+0000") + VALID[0] + "\n")
    picard, ours = run_picard([f"I={path}"]), run_ours([f"I={path}"])
    assert picard.returncode == 0 and ours.returncode == 0, ours.stderr
    without = lambda out: b"\n".join(  # noqa: E731
        b"\t".join(f for f in line.split(b"\t") if not (line.startswith(b"@RG") and f.startswith(b"DT:")))
        for line in out.split(b"\n"))
    assert without(picard.stdout) == without(ours.stdout)
    assert b"DT:2020-01-02T03:04:05+0000" in ours.stdout


def test_verb_forms(corpus, tmp_path):
    """`view-sam BAM` prints Picard's bytes, `-o FILE` writes them, `-` reads standard
    input, and `--format tsv` gives one row per record."""
    bam = corpus / "corpus.bam"
    expected = run_picard([f"I={bam}", SILENT]).stdout
    printed = run_ours([str(bam)], face=VERB)
    assert printed.returncode == 0 and printed.stdout == expected, printed.stderr
    out = tmp_path / "view.sam"
    to_file = run_ours([str(bam), "-o", str(out)], face=VERB)
    assert to_file.returncode == 0 and out.read_bytes() == expected, to_file.stderr
    with open(bam, "rb") as stdin:
        piped = run_ours(["-", "--pf-status", "PF"], face=VERB, stdin=stdin)
    assert piped.stdout == run_picard([f"I={bam}", "PF_STATUS=PF", SILENT]).stdout
    tsv = run_ours([str(bam), "--format", "tsv"], face=VERB)
    records = sum(1 for line in expected.split(b"\n") if line and not line.startswith(b"@"))
    assert tsv.returncode == 0 and tsv.stdout.count(b"\n") == records + 1, tsv.stderr
    assert tsv.stdout.startswith(b"QNAME\tFLAG\tRNAME\tPOS\tMAPQ\tCIGAR\tRNEXT\tPNEXT\tTLEN\tSEQ\tQUAL\tTAGS\n")
    refused = run_ours([str(bam), "-o", str(bam)], face=VERB)
    assert refused.returncode != 0 and bam.stat().st_size > 0
    lines = expected.count(b"\n")
    RESULTS.append((f"{VERB} verb forms", lines, lines, "stdout, -o and stdin identical"))


def chain_sam(source, out, reference=None):
    args = ["chain", str(source), "+", "write", str(out)]
    if reference:
        args[1:1] = ["-R", str(reference)]
    return subprocess.run([str(BAMMETRICS), *args], capture_output=True)


@pytest.mark.parametrize("source", ["markdup", "oxog", "foundation", "t1", "handwritten"])
def test_chain_writes_the_sam_file_picard_writes(corpus, handwritten, tmp_path, source):
    """`bammetrics chain X + write out.sam` against Picard `SamFormatConverter`: the
    whole file, byte for byte. The header is htsjdk's `SAMFileWriter` header (`VN:1.6`)."""
    path = {"markdup": FIXTURES / "markdup.bam", "oxog": FIXTURES / "oxog.bam",
            "foundation": corpus / "corpus.bam", "t1": CORPUS / "t1.bam", "handwritten": handwritten}[source]
    picard_out, ours_out = tmp_path / "picard.sam", tmp_path / "ours.sam"
    picard = run_picard([f"I={path}", f"O={picard_out}", SILENT], tool="SamFormatConverter")
    assert picard.returncode == 0, picard.stderr[-800:]
    ours = chain_sam(path, ours_out)
    assert ours.returncode == 0, ours.stderr
    a, b = picard_out.read_bytes(), ours_out.read_bytes()
    assert a == b, f"SamFormatConverter {source}: {first_difference(a, b)}"
    assert a.startswith(b"@HD\tVN:1.6")
    note_tags(a)
    lines = a.count(b"\n")
    RESULTS.append((f"chain write .sam {source}", lines, lines, "file identical to SamFormatConverter"))


@pytest.mark.parametrize("tool,args", [
    ("CollectInsertSizeMetrics", []),
    ("CollectAlignmentSummaryMetrics", []),
    ("CollectQualityYieldMetrics", []),
    ("MeanQualityByCycle", []),
    ("CollectDuplicateMetrics", []),
    ("EstimateLibraryComplexity", []),
    ("CollectWgsMetrics", ["R=REF"]),
    ("CollectGcBiasMetrics", ["R=REF", "S=OUT.summary"]),
    ("CollectHsMetrics", ["BI=BAITS", "TI=TARGETS"]),
    ("CheckDuplicateMarking", []),
])
def test_tools_read_sam_text_as_they_read_the_bam(corpus, foundation_sam, tmp_path, tool, args):
    """Each tool on the SAM text of the foundation corpus (SILENT, as its records need)
    writes the file it writes for the BAM."""
    bodies = []
    for side, source in (("sam", foundation_sam), ("bam", corpus / "corpus.bam")):
        out = tmp_path / f"{side}.txt"
        filled = [a.replace("REF", str(corpus / "ref.fa")).replace("OUT", str(tmp_path / side))
                  .replace("BAITS", str(corpus / "bait.interval_list"))
                  .replace("TARGETS", str(corpus / "target.interval_list")) for a in args]
        key = "M" if tool == "CollectDuplicateMetrics" else "O"
        run = subprocess.run([str(BAMMETRICS), tool, f"I={source}", f"{key}={out}", SILENT, *filled],
                             capture_output=True)
        assert run.returncode in (0, 1) and out.exists(), run.stderr
        bodies.append((run.returncode, [l for l in out.read_text().splitlines() if not l.startswith("#")]))
    assert bodies[0] == bodies[1], tool
    RESULTS.append((f"{tool} on SAM text", len(bodies[0][1]), len(bodies[0][1]), "same as on the BAM"))


def test_mark_duplicates_reads_sam_text(tmp_path):
    """MarkDuplicates on SAM text against Picard on the same text: the metrics and the
    records of the marked BAM."""
    sam = tmp_path / "markdup.sam"
    sam.write_bytes(sam_text(FIXTURES / "markdup.bam"))
    outputs = {}
    for side, run in (("picard", lambda a: run_picard(a, tool="MarkDuplicates")),
                      ("ours", lambda a: run_ours(a, face="MarkDuplicates"))):
        bam, metrics = tmp_path / f"{side}.bam", tmp_path / f"{side}.txt"
        done = run([f"I={sam}", f"O={bam}", f"M={metrics}"])
        assert done.returncode == 0, done.stderr[-800:]
        records = subprocess.run(["samtools", "view", str(bam)], check=True, capture_output=True).stdout
        body = [l for l in metrics.read_text().splitlines() if not l.startswith("#")]
        outputs[side] = (records, body)
    assert outputs["picard"] == outputs["ours"]
    RESULTS.append(("MarkDuplicates on SAM text", 36, 36, "records and metrics identical"))


def test_usage_errors_fail_in_both(corpus, tmp_path):
    bam = corpus / "corpus.bam"
    for label, args in [
        ("header and records only", [f"I={bam}", "HEADER_ONLY=true", "RECORDS_ONLY=true"]),
        ("OUTPUT", [f"I={bam}", f"O={tmp_path / 'x.sam'}"]),
        ("ASSUME_SORTED", [f"I={bam}", "ASSUME_SORTED=true"]),
        ("enum case", [f"I={bam}", "ALIGNMENT_STATUS=aligned"]),
        ("enum value", [f"I={bam}", "PF_STATUS=Maybe"]),
        ("stringency", [f"I={bam}", "VALIDATION_STRINGENCY=LOUD"]),
        ("missing input", [f"I={tmp_path / 'nope.bam'}"]),
        ("missing interval list", [f"I={bam}", f"INTERVAL_LIST={tmp_path / 'nope.interval_list'}"]),
        ("no input", ["HEADER_ONLY=true"]),
    ]:
        picard, ours = both_fail(label, args)
        assert ours.stdout == b"", label
    both = run_ours([f"I={bam}", "HEADER_ONLY=true", "RECORDS_ONLY=true"])
    assert b"Cannot specify both HEADER_ONLY=true and RECORDS_ONLY=true." in both.stderr
    text = tmp_path / "notes.txt"
    text.write_text("these are notes\nnot alignments\n")
    picard, ours = both_fail("a text file", [f"I={text}"])
    assert b"Not enough fields" in picard.stderr and b"not a BAM or CRAM file" in ours.stderr


def test_cases_are_covered():
    """Every tag type and every filter value was in compared output."""
    if not SEEN_TAGS:
        pytest.skip("no comparison ran")
    tags = {"A", "i", "f", "Z", "B:c", "B:C", "B:s", "B:S", "B:i", "B:I", "B:f"}
    assert tags <= SEEN_TAGS, tags - SEEN_TAGS
    options = {f"ALIGNMENT_STATUS={v}" for v in ("Aligned", "Unaligned", "All")}
    options |= {f"PF_STATUS={v}" for v in ("PF", "NonPF", "All")}
    options |= {"HEADER_ONLY=true", "RECORDS_ONLY=true"}
    options |= {f"VALIDATION_STRINGENCY={v}" for v in STRINGENCIES}
    assert options <= SEEN_OPTIONS, options - SEEN_OPTIONS

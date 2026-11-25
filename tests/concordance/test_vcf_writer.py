"""The VCF writer of wave W16 against Picard 3.5.0 (htsjdk 5.0.0).

The oracle is Picard VcfFormatConverter, a pass-through writer: it reads each record
with VCFCodec and writes it with VCFWriter, as `bammetrics chain X -o OUT` does.

- Round trips: every VCF of variant_vcfs.py, plain and .vcf.gz, and the corner files of
  CORNERS, written by both tools. compare_vcf's rule: every header line equal except the
  date and command lines and the lines the inventory (vcf_headers.json) marks volatile,
  then every record line equal. A .vcf.gz is BGZF with its end-of-file block, and its
  .tbi answers every tabix query as Picard's does (index bytes are not compared).
- The inventory: each case of vcf_headers.json whose output header is a pure rewrite
  (no line added or removed, the sample columns unchanged) is written by the chain from
  its input header, and the header lines must equal the recorded output header.
- Error parity: what VCFCodec or VCFEncoder refuses fails both tools, with no file left.
- Taps: CollectVariantCallingMetrics and GenotypeConcordance on a variant stream write
  the bodies the tools write alone (test_variant_tools.py holds those to Picard).
"""
import json
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS
from compare import BGZF_EOF, diff_vcf, vcf_header_lines, vcf_lines
from conftest import RESULTS
from run_picard import picard_env
from variant_vcfs import build_all

HERE = Path(__file__).resolve().parent
INVENTORY = json.loads((HERE / "vcf_headers.json").read_text())


def volatile_patterns():
    """The regexes of the header lines the inventory marks volatile in any case
    (a `##reference=file:` path), each once."""
    patterns = []
    for entry in INVENTORY:
        for case in [entry, *entry.get("variants", [])]:
            for line in case.get("volatile", []) or []:
                if line["regex"] not in patterns:
                    patterns.append(line["regex"])
    return [re.compile(p) for p in patterns]


VOLATILE = volatile_patterns()


def header_rule(lines):
    """The header lines compared: compare_vcf's (date and command lines aside),
    then without the lines an inventoried regex matches."""
    return [line for line in vcf_header_lines(lines)
            if not any(p.match(line) for p in VOLATILE)]


def diff_vcf_rule(picard, ours):
    """compare_vcf's differences with the header rule of the inventory."""
    problems = [p for p in diff_vcf(picard, ours) if not p.startswith("header line")]
    a, b = header_rule(vcf_lines(picard)), header_rule(vcf_lines(ours))
    if a != b:
        problems.append(f"header lines differ: {[x for x in a if x not in b]} / "
                        f"{[x for x in b if x not in a]}")
    return problems


HEADER = """##fileformat=VCFv4.2
##source=a caller
##contig=<ID=c2,length=3000,assembly=b37,md5=abc>
##contig=<ID=c1,length=5000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=FL,Number=0,Type=Flag,Description="A flag">
##INFO=<ID=FY,Number=.,Type=Flag,Description="A flag of any count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele freq",Source="gatk",Version="3">
##INFO=<ID=ZZ,Number=.,Type=String,Description="Has \\"quotes\\" and, commas">
##INFO=<ID=XT,Number=1,Type=String,Description="Extra attrs",Foo=bar>
##INFO=<ID=NODESC,Number=1,Type=Integer>
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Repaired to Float">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">
##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">
##FORMAT=<ID=XA,Number=A,Type=Float,Description="Per alternate">
##FORMAT=<ID=XS,Number=1,Type=String,Description="Text">
##ALT=<ID=NON_REF,Description="Any other allele">
##META=<ID=Assay,Type=String,Number=.,Values=[WholeGenome, Exome]>
##SAMPLE=<ID=S1,Assay=WholeGenome>
##noequals
##key=value with spaces
##key=value with spaces
"""

ROWS = [
    ("c2", 5, ".", "N", "A", "30.5", "PASS", "AF=1e-5", "GT", "0/1", "0/1"),
    ("c2", 6, ".", "C", "T,G", "1.005", "PASS", "AF=0.1,0.000001", "GT:AD:XA", "1/2:1,2,3:0.5", "0/1:.:."),
    ("c2", 7, ".", "C", "T", "2.675", "PASS", "FL;FL;FY=2", "GT", "0/1", "0/1"),
    ("c1", 10, ".", "a", "c", "30", "PASS", "DP=10;FL;AF=0.5", "GT:GQ:DP:AD:PL",
     "0/1:30.5:10:3,7:10,0,100", "1|1:.:.:.:."),
    ("c1", 11, "rs1;rs2", "AC", "A,ACT", "30.0", "q10;LowQual", "AF=0.25,0.5;DP=3", "GT:GQ:AD",
     "0/1:.:1,2,3", "./.:.:."),
    ("c1", 12, "rs3", "G", "T", "12.345", "LowQual;LowQual", ".", "GT:FT", "0/1:PASS", "0/1:q10;b;a"),
    ("c1", 13, ".", "G", ".", ".", ".", "FL=1;XT=ab;FL=0", "GT:FT", "0/0:PASS", "0/0:."),
    ("c1", 14, ".", "G", "<NON_REF>", "1e2", "PASS", "DP=7;NODESC=3;MQ=60", "GT:PL",
     "0/0:0,10,100", "0/0:0,10,100"),
    ("c1", 15, ".", "GTTT", "G", "0.001", "q10", "ZZ=x,y;AF=1.0", "GT:GQ:XS", "1/1:-1:x", "0|1:99.5"),
    ("c1", 16, ".", "C", "T", "1000000.5", "PASS", "DP=.", "GT:DP:AD", "0|1:-1:1,x", "1:5:2,3"),
    ("c1", 17, ".", "C", "T", "-1", "PASS", "DP=;DP", "GT:XS:XA", "0/1", "."),
    ("c1", 18, ".", "C", "T", "45.678", "PASS", "DP=10;DP=11", "XS:GQ", "a:20", "b"),
    ("c1", 19, ".", "A", "C,T,G", "-10", ".", ".", "GT:PL:XA", "0/3:.:1,2,3", "./.:.:."),
    ("c1", 20, ".", "A", "C", "inf", ".", ".", "GT", "0//1", "1/"),
]

COLUMNS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def vcf_text(version="VCFv4.2", samples=("S1", "S2"), rows=ROWS, header=HEADER):
    head = header.replace("##fileformat=VCFv4.2", f"##fileformat={version}")
    columns = COLUMNS + ("\tFORMAT\t" + "\t".join(samples) if samples else "")
    body = "".join("\t".join(map(str, row[:len(columns.split("\t"))])) + "\n" for row in rows)
    return head + columns + "\n" + body


def sorted_rows():
    order = {"c2": 0, "c1": 1}
    return sorted(ROWS, key=lambda r: (order[r[0]], r[1]))


# Corner files: name -> text. Their records are sorted by the contig lines' order.
CORNERS = {
    "corners_sorted_samples": vcf_text(rows=sorted_rows()),
    "corners_unsorted_samples": vcf_text(samples=("S2", "S1"), rows=sorted_rows()),
    "corners_v41": vcf_text("VCFv4.1", rows=sorted_rows()),
    "corners_v43": vcf_text("VCFv4.3", samples=("S2", "S1"), rows=sorted_rows()),
    "corners_sites_only": vcf_text(samples=(), rows=sorted_rows()),
    "corners_crlf": vcf_text(rows=sorted_rows()).replace("\n", "\r\n"),
    "corners_percent": vcf_text("VCFv4.3", rows=[
        ("c1", 5, ".", "A", "C", "30", ".", "XT=a%3Bb;ZZ=%2C,%41", "GT:XS", "0/1:p%2Cq", "1/1:%"),
    ]),
}

# The corners of the fixup of wave W16, each found against Picard: the header adds a
# key of Number=2 and GL. Samples out of sorted order, so htsjdk decodes the genotypes.
HEADER_W16 = HEADER.replace(
    "##FORMAT=<ID=XS,",
    '##FORMAT=<ID=X2,Number=2,Type=String,Description="Two">\n'
    '##FORMAT=<ID=XR,Number=R,Type=String,Description="Per allele">\n'
    '##FORMAT=<ID=GL,Number=G,Type=Float,Description="Likelihoods">\n'
    "##FORMAT=<ID=XS,")


def w16(rows, version="VCFv4.2", samples=("S2", "S1")):
    return vcf_text(version, samples=samples, rows=rows, header=HEADER_W16)


CORNERS.update({
    # GL becomes PL; a PL read earlier in the record wins over a later GL.
    "corners_gl": w16([
        ("c1", 1, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:-1.5,-0.1,-100", "1/1:0,-1,-2"),
        ("c1", 2, ".", "A", "C", ".", ".", ".", "GT:GL:PL", "0/1:-1.5,-0.1,-100:5,6,7", "1/1:0,-1,-2"),
        ("c1", 3, ".", "A", "C", ".", ".", ".", "GT:PL:GL", "0/1:5,6,7:-1,0,-2", "1/1:0,-1,-2"),
        ("c1", 4, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:.,.,.", "1/1:."),
        ("c1", 5, ".", "A", "C,G", ".", ".", ".", "GT:GL", "0/2:-0.04,-1.4,-2.5,-3,-4.44,-5.55",
         "1/1:-1e-5,-1,-2,-3,-4,-500"),
        ("c1", 6, ".", "A", "C", ".", ".", ".", "GT:GL", "0:-0.04,-1.4", "1:0,-2.05"),
        ("c1", 7, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:NaN,0,-1", "1/1:-Infinity,0,-1"),
        ("c1", 8, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:1.5,0.2,-3", "1/1:-0.05,-0.15,-0.25"),
    ]),
    # A missing value before a present one is one `.` whatever the key's Number.
    "corners_missing_values": w16([
        ("c1", 1, ".", "A", "T,G", ".", ".", ".", "GT:XA:DP", "0/1:.:5", "0/1:0.5,.:6"),
        ("c1", 2, ".", "A", "C,G,T", ".", ".", ".", "GT:XS:XA", "0/1:b", "0/1:c:1,2,3"),
        ("c1", 3, ".", "A", "C", ".", ".", ".", "GT:X2:XS:XR", "0/1:.:b:.", "0/1:.:.:z"),
        ("c1", 4, ".", "A", "C", ".", ".", ".", "GT:XR:X2:XS", "0/1:.:.:b", "0/1"),
        ("c1", 5, ".", "A", "C", ".", ".", ".", "GT:PL:XS", "0/1:.:b", "0/1:1,0,2"),
        ("c1", 6, ".", "A", "C", ".", ".", ".", "GT:AD:DP", "0/1:.:5", "0/1:1,2:6"),
    ]),
    # A GQ that rounds to -1 is no GQ; a QUAL of -0 prints 0.
    "corners_gq_qual": w16([
        ("c1", 1, ".", "A", "C", "-0", ".", ".", "GT:GQ", "0/1:-1.5", "0/1:-0.6"),
        ("c1", 2, ".", "A", "C", "-0.0", ".", ".", "GT:GQ", "0/1:-1.0", "0/1:-0.4"),
    ]),
    # An empty first sample cell (DP -1 and FF missing), and v4.3 text with %3A: both
    # written, and the .tbi built from the written lines.
    "corners_empty_cell": w16([
        ("c1", 1, ".", "A", "C", ".", ".", ".", "DP:XS", "-1:.", "5:x"),
    ]),
    "corners_v43_colon": w16([
        ("c1", 1, ".", "A", "C", ".", ".", ".", "GT:XS", "0/1:a%3Ab", "0/1:x"),
    ], version="VCFv4.3"),
    # Text is written as ISO-8859-1: up to U+00FF one byte, any other character `?`,
    # a byte that is not UTF-8 `?` too.
    "corners_latin1": (
        HEADER_W16.replace('Description="Text"', 'Description="DépĊth"')
        + "##comment=café 中 \U0001F600  \n"
        + COLUMNS + "\tFORMAT\tS2\tt中\n"
        + "c1\t1\tré\tA\tC\t.\t.\tXT=café中\U0001F600\tGT:XS\t0/1:é\t0/1:x\n"
    ).encode() .replace(b"XT=caf\xc3\xa9", b"XT=caf\xc3\xa9\xe9", 1),
    # A bare CR ends a line, as htsjdk's line reader ends it.
    "corners_cr_only": w16([
        ("c1", 1, ".", "A", "C", ".", ".", ".", "GT", "0/1", "0/1"),
        ("c1", 2, ".", "A", "C", ".", ".", ".", "GT", "0/1", "1/1"),
    ]).replace("\n", "\r"),
    # The version line need not be the first header line.
    "corners_source_first": "##source=first\n" + w16([
        ("c1", 1, ".", "A", "C", ".", ".", ".", "GT", "0/1", "0/1"),
    ]),
})

# Files both tools refuse: name -> text.
ERRORS = {
    "undeclared_info": vcf_text(rows=[("c1", 1, ".", "A", "C", ".", ".", "NOPE=1", "GT", "0", "0")]),
    "undeclared_filter": vcf_text(rows=[("c1", 1, ".", "A", "C", ".", "nofilter", ".", "GT", "0", "0")]),
    "whitespace_info": vcf_text(rows=[("c1", 1, ".", "A", "C", ".", ".", "XT=a b", "GT", "0", "0")]),
    "negative_quality": vcf_text(rows=[("c1", 1, ".", "A", "C", "-5", ".", ".", "GT", "0", "0")]),
    "too_many_values": vcf_text(samples=("S2", "S1"),
                                rows=[("c1", 1, ".", "A", "C", ".", ".", ".", "GT", "0/1", "0/1:5")]),
    "undeclared_format": vcf_text(samples=("S2", "S1"),
                                  rows=[("c1", 1, ".", "A", "C", ".", ".", ".", "GT:QQ", "0:1", "0:1")]),
    "bad_gt_index": vcf_text(samples=("S2", "S1"),
                             rows=[("c1", 1, ".", "A", "C", ".", ".", ".", "GT", "0/3", "0")]),
    "duplicate_allele": vcf_text(rows=[("c1", 1, ".", "A", "C,C", ".", ".", ".", "GT", "0", "0")]),
    "no_version": HEADER.replace("##fileformat=VCFv4.2\n", "") + COLUMNS + "\n",
    "gl_not_a_number": w16([("c1", 1, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:x,1,2", "0")]),
    "gl_partly_missing": w16([("c1", 1, ".", "A", "C", ".", ".", ".", "GT:GL", "0/1:.,1,2", "0")]),
    # A CR inside a field ends the line there: the record is cut short.
    "cr_in_id": w16([("c1", 1, "a\rb", "A", "C", ".", ".", ".", "GT", "0/1", "0/1")]),
    "cr_in_info": w16([("c1", 1, ".", "A", "C", ".", ".", "XT=a\rb", "GT", "0/1", "0/1")]),
}

# Inputs whose .vcf.gz index both tools refuse: a record past the 32,770 16 kb windows a
# .tbi holds, by its start on a long contig or by its END.
INDEX_LIMITS = {
    "past_tbi_start": HEADER_W16.replace("##contig=<ID=c1,length=5000>",
                                         "##contig=<ID=c1,length=600000000>")
    + COLUMNS + "\n" + "c1\t10\t.\tA\tC\t.\t.\t.\n" + "c1\t550000000\t.\tA\tC\t.\t.\t.\n",
    "past_tbi_end": HEADER_W16.replace(
        "##INFO=<ID=DP,", '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n##INFO=<ID=DP,')
    + COLUMNS + "\n" + "c1\t1\t.\tA\tC\t.\t.\tEND=2000000000\n",
}

# Inputs a .vcf.gz with an index refuses: records out of order (c1 12 before c1 11).
UNSORTED = vcf_text(rows=[ROWS[3], ROWS[5], ROWS[4]])


@pytest.fixture(scope="module")
def vcfs(tmp_path_factory):
    out = tmp_path_factory.mktemp("vcf_writer")
    files = build_all(out)
    for name, text in {**CORNERS, **ERRORS, **INDEX_LIMITS, "unsorted": UNSORTED}.items():
        path = out / f"{name}.vcf"
        path.write_bytes(text if isinstance(text, bytes) else text.encode())
        files[f"{name}.vcf"] = path
    files["dir"] = out
    return files


def run_picard(src, out, *args):
    return subprocess.run(["picard", "VcfFormatConverter", f"I={src}", f"O={out}",
                           "REQUIRE_INDEX=false", "VERBOSITY=ERROR", "QUIET=true", *args],
                          env=picard_env(), capture_output=True, text=True)


def run_chain(src, out, *args):
    return subprocess.run([str(BAMMETRICS), "chain", str(src), "-o", str(out), "-q", *args],
                          capture_output=True, text=True)


def tabix_regions(path, contigs):
    """Every record line `tabix` returns for whole contigs and 200-base windows."""
    regions = []
    for name, length in contigs:
        regions.append(name)
        regions.extend(f"{name}:{start}-{start + 199}" for start in range(1, length + 1, 200))
    out = []
    for region in regions:
        # Bytes: a record may hold ISO-8859-1 text, which is no UTF-8.
        done = subprocess.run(["tabix", str(path), region], capture_output=True)
        out.append((region, done.returncode, done.stdout))
    return out


def contigs_of(lines):
    found = []
    for line in lines:
        match = re.match(r"##contig=<ID=([^,>]+)(?:,length=(\d+))?", line)
        if match:
            found.append((match.group(1), int(match.group(2) or 5000)))
    return found


ROUND_TRIPS = ["truth.vcf", "call.vcf", "calls3.vcf", "dbsnp.vcf", "gvcf.vcf", "sites_only.vcf",
               "nocontigs.vcf", "calls3_nocontigs.vcf", "mismatch_truth.vcf", "mismatch_call.vcf",
               "haploid_call.vcf", *[f"{name}.vcf" for name in CORNERS]]


@pytest.mark.parametrize("name", ROUND_TRIPS)
@pytest.mark.parametrize("suffix", [".vcf", ".vcf.gz"])
def test_round_trip(vcfs, tmp_path, name, suffix):
    src = vcfs.get(name) or vcfs["dir"] / name
    if not Path(src).exists():
        src = vcfs["dir"] / (name + ".gz")
    picard, ours = tmp_path / f"picard{suffix}", tmp_path / f"ours{suffix}"
    indexed = suffix == ".vcf.gz"
    p = run_picard(src, picard, f"CREATE_INDEX={'true' if indexed else 'false'}")
    o = run_chain(src, ours)
    if p.returncode != 0:
        assert o.returncode != 0, f"{name}: Picard failed, bammetrics did not:\n{p.stderr[-800:]}"
        assert not ours.exists(), f"{name}: bammetrics failed and left {ours.name}"
        RESULTS.append((f"VcfWriter {name}{suffix} both fail", 1, 1, "both fail"))
        return
    assert o.returncode == 0, f"{name}: {o.stderr}"
    problems = diff_vcf_rule(picard, ours)
    assert not problems, f"{name}{suffix}:\n  " + "\n  ".join(problems[:20])
    lines = vcf_lines(ours)
    if indexed:
        assert ours.read_bytes().endswith(BGZF_EOF)
        ours_index = Path(f"{ours}.tbi")
        assert ours_index.exists() and Path(f"{picard}.tbi").exists()
        if shutil.which("tabix") and contigs_of(lines):
            assert tabix_regions(picard, contigs_of(lines)) == tabix_regions(ours, contigs_of(lines))
    else:
        assert not Path(f"{ours}.idx").exists()
    RESULTS.append((f"VcfWriter {name}{suffix}", len(lines), len(lines), "every line equal"))


@pytest.mark.parametrize("name", sorted(ERRORS))
def test_error_parity(vcfs, tmp_path, name):
    src = vcfs[f"{name}.vcf"]
    picard, ours = tmp_path / "picard.vcf", tmp_path / "ours.vcf"
    p = run_picard(src, picard, "CREATE_INDEX=false")
    o = run_chain(src, ours)
    assert p.returncode != 0, f"{name}: Picard succeeded"
    assert o.returncode != 0, f"{name}: bammetrics succeeded"
    assert not ours.exists(), f"{name}: bammetrics left {ours.name}"
    RESULTS.append((f"VcfWriter error {name}", 1, 1, "both fail"))


def test_unsorted_records_fail_the_index(vcfs, tmp_path):
    src = vcfs["unsorted.vcf"]
    picard, ours = tmp_path / "picard.vcf.gz", tmp_path / "ours.vcf.gz"
    p = run_picard(src, picard, "CREATE_INDEX=true")
    o = run_chain(src, ours)
    assert p.returncode != 0 and o.returncode != 0, (p.stderr[-500:], o.stderr)
    assert not ours.exists() and not Path(f"{ours}.tbi").exists()
    # Without the index both write every record in input order.
    p = run_picard(src, tmp_path / "p.vcf.gz", "CREATE_INDEX=false")
    o = subprocess.run([str(BAMMETRICS), "chain", str(src), "-q", "write", str(tmp_path / "o.vcf.gz")],
                       capture_output=True, text=True)
    assert p.returncode == 0 and o.returncode == 0, (p.stderr[-500:], o.stderr)
    assert not diff_vcf_rule(tmp_path / "p.vcf.gz", tmp_path / "o.vcf.gz")
    RESULTS.append(("VcfWriter unsorted", 2, 2, "index fails in both, copy equal"))


@pytest.mark.parametrize("name", sorted(INDEX_LIMITS))
def test_index_limits_fail_in_both(vcfs, tmp_path, name):
    src = vcfs[f"{name}.vcf"]
    picard, ours = tmp_path / "picard.vcf.gz", tmp_path / "ours.vcf.gz"
    p = run_picard(src, picard, "CREATE_INDEX=true")
    o = run_chain(src, ours)
    assert p.returncode != 0 and o.returncode != 0, (p.stderr[-500:], o.stderr)
    assert "ArrayIndexOutOfBounds" in p.stderr, p.stderr[-500:]
    assert ".tbi" in o.stderr, o.stderr
    assert not ours.exists() and not Path(f"{ours}.tbi").exists()
    # The plain copy holds every record in both.
    p = run_picard(src, tmp_path / "p.vcf", "CREATE_INDEX=false")
    o = run_chain(src, tmp_path / "o.vcf")
    assert p.returncode == 0 and o.returncode == 0, (p.stderr[-500:], o.stderr)
    assert not diff_vcf_rule(tmp_path / "p.vcf", tmp_path / "o.vcf")
    RESULTS.append((f"VcfWriter {name}", 2, 2, "index fails in both, copy equal"))


@pytest.mark.parametrize("suffix", [".gz", ".bgz", ".GZ"])
def test_compressed_names_are_bgzf(vcfs, tmp_path, suffix):
    """htsjdk writes BGZF and a .tbi for any name ending in .gz or .bgz."""
    src = vcfs["corners_unsorted_samples.vcf"]
    picard, ours = tmp_path / f"picard{suffix}", tmp_path / f"ours{suffix}"
    p = run_picard(src, picard, "CREATE_INDEX=true")
    o = run_chain(src, ours)
    assert p.returncode == 0 and o.returncode == 0, (p.stderr[-500:], o.stderr)
    assert ours.read_bytes()[:4] == picard.read_bytes()[:4] == b"\x1f\x8b\x08\x04"
    assert ours.read_bytes().endswith(BGZF_EOF)
    assert Path(f"{ours}.tbi").exists() and Path(f"{picard}.tbi").exists()
    assert not diff_vcf_rule(picard, ours)
    if shutil.which("tabix"):
        contigs = contigs_of(vcf_lines(ours))
        assert tabix_regions(picard, contigs) == tabix_regions(ours, contigs)
    RESULTS.append((f"VcfWriter name {suffix}", 1, 1, "BGZF and .tbi as Picard"))


def pure_rewrites():
    """(label, input header, output header) of every inventory case whose
    output header adds and removes no line and keeps the sample columns."""
    cases = []
    for entry in INVENTORY:
        if not entry.get("run"):
            continue
        for case in [entry, *entry.get("variants", [])]:
            if case.get("error") or case.get("added") or case.get("removed"):
                continue
            if "input_header" not in case or "output_header" not in case:
                continue
            rewritten = case.get("rewritten") or []
            if any(r["key"] != ["FORMAT", "FT"] for r in rewritten):
                continue
            label = f"{entry['tool']}:{case.get('label', 'default')}"
            cases.append((label, case["input_header"], case["output_header"]))
    return cases


@pytest.mark.parametrize("label,input_header,output_header", pure_rewrites(),
                         ids=[c[0] for c in pure_rewrites()])
def test_headers_match_the_inventory(tmp_path, label, input_header, output_header):
    src = tmp_path / "in.vcf"
    src.write_text("\n".join(input_header) + "\n")
    out = tmp_path / "out.vcf"
    done = run_chain(src, out)
    assert done.returncode == 0, done.stderr
    ours = header_rule(vcf_lines(out))
    recorded = header_rule(output_header)
    assert ours == recorded, f"{label}: {ours} != {recorded}"
    RESULTS.append((f"VcfWriter inventory {label}", len(ours), len(ours), "header lines equal"))


def bodies(paths):
    """Each file without its header comment lines."""
    return [[line for line in Path(p).read_text().splitlines() if not line.startswith("# ")]
            for p in paths]


def test_variant_calling_tap(vcfs, tmp_path):
    for label, args in [("default", []), ("targets", ["--ti", str(vcfs["dir"] / "targets.interval_list")])]:
        alone = tmp_path / f"alone_{label}"
        chained = tmp_path / f"chained_{label}"
        subprocess.run([str(BAMMETRICS), "variant-calling", str(vcfs["dir"] / "calls3.vcf.gz"),
                        "--dbsnp", str(vcfs["dir"] / "dbsnp.vcf.gz"), "-o", str(alone), *args],
                       check=True, capture_output=True)
        subprocess.run([str(BAMMETRICS), "chain", str(vcfs["dir"] / "calls3.vcf.gz"), "-q",
                        "-o", str(chained), "variant-calling", "--dbsnp",
                        str(vcfs["dir"] / "dbsnp.vcf.gz"), *args, "write", str(tmp_path / f"{label}.vcf.gz")],
                       check=True, capture_output=True)
        for ext in ("variant_calling_detail_metrics", "variant_calling_summary_metrics"):
            a, c = bodies([f"{alone}.{ext}", f"{chained}.{ext}"])
            assert a == c, f"{label} {ext}"
    RESULTS.append(("VcfWriter tap variant-calling", 2, 2, "bodies equal to the tool"))


def test_genotype_concordance_tap(vcfs, tmp_path):
    d = vcfs["dir"]
    for label, args in [("C", ["--call-sample", "C"]), ("C2", ["--call-sample", "C2", "--min-gq", "30"])]:
        alone, chained = tmp_path / f"alone_{label}", tmp_path / f"chained_{label}"
        subprocess.run([str(BAMMETRICS), "genotype-concordance", str(d / "truth.vcf.gz"),
                        str(d / "call.vcf.gz"), "--truth-sample", "T", "-o", str(alone), *args],
                       check=True, capture_output=True)
        subprocess.run([str(BAMMETRICS), "chain", str(d / "call.vcf.gz"), "-q", "-o", str(chained),
                        "genotype-concordance", "--truth-vcf", str(d / "truth.vcf.gz"),
                        "--truth-sample", "T", *args],
                       check=True, capture_output=True)
        for ext in ("summary", "detail", "contingency"):
            name = f"genotype_concordance_{ext}_metrics"
            a, c = bodies([f"{alone}.{name}", f"{chained}.{name}"])
            assert a == c, f"{label} {name}"
    RESULTS.append(("VcfWriter tap genotype-concordance", 2, 2, "bodies equal to the tool"))

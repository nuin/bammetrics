"""CollectExonCoverageMetrics against CollectHsMetrics PER_TARGET_COVERAGE (oracle tier A,
docs/exon_coverage.md): for the same intervals, each exon's MEAN_DEPTH and MINIMUM_DEPTH
print as the mean_coverage and min_coverage of Picard 3.5.0 and of bammetrics, to the
last printed digit.

Each input gets a GTF written here from its target list: genes of up to three
targets, every other gene on the minus strand, a second transcript whose exons are
the first's shifted by a third of their length, a one-base exon and a gene on a
contig the BAM lacks. Picard merges overlapping and abutting targets, so the exons
are split into layers in which no two exons touch, and CollectHsMetrics runs once
per layer with the layer as baits and targets: every exon is then one Picard target.

The cross-check is samtools bedcov (-Q at the MAPQ filter, -j without deletions):
its sum counts every base of a read the flags keep, so a mean from it is at least
MEAN_DEPTH; the gap is the base quality filter and the mate overlap counted once.
With MQ=0 Q=0 and CLIP_OVERLAPPING_READS=false only bases at quality 2 or less and
MAPQ 0 adapter reads differ, within BEDCOV_TOLERANCE.

The last test is the distinct-value gate: MEAN_DEPTH and MINIMUM_DEPTH took at least
3 distinct values over the exons compared with Picard.
"""
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS
from conftest import RESULTS
from hs_bams import build_all as build_hs
from picard_metrics import parse_picard_metrics
from run_picard import picard_env

TOOL = "CollectExonCoverageMetrics"
SILENT = "VALIDATION_STRINGENCY=SILENT"
# The Hs sensitivity sampler is not compared; a small sample keeps the runs short.
SMALL = "SAMPLE_SIZE=100"
SEEN = {"MEAN_DEPTH": set(), "MINIMUM_DEPTH": set()}
MISSING_CONTIG = "chrNotInTheBam"
# Mean depth from samtools bedcov over MEAN_DEPTH at MQ=0 Q=0 without clipping: at
# most this relative excess per exon (measured 0.55% on the exome slice, 0.48% on the
# foundation corpus, 0 on hs_deep; docs/exon_coverage.md).
BEDCOV_TOLERANCE = 0.01


def header(bam):
    return subprocess.run(["samtools", "view", "-H", str(bam)], check=True, capture_output=True,
                          text=True).stdout


def contigs(bam):
    """{name: length} from the @SQ lines."""
    out = {}
    for line in header(bam).splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            out[fields["SN"]] = int(fields["LN"])
    return out


def read_targets(path):
    """(contig, start, end) of an interval list or a BED, 1-based closed."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith(("@", "#", "track", "browser")):
            continue
        f = line.split("\t")
        if path.suffix == ".bed":
            rows.append((f[0], int(f[1]) + 1, int(f[2])))
        else:
            rows.append((f[0], int(f[1]), int(f[2])))
    return sorted(set(rows))


def annotation(targets, lengths):
    """GTF lines and the exons as (contig, start, end, transcript, strand)."""
    merged = []
    for contig, start, end in targets:
        if merged and merged[-1][0] == contig and start <= merged[-1][2] + 1:
            merged[-1] = (contig, merged[-1][1], max(end, merged[-1][2]))
        else:
            merged.append((contig, start, end))
    exons = []
    groups = [merged[i:i + 3] for i in range(0, len(merged), 3)]
    for g, group in enumerate(groups):
        strand = "-" if g % 2 else "+"
        previous = None
        for t, (contig, start, end) in enumerate(group):
            exons.append((contig, start, end, f"T{g}.1", strand, f"G{g}", t + 1))
            shift = max(1, (end - start + 1) // 3)
            shifted = (contig, min(start + shift, lengths[contig]), min(end + shift, lengths[contig]))
            # The tool merges overlapping exons of one transcript: keep them apart.
            if g % 2 == 0 and (previous is None or previous[0] != contig or shifted[1] > previous[2]):
                exons.append((*shifted, f"T{g}.2", strand, f"G{g}", t + 1))
                previous = shifted
        if g == 0:
            contig, start, end = group[0]
            middle = (start + end) // 2
            exons.append((contig, middle, middle, "T0.3", strand, "G0", 1))
    exons.append((MISSING_CONTIG, 100, 200, "TX.1", "+", "GX", 1))
    lines = [
        f"{c}\ttest\texon\t{s}\t{e}\t.\t{st}\t.\tgene_id \"{g}\"; gene_name \"N{g}\"; "
        f"transcript_id \"{t}\"; exon_number \"{n}\";\n"
        for c, s, e, t, st, g, n in exons
    ]
    return "".join(lines), exons


def layers(exons):
    """The distinct placed exon intervals in groups where no two touch on one contig."""
    distinct = sorted({(c, s, e) for c, s, e, *_ in exons if c != MISSING_CONTIG})
    out = []
    for interval in distinct:
        for layer in out:
            last = layer[-1]
            if last[0] != interval[0] or interval[1] > last[2] + 1:
                layer.append(interval)
                break
        else:
            out.append([interval])
    return out


def interval_list(bam, intervals, path):
    lines = [line for line in header(bam).splitlines() if line.startswith(("@HD", "@SQ"))]
    lines += [f"{c}\t{s}\t{e}\t+\tx{c}_{s}_{e}" for c, s, e in intervals]
    path.write_text("\n".join(lines) + "\n")
    return path


def per_target(path):
    """{(contig, start, end): (mean_coverage, min_coverage)} as printed."""
    lines = Path(path).read_text().splitlines()
    columns = lines[0].split("\t")
    out = {}
    for line in lines[1:]:
        row = dict(zip(columns, line.split("\t")))
        out[(row["chrom"], int(row["start"]), int(row["end"]))] = (row["mean_coverage"], row["min_coverage"])
    return out


def hs_per_target(command, bam, ref, layer_list, args, out):
    pt = out.with_suffix(".per_target")
    cmd = [command, "CollectHsMetrics", f"I={bam}", f"O={out}", f"BI={layer_list}", f"TI={layer_list}",
           f"R={ref}", f"PER_TARGET_COVERAGE={pt}", SILENT, SMALL, *args]
    if command == "picard":
        cmd += ["VERBOSITY=ERROR", "QUIET=true"]
        subprocess.run(cmd, check=True, env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    else:
        subprocess.run(cmd, check=True, capture_output=True)
    return per_target(pt)


def ours(bam, gtf, args, out):
    run = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"ANNOTATION={gtf}", f"O={out}", SILENT, *args],
                         check=True, capture_output=True, text=True)
    assert MISSING_CONTIG in run.stderr, "the contig the BAM lacks is reported"
    return parse_picard_metrics(out).rows


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("hs")
    bams = build_hs(out)
    return {name: (bam, out / "ref.fa", out / ("deep_panel.bed" if name == "deep" else "hs_panel.bed"))
            for name, bam in bams.items()}


def inputs(name, built, corpus, exome):
    """(bam, reference, target list)."""
    if name == "foundation":
        return corpus / "corpus.bam", corpus / "ref.fa", corpus / "target.interval_list"
    if name == "exome":
        return exome / "exome.bam", exome / "ref.fa", exome / "targets.interval_list"
    if name.startswith("hs_"):
        return built[name[3:]]
    return CORPUS / f"{name}.bam", CORPUS / "reference.fa", CORPUS / "targets.interval_list"


CASES = [("foundation", []), ("exome", []), ("exome", ["MQ=0", "Q=0"]), ("exome", ["CLIP_OVERLAPPING_READS=false"]),
         ("exome", ["INCLUDE_INDELS=true", "Q=30"]), ("hs_deep", []), ("hs_first_empty", []), ("hs_unknown", []),
         ("t1", []), ("t2", ["MQ=1"])]


@pytest.mark.parametrize("name, args", CASES)
def test_mean_and_minimum_equal_the_hs_per_target_coverage(built, corpus, exome, tmp_path, name, args):
    bam, ref, targets = inputs(name, built, corpus, exome)
    if not bam.exists():
        pytest.skip(f"{bam} not built (make corpus)")
    text, exons = annotation(read_targets(targets), contigs(bam))
    gtf = tmp_path / "exons.gtf"
    gtf.write_text(text)
    rows = ours(bam, gtf, args, tmp_path / "exons.txt")
    assert len(rows) == len(exons)
    expected = {}
    for i, layer in enumerate(layers(exons)):
        layer_list = interval_list(bam, layer, tmp_path / f"layer{i}.interval_list")
        picard = hs_per_target("picard", bam, ref, layer_list, args, tmp_path / f"picard{i}.txt")
        mine = hs_per_target(str(BAMMETRICS), bam, ref, layer_list, args, tmp_path / f"hs{i}.txt")
        assert picard == mine, f"{name} layer {i}: bammetrics CollectHsMetrics differs from Picard"
        assert set(picard) == set(layer), f"{name} layer {i}: Picard merged targets"
        expected.update(picard)
    compared = 0
    for row in rows:
        key = (row["CHROM"], int(row["START"]), int(row["END"]))
        if row["CHROM"] == MISSING_CONTIG:
            assert row["MEAN_DEPTH"] == row["MINIMUM_DEPTH"] == "", "an exon off the BAM is not counted"
            continue
        mean, minimum = expected[key]
        assert (row["MEAN_DEPTH"], row["MINIMUM_DEPTH"]) == (mean, minimum), f"{name} {key}"
        SEEN["MEAN_DEPTH"].add(mean)
        SEEN["MINIMUM_DEPTH"].add(minimum)
        compared += 1
    assert compared == len(exons) - 1
    RESULTS.append((f"{TOOL} {name} {' '.join(args)}".strip(), 2 * compared, 2 * compared,
                    "MEAN_DEPTH and MINIMUM_DEPTH equal Picard PER_TARGET_COVERAGE"))


def bedcov_means(bam, exons, tmp_path, mapq):
    """{(contig, start, end): mean depth} from samtools bedcov."""
    bed = tmp_path / "exons.bed"
    distinct = sorted({(c, s, e) for c, s, e, *_ in exons if c != MISSING_CONTIG})
    bed.write_text("".join(f"{c}\t{s - 1}\t{e}\n" for c, s, e in distinct))
    index = tmp_path / "input.bai"
    subprocess.run(["samtools", "index", "-o", str(index), str(bam)], check=True, capture_output=True)
    out = subprocess.run(["samtools", "bedcov", "-Q", str(mapq), "-j", "--max-depth", "100000000", "-X",
                          str(bed), str(bam), str(index)], check=True, capture_output=True, text=True).stdout
    means = {}
    for line in out.splitlines():
        c, s, e, total = line.split("\t")
        means[(c, int(s) + 1, int(e))] = int(total) / (int(e) - int(s))
    return means


BEDCOV_CASES = [("exome", ["MQ=0", "Q=0", "CLIP_OVERLAPPING_READS=false"], 0, True),
                ("foundation", ["MQ=0", "Q=0", "CLIP_OVERLAPPING_READS=false"], 0, True),
                ("hs_deep", ["MQ=0", "Q=0", "CLIP_OVERLAPPING_READS=false"], 0, True),
                ("exome", [], 20, False), ("foundation", [], 20, False)]


@pytest.mark.parametrize("name, args, mapq, close", BEDCOV_CASES)
def test_mean_depth_against_samtools_bedcov(built, corpus, exome, tmp_path, name, args, mapq, close):
    bam, _, targets = inputs(name, built, corpus, exome)
    text, exons = annotation(read_targets(targets), contigs(bam))
    gtf = tmp_path / "exons.gtf"
    gtf.write_text(text)
    rows = ours(bam, gtf, args, tmp_path / "exons.txt")
    means = bedcov_means(bam, exons, tmp_path, mapq)
    worst = 0.0
    for row in rows:
        if row["CHROM"] == MISSING_CONTIG:
            continue
        mine = float(row["MEAN_DEPTH"])
        theirs = means[(row["CHROM"], int(row["START"]), int(row["END"]))]
        # Six printed digits: half a unit of the last one.
        assert mine <= theirs + 5e-7 * max(1.0, theirs), f"{name}: {row['START']}: {mine} over bedcov {theirs}"
        if theirs > 0:
            worst = max(worst, (theirs - mine) / theirs)
    if close:
        assert worst <= BEDCOV_TOLERANCE, f"{name}: bedcov mean {worst:.2%} over MEAN_DEPTH"
    RESULTS.append((f"{TOOL} bedcov {name} {' '.join(args)}".strip(), len(rows) - 1, len(rows) - 1,
                    f"largest relative excess of bedcov {worst:.4%}"))


def test_mean_and_minimum_take_three_values():
    if not SEEN["MEAN_DEPTH"]:
        pytest.skip("no exon compared in this run")
    for column, values in SEEN.items():
        assert len(values) >= 3, f"{column}: {sorted(values)}"

"""CRAM input against Picard 3.5.0: every registry tool that reads a BAM, on
the same data written as CRAM.

Inputs: corpus t1 (tests/corpus/out), the foundation corpus and the exome
slice, each converted by `samtools view -C -T ref.fa` to CRAM 3.0 and to
CRAM 3.1. Both tools read the CRAM with `R=` the reference it was written
with. compare_metrics (compare_outputs for the tools with several files)
checks every cell and the file body byte for byte against Picard. Then
bammetrics on the BAM with the same arguments must write the same bodies:
the CRAM path gives the records the BAM path gives.

Error parity: a CRAM without `R=`, a reference with the same contig and
other bases (the slice MD5 fails), and a reference without the CRAM's
contigs. Each fails in both tools with no file written. STOP_AFTER runs on
the one-pass tools, and BedToIntervalList takes the CRAM as its
SEQUENCE_DICTIONARY.

Needs picard, samtools and pysam; run `make concordance`.
"""
import random
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS
from compare import body, compare_metrics, compare_outputs
from conftest import RESULTS
from run_picard import picard_env

SILENT = "VALIDATION_STRINGENCY=SILENT"
VERSIONS = ["3.0", "3.1"]
SOURCES = ["t1", "foundation", "exome"]
LEVELS = ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"]


def source(name, corpus, exome):
    """BAM, reference, baits, targets and refFlat (or None) of one input."""
    if name == "t1":
        return (CORPUS / "t1.bam", CORPUS / "reference.fa", CORPUS / "baits.interval_list",
                CORPUS / "targets.interval_list", None)
    if name == "foundation":
        return (corpus / "corpus.bam", corpus / "ref.fa", corpus / "bait.interval_list",
                corpus / "target.interval_list", corpus / "refflat.txt")
    return (exome / "exome.bam", exome / "ref.fa", exome / "baits.interval_list",
            exome / "targets.interval_list", None)


@pytest.fixture(scope="module")
def crams(tmp_path_factory, corpus, exome):
    """(name, version) -> (cram, bam, ref, baits, targets, refflat)."""
    directory = tmp_path_factory.mktemp("cram")
    made = {}
    for name in SOURCES:
        bam, ref, baits, targets, refflat = source(name, corpus, exome)
        for version in VERSIONS:
            cram = directory / f"{name}.{version}.cram"
            subprocess.run(["samtools", "view", "-C", "-T", str(ref), "--output-fmt-option",
                            f"version={version}", "-o", str(cram), str(bam)], check=True)
            assert cram.read_bytes()[:6] == b"CRAM" + bytes([3, int(version[-1])]), cram
            made[name, version] = (cram, bam, ref, baits, targets, refflat)
    return made


def run_ours(tool, bam, outputs, args, out_dir, side):
    files = [f"{option}={out_dir / (side + suffix)}" for option, suffix in outputs]
    subprocess.run([str(BAMMETRICS), tool, f"I={bam}", *files, *args], check=True,
                   capture_output=True, text=True)
    return {p.name[len(side):]: body(p) for p in out_dir.glob(side + "*")}


def same_as_bam(tool, cram, bam, outputs, args, tmp_path):
    """bammetrics writes the same bodies for the CRAM and for the BAM."""
    from_cram = run_ours(tool, cram, outputs, args, tmp_path, "cram")
    from_bam = run_ours(tool, bam, outputs, args, tmp_path, "bam")
    assert from_cram, f"{tool}: no file from the CRAM"
    assert from_cram.keys() == from_bam.keys(), (sorted(from_cram), sorted(from_bam))
    for suffix in from_cram:
        assert from_cram[suffix] == from_bam[suffix], f"{tool} file {suffix!r}: CRAM and BAM differ"


def cases(ref, baits, targets, refflat):
    """Label -> (tool, outputs, arguments) for every registry tool that reads
    a BAM. `outputs` None compares the one file `O` (or `M`) names."""
    listed = {
        "insert_size": ("CollectInsertSizeMetrics", None, ["H=/dev/null", *LEVELS]),
        "alignment_summary": ("CollectAlignmentSummaryMetrics", None, LEVELS),
        "quality_yield": ("CollectQualityYieldMetrics", None, []),
        "quality_distribution": ("QualityScoreDistribution", None, []),
        "quality_by_cycle": ("MeanQualityByCycle", None, []),
        "base_distribution": ("CollectBaseDistributionByCycle", None, []),
        "duplicates": ("CollectDuplicateMetrics", None, []),
        "gc_bias": ("CollectGcBiasMetrics", (("O", ".detail"), ("S", ".summary")),
                    ["CHART=/dev/null"]),
        "multiple": ("CollectMultipleMetrics", (("O", ""),),
                     ["PROGRAM=CollectGcBiasMetrics", "PROGRAM=CollectQualityYieldMetrics"]),
        "complexity": ("EstimateLibraryComplexity", None, []),
        "wgs": ("CollectWgsMetrics", None, []),
        "wgs_fast": ("CollectWgsMetrics", None, ["USE_FAST_ALGORITHM=true"]),
        "raw_wgs": ("CollectRawWgsMetrics", None, []),
        "non_zero_wgs": ("CollectWgsMetricsWithNonZeroCoverage", None, ["CHART=/dev/null"]),
        "hs": ("CollectHsMetrics", None, [f"BAIT_INTERVALS={baits}", f"TARGET_INTERVALS={targets}"]),
        "targeted_pcr": ("CollectTargetedPcrMetrics", None,
                         [f"AMPLICON_INTERVALS={baits}", f"TARGET_INTERVALS={targets}"]),
        "rna_seq": ("CollectRnaSeqMetrics", None,
                    [f"REF_FLAT={refflat}", "STRAND=SECOND_READ_TRANSCRIPTION_STRAND", *LEVELS]),
    }
    return {label: (tool, outputs, [SILENT, f"R={ref}", *args])
            for label, (tool, outputs, args) in listed.items()}


LABELS = list(cases("ref", "baits", "targets", "refflat"))
# Only the foundation corpus has a refFlat.
TOOL_CASES = [(name, label) for name in SOURCES for label in LABELS
              if label != "rna_seq" or name == "foundation"]


def output_key(tool):
    return "M" if tool == "CollectDuplicateMetrics" else "O"


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("name,label", TOOL_CASES)
def test_tool_on_cram(crams, tmp_path, name, label, version):
    cram, bam, ref, baits, targets, refflat = crams[name, version]
    tool, outputs, args = cases(ref, baits, targets, refflat)[label]
    if outputs is None:
        cells = compare_metrics(tool, cram, args, workdir=tmp_path, expect_file=True,
                                output_key=output_key(tool))
        outputs = ((output_key(tool), ""),)
    else:
        cells = compare_outputs(tool, cram, args, outputs, workdir=tmp_path)
    assert cells > 0, f"{tool} {name} CRAM {version}: no cell compared"
    same_as_bam(tool, cram, bam, outputs, args, tmp_path)
    RESULTS.append((f"{label} {name} CRAM {version}", cells, cells,
                    "Picard on the CRAM, bammetrics on the CRAM and the BAM"))


STOP_AFTER = [("CollectQualityYieldMetrics", []), ("CollectBaseDistributionByCycle", []),
              ("CollectWgsMetrics", []), ("CollectWgsMetrics", ["USE_FAST_ALGORITHM=true"]),
              ("CollectAlignmentSummaryMetrics", []), ("MeanQualityByCycle", [])]


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("stop_after", [1, 1000, 25000])
@pytest.mark.parametrize("tool,args", STOP_AFTER)
def test_stop_after_on_cram(crams, tmp_path, version, stop_after, tool, args):
    cram, bam, ref, *_ = crams["foundation", version]
    args = [SILENT, f"R={ref}", f"STOP_AFTER={stop_after}", *args]
    compare_metrics(tool, cram, args, workdir=tmp_path, expect_file=True)
    same_as_bam(tool, cram, bam, (("O", ""),), args, tmp_path)
    RESULTS.append((f"{tool} foundation CRAM {version} STOP_AFTER={stop_after}", 1, 1, "agrees"))


def test_multiple_metrics_stop_after_on_cram(crams, tmp_path):
    cram, bam, ref, *_ = crams["t1", "3.1"]
    args = [SILENT, f"R={ref}", "STOP_AFTER=5000"]
    assert compare_outputs("CollectMultipleMetrics", cram, args, (("O", ""),), workdir=tmp_path) > 0
    same_as_bam("CollectMultipleMetrics", cram, bam, (("O", ""),), args, tmp_path)


def test_bed_to_interval_list_takes_a_cram_dictionary(crams, tmp_path):
    cram, *_ = crams["foundation", "3.0"]
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t10\t200\ta\nchr2\t5\t50\tb\n")
    outputs = {}
    for side, cmd in (("picard", ["picard"]), ("ours", [str(BAMMETRICS)])):
        out = tmp_path / f"{side}.interval_list"
        extra = ["QUIET=true", "VERBOSITY=ERROR"] if side == "picard" else []
        subprocess.run([*cmd, "BedToIntervalList", f"I={bed}", f"SD={cram}", f"O={out}", *extra],
                       check=True, env=picard_env(), capture_output=True)
        outputs[side] = out.read_bytes()
    assert outputs["picard"] == outputs["ours"]


def wrong_bases(ref, out):
    """`ref` with every base drawn again: the same contigs, lengths and index."""
    rng = random.Random(1)
    lines = ref.read_text().split("\n")
    out.write_text("\n".join(line if line.startswith(">") else
                             "".join(rng.choice("ACGT") for _ in line) for line in lines))
    subprocess.run(["samtools", "faidx", str(out)], check=True)
    return out


def other_contig(out):
    out.write_text(">chrZ\nACGTACGTAC\n")
    subprocess.run(["samtools", "faidx", str(out)], check=True)
    return out


def both_fail(tool, cram, args, tmp_path, message):
    picard_out, ours_out = tmp_path / "picard.txt", tmp_path / "ours.txt"
    picard = subprocess.run(["picard", tool, f"I={cram}", f"O={picard_out}", SILENT, *args,
                             "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(), capture_output=True)
    assert picard.returncode != 0, f"Picard completed: {tool} {args}"
    ours = subprocess.run([str(BAMMETRICS), tool, f"I={cram}", f"O={ours_out}", *args],
                          capture_output=True, text=True)
    assert ours.returncode == 1, ours.stderr
    assert message in ours.stderr, ours.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{tool} {cram.name} {' '.join(args)[:40]}", 0, 0, "both exit non-zero, no file"))


ERROR_TOOLS = ["CollectQualityYieldMetrics", "CollectInsertSizeMetrics", "CollectMultipleMetrics"]


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("tool", ERROR_TOOLS)
def test_a_cram_without_reference_fails(crams, tmp_path, version, tool):
    cram, *_ = crams["t1", version]
    extra = ["H=/dev/null"] if tool == "CollectInsertSizeMetrics" else []
    both_fail(tool, cram, extra, tmp_path, "a CRAM input needs REFERENCE_SEQUENCE (R)")


@pytest.mark.parametrize("version", VERSIONS)
@pytest.mark.parametrize("tool", ERROR_TOOLS)
def test_a_cram_with_other_reference_bases_fails(crams, tmp_path, version, tool):
    cram, _, ref, *_ = crams["t1", version]
    wrong = wrong_bases(ref, tmp_path / "wrong.fa")
    extra = ["H=/dev/null"] if tool == "CollectInsertSizeMetrics" else []
    both_fail(tool, cram, [f"R={wrong}", *extra], tmp_path, "does not match the CRAM slice")


@pytest.mark.parametrize("tool", ERROR_TOOLS)
def test_a_cram_with_a_reference_without_its_contigs_fails(crams, tmp_path, tool):
    cram, *_ = crams["foundation", "3.1"]
    extra = ["H=/dev/null"] if tool == "CollectInsertSizeMetrics" else []
    both_fail(tool, cram, [f"R={other_contig(tmp_path / 'other.fa')}", *extra], tmp_path,
              "in the reference")


def test_a_cram_is_told_by_its_bytes_not_its_name(crams, tmp_path):
    cram, bam, ref, *_ = crams["foundation", "3.0"]
    renamed = tmp_path / "data.bam"
    renamed.write_bytes(Path(cram).read_bytes())
    args = [SILENT, f"R={ref}"]
    compare_metrics("CollectQualityYieldMetrics", renamed, args, workdir=tmp_path, expect_file=True)
    same_as_bam("CollectQualityYieldMetrics", renamed, bam, (("O", ""),), args, tmp_path)

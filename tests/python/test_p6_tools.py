"""CollectRnaSeqMetrics from Python through the registry, on the corpus the
concordance harness generates and on the committed fixture."""
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out" / "foundation"
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
COLUMNS = 30


@pytest.fixture(scope="module")
def corpus():
    if not (CORPUS / "corpus.bam").exists():
        subprocess.run(["python3", str(ROOT / "tests" / "data" / "make_corpus.py"), str(CORPUS)],
                       check=True, capture_output=True)
    return CORPUS


def write_refflat(directory):
    """One gene on chr1 with two transcripts: the fixture's reads lie in 1001..5500."""
    path = directory / "genes.refFlat"
    path.write_text("G1\tT1\tchr1\t+\t999\t60000\t1500\t59000\t3\t999,20000,40000,\t5000,30000,60000,\n"
                    "G1\tT2\tchr1\t+\t999\t60000\t1500\t59000\t2\t999,40000,\t5000,60000,\n")
    return str(path)


def test_tool_is_registered():
    assert "CollectRnaSeqMetrics" in bammetrics.tools()


def test_fixture_row_shape(tmp_path):
    rows = bammetrics.collect_rna_seq_metrics(FIXTURE, ref_flat=write_refflat(tmp_path),
                                              strand_specificity="NONE")
    assert len(rows) == 1
    row = rows[0]
    assert list(row)[:3] == ["PF_BASES", "PF_ALIGNED_BASES", "RIBOSOMAL_BASES"]
    assert list(row)[-3:] == ["SAMPLE", "LIBRARY", "READ_GROUP"]
    assert len(row) == COLUMNS
    assert row["PF_BASES"] == 1000 and row["PF_ALIGNED_BASES"] == 1000
    assert row["CODING_BASES"] == 600 and row["UTR_BASES"] == 200 and row["INTRONIC_BASES"] == 200
    assert row["RIBOSOMAL_BASES"] is None and row["PCT_RIBOSOMAL_BASES"] is None
    assert row["PCT_MRNA_BASES"] == 0.8 and row["NUM_R1_TRANSCRIPT_STRAND_READS"] == 4
    assert row["CORRECT_STRAND_READS"] == 0 and row["SAMPLE"] is None


def test_options_and_output_file(corpus, tmp_path):
    out = tmp_path / "rna.txt"
    rows = bammetrics.collect_rna_seq_metrics(
        str(corpus / "corpus.bam"), ref_flat=str(corpus / "refflat.txt"),
        strand_specificity="FIRST_READ_TRANSCRIPTION_STRAND",
        ribosomal_intervals=str(corpus / "target.interval_list"), minimum_length=100,
        end_bias_bases=20, rrna_fragment_percentage=0.5, ignore_sequence=["chrM"],
        metric_accumulation_level=["SAMPLE", "LIBRARY"], output=str(out))
    assert [r["SAMPLE"] for r in rows][:3] == [None, "S1", "S2"]
    assert {r["LIBRARY"] for r in rows[3:]} == {"LibA", "LibB", "LibC", "LibD"}
    first = rows[0]
    assert first["RIBOSOMAL_BASES"] > 0 and 0 < first["PCT_RIBOSOMAL_BASES"] < 1
    assert first["IGNORED_READS"] > 0
    assert first["CORRECT_STRAND_READS"] + first["INCORRECT_STRAND_READS"] > 0
    assert first["MEDIAN_CV_COVERAGE"] > 0 and first["MEDIAN_5PRIME_TO_3PRIME_BIAS"] > 0
    text = out.read_text()
    assert "## METRICS CLASS\tpicard.analysis.RnaSeqMetrics" in text
    assert "normalized_position\tAll_Reads.normalized_coverage\tS1.normalized_coverage" in text
    assert text.count("\n") > 100


def test_required_options_and_bad_values_raise_value_error(tmp_path):
    refflat = write_refflat(tmp_path)
    with pytest.raises(ValueError, match="STRAND_SPECIFICITY"):
        bammetrics.collect_rna_seq_metrics(FIXTURE, ref_flat=refflat)
    with pytest.raises(ValueError, match="REF_FLAT"):
        bammetrics.collect_rna_seq_metrics(FIXTURE, strand_specificity="NONE")
    with pytest.raises(ValueError, match="StrandSpecificity"):
        bammetrics.collect_rna_seq_metrics(FIXTURE, ref_flat=refflat, strand_specificity="BOTH")
    with pytest.raises(ValueError, match="IGNORE_SEQUENCE"):
        bammetrics.collect_rna_seq_metrics(FIXTURE, ref_flat=refflat, strand_specificity="NONE",
                                           ignore_sequence="chrX")
    with pytest.raises(ValueError, match="RRNA_FRAGMENT_PERCENTAGE"):
        bammetrics.collect_rna_seq_metrics(FIXTURE, ref_flat=refflat, strand_specificity="NONE",
                                           rrna_fragment_percentage=0)

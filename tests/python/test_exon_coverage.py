"""CollectExonCoverageMetrics from Python on the committed fixture: the exon rows,
the filters, both files read back by read_table with the glossary's dtypes (the
BASES_<N>X and PCT_<N>X rows stand for every MIN_DEPTH), and the usage errors.
Run: pytest tests/python after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
FIXTURE = str(ROOT / "data" / "insert_size.bam")
GTF = str(ROOT / "data" / "exons.gtf")
BED12 = str(ROOT / "data" / "exons.bed")


def test_exon_rows_of_the_fixture():
    rows = bammetrics.collect_exon_coverage_metrics(FIXTURE, annotation=GTF, min_depth=1)
    assert len(rows) == 9
    first = rows[0]
    assert (first["GENE_NAME"], first["TRANSCRIPT_ID"], first["EXON_NUMBER"]) == ("GENE_A", "TA1.2", 1)
    assert (first["START"], first["END"], first["LENGTH"]) == (1001, 1100, 100)
    assert (first["MEAN_DEPTH"], first["MINIMUM_DEPTH"], first["BASES_1X"], first["PCT_1X"]) == (1, 1, 100, 1)
    one_base = rows[6]
    assert (one_base["START"], one_base["END"], one_base["EXON_NUMBER"]) == (4050, 4050, 2)
    off_the_bam = rows[-1]
    assert off_the_bam["CHROM"] == "chrUn"
    assert off_the_bam["MEAN_DEPTH"] is None


def test_python_prints_no_note(capfd):
    rows = bammetrics.collect_exon_coverage_metrics(FIXTURE, annotation=GTF, min_depth=1)
    assert rows[-1]["CHROM"] == "chrUn", "the row of the contig the BAM lacks"
    assert capfd.readouterr().err == ""


def test_filters_choose_the_rows():
    def transcripts(**options):
        rows = bammetrics.collect_exon_coverage_metrics(FIXTURE, annotation=GTF, **options)
        return [row["TRANSCRIPT_ID"] for row in rows]

    assert transcripts(gene="GENE_B") == ["TB1.1"]
    assert transcripts(gene=["GENE_B", "ENSG_C"]) == ["TB1.1", "TC1.1", "TC1.1"]
    assert transcripts(transcript="TA2") == ["TA2.1", "TA2.1"]
    assert transcripts(gene="GENE_A", canonical=True) == ["TA1.2"] * 3
    assert transcripts(feature="CDS", gene="GENE_A") == ["TA1.2"] * 3 + ["TA2.1"] * 2
    rows = bammetrics.collect_exon_coverage_metrics(FIXTURE, annotation=BED12, min_depth=[2, 1])
    assert [row["CHROM"] for row in rows] == ["chr1"] * 5, "contig 1 read as chr1"
    assert "BASES_1X" in rows[0] and "PCT_2X" in rows[0]


def test_both_files_read_back_with_their_types(tmp_path):
    exons, genes = tmp_path / "exons.txt", tmp_path / "genes.txt"
    bammetrics.collect_exon_coverage_metrics(FIXTURE, annotation=GTF, min_depth=[1, 45], output=str(exons),
                                             summary_output=str(genes))
    table = bammetrics.read_table(str(exons))
    assert table.name == "bammetrics.qc.ExonCoverageMetrics"
    assert table.dtypes["MEAN_DEPTH"] == "float"
    assert table.dtypes["MINIMUM_DEPTH"] == "int"
    assert table.dtypes["BASES_45X"] == "int"
    # Every PCT_45X cell prints 0: a float by the glossary, an int by its text.
    assert table.dtypes["PCT_45X"] == "float"
    assert table.columns["PCT_45X"][0] == 0.0
    summary = bammetrics.read_table(str(genes))
    assert summary.name == "bammetrics.qc.GeneCoverageMetrics"
    assert summary.columns["LEVEL"][:3] == ["GENE", "TRANSCRIPT", "TRANSCRIPT"]
    assert summary.columns["TRANSCRIPT_ID"][0] is None
    assert summary.dtypes["EXON_COUNT"] == "int"
    assert summary.columns["LENGTH"][0] == 600


@pytest.mark.parametrize("options, text", [
    ({}, "ANNOTATION"),
    ({"annotation": GTF, "feature": "gene"}, "FEATURE=gene"),
    ({"annotation": GTF, "gene": "GENE_E"}, "did you mean GENE_A"),
    ({"annotation": GTF, "min_depth": 0}, "MIN_DEPTH=0"),
])
def test_usage_errors(options, text):
    with pytest.raises(ValueError, match=text):
        bammetrics.collect_exon_coverage_metrics(FIXTURE, **options)

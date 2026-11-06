"""CreateSequenceDictionary from Python: rows for each sequence, the .dict only with
output=, the Picard refusal of an existing file, and the ALT_NAMES column."""
import gzip
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
REF = ROOT / "tests" / "data" / "ref_100k.fa"
M5 = "45cf6d7380b910e1541ebad795fe3e49"


def test_rows_without_output_write_nothing(tmp_path):
    rows = bammetrics.create_sequence_dictionary(str(REF), genome_assembly="GRCh38")
    assert rows == [{
        "name": "chr1", "length": 100000, "md5": M5, "assembly": "GRCh38",
        "uri": f"file://{REF}", "species": None, "alternate_names": None,
    }]
    assert not (REF.parent / "ref_100k.dict").exists()


def test_output_is_written_once(tmp_path):
    out = tmp_path / "ref.dict"
    bammetrics.create_sequence_dictionary(str(REF), output=str(out), uri="ref.fa",
                                          create_md5_file=True)
    assert out.read_text() == f"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\tM5:{M5}\tUR:ref.fa\n"
    assert len((tmp_path / "ref.dict.md5").read_text()) == 32
    with pytest.raises(OSError, match="already exists"):
        bammetrics.create_sequence_dictionary(str(REF), output=str(out))


def test_gzip_fasta_and_alternative_names(tmp_path):
    fasta = tmp_path / "x.fa.gz"
    fasta.write_bytes(gzip.compress(b">chr1 desc\nacgt\n>chrM\nNN\n"))
    alt = tmp_path / "alt.tsv"
    alt.write_text("chr1\t1\nchrM\tMT\nchrM\tM\n")
    rows = bammetrics.create_sequence_dictionary(str(fasta), alt_names=str(alt), num_sequences=5)
    assert [(r["name"], r["length"], r["alternate_names"]) for r in rows] == [
        ("chr1", 4, "1"), ("chrM", 2, "MT,M")]
    assert rows[0]["md5"] == "f1f8f4bf413b16ad135722aa4591043e"
    with pytest.raises(OSError, match="doesn't match regex"):
        bammetrics.create_sequence_dictionary(str(fasta), truncate_names_at_whitespace=False)

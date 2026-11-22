"""The SQL dump (docs/formats.md "SQL dump"): qc(also="sql") and chain(also="sql")
write FILE.sql beside every metrics file, Python's sqlite3 loads each one with the
glossary types, and read_table refuses a dump."""
import sqlite3

import pytest

import bammetrics

BAM = "tests/data/insert_size.bam"


def load(path):
    db = sqlite3.connect(":memory:")
    with open(path, encoding="utf-8") as dump:
        db.executescript(dump.read())
    return db


def test_qc_also_sql_loads_in_sqlite(tmp_path):
    prefix = tmp_path / "s"
    run = bammetrics.qc(BAM, output=prefix, also="sql")
    dumps = sorted(path for path in run.files if path.endswith(".sql"))
    assert str(prefix) + ".insert_size_metrics.sql" in dumps
    db = load(str(prefix) + ".insert_size_metrics.sql")
    table = bammetrics.read_table(str(prefix) + ".insert_size_metrics")
    row = db.execute('SELECT "MEDIAN_INSERT_SIZE", "READ_PAIRS", "PAIR_ORIENTATION" '
                     'FROM "InsertSizeMetrics"').fetchall()
    assert row == [(table.columns["MEDIAN_INSERT_SIZE"][0], table.columns["READ_PAIRS"][0],
                    table.columns["PAIR_ORIENTATION"][0])]
    assert isinstance(row[0][0], float) and isinstance(row[0][1], int)
    types = dict(db.execute('SELECT "column", "type" FROM "columns" '
                            "WHERE \"table\" = 'InsertSizeMetrics'").fetchall())
    assert types == table.dtypes
    tables = db.execute('SELECT "name", "kind", "sample" FROM "tables"').fetchall()
    assert tables == [("InsertSizeMetrics", "metrics", run.sample),
                      ("histogram", "histogram", run.sample)]
    for dump in dumps:
        assert load(dump).execute('SELECT count(*) FROM "headers"').fetchone()[0] == 2


def test_chain_also_sql(tmp_path):
    run = bammetrics.chain(BAM, ["insert_size", "quality_yield"], output=str(tmp_path / "c"),
                           also="sql")
    dumps = [path for path in run.files if path.endswith(".sql")]
    assert len(dumps) == 2
    names = {load(path).execute('SELECT "name" FROM "tables"').fetchone()[0] for path in dumps}
    assert names == {"InsertSizeMetrics", "QualityYieldMetrics"}


def test_read_table_refuses_a_dump():
    with pytest.raises(bammetrics.InputError, match="SQL dump"):
        bammetrics.read_table("tests/data/insert_size_metrics.sql")

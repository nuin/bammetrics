"""Calls checked by mypy against bammetrics.pyi (`make wheel-check` runs
`mypy bammetrics.pyi tests/python/stub_types.py`); pytest never imports it.
Each assert_type is the type a caller gets for one form of the API."""
from typing import Any, Literal, assert_type

import bammetrics as bm


def forms(names: list[str]) -> None:
    assert_type(bm.collect("x.bam"), dict[str, Any])
    assert_type(bm.collect("x.bam", ["insert_size", "quality_yield"]), dict[str, Any])
    assert_type(bm.collect("x.bam", ["insert-size"]), bm.Run[None])
    assert_type(bm.collect("x.bam", metrics=["gc-bias"], reference="r.fa", output="o"),
                bm.Run[None])
    assert_type(bm.collect("x.bam", ["wgs"], reference="r.fa", output="o"), bm.Run[None])
    assert_type(bm.collect("x.bam", names), dict[str, Any] | bm.Run[None])
    run = bm.qc("x.bam", output="o", preset="clinical-panel")
    assert_type(run, bm.Run[bm.Verdict])
    assert_type(run.verdict.status, Literal["pass", "warn", "fail", "incomplete"])
    assert_type(bm.check("dir"), bm.Verdict | list[bm.Verdict])
    assert_type(bm.report("dir"), list[str])
    try:
        bm.qc("x.bam")
    except bm.InputError as e:
        assert_type(e.run, bm.Run[Any] | None)
    except bm.VerdictError as e:
        assert_type(e.run, bm.Run[bm.Verdict] | None)
    assert_type(bm.tools()[0].verb, str)
    assert_type(bm.presets()[0].hash, str)
    row = bm.bam_index_stats("x.bam")[0]
    assert_type(row["reference"], str)
    assert_type(row["aligned"], int)
    assert_type(bm.create_sequence_dictionary("r.fa")[0], dict[str, Any])

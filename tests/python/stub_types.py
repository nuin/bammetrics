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
    assert_type(bm.build_bam_index("x.bam", output="x.bai")[0]["unaligned"], int)
    assert_type(bm.create_sequence_dictionary("r.fa")[0], dict[str, Any])
    assert_type(bm.view_sam("x.sam", alignment_status="Aligned", header_only=False), list[dict[str, Any]])
    table = bm.read_table("x.insert_size_metrics")
    assert_type(table, bm.Table)
    assert_type(table.dtypes["READ_PAIRS"], Literal["int", "float", "string"])
    assert_type(bm.read_table("x.json", table="histogram").histogram, bm.Table | None)
    assert_type(bm.qc("x.bam", also=["json", "markdown"]), bm.Run[bm.Verdict])
    assert_type(bm.qc("x.bam", also="junit"), bm.Run[bm.Verdict])
    step = bm.step("mark_duplicates", remove_duplicates=True, as_="dedup")
    assert_type(step, bm.Step)
    assert_type(step.role, Literal["transform", "tap", "write"])
    assert_type(bm.Step("write", output="x.bam", index=True).args, list[str])
    chained = bm.chain("x.bam", [step, "insert_size", bm.step("write")], output="o",
                       reference="r.fa", also=["json", "markdown"], format="picard,json")
    assert_type(chained, bm.Run[None])
    assert_type(chained.verdict, None)
    assert_type(chained.stages, list[bm.Stage])
    assert_type(chained.stages[0].records_in, int | None)
    assert_type(chained.stages[0].status, Literal["ok", "failed", "stopped", "not_run"])
    assert_type(bm.chain("x.bam", names, threads=2), bm.Run[None])
    assert_type(bm.chain(plan="plan.toml", output="o"), bm.Run[None])
    assert_type(bm.chain("x.bam", "insert_size", dry_run=True), str)
    assert_type(bm.chain("x.bam", step, dry_run=bool(names)), bm.Run[None] | str)
    assert_type(run.stages, list[bm.Stage])

"""What chain, qc and collect refuse before their first read, from Python
(docs/chaining.md section 2.5; crates/cli/tests/chain_safety.rs holds the
same cases for the command): a thread count over 256, a value only the
tool can judge, an output through a directory that does not exist and
back, a plan file over 1 MiB. Each is a UsageError and writes nothing."""
from pathlib import Path

import pytest

import bammetrics as bm

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
BAM = DATA / "insert_size.bam"
MARKDUP = DATA / "markdup.bam"

OVER = "more than the maximum of 256 threads"


@pytest.mark.parametrize("threads", [257, 100_000_000])
@pytest.mark.parametrize("call", [
    lambda out, threads: bm.chain(BAM, ["insert_size"], output=out, threads=threads),
    lambda out, threads: bm.qc(BAM, output=out, threads=threads),
    lambda out, threads: bm.collect(BAM, metrics=["insert-size"], output=out, threads=threads),
    lambda out, threads: bm.collect(str(BAM), ["insert_size"], output_prefix=str(out),
                                    threads=threads),
], ids=["chain", "qc", "collect-verbs", "collect-dict"])
def test_a_thread_count_over_the_maximum_is_a_usage_error(tmp_path, call, threads):
    with pytest.raises(bm.UsageError, match=OVER) as caught:
        call(tmp_path / "s", threads)
    assert str(threads) in str(caught.value)
    assert list(tmp_path.iterdir()) == []


def test_the_maximum_thread_count_runs(tmp_path):
    run = bm.chain(BAM, ["quality_yield"], output=tmp_path / "s", threads=256)
    assert [stage.status for stage in run.stages] == ["ok"]


@pytest.mark.parametrize("steps, message", [
    ([bm.step("insert_size", metric_accumulation_level=["NOPE"])],
     r'step 1 \(insert-size\): METRIC_ACCUMULATION_LEVEL "NOPE": expected ALL_READS'),
    (["mark_duplicates", bm.step("alignment_summary", expected_pair_orientations=["UP"])],
     r'step 2 \(alignment-summary\): pair orientation "UP": expected FR, RF or TANDEM'),
])
def test_a_value_only_the_tool_can_judge_is_a_usage_error(tmp_path, steps, message):
    """In 0.22.0 the tap failed alone when it was built: an InputError
    after the other steps had run."""
    with pytest.raises(bm.UsageError, match=message):
        bm.chain(MARKDUP, steps, output=tmp_path / "made" / "s")
    assert list(tmp_path.iterdir()) == []


def test_an_output_through_a_missing_directory_and_back_is_refused(tmp_path):
    source = tmp_path / "in.bam"
    source.write_bytes(BAM.read_bytes())
    onto_input = tmp_path / "missing" / ".." / "in.bam"
    for step in (bm.step("write", output=onto_input), bm.step("insert_size", output=onto_input)):
        with pytest.raises(bm.UsageError, match=r"does not exist, so the \.\. behind it"):
            bm.chain(source, [step])
    assert [p.name for p in tmp_path.iterdir()] == ["in.bam"]
    assert source.read_bytes() == BAM.read_bytes()


def test_a_plan_file_over_one_mebibyte_is_refused(tmp_path):
    plan = tmp_path / "plan.toml"
    text = f'[chain]\nsource = "{BAM}"\n\n[[step]]\nverb = "quality-yield"\n# '
    plan.write_text(text + "x" * ((1 << 20) + 1 - len(text)))
    with pytest.raises(bm.UsageError, match="not a plan file: it is larger than 1 MiB"):
        bm.chain(plan=plan, output=tmp_path / "s")
    assert [p.name for p in tmp_path.iterdir()] == ["plan.toml"]

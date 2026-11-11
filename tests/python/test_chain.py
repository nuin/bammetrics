"""chain and step from Python (docs/python.md, docs/chaining.md) on the
committed fixtures: chained equals staged, the plan file equals the inline
steps, the bytes the command writes, Run.stages and Run.metrics, one error
per misuse, pickling and the GIL release. The byte comparisons with the
command need target/release/bammetrics (`cargo build --release`)."""
import ast
import concurrent.futures
import copy
import pickle
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

import bammetrics as bm
from test_mirror import ticks_during

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
MARKDUP = DATA / "markdup.bam"
BAM = DATA / "insert_size.bam"
QUERYNAME = DATA / "duplicate_marking.bam"
ARTIFACT = DATA / "artifact.bam"
REFERENCE = DATA / "ref_100k.fa"
STUB = ROOT / "bammetrics.pyi"
CLI = ROOT / "target" / "release" / "bammetrics"

needs_cli = pytest.mark.skipif(not CLI.exists(), reason="needs `cargo build --release`")

# The four collectors of both documented examples, as chain steps and as
# the extension qc gives each file.
TAPS = {
    "alignment_summary": ".alignment_summary_metrics",
    "insert_size": ".insert_size_metrics",
    "quality_yield": ".quality_yield_metrics",
    "duplicate_metrics": ".duplicate_metrics",
}


def body(path):
    """A metrics file without its StringHeader lines (command line, clock)."""
    lines = Path(path).read_text().splitlines()
    return [line for line in lines
            if not line.startswith("# ") and line != "## htsjdk.samtools.metrics.StringHeader"]


def masked(path):
    """A written file with the Started on line (the run's clock) masked."""
    lines = Path(path).read_bytes().split(b"\n")
    return b"\n".join(b"# Started on: <time>" if line.startswith(b"# Started on:") else line
                      for line in lines)


def bam(path):
    """The header lines of a BAM without the CL field of each @PG line, and
    every record as SAM text."""
    with pysam.AlignmentFile(str(path)) as reads:
        lines = [
            "\t".join(f for f in line.split("\t")
                      if not (line.startswith("@PG") and f.startswith("CL:")))
            for line in str(reads.header).splitlines()
        ]
        return lines, [record.to_string() for record in reads]


def names(directory):
    return sorted(p.name for p in Path(directory).iterdir())


def markdup_steps():
    """The first documented example: mark-duplicates, four collectors on the
    marked records, the marked BAM with its index."""
    return [bm.step("mark_duplicates"), *TAPS, bm.step("write", index=True)]


def taps_steps():
    """The second documented example: four collectors with options."""
    return [
        bm.step("quality_yield", use_original_qualities=False),
        bm.step("alignment_summary", metric_accumulation_level=["SAMPLE", "LIBRARY"],
                stop_after=8),
        bm.step("insert_size", minimum_pct=0.1),
        "duplicate_metrics",
    ]


def test_chained_equals_staged(tmp_path):
    staged = tmp_path / "staged.marked.bam"
    rows = bm.mark_duplicates(MARKDUP, output=staged,
                              metrics_file=tmp_path / "staged.marked.duplication_metrics")
    bm.collect_alignment_summary_metrics(staged, reference_sequence=REFERENCE,
                                         output=tmp_path / "staged.alignment_summary_metrics")
    bm.collect_insert_size_metrics(staged, output=tmp_path / "staged.insert_size_metrics")
    bm.collect_quality_yield_metrics(staged, output=tmp_path / "staged.quality_yield_metrics")
    bm.collect_duplicate_metrics(staged, output=tmp_path / "staged.duplicate_metrics")
    bm.collect_gc_bias_metrics(staged, reference_sequence=REFERENCE,
                               output=tmp_path / "staged.gc_bias.detail_metrics",
                               summary_output=tmp_path / "staged.gc_bias.summary_metrics")

    steps = [bm.step("mark_duplicates"), *TAPS, "gc_bias", "write"]
    run = bm.chain(MARKDUP, steps, reference=REFERENCE, output=tmp_path / "chain")
    assert isinstance(run, bm.Run) and run.exit_code == 0 and run.errors == {}
    staged_metrics = [n for n in names(tmp_path) if n.startswith("staged.") and "metrics" in n]
    assert len(staged_metrics) == 7
    for name in staged_metrics:
        chained = tmp_path / name.replace("staged.", "chain.", 1)
        assert body(chained) == body(tmp_path / name), name
    chain_header, chain_records = bam(tmp_path / "chain.marked.bam")
    staged_header, staged_records = bam(staged)
    assert chain_header == staged_header
    assert chain_records == staged_records and len(chain_records) == 36
    # The rows the staged function returned are the rows of the chain's file.
    marked = run.metrics["mark-duplicates"]["classes"][0]["rows"]
    assert [r["LIBRARY"] for r in marked] == [r["LIBRARY"] for r in rows]
    assert [r["READ_PAIR_DUPLICATES"] for r in marked] == [r["READ_PAIR_DUPLICATES"] for r in rows]
    assert run.metrics["insert-size"] == bm.read_metrics(tmp_path / "chain.insert_size_metrics")


def test_a_chain_without_a_transform_writes_the_files_of_qc(tmp_path):
    run = bm.chain(BAM, ["quality_yield", "alignment_summary", "insert_size",
                         "duplicate_metrics"], output=tmp_path / "c" / "s")
    bm.qc(BAM, output=tmp_path / "q" / "s", report=False)
    assert sorted(Path(f).name for f in run.files) == sorted("s" + ext for ext in TAPS.values())
    for extension in TAPS.values():
        chained = masked(tmp_path / "c" / ("s" + extension))
        assert chained.replace(b"/c/s.", b"/q/s.") == masked(tmp_path / "q" / ("s" + extension))
    assert set(run.metrics) == {"quality-yield", "alignment-summary", "insert-size",
                                "duplicate-metrics"}
    assert run.verdict is None and run.manifest is None and run._repr_html_() is None
    assert [s.status for s in run.stages] == ["ok"] * 4
    assert all(s.records_in is None and s.segment == 1 for s in run.stages)


@needs_cli
def test_chain_writes_the_bytes_the_command_writes(tmp_path, monkeypatch):
    a, b = tmp_path / "cli", tmp_path / "py"
    a.mkdir()
    b.mkdir()
    done = subprocess.run(
        [str(CLI), "chain", str(MARKDUP), "-R", str(REFERENCE), "-o", "s", "--also", "json",
         "+", "quality-yield", "+", "mark-duplicates", "--remove-duplicates",
         "+", "insert-size", "--minimum-pct", "0.1", "--level", "SAMPLE", "--level", "LIBRARY",
         "+", "write", "--index", "--md5"],
        cwd=a, capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    monkeypatch.chdir(b)
    run = bm.chain(MARKDUP, [
        "quality_yield",
        bm.step("mark_duplicates", remove_duplicates=True),
        bm.step("insert_size", minimum_pct=0.1,
                metric_accumulation_level=["SAMPLE", "LIBRARY"]),
        bm.step("write", index=True, md5=True),
    ], reference=REFERENCE, output="s", also="json")
    assert names(a) == names(b) == sorted(Path(f).name for f in run.files)
    assert "s.marked.bam.md5" in names(a) and "s.in.quality_yield_metrics.json" in names(a)
    for name in names(a):
        assert masked(a / name) == masked(b / name), name
    assert run.prefix == "s"


@pytest.mark.parametrize("plan, steps, keywords", [
    ("chain_markdup.toml", markdup_steps,
     {"source": "tests/data/markdup.bam", "reference": "tests/data/ref_100k.fa"}),
    ("chain_taps.toml", taps_steps,
     {"source": "tests/data/insert_size.bam", "format": ["picard", "json"], "sample": "S1"}),
])
def test_the_plan_file_equals_the_inline_steps(tmp_path, monkeypatch, plan, steps, keywords):
    # The plan files name their inputs from the repository root.
    monkeypatch.chdir(ROOT)
    a, b = tmp_path / "plan", tmp_path / "inline"
    from_plan = bm.chain(plan=DATA / plan, output=a / "s")
    inline = bm.chain(steps=steps(), output=b / "s", **keywords)
    assert names(a) == names(b) and len(names(a)) == len(from_plan.files)
    for name in names(a):
        if name.endswith(".bam"):
            # The @PG CL names the metrics file under each prefix.
            assert bam(a / name) == bam(b / name)
        elif not name.endswith(".bai"):
            assert masked(a / name).replace(b"/plan/s.", b"/inline/s.") == masked(b / name), name
    assert [(s.verb, s.role, s.status, s.reads, s.leaves, s.segment, s.records_in,
             s.records_out, s.changed) for s in from_plan.stages] == [
        (s.verb, s.role, s.status, s.reads, s.leaves, s.segment, s.records_in,
         s.records_out, s.changed) for s in inline.stages]
    assert list(from_plan.metrics) == list(inline.metrics)
    assert from_plan.sample == inline.sample
    text = bm.chain(plan=DATA / plan, output=a / "s", dry_run=True)
    assert text.replace("/plan/s", "/inline/s") == bm.chain(
        steps=steps(), output=b / "s", dry_run=True, **keywords)


def test_a_keyword_replaces_the_plan_files_value(tmp_path, monkeypatch):
    monkeypatch.chdir(ROOT)
    run = bm.chain(QUERYNAME, plan=DATA / "chain_markdup.toml", output=tmp_path / "s",
                   reference=None)
    assert run.stages[0].records_in == 13
    with pysam.AlignmentFile(str(tmp_path / "s.marked.bam")) as marked:
        assert marked.header["HD"]["SO"] == "queryname"
    with pytest.raises(bm.UsageError, match="--plan gives the steps"):
        bm.chain(plan=DATA / "chain_markdup.toml", steps=["insert_size"], output=tmp_path / "t")
    with pytest.raises(bm.InputError, match=r"nope.toml: .* \(from --plan\)"):
        bm.chain(plan=tmp_path / "nope.toml")
    bad = tmp_path / "bad.toml"
    bad.write_text('[chain]\nsorce = "x.bam"\n')
    with pytest.raises(bm.UsageError, match=r"bad.toml:2: \[chain\]: unknown key sorce\n"
                                            "  did you mean source"):
        bm.chain(plan=bad)
    # write --index writes no .bai of a stream that is not coordinate sorted.
    assert names(tmp_path) == ["bad.toml", "s.alignment_summary_metrics", "s.duplicate_metrics",
                               "s.insert_size_metrics", "s.marked.bam",
                               "s.marked.duplication_metrics", "s.quality_yield_metrics"]


def test_no_step_writes_the_plan_file(tmp_path):
    """The chain reads its plan file, so a step that names it as its output
    (by its name or through a link) is a UsageError and the plan stays whole;
    a plan that is not text is one too."""
    plan = tmp_path / "plan.toml"
    link = tmp_path / "link.toml"
    link.symlink_to(plan)
    for output in (plan, link):
        text = (f'[chain]\nsource = "{BAM}"\noutput-prefix = "{tmp_path / "s"}"\n'
                f'[[step]]\nverb = "quality-yield"\n[[step]]\nverb = "insert-size"\n'
                f'output = "{output}"\n')
        plan.write_text(text)
        for given in (plan, link):
            with pytest.raises(bm.UsageError, match=r"step 2 \(insert-size\): .*\.toml is an "
                                                    "input; writing it would destroy it"):
                bm.chain(plan=given)
            assert plan.read_text() == text
    binary = tmp_path / "binary.toml"
    binary.write_bytes(b"\xff\xfe\x00")
    with pytest.raises(bm.UsageError, match=r"binary.toml: not a plan file: it is not UTF-8 "
                                            r"text \(from --plan\)"):
        bm.chain(plan=binary)
    assert names(tmp_path) == ["binary.toml", "link.toml", "plan.toml"]


def test_panel_names_the_baits_and_the_targets_of_hs_metrics(tmp_path):
    """The command's error for missing intervals names --panel, so a step
    takes panel= too: one file as BAIT_INTERVALS and TARGET_INTERVALS."""
    targets = DATA / "targets.interval_list"
    with pytest.raises(bm.UsageError, match="Use --panel FILE when baits and targets are the "
                                            "same file"):
        bm.chain(BAM, ["hs_metrics"], output=tmp_path / "none")
    step = bm.step("hs_metrics", panel=targets, sample_size=1000)
    assert step.args == [f"--panel={targets}", "--sample-size=1000"]
    one = bm.chain(BAM, [step], reference=REFERENCE, output=tmp_path / "one")
    bm.chain(BAM, [bm.step("hs_metrics", bait_intervals=targets, target_intervals=targets,
                           sample_size=1000)],
             reference=REFERENCE, output=tmp_path / "two")
    assert list(one.metrics) == ["hs-metrics"]
    assert body(tmp_path / "one.hs_metrics") == body(tmp_path / "two.hs_metrics")
    assert all(copied == step for copied in round_trips(step))


def test_stages_hold_each_steps_status_files_and_record_counts(tmp_path):
    run = bm.chain(MARKDUP, ["quality_yield", bm.step("mark_duplicates", as_="dedup",
                                                      remove_duplicates=True),
                             "insert_size", bm.step("write", index=True)],
                   output=tmp_path / "s")
    assert [(s.number, s.verb, s.role, s.reads, s.leaves, s.segment) for s in run.stages] == [
        (1, "quality-yield", "tap", "in", None, 1),
        (2, "mark-duplicates", "transform", "in", "dedup", 1),
        (3, "insert-size", "tap", "dedup", None, 2),
        (4, "write", "write", "dedup", None, 2),
    ]
    assert all(s.status == "ok" and s.error is None for s in run.stages)
    marked = run.stages[1]
    assert (marked.records_in, marked.records_out, marked.changed) == (36, 19, 16)
    with pysam.AlignmentFile(str(tmp_path / "s.dedup.bam")) as written:
        assert sum(1 for _ in written) == marked.records_out
    assert [s.records_in for s in run.stages if s.role != "transform"] == [None] * 3
    assert [[Path(f).name for f in s.files] for s in run.stages] == [
        ["s.in.quality_yield_metrics"], ["s.dedup.duplication_metrics"],
        ["s.insert_size_metrics"], ["s.dedup.bam", "s.dedup.bai"]]
    assert run.files == [f for s in run.stages for f in s.files]
    assert repr(marked) == ('Stage(number=2, verb="mark-duplicates", role="transform", '
                            'status="ok", records_in=36, records_out=19, changed=16)')
    # A qc or collect run has no stage.
    assert bm.collect(BAM, metrics=["insert-size"], output=tmp_path / "c").stages == []


def test_metrics_are_keyed_by_verb_stream_and_file(tmp_path):
    run = bm.chain(MARKDUP, ["quality_yield", "mark_duplicates", "quality_yield", "gc_bias"],
                   reference=REFERENCE, output=tmp_path / "s")
    assert list(run.metrics) == ["in.quality-yield", "mark-duplicates", "quality-yield",
                                 "gc-bias", "gc-bias.summary"]
    for key, name in [("in.quality-yield", "s.in.quality_yield_metrics"),
                      ("mark-duplicates", "s.marked.duplication_metrics"),
                      ("gc-bias.summary", "s.gc_bias.summary_metrics")]:
        assert repr(run.metrics[key]) == repr(bm.read_metrics(tmp_path / name)), key
    # A tool that names its files when it finishes, and two taps of one tool
    # on one stream: no file is lost.
    run = bm.chain(ARTIFACT, [
        "sequencing_artifact", "quality_yield",
        bm.step("quality_yield", output=tmp_path / "second.txt", use_original_qualities=False),
    ], reference=REFERENCE, output=tmp_path / "a")
    assert list(run.metrics) == [
        "sequencing-artifact.pre_adapter_detail_metrics",
        "sequencing-artifact.pre_adapter_summary_metrics",
        "sequencing-artifact.bait_bias_detail_metrics",
        "sequencing-artifact.bait_bias_summary_metrics",
        "sequencing-artifact.error_summary_metrics",
        "quality-yield", "quality-yield#3"]
    assert run.metrics["quality-yield#3"] == bm.read_metrics(tmp_path / "second.txt")


def test_formats_and_also_add_files_beside_each_metrics_file(tmp_path):
    run = bm.chain(BAM, ["insert_size"], output=tmp_path / "s", format=["picard", "json"],
                   also=["tsv", "markdown"], sample="S1")
    assert [Path(f).name for f in run.files] == [
        "s.insert_size_metrics", "s.insert_size_metrics.json", "s.insert_size_metrics.tsv",
        "s.report.md"]
    assert run.sample == "S1" and list(run.metrics) == ["insert-size"]
    assert bm.read_table(tmp_path / "s.insert_size_metrics.json").sample == "S1"
    both = bm.chain(BAM, ["quality_yield", bm.step("insert_size", format="picard,csv")],
                    output=tmp_path / "t", format="picard")
    assert [Path(f).name for f in both.files] == [
        "t.quality_yield_metrics", "t.insert_size_metrics", "t.insert_size_metrics.csv"]


def test_dry_run_returns_the_plan_and_runs_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    text = bm.chain(MARKDUP, markdup_steps(), reference=REFERENCE, output="s", dry_run=True)
    assert isinstance(text, str)
    assert text.startswith(f"chain plan for {MARKDUP}\n")
    assert "segment 2  a second read of" in text and "s.marked.bam, s.marked.bai" in text
    assert text.endswith("nothing ran: this is a plan\n")
    assert names(tmp_path) == []
    if CLI.exists():
        done = subprocess.run(
            [str(CLI), "chain", str(MARKDUP), "-R", str(REFERENCE), "-o", "s", "--dry-run",
             "+", "mark-duplicates", "+", "alignment-summary", "+", "insert-size",
             "+", "quality-yield", "+", "duplicate-metrics", "+", "write", "--index"],
            cwd=tmp_path, capture_output=True, text=True)
        assert (done.returncode, done.stdout) == (0, text)


def test_a_step_is_checked_when_it_is_made():
    step = bm.step("mark_duplicates", remove_duplicates=True, comment=["a b", "c"], as_="dedup",
                   barcode_tag=None)
    assert isinstance(step, bm.Step)
    assert (step.verb, step.role) == ("mark-duplicates", "transform")
    assert step.options == {"remove_duplicates": True, "comment": ("a b", "c"), "as_": "dedup"}
    assert step.args == ["--remove-duplicates=true", "--comment=a b", "--comment=c",
                         "--as=dedup"]
    assert repr(step) == ("Step(\"mark-duplicates\", remove_duplicates=True, "
                          "comment=('a b', 'c'), as_='dedup')")
    # A verb reads in snake or kebab case, as its Picard class and as the name of
    # its function; the class and the function are other names of one step.
    forms = ["insert_size", "insert-size", "CollectInsertSizeMetrics",
             "collect_insert_size_metrics"]
    assert {bm.step(name, minimum_pct=0.1) for name in forms} == {bm.Step("insert_size",
                                                                           minimum_pct=0.1)}
    assert bm.step("insert_size") != bm.step("insert_size", minimum_pct=0.1)
    assert not bm.step("insert_size", minimum_pct=0.1) != bm.step("insert-size", minimum_pct=0.1)
    assert bm.step("insert_size") != "insert_size"
    write = bm.step("write", output=Path("out") / "x.bam", index=True, compression_level=9)
    assert (write.role, write.args) == ("write", ["--output=out/x.bam", "--index=true",
                                                  "--compression-level=9"])
    assert bm.step("quality_yield", use_original_qualities=False).args == [
        "--use-original-qualities=false"]


@pytest.mark.parametrize("call, error, message", [
    (lambda: bm.step("inser_size"), bm.UsageError,
     r'chain: unknown step "inser_size"\n  did you mean insert_size\?\n  steps: insert_size, '),
    (lambda: bm.step("inser-size"), bm.UsageError, r"did you mean insert-size\?"),
    (lambda: bm.Step("mark_duplicats"), bm.UsageError, r"did you mean mark_duplicates\?"),
    (lambda: bm.step("wgs"), bm.UsageError,
     "step wgs: CollectWgsMetrics has no stage and does not read a BAM in one pass"),
    (lambda: bm.step("insert_size", minmum_pct=0.1), bm.UsageError,
     r'step insert-size: unknown keyword "minmum_pct"\n  did you mean minimum_pct\?\n'
     "  keywords: deviations, "),
    (lambda: bm.step("write", indx=True), bm.UsageError, r"did you mean index\?"),
    (lambda: bm.step("insert_size", input="x.bam"), bm.UsageError,
     "input is not a keyword of a step: the chain's source feeds every step"),
    (lambda: bm.step("mark_duplicates", output="x.bam"), bm.UsageError,
     r'mark-duplicates writes no BAM in a chain, add step\("write", output=...\) after it'),
    (lambda: bm.step("insert_size", as_="x"), bm.UsageError, 'unknown keyword "as_"'),
    (lambda: bm.step("insert_size", minimum_pct="0.1"), bm.UsageError,
     "step insert-size: minimum_pct takes a number, not str"),
    (lambda: bm.step("insert_size", histogram_width=1.5), bm.UsageError,
     "histogram_width takes an int, not float"),
    (lambda: bm.step("insert_size", include_duplicates=1), bm.UsageError,
     "include_duplicates takes True or False, not int"),
    (lambda: bm.step("insert_size", deviations=[1, 2]), bm.UsageError,
     "deviations takes one value, not a list"),
    (lambda: bm.step("insert_size", output=3), TypeError,
     "step insert-size: output takes a str or os.PathLike, not int"),
    (lambda: bm.step("mark_duplicates", comment="a\0b"), bm.UsageError,
     "comment cannot hold a NUL byte"),
    # A float no option can hold, and a list of texts with another type in it:
    # found here, not when the tap is built.
    (lambda: bm.step("insert_size", minimum_pct=float("nan")), bm.UsageError,
     "step insert-size: minimum_pct takes a finite number, not nan"),
    (lambda: bm.step("insert_size", minimum_pct=float("-inf")), bm.UsageError,
     "minimum_pct takes a finite number, not -inf"),
    (lambda: bm.step("mark_duplicates", program_record_id=float("inf")), bm.UsageError,
     "program_record_id takes a finite number, not inf"),
    (lambda: bm.step("insert_size", metric_accumulation_level=["SAMPLE", 1]), bm.UsageError,
     "step insert-size: metric_accumulation_level takes a list of str, not one holding int"),
    (lambda: bm.step("mark_duplicates", comment=("a", 2.5)), bm.UsageError,
     "comment takes a list of str, not one holding float"),
    (lambda: bm.step("insert_size", panel="targets.interval_list"), bm.UsageError,
     'step insert-size: unknown keyword "panel"'),
    (lambda: bm.step(3), TypeError, "verb"),
])
def test_a_misused_step_raises_when_it_is_made(call, error, message):
    with pytest.raises(error, match=message):
        call()


def stub_function(name):
    tree = ast.parse(STUB.read_text())
    return next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name)


def test_a_step_takes_the_keywords_of_its_tools_function():
    """Every keyword the generated registry stub gives a tool's function is a
    keyword of its step, apart from the BAM a transform leaves to write."""
    steps = {tool.verb.replace("-", "_"): tool for tool in bm.tools()}
    checked = 0
    for verb in steps:
        try:
            bm.step(verb)
        except bm.UsageError:
            continue
        function = stub_function("".join(
            "_" + c.lower() if c.isupper() else c for c in steps[verb]).lstrip("_"))
        for keyword in (a.arg for a in function.args.kwonlyargs):
            if (verb, keyword) == ("mark_duplicates", "output"):
                with pytest.raises(bm.UsageError, match="writes no BAM in a chain"):
                    bm.step(verb, output=None)
            else:
                assert bm.step(verb, **{keyword: None}).args == []
            checked += 1
    assert checked > 150


@pytest.mark.parametrize("call, error, message", [
    (lambda out: bm.chain(BAM, ["insert_size"], outptu=out), bm.UsageError,
     r'chain: unknown keyword "outptu"\n  did you mean output\?'),
    (lambda out: bm.chain(BAM, output=out), bm.UsageError, "no step: a chain needs at least one"),
    (lambda out: bm.chain(steps=["insert_size"], output=out), bm.UsageError, "missing BAM"),
    (lambda out: bm.chain(out + ".nope.bam", ["insert_size"], output=out), bm.InputError,
     r"nope.bam not found \(from BAM\)"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, reference=out + ".fa"),
     bm.InputError, r"\.fa not found \(from -R\)"),
    (lambda out: bm.chain(BAM, ["inser_size"], output=out), bm.UsageError,
     "did you mean insert_size"),
    (lambda out: bm.chain(BAM, [3], output=out), TypeError,
     r"steps takes verb names and step\(...\) objects, not int"),
    (lambda out: bm.chain(5, ["insert_size"], output=out), TypeError,
     "chain: source takes a str or os.PathLike, not int"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, threads="2"), bm.UsageError,
     "chain: threads takes an int, not str"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, dry_run="yes"), bm.UsageError,
     "dry_run takes True or False, not str"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, also="junit"), bm.UsageError,
     "--also junit prints a verdict, and a chain has none"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, format="json"), bm.UsageError,
     "--format json starts with json: a chain names Picard files"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, sample="a\0b"), bm.UsageError,
     "sample cannot hold a NUL byte"),
    (lambda out: bm.chain(MARKDUP, ["mark_duplicates", "mark_duplicates"], output=out),
     bm.UsageError, r"step 2 \(mark-duplicates\): it reads its input twice and cannot follow"),
    (lambda out: bm.chain(MARKDUP, ["insert_size", "insert_size"], output=out), bm.UsageError,
     "insert_size_metrics"),
    (lambda out: bm.chain(MARKDUP, [bm.step("write", output=out + ".sam")]), bm.UsageError,
     "W05"),
    (lambda out: bm.chain(MARKDUP, [bm.step("write", output=MARKDUP)]), bm.UsageError,
     "is an input; writing it would destroy it"),
    (lambda out: bm.chain(MARKDUP, ["hs_metrics"], output=out), bm.UsageError,
     r"step 1 \(hs-metrics\): missing --target-intervals \(Picard: TARGET_INTERVALS\)"),
    (lambda out: bm.chain(MARKDUP, ["insert_size", "gc_bias"], output=out), bm.UsageError,
     r"step 2 \(gc-bias\): missing --reference-sequence"),
    (lambda out: bm.chain(MARKDUP, [bm.step("mark_duplicates",
                                            duplicate_scoring_strategy="NOPE")], output=out),
     bm.UsageError, r"step 1 \(mark-duplicates\): DUPLICATE_SCORING_STRATEGY=NOPE: expected "),
    (lambda out: bm.chain(MARKDUP, [bm.step("mark_duplicates", as_="in")], output=out),
     bm.UsageError, "the stream name in is taken: it names the source"),
    (lambda out: bm.chain(MARKDUP, [bm.step("write", compression_level=11)], output=out),
     bm.UsageError, r'step 1 \(write\): expected 0 to 9, got "11"'),
    (lambda out: bm.chain(MARKDUP, [bm.step("insert_size", output="-"),
                                    bm.step("write", output="-")], output=out),
     bm.UsageError, "both write standard output; it takes one file"),
    (lambda out: bm.chain(DATA / "insert_size.cram", ["insert_size"], output=out), bm.UsageError,
     "is a CRAM: give the reference"),
    # An empty prefix would name hidden files (.insert_size_metrics) in the
    # working directory; an empty file name names no file.
    (lambda out: bm.chain(BAM, ["insert_size"], output=""), bm.UsageError,
     r"chain: the output prefix is empty \(from -o\)"),
    (lambda out: bm.chain(BAM, [bm.step("insert_size", output="")], output=out), bm.UsageError,
     r"step 1 \(insert-size\): --output names no file: its value is empty"),
    (lambda out: bm.chain(BAM, [bm.step("write", output="")], output=out), bm.UsageError,
     r"step 1 \(write\): --output names no file: its value is empty"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, also=["json", 1]), bm.UsageError,
     "chain: also takes a list of str, not one holding int"),
    (lambda out: bm.chain(BAM, ["insert_size"], output=out, sample=float("nan")),
     bm.UsageError, "chain: sample takes a finite number, not nan"),
])
def test_a_misused_chain_raises_before_anything_runs(tmp_path, call, error, message):
    with pytest.raises(error, match=message) as caught:
        call(str(tmp_path / "s"))
    assert names(tmp_path) == []
    assert getattr(caught.value, "run", None) is None


def test_a_failing_tap_fails_alone(tmp_path):
    with pytest.raises(bm.InputError, match=r"chain: 1 step\(s\) failed: step 1 "
                                            r"\(quality-yield\): .* should be coordinate "
                                            "sorted") as caught:
        bm.chain(QUERYNAME, [bm.step("quality_yield", assume_sorted=False), "mark_duplicates",
                             "insert_size", "write"], output=tmp_path / "s")
    run = caught.value.run
    assert isinstance(run, bm.Run) and run.exit_code == 1
    assert [s.status for s in run.stages] == ["failed", "ok", "ok", "ok"]
    assert list(run.errors) == ["step 1 (quality-yield)"]
    assert run.errors["step 1 (quality-yield)"] == run.stages[0].error
    assert run.stages[0].files == [] and list(run.metrics) == ["mark-duplicates", "insert-size"]
    assert names(tmp_path) == ["s.insert_size_metrics", "s.marked.bam",
                               "s.marked.duplication_metrics"]
    assert sorted(Path(f).name for f in run.files) == names(tmp_path)


def test_a_failing_transform_stops_its_segment(tmp_path):
    # insert_size.bam has no read group, which MarkDuplicates needs.
    with pytest.raises(OSError, match=r"chain: 1 step\(s\) failed: step 2 \(mark-duplicates\): "
                                      r".*; 1 stopped; 2 not run$") as caught:
        bm.chain(BAM, ["quality_yield", "mark_duplicates", "insert_size", "write"],
                 output=tmp_path / "s")
    run = caught.value.run
    assert [s.status for s in run.stages] == ["stopped", "failed", "not_run", "not_run"]
    assert run.stages[0].error.startswith("stopped: step 2 (mark-duplicates) failed: ")
    assert run.stages[2].error.startswith("not run: step 2 (mark-duplicates) failed: ")
    assert len(run.errors) == 4 and run.files == [] and run.metrics == {}
    assert names(tmp_path) == []


def round_trips(value):
    """copy, deepcopy and pickle of `value`. The pickled bytes are made here,
    from objects this test built: nothing untrusted is loaded."""
    return [copy.copy(value), copy.deepcopy(value), pickle.loads(pickle.dumps(value))]


STAGE_FIELDS = ("number", "verb", "role", "status", "error", "reads", "leaves", "segment",
                "files", "records_in", "records_out", "changed")


def test_steps_stages_and_chain_runs_copy_and_pickle(tmp_path):
    step = bm.step("mark_duplicates", remove_duplicates=True, comment=["a", "b"],
                   metrics_file=tmp_path / "m.txt")
    for other in round_trips(step):
        assert isinstance(other, bm.Step) and other == step and hash(other) == hash(step)
        assert (other.verb, other.role, other.args, other.options) == (
            step.verb, step.role, step.args, step.options)
    run = bm.chain(MARKDUP, ["quality_yield", step, "insert_size", "write"],
                   output=tmp_path / "s")
    for other in round_trips(run):
        assert isinstance(other, bm.Run)
        assert (other.sample, other.prefix, other.files, other.manifest, other.exit_code,
                other.errors, other.verdict) == (run.sample, run.prefix, run.files, None, 0,
                                                 {}, None)
        assert repr(other.metrics) == repr(run.metrics)  # NaN cells: repr, not ==
        assert [[getattr(s, f) for f in STAGE_FIELDS] for s in other.stages] == [
            [getattr(s, f) for f in STAGE_FIELDS] for s in run.stages]
    # A qc run keeps pickling, with no stage.
    qc = bm.qc(BAM, output=tmp_path / "q")
    assert pickle.loads(pickle.dumps(qc)).stages == []


def test_a_failed_chain_pickles_with_its_run(tmp_path):
    with pytest.raises(bm.InputError) as caught:
        bm.chain(BAM, ["mark_duplicates", "write"], output=tmp_path / "s")
    run = pickle.loads(pickle.dumps(caught.value.run))
    assert [s.status for s in run.stages] == ["failed", "not_run"]
    assert run.errors == caught.value.run.errors


def test_steps_and_a_chain_run_cross_a_process_pool(tmp_path):
    steps = [bm.step("mark_duplicates", remove_duplicates=True),
             bm.step("insert_size", minimum_pct=0.1), "write"]
    with concurrent.futures.ProcessPoolExecutor(1) as pool:
        run = pool.submit(bm.chain, MARKDUP, steps, output=tmp_path / "s").result()
    assert run.stages[0].records_out == 19
    assert (tmp_path / "s.marked.bam").exists()


def test_chain_releases_the_gil(tmp_path):
    runs = iter(range(5))
    inside = ticks_during(lambda: bm.chain(MARKDUP, markdup_steps(), reference=REFERENCE,
                                           output=tmp_path / f"s{next(runs)}", threads=1), 5)
    assert inside > 10, f"{inside} ticks while chain ran"


def test_chain_reads_standard_input_and_a_cram(tmp_path):
    cram = bm.chain(DATA / "insert_size.cram", ["insert_size"], reference=REFERENCE,
                    output=tmp_path / "c")
    plain = bm.chain(BAM, ["insert_size"], output=tmp_path / "b")
    assert cram.metrics["insert-size"]["classes"] == plain.metrics["insert-size"]["classes"]
    # `-` is the standard input of the process, as on the command line.
    code = ("import sys, bammetrics as bm\n"
            "run = bm.chain('-', ['mark_duplicates', 'insert_size', 'write'], output=sys.argv[1])\n"
            "print(run.sample, run.stages[0].records_in)\n")
    with open(MARKDUP, "rb") as source:
        done = subprocess.run([sys.executable, "-c", code, str(tmp_path / "p")], stdin=source,
                              capture_output=True, text=True)
    assert (done.returncode, done.stdout) == (0, "stdin 36\n"), done.stderr
    piped = bam(tmp_path / "p.marked.bam")
    bm.chain(MARKDUP, ["mark_duplicates", "write"], output=tmp_path / "f")
    assert piped[1] == bam(tmp_path / "f.marked.bam")[1]

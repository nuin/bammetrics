"""The Python mirror of the verbs (docs/05 section 5 and section 7 row 10):
qc, check, report, collect(metrics=[...]), tools and presets against the
command on the committed fixtures, the stub against the module both ways
(chain and step included; tests/python/test_chain.py runs them), the error
classes and the GIL release. Needs target/release/bammetrics
(`cargo build --release`) for the byte comparisons."""
import ast
import concurrent.futures
import copy
import inspect
import json
import os
import pickle
import shutil
import subprocess
import sys
import threading
import time
import types
from pathlib import Path

import pytest

import bammetrics as bm

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests" / "data"
BAM = DATA / "insert_size.bam"
REFERENCE = DATA / "ref_100k.fa"
PANEL = DATA / "targets.bed"
PICARD = DATA / "picard"
STUB = ROOT / "bammetrics.pyi"
CLI = ROOT / "target" / "release" / "bammetrics"
DATE = "2026-01-02 03:04:05 UTC"

needs_cli = pytest.mark.skipif(not CLI.exists(), reason="needs `cargo build --release`")


def cli(*args, cwd=None):
    return subprocess.run([str(CLI), *map(str, args)], cwd=cwd, capture_output=True, text=True)


def stub_tree():
    return ast.parse(STUB.read_text())


def stub_names():
    """Public top-level names the stub defines: functions, classes, variables."""
    names = set()
    for node in stub_tree().body:
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return {n for n in names if not n.startswith("_")}


def module_names():
    """Public names of the module; the maturin package shim adds the native
    submodule, which is no API."""
    return {n for n in dir(bm)
            if not n.startswith("_") and not isinstance(getattr(bm, n), types.ModuleType)}


def test_every_stub_name_exists_in_the_module():
    missing = stub_names() - module_names()
    assert not missing, f"in bammetrics.pyi, not in the module: {sorted(missing)}"


def test_every_module_name_is_in_the_stub():
    missing = module_names() - stub_names()
    assert not missing, f"in the module, not in bammetrics.pyi: {sorted(missing)}"


def test_all_lists_every_public_name():
    native = sys.modules.get("bammetrics.bammetrics", bm)
    assert set(native.__all__) == module_names()


@pytest.mark.parametrize("name", ["Run", "Verdict", "Check", "Step", "Stage"])
def test_stub_classes_match_the_runtime_classes_both_ways(name):
    node = next(n for n in stub_tree().body if isinstance(n, ast.ClassDef) and n.name == name)
    declared = {n.target.id for n in node.body if isinstance(n, ast.AnnAssign)}
    declared |= {n.name for n in node.body if isinstance(n, ast.FunctionDef)}
    declared = {n for n in declared if not n.startswith("_")}
    runtime = {n for n in dir(getattr(bm, name)) if not n.startswith("_")}
    assert declared == runtime


def test_error_classes_keep_value_error_and_os_error():
    assert issubclass(bm.UsageError, ValueError)
    assert issubclass(bm.InputError, OSError)
    for cls in (bm.UsageError, bm.InputError, bm.VerdictError):
        assert issubclass(cls, bm.BammetricsError)
        assert cls.__module__ == "bammetrics"


def test_old_functions_raise_the_new_classes(tmp_path):
    with pytest.raises(bm.InputError, match="nope.bam"):
        bm.insert_size(str(tmp_path / "nope.bam"))
    with pytest.raises(bm.UsageError, match="XX"):
        bm.alignment_summary(str(BAM), expected_pair_orientations=["XX"])
    with pytest.raises(bm.UsageError):
        bm.collect_insert_size_metrics(str(BAM), no_such_option=1)
    # Code written against 0.14 keeps catching them.
    try:
        bm.quality_yield(str(tmp_path / "nope.bam"))
    except OSError as e:
        assert "nope.bam" in str(e)


def test_qc_usage_and_input_errors(tmp_path):
    with pytest.raises(bm.UsageError, match="did you mean output"):
        bm.qc(BAM, outptu=tmp_path / "s")
    with pytest.raises(ValueError, match=r"invalid value \"maybe\" \(from --fail-on\)"):
        bm.qc(BAM, output=tmp_path / "s", fail_on="maybe")
    with pytest.raises(bm.UsageError, match="report takes True or False"):
        bm.qc(BAM, output=tmp_path / "s", report="no")
    with pytest.raises(OSError, match="nope.bam not found"):
        bm.qc(tmp_path / "nope.bam", output=tmp_path / "s")
    with pytest.raises(bm.InputError, match="absent.fa not found"):
        bm.qc(BAM, output=tmp_path / "s", reference=tmp_path / "absent.fa")
    bm.qc(BAM, output=tmp_path / "s")
    with pytest.raises(bm.UsageError, match="holds an earlier run"):
        bm.qc(BAM, output=tmp_path / "s")


def masked(path):
    """A written file with the Started on line (the run's clock) masked."""
    lines = path.read_bytes().split(b"\n")
    return b"\n".join(b"# Started on: <time>" if l.startswith(b"# Started on:") else l for l in lines)


@needs_cli
@pytest.mark.parametrize("extra", [
    [],
    ["--preset", "clinical-panel", "--panel", PANEL, "-R", REFERENCE],
])
def test_qc_writes_the_bytes_the_command_writes(tmp_path, monkeypatch, extra):
    a, b = tmp_path / "cli", tmp_path / "py"
    a.mkdir()
    b.mkdir()
    done = cli("qc", BAM, "-o", "s", "--date", DATE, *extra, cwd=a)
    assert done.returncode == 0, done.stderr
    keywords = {"output": "s", "date": DATE}
    flags = {"--preset": "preset", "--panel": "panel", "-R": "reference"}
    for flag, value in zip(extra[::2], extra[1::2]):
        keywords[flags[flag]] = value
    monkeypatch.chdir(b)
    run = bm.qc(BAM, **keywords)
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir())
    assert sorted(Path(f).name for f in run.files) == names
    for name in names:
        assert masked(a / name) == masked(b / name), name
    assert run.manifest == "s.manifest.json"
    assert run.exit_code == 0


@pytest.fixture(scope="module")
def failing_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("qc") / "s"
    return bm.qc(BAM, output=out, preset="clinical-panel", panel=PANEL, reference=REFERENCE,
                 date=DATE)


def test_a_failed_verdict_is_a_value(failing_run):
    verdict = failing_run.verdict
    assert verdict.status == "fail"
    statuses = {c.status for c in verdict.checks}
    assert statuses <= {"pass", "warn", "fail", "not_measured"}
    fail = next(c for c in verdict.checks if c.status == "fail")
    assert fail.limit == fail.fail and fail.why
    assert verdict.counts["fail"] == sum(c.status == "fail" for c in verdict.checks)
    assert json.loads(verdict.to_json()) == verdict.to_dict()
    assert verdict.to_dict() == json.loads(Path(failing_run.prefix + ".verdict.json").read_text())
    assert verdict.exit_code("fail") == 10


def test_run_metrics_are_parsed_files_by_verb(failing_run):
    assert {"quality-yield", "alignment-summary", "insert-size", "hs-metrics"} <= set(failing_run.metrics)
    insert = failing_run.metrics["insert-size"]
    assert insert == bm.read_metrics(failing_run.prefix + ".insert_size_metrics")
    assert insert["classes"][0]["rows"][0]["MEDIAN_INSERT_SIZE"] == 310


def test_fail_on_raises_with_the_run(tmp_path):
    with pytest.raises(bm.VerdictError, match="reaches fail_on") as caught:
        bm.qc(BAM, output=tmp_path / "s", preset="clinical-panel", panel=PANEL, fail_on="fail")
    assert caught.value.code == 10
    assert caught.value.run.verdict.status == "fail"
    assert [v.status for v in caught.value.verdicts] == ["fail"]
    assert Path(caught.value.run.prefix + ".report.html").exists()


def test_report_false_leaves_out_the_page(tmp_path):
    run = bm.qc(BAM, output=tmp_path / "s", report=False)
    assert not (tmp_path / "s.report.html").exists()
    assert (tmp_path / "s.verdict.json").exists()
    assert not any(f.endswith(".report.html") for f in run.files)


def test_repr_html_is_the_scoped_body_without_script(tmp_path):
    hostile = "<script>alert(1)</script>"
    run = bm.qc(BAM, output=tmp_path / "s", sample=hostile, preset="clinical-panel",
                panel=PANEL, all_plots=True)
    html = run._repr_html_()
    assert html.startswith('<div class="bm">\n<style>\n.bm {')
    assert "<script" not in html.lower()
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in html
    assert "<svg" in html
    assert run.sample == hostile
    page = (tmp_path / "s.report.html").read_text()
    assert "<script" not in page.lower()


@needs_cli
def test_check_on_picard_output_matches_the_command():
    for preset in ("generic", "clinical-panel"):
        done = cli("check", PICARD, "--preset", preset, "--format", "json")
        expected = json.loads(done.stdout)
        verdict = bm.check(PICARD, preset=preset, fail_on=None)
        assert isinstance(verdict, bm.Verdict)
        assert [verdict.to_dict()] == expected
        many = bm.check([PICARD], preset=preset, fail_on=None)
        assert [v.to_dict() for v in many] == expected
        if done.returncode == 10:
            with pytest.raises(bm.VerdictError) as caught:
                bm.check(PICARD, preset=preset)
            assert [v.sample for v in caught.value.verdicts] == [e["sample"] for e in expected]
        else:
            assert bm.check(PICARD, preset=preset).status == expected[0]["status"]


def test_check_reads_a_qc_prefix(failing_run):
    verdict = bm.check(failing_run.prefix, fail_on=None)
    assert verdict.to_dict() == failing_run.verdict.to_dict()


@needs_cli
def test_report_writes_the_files_the_command_writes(tmp_path):
    done = cli("report", PICARD, "-o", tmp_path / "cli", "--date", DATE)
    assert done.returncode == 0, done.stderr
    written = bm.report(PICARD, output=tmp_path / "py", date=DATE)
    names = sorted(p.name for p in (tmp_path / "cli").iterdir())
    assert sorted(Path(p).name for p in written) == names
    for name in names:
        assert (tmp_path / "cli" / name).read_bytes() == (tmp_path / "py" / name).read_bytes()
    with pytest.raises(bm.UsageError, match="several directories"):
        bm.report([PICARD, PICARD])


def test_collect_keeps_its_old_forms(tmp_path):
    positional = bm.collect(str(BAM), ["quality_yield", "insert_size"])
    assert sorted(positional) == ["insert_size", "quality_yield"]
    keyword = bm.collect(str(BAM), metrics=["quality_yield"], output_prefix=str(tmp_path / "o"),
                         threads=2)
    assert keyword["quality_yield"].total_reads > 0
    assert (tmp_path / "o.quality_yield_metrics").exists()
    assert sorted(bm.collect(str(BAM))) == ["alignment_summary", "insert_size", "quality_yield"]
    with pytest.raises(ValueError, match="wgs"):
        bm.collect(str(BAM), metrics=["wgs"])


def test_collect_with_verb_names_returns_a_run(tmp_path):
    run = bm.collect(BAM, metrics=["quality-yield", "insert-size", "gc-bias"],
                     reference=REFERENCE, output=tmp_path / "s", threads=2)
    assert isinstance(run, bm.Run)
    assert run.verdict is None and run.manifest is None
    assert set(run.metrics) == {"quality-yield", "insert-size", "gc-bias", "gc-bias.summary"}
    assert sorted(Path(f).name for f in run.files) == [
        "s.gc_bias.detail_metrics", "s.gc_bias.summary_metrics",
        "s.insert_size_metrics", "s.quality_yield_metrics"]
    assert run._repr_html_() is None
    # The same files qc writes for the same tools.
    qc = bm.qc(BAM, output=tmp_path / "q", preset="generic", report=False)
    assert masked(tmp_path / "s.insert_size_metrics").replace(b"/s.", b"/q.") == masked(
        tmp_path / "q.insert_size_metrics")
    assert qc.metrics["insert-size"]["classes"] == run.metrics["insert-size"]["classes"]
    with pytest.raises(bm.UsageError, match="needs reference"):
        bm.collect(BAM, metrics=["gc-bias"], output=tmp_path / "t")
    with pytest.raises(bm.UsageError, match="did you mean insert-size"):
        bm.collect(BAM, metrics=["insert-sise"], output=tmp_path / "t")
    with pytest.raises(bm.UsageError, match="output_prefix"):
        bm.collect(BAM, metrics=["insert-size"], output_prefix=tmp_path / "t")


@needs_cli
def test_tools_and_presets_are_the_info_verbs_data():
    table = [line.split() for line in cli("tools").stdout.split("\n\n")[0].splitlines()[1:]]
    tools = bm.tools()
    assert [[t.verb, t.picard_class] for t in tools] == table
    assert all(isinstance(t, str) and t == t.picard_class and t.summary for t in tools)
    printed = [line.split() for line in cli("presets").stdout.splitlines()]
    presets = bm.presets()
    assert presets == ["generic", "clinical-panel", "wgs", "rna-seq"]
    for row, preset in zip(printed, presets):
        assert row[:4] == [preset.name, f"v{preset.version}", preset.hash[:6], str(preset.rules)]


def test_qc_releases_the_gil(tmp_path):
    """A second Python thread runs while qc reads the BAM: with a switch
    interval of 1 s it could not, had qc kept the GIL."""
    runs = iter(range(5))
    inside = ticks_during(lambda: bm.qc(BAM, output=tmp_path / f"s{next(runs)}",
                                        preset="clinical-panel", panel=PANEL,
                                        reference=REFERENCE, threads=1), 5)
    assert inside > 10, f"{inside} ticks while qc ran"


def ticks_during(call, runs):
    """Ticks a second Python thread makes while `call` runs `runs` times, with
    a switch interval of 1 s: none had `call` kept the GIL."""
    ticks = []
    stop = threading.Event()

    def ticker():
        while not stop.is_set():
            ticks.append(time.perf_counter())
            time.sleep(0)

    interval = sys.getswitchinterval()
    sys.setswitchinterval(1.0)
    thread = threading.Thread(target=ticker)
    thread.start()
    try:
        time.sleep(0.01)
        inside = 0
        for _ in range(runs):
            start = time.perf_counter()
            call()
            end = time.perf_counter()
            inside += sum(start < t < end for t in ticks)
    finally:
        stop.set()
        thread.join()
        sys.setswitchinterval(interval)
    return inside


@pytest.mark.parametrize("name", ["insert_size", "alignment_summary", "quality_yield",
                                  "hs_metrics"])
def test_typed_functions_release_the_gil(name):
    targets = str(DATA / "targets.interval_list")
    calls = {
        "insert_size": lambda: bm.insert_size(str(BAM)),
        "alignment_summary": lambda: bm.alignment_summary(str(BAM)),
        "quality_yield": lambda: bm.quality_yield(str(BAM)),
        "hs_metrics": lambda: bm.hs_metrics(str(BAM), targets, targets),
    }
    inside = ticks_during(calls[name], 20)
    assert inside > 10, f"{inside} ticks while {name} ran"


def round_trips(value):
    """copy, deepcopy and pickle of `value`."""
    return [copy.copy(value), copy.deepcopy(value), pickle.loads(pickle.dumps(value))]


def test_tools_and_presets_rows_copy_and_pickle():
    for row in (bm.tools()[0], bm.presets()[1]):
        for other in round_trips(row):
            assert type(other) is type(row)
            assert other == row and vars(other) == vars(row)


CHECK_FIELDS = ("rule", "status", "value", "limit", "why", "label", "metric", "level", "unit",
                "warn", "fail")


def test_verdicts_checks_and_runs_copy_and_pickle(failing_run):
    verdict = failing_run.verdict
    for other in round_trips(verdict):
        assert isinstance(other, bm.Verdict) and other.to_json() == verdict.to_json()
    for check in verdict.checks:
        for other in round_trips(check):
            assert [getattr(other, f) for f in CHECK_FIELDS] == [
                getattr(check, f) for f in CHECK_FIELDS]
    for other in round_trips(failing_run):
        assert isinstance(other, bm.Run)
        assert (other.sample, other.prefix, other.files, other.manifest, other.exit_code) == (
            failing_run.sample, failing_run.prefix, failing_run.files, failing_run.manifest,
            failing_run.exit_code)
        assert repr(other.metrics) == repr(failing_run.metrics)  # NaN cells: repr, not ==
        assert other.verdict.to_json() == verdict.to_json()
        assert other._repr_html_() == failing_run._repr_html_()


def test_errors_copy_and_pickle(failing_run):
    error = bm.VerdictError("m", [failing_run.verdict], failing_run)
    for other in round_trips(error):
        assert str(other) == "m" and other.code == 10
        assert other.verdicts[0].to_json() == failing_run.verdict.to_json()
        assert other.run.prefix == failing_run.prefix
    for other in round_trips(bm.VerdictError("m", [])):
        assert other.verdicts == [] and other.run is None


def test_a_verdict_error_crosses_a_process_pool(failing_run):
    with concurrent.futures.ProcessPoolExecutor(1) as pool:
        future = pool.submit(bm.check, failing_run.prefix, preset="clinical-panel")
        with pytest.raises(bm.VerdictError) as caught:
            future.result()
    assert caught.value.verdicts[0].status == "fail"


def test_input_error_run_is_none_without_a_failed_collector(tmp_path):
    with pytest.raises(bm.InputError) as caught:
        bm.qc(tmp_path / "nope.bam", output=tmp_path / "s")
    assert caught.value.run is None


def test_path_keywords_take_a_str_or_a_path_only(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(TypeError, match="qc: output takes a str or os.PathLike, not int"):
        bm.qc(BAM, output=5)
    with pytest.raises(TypeError, match="qc: bam takes a str or os.PathLike, not NoneType"):
        bm.qc(None)
    with pytest.raises(TypeError, match="check: inputs takes"):
        bm.check(None)
    with pytest.raises(TypeError, match="qc: hs takes"):
        bm.qc(BAM, hs=(PANEL, 3))
    assert list(tmp_path.iterdir()) == []
    assert bm.qc(BAM, output="s", threads=2).exit_code == 0


def test_a_nul_byte_is_a_usage_error(tmp_path):
    with pytest.raises(bm.UsageError, match="cannot hold a NUL byte"):
        bm.qc(BAM, output=str(tmp_path / "s\0x"))
    with pytest.raises(bm.UsageError, match="cannot hold a NUL byte"):
        bm.qc(BAM, output=tmp_path / "s", sample="a\0b")


def test_an_unreadable_manifest_is_named_with_its_reason(tmp_path):
    (tmp_path / "s.manifest.json").write_text("{")
    with pytest.raises(bm.UsageError, match="is not a whole manifest: --resume or --force"):
        bm.qc(BAM, output=tmp_path / "s")
    assert bm.qc(BAM, output=tmp_path / "s", force=True).exit_code == 0


def test_report_false_keeps_a_page_no_run_wrote(tmp_path):
    page = tmp_path / "s.report.html"
    page.write_text("USER NOTES")
    bm.qc(BAM, output=tmp_path / "s", report=False)
    assert page.read_text() == "USER NOTES"
    bm.qc(BAM, output=tmp_path / "t")
    assert (tmp_path / "t.report.html").exists()
    bm.qc(BAM, output=tmp_path / "t", report=False, force=True)
    assert not (tmp_path / "t.report.html").exists()


def test_collect_refuses_the_prefix_of_a_qc_run(tmp_path):
    bm.qc(BAM, output=tmp_path / "s")
    with pytest.raises(bm.UsageError, match="holds a qc run"):
        bm.collect(BAM, metrics=["insert-size"], output=tmp_path / "s")


def test_old_collect_type_errors_name_the_argument():
    with pytest.raises(TypeError, match="argument 'path'"):
        bm.collect(5)
    with pytest.raises(TypeError, match="argument 'hs'"):
        bm.collect(str(BAM), metrics=["hs_metrics"], hs="x")


def test_check_returns_an_incomplete_verdict_under_fail_on(tmp_path):
    shutil.copy(PICARD / "NA12878" / "NA12878.insert_size_metrics",
                tmp_path / "NA12878.insert_size_metrics")
    for fail_on in ("fail", "warn", None):
        verdict = bm.check(tmp_path, preset="clinical-panel", fail_on=fail_on)
        assert verdict.status == "incomplete"
    assert verdict.exit_code("fail") == 11


def stub_function(name):
    return next(n for n in stub_tree().body if isinstance(n, ast.FunctionDef) and n.name == name)


@pytest.mark.parametrize("name", ["qc", "check", "report", "chain", "step"])
def test_help_shows_the_stub_signature_and_docstring(name):
    # The first entry of an overloaded function (chain) carries the docstring
    # and the defaults help() shows.
    node = stub_function(name)
    function = getattr(bm, name)
    assert " ".join(function.__doc__.split()) == " ".join(ast.get_docstring(node).split())
    args = node.args
    empty = inspect.Parameter.empty
    expected = [(a.arg, inspect.Parameter.POSITIONAL_ONLY, empty) for a in args.posonlyargs]
    defaults = [empty] * (len(args.args) - len(args.defaults)) + [
        ast.literal_eval(d) for d in args.defaults]
    expected += [(a.arg, inspect.Parameter.POSITIONAL_OR_KEYWORD, d)
                 for a, d in zip(args.args, defaults)]
    expected += [(a.arg, inspect.Parameter.KEYWORD_ONLY, ast.literal_eval(d))
                 for a, d in zip(args.kwonlyargs, args.kw_defaults)]
    if args.kwarg:
        expected.append((args.kwarg.arg, inspect.Parameter.VAR_KEYWORD, empty))
    parameters = inspect.signature(function).parameters.values()
    assert [(p.name, p.kind, p.default) for p in parameters] == expected


def test_stub_collector_verbs_are_the_runtime_names():
    node = next(n for n in stub_tree().body if isinstance(n, ast.Assign)
                and n.targets[0].id == "_CollectorVerb")
    stub = [e.value for e in node.value.slice.elts]
    with pytest.raises(bm.UsageError) as caught:
        bm.collect(BAM, metrics=["no-such"])
    runtime = str(caught.value).split("one of: ")[1].split(", ")
    assert stub == [n for n in runtime if "-" in n]


def test_stub_step_verbs_are_the_runtime_names():
    node = next(n for n in stub_tree().body if isinstance(n, ast.Assign)
                and n.targets[0].id == "_StepVerb")
    stub = [e.value for e in node.value.slice.elts]
    with pytest.raises(bm.UsageError) as caught:
        bm.step("no_such")
    assert stub == str(caught.value).split("steps: ")[1].split(", ")
    assert stub[-1] == "write"
    for name in stub:
        assert bm.step(name).verb == name.replace("_", "-")


def test_stub_chain_overloads_differ_in_dry_run_only():
    overloads = [n for n in stub_tree().body
                 if isinstance(n, ast.FunctionDef) and n.name == "chain"]
    returns = [ast.unparse(n.returns) for n in overloads]
    assert returns == ["Run[None]", "str", "Run[None] | str"]
    keywords = [[a.arg for a in n.args.kwonlyargs] for n in overloads]
    assert keywords[0] == keywords[1] == keywords[2] and keywords[0][-1] == "dry_run"
    dry_run = [ast.unparse(n.args.kwonlyargs[-1].annotation) for n in overloads]
    assert dry_run == ["Literal[False]", "Literal[True]", "bool"]

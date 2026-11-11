"""The Picard queue (tests/concordance/picard_lock.py): one JVM at a time,
JAVA_OPTS and the exit code passed through, and every harness call routed
through it. A shell script stands in for `picard`."""
import ast
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pytest

import picard_lock
from conftest import ROOT
from run_picard import picard_env

# Writes one `start` and one `end` line around a pause, then echoes its input.
SLOW_PICARD = '''#!/bin/sh
echo "start $$" >> "$FAKE_LOG"
sleep "${FAKE_SLEEP:-0}"
echo "end $$" >> "$FAKE_LOG"
echo "opts=$JAVA_OPTS args=$*"
exit "${FAKE_EXIT:-0}"
'''


@pytest.fixture
def slow_picard(tmp_path, monkeypatch, queue_env):
    """An environment whose `picard` is the queue shim in front of SLOW_PICARD."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    picard_lock.write_script(bin_dir / "picard", SLOW_PICARD)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("FAKE_LOG", str(tmp_path / "log"))
    return tmp_path / "log"


def test_java_opts_arguments_and_exit_code_pass_through(slow_picard, monkeypatch):
    monkeypatch.setenv("JAVA_OPTS", "-Xmx1g")
    monkeypatch.setenv("FAKE_EXIT", "3")
    done = subprocess.run(["picard", "MarkDuplicates", "I=a b.bam", "--help"], env=picard_env(),
                          capture_output=True, text=True)
    assert done.stdout == "opts=-Xmx1g args=MarkDuplicates I=a b.bam --help\n"
    assert done.returncode == 3 and done.stderr == ""


def test_picard_env_puts_the_shim_first_and_keeps_the_real_command(slow_picard):
    env = picard_env()
    first = env["PATH"].split(os.pathsep)[0]
    assert picard_lock.is_shim(Path(first) / "picard")
    assert picard_lock.real_picard(env) == str(slow_picard.parent / "bin" / "picard")


def test_concurrent_calls_never_overlap(slow_picard, monkeypatch):
    monkeypatch.setenv("FAKE_SLEEP", "0.3")
    env = picard_env()
    running = [subprocess.Popen(["picard", "X"], env=env, stdout=subprocess.DEVNULL) for _ in range(4)]
    assert [p.wait() for p in running] == [0, 0, 0, 0]
    events = [line.split() for line in slow_picard.read_text().splitlines()]
    assert len(events) == 8
    for start, end in zip(events[0::2], events[1::2]):
        assert (start[0], end[0]) == ("start", "end") and start[1] == end[1], events


def test_without_the_queue_the_same_calls_overlap(slow_picard, monkeypatch):
    """The control for the test above: the stand-in started directly runs four at once."""
    monkeypatch.setenv("FAKE_SLEEP", "1")
    running = [subprocess.Popen(["picard", "X"], stdout=subprocess.DEVNULL) for _ in range(4)]
    assert [p.wait() for p in running] == [0, 0, 0, 0]
    words = [line.split()[0] for line in slow_picard.read_text().splitlines()]
    assert words[:2] == ["start", "start"], words


def test_the_lock_lasts_as_long_as_the_process(slow_picard, queue_env, monkeypatch):
    monkeypatch.setenv("FAKE_SLEEP", "1")
    assert picard_lock.status().startswith("free")
    child = subprocess.Popen(["picard", "CollectWgsMetrics", "I=x.bam"], env=picard_env(),
                             stdout=subprocess.DEVNULL)
    deadline = time.monotonic() + 10
    while not slow_picard.exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    held = picard_lock.status()
    assert held.startswith("held: pid ") and "picard CollectWgsMetrics I=x.bam" in held
    assert child.wait() == 0
    assert picard_lock.status().startswith("free: last holder pid ")


def test_timeout_exits_75_and_names_the_holder(slow_picard, queue_env, monkeypatch):
    fd = picard_lock.acquire(queue_env)
    picard_lock.note_holder(fd, ["picard", "LongRun"])
    monkeypatch.setenv("PICARD_LOCK_TIMEOUT", "0.2")
    done = subprocess.run(["picard", "X"], env=picard_env(), capture_output=True, text=True)
    os.close(fd)
    assert done.returncode == picard_lock.TIMED_OUT
    assert "picard LongRun" in done.stderr and not slow_picard.exists()


def test_held_lets_its_own_children_through_and_blocks_the_rest(slow_picard, queue_env, monkeypatch):
    monkeypatch.setenv("PICARD_LOCK_TIMEOUT", "0.2")
    with picard_lock.held(["picard", "bench"]):
        inside = subprocess.run(["picard", "X"], env=picard_env(), capture_output=True)
        outsider = {k: v for k, v in picard_env().items() if k != picard_lock.HELD}
        blocked = subprocess.run(["picard", "X"], env=outsider, capture_output=True)
    after = subprocess.run(["picard", "X"], env=picard_env(), capture_output=True)
    assert (inside.returncode, blocked.returncode, after.returncode) == (0, picard_lock.TIMED_OUT, 0)
    assert picard_lock.HELD not in os.environ


def test_lock_file_is_private_and_not_followed_through_a_link(tmp_path):
    path = tmp_path / "cache" / "picard.lock"
    os.close(picard_lock.acquire(path))
    assert path.stat().st_mode & 0o777 == 0o600
    assert path.parent.stat().st_mode & 0o777 == 0o700
    link = tmp_path / "link.lock"
    link.symlink_to(tmp_path / "elsewhere")
    with pytest.raises(OSError):
        picard_lock.acquire(link)
    assert not (tmp_path / "elsewhere").exists()


PS = """\
  501 /opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home/bin/java -Xmx8g -jar /opt/picard/picard.jar MarkDuplicates I=a.bam
  502 java -jar /home/runner/picard/picard.jar CollectWgsMetrics
  503 /usr/bin/python3 tests/concordance/picard_lock.py --status
  504 vim /opt/picard/picard.jar.notes
  505 java -jar gatk.jar HaplotypeCaller
  506 /bin/sh /opt/homebrew/bin/picard MarkDuplicates
"""


def test_picard_jvms_are_java_processes_running_picard_jar():
    assert picard_lock.picard_jvms(PS.splitlines()) == [501, 502]


def test_wait_for_foreign_polls_until_no_jvm_is_listed():
    tables = [PS.splitlines(), PS.splitlines(), ["503 python3 x"]]
    naps = []
    picard_lock.wait_for_foreign(table=lambda: tables.pop(0), sleep=naps.append)
    assert naps == [picard_lock.FOREIGN_POLL] * 2 and not tables
    with pytest.raises(picard_lock.Timeout, match="pid 501"):
        picard_lock.wait_for_foreign(deadline=time.monotonic() - 1, table=PS.splitlines,
                                     sleep=naps.append)


def test_each_interpreter_has_a_shim_that_is_never_rewritten(tmp_path):
    """Two harness processes started from different interpreters must not
    replace each other's script: a queued call would start an interpreter
    its own process never named."""
    script = picard_lock.shim(tmp_path / "queue") / "picard"
    assert os.access(script, os.X_OK) and picard_lock.is_shim(script)
    assert sys.executable in script.read_text()
    before = (script.stat().st_ino, script.read_text())
    other = picard_lock.shim(tmp_path / "queue", python="/another/python") / "picard"
    assert other != script and "/another/python" in other.read_text()
    assert picard_lock.shim(tmp_path / "queue") / "picard" == script
    assert (script.stat().st_ino, script.read_text()) == before
    homes = sorted(p.name for p in (tmp_path / "queue").iterdir())
    assert homes == sorted([script.parent.name, other.parent.name])
    assert all(re.fullmatch("[0-9a-f]{12}", name) for name in homes)
    assert [p.name for p in script.parent.iterdir()] == ["picard"]


@pytest.mark.parametrize("value", ["abc", "nan", "-1", "inf", "1s"])
def test_a_timeout_that_is_not_seconds_is_a_usage_error(slow_picard, monkeypatch, value):
    monkeypatch.setenv("PICARD_LOCK_TIMEOUT", value)
    done = subprocess.run(["picard", "X"], env=picard_env(), capture_output=True, text=True)
    assert done.returncode == picard_lock.USAGE == 2
    assert done.stderr == f"picard_lock: PICARD_LOCK_TIMEOUT must be seconds, not {value!r}\n"
    assert not slow_picard.exists()
    with pytest.raises(ValueError, match="PICARD_LOCK_TIMEOUT must be seconds"):
        with picard_lock.held():
            pass


def test_timeout_values():
    assert picard_lock.timeout({}) is None and picard_lock.timeout({"PICARD_LOCK_TIMEOUT": ""}) is None
    assert picard_lock.timeout({"PICARD_LOCK_TIMEOUT": "0"}) == 0
    assert picard_lock.timeout({"PICARD_LOCK_TIMEOUT": "1.5"}) == 1.5


# --- every harness call goes through the queue --------------------------------------

HARNESS = [*sorted((ROOT / "tests").rglob("*.py")), *sorted((ROOT / "benches").glob("*.py")),
           *sorted((ROOT / "tools").glob("*.py"))]
STARTERS = {"run", "Popen", "call", "check_call", "check_output"}
# `run` of these modules gives a `picard` command the queue environment itself.
QUEUED_RUN = re.compile(r"from (asm_bams|make_corpus) import [^)]*\brun\b")


# The output formats of `--format` and `--also` (docs/formats.md).
FORMATS = {"picard", "json", "jsonl", "tsv", "csv", "human", "markdown", "junit", "multiqc"}


def format_list(node):
    """A list of two or more format names, such as `format=["picard", "json"]`
    of the Python `chain`: the names of files to write, not a command."""
    return len(node.elts) > 1 and all(
        isinstance(item, ast.Constant) and item.value in FORMATS for item in node.elts)


def starts_picard(node):
    """A list literal whose first element is the string "picard", unless it
    is a `format_list`."""
    return (isinstance(node, ast.List) and node.elts and isinstance(node.elts[0], ast.Constant)
            and node.elts[0].value == "picard" and not format_list(node))


def unqueued_calls(path, name):
    """Problems in one file, each starting with `name`. A literal
    `["picard", ...]` handed to a `subprocess` call needs `env=`. One handed
    to a bare `run(...)` needs that `run` to come from a QUEUED_RUN module.
    Any other use needs the file to name `picard_env` or `picard_lock`."""
    text = path.read_text()
    tree = ast.parse(text)
    names_queue = "picard_env" in text or "picard_lock" in text
    handed = set()
    problems = []
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and node.args and starts_picard(node.args[0])):
            continue
        handed.add(id(node.args[0]))
        if isinstance(node.func, ast.Attribute) and node.func.attr in STARTERS:
            if not any(k.arg == "env" for k in node.keywords):
                problems.append(f"{name}:{node.lineno}: picard started without env=")
        elif not (names_queue or (isinstance(node.func, ast.Name) and node.func.id == "run"
                                  and QUEUED_RUN.search(text))):
            problems.append(f"{name}:{node.lineno}: picard command handed to an unknown runner")
    loose = [n for n in ast.walk(tree) if starts_picard(n) and id(n) not in handed]
    if loose and not names_queue:
        problems.append(f"{name}:{loose[0].lineno}: builds a picard command without picard_env()")
    return problems


def test_every_picard_call_of_the_harness_takes_the_queue_environment():
    assert len(HARNESS) > 40
    problems = [problem for path in HARNESS if path.parent.name != "harness"
                for problem in unqueued_calls(path, path.relative_to(ROOT))]
    assert not problems, "\n".join(problems)


def test_the_check_catches_a_direct_call(tmp_path):
    new = tmp_path / "test_new_tool.py"
    call = 'subprocess.run(["picard", "SortSam", "I=a.bam"]'
    for source, expected in (
            (f"import subprocess\n{call}, check=True)\n", ["new:2: picard started without env="]),
            (f"import subprocess\nfrom run_picard import picard_env\n{call}, env=picard_env())\n", []),
            ('cmd = ["picard", "SortSam"]\n', ["new:1: builds a picard command without picard_env()"]),
            ('launch(["picard", "SortSam"])\n', ["new:1: picard command handed to an unknown runner"]),
            # A `--format` list names output formats and starts no JVM.
            ('chain(bam, ["insert_size"], format=["picard", "json"])\n', []),
            ('plan = {"format": ["picard", "json", "csv"]}\n', []),
            ('cmd = ["picard"]\n', ["new:1: builds a picard command without picard_env()"]),
            ('cmd = ["picard", "json", "SortSam"]\n',
             ["new:1: builds a picard command without picard_env()"]),
            ('from asm_bams import (Sam,\n    run)\nrun(["picard", "SortSam"])\n', [])):
        new.write_text(source)
        assert unqueued_calls(new, "new") == expected, source

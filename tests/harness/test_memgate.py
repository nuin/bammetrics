"""benches/memgate.py on stand-in commands: the table reader, the two `time` report
formats, a row within its limits, a row over each limit, a failing command, a missing
input in each tier, and that a row's directory is gone after the measurement. No
release binary and no real data."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "benches"))

import memgate  # noqa: E402

MACOS = """        1.25 real         0.90 user         0.10 sys
            52428800  maximum resident set size
                   0  average shared memory size
"""
GNU = """\tCommand being timed: "x"
\tElapsed (wall clock) time (h:mm:ss or m:ss): 1:02.50
\tMaximum resident set size (kbytes): 51200
"""
# A command that holds `mib` MiB, optionally writes `spill` MiB into `{out}/tmp`.
HOLD = ("import sys, os, time; mib, spill, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; "
        "held = b'x' * (mib << 20); "
        "os.makedirs(out + '/tmp', exist_ok=True); "
        "open(out + '/tmp/run', 'wb').write(b'x' * (spill << 20)); "
        "time.sleep(0.5 if spill else 0); "
        "open(out + '/product', 'w').write('done')")


def row(name="hold", mib=1, spill=0, **keys):
    argv = [sys.executable, "-c", HOLD, str(mib), str(spill), "{out}"]
    return memgate.Row(name=name, argv=argv, **{"max_rss_mib": 200, **keys})


def test_both_time_reports_are_read():
    assert memgate.parse_time(MACOS, 1) == (1.25, 52428800)
    assert memgate.parse_time(GNU, 1024) == (62.5, 51200 * 1024)
    with pytest.raises(ValueError):
        memgate.parse_time("no report here", 1)


def test_the_committed_table_parses_and_names_both_tiers():
    rows = memgate.parse_table(memgate.TABLE.read_text())
    assert {r.tier for r in rows} == {"fixture", "scale"}
    fixture = [r.name for r in rows if r.tier == "fixture"]
    assert fixture == ["build-bam-index", "mark-duplicates", "chain write --index"]
    for r in rows:
        assert r.argv[0] == "{bin}" and r.max_rss_mib > 0
        if r.tier == "fixture":
            assert all((ROOT / need).is_file() for need in r.needs), r.name
        else:
            assert r.needs and all(need.startswith("~/data/") for need in r.needs), r.name


@pytest.mark.parametrize("text, message", [
    ("", "no [[command]]"),
    ("[[command]]\nargv = ['x']\nmax_rss_mib = 1\n", "no name"),
    ("[[command]]\nname = 'a'\nargv = ['x']\n", "a: max_rss_mib must be a number above 0"),
    ("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 0\n", "max_rss_mib must be a number above 0"),
    ("[[command]]\nname = 'a'\nargv = []\nmax_rss_mib = 1\n", "a: argv must be a list of strings"),
    ("[[command]]\nname = 'a'\nargv = 'x'\nmax_rss_mib = 1\n", "a: argv must be a list of strings"),
    ("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 1\nrss = 2\n", "a: unknown key 'rss'"),
    ("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 1\ntier = 'big'\n", "a: tier must be one of"),
    ("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 1\nmax_spill_mib = 1\n", "needs spill_dir"),
    ("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 1\n" * 2, "a: the name is used twice"),
    ("[[commands]]\nname = 'a'\n", "unknown table 'commands'"),
    ("name = ", "not TOML"),
])
def test_a_bad_table_is_refused_with_the_row_and_the_key(text, message):
    with pytest.raises(memgate.TableError, match=message.replace("[", r"\[")):
        memgate.parse_table(text)


def test_a_row_within_its_limits_passes_and_leaves_no_file(tmp_path):
    measured = memgate.measure(row(mib=20, max_seconds=60), Path("unused"), tmp_path)
    assert measured.status == "ok", measured.note
    assert 20 <= measured.rss_mib < 200 and 0 <= measured.seconds < 60
    assert measured.spill_mib is None
    assert list(tmp_path.iterdir()) == []


def test_a_row_over_its_rss_limit_fails(tmp_path):
    measured = memgate.measure(row(mib=120, max_rss_mib=60), Path("unused"), tmp_path)
    assert measured.status == "over" and "MiB over 60" in measured.note
    assert measured.rss_mib >= 120
    assert list(tmp_path.iterdir()) == []


def test_a_row_over_its_time_or_spill_limit_fails(tmp_path):
    slow = memgate.measure(row(spill=1, max_seconds=0.1, spill_dir="{out}/tmp"), Path("unused"), tmp_path)
    assert slow.status == "over" and "s over 0" in slow.note
    spilled = memgate.measure(row(spill=8, spill_dir="{out}/tmp", max_spill_mib=4), Path("unused"), tmp_path)
    assert spilled.status == "over" and "spill 8 MiB over 4" in spilled.note
    within = memgate.measure(row(spill=8, spill_dir="{out}/tmp", max_spill_mib=16), Path("unused"), tmp_path)
    assert within.status == "ok" and 8 <= within.spill_mib < 16


def test_a_failing_command_fails_the_row_with_its_first_line(tmp_path):
    bad = memgate.Row(name="bad", argv=[sys.executable, "-c", "import sys; sys.exit('went wrong')"],
                      max_rss_mib=100)
    measured = memgate.measure(bad, Path("unused"), tmp_path)
    assert measured.status == "failed" and "exit 1: went wrong" in measured.note


def test_a_missing_input_fails_a_fixture_row_and_skips_a_scale_row(tmp_path):
    missing = str(tmp_path / "none.bam")
    assert memgate.measure(row(needs=[missing]), Path("unused"), tmp_path).status == "failed"
    scale = row(needs=[missing], tier="scale")
    skipped = memgate.measure(scale, Path("unused"), tmp_path)
    assert skipped.status == "skipped" and "missing" in skipped.note
    assert memgate.measure(scale, Path("unused"), tmp_path, strict=True).status == "failed"


def test_the_report_names_each_row_its_numbers_and_its_verdict(tmp_path):
    results = [memgate.measure(row(mib=10, max_seconds=60), Path("unused"), tmp_path),
               memgate.measure(row(name="big", mib=120, max_rss_mib=60), Path("unused"), tmp_path)]
    text = memgate.render(results)
    lines = text.splitlines()
    assert lines[0].split()[:2] == ["command", "tier"]
    assert lines[1].startswith("hold") and lines[1].rstrip().endswith("ok")
    assert "over: RSS" in lines[2]
    markdown = memgate.render(results, markdown=True)
    assert markdown.splitlines()[1].startswith("|---|") and "| big | fixture |" in markdown


def test_main_lists_filters_and_reports_usage_errors(tmp_path, capsys):
    table = tmp_path / "table.toml"
    table.write_text("[[command]]\nname = 'a'\nargv = ['x']\nmax_rss_mib = 1\n"
                     "[[command]]\nname = 'b'\ntier = 'scale'\nargv = ['y']\nmax_rss_mib = 2\n")
    assert memgate.main(["--table", str(table), "--list"]) == 0
    assert capsys.readouterr().out == "a\tfixture\t1 MiB\tx\n"
    assert memgate.main(["--table", str(table), "--list", "--tier", "all"]) == 0
    assert capsys.readouterr().out.count("\n") == 2
    assert memgate.main(["--table", str(table), "--list", "--only", "b"]) == 0
    assert capsys.readouterr().out.startswith("b\tscale")
    assert memgate.main(["--table", str(table), "--only", "c"]) == 2
    assert "no row named 'c'" in capsys.readouterr().err
    assert memgate.main(["--table", str(tmp_path / "none.toml")]) == 2
    table.write_text("[[command]]\nname = 'a'\n")
    assert memgate.main(["--table", str(table)]) == 2
    assert "a: argv must be a list of strings" in capsys.readouterr().err

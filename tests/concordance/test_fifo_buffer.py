"""FifoBuffer against Picard 3.5.0: the bytes on standard output, the exit code and the
log lines of DEBUG_FREQUENCY without their timestamp.

Six input shapes go through both tools under every option: no byte, one byte, FASTQ text
without a final newline, every byte value, 3 MiB of seeded random bytes (larger than the
small BUFFER_SIZE values) and the foundation corpus BAM. Options: BUFFER_SIZE (1, 1000,
the default), IO_SIZE (1, 7, one larger than the buffer), DEBUG_FREQUENCY (-1, 0, 1),
NAME (a value, `null`), VERBOSITY (each level), the dashed spellings and the common
options Picard accepts. The `fifo-buffer` verb must write the bytes Picard writes.

Log lines: Picard starts its log thread beside the copy, so its first line holds whatever
arrived by then. Two shapes make every compared line the same in both tools: an input
that arrives after the second line (`used 0` twice) and a standard output nobody reads
(the buffer is full from the second line on). bammetrics prints its first line before
the first read, so that line always reads `used 0` (KNOWN_DIFFERENCES.md).

Error parity: a value that is no Java int, a negative size, `null`, an option given
twice, INPUT, OUTPUT, ASSUME_SORTED, GA4GH_CLIENT_SECRETS, an unknown option and a bad
VERBOSITY exit non-zero in both tools with nothing on standard output. A size of 0 never
ends in Picard (checked with a timeout) and exits 1 here. A reader that closes standard
output ends both tools with exit 0; Picard reads the rest of its input, bammetrics stops.

Picard's own QUIET default for this tool is true, so its standard error carries the JVM's
library warnings, the syntax banner and the tool's log lines; `LOG` picks the log lines.
JAVA_OPTS=-Xmx8g. Needs picard; run `make concordance`.
"""
import random
import re
import subprocess
import threading

import pytest

from cases import BAMMETRICS
from conftest import RESULTS
from run_picard import picard_env

TOOL = "FifoBuffer"
VERB = "fifo-buffer"
# `(.+)`: Picard prints an empty INFO line (its banner) before the tool's lines.
LOG = re.compile(r"^(INFO|WARNING|ERROR)\t\S+ \S+\tFifoBuffer\t(.+)$")
OURS_LOG = re.compile(r"^(INFO|WARNING|ERROR)\tFifoBuffer\t(.+)$")
VERB_LOG = re.compile(r"^bammetrics fifo-buffer: (Fifo buffer .+)$")
CHAIN_NOTE = ("bammetrics fifo-buffer: bammetrics chain runs its steps in one process and "
              "removes the pipe this buffer serves (-q drops this line)")
MIB = 1024 * 1024

# Each option set in Picard syntax and in the verb's.
OPTION_SETS = [
    ([], []),
    (["BUFFER_SIZE=1000"], ["--buffer-size", "1000"]),
    (["BUFFER_SIZE=1000", "IO_SIZE=7"], ["--buffer-size=1000", "--io-size=7"]),
    (["BUFFER_SIZE=10", "IO_SIZE=4096"], ["BUFFER_SIZE=10", "IO_SIZE=4096"]),
    (["IO_SIZE=100000"], ["--io-size", "100000"]),
    (["--BUFFER_SIZE", "65536", "--IO_SIZE", "65536"], ["--BUFFER_SIZE", "65536", "--IO_SIZE", "65536"]),
    (["DEBUG_FREQUENCY=0", "NAME=lane1"], ["--debug-frequency", "0", "--name", "lane1"]),
    (["DEBUG_FREQUENCY=-1", "NAME=null"], ["--debug-frequency=-1", "--name", "null"]),
    (["DEBUG_FREQUENCY=1", "VERBOSITY=WARNING"], ["--debug-frequency", "1", "--verbosity", "WARNING"]),
    (["DEBUG_FREQUENCY=1", "VERBOSITY=ERROR", "NAME=x"],
     ["--debug-frequency", "1", "--verbosity", "ERROR", "--name", "x"]),
    (["VERBOSITY=DEBUG"], ["--verbosity", "DEBUG"]),
    (["VERBOSITY=INFO", "VALIDATION_STRINGENCY=SILENT", "COMPRESSION_LEVEL=1",
      "MAX_RECORDS_IN_RAM=5", "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "TMP_DIR=/tmp",
      "USE_JDK_DEFLATER=true", "R=no_such.fa"],
     ["--verbosity", "INFO", "--validation-stringency", "SILENT", "--compression-level", "1",
      "--max-records-in-ram", "5", "--create-index", "--create-md5-file", "--tmp-dir", "/tmp",
      "--use-jdk-deflater", "-R", "no_such.fa"]),
]
# One byte per read and per buffer: two thread switches a byte, so the small shapes only.
ONE_BYTE = (["BUFFER_SIZE=1", "IO_SIZE=1"], ["--buffer-size", "1", "--io-size", "1"])
SMALL = 5000

# What each bad argument list must say in both tools; None where Picard throws a Java
# exception (NegativeArraySizeException, IllegalArgumentException) bammetrics rewords.
BAD = [
    (["BUFFER_SIZE=abc"], "Problem constructing Integer from the string 'abc'."),
    (["BUFFER_SIZE=1.5"], "Problem constructing Integer from the string '1.5'."),
    (["BUFFER_SIZE=3000000000"], "Problem constructing Integer from the string '3000000000'."),
    (["BUFFER_SIZE=-1"], None),
    (["BUFFER_SIZE=null"], None),
    (["BUFFER_SIZE=2147483647"], None),
    (["IO_SIZE=abc"], "Problem constructing Integer from the string 'abc'."),
    (["IO_SIZE=-1"], None),
    (["IO_SIZE=null"], None),
    (["DEBUG_FREQUENCY=abc"], "Problem constructing Integer from the string 'abc'."),
    (["DEBUG_FREQUENCY=null"], None),
    (["VERBOSITY=LOUD"], None),
    (["BUFFER_SIZE=10", "BUFFER_SIZE=20"], "Option 'BUFFER_SIZE' cannot be specified more than once."),
    (["INPUT=x"], "Unrecognized option: INPUT"),
    (["I=x"], "Unrecognized option: I"),
    (["OUTPUT=x"], "Unrecognized option: OUTPUT"),
    (["O=x"], "Unrecognized option: O"),
    (["ASSUME_SORTED=true"], "Unrecognized option: ASSUME_SORTED"),
    (["GA4GH_CLIENT_SECRETS=x"], "Unrecognized option: GA4GH_CLIENT_SECRETS"),
    (["STOP_AFTER=1"], None),
    (["NO_SUCH_OPTION=1"], None),
]


@pytest.fixture(scope="module")
def shapes(corpus):
    """The six input shapes by name."""
    rng = random.Random(20261110)
    return {
        "empty": b"",
        "one_byte": b"\x00",
        "text": (b"@r1\nACGT\n+\nIIII\n" * 200)[:-1],
        "every_byte_value": bytes(range(256)) * 16,
        "random_3mib": rng.randbytes(3 * MIB + 5),
        "bam": (corpus / "corpus.bam").read_bytes(),
    }


def command(side, args):
    """The argv and environment of one tool: `picard`, `ours` (the Picard face) or `verb`."""
    if side == "picard":
        env = picard_env()
        env["JAVA_OPTS"] = "-Xmx8g"
        return ["picard", TOOL, *args], env
    return [str(BAMMETRICS), VERB if side == "verb" else TOOL, *args], None


def run(side, args, data, timeout=300):
    argv, env = command(side, args)
    return subprocess.run(argv, input=data, capture_output=True, env=env, timeout=timeout)


def log_lines(side, stderr):
    pattern = {"picard": LOG, "ours": OURS_LOG, "verb": VERB_LOG}[side]
    lines = stderr.decode(errors="replace").splitlines()
    return [m.groups()[-1] for m in map(pattern.match, lines) if m]


def test_bytes_and_exit_code_equal_for_every_shape_and_option(shapes):
    assert len(shapes) == 6
    runs = 0
    for name, data in shapes.items():
        option_sets = OPTION_SETS + ([ONE_BYTE] if len(data) <= SMALL else [])
        for picard_args, verb_args in option_sets:
            label = f"{name} {picard_args}"
            picard = run("picard", picard_args, data)
            ours = run("ours", picard_args, data)
            verb = run("verb", verb_args, data)
            assert picard.returncode == 0, f"{label}: {picard.stderr[-400:]!r}"
            assert picard.stdout == data, f"{label}: Picard changed the bytes"
            assert (ours.returncode, verb.returncode) == (0, 0), f"{label}: {ours.stderr!r} {verb.stderr!r}"
            assert ours.stdout == picard.stdout, f"{label}: the Picard face differs"
            assert verb.stdout == picard.stdout, f"{label}: the verb differs"
            # No log line in either: DEBUG_FREQUENCY is 0 or less, or VERBOSITY drops them.
            assert log_lines("picard", picard.stderr) == [], label
            assert log_lines("ours", ours.stderr) == [], label
            assert log_lines("verb", verb.stderr) == [], label
            assert CHAIN_NOTE in verb.stderr.decode(), label
            runs += 1
    total = sum(len(data) for data in shapes.values())
    RESULTS.append((f"{TOOL} bytes", runs, runs, f"6 shapes, {total} bytes each pass, stdout identical"))


class Running:
    """One tool running on pipes: its log lines read one at a time, bytes fed to its
    standard input on a thread, standard output read only by `finish`."""

    def __init__(self, side, args):
        argv, env = command(side, args)
        self.label = f"{side} {args}"
        self.pattern = {"picard": LOG, "ours": OURS_LOG, "verb": VERB_LOG}[side]
        self.proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, env=env)
        self.feeders = []

    def feed(self, data, close=False):
        """Write `data` to standard input on a thread, after every earlier `feed`."""
        earlier = list(self.feeders)

        def write():
            for thread in earlier:
                thread.join()
            try:
                self.proc.stdin.write(data)
                if close:
                    self.proc.stdin.close()
                else:
                    self.proc.stdin.flush()
            except BrokenPipeError:
                pass

        thread = threading.Thread(target=write)
        thread.start()
        self.feeders.append(thread)

    def next_line(self):
        """The message of the next log line on standard error."""
        while True:
            line = self.proc.stderr.readline()
            assert line, f"{self.label}: standard error ended"
            match = self.pattern.match(line.decode(errors="replace").rstrip("\n"))
            if match:
                return match.groups()[-1]

    def finish(self):
        """Close standard input after what was fed, read standard output to its end
        (standard error beside it, so neither pipe fills) and return it with the exit
        code."""
        self.feed(b"", close=True)
        rest = []
        reader = threading.Thread(target=lambda: rest.append(self.proc.stderr.read()))
        reader.start()
        out = self.proc.stdout.read()
        code = self.proc.wait(timeout=300)
        for thread in self.feeders:
            thread.join()
        reader.join()
        return out, code


def staged(side, args, data, lines, drain_first):
    """Start a tool and wait for `lines` log lines on standard error. With `drain_first`
    the input is fed from the start and standard output is read only after those lines
    (the buffer fills); without it the input is fed after them (the buffer is idle).
    Returns the log lines, standard output and the exit code."""
    tool = Running(side, args)
    if drain_first:
        tool.feed(data)
    seen = [tool.next_line() for _ in range(lines)]
    if not drain_first:
        tool.feed(data)
    out, code = tool.finish()
    return seen, out, code


def test_log_lines_of_an_idle_buffer_equal():
    """The input arrives after two lines: both read `used 0` in both tools."""
    cases = [
        (["DEBUG_FREQUENCY=1"], "Fifo buffer used 0 / 536,870,912 (0%)."),
        (["DEBUG_FREQUENCY=1", "BUFFER_SIZE=1000", "NAME=lane1"], "Fifo buffer lane1 used 0 / 1,000 (0%)."),
        (["DEBUG_FREQUENCY=1", "BUFFER_SIZE=1234567", "NAME=null", "VERBOSITY=DEBUG"],
         "Fifo buffer used 0 / 1,234,567 (0%)."),
        (["--DEBUG_FREQUENCY", "1", "--BUFFER_SIZE", "1", "--NAME", "a b"], "Fifo buffer a b used 0 / 1 (0%)."),
    ]
    for args, expected in cases:
        results = {side: staged(side, args, b"late bytes\n", 2, drain_first=False)
                   for side in ("picard", "ours")}
        for side, (seen, out, code) in results.items():
            assert (code, out) == (0, b"late bytes\n"), f"{side} {args}"
            assert seen == [expected, expected], f"{side} {args}"
    verb = staged("verb", ["--debug-frequency", "1", "--buffer-size", "1000", "--name", "lane1", "-q"],
                  b"late bytes\n", 2, drain_first=False)
    assert verb == ([cases[1][1]] * 2, b"late bytes\n", 0)
    RESULTS.append((f"{TOOL} idle log lines", 2 * len(cases), 2 * len(cases), "message text identical"))


def test_log_lines_of_a_full_buffer_equal(shapes):
    """Nobody reads standard output until the third line: the buffer is full in the
    second and third lines of both tools, and every byte still arrives."""
    data = shapes["random_3mib"]
    cases = [
        (["DEBUG_FREQUENCY=1", "BUFFER_SIZE=1000", "IO_SIZE=100"],
         "Fifo buffer used 0 / 1,000 (0%).", "Fifo buffer used 1,000 / 1,000 (100%)."),
        (["DEBUG_FREQUENCY=1", "BUFFER_SIZE=131072", "IO_SIZE=4096", "NAME=full"],
         "Fifo buffer full used 0 / 131,072 (0%).", "Fifo buffer full used 131,072 / 131,072 (100%)."),
    ]
    for args, first, full in cases:
        for side in ("picard", "ours"):
            seen, out, code = staged(side, args, data, 3, drain_first=True)
            assert (code, out == data) == (0, True), f"{side} {args}"
            assert seen[1:] == [full, full], f"{side} {args}: {seen}"
            # Picard's first line holds what its reader thread put before the log thread ran.
            if side == "ours":
                assert seen[0] == first
    RESULTS.append((f"{TOOL} full log lines", 2 * len(cases), 2 * len(cases), "message text identical"))


USED = re.compile(r"^Fifo buffer used ([\d,]+) / ([\d,]+) \((\d+)%\)\.$")


def counts(message):
    """`used`, `capacity` and the printed percent of one log message."""
    used, capacity, percent = USED.match(message).groups()
    return int(used.replace(",", "")), int(capacity.replace(",", "")), int(percent)


def passed_on(side):
    """The bytes a tool passes on before its buffer starts to fill when nobody reads
    standard output and IO_SIZE is 1: what the pipe takes, plus in Picard the 8192 bytes
    of `System.out`. Found as the bytes fed less the bytes a large buffer holds."""
    fed = 200_000
    tool = Running(side, ["DEBUG_FREQUENCY=1", "BUFFER_SIZE=400000", "IO_SIZE=1"])
    tool.feed(bytes(fed))
    tool.next_line()
    held = [counts(tool.next_line())[0]]
    while len(held) < 2 or held[-1] != held[-2]:
        held.append(counts(tool.next_line())[0])
    out, code = tool.finish()
    assert (code, len(out)) == (0, fed), tool.label
    return fed - held[-1]


def ramp(side, steps, step):
    """Log messages of a 200-byte buffer filled `step` bytes a second, by bytes held.
    Standard output is not read meanwhile, so each byte past `passed_on` stays held."""
    before = passed_on(side)
    tool = Running(side, ["DEBUG_FREQUENCY=1", "BUFFER_SIZE=200", "IO_SIZE=1"])
    tool.next_line()
    tool.feed(bytes(before + 3))
    seen = {}
    for _ in range(steps):
        message = tool.next_line()
        seen[counts(message)[0]] = message
        tool.feed(bytes(step))
    out, code = tool.finish()
    assert (code, len(out)) == (0, before + 3 + steps * step), tool.label
    return seen


def test_log_lines_of_a_part_full_buffer_round_as_the_jdk_does():
    """A 200-byte buffer holds 3, 40, 77, 114, 151 and 188 bytes in turn: 1.5%, 38.5%
    and 75.5% are ties, printed as 2%, 38% and 76% (to the even neighbour) by both
    tools. Every message either tool prints must follow that rule, and the messages for
    the same count must be equal. The pipe's capacity is measured per run; a run where
    it changed between the two measurements is repeated."""
    expected = {3: 2, 40: 20, 77: 38, 114: 57, 151: 76, 188: 94}
    for _attempt in range(3):
        picard = ramp("picard", 6, 37)
        ours = ramp("ours", 6, 37)
        if all(set(expected) <= set(seen) for seen in (picard, ours)):
            break
    for side, seen in (("picard", picard), ("ours", ours)):
        for message in seen.values():
            used, capacity, percent = counts(message)
            assert capacity == 200, f"{side}: {message}"
            assert f"{percent}" == f"{used / capacity * 100:.0f}", f"{side}: {message}"
    for used, percent in expected.items():
        message = f"Fifo buffer used {used} / 200 ({percent}%)."
        assert picard.get(used) == message, f"Picard: {sorted(picard.items())}"
        assert ours.get(used) == message, f"bammetrics: {sorted(ours.items())}"
    RESULTS.append((f"{TOOL} part-full log lines", len(expected), len(expected),
                    "message text identical, 3 ties to even"))


def test_bad_options_fail_in_both_with_nothing_copied():
    for args, message in BAD:
        picard = run("picard", args, b"never copied")
        ours = run("ours", args, b"never copied")
        assert picard.returncode == 1, f"{args}: Picard exit {picard.returncode}"
        assert ours.returncode == 1, f"{args}: {ours.stderr!r}"
        assert (picard.stdout, ours.stdout) == (b"", b""), args
        if message is not None:
            assert message in picard.stderr.decode(errors="replace"), args
            assert message in ours.stderr.decode(), f"{args}: {ours.stderr!r}"
        # The verb exits 2 for a usage error. It keeps the last value of an option
        # given twice, as every verb does.
        verb = run("verb", args, b"never copied")
        twice = len(args) == 2
        assert (verb.returncode, verb.stdout) == ((0, b"never copied") if twice else (2, b"")), \
            f"{args}: {verb.stderr!r}"
    RESULTS.append((f"{TOOL} bad options", len(BAD), len(BAD), "both exit 1, standard output empty"))


def still_running(args, data, seconds):
    """True when Picard, once its JVM has printed the banner line, has not ended
    `seconds` after its input was written and closed. The process is then killed."""
    argv, env = command("picard", args)
    proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, env=env)
    try:
        while True:
            line = proc.stderr.readline()
            assert line, f"{args}: Picard ended before its banner"
            if line.startswith(b"INFO\t") and b"\tFifoBuffer\t" in line:
                break
        proc.stdin.write(data)
        proc.stdin.close()
        try:
            proc.wait(timeout=seconds)
        except subprocess.TimeoutExpired:
            return True
        return False
    finally:
        proc.kill()
        proc.wait()


def test_a_size_of_zero_never_ends_in_picard_and_is_refused_here():
    for args in (["BUFFER_SIZE=0"], ["IO_SIZE=0"]):
        assert still_running(args, b"waits forever", 6), args
        ours = run("ours", args, b"waits forever")
        assert (ours.returncode, ours.stdout) == (1, b"")
        assert b"a size from 1 to 2147483645 bytes is needed" in ours.stderr
    RESULTS.append((f"{TOOL} size 0", 2, 2, "Picard never ends; bammetrics exits 1"))


def test_a_closed_standard_output_ends_both_with_exit_0(shapes):
    """The reader takes 10 bytes and closes. Picard's PrintStream hides the error and
    the tool reads the rest of its input; bammetrics stops at the failed write."""
    data = shapes["random_3mib"] * 8
    fed = {}
    # A buffer smaller than the input: the reader thread cannot hold all of it.
    for side in ("picard", "ours", "verb"):
        argv, env = command(side, ["BUFFER_SIZE=65536"] + (["-q"] if side == "verb" else []))
        proc = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, env=env)
        errors = []
        reader = threading.Thread(target=lambda: errors.append(proc.stderr.read()))
        reader.start()

        def feed():
            try:
                proc.stdin.write(data)
                proc.stdin.close()
                fed[side] = True
            except BrokenPipeError:
                fed[side] = False

        feeder = threading.Thread(target=feed)
        feeder.start()
        first = proc.stdout.read(10)
        proc.stdout.close()
        code = proc.wait(timeout=300)
        feeder.join()
        reader.join()
        assert first == data[:10], side
        assert code == 0, f"{side}: {errors[0][-300:]!r}"
        assert log_lines(side, errors[0]) == [], side
        if side != "picard":
            assert errors[0] == b"", f"{side}: {errors[0]!r}"
    assert fed == {"picard": True, "ours": False, "verb": False}
    RESULTS.append((f"{TOOL} closed output", 3, 3, "exit 0 in both; bammetrics stops reading"))

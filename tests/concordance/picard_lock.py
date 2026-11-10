"""One Picard JVM at a time for every harness of one user on this machine.

    python3 tests/concordance/picard_lock.py CollectInsertSizeMetrics I=in.bam O=out.txt
    python3 tests/concordance/picard_lock.py --status

The command form takes an exclusive `flock` on `~/.cache/bammetrics/picard.lock`
(`PICARD_LOCK_FILE` names another file), waits until `/bin/ps` lists no JVM
running `picard.jar`, then replaces itself with the real `picard`. The JVM
inherits the lock's descriptor, so the lock lasts as long as the JVM and is
released by the kernel when it exits, however it exits. Standard input,
output, error, the exit code and the environment are the real command's:
`JAVA_OPTS` reaches the Homebrew `picard` wrapper unchanged.

No test calls this file by name. `run_picard.picard_env()` puts `shim()`, a
two-line `picard` script in `target/picard-queue`, first on `PATH`, so every
`subprocess.run(["picard", ...], env=picard_env())` of the harness queues here.
A cached Picard result (`target/fuzz/cache`, the `run_picard.py` manifest) starts
no process and takes no lock. `held()` is the same lock for a caller that times
the JVM itself (`benches/bench.py`).

Environment:

  PICARD_LOCK_FILE          the lock file (default ~/.cache/bammetrics/picard.lock)
  PICARD_LOCK_TIMEOUT       seconds to wait for the lock and for other JVMs; unset
                            waits without limit; exit 75 when it runs out
  PICARD_LOCK_WAIT_FOREIGN  0 skips the wait for JVMs started outside the lock

The wait for other JVMs covers a harness of a branch older than this file and a
Picard run started by hand. It is one check before the start: a JVM started
outside the lock during a queued run is not seen.
"""
import contextlib
import fcntl
import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_LOCK = Path.home() / ".cache" / "bammetrics" / "picard.lock"
# Exit code when PICARD_LOCK_TIMEOUT runs out: EX_TEMPFAIL of sysexits.h.
TIMED_OUT = 75
# Set for the JVM and its children: a nested call must not wait for its own lock.
HELD = "PICARD_LOCK_HELD"
LOCK_POLL = 0.05
FOREIGN_POLL = 1.0


class Timeout(Exception):
    """The wait passed PICARD_LOCK_TIMEOUT."""


def lock_path(env=os.environ):
    return Path(env.get("PICARD_LOCK_FILE") or DEFAULT_LOCK)


def timeout(env=os.environ):
    """PICARD_LOCK_TIMEOUT as seconds, None when unset or empty."""
    text = env.get("PICARD_LOCK_TIMEOUT", "")
    return float(text) if text else None


def open_lock(path):
    """The lock file's descriptor: created 0600 in a 0700 directory, never
    through a symbolic link."""
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    return os.open(path, os.O_RDWR | os.O_CREAT | os.O_NOFOLLOW, 0o600)


def acquire(path, deadline=None, sleep=time.sleep):
    """Take the exclusive lock, polling until `deadline` (a time.monotonic
    value, None for no limit). Returns the descriptor that holds it."""
    fd = open_lock(path)
    while True:
        try:
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return fd
        except BlockingIOError:
            if deadline is not None and time.monotonic() >= deadline:
                os.close(fd)
                raise Timeout(f"{path} is held: {holder(path)}") from None
            sleep(LOCK_POLL)


def note_holder(fd, argv):
    """Write who holds the lock into the file, for `--status` and `cat`."""
    line = f"pid {os.getpid()}\tsince {time.strftime('%Y-%m-%d %H:%M:%S')}\t{shlex.join(argv)}\n"
    os.ftruncate(fd, 0)
    os.pwrite(fd, line.encode(), 0)


def holder(path):
    """The line the last holder wrote."""
    try:
        return path.read_text(errors="replace").strip() or "no holder line"
    except OSError as error:
        return f"unreadable: {error}"


def process_table():
    """`pid args` lines of every process, from /bin/ps."""
    out = subprocess.run(["/bin/ps", "-eo", "pid=,args="], check=True, capture_output=True)
    return out.stdout.decode(errors="replace").splitlines()


def picard_jvms(lines):
    """Process ids of the JVMs running picard.jar among `pid args` lines."""
    pids = []
    for line in lines:
        fields = line.split()
        if len(fields) < 2 or not fields[0].isdigit():
            continue
        program = fields[1]
        if os.path.basename(program) == "java" and "picard.jar" in line:
            pids.append(int(fields[0]))
    return pids


def wait_for_foreign(deadline=None, table=process_table, sleep=time.sleep):
    """Return when no JVM runs picard.jar. This process holds the lock, so
    any such JVM was started outside it."""
    while True:
        pids = picard_jvms(table())
        if not pids:
            return
        if deadline is not None and time.monotonic() >= deadline:
            raise Timeout(f"another Picard JVM is running outside the lock: pid {pids[0]}")
        sleep(FOREIGN_POLL)


def waits_for_foreign(env=os.environ):
    return env.get("PICARD_LOCK_WAIT_FOREIGN", "1") != "0"


def deadline_from(env=os.environ):
    seconds = timeout(env)
    return None if seconds is None else time.monotonic() + seconds


@contextlib.contextmanager
def held(argv=("held",), env=os.environ):
    """The lock for the length of a `with` block, for a caller that starts
    the real `picard` itself. Nested use in one process tree takes it once."""
    if env.get(HELD):
        yield
        return
    deadline = deadline_from(env)
    fd = acquire(lock_path(env), deadline)
    try:
        note_holder(fd, list(argv))
        if waits_for_foreign(env):
            wait_for_foreign(deadline)
        os.environ[HELD] = str(os.getpid())
        try:
            yield
        finally:
            del os.environ[HELD]
    finally:
        os.close(fd)


def shim(directory, python=sys.executable):
    """Write `directory/picard`, a script that runs this file with its
    arguments, and return the directory. The script is replaced by rename, and
    only when its text changes, so a process starting it never reads half."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    return write_script(directory / "picard",
                        f'#!/bin/sh\nexec {shlex.quote(str(python))} '
                        f'{shlex.quote(str(Path(__file__).resolve()))} "$@"\n').parent


def write_script(path, text):
    """An executable script at `path` holding `text`, replaced by rename when
    the text differs."""
    if not (path.exists() and path.read_text() == text and os.access(path, os.X_OK)):
        fresh = path.with_name(f".{path.name}.{os.getpid()}")
        fresh.write_text(text)
        fresh.chmod(0o755)
        os.replace(fresh, path)
    return path


def real_picard(env=os.environ):
    """The first `picard` on PATH that is not a shim of this file."""
    for directory in env.get("PATH", "").split(os.pathsep):
        found = shutil.which("picard", path=directory) if directory else None
        if found and not is_shim(found):
            return found
    raise FileNotFoundError("no picard on PATH beside the queue shim")


def is_shim(path):
    try:
        with open(path, "rb") as handle:
            return b"picard_lock.py" in handle.read(4096)
    except OSError:
        return False


def status(env=os.environ):
    """One line: whether the lock is free, and the last holder's line."""
    path = lock_path(env)
    if not path.exists():
        return f"free: {path} does not exist yet"
    fd = open_lock(path)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        return f"held: {holder(path)}"
    finally:
        os.close(fd)
    return f"free: last holder {holder(path)}"


def main(argv):
    if argv == ["--status"]:
        print(status())
        return 0
    picard = real_picard()
    if not os.environ.get(HELD):
        deadline = deadline_from()
        try:
            fd = acquire(lock_path(), deadline)
            note_holder(fd, ["picard", *argv])
            if waits_for_foreign():
                wait_for_foreign(deadline)
        except Timeout as error:
            print(f"picard_lock: {error}", file=sys.stderr)
            return TIMED_OUT
        os.set_inheritable(fd, True)
        os.environ[HELD] = str(os.getpid())
    os.execv(picard, ["picard", *argv])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

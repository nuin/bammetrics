//! `FifoBuffer`, mirroring Picard 3.5.0 `picard.util.FifoBuffer` and its
//! `CircularByteBuffer`.
//!
//! Picard copies standard input to standard output through a ring of
//! `BUFFER_SIZE` bytes (default 512 MiB). One thread reads up to `IO_SIZE`
//! bytes (default 64 KiB) at a time and puts them into the ring, waiting
//! while it is full. A second thread takes up to `IO_SIZE` bytes and
//! writes them. A producer that writes in bursts therefore never waits
//! for a consumer that reads in bursts, up to `BUFFER_SIZE` bytes.
//!
//! With `DEBUG_FREQUENCY=N` above 0 Picard logs one `INFO` line when it
//! starts and one every N seconds:
//! `Fifo buffer [NAME ]used 1,024 / 536,870,912 (0%).` The counts are
//! `DecimalFormat("#,##0")`, the percent is `NumberFormat.getPercentInstance`
//! (no fraction digit, `HALF_EVEN`). `VERBOSITY=WARNING` or `ERROR` drops
//! the lines. The tool takes no `INPUT`, no `OUTPUT`, no `ASSUME_SORTED`
//! and no `GA4GH_CLIENT_SECRETS` (`Unrecognized option`).
//!
//! What bammetrics does differently (KNOWN_DIFFERENCES.md):
//!
//! - The first log line is printed before the first read, so it always
//!   reads `used 0`. Picard starts its log thread beside the copy, so its
//!   first line holds whatever arrived by then.
//! - `BUFFER_SIZE=0` and `IO_SIZE=0` are usage errors. Picard waits
//!   forever on the first byte.
//! - When the reader of standard output closes it, the copy stops and
//!   exits 0. Picard reads the rest of its input and drops it.
//! - A failed write of another kind (a full disk) exits 1. Picard's
//!   `PrintStream` hides it and exits 0.
//! - The ring starts again at its first byte whenever it is empty, so the
//!   memory touched follows the most bytes ever held.
//!
//! The registry row gives the Picard face, the `fifo-buffer` verb and the
//! Python keywords. [`TOOL`]'s `run` refuses: a call that returns rows has
//! no pipe to buffer. Both command faces call [`copy_standard_streams`].

use std::fs::File;
use std::io::{self, Read, Write};
use std::os::fd::AsFd;
use std::path::PathBuf;
use std::sync::{Arc, Condvar, Mutex, MutexGuard, PoisonError};
use std::thread;
use std::time::Duration;

use crate::error::{Error, Result};
use crate::java;
use crate::metric::ToolOutput;
use crate::registry::{Invocation, Opt, REFERENCE_NOT_DECODED, Tool};

/// Picard's default `BUFFER_SIZE`: 512 MiB.
pub const DEFAULT_BUFFER_SIZE: usize = 512 * 1024 * 1024;
/// Picard's default `IO_SIZE`: 64 KiB, the usual capacity of a pipe.
pub const DEFAULT_IO_SIZE: usize = 64 * 1024;
/// The largest array the JVM allocates. Picard fails above it with
/// `Requested array size exceeds VM limit`.
pub const MAX_SIZE: i64 = 2_147_483_645;

/// What the Python function and every other caller of [`TOOL`]'s `run` get.
pub const NO_ROWS: &str = "FifoBuffer copies the standard input of the command to its standard \
                           output and returns no rows; run `bammetrics fifo-buffer` (or \
                           `bammetrics FifoBuffer`) between two commands of a pipe";

pub const TOOL: Tool = Tool {
    name: "FifoBuffer",
    summary: "standard input copied to standard output through a memory buffer of BUFFER_SIZE \
              bytes, so a producer and a consumer that work in bursts do not wait for each \
              other. A command only: it returns no rows, and the Python function raises \
              UsageError, because a Python call has no pipe to buffer.",
    options: &[
        Opt::int(
            &["BUFFER_SIZE"],
            "The size of the memory buffer in bytes, 1 to 2147483645 (default 536870912).",
        ),
        Opt::int(
            &["IO_SIZE"],
            "The most bytes one read or one write moves, 1 to 2147483645 (default 65536).",
        ),
        Opt::int(
            &["DEBUG_FREQUENCY"],
            "Seconds between two log lines of the bytes held; 0 or less logs none (default 0).",
        ),
        Opt::text(&["NAME"], "A name for the buffer in the log lines."),
        Opt::text(
            &["VERBOSITY"],
            "ERROR or WARNING drops the log lines of DEBUG_FREQUENCY (default INFO).",
        ),
        Opt::rejected(&["INPUT", "I"]),
        Opt::rejected(&["OUTPUT", "O"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
        REFERENCE_NOT_DECODED,
    ],
    run: run_tool,
    collector: None,
};

/// The options of one copy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    /// The ring's capacity in bytes, at least 1.
    pub buffer_size: usize,
    /// The most bytes of one read or one write, at least 1.
    pub io_size: usize,
    /// Seconds between two log lines, `None` for no line.
    pub debug_frequency: Option<u64>,
    /// The buffer's name in the log lines.
    pub name: Option<String>,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            buffer_size: DEFAULT_BUFFER_SIZE,
            io_size: DEFAULT_IO_SIZE,
            debug_frequency: None,
            name: None,
        }
    }
}

/// A Picard `int` option: `Integer.parseInt` range, `null` refused as
/// Picard's reflection refuses it.
fn int_option(inv: &Invocation, name: &str) -> Result<Option<i64>> {
    let Some(text) = inv.get(name) else {
        return Ok(None);
    };
    if text == "null" {
        return Err(Error::Usage(format!(
            "{name}=null: an int option cannot be null"
        )));
    }
    java::parse_integer(text, 32).map(Some).map_err(|_| {
        Error::Usage(format!(
            "{name}={text}: Problem constructing Integer from the string '{text}'."
        ))
    })
}

/// A size option: 1 to [`MAX_SIZE`]. Picard throws
/// `NegativeArraySizeException` below 0, never ends at 0 and fails with
/// `Requested array size exceeds VM limit` above [`MAX_SIZE`].
fn size_option(inv: &Invocation, name: &str, default: usize) -> Result<usize> {
    let Some(value) = int_option(inv, name)? else {
        return Ok(default);
    };
    if !(1..=MAX_SIZE).contains(&value) {
        return Err(Error::Usage(format!(
            "{name}={value}: a size from 1 to {MAX_SIZE} bytes is needed"
        )));
    }
    usize::try_from(value)
        .map_err(|_| Error::Usage(format!("{name}={value}: too large for this machine")))
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let shows_info = match inv.get("VERBOSITY") {
            None | Some("INFO" | "DEBUG") => true,
            Some("WARNING" | "ERROR") => false,
            Some(other) => {
                return Err(Error::Usage(format!(
                    "VERBOSITY={other}: one of ERROR, WARNING, INFO, DEBUG"
                )));
            }
        };
        let frequency = int_option(inv, "DEBUG_FREQUENCY")?.unwrap_or(0);
        Ok(Self {
            buffer_size: size_option(inv, "BUFFER_SIZE", DEFAULT_BUFFER_SIZE)?,
            io_size: size_option(inv, "IO_SIZE", DEFAULT_IO_SIZE)?,
            debug_frequency: u64::try_from(frequency)
                .ok()
                .filter(|seconds| *seconds > 0 && shows_info),
            name: inv
                .get("NAME")
                .filter(|name| *name != "null")
                .map(str::to_owned),
        })
    }
}

/// `n` as `DecimalFormat("#,##0")` prints it: `536,870,912`.
fn grouped(n: usize) -> String {
    let digits = n.to_string();
    let mut text = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, digit) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            text.push(',');
        }
        text.push(digit);
    }
    text
}

/// Picard's log message for `used` bytes held of `capacity`. The percent
/// is the double `used / capacity * 100` with no fraction digit, a tie to
/// the even neighbour, as `NumberFormat.getPercentInstance` prints it.
pub fn debug_line(name: Option<&str>, used: usize, capacity: usize) -> String {
    let name = name.map_or_else(String::new, |name| format!("{name} "));
    let percent = used as f64 / capacity as f64 * 100.0;
    format!(
        "Fifo buffer {name}used {} / {} ({percent:.0}%).",
        grouped(used),
        grouped(capacity)
    )
}

/// The bytes held between the reader and the writer.
struct Ring {
    bytes: Vec<u8>,
    /// Index of the oldest byte.
    start: usize,
    /// Bytes held.
    len: usize,
    /// The input ended or failed: nothing more is put.
    input_done: bool,
    /// The output ended or failed: nothing more is taken.
    output_done: bool,
    input_error: Option<io::Error>,
    /// Bytes put so far.
    read: u64,
}

impl Ring {
    /// Copy as much of `data` as fits behind the newest byte. Returns the
    /// count.
    fn put(&mut self, data: &[u8]) -> usize {
        let capacity = self.bytes.len();
        let count = data.len().min(capacity - self.len);
        let end = (self.start + self.len) % capacity;
        let first = count.min(capacity - end);
        self.bytes[end..end + first].copy_from_slice(&data[..first]);
        self.bytes[..count - first].copy_from_slice(&data[first..count]);
        self.len += count;
        self.read += count as u64;
        count
    }

    /// Move the oldest bytes into `out`. Returns the count. An empty ring
    /// starts again at its first byte.
    fn take(&mut self, out: &mut [u8]) -> usize {
        let capacity = self.bytes.len();
        let count = out.len().min(self.len);
        let first = count.min(capacity - self.start);
        out[..first].copy_from_slice(&self.bytes[self.start..self.start + first]);
        out[first..count].copy_from_slice(&self.bytes[..count - first]);
        self.len -= count;
        self.start = if self.len == 0 {
            0
        } else {
            (self.start + count) % capacity
        };
        count
    }
}

/// The ring with its two wake-ups: `filled` after a put, `drained` after a
/// take.
struct Fifo {
    ring: Mutex<Ring>,
    filled: Condvar,
    drained: Condvar,
}

impl Fifo {
    fn new(capacity: usize) -> Self {
        Self {
            ring: Mutex::new(Ring {
                bytes: vec![0; capacity],
                start: 0,
                len: 0,
                input_done: false,
                output_done: false,
                input_error: None,
                read: 0,
            }),
            filled: Condvar::new(),
            drained: Condvar::new(),
        }
    }

    /// The ring. A thread that panicked while holding it left counters
    /// that still describe the bytes held, so the lock is taken anyway.
    fn lock(&self) -> MutexGuard<'_, Ring> {
        self.ring.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// Put all of `data`, waiting while the ring is full. False when the
    /// output ended first.
    fn put(&self, mut data: &[u8]) -> bool {
        let mut ring = self.lock();
        while !data.is_empty() {
            while ring.len == ring.bytes.len() && !ring.output_done {
                ring = self
                    .drained
                    .wait(ring)
                    .unwrap_or_else(PoisonError::into_inner);
            }
            if ring.output_done {
                return false;
            }
            let count = ring.put(data);
            data = &data[count..];
            self.filled.notify_one();
        }
        true
    }

    /// Take up to `out.len()` bytes, waiting while the ring is empty and
    /// the input goes on. 0 when the input ended and nothing is held.
    fn take(&self, out: &mut [u8]) -> usize {
        let mut ring = self.lock();
        while ring.len == 0 && !ring.input_done {
            ring = self
                .filled
                .wait(ring)
                .unwrap_or_else(PoisonError::into_inner);
        }
        let count = ring.take(out);
        if count > 0 {
            self.drained.notify_one();
        }
        count
    }

    fn end_input(&self, error: Option<io::Error>) {
        let mut ring = self.lock();
        ring.input_done = true;
        ring.input_error = error;
        self.filled.notify_one();
    }

    fn end_output(&self) {
        self.lock().output_done = true;
        self.drained.notify_one();
    }
}

/// The reader thread: `input` into the ring, `io_size` bytes at a time,
/// until it ends, fails or the output ends.
fn fill(fifo: &Fifo, mut input: impl Read, io_size: usize) {
    let mut chunk = vec![0; io_size];
    let error = loop {
        match input.read(&mut chunk) {
            Ok(0) => break None,
            Ok(count) => {
                if !fifo.put(&chunk[..count]) {
                    break None;
                }
            }
            Err(error) if error.kind() == io::ErrorKind::Interrupted => {}
            Err(error) => break Some(error),
        }
    };
    fifo.end_input(error);
}

/// Tells the log thread that the copy ended.
struct Stop {
    stopped: Mutex<bool>,
    changed: Condvar,
}

/// The log thread: one line every `period` until `stop`.
fn report(fifo: &Fifo, stop: &Stop, period: Duration, name: Option<&str>, log: impl Fn(&str)) {
    let mut stopped = stop.stopped.lock().unwrap_or_else(PoisonError::into_inner);
    loop {
        let (guard, wait) = stop
            .changed
            .wait_timeout(stopped, period)
            .unwrap_or_else(PoisonError::into_inner);
        stopped = guard;
        if *stopped {
            return;
        }
        if wait.timed_out() {
            let ring = fifo.lock();
            let (used, capacity) = (ring.len, ring.bytes.len());
            drop(ring);
            log(&debug_line(name, used, capacity));
        }
    }
}

/// What one [`copy`] moved.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Copied {
    /// Bytes read from the input.
    pub read: u64,
    /// Bytes written to the output.
    pub written: u64,
    /// The reader of the output closed it before the input ended.
    pub output_closed: bool,
}

/// The names the two streams have in an error.
const STANDARD_INPUT: &str = "standard input";
const STANDARD_OUTPUT: &str = "standard output";

/// Copy `input` to `output` through a ring of `options.buffer_size` bytes.
/// `input` is read on a thread of its own, `output` is written by the
/// caller's thread. `log` gets each message of `DEBUG_FREQUENCY`: one
/// before the first read, then one per period on a third thread.
///
/// An output closed by its reader (`BrokenPipe`) ends the copy without an
/// error; the reader thread is then left to end at its next read. Any
/// other write error, and a read error once the bytes read before it are
/// written, is the error returned.
pub fn copy<R, W, L>(input: R, mut output: W, options: &Options, log: L) -> Result<Copied>
where
    R: Read + Send + 'static,
    W: Write,
    L: Fn(&str) + Send + 'static,
{
    let capacity = options.buffer_size;
    if capacity == 0 || options.io_size == 0 {
        return Err(Error::Usage(
            "BUFFER_SIZE and IO_SIZE must be 1 or more".to_owned(),
        ));
    }
    // A size the allocator refuses is an error here, not an abort below.
    Vec::<u8>::new().try_reserve_exact(capacity).map_err(|_| {
        Error::Usage(format!(
            "BUFFER_SIZE={capacity}: the memory cannot be allocated"
        ))
    })?;
    let fifo = Arc::new(Fifo::new(capacity));
    let stop = Arc::new(Stop {
        stopped: Mutex::new(false),
        changed: Condvar::new(),
    });
    let reporter = options.debug_frequency.map(|seconds| {
        let name = options.name.clone();
        log(&debug_line(name.as_deref(), 0, capacity));
        let (fifo, stop) = (Arc::clone(&fifo), Arc::clone(&stop));
        thread::spawn(move || {
            report(
                &fifo,
                &stop,
                Duration::from_secs(seconds),
                name.as_deref(),
                log,
            );
        })
    });
    let io_size = options.io_size;
    let reader = {
        let fifo = Arc::clone(&fifo);
        thread::spawn(move || fill(&fifo, input, io_size))
    };

    let mut chunk = vec![0; io_size];
    let mut written = 0u64;
    let outcome = loop {
        let count = fifo.take(&mut chunk);
        if count == 0 {
            break output.flush();
        }
        match output.write_all(&chunk[..count]) {
            Ok(()) => written += count as u64,
            Err(error) => break Err(error),
        }
    };
    fifo.end_output();
    *stop.stopped.lock().unwrap_or_else(PoisonError::into_inner) = true;
    stop.changed.notify_one();
    if let Some(reporter) = reporter {
        // The log thread waits on `stop` alone, so it ends at once.
        let _ = reporter.join();
    }
    if outcome.is_ok() {
        // The ring was empty with the input ended: the reader has returned.
        let _ = reader.join();
    }
    let mut ring = fifo.lock();
    let copied = Copied {
        read: ring.read,
        written,
        output_closed: false,
    };
    match outcome {
        Err(error) if error.kind() == io::ErrorKind::BrokenPipe => Ok(Copied {
            output_closed: true,
            ..copied
        }),
        Err(source) => Err(Error::Output {
            path: PathBuf::from(STANDARD_OUTPUT),
            source,
        }),
        Ok(()) => match ring.input_error.take() {
            Some(source) => Err(Error::Input {
                path: PathBuf::from(STANDARD_INPUT),
                source,
            }),
            None => Ok(copied),
        },
    }
}

/// [`copy`] from the process's standard input to its standard output, each
/// through a duplicate of its descriptor, so no buffer of the standard
/// library sits in between.
pub fn copy_standard_streams<L>(options: &Options, log: L) -> Result<Copied>
where
    L: Fn(&str) + Send + 'static,
{
    let input = io::stdin()
        .as_fd()
        .try_clone_to_owned()
        .map_err(|source| Error::Input {
            path: PathBuf::from(STANDARD_INPUT),
            source,
        })?;
    let output = io::stdout()
        .as_fd()
        .try_clone_to_owned()
        .map_err(|source| Error::Output {
            path: PathBuf::from(STANDARD_OUTPUT),
            source,
        })?;
    copy(File::from(input), File::from(output), options, log)
}

/// The registry's `run`: the options are checked, then the call is
/// refused, because it would return rows and the tool has none.
fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    Options::from_invocation(inv)?;
    Err(Error::Usage(NO_ROWS.to_owned()))
}

/// The Picard face: parse, copy standard input to standard output, print
/// each log line on standard error as Picard does without its timestamp
/// (`INFO\tFifoBuffer\t<message>`). The status is Picard's exit code, 0.
pub fn run_cli(args: &[String]) -> Result<u8> {
    let inv = TOOL.parse(args)?;
    let options = Options::from_invocation(&inv)?;
    for note in &inv.notes {
        eprintln!("bammetrics: {note}");
    }
    copy_standard_streams(&options, |line| eprintln!("INFO\tFifoBuffer\t{line}"))?;
    Ok(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::mpsc;

    fn options(buffer_size: usize, io_size: usize) -> Options {
        Options {
            buffer_size,
            io_size,
            ..Options::default()
        }
    }

    /// Bytes that differ at every position modulo 251, zero and 0xFF
    /// included.
    fn bytes(count: usize) -> Vec<u8> {
        (0..count).map(|i| (i * 7 % 251) as u8).collect()
    }

    fn copied(data: &[u8], options: &Options) -> (Vec<u8>, Copied) {
        let mut out = Vec::new();
        let result = copy(io::Cursor::new(data.to_vec()), &mut out, options, |_| {}).expect("copy");
        (out, result)
    }

    #[test]
    fn output_equals_input_for_every_size_of_buffer_and_io() {
        for (len, buffer_size, io_size) in [
            (0, 16, 4),
            (1, 1, 1),
            (10, 16, 4),
            (1000, 16, 4),
            (1000, 1, 1),
            (1000, 10, 100),
            (1000, 7, 3),
            (100_000, 4096, 65_536),
            (300_000, DEFAULT_BUFFER_SIZE, DEFAULT_IO_SIZE),
        ] {
            let data = bytes(len);
            let (out, result) = copied(&data, &options(buffer_size, io_size));
            assert_eq!(out, data, "{len} bytes, buffer {buffer_size}, io {io_size}");
            assert_eq!(
                result,
                Copied {
                    read: len as u64,
                    written: len as u64,
                    output_closed: false
                }
            );
        }
    }

    #[test]
    fn the_ring_wraps_and_starts_again_when_empty() {
        let mut ring = Fifo::new(8).ring.into_inner().expect("ring");
        assert_eq!(ring.put(b"abcdef"), 6);
        let mut out = [0; 4];
        assert_eq!(ring.take(&mut out), 4);
        assert_eq!(&out, b"abcd");
        assert_eq!(
            ring.put(b"ghijklmn"),
            6,
            "two free at the end, four in front"
        );
        assert_eq!((ring.start, ring.len), (4, 8));
        let mut out = [0; 16];
        assert_eq!(ring.take(&mut out), 8);
        assert_eq!(&out[..8], b"efghijkl");
        assert_eq!((ring.start, ring.len), (0, 0), "empty: back at the start");
        assert_eq!(ring.take(&mut out), 0);
    }

    /// A writer whose reader went away after `room` bytes.
    struct Closing {
        room: usize,
        kind: io::ErrorKind,
    }

    impl Write for Closing {
        fn write(&mut self, data: &[u8]) -> io::Result<usize> {
            if self.room == 0 {
                return Err(io::Error::from(self.kind));
            }
            let count = data.len().min(self.room);
            self.room -= count;
            Ok(count)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_closed_output_ends_the_copy_without_an_error() {
        let closing = Closing {
            room: 10,
            kind: io::ErrorKind::BrokenPipe,
        };
        let result = copy(
            io::Cursor::new(bytes(100_000)),
            closing,
            &options(64, 8),
            |_| {},
        )
        .expect("quiet");
        assert!(result.output_closed);
        assert_eq!(result.written, 8, "whole chunks before the failed one");
    }

    #[test]
    fn another_write_error_is_an_output_error() {
        let full = Closing {
            room: 10,
            kind: io::ErrorKind::StorageFull,
        };
        let error =
            copy(io::Cursor::new(bytes(1000)), full, &options(64, 8), |_| {}).expect_err("full");
        assert!(matches!(error, Error::Output { .. }), "{error}");
        assert!(
            error.to_string().starts_with("standard output: "),
            "{error}"
        );
    }

    /// A reader that gives `first`, then waits for `gate` and fails.
    struct Gated {
        first: Option<Vec<u8>>,
        gate: mpsc::Receiver<()>,
    }

    impl Read for Gated {
        fn read(&mut self, out: &mut [u8]) -> io::Result<usize> {
            if let Some(first) = self.first.take() {
                out[..first.len()].copy_from_slice(&first);
                return Ok(first.len());
            }
            let _ = self.gate.recv();
            Err(io::Error::other("the input broke"))
        }
    }

    #[test]
    fn a_read_error_comes_after_the_bytes_read_before_it() {
        let (open, gate) = mpsc::channel();
        drop(open);
        let input = Gated {
            first: Some(b"kept".to_vec()),
            gate,
        };
        let mut out = Vec::new();
        let error = copy(input, &mut out, &options(64, 8), |_| {}).expect_err("broken input");
        assert_eq!(out, b"kept");
        assert!(matches!(error, Error::Input { .. }), "{error}");
        assert_eq!(error.to_string(), "standard input: the input broke");
    }

    #[test]
    fn a_closed_output_does_not_wait_for_an_idle_input() {
        let (open, gate) = mpsc::channel();
        let input = Gated {
            first: Some(b"some".to_vec()),
            gate,
        };
        let closed = Closing {
            room: 0,
            kind: io::ErrorKind::BrokenPipe,
        };
        // The input thread is waiting on the gate while the copy returns.
        let result = copy(input, closed, &options(64, 8), |_| {}).expect("quiet");
        assert!(result.output_closed);
        drop(open);
    }

    #[test]
    fn the_first_log_line_comes_before_the_first_read() {
        let lines = Arc::new(Mutex::new(Vec::new()));
        let sink = Arc::clone(&lines);
        let options = Options {
            debug_frequency: Some(3600),
            name: Some("lane 1".to_owned()),
            ..options(1024, 16)
        };
        let mut out = Vec::new();
        copy(io::Cursor::new(bytes(5000)), &mut out, &options, move |l| {
            sink.lock().expect("lines").push(l.to_owned());
        })
        .expect("copy");
        assert_eq!(out, bytes(5000));
        assert_eq!(
            *lines.lock().expect("lines"),
            ["Fifo buffer lane 1 used 0 / 1,024 (0%)."]
        );
    }

    #[test]
    fn a_full_buffer_is_logged_every_period() {
        // The output takes nothing until the gate opens, so the ring fills.
        struct Held {
            gate: mpsc::Receiver<()>,
            open: bool,
            bytes: Vec<u8>,
        }
        impl Write for Held {
            fn write(&mut self, data: &[u8]) -> io::Result<usize> {
                if !self.open {
                    let _ = self.gate.recv();
                    self.open = true;
                }
                self.bytes.extend_from_slice(data);
                Ok(data.len())
            }
            fn flush(&mut self) -> io::Result<()> {
                Ok(())
            }
        }
        let (open, gate) = mpsc::channel();
        let (seen, lines) = mpsc::channel();
        let options = Options {
            debug_frequency: Some(1),
            ..options(1000, 100)
        };
        let mut held = Held {
            gate,
            open: false,
            bytes: Vec::new(),
        };
        let data = bytes(5000);
        let opener = thread::spawn(move || {
            let first: String = lines.recv().expect("first line");
            let second: String = lines.recv().expect("second line");
            open.send(()).expect("open the output");
            (first, second)
        });
        copy(
            io::Cursor::new(data.clone()),
            &mut held,
            &options,
            move |l| {
                let _ = seen.send(l.to_owned());
            },
        )
        .expect("copy");
        let (first, second) = opener.join().expect("opener");
        assert_eq!(first, "Fifo buffer used 0 / 1,000 (0%).");
        assert_eq!(second, "Fifo buffer used 1,000 / 1,000 (100%).");
        assert_eq!(held.bytes, data);
    }

    /// Each expected text is what JDK 27 printed for the same pair
    /// (`DecimalFormat("#,##0")` and `NumberFormat.getPercentInstance()`).
    #[test]
    fn log_lines_are_formatted_as_the_jdk_formats_them() {
        for (used, capacity, expected) in [
            (0, 1, "0 / 1 (0%)"),
            (1, 1, "1 / 1 (100%)"),
            (1, 2, "1 / 2 (50%)"),
            (1, 8, "1 / 8 (12%)"),
            (3, 8, "3 / 8 (38%)"),
            (5, 8, "5 / 8 (62%)"),
            (7, 8, "7 / 8 (88%)"),
            (1, 3, "1 / 3 (33%)"),
            (2, 3, "2 / 3 (67%)"),
            (1, 200, "1 / 200 (0%)"),
            (3, 200, "3 / 200 (2%)"),
            (5, 200, "5 / 200 (2%)"),
            (199, 200, "199 / 200 (100%)"),
            (15, 1000, "15 / 1,000 (2%)"),
            (25, 1000, "25 / 1,000 (2%)"),
            (35, 1000, "35 / 1,000 (4%)"),
            (45, 1000, "45 / 1,000 (4%)"),
            (145, 1000, "145 / 1,000 (14%)"),
            (155, 1000, "155 / 1,000 (16%)"),
            (995, 1000, "995 / 1,000 (100%)"),
            (7, 40, "7 / 40 (18%)"),
            (9, 40, "9 / 40 (22%)"),
            (11, 40, "11 / 40 (28%)"),
            (81_920, 536_870_912, "81,920 / 536,870,912 (0%)"),
            (2_684_354, 536_870_912, "2,684,354 / 536,870,912 (0%)"),
            (2_684_355, 536_870_912, "2,684,355 / 536,870,912 (1%)"),
            (268_435_456, 536_870_912, "268,435,456 / 536,870,912 (50%)"),
            (536_870_911, 536_870_912, "536,870,911 / 536,870,912 (100%)"),
            (1_234_567, 2_147_483_645, "1,234,567 / 2,147,483,645 (0%)"),
        ] {
            assert_eq!(
                debug_line(None, used, capacity),
                format!("Fifo buffer used {expected}.")
            );
        }
        assert_eq!(
            debug_line(Some("a b"), 6, 536_870_912),
            "Fifo buffer a b used 6 / 536,870,912 (0%)."
        );
    }

    fn parsed(args: &[&str]) -> Result<Options> {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        Options::from_invocation(&TOOL.parse(&args)?)
    }

    #[test]
    fn options_follow_picard() {
        assert_eq!(parsed(&[]).expect("defaults"), Options::default());
        let given = parsed(&[
            "BUFFER_SIZE=1000",
            "IO_SIZE=+7",
            "DEBUG_FREQUENCY=2",
            "NAME=x",
        ])
        .expect("values");
        assert_eq!(
            given,
            Options {
                buffer_size: 1000,
                io_size: 7,
                debug_frequency: Some(2),
                name: Some("x".to_owned()),
            }
        );
        // Picard's `null` clears NAME; a frequency of 0 or less logs nothing.
        let quiet = parsed(&["NAME=null", "DEBUG_FREQUENCY=-1"]).expect("values");
        assert_eq!((quiet.name, quiet.debug_frequency), (None, None));
        for verbosity in ["WARNING", "ERROR"] {
            let options = parsed(&["DEBUG_FREQUENCY=1", &format!("VERBOSITY={verbosity}")]);
            assert_eq!(options.expect("values").debug_frequency, None);
        }
        let options = parsed(&["DEBUG_FREQUENCY=1", "VERBOSITY=DEBUG"]);
        assert_eq!(options.expect("values").debug_frequency, Some(1));
        let top = parsed(&["BUFFER_SIZE=2147483645"]).expect("the JVM limit");
        assert_eq!(top.buffer_size, 2_147_483_645);
    }

    #[test]
    fn bad_values_and_options_picard_lacks_are_usage_errors() {
        for (args, message) in [
            (vec!["BUFFER_SIZE=abc"], "Problem constructing Integer"),
            (vec!["BUFFER_SIZE=1.5"], "Problem constructing Integer"),
            (
                vec!["BUFFER_SIZE=3000000000"],
                "Problem constructing Integer",
            ),
            (vec!["BUFFER_SIZE=null"], "cannot be null"),
            (vec!["BUFFER_SIZE=-1"], "from 1 to 2147483645"),
            (vec!["BUFFER_SIZE=0"], "from 1 to 2147483645"),
            (vec!["BUFFER_SIZE=2147483646"], "from 1 to 2147483645"),
            (vec!["IO_SIZE=abc"], "Problem constructing Integer"),
            (vec!["IO_SIZE=0"], "from 1 to 2147483645"),
            (vec!["IO_SIZE=-1"], "from 1 to 2147483645"),
            (vec!["IO_SIZE=null"], "cannot be null"),
            (vec!["DEBUG_FREQUENCY=abc"], "Problem constructing Integer"),
            (vec!["DEBUG_FREQUENCY=null"], "cannot be null"),
            (
                vec!["VERBOSITY=BOGUS"],
                "one of ERROR, WARNING, INFO, DEBUG",
            ),
            (vec!["INPUT=x"], "Unrecognized option: INPUT"),
            (vec!["I=x"], "Unrecognized option: I"),
            (vec!["OUTPUT=x"], "Unrecognized option: OUTPUT"),
            (vec!["O=x"], "Unrecognized option: O"),
            (
                vec!["ASSUME_SORTED=true"],
                "Unrecognized option: ASSUME_SORTED",
            ),
            (vec!["GA4GH_CLIENT_SECRETS=x"], "Unrecognized option"),
            (vec!["STOP_AFTER=1"], "STOP_AFTER"),
            (vec!["BUFFER_SIZE=10", "BUFFER_SIZE=20"], "more than once"),
        ] {
            let error = parsed(&args).expect_err("refused");
            assert!(matches!(error, Error::Usage(_)), "{args:?}: {error}");
            assert!(error.to_string().contains(message), "{args:?}: {error}");
        }
    }

    #[test]
    fn the_registry_run_checks_the_options_then_refuses() {
        let args =
            |list: &[&str]| -> Vec<String> { list.iter().map(|a| (*a).to_owned()).collect() };
        let error = TOOL
            .execute(&args(&["BUFFER_SIZE=10"]))
            .expect_err("no rows");
        assert_eq!(error.to_string(), NO_ROWS);
        let error = TOOL
            .execute(&args(&["BUFFER_SIZE=0"]))
            .expect_err("bad size");
        assert!(error.to_string().starts_with("BUFFER_SIZE=0: "), "{error}");
    }
}

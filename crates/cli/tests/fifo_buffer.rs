//! `fifo-buffer` and `FifoBuffer` end to end, on real pipes: the bytes
//! written equal the bytes read for every input shape and option, both
//! faces agree, the three spellings of an option agree, a closed standard
//! output ends the copy with exit 0 and no message, the verb's note about
//! `bammetrics chain` follows `-q`, the log lines follow
//! `DEBUG_FREQUENCY`, `NAME` and `VERBOSITY`, and a bad option exits
//! before one byte is written. No file is used: every input is a pipe.

use std::io::{BufRead, BufReader, Read, Write};
use std::process::{Command, Output, Stdio};
use std::thread;

const CHAIN_NOTE: &str = "bammetrics fifo-buffer: bammetrics chain runs its steps in one process \
                          and removes the pipe this buffer serves (-q drops this line)\n";

fn command(args: &[&str]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_bammetrics"));
    command
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped());
    command
}

/// Run with `input` on standard input.
fn run(args: &[&str], input: &[u8]) -> Output {
    let mut child = command(args).spawn().expect("run bammetrics");
    let mut stdin = child.stdin.take().expect("standard input");
    let bytes = input.to_vec();
    let feeder = thread::spawn(move || {
        // A child that refused its options has closed the pipe.
        let _ = stdin.write_all(&bytes);
    });
    let output = child.wait_with_output().expect("wait");
    feeder.join().expect("feeder");
    output
}

fn text(bytes: &[u8]) -> String {
    String::from_utf8_lossy(bytes).into_owned()
}

/// Bytes with no period a buffer size divides: every value from 0 to 250.
fn pattern(count: usize) -> Vec<u8> {
    (0..count).map(|i| (i * 7 % 251) as u8).collect()
}

/// The six input shapes: nothing, one byte, text without a final newline,
/// every byte value, more than the default `IO_SIZE`, and 3 MiB.
fn shapes() -> Vec<(&'static str, Vec<u8>)> {
    vec![
        ("empty", Vec::new()),
        ("one byte", vec![0]),
        ("text", b"@r1\nACGT\n+\nIIII\n@r2\nAC".to_vec()),
        ("every byte value", (0..=255).cycle().take(1024).collect()),
        ("above IO_SIZE", pattern(70_001)),
        ("3 MiB", pattern(3 * 1024 * 1024 + 5)),
    ]
}

#[test]
fn both_faces_write_the_bytes_they_read() {
    let option_sets: [(&[&str], &[&str]); 5] = [
        (&[], &[]),
        // The input is larger than the buffer.
        (
            &["BUFFER_SIZE=1024", "IO_SIZE=100"],
            &["--buffer-size", "1024", "--io-size", "100"],
        ),
        (
            &["BUFFER_SIZE=1", "IO_SIZE=1"],
            &["--buffer-size=1", "--io-size=1"],
        ),
        // One read is larger than the buffer.
        (
            &["BUFFER_SIZE=10", "IO_SIZE=4096"],
            &["BUFFER_SIZE=10", "IO_SIZE=4096"],
        ),
        (
            &["--BUFFER_SIZE", "2147483645", "--IO_SIZE", "7"],
            &["--BUFFER_SIZE", "2147483645", "--IO_SIZE", "7"],
        ),
    ];
    for (name, bytes) in shapes() {
        for (picard, verb) in option_sets {
            // A buffer of 1 or 10 bytes costs two thread switches per fill:
            // the small shapes only.
            let tiny = picard.contains(&"BUFFER_SIZE=1") || picard.contains(&"BUFFER_SIZE=10");
            if bytes.len() > 5_000 && tiny {
                continue;
            }
            let mut args = vec!["FifoBuffer"];
            args.extend(picard);
            let face = run(&args, &bytes);
            assert_eq!(face.status.code(), Some(0), "{name} {args:?}");
            assert!(face.stdout == bytes, "{name} {args:?}: bytes differ");
            assert_eq!(text(&face.stderr), "", "{name} {args:?}");

            let mut args = vec!["fifo-buffer"];
            args.extend(verb);
            let run = run(&args, &bytes);
            assert_eq!(run.status.code(), Some(0), "{name} {args:?}");
            assert!(run.stdout == bytes, "{name} {args:?}: bytes differ");
            assert_eq!(text(&run.stderr), CHAIN_NOTE, "{name} {args:?}");
        }
    }
}

#[test]
fn the_note_about_chain_is_the_verbs_and_quiet_drops_it() {
    for quiet in [
        &["-q"][..],
        &["--quiet"],
        &["QUIET=true"],
        &["--QUIET", "true"],
    ] {
        let mut args = vec!["fifo-buffer"];
        args.extend(quiet);
        let run = run(&args, b"bytes");
        assert_eq!(run.status.code(), Some(0), "{args:?}");
        assert_eq!(run.stdout, b"bytes");
        assert_eq!(text(&run.stderr), "", "{args:?}");
    }
    let run = run(&["fifo-buffer", "--no-quiet"], b"bytes");
    assert_eq!(text(&run.stderr), CHAIN_NOTE);
    // The Picard face prints what Picard logs: nothing without DEBUG_FREQUENCY.
    let face = self::run(&["FifoBuffer"], b"bytes");
    assert_eq!(text(&face.stderr), "");
    // -v adds the Picard command line and the byte count.
    let verbose = self::run(&["fifo-buffer", "-v", "--io-size", "2"], b"bytes");
    assert_eq!(
        text(&verbose.stderr),
        format!(
            "bammetrics fifo-buffer: FifoBuffer IO_SIZE=2\n{CHAIN_NOTE}bammetrics fifo-buffer: \
             copied 5 of 5 bytes read\n"
        )
    );
}

/// Run, wait for `lines` log lines on standard error, then write `input`.
/// Returns those lines, the rest of standard error and standard output.
fn run_after_lines(args: &[&str], input: &[u8], lines: usize) -> (Vec<String>, String, Vec<u8>) {
    let mut child = command(args).spawn().expect("run bammetrics");
    let mut stderr = BufReader::new(child.stderr.take().expect("standard error"));
    let seen: Vec<String> = (0..lines)
        .map(|_| {
            let mut line = String::new();
            stderr.read_line(&mut line).expect("a log line");
            line
        })
        .collect();
    let mut stdin = child.stdin.take().expect("standard input");
    stdin.write_all(input).expect("feed");
    drop(stdin);
    let mut rest = String::new();
    stderr.read_to_string(&mut rest).expect("standard error");
    let output = child.wait_with_output().expect("wait");
    assert_eq!(output.status.code(), Some(0), "{args:?}");
    (seen, rest, output.stdout)
}

#[test]
fn log_lines_follow_debug_frequency_name_and_verbosity() {
    // The input arrives after the second line, so both read `used 0`.
    let (seen, rest, stdout) = run_after_lines(
        &[
            "FifoBuffer",
            "DEBUG_FREQUENCY=1",
            "BUFFER_SIZE=1000",
            "NAME=lane 1",
        ],
        b"late",
        2,
    );
    let line = "INFO\tFifoBuffer\tFifo buffer lane 1 used 0 / 1,000 (0%).\n";
    assert_eq!(seen, [line, line]);
    assert_eq!((rest.as_str(), stdout.as_slice()), ("", &b"late"[..]));

    let (seen, rest, stdout) =
        run_after_lines(&["fifo-buffer", "-q", "--debug-frequency", "1"], b"late", 2);
    let line = "bammetrics fifo-buffer: Fifo buffer used 0 / 536,870,912 (0%).\n";
    assert_eq!(seen, [line, line]);
    assert_eq!((rest.as_str(), stdout.as_slice()), ("", &b"late"[..]));

    for silent in [
        &["DEBUG_FREQUENCY=1", "VERBOSITY=WARNING"][..],
        &["DEBUG_FREQUENCY=1", "VERBOSITY=ERROR"],
        &["DEBUG_FREQUENCY=0"],
        &["DEBUG_FREQUENCY=-5"],
    ] {
        let mut args = vec!["FifoBuffer"];
        args.extend(silent);
        let face = run(&args, b"x");
        assert_eq!(face.status.code(), Some(0), "{args:?}");
        assert_eq!(text(&face.stderr), "", "{args:?}");
    }
    let named = run(
        &[
            "FifoBuffer",
            "DEBUG_FREQUENCY=1",
            "VERBOSITY=DEBUG",
            "NAME=null",
        ],
        b"x",
    );
    assert_eq!(
        text(&named.stderr),
        "INFO\tFifoBuffer\tFifo buffer used 0 / 536,870,912 (0%).\n"
    );
}

/// A full buffer is logged as full: nothing reads standard output until
/// the third log line, then all of it is read.
#[test]
fn a_reader_that_waits_sees_a_full_buffer_logged_then_every_byte() {
    let bytes = pattern(2_000_000);
    let mut child = command(&[
        "FifoBuffer",
        "BUFFER_SIZE=1000",
        "IO_SIZE=100",
        "DEBUG_FREQUENCY=1",
    ])
    .spawn()
    .expect("run bammetrics");
    let mut stdin = child.stdin.take().expect("standard input");
    let input = bytes.clone();
    let feeder = thread::spawn(move || stdin.write_all(&input).expect("feed"));
    let mut stderr = BufReader::new(child.stderr.take().expect("standard error"));
    let seen: Vec<String> = (0..3)
        .map(|_| {
            let mut line = String::new();
            stderr.read_line(&mut line).expect("a log line");
            line
        })
        .collect();
    let output = child.wait_with_output().expect("wait");
    feeder.join().expect("feeder");
    assert_eq!(output.status.code(), Some(0));
    assert!(output.stdout == bytes, "bytes differ");
    let full = "INFO\tFifoBuffer\tFifo buffer used 1,000 / 1,000 (100%).\n";
    assert_eq!(
        seen,
        [
            "INFO\tFifoBuffer\tFifo buffer used 0 / 1,000 (0%).\n",
            full,
            full
        ]
    );
}

/// The reader of standard output takes 10 bytes and closes: both faces
/// exit 0 with nothing on standard error, and the writer of standard
/// input sees its pipe close instead of writing 64 MiB.
#[test]
fn a_closed_standard_output_ends_the_copy_quietly() {
    // The buffer is smaller than the input, so the reader thread cannot
    // take all of it before the writer meets the closed pipe.
    for args in [
        &["FifoBuffer", "BUFFER_SIZE=65536"][..],
        &["fifo-buffer", "-q", "--buffer-size", "65536"],
    ] {
        let mut child = command(args).spawn().expect("run bammetrics");
        let mut stdin = child.stdin.take().expect("standard input");
        let feeder = thread::spawn(move || {
            let chunk = vec![7u8; 65_536];
            (0..1024)
                .take_while(|_| stdin.write_all(&chunk).is_ok())
                .count()
        });
        let mut stdout = child.stdout.take().expect("standard output");
        let mut first = [0u8; 10];
        stdout.read_exact(&mut first).expect("ten bytes");
        drop(stdout);
        let output = child.wait_with_output().expect("wait");
        let chunks = feeder.join().expect("feeder");
        assert_eq!(first, [7u8; 10]);
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        assert_eq!(text(&output.stderr), "", "{args:?}");
        assert!(chunks < 1024, "{args:?}: the whole input was read");
    }
}

#[test]
fn three_spellings_of_an_option_copy_the_same_bytes() {
    let bytes = pattern(50_000);
    for spelling in [
        &["--buffer-size", "333"][..],
        &["--BUFFER_SIZE", "333"],
        &["BUFFER_SIZE=333"],
    ] {
        let mut args = vec!["fifo-buffer", "-v"];
        args.extend(spelling);
        let run = run(&args, &bytes);
        assert_eq!(run.status.code(), Some(0), "{args:?}");
        assert!(run.stdout == bytes, "{args:?}: bytes differ");
        assert!(
            text(&run.stderr).starts_with("bammetrics fifo-buffer: FifoBuffer BUFFER_SIZE=333\n"),
            "{args:?}: {}",
            text(&run.stderr)
        );
    }
}

/// Picard exits 1 for each of these (it waits forever on a size of 0).
/// The Picard face exits 1 and the verb 2, with no byte on standard
/// output.
#[test]
fn bad_options_exit_before_a_byte_is_copied() {
    let cases: [(&[&str], &[&str], &str); 12] = [
        (
            &["BUFFER_SIZE=abc"],
            &["--buffer-size", "abc"],
            "Problem constructing Integer from the string 'abc'.",
        ),
        (
            &["BUFFER_SIZE=3000000000"],
            &["--buffer-size", "3000000000"],
            "Problem constructing Integer",
        ),
        (
            &["BUFFER_SIZE=-1"],
            &["--buffer-size=-1"],
            "a size from 1 to 2147483645 bytes is needed",
        ),
        (
            &["BUFFER_SIZE=0"],
            &["--buffer-size", "0"],
            "a size from 1 to 2147483645 bytes is needed",
        ),
        (
            &["BUFFER_SIZE=2147483647"],
            &["--buffer-size", "2147483647"],
            "a size from 1 to 2147483645 bytes is needed",
        ),
        (
            &["BUFFER_SIZE=null"],
            &["--buffer-size", "null"],
            "an int option cannot be null",
        ),
        (
            &["IO_SIZE=0"],
            &["--io-size", "0"],
            "a size from 1 to 2147483645 bytes is needed",
        ),
        (
            &["IO_SIZE=1.5"],
            &["--io-size", "1.5"],
            "Problem constructing Integer",
        ),
        (
            &["DEBUG_FREQUENCY=soon"],
            &["--debug-frequency", "soon"],
            "Problem constructing Integer",
        ),
        (
            &["VERBOSITY=LOUD"],
            &["--verbosity", "LOUD"],
            "one of ERROR, WARNING, INFO, DEBUG",
        ),
        (&["OUTPUT=x"], &["-o", "x"], ""),
        (&["INPUT=x"], &["x"], ""),
    ];
    for (picard, verb, message) in cases {
        let mut args = vec!["FifoBuffer"];
        args.extend(picard);
        let face = run(&args, b"never copied");
        assert_eq!(face.status.code(), Some(1), "{args:?}");
        assert_eq!(face.stdout, b"", "{args:?}");
        assert!(text(&face.stderr).contains(message), "{args:?}: {face:?}");

        let mut args = vec!["fifo-buffer"];
        args.extend(verb);
        let run = run(&args, b"never copied");
        assert_eq!(run.status.code(), Some(2), "{args:?}");
        assert_eq!(run.stdout, b"", "{args:?}");
        assert!(text(&run.stderr).contains(message), "{args:?}: {run:?}");
        assert!(
            !text(&run.stderr).contains("bammetrics chain"),
            "{args:?}: the note comes after the options are read"
        );
    }
    let face = run(&["FifoBuffer", "O=x"], b"");
    assert_eq!(
        text(&face.stderr),
        "bammetrics: FifoBuffer: Unrecognized option: O\n"
    );
    let twice = run(&["FifoBuffer", "IO_SIZE=1", "IO_SIZE=2"], b"");
    assert_eq!(twice.status.code(), Some(1));
    assert_eq!(
        text(&twice.stderr),
        "bammetrics: Option 'IO_SIZE' cannot be specified more than once.\n"
    );
    // The verb keeps the last value of an option given twice, as every verb does.
    let last = run(
        &["fifo-buffer", "-q", "--io-size", "1", "--io-size", "2"],
        b"ab",
    );
    assert_eq!(
        (last.status.code(), last.stdout.as_slice()),
        (Some(0), &b"ab"[..])
    );
}

/// The options Picard's FifoBuffer takes from `CommandLineProgram` are
/// accepted and noted; the copy is the same.
#[test]
fn common_options_are_accepted_and_change_nothing() {
    let face = run(
        &[
            "FifoBuffer",
            "VALIDATION_STRINGENCY=SILENT",
            "TMP_DIR=/tmp",
            "COMPRESSION_LEVEL=1",
            "MAX_RECORDS_IN_RAM=5",
            "CREATE_INDEX=true",
            "CREATE_MD5_FILE=true",
            "R=no_such.fa",
            "USE_JDK_DEFLATER=true",
            "QUIET=false",
        ],
        b"same",
    );
    assert_eq!(face.status.code(), Some(0), "{face:?}");
    assert_eq!(face.stdout, b"same");
    let notes = text(&face.stderr);
    assert_eq!(notes.lines().count(), 9, "{notes}");
    assert!(
        notes
            .lines()
            .all(|line| line.starts_with("bammetrics: ignoring ")),
        "{notes}"
    );
}

#[test]
fn a_closed_standard_input_is_an_empty_input() {
    for args in [&["FifoBuffer"][..], &["fifo-buffer", "-q"]] {
        let output = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args(args)
            .stdin(Stdio::null())
            .output()
            .expect("run bammetrics");
        assert_eq!(output.status.code(), Some(0), "{args:?}");
        assert_eq!(output.stdout, b"", "{args:?}");
        assert_eq!(text(&output.stderr), "", "{args:?}");
    }
}

#[test]
fn help_names_the_streams_and_no_output_option() {
    let help = run(&["fifo-buffer", "--help"], b"");
    let page = text(&help.stdout);
    assert!(
        page.contains("Usage: bammetrics fifo-buffer [OPTIONS] < IN > OUT\n"),
        "{page}"
    );
    for flag in [
        "--buffer-size INT",
        "--io-size INT",
        "--debug-frequency INT",
        "--name TEXT",
    ] {
        assert!(page.contains(flag), "{flag}: {page}");
    }
    for absent in ["--output", "--format", "--sample", "--input"] {
        assert!(!page.contains(absent), "{absent}: {page}");
    }
    let face = run(&["FifoBuffer", "--help"], b"");
    assert_eq!(text(&face.stdout), page);
    let by_help = run(&["help", "FifoBuffer"], b"");
    assert_eq!(text(&by_help.stdout), page);
}

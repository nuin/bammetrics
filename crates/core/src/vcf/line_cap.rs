//! A `BufRead` that fails once a line grows past a byte cap or the lines
//! past a count cap, so a VCF reader never holds more than one capped line
//! and one buffer.

use std::io::{self, BufRead, Read};

pub(super) struct LineCap {
    inner: Box<dyn BufRead>,
    max_line_bytes: usize,
    max_lines: u64,
    /// Bytes handed out since the last line feed.
    line_bytes: usize,
    /// Line feeds handed out.
    lines: u64,
    /// The longest line handed out, its line feed not counted.
    longest: usize,
    max_header_lines: u64,
    /// Lines whose first byte is `#`.
    header_lines: u64,
    /// True when the next byte handed out starts a line.
    at_line_start: bool,
}

impl LineCap {
    pub(super) fn new(inner: Box<dyn BufRead>, max_line_bytes: usize, max_lines: u64) -> Self {
        Self {
            inner,
            max_line_bytes,
            max_lines,
            line_bytes: 0,
            lines: 0,
            longest: 0,
            max_header_lines: u64::MAX,
            header_lines: 0,
            at_line_start: true,
        }
    }

    /// Fail too once more than `max` lines start with `#`.
    pub(super) fn with_header_lines(mut self, max: u64) -> Self {
        self.max_header_lines = max;
        self
    }

    fn check(&self) -> io::Result<()> {
        if self.longest > self.max_line_bytes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("a line is over {} bytes", self.max_line_bytes),
            ));
        }
        if self.lines > self.max_lines {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("over {} lines", self.max_lines),
            ));
        }
        if self.header_lines > self.max_header_lines {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("over {} header lines", self.max_header_lines),
            ));
        }
        Ok(())
    }
}

impl Read for LineCap {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let n = available.len().min(buf.len());
        buf[..n].copy_from_slice(&available[..n]);
        self.consume(n);
        Ok(n)
    }
}

impl BufRead for LineCap {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.check()?;
        self.inner.fill_buf()
    }

    fn consume(&mut self, amount: usize) {
        if let Ok(buffer) = self.inner.fill_buf() {
            let taken = &buffer[..amount.min(buffer.len())];
            for (at, segment) in taken.split(|&b| b == b'\n').enumerate() {
                if at == 0 {
                    self.line_bytes = self.line_bytes.saturating_add(segment.len());
                } else {
                    self.lines = self.lines.saturating_add(1);
                    self.line_bytes = segment.len();
                    self.at_line_start = true;
                }
                self.longest = self.longest.max(self.line_bytes);
                if self.at_line_start
                    && let Some(&first) = segment.first()
                {
                    self.at_line_start = false;
                    if first == b'#' {
                        self.header_lines = self.header_lines.saturating_add(1);
                    }
                }
            }
        }
        self.inner.consume(amount);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn capped(text: &'static [u8], max_line_bytes: usize, max_lines: u64) -> LineCap {
        LineCap::new(Box::new(text), max_line_bytes, max_lines)
    }

    #[test]
    fn lines_within_the_caps_pass_and_others_fail() {
        let mut out = String::new();
        capped(b"ab\ncd\n", 2, 2).read_to_string(&mut out).unwrap();
        assert_eq!(out, "ab\ncd\n");
        let mut lines = capped(b"abc\nd\n", 2, 9).lines();
        assert!(lines.next().unwrap().is_ok());
        assert!(lines.next().unwrap().is_err());
        let mut out = String::new();
        let over = capped(b"a\nb\nc\nd", 9, 2).read_to_string(&mut out);
        assert!(over.is_err());
    }

    #[test]
    fn header_lines_are_counted_across_reads() {
        let text: &'static [u8] = b"##a\n##b\n#c\nx#\ny\n";
        let mut out = String::new();
        capped(text, 9, 9)
            .with_header_lines(3)
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out.as_bytes(), text);
        let mut out = String::new();
        let over = capped(text, 9, 9)
            .with_header_lines(2)
            .read_to_string(&mut out);
        assert!(over.is_err());
        // One byte at a time: a line start split from its first byte.
        let mut capped = capped(text, 9, 9).with_header_lines(2);
        let mut byte = [0u8; 1];
        let mut failed = false;
        for _ in 0..text.len() {
            if capped.read(&mut byte).is_err() {
                failed = true;
                break;
            }
        }
        assert!(failed);
    }
}

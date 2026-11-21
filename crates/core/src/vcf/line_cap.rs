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
        }
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
            let mut segments = taken.split(|&b| b == b'\n');
            if let Some(first) = segments.next() {
                self.line_bytes = self.line_bytes.saturating_add(first.len());
                self.longest = self.longest.max(self.line_bytes);
            }
            for segment in segments {
                self.lines = self.lines.saturating_add(1);
                self.line_bytes = segment.len();
                self.longest = self.longest.max(self.line_bytes);
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
}

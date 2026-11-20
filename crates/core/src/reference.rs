//! Reference bases one contig at a time, as Picard's
//! `ReferenceSequenceFileWalker` hands them to `SinglePassSamProgram`.

use std::fs::File;
use std::io::{self, BufRead, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam;

use crate::error::{Error, Result};

type Reader = fasta::io::IndexedReader<fasta::io::BufReader<File>>;

/// The files a reader of the FASTA `fasta` opens: the FASTA and its
/// `.fai`. The CRAM decoder and [`Reference::open`] read no `.dict` and no
/// `.gzi`. A tool refuses to write any of them.
pub fn files_read(fasta: &Path) -> [PathBuf; 2] {
    let mut fai = fasta.as_os_str().to_owned();
    fai.push(".fai");
    [fasta.to_path_buf(), PathBuf::from(fai)]
}

/// A FASTA with its `.fai`, read by random access. Only the contig of the
/// current record is held in memory.
pub struct Reference {
    path: PathBuf,
    reader: Reader,
    /// The FASTA's size in bytes: the most a contig buffer reserves.
    size: u64,
    current: Option<(usize, Vec<u8>)>,
}

impl Reference {
    /// Open `path` and its `path.fai`. When the BAM header has `@SQ` lines,
    /// they must name the same contigs with the same lengths in the same
    /// order (`SequenceDictionaryUtils.assertSequenceDictionariesEqual`).
    pub fn open(path: &Path, header: &sam::Header) -> Result<Self> {
        let error = |message: String| Error::Reference {
            path: path.to_path_buf(),
            message,
        };
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| error(format!("cannot open with its .fai index: {e}")))?;
        let contigs = header.reference_sequences();
        if !contigs.is_empty() {
            let index = reader.index().as_ref();
            if index.len() != contigs.len() {
                return Err(error(format!(
                    "{} contigs in the reference, {} in the BAM header",
                    index.len(),
                    contigs.len()
                )));
            }
            for (record, (name, contig)) in index.iter().zip(contigs) {
                let length = contig.length().get();
                if record.name() != name.as_slice() || record.length() != length as u64 {
                    return Err(error(format!(
                        "contig {} ({} bp) in the reference, {} ({length} bp) in the BAM header",
                        record.name(),
                        record.length(),
                        name
                    )));
                }
            }
        }
        Ok(Self {
            path: path.to_path_buf(),
            reader,
            size: std::fs::metadata(path).map_or(0, |m| m.len()),
            current: None,
        })
    }

    /// Bases of the contig at `id` in the reference dictionary.
    /// `ReferenceSequenceFileWalker.get` only moves forward, so a record on
    /// an earlier contig than the one before it is an error, as in Picard.
    pub fn contig(&mut self, id: usize) -> Result<&[u8]> {
        match self.current.as_ref().map(|(i, _)| *i) {
            Some(previous) if id < previous => {
                return Err(self.error(format!(
                    "a record on contig {id} follows one on contig {previous}; the input must be coordinate sorted"
                )));
            }
            Some(previous) if id == previous => {}
            _ => {
                // The buffer of the contig left is refilled: one allocation
                // for the walk, at the largest contig's size.
                let mut bases = self.current.take().map(|(_, b)| b).unwrap_or_default();
                self.read_contig(id, &mut bases)?;
                self.current = Some((id, bases));
            }
        }
        Ok(self.current.as_ref().map_or(&[], |(_, bases)| bases))
    }

    /// Bases `start..=end` (1-based) of the contig at `id`, read through the
    /// `.fai` without loading the contig (`ReferenceSequenceFile.getSubsequenceAt`).
    /// It does not move the forward-only `contig` walker.
    pub fn subsequence(&mut self, id: usize, start: i64, end: i64) -> Result<Vec<u8>> {
        if end < start {
            return Ok(Vec::new());
        }
        let position = |p: i64| {
            usize::try_from(p)
                .ok()
                .and_then(Position::new)
                .ok_or_else(|| self.error(format!("bad position {p}")))
        };
        let region = Region::new(self.name(id)?, position(start)?..=position(end)?);
        self.query(&region)
    }

    /// Every base of the contig at `id`, read fresh through the `.fai`.
    /// It does not move the forward-only `contig` walker.
    pub fn load(&mut self, id: usize) -> Result<Vec<u8>> {
        let mut bases = Vec::new();
        self.read_contig(id, &mut bases)?;
        Ok(bases)
    }

    /// Every base of the contig at `id` into `bases`, cleared first: the
    /// bases noodles' query of the whole contig returns, read into a buffer
    /// reserved at the `.fai` length (at most the FASTA's size) with no
    /// copy. Picard's walker holds one contig too.
    fn read_contig(&mut self, id: usize, bases: &mut Vec<u8>) -> Result<()> {
        let region = Region::new(self.name(id)?, ..);
        let index = self.reader.index();
        let length = index.as_ref().get(id).map_or(0, |r| r.length());
        let offset = index
            .query(&region)
            .map_err(|e| self.error(e.to_string()))?;
        bases.clear();
        bases.reserve(usize::try_from(length.min(self.size)).unwrap_or(0));
        read_bases(self.reader.get_mut(), offset, bases).map_err(|e| self.error(e.to_string()))
    }

    /// Number of contigs in the `.fai`.
    pub fn len(&self) -> usize {
        self.reader.index().as_ref().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn name(&self, id: usize) -> Result<Vec<u8>> {
        self.reader
            .index()
            .as_ref()
            .get(id)
            .map(|r| r.name().to_vec())
            .ok_or_else(|| self.error(format!("no contig {id} in the .fai index")))
    }

    fn query(&mut self, region: &Region) -> Result<Vec<u8>> {
        let record = self
            .reader
            .query(region)
            .map_err(|e: io::Error| self.error(e.to_string()))?;
        Ok(record.sequence().as_ref().to_vec())
    }

    fn error(&self, message: String) -> Error {
        Error::Reference {
            path: self.path.clone(),
            message,
        }
    }
}

/// The sequence lines from `offset` to the next `>` or the end of the
/// file, as noodles' sequence reader takes them: empty lines skipped, the
/// line end and a `\r` before it dropped.
fn read_bases<R: BufRead + Seek>(
    reader: &mut R,
    offset: u64,
    bases: &mut Vec<u8>,
) -> io::Result<()> {
    reader.seek(SeekFrom::Start(offset))?;
    loop {
        loop {
            let mut newline = false;
            if reader.fill_buf()?.starts_with(b"\r") {
                reader.consume(1);
                newline = true;
            }
            if reader.fill_buf()?.starts_with(b"\n") {
                reader.consume(1);
                newline = true;
            }
            if !newline {
                break;
            }
        }
        let src = reader.fill_buf()?;
        if src.first().is_none_or(|&b| b == b'>') {
            return Ok(());
        }
        let line = src
            .iter()
            .position(|&b| b == b'\n')
            .map_or(src, |end| &src[..end]);
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        bases.extend_from_slice(line);
        let used = line.len();
        reader.consume(used);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// chr1 = ACGTACGTAC, chr2 = GGGG, one line each, with a `.fai`, in a
    /// directory of its own per test so parallel tests do not race.
    fn two_contigs(test: &str) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("bammetrics-ref-{}-{test}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("ref.fa");
        std::fs::write(&path, ">chr1\nACGTACGTAC\n>chr2\nGGGG\n").expect("write FASTA");
        std::fs::write(
            dir.join("ref.fa.fai"),
            "chr1\t10\t6\t10\t11\nchr2\t4\t23\t4\t5\n",
        )
        .expect("write .fai");
        path
    }

    #[test]
    fn contigs_load_forward_only() {
        let mut reference =
            Reference::open(&two_contigs("walk"), &sam::Header::default()).expect("opens");
        assert_eq!(reference.contig(0).expect("chr1"), b"ACGTACGTAC");
        assert_eq!(reference.contig(1).expect("chr2"), b"GGGG");
        assert_eq!(reference.contig(1).expect("chr2 again"), b"GGGG");
        assert!(matches!(
            reference.contig(0),
            Err(Error::Reference { message, .. }) if message.contains("coordinate sorted")
        ));
        assert!(reference.contig(5).is_err());
    }

    #[test]
    fn subsequences_are_one_based_and_closed() {
        let mut reference =
            Reference::open(&two_contigs("region"), &sam::Header::default()).expect("opens");
        assert_eq!(reference.subsequence(0, 2, 5).expect("chr1"), b"CGTA");
        assert_eq!(reference.subsequence(1, 4, 4).expect("chr2"), b"G");
        assert!(reference.subsequence(1, 3, 2).expect("empty").is_empty());
        assert_eq!(
            reference.contig(0).expect("walker untouched"),
            b"ACGTACGTAC"
        );
    }

    /// The walker's bases are noodles' query of the whole contig on every
    /// line form: CRLF, empty lines, a lone `\r`, lower case, a `.fai`
    /// length shorter or longer than the lines, no newline at the end.
    #[test]
    fn contigs_read_as_noodles_reads_them() {
        let dir = std::env::temp_dir().join(format!("bammetrics-ref-{}-lines", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("ref.fa");
        let fasta = ">a\r\nACgt\r\n\r\nNNac\r\n>b desc\nAC\rGT\n\n\nTT\n>c\nGGG";
        std::fs::write(&path, fasta).expect("write FASTA");
        let offset = |name: &str| fasta.find(name).expect("name") + name.len() + 1;
        let fai = format!(
            "a\t5\t{}\t4\t6\nb\t20\t{}\t5\t6\nc\t3\t{}\t3\t4\n",
            offset(">a\r"),
            offset(">b desc"),
            offset(">c"),
        );
        std::fs::write(dir.join("ref.fa.fai"), fai).expect("write .fai");
        let mut reference = Reference::open(&path, &sam::Header::default()).expect("opens");
        let mut noodles = fasta::io::indexed_reader::Builder::default()
            .build_from_path(&path)
            .expect("noodles");
        for (id, name) in ["a", "b", "c"].into_iter().enumerate() {
            let region = Region::new(name, ..);
            let expected = noodles.query(&region).expect("query");
            let expected = expected.sequence().as_ref().to_vec();
            assert_eq!(reference.contig(id).expect("walk"), expected, "{name}");
            assert_eq!(reference.load(id).expect("load"), expected, "{name}");
        }
        assert_eq!(reference.load(0).expect("a"), b"ACgtNNac");
        assert_eq!(reference.load(1).expect("b"), b"AC\rGTTT");
        assert_eq!(reference.load(2).expect("c"), b"GGG");
        std::fs::remove_dir_all(&dir).expect("cleanup");
    }
}

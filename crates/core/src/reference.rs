//! Reference bases one contig at a time, as Picard's
//! `ReferenceSequenceFileWalker` hands them to `SinglePassSamProgram`.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use noodles::core::Position;
use noodles::core::region::Interval;
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
        let reference = Self::open_unchecked(path)?;
        let error = |message: String| Error::Reference {
            path: path.to_path_buf(),
            message,
        };
        let contigs = header.reference_sequences();
        if !contigs.is_empty() {
            let index = reference.reader.index().as_ref();
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
        Ok(reference)
    }

    /// Open `path` and its `path.fai` with no check against a header: the
    /// caller finds each contig it needs ([`Reference::places`]).
    pub fn open_unchecked(path: &Path) -> Result<Self> {
        let reader = fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| Error::Reference {
                path: path.to_path_buf(),
                message: format!("cannot open with its .fai index: {e}"),
            })?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
            size: std::fs::metadata(path).map_or(0, |m| m.len()),
            current: None,
        })
    }

    /// Where the `.fai` has each reference sequence of `header`, found by
    /// name: its index there and its length there, `None` for a name the
    /// `.fai` lacks. The first line of a name counts. One pass over the
    /// `.fai` and one over the header.
    pub fn places(&self, header: &sam::Header) -> Vec<Option<(usize, u64)>> {
        let index = self.reader.index().as_ref();
        let mut by_name: HashMap<&[u8], usize> = HashMap::with_capacity(index.len());
        for (id, record) in index.iter().enumerate() {
            by_name.entry(record.name().as_ref()).or_insert(id);
        }
        header
            .reference_sequences()
            .keys()
            .map(|name| {
                let id = *by_name.get(AsRef::<[u8]>::as_ref(name))?;
                Some((id, index[id].length()))
            })
            .collect()
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
    /// `.fai` without loading the contig (`ReferenceSequenceFile.getSubsequenceAt`):
    /// fewer where the contig's lines end before `end`. The `.fai` line is
    /// taken by `id`, with no search by name. It does not move the
    /// forward-only `contig` walker.
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
        let (start, end) = (position(start)?, position(end)?);
        let record = self.record(id)?;
        if usize::from(start) as u64 > record.length() {
            // The offset of such a start lies in the lines of a later contig.
            return Ok(Vec::new());
        }
        let offset = record
            .query(Interval::from(start..=end))
            .map_err(|e| self.error(e.to_string()))?;
        let limit = usize::from(end) - usize::from(start) + 1;
        // Nothing is reserved: `limit` comes from a record's CIGAR.
        let mut bases = Vec::new();
        read_bases(self.reader.get_mut(), offset, limit, &mut bases)
            .map_err(|e| self.error(e.to_string()))?;
        Ok(bases)
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
    /// copy. Picard's walker holds one contig too. The `.fai` line is taken
    /// by `id`: noodles' query by name scans the index, which made a walk
    /// over every contig quadratic in their number.
    fn read_contig(&mut self, id: usize, bases: &mut Vec<u8>) -> Result<()> {
        let record = self.record(id)?;
        let (length, offset) = (record.length(), record.position());
        bases.clear();
        bases.reserve(usize::try_from(length.min(self.size)).unwrap_or(0));
        read_bases(self.reader.get_mut(), offset, usize::MAX, bases)
            .map_err(|e| self.error(e.to_string()))
    }

    /// Number of contigs in the `.fai`.
    pub fn len(&self) -> usize {
        self.reader.index().as_ref().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// The `.fai` line at `id`.
    fn record(&self, id: usize) -> Result<&fasta::fai::Record> {
        self.reader
            .index()
            .as_ref()
            .get(id)
            .ok_or_else(|| self.error(format!("no contig {id} in the .fai index")))
    }

    /// An error that names this FASTA.
    pub(crate) fn error(&self, message: String) -> Error {
        Error::Reference {
            path: self.path.clone(),
            message,
        }
    }
}

/// The sequence lines from `offset` to the next `>` or the end of the
/// file, as noodles' sequence reader takes them: empty lines skipped, the
/// line end and a `\r` before it dropped. At most `limit` bases are added
/// to `bases`.
fn read_bases<R: BufRead + Seek>(
    reader: &mut R,
    offset: u64,
    limit: usize,
    bases: &mut Vec<u8>,
) -> io::Result<()> {
    reader.seek(SeekFrom::Start(offset))?;
    let end = bases.len().saturating_add(limit);
    loop {
        if bases.len() >= end {
            return Ok(());
        }
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
        let line = &line[..line.len().min(end - bases.len())];
        bases.extend_from_slice(line);
        let used = line.len();
        reader.consume(used);
    }
}

#[cfg(test)]
mod tests {
    use noodles::core::Region;

    use super::*;

    /// A contig is read through its own `.fai` line, with no search by
    /// name. 100,000 contigs, each loaded whole and read in part: noodles'
    /// query by name took 176 s in a debug build on the machine this was
    /// written on, and reading by index 0.3 s.
    #[test]
    fn a_walk_over_100_000_contigs_is_linear() {
        const CONTIGS: usize = 100_000;
        let dir = std::env::temp_dir().join(format!("bammetrics-ref-{}-many", std::process::id()));
        std::fs::create_dir_all(&dir).expect("temp dir");
        let path = dir.join("ref.fa");
        let (mut fasta, mut fai) = (String::new(), String::new());
        for id in 0..CONTIGS {
            fasta.push_str(&format!(">c{id}\n"));
            fai.push_str(&format!("c{id}\t8\t{}\t8\t9\n", fasta.len()));
            fasta.push_str("ACGTACGT\n");
        }
        std::fs::write(&path, fasta).expect("write FASTA");
        std::fs::write(dir.join("ref.fa.fai"), fai).expect("write .fai");
        let mut reference = Reference::open_unchecked(&path).expect("opens");
        let started = std::time::Instant::now();
        for id in 0..CONTIGS {
            assert_eq!(reference.load(id).expect("load"), b"ACGTACGT");
            assert_eq!(reference.subsequence(id, 3, 6).expect("part"), b"GTAC");
        }
        let took = started.elapsed();
        std::fs::remove_dir_all(&dir).expect("cleanup");
        assert!(took.as_secs() < 10, "{took:?} for {CONTIGS} contigs");
    }

    /// The contigs of a header are found by name: a FASTA with more
    /// contigs, in another order, of another length, or without one.
    #[test]
    fn places_are_found_by_name() {
        let reference = Reference::open_unchecked(&two_contigs("places")).expect("opens");
        let header: sam::Header = "@SQ\tSN:chr2\tLN:9\n@SQ\tSN:chrX\tLN:5\n@SQ\tSN:chr1\tLN:10\n"
            .parse()
            .expect("header");
        assert_eq!(
            reference.places(&header),
            [Some((1, 4)), None, Some((0, 10))]
        );
        assert!(Reference::open(&two_contigs("strict"), &header).is_err());
    }

    /// A part of a contig stops at the contig's last line, and at the end
    /// asked for inside a line.
    #[test]
    fn a_subsequence_ends_with_its_contig() {
        let mut reference = Reference::open_unchecked(&two_contigs("ends")).expect("opens");
        assert_eq!(reference.subsequence(0, 9, 40).expect("chr1"), b"AC");
        assert_eq!(reference.subsequence(1, 2, 3).expect("chr2"), b"GG");
        assert!(reference.subsequence(0, 11, 14).expect("past").is_empty());
        assert!(reference.subsequence(2, 1, 4).is_err());
    }

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

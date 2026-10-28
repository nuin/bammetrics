//! Reference bases one contig at a time, as Picard's
//! `ReferenceSequenceFileWalker` hands them to `SinglePassSamProgram`.

use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

use noodles::core::{Position, Region};
use noodles::fasta;
use noodles::sam;

use crate::error::{Error, Result};

type Reader = fasta::io::IndexedReader<fasta::io::BufReader<File>>;

/// A FASTA with its `.fai`, read by random access. Only the contig of the
/// current record is held in memory.
pub struct Reference {
    path: PathBuf,
    reader: Reader,
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
                let bases = self.load(id)?;
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

    fn load(&mut self, id: usize) -> Result<Vec<u8>> {
        let region = Region::new(self.name(id)?, ..);
        self.query(&region)
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
}

//! The record edit of Picard 3.5.0 `MakeSitesOnlyVcf`: each record read
//! as `VCFCodec` reads it, then cut to the samples of `SAMPLE`
//! (`subsetToSamplesWithOriginalAnnotations`: the alleles and `INFO` as
//! read, `GenotypesContext.subsetToSamples` over the genotypes).
//!
//! - No sample (the default): the eight site columns, the genotype text
//!   left unread unless the input header's samples are out of sorted
//!   order, when the codec decodes it as it reads the record.
//! - Samples: the header names them in `TreeSet` order (`String.compareTo`,
//!   each once). Every genotype of the record is decoded (the lazy context
//!   decodes all of them on the first lookup), and the genotypes of the
//!   named samples are written anew as `VCFEncoder` writes them; a name
//!   the input lacks gets a no-call at the record's largest ploidy.
//!
//! The site columns are handed on as read: the writer encodes them.

use std::collections::BTreeMap;
use std::sync::Arc;

use super::encode_genotypes::Subset;
use super::full_header::{VcfHeader, java_compare};
use super::site::{self, Rules, split_condensed};

/// A name ordered as Java's `TreeSet<String>` orders it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct JavaName(String);

impl PartialOrd for JavaName {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for JavaName {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        java_compare(&self.0, &other.0)
    }
}

/// `SAMPLE` as Picard's `TreeSet<String>` holds it: each name once, in
/// `String.compareTo` order.
pub fn tree_set(names: &[String]) -> Vec<String> {
    let set: BTreeMap<JavaName, ()> = names.iter().map(|n| (JavaName(n.clone()), ())).collect();
    set.into_keys().map(|JavaName(name)| name).collect()
}

/// The edit of every record of one header.
#[derive(Debug)]
pub struct SampleSubset {
    rules: Rules,
    percent: bool,
    output: Arc<VcfHeader>,
    /// Per output sample, its column in the input; `None` when the input
    /// lacks it.
    picks: Vec<Option<usize>>,
}

impl SampleSubset {
    /// The edit of records read under `input`, keeping `samples` (any
    /// order, repeats allowed: they are read as a `TreeSet`).
    pub fn new(input: &VcfHeader, samples: &[String]) -> Self {
        let kept = tree_set(samples);
        let column: std::collections::HashMap<&str, usize> = input
            .samples()
            .iter()
            .enumerate()
            .map(|(i, name)| (name.as_str(), i))
            .collect();
        let picks = kept
            .iter()
            .map(|name| column.get(name.as_str()).copied())
            .collect();
        Self {
            rules: Rules::of(input.site()),
            percent: input.minor() >= 3,
            output: Arc::new(input.with_samples(kept)),
            picks,
        }
    }

    /// The header the edited records are written under.
    pub fn header(&self) -> &VcfHeader {
        &self.output
    }

    /// [`SampleSubset::header`], one copy shared by every holder.
    pub fn shared_header(&self) -> Arc<VcfHeader> {
        Arc::clone(&self.output)
    }

    /// Replace `out` by the edited `line`: false for a line that starts
    /// with `#`, which the codec skips. An error is the codec's or the
    /// encoder's message.
    pub fn edit(&self, line: &str, out: &mut String) -> Result<bool, String> {
        out.clear();
        let Some((site, alleles, genotypes)) = site::parse_fixed(line, &self.rules)? else {
            return Ok(false);
        };
        if self.picks.is_empty() && self.rules.eager_genotypes {
            super::encode_genotypes::check(
                genotypes,
                self.rules.samples,
                alleles,
                self.rules.percent_decoded,
            )?;
        }
        for (i, part) in split_condensed(line, 9).iter().take(8).enumerate() {
            if i > 0 {
                out.push('\t');
            }
            out.push_str(part);
        }
        if self.picks.is_empty() {
            return Ok(true);
        }
        out.push('\t');
        Subset {
            samples: self.rules.samples,
            percent: self.percent,
            output: &self.output,
            picks: &self.picks,
        }
        .write(genotypes, alleles, &site, out)?;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(samples: &[&str]) -> VcfHeader {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"q\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        );
        if !samples.is_empty() {
            text.push_str("\tFORMAT");
            for s in samples {
                text.push('\t');
                text.push_str(s);
            }
        }
        text.push('\n');
        VcfHeader::read(&mut text.as_bytes()).expect("a header")
    }

    fn edit(input: &[&str], keep: &[&str], line: &str) -> Result<(Vec<String>, String), String> {
        let h = header(input);
        let keep: Vec<String> = keep.iter().map(|s| (*s).to_owned()).collect();
        let subset = SampleSubset::new(&h, &keep);
        let mut out = String::new();
        subset.edit(line, &mut out)?;
        Ok((subset.header().samples().to_vec(), out))
    }

    #[test]
    fn samples_are_kept_in_tree_set_order() {
        assert_eq!(
            tree_set(&["b".into(), "B".into(), "a".into(), "b".into()]),
            ["B", "a", "b"]
        );
        let line = "c1\t1\t.\tA\tC\t.\t.\t.\tGT:GQ:DP\t0/1:30:5\t1/1:.:7\t0/0:99:3";
        let (samples, out) = edit(&["S1", "S2", "S3"], &[], line).expect("edits");
        assert!(samples.is_empty());
        assert_eq!(out, "c1\t1\t.\tA\tC\t.\t.\t.");
        let (samples, out) = edit(&["S1", "S2", "S3"], &["S3", "S1"], line).expect("edits");
        assert_eq!(samples, ["S1", "S3"]);
        assert_eq!(out, "c1\t1\t.\tA\tC\t.\t.\t.\tGT:DP:GQ\t0/1:5:30\t0/0:3:99");
        // A name the input lacks is a no-call at the record's ploidy.
        let (samples, out) = edit(
            &["S1", "S2"],
            &["S1", "NOPE"],
            line.rsplit_once('\t').expect("cut").0,
        )
        .expect("edits");
        assert_eq!(samples, ["NOPE", "S1"]);
        assert_eq!(out, "c1\t1\t.\tA\tC\t.\t.\t.\tGT:DP:GQ\t./.\t0/1:5:30");
        let (_, out) = edit(&["S1"], &["NOPE"], "c1\t1\t.\tA\tC\t.\t.\t.\tGT\t1").expect("edits");
        assert_eq!(out, "c1\t1\t.\tA\tC\t.\t.\t.\tGT\t./.");
        let (_, out) = edit(&[], &["X"], "c1\t1\t.\tA\tC\t.\t.\t.").expect("edits");
        assert_eq!(out, "c1\t1\t.\tA\tC\t.\t.\t.\tGT\t./.");
    }

    #[test]
    fn every_genotype_is_decoded_when_samples_are_kept() {
        let line = "c1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/5";
        assert!(edit(&["S1", "S2"], &["S1"], line).is_err());
        // With no sample kept the genotypes stay unread under sorted names.
        assert!(edit(&["S1", "S2"], &[], line).is_ok());
        assert!(edit(&["S2", "S1"], &[], line).is_err());
        let mut out = String::from("old");
        let subset = SampleSubset::new(&header(&["S1"]), &[]);
        assert!(!subset.edit("#comment", &mut out).expect("skips"));
        assert!(out.is_empty());
    }
}

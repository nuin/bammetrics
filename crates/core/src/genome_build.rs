//! The human genome build of a sequence dictionary, told from the `@SQ`
//! names and lengths alone (docs/doctor.md, "build"). No file and no
//! network is read: the primary chromosome lengths of GRCh37, GRCh38 and
//! T2T-CHM13 v2.0 are compiled in ([`PRIMARY`], from the UCSC
//! `chrom.sizes` of hg19, hg38 and hs1).
//!
//! A name is compared without a leading `chr`, and `MT` is `M`. Each
//! primary chromosome (1 to 22, X, Y) present names the assemblies whose
//! length it has. The dictionary is of an assembly when every primary
//! chromosome present has that assembly's length; a primary chromosome
//! of another assembly's length makes it mixed. The assembly with `chr`
//! names is reported by its UCSC name (`hg19`, `hg38`), without them by
//! its GRC name (`GRCh37`, `GRCh38`); T2T-CHM13 has one name. The other
//! contigs are counted by kind: alt (`_alt`), decoy (`_decoy`,
//! `hs37d5`), HLA (`HLA-`), and the rest (unplaced, random, EBV).

use std::fmt;

/// An assembly the table holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Assembly {
    Grch37,
    Grch38,
    T2t,
}

impl Assembly {
    pub const ALL: [Assembly; 3] = [Assembly::Grch37, Assembly::Grch38, Assembly::T2t];

    /// The name with or without `chr` names.
    pub fn name(self, chr: bool) -> &'static str {
        match (self, chr) {
            (Assembly::Grch37, false) => "GRCh37",
            (Assembly::Grch37, true) => "hg19",
            (Assembly::Grch38, false) => "GRCh38",
            (Assembly::Grch38, true) => "hg38",
            (Assembly::T2t, _) => "T2T-CHM13",
        }
    }

    fn column(self) -> usize {
        match self {
            Assembly::Grch37 => 0,
            Assembly::Grch38 => 1,
            Assembly::T2t => 2,
        }
    }
}

/// A build as a user names it: the assembly, and whether its contigs carry
/// `chr` (`None` for a name that says neither, as `T2T-CHM13`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Build {
    pub assembly: Assembly,
    pub chr: Option<bool>,
}

/// The names [`Build::parse`] takes, for an error message.
pub const BUILD_NAMES: &str = "GRCh37, b37, hg19, GRCh38, hg38, T2T-CHM13, chm13, hs1";

impl Build {
    /// A build name, any case: `GRCh37` and `b37` without `chr`, `hg19`
    /// with it; `GRCh38` without, `hg38` with; `T2T-CHM13`, `CHM13`,
    /// `T2T` and `hs1` either way.
    pub fn parse(text: &str) -> Option<Self> {
        let build = |assembly, chr| Some(Build { assembly, chr });
        match text.trim().to_ascii_lowercase().as_str() {
            "grch37" | "b37" | "hs37d5" => build(Assembly::Grch37, Some(false)),
            "hg19" => build(Assembly::Grch37, Some(true)),
            "grch38" => build(Assembly::Grch38, Some(false)),
            "hg38" => build(Assembly::Grch38, Some(true)),
            "t2t-chm13" | "t2t" | "chm13" | "hs1" | "chm13v2" | "t2t-chm13v2" => {
                build(Assembly::T2t, None)
            }
            _ => None,
        }
    }
}

impl fmt::Display for Build {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.assembly.name(self.chr.unwrap_or(false)))
    }
}

/// The primary chromosomes and their lengths in GRCh37, GRCh38 and
/// T2T-CHM13 v2.0, in that column order.
pub const PRIMARY: [(&str, [u64; 3]); 24] = [
    ("1", [249_250_621, 248_956_422, 248_387_328]),
    ("2", [243_199_373, 242_193_529, 242_696_752]),
    ("3", [198_022_430, 198_295_559, 201_105_948]),
    ("4", [191_154_276, 190_214_555, 193_574_945]),
    ("5", [180_915_260, 181_538_259, 182_045_439]),
    ("6", [171_115_067, 170_805_979, 172_126_628]),
    ("7", [159_138_663, 159_345_973, 160_567_428]),
    ("8", [146_364_022, 145_138_636, 146_259_331]),
    ("9", [141_213_431, 138_394_717, 150_617_247]),
    ("10", [135_534_747, 133_797_422, 134_758_134]),
    ("11", [135_006_516, 135_086_622, 135_127_769]),
    ("12", [133_851_895, 133_275_309, 133_324_548]),
    ("13", [115_169_878, 114_364_328, 113_566_686]),
    ("14", [107_349_540, 107_043_718, 101_161_492]),
    ("15", [102_531_392, 101_991_189, 99_753_195]),
    ("16", [90_354_753, 90_338_345, 96_330_374]),
    ("17", [81_195_210, 83_257_441, 84_276_897]),
    ("18", [78_077_248, 80_373_285, 80_542_538]),
    ("19", [59_128_983, 58_617_616, 61_707_364]),
    ("20", [63_025_520, 64_444_167, 66_210_255]),
    ("21", [48_129_895, 46_709_983, 45_090_682]),
    ("22", [51_304_566, 50_818_468, 51_324_926]),
    ("X", [155_270_560, 156_040_895, 154_259_566]),
    ("Y", [59_373_566, 57_227_415, 62_460_029]),
];

/// The mitochondrial length of the rCRS (GRCh37 `MT`, GRCh38, T2T) and of
/// the older UCSC hg19 `chrM`.
pub const MITO_RCRS: u64 = 16_569;
pub const MITO_HG19: u64 = 16_571;

/// The contigs of a dictionary by kind.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Contigs {
    /// 1 to 22, X and Y.
    pub primary: usize,
    pub mito: usize,
    pub alt: usize,
    pub decoy: usize,
    pub hla: usize,
    /// Unplaced, random, EBV and every other name.
    pub other: usize,
}

/// What [`detect`] found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Detection {
    /// Every primary chromosome present has the lengths of `assembly`.
    Known {
        assembly: Assembly,
        /// The primary chromosomes carry `chr`.
        chr: bool,
        /// Primary chromosomes present, of 24.
        matched: usize,
    },
    /// Primary chromosomes of two or more assemblies, or some of none: how
    /// many agree with the likeliest assembly, then each name that does
    /// not, with what its length says.
    Mixed(Vec<String>),
    /// No primary chromosome name, or none with a length the table knows.
    Unknown,
}

impl Detection {
    /// The build name: `hg38`, `GRCh37`, `T2T-CHM13`; `None` otherwise.
    pub fn name(&self) -> Option<&'static str> {
        match self {
            Detection::Known { assembly, chr, .. } => Some(assembly.name(*chr)),
            _ => None,
        }
    }
}

/// The primary chromosome `name` stands for, without `chr`; `M` for `MT`.
fn primary_key(name: &str) -> (&str, bool) {
    let (bare, chr) = match name.strip_prefix("chr") {
        Some(rest) => (rest, true),
        None => (name, false),
    };
    (if bare == "MT" { "M" } else { bare }, chr)
}

/// The mitochondrial contig names, in the order a tool looks for them:
/// every name [`chromosome`] reads as `M`.
pub const MITO_NAMES: [&str; 4] = ["chrM", "MT", "M", "chrMT"];

/// The human chromosome `name` stands for, without `chr`: `1` to `22`,
/// `X`, `Y`, or `M` for any of [`MITO_NAMES`]; `None` for every other name.
/// The one table of chromosome names: doctor, `CollectSexCheckMetrics`,
/// `CollectMtCopyNumberMetrics` and the annotation aliases of
/// `CollectExonCoverageMetrics` read it.
pub fn chromosome(name: &str) -> Option<&str> {
    let (key, _) = primary_key(name);
    (key == "M" || PRIMARY.iter().any(|(n, _)| *n == key)).then_some(key)
}

/// True for `chr1` to `chr22` and `1` to `22`.
pub fn is_autosome(name: &str) -> bool {
    chromosome(name).is_some_and(|key| !matches!(key, "X" | "Y" | "M"))
}

/// True for each of [`MITO_NAMES`].
pub fn is_mito(name: &str) -> bool {
    chromosome(name) == Some("M")
}

/// The other spellings of a contig name, in the order to try them: `chr`
/// added or removed, then for a mitochondrial name the rest of
/// [`MITO_NAMES`].
pub fn aliases(name: &str) -> Vec<String> {
    let mut names = match name.strip_prefix("chr") {
        Some(rest) if !rest.is_empty() => vec![rest.to_owned()],
        _ => vec![format!("chr{name}")],
    };
    if is_mito(name) {
        for other in MITO_NAMES {
            if other != name && !names.iter().any(|n| n == other) {
                names.push(other.to_owned());
            }
        }
    }
    names
}

/// The kind of contig `name` is.
fn kind(name: &str) -> Kind {
    let (key, _) = primary_key(name);
    if PRIMARY.iter().any(|(n, _)| *n == key) {
        Kind::Primary
    } else if key == "M" {
        Kind::Mito
    } else if name.ends_with("_alt") {
        Kind::Alt
    } else if name.ends_with("_decoy") || name == "hs37d5" {
        Kind::Decoy
    } else if name.starts_with("HLA-") {
        Kind::Hla
    } else {
        Kind::Other
    }
}

enum Kind {
    Primary,
    Mito,
    Alt,
    Decoy,
    Hla,
    Other,
}

/// The contigs of `dictionary` (name and length) by kind.
pub fn contigs<'a>(dictionary: impl IntoIterator<Item = (&'a str, u64)>) -> Contigs {
    let mut counts = Contigs::default();
    for (name, _) in dictionary {
        match kind(name) {
            Kind::Primary => counts.primary += 1,
            Kind::Mito => counts.mito += 1,
            Kind::Alt => counts.alt += 1,
            Kind::Decoy => counts.decoy += 1,
            Kind::Hla => counts.hla += 1,
            Kind::Other => counts.other += 1,
        }
    }
    counts
}

/// The assembly of `dictionary` (name and length per `@SQ`).
pub fn detect<'a>(dictionary: impl IntoIterator<Item = (&'a str, u64)>) -> Detection {
    let mut votes = [0usize; 3];
    let mut chr_names = 0usize;
    let mut seen: Vec<(&'a str, u64, Vec<Assembly>)> = Vec::new();
    for (name, length) in dictionary {
        let (key, chr) = primary_key(name);
        let Some((_, lengths)) = PRIMARY.iter().find(|(n, _)| *n == key) else {
            continue;
        };
        chr_names += usize::from(chr);
        let matching: Vec<Assembly> = Assembly::ALL
            .into_iter()
            .filter(|a| lengths[a.column()] == length)
            .collect();
        for assembly in &matching {
            votes[assembly.column()] += 1;
        }
        seen.push((name, length, matching));
    }
    let present = seen.len();
    if present == 0 || votes.iter().all(|&v| v == 0) {
        return Detection::Unknown;
    }
    let chr = chr_names * 2 > present;
    if let Some(assembly) = Assembly::ALL
        .into_iter()
        .find(|a| votes[a.column()] == present)
    {
        return Detection::Known {
            assembly,
            chr,
            matched: present,
        };
    }
    // The assembly most chromosomes agree with first, then each one that
    // does not, with what its length says.
    let best = Assembly::ALL
        .into_iter()
        .max_by_key(|a| (votes[a.column()], std::cmp::Reverse(a.column())))
        .unwrap_or(Assembly::Grch38);
    let mut said = vec![format!(
        "{} of {present} have {} lengths",
        votes[best.column()],
        best.name(false)
    )];
    for (name, length, matching) in seen {
        if matching.contains(&best) {
            continue;
        }
        let what = match matching.as_slice() {
            [] => "no build".to_owned(),
            several => several
                .iter()
                .map(|a| a.name(false))
                .collect::<Vec<_>>()
                .join(" or "),
        };
        said.push(format!("{name} has {length} ({what})"));
    }
    Detection::Mixed(said)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dictionary(assembly: Assembly, chr: bool) -> Vec<(String, u64)> {
        PRIMARY
            .iter()
            .map(|(name, lengths)| {
                let name = if chr {
                    format!("chr{name}")
                } else {
                    (*name).to_owned()
                };
                (name, lengths[assembly.column()])
            })
            .collect()
    }

    fn detected(dict: &[(String, u64)]) -> Detection {
        detect(dict.iter().map(|(n, l)| (n.as_str(), *l)))
    }

    #[test]
    fn each_assembly_is_found_with_and_without_chr() {
        for assembly in Assembly::ALL {
            for chr in [false, true] {
                let found = detected(&dictionary(assembly, chr));
                assert_eq!(
                    found,
                    Detection::Known {
                        assembly,
                        chr,
                        matched: 24
                    }
                );
            }
        }
        let hg19 = detected(&dictionary(Assembly::Grch37, true));
        assert_eq!(hg19.name(), Some("hg19"));
        assert_eq!(
            detected(&dictionary(Assembly::Grch38, false)).name(),
            Some("GRCh38")
        );
    }

    #[test]
    fn no_two_assemblies_share_a_primary_length() {
        for (name, lengths) in PRIMARY {
            assert_ne!(lengths[0], lengths[1], "{name}");
            assert_ne!(lengths[0], lengths[2], "{name}");
            assert_ne!(lengths[1], lengths[2], "{name}");
        }
    }

    #[test]
    fn a_subset_is_enough_and_other_contigs_are_ignored() {
        let mut dict = vec![("chr21".to_owned(), 46_709_983)];
        dict.push(("chr1_KI270706v1_random".to_owned(), 175_055));
        dict.push(("chrUn_JTFH01000001v1_decoy".to_owned(), 25_139));
        dict.push(("HLA-A*01:01:01:01".to_owned(), 3_503));
        dict.push(("chr6_GL000250v2_alt".to_owned(), 4_672_374));
        dict.push(("chrM".to_owned(), MITO_RCRS));
        assert_eq!(
            detected(&dict),
            Detection::Known {
                assembly: Assembly::Grch38,
                chr: true,
                matched: 1
            }
        );
        let counts = contigs(dict.iter().map(|(n, l)| (n.as_str(), *l)));
        assert_eq!(
            counts,
            Contigs {
                primary: 1,
                mito: 1,
                alt: 1,
                decoy: 1,
                hla: 1,
                other: 1
            }
        );
    }

    #[test]
    fn mixed_and_unknown_dictionaries_are_told_apart() {
        let mut dict = dictionary(Assembly::Grch38, true);
        dict[1].1 = PRIMARY[1].1[0];
        assert_eq!(
            detected(&dict),
            Detection::Mixed(vec![
                "23 of 24 have GRCh38 lengths".to_owned(),
                "chr2 has 243199373 (GRCh37)".to_owned()
            ])
        );
        assert_eq!(
            detected(&[("chr1".to_owned(), 100_000)]),
            Detection::Unknown
        );
        assert_eq!(detected(&[("contig7".to_owned(), 5)]), Detection::Unknown);
        assert_eq!(detected(&[]), Detection::Unknown);
    }

    #[test]
    fn chromosome_names_are_read_with_and_without_chr() {
        for name in MITO_NAMES {
            assert!(is_mito(name), "{name}");
            assert!(!is_autosome(name), "{name}");
        }
        for name in ["chr1", "22", "chr22"] {
            assert!(is_autosome(name), "{name}");
        }
        for name in [
            "chrX",
            "Y",
            "chr23",
            "01",
            "chr",
            "",
            "chrchr1",
            "MT1",
            "chrUn_KI270302v1",
        ] {
            assert!(!is_autosome(name), "{name}");
            assert!(!is_mito(name), "{name}");
        }
        assert_eq!(chromosome("chrX"), Some("X"));
        assert_eq!(chromosome("MT"), Some("M"));
        assert_eq!(chromosome("chr1_KI270706v1_random"), None);
        assert_eq!(aliases("chr7"), ["7"]);
        assert_eq!(aliases("7"), ["chr7"]);
        assert_eq!(aliases("chrM"), ["M", "MT", "chrMT"]);
        assert_eq!(aliases("MT"), ["chrMT", "chrM", "M"]);
        assert_eq!(aliases("chr"), ["chrchr"]);
    }

    #[test]
    fn build_names_parse_in_any_case() {
        assert_eq!(
            Build::parse("hg19"),
            Some(Build {
                assembly: Assembly::Grch37,
                chr: Some(true)
            })
        );
        assert_eq!(
            Build::parse("b37").map(|b| b.assembly),
            Some(Assembly::Grch37)
        );
        assert_eq!(Build::parse("GRCH38").map(|b| b.chr), Some(Some(false)));
        assert_eq!(Build::parse("hs1").map(|b| b.chr), Some(None));
        assert_eq!(Build::parse("hg18"), None);
    }
}

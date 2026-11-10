//! The `Coverage gaps` table (docs/06 section 1.7): the BED
//! `CollectCoverageGaps` writes as `GAPS_FILE`, read from beside the
//! `CoverageGapMetrics` file (`<metrics path>.bed`, the name `qc` gives
//! it), its longest gaps first.

/// The short class name of the gap metrics.
pub const CLASS: &str = "CoverageGapMetrics";

/// The most gaps the table shows.
pub const SHOWN: usize = 200;

/// One line of the gaps BED. Depths are kept as printed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gap {
    pub chrom: String,
    /// 0-based, half-open, as BED has it.
    pub start: i64,
    pub end: i64,
    pub target: String,
    pub mean_depth: String,
    pub min_depth: String,
}

/// The gaps of a BED text: `#` lines are skipped, a line with fewer than
/// six fields or a start or end that is not a number is skipped too.
pub fn parse(text: &str) -> Vec<Gap> {
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            let [chrom, start, end, target, mean_depth, min_depth, ..] = fields[..] else {
                return None;
            };
            Some(Gap {
                chrom: chrom.to_owned(),
                start: start.parse().ok()?,
                end: end.parse().ok()?,
                target: target.to_owned(),
                mean_depth: mean_depth.to_owned(),
                min_depth: min_depth.to_owned(),
            })
        })
        .collect()
}

/// The `n` longest gaps, longest first (file order between equal lengths),
/// and the total count.
pub fn longest(gaps: &[Gap], n: usize) -> (Vec<&Gap>, usize) {
    let mut sorted: Vec<&Gap> = gaps.iter().collect();
    sorted.sort_by_key(|g| std::cmp::Reverse(g.end - g.start));
    sorted.truncate(n);
    (sorted, gaps.len())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bed_lines_parse_and_sort_by_length() {
        let text = "#chrom\tstart\tend\ttarget\tmean_depth\tmin_depth\n\
                    chr1\t1000\t1200\tt1\t0.5\t0\n\
                    chr1\t2000\t2400\tt2\t3.25\t1\n\
                    bad line\n\
                    chr2\t10\t20\tt3\t0\t0\n";
        let gaps = parse(text);
        assert_eq!(gaps.len(), 3);
        assert_eq!(gaps[1].mean_depth, "3.25");
        let (shown, total) = longest(&gaps, 2);
        assert_eq!(total, 3);
        assert_eq!(
            shown.iter().map(|g| g.target.as_str()).collect::<Vec<_>>(),
            ["t2", "t1"]
        );
        assert_eq!(longest(&[], 5), (Vec::new(), 0));
    }
}

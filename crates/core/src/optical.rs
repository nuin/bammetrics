//! Read-name tile parser and optical duplicate finder (plan item M18),
//! mirroring Picard 3.5.0 `picard.sam.util.ReadNameParser` and
//! `picard.sam.markduplicates.util.OpticalDuplicateFinder`.
//!
//! The default `READ_NAME_REGEX` is not a regex: the name is split on `:`
//! and the last three fields are parsed as tile, x and y when the name has
//! 5 or 7 fields. Each field is read with `rapidParseInt`: an optional `-`,
//! then the leading digits, with 32-bit wrap; no digit at all means the
//! name has no location. The tile is cut to 16 bits (Java's `(short)`).
//! Any other value is a regex matched against the whole name with three
//! capture groups; `null` turns location parsing off. Picard compiles it as
//! a Java regex, bammetrics with the Rust `regex` crate: the common syntax
//! (classes, groups, named groups, `(?i)`, `\b`, `\p{..}`, lazy quantifiers)
//! means the same, while look-around and `\Q...\E` are refused with an
//! error.
//!
//! Two reads are an optical duplicate pair when both have a location, share
//! the read group index and the tile, and are within the pixel distance on
//! both axes. For a set of 2 reads, or 3 with a keeper, the pairs are
//! tested directly; larger sets are clustered through a graph, one read per
//! cluster kept (the keeper, or the smallest x then y).

use regex::Regex;

use crate::error::{Error, Result};

/// `ReadNameParser.DEFAULT_READ_NAME_REGEX`.
pub const DEFAULT_READ_NAME_REGEX: &str =
    "<optimized capture of last three ':' separated fields as numeric values>";
/// `OpticalDuplicateFinder.DEFAULT_OPTICAL_DUPLICATE_DISTANCE`.
pub const DEFAULT_OPTICAL_DUPLICATE_DISTANCE: i32 = 100;
/// `OpticalDuplicateFinder.DEFAULT_MAX_DUPLICATE_SET_SIZE`.
pub const DEFAULT_MAX_DUPLICATE_SET_SIZE: i64 = 300_000;

/// `PhysicalLocation`: tile, x and y from the read name, plus the read
/// group index optical duplicates must share. `tile` -1 means none.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Location {
    pub read_group: i16,
    pub tile: i16,
    pub x: i32,
    pub y: i32,
}

impl Default for Location {
    fn default() -> Self {
        Self {
            read_group: -1,
            tile: -1,
            x: -1,
            y: -1,
        }
    }
}

impl Location {
    pub fn has_location(&self) -> bool {
        self.tile != -1
    }
}

/// How read names are parsed.
#[derive(Clone, Debug)]
pub enum ReadNameParser {
    /// The 5 or 7 colon-separated field forms.
    Default,
    /// `READ_NAME_REGEX=null`: no location for any read.
    None,
    Regex(Regex),
}

impl ReadNameParser {
    /// The parser for a `READ_NAME_REGEX` value; `None` means the option
    /// was absent.
    pub fn new(read_name_regex: Option<&str>) -> Result<Self> {
        match read_name_regex {
            None | Some(DEFAULT_READ_NAME_REGEX) => Ok(Self::Default),
            Some("null") => Ok(Self::None),
            Some(pattern) => Regex::new(&format!("^(?:{pattern})$"))
                .map(Self::Regex)
                .map_err(|e| Error::Usage(format!("READ_NAME_REGEX {pattern:?}: {e}"))),
        }
    }

    /// `(tile, x, y)` of a read name, or `None` when the name does not
    /// carry a location.
    pub fn location(&self, name: &str) -> Option<(i16, i32, i32)> {
        match self {
            Self::None => None,
            Self::Default => {
                let (fields, tokens) = last_three_fields(name)?;
                (fields == 5 || fields == 7).then_some((tokens[0] as i16, tokens[1], tokens[2]))
            }
            Self::Regex(regex) => {
                let captures = regex.captures(name)?;
                let group = |i: usize| captures.get(i)?.as_str().parse::<i32>().ok();
                Some((group(1)? as i16, group(2)?, group(3)?))
            }
        }
    }
}

impl ReadNameParser {
    /// `readLocationInformation` into `loc`: true when the name gave tile,
    /// x and y. With a regex the fields are set one by one, so a name whose
    /// x or y is not a Java `int` leaves the tile (and x) set and returns
    /// false, as Picard does.
    pub fn read_location(&self, name: &str, loc: &mut Location) -> bool {
        match self {
            Self::None => false,
            Self::Default => match self.location(name) {
                Some((tile, x, y)) => {
                    (loc.tile, loc.x, loc.y) = (tile, x, y);
                    true
                }
                None => false,
            },
            Self::Regex(regex) => {
                let Some(captures) = regex.captures(name) else {
                    return false;
                };
                let group = |i: usize| captures.get(i)?.as_str().parse::<i32>().ok();
                let Some(tile) = group(1) else { return false };
                loc.tile = tile as i16;
                let Some(x) = group(2) else { return false };
                loc.x = x;
                let Some(y) = group(3) else { return false };
                loc.y = y;
                true
            }
        }
    }
}

/// `ReadNameParser.getLastThreeFields`: the field count and the last three
/// fields as numbers. `None` when a field is not a number
/// (`NumberFormatException`) or there are fewer than three.
fn last_three_fields(name: &str) -> Option<(usize, [i32; 3])> {
    let bytes = name.as_bytes();
    let mut tokens = [0i32; 3];
    let mut token = 2usize;
    let mut fields = 0usize;
    let mut end = bytes.len();
    let mut i = bytes.len() as isize - 1;
    while i >= 0 && token < 3 {
        let at = i as usize;
        if bytes[at] == b':' || at == 0 {
            fields += 1;
            let start = if at == 0 { 0 } else { at + 1 };
            tokens[token] = rapid_parse_int(&name[start..end])?;
            token = token.wrapping_sub(1);
            end = at;
        }
        i -= 1;
    }
    while i >= 0 {
        let at = i as usize;
        if bytes[at] == b':' || at == 0 {
            fields += 1;
        }
        i -= 1;
    }
    (fields >= 3).then_some((fields, tokens))
}

/// `ReadNameParser.rapidParseInt`: an optional `-` and the leading digits,
/// wrapping as a Java `int`; `None` without a digit.
fn rapid_parse_int(text: &str) -> Option<i32> {
    let bytes = text.as_bytes();
    let (negative, mut i) = match bytes.first() {
        Some(b'-') => (true, 1),
        _ => (false, 0),
    };
    let mut value = 0i32;
    let mut digits = false;
    while let Some(&c) = bytes.get(i) {
        if !c.is_ascii_digit() {
            break;
        }
        value = value.wrapping_mul(10).wrapping_add(i32::from(c - b'0'));
        digits = true;
        i += 1;
    }
    if !digits {
        return None;
    }
    Some(if negative {
        value.wrapping_neg()
    } else {
        value
    })
}

/// `OpticalDuplicateFinder`.
#[derive(Clone, Debug)]
pub struct OpticalDuplicateFinder {
    pub parser: ReadNameParser,
    pub pixel_distance: i32,
    /// `MAX_OPTICAL_DUPLICATE_SET_SIZE`: a larger set has no optical
    /// duplicates. Picard's constructor stores a negative value as given,
    /// so -1 turns detection off rather than the limit.
    pub max_set_size: i64,
}

impl OpticalDuplicateFinder {
    /// `closeEnough`: both located, same read group and tile, within the
    /// distance on both axes.
    fn close_enough(&self, lhs: &Location, rhs: &Location) -> bool {
        lhs.has_location()
            && rhs.has_location()
            && lhs.read_group == rhs.read_group
            && lhs.tile == rhs.tile
            && self.close_enough_short(lhs, rhs)
    }

    /// `closeEnoughShort`: the distance test alone (`Math.abs` of a Java
    /// `int`, so a wrapped difference keeps its sign).
    fn close_enough_short(&self, lhs: &Location, rhs: &Location) -> bool {
        lhs.x.wrapping_sub(rhs.x).wrapping_abs() <= self.pixel_distance
            && lhs.y.wrapping_sub(rhs.y).wrapping_abs() <= self.pixel_distance
    }

    /// `findOpticalDuplicates`: one flag per entry of `list`, true for a
    /// read that is an optical duplicate of another. `keeper` is the index
    /// of the read that must stay unflagged.
    pub fn find_optical_duplicates(&self, list: &[Location], keeper: Option<usize>) -> Vec<bool> {
        self.find_with_keeper(list, keeper, keeper.is_some())
    }

    /// `findOpticalDuplicates` with a keeper that is not in `list`
    /// (`MarkDuplicates` splits a set by strand and passes the set's keeper
    /// to both halves): no read is kept, and the graph starts at 4 reads,
    /// as with a keeper.
    pub fn find_optical_duplicates_keeper_elsewhere(&self, list: &[Location]) -> Vec<bool> {
        self.find_with_keeper(list, None, true)
    }

    fn find_with_keeper(
        &self,
        list: &[Location],
        keeper: Option<usize>,
        keeper_given: bool,
    ) -> Vec<bool> {
        let length = list.len();
        let mut flags = vec![false; length];
        if matches!(self.parser, ReadNameParser::None)
            || length < 2
            || length as i64 > self.max_set_size
        {
            return flags;
        }
        let actual_keeper = keeper.filter(|&k| list[k].has_location());
        if length >= if keeper_given { 4 } else { 3 } {
            self.flags_with_graph(list, actual_keeper, &mut flags);
        } else {
            self.flags_fast(list, actual_keeper, &mut flags);
        }
        flags
    }

    /// `getOpticalDuplicatesFlagFast`: every read close to the keeper is
    /// flagged; then each pair of the rest, the keeper left out, flags its
    /// second read, or its first when the second is flagged already.
    fn flags_fast(&self, list: &[Location], keeper: Option<usize>, flags: &mut [bool]) {
        if let Some(k) = keeper {
            for (i, other) in list.iter().enumerate() {
                flags[i] = i != k && self.close_enough(&list[k], other);
            }
        }
        for i in 0..list.len() {
            if Some(i) == keeper {
                continue;
            }
            for j in i + 1..list.len() {
                if Some(j) == keeper || (flags[i] && flags[j]) {
                    continue;
                }
                if self.close_enough(&list[i], &list[j]) {
                    let index = if flags[j] { i } else { j };
                    flags[index] = true;
                }
            }
        }
    }

    /// `getOpticalDuplicatesFlagWithGraph`: reads of one read group and tile
    /// are joined when close enough, clusters found by union-find, and in
    /// each cluster one read stays: the keeper's cluster keeps the keeper,
    /// any other keeps the read with the smallest x, then y.
    fn flags_with_graph(&self, list: &[Location], keeper: Option<usize>, flags: &mut [bool]) {
        let mut parent: Vec<usize> = (0..list.len()).collect();
        fn find(parent: &mut [usize], mut i: usize) -> usize {
            while parent[i] != i {
                parent[i] = parent[parent[i]];
                i = parent[i];
            }
            i
        }
        let mut by_tile: Vec<(i32, Vec<usize>)> = Vec::new();
        for (i, location) in list.iter().enumerate() {
            if !location.has_location() {
                continue;
            }
            let key = (i32::from(location.read_group) << 16).wrapping_add(i32::from(location.tile));
            match by_tile.iter_mut().find(|(k, _)| *k == key) {
                Some((_, group)) => group.push(i),
                None => by_tile.push((key, vec![i])),
            }
        }
        for (_, group) in &by_tile {
            for (a, &i) in group.iter().enumerate() {
                for &j in &group[a + 1..] {
                    if self.close_enough_short(&list[i], &list[j]) {
                        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                        if ri != rj {
                            parent[ri] = rj;
                        }
                    }
                }
            }
        }
        let clusters: Vec<usize> = (0..list.len()).map(|i| find(&mut parent, i)).collect();
        let mut representative: Vec<Option<usize>> = vec![None; list.len()];
        let keeper_cluster = keeper.map(|k| {
            representative[clusters[k]] = Some(k);
            clusters[k]
        });
        // `GraphUtils.cluster` returns a `HashMap<Integer, Integer>`, whose
        // keys come out by `i ^ (i >>> 16)`: the natural order below 65536.
        let mut order: Vec<usize> = (0..list.len()).collect();
        if list.len() > 1 << 16 {
            order.sort_unstable_by_key(|&i| i ^ (i >> 16));
        }
        for i in order {
            let cluster = clusters[i];
            match representative[cluster] {
                Some(rep) if Some(i) != keeper => {
                    let (current, best) = (&list[i], &list[rep]);
                    let smaller = current.x < best.x || (current.x == best.x && current.y < best.y);
                    if Some(cluster) != keeper_cluster && smaller {
                        flags[rep] = true;
                        representative[cluster] = Some(i);
                    } else {
                        flags[i] = true;
                    }
                }
                Some(_) => {}
                None => representative[cluster] = Some(i),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_names_take_the_last_three_of_five_or_seven_fields() {
        let parser = ReadNameParser::Default;
        assert_eq!(
            parser.location("A00123:45:HXXX:1:1101:1000:2000"),
            Some((1101, 1000, 2000))
        );
        assert_eq!(
            parser.location("M01234:12:000000000-ABCDE:1:1101"),
            Some((0, 1, 1101))
        );
        assert_eq!(
            parser.location("run:1:1101:1000:2000"),
            Some((1101, 1000, 2000))
        );
        assert_eq!(parser.location("1:2:3:4:5:6"), None);
        assert_eq!(parser.location("p005469"), None);
        assert_eq!(parser.location("a:b:c:d:e"), None);
        assert_eq!(parser.location("r:1:1101:10x:20"), Some((1101, 10, 20)));
        assert_eq!(parser.location("r:1:70000:1:-5"), Some((4464, 1, -5)));
        assert_eq!(
            parser.location("r:1:1:99999999999:1"),
            Some((1, 1_215_752_191, 1))
        );
        assert_eq!(ReadNameParser::None.location("r:1:2:3:4"), None);
    }

    #[test]
    fn regex_names_use_three_groups_over_the_whole_name() {
        let parser =
            ReadNameParser::new(Some(r"(?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*"))
                .expect("compiles");
        assert_eq!(parser.location("x;7a;300b;400c"), Some((7, 300, 400)));
        assert_eq!(parser.location("x;7a;300b;400c extra"), Some((7, 300, 400)));
        assert_eq!(parser.location("nothing"), None);
        let strict = ReadNameParser::new(Some(r"(\d+):(\d+):(\d+)")).expect("compiles");
        assert_eq!(strict.location("1:2:3"), Some((1, 2, 3)));
        assert_eq!(strict.location("1:2:3:4"), None);
        assert_eq!(strict.location("1:2:99999999999"), None);
        assert!(ReadNameParser::new(Some("(")).is_err());
    }

    fn finder() -> OpticalDuplicateFinder {
        OpticalDuplicateFinder {
            parser: ReadNameParser::Default,
            pixel_distance: 100,
            max_set_size: DEFAULT_MAX_DUPLICATE_SET_SIZE,
        }
    }

    fn at(tile: i16, x: i32, y: i32) -> Location {
        Location {
            read_group: 0,
            tile,
            x,
            y,
        }
    }

    #[test]
    fn small_sets_flag_pairs_directly() {
        let f = finder();
        assert_eq!(
            f.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50)], Some(1)),
            [true, false]
        );
        assert_eq!(
            f.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50)], Some(0)),
            [false, true]
        );
        assert_eq!(
            f.find_optical_duplicates(&[at(1, 10, 10), at(2, 50, 50)], Some(0)),
            [false, false]
        );
        // Three with a keeper at index 2: 0 and 1 are close to each other
        // but not to the keeper, so the second of the pair is flagged.
        assert_eq!(
            f.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50), at(1, 900, 900)], Some(2)),
            [false, true, false]
        );
        // A keeper without a location is no keeper.
        let none = Location::default();
        assert_eq!(
            f.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50), none], Some(2)),
            [false, true, false]
        );
        assert_eq!(f.find_optical_duplicates(&[at(1, 1, 1)], Some(0)), [false]);
        let off = OpticalDuplicateFinder {
            parser: ReadNameParser::None,
            ..finder()
        };
        assert_eq!(
            off.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50)], Some(1)),
            [false, false]
        );
        let tiny = OpticalDuplicateFinder {
            max_set_size: 1,
            ..finder()
        };
        assert_eq!(
            tiny.find_optical_duplicates(&[at(1, 10, 10), at(1, 50, 50)], Some(1)),
            [false, false]
        );
    }

    #[test]
    fn large_sets_keep_one_read_per_cluster() {
        let f = finder();
        // Chain 0-1-2 within 100 each, 3 far away, keeper 3 untouched.
        let list = [at(1, 0, 0), at(1, 90, 0), at(1, 180, 0), at(1, 5000, 0)];
        assert_eq!(
            f.find_optical_duplicates(&list, Some(3)),
            [false, true, true, false]
        );
        // The keeper's cluster keeps the keeper even when another read is
        // smaller; another cluster keeps its smallest x then y.
        let list = [
            at(1, 50, 50),
            at(1, 10, 10),
            at(1, 900, 5),
            at(1, 900, 1),
            at(1, 900, 3),
        ];
        assert_eq!(
            f.find_optical_duplicates(&list, Some(0)),
            [false, true, true, false, true]
        );
        // Different read groups never join.
        let mut other = at(1, 60, 60);
        other.read_group = 1;
        let list = [
            at(1, 50, 50),
            other,
            at(1, 55, 55),
            at(1, 900, 900),
            Location::default(),
        ];
        assert_eq!(
            f.find_optical_duplicates(&list, Some(3)),
            [false, false, true, false, false]
        );
    }
}

//! `ReadEndsForMarkDuplicates` (with the barcode fields of
//! `ReadEndsForMarkDuplicatesWithBarcodes`), its sort order and the record
//! a duplicate set's representative leaves for `TAG_DUPLICATE_SET_MEMBERS`.

use std::cmp::Ordering;
use std::io::{self, BufRead};

use crate::optical::Location;
use crate::spill::{Spill, get_u64, need_i64, put_i64};

/// `ReadEnds` orientation bytes.
pub const F: i8 = 0;
pub const R: i8 = 1;
pub const FF: i8 = 2;
pub const FR: i8 = 3;
pub const RR: i8 = 4;
pub const RF: i8 = 5;

/// `ReadEnds.getOrientationByte`.
pub fn orientation_byte(read1_negative: bool, read2_negative: bool) -> i8 {
    match (read1_negative, read2_negative) {
        (true, true) => RR,
        (true, false) => RF,
        (false, true) => FR,
        (false, false) => FF,
    }
}

/// One fragment end, or both ends of a pair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReadEnds {
    pub library_id: i16,
    /// `Objects.hash` of the top-strand UMI; 0 without barcode tags.
    pub barcode: i32,
    pub read_one_barcode: i32,
    pub read_two_barcode: i32,
    pub read1_reference: i32,
    pub read1_coordinate: i32,
    pub orientation: i8,
    pub read2_reference: i32,
    pub read2_coordinate: i32,
    pub tile: i16,
    pub x: i32,
    pub y: i32,
    pub read_group: i16,
    pub read1_index: i64,
    pub read2_index: i64,
    pub score: i16,
    pub orientation_for_optical: i8,
}

impl Default for ReadEnds {
    fn default() -> Self {
        Self {
            library_id: 0,
            barcode: 0,
            read_one_barcode: 0,
            read_two_barcode: 0,
            read1_reference: -1,
            read1_coordinate: -1,
            orientation: 0,
            read2_reference: -1,
            read2_coordinate: -1,
            tile: -1,
            x: -1,
            y: -1,
            read_group: -1,
            read1_index: -1,
            read2_index: -1,
            score: 0,
            orientation_for_optical: -1,
        }
    }
}

impl ReadEnds {
    /// `isPaired`: the mate is mapped.
    pub fn is_paired(&self) -> bool {
        self.read2_reference != -1
    }

    pub fn location(&self) -> Location {
        Location {
            read_group: self.read_group,
            tile: self.tile,
            x: self.x,
            y: self.y,
        }
    }

    /// `areComparableForDuplicates`: one duplicate set. The barcode fields
    /// are 0 when no barcode tag is set, so comparing them always is the
    /// same as Picard's `useBarcodes` test.
    pub fn same_set(&self, other: &Self, compare_read2: bool) -> bool {
        self.library_id == other.library_id
            && self.barcode == other.barcode
            && self.read_one_barcode == other.read_one_barcode
            && self.read_two_barcode == other.read_two_barcode
            && self.read1_reference == other.read1_reference
            && self.read1_coordinate == other.read1_coordinate
            && self.orientation == other.orientation
            && (!compare_read2
                || (self.read2_reference == other.read2_reference
                    && self.read2_coordinate == other.read2_coordinate))
    }
}

/// The sign of a Java `int` subtraction, which wraps.
fn java_difference(lhs: i32, rhs: i32) -> Ordering {
    lhs.wrapping_sub(rhs).cmp(&0)
}

impl Ord for ReadEnds {
    /// `ReadEndsMDComparator`.
    fn cmp(&self, other: &Self) -> Ordering {
        i32::from(self.library_id)
            .cmp(&i32::from(other.library_id))
            .then(self.barcode.cmp(&other.barcode))
            .then(self.read_one_barcode.cmp(&other.read_one_barcode))
            .then(self.read_two_barcode.cmp(&other.read_two_barcode))
            .then_with(|| java_difference(self.read1_reference, other.read1_reference))
            .then_with(|| java_difference(self.read1_coordinate, other.read1_coordinate))
            .then(self.orientation.cmp(&other.orientation))
            .then_with(|| java_difference(self.read2_reference, other.read2_reference))
            .then_with(|| java_difference(self.read2_coordinate, other.read2_coordinate))
            .then(self.tile.cmp(&other.tile))
            .then_with(|| java_difference(self.x, other.x))
            .then_with(|| java_difference(self.y, other.y))
            .then(self.read1_index.cmp(&other.read1_index))
            .then(self.read2_index.cmp(&other.read2_index))
    }
}

impl PartialOrd for ReadEnds {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Spill for ReadEnds {
    fn encode(&self, out: &mut Vec<u8>) {
        for v in [
            i64::from(self.library_id),
            i64::from(self.barcode),
            i64::from(self.read_one_barcode),
            i64::from(self.read_two_barcode),
            i64::from(self.read1_reference),
            i64::from(self.read1_coordinate),
            i64::from(self.orientation),
            i64::from(self.read2_reference),
            i64::from(self.read2_coordinate),
            i64::from(self.tile),
            i64::from(self.x),
            i64::from(self.y),
            i64::from(self.read_group),
            self.read1_index,
            self.read2_index,
            i64::from(self.score),
            i64::from(self.orientation_for_optical),
        ] {
            put_i64(out, v);
        }
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        let Some(first) = get_u64(input)? else {
            return Ok(None);
        };
        let library_id = ((first >> 1) as i64 ^ -((first & 1) as i64)) as i16;
        let mut next = || need_i64(input);
        Ok(Some(Self {
            library_id,
            barcode: next()? as i32,
            read_one_barcode: next()? as i32,
            read_two_barcode: next()? as i32,
            read1_reference: next()? as i32,
            read1_coordinate: next()? as i32,
            orientation: next()? as i8,
            read2_reference: next()? as i32,
            read2_coordinate: next()? as i32,
            tile: next()? as i16,
            x: next()? as i32,
            y: next()? as i32,
            read_group: next()? as i16,
            read1_index: next()?,
            read2_index: next()?,
            score: next()? as i16,
            orientation_for_optical: next()? as i8,
        }))
    }
}

/// `RepresentativeReadIndexer`: a read of a duplicate set, the set's
/// representative and its size, as Java `int`s. Sorted by the read only
/// (`Comparator.comparing(read -> read.readIndexInFile)`).
#[derive(Clone, Copy, Debug)]
pub struct Representative {
    pub read_index: i32,
    pub representative_index: i32,
    pub set_size: i32,
}

impl PartialEq for Representative {
    fn eq(&self, other: &Self) -> bool {
        self.read_index == other.read_index
    }
}

impl Eq for Representative {}

impl Ord for Representative {
    fn cmp(&self, other: &Self) -> Ordering {
        self.read_index.cmp(&other.read_index)
    }
}

impl PartialOrd for Representative {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Spill for Representative {
    fn encode(&self, out: &mut Vec<u8>) {
        put_i64(out, i64::from(self.read_index));
        put_i64(out, i64::from(self.representative_index));
        put_i64(out, i64::from(self.set_size));
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        let Some(first) = get_u64(input)? else {
            return Ok(None);
        };
        let read_index = ((first >> 1) as i64 ^ -((first & 1) as i64)) as i32;
        Ok(Some(Self {
            read_index,
            representative_index: need_i64(input)? as i32,
            set_size: need_i64(input)? as i32,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ends(library_id: i16, coordinate: i32, index: i64) -> ReadEnds {
        ReadEnds {
            library_id,
            read1_reference: 0,
            read1_coordinate: coordinate,
            read1_index: index,
            ..ReadEnds::default()
        }
    }

    #[test]
    fn library_comes_before_position_and_file_order_breaks_ties() {
        let mut list = [ends(2, 5, 0), ends(1, 9, 1), ends(1, 5, 3), ends(1, 5, 2)];
        list.sort();
        let order: Vec<i64> = list.iter().map(|e| e.read1_index).collect();
        assert_eq!(order, [2, 3, 1, 0]);
    }

    #[test]
    fn spilled_ends_come_back_unchanged() {
        let original = ReadEnds {
            library_id: 3,
            barcode: -7,
            read_one_barcode: i32::MIN,
            read_two_barcode: i32::MAX,
            read1_reference: 4,
            read1_coordinate: -12,
            orientation: FR,
            read2_reference: 5,
            read2_coordinate: 1_000_000,
            tile: 1101,
            x: 20_000,
            y: -1,
            read_group: 2,
            read1_index: 77,
            read2_index: 1 << 40,
            score: -32_000,
            orientation_for_optical: RF,
        };
        let mut bytes = Vec::new();
        original.encode(&mut bytes);
        original.encode(&mut bytes);
        let mut input = bytes.as_slice();
        assert_eq!(
            ReadEnds::decode(&mut input).expect("one"),
            Some(original.clone())
        );
        assert_eq!(ReadEnds::decode(&mut input).expect("two"), Some(original));
        assert_eq!(ReadEnds::decode(&mut input).expect("end"), None);
        assert_eq!(orientation_byte(true, false), RF);
    }
}

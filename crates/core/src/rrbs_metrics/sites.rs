//! The CpG sites of a contig a unit has left, packed until the detail file
//! is written: sorted by position, each site as three LEB128 numbers (the
//! distance from the site before, the reads, the converted reads). A site
//! of an exome takes 3 to 4 bytes this way, 12 as three `u32`.

/// A CpG site: the 0-based position of its C, its reads and its converted
/// reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Site {
    pub position: u32,
    pub total: u32,
    pub converted: u32,
}

/// The sites of one contig, packed.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Sites {
    bytes: Vec<u8>,
    len: usize,
}

impl Sites {
    /// `sites` sorted by position and packed. Two sites never share a
    /// position: a unit counts each position once.
    pub fn pack(sites: &mut [Site]) -> Self {
        sites.sort_unstable_by_key(|site| site.position);
        let mut bytes = Vec::with_capacity(sites.len() * 4);
        let mut previous = 0;
        for site in sites.iter() {
            put(&mut bytes, site.position - previous);
            put(&mut bytes, site.total);
            put(&mut bytes, site.converted);
            previous = site.position;
        }
        bytes.shrink_to_fit();
        Self {
            bytes,
            len: sites.len(),
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    /// The sites in position order.
    pub fn iter(&self) -> impl Iterator<Item = Site> + '_ {
        let mut bytes = self.bytes.as_slice();
        let mut position = 0;
        (0..self.len).map(move |_| {
            position += take(&mut bytes);
            Site {
                position,
                total: take(&mut bytes),
                converted: take(&mut bytes),
            }
        })
    }
}

/// `value` as LEB128: 7 bits a byte, the high bit set on every byte but
/// the last.
fn put(bytes: &mut Vec<u8>, mut value: u32) {
    while value >= 0x80 {
        bytes.push((value as u8) | 0x80);
        value >>= 7;
    }
    bytes.push(value as u8);
}

/// The next LEB128 number of `bytes`, which [`put`] wrote.
fn take(bytes: &mut &[u8]) -> u32 {
    let mut value = 0;
    let mut shift = 0;
    while let Some((&byte, rest)) = bytes.split_first() {
        *bytes = rest;
        value |= u32::from(byte & 0x7f).checked_shl(shift).unwrap_or(0);
        if byte < 0x80 {
            break;
        }
        shift += 7;
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sites_come_back_sorted_and_whole() {
        let site = |position, total, converted| Site {
            position,
            total,
            converted,
        };
        let mut sites = vec![
            site(u32::MAX, u32::MAX, 0),
            site(5, 1, 1),
            site(0, 128, 127),
            site(300_000, 16_384, 2),
        ];
        let packed = Sites::pack(&mut sites);
        assert_eq!(packed.len(), 4);
        assert_eq!(packed.iter().collect::<Vec<_>>(), sites);
        assert_eq!(sites[0].position, 0);
        assert!(Sites::pack(&mut []).iter().next().is_none());
    }
}

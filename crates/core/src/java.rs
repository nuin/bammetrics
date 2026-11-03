//! Java behaviour some Picard output depends on (plan item M22):
//! `String.hashCode` and the iteration order of a `HashMap<String, _>`.
//!
//! `EstimateLibraryComplexity` keeps its per-library histograms in a
//! `HashMap` and writes rows in its iteration order. That order is the
//! bucket index of each key's spread hash in the table, ties in insertion
//! order (Java 8 keeps a bucket's list order through a resize). The table
//! starts at 16 slots and doubles whenever the entry count passes three
//! quarters of it. A bucket turns into a tree only past 8 entries in one
//! bucket with 64 slots or more, which no library list reaches.

/// `String.hashCode`: `31 * h + c` over UTF-16 code units, wrapping.
pub fn string_hash_code(text: &str) -> i32 {
    text.encode_utf16()
        .fold(0i32, |h, c| h.wrapping_mul(31).wrapping_add(i32::from(c)))
}

/// `HashMap.hash`: the hash code with its upper half folded in.
fn spread(hash: i32) -> i32 {
    hash ^ ((hash as u32) >> 16) as i32
}

/// Slots of a `HashMap` holding `size` entries, grown from 16 by doubling
/// while `size` exceeds the load threshold.
fn capacity(size: usize) -> usize {
    let mut capacity = 16usize;
    while size > capacity * 3 / 4 {
        capacity *= 2;
    }
    capacity
}

/// Indices of `keys` (given in insertion order, no repeats) in the order a
/// `HashMap<String, _>` iterates them.
pub fn hash_map_order(keys: &[&str]) -> Vec<usize> {
    let hashes: Vec<i32> = keys.iter().map(|k| string_hash_code(k)).collect();
    hash_map_order_by_hash(&hashes)
}

/// [`hash_map_order`] for keys of any class: `hashes` are their `hashCode`
/// values in insertion order. A `HashSet` iterates like the `HashMap` it
/// wraps.
pub fn hash_map_order_by_hash(hashes: &[i32]) -> Vec<usize> {
    let mask = (capacity(hashes.len()) - 1) as i32;
    let mut order: Vec<usize> = (0..hashes.len()).collect();
    order.sort_by_key(|&i| spread(hashes[i]) & mask);
    order
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hash_code_matches_java() {
        assert_eq!(string_hash_code(""), 0);
        assert_eq!(string_hash_code("Unknown"), 1_379_812_394);
        assert_eq!(string_hash_code("lib1"), 3_321_420);
        assert_eq!(string_hash_code("LibA"), 2_368_124);
        assert_eq!(string_hash_code("polygenelubricants"), i32::MIN);
    }

    /// Each expectation is the `keySet()` order a JVM (OpenJDK 27) printed
    /// for a `HashMap` filled in the given order.
    #[test]
    fn iteration_order_is_by_bucket_then_insertion() {
        let order = hash_map_order(&["L2", "L1", "Unknown", "L3"]);
        assert_eq!(order, [2, 1, 0, 3]);
        let order = hash_map_order(&["LibD", "LibA", "Unknown", "LibC", "LibB"]);
        assert_eq!(order, [2, 1, 4, 3, 0]);
        // "Aa" and "BB" share a hash code: insertion order within the bucket.
        assert_eq!(hash_map_order(&["Aa", "BB"]), [0, 1]);
        assert_eq!(hash_map_order(&["BB", "Aa"]), [0, 1]);
        // 14 keys: the table has grown to 32 slots.
        let keys: Vec<String> = (1..=14).map(|i| format!("a{i}")).collect();
        let keys: Vec<&str> = keys.iter().map(String::as_str).collect();
        let order: Vec<&str> = hash_map_order(&keys).into_iter().map(|i| keys[i]).collect();
        assert_eq!(
            order,
            [
                "a11", "a10", "a13", "a12", "a14", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8",
                "a9"
            ]
        );
        assert_eq!(capacity(12), 16);
        assert_eq!(capacity(13), 32);
        assert_eq!(capacity(25), 64);
    }
}

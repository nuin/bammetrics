//! Orders a JVM (OpenJDK 27) printed: `new ArrayList<>(set)` of a
//! `HashSet<String>` filled in the order given, as insertion numbers. The
//! long orders are kept as a checksum ([`checksum`]).

use super::*;

/// Strings of `k` pairs of `Aa` and `BB`, which all share one hash code,
/// in binary order (`Aa` for 0).
fn combos(k: u32) -> Vec<String> {
    (0..1u32 << k)
        .map(|m| {
            (0..k)
                .rev()
                .map(|i| if (m >> i) & 1 == 0 { "Aa" } else { "BB" })
                .collect()
        })
        .collect()
}

fn order(keys: &[String]) -> Vec<usize> {
    let mut set = HashSetOrder::new();
    for key in keys {
        set.insert(key.as_bytes());
    }
    set.order()
}

fn checksum(order: &[usize]) -> u64 {
    order.iter().fold(0u64, |s, &x| {
        (s.wrapping_mul(1_000_003).wrapping_add(x as u64)) % (1 << 61)
    })
}

#[test]
fn one_hash_code_makes_a_tree() {
    assert_eq!(
        order(&combos(5)),
        [
            7, 3, 0, 1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31
        ]
    );
}

#[test]
fn a_tree_among_lists() {
    let six = combos(6);
    let mut keys = Vec::new();
    for i in 0..2000 {
        keys.push(format!("x{i}"));
        if i % 30 == 0 && i / 30 < six.len() {
            keys.push(six[i / 30].clone());
        }
    }
    let order = order(&keys);
    assert_eq!(order.len(), 2064);
    assert_eq!(
        &order[..12],
        [
            1974, 1966, 1967, 1964, 1965, 1970, 1971, 1968, 1969, 1972, 1973, 1995
        ]
    );
    assert_eq!(checksum(&order), 775836436710489108);
}

#[test]
fn illumina_names() {
    let keys: Vec<String> = (0..20000)
        .map(|i| {
            format!(
                "M1:7:FC1:1:{}:{}:{}",
                1101 + i % 7,
                1000 + (i * 37) % 30000,
                2000 + i / 7
            )
        })
        .collect();
    let order = order(&keys);
    assert_eq!(
        &order[..12],
        [
            11063, 11855, 4108, 10557, 3497, 9218, 5334, 4904, 18221, 14223, 5557, 5794
        ]
    );
    assert_eq!(checksum(&order), 619249332810128294);
}

/// Two groups of 16 whose spread hashes share the low 9 bits: one tree of
/// 32 until the table reaches 1,024 slots, then two trees made again.
#[test]
fn a_tree_splits_into_two_trees() {
    let four = combos(4);
    let mut keys = Vec::new();
    let mut group = 0;
    for i in 0..700 {
        keys.push(format!("f{i}"));
        if i % 20 == 0 && group < 16 {
            keys.push(format!("p0{}", four[group]));
            keys.push(format!("p215{}", four[group]));
            group += 1;
        }
    }
    let order = order(&keys);
    assert_eq!(order.len(), 732);
    assert_eq!(checksum(&order), 1900439107876390308);
}

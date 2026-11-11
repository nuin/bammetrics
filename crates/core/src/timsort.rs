//! `java.util.TimSort` (JDK 11 to 27, unchanged since the `mergeCollapse`
//! fix of JDK 11): the sort behind `Arrays.sort(T[], Comparator)`, which
//! htsjdk's `SortingCollection` calls on every run it holds in memory.
//!
//! The port gives the order Java gives for any comparator. For a total
//! order that is the order of any stable sort. For a comparator that is not
//! one (`SAMRecordQueryNameComparator` on a record flagged both first and
//! second of its pair) it is still Java's order, where Rust's own sort may
//! panic. Java throws `IllegalArgumentException("Comparison method violates
//! its general contract!")` where TimSort notices; [`sort_by`] returns
//! [`ContractViolation`] there and never panics.
//!
//! The sort moves indices, not items: the comparator reads `items[i]`, and
//! the items are put in order once at the end.

use std::cmp::Ordering;

/// Runs shorter than this are sorted by binary insertion (`MIN_MERGE`).
const MIN_MERGE: usize = 32;
/// `MIN_GALLOP`.
const MIN_GALLOP: usize = 7;

/// The message Java throws when TimSort sees an inconsistent comparator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ContractViolation;

impl std::fmt::Display for ContractViolation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Comparison method violates its general contract!")
    }
}

/// Sort `items` as `Arrays.sort(items, compare)` sorts them.
pub fn sort_by<T>(
    items: &mut Vec<T>,
    mut compare: impl FnMut(&T, &T) -> Ordering,
) -> Result<(), ContractViolation> {
    if items.len() < 2 {
        return Ok(());
    }
    let order = {
        let view: &[T] = items;
        let mut c = |i: usize, j: usize| compare(&view[i], &view[j]);
        sort_indices(view.len(), &mut c)?
    };
    let mut slots: Vec<Option<T>> = std::mem::take(items).into_iter().map(Some).collect();
    items.extend(
        order
            .into_iter()
            .map(|i| slots[i].take().expect("TimSort moves a permutation")),
    );
    Ok(())
}

/// `0..n` sorted by `c`, which compares two indices.
fn sort_indices(
    n: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> Result<Vec<usize>, ContractViolation> {
    let mut a: Vec<usize> = (0..n).collect();
    let mut lo = 0;
    let hi = n;
    let mut remaining = n;
    if remaining < MIN_MERGE {
        let run = count_run_and_make_ascending(&mut a, lo, hi, c);
        binary_sort(&mut a, lo, hi, lo + run, c);
        return Ok(a);
    }
    let mut ts = TimSort {
        a,
        tmp: Vec::new(),
        runs: Vec::new(),
        min_gallop: MIN_GALLOP,
    };
    let min_run = min_run_length(remaining);
    loop {
        let mut run = count_run_and_make_ascending(&mut ts.a, lo, hi, c);
        if run < min_run {
            let force = remaining.min(min_run);
            binary_sort(&mut ts.a, lo, lo + force, lo + run, c);
            run = force;
        }
        ts.runs.push((lo, run));
        ts.merge_collapse(c)?;
        lo += run;
        remaining -= run;
        if remaining == 0 {
            break;
        }
    }
    ts.merge_force_collapse(c)?;
    Ok(ts.a)
}

/// `binarySort`: `a[lo..start]` is sorted; insert the rest.
fn binary_sort(
    a: &mut [usize],
    lo: usize,
    hi: usize,
    mut start: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) {
    if start == lo {
        start += 1;
    }
    while start < hi {
        let pivot = a[start];
        let (mut left, mut right) = (lo, start);
        while left < right {
            let mid = (left + right) >> 1;
            if c(pivot, a[mid]) == Ordering::Less {
                right = mid;
            } else {
                left = mid + 1;
            }
        }
        a.copy_within(left..start, left + 1);
        a[left] = pivot;
        start += 1;
    }
}

/// `countRunAndMakeAscending`.
fn count_run_and_make_ascending(
    a: &mut [usize],
    lo: usize,
    hi: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> usize {
    let mut run_hi = lo + 1;
    if run_hi == hi {
        return 1;
    }
    let descending = c(a[run_hi], a[lo]) == Ordering::Less;
    run_hi += 1;
    if descending {
        while run_hi < hi && c(a[run_hi], a[run_hi - 1]) == Ordering::Less {
            run_hi += 1;
        }
        a[lo..run_hi].reverse();
    } else {
        while run_hi < hi && c(a[run_hi], a[run_hi - 1]) != Ordering::Less {
            run_hi += 1;
        }
    }
    run_hi - lo
}

/// `minRunLength`.
fn min_run_length(mut n: usize) -> usize {
    let mut r = 0;
    while n >= MIN_MERGE {
        r |= n & 1;
        n >>= 1;
    }
    n + r
}

/// `gallopLeft`: the leftmost place for `key` in `run[base..base + len]`.
fn gallop_left(
    key: usize,
    run: &[usize],
    base: usize,
    len: usize,
    hint: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> usize {
    let mut last_ofs = 0usize;
    let mut ofs = 1usize;
    if c(key, run[base + hint]) == Ordering::Greater {
        let max_ofs = len - hint;
        while ofs < max_ofs && c(key, run[base + hint + ofs]) == Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        last_ofs += hint;
        ofs += hint;
    } else {
        let max_ofs = hint + 1;
        while ofs < max_ofs && c(key, run[base + hint - ofs]) != Ordering::Greater {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        // Java keeps `lastOfs` as `hint - ofs`, which may be -1; it is
        // incremented before use, so `left` holds the incremented value.
        let tmp = last_ofs;
        let left = hint + 1 - ofs;
        ofs = hint - tmp;
        return search_left(key, run, base, left, ofs, c);
    }
    search_left(key, run, base, last_ofs + 1, ofs, c)
}

/// The binary search that ends `gallopLeft`, over `[left, right)`.
fn search_left(
    key: usize,
    run: &[usize],
    base: usize,
    mut left: usize,
    mut right: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> usize {
    while left < right {
        let m = left + ((right - left) >> 1);
        if c(key, run[base + m]) == Ordering::Greater {
            left = m + 1;
        } else {
            right = m;
        }
    }
    right
}

/// `gallopRight`: the rightmost place for `key` in `run[base..base + len]`.
fn gallop_right(
    key: usize,
    run: &[usize],
    base: usize,
    len: usize,
    hint: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> usize {
    let mut ofs = 1usize;
    let mut last_ofs = 0usize;
    if c(key, run[base + hint]) == Ordering::Less {
        let max_ofs = hint + 1;
        while ofs < max_ofs && c(key, run[base + hint - ofs]) == Ordering::Less {
            last_ofs = ofs;
            ofs = (ofs << 1) + 1;
        }
        ofs = ofs.min(max_ofs);
        let tmp = last_ofs;
        let left = hint + 1 - ofs;
        ofs = hint - tmp;
        return search_right(key, run, base, left, ofs, c);
    }
    let max_ofs = len - hint;
    while ofs < max_ofs && c(key, run[base + hint + ofs]) != Ordering::Less {
        last_ofs = ofs;
        ofs = (ofs << 1) + 1;
    }
    ofs = ofs.min(max_ofs);
    last_ofs += hint;
    ofs += hint;
    search_right(key, run, base, last_ofs + 1, ofs, c)
}

/// The binary search that ends `gallopRight`, over `[left, right)`.
fn search_right(
    key: usize,
    run: &[usize],
    base: usize,
    mut left: usize,
    mut right: usize,
    c: &mut dyn FnMut(usize, usize) -> Ordering,
) -> usize {
    while left < right {
        let m = left + ((right - left) >> 1);
        if c(key, run[base + m]) == Ordering::Less {
            right = m;
        } else {
            left = m + 1;
        }
    }
    right
}

/// The state of one sort: the indices, the merge buffer and the run stack.
struct TimSort {
    a: Vec<usize>,
    tmp: Vec<usize>,
    /// `(runBase, runLen)` per pending run.
    runs: Vec<(usize, usize)>,
    min_gallop: usize,
}

impl TimSort {
    fn len(&self, i: usize) -> usize {
        self.runs[i].1
    }

    /// `mergeCollapse`, as fixed in JDK 11.
    fn merge_collapse(
        &mut self,
        c: &mut dyn FnMut(usize, usize) -> Ordering,
    ) -> Result<(), ContractViolation> {
        while self.runs.len() > 1 {
            let mut n = self.runs.len() - 2;
            if (n > 0 && self.len(n - 1) <= self.len(n) + self.len(n + 1))
                || (n > 1 && self.len(n - 2) <= self.len(n) + self.len(n - 1))
            {
                if self.len(n - 1) < self.len(n + 1) {
                    n -= 1;
                }
            } else if self.len(n) > self.len(n + 1) {
                break;
            }
            self.merge_at(n, c)?;
        }
        Ok(())
    }

    /// `mergeForceCollapse`.
    fn merge_force_collapse(
        &mut self,
        c: &mut dyn FnMut(usize, usize) -> Ordering,
    ) -> Result<(), ContractViolation> {
        while self.runs.len() > 1 {
            let mut n = self.runs.len() - 2;
            if n > 0 && self.len(n - 1) < self.len(n + 1) {
                n -= 1;
            }
            self.merge_at(n, c)?;
        }
        Ok(())
    }

    /// `mergeAt`: merge runs `i` and `i + 1`.
    fn merge_at(
        &mut self,
        i: usize,
        c: &mut dyn FnMut(usize, usize) -> Ordering,
    ) -> Result<(), ContractViolation> {
        let (mut base1, mut len1) = self.runs[i];
        let (base2, mut len2) = self.runs[i + 1];
        self.runs[i].1 = len1 + len2;
        self.runs.remove(i + 1);
        let k = gallop_right(self.a[base2], &self.a, base1, len1, 0, c);
        base1 += k;
        len1 -= k;
        if len1 == 0 {
            return Ok(());
        }
        len2 = gallop_left(self.a[base1 + len1 - 1], &self.a, base2, len2, len2 - 1, c);
        if len2 == 0 {
            return Ok(());
        }
        if len1 <= len2 {
            self.merge_lo(base1, len1, base2, len2, c)
        } else {
            self.merge_hi(base1, len1, base2, len2, c)
        }
    }

    /// `mergeLo`: the first run is copied out and merged from the left.
    fn merge_lo(
        &mut self,
        base1: usize,
        mut len1: usize,
        base2: usize,
        mut len2: usize,
        c: &mut dyn FnMut(usize, usize) -> Ordering,
    ) -> Result<(), ContractViolation> {
        self.tmp.clear();
        self.tmp.extend_from_slice(&self.a[base1..base1 + len1]);
        let (a, tmp) = (&mut self.a, &self.tmp);
        let mut cursor1 = 0;
        let mut cursor2 = base2;
        let mut dest = base1;
        a[dest] = a[cursor2];
        dest += 1;
        cursor2 += 1;
        len2 -= 1;
        if len2 == 0 {
            a[dest..dest + len1].copy_from_slice(&tmp[cursor1..cursor1 + len1]);
            return Ok(());
        }
        if len1 == 1 {
            a.copy_within(cursor2..cursor2 + len2, dest);
            a[dest + len2] = tmp[cursor1];
            return Ok(());
        }
        let mut min_gallop = self.min_gallop;
        'outer: loop {
            let mut count1 = 0usize;
            let mut count2 = 0usize;
            loop {
                if c(a[cursor2], tmp[cursor1]) == Ordering::Less {
                    a[dest] = a[cursor2];
                    dest += 1;
                    cursor2 += 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 0 {
                        break 'outer;
                    }
                } else {
                    a[dest] = tmp[cursor1];
                    dest += 1;
                    cursor1 += 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 1 {
                        break 'outer;
                    }
                }
                if (count1 | count2) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 = gallop_right(a[cursor2], tmp, cursor1, len1, 0, c);
                if count1 != 0 {
                    a[dest..dest + count1].copy_from_slice(&tmp[cursor1..cursor1 + count1]);
                    dest += count1;
                    cursor1 += count1;
                    len1 -= count1;
                    if len1 <= 1 {
                        break 'outer;
                    }
                }
                a[dest] = a[cursor2];
                dest += 1;
                cursor2 += 1;
                len2 -= 1;
                if len2 == 0 {
                    break 'outer;
                }
                count2 = gallop_left(tmp[cursor1], a, cursor2, len2, 0, c);
                if count2 != 0 {
                    a.copy_within(cursor2..cursor2 + count2, dest);
                    dest += count2;
                    cursor2 += count2;
                    len2 -= count2;
                    if len2 == 0 {
                        break 'outer;
                    }
                }
                a[dest] = tmp[cursor1];
                dest += 1;
                cursor1 += 1;
                len1 -= 1;
                if len1 == 1 {
                    break 'outer;
                }
                min_gallop = min_gallop.saturating_sub(1);
                if count1 < MIN_GALLOP && count2 < MIN_GALLOP {
                    break;
                }
            }
            min_gallop += 2;
        }
        self.min_gallop = min_gallop.max(1);
        match len1 {
            1 => {
                a.copy_within(cursor2..cursor2 + len2, dest);
                a[dest + len2] = tmp[cursor1];
                Ok(())
            }
            0 => Err(ContractViolation),
            _ => {
                a[dest..dest + len1].copy_from_slice(&tmp[cursor1..cursor1 + len1]);
                Ok(())
            }
        }
    }

    /// `mergeHi`: the second run is copied out and merged from the right.
    /// Cursors are kept one past the element Java's cursor names, so none
    /// goes below zero.
    fn merge_hi(
        &mut self,
        base1: usize,
        mut len1: usize,
        base2: usize,
        mut len2: usize,
        c: &mut dyn FnMut(usize, usize) -> Ordering,
    ) -> Result<(), ContractViolation> {
        self.tmp.clear();
        self.tmp.extend_from_slice(&self.a[base2..base2 + len2]);
        let (a, tmp) = (&mut self.a, &self.tmp);
        // `end1` is Java's `cursor1 + 1`, `end2` its `cursor2 + 1`, `top`
        // its `dest + 1`.
        let mut end1 = base1 + len1;
        let mut end2 = len2;
        let mut top = base2 + len2;
        a[top - 1] = a[end1 - 1];
        top -= 1;
        end1 -= 1;
        len1 -= 1;
        if len1 == 0 {
            a[top - len2..top].copy_from_slice(&tmp[..len2]);
            return Ok(());
        }
        if len2 == 1 {
            top -= len1;
            end1 -= len1;
            a.copy_within(end1..end1 + len1, top);
            a[top - 1] = tmp[end2 - 1];
            return Ok(());
        }
        let mut min_gallop = self.min_gallop;
        'outer: loop {
            let mut count1 = 0usize;
            let mut count2 = 0usize;
            loop {
                if c(tmp[end2 - 1], a[end1 - 1]) == Ordering::Less {
                    a[top - 1] = a[end1 - 1];
                    top -= 1;
                    end1 -= 1;
                    count1 += 1;
                    count2 = 0;
                    len1 -= 1;
                    if len1 == 0 {
                        break 'outer;
                    }
                } else {
                    a[top - 1] = tmp[end2 - 1];
                    top -= 1;
                    end2 -= 1;
                    count2 += 1;
                    count1 = 0;
                    len2 -= 1;
                    if len2 == 1 {
                        break 'outer;
                    }
                }
                if (count1 | count2) >= min_gallop {
                    break;
                }
            }
            loop {
                count1 = len1 - gallop_right(tmp[end2 - 1], a, base1, len1, len1 - 1, c);
                if count1 != 0 {
                    top -= count1;
                    end1 -= count1;
                    len1 -= count1;
                    a.copy_within(end1..end1 + count1, top);
                    if len1 == 0 {
                        break 'outer;
                    }
                }
                a[top - 1] = tmp[end2 - 1];
                top -= 1;
                end2 -= 1;
                len2 -= 1;
                if len2 == 1 {
                    break 'outer;
                }
                count2 = len2 - gallop_left(a[end1 - 1], tmp, 0, len2, len2 - 1, c);
                if count2 != 0 {
                    top -= count2;
                    end2 -= count2;
                    len2 -= count2;
                    a[top..top + count2].copy_from_slice(&tmp[end2..end2 + count2]);
                    if len2 <= 1 {
                        break 'outer;
                    }
                }
                a[top - 1] = a[end1 - 1];
                top -= 1;
                end1 -= 1;
                len1 -= 1;
                if len1 == 0 {
                    break 'outer;
                }
                min_gallop = min_gallop.saturating_sub(1);
                if count1 < MIN_GALLOP && count2 < MIN_GALLOP {
                    break;
                }
            }
            min_gallop += 2;
        }
        self.min_gallop = min_gallop.max(1);
        match len2 {
            1 => {
                top -= len1;
                end1 -= len1;
                a.copy_within(end1..end1 + len1, top);
                a[top - 1] = tmp[end2 - 1];
                Ok(())
            }
            0 => Err(ContractViolation),
            _ => {
                a[top - len2..top].copy_from_slice(&tmp[..len2]);
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small xorshift generator: the tests need no dependency.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    #[test]
    fn a_total_order_sorts_stably_at_every_size() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for n in [0, 1, 2, 31, 32, 33, 64, 65, 500, 4_097, 20_000] {
            for keys in [3u64, 100, 1 << 40] {
                let mut items: Vec<(u64, usize)> = (0..n).map(|i| (rng.next() % keys, i)).collect();
                // Long runs in both directions exercise galloping.
                if n > 1000 {
                    items[100..900].sort_by_key(|x| x.0);
                    items[1000..1800].sort_by_key(|x| std::cmp::Reverse(x.0));
                }
                let mut expected = items.clone();
                expected.sort_by_key(|x| x.0);
                sort_by(&mut items, |a, b| a.0.cmp(&b.0)).expect("a total order");
                assert_eq!(items, expected, "n {n}, keys {keys}");
            }
        }
    }

    /// Any comparator, however wrong, ends in an order or the violation
    /// Java throws: never a panic, and every item is kept.
    #[test]
    fn a_comparator_that_is_not_an_order_never_panics() {
        let mut rng = Rng(42);
        let mut violations = 0;
        for round in 0..400 {
            let n = (rng.next() % 3000) as usize;
            let mut items: Vec<usize> = (0..n).collect();
            let mut answers = Rng(round + 7);
            let result = sort_by(&mut items, |_, _| match answers.next() % 3 {
                0 => Ordering::Less,
                1 => Ordering::Equal,
                _ => Ordering::Greater,
            });
            violations += usize::from(result.is_err());
            if result.is_ok() {
                let mut seen = items.clone();
                seen.sort_unstable();
                assert_eq!(seen, (0..n).collect::<Vec<_>>());
            }
        }
        assert!(
            violations > 0,
            "random answers must hit the violation check"
        );
    }

    /// Java 27's `Arrays.sort` on the same records with htsjdk's queryname
    /// comparator, whose order is not total for a record flagged both
    /// first and second: the checksum of the order Java gave, or the
    /// violation Java threw (`TimCheck.java`, run once to make this table).
    #[test]
    fn queryname_records_sort_as_java_sorts_them() {
        use crate::sam_order::QuerynameKey;
        let both = 0x1 | 0x40 | 0x80;
        let check = |names: Vec<String>, flags: Vec<u16>| -> Option<u64> {
            let mut order: Vec<usize> = (0..names.len()).collect();
            let key = |i: usize| QuerynameKey::new(names[i].as_bytes(), flags[i], None);
            sort_by(&mut order, |&a, &b| key(a).compare(&key(b))).ok()?;
            Some(
                order
                    .iter()
                    .fold(0u64, |sum, &i| (sum * 31 + i as u64) % 1_000_000_007),
            )
        };
        let a = check(
            (0..200).map(|i| format!("r{}", i % 7)).collect(),
            (0..200).map(|i| both | ((i % 2) * 0x10)).collect(),
        );
        assert_eq!(a, None);
        let b = check(
            (0..50).map(|i| format!("r{}", i % 3)).collect(),
            vec![both; 50],
        );
        assert_eq!(b, None);
        let c = check(
            (0..20).map(|i| format!("q{}", i % 2)).collect(),
            (0..20).map(|i| both | ((i % 3) * 0x10)).collect(),
        );
        assert_eq!(c, Some(599_806_246));
        let choices: [u16; 8] = [0x41, 0x81, both, 0x1, 0x0, 0x51, 0x141, 0x841];
        let expected = [
            Some(456_983_985),
            Some(663_181_677),
            Some(876_757_764),
            Some(20_839_592),
            None,
            None,
        ];
        let mut seed: u64 = 12345;
        for (round, want) in expected.iter().enumerate() {
            let n = 40 + round * 300;
            let (mut names, mut flags) = (Vec::new(), Vec::new());
            for _ in 0..n {
                seed = seed
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let r = ((seed >> 33) & 0x7fff_ffff) as usize;
                names.push(format!("n{}", r % 23));
                flags.push(choices[(r / 23) % choices.len()]);
            }
            assert_eq!(check(names, flags), *want, "round {round}");
        }
    }

    /// Java's answer for a comparator that says "less" for every pair, on
    /// fewer than 32 items: binary insertion puts each new item first, so
    /// the order is reversed after the first descending run.
    #[test]
    fn always_less_reverses_a_short_array_as_java_does() {
        let mut items: Vec<u32> = (0..10).collect();
        sort_by(&mut items, |_, _| Ordering::Less).expect("no merge");
        assert_eq!(items, (0..10).rev().collect::<Vec<_>>());
    }
}

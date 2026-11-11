//! `java.util.TreeSet` as `java.util.TreeMap` (JDK 8 to 27) builds it: a
//! red-black tree whose insertion walks from the root calling
//! `compare(new, node)`, and whose `pollFirst` takes the leftmost node.
//!
//! htsjdk's `SortingCollection` merges its spilled runs through a
//! `TreeSet` of run heads. For a total order any stable merge gives the same
//! records. For a comparator that is not one (`SAMRecordQueryNameComparator`
//! on a record flagged both first and second of its pair) the record taken
//! first depends on the shape of the tree, so the port keeps Java's
//! insertion, rotations and deletion exactly.

use std::cmp::Ordering;

struct Node<T> {
    item: Option<T>,
    left: Option<usize>,
    right: Option<usize>,
    parent: Option<usize>,
    black: bool,
}

/// The set; nodes live in one vector and freed slots are reused.
pub struct TreeSet<T> {
    nodes: Vec<Node<T>>,
    free: Vec<usize>,
    root: Option<usize>,
    len: usize,
}

impl<T> TreeSet<T> {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(capacity),
            free: Vec::new(),
            root: None,
            len: 0,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// `TreeSet.add`: false, with the set unchanged, when `compare` finds an
    /// equal item on the way down.
    pub fn add(&mut self, item: T, mut compare: impl FnMut(&T, &T) -> Ordering) -> bool {
        let Some(mut t) = self.root else {
            let node = self.alloc(item, None);
            self.nodes[node].black = true;
            self.root = Some(node);
            self.len = 1;
            return true;
        };
        let (parent, to_left) = loop {
            let ordering = compare(&item, self.item(t));
            let next = match ordering {
                Ordering::Less => self.nodes[t].left,
                Ordering::Greater => self.nodes[t].right,
                Ordering::Equal => return false,
            };
            match next {
                Some(child) => t = child,
                None => break (t, ordering == Ordering::Less),
            }
        };
        let node = self.alloc(item, Some(parent));
        if to_left {
            self.nodes[parent].left = Some(node);
        } else {
            self.nodes[parent].right = Some(node);
        }
        self.fix_after_insertion(node);
        self.len += 1;
        true
    }

    /// `TreeSet.pollFirst`: the leftmost item, removed.
    pub fn poll_first(&mut self) -> Option<T> {
        let mut p = self.root?;
        while let Some(left) = self.nodes[p].left {
            p = left;
        }
        self.delete_entry(p);
        self.len -= 1;
        self.free.push(p);
        self.nodes[p].item.take()
    }

    fn item(&self, node: usize) -> &T {
        match &self.nodes[node].item {
            Some(item) => item,
            None => unreachable!("a node in the tree holds an item"),
        }
    }

    fn alloc(&mut self, item: T, parent: Option<usize>) -> usize {
        let node = Node {
            item: Some(item),
            left: None,
            right: None,
            parent,
            black: false,
        };
        match self.free.pop() {
            Some(slot) => {
                self.nodes[slot] = node;
                slot
            }
            None => {
                self.nodes.push(node);
                self.nodes.len() - 1
            }
        }
    }

    fn parent_of(&self, x: Option<usize>) -> Option<usize> {
        x.and_then(|x| self.nodes[x].parent)
    }

    fn left_of(&self, x: Option<usize>) -> Option<usize> {
        x.and_then(|x| self.nodes[x].left)
    }

    fn right_of(&self, x: Option<usize>) -> Option<usize> {
        x.and_then(|x| self.nodes[x].right)
    }

    /// A missing node is black.
    fn is_black(&self, x: Option<usize>) -> bool {
        x.is_none_or(|x| self.nodes[x].black)
    }

    fn set_black(&mut self, x: Option<usize>, black: bool) {
        if let Some(x) = x {
            self.nodes[x].black = black;
        }
    }

    fn rotate_left(&mut self, p: Option<usize>) {
        let Some(p) = p else { return };
        let Some(r) = self.nodes[p].right else { return };
        let r_left = self.nodes[r].left;
        self.nodes[p].right = r_left;
        if let Some(rl) = r_left {
            self.nodes[rl].parent = Some(p);
        }
        let pp = self.nodes[p].parent;
        self.nodes[r].parent = pp;
        match pp {
            None => self.root = Some(r),
            Some(pp) if self.nodes[pp].left == Some(p) => self.nodes[pp].left = Some(r),
            Some(pp) => self.nodes[pp].right = Some(r),
        }
        self.nodes[r].left = Some(p);
        self.nodes[p].parent = Some(r);
    }

    fn rotate_right(&mut self, p: Option<usize>) {
        let Some(p) = p else { return };
        let Some(l) = self.nodes[p].left else { return };
        let l_right = self.nodes[l].right;
        self.nodes[p].left = l_right;
        if let Some(lr) = l_right {
            self.nodes[lr].parent = Some(p);
        }
        let pp = self.nodes[p].parent;
        self.nodes[l].parent = pp;
        match pp {
            None => self.root = Some(l),
            Some(pp) if self.nodes[pp].right == Some(p) => self.nodes[pp].right = Some(l),
            Some(pp) => self.nodes[pp].left = Some(l),
        }
        self.nodes[l].right = Some(p);
        self.nodes[p].parent = Some(l);
    }

    fn fix_after_insertion(&mut self, node: usize) {
        let mut x = Some(node);
        self.set_black(x, false);
        while x.is_some() && x != self.root && !self.is_black(self.parent_of(x)) {
            let parent = self.parent_of(x);
            let grand = self.parent_of(parent);
            if parent == self.left_of(grand) {
                let y = self.right_of(grand);
                if !self.is_black(y) {
                    self.set_black(parent, true);
                    self.set_black(y, true);
                    self.set_black(grand, false);
                    x = grand;
                } else {
                    if x == self.right_of(parent) {
                        x = parent;
                        self.rotate_left(x);
                    }
                    let parent = self.parent_of(x);
                    self.set_black(parent, true);
                    self.set_black(self.parent_of(parent), false);
                    self.rotate_right(self.parent_of(parent));
                }
            } else {
                let y = self.left_of(grand);
                if !self.is_black(y) {
                    self.set_black(parent, true);
                    self.set_black(y, true);
                    self.set_black(grand, false);
                    x = grand;
                } else {
                    if x == self.left_of(parent) {
                        x = parent;
                        self.rotate_right(x);
                    }
                    let parent = self.parent_of(x);
                    self.set_black(parent, true);
                    self.set_black(self.parent_of(parent), false);
                    self.rotate_left(self.parent_of(parent));
                }
            }
        }
        self.set_black(self.root, true);
    }

    /// `TreeMap.deleteEntry` of a node without a left child (the first).
    fn delete_entry(&mut self, p: usize) {
        let replacement = self.nodes[p].left.or(self.nodes[p].right);
        let parent = self.nodes[p].parent;
        if let Some(r) = replacement {
            self.nodes[r].parent = parent;
            match parent {
                None => self.root = Some(r),
                Some(pp) if self.nodes[pp].left == Some(p) => self.nodes[pp].left = Some(r),
                Some(pp) => self.nodes[pp].right = Some(r),
            }
            self.nodes[p].left = None;
            self.nodes[p].right = None;
            self.nodes[p].parent = None;
            if self.nodes[p].black {
                self.fix_after_deletion(r);
            }
        } else if parent.is_none() {
            self.root = None;
        } else {
            if self.nodes[p].black {
                self.fix_after_deletion(p);
            }
            if let Some(pp) = self.nodes[p].parent {
                if self.nodes[pp].left == Some(p) {
                    self.nodes[pp].left = None;
                } else if self.nodes[pp].right == Some(p) {
                    self.nodes[pp].right = None;
                }
                self.nodes[p].parent = None;
            }
        }
    }

    fn fix_after_deletion(&mut self, node: usize) {
        let mut x = Some(node);
        while x != self.root && self.is_black(x) {
            let parent = self.parent_of(x);
            if x == self.left_of(parent) {
                let mut sib = self.right_of(parent);
                if !self.is_black(sib) {
                    self.set_black(sib, true);
                    self.set_black(parent, false);
                    self.rotate_left(parent);
                    sib = self.right_of(self.parent_of(x));
                }
                if self.is_black(self.left_of(sib)) && self.is_black(self.right_of(sib)) {
                    self.set_black(sib, false);
                    x = self.parent_of(x);
                } else {
                    if self.is_black(self.right_of(sib)) {
                        self.set_black(self.left_of(sib), true);
                        self.set_black(sib, false);
                        self.rotate_right(sib);
                        sib = self.right_of(self.parent_of(x));
                    }
                    let parent = self.parent_of(x);
                    self.set_black(sib, self.is_black(parent));
                    self.set_black(parent, true);
                    self.set_black(self.right_of(sib), true);
                    self.rotate_left(parent);
                    x = self.root;
                }
            } else {
                let mut sib = self.left_of(parent);
                if !self.is_black(sib) {
                    self.set_black(sib, true);
                    self.set_black(parent, false);
                    self.rotate_right(parent);
                    sib = self.left_of(self.parent_of(x));
                }
                if self.is_black(self.right_of(sib)) && self.is_black(self.left_of(sib)) {
                    self.set_black(sib, false);
                    x = self.parent_of(x);
                } else {
                    if self.is_black(self.left_of(sib)) {
                        self.set_black(self.right_of(sib), true);
                        self.set_black(sib, false);
                        self.rotate_left(sib);
                        sib = self.left_of(self.parent_of(x));
                    }
                    let parent = self.parent_of(x);
                    self.set_black(sib, self.is_black(parent));
                    self.set_black(parent, true);
                    self.set_black(self.left_of(sib), true);
                    self.rotate_right(parent);
                    x = self.root;
                }
            }
        }
        self.set_black(x, true);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn by_value(a: &(u32, u32), b: &(u32, u32)) -> Ordering {
        a.0.cmp(&b.0).then(a.1.cmp(&b.1))
    }

    #[test]
    fn a_total_order_comes_out_sorted() {
        let mut set = TreeSet::with_capacity(8);
        let mut seed = 11u64;
        let mut out = Vec::new();
        for i in 0..2000u32 {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            assert!(set.add(((seed >> 40) as u32 % 97, i), by_value));
            if i % 5 == 4 {
                out.push(set.poll_first().expect("an item"));
            }
        }
        let mut rest = Vec::new();
        while let Some(item) = set.poll_first() {
            rest.push(item);
        }
        assert!(
            rest.windows(2)
                .all(|w| by_value(&w[0], &w[1]) == Ordering::Less)
        );
        assert!(set.is_empty());
    }

    #[test]
    fn an_equal_item_is_refused() {
        let mut set = TreeSet::with_capacity(2);
        assert!(set.add((1, 1), by_value));
        assert!(!set.add((1, 1), by_value));
        assert_eq!(set.len(), 1);
    }

    /// Two multiples of 3 each sort before the other. The expected line is
    /// what `java.util.TreeSet` printed for the same comparator, adds and
    /// polls written in Java and run on openjdk 27.0.1.
    #[test]
    fn a_comparator_that_is_not_an_order_gives_java_order() {
        fn odd(a: &(i64, i64), b: &(i64, i64)) -> Ordering {
            let r = if a.0 % 3 == 0 && b.0 % 3 == 0 {
                Ordering::Less
            } else {
                (a.0 / 2).cmp(&(b.0 / 2))
            };
            r.then(a.1.cmp(&b.1))
        }
        let mut set = TreeSet::with_capacity(8);
        let mut seed = 7i64;
        let mut out = Vec::new();
        for i in 0..40i64 {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let v = ((seed as u64) >> 33) as i64 % 30;
            set.add((v, i), odd);
            if i % 3 == 2 {
                out.push(set.poll_first().expect("an item"));
            }
        }
        while let Some(item) = set.poll_first() {
            out.push(item);
        }
        let text: Vec<String> = out.iter().map(|(v, r)| format!("{v}:{r}")).collect();
        assert_eq!(
            text.join(" "),
            "3:2 8:0 4:6 11:1 2:12 9:17 6:18 4:22 11:19 3:28 5:27 7:31 7:37 12:39 9:36 \
             13:21 15:16 9:14 15:11 16:9 16:20 19:8 18:23 20:13 20:26 20:32 23:3 22:15 \
             23:24 22:29 22:30 23:34 22:38 25:4 24:7 25:10 24:25 25:33 29:5 29:35"
        );
    }
}

//! The iteration order of a `java.util.HashSet<String>` filled by `add`
//! alone, tree bins included.
//!
//! `HighAccuracyDownsamplingIterator` copies a `HashSet` of read names into
//! an `ArrayList` before it shuffles it, so the names it keeps depend on the
//! set's order. `java.rs` models that order for tables whose bins stay
//! lists; a chunk of 10,000 names fills 16,384 slots, and a bin past 8
//! entries then turns into a red-black tree, which reorders it. This module
//! follows `HashMap.putVal` step by step: lists append, `treeifyBin`
//! doubles a table under 64 slots or builds the tree, `putTreeVal` links a
//! new node after its tree parent, `moveRootToFront` puts the root first,
//! and `resize` splits each bin in order into its low and high halves
//! (`TreeNode.split`: a half of at most 6 nodes becomes a list, a larger
//! half of a bin that split is treeified again). Keys are strings of chars
//! up to U+00FF given as bytes; a tree orders them by spread hash, then by
//! `String.compareTo`, which on such strings is the byte order.

use std::cmp::Ordering;

const TREEIFY_THRESHOLD: usize = 8;
const UNTREEIFY_THRESHOLD: usize = 6;
const MIN_TREEIFY_CAPACITY: usize = 64;

/// `String.hashCode` of chars given as bytes.
pub fn bytes_hash_code(bytes: &[u8]) -> i32 {
    bytes
        .iter()
        .fold(0i32, |h, &c| h.wrapping_mul(31).wrapping_add(i32::from(c)))
}

/// `HashMap.hash`.
fn spread(hash: i32) -> i32 {
    hash ^ ((hash as u32) >> 16) as i32
}

/// One entry: a `TreeNode` once its bin is a tree.
#[derive(Clone, Copy, Debug, Default)]
struct Node {
    hash: i32,
    parent: Option<usize>,
    left: Option<usize>,
    right: Option<usize>,
    red: bool,
}

/// One table slot: its nodes in `next` order and, for a tree bin, the root.
#[derive(Clone, Debug, Default)]
struct Bin {
    list: Vec<usize>,
    root: Option<usize>,
}

/// A `HashSet<String>` built by `add` in the order given, reduced to its
/// layout. Each key is known by its insertion number.
pub struct HashSetOrder<'k> {
    keys: Vec<&'k [u8]>,
    nodes: Vec<Node>,
    table: Vec<Bin>,
}

impl<'k> Default for HashSetOrder<'k> {
    fn default() -> Self {
        Self::new()
    }
}

impl<'k> HashSetOrder<'k> {
    /// `new HashSet<>()`: no table until the first add, then 16 slots.
    pub fn new() -> Self {
        Self {
            keys: Vec::new(),
            nodes: Vec::new(),
            table: vec![Bin::default(); 16],
        }
    }

    pub fn len(&self) -> usize {
        self.keys.len()
    }

    pub fn is_empty(&self) -> bool {
        self.keys.is_empty()
    }

    /// `add` of a key not in the set yet; the caller keeps keys distinct.
    pub fn insert(&mut self, key: &'k [u8]) {
        let hash = spread(bytes_hash_code(key));
        let id = self.keys.len();
        self.keys.push(key);
        self.nodes.push(Node {
            hash,
            ..Node::default()
        });
        let at = self.slot(hash);
        if self.table[at].root.is_some() {
            self.put_tree(at, id);
        } else {
            self.table[at].list.push(id);
            if self.table[at].list.len() > TREEIFY_THRESHOLD {
                if self.table.len() < MIN_TREEIFY_CAPACITY {
                    self.resize();
                } else {
                    self.treeify(at);
                }
            }
        }
        if self.keys.len() > self.table.len() / 4 * 3 {
            self.resize();
        }
    }

    /// The insertion numbers in iteration order.
    pub fn order(&self) -> Vec<usize> {
        self.table
            .iter()
            .flat_map(|bin| bin.list.iter().copied())
            .collect()
    }

    fn slot(&self, hash: i32) -> usize {
        (hash as u32 as usize) & (self.table.len() - 1)
    }

    /// The tree's order of node `a` against node `b`: hash, then the
    /// string (`compareComparables`).
    fn compare(&self, a: usize, b: usize) -> Ordering {
        self.nodes[a]
            .hash
            .cmp(&self.nodes[b].hash)
            .then_with(|| self.keys[a].cmp(self.keys[b]))
    }

    /// Walk from `root` to the leaf place of `x` and link it there.
    /// Returns its parent.
    fn attach(&mut self, root: usize, x: usize) -> usize {
        let mut p = root;
        loop {
            let left = self.compare(x, p) == Ordering::Less;
            let next = if left {
                self.nodes[p].left
            } else {
                self.nodes[p].right
            };
            match next {
                Some(child) => p = child,
                None => {
                    if left {
                        self.nodes[p].left = Some(x);
                    } else {
                        self.nodes[p].right = Some(x);
                    }
                    self.nodes[x].parent = Some(p);
                    return p;
                }
            }
        }
    }

    /// `putTreeVal`: the new node goes into the list right after its tree
    /// parent, then the balanced root moves to the front.
    fn put_tree(&mut self, at: usize, x: usize) {
        let root = self.table[at].root.expect("a tree bin");
        let parent = self.attach(root, x);
        let list = &mut self.table[at].list;
        let place = list.iter().position(|&n| n == parent).map_or(0, |p| p + 1);
        list.insert(place, x);
        let root = self.balance_insertion(root, x);
        self.move_root_to_front(at, root);
    }

    /// `treeifyBin` past the table check, and `TreeNode.treeify`: the list's
    /// nodes go into a new tree in list order, then the root moves first.
    fn treeify(&mut self, at: usize) {
        let list = self.table[at].list.clone();
        let mut root: Option<usize> = None;
        for &x in &list {
            self.nodes[x].left = None;
            self.nodes[x].right = None;
            self.nodes[x].parent = None;
            root = Some(match root {
                None => {
                    self.nodes[x].red = false;
                    x
                }
                Some(r) => {
                    self.nodes[x].red = true;
                    self.attach(r, x);
                    self.balance_insertion(r, x)
                }
            });
        }
        if let Some(root) = root {
            self.table[at].root = Some(root);
            self.move_root_to_front(at, root);
        }
    }

    fn move_root_to_front(&mut self, at: usize, root: usize) {
        self.table[at].root = Some(root);
        let list = &mut self.table[at].list;
        if list.first() != Some(&root)
            && let Some(place) = list.iter().position(|&n| n == root)
        {
            list.remove(place);
            list.insert(0, root);
        }
    }

    /// `resize`: double the table; each bin's nodes go, in order, to the
    /// same index or the index plus the old size.
    fn resize(&mut self) {
        let old = std::mem::take(&mut self.table);
        let old_cap = old.len();
        self.table = vec![Bin::default(); old_cap * 2];
        for (j, bin) in old.into_iter().enumerate() {
            let (lo, hi): (Vec<usize>, Vec<usize>) = bin
                .list
                .iter()
                .partition(|&&n| (self.nodes[n].hash as u32 as usize) & old_cap == 0);
            let tree = bin.root.is_some();
            for (place, half, other_empty) in [
                (j, lo.clone(), hi.is_empty()),
                (j + old_cap, hi, lo.is_empty()),
            ] {
                if half.is_empty() {
                    continue;
                }
                self.table[place].list = half;
                if !tree {
                    continue;
                }
                if self.table[place].list.len() <= UNTREEIFY_THRESHOLD {
                    self.table[place].root = None;
                } else if other_empty {
                    self.table[place].root = bin.root;
                } else {
                    self.treeify(place);
                }
            }
        }
    }

    fn rotate_left(&mut self, mut root: usize, p: usize) -> usize {
        if let Some(r) = self.nodes[p].right {
            let rl = self.nodes[r].left;
            self.nodes[p].right = rl;
            if let Some(rl) = rl {
                self.nodes[rl].parent = Some(p);
            }
            let pp = self.nodes[p].parent;
            self.nodes[r].parent = pp;
            match pp {
                None => {
                    root = r;
                    self.nodes[r].red = false;
                }
                Some(pp) if self.nodes[pp].left == Some(p) => self.nodes[pp].left = Some(r),
                Some(pp) => self.nodes[pp].right = Some(r),
            }
            self.nodes[r].left = Some(p);
            self.nodes[p].parent = Some(r);
        }
        root
    }

    fn rotate_right(&mut self, mut root: usize, p: usize) -> usize {
        if let Some(l) = self.nodes[p].left {
            let lr = self.nodes[l].right;
            self.nodes[p].left = lr;
            if let Some(lr) = lr {
                self.nodes[lr].parent = Some(p);
            }
            let pp = self.nodes[p].parent;
            self.nodes[l].parent = pp;
            match pp {
                None => {
                    root = l;
                    self.nodes[l].red = false;
                }
                Some(pp) if self.nodes[pp].right == Some(p) => self.nodes[pp].right = Some(l),
                Some(pp) => self.nodes[pp].left = Some(l),
            }
            self.nodes[l].right = Some(p);
            self.nodes[p].parent = Some(l);
        }
        root
    }

    /// `TreeNode.balanceInsertion`.
    fn balance_insertion(&mut self, mut root: usize, mut x: usize) -> usize {
        self.nodes[x].red = true;
        loop {
            let Some(xp) = self.nodes[x].parent else {
                self.nodes[x].red = false;
                return x;
            };
            let Some(xpp) = self.nodes[xp].parent.filter(|_| self.nodes[xp].red) else {
                return root;
            };
            let xppl = self.nodes[xpp].left;
            if Some(xp) == xppl {
                match self.nodes[xpp].right.filter(|&u| self.nodes[u].red) {
                    Some(xppr) => {
                        self.nodes[xppr].red = false;
                        self.nodes[xp].red = false;
                        self.nodes[xpp].red = true;
                        x = xpp;
                    }
                    None => {
                        let mut xp = xp;
                        if self.nodes[xp].right == Some(x) {
                            x = xp;
                            root = self.rotate_left(root, x);
                            xp = self.nodes[x].parent.expect("rotated under a parent");
                        }
                        self.nodes[xp].red = false;
                        if let Some(xpp) = self.nodes[xp].parent {
                            self.nodes[xpp].red = true;
                            root = self.rotate_right(root, xpp);
                        }
                    }
                }
            } else {
                match xppl.filter(|&u| self.nodes[u].red) {
                    Some(xppl) => {
                        self.nodes[xppl].red = false;
                        self.nodes[xp].red = false;
                        self.nodes[xpp].red = true;
                        x = xpp;
                    }
                    None => {
                        let mut xp = xp;
                        if self.nodes[xp].left == Some(x) {
                            x = xp;
                            root = self.rotate_right(root, x);
                            xp = self.nodes[x].parent.expect("rotated under a parent");
                        }
                        self.nodes[xp].red = false;
                        if let Some(xpp) = self.nodes[xp].parent {
                            self.nodes[xpp].red = true;
                            root = self.rotate_left(root, xpp);
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;

//! `SamFileHeaderMerger.mergeSequenceDictionaries`: the dictionaries
//! folded one by one through `mergeSequences` from an empty one.
//!
//! Each sequence of the next dictionary that the merged one lacks waits
//! until the next sequence they share, and goes in before it; the ones
//! after the last shared sequence go at the end. A shared sequence before
//! the one shared last is an error. A name the merged dictionary has keeps
//! its first record, whatever its length.
//!
//! The merged list is a linked list over borrowed records: a fold costs
//! its own dictionary plus the sequences it inserts, never the whole list.
//! Each node carries an order label, so the order test of two shared
//! sequences compares labels; the labels are spread again only when an
//! insert finds no room between two of them.

use std::collections::HashMap;

use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReferenceSequence;

use crate::error::{Error, Result};

const NIL: usize = usize::MAX;

/// The label distance a relabel leaves between two nodes.
const SPACING: u64 = 1 << 32;

struct Node<'h> {
    name: &'h [u8],
    map: &'h Map<ReferenceSequence>,
    prev: usize,
    next: usize,
    label: u64,
}

type Sequence<'h> = (&'h [u8], &'h Map<ReferenceSequence>);

struct Merged<'h> {
    nodes: Vec<Node<'h>>,
    by_name: HashMap<&'h [u8], usize>,
    head: usize,
    tail: usize,
}

/// The merged dictionary of `headers`, in input order.
pub fn merge_dictionaries(headers: &[sam::Header]) -> Result<sam::header::ReferenceSequences> {
    let mut merged = Merged {
        nodes: Vec::new(),
        by_name: HashMap::new(),
        head: NIL,
        tail: NIL,
    };
    for header in headers {
        merged.fold(header.reference_sequences())?;
    }
    Ok(merged.into_sequences())
}

impl<'h> Merged<'h> {
    /// `mergeSequences(merged, next)`. Positions in the list Java grows
    /// keep the order of the merged list before the fold, so the order test
    /// compares labels set before any insert of this fold.
    fn fold(&mut self, next: &'h sam::header::ReferenceSequences) -> Result<()> {
        let mut holder: Vec<Sequence<'h>> = Vec::new();
        let mut groups: Vec<(usize, Vec<Sequence<'h>>)> = Vec::new();
        let mut previous: Option<(usize, &[u8])> = None;
        for (name, map) in next {
            let Some(&at) = self.by_name.get(name.as_slice()) else {
                holder.push((name.as_slice(), map));
                continue;
            };
            if let Some((last, last_name)) = previous
                && self.nodes[last].label > self.nodes[at].label
            {
                return Err(Error::Metric(format!(
                    "Cannot merge sequence dictionaries because sequence {} and {} are in \
                     different orders in two input sequence dictionaries.",
                    String::from_utf8_lossy(name),
                    String::from_utf8_lossy(last_name)
                )));
            }
            if !holder.is_empty() {
                groups.push((at, std::mem::take(&mut holder)));
            }
            previous = Some((at, name.as_slice()));
        }
        for (before, group) in groups {
            self.insert_before(before, &group);
        }
        for sequence in holder {
            self.push_back(sequence);
        }
        Ok(())
    }

    fn add(&mut self, (name, map): Sequence<'h>, prev: usize, next: usize, label: u64) -> usize {
        let at = self.nodes.len();
        self.nodes.push(Node {
            name,
            map,
            prev,
            next,
            label,
        });
        self.by_name.insert(name, at);
        match prev {
            NIL => self.head = at,
            p => self.nodes[p].next = at,
        }
        match next {
            NIL => self.tail = at,
            n => self.nodes[n].prev = at,
        }
        at
    }

    fn push_back(&mut self, sequence: Sequence<'h>) {
        let label = match self.tail {
            NIL => Some(SPACING),
            tail => self.nodes[tail].label.checked_add(SPACING),
        };
        let label = label.unwrap_or_else(|| {
            self.relabel();
            self.nodes[self.tail].label + SPACING
        });
        self.add(sequence, self.tail, NIL, label);
    }

    /// The labels between the node before `at` and `at`.
    fn gap(&self, at: usize) -> (u64, u64) {
        let low = match self.nodes[at].prev {
            NIL => 0,
            p => self.nodes[p].label,
        };
        (low, self.nodes[at].label)
    }

    fn insert_before(&mut self, at: usize, group: &[Sequence<'h>]) {
        let count = group.len() as u64;
        let (mut low, mut high) = self.gap(at);
        if high - low <= count {
            self.relabel();
            (low, high) = self.gap(at);
        }
        let mut prev = self.nodes[at].prev;
        for (i, sequence) in group.iter().enumerate() {
            let step = u128::from(high - low) * (i as u128 + 1) / (u128::from(count) + 1);
            let label = low + step as u64;
            prev = self.add(*sequence, prev, at, label);
        }
    }

    /// Labels `SPACING`, `2 x SPACING`, ... in list order.
    fn relabel(&mut self) {
        let mut at = self.head;
        let mut label = SPACING;
        while at != NIL {
            self.nodes[at].label = label;
            label += SPACING;
            at = self.nodes[at].next;
        }
    }

    fn into_sequences(self) -> sam::header::ReferenceSequences {
        let mut sequences = sam::header::ReferenceSequences::default();
        let mut at = self.head;
        while at != NIL {
            let node = &self.nodes[at];
            sequences.insert(node.name.into(), node.map.clone());
            at = node.next;
        }
        sequences
    }
}

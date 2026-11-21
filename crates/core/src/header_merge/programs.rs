//! The programs `mergeProgramGroups` has yet to merge, indexed so each
//! level of the `PP` walk touches only the programs it changes.
//!
//! htsjdk translates the ID and the `PP` of every program left after each
//! level through its input's table (`translateIds`), then moves the ones
//! whose `PP` names a program just merged from the same input. A program
//! changes at a level only when its ID or `PP` is a key the table maps to
//! another ID, so [`Left`] keeps the programs by input and ID and by input
//! and `PP`, looks up the keys a level added, and carries a program whose
//! new value is again such a key to the next level. Each program is
//! translated once per change and moved once, in the list's order.

use std::collections::HashMap;

use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;

use super::{Entry, Id};

/// A program's ID or `PP` within its input.
type Key = (usize, Id);

/// The programs with a `PP` not yet merged, in htsjdk's list order.
pub(super) struct Left {
    entries: Vec<Entry<Map<Program>>>,
    /// Each program's `PP`, as translated so far.
    previous: Vec<Id>,
    alive: Vec<bool>,
    remaining: usize,
    /// Positions by (input, ID) and by (input, `PP`); a position may stay
    /// under a value it no longer has, and every lookup checks it.
    by_id: HashMap<Key, Vec<usize>>,
    by_previous: HashMap<Key, Vec<usize>>,
    /// Programs whose ID or `PP` the table already maps to another ID.
    carry: Vec<usize>,
}

/// A program's `PP`.
pub(super) fn previous_of(map: &Map<Program>) -> Option<Id> {
    map.other_fields()
        .get(&PREVIOUS_PROGRAM_ID)
        .map(|v| v.to_vec())
}

/// True when `table` maps `value` to another ID.
fn renamed(table: &HashMap<Id, Id>, value: &[u8]) -> bool {
    table.get(value).is_some_and(|new| new.as_slice() != value)
}

impl Left {
    /// The programs of `entries`, each with a `PP`, in order.
    pub(super) fn new(entries: Vec<Entry<Map<Program>>>) -> Self {
        let previous: Vec<Id> = entries
            .iter()
            .map(|entry| previous_of(&entry.map).expect("only programs with a PP are left"))
            .collect();
        let mut left = Self {
            alive: vec![true; entries.len()],
            remaining: entries.len(),
            previous,
            entries,
            by_id: HashMap::new(),
            by_previous: HashMap::new(),
            carry: Vec::new(),
        };
        for at in 0..left.entries.len() {
            left.index(at);
        }
        left
    }

    fn index(&mut self, at: usize) {
        let input = self.entries[at].input;
        self.by_id
            .entry((input, self.entries[at].id.clone()))
            .or_default()
            .push(at);
        self.by_previous
            .entry((input, self.previous[at].clone()))
            .or_default()
            .push(at);
    }

    /// `translateIds(programGroupsLeftToProcess, table, true)` after a
    /// level that mapped `changed` (input and old ID) to new IDs.
    pub(super) fn translate(&mut self, ids: &[HashMap<Id, Id>], changed: &[Key]) {
        let mut candidates = std::mem::take(&mut self.carry);
        for key in changed {
            candidates.extend(self.by_id.remove(key).unwrap_or_default());
            candidates.extend(self.by_previous.remove(key).unwrap_or_default());
        }
        candidates.sort_unstable();
        candidates.dedup();
        for at in candidates {
            if !self.alive[at] {
                continue;
            }
            let table = &ids[self.entries[at].input];
            let entry = &mut self.entries[at];
            let new_id = table
                .get(&entry.id)
                .filter(|new| **new != entry.id)
                .cloned();
            let new_previous = table
                .get(&self.previous[at])
                .filter(|new| **new != self.previous[at])
                .cloned();
            if let Some(new_id) = new_id {
                entry.id = new_id;
            }
            if let Some(new_previous) = new_previous {
                entry
                    .map
                    .other_fields_mut()
                    .insert(PREVIOUS_PROGRAM_ID, new_previous.clone().into());
                self.previous[at] = new_previous;
            }
            if renamed(table, &self.entries[at].id) || renamed(table, &self.previous[at]) {
                self.carry.push(at);
            }
            self.index(at);
        }
    }

    /// The programs of `input` whose `PP` is `id`, removed, as positions.
    pub(super) fn take_children(&mut self, input: usize, id: &[u8]) -> Vec<usize> {
        let Some(positions) = self.by_previous.remove(&(input, id.to_vec())) else {
            return Vec::new();
        };
        let mut taken = Vec::new();
        for at in positions {
            if self.alive[at] && self.previous[at] == id {
                self.alive[at] = false;
                self.remaining -= 1;
                taken.push(at);
            }
        }
        taken
    }

    pub(super) fn entry(&self, at: usize) -> &Entry<Map<Program>> {
        &self.entries[at]
    }

    /// The programs never moved: their `PP` names no program merged.
    pub(super) fn remaining(&self) -> usize {
        self.remaining
    }
}

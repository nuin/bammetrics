use std::time::{Duration, Instant};

use noodles::sam::header::record::value::map::program::tag::PREVIOUS_PROGRAM_ID;

use super::*;

#[test]
fn base36_counts_like_htsjdk() {
    assert_eq!(base36(0), "0");
    assert_eq!(base36(35), "Z");
    assert_eq!(base36(36), "10");
    assert_eq!(base36(1297), "101");
}

/// xorshift64: the random inputs of the tests below, from a seed.
struct Draw(u64);

impl Draw {
    fn below(&mut self, n: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % n as u64) as usize
    }
}

fn header(text: &str) -> sam::Header {
    text.parse().expect("header")
}

/// The merged programs and each input's ID table.
type Programs = (Vec<(Id, Map<Program>)>, Vec<HashMap<Id, Id>>);

/// htsjdk's `mergeProgramGroups` as written: every level translates and
/// scans every program left. The indexed walk must give its result.
fn reference_programs(headers: &[sam::Header]) -> std::result::Result<Programs, String> {
    let previous = |e: &Entry<Map<Program>>| {
        e.map
            .other_fields()
            .get(&PREVIOUS_PROGRAM_ID)
            .map(|v| v.to_vec())
    };
    let mut left: Vec<Entry<Map<Program>>> = headers
        .iter()
        .enumerate()
        .flat_map(|(input, h)| {
            h.programs().as_ref().iter().map(move |(id, map)| Entry {
                id: id.to_vec(),
                map: map.clone(),
                input,
            })
        })
        .collect();
    let mut current: Vec<_> = Vec::new();
    left.retain(|e| {
        let root = previous(e).is_none();
        if root {
            current.push(e.clone());
        }
        !root
    });
    let mut counter = 0;
    let mut taken: HashSet<Id> = HashSet::new();
    let mut ids = vec![HashMap::new(); headers.len()];
    let mut overall = Vec::new();
    let translate_all = |entries: Vec<Entry<Map<Program>>>, ids: &[HashMap<Id, Id>], pp: bool| {
        entries
            .into_iter()
            .map(|mut entry| {
                let table = &ids[entry.input];
                if let Some(new_id) = table.get(&entry.id) {
                    entry.id = new_id.clone();
                }
                if pp && let Some(new_pp) = previous(&entry).and_then(|pp| table.get(&pp).cloned())
                {
                    entry
                        .map
                        .other_fields_mut()
                        .insert(PREVIOUS_PROGRAM_ID, new_pp.into());
                }
                entry
            })
            .collect::<Vec<_>>()
    };
    while !current.is_empty() {
        let mut by_id: Vec<(Id, Distinct<Map<Program>>)> = Vec::new();
        for entry in &current {
            let records = match by_id.iter_mut().find(|(id, _)| *id == entry.id) {
                Some((_, records)) => records,
                None => {
                    by_id.push((entry.id.clone(), Vec::new()));
                    &mut by_id.last_mut().expect("pushed").1
                }
            };
            match records.iter_mut().find(|(m, _)| *m == entry.map) {
                Some((_, inputs)) => inputs.push(entry.input),
                None => records.push((entry.map.clone(), vec![entry.input])),
            }
        }
        for (id, records) in by_id {
            for (map, inputs) in records {
                let new_id = if taken.insert(id.clone()) {
                    counter += 1;
                    id.clone()
                } else {
                    loop {
                        let candidate =
                            format!("{}.{}", String::from_utf8_lossy(&id), base36(counter))
                                .into_bytes();
                        counter += 1;
                        if taken.insert(candidate.clone()) {
                            break candidate;
                        }
                    }
                };
                for input in inputs {
                    ids[input].insert(id.clone(), new_id.clone());
                }
                overall.push((new_id, map));
            }
        }
        current = translate_all(current, &ids, false);
        left = translate_all(left, &ids, true);
        let mut next = Vec::new();
        let mut i = 0;
        while i < left.len() {
            let pp = previous(&left[i]).expect("PP");
            if current
                .iter()
                .any(|done| done.input == left[i].input && done.id == pp)
            {
                next.push(left.remove(i));
            } else {
                i += 1;
            }
        }
        current = next;
    }
    if !left.is_empty() {
        return Err(format!("{} left", left.len()));
    }
    overall.sort_by(|a, b| a.0.cmp(&b.0));
    Ok((overall, ids))
}

/// Programs whose IDs collide across inputs and with the renamed forms
/// (`A.1`), with `PP` chains, forks and dangling `PP`s.
fn random_programs(draw: &mut Draw) -> Vec<sam::Header> {
    const IDS: [&str; 7] = ["A", "B", "A.1", "A.2", "A.1.3", "B.4", "C"];
    (0..1 + draw.below(4))
        .map(|_| {
            let mut text = String::from("@HD\tVN:1.6\n");
            let mut ids: Vec<&str> = IDS.to_vec();
            let count = 1 + draw.below(IDS.len());
            let mut chosen: Vec<&str> = Vec::new();
            for _ in 0..count {
                let id = ids.swap_remove(draw.below(ids.len()));
                text.push_str(&format!("@PG\tID:{id}\tPN:p{}", draw.below(2)));
                // Mostly a program listed before it, so chains end at a root;
                // now and then any ID, which may dangle or loop.
                let pp = match draw.below(20) {
                    0 => Some(IDS[draw.below(IDS.len())]),
                    n if n < 13 && !chosen.is_empty() => Some(chosen[draw.below(chosen.len())]),
                    _ => None,
                };
                if let Some(pp) = pp {
                    text.push_str(&format!("\tPP:{pp}"));
                }
                text.push('\n');
                chosen.push(id);
            }
            header(&text)
        })
        .collect()
}

#[test]
fn the_indexed_program_walk_equals_htsjdks_scan() {
    let mut merged_runs = 0;
    for seed in 1..=3_000u64 {
        let mut draw = Draw(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        let headers = random_programs(&mut draw);
        let expected = reference_programs(&headers);
        let got = HeaderMerger::merge(&headers);
        match (expected, got) {
            (Ok((programs, ids)), Ok(merger)) => {
                merged_runs += 1;
                let got: Vec<(Id, Map<Program>)> = merger
                    .header
                    .programs()
                    .as_ref()
                    .iter()
                    .map(|(id, map)| (id.to_vec(), map.clone()))
                    .collect();
                assert_eq!(got, programs, "seed {seed}");
                assert_eq!(merger.program_ids, ids, "seed {seed}");
            }
            (Err(_), Err(_)) => {}
            (expected, got) => panic!(
                "seed {seed}: htsjdk's walk {:?}, the indexed walk {:?}",
                expected.is_ok(),
                got.is_ok()
            ),
        }
    }
    assert!(merged_runs > 2_000, "{merged_runs} seeds merged");
}

/// `mergeSequences` as htsjdk writes it, over a list it grows.
fn reference_dictionaries(headers: &[sam::Header]) -> std::result::Result<Vec<Vec<u8>>, String> {
    let mut merged: Vec<Vec<u8>> = Vec::new();
    for header in headers {
        let mut result = merged.clone();
        let mut previous: Option<usize> = None;
        let mut holder: Vec<Vec<u8>> = Vec::new();
        for name in header.reference_sequences().keys() {
            match result.iter().position(|n| n.as_slice() == name.as_slice()) {
                None => holder.push(name.to_vec()),
                Some(at) if previous.is_some_and(|p| p > at) => return Err("order".to_owned()),
                Some(at) => {
                    let count = holder.len();
                    result.splice(at..at, holder.drain(..));
                    previous = Some(at + count);
                }
            }
        }
        result.append(&mut holder);
        merged = result;
    }
    Ok(merged)
}

#[test]
fn the_linked_dictionary_merge_equals_htsjdks_list() {
    for seed in 1..=3_000u64 {
        let mut draw = Draw(seed.wrapping_mul(0xD1B5_4A32_D192_ED03));
        let headers: Vec<sam::Header> = (0..1 + draw.below(6))
            .map(|_| {
                let mut names: Vec<usize> = (0..12).collect();
                let mut text = String::new();
                for _ in 0..draw.below(8) {
                    let name = names.swap_remove(draw.below(names.len()));
                    text.push_str(&format!("@SQ\tSN:c{name}\tLN:{}\n", 100 + name));
                }
                header(&text)
            })
            .collect();
        let expected = reference_dictionaries(&headers);
        let got = merge_dictionaries(&headers)
            .map(|merged| merged.keys().map(|name| name.to_vec()).collect::<Vec<_>>());
        match (expected, got) {
            (Ok(expected), Ok(got)) => assert_eq!(got, expected, "seed {seed}"),
            (Err(_), Err(_)) => {}
            (expected, got) => panic!("seed {seed}: {expected:?} against {:?}", got.is_ok()),
        }
    }
}

#[test]
fn inserts_at_one_place_spread_the_labels_again() {
    // Each input puts one new contig before the shared `s`, so every fold
    // halves the same gap until a relabel.
    let headers: Vec<sam::Header> = (0..200)
        .map(|i| header(&format!("@SQ\tSN:n{i}\tLN:10\n@SQ\tSN:s\tLN:10\n")))
        .collect();
    let merged = merge_dictionaries(&headers).expect("merge");
    let names: Vec<String> = merged
        .keys()
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();
    let mut expected: Vec<String> = (0..200).map(|i| format!("n{i}")).collect();
    expected.push("s".to_owned());
    assert_eq!(names, expected);
}

/// Time bounds far above a linear merge and far below a quadratic one in
/// a debug build.
const BOUND: Duration = Duration::from_secs(20);

#[test]
fn a_hundred_thousand_read_groups_merge_in_linear_time() {
    let mut text = String::from("@HD\tVN:1.6\n");
    for i in 0..100_000 {
        text.push_str(&format!("@RG\tID:rg{i}\tSM:s\tLB:l{}\n", i % 7));
    }
    let headers = [header(&text), header(&text)];
    let start = Instant::now();
    let merger = HeaderMerger::merge(&headers).expect("merge");
    assert!(start.elapsed() < BOUND, "{:?}", start.elapsed());
    assert_eq!(merger.header.read_groups().len(), 100_000);
}

#[test]
fn a_chain_of_twenty_thousand_programs_merges_in_linear_time() {
    let chain = |pn: &str| {
        let mut text = String::from("@HD\tVN:1.6\n@PG\tID:p0\tPN:x\n");
        for i in 1..20_000 {
            text.push_str(&format!("@PG\tID:p{i}\tPN:{pn}\tPP:p{}\n", i - 1));
        }
        header(&text)
    };
    // The second chain differs, so each of its programs is renamed.
    let headers = [chain("x"), chain("y")];
    let start = Instant::now();
    let merger = HeaderMerger::merge(&headers).expect("merge");
    assert!(start.elapsed() < BOUND, "{:?}", start.elapsed());
    assert_eq!(merger.header.programs().as_ref().len(), 39_999);
}

#[test]
fn many_inputs_of_new_contigs_merge_in_linear_time() {
    let headers: Vec<sam::Header> = (0..300)
        .map(|i| {
            let mut text = String::new();
            for c in 0..300 {
                text.push_str(&format!("@SQ\tSN:i{i}c{c}\tLN:10\n"));
            }
            text.push_str("@SQ\tSN:shared\tLN:10\n");
            header(&text)
        })
        .collect();
    let start = Instant::now();
    let merged = merge_dictionaries(&headers).expect("merge");
    assert!(start.elapsed() < BOUND, "{:?}", start.elapsed());
    assert_eq!(merged.len(), 90_001);
    assert_eq!(
        merged.get_index(300).map(|(name, _)| name.to_vec()),
        Some(b"i1c0".to_vec())
    );
}

use std::io;

use super::indexer::{region_to_bin, window};
use super::*;

fn offset(block: u64, within: u64) -> u64 {
    block << 16 | within
}

fn mapped(reference: i32, start: i32, length: i32, chunk: Chunk) -> Alignment {
    Alignment {
        reference,
        start,
        end: start + length - 1,
        unmapped: false,
        chunk,
    }
}

fn build(lengths: &[usize], alignments: &[Alignment]) -> (Index, Counts) {
    let mut indexer = Indexer::new(Vec::new(), lengths.to_vec()).expect("memory");
    for alignment in alignments {
        indexer
            .add(alignment)
            .expect("memory")
            .unwrap_or_else(|refusal| panic!("{alignment:?}: {refusal}"));
    }
    let (bytes, counts) = indexer.finish().expect("memory");
    let index = Index::read(&mut bytes.as_slice()).expect("an index");
    let mut again = Vec::new();
    index.write(&mut again).expect("memory");
    assert_eq!(again, bytes, "read then written is the same file");
    (index, counts)
}

fn refusal(lengths: &[usize], alignments: &[Alignment]) -> Refusal {
    let mut indexer = Indexer::new(Vec::new(), lengths.to_vec()).expect("memory");
    alignments
        .iter()
        .find_map(|a| indexer.add(a).expect("memory").err())
        .expect("a refusal")
}

#[test]
fn bins_follow_the_sam_specification() {
    assert_eq!(region_to_bin(0, 1), 4681);
    assert_eq!(region_to_bin(16383, 16384), 4681);
    assert_eq!(region_to_bin(16383, 16385), 585);
    assert_eq!(region_to_bin(16384, 16385), 4682);
    assert_eq!(region_to_bin(0, 1 << 17), 585);
    assert_eq!(region_to_bin(0, (1 << 17) + 1), 73);
    assert_eq!(region_to_bin(0, (1 << 26) + 1), 0);
    assert_eq!(region_to_bin((1 << 29) - 1, 1 << 29), 37448);
    assert_eq!(window(0), 0);
    assert_eq!(window(1), 0);
    assert_eq!(window(16384), 0);
    assert_eq!(window(16385), 1);
}

#[test]
fn an_index_without_records_is_zero_bytes_per_reference() {
    let mut indexer = Indexer::new(Vec::new(), vec![1000, 2000]).expect("memory");
    let unplaced = Alignment {
        reference: -1,
        start: 0,
        end: 0,
        unmapped: true,
        chunk: Chunk { start: 1, end: 2 },
    };
    indexer.add(&unplaced).expect("memory").expect("indexed");
    let (bytes, counts) = indexer.finish().expect("memory");
    let mut expected = b"BAI\x01\x02\0\0\0".to_vec();
    expected.extend([0u8; 16]);
    expected.extend(1u64.to_le_bytes());
    assert_eq!(bytes, expected);
    assert_eq!(counts.references, [(0, 0), (0, 0)]);
    assert_eq!((counts.no_coordinate, counts.records), (1, 1));
    let none = Indexer::new(Vec::new(), Vec::new())
        .and_then(Indexer::finish)
        .expect("memory");
    assert_eq!(none.0, b"BAI\x01\0\0\0\0\0\0\0\0\0\0\0\0");
}

#[test]
fn a_chunk_grows_within_a_block_and_a_new_block_starts_one() {
    let chunk = |a: (u64, u64), b: (u64, u64)| Chunk {
        start: offset(a.0, a.1),
        end: offset(b.0, b.1),
    };
    let (index, counts) = build(
        &[100_000],
        &[
            mapped(0, 1, 50, chunk((100, 0), (100, 90))),
            mapped(0, 10, 50, chunk((100, 90), (100, 180))),
            // The record ends block 100, so it ends at block 300.
            mapped(0, 20, 50, chunk((100, 180), (300, 0))),
            mapped(0, 30, 50, chunk((300, 0), (300, 90))),
            // The record runs from block 300 into block 500.
            mapped(0, 40, 50, chunk((300, 90), (500, 10))),
            mapped(0, 50, 50, chunk((500, 10), (500, 100))),
            mapped(0, 20_000, 50, chunk((500, 100), (700, 0))),
            mapped(0, 20_001, 50, chunk((700, 0), (900, 0))),
        ],
    );
    let reference = &index.references[0];
    assert_eq!(
        reference.bins,
        [
            (4681, vec![chunk((100, 0), (500, 100))]),
            (4682, vec![chunk((500, 100), (900, 0))]),
        ]
    );
    assert_eq!(
        reference.metadata,
        Some(Metadata {
            first_offset: offset(100, 0),
            last_offset: offset(900, 0),
            aligned: 8,
            unaligned: 0,
        })
    );
    assert_eq!(reference.intervals, [offset(100, 0), offset(500, 100)]);
    assert_eq!(counts.references, [(8, 0)]);
    assert_eq!(index.no_coordinate, Some(0));
}

#[test]
fn a_bin_met_again_joins_its_last_chunk_only_in_the_same_block() {
    let chunk = |a: (u64, u64), b: (u64, u64)| Chunk {
        start: offset(a.0, a.1),
        end: offset(b.0, b.1),
    };
    let spanning = |start: i32, chunk: Chunk| mapped(0, start, 20_000, chunk);
    let (index, _) = build(
        &[100_000],
        &[
            spanning(1, chunk((100, 0), (100, 90))),
            mapped(0, 2, 50, chunk((100, 90), (100, 180))),
            // Bin 585 again, in the block its chunk ended in: one chunk,
            // which now also covers the record of bin 4681 between.
            spanning(3, chunk((100, 180), (100, 270))),
            mapped(0, 4, 50, chunk((100, 270), (300, 0))),
            mapped(0, 5, 50, chunk((300, 0), (500, 0))),
            // Bin 585 again, two blocks after its chunk ended.
            spanning(6, chunk((500, 0), (500, 90))),
            // A block one byte after the last one counts as the same.
            spanning(7, chunk((501, 0), (501, 90))),
        ],
    );
    assert_eq!(
        index.references[0].bins,
        [
            (
                585,
                vec![chunk((100, 0), (100, 270)), chunk((500, 0), (501, 90))]
            ),
            (4681, vec![chunk((100, 90), (500, 0))]),
        ]
    );
}

#[test]
fn windows_between_records_take_the_offset_before_them() {
    let at = |within: u64| Chunk {
        start: offset(7, within),
        end: offset(7, within + 10),
    };
    let unmapped = |start: i32, chunk: Chunk| Alignment {
        reference: 0,
        start,
        end: 0,
        unmapped: true,
        chunk,
    };
    let (index, counts) = build(
        &[200_000],
        &[
            // Window 1 only; window 0 stays 0.
            mapped(0, 16_385, 50, at(0)),
            // Windows 3 and 4; window 2 takes window 1's offset.
            mapped(0, 65_500, 100, at(10)),
            // An unmapped record on a window edge counts in the window
            // of the base before it: 5, not 6.
            unmapped(98_305, at(20)),
            unmapped(98_306, at(30)),
        ],
    );
    assert_eq!(
        index.references[0].intervals,
        [
            0,
            offset(7, 0),
            offset(7, 0),
            offset(7, 10),
            offset(7, 10),
            offset(7, 20),
            offset(7, 30),
        ]
    );
    assert_eq!(counts.references, [(2, 2)]);
    let bins: Vec<u32> = index.references[0].bins.iter().map(|b| b.0).collect();
    assert_eq!(bins, [585, 4682, 4687]);
}

#[test]
fn a_reference_is_written_when_a_later_one_starts() {
    let chunk = Chunk {
        start: offset(9, 0),
        end: offset(9, 50),
    };
    let no_position = Alignment {
        reference: 2,
        start: 0,
        end: 0,
        unmapped: true,
        chunk,
    };
    let (index, counts) = build(
        &[1000, 1000, 1000, 1000],
        &[mapped(1, 5, 10, chunk), no_position],
    );
    assert_eq!(index.references.len(), 4);
    assert!(index.references[0].bins.is_empty());
    assert_eq!(index.references[1].bins.len(), 1);
    assert!(index.references[2].metadata.is_none());
    assert_eq!(counts.references, [(0, 0), (1, 0), (0, 0), (0, 0)]);
    assert_eq!(counts.no_coordinate, 1);
}

#[test]
fn records_picard_refuses_are_refused() {
    let chunk = Chunk { start: 1, end: 2 };
    assert_eq!(
        refusal(
            &[1000, 1000],
            &[mapped(1, 5, 10, chunk), mapped(0, 5, 10, chunk)]
        ),
        Refusal::UnexpectedReference {
            reference: 0,
            indexing: 1
        }
    );
    assert_eq!(
        refusal(&[1000], &[mapped(-1, 5, 10, chunk)]).to_string(),
        "Unexpected reference -1 when constructing index for 0"
    );
    assert_eq!(
        refusal(&[1000], &[mapped(3, 5, 10, chunk)]),
        Refusal::UnknownReference(3)
    );
    let top = 1 << 29;
    let long = [(1 << 29) + 40_000];
    build(
        &long,
        &[mapped(0, top - 49, 50, chunk), mapped(0, top, 1, chunk)],
    );
    assert_eq!(
        refusal(&long, &[mapped(0, top - 48, 50, chunk)]),
        Refusal::PositionTooHigh
    );
    assert_eq!(
        refusal(&[2000], &[mapped(0, 40_000, 50, chunk)]),
        Refusal::PastReference {
            bin: 4683,
            last: 4681
        }
    );
}

#[test]
fn a_summary_reads_the_counts_and_no_bin() {
    let chunk = Chunk {
        start: offset(3, 0),
        end: offset(3, 9),
    };
    let mut indexer = Indexer::new(Vec::new(), vec![1000, 1000, 1000]).expect("memory");
    for alignment in [mapped(0, 1, 5, chunk), mapped(2, 1, 5, chunk)] {
        indexer.add(&alignment).expect("memory").expect("indexed");
    }
    let (bytes, _) = indexer.finish().expect("memory");
    let summary = Summary::read(&mut bytes.as_slice(), 2).expect("a summary");
    assert_eq!(summary.references, 3);
    assert_eq!(summary.metadata.len(), 2);
    assert_eq!(summary.metadata[0].map(|m| m.aligned), Some(1));
    assert_eq!(summary.metadata[1], None);
    assert_eq!(summary.no_coordinate, Some(0));

    // Without the last 8 bytes the count is unknown; cut further, the
    // index is an error.
    let cut = &bytes[..bytes.len() - 8];
    assert_eq!(
        Summary::read(&mut &cut[..], 3)
            .expect("a summary")
            .no_coordinate,
        None
    );
    let short = &bytes[..bytes.len() - 20];
    let error = Summary::read(&mut &short[..], 3).expect_err("cut inside a reference");
    assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
    assert!(Index::read(&mut &short[..]).is_err());
    assert!(Index::read(&mut &b"BAM\x01\0\0\0\0"[..]).is_err());
}

#[test]
fn a_count_the_file_does_not_back_is_an_early_end() {
    // 2^31 - 1 references, bins and chunks, and no byte behind them.
    let mut bytes = b"BAI\x01".to_vec();
    bytes.extend(i32::MAX.to_le_bytes());
    bytes.extend(i32::MAX.to_le_bytes());
    bytes.extend(0u32.to_le_bytes());
    bytes.extend(i32::MAX.to_le_bytes());
    for keep in [true, false] {
        let error = if keep {
            Index::read(&mut bytes.as_slice()).map(|_| ())
        } else {
            Summary::read(&mut bytes.as_slice(), usize::MAX).map(|_| ())
        };
        assert_eq!(
            error.expect_err("no chunk follows").kind(),
            io::ErrorKind::UnexpectedEof
        );
    }
    let mut negative = b"BAI\x01".to_vec();
    negative.extend((-1i32).to_le_bytes());
    assert!(Index::read(&mut negative.as_slice()).is_err());
}

#[test]
fn a_reference_starts_without_the_chunks_of_the_one_before() {
    let chunk = |block: u64| Chunk {
        start: offset(block, 0),
        end: offset(block, 50),
    };
    let (index, counts) = build(
        &vec![100_000; 1000],
        &[
            mapped(0, 1, 10, chunk(5)),
            mapped(0, 20_000, 10, chunk(7)),
            mapped(999, 1, 10, chunk(9)),
        ],
    );
    assert_eq!(index.references.len(), 1000);
    assert_eq!(index.references[0].bins.len(), 2);
    assert!(
        index.references[1..999]
            .iter()
            .all(|r| *r == Reference::default())
    );
    assert_eq!(index.references[999].bins, [(4681, vec![chunk(9)])]);
    assert_eq!(index.references[999].intervals, [offset(9, 0)]);
    assert_eq!(counts.references[999], (1, 0));
    assert_eq!(counts.references.iter().filter(|c| c.0 > 0).count(), 2);
}

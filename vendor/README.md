# Vendored crates

`[patch.crates-io]` in the root `Cargo.toml` builds these in place of the
crates.io releases. `exclude = ["vendor"]` keeps them out of the workspace, so
`cargo fmt`, `cargo clippy` and `cargo test --workspace` do not run on them.

## noodles-cram 0.102.0

The crates.io tarball (`noodles-cram-0.102.0.crate`), unpacked, with
`Cargo.lock` dropped, an empty `[workspace]` table added to `Cargo.toml`, and
ten changes marked `bammetrics patch`: six in the reader, three in the writer,
one in the rANS codec. The first is in
`src/io/reader/container/slice.rs`:

- `calculate_template_length` is htsjdk 4's
  `CRAMCompressionRecord.computeInsertSize`, which Picard 3.5.0 decodes with:
  0 when either read is unmapped or the two lie on different references,
  otherwise the span of both reads, positive on the read that starts first
  and, at equal starts, on the first of pair. noodles gave the plus sign to
  the read that comes first in the file and ignored unmapped reads and
  references.
- `resolve_mates` sets the template length on the first and the last
  segment of a chain of attached mates and leaves a middle segment at 0, as
  htsjdk's `restoreMateInfo` does. noodles gave every later segment the
  negated length.

The tests of that file cover both. `cargo test --manifest-path
vendor/noodles-cram/Cargo.toml` runs the crate's tests.

Without the change, a mate pair samtools stored attached at equal starts and
read without an `MC` tag gets the opposite orientation from Picard
(`tests/concordance/test_cram.py::test_attached_mates_without_mc_tag`).

The second is in `src/io/reader/container.rs`: `read_container` grows the
container buffer one MiB at a time up to the container length, in place of one
`resize` to that length. A container header (its CRC32 valid) that claims
0x7fffff00 bytes in a file of 343 bytes cost 1,578 MiB peak RSS in `insert-size`;
it now ends as `truncated CRAM` at 5 MiB
(`crates/cli/tests/forged_block_size.rs`).

The third is in `src/record/sequence/iter.rs`: a read that runs past the end of
its reference sequence reads each reference base there as `N`, as htsjdk 5.0.0's
`CRAMRecordReadFeatures.restoreReadBases` does. noodles panicked on the slice
index (`range start index 701 out of range for slice of length 700`).

The fourth is in `src/io/reader/container/slice/records.rs`: `read_sequence`
and `read_quality_scores` return no bytes for a read of length 0 before they
look up the `BA` or `QS` encoding and block. A container of unmapped reads with
SEQ `*` alone has neither, and noodles failed with `missing data series
encoding: Bases`. Picard reads both shapes; `ViewSam` on each prints Picard's
text (`tests/concordance/test_replace_sam_header.py::test_cram_shapes_as_htsjdk_reads_them`).

The fifth is in `src/codecs/name_tokenizer/decode.rs`: `TokenReader::get`
returns the `DZLen` stream, so a token stream that duplicates one decodes.
noodles had `DZLen` in `set` and `get_mut` but not in `get`, and panicked with
`unhandled ty: DZLen` on a CRAM 3.1 container samtools 1.25 wrote. An unknown
stream type in `get` or `set`, a stream before the first token and a duplicate
of a token not read yet are now `InvalidData` errors in place of panics.

The sixth is the writer in pieces, `src/io/writer/raw.rs` (new, the whole file):
`write_file_definition`, `write_header_text`, `encode_container` and
`write_eof`, with `Record`, `Feature`, `Flags`, `MateFlags`, `Base`,
`SliceContext` and `Options`. `crates/core/src/cram_writer.rs` is its one
caller. noodles' own `Writer` converts each record against a
`fasta::Repository` that keeps every contig it has read, attaches mates by read
name alone and writes no mate flags for a detached record, cuts a container
every 10,240 records whatever their reference, panics on a read past the end of
its contig, and asks an `M5` of every `@SQ` line. htsjdk 5.0.0 writes none and
cuts its slices at a change of reference. The module takes records the caller
converted, one slice per container with the caller's reference context and
slice MD5, and the header text as given. It touches these files beside its own:

- `src/io/writer.rs`: `pub mod raw`.
- `src/io/writer/record.rs`: the fields of `Record` are public.
- `src/record.rs`: `Flags` and `MateFlags` are exported.
- `src/io/writer/container.rs`: `encode_prepared_container`.
- `src/io/writer/container/slice.rs`: `build_prepared_slice`.
- `src/io/writer/header.rs` and `header/container.rs`: `write_text_container`,
  the header container without the `M5` check.
- `#![allow(missing_docs)]` at the top of the six files whose items the module
  exports (`io/writer/record.rs`, `record/feature.rs`, `record/convert.rs`,
  `record/flags.rs`, `record/mate_flags.rs`,
  `container/compression_header/preservation_map/substitution_matrix/base.rs`).

The seventh makes the writer's bytes a function of its records, and its tag
lookups linear:

- `src/io/writer/container/slice.rs` `write_records` sorts the external blocks
  by content ID, and `compression_header/tag_encodings.rs` writes the tag
  encodings in that order. noodles wrote both in `HashMap` order, so two runs
  on the same records gave different bytes.
- `compression_header/preservation_map/tag_sets.rs` `build_tag_sets` keeps the
  sets seen in a `HashSet`, and `slice/records.rs` `write_data` finds a
  record's set in a `HashMap`. noodles compared each record's tag set with
  every set of the container, twice.

The eighth is in `src/io/writer/container/slice/records.rs`: `write_data`
stores a `Z` or `H` tag value as the bytes given and a NUL, as a BAM stores it.
noodles' BAM value encoder refused a `Z` byte outside ` ` to `~` and an `H`
text that is not pairs of `0-9A-F` (`invalid string`, `invalid hex`), so one
record with `XZ:Z:café` in UTF-8 failed the whole file. A BAM holds such texts
and Picard and samtools write a CRAM of them. A text with a NUL is an error.

The ninth is in `src/codecs/rans_4x8/encode/order_0.rs`, `normalize_frequencies`,
which the order-1 encoder of the quality scores calls per context:

- The product `count * 4095` is taken in 64 bits. In 32 bits it wrapped for a
  symbol seen more than 1,048,832 times in a block, a panic in a debug build
  and a table that no longer followed the counts in a release build: three
  symbols seen 1,400,000 times each got 3411, 342 and 342 of 4095 in place of
  1365 each. The file still decoded. A slice of 10,000 reads of 20,000 bases
  is such a block; `cram_writer.rs` now also ends a slice at 8 MiB of BAM
  records, which keeps most blocks under that count.
- What the normalized sum is over 4095 is taken from the most frequent symbol
  and, when that leaves it at 1, from the next. noodles took it all from the
  most frequent symbol, whose frequency went below 0 for a block of 60 common
  and 196 rare symbols (181 over, the largest frequency 68).

The tenth is in `src/record.rs`: `Record::mask_texts`, for the reader in
`crates/core/src/cram.rs`, which turns each CRAM record into a BAM record
through noodles-bam. Its encoder refuses a read name with a byte outside `!` to
`~` or an `@`, and the tag texts of the eighth change. `mask_texts` replaces
each such text with as many `x` bytes and returns the bytes, and
`bam_raw::restore_texts` puts them back into the BAM record. Without it
bammetrics read no CRAM that holds such a name or tag, its own, Picard's or
samtools'.

`cargo test --manifest-path vendor/noodles-cram/Cargo.toml` runs 317 tests, one
of them for the third change, one for the fifth, two for the sixth, one for the
eighth and two for the ninth. `cargo test -p bammetrics-core cram_writer` holds
the tenth.

To move to a later noodles-cram: unpack its tarball here, apply the
`bammetrics patch` changes again if upstream has not, and run that test module,
`crates/cli/tests/forged_block_size.rs` and `cargo test -p bammetrics-core
cram_writer`.

# Vendored crates

`[patch.crates-io]` in the root `Cargo.toml` builds these in place of the
crates.io releases. `exclude = ["vendor"]` keeps them out of the workspace, so
`cargo fmt`, `cargo clippy` and `cargo test --workspace` do not run on them.

## noodles-cram 0.102.0

The crates.io tarball (`noodles-cram-0.102.0.crate`), unpacked, with
`Cargo.lock` dropped, an empty `[workspace]` table added to `Cargo.toml`, and
five changes marked `bammetrics patch`. The first is in
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

`cargo test --manifest-path vendor/noodles-cram/Cargo.toml` runs 312 tests, one
of them for the third change and one for the fifth.

To move to a later noodles-cram: unpack its tarball here, apply the
`bammetrics patch` changes again if upstream has not, and run that test module
and `crates/cli/tests/forged_block_size.rs`.

# Vendored crates

`[patch.crates-io]` in the root `Cargo.toml` builds these in place of the
crates.io releases. `exclude = ["vendor"]` keeps them out of the workspace, so
`cargo fmt`, `cargo clippy` and `cargo test --workspace` do not run on them.

## noodles-cram 0.102.0

The crates.io tarball (`noodles-cram-0.102.0.crate`), unpacked, with
`Cargo.lock` dropped, an empty `[workspace]` table added to `Cargo.toml`, and
one change in `src/io/reader/container/slice.rs`, marked `bammetrics patch`:

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
vendor/noodles-cram/Cargo.toml` runs the crate's 310 tests.

Without the change, a mate pair samtools stored attached at equal starts and
read without an `MC` tag gets the opposite orientation from Picard
(`tests/concordance/test_cram.py::test_attached_mates_without_mc_tag`).

To move to a later noodles-cram: unpack its tarball here, apply the two
changes again if upstream has not, and run that test module.

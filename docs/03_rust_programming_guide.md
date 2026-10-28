# Rust for the Experienced Programmer — a Detailed Guide

Written for someone fluent in Python/TypeScript/Swift who needs to *own* a Rust codebase, not just paste snippets. Every concept lands on the BAM/genomics domain you're building (`bammetrics`). All non-trivial examples here were compiled on rustc 1.75; the two compiler-error sections quote real output.

Read it once top-to-bottom for the mental model, then use it as a reference. The single most important section is **2–4 (ownership/borrowing/lifetimes)** — everything else is comparatively familiar.

---

## 0. The mental-model shift

You already know variables, functions, types, generics, closures, and iterators from Python/TS/Swift. Three things are genuinely new:

1. **Ownership.** Every value has exactly one owner. When the owner goes out of scope, the value is freed — no GC, no manual `free`. Assigning or passing a value usually *moves* ownership; the old binding becomes unusable. This is the thing that feels alien for a week and then disappears into reflex.
2. **Borrowing.** Instead of moving, you lend a reference: `&T` (shared, many readers) or `&mut T` (exclusive, one writer). The compiler enforces that you never have an exclusive borrow alongside any other borrow. This is what makes data races a compile error.
3. **The compiler is a pair-programmer, not a gate.** Rust's errors are unusually specific and almost always tell you the fix. You will fight the borrow checker early; the fight *is* the learning. Closest analogies: Swift's value semantics and ARC give you intuition for moves; TypeScript's structural typing is the opposite of Rust's nominal `trait` system, so unlearn that.

What you can carry over directly: expressions-as-values (like Rust, the last expression in a block is its value — no `return` needed), `match` ≈ Swift's `switch` with exhaustiveness, `Option`/`Result` ≈ Swift's optionals and `throws`, iterators ≈ Python generators but lazy and zero-cost.

---

## 1. Project mechanics (cargo, modules, crates)

`cargo` is npm + pip + the build system in one. A *crate* is a compilation unit (a library or binary); a *package* (one `Cargo.toml`) ships one or more crates; a *workspace* groups packages.

```
cargo new --lib foo        # library crate
cargo new foo              # binary crate (has src/main.rs)
cargo check                # typecheck without codegen — your fastest feedback loop
cargo build [--release]    # debug is fast-compile/slow-run; release is the opposite
cargo run -- arg1 arg2     # build+run a binary, args after --
cargo test                 # run #[test] functions
cargo doc --open           # build & view API docs for your deps (invaluable for noodles)
cargo add noodles --features bam,sam    # edit Cargo.toml for you
```

**Modules** are namespaces, declared in code, not by file convention alone. `mod histogram;` in `lib.rs` tells the compiler to pull in `src/histogram.rs`. Items are private by default; `pub` exposes them. Paths use `::`. `use` imports a path into scope.

```rust
// src/lib.rs
mod histogram;             // -> src/histogram.rs
mod insert_size;           // -> src/insert_size.rs
pub use insert_size::Row;  // re-export for downstream users

// inside insert_size.rs
use crate::histogram::Histogram;   // `crate::` = this crate's root
```

**Editions** (`edition = "2021"`) are opt-in language epochs. Relevant gotcha you've already met: some recent dependency crates require `edition2024`, which needs rustc ≥ 1.85; on older compilers, pin those deps down.

---

## 2. Ownership — the core idea

```rust
let a = vec![300, 310, 320];   // `a` owns the Vec (heap buffer)
let b = a;                     // MOVE: ownership transfers to `b`; `a` is now invalid
// println!("{:?}", a);        // compile error: borrow of moved value
```

A `Vec` owns a heap allocation; copying the pointer and letting both `a` and `b` think they own it would mean a double-free at scope end. Rust forbids it: the move invalidates `a`. This is unlike Python (where `b = a` aliases the same list) and unlike a deep copy.

**`Copy` types are the exception.** Small, stack-only types (`i32`, `u8`, `f64`, `bool`, `char`, and tuples/arrays of them) implement the `Copy` trait, so `let y = x;` duplicates the bits and leaves `x` usable. That's why integer-heavy code (TLEN, MAPQ, insert sizes) rarely fights the borrow checker, while `Vec`/`String` code does.

```rust
let t: i32 = -300;
let u = t;          // copied; t still valid (i32 is Copy)
let s = String::from("chr1");
let r = s;          // moved; s invalid (String is not Copy)
```

**Three ways to give a function a value:**

```rust
fn takes_ownership(v: Vec<i64>) -> usize { v.len() }   // consumes — caller loses v
fn borrows(v: &[i64]) -> i64 { v.iter().copied().sum() } // reads — caller keeps v
fn mutates(v: &mut Vec<i64>) { v.push(0); }            // exclusive write borrow

let mut h = vec![300, 310, 310];
let s = borrows(&h);        // shared borrow, ends when the call returns
mutates(&mut h);            // exclusive borrow, allowed now that the shared one is gone
let n = takes_ownership(h); // h moved into the function; h unusable after this line
```

Default to **borrowing** (`&[T]`, `&str`) for function parameters; take ownership only when you genuinely need to store or consume the value. `clone()` makes an explicit deep copy when you really need an independent owner — it's never implicit, so a `.clone()` in a hot loop is a visible cost you chose.

---

## 3. Borrowing rules and the borrow checker

Two rules, enforced at compile time:

1. Any number of shared borrows `&T` **or** exactly one exclusive borrow `&mut T` — never both at once.
2. A borrow must not outlive the value it points to.

This is exactly why iterating a collection while mutating it is a *compile error* in Rust rather than a runtime surprise. The canonical mistakes and their real messages:

**Use after move (E0382):**

```rust
let h: Vec<i64> = vec![300, 310];
consume(h);                       // h moved
println!("{}", h.len());          // error
```
```
error[E0382]: borrow of moved value: `h`
  | move occurs because `h` has type `Vec<i64>`, which does not implement the `Copy` trait
```
Fix: borrow instead (`consume(&h)`), or `clone()` if you need two owners, or restructure so you don't use `h` after.

**Aliasing a mutable borrow (E0502):**

```rust
let mut h = vec![300, 310];
let first = &h[0];   // shared borrow
h.push(320);         // mutating while `first` is still live
println!("{}", first);
```
```
error[E0502]: cannot borrow `h` as mutable because it is also borrowed as immutable
  |  - immutable borrow occurs here
  |  ^^^^^^^^^^^ mutable borrow occurs here
  |  ----- immutable borrow later used here
```
Fix: finish using `first` before mutating, or copy the value out (`let first = h[0];` — `i64` is `Copy`) so no borrow is held. **Non-lexical lifetimes** mean a borrow ends at its last *use*, not at the end of the block — so often just reordering lines fixes it.

---

## 4. Lifetimes

A lifetime is the compiler's name for "how long this reference is valid." Most of the time they're inferred and invisible. You only write them when a function or struct *holds or returns* a reference and the compiler can't tell which input it's tied to.

```rust
// A zero-copy view that borrows a buffer it does not own.
pub struct ReadView<'a> { bases: &'a [u8] }

impl<'a> ReadView<'a> {
    pub fn new(bases: &'a [u8]) -> Self { Self { bases } }
    // the returned slice is valid as long as the borrowed buffer is
    pub fn first_n(&self, n: usize) -> &'a [u8] {
        &self.bases[..n.min(self.bases.len())]
    }
}

// "the result lives as long as both inputs"
pub fn longer<'a>(a: &'a [u8], b: &'a [u8]) -> &'a [u8] {
    if a.len() >= b.len() { a } else { b }
}
```

`'a` is just a label tying outputs to inputs so the compiler can prove the reference never dangles. This is the mechanism that lets `bammetrics` wrap a noodles record in a `ReadView` and hand out slices into its sequence/qualities **without copying** — exactly what you want when iterating tens of millions of reads. If lifetimes get hairy early on, the escape hatch is to own the data (`Vec<u8>`, `String`) instead of borrowing; optimize to borrows later.

---

## 5. The type system: structs, enums, pattern matching

**Structs** are product types (records).

```rust
#[derive(Debug, Clone, Default)]   // auto-generate these trait impls
pub struct InsertSizeRow {
    pub orientation: Orientation,
    pub read_pairs: i64,
    pub median_insert_size: f64,
}
```

**Enums** are sum types — far more powerful than C enums; each variant can carry data. This is the workhorse of idiomatic Rust.

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Orientation { #[default] Fr, Rf, Tandem }
```

**`Option<T>` and `Result<T, E>` are just enums** from the standard library:

```rust
enum Option<T> { Some(T), None }            // a value, or nothing (no null)
enum Result<T, E> { Ok(T), Err(E) }         // success or typed failure
```

There is no `null` and no untyped exception. Absence is `Option`; failure is `Result`. You handle them with `match` or the combinators (`map`, `unwrap_or`, `?`).

**`match` is exhaustive** — the compiler forces you to handle every variant, which is how adding an enum variant later surfaces every place that must change:

```rust
fn describe(o: Orientation) -> &'static str {
    match o {
        Orientation::Fr => "forward-reverse",
        Orientation::Rf => "reverse-forward",
        Orientation::Tandem => "same strand",
    }
}

// `if let` for the single-case shortcut:
if let Some(start) = record.alignment_start() { /* ... */ }
```

**The `Option<Result<T>>` pattern you'll hit constantly in noodles.** Accessors like `alignment_start()` return `Option<io::Result<Position>>` — "the field may be absent, and decoding it may fail." Two layers, peeled explicitly:

```rust
let rid: usize = match record.reference_sequence_id() {
    Some(Ok(id)) => id,
    Some(Err(e)) => return Err(e.into()),
    None => continue,            // unmapped: no reference id
};
// or terser, when you want to bail on error and skip on absence:
let rid = match record.reference_sequence_id() { Some(x) => x?, None => continue };
```

That `x?` propagates the inner error; the outer `Some/None` you decide per the domain (here, skip unmapped reads).

---

## 6. Error handling

No exceptions. Fallible functions return `Result<T, E>`. The `?` operator is the ergonomic core: on `Ok(v)` it unwraps to `v`; on `Err(e)` it returns early, converting `e` into the function's error type via the `From` trait.

```rust
use std::fs::File;
use std::io::{self, Read};
use std::num::ParseIntError;

#[derive(Debug)]
pub enum MetricsError {
    Io(io::Error),
    Parse(ParseIntError),
    BadRecord { name: String, reason: String },
}
impl From<io::Error> for MetricsError {        // enables `?` on io::Result
    fn from(e: io::Error) -> Self { MetricsError::Io(e) }
}
impl From<ParseIntError> for MetricsError {    // enables `?` on parse()
    fn from(e: ParseIntError) -> Self { MetricsError::Parse(e) }
}

pub fn read_tlen(path: &str) -> Result<i64, MetricsError> {
    let mut s = String::new();
    File::open(path)?.read_to_string(&mut s)?;   // io::Error auto-converted
    let n: i64 = s.trim().parse()?;              // ParseIntError auto-converted
    if n == 0 {
        return Err(MetricsError::BadRecord { name: path.into(), reason: "TLEN zero".into() });
    }
    Ok(n.abs())
}
```

Writing `From` by hand is boilerplate. **`thiserror`** (for libraries) derives `Display` and `From` for you:

```rust
#[derive(Debug, thiserror::Error)]
pub enum BamError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),                        // #[from] generates the From impl
    #[error("record {name} invalid: {reason}")]
    BadRecord { name: String, reason: String },
}
```

Rule of thumb: **`thiserror` for libraries** (callers want to match on specific error variants), **`anyhow` for applications/CLIs** (you just want `Result<()>` and a good error message; `anyhow::Result<T>` boxes any error). `bammetrics`-core is a library → `thiserror`; a CLI on top → `anyhow`.

Avoid `.unwrap()`/`.expect()` outside tests and truly-impossible cases — they panic (abort the thread). `panic!` is for bugs, not for expected failure; expected failure is `Result`.

---

## 7. Traits and generics

A `trait` is an interface — closer to a Swift protocol or a Haskell typeclass than a TS interface, because it's nominal (you explicitly `impl Trait for Type`) and can carry default methods.

```rust
pub trait Metric {
    fn name(&self) -> &str;
    fn value(&self) -> f64;
    fn formatted(&self) -> String {                 // default method
        format!("{}={:.4}", self.name(), self.value())
    }
}

pub struct Median(pub f64);
impl Metric for Median {
    fn name(&self) -> &str { "MEDIAN_INSERT_SIZE" }
    fn value(&self) -> f64 { self.0 }
    // formatted() inherited from the default
}
```

**Two dispatch styles:**

```rust
// Static dispatch: generic, monomorphized at compile time — zero runtime cost.
pub fn print_static<M: Metric>(m: &M) -> String { m.formatted() }
// equivalent sugar:
pub fn print_static2(m: &impl Metric) -> String { m.formatted() }

// Dynamic dispatch: one vtable, lets you store mixed types together.
pub fn print_all(ms: &[Box<dyn Metric>]) -> Vec<String> {
    ms.iter().map(|m| m.formatted()).collect()
}
```

Use generics (`<M: Metric>`) by default — no overhead. Reach for `dyn Trait` (a "trait object") only when you need a heterogeneous collection, like `Vec<Box<dyn Metric>>` holding several different metric types.

**The standard traits you'll use everywhere**, almost always via `#[derive(...)]`:

| Trait | What it gives you | Derive? |
|---|---|---|
| `Debug` | `{:?}` formatting | yes |
| `Clone` | explicit deep copy `.clone()` | yes |
| `Copy` | implicit bit-copy (small types) | yes (with Clone) |
| `Default` | `T::default()` | yes |
| `PartialEq`/`Eq` | `==` | yes |
| `PartialOrd`/`Ord` | `<`, sorting, `BTreeMap` keys | yes |
| `Hash` | `HashMap`/`HashSet` keys | yes |
| `From`/`Into` | conversions; powers `?` | by hand or `thiserror` |
| `Iterator` | `for` loops, adapters | by hand (see §8) |

`From<A> for B` automatically gives you `Into<B> for A` — implement `From`, get `.into()` free.

---

## 8. Iterators and closures — the workhorse

Rust iterators are **lazy** (nothing happens until consumed) and **zero-cost** (they compile down to the same machine code as a hand-written loop). This is where idiomatic Rust reads like a clean data pipeline:

```rust
pub fn mean_insert_size(tlens: &[i32]) -> Option<f64> {
    let (sum, count) = tlens.iter()
        .filter(|&&t| t != 0)                  // skip TLEN==0
        .map(|&t| t.unsigned_abs() as u64)     // abs, widen
        .fold((0u64, 0u64), |(s, c), x| (s + x, c + 1));
    if count == 0 { None } else { Some(sum as f64 / count as f64) }
}
```

`|&&t| t != 0` is a closure; the `&&` destructures the double reference that `.iter()` over a slice produces (`&i32`, and the filter closure gets `&&i32`). Closures capture their environment like in Swift/JS; `move |x| ...` forces capture-by-value (needed when handing a closure to another thread).

Key adapters: `map`, `filter`, `filter_map`, `flat_map`, `take`, `skip`, `enumerate`, `zip`, `chain`, `scan`, `inspect`. Consumers (force evaluation): `collect`, `sum`, `count`, `fold`, `for_each`, `any`, `all`, `find`, `max`, `min`, `position`. `collect()` is type-directed — annotate the target:

```rust
let v: Vec<u64> = it.collect();
let m: std::collections::BTreeMap<i64, f64> = pairs.collect();
```

**Writing your own iterator** is one trait, one method:

```rust
pub struct RunningMax<I: Iterator<Item = i64>> { inner: I, cur: i64 }
impl<I: Iterator<Item = i64>> Iterator for RunningMax<I> {
    type Item = i64;
    fn next(&mut self) -> Option<i64> {
        let x = self.inner.next()?;          // None ends the iteration
        self.cur = self.cur.max(x);
        Some(self.cur)
    }
}
```

Once you `impl Iterator`, you get all the adapters for free. noodles' `reader.records()` is exactly this pattern — an iterator of `io::Result<Record>` you can `filter`/`map` over.

---

## 9. Collections — pick the right one

- **`Vec<T>`** — growable array; your default sequence. `with_capacity(n)` preallocates when you know the size (avoids reallocation in hot paths).
- **`BTreeMap<K, V>`** — sorted map. **This is what the Picard histogram must use**: medians, modes, and min/max all depend on ascending key iteration order. `HashMap` would give wrong answers there.
- **`HashMap<K, V>`** — unordered, faster lookups; for read-group → accumulator tables where order doesn't matter.
- **`HashSet`/`BTreeSet`** — membership.
- **`VecDeque`** — ring buffer (e.g., a sliding window over positions).

```rust
use std::collections::BTreeMap;
let mut hist: BTreeMap<i64, f64> = BTreeMap::new();
*hist.entry(insert_size).or_insert(0.0) += 1.0;   // the increment idiom
for (&size, &count) in &hist { /* ascending by size */ }
```

The `entry(k).or_insert(default)` pattern is the canonical "get-or-create-then-update" — no double lookup.

---

## 10. Strings and bytes (matters for sequence data)

- **`String`** owns a heap UTF-8 buffer; **`&str`** borrows one. Take `&str` as parameters, return `String` when you produce text. They are **guaranteed valid UTF-8**.
- **Sequence bases and qualities are not text — they're raw bytes.** Use `Vec<u8>` / `&[u8]`, not `String`. A BAM read's `A/C/G/T/N` are ASCII byte values; comparing read base to reference base is a `u8 == u8` after case-folding:

```rust
fn bases_equal(a: u8, b: u8) -> bool { a.to_ascii_uppercase() == b.to_ascii_uppercase() }
```

This is why noodles' `sequence().iter()` yields `u8` and `quality_scores().as_ref()` is `&[u8]`. Don't reach for `String` here; you'd pay UTF-8 validation for no reason and complicate indexing.

---

## 11. Concurrency and parallelism

Rust's ownership rules make data races a *compile error*: a value can't be shared mutably across threads unless the types opt in via the `Send` (movable to another thread) and `Sync` (shareable by reference) marker traits, which the compiler tracks for you.

For data-parallel work — the natural shape of "process millions of reads" — **`rayon`** turns a sequential iterator into a parallel one by changing `iter()` to `par_iter()`:

```rust
use rayon::prelude::*;

pub fn total_aligned_bases(per_read: &[u64]) -> u64 {
    per_read.par_iter().copied().sum()              // splits across all cores
}
pub fn count_high_mapq(mapqs: &[u8]) -> usize {
    mapqs.par_iter().filter(|&&q| q >= 20).count()
}
```

For BAM specifically, the clean parallelization is **per-contig or per-chunk**: decode regions independently, accumulate a `Histogram` per chunk, then merge the histograms at the end (addition is associative, so the merge is trivial and order-independent). That keeps the parallel section lock-free. Pin `rayon-core` to a version your compiler supports if you're on older rustc.

When you do need shared mutable state: `Arc<T>` (atomic reference count, shared ownership across threads) wrapping `Mutex<T>`/`RwLock<T>`. Reach for these last — message-passing (`std::sync::mpsc` channels) or the merge-at-the-end pattern is usually cleaner and faster.

---

## 12. Performance idioms

Rust gives you C-level performance *if* you don't accidentally allocate. In a loop over tens of millions of records:

- **Borrow, don't clone.** `&[u8]` over the record's buffer beats `to_vec()`. Use the `ReadView<'a>` lifetime pattern (§4) to hand out slices without copying.
- **Iterators over intermediate `Vec`s.** `xs.iter().filter(...).sum()` allocates nothing; `xs.iter().filter(...).collect::<Vec<_>>().iter().sum()` allocates a throwaway vector.
- **Preallocate.** `Vec::with_capacity(n)`, `String::with_capacity(n)` when the size is known.
- **Avoid `format!`/`String` in hot loops.** Build output once, at the end.
- **`#[inline]`** on tiny hot functions if profiling shows call overhead (usually unnecessary — LLVM inlines aggressively in release).
- **Always profile in `--release`.** Debug builds are 10–50× slower and will mislead you. Use `criterion` for stable microbenchmarks and `perf`/`cargo flamegraph` for hotspots.
- **`u32`/`u64` choice matters** for cache footprint when you have billions of values; pick the smallest correct integer.

The discipline is the same one you apply in production pipelines: measure first, and treat every allocation in the inner loop as suspect.

---

## 13. Testing and benchmarking

Tests live next to the code; `#[cfg(test)]` keeps them out of release builds.

```rust
pub fn median_of_two(a: i64, b: i64) -> f64 { (a + b) as f64 / 2.0 }

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn median_even() {
        assert_eq!(median_of_two(300, 320), 310.0);
    }
    #[test]
    fn floats_close() {
        let got = 8.16496_f64;
        assert!((got - 8.16497).abs() < 1e-3);   // float comparisons need tolerance
    }
}
```

`cargo test` runs them. **Integration tests** go in a top-level `tests/` directory and test your crate as an external user would (this is where the BAM-concordance harness lives). `assert_eq!`, `assert!`, `#[should_panic]`, and `-- --nocapture` to see `println!` output. For benchmarks use **`criterion`** (stable, statistically rigorous) rather than the nightly-only built-in `bench`.

For `bammetrics`, the highest-value tests are the hand-computed histogram cases (even/odd counts, half-integer median for the MAD edge) — small, deterministic, and they lock the concordance-critical math.

---

## 14. Decoding the compiler (a field guide)

Rust errors are verbose but precise. The frequent ones, with the fix:

- **E0382 "borrow of moved value"** — you used a value after moving it. Fix: borrow (`&x`), `clone()`, or reorder so the last use precedes the move. (Real message in §3.)
- **E0502 "cannot borrow as mutable because also borrowed as immutable"** — an outstanding shared borrow overlaps a mutation. Fix: shorten the shared borrow (copy the value out if `Copy`), or reorder so it ends first. (Real message in §3.)
- **E0499 "cannot borrow as mutable more than once"** — two `&mut` to the same value alive together. Fix: scope them apart, or use indices instead of holding references.
- **"does not live long enough"** — a reference outlives its referent. Fix: own the data, or tie the lifetimes with `'a`, or extend the owner's scope.
- **"the trait bound `T: Trait` is not satisfied"** — you called a method/used a feature the type doesn't implement. Fix: add the `#[derive(...)]` or `impl`, or add the bound to your generic (`<T: Trait>`).
- **"mismatched types: expected `&str`, found `String`"** (and friends) — usually fixed with `&`, `.as_str()`, `.to_string()`, `.into()`, or `.as_slice()`.
- **"cannot move out of borrowed content"** — you tried to take ownership through a `&`. Fix: `clone()`, or work with the reference, or `std::mem::take`/`replace`.

Run `rustc --explain E0502` for a full worked explanation of any code. `cargo clippy` is the linter — run it; it teaches idiom and catches real bugs (needless clones, `.unwrap()` smells, inefficient patterns).

---

## 15. A learning path

In order of leverage for actually shipping `bammetrics`:

1. **Internalize move vs borrow vs clone** (§2–3). Write small functions that take `&[T]`, `&mut Vec<T>`, and `T`; feel when each is needed. This removes 80% of early friction.
2. **`Option`/`Result`/`?` and `match`** (§5–6). Build the noodles `Option<Result<T>>` unwrap reflex.
3. **Iterators** (§8). Rewrite a manual loop as an iterator pipeline; write one `impl Iterator`.
4. **Traits + derive** (§7). Define `Metric`; derive `Debug/Clone/Default` on your structs.
5. **Lifetimes only when needed** (§4). Don't pre-optimize to borrows; own first, borrow when profiling demands it.
6. **rayon + chunked merge** (§11) once the sequential version is correct and tested.

Resources worth your time: *The Rust Programming Language* ("the book", free online) for depth, *Rust by Example* for runnable snippets, `cargo clippy` as a continuous tutor, and the std docs (`cargo doc --open`) which are excellent. Your Python/SQL/Swift instincts transfer; the borrow checker is the one genuinely new muscle, and it's built by losing arguments to it for a week.

---

*Non-trivial examples (errors, iterators, traits, lifetimes, rayon, ownership) were compiled on rustc 1.75.0; the E0382/E0502 messages in §3 are real compiler output. The iterator/trait/error patterns are stable across modern Rust; only the dependency-version pins (rayon-core, hashbrown, indexmap for edition2024) are toolchain-sensitive.*
